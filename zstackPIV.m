function [u_filt, v_filt, x,y] =zstackPIV(imagesbefore,imagesafter,parpiv,par2dpostprocess)


xsize = length(imagesbefore(:,1,1));
ysize = length(imagesbefore(1,:,1));
zsize = length(imagesbefore(1,1,:));


amount = zsize;

%% Standard PIV Settings
s = cell(10,2); % To make it more readable, let's create a "settings table"
%Parameter                       %Setting           %Options
s{1,1}= 'Int. area 1';           s{1,2}=32;         % window size of first pass
s{2,1}= 'Step size 1';           s{2,2}=16;         % step of first pass
s{3,1}= 'Subpix. finder';        s{3,2}=1;          % 1 = 3point Gauss, 2 = 2D Gauss
s{4,1}= 'Mask';                  s{4,2}=[];         % If needed, generate via: imagesc(image); [temp,Mask{1,1},Mask{1,2}]=roipoly;
s{5,1}= 'ROI';                   s{5,2}=[];         % Region of interest: [x,y,width,height] in pixels, may be left empty
s{6,1}= 'Nr. of passes';         s{6,2}=1;          % 1-4 nr. of passes
s{7,1}= 'Int. area 2';           s{7,2}=16;         % second pass window size
s{8,1}= 'Int. area 3';           s{8,2}=8;         % third pass window size
s{9,1}= 'Int. area 4';           s{9,2}=16;         % fourth pass window size
s{10,1}='Window deformation';    s{10,2}='*linear'; % '*spline' is more accurate, but slower
if(nargin>2)
    s = parpiv;
end

counter=0;
%% PIV analysis loop:


x=cell(amount,1);
y=x;
u=x;
v=x;
typevector=x; %typevector will be 1 for regular vectors, 0 for masked areas
for i=1:amount       % Here "i" is defined as a zstack number
    counter=counter+1;
    image1= single(squeeze(imagesbefore(:,:,i)));
    image2= single(squeeze(imagesafter(:,:,i)));
    
    [x{counter} y{counter} u{counter} v{counter} typevector{counter}] = piv_FFTmulti (image1,image2,s{1,2},s{2,2},s{3,2},s{4,2},s{5,2},s{6,2},s{7,2},s{8,2},s{9,2},s{10,2});
    
    sprintf('2D PIV:frame %d out of %d', i, amount)
end
%% PIV postprocessing loop
% Settings
if(nargin>3)
    umin = par2dpostprocess.umin;
    umax = par2dpostprocess.umax;
    vmin = par2dpostprocess.vmin;
    vmax = par2dpostprocess.vmax;
    stdthresh = par2dpostprocess.stdthresh;
    epsilon = par2dpostprocess.epsilon;
    thresh = par2dpostprocess.thresh;
    umargin = par2dpostprocess.umargin;
    vmargin = par2dpostprocess.vmargin;
else
    umin = -25; % minimum allowed u velocity
    umax = 25; % maximum allowed u velocity
    vmin = -25; % minimum allowed v velocity
    vmax = 25; % maximum allowed v velocity
    stdthresh=4; % threshold for standard deviation check
    epsilon=0.1; % epsilon for normalized median test
    thresh=0.5; % threshold for normalized median test
    umargin = []; %a corner patch to calculate shift in u;
    vmargin = []; %a corner patch to calculate shift in v;
    
end


u_filt=cell(amount,1);
v_filt=u_filt;
typevector_filt=u_filt;
for PIVresult=1:size(x,1)
    if(isempty(umargin))
        shiftx = 0;
        shifty = 0;
    else
        shiftx = nanmean(nanmean(u{PIVresult,1}(umargin(1):umargin(2),...
            umargin(3):umargin(4))));
        shifty = nanmean(nanmean(v{PIVresult,1}(vmargin(1):vmargin(2),...
            vmargin(3):vmargin(4))));
    end
    u_filtered=u{PIVresult,1}-shiftx;
    v_filtered=v{PIVresult,1}-shifty;
    
    typevector_filtered=typevector{PIVresult,1};
    %vellimit check
    u_filtered(u_filtered<umin)=NaN;
    u_filtered(u_filtered>umax)=NaN;
    v_filtered(v_filtered<vmin)=NaN;
    v_filtered(v_filtered>vmax)=NaN;
    if(isnan(max(abs(u_filtered(:)+v_filtered(:)))))
        sprintf('velolcity in frame %d may get upper/lower limit correction',PIVresult)
    end
    % stddev check
    meanu=nanmean(nanmean(u_filtered));
    meanv=nanmean(nanmean(v_filtered));
    std2u=nanstd(reshape(u_filtered,size(u_filtered,1)*size(u_filtered,2),1));
    std2v=nanstd(reshape(v_filtered,size(v_filtered,1)*size(v_filtered,2),1));
    minvalu=meanu-stdthresh*std2u;
    maxvalu=meanu+stdthresh*std2u;
    minvalv=meanv-stdthresh*std2v;
    maxvalv=meanv+stdthresh*std2v;
    u_filtered(u_filtered<minvalu)=NaN;
    u_filtered(u_filtered>maxvalu)=NaN;
    v_filtered(v_filtered<minvalv)=NaN;
    v_filtered(v_filtered>maxvalv)=NaN;
    if(isnan(max(abs(u_filtered(:)+v_filtered(:)))))
        sprintf('velolcity in frame %d may get stddev correction',PIVresult)
    end
    
    % normalized median check
    %Westerweel & Scarano (2005): Universal Outlier detection for PIV data
    [J,I]=size(u_filtered);
    medianres=zeros(J,I);
    normfluct=zeros(J,I,2);
    b=1;
    for c=1:2
        if c==1; velcomp=u_filtered;else;velcomp=v_filtered;end %#ok<*NOSEM>
        for k=1+b:I-b
            for j=1+b:J-b
                neigh=velcomp(j-b:j+b,k-b:k+b);
                neighcol=neigh(:);
                neighcol2=[neighcol(1:(2*b+1)*b+b);neighcol((2*b+1)*b+b+2:end)];
                med=median(neighcol2);
                fluct=velcomp(j,k)-med;
                res=neighcol2-med;
                medianres=median(abs(res));
                normfluct(j,k,c)=abs(fluct/(medianres+epsilon));
            end
        end
    end
    info1=(sqrt(normfluct(:,:,1).^2+normfluct(:,:,2).^2)>thresh);
    u_filtered(info1==1)=NaN;
    v_filtered(info1==1)=NaN;
    if(isnan(max(abs(u_filtered(:)+v_filtered(:)))))
        sprintf('velolcity in frame %d may get median correction',PIVresult)
    end
    
    typevector_filtered(isnan(u_filtered))=2;
    typevector_filtered(isnan(v_filtered))=2;
    typevector_filtered(typevector{PIVresult,1}==0)=0; %restores typevector for mask
    
    %Interpolate missing data
    u_filtered=inpaint_nans(u_filtered,4);
    v_filtered=inpaint_nans(v_filtered,4);
    if(isempty(umargin))
        shiftx = 0;
        shifty = 0;
    else
        shiftx = nanmean(nanmean(u_filtered(umargin(1):umargin(2),...
            umargin(3):umargin(4))));
        shifty = nanmean(nanmean(v_filtered(vmargin(1):vmargin(2),...
            vmargin(3):vmargin(4))));
    end
    u_filtered = u_filtered+shiftx;
    v_filtered = v_filtered+shifty;
        
    u_filt{PIVresult,1}=smoothn(u_filtered,'robust');
    v_filt{PIVresult,1}=smoothn(v_filtered,'robust');
                
    typevector_filt{PIVresult,1}=typevector_filtered;
end



end