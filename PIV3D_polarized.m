function [udef,vdef,wdef,wdefraw,xdef,ydef,zdef,badfit,outliers,direcvolcorr,pargrid] = PIV3D_polarized(imagesbefore,imagesafter,pars)
%The deformation field tells how to deform imagesbefore into
%imagesafter. The coordinates are based on imagesbefore.

%This 3D PIV assumes 2D deformation can be obtained slice by
%slice first, and z deformation is a small perturbation. This is the case
%for confocal reflection imaging. 

%The 2D PIV requires PIV_Lab has been installed, and is in the search
%directory of MATLAB. PIV_Lab can be found on MATLAB central. However, you
%may replace the PIV_Lab by other algorithms, simply rewrite the container
%function zstackPIV().

%output: 

%deformation field is described by[udef,vdef,wdef,xdef,ydef,zdef]. Notice
%that here xdef and udef
%correspond to the second dimension of imagesbefore, or imagesafter, as a
%funny convention of MATLAB.

%badfit is a matrix the same size of udef,vdef..... If badfit is 1 or -1 that
%means the corresponding zdisplacement wdef was obtained from wdefraw by
%nearest-neighbor interpolations.

%direcvolcorr is a cell array which is the correlation curve of each
%subvolume.

%pargrid collects the input parameters.

%input:

%imagesbefore, imagesafter are both 3D matrices.

%pars: a structure containing the parameters of PVC.
%see the comments below for interpretations.


%Bo Sun March 19th 2015
imagesbefore = single(imagesbefore);
imagesafter = single(imagesafter);

%similar to window size in 2D PIV
if(isfield(pars,'zdepth'))
    zdepth = pars.zdepth;
else
    zdepth = 4;
end

%similar to step size in 2D PIV
if(isfield(pars,'zstep'))
    zstep = pars.zstep;
else
    zstep = 2;
end

%maximum displacement in z direction
if(isfield(pars,'maxzdis'))
    maxzdis = pars.maxzdis;
else
    maxzdis= 10;
end

%The length of curve to fit for the z displacement. Typically no need to
%change.
if(isfield(pars,'fitwidth'))
    fitwidth = pars.fitwidth;
else
    fitwidth = 3; %this gives a 7 point gaussian fit
end

%The radius of local standard deviation thresholding. 
if(isfield(pars,'radius'))
   radius  = pars.radius;
else
   radius  = 3; 
end

%The threshold of local standard deviation thresholding. 
if(isfield(pars,'threshold'))
   threshold  = pars.threshold;
else
   threshold  = 4; 
end


%First do 2D PIV

%par2dpostprocess is the parameter for postprocessing 2D PIV
par2dpostprocess.umin = -20;
par2dpostprocess.umax = 20;
par2dpostprocess.vmin = -20;
par2dpostprocess.vmax = 20;
par2dpostprocess.stdthresh = 10;
par2dpostprocess.epsilon = 0.5;
par2dpostprocess.thresh = 5;
par2dpostprocess.umargin = [];
par2dpostprocess.vmargin = [];


if(isfield(pars,'par2dpostprocess'))
    par2dpostprocess = pars.par2dpostprocess;
end

%parpiv is the parameter for PIV
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
parpiv = s;
if(isfield(pars,'parpiv'))
    parpiv = pars.parpiv;
end

%u_xy,v_xy will tell how to deform imagesbefore into images after
[u_xy, v_xy, x_xy,y_xy] = zstackPIV(imagesbefore,imagesafter,parpiv,par2dpostprocess);
if(parpiv{6,2} == 1)
    xywindowsz = parpiv{1,2};
else
    xywindowsz = parpiv{s{6,2}+5,2};
end

xystep = x_xy{1}(1,2)-x_xy{1}(1,1);


maxx = x_xy{1}(1,end);
maxy = y_xy{1}(end,1);



pargrid.xywindowsz = xywindowsz;
pargrid.xystep = xystep;
pargrid.maxx = maxx;
pargrid.maxy = maxy;
pargrid.zdepth = zdepth;
pargrid.zstep = zstep;


indz = 0;
direcvolcorr = {0};
udef = zeros(length(xywindowsz/2+1:xystep:maxy),...
    length(xywindowsz/2+1:xystep:maxx),...
    length(1:zstep:(length(u_xy)-zdepth)));
vdef = udef;
wdef = udef;
xdef = udef;
ydef = udef;
zdef = udef;
badfit = udef;
for zz = 1:zstep:(length(u_xy)-zdepth)
    indx = 0;
    indz = indz+1;
    %take average xy motion in the depth
    uu = 0;
    vv = 0;
    for zt = 0:zdepth-1
        uu = uu + u_xy{zz+zt};
        vv = vv + v_xy{zz+zt};
    end
    uavg = uu/zdepth;
    vavg = vv/zdepth;
    uu = round(uu/zdepth);
    vv = round(vv/zdepth);
    
    for xx =  xywindowsz/2+1:xystep:maxx
        indx = indx + 1;
        indy = 0;
        for yy = xywindowsz/2+1:xystep:maxy
            indy = indy + 1;
            ystart = yy-xywindowsz/2;
            yend = yy+xywindowsz/2-1;
            xstart = xx-xywindowsz/2;
            xend = xx+xywindowsz/2-1;
            zstart = zz;
            zend = zz+zdepth-1;
            subvolbefore = imagesbefore(yend:-1:ystart,xend:-1:xstart,zend:-1:zstart);
            subvolbefore = subvolbefore-mean(subvolbefore(:));
            ystart2 = ystart+vv(indy,indx);
            yend2 = yend+vv(indy,indx);
            xstart2 = xstart+uu(indy,indx);
            xend2 = xend+uu(indy,indx);
            zstart2 = zz-maxzdis;
            zend2 = zz+zdepth-1+maxzdis;
            subvolafter = matrixnolimit(imagesafter,[[ystart2,yend2];[xstart2,xend2];[zstart2,zend2]],0);
            similarity = squeeze(convn(subvolafter,subvolbefore,'valid'));
            direcvolcorr{indy,indx,indz} = similarity;
            udef(indy,indx,indz) = uavg(indy,indx);
            vdef(indy,indx,indz) = vavg(indy,indx);
            xdef(indy,indx,indz) = xx;
            ydef(indy,indx,indz) = yy;
            zdef(indy,indx,indz) = (zstart+zend)/2;
            [tempm,tempc] = max(similarity);
            simind1 = max([1,tempc-fitwidth]);
            simind2 = min([length(similarity),tempc+fitwidth]);
            
            if(tempc == 1)
                wdef(indy,indx,indz) = -maxzdis;
                badfit(indy,indx,indz) = -1;
            elseif( tempc == length(similarity))
                wdef(indy,indx,indz) = length(similarity)-maxzdis-1;
                badfit(indy,indx,indz) = 1;
            else
                simtofit = similarity(simind1:simind2);
                simtofit = simtofit - min(simtofit);
                fitx = simind1:simind2;
            try %fit may cause error
                f = fit(fitx',simtofit,'gauss1');
            catch %first solution is to reduce fitwidth to 1
                display(sprintf('fit with shorter width at x-y-z index %d, %d, %d',...
                          indx,indy,indz));
                 simind1 = max([1,tempc-1]);
                 simind2 = min([length(similarity),tempc+1]);  
           
                 simtofit = similarity(simind1:simind2);
                 simtofit = simtofit - min(simtofit);
                 fitx = simind1:simind2;
                 try
                     f = fit(fitx',simtofit,'gauss1');
                 catch
                      f.b1 = NaN;
                      display(sprintf('fit failed at x-y-z index %d, %d, %d',...
                          indx,indy,indz));
                 end
            end                                                                                
                wdef(indy,indx,indz) = f.b1-maxzdis-1;
            end
            
        end
    end
    display(sprintf('finished z substack from z = %d to %d',zstart, zend));
end


wdefraw = wdef;
%Now interpolate the missing data and also take away the outliers
display('filter out outliers in z displacement');
wdef(abs(badfit)>0) = NaN;
wdef(abs(wdef)>maxzdis) = NaN;
[data,outliers] = outlierfilter3d_mstd(wdef,radius,threshold);
wdef = data;
wdef = inpaintn(wdef,400);

end
