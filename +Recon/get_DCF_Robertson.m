function DCF = get_DCF_Robertson(traj,ImSize,nIter)

%Do lots of things in this function - basically up to where we need the
%actual data - The DCF being passed back is not strictly the Density
%Compensation Function, but rather the recon object used by this
%reconstruction. traj need to be a column vector

verbose = false;%true();

%Some decent values for reconstruction:
overgrid_factor = 2;
kernel.sharpness = .3;
kernel.extent = 9*kernel.sharpness;
cropOvergriddedImage = true;
deapodizeImage = false;

%Set up kernel
kernelObj = Recon.SysModel.Kernel.Gaussian(kernel.sharpness, kernel.extent, verbose);

%set up L2 proximity
proxObj = Recon.SysModel.Proximity.L2Proximity(kernelObj, verbose);

%set up system
systemObj = Recon.SysModel.MatrixSystemModel(traj, overgrid_factor, ...
    ImSize, proxObj, verbose);

dcfObj = Recon.DCF.Iterative(systemObj, nIter, verbose);

reconObj = Recon.ReconModel.LSQGridded(systemObj, dcfObj, verbose);
clear modelObj;
clear dcfObj;
reconObj.crop = cropOvergriddedImage;
reconObj.deapodize = deapodizeImage;

DCF = reconObj;