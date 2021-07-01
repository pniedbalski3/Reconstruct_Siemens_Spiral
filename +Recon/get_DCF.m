function DCF = get_DCF(traj,ImSize,nIter)

if nargin < 3
    nIter = 10;
end

%Default to no verbose:
verbose = 0;
osf = 2.1; %evidently this is the optimal?

DCF = Recon.DC.sdc3_MAT(traj,nIter,ImSize,verbose,osf);
DCF = double(DCF);