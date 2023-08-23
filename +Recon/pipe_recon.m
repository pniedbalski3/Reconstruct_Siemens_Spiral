function Image_Out = pipe_recon(ImageSize,data,traj,DCF,cur_it,tot_coils,PixelShift)
%% A Function written to reconstruct Images when K-space data and trajectories are passed to it
% Uses Pipe's Group DCF and gridding code. This is for 3D data
% 
% ImageSize - Scalar: Output image matrix size
%
% data - KSpace Data in (N_ro x N_Proj) or (N_ro x N_Proj x N_Chan)
%
% traj - point in kspace corresponding to the data vector - columns for
% x,y, and z. (3 x N_ro x N_Proj)
%
% PixelShift - array of size 3 with pixels to shift by

%% Settings
numIter = 5;
effMtx  = ImageSize;
alpha = 2; %Image overgrid factor
osf = 2.1; %grid oversample factor on top of effMtx, at least 2, 2.1 optimal, >2.6 conservative
verbose = 0; %c verbose is not real time so not very useful
numThread = 1; %unthreaded since not Posix
if exist('PixelShift','var')==0%if not passed, set to 0's
    PixelShift = [0, 0, 0];
end
if exist('tot_coils','var')==0%if not passed, set to 0's
    tot_coils = size(data,3);
end
if exist('cur_it','var')==0%if not passed, set to 0's
    cur_it = nan;
end

%% Prep inputs
data_real = reshape(real(data), [1 size(data)]);
data_imag = reshape(imag(data), [1 size(data)]);
data_comb = cat(1, data_real, data_imag);
data_comb = double(data_comb);
traj = double(traj);

%Free up memory
clear data_real;
clear data_imag;
clear data;
chan = size(data_comb,4);

% %% DCF
% disp('Calculatuing DCF...');
% DCF = recon.DC.sdc3_MAT(traj,numIter,effMtx,verbose,osf);
% DCF = double(DCF);

%% Grid
gdata = zeros(2,effMtx*alpha,effMtx*alpha,effMtx*alpha,size(data_comb,4));
for coil = 1:chan
    if isnan(cur_it)
        n = coil;
    else
        n = cur_it;
    end
    disp(['Gridding channel ', num2str(n), ' of ', num2str(tot_coils), '...']);
    gdata(:,:,:,:,coil) = Recon.Grid.grid3_MAT(squeeze(data_comb(:,:,:,coil)),traj,DCF,effMtx*alpha,numThread);
end
%Free memory now that done with data_comb and traj
clear data_comb;
clear traj;
%%
gdata = squeeze(gdata(1,:,:,:,:) + 1j*gdata(2,:,:,:,:));

%% FFT to image
disp('FFTing gridded data...');
Image = fft(gdata,[],1);
Image = fft(Image,[],2);
Image = fft(Image,[],3);
Image = fftshift(Image,1);
Image = fftshift(Image,2);
Image = fftshift(Image,3);

%% Calculated rolloff kernel
disp('Calculating rolloff kernel...');
delta = [1.0, 0.0];
k_not = [0.0, 0.0, 0.0];
DCF_not = 1.0;
rokern = Recon.Grid.grid3_MAT(delta',k_not',DCF_not,effMtx*alpha,numThread);

%% FFT to rolloff image
disp('FFTing rolloff kernel to image space...');
rokern = squeeze(rokern(1,:,:,:) + 1j*rokern(2,:,:,:));
rokern = fftn(rokern);
rokern = fftshift(rokern,1);
rokern = fftshift(rokern,2);
rokern = fftshift(rokern,3);

%% Apply rolloff kernel image
for coil = 1:chan
    if isnan(cur_it)
        n = coil;
    else
        n = cur_it;
    end
    disp(['Applying rolloff for channel ', num2str(n), ' of ', num2str(tot_coils), '...']);
    Image(:,:,:,coil) = Image(:,:,:,coil) ./ rokern;
end

%% Shift pixels
disp('Applying shift and cropping image...');
Image = circshift(circshift(circshift(Image,round(PixelShift(1)),1),round(PixelShift(2)),2),round(PixelShift(3)),3);

%% Crop
xs = floor(effMtx - effMtx/alpha)+1;
xe = floor(effMtx + effMtx/alpha);
Image_Out = Image(xs:xe,xs:xe,xs:xe,:);

disp('Reconstruction complete.');

end





