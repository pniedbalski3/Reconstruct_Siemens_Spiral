function reconstruct_spiral(file)

%Function to Reconstruct Spiral Data from SPIRAL sequence created by Peter
%Niedbalski for Siemens Skyra scanner.

%% Get some path information - where to save images, where to pull some ancillary info from
idcs = strfind(file,filesep);%determine location of file separators
write_path = file(1:idcs(end)-1);%remove file

parent_path = which('reconstruct_spiral');
idcs = strfind(parent_path,filesep);%determine location of file separators
parent_path = parent_path(1:idcs(end)-1);%remove file


%% Read in Data

disp('Reading twix file')
twix_obj = Import.mapVBVD(file,'ignoreSeg');
%For now, throw out coil sensitivity information... may want to fix in the
%future
if length(twix_obj) > 1
    twix_obj = twix_obj{end};
end
%When I coded this sequence, I made it possible to have more lines of
%acquisition than is planned for in the underlying data architecture. The
%data is still stored, just a little harder to get to - check for that
%here.
%I think this should correctly identify this situation 
if twix_obj.image.NLin ~= length(twix_obj.image.Lin) && mod(length(twix_obj.image.Lin),twix_obj.image.NLin) ~=0
    disp('Data Overload, perorming alternate twix read')
    twix_obj = Import.mapVBVD_pjn(file,'ignoreSeg');
end

raw = squeeze(double(twix_obj.image()));

slices = twix_obj.hdr.MeasYaps.sSliceArray.lSize;
%Handle 2D vs 3D data
if slices > 1
    ImSize = [twix_obj.hdr.MeasYaps.sKSpace.lBaseResolution twix_obj.hdr.MeasYaps.sKSpace.lBaseResolution 1]; %Desired Output Image Size
    Dim = 0; %Dimension - 0 means 2D
else
    ImSize = twix_obj.hdr.MeasYaps.sKSpace.lBaseResolution; %Desired Output Image Size
    Dim = 1; %Dimension - 1 means 3D
end

%Need to reshape data intelligently:
%Order of loops is points, coils, projections, slices, bvalues (I think...)
%Basically, put coils last - if doing xenon, this should be fine. If doing
%1H, need to move coils.
%Syntax in try works for KUMC data (VE11C) - ResonantNucleus is no longer a
%field in XA20A, so need something else. Haven't found nucleus yet, so use
%coil
try
    if strcmp(twix_obj.hdr.Spice.ResonantNucleus,'1H') && Dim == 1 
        raw = permute(raw,[1,3,2]);
    elseif strcmp(twix_obj.hdr.Spice.ResonantNucleus,'1H') && Dim == 0 
        raw = permute(raw,[1 3 4 2]);
    end
catch
    if ~strcmp(twix_obj.hdr.Dicom.TransmittingCoil,'129Xe_Vest') && Dim == 1 
        raw = permute(raw,[1,3,2]);
    elseif ~strcmp(twix_obj.hdr.Dicom.TransmittingCoil,'129Xe_Vest') && Dim == 0 
        raw = permute(raw,[1 3 4 2]);
    end
end
%% Figure out Trajectories:
traj = Tools.seek_spiral_traj(twix_obj);

%% Optional Trajectory Correction
% dwell = twix_obj.hdr.MeasYaps.sRXSPEC.alDwellTime{1,1}/1000;
% traj = traj_delay_correction(traj,dwell,1*dwell);

%% Density Compensation
%If doing 3D, we can use Jim Pipe's really good, fast recon. Otherwise,
%have to Use Scott Robertson's (still good) recon
if Dim == 1
    nIter = 10; %10 generally seems a good number of iterations for density compensation
    DCF = Recon.get_DCF(traj,ImSize,nIter);
else
    traj = Tools.column_traj(traj);
    DCF = Tools.get_DCF_Robertson(traj,ImSize,nIter);
end

%% Reconstruction

%I think this should take into account most situations.
if Dim == 1
    Image = zeros([ImSize ImSize ImSize size(raw,3)]);
    for i = 1:size(raw,3)
        Image(:,:,:,i) = Recon.pipe_recon(ImSize,raw(:,:,i),traj,DCF,i,size(raw,3));
        %If this is xenon, each image matters, so write out as we go:
        try
            if ~strcmp(twix_obj.hdr.Spice.ResonantNucleus,'1H')
                niftiwrite(squeeze(abs(Image(:,:,:,i))),fullfile(write_path,['Reconstructed_Image_for_Index_' num2str(i)]),'Compressed',true);
            end
        catch
            if strcmp(twix_obj.hdr.Dicom.TransmittingCoil,'129Xe_Vest')
                niftiwrite(squeeze(abs(Image(:,:,:,i))),fullfile(write_path,['Reconstructed_Image_for_Index_' num2str(i)]),'Compressed',true);
            end
        end
    end
    %If 1H, coil combine at the end and write out:
    try
        if strcmp(twix_obj.hdr.Spice.ResonantNucleus,'1H')
            Image = Tools.soscoilcombine(Image);
            niftiwrite(Image,fullfile(write_path,'Reconstructed_Image'),'Compressed',true);
        end
    catch
        if strcmp(twix_obj.hdr.Dicom.TransmittingCoil,'129Xe_Vest')
            Image = Tools.soscoilcombine(Image);
            niftiwrite(Image,fullfile(write_path,'Reconstructed_Image'),'Compressed',true);
        end
    end
else
    %otherwise, loop through coils and slices, doing Robertson
    %Reconstruction - This doesn't take into account a 2D Diffusion Image
    Image = zeros([ImSize(1) ImSize(2) size(raw,3) size(raw,4)]);
    for i = 1:size(raw,3)
        for j = 1:size(raw,4)
            fid_recon = reshape(raw(:,:,i,j),1,[])';
            [Image(:,:,i,j),~] = DCF.reconstruct(fid_recon, traj);
        end
    end
    Image = Tools.soscoilcombine(Image);
    niftiwrite(Image,fullfile(write_path,'Reconstructed_Image'),'Compressed',true);
end
%Clean up
fclose all;