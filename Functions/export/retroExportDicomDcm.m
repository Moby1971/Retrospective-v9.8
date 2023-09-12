% ---------------------------------------------------------------
%           DICOM Export for Retrospective app
% ---------------------------------------------------------------
function folderName = retroExportDicomDcm(app,dcmdir,dcmExportDir,acqDur,recoType)

movie = app.retroRecoPars.movieExp;
      
if ~app.retroDataPars.PHASE_ORIENTATION
    movie = permute(rot90(permute(movie,[2,3,4,1,5]),1),[4,1,2,3,5]);
    movie = flip(movie,3);
end

% Data info
nrFrames = size(movie,1);
dimx = size(movie,2);
dimy = size(movie,3);
dimz = size(movie,4);
nrDynamics = size(movie,5);
heartRate = app.retroNavPars.meanHeartRate;
respRate = app.retroNavPars.meanRespRate;
slope = double(app.retroRecoPars.rescaleSlope);
intercept = double(app.retroRecoPars.rescaleIntercept);

% Reading in the DICOM header information
listing = dir(fullfile(dcmdir, '*.dcm'));
dcmFilename = [listing(1).folder,filesep,listing(1).name];
dcmHead = dicominfo(dcmFilename);
TextMessage(app,strcat('Reading DICOM info from',{' '},dcmFilename));

% Create folder if not exist, and delete folder content
dir1 = 'DICOM';
dir2 = strcat(num2str(dcmHead.SeriesNumber),'R');
dir3 = '1';
folderName = strcat(dcmExportDir,dir1,filesep,dir2,filesep,dir3,filesep);
if (~exist(folderName, 'dir')) 
    mkdir(folderName); 
end
delete([folderName,filesep,'*']);

% Variable flip-angle
if app.retroDataPars.VFA_size > 1
    dynamicLabel = '_flipangle_';
else
    dynamicLabel = '_dynamic_';
end

% Export the dicom images
cnt = 1;

for frame=1:nrFrames

    for dyn = 1:nrDynamics

        for slice = 1:dimz

            % Read header for specific slice location
            loc = dimz+1-slice;
            dcmFilename = [listing(loc).folder,filesep,listing(loc).name];
            dcmHead = dicominfo(dcmFilename);

            % Filename
            fn = ['0000',num2str(cnt)];
            fn = fn(size(fn,2)-4:size(fn,2));
            fname = strcat(folderName,filesep,fn,'_frame_',num2str(frame),'_slice_',num2str(slice),dynamicLabel,num2str(dyn),'.dcm');

            % Generate a dicom header
            dcmHeader = generate_dicomheader_dcm;

            % The image
            im = squeeze(cast(round(movie(frame,:,:,slice,dyn)),'uint16'));

            % Write dicom file
            dicomwrite(im, fname, dcmHeader, "CreateMode","copy");

            cnt = cnt + 1;

        end

    end

end


    function dicomHeader = generate_dicomheader_dcm
    % ---------------------------------------------------------------
    %           Generate DICOM header
    % ---------------------------------------------------------------

        % flip angle
        if app.retroDataPars.VFA_size>1
            dcmHead.FlipAngle = app.retroDataPars.VFA_angles(dyn);
        end

        % Aspectratio and pixel x and y dimensions
        aspectRatio = app.retroDataPars.FOVf/8;
        if app.retroDataPars.PHASE_ORIENTATION
            pixely = app.retroDataPars.FOV/dimy;
            pixelx = aspectRatio*app.retroDataPars.FOV/dimx;
        else
            pixely = aspectRatio*app.retroDataPars.FOV/dimy;
            pixelx = app.retroDataPars.FOV/dimx;
        end

        if strcmp(recoType,'realtime')

            % ---------------------------------------------------------------
            %          realtime dataset
            % ---------------------------------------------------------------

            TR = 1000*(60/heartRate)/nrFrames;            % time between cardiac frames in ms
            TD = 1000*acqDur/nrDynamics;                  % time between dynamics

            dcmHead.Filename = fname;
            dcmHead.FileModDate = app.retroDataPars.date;
            dcmHead.FileSize = dimy*dimx*2;
            dcmHead.Width = dimx;
            dcmHead.Height = dimy;
            dcmHead.BitDepth = 15;
            dcmHead.InstitutionName = 'Amsterdam UMC';
            dcmHead.ReferringPhysicianName.FamilyName = 'AMC preclinical MRI';
            dcmHead.StudyDescription = 'CINE';
            dcmHead.InstitutionalDepartmentName = 'Amsterdam UMC preclinical MRI';
            dcmHead.PhysicianOfRecord.FamilyName = 'Amsterdam UMC preclinical MRI';
            dcmHead.PerformingPhysicianName.FamilyName = 'Amsterdam UMC preclinical MRI';
            dcmHead.PhysicianReadingStudy.FamilyName = 'Amsterdam UMC preclinical MRI';
            dcmHead.OperatorName.FamilyName = 'manager';
            dcmHead.ManufacturerModelName = 'MRS7024';
            dcmHead.ReferencedFrameNumber = [];
            dcmHead.RepetitionTime = TR;     % time between frames
            dcmHead.EchoTime = 2;         % approximately correct, unknown because of shortest TE option
            dcmHead.NumberOfAverages = app.retroDataPars.NO_AVERAGES;
            dcmHead.InversionTime = 0;
            dcmHead.ImagedNucleus = '1H';
            dcmHead.MagneticFieldStrength = 7;
            dcmHead.TriggerTime = (dyn-1)*TD + (frame-1)*TR;
            dcmHead.AcquisitionMatrix = uint16([dimy 0 0 dimx])';
            dcmHead.AcquisitionDeviceProcessingDescription = '';
            dcmHead.AcquisitionDuration = acqDur;
            dcmHead.InstanceNumber = (slice-1)*nrDynamics*nrFrames + (frame-1)*nrDynamics + dyn;          % instance number
            dcmHead.TemporalPositionIdentifier = (dyn-1)*nrFrames + frame;
            dcmHead.NumberOfTemporalPositions = nrDynamics * nrFrames;
            dcmHead.TemporalResolution = TR;
            dcmHead.ImagesInAcquisition = nrDynamics * nrFrames * dimz;
            dcmHead.TemporalPositionIndex = uint32([]);
            dcmHead.Rows = dimy;
            dcmHead.Columns = dimx;
            dcmHead.PixelSpacing = [pixely pixelx]';
            dcmHead.PixelAspectRatio = [1 1]';
            dcmHead.BitsAllocated = 16;
            dcmHead.BitsStored = 15;
            dcmHead.HighBit = 14;
            dcmHead.PixelRepresentation = 0;
            dcmHead.PixelPaddingValue = 0;
            dcmHead.RescaleIntercept = intercept;
            dcmHead.RescaleSlope = slope;
            dcmHead.HeartRate = heartRate;
            dcmHead.NumberOfSlices = dimz;
            dcmHead.CardiacNumberOfImages = nrFrames;

        else

            % ---------------------------------------------------------------
            %          CINE dataset
            % ---------------------------------------------------------------

            if strcmp(recoType,'respiratory')
                TR = 1000*(60/respRate)/nrFrames;       % time between frames in ms
            else
                TR = 1000*(60/heartRate)/nrFrames;      % time between frames in ms
            end

            dcmHead.Filename = fname;
            dcmHead.FileModDate = app.retroDataPars.date;
            dcmHead.FileSize = dimy*dimx*2;
            dcmHead.Width = dimx;
            dcmHead.Height = dimy;
            dcmHead.BitDepth = 15;
            dcmHead.InstitutionName = 'Amsterdam UMC';
            dcmHead.ReferringPhysicianName.FamilyName = 'AMC preclinical MRI';
            dcmHead.StudyDescription = 'CINE';
            dcmHead.InstitutionalDepartmentName = 'Amsterdam UMC preclinical MRI';
            dcmHead.PhysicianOfRecord.FamilyName = 'Amsterdam UMC preclinical MRI';
            dcmHead.PerformingPhysicianName.FamilyName = 'Amsterdam UMC preclinical MRI';
            dcmHead.PhysicianReadingStudy.FamilyName = 'Amsterdam UMC preclinical MRI';
            dcmHead.OperatorName.FamilyName = 'manager';
            dcmHead.ManufacturerModelName = 'MRS7024';
            dcmHead.ReferencedFrameNumber = [];
            dcmHead.RepetitionTime = TR;     % time between frames
            dcmHead.EchoTime = 2;         % approximately correct, unknown because of shortest TE option
            dcmHead.NumberOfAverages = app.retroDataPars.NO_AVERAGES;
            dcmHead.InversionTime = 0;
            dcmHead.ImagedNucleus = '1H';
            dcmHead.MagneticFieldStrength = 7;
            dcmHead.TriggerTime = (frame-1)*TR;    % frame time = frame number times calculated TR
            dcmHead.AcquisitionMatrix = uint16([dimy 0 0 dimx])';
            dcmHead.AcquisitionDeviceProcessingDescription = '';
            dcmHead.AcquisitionDuration = acqDur/nrDynamics;
            dcmHead.InstanceNumber = (slice-1)*nrFrames + frame;          % instance number
            dcmHead.TemporalPositionIdentifier = frame;     % frame number
            dcmHead.NumberOfTemporalPositions = nrFrames;
            dcmHead.TemporalResolution = TR;
            dcmHead.ImagesInAcquisition = nrFrames*dimz;
            dcmHead.TemporalPositionIndex = uint32([]);
            dcmHead.Rows = dimy;
            dcmHead.Columns = dimx;
            dcmHead.PixelSpacing = [pixely pixelx]';
            dcmHead.PixelAspectRatio = [1 1]';
            dcmHead.BitsAllocated = 16;
            dcmHead.BitsStored = 15;
            dcmHead.HighBit = 14;
            dcmHead.PixelRepresentation = 0;
            dcmHead.PixelPaddingValue = 0;
            dcmHead.(dicomlookup('0028','1052')) = intercept;
            dcmHead.RescaleSlope = slope;
            dcmHead.HeartRate = heartRate;
            dcmHead.NumberOfSlices = dimz;
            dcmHead.CardiacNumberOfImages = nrFrames;

        end

        dicomHeader = dcmHead;

    end



end