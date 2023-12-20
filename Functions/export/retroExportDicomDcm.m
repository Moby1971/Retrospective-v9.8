function folderName = retroExportDicomDcm(app, dcmdir)

% ---------------------------------------------------------------
% DICOM Export for Retrospective app
% Gustav Strijkers
% December 2023
%
% ---------------------------------------------------------------


% Movie
movie = app.r.movieExp;

% Phase orientation
if ~app.r.PHASE_ORIENTATION
    movie = permute(rot90(permute(movie,[2,3,4,1,5]),1),[4,1,2,3,5]);
    movie = flip(movie,3);
end

% Data info
nrFrames = size(movie,1);
dimy = size(movie,2);
dimx = size(movie,3);
dimz = size(movie,4);
nrDynamics = size(movie,5);
heartRate = app.r.meanHeartRate;
respRate = app.r.meanRespRate;
slope = double(app.r.rescaleSlope);
intercept = double(app.r.rescaleIntercept);

% Square pixel correction
aspectRatio = app.r.FOVf/8;
if app.r.PHASE_ORIENTATION
    pixely = app.r.FOV/dimy;
    pixelx = aspectRatio*app.r.FOV/dimx;
else
    pixely = aspectRatio*app.r.FOV/dimy;
    pixelx = app.r.FOV/dimx;
end
if pixely > pixelx
    ratio = pixely/pixelx;
    movie = matrixInterpolate(movie,[1 ratio 1 1 1],'cubic');
end
if pixelx > pixely
    ratio = pixelx/pixely;
    movie = matrixInterpolate(movie,[1 1 ratio 1 1],'cubic');
end
dimy = size(movie,2);
dimx = size(movie,3);

% Reading in the DICOM header information
listing = dir(fullfile(dcmdir, '*.dcm'));
dcmFilename = strcat(listing(1).folder,filesep,listing(1).name);
dcmHead = dicominfo(dcmFilename);
TextMessage(app,strcat("Reading DICOM info from ",dcmFilename));

% Create new directory
ready = false;
cnt = 1;
while ~ready
    folderName = strcat(app.dicomExportPath,'DICOM',filesep,num2str(dcmHead.SeriesNumber),'R',filesep,num2str(cnt),filesep);
    if ~exist(folderName, 'dir')
        mkdir(folderName);
        ready = true;
    end
    cnt = cnt + 1;
end

% Message export folder
app.TextMessage(strcat("DICOM export folder = ",folderName));

% Variable flip-angle
if app.r.VFA_size > 1
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
            dcmFilename = strcat(listing(loc).folder,filesep,listing(loc).name);
            dcmHead = dicominfo(dcmFilename);

            % Filename
            fn = strcat('0000',num2str(cnt));
            fn = fn(size(fn,2)-4:size(fn,2));
            fname = strcat(folderName,filesep,fn,'_frame_',num2str(frame),'_slice_',num2str(slice),dynamicLabel,num2str(dyn),'.dcm');

            % Generate a dicom header
            dcmHeader = generate_dicomheader_dcm;

            % The image
            im = squeeze(cast(round(movie(frame,:,:,slice,dyn)),'uint16'));

            % Write dicom file
            dicomwrite(im, fname, dcmHeader); %, "CreateMode","copy");

            cnt = cnt + 1;

        end

    end

end


    function dicomHeader = generate_dicomheader_dcm
        % ---------------------------------------------------------------
        %           Generate DICOM header
        % ---------------------------------------------------------------

        % flip angle
        if app.r.VFA_size>1
            dcmHead.FlipAngle = app.r.VFA_angles(dyn);
        end

        % Aspectratio and pixel x and y dimensions
        aspectRatio = app.r.FOVf/8;
        if app.r.PHASE_ORIENTATION
            pixely = app.r.FOV/dimy;
            pixelx = aspectRatio*app.r.FOV/dimx;
        else
            pixely = aspectRatio*app.r.FOV/dimy;
            pixelx = app.r.FOV/dimx;
        end

        if strcmp(app.RecoTypeDropDown.Value,'realtime')

            % ---------------------------------------------------------------
            %          realtime dataset
            % ---------------------------------------------------------------

            TR = 1000*(60/heartRate)/nrFrames;            % time between cardiac frames in ms
            TD = 1000*app.acqDur/nrDynamics;                  % time between dynamics

            dcmHead.ImageType = 'ORIGINAL\PRIMARY\M_FFE\M\FFE';
            dcmHead.ScanningSequence = 'GR';
            dcmHead.Filename = fname;
            dcmHead.FileModDate = app.r.date;
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
            dcmHead.NumberOfAverages = app.r.NO_AVERAGES;
            dcmHead.InversionTime = 0;
            dcmHead.ImagedNucleus = '1H';
            dcmHead.MagneticFieldStrength = 7;
            dcmHead.TriggerTime = (dyn-1)*TD + (frame-1)*TR;
            dcmHead.AcquisitionMatrix = uint16([dimy 0 0 dimx])';
            dcmHead.AcquisitionDeviceProcessingDescription = '';
            dcmHead.AcquisitionDuration = app.acqDur;
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

            if strcmp(app.RecoTypeDropDown.Value,'respiratory')
                TR = 1000*(60/respRate)/nrFrames;       % time between frames in ms
                heartRespRate = respRate;
            else
                TR = 1000*(60/heartRate)/nrFrames;      % time between frames in ms
                heartRespRate = heartRate;
            end

            dcmHead.ImageType = 'ORIGINAL\PRIMARY\M_FFE\M\FFE';
            dcmHead.ScanningSequence = 'GR';
            dcmHead.Filename = fname;
            dcmHead.FileModDate = app.r.date;
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
            dcmHead.NumberOfAverages = app.r.NO_AVERAGES;
            dcmHead.InversionTime = 0;
            dcmHead.ImagedNucleus = '1H';
            dcmHead.MagneticFieldStrength = 7;
            dcmHead.TriggerTime = (frame-1)*TR;    % frame time = frame number times calculated TR
            dcmHead.AcquisitionMatrix = uint16([dimy 0 0 dimx])';
            dcmHead.AcquisitionDeviceProcessingDescription = '';
            dcmHead.AcquisitionDuration = app.acqDur/nrDynamics;
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
            dcmHead.HeartRate = heartRespRate;
            dcmHead.NumberOfSlices = dimz;
            dcmHead.CardiacNumberOfImages = nrFrames;

        end

        dicomHeader = dcmHead;

    end

    function outputMatrix = matrixInterpolate(inputMatrix, scaling, varargin)

        % inputMatrx = n-dimensional matrix
        % scaling = scaling factor
        % varargin = 'linear', 'cubic' interpolation method

        N = ndims(inputMatrix);
        scaling = scaling(1:N);
        scaling(1,1:N) = scaling(:).';
        sz = size(inputMatrix);
        xvec = cell(1,N);
        yvec = cell(1,N);
        szy = nan(1,N);
        nonsing = true(1,N);

        for i = 1:N

            n = sz(i);

            if n==1 %for vector input
                nonsing(i) = 0;
                szy(i) = 1;
                continue
            end

            szy(i) = round(sz(i)*scaling(i));
            m = szy(i);

            xax = linspace(1/n/2, 1-1/n/2 ,n);
            xax = xax-.5;

            yax = linspace(1/m/2, 1-1/m/2 ,m);
            yax = yax-.5;

            xvec{i} = xax;
            yvec{i} = yax;

        end

        xvec = xvec(nonsing);
        yvec = yvec(nonsing);
        F = griddedInterpolant(xvec,squeeze(inputMatrix),varargin{:});
        outputMatrix = reshape(F(yvec),szy);

    end % matrixInterpolate



    % Find parameter in RPR file
    function outputString = rprFind(rprFileName,searchString) %#ok<DEFNU>

        try
            loc1 = strfind(rprFileName, searchString);
            loc2 = strfind(rprFileName(loc1:end),newline);
            loc3 = strfind(rprFileName(loc1:loc1+loc2(1)),'"');
            outputString = rprFileName(loc1+loc3(1):loc1+loc3(2)-2);
        catch
            outputString = '';
        end

    end


end