% ---------------------------------------------------------------
%  DICOM Export for Retrospective app, DICOM info not available
% ---------------------------------------------------------------
function folderName = retroExportDicomMat(app,dcmExportDir,acqDur,tag,recoType)

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

dcmid = dicomuid;   % unique identifier
dcmid = dcmid(1:50);

% Create folder if not exist, and clear
folderName = [dcmExportDir,[filesep,'RETRO_DICOM_',num2str(nrFrames),'_',num2str(dimz),'_',num2str(nrDynamics),'_',tag]];
if (~exist(folderName, 'dir')); mkdir(folderName); end
delete([folderName,filesep,'*']);

% Variable flip-angle
if app.retroDataPars.VFA_size > 1
    dynamicLabel = '_flipangle_';
else
    dynamicLabel = '_dynamic_';
end


% export to Dicom images
cnt = 1;

for frame = 1:nrFrames

    for dyn = 1:nrDynamics

        for slice = 1:dimz

            fn = ['0000',num2str(cnt)];
            fn = fn(size(fn,2)-4:size(fn,2));
            fname = [folderName,filesep,fn,'_frame_',num2str(frame),'_slice_',num2str(slice),dynamicLabel,num2str(dyn),'.dcm'];
  
            dcmHeader = generate_dicomheader_mat;

            im = squeeze(cast(round(movie(frame,:,:,slice,dyn)),'uint16'));
      
            dicomwrite(im, fname, dcmHeader); %, "CreateMode","copy");

            cnt = cnt + 1;

        end

    end

end



    function dicomHeader = generate_dicomheader_mat
    % ---------------------------------------------------------------
    %           Generate DICOM header
    % ---------------------------------------------------------------

    try
        studyName = str2double(app.retroDataPars.filename(end-9:end-6));
    catch
        studyName = 1111;
    end

    % Flip angle
    if app.retroDataPars.VFA_size>1
        dcmHead.FlipAngle = app.retroDataPars.VFA_angles(dyn);
    else
        dcmHead.FlipAngle = app.retroDataPars.alpha;
    end

    % Aspect ratio and pixel x and y dimensions
    aspectRatio = app.retroDataPars.FOVf/8;
    if app.retroDataPars.PHASE_ORIENTATION
        pixely = app.retroDataPars.FOV/dimy;
        pixelx = aspectRatio*app.retroDataPars.FOV/dimx;
    else
        pixely = aspectRatio*app.retroDataPars.FOV/dimy;
        pixelx = app.retroDataPars.FOV/dimx;
    end

    % Date and time
    dt = datetime(app.retroDataPars.date,'InputFormat','dd-MMM-yyyy HH:mm:ss');
    year = num2str(dt.Year);
    month = ['0',num2str(dt.Month)]; month = month(end-1:end);
    day = ['0',num2str(dt.Day)]; day = day(end-1:end);
    date = [year,month,day];

    hour = ['0',num2str(dt.Hour)]; hour = hour(end-1:end);
    minute = ['0',num2str(dt.Minute)]; minute = minute(end-1:end);
    seconds = ['0',num2str(dt.Second)]; seconds = seconds(end-1:end);
    time = [hour,minute,seconds];

    if strcmp(recoType,'realtime')

        % ---------------------------------------------------------------
        %          realtime dataset
        % ---------------------------------------------------------------

        TR = 1000*(60/heartRate)/nrFrames;    % time between cardiac frames in ms
        TD = 1000*acqDur/nrDynamics;                % time between dynamics

        dcmHead.Filename = fname;
        dcmHead.FileModDate = app.retroDataPars.date;
        dcmHead.FileSize = dimy*dimx*2;
        dcmHead.Format = 'DICOM';
        dcmHead.FormatVersion = 3;
        dcmHead.Width = dimx;
        dcmHead.Height = dimy;
        dcmHead.BitDepth = 15;
        dcmHead.ColorType = 'grayscale';
        dcmHead.FileMetaInformationGroupLength = 178;
        dcmHead.FileMetaInformationVersion = uint8([0, 1])';
        dcmHead.MediaStorageSOPClassUID = '1.2.840.10008.5.1.4.1.1.4';
        dcmHead.TransferSyntaxUID = '1.2.840.10008.1.2.1';
        dcmHead.ImplementationClassUID = '1.2.826.0.9717382.3.0.3.6.0';
        dcmHead.ImplementationVersionName = 'OFFIS_DCMTK_360';
        dcmHead.SpecificCharacterSet = 'ISO_IR 100';
        dcmHead.ImageType = 'ORIGINAL\PRIMARY\';
        dcmHead.SOPClassUID = '1.2.840.10008.5.1.4.1.1.4';
        dcmHead.StudyDate = date;
        dcmHead.SeriesDate = date;
        dcmHead.AcquisitionDate = date;
        dcmHead.StudyTime = time;
        dcmHead.SeriesTime = time;
        dcmHead.AcquisitionTime = time;
        dcmHead.ContentTime = time;
        dcmHead.Modality = 'MR';
        dcmHead.Manufacturer = 'MR Solutions Ltd';
        dcmHead.InstitutionName = 'Amsterdam UMC';
        dcmHead.InstitutionAddress = 'Amsterdam, Netherlands';
        dcmHead.ReferringPhysicianName.FamilyName = 'Amsterdam UMC preclinical MRI';
        dcmHead.ReferringPhysicianName.GivenName = '';
        dcmHead.ReferringPhysicianName.MiddleName = '';
        dcmHead.ReferringPhysicianName.NamePrefix = '';
        dcmHead.ReferringPhysicianName.NameSuffix = '';
        dcmHead.StationName = 'MRI Scanner';
        dcmHead.StudyDescription = 'CINE';
        dcmHead.SeriesDescription = '';
        dcmHead.InstitutionalDepartmentName = 'Amsterdam UMC preclinical MRI';
        dcmHead.PhysicianOfRecord.FamilyName = 'Amsterdam UMC preclinical MRI';
        dcmHead.PhysicianOfRecord.GivenName = '';
        dcmHead.PhysicianOfRecord.MiddleName = '';
        dcmHead.PhysicianOfRecord.NamePrefix = '';
        dcmHead.PhysicianOfRecord.NameSuffix = '';
        dcmHead.PerformingPhysicianName.FamilyName = 'Amsterdam UMC preclinical MRI';
        dcmHead.PerformingPhysicianName.GivenName = '';
        dcmHead.PerformingPhysicianName.MiddleName = '';
        dcmHead.PerformingPhysicianName.NamePrefix = '';
        dcmHead.PerformingPhysicianName.NameSuffix = '';
        dcmHead.PhysicianReadingStudy.FamilyName = 'AMC preclinical MRI';
        dcmHead.PhysicianReadingStudy.GivenName = '';
        dcmHead.PhysicianReadingStudy.MiddleName = '';
        dcmHead.PhysicianReadingStudy.NamePrefix = '';
        dcmHead.PhysicianReadingStudy.NameSuffix = '';
        dcmHead.OperatorName.FamilyName = 'manager';
        dcmHead.AdmittingDiagnosesDescription = '';
        dcmHead.ManufacturerModelName = 'MRS7024';
        dcmHead.ReferencedSOPClassUID = '';
        dcmHead.ReferencedSOPInstanceUID = '';
        dcmHead.ReferencedFrameNumber = [];
        dcmHead.DerivationDescription = '';
        dcmHead.FrameType = '';
        dcmHead.PatientName.FamilyName = 'Amsterdam UMC preclinical MRI';
        dcmHead.PatientID = '01';
        dcmHead.PatientBirthDate = date;
        dcmHead.PatientBirthTime = '';
        dcmHead.PatientSex = 'F';
        dcmHead.OtherPatientID = '';
        dcmHead.OtherPatientName.FamilyName = 'Amsterdam UMC preclinical MRI';
        dcmHead.OtherPatientName.GivenName = '';
        dcmHead.OtherPatientName.MiddleName = '';
        dcmHead.OtherPatientName.NamePrefix = '';
        dcmHead.OtherPatientName.NameSuffix = '';
        dcmHead.PatientAge = '1';
        dcmHead.PatientSize = [];
        dcmHead.PatientWeight = 0.0300;
        dcmHead.Occupation = '';
        dcmHead.AdditionalPatientHistory = '';
        dcmHead.PatientComments = '';
        dcmHead.BodyPartExamined = '';
        dcmHead.SequenceName = 'CINE';
        dcmHead.SliceThickness = app.retroDataPars.SLICE_THICKNESS/dimz;
        dcmHead.KVP = 0;
        dcmHead.RepetitionTime = TR;     % time between frames
        dcmHead.EchoTime = 2;         % approximately correct
        dcmHead.InversionTime = 0;
        dcmHead.NumberOfAverages = app.retroDataPars.NO_AVERAGES;
        dcmHead.ImagedNucleus = '1H';
        dcmHead.MagneticFieldStrength = 7;
        dcmHead.SpacingBetweenSlices = 0;
        dcmHead.EchoTrainLength = 1;
        dcmHead.DeviceSerialNumber = '0034';
        dcmHead.PlateID = '';
        dcmHead.SoftwareVersion = '1.0.0.0';
        dcmHead.ProtocolName = 'CINE';
        dcmHead.SpatialResolution = [];
        dcmHead.TriggerTime = (frame-1)*TR;    % frame time = frame number times calculated TR
        dcmHead.DistanceSourceToDetector = [];
        dcmHead.DistanceSourceToPatient = [];
        dcmHead.FieldOfViewDimensions = [app.retroDataPars.FOV app.retroDataPars.FOV app.retroDataPars.SLICE_THICKNESS/dimz];
        dcmHead.ExposureTime = [];
        dcmHead.XrayTubeCurrent = [];
        dcmHead.Exposure = [];
        dcmHead.ExposureInuAs = [];
        dcmHead.FilterType = '';
        dcmHead.GeneratorPower = [];
        dcmHead.CollimatorGridName = '';
        dcmHead.FocalSpot = [];
        dcmHead.DateOfLastCalibration = '';
        dcmHead.TimeOfLastCalibration = '';
        dcmHead.PlateType = '';
        dcmHead.PhosphorType = '';
        dcmHead.AcquisitionMatrix = uint16([dimy 0 0 dimx])';

        dcmHead.AcquisitionDeviceProcessingDescription = '';
        dcmHead.CassetteOrientation = 'PORTRAIT';
        dcmHead.CassetteSize = '25CMX25CM';
        dcmHead.ExposuresOnPlate = 0;
        dcmHead.RelativeXrayExposure = [];
        dcmHead.AcquisitionComments = '';
        dcmHead.PatientPosition = 'HFS';
        dcmHead.Sensitivity = [];
        dcmHead.FieldOfViewOrigin = [];
        dcmHead.FieldOfViewRotation = [];
        dcmHead.AcquisitionDuration = acqDur;
        dcmHead.StudyInstanceUID = dcmid(1:18);
        dcmHead.SeriesInstanceUID = [dcmid(1:18),'.',num2str(studyName)];
        dcmHead.StudyID = '01';
        dcmHead.SeriesNumber = studyName;
        dcmHead.AcquisitionNumber = 1;

        dcmHead.TriggerTime = (dyn-1)*TD + (frame-1)*TR;
        dcmHead.AcquisitionDuration = acqDur;
        dcmHead.InstanceNumber = (slice-1)*nrDynamics*nrFrames + (frame-1)*nrDynamics + dyn;          % instance number
        dcmHead.TemporalPositionIdentifier = (dyn-1)*nrFrames + frame;
        dcmHead.NumberOfTemporalPositions = nrDynamics * nrFrames;
        dcmHead.TemporalResolution = TR;
        dcmHead.ImagesInAcquisition = nrDynamics * nrFrames * dimz;

        dcmHead.ImagePositionPatient = [(pixely/2)-app.retroDataPars.FOV/2, (pixelx/2)-app.retroDataPars.FOV/2, 0.0]';
        dcmHead.ImageOrientationPatient = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0]';
        dcmHead.FrameOfReferenceUID = '';
        dcmHead.SliceLocation = (app.retroDataPars.SLICE_THICKNESS/dimz)*((slice-1)-(round(dimz/2)));
        dcmHead.ImageComments = '';
        dcmHead.TemporalPositionIndex = uint32([]);
        dcmHead.SamplesPerPixel = 1;
        dcmHead.PhotometricInterpretation = 'MONOCHROME2';
        dcmHead.PlanarConfiguration = 0;
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
        dcmHead.Format = 'DICOM';
        dcmHead.FormatVersion = 3;
        dcmHead.Width = dimx;
        dcmHead.Height = dimy;
        dcmHead.BitDepth = 15;
        dcmHead.ColorType = 'grayscale';
        dcmHead.FileMetaInformationGroupLength = 178;
        dcmHead.FileMetaInformationVersion = uint8([0, 1])';
        dcmHead.MediaStorageSOPClassUID = '1.2.840.10008.5.1.4.1.1.4';
        dcmHead.TransferSyntaxUID = '1.2.840.10008.1.2.1';
        dcmHead.ImplementationClassUID = '1.2.826.0.9717382.3.0.3.6.0';
        dcmHead.ImplementationVersionName = 'OFFIS_DCMTK_360';
        dcmHead.SpecificCharacterSet = 'ISO_IR 100';
        dcmHead.ImageType = 'ORIGINAL\PRIMARY\';
        dcmHead.SOPClassUID = '1.2.840.10008.5.1.4.1.1.4';
        dcmHead.StudyDate = date;
        dcmHead.SeriesDate = date;
        dcmHead.AcquisitionDate = date;
        dcmHead.StudyTime = time;
        dcmHead.SeriesTime = time;
        dcmHead.AcquisitionTime = time;
        dcmHead.ContentTime = time;
        dcmHead.Modality = 'MR';
        dcmHead.Manufacturer = 'MR Solutions Ltd';
        dcmHead.InstitutionName = 'Amsterdam UMC';
        dcmHead.InstitutionAddress = 'Amsterdam, Netherlands';
        dcmHead.ReferringPhysicianName.FamilyName = 'Amsterdam UMC preclinical MRI';
        dcmHead.ReferringPhysicianName.GivenName = '';
        dcmHead.ReferringPhysicianName.MiddleName = '';
        dcmHead.ReferringPhysicianName.NamePrefix = '';
        dcmHead.ReferringPhysicianName.NameSuffix = '';
        dcmHead.StationName = 'MRI Scanner';
        dcmHead.StudyDescription = 'CINE';
        dcmHead.SeriesDescription = '';
        dcmHead.InstitutionalDepartmentName = 'Amsterdam UMC preclinical MRI';
        dcmHead.PhysicianOfRecord.FamilyName = 'Amsterdam UMC preclinical MRI';
        dcmHead.PhysicianOfRecord.GivenName = '';
        dcmHead.PhysicianOfRecord.MiddleName = '';
        dcmHead.PhysicianOfRecord.NamePrefix = '';
        dcmHead.PhysicianOfRecord.NameSuffix = '';
        dcmHead.PerformingPhysicianName.FamilyName = 'Amsterdam UMC preclinical MRI';
        dcmHead.PerformingPhysicianName.GivenName = '';
        dcmHead.PerformingPhysicianName.MiddleName = '';
        dcmHead.PerformingPhysicianName.NamePrefix = '';
        dcmHead.PerformingPhysicianName.NameSuffix = '';
        dcmHead.PhysicianReadingStudy.FamilyName = 'AMC preclinical MRI';
        dcmHead.PhysicianReadingStudy.GivenName = '';
        dcmHead.PhysicianReadingStudy.MiddleName = '';
        dcmHead.PhysicianReadingStudy.NamePrefix = '';
        dcmHead.PhysicianReadingStudy.NameSuffix = '';
        dcmHead.OperatorName.FamilyName = 'manager';
        dcmHead.AdmittingDiagnosesDescription = '';
        dcmHead.ManufacturerModelName = 'MRS7024';
        dcmHead.ReferencedSOPClassUID = '';
        dcmHead.ReferencedSOPInstanceUID = '';
        dcmHead.ReferencedFrameNumber = [];
        dcmHead.DerivationDescription = '';
        dcmHead.FrameType = '';
        dcmHead.PatientName.FamilyName = 'Amsterdam UMC preclinical MRI';
        dcmHead.PatientID = '01';
        dcmHead.PatientBirthDate = date;
        dcmHead.PatientBirthTime = '';
        dcmHead.PatientSex = 'F';
        dcmHead.OtherPatientID = '';
        dcmHead.OtherPatientName.FamilyName = 'Amsterdam UMC preclinical MRI';
        dcmHead.OtherPatientName.GivenName = '';
        dcmHead.OtherPatientName.MiddleName = '';
        dcmHead.OtherPatientName.NamePrefix = '';
        dcmHead.OtherPatientName.NameSuffix = '';
        dcmHead.PatientAge = '1';
        dcmHead.PatientSize = [];
        dcmHead.PatientWeight = 0.0300;
        dcmHead.Occupation = '';
        dcmHead.AdditionalPatientHistory = '';
        dcmHead.PatientComments = '';
        dcmHead.BodyPartExamined = '';
        dcmHead.SequenceName = 'CINE';
        dcmHead.SliceThickness = app.retroDataPars.SLICE_THICKNESS/dimz;
        dcmHead.KVP = 0;
        dcmHead.RepetitionTime = TR;     % time between frames
        dcmHead.EchoTime = 2;         % approximately correct
        dcmHead.InversionTime = 0;
        dcmHead.NumberOfAverages = app.retroDataPars.NO_AVERAGES;
        dcmHead.ImagedNucleus = '1H';
        dcmHead.MagneticFieldStrength = 7;
        dcmHead.SpacingBetweenSlices = 0;
        dcmHead.EchoTrainLength = 1;
        dcmHead.DeviceSerialNumber = '0034';
        dcmHead.PlateID = '';
        dcmHead.SoftwareVersion = '1.0.0.0';
        dcmHead.ProtocolName = 'CINE';
        dcmHead.SpatialResolution = [];
        dcmHead.TriggerTime = (frame-1)*TR;    % frame time = frame number times calculated TR
        dcmHead.DistanceSourceToDetector = [];
        dcmHead.DistanceSourceToPatient = [];
        dcmHead.FieldOfViewDimensions = [app.retroDataPars.FOV app.retroDataPars.FOV app.retroDataPars.SLICE_THICKNESS/dimz];
        dcmHead.ExposureTime = [];
        dcmHead.XrayTubeCurrent = [];
        dcmHead.Exposure = [];
        dcmHead.ExposureInuAs = [];
        dcmHead.FilterType = '';
        dcmHead.GeneratorPower = [];
        dcmHead.CollimatorGridName = '';
        dcmHead.FocalSpot = [];
        dcmHead.DateOfLastCalibration = '';
        dcmHead.TimeOfLastCalibration = '';
        dcmHead.PlateType = '';
        dcmHead.PhosphorType = '';
        dcmHead.AcquisitionMatrix = uint16([dimy 0 0 dimx])';
        dcmHead.AcquisitionDeviceProcessingDescription = '';
        dcmHead.CassetteOrientation = 'PORTRAIT';
        dcmHead.CassetteSize = '25CMX25CM';
        dcmHead.ExposuresOnPlate = 0;
        dcmHead.RelativeXrayExposure = [];
        dcmHead.AcquisitionComments = '';
        dcmHead.PatientPosition = 'HFS';
        dcmHead.Sensitivity = [];
        dcmHead.FieldOfViewOrigin = [];
        dcmHead.FieldOfViewRotation = [];
        dcmHead.AcquisitionDuration = acqDur/nrDynamics;
        dcmHead.StudyInstanceUID = dcmid(1:18);
        dcmHead.SeriesInstanceUID = [dcmid(1:18),'.',num2str(studyName)];
        dcmHead.StudyID = '01';
        dcmHead.SeriesNumber = studyName;
        dcmHead.AcquisitionNumber = 1;
        dcmHead.InstanceNumber = (slice-1)*nrFrames + frame;          % instance number
        dcmHead.ImagePositionPatient = [(pixely/2)-app.retroDataPars.FOV/2, (pixelx/2)-app.retroDataPars.FOV/2, 0.0]';
        dcmHead.ImageOrientationPatient = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0]';
        dcmHead.FrameOfReferenceUID = '';
        dcmHead.TemporalPositionIdentifier = frame;
        dcmHead.NumberOfTemporalPositions = nrFrames;
        dcmHead.TemporalResolution = TR;
        dcmHead.ImagesInAcquisition = nrFrames*dimz;
        dcmHead.SliceLocation = (app.retroDataPars.SLICE_THICKNESS/dimz)*((slice-1)-(round(dimz/2)));
        dcmHead.ImageComments = '';
        dcmHead.TemporalPositionIndex = uint32([]);
        dcmHead.SamplesPerPixel = 1;
        dcmHead.PhotometricInterpretation = 'MONOCHROME2';
        dcmHead.PlanarConfiguration = 0;
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

    end

    dicomHeader = dcmHead;

    end



end