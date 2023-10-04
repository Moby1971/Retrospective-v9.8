classdef retro 

    % ---------------------------------------------------------
    %
    % Data and parameter class for retrospective app
    %
    % Gustav Strijkers
    % g.j.strijkers@amsterdamumc.nl
    % Okt 2023
    %
    % ---------------------------------------------------------


    %%

    % Methods
    % -------

    % retro(parameter, mridata)
    % checkNumberOfExperiments(obj, app)
    % checkNumberOfAverages(obj, app)
    % checkForMultiSlab(obj, app)
    % checkTR(obj, app)
    % checkForVFA(obj, app)
    % setVariableFlipAngles(obj)
    % acquisitionType(obj, app)
    % guessRecoType(obj)
    % importMRD(filename, reordering1, reordering2)
    % importB(obj, app)
    % readMrdFooter(obj, app, mrdfile)
    % writeDataToMrd(obj, filename, parameters)
    % makeMrdFooter(obj, par)
    % readRPRfile(obj, app, filename)
    % writeToRprFile(obj, filename)
    % makeRprFile(obj, par)
    % get3DimageShift(obj, image, app)
    % readSQLfile(obj, app, filename)
    % sqlParameters(obj, app)
    % imageOrientLabels(obj, app)
    % writePhysLog(obj, exportPath)
    % ExportRecoParametersFcn(obj, app)
    % extractData(obj)
    % assignBinTimes(obj, app)
    % assignBinTimesAbs(obj, app)
    % assignBinFrames(obj, app)
    % assignBinFramesAbs(obj, app)
    % fillKspace2D(obj, app)
    % fillKspace3D(obj, app)
    % fillKspace2Drealtime(obj, app)
    % fillKspace3Dp2roud(obj, app)
    % fillKspace2Dradial(obj, app)
    % fillKspace2DradialReg(obj, app)
    % fillKspace3Dute(obj, app)
    % reshapeKspace(obj)
    % kSpaceStats(obj, app)
    % getKspaceTrajectory(obj, app)
    % backToKspace(obj)
    % readFlipSwitch(obj, app)
    % autoDiscardNav(obj, app)
    % extractNavigator(obj)
    % determinePowerSpectrumPCA(obj)
    % filterNavPCA(obj)
    % filterShape(obj, app, type)
    % trigPointsHeart(obj)
    % trigPointsResp(obj)
    % calcHeartRate(obj)
    % calcRespRate(obj)
    % makeRespWindow(obj, app)
    % reco2D(obj, app)
    % reco3D(obj, app)
    % reco2Dradial(obj, app)
    % reco3Dute(obj, app)
    % PCAdenoise(obj, app)
    % unRing(obj, app)
    % normImages(obj)
    % shiftImages2D(obj, app)
    % shiftImages3D(obj, app)
    % determineMultiDimensions(obj)
    % recoSurFiles(obj, surpath, suffix, mrdfilename, rprfilename)
    % recoLcurve(obj, app)


    % Static Methods
    % --------------

    % circtukey2D(dimy, dimx, row, col, filterWidth)
    % circtukey3D(dimz, dimy, dimx, lev, row, col, filterWidth)
    % centricFilling(noViews2)
    % gauss(x, s, m)
    % fracCircShift(input, shiftsize)
    % peakFinder(x0, sel, thresh, extrema, includeEndpoints, interpolate)
    % partialDerivative3D(app, kTraj, xNew, calibSize)
    % partialDerivative2D(app, kTraj, Xnew, calibSize)
    % lowRankThresh3D(Xold, kSize, thresh)
    % lowRankThresh2D(Xold, kSize, thresh)
    % trajInterpolation(kSpaceOld, dShift)
    % im2row3D(im, winSize)
    % im2row2D(im, winSize)
    % row2im3D(mtx, imSize, winSize)
    % row2im2D(mtx, imSize, winSize)
    % vec(x)
    % fft3Dmri(x)
    % ifft3Dmri(X)
    % fft2Dmri(x)
    % ifft2Dmri(X)
    % image2Dshift(im, xShift, yShift)




    %%



    properties

        % Data
        data                                                    % raw k-space data
        mrdFooter                                               % original MRD file footer
        newMrdFooter                                            % new MRD file footer
        rprFile = []                                            % original RPR file data
        sqlFile = []                                            % SQL file data
        newRprFile                                              % new RPR file data
        filename                                                % MRD file name
        exportDir                                               % export directory
        includeWindow                                           % data include window
        excludeWindow                                           % data exclude window
        nrKlines                                                % number of k-lines

        % K-space data
        raw                                                     % raw k-space data (without navigator)
        kSpace                                                  % sorted k-space data
        kSpaceMrd                                               % k-space for new MRD file
        kSpaceAvg                                               % k-space averages
        kSpaceTraj                                              % trajectory points for entire k-space
        trajectory                                              % k-space trajectory
        trajectoryFileName = "none"                             % trajectory filename
        gradTrajectory                                          % gradient calibration trajectory

        % Sequence parameters mostly from MRD file header
        PPL                                                     % PPL name
        NO_SAMPLES = 1                                          % number of samples
        NO_VIEWS = 1                                            % number of views / phase encoding 1
        NO_VIEWS_2 = 1                                          % number of views 2 / phase encoding 2
        oversample = 0                                          % oversample factor in readout direction
        slab_ratio = 100                                        % slab ratio (oversampling)
        EXPERIMENT_ARRAY = 1                                    % number of experiments
        nr_repetitions = 1                                      % number of repetitions / experiments
        NO_AVERAGES = 1                                         % number of averages
        NO_SLICES = 1                                           % number of slices
        SLICE_THICKNESS = 1                                     % slice thickness (mm)
        SLICE_SEPARATION = 1                                    % slice separations (mm)
        SLICE_INTERLEAVE = 1                                    % slice interleave value
        r_angle_var                                             % read angle
        p_angle_var                                             % phase encoding angle
        s_angle_var                                             % slice angle
        nr_coils = 1                                            % number of coils
        FOV = 30                                                % field of view (mm)
        PHASE_ORIENTATION = 0                                   % phase orientation 1 = hor. 0 = vert.
        FOVf = 8                                                % field of view factor/8 = aspect ratio
        aspectratio = 1                                         % aspectratio
        alpha = 20                                              % flip angle
        te = 2                                                  % echo time (ms)
        te_us = 0                                               % additional echo time (us)
        TE                                                      % echo time (ms) = te + te_us
        tr = 10                                                 % repetition time (ms)
        tr_extra_us = 0                                         % additional repetition time (us)
        TR                                                      % repetition time (ms) = tr + tr_extra_us
        ti = 1000                                               % inversion time (ms)
        VFA_angles = []                                         % flip angle array
        VFA_size = 0                                            % flip angle array size (0 = 1 flip angle)
        frame_loop_on                                           % CINE loop on (1) or off (0)
        radial_on = 0                                           % radial on (1) or off (0)
        spoke_increment = 1                                     % radial spoke angle increment (not always used)
        slice_nav = 0                                           % slice navigator on (1) or off (0)
        date                                                    % scan date
        pixelshift1 = 0                                         % image shift in views direction
        pixelshift2 = 0                                         % image shift in views_2 direction
        coil_scaling = 1                                        % coil intensity scaling parameter
        scanner = 'MRS'                                         % scanner type
        no_samples_nav = 10                                     % number of navigator samples
        fov_read_off = 0                                        % read-offset from MRD file, relative offset = value/4000
        fov_phase_off = 0                                       % phase-offset from MRD file, relative offset = value/4000
        fov_slice_off = 0                                       % slice-offset idem
        fov_offsets = [0 0 0]                                   % fov offsets
        SAMPLE_PERIOD                                           % sample period

        % K-space trajectory related
        pe1_order = 3                                           % phase-encoding order
        pe2_centric_on = 1                                      % phase-encoding 2 centric order on (1) or off (0)
        pe2_traj = 0                                            % phase-encoding 2 trajectory type
        gp_var_mul = []                                         % trajectory array
        interpFactor = 4                                        % k-space data interpolation factor

        % Navigator related
        primaryNavigatorPoint = 10                              % primary navigator point
        nrNavPointsDiscarded = 35                               % number of data points discarded after the navigator
        nrNavPointsUsed = 5                                     % number of navigator points used

        % Final dimensions
        dimx                                                    % image x dimension
        dimy                                                    % image y dimension
        dimz                                                    % image z dimension
        nrKsteps                                                % number of k-space steps / trajectory steps

        % Flags
        rprFlag = false                                         % RPR file available true / false
        sqlFlag = false                                         % SQL file available true / false
        validDataFlag = false                                   % valid data true / false
        multiCoilFlag = false                                   % multi coil data true / false
        multi2DFlag = false                                     % multi slice 2D data true / false
        vfaDataFlag = false                                     % variable flip angle data true / false

        % Data and reconstruction type
        dataType = '2D'                                         % type of data
        recoGuess = 'systolic function'                         % type of reconstruction

        % Parameters from SQL file
        SQLnumberOfSlices = 1                                   % number of slices
        SQLsliceGap = 0                                         % slice gap
        SQLangleX = 0                                           % angle X
        SQLangleY = 0                                           % angle Y
        SQLangleZ = 0                                           % angle Z
        SQLoffsetX = 0                                          % offset X
        SQLoffsetY = 0                                          % offset Y
        SQLoffsetZ = 0                                          % offset Z

        % Image shifts & orientations
        xShift = 0                                              % image shift in X direction
        yShift = 0                                              % image shift in Y direction
        zShift = 0                                              % image shift in Z direction
        LRvec = [1 0 0]'                                        % left-right orientation vector
        APvec = [0 1 0]'                                        % anterior-posterior orientation vector
        HFvec = [0 0 1]'                                        % head-feet orientation vector
        orientationLabels = [' ',' ',' ',' '];                  % orientation labels

        % Raw navigator
        navAmplitude                                            % navigator amplitude
        navPhase                                                % navigator phase (not used)
        upDown = 1                                              % flip the navigator up or down
        minDiscard = 15                                         % minimum number of discarded samples between navigator and echo
        maxDiscard = 60                                         % maximum number of discarded samples between navigator and echo
        defaultDiscard = 35                                     % default number of discarded samples between navigator and echo

        % Filtering
        physioFilterSettings                                    % filter settings for navigator
        detectedHR                                              % heart rate from navigator
        detectedRR                                              % respiratory rate from navigator
        powerSpectrum                                           % navigator power spectrum
        frequency                                               % navigator frequency range
        bandwidthHR                                             % navigator heart rate filter bandwidth
        bandwidthRR                                             % navigator respiration rate filter bandwidth
        heartNav                                                % heart rate navigator
        respNav                                                 % respiration rate navigator
        PCANav                                                  % principal component analysis of navigator

        % Triggering assignments
        splineFactor = 60                                       % data interpolation factor to prevent navigator discretization by TR
        heartTrigPoints                                         % heart trigger points in units of samples
        respTrigPoints                                          % respiratory trigger points in units of samples
        heartTrigTime                                           % heart trigger points in units of ms
        respTrigTime                                            % respiratory trigger points in units of ms
        respPercentage = 30                                     % percentage of data discarded during respiration
        respWindow                                              % Resp window in units of samples
        respWindowTime                                          % Resp start and end values in units of ms

        % Rates
        heartRateTime                                           % heart rate as function of time
        heartRateTimeFiltered                                   % filtered/smoothed heart rate as function of time
        respRateTime                                            % respiration rate as function of time
        respRateTimeFiltered                                    % filtered/smoothed respiration rate as function of time
        meanHeartRate                                           % average heart rate
        meanRespRate                                            % average respiration rate

        % Cardiac and respiratory binning
        cardBins                                                % time stamps of the cardiac bins in unit samples
        respBins                                                % time stamps of the respiratory bins in unit samples
        cardBinNrs                                              % array with bin assignments of data to cardiac frames
        respBinNrs                                              % array with bin assignments of data to respiratory frames

        % Reconstruction
        movieExp                                            % Movie for movie export
        movieApp                                            % Movie for viewing in the app
        senseMap                                            % Sense map
        rescaleSlope                                        % Dicom info RescaleSlope for image scaling
        rescaleIntercept                                    % Dicom info RescaleIntercept for image scaling
        multiSliceFlag = false                              % Multi-slice true or false
        multiDynamicFlag = false                            % Mutli-dynamic true or false
        totalVariation = 'T'                                % Total variation (T) or total generalized variation (G)
        maxImageValue = 32767                               % Maximum image value


    end % properties



    methods

        % ---------------------------------------------------------------------------------
        % Object constructor
        % ---------------------------------------------------------------------------------
        function obj = retro(parameter, mridata)

            if nargin == 2

                obj.data = mridata;

                if isfield(parameter,'filename')
                    obj.filename = parameter.filename;
                end

                if isfield(parameter,'PPL')
                    obj.PPL = parameter.PPL;
                end

                if isfield(parameter,'NO_SAMPLES')
                    obj.NO_SAMPLES = parameter.NO_SAMPLES;
                end

                if isfield(parameter,'NO_VIEWS')
                    obj.NO_VIEWS = parameter.NO_VIEWS;
                end

                if isfield(parameter,'NO_VIEWS_2')
                    obj.NO_VIEWS_2 = parameter.NO_VIEWS_2;
                end

                if isfield(parameter,'EXPERIMENT_ARRAY')
                    obj.EXPERIMENT_ARRAY = parameter.EXPERIMENT_ARRAY;
                    obj.nr_repetitions = parameter.EXPERIMENT_ARRAY;
                end

                if isfield(parameter,'NO_AVERAGES')
                    obj.NO_AVERAGES = parameter.NO_AVERAGES;
                end

                if isfield(parameter,'NO_SLICES')
                    obj.NO_SLICES = parameter.NO_SLICES;
                    obj.SQLnumberOfSlices = parameter.NO_SLICES;
                end

                if isfield(parameter,'SLICE_THICKNESS')
                    obj.SLICE_THICKNESS= parameter.SLICE_THICKNESS;
                end

                if isfield(parameter,'SLICE_SEPARATION')
                    obj.SLICE_SEPARATION = parameter.SLICE_SEPARATION;
                end

                if isfield(parameter,'SLICE_INTERLEAVE')
                    obj.SLICE_INTERLEAVE = parameter.SLICE_INTERLEAVE;
                end

                if isfield(parameter,'r_angle_var')
                    obj.r_angle_var = parameter.r_angle_var;
                end

                if isfield(parameter,'p_angle_var')
                    obj.p_angle_var = parameter.p_angle_var;
                end

                if isfield(parameter,'s_angle_var')
                    obj.s_angle_var = parameter.s_angle_var;
                end

                if isfield(parameter,'nr_coils')
                    obj.nr_coils = parameter.nr_coils;
                    if obj.nr_coils > 1
                        obj.multiCoilFlag = true;
                    else
                        obj.multiCoilFlag = false;
                    end
                end

                if isfield(parameter,'FOV')
                    obj.FOV = parameter.FOV;
                end

                if isfield(parameter,'PHASE_ORIENTATION')
                    obj.PHASE_ORIENTATION = parameter.PHASE_ORIENTATION;
                end

                if isfield(parameter,'FOVf')
                    obj.FOVf = parameter.FOVf;
                end

                obj.aspectratio = obj.FOVf/8;

                if isfield(parameter,'alpha')
                    obj.alpha = parameter.alpha;
                end

                if isfield(parameter,'te')
                    obj.te = parameter.te;
                end

                if isfield(parameter,'te_us')
                    obj.te_us = parameter.te_us;
                end

                obj.TE = obj.te + obj.te_us/1000;

                if isfield(parameter,'tr')
                    obj.tr = parameter.tr;
                end

                if isfield(parameter,'tr_extra_us')
                    obj.tr_extra_us = parameter.tr_extra_us;
                end

                obj.TR = obj.tr + obj.tr_extra_us/1000;

                if isfield(parameter,'pe1_order')
                    obj.pe1_order = parameter.pe1_order;
                end

                if isfield(parameter,'pe2_centric_on')
                    obj.pe2_centric_on = parameter.pe2_centric_on;
                end

                if isfield(parameter,'ti')
                    obj.ti = parameter.ti;
                end

                if isfield(parameter,'VFA_angles')
                    obj.VFA_angles = parameter.VFA_angles;
                end

                if isfield(parameter,'VFA_size')
                    obj.VFA_size = parameter.VFA_size;
                end

                if isfield(parameter,'frame_loop_on')
                    obj.frame_loop_on = parameter.frame_loop_on;
                end

                if isfield(parameter,'radial_on')
                    obj.radial_on = parameter.radial_on;
                end

                if isfield(parameter,'spoke_increment')
                    obj.spoke_increment = parameter.spoke_increment/10;
                end

                if isfield(parameter,'slice_nav')
                    obj.slice_nav = parameter.slice_nav;
                end

                if isfield(parameter,'no_samples_nav')
                    obj.no_samples_nav = parameter.no_samples_nav;
                end

                if isfield(parameter,'gp_var_mul')
                    obj.gp_var_mul = parameter.gp_var_mul;
                end

                if isfield(parameter,'date')
                    obj.date = parameter.date;
                end

                if isfield(parameter,'coil_scaling')
                    obj.coil_scaling = parameter.coil_scaling;
                end

                if isfield(parameter,'pixelshift1')
                    obj.pixelshift1 = parameter.pixelshift1;
                    obj.SQLoffsetY = parameter.pixelshift1;
                end

                if isfield(parameter,'pixelshift2')
                    obj.pixelshift2 = parameter.pixelshift2;
                    obj.SQLoffsetZ = parameter.pixelshift2;
                end

                if isfield(parameter,'pe2_traj')
                    obj.pe2_traj = parameter.pe2_traj;
                end

                if isfield(parameter,'scanner')
                    obj.scanner = parameter.scanner;
                end

                if isfield(parameter,'fov_read_off')
                    obj.fov_read_off = parameter.fov_read_off;
                end

                if isfield(parameter,'fov_phase_off')
                    obj.fov_phase_off = parameter.fov_phase_off;
                end

                if isfield(parameter,'fov_slice_off')
                    obj.fov_slice_off = parameter.fov_slice_off;
                end

                if isfield(parameter,'FOV_OFFSETS')
                    obj.fov_offsets = parameter.FOV_OFFSETS;
                end

                if isfield(parameter,'SAMPLE_PERIOD')
                    obj.SAMPLE_PERIOD = parameter.SAMPLE_PERIOD;
                end

                if isfield(parameter,'oversample')
                    obj.oversample = parameter.oversample;
                end

                if isfield(parameter,'slab_ratio')
                    obj.slab_ratio = parameter.slab_ratio;
                end

            end

        end % retro




        % ---------------------------------------------------------------------------------
        % Check whether there are sufficient experiments to peform a valid reconstruction
        % ---------------------------------------------------------------------------------
        function obj = checkNumberOfExperiments(obj, app)

            % Check if there is enough data to do some meaningful reconstruction
            % Kind of random number
            if (obj.EXPERIMENT_ARRAY * obj.NO_VIEWS * obj.NO_VIEWS_2 * obj.NO_SLICES) < 1000
                obj.validDataFlag = false;
                app.TextMessage('ERROR: Not enough k-lines for any reconstruction ...');
            end

        end % checkNumberOfExperiments




        % ---------------------------------------------------------------------------------
        % Check the number of averages
        % ---------------------------------------------------------------------------------
        function obj = checkNumberOfAverages(obj, app)

            % For retrospective imaging number of averages should be 1
            if obj.NO_AVERAGES > 1
                obj.validDataFlag = false;
                app.TextMessage('ERROR: Number of averages should be 1 ...');
            end

        end % checkNumberOfAverages




        % ---------------------------------------------------------------------------------
        % Check for multi-slab data
        % ---------------------------------------------------------------------------------
        function obj = checkForMultiSlab(obj, app)

            % 3D multi-slab data is not supported
            if contains(obj.dataType,{'3D','3Dp2roud'}) && obj.NO_SLICES > 1
                obj.validDataFlag = false;
                app.TextMessage('ERROR: Only 3D single-slab data supported ...');
            end

        end % checkForMultiSlab



        % ---------------------------------------------------------------------------------
        % Check for TR > 0
        % ---------------------------------------------------------------------------------
        function obj = checkTR(obj, app)

            % TR should not be zero (minimum in MR Solutions software)
            % A positive TR is needed to calculate heart and respiratory rates
            if obj.tr == 0
                obj.validDataFlag = false;
                app.TextMessage('ERROR: TR = 0 is not allowed ...');
            end

        end % checkTR



        % ---------------------------------------------------------------------------------
        % Check for variable flip angle data
        % ---------------------------------------------------------------------------------
        function obj = checkForVFA(obj, app)

            % Variable flip-angle data in 3D is supported
            if contains(obj.dataType,{'3D','3Dp2roud'}) && obj.VFA_size > 0
                obj.vfaDataFlag = true;
                obj = setVariableFlipAngles(obj);
                app.TextMessage(strcat('INFO:',{' '},num2str(obj.VFA_size),{' '},'flip angles detected ...'));
            end

            % Set number of flip angles to 1 if obj.VFA_size = 0
            if obj.VFA_size == 0
                app.NrFlipAnglesViewField.Value = 1;
            else
                app.NrFlipAnglesViewField.Value = obj.VFA_size;
            end

        end % checkForVFA




        % ---------------------------------------------------------------------------------
        % Set the variable flip angles, sort the angles in groups
        % ---------------------------------------------------------------------------------
        function obj = setVariableFlipAngles(obj)

            % The different flip angles in a variable flip-angle experiment can be ordered in different ways
            % Here they are sorted, such that they can be combined
            a = obj.VFA_size;
            b = unique(obj.VFA_angles(1:obj.VFA_size),'Stable');
            c = length(b);
            d = obj.VFA_angles(1:obj.VFA_size);
            if a == c
                nrFlipAngles = a;     % FA1, FA2, FA3, FA4, .... = dyn1, dyn2, dyn3, dyn4, ...
                lsFlipAngles = b;
            elseif mod(a,c) == 0
                nrFlipAngles = c;     % FA1, FA1, ..., FA2, FA2, ..., FA3, FA3, ... = dyn1, dyn2, dyn3, dyn4, ...
                lsFlipAngles = b;
            else
                nrFlipAngles = a;     % Each dynamic has its own flip-angle
                lsFlipAngles = d;
            end
            obj.VFA_size = nrFlipAngles;
            obj.VFA_angles = lsFlipAngles;

        end % setVariableFlipAngles




        % ---------------------------------------------------------------------------------
        % Check if the data is a valid file with a motion navigator acquisition
        % ---------------------------------------------------------------------------------
        function obj = acquisitionType(obj, app)

            % Determine acquisition type and put the data in the right order

            if obj.NO_VIEWS == 1 && obj.NO_VIEWS_2 == 1 && obj.EXPERIMENT_ARRAY > 1000

                % 3D UTE data

                obj.validDataFlag = true;

                obj.dataType = '3Dute';

                for i=1:obj.nr_coils
                    %                             spokes,1,1,X
                    obj.data{i} = permute(obj.data{i},[1,4,3,2]);
                end

                obj.primaryNavigatorPoint = 1;
                obj.nr_repetitions = size(obj.data{1},1);

            elseif obj.radial_on == 1

                % 2D radial data

                obj.validDataFlag = true;

                obj.dataType = '2Dradial';

                for i=1:obj.nr_coils
                    if ismatrix(obj.data{i})
                        %                                  1 1 Y X
                        obj.data{i} = permute(obj.data{i},[3,4,1,2]);
                    end
                    if ndims(obj.data{i}) == 3 && obj.NO_SLICES == 1
                        %                                 NR 1 Y X
                        obj.data{i} = permute(obj.data{i},[1,4,2,3]);
                    end
                    if ndims(obj.data{i}) == 3 && obj.NO_SLICES > 1
                        %                                  1 Z Y X
                        obj.data{i} = permute(obj.data{i},[4,1,2,3]);
                        obj.dataType = '2Dradialms';
                    end
                    if ndims(obj.data{i}) == 4 && obj.NO_SLICES > 1
                        %                                 NR Z Y X
                        obj.data{i} = permute(obj.data{i},[1,2,3,4]);
                        obj.dataType = '2Dradialms';
                    end

                end

                % Center the echo for the navigator
                % This is needed to know that the navigator is the center k-space point
                % For out-of-center points the navigator is disturbed by the radial frequency
                % Loop through each dynamic (time point) in the data

                pcorr = true;

                if pcorr

                    dims = size(obj.data{1},4);

                    for dynamic = 1:size(obj.data{1},1)

                        % Loop through each slice in the data
                        for slice = 1:size(obj.data{1},2)

                            % Loop through each spoke (k-space line) in the data
                            for spoke = 1:size(obj.data{1},3)

                                % Extract the current k-line (spoke) from the data for the current coil, dynamic, slice
                                tmpKline1 = squeeze(obj.data{i}(dynamic,slice,spoke,:));

                                % Interpolate the k-line to increase the resolution
                                tmpKline2 = interp(tmpKline1,obj.interpFactor);

                                % Find the index of the maximum absolute value of the interpolated k-line, which represents the center of k-space
                                [~,kCenter] = max(abs(tmpKline2));

                                % Calculate the amount of shift needed to center the k-line by subtracting the k-center from the midpoint of k-space
                                kShift = floor(dims/2)-kCenter/obj.interpFactor;

                                % Shift the k-line by the calculated amount, using a circular shift to avoid edge effects
                                tmpKline1 = retro.fracCircShift(tmpKline1,kShift);

                                % Replace the original k-line with the shifted one in the data for the current coil, dynamic, slice, spoke
                                obj.data{i}(dynamic,slice,spoke,:) = tmpKline1;

                            end

                        end

                    end

                end

                kSpaceSum = squeeze(sum(abs(obj.data{1}),[1,2,3]));     % Sum over all dimensions
                [~,kCenter] = max(kSpaceSum);                           % Determine the center
                obj.primaryNavigatorPoint = kCenter;                    % Navigator = k-space center
                obj.nrNavPointsUsed = 1;                                % Number of navigator points = 1
                obj.nr_repetitions = size(obj.data{1},1);               % Number of k-space repetitions

            elseif (obj.slice_nav == 1) && (obj.no_samples_nav > 0)

                % 3D Cartesian data

                obj.validDataFlag = true;

                if obj.NO_VIEWS_2 > 1

                    obj.dataType = '3D';

                    for i=1:obj.nr_coils
                        if obj.EXPERIMENT_ARRAY > 1
                            %                                 NR Z Y X
                            obj.data{i} = permute(obj.data{i},[1,3,2,4]);
                        else
                            %                                  1 Z Y X
                            obj.data{i} = permute(obj.data{i},[4,2,1,3]);
                        end
                    end

                    % Pseudo Cartesion k-space filling / P2ROUD
                    if obj.pe1_order == 4
                        obj.dataType = '3Dp2roud';
                    end

                else

                    % 2D single-slice Cartesian data

                    obj.dataType = '2D';

                    for i=1:obj.nr_coils
                        if ismatrix(obj.data{i})
                            %                                  1 1 Y X
                            obj.data{i} = permute(obj.data{i},[3,4,1,2]);
                        end
                        if ndims(obj.data{i}) == 3 && obj.NO_SLICES == 1
                            %                                 NR 1 Y X
                            obj.data{i} = permute(obj.data{i},[1,4,2,3]);
                        end
                        if ndims(obj.data{i}) == 3 && obj.NO_SLICES > 1
                            %                                  1 Z Y X
                            obj.data{i} = permute(obj.data{i},[4,1,2,3]);
                        end

                    end

                    % 2D multi-slice data

                    if obj.NO_SLICES > 1
                        obj.dataType = '2Dms';
                        obj.multi2DFlag = true;
                    else
                        obj.multi2DFlag = false;
                    end

                end

                % Set the navigator points
                obj.primaryNavigatorPoint = obj.no_samples_nav;
                if obj.nrNavPointsUsed > obj.primaryNavigatorPoint
                    obj.nrNavPointsUsed = obj.primaryNavigatorPoint;
                end

                % Number of k-space repetitions
                obj.nr_repetitions = size(obj.data{1},1);

            else

                % If not one of the above, data cannot be used
                obj.validDataFlag = false;

            end

            % Message the user on the type of data
            switch obj.dataType

                case '2D'
                    app.TextMessage('2D single-slice data ...');
                case '3D'
                    app.TextMessage('3D data ...');
                case '3Dp2roud'
                    app.TextMessage('3D P2ROUD data ...');
                case '2Dms'
                    app.TextMessage('2D multi-slice data ...');
                case '2Dradial'
                    app.TextMessage('2D radial data ...');
                case '2Dradialms'
                    app.TextMessage('2D multi-slice radial data ...');
                case '3Dute'
                    app.TextMessage('3D UTE data ...');

            end

        end % acquisitionType




        % ---------------------------------------------------------------------------------
        % Guess which reconstruction type (systolic function, diastolic
        % function, scout, VFA
        % ---------------------------------------------------------------------------------
        function obj = guessRecoType(obj)

            % Guess which type of reconstruction would be suitable for the data

            switch obj.dataType

                case {'2D','2Dradial'}

                    if obj.nr_repetitions > 200
                        obj.recoGuess = 'diastolic function';
                    else
                        obj.recoGuess = 'systolic function';
                    end

                case {'2Dms','2Dradialms'}

                    if obj.nr_repetitions > 200
                        obj.recoGuess = 'diastolic function';
                    else
                        obj.recoGuess = 'systolic function';
                    end

                    % Check if there is more than one slice orientation
                    % Then it is probably a scout measurement
                    slices = [obj.s_angle_var ; obj.p_angle_var ; obj.r_angle_var]';
                    nr_unique_slices = size(unique(slices,'rows'),1);
                    if nr_unique_slices > 1
                        obj.recoGuess = 'scout';
                    end

                case {'3D','3Dute','3Dp2roud'}

                    obj.recoGuess = 'systolic function';

                    if obj.vfaDataFlag
                        obj.recoGuess = 'variable flip-angle';
                    end

            end

        end % guessRecoType




        % ---------------------------------------------------------------------------------
        % Read MRD file
        % ---------------------------------------------------------------------------------
        function [im, dim, par, unsortedKspace] = importMRD(obj, filename, reordering1, reordering2) %#ok<INUSD> 

            % Description: Function to open multidimensional MRD/SUR files given a filename with PPR-parsing
            % Read in MRD and SUR file formats
            % Inputs: string filename, reordering1, reordering2
            % Reordering1, 2 is 'seq' or 'cen'
            % Reordering1 is for 2D (views)
            % Reordering2 is for 3D (views2)
            % Outputs: complex data, raw dimension [no_expts,no_echoes,no_slices,no_views,no_views_2,no_samples], MRD/PPR parameters
            % Author: Ruslan Garipov
            % Date: 01/03/2014 - swapped views and views2 dimension - now correct
            % 30 April 2014 - support for reading orientations added
            % 11 September 2014 - swapped views and views2 in the array (otherwise the images are rotated)
            % 13 October 2015 - scaling added as a parameter

            fid = fopen(filename,'r');      % Define the file id
            val = fread(fid,4,'int32');
            xdim = val(1);
            ydim = val(2);
            zdim = val(3);
            dim4 = val(4);
            fseek(fid,18,'bof');
            datatype=fread(fid,1, 'uint16');
            datatype = dec2hex(datatype);
            fseek(fid,48,'bof');
            scaling = fread(fid,1, 'float32');
            bitsperpixel = fread(fid,1, 'uchar');
            fseek(fid,152,'bof');
            val = fread(fid,2, 'int32');
            dim5 = val(1);
            dim6 = val(2);
            fseek(fid,256,'bof');
            text = fread(fid,256);
            no_samples = xdim;
            no_views = ydim;
            no_views_2 = zdim;
            no_slices = dim4;
            no_echoes = dim5;
            no_expts = dim6;

            % Read in the complex image data
            dim = [no_expts,no_echoes,no_slices,no_views_2,no_views,no_samples];

            if size(datatype,2)>1
                onlydatatype = datatype(2);
                iscomplex = 2;
            else
                onlydatatype = datatype(1);
                iscomplex = 1;
            end

            switch onlydatatype
                case '0'
                    dataformat = 'uchar';
                case '1'
                    dataformat = 'schar';
                case '2'
                    dataformat = 'short';
                case '3'
                    dataformat = 'int16';
                case '4'
                    dataformat = 'int32';
                case '5'
                    dataformat = 'float32';
                case '6'
                    dataformat = 'double';
                otherwise
                    dataformat = 'int32';
            end

            % Read the data
            num2Read = no_expts*no_echoes*no_slices*no_views_2*no_views*no_samples*iscomplex; %*datasize;
            [m_total, count] = fread(fid,num2Read,dataformat); % Reading all the data at once

            % Check if expected size of data was read
            % If not, this means that the acquisition was prematurely stopped
            % and only part of the data is available
            if count < num2Read

                % Find the end of the data by looking for :PPL string
                textData = fileread(filename);
                targetText = ":PPL";
                amountOfData = strfind(textData,targetText);

                % Number of floats to read
                newNum2Read = (amountOfData-4)/4 - 512;

                % Reset the file position indicator to beginning of the data
                fseek(fid,512,'bof');

                % Read the data again
                [m_total, count] = fread(fid,newNum2Read ,dataFormat);

            end

            if iscomplex == 2
                a=1:count/2;
                m_real = m_total(2*a-1);
                m_imag = m_total(2*a);
                clear m_total;
                m_C_tmp = m_real+m_imag*1i;
                clear m_real m_imag;
            else
                m_C_tmp = m_total;
                clear m_total;
            end

            % Pre-allocate the expected size of m_C, in case of missing data
            m_C = zeros(num2Read,1);
            m_C(1:length(m_C_tmp)) = m_C_tmp;

            % The unsorted k-space
            unsortedKspace = m_C;

            % Centric k-space ordering views
            ord=1:no_views;
            if strcmp(reordering1,'cen')
                for g=1:no_views/2
                    ord(2*g-1)=no_views/2+g;
                    ord(2*g)=no_views/2-g+1;
                end
            end

            % Centric k-space ordering views2
            ord1 = 1:no_views_2;
            ord2 = ord1;
            if strcmp(reordering2,'cen')
                for g=1:no_views_2/2
                    ord2(2*g-1)=no_views_2/2+g;
                    ord2(2*g)=no_views_2/2-g+1;
                end
            end

            % Pre-allocating the data matrix speeds up this function significantly
            m_C_1=zeros(no_expts,no_echoes,no_slices,max(ord(:)),max(ord2(:)),no_samples);
            n = 0;
            for a=1:no_expts
                for b=1:no_echoes
                    for c=1:no_slices
                        for d=1:no_views
                            for e=1:no_views_2
                                m_C_1(a,b,c,ord(d),ord2(e),:) = m_C(1+n:no_samples+n); % sequential ordering
                                n=n+no_samples;
                            end
                        end
                    end
                end
            end

            clear ord;
            clear ord2;
            m_C = squeeze(m_C_1);
            clear m_C_1;
            im=m_C;
            clear m_C;
            sample_filename = char(fread(fid,120,'uchar')');
            ppr_text = char(fread(fid,Inf,'uchar')');
            fclose(fid);

            % Parse fields in ppr section of the MRD file
            if numel(ppr_text)>0

                cell_text = textscan(ppr_text,'%s','delimiter',char(13));
                PPR_keywords = {'BUFFER_SIZE','DATA_TYPE','DECOUPLE_FREQUENCY','DISCARD','DSP_ROUTINE','EDITTEXT','EXPERIMENT_ARRAY','FOV','FOV_READ_OFF','FOV_PHASE_OFF','FOV_SLICE_OFF','GRADIENT_STRENGTH','MULTI_ORIENTATION','Multiple Receivers','NO_AVERAGES','NO_ECHOES','NO_RECEIVERS','NO_SAMPLES','NO_SLICES','NO_VIEWS','NO_VIEWS_2','OBLIQUE_ORIENTATION','OBSERVE_FREQUENCY','ORIENTATION','PHASE_CYCLE','READ/PHASE/SLICE_SELECTION','RECEIVER_FILTER','SAMPLE_PERIOD','SAMPLE_PERIOD_2','SCROLLBAR','SLICE_BLOCK','SLICE_FOV','SLICE_INTERLEAVE','SLICE_THICKNESS','SLICE_SEPARATION','SPECTRAL_WIDTH','SWEEP_WIDTH','SWEEP_WIDTH_2','VAR_ARRAY','VIEW_BLOCK','VIEWS_PER_SEGMENT','SMX','SMY','SWX','SWY','SMZ','SWZ','VAR','PHASE_ORIENTATION','X_ANGLE','Y_ANGLE','Z_ANGLE','PPL','IM_ORIENTATION','IM_OFFSETS','FOV_OFFSETS'};
                %PPR_type_0 keywords have text fields only, e.g. ":PPL C:\ppl\smisim\1ge_tagging2_1.PPL"
                PPR_type_0 = [23 53];
                %PPR_type_1 keywords have single value, e.g. ":FOV 300"
                PPR_type_1 = [8 42:47];
                %PPR_type_2 keywords have single variable and single value, e.g. ":NO_SAMPLES no_samples, 16"
                PPR_type_2 = [4 7 15:21 25 31 33 41 49];
                PPR_type_3 = 48; % VAR keyword only (syntax same as above)
                PPR_type_4 = [28 29]; % :SAMPLE_PERIOD sample_period, 300, 19, "33.3 KHz  30 ?s" and SAMPLE_PERIOD_2 - read the first number=timeincrement in 100ns
                %PPR_type_5 keywords have single variable and two values, e.g. ":SLICE_THICKNESS gs_var, -799, 100"
                PPR_type_5 = [34 35];
                % KEYWORD [pre-prompt,] [post-prompt,] [min,] [max,] default, variable [,scale] [,further parameters ...];
                PPR_type_6 = [9:11 39 50:52]; % VAR_ARRAY and angles keywords
                PPR_type_7 = [54 55 56]; % IM_ORIENTATION and IM_OFFSETS (SUR only)

                par = struct('filename',filename);
                for j=1:size(cell_text{1},1)

                    char1 = char(cell_text{1}(j,:));
                    field_ = '';
                    if ~isempty(char1)
                        C = textscan(char1, '%*c%s %s', 1);
                        field_ = char(C{1});
                    end

                    % Find matching number in PPR_keyword array:
                    num = find(strcmp(field_,PPR_keywords));

                    if num>0

                        if find(PPR_type_3==num) % :VAR keyword
                            C = textscan(char1, '%*s %s %f');
                            field_title = char(C{1}); field_title(numel(field_title)) = [];
                            numeric_field = C{2};
                            par = setfield(par, field_title, numeric_field); %#ok<*SFLD>

                        elseif find(PPR_type_1==num)
                            C = textscan(char1, '%*s %f');
                            numeric_field = C{1};
                            par = setfield(par, field_, numeric_field);

                        elseif find(PPR_type_2==num)
                            C = textscan(char1, '%*s %s %f');
                            numeric_field = C{2};
                            par = setfield(par, field_, numeric_field);

                        elseif find(PPR_type_4==num)
                            C = textscan(char1, '%*s %s %n %n %s');
                            field_title = char(C{1}); field_title(numel(field_title)) = []; %#ok<*NASGU>
                            numeric_field = C{2};
                            par = setfield(par, field_, numeric_field);

                        elseif find(PPR_type_0==num)
                            C = textscan(char1, '%*s %[^\n]');
                            text_field = char(C{1}); %text_field = reshape(text_field,1,[]);
                            par = setfield(par, field_, text_field);

                        elseif  find(PPR_type_5==num)
                            C = textscan(char1, '%*s %s %f %c %f');
                            numeric_field = C{4};
                            par = setfield(par, field_, numeric_field);

                        elseif  find(PPR_type_6==num)
                            C = textscan(char1, '%*s %s %f %c %f', 100);
                            field_ = char(C{1}); field_(end) = [];% the name of the array
                            num_elements = C{2}; % the number of elements of the array
                            numeric_field = C{4};
                            multiplier = [];
                            for l=4:numel(C)
                                multiplier = [multiplier C{l}];
                            end
                            pattern = ':';
                            k=1;
                            tline = char(cell_text{1}(j+k,:));
                            while (~contains(tline, pattern))
                                tline = char(cell_text{1}(j+k,:));
                                arr = textscan(tline, '%*s %f', num_elements);
                                multiplier = [multiplier, arr{1}']; %#ok<*AGROW>
                                k = k+1;
                                tline = char(cell_text{1}(j+k,:));
                            end
                            par = setfield(par, field_, multiplier);

                        elseif find(PPR_type_7==num) % :IM_ORIENTATION keyword
                            C = textscan(char1, '%s %f %f %f');
                            field_title = char(C{1}); field_title(1) = [];
                            char2 = char(cell_text{1}(j+1,:));
                            C = textscan(char2, ',%f, %f, %f');
                            numeric_field = [C{1}, C{2}, C{3}];
                            par = setfield(par, field_title, numeric_field);

                        end

                    end

                end

                if isfield('OBSERVE_FREQUENCY','par')
                    C = textscan(par.OBSERVE_FREQUENCY, '%q');
                    text_field = char(C{1});
                    par.Nucleus = text_field(1,:);
                else
                    par.Nucleus = 'Unspecified';
                end

                par.datatype = datatype;
                file_pars = dir(filename);
                par.date = file_pars.date;

            else
                par = [];
            end

            par.scaling = scaling;

        end % ImportMRD




        % ---------------------------------------------------------------------------------
        % Import B-type scanner data
        % ---------------------------------------------------------------------------------
        function [obj, rawData, parameters] = importB(obj, app)

            % Import path
            importPath = app.mrdImportPath;

            % Parameters
            info1 = jCampRead(strcat(importPath,'acqp'));
            info2 = jCampRead(strcat(importPath,'method'));

            % Scanner type
            parameters.scanner = 'B-type';

            % Slices
            parameters.NO_SLICES = str2num(info1.NSLICES);
            parameters.SLICE_THICKNESS = str2num(info2.pvm.slicethick) * parameters.NO_SLICES;

            % Matrix in readout direction
            parameters.NO_SAMPLES = info1.acq.size(1) / 2;
            if isfield(info2.pvm,"matrix")
                parameters.NO_SAMPLES = info2.pvm.encmatrix(1);
            end

            % Matrix in phase encoding direction
            parameters.NO_VIEWS = info1.acq.size(2);
            if isfield(info2.pvm,"matrix")
                parameters.NO_VIEWS = info2.pvm.encmatrix(2);
            end

            % Phase encoding orientation
            parameters.PHASE_ORIENTATION = 1;
            pm1 = -1;
            pm2 = -1;

            % Determine how to flip the data for different orientations
            if isfield(info2.pvm,'spackarrreadorient')
                if strcmp(info2.pvm.spackarrreadorient(1:3),'L_R')
                    parameters.PHASE_ORIENTATION = 0;
                    flr =  1;
                    pm1 = +1;
                    pm2 = -1;
                end
                if strcmp(info2.pvm.spackarrreadorient(1:3),'A_P')
                    parameters.PHASE_ORIENTATION = 0;
                    flr =  1;
                    pm1 = +1;
                    pm2 = -1;
                end
                if strcmp(info2.pvm.spackarrreadorient(1:3),'H_F')
                    parameters.PHASE_ORIENTATION = 1;
                    flr =  0;
                    pm1 = -1;
                    pm2 = -1;
                end
            end

            % Matrix in 2nd phase encoding direction
            parameters.NO_VIEWS_2 = 1;
            parameters.pe2_centric_on = 0;

            % FOV
            parameters.FOV = info1.acq.fov(1)*10;
            parameters.FOV2 = info1.acq.fov(2)*10;
            parameters.FOVf = round(8*parameters.FOV2/parameters.FOV);

            % Sequence parameters
            parameters.tr = info1.acq.repetition_time;
            parameters.te = info1.acq.echo_time;
            parameters.alpha = str2num(info1.acq.flip_angle);
            parameters.NO_ECHOES = 1;
            parameters.NO_AVERAGES = str2num(info1.NA);

            % Other parameters
            parameters.date = datetime;
            parameters.nucleus = '1H';
            parameters.PPL = 'Navigator Sequence';
            parameters.filename = 'Proton';
            parameters.field_strength = str2num(info1.BF1)/42.58; %#ok<*ST2NM>
            parameters.filename = 111;
            parameters.pe1_order = 2;
            parameters.radial_on = 0;
            parameters.slice_nav = 1;

            % Number of navigator points
            if isfield(info2.pvm,"navpoints")
                parameters.no_samples_nav = str2num(info2.pvm.navpoints);
            else
                parameters.no_samples_nav = str2num(info2.NavSize) / 2;
            end

            % Number of receiver coils
            parameters.nr_coils = str2num(info2.pvm.encnreceivers);

            % Trajectory 1st phase encoding direction
            if isfield(info2.pvm,'ppggradamparray1')
                parameters.gp_var_mul = round(pm1 * info2.pvm.ppggradamparray1 * (parameters.NO_VIEWS / 2 - 0.5));
                if isfield(info2.pvm,'enczfaccel1') && isfield(info2.pvm,'encpftaccel1')
                    % * str2num(info2.pvm.enczfaccel1)
                    parameters.gp_var_mul = round(pm1 * info2.pvm.ppggradamparray1 * str2num(info2.pvm.encpftaccel1) * (parameters.NO_VIEWS / 2 - 0.5));
                else
                    parameters.gp_var_mul = round(pm1 * info2.pvm.ppggradamparray1 * (parameters.NO_VIEWS / 2 - 0.5));
                end
                parameters.pe1_order = 3;
            elseif isfield(info2.pvm,'encvalues1')
                parameters.gp_var_mul = round(pm1 * info2.pvm.encvalues1 * (parameters.NO_VIEWS / 2 - 0.5));
                if isfield(info2.pvm,'enczf') && isfield(info2.pvm,'encpft')
                    % * info2.pvm.enczf(2)
                    parameters.gp_var_mul = round(pm1 * info2.pvm.encvalues1 * info2.pvm.encpft(2) * (parameters.NO_VIEWS / 2 - 0.5));
                else
                    parameters.gp_var_mul = round(pm1 * info2.pvm.encvalues1 * (parameters.NO_VIEWS / 2 - 0.5));
                end
                parameters.pe1_order = 3;
            else
                % Assume linear
                parameters.pe1_order = 1;
            end

            % Data type
            datatype = 'int32';
            if isfield(info1.acq,'word_size')
                if strcmp(info1.acq.word_size,'_32_BIT')
                    datatype = 'int32';
                end
                if strcmp(info1.acq.word_size,'_16_BIT')
                    datatype = 'int16';
                end
            end

            % Read data
            if isfile(strcat(importPath,'fid.orig'))
                fileID = fopen(strcat(importPath,'fid.orig'));
            else
                fileID = fopen(strcat(importPath,'rawdata.job0'));
            end
            dataRaw = fread(fileID,datatype);
            fclose(fileID);
            kReal = dataRaw(1:2:end);
            kIm = dataRaw(2:2:end);
            kSpaceB = kReal + 1j*kIm;

            % Read navigator
            if isfile(strcat(importPath,'fid.NavFid'))
                fileID = fopen(strcat(importPath,'fid.NavFid'));
            else
                fileID = fopen(strcat(importPath,'rawdata.job1'));
            end
            navdata = fread(fileID,datatype);
            fclose(fileID);
            kReal = navdata(1:2:end);
            kIm = navdata(2:2:end);
            navKspace = kReal + 1j*kIm;

            % Phase offset
            if isfield(info1.acq,'phase1_offset')
                parameters.pixelshift1 = round(pm1 * parameters.NO_VIEWS * info1.acq.phase1_offset / parameters.FOV);
            end

            % 2D data
            if strcmp(info2.pvm.spatdimenum,"2D") || strcmp(info2.pvm.spatdimenum,"<2D>")

                % In case k-space trajectory was hacked by macro
                if isfield(info2.pvm,'ppggradamparray1')
                    parameters.NO_VIEWS = length(info2.pvm.ppggradamparray1);
                end

                % Imaging k-space
                singleRep = parameters.NO_SLICES*parameters.NO_SAMPLES*parameters.nr_coils*parameters.NO_VIEWS;
                parameters.EXPERIMENT_ARRAY = floor(length(kSpaceB)/singleRep);
                kSpaceB = kSpaceB(1:singleRep*parameters.EXPERIMENT_ARRAY,:);
                kSpaceB = reshape(kSpaceB,parameters.NO_SLICES,parameters.NO_SAMPLES,parameters.nr_coils,parameters.NO_VIEWS,parameters.EXPERIMENT_ARRAY);

                parameters.EXPERIMENT_ARRAY = size(kSpaceB,5);
                kSpaceB = permute(kSpaceB,[3,5,1,4,2]); % nc, nr, ns, np, nf

                % Flip readout if needed
                if flr
                    kSpaceB = flip(kSpaceB,5);
                end

                % Coil intensity scaling
                if isfield(info2.pvm,'encchanscaling')
                    for i = 1:parameters.nr_coils
                        kSpaceB(i,:) = kSpaceB(i,:) * info2.pvm.encchanscaling(i);
                    end
                end

                % Navigator
                navKspace = navKspace(1:parameters.NO_SLICES*parameters.no_samples_nav*parameters.nr_coils*parameters.NO_VIEWS*parameters.EXPERIMENT_ARRAY);
                navKspace = reshape(navKspace,parameters.NO_SLICES,parameters.no_samples_nav,parameters.nr_coils,parameters.NO_VIEWS,parameters.EXPERIMENT_ARRAY);
                navKspace = permute(navKspace,[3,5,1,4,2]);

                % Insert 34 point spacer to make the data compatible with MR Solutions data
                kSpacer = zeros(parameters.nr_coils,parameters.EXPERIMENT_ARRAY,parameters.NO_SLICES,parameters.NO_VIEWS,35);

                % Combine navigator + spacer + k-space
                rawB = cat(5,navKspace,kSpacer,kSpaceB);
                rawData = cell(parameters.nr_coils);
                for i = 1:parameters.nr_coils
                    rawData{i} = squeeze(rawB(i,:,:,:,:));
                    rawData{i} = reshape(rawData{i},parameters.EXPERIMENT_ARRAY,parameters.NO_SLICES,parameters.NO_VIEWS,parameters.NO_SAMPLES+35+parameters.no_samples_nav);
                end

            end

            % 3D data
            if strcmp(info2.pvm.spatdimenum,"3D") || strcmp(info2.pvm.spatdimenum,"<3D>")

                % 2nd phase encoding direction
                parameters.NO_VIEWS_2 = info1.acq.size(3);
                if isfield(info2.pvm,"matrix")
                    parameters.NO_VIEWS_2 = info2.pvm.encmatrix(3);
                end

                % Phase offset 2
                if isfield(info1.acq,'phase2_offset')
                    parameters.pixelshift2 = round(pm2 * parameters.NO_VIEWS_2 * info1.acq.phase2_offset/parameters.FOV);
                end

                % Slice thickness
                parameters.SLICE_THICKNESS = str2num(info2.pvm.slicethick);

                % 2nd phase encoding trajectory
                parameters.pe2_centric_on = 0;
                if isfield(info2.pvm,"encsteps2")
                    parameters.pe2_traj = info2.pvm.encsteps2;
                    parameters.pe2_centric_on = 2;
                end
                if isfield(info2.pvm,'encvalues2')
                    parameters.pe2_traj = round(info2.pvm.encvalues2 * (parameters.NO_VIEWS_2/2-0.5));
                    parameters.pe2_centric_on = 2;
                end

                % K-space
                kSpaceB = reshape(kSpaceB,parameters.nr_coils,parameters.NO_SAMPLES,parameters.NO_VIEWS,parameters.NO_VIEWS_2,[]);
                parameters.EXPERIMENT_ARRAY = size(kSpaceB,5);
                kSpaceB = permute(kSpaceB,[1,5,4,3,2]);

                % Flip readout if needed
                if flr
                    kSpaceB = flip(kSpaceB,5);
                end

                % Coil intesnity scaling
                if isfield(info2.pvm,'encchanscaling')
                    for i = 1:parameters.nr_coils
                        kSpaceB(i,:) = kSpaceB(i,:) * info2.pvm.encchanscaling(i);
                    end
                end

                % Navigator
                navKspace = reshape(navKspace,parameters.nr_coils,parameters.no_samples_nav,parameters.NO_VIEWS,parameters.NO_VIEWS_2,parameters.EXPERIMENT_ARRAY);
                navKspace = permute(navKspace,[1,5,4,3,2]);

                % Insert 34 point spacer to make the data compatible with MR Solutions data
                kSpacer = zeros(parameters.nr_coils,parameters.EXPERIMENT_ARRAY,parameters.NO_VIEWS_2,parameters.NO_VIEWS,35);

                % Combine navigator + spacer + k-space
                rawB = cat(5,navKspace,kSpacer,kSpaceB);
                for i = 1:parameters.nr_coils
                    rawData{i} = squeeze(rawB(i,:,:,:,:));
                    rawData{i} = reshape(rawData{i},parameters.EXPERIMENT_ARRAY,parameters.NO_VIEWS_2,parameters.NO_VIEWS,parameters.NO_SAMPLES+35+parameters.no_samples_nav);
                end

            end


            % Read reco files to a structure
            function struct = jCampRead(filename) %#ok<STOUT>

                % Open file read-only big-endian
                fid = fopen(filename,'r','b');
                skipLine = 0;

                % Loop through separate lines
                if fid~=-1

                    while 1

                        if skipLine
                            line = nextLine;
                            skipLine = 0;
                        else
                            line = fgetl(fid);
                        end

                        % Testing the text lines
                        while length(line) < 2
                            line = fgetl(fid);
                        end

                        % Parameters and optional size of parameter are on lines starting with '##'
                        if line(1:2) == '##' %#ok<*BDSCA>

                            % Parameter extracting and formatting
                            % Read parameter name
                            paramName = fliplr(strtok(fliplr(strtok(line,'=')),'#'));

                            % Check for illegal parameter names starting with '$' and correct (Matlab does not accepts variable names starting with $)
                            if paramName(1) == '$'
                                paramName = paramName(2:length(paramName));
                                % Check if EOF, if true return
                            elseif paramName(1:3) == 'END'
                                break
                            end

                            % Parameter value formatting
                            paramValue = fliplr(strtok(fliplr(line),'='));

                            % Check if parameter values are in a matrix and read the next line
                            if paramValue(1) == '('

                                paramValueSize = str2num(fliplr(strtok(fliplr(strtok(paramValue,')')),'(')));

                                % Create an empty matrix with size 'paramvaluesize' check if only one dimension
                                if ~isempty(paramValueSize)

                                    if size(paramValueSize,2) == 1
                                        paramValueSize = [paramValueSize,1];
                                    end

                                    % Read the next line
                                    nextLine = fgetl(fid);

                                    % See whether next line contains a character array
                                    if nextLine(1) == '<'
                                        paramValue = fliplr(strtok(fliplr(strtok(nextLine,'>')),'<')); %#ok<*NASGU>
                                    elseif strcmp(nextLine(1),'L') || strcmp(nextLine(1),'A') || strcmp(nextLine(1),'H')
                                        paramValue = nextLine;
                                    else

                                        % Check if matrix has more then one dimension
                                        if paramValueSize(2) ~= 1

                                            paramValueLong = str2num(nextLine);
                                            while (length(paramValueLong)<(paramValueSize(1)*paramValueSize(2))) & (nextLine(1:2) ~= '##') %#ok<*AND2>
                                                nextLine = fgetl(fid);
                                                paramValueLong = [paramValueLong str2num(nextLine)];
                                            end

                                            if (length(paramValueLong) == (paramValueSize(1)*paramValueSize(2))) & (~isempty(paramValueLong))
                                                paramValue=reshape(paramValueLong,paramValueSize(1),paramValueSize(2));
                                            else
                                                paramValue = paramValueLong;
                                            end

                                            if length(nextLine) > 1
                                                if (nextLine(1:2) ~= '##')
                                                    skipLine = 1;
                                                end
                                            end

                                        else

                                            % If only 1 dimension just assign whole line to paramvalue
                                            paramValue = str2num(nextLine);
                                            if ~isempty(str2num(nextLine))
                                                while length(paramValue)<paramValueSize(1)
                                                    line = fgetl(fid);
                                                    paramValue = [paramValue str2num(line)];
                                                end
                                            end

                                        end

                                    end

                                else
                                    paramValue = '';
                                end

                            end

                            % Add paramvalue to structure.paramname
                            if isempty(findstr(paramName,'_'))
                                eval(['struct.' paramName '= paramValue;']); %#ok<*EVLDOT>
                            else
                                try
                                    eval(['struct.' lower(paramName(1:findstr(paramName,'_')-1)) '.' lower(paramName(findstr(paramName,'_')+1:length(paramName))) '= paramValue;']);
                                catch
                                    eval(['struct.' lower(paramName(1:findstr(paramName,'_')-1)) '.' datestr(str2num(paramName(findstr(paramName,'_')+1:findstr(paramName,'_')+2)),9) ...
                                        paramName(findstr(paramName,'_')+2:length(paramName)) '= paramValue;']); %#ok<*DATST,*FSTR>
                                end
                            end

                        elseif line(1:2) == '$$'
                            % The two $$ lines are not parsed for now
                        end

                    end

                    % Close file
                    fclose(fid);

                end

            end

        end % importB




        % ---------------------------------------------------------------------------------
        % Read the MRD footer
        % ---------------------------------------------------------------------------------
        function obj = readMrdFooter(obj, app, mrdfile)

            try

                % Read information from the header and footer first
                fid = fopen(mrdfile,'r');
                val = fread(fid,4,'int32');
                xdim = val(1);
                ydim = val(2);
                zdim = val(3);
                dim4 = val(4);
                fseek(fid,18,'bof');
                data_type=fread(fid,1, 'uint16');
                data_type = dec2hex(data_type);
                fseek(fid,152,'bof');
                val = fread(fid,2, 'int32');
                dim5 = val(1);
                dim6 = val(2);
                no_samples = xdim;
                no_views = ydim;
                no_views_2 = zdim;
                no_slices = dim4;
                no_echoes = dim5;
                no_expts = dim6;

                % Determine datatype
                if size(data_type,2)>1
                    onlydatatype = data_type(2);
                    iscomplex = 2;
                else
                    onlydatatype = data_type(1);
                    iscomplex = 1;
                end
                switch onlydatatype
                    case '0'
                        datasize = 1; % size in bytes
                    case '1'
                        datasize = 1; % size in bytes
                    case '2'
                        datasize = 2; % size in bytes
                    case '3'
                        datasize = 2; % size in bytes
                    case '4'
                        datasize = 4; % size in bytes
                    case '5'
                        datasize = 4; % size in bytes
                    case '6'
                        datasize = 8; % size in bytes
                    otherwise
                        datasize = 4; % size in bytes
                end

                % Fast forward to beginning of footer
                fseek(fid,512,'bof');
                num2read = no_expts*no_echoes*no_slices*no_views_2*no_views*no_samples*iscomplex;
                fseek(fid,num2read*datasize,'cof');

                % Read the footer
                obj.mrdFooter = char(fread(fid,Inf,'uchar')');
                fclose(fid);

            catch ME

                % If unsuccesful data is invalid
                obj.validDataFlag = false;
                app.TextMessage(ME.message);

            end

        end % readMrdFooter




        % ---------------------------------------------------------------------------------
        % Write MRD file
        % ---------------------------------------------------------------------------------
        function obj = writeDataToMrd(obj, filename, parameters)

            % Description: Function to convert multidimensional complex data to MRD format file
            % Author: Ruslan Garipov / MR Solutions Ltd
            % Date: 17/04/2020
            % Inputs: string filename, N-dimensional data matrix, dimensions structure with the following fields:
            % .NoExperiments
            % .NoEchoes
            % .NoSlices
            % .NoSamples
            % .NoViews
            % .NoViews2
            % Footer - int8 data type footer as copied from an MRD containing a copy of
            % The PPR file, including preceeding 120-byte zeros
            % Output: 1 if write was successful, 0 if not, -1 if failed early (e.g. the dimension checks)

            % Data: multidimensional, complex float, with dimensions arranged
            % Dimensions: structure

            kSpaceMRDdata = obj.kSpaceMrd;
            footer = obj.newMrdFooter;

            header1 = zeros(128,1);
            header1(1)  = parameters.NoSamples;
            header1(2)  = parameters.NoViews;
            header1(3)  = parameters.NoViews2;
            header1(4)  = parameters.NoSlices;
            header1(39) = parameters.NoEchoes;
            header1(40) = parameters.NoExperiments;

            % Set datatype - 'complex float'
            header1(5)  = hex2dec('150000');

            % Open new file for writing
            fid1 = fopen(filename,'wb');

            % Write 512 byte header
            fwrite(fid1,header1,'int32');

            switch obj.dataType

                case {'3D','3Dp2roud'}
                    kSpaceMRDdata = flip(permute(kSpaceMRDdata,[1,3,2,4,5,6,7]),1);

                otherwise
                    kSpaceMRDdata = flip(kSpaceMRDdata,1);
            end

            % Convert to 1D array with alternating real and imag part of the data
            temp = kSpaceMRDdata;
            temp = temp(:);
            a = real(temp);
            b = imag(temp);
            temp = transpose([a b]);
            temp = temp(:);

            % Write data at once
            fwrite(fid1,temp,'float32');

            % write the footer
            fwrite(fid1,footer,'int8');

            % close file
            fclose(fid1);

        end % writeDataToMrd




        % ---------------------------------------------------------------------------------
        % Make MRD footer
        % ---------------------------------------------------------------------------------
        function obj = makeMrdFooter(obj, par)

            inputFooter = obj.mrdFooter;

            % Parameter names
            parameters = {':NO_SAMPLES no_samples, ',':NO_VIEWS no_views, ',':NO_VIEWS_2 no_views_2, ', ...
                ':NO_ECHOES no_echoes, ',':EXPERIMENT_ARRAY no_experiments, ',':NO_AVERAGES no_averages, ', ...
                ':VAR pe1_order, ',':VAR slice_nav, ',':VAR radial_on, ', ...
                ':VAR frame_loop_on, ',':VAR tr, ',':VAR te, ', ...
                ':BATCH_SLICES batch_slices, ',':NO_SLICES no_slices, ', ...
                ':VAR VFA_size, ',':VAR ti, ',':VAR pe2_centric_on, ', ...
                ':VAR tr_extra_us, '
                };

            % Parameter replacement values
            replacePars = {par.NoSamples,par.NoViews,par.NoViews2, ...
                par.NoEchoes,par.NoExperiments,par.NoAverages, ...
                par.peorder,par.slicenav,par.radialon, ...
                par.frameloopon,par.tr,par.te, ...
                par.batchslices,par.NoSlices, ...
                par.vfasize, par.ti, par.pe2_centric_on, ...
                par.tr_extra_us
                };

            % Replace all simple valued parameters
            for i = 1:length(parameters)

                txt = parameters{i}; % parameter name
                var = replacePars{i}; % parameter value
                pos = strfind(inputFooter,txt); % find position of parameter name in the footer

                if ~isempty(pos) % if parameter name is found
                    try
                        oldTextLength = strfind(inputFooter(pos+length(txt):pos+length(txt)+6),char(13))-1; % find length of the old parameter value
                        if isempty(oldTextLength)
                            oldTextLength = strfind(inputFooter(pos+length(txt):pos+length(txt)+6),newline)-1;
                        end

                        newText = [num2str(var),'     '];
                        newText = newText(1:6);
                        inputFooter = replaceBetween(inputFooter,pos+length(txt),pos+length(txt)+oldTextLength-1,newText);
                    catch
                    end
                end

            end


            % VFA array replace
            newText = '';
            for i = 1:length(par.vfaangles)
                newText = newText + ", " + num2str(par.vfaangles(i));
                if mod(i,8) == 2
                    newText = newText + newline;
                end
            end
            txt = ':VAR_ARRAY VFA_angles, ';
            pos1 = strfind(inputFooter,txt);
            if ~isempty(pos1)
                newStr = extractAfter(inputFooter,pos1);
                pos2 = strfind(newStr,':');
                pos3 = pos1+pos2(1)-2;
                inputFooter = replaceBetween(inputFooter,pos1+length(txt)-2,pos3,newText);
            end

            % Return the object
            obj.newMrdFooter  = inputFooter;

        end % makeMrdFooter




        % ---------------------------------------------------------------------------------
        % Read RPR file
        % ---------------------------------------------------------------------------------
        function obj = readRPRfile(obj, app, filename)

            try
                fid = fopen(filename,'r');
                obj.rprFile = char(fread(fid,Inf,'uchar')');
                fclose(fid);
                obj.rprFlag = true;
            catch
                obj.rprFile = '';
                obj.rprFlag = false;
                app.TextMessage('WARNING: RPR file not found ...');
            end

        end % readRPRfile




        % ---------------------------------------------------------------------------------
        % Write RPR file
        % ---------------------------------------------------------------------------------
        function obj = writeToRprFile(obj, filename)

            fid = fopen(filename,'wb');
            fwrite(fid,obj.newRprFile,'int8');
            fclose(fid);

        end % writeToRprFile




        % ---------------------------------------------------------------------------------
        % Make RPR file
        % ---------------------------------------------------------------------------------
        function obj = makeRprFile(obj, par)
            % This function makes a new RPR file from an existing RPR file
            % and a parameter structure.  The parameter structure can contain
            % any fields that the RPR file can contain.  If a field is not
            % specified, the existing value is used.  If a field is specified,
            % the new value is used.  If a field is specified as [], the field
            % is removed from the RPR file.

            inputRpr = obj.rprFile;

            % Parameter names
            parameters = {
                ':EDITTEXT LAST_ECHO ',':EDITTEXT MAX_ECHO ', ...
                ':EDITTEXT LAST_EXPT ',':EDITTEXT MAX_EXPT ', ...
                ':EDITTEXT SAMPLES_DIM1 ',':EDITTEXT DATA_LENGTH1 ', ':EDITTEXT OUTPUT_SIZE1 ', ...
                ':EDITTEXT SAMPLES_DIM2 ',':EDITTEXT DATA_LENGTH2 ', ':EDITTEXT OUTPUT_SIZE2 ', ...
                ':EDITTEXT SAMPLES_DIM3 ',':EDITTEXT DATA_LENGTH3 ', ':EDITTEXT OUTPUT_SIZE3 ', ...
                ':EDITTEXT LAST_SLICE ',':EDITTEXT MAX_SLICE ', ...
                ':COMBOBOX FFT_DIM1 ',':COMBOBOX FFT_DIM2 ',':COMBOBOX FFT_DIM3 ', ...
                ':RADIOBUTTON VIEW_ORDER_1',':RADIOBUTTON VIEW_ORDER_2'
                };

            % Parameter replacement values
            replacePars = {par.NoEchoes,par.NoEchoes, ...
                par.NoExperiments, par.NoExperiments, ...
                par.NoSamples, par.NoSamples, par.NoSamples, ...
                par.NoViews, par.NoViews, par.NoViews, ...
                par.NoViews2, par.NoViews2, par.NoViews2, ...
                par.NoSlices, par.NoSlices, ...
                par.NoSamples, par.NoViews, par.NoViews2, ...
                par.View1order, par.View2order
                };

            % Search for all parameters and replace values
            for i = 1:length(parameters)

                txt = parameters{i};
                var = replacePars{i};

                pos = strfind(inputRpr,txt);

                if ~isempty(pos)

                    if ~isstring(var)

                        try
                            oldTxtLength = strfind(inputRpr(pos+length(txt):pos+length(txt)+15),char(13))-1;
                            if isempty(oldTxtLength)
                                oldTxtLength = strfind(inputRpr(pos+length(txt):pos+length(txt)+15),newline)-1;
                            end
                            newText = [num2str(var),'     '];
                            newText = newText(1:6);
                            inputRpr = replaceBetween(inputRpr,pos+length(txt),pos+length(txt)+oldTxtLength-1,newText);
                        catch
                        end

                    else

                        try
                            oldTxtLength = strfind(inputRpr(pos+length(txt):pos+length(txt)+15),char(13))-1;
                            if isempty(oldTxtLength)
                                oldTxtLength = strfind(inputRpr(pos+length(txt):pos+length(txt)+15),newline)-1;
                            end
                            newText = strcat(" ",var,"           ");
                            newText = extractBefore(newText,12);
                            inputRpr = replaceBetween(inputRpr,pos+length(txt),pos+length(txt)+oldTxtLength-1,newText);
                        catch
                        end

                    end

                end

            end

            obj.newRprFile = inputRpr;

        end % makeRprFile




        % ---------------------------------------------------------------------------------
        % Retrieve the 3D image shift for off-center and oblique Radial and P2ROUD sequences
        % ---------------------------------------------------------------------------------
        function obj = get3DimageShift(obj, image, app)

            % Image dimensions in pixels
            dx = size(image,2);
            dy = size(image,3);
            dz = size(image,4);

            % Calculate the shift
            for i = 1:length(obj.fov_read_off)
                relShiftX = dx*obj.fov_read_off(i)/4000;      % Relative offset, scaling from PPL file
                relShiftY = dy*obj.fov_phase_off(i)/4000;
                relShiftZ = dz*obj.fov_slice_off(i)/400;
            end

            % Different readout / phase depending on phase_orientation value
            if obj.PHASE_ORIENTATION

                shiftInX = +relShiftX;
                shiftInY = -relShiftY;
                shiftInZ = -relShiftZ;

            else

                shiftInX = -relShiftX;
                shiftInY = -relShiftY;
                shiftInZ = -relShiftZ;

            end

            % Report the values back / return the object
            obj.xShift = shiftInX;
            obj.yShift = shiftInY;
            obj.zShift = shiftInZ;

            % Readout shift already taken care of in PPL sequence
            if strcmp(obj.dataType,'3Dp2roud')
                obj.xShift = 0;
            end

            % Textmessage
            app.TextMessage(sprintf('Image shift ΔX = %.2f, ΔY = %.2f pixels, ΔZ = %.2f pixels ...',obj.xShift,obj.yShift,obj.zShift));

        end % get3DimageShift




        % ---------------------------------------------------------------------------------
        % Read SQL file
        % ---------------------------------------------------------------------------------
        function obj = readSQLfile(obj, app, filename)

            try
                fid = fopen(filename,'r');
                obj.sqlFile = char(fread(fid,Inf,'uchar')');
                fclose(fid);
                obj.sqlFlag = true;
            catch
                obj.sqlFile = '';
                obj.sqlFlag = false;
                app.TextMessage('WARNING: SQL file not found ...');
            end

        end % readSQLfile




        % ---------------------------------------------------------------------------------
        % Get some parameters from SQL file
        % ---------------------------------------------------------------------------------
        function obj = sqlParameters(obj, app)

            if obj.sqlFlag

                if ~isempty(obj.sqlFile)

                    sqlData = obj.sqlFile;

                    try

                        % Group settings
                        groupSettings = strfind(sqlData,'[GROUP SETTINGS]');
                        posStart = strfind(sqlData(groupSettings:end),'VALUES(');
                        posStart = posStart(1)+groupSettings+6;
                        posEnd = strfind(sqlData(posStart:end),')');
                        posEnd = posEnd(1)+posStart-2;
                        groupData = sqlData(posStart:posEnd);
                        values = textscan(groupData, '%f %s %f %f %f %f %f %f %f %f %f %f','Delimiter',',');

                        obj.SQLnumberOfSlices = values{4};
                        obj.SQLsliceGap = values{5};
                        obj.SQLangleX = values{6};
                        obj.SQLangleY = values{7};
                        obj.SQLangleZ = values{8};
                        obj.SQLoffsetX = values{9};
                        obj.SQLoffsetY = values{10};
                        obj.SQLoffsetZ = values{11};

                    catch ME

                        app.TextMessage(ME.message);
                        app.TextMessage('WARNING: Something went wrong analyzing SQL file group settings ...');
                        app.SetStatus(1);

                    end

                end

            end

        end % sqlParameters



        % ---------------------------------------------------------------------------------
        % Calculate image orientation labels
        % ---------------------------------------------------------------------------------
        function obj = imageOrientLabels(obj, app)

            obj.orientationLabels = [' ',' ',' ',' '];

            try

                if obj.sqlFlag

                    % Start from axial orientation, head first, supine
                    obj.LRvec = [ 1  0  0 ]';
                    obj.APvec = [ 0  1  0 ]';
                    obj.HFvec = [ 0  0  1 ]';

                    % Rotate the vectors according to the angle values
                    % Add a tiny angle to make the chance of hitting 45 degree angles for which orientation is indetermined very unlikely
                    tinyAngle = 0.00001;
                    obj.LRvec = rotz(obj.SQLangleZ+tinyAngle)*roty(-obj.SQLangleY+tinyAngle)*rotx(-obj.SQLangleX+tinyAngle)*obj.LRvec;
                    obj.APvec = rotz(obj.SQLangleZ+tinyAngle)*roty(-obj.SQLangleY+tinyAngle)*rotx(-obj.SQLangleX+tinyAngle)*obj.APvec;
                    obj.HFvec = rotz(obj.SQLangleZ+tinyAngle)*roty(-obj.SQLangleY+tinyAngle)*rotx(-obj.SQLangleX+tinyAngle)*obj.HFvec;

                    % Determine the orientation combination
                    % This is done by determining the main direction of the vectors
                    [~, indxLR1] = max(abs(obj.LRvec(:)));
                    [~, indxAP1] = max(abs(obj.APvec(:)));
                    indxLR2 = sign(obj.LRvec(indxLR1));
                    indxAP2 = sign(obj.APvec(indxAP1));
                    indxLR2(indxLR2 == -1) = 2;
                    indxAP2(indxAP2 == -1) = 2;

                    labelsPrimary   = [ 'L','R' ; 'A','P' ; 'F','H'];
                    labelsSecondary = [ 'R','L' ; 'P','A' ; 'H','F'];

                    % Sort the labels according to the starting orientation
                    labelsPrimary   = [labelsPrimary(indxLR1,indxLR2),labelsPrimary(indxAP1,indxAP2)];
                    labelsSecondary = [labelsSecondary(indxLR1,indxLR2),labelsSecondary(indxAP1,indxAP2)];

                    % Assign the labels
                    obj.orientationLabels = [labelsPrimary(1),labelsPrimary(2),labelsSecondary(1),labelsSecondary(2)];

                end

            catch ME

                app.TextMessage(ME.message);

            end

        end % imageOrientLabels




        % ---------------------------------------------------------------------------------
        % Write physiological data to files
        % ---------------------------------------------------------------------------------
        function obj = writePhysLog(obj,exportPath)

            % Write the cardiac triggering, respiratory triggering and window values to files

            if ~isempty(obj.heartTrigTime)
                writematrix(obj.heartTrigTime',strcat(exportPath,filesep,'cardtrigpoints.txt'),'Delimiter','tab');
            end

            if ~isempty(obj.heartRateTime)
                heartRate = [obj.TR*obj.heartTrigPoints(1:end-1)/1000 ; obj.heartRateTime]';
                writematrix(heartRate,strcat(exportPath,filesep,'heartrate.txt'),'Delimiter','tab');
            end

            if ~isempty(obj.respTrigTime)
                writematrix(obj.respTrigTime',strcat(exportPath,filesep,'resptrigpoints.txt'),'Delimiter','tab');
            end

            if ~isempty(obj.respRateTime)
                respRate = [obj.TR*obj.respTrigPoints(1:end-1)/1000 ; obj.respRateTime]';
                writematrix(respRate,strcat(exportPath,filesep,'resprate.txt'),'Delimiter','tab');
            end

            if ~isempty(obj.respWindowTime)
                writematrix(obj.respWindowTime,strcat(exportPath,filesep,'respwindow.txt'),'Delimiter','tab');
            end

        end % writePhysLog



        % Export the reconstruction settings
        function obj = ExportRecoParametersFcn(obj, app)

            % Write the reconstruction information to a txt file

            pars = strcat(...
                "------------------------- \n\n", ...
                "RETROSPECTIVE ", app.appVersion,"\n\n", ...
                "Gustav Strijkers\n", ...
                "Amsterdam UMC\n", ...
                "g.j.strijkers@amsterdamumc.nl\n\n", ...
                "------------------------- \n", ...
                "\nDATA \n\n", ...
                "file = ", app.mrdFileName, "\n", ...
                "\nNAVIGATOR \n\n", ...
                "primary = ", num2str(app.NavigatorEndEditField.Value), "\n", ...
                "#points = ", num2str(app.NavigatorNRPointsEditField.Value), "\n", ...
                "switch = ",app.NavigatorFlipSwitch.Value , "\n", ...
                "\nFILTERS \n\n", ...
                "heart = ", num2str(app.HeartEditField.Value), " bpm\n", ...
                "heart width = ", num2str(app.HeartWidthEditField.Value), " bpm\n", ...
                "respiration = ", num2str(app.RespirationEditField.Value), " bpm\n", ...
                "respiration width = ", num2str(app.RespirationWidthEditField.Value), " bpm\n", ...
                "respiration window = ", num2str(app.RespirationWindowEditField.Value), " %%\n", ...
                "expiration (0) / inspiration (1) = ",num2str(app.RespirationToggleCheckBox.Value), " \n", ...
                "ignore respiration = ",num2str(app.IgnoreRespirationToggleCheckBox.Value), " \n", ...
                "\nCINE\n\n", ...
                "type = ",app.RecoTypeDropDown.Value, "\n", ...
                "#frames = ",num2str(app.FramesEditField.Value), "\n", ...
                "#dynamics = ",num2str(app.DynamicsEditField.Value), "\n", ...
                "\nRECO \n\n", ...
                "Bart version = ",app.bartVersion, "\n", ...
                "WVxyz = ",num2str(app.WVxyzEditField.Value), "\n", ...
                "TVxyz = ",num2str(app.TVxyzEditField.Value), "\n", ...
                "LLRxyz = ",num2str(app.LLRxyzEditField.Value), "\n", ...
                "TVcine = ",num2str(app.TVcineEditField.Value), "\n", ...
                "TVdyn = ",num2str(app.TVdynEditField.Value), "\n", ...
                "ESPIRiT = ",num2str(app.ESPIRiTCheckBox.Value), "\n", ...
                "PCA filtering = ",num2str(app.PCACheckBox.Value), "\n", ...
                "PCA window = ",num2str(app.PCAwindowEditField.Value), "\n", ...
                "ringing filter = ",num2str(app.RingRmCheckBox.Value), "\n", ...
                "Matlab reco = ",num2str(app.MatlabRecoCheckBox.Value), "\n", ...
                "sharing = ",num2str(app.SharingEditField.Value), "\n", ...
                "\nINCLUDE & EXCLUDE WINDOWS\n\n", ...
                "include = ",num2str(app.IncludeCheckBox.Value), "\n", ...
                "include start = ",num2str(app.IncludeStartSlider.Value), " s\n", ...
                "include end = ",num2str(app.IncludeEndSlider.Value), " s\n", ...
                "exclude = ",num2str(app.ExcludeCheckBox.Value), "\n", ...
                "exclude start = ",num2str(app.ExcludeStartSlider.Value), " s\n", ...
                "exclude end = ",num2str(app.ExcludeEndSlider.Value), " s\n", ...
                "\nRATES\n\n", ...
                "respiration rate = ",num2str(app.FinalRespirationViewField.Value), " bpm\n", ...
                "heart rate = ",num2str(app.FinalHeartViewField.Value), " bpm\n" ...
                );

            if strcmp(obj.dataType,'3Dute')
                pars = strcat(pars, ...
                    "\n3D UTE\n\n", ...
                    "Gx delay = ",num2str(app.GxDelayEditField.Value),"\n",...
                    "Gy delay = ",num2str(app.GyDelayEditField.Value),"\n",...
                    "Gz delay = ",num2str(app.GzDelayEditField.Value),"\n",...
                    "offset = ",num2str(app.DataOffsetRadialEditField.Value),"\n" ...
                    );
            end

            if strcmp(obj.dataType,'3Dp2roud')
                pars = strcat(pars, ...
                    "\n3D P2ROUD\n\n", ...
                    "trajectory = ",obj.trajectoryFileName,"\n" ...
                    );
            end

            if strcmp(obj.dataType,'2Dradial') || strcmp(obj.dataType,'2Dradialms')

                if app.HalfCircleRadialButton.Value == 1    trajType = 1; end
                if app.FullCircleRadialButton.Value == 1    trajType = 2; end
                if app.UserDefinedRadialButton.Value == 1   trajType = 3; end

                pars = strcat(pars,...
                    "\n2D radial\n\n", ...
                    "Gx delay = ",num2str(app.GxDelayEditField.Value),"\n",...
                    "Gy delay = ",num2str(app.GyDelayEditField.Value),"\n",...
                    "Gz delay = ",num2str(app.GzDelayEditField.Value),"\n",...
                    "trajectory = ",num2str(trajType)...
                    );
            end

            fid = fopen(strcat(obj.exportDir,filesep,'recoparameters_',app.tag,'.txt'),'wt');
            fprintf(fid,pars);
            fclose(fid);

        end % ExportRecoParametersFcn



        % ---------------------------------------------------------------------------------
        % Extract unsorted k-data from raw data
        % ---------------------------------------------------------------------------------
        function obj = extractData(obj)

            obj.raw = cell(obj.nr_coils); % Initialize cell array with raw k-space data

            switch obj.dataType

                case {'2D','2Dms','3D','3Dp2roud'} % 2D and 3D data
                    for i = 1:obj.nr_coils
                        % Cut off the navigator and spacer
                        obj.raw{i} = obj.data{i}(:,:,:,obj.primaryNavigatorPoint+obj.nrNavPointsDiscarded+1:end);
                    end

                case {'2Dradial','2Dradialms','3Dute'} % 2D radial and 3D UTE data
                    for i = 1:obj.nr_coils
                        obj.raw{i} = obj.data{i};
                    end

            end

        end % extractData



        % ---------------------------------------------------------------------------------
        % Assign bin times, relative binning
        % ---------------------------------------------------------------------------------
        function obj = assignBinTimes(obj, app)

            cardLocations = obj.heartTrigPoints;     % Cardiac trigger points in units of samples
            respLocations = obj.respTrigPoints;      % Respiratory trigger points in units of samples
            nrCardFrames = app.nrCardFrames;            % Number of cardiac frames
            nrRespFrames = app.nrRespFrames;            % Number of respiratory frames

            nrCard = length(cardLocations);             % Number of cardiac time points
            nrResp = length(respLocations);             % Number of respiratory time points

            % Cardiac binning
            %
            %  INPUT: card_locations = array with the cardiac trigger points in units of samples
            %         nr_card_frames = number of desired cardiac bins
            %
            %  card_locations(i)                                card_locations(i+1)
            %       ||  j=1     j=2     j=3     j=4     j=5     j=nr_of_card_frames
            %       ||       |       |       |       |       |       ||
            %       ||       |       |       |       |       |       ||
            %       ||       |       |       |       |       |       ||
            %       cnt=1   cnt=2   cnt=3   cnt=4   cnt=5   cnt=6   cnt=7   cnt=...
            %      cbins(1) cbins(2) .....
            %
            %  RESULT: array with time-stamp of all the cardiac bins for all heartbeats in the measurement in units of samples
            %

            cBins = zeros(1,(nrCard-1)*nrCardFrames); % Array with time-stamp of all the cardiac bins for all heartbeats in the measurement in units of samples

            cnt = 1;                        % Counter for the cardiac bins
            for i=1:nrCard-1                % Loop over all cardiac trigger points
                for j=1:nrCardFrames        % Loop over all cardiac frames
                    cBins(cnt)=cardLocations(i)+(j-1)*(cardLocations(i+1)-cardLocations(i))/nrCardFrames; % Time-stamp of the cardiac bin
                    cnt = cnt + 1;          % Increase counter
                end
            end

            % Respiratory binning

            rBins = zeros(1,(nrResp-1)*nrRespFrames); % Array with time-stamp of all the respiratory bins for all respiratory cycles in the measurement in units of samples

            cnt = 1;                        % Counter for the respiratory bins
            for i=1:nrResp-1                % Loop over all respiratory trigger points
                for j=1:nrRespFrames        % Loop over all respiratory frames
                    rBins(cnt)=respLocations(i)+(j-1)*(respLocations(i+1)-respLocations(i))/nrRespFrames; % Time-stamp of the respiratory bin
                    cnt = cnt + 1;          % Increase counter
                end
            end

            obj.cardBins = cBins;     % Time-stamp of the cardiac bins for all heartbeats in the measurement in units of samples
            obj.respBins = rBins;     % Time-stamp of the respiratory bins for all respiratory cycles in the measurement in units of samples

        end % assignBinTimes



        % ---------------------------------------------------------------------------------
        % Assign bin times, absolute binning
        % ---------------------------------------------------------------------------------
        function obj = assignBinTimesAbs(obj, app)

            cardLocations = obj.heartTrigPoints;
            respLocations = obj.respTrigPoints;
            heartRate = obj.meanHeartRate;
            repetitionTime = obj.TR;
            nrCardFrames = app.nrCardFrames;
            nrRespFrames = app.nrRespFrames;

            nrCard = length(cardLocations); % Number of cardiac time points
            nrResp = length(respLocations); % Number of respiratory time points

            %  Cardiac binning = ABSOLUTE BINNING
            %
            %  INPUT: card_locations = array with the cardiac trigger points in units of samples
            %         nr_card_frames = number of desired cardiac bins
            %
            %  card_locations(i)                                card_locations(i+1)
            %       ||  j=1     j=2     j=3     j=4     j=5     j=nr_of_card_frames
            %       ||       |       |       |       |       |       ||
            %       ||       |       |       |       |       |       ||
            %       ||       |       |       |       |       |       ||
            %       cnt=1   cnt=2   cnt=3   cnt=4   cnt=5   cnt=6   cnt=7   cnt=...
            %      cbins(1) cbins(2) .....
            %
            %  RESULT: cbins = array with time-stamp of all the cardiac bins for all heartbeats in the measurement in units of samples
            %

            % Frame duration in sample tijd
            frameDuration = (60/heartRate)*(1000/repetitionTime)/nrCardFrames;

            cBins = zeros(2,(nrCard-1)*nrCardFrames);

            cnt = 1;
            for i=1:nrCard-1

                nrframes = floor((cardLocations(i+1)-cardLocations(i))/frameDuration); % number of frames for this particular heartbeat
                if nrframes>nrCardFrames
                    nrframes = nrCardFrames;
                end

                for j=1:nrframes
                    cBins(1,cnt) = cardLocations(i)+(j-1)*frameDuration;    % divide heart-beat in nrframes equal bins of duration framedur
                    cBins(2,cnt) = j;                                       % immediately assign frame number
                    cnt = cnt + 1;
                end

            end

            cBins = cBins(:,1:cnt-1);

            % Respiratory binning = regular relative binning

            rBins = zeros(1,(nrResp-1)*nrRespFrames);

            cnt = 1;
            for i=1:nrResp-1
                for j=1:nrRespFrames
                    rBins(cnt)=respLocations(i)+(j-1)*(respLocations(i+1)-respLocations(i))/nrRespFrames;
                    cnt = cnt + 1;
                end
            end

            obj.cardBins = cBins;
            obj.respBins = rBins;

        end % assignBinTimesAbs



        % ---------------------------------------------------------------------------------
        % Assign bin frames, relative binning
        % ---------------------------------------------------------------------------------
        function obj = assignBinFrames(obj,  app)

            % Assigns all the measured k-lines to a specific cardiac phase and respiratory phase bin

            binTimesCard = obj.cardBins';
            binTimesResp = obj.respBins';
            respiratoryWindow = obj.respWindow;
            nrOfKlines = obj.nrKlines;
            nrCardFrames = app.nrCardFrames;
            nrRespFrames = app.nrRespFrames;

            % Assignmnent of all k-line acquisitions to a specific bin
            %
            % INPUT : bin_times_card = array of all time-stamps (units of samples) of the whole acquisition
            %
            %                       >    <--- decreases ------ bin_times_card(j)
            %                      ||            |
            %                      ||  j-th bin  |
            %                      ||            |
            %                      ||            |
            %       - increases -----> i-th time point / sample
            %
            % OUTPUT: For each time point / k-line sample assignment to cardiac bin (saw-tooth pattern)
            %
            %         /|      /|    /|        /|
            %        / |     / |   / |      /  |
            %       /  |   /   |  /  |    /    |  /     = PHASE BINNING, each heart-beat has equal number of bins
            %      /   |  /    | /   |  /      | /                       despite differences in heart beat duration
            %     /    |/      |/    |/        |/
            %     12345 1 23 45 12345 1 2 3 4 5  <- bins
            %

            startPoint = round(binTimesResp(2)+1);                          % First measurement to be considered
            endPoint = round(binTimesResp(length(binTimesResp))-1);         % Last measurement to be considered
            locMaxCard = length(binTimesCard);                              % Last bin border
            locMaxResp = length(binTimesResp);

            cardAssignments = zeros(nrOfKlines,1);                            % Zero = no assignment (breathing, begin/end of data)
            parfor i=startPoint:endPoint                                    % Start search to which heartbeat the measurement belongs
                j=locMaxCard;
                while j>1 && i<binTimesCard(j)  %#ok<*PFBNS>
                    j=j-1;
                end
                cardAssignments(i) = mod(j-1,nrCardFrames)+1;            % Assign to bin frame number = j modulus nr_frames
                if nrRespFrames==1 && respiratoryWindow(i)==1                   % If measurement is during respiration and only 1 resp state, put back to 0 to discard this k-line
                    cardAssignments(i) = 0;
                end
            end

            respAssignments = zeros(nrOfKlines,1);                        % Zero = no assignment (breathing, begin/end of data)
            parfor i=startPoint:endPoint                                % Start search to which respiration the measurement belongs
                j=locMaxResp;
                while j>1 && i<binTimesResp(j)
                    j=j-1;
                end
                respAssignments(i) = mod(j-1,nrRespFrames)+1;        % Assign to bin frame number = j modulus nr_frames
            end

            % Report back
            obj.cardBinNrs = cardAssignments;
            obj.respBinNrs = respAssignments;

        end % assignBinFrames



        % ---------------------------------------------------------------------------------
        % Assign bin frames, absolute binning
        % ---------------------------------------------------------------------------------
        function obj = assignBinFramesAbs(obj,  app)

            binTimesCard = obj.cardBins';
            binTimesResp = obj.respBins';
            respiratoryWindow = obj.respWindow;
            nrOfKlines = obj.nrKlines;
            nrRespFrames = app.nrRespFrames;

            % Assigns all the measured k-lines to a specific cardiac phase and respiratory phase bin
            %
            % ABSOLUTE BINNING

            startPoint = round(binTimesResp(2)+1);                          % First measurement to be considered
            endPoint = round(binTimesResp(length(binTimesResp))-1);         % Last measurement to be considered
            locMaxCard = length(binTimesCard);                              % Last bin border
            locMaxResp = length(binTimesResp);

            % Assignmnent of all k-line acquisitions to a specific bin
            %
            % INPUT : bin_times_card = array of all time-stamps (units of samples) of the whole acquisition
            %
            %                       >    <--- decreases ------ bin_times_card(j)
            %                      ||            |
            %                      ||  j-th bin  |
            %                      ||            |
            %                      ||            |
            %       - increases -----> i-th time point / sample
            %
            % OUTPUT: For each time point / k-line sample assignment to cardiac bin (saw-tooth pattern)
            %
            %
            %
            %         /|      /|    /|        /|
            %        / |     / |   / |      /  |
            %       /  |   /   |  /  |    /    |  /     = ABSOLUTE BINNING, time between frames is equal
            %      /   |  /    | /   |  /      | /                          despite differences in heart beat duration
            %     /    |/      |/    |/        |/                           -> later frames receive less data
            %     12345 1234567 12345 123456789  <- frames/bins
            %

            cardAssignments = zeros(nrOfKlines,1);                            % Zero = no assignment (breathing, begin/end of data)
            parfor i=startPoint:endPoint                                    % Start search to which heartbeat the measurement belongs
                j=locMaxCard;

                while j>1 && i<binTimesCard(j,1)
                    j=j-1;
                end

                cardAssignments(i) = binTimesCard(j,2);                     % Assign to bin frame number

                if nrRespFrames==1 && respiratoryWindow(i)==1
                    cardAssignments(i) = 0;                                 % If measurement is during respiration and only 1 resp state, put back to 0 to discard this k-line
                end

            end

            % Respiratory binning is still done using relative binning

            respAssignments = zeros(nrOfKlines,1);                            % Zero = no assignment (breathing, begin/end of data)
            parfor i=startPoint:endPoint                                    % Start search to which respiration the measurement belongs
                j=locMaxResp;
                while j>1 && i<binTimesResp(j)
                    j=j-1;
                end
                respAssignments(i) = mod(j-1,nrRespFrames)+1;               % Assign to bin frame number = j modulus nr_frames
            end

            obj.cardBinNrs = cardAssignments;
            obj.respBinNrs = respAssignments;

        end % assignBinFramesAbs



        % ---------------------------------------------------------------------------------
        % Fill K-space 2D
        % ---------------------------------------------------------------------------------
        function obj = fillKspace2D(obj, app)

            % This function creates 2 arrays
            % (1) The kspace data sorted into the correct cardiac frames and phase-encoding positions
            % (2) An array with the same size that keeps track of the number of averages per k-space point for normalization/statistics purposes

            % Required input:
            %
            % rawdata               = unsorted k-space data
            % navheart              = navigator signal, used to construct an average heartbeat
            % nrCardFrames          = number of desired cardiac frames, must be consisted with bin assignments (bin_ass)
            % nrRespFrames          = number of desired respiratory frames
            % nrDynamics            = number of desired dynamics
            % dimz                  = number of slices
            % dimy                  = dimensions of the images: dimy (phase encoding)
            % nrKsteps              = number of k-lines in 1 repetition
            % dimx                  = dimensions of the images: dimx (readout)
            % cardBinAss            = the cardiac bin assignment array for all measured k-lines
            % respBinAss            = the respiratory bin assignment array for all measured k-lines
            % trajectory            = the k-space trajectory
            % includeWindow         = data which should be include: 1 = yes, 0 = no
            % share                 = number of k-space sharing points

            share = app.SharingEditField.Value;
            obj.kSpace = cell(obj.nr_coils);
            obj.kSpaceAvg = [];
            includeAcqWindow = obj.includeWindow.*obj.excludeWindow;
            cardBinAss = obj.cardBinNrs;
            respBinAss = obj.respBinNrs;
            dimX = obj.dimx;
            dimY = obj.dimy;
            dimZ = obj.dimz;
            nrOfKsteps = obj.nrKsteps;
            traj = obj.trajectory;
            nrCardFrames = app.nrCardFrames;
            nrRespFrames = app.nrRespFrames;
            nrDynamics = app.nrDynamics;

            for coilnr = 1:obj.nr_coils

                rawdata = obj.raw{coilnr};

                nrReps = size(rawdata,1);                                       % Number of k-space repetitions
                unsortedKspace = reshape(rawdata,[1,size(rawdata),1]);

                % Dynamics assignment
                totalk = nrReps * nrOfKsteps * dimZ;
                dynBinAss = round(linspace(0.5, nrDynamics+0.49, totalk));      % List of increasing integer number 1 .. nr_dynamics evenly spaced over the entire acquistion time

                % Sorting
                sortedKspace = complex(zeros(nrRespFrames,nrCardFrames,dimZ,dimY,dimX,nrDynamics));   % Fill temp k-space with zeros
                sortedAverages = zeros(nrRespFrames,nrCardFrames,dimZ,dimY,dimX,nrDynamics);          % Fill temp nr averages array with zeros
                cnt = 0;
                for slice=1:dimZ                    % Loop over slices
                    for i=1:nrReps                  % Loop through all repetitions
                        for j=1:nrOfKsteps            % Loop through all the phase-encoding steps
                            cnt = cnt + 1;
                            if (cardBinAss(cnt) > 0) && (includeAcqWindow(cnt) == 1)       % If assigment = 0, this acquisition is discarded
                                kline = traj(mod(cnt - 1,nrOfKsteps) + 1);                % The phase-encoding step using the trajectory info
                                sortedKspace(respBinAss(cnt),cardBinAss(cnt),slice,kline,:,dynBinAss(cnt)) = sortedKspace(respBinAss(cnt),cardBinAss(cnt),slice,kline,:,dynBinAss(cnt)) + unsortedKspace(1,i,slice,j,:,1);   % add the data to the correct k-position
                                sortedAverages(respBinAss(cnt),cardBinAss(cnt),slice,kline,:,dynBinAss(cnt)) = sortedAverages(respBinAss(cnt),cardBinAss(cnt),slice,kline,:,dynBinAss(cnt)) + 1;        % increase the number of averages with 1
                            end
                        end
                    end
                end

                % Temp k-space
                tmpKspace = sortedKspace;
                tmpAverages = sortedAverages;

                % Find center of k-space
                kSpaceSum = squeeze(sum(sortedKspace,[1 2 3 6]));                   % Sum over all slices frames and dynamics
                [row, col] = find(ismember(kSpaceSum, max(kSpaceSum(:))));          % Coordinate of k-space maximum = center of k-space

                % Weighted view sharing
                if (share > 0) && (nrCardFrames > 1 || nrRespFrames > 1)

                    % Respiratory or cardiac frames
                    nrFrames = nrCardFrames;
                    if nrRespFrames > 1
                        nrFrames = nrRespFrames;
                        tmpKspace = permute(tmpKspace,[2,1,3,4,5,6]);
                        tmpAverages = permute(tmpAverages,[2,1,3,4,5,6]);
                        sortedKspace = permute(sortedKspace,[2,1,3,4,5,6]);
                        sortedAverages = permute(sortedAverages,[2,1,3,4,5,6]);
                    end

                    % determine share range
                    maxShare = round(max([nrCardFrames nrRespFrames])/2);     % Maximum number of shares
                    share(share > maxShare) = maxShare;

                    % define ellipsoid regions
                    Ry = round(dimY/share/2);
                    Rx = round(dimX/share/2);
                    [Y,X] = ndgrid(1:dimY,1:dimX);
                    L = zeros(share,dimY,dimX);
                    for i = 1:share
                        L(i,:,:) = sqrt( ((row-Y)/(Ry*i)).^2 + ((col-X)/(Rx*i)).^2 ) <= 1;
                    end
                    C(1,:,:) = L(1,:,:);
                    if share > 1
                        for i = 2:share
                            C(i,:,:) = L(i,:,:) - L(i-1,:,:);
                        end
                    end

                    % Weights
                    weights = zeros(share);
                    for i = 1:share
                        for j = 1:share
                            weights(i,j) = retro.gauss(i+j-1,share,0);
                        end
                    end
                    weights = 0.5*weights/max(weights(:));

                    % Apply sharing to k-space
                    for frame = 1:nrFrames

                        for i = -share:share

                            sharedFrame = frame + i;
                            sharedFrame(sharedFrame < 1) = nrFrames - sharedFrame - 1;
                            sharedFrame(sharedFrame > nrFrames) = sharedFrame - nrFrames;

                            if i~=0

                                for j = 1:share

                                    ROI = reshape(squeeze(C(j,:,:)),[1 1 1 dimY dimX 1])*weights(j,abs(i));
                                    tmpKspace(:,frame,:,:,:,:)   = tmpKspace(:,frame,:,:,:,:)   + sortedKspace(:,sharedFrame,:,:,:,:)   .* ROI;
                                    tmpAverages(:,frame,:,:,:,:) = tmpAverages(:,frame,:,:,:,:) + sortedAverages(:,sharedFrame,:,:,:,:) .* ROI;

                                end

                            end

                        end

                    end

                    % Respiratory or cardiac frames
                    if nrRespFrames > 1
                        tmpKspace = permute(tmpKspace,[2,1,3,4,5,6]);
                        tmpAverages = permute(tmpAverages,[2,1,3,4,5,6]);
                    end

                end

                % Normalize by number of averages
                tmpKspace = tmpKspace./tmpAverages;
                tmpKspace(isnan(tmpKspace)) = complex(0); % correct for NaN or Inf because of division by zero in case of missing k-lines
                tmpKspace(isinf(tmpKspace)) = complex(0);

                % Apply a circular Tukey filter
                filterWidth = 0.1;
                tukeyFilter(1,1,1,:,:,1) = retro.circtukey2D(dimY,dimX,row,col,filterWidth);
                tmpKspace = tmpKspace.*tukeyFilter;

                % Report back k-space per coil
                obj.kSpace{coilnr} = tmpKspace;

            end

            % Report back averages
            obj.kSpaceAvg = tmpAverages;

        end % fillKspace2D



        % ---------------------------------------------------------------------------------
        % Fill K-space 3D
        % ---------------------------------------------------------------------------------
        function obj = fillKspace3D(obj, app)

            % This function creates 2 arrays
            % (1) The 3D kspace data sorted into the correct cardiac frames and phase-encoding positions
            % (2) An array with the same size that keeps track of the number of averages per k-space point for normalization/statistics purposes

            % Required input:
            %
            % raw                   = unsorted k-space data
            % repstart              = first k-space repetition that is considered, used for partial reconstruction of the data
            % navheart              = navigator signal, used to construct an average heartbeat
            % nr_of_card_frames     = number of desired cardiac frames, must be consisted with bin assignments (bin_ass)
            % dimz                  = 2nd phase encoding dimension
            % dimy                  = dimensions of the images: dimy (phase encoding)
            % nrKsteps              = number of k-lines in 1 repetition
            % dimx                  = dimensions of the images: dimx (readout)
            % card_bin_ass          = the cardiac bin assignment array for all measured k-lines
            % resp_bin_ass          = the respiratory bin assignment array for all measured k-lines
            % traj                  = the k-space trajectory
            % includewindow         = data which should be include: 1 = yes, 0 = no
            % share                 = number of weighted view sharing points

            share = app.SharingEditField.Value;
            obj.kSpace = cell(obj.nr_coils);
            obj.kSpaceAvg = [];
            includeAcqWindow = obj.includeWindow.*obj.excludeWindow;
            cardBinAss = obj.cardBinNrs;
            respBinAss = obj.respBinNrs;
            dimX = obj.dimx;
            dimY = obj.dimy;
            dimZ = obj.dimz;
            nrOfKsteps = obj.nrKsteps;
            traj = obj.trajectory;
            nrCardFrames = app.nrCardFrames;
            nrRespFrames = app.nrRespFrames;
            nrDynamics = app.nrDynamics;

            for coilnr = 1:obj.nr_coils

                rawdata = obj.raw{coilnr};

                sortedKspace = complex(zeros(nrRespFrames,nrCardFrames,dimZ,dimY,dimX,nrDynamics));     % Fill temp k-space with zeros
                sortedAverages = zeros(nrRespFrames,nrCardFrames,dimZ,dimY,dimX,nrDynamics);            % Fill temp nr averages array with zeros
                nrReps = size(rawdata,1);                                                               % Number of k-space repetitions
                unsortedKspace = reshape(rawdata,[1,size(rawdata),1]);

                % Dynamics assignment
                totalk = nrReps * nrOfKsteps * dimZ;
                dynBinAss = round(linspace(0.5, nrDynamics+0.49, totalk));       % List of increasing integer number 1 .. nr_dynamics evenly spaced over the entire acquistion time

                % Adapt trajectory for 3D acqusition

                % The y-dimension
                traj3Dy = zeros(nrReps * nrOfKsteps * dimZ,1);
                cnt = 1;
                for i = 1:nrReps
                    for j = 1:nrOfKsteps
                        for k = 1:dimZ
                            traj3Dy(cnt) = traj(j);
                            cnt = cnt + 1;
                        end
                    end
                end

                % The z-dimension
                traj3Dz = zeros(nrReps * nrOfKsteps * dimZ,1);

                switch obj.pe2_centric_on

                    case 0

                        % Linear in the 3rd dimension
                        cnt = 1;
                        for i = 1:nrReps
                            for j = 1:nrOfKsteps
                                for k = 1:dimZ
                                    traj3Dz(cnt) = k;  % Linear
                                    cnt = cnt + 1;
                                end
                            end
                        end

                    case 1

                        % Centric in the 3rd dimension
                        cnt = 1;
                        cf = retro.centricFilling(dimZ);
                        for i = 1:nrReps
                            for j = 1:nrOfKsteps
                                for k = 1:dimZ
                                    traj3Dz(cnt) = cf(k);  % Centric
                                    cnt = cnt + 1;
                                end
                            end
                        end

                    case 2

                        % Special case, trajectory
                        cnt = 1;
                        for i = 1:nrReps
                            for j = 1:nrOfKsteps
                                for k = 1:dimZ
                                    traj3Dz(cnt) = obj.pe2_traj(k) + round(dimZ/2) + 1;
                                    cnt = cnt + 1;
                                end
                            end
                        end

                end

                % Do the filling of k-space
                cnt = 0;

                for i = 1:nrReps                    % Loop through all repetitions

                    for j = 1:nrOfKsteps              % Loop through all the phase-encoding steps

                        for k = 1:dimZ              % Loop through phase-encoding 3rd dimension

                            cnt = cnt + 1;

                            if (cardBinAss(cnt) > 0) && (includeAcqWindow(cnt) == 1)     % If assigment == 0, this acquisition is discarded

                                kline_y = traj3Dy(cnt);            % The phase-encoding step using the 3D trajectory info
                                kline_z = traj3Dz(cnt);            % The 2nd phase-encoding

                                sortedKspace(respBinAss(cnt),cardBinAss(cnt),kline_z,kline_y,:,dynBinAss(cnt))   = sortedKspace(respBinAss(cnt),cardBinAss(cnt),kline_z,kline_y,:,dynBinAss(cnt))   + unsortedKspace(1,i,k,j,:,1);       % add the data to the correct k-position
                                sortedAverages(respBinAss(cnt),cardBinAss(cnt),kline_z,kline_y,:,dynBinAss(cnt)) = sortedAverages(respBinAss(cnt),cardBinAss(cnt),kline_z,kline_y,:,dynBinAss(cnt)) + 1;                                  % increase the number of averages with 1

                            end

                        end

                    end

                end

                % Temp new k-space copy
                tmpKspace = sortedKspace;
                tmpAverages = sortedAverages;

                % Find center of k-space
                kSpaceSum = squeeze(sum(sortedKspace,[1 2 6]));     % Sum over all slices frames and dynamics
                [~,idx] = max(kSpaceSum(:));
                [lev, row, col] = ind2sub(size(kSpaceSum),idx);     % Coordinate of k-space maximum = center of k-space

                % Weighted view sharing
                if (share > 0) && (nrCardFrames > 1 || nrRespFrames > 1)

                    % Respiratory or cardiac frames
                    nrFrames = nrCardFrames;
                    if nrRespFrames > 1
                        nrFrames = nrRespFrames;
                        tmpKspace = permute(tmpKspace,[2,1,3,4,5,6]);
                        tmpAverages = permute(tmpAverages,[2,1,3,4,5,6]);
                        sortedKspace = permute(sortedKspace,[2,1,3,4,5,6]);
                        sortedAverages = permute(sortedAverages,[2,1,3,4,5,6]);
                    end

                    % Determine share range
                    maxShare = round(max([nrCardFrames nrRespFrames])/2); % Maximum number of shares
                    share(share > maxShare) = maxShare;
                    weights = retro.gauss(1:share+1,share,0);
                    weights = weights/max(weights);

                    % Define ellipsoid regions
                    Rz = round(dimZ/share/2);
                    Ry = round(dimY/share/2);
                    Rx = round(dimX/share/2);
                    [Z,Y,X] = ndgrid(1:dimZ,1:dimY,1:dimX);
                    L = zeros(share,dimZ,dimY,dimX);
                    for i = 1:share
                        L(i,:,:,:) = sqrt( ((lev-Z)/(Rz*i)).^2 + ((row-Y)/(Ry*i)).^2 + ((col-X)/(Rx*i)).^2 ) <= 1;
                    end
                    C(1,:,:,:) = L(1,:,:,:);
                    if share > 1
                        for i = 2:share
                            C(i,:,:,:) = L(i,:,:,:) - L(i-1,:,:,:);
                        end
                    end

                    % Weights
                    for i = 1:share
                        for j = 1:share
                            weights(i,j) = retro.gauss(i+j-1,share,0);
                        end
                    end
                    weights = 0.5*weights/max(weights(:));

                    % Apply sharing to k-space
                    for frame = 1:nrFrames

                        for i = -share:share

                            sharedframe = frame + i;
                            sharedframe(sharedframe < 1) = nrFrames - sharedframe - 1;
                            sharedframe(sharedframe > nrFrames) = sharedframe - nrFrames;

                            if i~=0

                                for j = 1:share

                                    ROI = reshape(squeeze(C(j,:,:,:)),[1 1 dimZ dimY dimX 1])*weights(j,abs(i));
                                    tmpKspace(:,frame,:,:,:,:)   = tmpKspace(:,frame,:,:,:,:)   + sortedKspace(:,sharedframe,:,:,:,:)   .* ROI;
                                    tmpAverages(:,frame,:,:,:,:) = tmpAverages(:,frame,:,:,:,:) + sortedAverages(:,sharedframe,:,:,:,:) .* ROI;

                                end

                            end

                        end

                    end

                    % Respiratory of cardiac frames
                    if nrRespFrames > 1
                        tmpKspace = permute(tmpKspace,[2,1,3,4,5,6]);
                        tmpAverages = permute(tmpAverages,[2,1,3,4,5,6]);
                    end

                end

                % Normalize by number of averages
                tmpKspace = tmpKspace./tmpAverages;
                tmpKspace(isnan(tmpKspace)) = complex(0);     % Correct for NaN because of division by zero in case of missing k-lines
                tmpKspace(isinf(tmpKspace)) = complex(0);

                % Apply a circular Tukey filter
                filterWidth = 0.1;
                tukeyFilter(1,1,:,:,:,1) = retro.circtukey3D(dimZ,dimY,dimX,lev,row,col,filterWidth);
                tmpKspace = tmpKspace.*tukeyFilter;

                % Report back
                obj.kSpace{coilnr} = tmpKspace;

            end

            obj.kSpaceAvg = tmpAverages;

        end % fillKspace3D



        % ---------------------------------------------------------------------------------
        % Fill K-space 2D real-time
        % ---------------------------------------------------------------------------------
        function obj = fillKspace2Drealtime(obj, app)

            % This function creates 3 arrays
            % (1) The kspace data sorted into the correct cardiac frames,  and phase-encoding positions
            % (2) An array with the same size that keeps track of the number of averages per k-space point for normalization/statistics purposes
            % (3) The number of dynamics = number of heart-beats in the datasets

            % Required input:
            %
            % rawdata               = unsorted k-space data
            % navheart              = navigator signal, used to construct an average heartbeat
            % nr_of_card_frames     = number of desired cardiac frames, must be consisted with bin assignments (bin_ass)
            % nr_of_resp_frames     = number of desired respiratory frames
            % dimz                  = number of slices
            % dimy                  = dimensions of the images: dimy (phase encoding)
            % nrKsteps              = number of k-lines in 1 repetition
            % dimx                  = dimensions of the images: dimx (readout)
            % card_bin_ass          = the cardiac bin assignment array for all measured k-lines
            % traj                  = the k-space trajectory
            % share                 = weighted view sharing of neighboring data

            share = app.SharingEditField.Value;
            obj.kSpace = cell(obj.nr_coils);
            obj.kSpaceAvg = [];
            binTimesCard = obj.heartTrigPoints;
            dimX = obj.dimx;
            dimY = obj.dimy;
            dimZ = obj.dimz;
            nrOfKsteps = obj.nrKsteps;
            traj = obj.trajectory;
            nrCardFrames = app.nrCardFrames;
            nrRespFrames = app.nrRespFrames;

            for coilnr = 1:obj.nr_coils

                rawData = obj.raw{coilnr};

                nrDynamics = round(length(binTimesCard)/nrCardFrames);                                     % number of dynamics equals the number of heartbeats in the acquisition
                sortedKspace = complex(zeros(nrRespFrames,nrCardFrames,dimZ,dimY,dimX,nrDynamics));       % fill temp k-space with zeros
                sortedAverages = zeros(nrRespFrames,nrCardFrames,dimZ,dimY,dimX,nrDynamics);              % fill temp averages array with zeros
                nrReps = size(rawData,1);                                                                          % number of k-space repetitions
                nrOfKlines = nrReps * nrOfKsteps * dimZ;                                                         % total number of k-lines

                % Unsorted k-lines
                rawData = permute(rawData,[3,1,2,4]);
                unsortedKlines = reshape(rawData,[1,1,1,nrOfKlines,dimX]);

                % Fill k-space
                fcnt = 1;

                for i = 1 : nrDynamics

                    for j = 1 : nrCardFrames

                        for k = 1:nrOfKlines

                            if fcnt < length(binTimesCard)

                                if k > binTimesCard(1,fcnt) && k < binTimesCard(1,fcnt+1)

                                    % The phase-encoding step using the trajectory info
                                    kline = traj(mod(k - 1,nrOfKsteps) + 1);

                                    % Fill k-line
                                    sortedKspace(1,j,1,kline,:,i) = sortedKspace(1,j,1,kline,:,i) + unsortedKlines(1,1,1,k,:);   % add the data to the correct k-position
                                    sortedAverages(1,j,1,kline,:,i) = sortedAverages(1,j,1,kline,:,i) + 1;

                                end

                            end

                        end

                        fcnt = fcnt + 1;

                    end

                end

                % Find center of k-space
                kSpaceSum = squeeze(sum(sortedKspace,[1 2 3 6]));               % Ssum over all slices frames and dynamics
                [row, col] = find(ismember(kSpaceSum, max(kSpaceSum(:))));      % Coordinate of k-space maximum = center of k-space

                % Temp k-space
                tmpKspace = sortedKspace;
                tmpAverages = sortedAverages;

                % Weighted view sharing
                if (share > 0) && (nrDynamics > 1)

                    % Determine share range
                    maxShare = 20;    % Maximum number of shares
                    share(share > maxShare) = maxShare;

                    % Define ellipsoid regions
                    Ry = round(dimY/share/2);
                    Rx = round(dimX/share/2);
                    [Y,X] = ndgrid(1:dimY,1:dimX);
                    L = zeros(share,dimY,dimX);
                    for i = 1:share
                        L(i,:,:) = sqrt( ((row-Y)/(Ry*i)).^2 + ((col-X)/(Rx*i)).^2 ) <= 1;
                    end
                    C(1,:,:) = L(1,:,:);
                    if share > 1
                        for i = 2:share
                            C(i,:,:) = L(i,:,:) - L(i-1,:,:);
                        end
                    end

                    % Weights
                    weights = zeros(share);
                    for i = 1:share
                        for j = 1:share
                            weights(i,j) = retro.gauss(i+j-1,share,0);
                        end
                    end
                    weights = 0.5*weights/max(weights(:));

                    % Apply sharing to k-space
                    for frame = 1:nrDynamics

                        for i = -share:share

                            sharedframe = frame + i;

                            if sharedframe > 0 && sharedframe < nrDynamics

                                if i~=0

                                    for j = 1:share

                                        ROI = reshape(squeeze(C(j,:,:)),[1 1 1 dimY dimX 1])*weights(j,abs(i));
                                        tmpKspace(:,:,:,:,:,frame)   = tmpKspace(:,:,:,:,:,frame)   + sortedKspace(:,:,:,:,:,sharedframe)   .* ROI;
                                        tmpAverages(:,:,:,:,:,frame) = tmpAverages(:,:,:,:,:,frame) + sortedAverages(:,:,:,:,:,sharedframe) .* ROI;

                                    end

                                end

                            end

                        end

                    end

                end

                % Normalize by number of averages
                tmpKspace = tmpKspace./tmpAverages;
                tmpKspace(isnan(tmpKspace)) = complex(0);
                tmpKspace(isinf(tmpKspace)) = complex(0);

                % Apply a circular Tukey filter
                filterWidth = 0.1;
                tukeyFilter(1,1,1,:,:,1) = retro.circtukey2D(dimY,dimX,row,col,filterWidth);
                tmpKspace = tmpKspace.*tukeyFilter;

                % Report back
                obj.kSpace{coilnr} = tmpKspace;

            end

            obj.kSpaceAvg = tmpAverages;
            app.nrDynamics = nrDynamics;

        end % fillKspace2Drealtime



        % ---------------------------------------------------------------------------------
        % Fill K-space 3D P2ROUD
        % ---------------------------------------------------------------------------------
        function obj = fillKspace3Dp2roud(obj, app)

            % This function creates 2 arrays
            % (1) the 3D kspace data sorted into the correct cardiac frames and phase-encoding positions
            % (2) an array with the same size that keeps track of the number of averages per k-space point for normalization/statistics purposes

            % Required input:
            %
            % raw                   = unsorted k-space data
            % repstart              = first k-space repetition that is considered, used for partial reconstruction of the data
            % navheart              = navigator signal, used to construct an average heartbeat
            % nr_of_card_frames     = number of desired cardiac frames, must be consisted with bin assignments (bin_ass)
            % dimz                  = 2nd phase encoding dimension
            % dimy                  = dimensions of the images: dimy (phase encoding)
            % nrKsteps              = number of k-lines in 1 repetition
            % dimx                  = dimensions of the images: dimx (readout)
            % card_bin_ass          = the cardiac bin assignment array for all measured k-lines
            % resp_bin_ass          = the respiratory bin assignment array for all measured k-lines
            % traj                  = the k-space trajectory
            % includewindow         = data which should be include: 1 = yes, 0 = no

            share = app.SharingEditField.Value;
            obj.kSpace = cell(obj.nr_coils);
            obj.kSpaceAvg = [];
            includeAcqWindow = obj.includeWindow.*obj.excludeWindow;
            cardBinAss = obj.cardBinNrs;
            respBinAss = obj.respBinNrs;
            dimX = obj.dimx;
            dimY = obj.dimy;
            dimZ = obj.dimz;
            nrOfKsteps = obj.nrKsteps;
            traj = obj.trajectory;
            nrCardFrames = app.nrCardFrames;
            nrRespFrames = app.nrRespFrames;
            nrDynamics = app.nrDynamics;

            for coilnr = 1:obj.nr_coils

                rawData = obj.raw{coilnr};

                sortedKspace = complex(zeros(nrRespFrames,nrCardFrames,dimZ,dimY,dimX,nrDynamics));       % fill temp k-space with zeros
                sortedAverages = zeros(nrRespFrames,nrCardFrames,dimZ,dimY,dimX,nrDynamics);              % fill temp nr averages array with zeros
                nrReps = size(rawData,1);                                                                      % number of k-space repetitions
                unsortedKspace = reshape(rawData,[1,size(rawData),1]);

                % Dynamics assignment
                totalk = nrReps * nrOfKsteps * dimZ;
                dynBinAss = round(linspace(0.5, nrDynamics+0.49, totalk));       % List of increasing integer number 1 .. nr_dynamics evenly spaced over the entire acquistion time

                % Trajectory for 3D p2roud acquisition
                offset1 = floor(dimY/2 + 1);
                offset2 = floor(dimZ/2 + 1);
                cnt1 = 1;
                for i = 1:nrReps
                    cnt2 = 1;
                    for j = 1:nrOfKsteps*dimZ
                        traj3Dy(cnt1) = traj(cnt2)   + offset1; %#ok<*AGROW>
                        traj3Dz(cnt1) = traj(cnt2+1) + offset2;
                        cnt1 = cnt1 + 1;
                        cnt2 = cnt2 + 2;
                    end
                end

                % Do the filling of k-space
                cnt = 0;

                for i = 1:nrReps                  % Loop through all repetitions

                    for j = 1:nrOfKsteps            % Loop through all the phase-encoding steps

                        for k = 1:dimZ            % Loop through phase-encoding 3rd dimension

                            cnt = cnt + 1;

                            if (cardBinAss(cnt) > 0) && (includeAcqWindow(cnt) == 1)     % If assigment == 0, this acquisition is discarded

                                kLineY = traj3Dy(cnt);            % The phase-encoding step using the 3D trajectory info
                                kLineZ = traj3Dz(cnt);            % The 2nd phase-encoding

                                sortedKspace(respBinAss(cnt),cardBinAss(cnt),kLineZ,kLineY,:,dynBinAss(cnt)) = sortedKspace(respBinAss(cnt),cardBinAss(cnt),kLineZ,kLineY,:,dynBinAss(cnt)) + unsortedKspace(1,i,k,j,:,1);     % add the data to the correct k-position
                                sortedAverages(respBinAss(cnt),cardBinAss(cnt),kLineZ,kLineY,:,dynBinAss(cnt)) = sortedAverages(respBinAss(cnt),cardBinAss(cnt),kLineZ,kLineY,:,dynBinAss(cnt)) + 1;                           % increase the number of averages with 1

                            end

                        end

                    end

                end

                % Temp new k-space copy
                tmpKspace = sortedKspace;
                tmpAverages = sortedAverages;

                % Find center of k-space
                kSpaceSum = squeeze(sum(sortedKspace,[1 2 6]));     % Sum over all slices frames and dynamics
                [~,idx] = max(kSpaceSum(:));
                [lev, row, col] = ind2sub(size(kSpaceSum),idx);     % Coordinate of k-space maximum = center of k-space

                % Weighted view sharing
                if (share > 0) && (nrCardFrames > 1 || nrRespFrames > 1)

                    % Respiratory or cardiac frames
                    nrFrames = nrCardFrames;
                    if nrRespFrames > 1
                        nrFrames = nrRespFrames;
                        tmpKspace = permute(tmpKspace,[2,1,3,4,5,6]);
                        tmpAverages = permute(tmpAverages,[2,1,3,4,5,6]);
                        sortedKspace = permute(sortedKspace,[2,1,3,4,5,6]);
                        sortedAverages = permute(sortedAverages,[2,1,3,4,5,6]);
                    end

                    % Determine share range
                    maxShare = round(max([nrCardFrames nrRespFrames])/2); % Maximum number of shares
                    share(share > maxShare) = maxShare;
                    weights = retro.gauss(1:share+1,share,0);
                    weights = weights/max(weights);

                    % Define ellipsoid regions
                    Rz = round(dimZ/share/2);
                    Ry = round(dimY/share/2);
                    Rx = round(dimX/share/2);
                    [Z,Y,X] = ndgrid(1:dimZ,1:dimY,1:dimX);
                    L = zeros(share,dimZ,dimY,dimX);
                    for i = 1:share
                        L(i,:,:,:) = sqrt( ((lev-Z)/(Rz*i)).^2 + ((row-Y)/(Ry*i)).^2 + ((col-X)/(Rx*i)).^2 ) <= 1;
                    end
                    C(1,:,:,:) = L(1,:,:,:);
                    if share > 1
                        for i = 2:share
                            C(i,:,:,:) = L(i,:,:,:) - L(i-1,:,:,:);
                        end
                    end

                    % Weights
                    for i = 1:share
                        for j = 1:share
                            weights(i,j) = retro.gauss(i+j-1,share,0);
                        end
                    end
                    weights = 0.5*weights/max(weights(:));

                    % Apply sharing to k-space
                    for frame = 1:nrFrames

                        for i = -share:share

                            sharedframe = frame + i;
                            sharedframe(sharedframe < 1) = nrFrames - sharedframe - 1;
                            sharedframe(sharedframe > nrFrames) = sharedframe - nrFrames;

                            if i~=0

                                for j = 1:share

                                    ROI = reshape(squeeze(C(j,:,:,:)),[1 1 dimZ dimY dimX 1])*weights(j,abs(i));
                                    tmpKspace(:,frame,:,:,:,:)   = tmpKspace(:,frame,:,:,:,:)   + sortedKspace(:,sharedframe,:,:,:,:)   .* ROI;
                                    tmpAverages(:,frame,:,:,:,:) = tmpAverages(:,frame,:,:,:,:) + sortedAverages(:,sharedframe,:,:,:,:) .* ROI;

                                end

                            end

                        end

                    end

                    % Respiratory of cardiac frames
                    if nrRespFrames > 1
                        tmpKspace = permute(tmpKspace,[2,1,3,4,5,6]);
                        tmpAverages = permute(tmpAverages,[2,1,3,4,5,6]);
                    end

                end

                % Normalize by number of averages
                tmpKspace = tmpKspace./tmpAverages;
                tmpKspace(isnan(tmpKspace)) = complex(0);     % Correct for NaN because of division by zero in case of missing k-lines
                tmpKspace(isinf(tmpKspace)) = complex(0);

                % Apply a circular Tukey filter
                filterWidth = 0.1;
                tukeyFilter(1,1,:,:,:,1) = retro.circtukey3D(dimZ,dimY,dimX,lev,row,col,filterWidth);
                tmpKspace = tmpKspace.*tukeyFilter;

                % Report back
                obj.kSpace{coilnr} = tmpKspace;

            end

            % Report back
            obj.kSpaceAvg = tmpAverages;

        end % fillKspace3Dp2roud




        % ---------------------------------------------------------------------------------
        % Fill k-space 2D RADIAL
        % Unique spokes (e.g. golden angle)
        % ---------------------------------------------------------------------------------
        function obj = fillKspace2Dradial(obj, app)

            obj.kSpace = cell(obj.nr_coils);
            obj.kSpaceTraj = [];
            obj.kSpaceAvg = [];

            includeAcqWindow = obj.includeWindow.*obj.excludeWindow;
            cardBinAss = obj.cardBinNrs;
            respBinAss = obj.respBinNrs;

            dimX = obj.dimx;
            dimZ = obj.dimz;

            nrCardFrames = app.nrCardFrames;
            nrRespFrames = app.nrRespFrames;
            nrDynamics = app.nrDynamics;
            interpolationFactor = 16;

            % Create spokes from trajectory list of angles
            spoke = zeros(1,1,obj.nrKlines,dimZ,dimX,1,3);

            for z = 1:dimZ
                for cnt = 1:obj.nrKlines
                    spoke(1,1,cnt,z,:,1,1) = (-floor(dimX/2)+0.5:floor(dimX/2)-0.5)*cos(obj.trajectory(cnt)*pi/180);
                    spoke(1,1,cnt,z,:,1,2) = (-floor(dimX/2)+0.5:floor(dimX/2)-0.5)*sin(obj.trajectory(cnt)*pi/180);
                end
            end

            for coilnr = 1:obj.nr_coils

                rawData = permute(obj.raw{coilnr},[3 2 1 4]); % needed to be able to rearrange according to spokes, readout
                sortedKspace = complex(zeros(nrRespFrames,nrCardFrames,obj.nrKlines,dimZ,dimX,nrDynamics));   % fill temp k-space with zeros
                sortedAverages = zeros(nrRespFrames,nrCardFrames,obj.nrKlines,dimZ,dimX,nrDynamics);
                sortedTraj = zeros(nrRespFrames,nrCardFrames,obj.nrKlines,dimZ,dimX,nrDynamics,3);
                unsortedKspace = reshape(rawData,[1,1,obj.nrKlines,1,dimX]);


                % Dynamics and slice assignment
                if dimZ > 1
                    dynBinAss = ones(obj.nrKlines);
                    sliceAss = round(linspace(0.5, dimZ+0.49, obj.nrKlines));

                else
                    dynBinAss = round(linspace(0.5, nrDynamics+0.49, obj.nrKlines));       % list of increasing integer number 1 .. nr_dynamics evenly spaced over the entire acquistion time
                    sliceAss = ones(obj.nrKlines);
                end

                % Loop over all acquired 2D spokes
                for cnt = 1:obj.nrKlines

                    if (cardBinAss(cnt) > 0) && (includeAcqWindow(cnt) == 1)

                        tmpKline1 = squeeze(unsortedKspace(1,1,cnt,1,:));

                        % Center echo
                        if app.CenterEchoCheckBox.Value
                            tmpKline2 = interp(tmpKline1,interpolationFactor);
                            [~,kCenter] = max(abs(tmpKline2));
                            kShift = floor(dimX/2)-kCenter/interpolationFactor;
                            tmpKline1 = retro.fracCircShift(tmpKline1,kShift);
                        end

                        % Phase correction for k-space center
                        if app.PhaseCorrectCheckBox.Value
                            kCenterPhase = mean(angle(tmpKline1((floor(dimX/2)+1)-1:(floor(dimX/2)+1)+1)));
                            tmpKline1 = tmpKline1.*exp(-1j.*kCenterPhase);
                        end

                        % Assign the k-line to respiratory, cardiac, slice, and dynamic bins
                        sortedKspace(respBinAss(cnt),cardBinAss(cnt),cnt,sliceAss(cnt),:,dynBinAss(cnt)) = tmpKline1;
                        sortedAverages(respBinAss(cnt),cardBinAss(cnt),cnt,sliceAss(cnt),:,dynBinAss(cnt)) = sortedAverages(respBinAss(cnt),cardBinAss(cnt),cnt,sliceAss(cnt),:,dynBinAss(cnt)) + 1;
                        sortedTraj(respBinAss(cnt),cardBinAss(cnt),cnt,sliceAss(cnt),:,dynBinAss(cnt),:) = spoke(1,1,cnt,1,:,1,:);

                    end

                end

                % Apply 1-D Tukey filter
                filterWidth = 0.25;
                tukeyWindow(1,1,1,1,:,1) = tukeywin(dimX,filterWidth);
                sortedKspace = sortedKspace.*tukeyWindow;

                % Report back k-space
                obj.kSpace{coilnr} = sortedKspace;

            end

            % Show typical radial filling for frame 1, dynamic 1, slice 1
            xc = squeeze(sortedTraj(1,1,:,1,:,1,1));
            yc = squeeze(sortedTraj(1,1,:,1,:,1,2));

            % Remove all zero entries
            ze = ~((xc(:,end)==0) & (yc(:,end)==0));
            xc = xc(ze,:);
            yc = yc(ze,:);

            % Restrict number of spokes to 1000
            trajSize = size(xc,1);
            skip = ceil(trajSize/1000);
            xc = xc(1:skip:end,[1 end])';
            yc = yc(1:skip:end,[1 end])';

            % Plot the initial trajectory
            app.PlotTrajectoryMovieFrameFcn(xc,yc,[],dimX,true);

            % Report back averages and trajectory
            obj.kSpaceTraj = sortedTraj;
            obj.kSpaceAvg = sortedAverages;

        end % fillKspace2Dradial




        % ---------------------------------------------------------------------------------
        % Fill k-space 2D RADIAL
        % Non-unique spokes
        % ---------------------------------------------------------------------------------
        function obj = fillKspace2DradialReg(obj, app)

            obj.kSpace = cell(obj.nr_coils);
            obj.kSpaceTraj = [];
            obj.kSpaceAvg = [];

            includeAcqWindow = obj.includeWindow.*obj.excludeWindow;
            cardBinAss = obj.cardBinNrs;
            respBinAss = obj.respBinNrs;

            dimX = obj.dimx;
            dimZ = obj.dimz;

            nrCardFrames = app.nrCardFrames;
            nrRespFrames = app.nrRespFrames;
            nrDynamics = app.nrDynamics;
            interpolationFactor = 16;

            % Create spokes from trajectory list of angles
            spoke = zeros(1,1,obj.nrKsteps,dimZ,dimX,1,3);
            srange = floor(dimX/2);

            for z = 1:dimZ
                for cnt = 1:obj.nrKsteps
                    spoke(1,1,cnt,z,:,1,1) = (-srange:srange-1)*cos(obj.trajectory(cnt)*pi/180);
                    spoke(1,1,cnt,z,:,1,2) = (-srange:srange-1)*sin(obj.trajectory(cnt)*pi/180);
                end
            end

            for coilnr = 1:obj.nr_coils

                rawData = permute(obj.raw{coilnr},[3 2 1 4]); % needed to be able to rearrange according to spokes, readout
                sortedKspace = complex(zeros(nrRespFrames,nrCardFrames,obj.nrKsteps,dimZ,dimX,nrDynamics));   % fill temp k-space with zeros
                sortedAverages = zeros(nrRespFrames,nrCardFrames,obj.nrKsteps,dimZ,dimX,nrDynamics);
                sortedTraj = zeros(nrRespFrames,nrCardFrames,obj.nrKsteps,dimZ,dimX,nrDynamics,3);
                unsortedKspace = reshape(rawData,[1,1,obj.nrKlines,1,dimX]);

                % Dynamics and slice assignment
                if dimZ > 1
                    dynBinAss = ones(obj.nrKlines);
                    sliceAss = round(linspace(0.5, dimZ+0.49, obj.nrKlines));
                else
                    dynBinAss = round(linspace(0.5, nrDynamics+0.49, obj.nrKlines));       % list of increasing integer number 1 .. nr_dynamics evenly spaced over the entire acquistion time
                    sliceAss = ones(obj.nrKlines);
                end

                % Spoke assignment
                spokeAss = repmat(1:obj.nrKsteps,[1 dimZ*obj.nr_repetitions]);

                % Loop over all acquired 2D spokes
                for cnt = 1:obj.nrKlines

                    if (cardBinAss(cnt) > 0) && (includeAcqWindow(cnt) == 1)

                        tmpKline1 = squeeze(unsortedKspace(1,1,cnt,1,:));

                        % Center echo
                        if app.CenterEchoCheckBox.Value
                            tmpKline2 = interp(tmpKline1,interpolationFactor);
                            [~,kCenter] = max(abs(tmpKline2));
                            kShift = floor(dimX/2)-kCenter/interpolationFactor;
                            tmpKline1 = retro.fracCircShift(tmpKline1,kShift);
                        end

                        % Phase correction for k-space center
                        if app.PhaseCorrectCheckBox.Value
                            kCenterPhase = mean(angle(tmpKline1((floor(dimX/2)+1)-1:(floor(dimX/2)+1)+1)));
                            tmpKline1 = tmpKline1.*exp(-1j.*kCenterPhase);
                            % kCenterAmp = mean(abs(tmpKline1((floor(dimx/2)+1)-1:(floor(dimx/2)+1)+1)));
                            % tmpKline1 = tmpKline1./kCenterAmp;
                        end

                        % Reshape for adding to sortedKspace
                        tmpKline1 = reshape(tmpKline1,[1 1 1 1 dimX 1]);

                        % Assign the k-line to respiratory, cardiac, slice, and dynamic bins
                        sortedKspace(respBinAss(cnt),cardBinAss(cnt),spokeAss(cnt),sliceAss(cnt),:,dynBinAss(cnt)) = sortedKspace(respBinAss(cnt),cardBinAss(cnt),spokeAss(cnt),sliceAss(cnt),:,dynBinAss(cnt)) + tmpKline1;
                        sortedAverages(respBinAss(cnt),cardBinAss(cnt),spokeAss(cnt),sliceAss(cnt),:,dynBinAss(cnt)) = sortedAverages(respBinAss(cnt),cardBinAss(cnt),spokeAss(cnt),sliceAss(cnt),:,dynBinAss(cnt)) + 1;
                        sortedTraj(respBinAss(cnt),cardBinAss(cnt),spokeAss(cnt),sliceAss(cnt),:,dynBinAss(cnt),:) = spoke(1,1,spokeAss(cnt),1,:,1,:);

                    end

                end

                % Normalize by number of averages
                sortedKspace = sortedKspace./sortedAverages;
                sortedKspace(isnan(sortedKspace)) = complex(0);
                sortedKspace(isinf(sortedKspace)) = complex(0);

                % Apply 1-D Tukey filter
                filterWidth = 0.25;
                tukeyWindow(1,1,1,1,:,1) = tukeywin(dimX,filterWidth);
                sortedKspace = sortedKspace.*tukeyWindow;

                % Report back k-space
                obj.kSpace{coilnr} = sortedKspace;

            end

            % Show typical radial filling for frame 1, dynamic 1, slice 1
            xc = squeeze(sortedTraj(1,1,:,1,:,1,1));
            yc = squeeze(sortedTraj(1,1,:,1,:,1,2));

            % Remove all zero entries
            ze = ~((xc(:,end)==0) & (yc(:,end)==0));
            xc = xc(ze,:);
            yc = yc(ze,:);

            % Restrict number of spokes to 1000
            trajSize = size(xc,1);
            skip = ceil(trajSize/1000);
            xc = xc(1:skip:end,[1 end])';
            yc = yc(1:skip:end,[1 end])';

            % Plot the initial trajectory
            app.PlotTrajectoryMovieFrameFcn(xc,yc,[],dimX,true);

            % Report back averages and trajectory
            obj.kSpaceTraj = sortedTraj;
            obj.kSpaceAvg = sortedAverages;

        end % fillKspace2DradialReg




        % ---------------------------------------------------------------------------------
        % Fill k-space 3D UTE
        % ---------------------------------------------------------------------------------
        function obj = fillKspace3Dute(obj, app)

            % This function creates 2 arrays
            % (1) the kspace data sorted into the correct cardiac frames and phase-encoding positions
            % (2) the corresponding trajectories

            obj.kSpace = cell(obj.nr_coils);
            obj.kSpaceTraj = [];
            obj.kSpaceAvg = [];
            includeAcqWindow = obj.includeWindow.*obj.excludeWindow;
            cardBinAss = obj.cardBinNrs;
            respBinAss = obj.respBinNrs;
            offset = app.DataOffsetRadialEditField.Value;
            dimX = length(obj.gradTrajectory);

            for coilnr = 1:obj.nr_coils

                % Unsorted data
                rawData = obj.raw{coilnr};

                % Dynamics assignment
                dynBinAss = round(linspace(0.5, app.nrDynamics+0.49, obj.nrKlines));       % list of increasing integer number 1 .. nr_dynamics evenly spaced over the entire acquistion time

                % K-space radial spokes
                for cnt = 1:dimX
                    spoke(1,1,:,1,cnt,1,1) = dimX*(obj.trajectory(1,:)/32767)*obj.gradTrajectory(cnt);
                    spoke(1,1,:,1,cnt,1,2) = dimX*(obj.trajectory(2,:)/32767)*obj.gradTrajectory(cnt);
                    spoke(1,1,:,1,cnt,1,3) = dimX*(obj.trajectory(3,:)/32767)*obj.gradTrajectory(cnt);
                end

                % Sort
                sortedKspace = complex(zeros(app.nrRespFrames,app.nrCardFrames,obj.nrKlines,1,dimX,app.nrDynamics));   % fill temp k-space with zeros
                sortedAverages = zeros(app.nrRespFrames,app.nrCardFrames,obj.nrKlines,1,dimX,app.nrDynamics);
                sortedTraj = zeros(app.nrRespFrames,app.nrCardFrames,obj.nrKlines,1,dimX,app.nrDynamics,3);            % fill temp trajectory with zeros
                unsortedKspace = reshape(rawData,[1,size(rawData),1]);

                % Check if offset does not lead to index larger than available
                offset(offset < 0) = 0;
                offset(offset > (size(unsortedKspace,5)-dimX)) = size(unsortedKspace,5)-dimX;
                app.DataOffsetRadialEditField.Value = offset;

                % Loop over acquired 3D spokes
                for cnt = 1:obj.nrKlines

                    if (cardBinAss(cnt) > 0) && (includeAcqWindow(cnt) == 1)

                        tmpKline = squeeze(unsortedKspace(1,cnt,1,1,1+offset:dimX+offset,1));

                        % Phase correction for k-space center
                        if app.PhaseCorrectCheckBox.Value
                            kCenterPhase = angle(tmpKline(1));
                            tmpKline = tmpKline.*exp(-1j.*kCenterPhase);
                        end

                        % Assign the k-line to respiratory, cardiac, slice, and dynamic bins
                        sortedKspace(respBinAss(cnt),cardBinAss(cnt),cnt,1,:,dynBinAss(cnt)) = tmpKline;
                        sortedAverages(respBinAss(cnt),cardBinAss(cnt),cnt,1,:,dynBinAss(cnt)) = sortedAverages(respBinAss(cnt),cardBinAss(cnt),cnt,1,:,dynBinAss(cnt)) + 1;
                        sortedTraj(respBinAss(cnt),cardBinAss(cnt),cnt,1,:,dynBinAss(cnt),:) = spoke(1,1,cnt,1,:,1,:);

                    end

                end

                % Apply 1-D Tukey filter
                filterWidth = 0.25/2;
                tmpFilter = tukeywin(dimX*2,filterWidth);
                tmpFilter = tmpFilter(dimX+1:dimX*2);
                tukeyWindow(1,1,1,1,:,1) = tmpFilter;
                sortedKspace = sortedKspace.*tukeyWindow;

                % Report back k-space per coil
                obj.kSpace{coilnr} = sortedKspace;

            end

            % Show typical radial filling for frame 1, dynamic 1, slice 1
            xc = squeeze(sortedTraj(1,1,:,1,:,1,1));
            yc = squeeze(sortedTraj(1,1,:,1,:,1,2));
            zc = squeeze(sortedTraj(1,1,:,1,:,1,3));

            % Remove all zero entries
            ze = ~((xc(:,end)==0) & (yc(:,end)==0) & (zc(:,end)==0));
            xc = xc(ze,:);
            yc = yc(ze,:);
            zc = zc(ze,:);

            % Restrict number of spokes to 1000
            trajSize = size(xc,1);
            skip = ceil(trajSize/1000);
            xc = xc(1:skip:end,[1 end])';
            yc = yc(1:skip:end,[1 end])';
            zc = zc(1:skip:end,[1 end])';

            % Plot the initial trajectory
            app.PlotTrajectoryMovieFrameFcn(xc,yc,zc,dimX,true);

            % Report back k-space trajectory and averages
            obj.kSpaceTraj = sortedTraj;
            obj.kSpaceAvg = sortedAverages;


        end % fillKspace3Dute



        % ---------------------------------------------------------------------------------
        % Reshape k-space
        % ---------------------------------------------------------------------------------
        function obj = reshapeKspace(obj)

            switch obj.dataType

                case {'2D','2Dms','3D','3Dp2roud'}

                    % Reshape to either cardiac or respiratory CINE

                    [s1,s2,s3,s4,s5,s6] = size(obj.kSpaceAvg);
                    s = max([s1,s2]);
                    nrCoils = length(obj.kSpace);
                    for i = 1:nrCoils
                        obj.kSpace{i} = reshape(obj.kSpace{i},[s,s3,s4,s5,s6]);
                        obj.kSpace{i} = permute(obj.kSpace{i},[1,4,3,2,5,6]);
                    end
                    obj.kSpaceAvg = reshape(obj.kSpaceAvg,[s,s3,s4,s5,s6]);
                    obj.kSpaceAvg = permute(obj.kSpaceAvg,[1,4,3,2,5,6]);

                    % kSpace = frames, X, Y, Z(/slices), dynamics

                case {'2Dradial','2Dradialms','3Dute'}

                    % Reshape to either cardiac or respiratory CINE

                    [s1,s2,s3,s4,s5,s6] = size(obj.kSpace{1});
                    s = max([s1,s2]);
                    nrCoils = length(obj.kSpace);
                    for i = 1:nrCoils
                        obj.kSpace{i} = reshape(obj.kSpace{i},[s,s3,s4,s5,s6]);
                        obj.kSpace{i} = permute(obj.kSpace{i},[1,4,2,3,5,6]);
                    end

                    obj.kSpaceTraj = reshape(obj.kSpaceTraj,[s,s3,s4,s5,s6,3]);
                    obj.kSpaceTraj = permute(obj.kSpaceTraj,[1,4,2,3,5,6]);
                    obj.kSpaceAvg = reshape(obj.kSpaceAvg,[s,s3,s4,s5,s6]);
                    obj.kSpaceAvg = permute(obj.kSpaceAvg,[1,4,2,3,5]);

                    % kSpace     = frames, X(readout), Y(spokes), Z(1 or slices), dynamics
                    % kSpaceTraj = frames, X(readout), Y(spokes), Z(1 or slices), dynamics, 3(trajectory coordinates)

                    % Set multi-slice and multi-dynamic flags for
                    % visualization of k-space trajectory
                    if size(obj.kSpace{1},4) > 1
                        obj.multiSliceFlag = true;
                    else
                        obj.multiSliceFlag = false;
                    end

                    if size(obj.kSpace{1},5) > 1
                        obj.multiDynamicFlag = true;
                    else
                        obj.multiDynamicFlag = false;
                    end

            end

        end % reshapeKspace



        % ---------------------------------------------------------------------------------
        % K-space statistics
        % ---------------------------------------------------------------------------------
        function [obj, app] = kSpaceStats(obj, app)

            switch app.r.dataType

                case {'2D','2Dms','3D','3Dp2roud'}

                    % Calculate effective number of averages
                    app.AveragesViewField.Value = mean2(obj.kSpaceAvg);
                    obj.NO_AVERAGES = round(app.AveragesViewField.Value);

                    % Filling
                    app.FillingViewField.Value = round(100*nnz(obj.kSpaceAvg)/numel(obj.kSpaceAvg));
                    if app.FillingViewField.Value < 20
                        app.SetStatus(1);
                        app.TextMessage('WARNING: Low k-space filling, decrease number of frames or dynamics ...');
                    end

                    % Warning in case of large dataset
                    if app.FramesEditField.Value*app.DynamicsEditField.Value > 2000
                        app.SetStatus(1);
                        app.TextMessage('WARNING: large dataset, reconstruction may take a very long time ...');
                    end

                case '3Dute'

                    % Kspace radial spokes
                    nrKpoints = length(obj.gradTrajectory);

                    % Kspace
                    sortedKspace = obj.kSpace{1};

                    % Number of averages
                    nrFullNyquist = (pi/2)*nrKpoints^2;
                    nrNonZeros = nnz(sortedKspace(:))/nrKpoints;
                    app.AveragesViewField.Value = (nrNonZeros/nrFullNyquist)/(app.nrCardFrames*app.nrRespFrames*app.nrDynamics);
                    obj.NO_AVERAGES = round(app.AveragesViewField.Value);

                    % Filling
                    filling = round(100*(nrNonZeros/nrFullNyquist)/(app.nrCardFrames*app.nrRespFrames*app.nrDynamics));
                    filling(filling>100) = 100;
                    app.FillingViewField.Value = filling;

                case {'2Dradial','2Dradialms'}

                    % Kspace radial spokes
                    nrKpoints = obj.dimx;

                    % Kspace
                    sortedKspace = obj.kSpace{1};

                    % Number of averages
                    nrFullNyquist = (pi/2)*nrKpoints/2;
                    nrNonZeros = nnz(sortedKspace(:))/nrKpoints;
                    app.AveragesViewField.Value = (nrNonZeros/nrFullNyquist)/(app.nrCardFrames*app.nrRespFrames*app.nrDynamics);
                    obj.NO_AVERAGES = round(app.AveragesViewField.Value);

                    % Filling
                    filling = round(100*(nrNonZeros/nrFullNyquist)/(app.nrCardFrames*app.nrRespFrames*app.nrDynamics));
                    filling(filling>100) = 100;
                    app.FillingViewField.Value = filling;

            end

        end % kSpaceStats



        % ---------------------------------------------------------------------------------
        % K-space trajectory
        % ---------------------------------------------------------------------------------
        function obj = getKspaceTrajectory(obj, app)

            % Trajectory (linear, zigzag, user-defined, radial)
            switch obj.pe1_order

                case 0
                    obj.trajectory = linearTrajectory(obj.nrKsteps);

                case 2
                    obj.trajectory = zigzagTrajectory(obj.nrKsteps);

                case 3
                    obj.trajectory = arrayTrajectory(obj.NO_VIEWS,obj.gp_var_mul);

                case {4, 5} % reserved for 3D P2ROUD acquisition
                    flist = dir(fullfile(app.mrdImportPath,'*.txt'));
                    if ~isempty(flist)
                        obj.trajectoryFileName = flist(1).name;
                        obj.trajectory = load([flist(1).folder,filesep,flist(1).name]);
                        obj.dataType = '3Dp2roud';
                        app.TextMessage(strcat('P2ROUD trajectory: ',{' '},flist(1).name));
                    else
                        obj.dataType = '3D';
                        obj.trajectory = linearTrajectory(obj.nrKsteps);
                        app.TextMessage('WARNING: P2ROUD trajectory file not found, assuming linear k-space filling ...');
                        app.SetStatus(1);
                    end

                case 8 % Reserved for 2D radial
                    if app.HalfCircleRadialButton.Value == 1              trajType = 1; end
                    if app.FullCircleRadialButton.Value == 1              trajType = 2; end
                    if app.UserDefinedRadialButton.Value == 1             trajType = 3; end
                    obj.trajectory = radialTrajectory(obj.nrKsteps, obj.nr_repetitions, obj.NO_SLICES, app.RadialAngleEditField.Value,trajType);

                case 9 % 3D UTE
                    flist = dir(fullfile(app.mrdImportPath,'lut*.txt'));
                    try
                        obj.trajectory = load([flist(1).folder,filesep,flist(1).name]);
                        obj.gradTrajectory = load('ktraj.txt');
                        obj.dataType = '3Dute';
                        app.TextMessage(strcat('3D UTE trajectory: ',{' '},flist(1).name));
                        obj.trajectory = reshape(obj.trajectory,[3,length(obj.trajectory)/3]);
                        obj.trajectory = obj.trajectory(:,1:obj.nr_repetitions);
                    catch ME
                        app.TextMessage(ME.message);
                        obj.trajectory = linearTrajectory(obj.nrKsteps);
                        app.TextMessage('ERROR: 3D UTE trajectory not found or invalid ...');
                        app.SetStatus(2);
                        app.obj.validDataFlag = false;
                    end

            end


            % ---------------------------------------------------------------------------------
            % Linear trajectory
            % ---------------------------------------------------------------------------------
            function traject = linearTrajectory(nrPE)

                for i = 1:nrPE
                    traj(i) = i;
                end

                traject = traj;

            end % linearTrajectory


            % ---------------------------------------------------------------------------------
            % Zig-zag trajectory
            % ---------------------------------------------------------------------------------
            function traject = zigzagTrajectory(nrPE)

                for i = 1:nrPE/2
                    traj(i) = 2*(i-1)+1;
                    traj(i+nrPE/2) = nrPE-2*i+2;
                end

                traject = traj;

            end % zigzagTrajectory


            % ---------------------------------------------------------------------------------
            % Array trajectory
            % ---------------------------------------------------------------------------------
            function traject = arrayTrajectory(nrPE,peArray)

                traject = round(peArray(1:nrPE) - min(peArray) + 1);

            end % arrayTrajectory


            % ---------------------------------------------------------------------------------
            % Radial trajectory
            % ---------------------------------------------------------------------------------
            function traject = radialTrajectory(nrPE,nrReps,nrSlices,radialAngle,trajType)

                % Radial k-space trajectory
                % The actual trajectory is later incorporated in the trajectory that is used by the Bart toolbox
                % Therefore the trajectory is simply linear here

                switch trajType

                    case 1

                        % Half circle linear
                        fullAngle = 180;

                        cnt = 1;
                        for k = 1:nrSlices
                            for i = 1:nrReps
                                for j = 1:nrPE
                                    traject(cnt) = (j-1)*fullAngle/nrPE;
                                    cnt = cnt + 1;
                                end
                            end
                        end

                    case 2

                        % Full circle linear
                        fullAngle = 360;

                        cnt = 1;
                        for k = 1:nrSlices
                            for i = 1:nrReps
                                for j = 1:nrPE
                                    traject(cnt) = (j-1)*fullAngle/nrPE;
                                    cnt = cnt + 1;
                                end
                            end
                        end

                    case 3

                        % User defined trajectory
                        cnt = 1;
                        for k = 1:nrSlices
                            for i = 1:nrReps
                                for j = 1:nrPE
                                    traject(cnt) = mod(j*radialAngle,360);
                                    if traject(cnt) > 180
                                        traject(cnt) = traject(cnt) - 180;
                                    end
                                    cnt = cnt + 1;
                                end
                            end
                        end

                end

            end % radialTrajectory

        end % getKspaceTrajectory



        % ---------------------------------------------------------------------------------
        % Back to K-space
        % ---------------------------------------------------------------------------------
        function obj = backToKspace(obj)

            obj.kSpaceMrd = [];
            [nf, ~, ~, dimZ, nr] = size(obj.movieExp);

            switch obj.dataType

                case {'2D','2Dms','2Dradial','2Dradialms'}

                    for i = 1:nf
                        for j = 1:nr
                            for k = 1:dimZ
                                obj.kSpaceMrd(i,:,:,k,j) = retro.fft2Dmri(squeeze(obj.movieExp(i,:,:,k,j)));
                            end
                        end
                    end

                    % Samples, views, views2, slices, echoes (frames), experiments
                    obj.kSpaceMrd = permute(obj.kSpaceMrd,[2,3,6,4,1,5]);

                case {'3D','3Dp2roud','3Dute'}

                    for i = 1:nf
                        for j = 1:nr
                            obj.kSpaceMrd(i,:,:,:,j) = retro.fft3Dmri(squeeze(obj.movieExp(i,:,:,:,j)));
                        end
                    end

                    % Samples, views, views2, slices, echoes (frames), experiments
                    obj.kSpaceMrd = permute(obj.kSpaceMrd,[2,3,4,6,1,5]);

            end

        end





        % ---------------------------------------------------------------------------------
        % Read navigator up/down flip switch value
        % ---------------------------------------------------------------------------------
        function obj = readFlipSwitch(obj, app)

            switch app.NavigatorFlipSwitch.Value

                case 'Up'
                    obj.upDown = 1;

                case 'Down'
                    obj.upDown = -1;

            end

        end % readFlipSwitch



        % ---------------------------------------------------------------------------------
        % Auto determine the number of discarded points between navigator and imaging data
        % ---------------------------------------------------------------------------------
        function obj = autoDiscardNav(obj, app)

            switch obj.dataType

                case {'2D','2Dms','3D','3Dp2roud'}

                    % Sum data over all coils
                    dataNav = zeros(size(obj.data{1}));
                    for cnt = 1:obj.nr_coils
                        dataNav = dataNav + abs(obj.data{cnt});
                    end

                    % Sum all the data into 1 echo
                    dataNav = squeeze(sum(dataNav,[1,2,3]));
                    dataNav = dataNav(obj.primaryNavigatorPoint:end);
                    dataNav = smooth(dataNav,5);

                    % Estimate the location of the echo
                    echoLocation = round(length(dataNav)/2);

                    % Find the location of the minimum between navigator and echo
                    [~,idx] = min(dataNav(1:echoLocation));

                    % Return a reasonable value for the number of discarded points
                    if idx<obj.minDiscard || idx>obj.maxDiscard
                        obj.nrNavPointsDiscarded = obj.defaultDiscard; % default value
                    else
                        obj.nrNavPointsDiscarded = idx;
                    end

                    % Set the value
                    app.DiscardEditField.Value = obj.nrNavPointsDiscarded;

            end

        end % autoDiscardNav



        % ---------------------------------------------------------------------------------
        % Extract the navigator signals
        % ---------------------------------------------------------------------------------
        function obj = extractNavigator(obj)

            switch obj.dataType

                case '2D'
                    extractNavigator2D;

                case '2Dms'
                    extractNavigator2Dms;

                case {'3D','3Dp2roud'}
                    extractNavigator3D;

                case '2Dradial'
                    extractNavigator2Dradial;

                case '2Dradialms'
                    extractNavigator2Dradialms;

                case '3Dute'
                    extractNavigator3Dute;

            end


            % ---------------------------------------------------------------------------------
            % ----- 2D single-slice data -------------------------------------------------------
            % ---------------------------------------------------------------------------------
            function amplitude = extractNavigator2D

                % Extracts the navigator data from the raw k-space data
                % Outputs a 1D array of doubles

                obj.navAmplitude = cell(obj.nr_coils);

                for coilnr = 1:obj.nr_coils

                    % Size of the input data, 4th dimension is the readout direction which contains the navigator
                    [nrEperiments,dimZ,dimY,dimX] = size(obj.data{coilnr});

                    % Extract the navigator and put it in a long array
                    navdataAmplitude = reshape(permute(obj.data{coilnr},[3,2,1,4]),nrEperiments*dimY*dimZ,dimX);

                    if obj.nrNavPointsUsed > 1

                        % Take the principal component of the data
                        dataNav = navdataAmplitude(:,obj.primaryNavigatorPoint-obj.nrNavPointsUsed+1:obj.primaryNavigatorPoint);
                        [coeff,~,~] = pca(dataNav);
                        dataPCA = dataNav*coeff;

                        % Take the principal component of the data
                        amplitude = abs(dataPCA(:,1))';

                    else

                        % Single nav point
                        amplitude = abs(navdataAmplitude(:,obj.primaryNavigatorPoint))';

                    end

                    % Detrend
                    amplitude(1,:) = detrend(amplitude(1,:));

                    % Make a guess whether the respiration peaks are positive or negative
                    nrElements = length(amplitude);
                    firstElement = round(0.4*nrElements);
                    lastElement = round(0.8*nrElements);
                    maxAmplitude = abs(max(detrend(amplitude(1,firstElement:lastElement))));
                    minAmplitude = abs(min(detrend(amplitude(1,firstElement:lastElement))));
                    if minAmplitude > maxAmplitude
                        amplitude(1,:) = -amplitude(1,:);
                    end

                    % Multiple with +1 or -1 depending on switch
                    amplitude = amplitude * obj.upDown;

                    % Return the final nav amplitude array
                    obj.navAmplitude{coilnr} = amplitude;

                end

            end % extractNavigator2D




            % ---------------------------------------------------------------------------------
            % ----- 2D multi-slice data -------------------------------------------------------
            % ---------------------------------------------------------------------------------
            function amplitude = extractNavigator2Dms

                % Extracts the navigator data from the raw k-space data of multi-slice 2D data
                % Outputs a 1D array of doubles

                obj.navAmplitude = cell(obj.nr_coils);

                for coilnr = 1:obj.nr_coils

                    % Size of the input data, 4th dimension is the readout direction which contains the navigator
                    [nrExperiments,dimZ,dimY,dimX] = size(obj.data{coilnr});

                    % Extract the navigator and put it in a long array
                    % Y-dimension, repetitions, slice, readout
                    navDataAmplitude = reshape(permute(obj.data{coilnr},[3,1,2,4]),nrExperiments*dimY*dimZ,dimX);

                    if obj.nrNavPointsUsed > 1

                        % Take the principal component of the data
                        dataNav = navDataAmplitude(:,obj.primaryNavigatorPoint-obj.nrNavPointsUsed+1:obj.primaryNavigatorPoint);
                        [coeff,~,~] = pca(dataNav);
                        dataPCA = dataNav*coeff;

                        % Take the principal component of the data
                        amplitude = abs(dataPCA(:,1))';

                    else

                        % Single nav point
                        amplitude = abs(navDataAmplitude(:,obj.primaryNavigatorPoint))';

                    end

                    % Detrend
                    amplitude(1,:) = detrend(amplitude(1,:));

                    % Make a guess whether the respiration peaks are positive or negative in the different slices
                    % This will not be needed with out-of-slice navigator

                    nrElements = length(amplitude);
                    nrElementsPerSlice = round(nrElements/dimZ);

                    for i = 1:dimZ

                        % First and last navigator point for each slice
                        firstElement0 = (i-1)*nrElementsPerSlice + 1;
                        lastElement0 = i*nrElementsPerSlice;

                        % Only look at part of that data away from the start to prevent transient effects
                        firstElement1 = (i-1)*nrElementsPerSlice + 1 + round(0.4*nrElementsPerSlice);
                        lastElement1 = i*nrElementsPerSlice - round(0.1*nrElementsPerSlice);

                        % Min/max of navigator
                        maxAmplitude = abs(max(detrend(amplitude(1,firstElement1:lastElement1))));
                        minAmplitude = abs(min(detrend(amplitude(1,firstElement1:lastElement1))));

                        if minAmplitude > maxAmplitude
                            amplitude(1,firstElement0:lastElement0) = -amplitude(1,firstElement0:lastElement0);
                        end

                        amplitude(1,firstElement0:lastElement0) = detrend(amplitude(1,firstElement0:lastElement0));

                    end

                    % Multiple with +1 or -1 depending on global switch
                    amplitude = amplitude * obj.upDown;

                    % Return the final nav amplitude
                    obj.navAmplitude{coilnr} = amplitude;

                end

            end % extractNavigator2Dms




            % ---------------------------------------------------------------------------------
            % ----- 3D data -------------------------------------------------------------------
            % ---------------------------------------------------------------------------------
            function amplitude = extractNavigator3D

                % Extracts the navigator data from the raw k-space data
                % Outputs a 1D array of doubles

                obj.navAmplitude = cell(obj.nr_coils);

                for coilnr = 1:obj.nr_coils

                    % Size of the input data, 4th dimension is the readout direction which contains the navigator
                    [nrRepetitions,dimZ,dimY,dimX] = size(obj.data{coilnr});

                    % Extract the navigator and put it in a long array
                    navDataAmplitude = reshape(permute(obj.data{coilnr},[2,3,1,4]),nrRepetitions*dimY*dimZ,dimX);

                    if obj.nrNavPointsUsed > 1

                        % Take the principal component of the data
                        dataNav = navDataAmplitude(:,obj.primaryNavigatorPoint-obj.nrNavPointsUsed+1:obj.primaryNavigatorPoint);
                        [coeff,~,~] = pca(dataNav);
                        dataPCA = dataNav*coeff;

                        % Take the principal component of the data
                        amplitude = abs(dataPCA(:,1))';

                    else

                        % Single nav point
                        amplitude = abs(navDataAmplitude(:,obj.primaryNavigatorPoint))';

                    end

                    % Detrend
                    amplitude(1,:) = detrend(amplitude(1,:));

                    % Make a guess whether the respiration peaks are positive or negative
                    nrElements = length(amplitude);
                    firstElement = round(0.4*nrElements);
                    lastElement = round(0.6*nrElements);
                    maxAmplitude = abs(max(amplitude(1,firstElement:lastElement)));
                    minAmplitude = abs(min(amplitude(1,firstElement:lastElement)));
                    if minAmplitude > maxAmplitude
                        amplitude = -amplitude;
                    end

                    % Multiple with +1 or -1 depending on switch
                    amplitude = amplitude * obj.upDown;

                    % Return the final nav amplitude
                    obj.navAmplitude{coilnr} = amplitude;

                end

            end % extractNavigator3D



            % ---------------------------------------------------------------------------------
            % ----- 2D radial data ------------------------------------------------------------
            % ---------------------------------------------------------------------------------
            function extractNavigator2Dradial

                % Extracts the navigator data from the raw k-space data
                % Outputs a 1D array of doubles

                obj.navAmplitude = cell(obj.nr_coils);

                for coilnr = 1:obj.nr_coils

                    % Size of the input data, 4th dimension is the readout direction which contains the navigator
                    [nrRepetitions,dimZ,dimY,dimX] = size(obj.data{coilnr});

                    % Extract the navigator and put it in a long array:      y repetions slices, x
                    navDataAmplitude = reshape(permute(obj.data{coilnr},[3,1,2,4]),nrRepetitions*dimY*dimZ,dimX);

                    if obj.nrNavPointsUsed > 1

                        range = round(obj.nrNavPointsUsed/2);

                        % Take the principal component of the data
                        dataNav = navDataAmplitude(:,obj.primaryNavigatorPoint-range:obj.primaryNavigatorPoint+range);
                        [coeff,~,~] = pca(dataNav);
                        dataPCA = dataNav*coeff;

                        % Take the principal component of the data
                        amplitude = abs(dataPCA(:,1))';
                        phase = angle(dataPCA(:,1))';

                    else

                        % Single nav point
                        amplitude = abs(navDataAmplitude(:,obj.primaryNavigatorPoint))';
                        phase = angle(navDataAmplitude(:,obj.primaryNavigatorPoint))';

                    end

                    % Detrend
                    amplitude(1,:) = detrend(amplitude(1,:));

                    % Make a guess whether the respiration peaks are positive or negative
                    nrElements = length(amplitude);
                    firstElement = round(0.4*nrElements);
                    lastElement = round(0.6*nrElements);
                    maxAmplitude = abs(max(amplitude(1,firstElement:lastElement)));
                    minAmplitude = abs(min(amplitude(1,firstElement:lastElement)));
                    if minAmplitude > maxAmplitude
                        amplitude = -amplitude;
                    end

                    % Multiple with +1 or -1 depending on switch
                    amplitude = amplitude * obj.upDown;

                    % Return the final nav amplitude
                    obj.navAmplitude{coilnr} = amplitude;
                    obj.navPhase{coilnr} = phase;

                end

            end % extractNavigator2Dradial



            % ---------------------------------------------------------------------------------
            % ----- 2D muslti-slice radial data ------------------------------------------------------------
            % ---------------------------------------------------------------------------------
            function extractNavigator2Dradialms

                % Extracts the navigator data from the raw k-space data
                % Outputs a 1D array of doubles

                obj.navAmplitude = cell(obj.nr_coils);

                for coilnr = 1:obj.nr_coils

                    % Size of the input data, 4th dimension is the readout direction which contains the navigator
                    [nrExperiments,dimZ,dimY,dimX] = size(obj.data{coilnr});

                    % Extract the navigator and put it in a long array
                    % Y-dimension, repetitions, slice, readout
                    navDataAmplitude = reshape(permute(obj.data{coilnr},[3,1,2,4]),nrExperiments*dimY*dimZ,dimX);

                    if obj.nrNavPointsUsed > 1

                        % Take the principal component of the data
                        dataNav = navDataAmplitude(:,obj.primaryNavigatorPoint-obj.nrNavPointsUsed+1:obj.primaryNavigatorPoint);
                        [coeff,~,~] = pca(dataNav);
                        dataPCA = dataNav*coeff;

                        % Take the principal component of the data
                        amplitude = abs(dataPCA(:,1))';

                    else

                        % Single nav point
                        amplitude = abs(navDataAmplitude(:,obj.primaryNavigatorPoint))';

                    end

                    % Detrend
                    amplitude(1,:) = detrend(amplitude(1,:));

                    % Make a guess whether the respiration peaks are positive or negative in the different slices
                    % This will not be needed with out-of-slice navigator

                    nrElements = length(amplitude);
                    nrElementsPerSlice = round(nrElements/dimZ);

                    for i = 1:dimZ

                        % First and last navigator point for each slice
                        firstElement0 = (i-1)*nrElementsPerSlice + 1;
                        lastElement0 = i*nrElementsPerSlice;

                        % Only look at part of that data away from the start to prevent transient effects
                        firstElement1 = (i-1)*nrElementsPerSlice + 1 + round(0.4*nrElementsPerSlice);
                        lastElement1 = i*nrElementsPerSlice - round(0.1*nrElementsPerSlice);

                        % Min/max of navigator
                        maxAmplitude = abs(max(detrend(amplitude(1,firstElement1:lastElement1))));
                        minAmplitude = abs(min(detrend(amplitude(1,firstElement1:lastElement1))));

                        if minAmplitude > maxAmplitude
                            amplitude(1,firstElement0:lastElement0) = -amplitude(1,firstElement0:lastElement0);
                        end

                        amplitude(1,firstElement0:lastElement0) = detrend(amplitude(1,firstElement0:lastElement0));

                    end

                    % Multiple with +1 or -1 depending on global switch
                    amplitude = amplitude * obj.upDown;

                    % Return the final nav amplitude
                    obj.navAmplitude{coilnr} = amplitude;

                end

            end % extractNavigator2Dradialms



            % ---------------------------------------------------------------------------------
            % ----- 3D UTE data ---------------------------------------------------------------
            % ---------------------------------------------------------------------------------
            function extractNavigator3Dute

                % Extracts the navigator data from the raw k-space data
                % Outputs a 1D array of doubles

                obj.navAmplitude = cell(obj.nr_coils);

                for coilnr = 1:obj.nr_coils

                    % Size of the input data, 4th dimension is the readout direction which contains the navigator
                    [nrRepetitions,~,~,dimX] = size(obj.data{coilnr});

                    % Extract the navigator and put it in a long array
                    navDataAmplitude = reshape(permute(obj.data{coilnr},[3,2,1,4]),nrRepetitions,dimX);

                    if obj.nrNavPointsUsed > 1

                        range = obj.nrNavPointsUsed;

                        % Take the principal component of the data
                        dataNav = navDataAmplitude(:,obj.primaryNavigatorPoint:obj.primaryNavigatorPoint+range);
                        [coeff,~,~] = pca(dataNav);
                        dataPCA = dataNav*coeff;

                        % Take the principal component of the data
                        amplitude = abs(dataPCA(:,1))';
                        phase = angle(dataPCA(:,1))';

                    else

                        % Single nav point
                        amplitude = abs(navDataAmplitude(:,obj.primaryNavigatorPoint))';
                        phase = angle(navDataAmplitude(:,obj.primaryNavigatorPoint))';

                    end

                    % Detrend
                    amplitude(1,:) = detrend(amplitude(1,:));

                    % Make a guess whether the respiration peaks are positive or negative
                    nrElements = length(amplitude);
                    firstElement = round(0.4*nrElements);
                    lastElement = round(0.6*nrElements);
                    maxAmplitude = abs(max(amplitude(1,firstElement:lastElement)));
                    minAmplitude = abs(min(amplitude(1,firstElement:lastElement)));
                    if minAmplitude > maxAmplitude
                        amplitude = -amplitude;
                    end

                    % Multiple with +1 or -1 depending on switch
                    amplitude = amplitude * obj.upDown;

                    % Return the final nav amplitude
                    obj.navAmplitude{coilnr} = amplitude;
                    obj.navPhase{coilnr} = phase;

                end

            end % extractNavigator3Dute



        end % extractNavigator



        % ---------------------------------------------------------------------------------
        % Determine the power-frequency spectrum of the navigator
        % ---------------------------------------------------------------------------------
        function  obj = determinePowerSpectrumPCA(obj)

            if obj.nr_coils > 1

                dataNav = zeros([length(obj.navAmplitude{1}),obj.nr_coils]);
                for i = 1:obj.nr_coils
                    dataNav(:,i) = obj.navAmplitude{i};
                end

                % Take the principal component of the data
                [coeff,~,~] = pca(dataNav);
                dataPCA = dataNav*coeff;
                amplitude = dataPCA(:,1);

            else

                amplitude = obj.navAmplitude{1}';

            end

            % Include only those navigators when includewindow == 1 and excludewindow == 1
            amplitude = amplitude.*obj.includeWindow.*obj.excludeWindow;

            % Determine the frequency power spectrum
            y = fft(amplitude);
            fs = 1000/obj.TR;                       % Sample frequency in Hz
            n = length(amplitude);                      % Number of samples
            obj.frequency = (0:n-1)*(fs/n)*60;       % Frequency range in bpm
            power = abs(y).^2/n;                        % Power of the DFT

            % Determine frequency and harmonics of k-space trajectory and set those to zero
            if obj.NO_VIEWS > 1

                kfreq = 60/(0.001*obj.NO_VIEWS*obj.TR);
                ifreq = (fs/n)*60;

                for i = 1:10
                    try
                        power(round(i*kfreq/ifreq)+1) = 0;
                        power(round(i*kfreq/ifreq))   = 0;
                        power(round(i*kfreq/ifreq)-1) = 0;
                    catch
                    end
                end

            end

            % Smooth the power spectrum with moving average
            power = medfilt1(power,6);
            obj.powerSpectrum = power;

            % Detect heart rate
            minheartbpm = obj.physioFilterSettings(1);
            maxheartbpm = obj.physioFilterSettings(2);
            minidx = round(minheartbpm*n/(fs*60));
            maxidx = round(maxheartbpm*n/(fs*60));
            [~, idx] = max(power(minidx:maxidx));
            obj.detectedHR = round(idx*fs*60/n + minheartbpm);

            % Detect respiratory rate
            minRRbpm = obj.physioFilterSettings(3);
            maxRRbpm = obj.physioFilterSettings(4);
            minidx = round(minRRbpm*n/(fs*60));
            maxidx = round(maxRRbpm*n/(fs*60));
            [~, idx] = max(power(minidx:maxidx));
            obj.detectedRR = round(idx*fs*60/n + minRRbpm);

        end % determinePowerSpectrumPCA



        % ---------------------------------------------------------------------------------
        % Filter the navigator
        % ---------------------------------------------------------------------------------
        function obj = filterNavPCA(obj)

            % Applies a bandwidth filter on the navigator data

            sf = 1000/obj.TR;                   % Sampling frequency in Hz = 1/TR[ms]
            respHarmonics = 2;                      % Number of higher order harmonics for respiratory frequency, 2 = main + 1 harmonic
            order = obj.physioFilterSettings(5); % Butterworth filter order

            if obj.nr_coils > 1

                dataNav = zeros([length(obj.navAmplitude{1}),obj.nr_coils]);
                for i = 1:obj.nr_coils
                    dataNav(:,i) = obj.navAmplitude{i};
                end

                % Take the principal component of the data
                [coeff,~,~] = pca(dataNav);
                dataPCA = dataNav*coeff;

                % Take the principal component of the data
                amplitude = dataPCA(:,1);

            else

                amplitude = obj.navAmplitude{1}';

            end

            % Filter for heart motion
            hrf = obj.detectedHR/60;       % Expected heartrate in Hz = hr[bpm]/60
            bwh = obj.bandwidthHR/60;      % Bandwidth heartrate in Hz = [bpm]/60
            [b,a] = butter(order,[hrf-0.5*bwh,hrf+0.5*bwh]/(sf/2),'bandpass');     % Butterworth bandpass filter
            heartOutputData = filtfilt(b,a,amplitude);                             % Apply zero-phase shift filtering

            % Detrend
            heartOutputData = detrend(heartOutputData);

            % Normalize envelope
            factor = round(100/hrf);      % Adapt the envelope setting to match the expected heart rate frequency
            [env,~] = envelope(heartOutputData,factor,'peak');
            obj.heartNav = heartOutputData./abs(env);

            % Filter for respiration motion
            while true

                rrf = obj.detectedRR/60;       % Expected resprate in Hz = rr[bpm]/60
                bwr = obj.bandwidthRR/60;      % Bandwidth resprate in Hz = [bpm]/60

                respOutputData = zeros(size(amplitude));

                if obj.detectedRR<45
                    [b, a] = butter(order,(rrf+0.5*bwr)/(sf/2),'low');    % Butterworth lowpass filter for low frequencies
                    respOutputData = filtfilt(b,a,amplitude);
                else
                    for i = 1:respHarmonics
                        [b, a] = butter(order,[i*rrf-0.5*bwr,i*rrf+0.5*bwr]/(sf/2),'bandpass');      % Butterworth bandpass filter
                        respOutputData = respOutputData + (1/i^2.75)*filtfilt(b,a,amplitude);        % Apply zero-phase shift filtering
                    end
                end

                % In some cases the filter produces NaN when filtering with too low
                % frequency. In those cases the respirate and bandwidth will be
                % increased until there are no more NaN

                if sum(isnan(respOutputData))==0
                    break;
                end

                obj.detectedRR = obj.detectedRR+1;
                obj.bandwidthRR = obj.bandwidthRR+1;

            end

            % Detrend and normalize envelope
            respOutputData = detrend(respOutputData);

            % Normalize envelope
            factor = round(150/rrf);  % Adapt the envelope setting to match the expected respiration rate frequency
            [env,~] = envelope(respOutputData,factor,'peak');
            obj.respNav = respOutputData./abs(env);

            % Return the principal component, or in case 1 coil the orignal navigator data, and detrend
            obj.PCANav = detrend(amplitude);

        end % filterNavPCA



        % ---------------------------------------------------------------------------------
        % Heart filter shape for plotting
        % ---------------------------------------------------------------------------------
        function output = filterShape(obj, app, type)

            % Determines the filter shape for display purposes
            % Order = order of the filter
            % Type = 'heart' or 'resp'

            order = app.FilterOrderEditField.Value;
            sf = 1000/obj.TR;           % Sampling frequency in Hz = 1/TR[ms]

            if strcmp(type,'resp')
                hr = obj.detectedRR;
                bw = obj.bandwidthRR;
            else
                hr = obj.detectedHR;
                bw = obj.bandwidthHR;
            end

            cf = hr/60;       % Expected heartrate in Hz = hr[bpm]/60
            bw1 = bw/60;      % Bandwidth in Hz = [bpm]/60

            % Check that bw1 is smaller than center-frequency
            if cf <= bw1
                bw1 = 0.98*cf;
            end

            if hr < 45
                % Butterworth lowpass filter for low frequencies
                [b, a] = butter(order,(cf+0.5*bw1)/(sf/2),'low');
            else
                % Butterworth bandpass filter
                [b,a] = butter(order,[cf-0.5*bw1,cf+0.5*bw1]/(sf/2),'bandpass');
            end

            [h,w] = freqz(b,a,2000,sf);
            output = [w*60, abs(h)];

        end



        % ---------------------------------------------------------------------------------
        % Determine the cardiac trigger points
        % ---------------------------------------------------------------------------------
        function obj = trigPointsHeart(obj)

            % Extracts the ECG trigger points from the navigators

            % Find the peaks and locations = fast
            try
                obj.heartTrigPoints = retro.peakFinder(obj.heartNav',[],[],[],false,true);
            catch
                obj.heartTrigPoints = []; % Failure
            end

            % Trigger points are in units of samples (actual time is thus heartTrigPoints*TR)

            % Backup plan in case peakFinder fails = slower
            if length(obj.heartTrigPoints)<20

                % Minimal distance 50% of expected heart rate [in points]
                if(isnan(obj.meanHeartRate))
                    obj.meanHeartRate = 500;
                end
                dist = 0.50*(60/obj.meanHeartRate)/(obj.TR/1000);
                interPolationFactor = obj.splineFactor;

                % Cubic spline interpolation of the data
                nrl = length(obj.heartNav);
                navi = interp1(1:nrl,obj.heartNav(1:nrl),1:1/interPolationFactor:nrl,'spline');

                % Find the peaks and locations
                [~,locs]=findpeaks(navi,'MinPeakDistance',dist*interPolationFactor);
                locs = locs + interPolationFactor/2;

                % Recalculate orginal fractional time point peak positions
                obj.heartTrigPoints = locs/interPolationFactor;

            end

            obj.heartTrigTime = obj.heartTrigPoints*obj.TR - obj.TR; % Time starts at zero

        end % trigPointsHeart



        % ---------------------------------------------------------------------------------
        % Determine the respiration trigger points
        % ---------------------------------------------------------------------------------
        function obj = trigPointsResp(obj)

            % Extracts the ECG trigger points from the navigators

            % Find the peaks and locations
            try
                obj.respTrigPoints = retro.peakFinder(obj.respNav',[],[],[],false,true);
            catch
                obj.respTrigPoints = []; % Failure
            end

            % Backup plan in case peakFinder fails
            if length(obj.respTrigPoints)<10

                interPolationFactor = obj.splineFactor;
                nrl = size(obj.respNav,1);
                navi = interp1(1:nrl,obj.respNav,1:1/interPolationFactor:nrl,'spline');

                % Find the peaks and locations
                [~,locs] = findpeaks(navi,'MinPeakProminence',0.1);
                locs = locs + interPolationFactor/2;

                % Recalculate orginal fractional time point peak positions
                obj.respTrigPoints = locs/interPolationFactor;

            end

            % Trigger points are in units of samples (actual time is thus heartTrigPoints*TR - TR)
            obj.respTrigTime = obj.respTrigPoints*obj.TR  - obj.TR; % Time starts at zero

        end % trigPointsResp



        % ---------------------------------------------------------------------------------
        % Determine average heart rate
        % ---------------------------------------------------------------------------------
        function obj = calcHeartRate(obj)
            % This function calculates heart rate as a function of time in bpm

            nrLocs = size(obj.heartTrigPoints,2); % Determine the number of trigger points

            rate = zeros(nrLocs-1,1); % Initialize the heart rate array
            for i=1:nrLocs-1
                % Calculate heart rate using adjacent trigger points and TR (repetition time)
                rate(i) = 60/((obj.heartTrigPoints(i+1)-obj.heartTrigPoints(i))*(obj.TR/1000));
            end

            hr = rate'; % Transpose the heart rate array
            hrf = movmedian(hr,32); % Smooth the trendline using a moving median filter

            % Determine the mean cardiac rate
            includeData = obj.includeWindow.*obj.excludeWindow; % Determine the included data window
            try
                includeData = round(resample(includeData,2*length(hr),length(includeData)));
            catch
                includeData = round(resample(includeData,3*length(hr),length(includeData))); % Resample the data-window to match the number of samples in the heart rate
            end
            includeData = round(resample(includeData,length(hr),length(includeData))); % Resample again to prevent overflow error

            obj.meanHeartRate = round(median(nonzeros(includeData'.*hr))); % Calculate the median heart rate using the included data
            obj.heartRateTime = hr; % Store the heart rate array
            obj.heartRateTimeFiltered = hrf; % Store the filtered heart rate array

        end % calcHeartRate



        % ---------------------------------------------------------------------------------
        % Determine average respiration rate
        % ---------------------------------------------------------------------------------
        function obj = calcRespRate(obj)

            % Determine respiration rate as function of time in bpm

            nrLocs = size(obj.respTrigPoints,2);

            rate = zeros(nrLocs-1,1);
            for i=1:nrLocs-1
                rate(i) = 60/((obj.respTrigPoints(i+1)-obj.respTrigPoints(i))*(obj.TR/1000));
            end

            resp = rate';
            respf = movmedian(resp,32);     % Smooth trendline

            % Determine the mean respiration rate
            includeData = obj.includeWindow.*obj.excludeWindow;                         % Data-window which is included
            try
                includeData = round(resample(includeData,2*length(resp),length(includeData)));  % Resample the data-window to nr samples respiration rate
            catch
                includeData = round(resample(includeData,3*length(resp),length(includeData)));  % Resample the data-window to nr samples respiration rate
            end
            includeData = round(resample(includeData,length(resp),length(includeData)));        % In 2 steps, to prevent an overflow error

            obj.meanRespRate = round(median(nonzeros(includeData'.*resp)));                  % Take the median of the respirationrate
            obj.respRateTime = resp;
            obj.respRateTimeFiltered = respf;

        end % calcRespRate



        % ---------------------------------------------------------------------------------
        % Determine respiration windows
        % ---------------------------------------------------------------------------------
        function obj = makeRespWindow(obj, app)

            rLocs = obj.respTrigPoints;
            meanResp = mean(obj.respRateTimeFiltered);
            respPercent = obj.respPercentage;
            nrOfKlines = obj.nrKlines;
            obj.respWindowTime = [];

            % This function creates an array (time line) of rectangular boxes of 0's and 1's around the detected respiratory signals
            % Later on 1 means that there is a respiration, for which the k-lines will be discarded

            respWin = 0.5*(respPercent/100)*(60/meanResp)*1000/obj.TR;   % 1/2 window width around respiratory peak locations

            window = zeros(nrOfKlines,1);   % Fill array with zeros

            if ~app.IgnoreRespirationToggleCheckBox.Value

                for i = 1 : size(rLocs,2)

                    % Center of the respiration window
                    center = rLocs(1,i);

                    % Beginning and end of the respiration window in ms, time starts at zero
                    obj.respWindowTime(i,1) = (center - respWin)*obj.TR  - obj.TR;
                    obj.respWindowTime(i,2) = (center + respWin)*obj.TR  - obj.TR;

                    % Set respiration window mask to 1 during respiration
                    for j = round(center - respWin) : round(center + respWin)
                        if (j>0) && (j<=nrOfKlines)
                            window(j) = 1;
                        end
                    end

                end

                % Set window to expiration (0) or inspiration (1)
                if app.RespirationToggleCheckBox.Value
                    window = 1 - window;
                end

            end

            obj.respWindow = window;

        end % makeRespWindow





        % ---------------------------------------------------------------------------------
        % 2D reconstruction
        % ---------------------------------------------------------------------------------
        function obj = reco2D(obj, app)

            if app.bartDetectedFlag && ~app.MatlabRecoCheckBox.Value

                % Perform CS reconstruction with the Bart toolbox, preferred option
                app.TextMessage('Reconstructing the data with the BART toolbox ...');
                app.ProgressGauge.Value = 0;
                drawnow;
                csReco2Dmc;
                app.ProgressGauge.Value = 100;
                app.validRecoFlag = true;
                app.validSenseMapFlag = true;
                drawnow;

            else

                % Perform CS reconstruction in matlab, slower, but works
                app.TextMessage('WARNING: Reconstructing the data with MATLAB ...');
                app.TextMessage('Slow reconstruction, please be patient ...');
                app.ProgressGauge.Value = 0;
                app.SetStatus(1);
                drawnow;
                csReco2Dmatmc;
                app.ProgressGauge.Value = 100;
                app.validRecoFlag = true;
                app.validSenseMapFlag = true;
                drawnow;

            end


            % ---------------------------------------------------------------------------------
            % 2D compressed sensing reconstruction with Matlab
            % ---------------------------------------------------------------------------------
            function csReco2Dmatmc

                kSpaceIn = obj.kSpace;
                averages = obj.kSpaceAvg;
                lambdaTV = app.TVcineEditField.Value;

                % K-space = {coil}[frames, x, y, slice, dynamic]
                %                    1     2  3    4       5
                [dimF,dimX,dimY,dimZ,dimD] = size(kSpaceIn{1});
                nrCoils = obj.nr_coils;

                % In case of 1 frame, duplicate that frame to facilitate reconstruction
                if dimF == 1
                    for i = 1:nrCoils
                        kSpaceIn{i}(2,:,:,:,:) = kSpaceIn{i}(1,:,:,:,:);
                    end
                end

                % K-space data: x,y,frames,slices,dynamics
                for i = 1:nrCoils
                    kSpaceIn{i} = permute(kSpaceIn{i},[2,3,1,4,5]);
                end

                % K-space data: x,y,frames,slices,dynamics,coils
                kSpaceData = zeros(dimX,dimY,dimF,dimZ,dimD,nrCoils);
                for i = 1:nrCoils
                    kSpaceData(:,:,:,:,:,i) = kSpaceIn{i}(:,:,:,:,:);
                end

                % Averages data: x,y,frames,slices,dynamics
                averages = permute(averages,[2,3,1,4,5]);

                % Reset progress counter
                param.iteration = 0;

                % For convenience make rectangular matrix size of power of 2
                mdimxy = max([dimX,dimY]);
                dimX = 2^nextpow2(mdimxy);
                dimY = dimX;
                dimF = size(kSpaceData,3);
                app.TextMessage(strcat('Reconstruction matrix =',{' '},num2str(dimX),'x',num2str(dimY)));

                % Pre-allocate memory for image_out
                imageOut = zeros(dimF,dimX,dimY,dimZ,dimD);

                % Slice and dynamic loop
                for slice = 1:dimZ

                    for dynamic = 1:dimD

                        % K-space of slice and dynamic
                        kData = squeeze(kSpaceData(:,:,:,slice,dynamic,:));

                        % Zero padding
                        padSizex = round((dimX - size(kData,1))/2);
                        padSizey = round((dimY - size(kData,2))/2);
                        kData = padarray(kData,[padSizex,padSizey,0],'both');
                        kData = kData(1:dimX,1:dimY,:,:);

                        % Size of the data
                        [nx,ny,~,nrCoils] = size(kData);

                        % Normalize the data in the range of approx 0 - 1 for better numerical stability
                        kData = kData/max(abs(kData(:)));

                        % K-space mask: 0 = no data, 1 = data, zero-pad to same size as k-space
                        mask = squeeze(averages(:,:,:,slice,dynamic));
                        mask = padarray(mask,[padSizex,padSizey,0],'both');
                        mask = mask(1:dimX,1:dimY,:);
                        mask = mask./mask;
                        mask(isnan(mask)) = 1;
                        mask = logical(mask);

                        % Coil sensitivity map
                        b1 = ones(nx,ny,nrCoils);

                        % Data
                        param.y = kData;

                        % Reconstruction design matrix
                        param.E = Emat_yxt(mask,b1);

                        % Total variation (TV) constraint in the temporal domain
                        % for 'consistency' with Bart reconstruction, TV seems to be scale empirically by a factor of 8
                        % TV only in the time domain
                        param.TV = TVOP;
                        param.TVWeight = lambdaTV/8;

                        % Number of iterations, 2 x 10 iterations
                        param.nite = 10;
                        param.nouter = 2;
                        param.totaliterations = dimZ * dimD * param.nouter * param.nite;

                        % Add a little bit of randomness, such that linear reco is not exactly right to start with
                        kData1 = randn(size(kData))/2000 + kData;

                        % Linear reconstruction
                        reconDft = param.E'*kData1;

                        % Iterative reconstruction
                        reconCS = reconDft;
                        for n = 1:param.nouter
                            [reconCS,param.iteration] = CSL1NlCg(app,reconCS,param);
                        end

                        % Rearrange to correct orientation: frames, x, y
                        imageTmp = flip(permute(squeeze(abs(reconCS)),[3, 1, 2]),2);

                        % Output reconstructed image
                        imageOut(:,:,:,slice,dynamic) = imageTmp;

                    end

                end

                % Correct back to 1 frame reconstruction
                if dimF == 1
                    imageOut = imageOut(1,:,:,:,:);
                end

                % Shift image in phase-encoding direction if needed
                obj.movieExp = circshift(imageOut,-obj.pixelshift1,3);
                obj.senseMap = ones(size(obj.movieExp));

            end % csReco2Dmatmc


            % ---------------------------------------------------------------------------------
            % 2D compressed sensing reconstruction with Bart toolbox
            % ---------------------------------------------------------------------------------
            function csReco2Dmc

                % app = matlab app
                % kSpaceIn = sorted k-space
                % nrCoils = number of RF receiver coils
                % Wavelet = wavelet L1-norm regularization factor
                % LR = low rank regularization
                % TVt = total variation in CINE dimension
                % TVxy = total variation in xy-dimension regularization
                % TVd = total variation in dynamic dimension
                % ESPIRiT = reconstruction of multi-coil data with ESPIRiT true or false

                kSpaceIn = obj.kSpace;
                Wavelet = app.WVxyzEditField.Value;
                TVxy = app.TVxyzEditField.Value;
                LR = app.LLRxyzEditField.Value;
                TVt = app.TVcineEditField.Value;
                TVd = app.TVdynEditField.Value;
                ESPIRiT = app.ESPIRiTCheckBox.Value;

                % kSpaceIn = {coil}[CINE, x, y, slice, dynamic]
                %                    1    2  3    4       5
                [dimF,dimX,dimY,dimZ,dimD] = size(kSpaceIn{1});
                nrCoils = obj.nr_coils;

                % For convenience make rectangular matrix size of power of 2
                mdimxy = max([dimX,dimY]);
                dimX = 2^nextpow2(mdimxy);
                dimY = dimX;
                app.TextMessage(strcat('Reconstruction matrix =',{' '},num2str(dimX),'x',num2str(dimY)));

                % Resize k-space to next power of 2
                for i = 1:nrCoils
                    kSpaceIn{i} = bart(app,['resize -c 1 ',num2str(dimX),' 2 ',num2str(dimY)],kSpaceIn{i});
                end

                % Put all data in a normal matrix
                kSpaceData = zeros(dimF,dimX,dimY,dimZ,dimD);
                for i = 1:nrCoils
                    kSpaceData(:,:,:,:,:,i) = kSpaceIn{i}(:,:,:,:,:);
                end

                % Bart dimensions  Bart   Matlab
                % 	READ_DIM,       0       1   z
                % 	PHS1_DIM,       1       2   y
                % 	PHS2_DIM,       2       3   x
                % 	COIL_DIM,       3       4   coils
                % 	MAPS_DIM,       4       5   sense maps
                % 	TE_DIM,         5       6
                % 	COEFF_DIM,      6       7
                % 	COEFF2_DIM,     7       8
                % 	ITER_DIM,       8       9
                % 	CSHIFT_DIM,     9       10
                % 	TIME_DIM,       10      11  cardiac / respiratory CINE frames
                % 	TIME2_DIM,      11      12  dynamics
                % 	LEVEL_DIM,      12      13
                % 	SLICE_DIM,      13      14  slices
                % 	AVG_DIM,        14      15

                kSpacePics = permute(kSpaceData,[7,3,2,6,8,9,10,11,12,13,1,5,14,4]);

                % Wavelet in spatial dimensions 2^1+2^2=6
                % Total variation in spatial dimensions 2^1+2^2=6
                % Total variation in cine dimension 2^10 = 1024
                % Total variation in dynamic dimension 2^11 = 2048

                if ESPIRiT && nrCoils>1

                    % ESPIRiT reconstruction
                    app.TextMessage('ESPIRiT reconstruction ...');

                    % Calculate coil sensitivity maps with ecalib bart function
                    kSpacePicsSum = sum(kSpacePics,[11,12]);
                    sensitivities = bart(app,'ecalib -S -I -a', kSpacePicsSum);      % Ecalib with softsense

                else

                    % Reconstruction without sensitivity correction
                    sensitivities = ones(1,dimY,dimX,nrCoils,1,1,1,1,1,1,1,1,1,dimZ);

                end

                % PICS command
                picsCommand = 'pics ';
                if Wavelet>0
                    picsCommand = [picsCommand,' -RW:6:0:',num2str(Wavelet)];
                end
                if TVxy>0
                    picsCommand = [picsCommand,' -R',obj.totalVariation,':6:0:',num2str(TVxy)];
                end
                if LR>0
                    % Locally low-rank in the spatial domain
                    blocksize = round(dimX/16);  % Block size
                    blocksize(blocksize < 8) = 8;
                    picsCommand = [picsCommand,' -RL:6:6:',num2str(LR),' -b',num2str(blocksize)];
                end
                if TVt>0
                    picsCommand = [picsCommand,' -R',obj.totalVariation,':1024:0:',num2str(TVt)];
                end
                if TVd>0
                    picsCommand = [picsCommand,' -R',obj.totalVariation,':2048:0:',num2str(TVd)];
                end

                % BART reconstruction
                imageReco = bart(app, picsCommand, kSpacePics, sensitivities);

                % Sum of squares in the coil dimension
                imageReco = bart(app,'rss 16', imageReco);
                imageReco = abs(imageReco);

                % Reshape the image :           y    x  frames dynamic slices
                imageReco = reshape(imageReco,[dimY,dimX,dimF,dimD,dimZ]);

                % Rearrange to correct orientation: frames, x, y, slices, dynamics
                imageOut = flip(permute(imageReco,[3,2,1,5,4]),3);

                % Sense map orientations: x, y, slices, map1, map2
                senseMap1 = flip(permute(abs(sensitivities),[3,2,14,4,5,1,6,7,8,9,10,11,12,13]),2);

                % Normalize sense map to reasonable value range
                senseMap1 = senseMap1*4095/max(senseMap1(:));

                % Shift image in phase-encoding direction if needed
                obj.movieExp = circshift(imageOut,-obj.pixelshift1,3);
                obj.senseMap = circshift(senseMap1,-obj.pixelshift1,3);

            end % csReco2Dmc

        end % reco2D



        % ---------------------------------------------------------------------------------
        % 3D reconstruction
        % ---------------------------------------------------------------------------------
        function obj = reco3D(obj, app)

            if app.bartDetectedFlag && ~app.MatlabRecoCheckBox.Value

                % Perform CS reconstruction with the Bart toolbox, preferred option
                app.TextMessage('Reconstructing the 3D data with the BART toolbox ...');
                app.ProgressGauge.Value = 0;
                drawnow;
                csReco3Dmc;
                app.ProgressGauge.Value = 100;
                app.validRecoFlag = true;
                app.validSenseMapFlag = true;
                drawnow;

            else

                % Perform cs reconstruction in matlab, slower, but works
                app.TextMessage('WARNING: Reconstructing the data with MATLAB ...');
                app.TextMessage('Slow reconstruction, please be patient ...');
                app.ProgressGauge.Value = 0;
                app.SetStatus(1);
                csReco3Dmatmc;
                app.ProgressGauge.Value = 100;
                app.validRecoFlag = true;
                app.validSenseMapFlag = true;
                drawnow;

            end


            % ---------------------------------------------------------------------------------
            % 3D compressed sensing reconstruction with Matlab
            % ---------------------------------------------------------------------------------
            function csReco3Dmatmc

                kSpaceIn = obj.kSpace;
                averages = obj.kSpaceAvg;
                lambdaTV = app.TVcineEditField.Value;

                % kSpaceIn = {coil}[frames, x, y, z, dynamics]
                %                      1    2  3  4     5
                dimF = size(kSpaceIn{1},1);
                dimD = size(kSpaceIn{1},5);
                nrCoils = obj.nr_coils;

                % In case of 1 frame, duplicate that frame to facilitate reconstruction
                if dimF == 1
                    for i = 1:nrCoils
                        kSpaceIn{i}(2,:,:,:,:) = kSpaceIn{i}(1,:,:,:,:);
                    end
                end

                % K-space data: x,y,z,frames,dynamics
                for i = 1:nrCoils
                    kSpaceIn{i} = permute(kSpaceIn{i},[2,3,4,1,5]);
                end

                % K-space data: x,y,z,frames,dynamics,coils
                kSpaceData = zeros([size(kSpaceIn{1}),nrCoils]);
                for i = 1:nrCoils
                    kSpaceData(:,:,:,:,:,i) = kSpaceIn{i};
                end

                % Averages data: x,y,z,frames,dynamics
                averages = permute(averages,[2,3,4,1,5]);

                % Reset progress counter
                param.iteration = 0;

                % Pad to next power of 2
                dimX = 2^nextpow2(size(kSpaceData,1));
                dimY = 2^nextpow2(size(kSpaceData,2));
                dimZ = size(kSpaceData,3);

                % For convenience make rectangular matrix size
                mdimxy = max([dimX,dimY]);
                dimX = mdimxy;
                dimY = mdimxy;
                dimF = size(kSpaceData,4);

                % Pre-allocate memory for imageOut
                imageOut = zeros(dimF,dimX,dimY,dimZ,dimD);

                for dynamic = 1:dimD

                    % K-space of slice and dynamic
                    kData = squeeze(kSpaceData(:,:,:,:,dynamic,:));

                    % Zero padding
                    padSizex = round((dimX - size(kData,1))/2);
                    padSizey = round((dimY - size(kData,2))/2);
                    padSizez = round((dimZ - size(kData,3))/2);
                    kData = padarray(kData,[padSizex,padSizey,padSizez,0],'both');
                    kData = kData(1:dimX,1:dimY,1:dimZ,:,:);

                    % Size of the data: z,y,x,frames,coils
                    [nx,ny,nz,~,nrCoils] = size(kData);

                    % Normalize the data in the range of approx 0 - 1 for better numerical stability
                    kData = kData/max(abs(kData(:)));

                    % K-space mask: 0 = nodata, 1 = data, zero-pad to same size as k-space
                    mask = squeeze(averages(:,:,:,:,dynamic));
                    mask = padarray(mask,[padSizex,padSizey,padSizez,0],'both');
                    mask = mask(1:dimX,1:dimY,1:dimZ,:);
                    mask = mask./mask;
                    mask(isnan(mask)) = 1;
                    mask = logical(mask);

                    % Coil sensitivity map
                    b1 = ones(nx,ny,nz,nrCoils);

                    % Data
                    param.y = kData;

                    % Reconstruction design matrix
                    param.E = Emat_zyxt(mask,b1);

                    % Total variation (TV) constraint in the temporal domain
                    % for 'consistency' with Bart reconstruction, TV seems to be scale empirically by a factor of 8
                    % TV only in the time domain
                    param.TV = TVOP3D;
                    param.TVWeight = lambdaTV/8;

                    % Number of iterations, 2 x 10 iterations
                    param.nite = 10;
                    param.nouter = 2;
                    param.totaliterations = dimD * param.nouter * param.nite;

                    % Linear reconstruction
                    kData1 = randn(size(kData))/2000 + kData;  % add a little bit of randomness, such that linear reco is not exactly right
                    reconDft = param.E'*kData1;

                    % Iterative reconstruction
                    reconCs=reconDft;
                    for n=1:param.nouter
                        [reconCs,param.iteration] = CSL1NlCg(app,reconCs,param);
                    end

                    % Rearrange to correct orientation: frames, x, y, z
                    imageTmp = flip(permute(abs(reconCs),[4, 1, 2, 3]),2);

                    % Output reconstructed image
                    imageOut(:,:,:,:,dynamic) = imageTmp;

                end

                % Correct back to 1 frame reconstruction
                if dimF == 1
                    imageOut = imageOut(1,:,:,:,:);
                end

                % Shift image in phase-encoding directions if needed
                obj.movieExp = circshift(imageOut,-2-obj.pixelshift1,3);
                obj.movieExp = circshift(imageOut,1-obj.pixelshift2,4);
                obj.senseMap = ones(size(obj.movieExp));

            end % csReco3Dmatmc




            % ---------------------------------------------------------------------------------
            % 3D compressed sensing reconstruction with Bart toolbox
            % ---------------------------------------------------------------------------------
            function csReco3Dmc

                % app = matlab app
                % kSpaceIn = sorted k-space
                % nrCoils = number of RF receiver coils
                % averages_in = number of averages per k-line
                % Wavelet = wavelet L1-norm regularization factor
                % LR = low rank regularization
                % TVt = total variation in time regularization
                % TVxyz = total variation in xyz-dimension regularization
                % ESPIRiT = reconstruction of multi-coil data with ESPIRiT true or false

                kSpaceIn = obj.kSpace;
                Wavelet = app.WVxyzEditField.Value;
                TVxyz = app.TVxyzEditField.Value;
                LR = app.LLRxyzEditField.Value;
                TVt = app.TVcineEditField.Value;
                TVd = app.TVdynEditField.Value;
                ESPIRiT = app.ESPIRiTCheckBox.Value;
                obj.totalVariation = 'T';

                % kSpaceIn = {coil}[frames, x, y, z, dynamics]
                %                    1      2  3  4     5
                dimF = size(kSpaceIn{1},1);
                dimX = 2^nextpow2(size(kSpaceIn{1},2));
                dimY = 2^nextpow2(size(kSpaceIn{1},3));
                dimZ = size(kSpaceIn{1},4);
                dimD = size(kSpaceIn{1},5);
                nrCoils = obj.nr_coils;

                % For convenience make rectangular matrix size
                mdimxy = max([dimX,dimY]);
                dimX = mdimxy;
                dimY = mdimxy;

                % Resize to next power of 2
                for i = 1:nrCoils
                    % kSpaceIn{i} = bart(app,['resize -c 1 ',num2str(dimx),' 2 ',num2str(dimy),' 3 ',num2str(dimz)],kSpaceIn{i});
                    kSpaceIn{i} = bart(app,['resize -c 1 ',num2str(dimX),' 2 ',num2str(dimY)],kSpaceIn{i});
                end

                % K-space suitable for bart
                kSpaceData = zeros(dimF,dimX,dimY,dimZ,dimD,nrCoils);
                for i = 1:nrCoils
                    kSpaceData(:,:,:,:,:,i) = kSpaceIn{i}(:,:,:,:,:);
                end

                % Bart dimensions
                % 	READ_DIM,       1   z
                % 	PHS1_DIM,       2   y
                % 	PHS2_DIM,       3   x
                % 	COIL_DIM,       4   coils
                % 	MAPS_DIM,       5   sense maps
                % 	TE_DIM,         6
                % 	COEFF_DIM,      7
                % 	COEFF2_DIM,     8
                % 	ITER_DIM,       9
                % 	CSHIFT_DIM,     10
                % 	TIME_DIM,       11  cardiac / respiratory CINE
                % 	TIME2_DIM,      12  dynamics
                % 	LEVEL_DIM,      13
                % 	SLICE_DIM,      14  slices
                % 	AVG_DIM,        15

                kSpacePics = permute(kSpaceData,[4,3,2,6,7,8,9,10,11,12,1,5,13,14]);

                % wavelet in spatial dimensions 2^0+2^1+2^2=7
                % total variation in spatial dimensions 2^0+2^1+2^2=7
                % total variation in time 2^10 = 1024
                % total variation in dynamic dimension 2^11 = 2048

                if ESPIRiT && nrCoils>1

                    app.TextMessage('ESPIRiT reconstruction ...');
                    kspace_pics_sum = sum(kSpacePics,[11,12]);
                    sensitivities = bart(app,'ecalib -I -S -a', kspace_pics_sum);

                    app.ProgressGauge.Value = 25;
                    drawnow;

                else

                    % Reconstruction without sensitivity correction
                    sensitivities = ones(dimZ,dimY,dimX,nrCoils,1,1,1,1,1,1,1,1,1,1);

                end

                picsCommand = 'pics -S ';
                if Wavelet>0
                    picsCommand = [picsCommand,' -RW:7:0:',num2str(Wavelet)];
                end
                if TVxyz>0
                    picsCommand = [picsCommand,' -R',obj.totalVariation,':7:0:',num2str(TVxyz)];
                end
                if LR>0
                    % Locally low-rank in the spatial domain
                    blocksize = round(dimX/16);  % Block size
                    blocksize(blocksize < 8) = 8;
                    picsCommand = [picsCommand,' -RL:7:7:',num2str(LR),' -b',num2str(blocksize)];
                end
                if TVt>0
                    picsCommand = [picsCommand,' -R',obj.totalVariation,':1024:0:',num2str(TVt)];
                end
                if TVd>0
                    picsCommand = [picsCommand,' -R',obj.totalVariation,':2048:0:',num2str(TVd)];
                end
                imageReg = bart(app,picsCommand,kSpacePics,sensitivities);

                app.ProgressGauge.Value = 95;
                drawnow;

                % Sum of squares over the coil dimension
                imageReg = bart(app,'rss 16', imageReg);

                % Absolute value
                imageReg = abs(imageReg);

                % Reshape to proper dimensions
                imageReg = reshape(imageReg,[dimZ,dimY,dimX,dimF,dimD]);

                % Rearrange to correct orientation: frames, x, y, z, dynamics
                imageOut = flip(flip(permute(imageReg,[4,3,2,1,5]),3),4);
                imageOut = circshift(imageOut,1,4);

                % Sense map orientations: x, y, z, map1, map2
                senseMap1 = flip(flip(permute(abs(sensitivities),[3,2,1,4,5,6,7,8,9,10,11,12,13,14]),2),3);
                senseMap1 = circshift(senseMap1,1,4);

                % Normalize sense map to reasonable value range
                senseMap1 = senseMap1*4095/max(senseMap1(:));

                % Shift image in phase-encoding directions with integer if needed
                obj.movieExp = circshift(imageOut,-obj.pixelshift1,3);
                obj.movieExp = circshift(imageOut,-obj.pixelshift2,4);
                obj.senseMap = circshift(senseMap1,-obj.pixelshift1,3);
                obj.senseMap = circshift(senseMap1,-obj.pixelshift2,4);

                app.ProgressGauge.Value = 100;
                drawnow;

            end % csReco3Dmc

        end % reco3D



        % ---------------------------------------------------------------------------------
        % 2D radial reconstruction with the Bart toolbox or Matlab
        % ---------------------------------------------------------------------------------
        function obj = reco2Dradial(obj, app)

            if app.bartDetectedFlag && ~app.MatlabRecoCheckBox.Value

                app.TextMessage('Reconstructing 2D radial data with the BART toolbox ...');
                app.ProgressGauge.Value = 0;
                drawnow;

                % Wavelet = wavelet L1-norm regularization factor
                % TVxyz = total variation in xyz-dimension regularization
                % LR = low rank regularization
                % TVt = total variation in time regularization
                % TVd = total variation in dynamic dimension regularization
                % nrCoils = number of RF receiver coils

                Wavelet = app.WVxyzEditField.Value;
                TVxyz = app.TVxyzEditField.Value;
                LR = app.LLRxyzEditField.Value;
                TVt = app.TVcineEditField.Value;
                TVd = app.TVdynEditField.Value;
                nrCoils = obj.nr_coils;
                ESPIRiT = app.ESPIRiTCheckBox.Value;

                % Get the data dimensions from the kSpace object
                dimF = size(obj.kSpace{1},1);
                dimX = size(obj.kSpace{1},2);
                dimY = dimX;
                dimZ = size(obj.kSpace{1},4);
                dimD = size(obj.kSpace{1},5);

                % Retrieve the k-Space, trajectory, and averages from the object
                kSpaceData = zeros(size(obj.kSpace{1}));
                for i = 1:nrCoils
                    kSpaceData(:,:,:,:,:,i) = obj.kSpace{i};      % K-space
                end
                traj = obj.kSpaceTraj;                        % Trajectory
                averages = obj.kSpaceAvg;                     % Averages

                % Bart dimensions  Bart   Matlab
                % 	READ_DIM,       0       1   z
                % 	PHS1_DIM,       1       2   y
                % 	PHS2_DIM,       2       3   x
                % 	COIL_DIM,       3       4   coils
                % 	MAPS_DIM,       4       5   sense maps
                % 	TE_DIM,         5       6
                % 	COEFF_DIM,      6       7
                % 	COEFF2_DIM,     7       8
                % 	ITER_DIM,       8       9
                % 	CSHIFT_DIM,     9       10
                % 	TIME_DIM,       10      11  cardiac / respiratory CINE frames
                % 	TIME2_DIM,      11      12  dynamics
                % 	LEVEL_DIM,      12      13
                % 	SLICE_DIM,      13      14  slices
                % 	AVG_DIM,        14      15

                % Rearrange for BART         1  2  3  4  5  6  7  8  9 10 11 12 13 14
                kSpacePics = permute(kSpaceData,[7, 2, 3, 6,14, 8, 9,10,11,12, 1, 5,13, 4]);

                % Rearrange for BART        1  2  3  4  5  6  7  8  9 10 11 12 13 14
                avgPics = permute(averages,[7, 2, 3, 6,14, 8, 9,10,11,12, 1, 5,13, 4]);

                % Rearrange for BART     1  2  3  4  5  6  7  8  9 10 11 12 13 14
                trajPics = permute(traj,[6, 2, 3,14, 7, 8, 9,10,11,12, 1, 5,13, 4]);

                % Calibration and density correction size
                kdim = round(dimX/2);
                if mod(kdim,2) == 1
                    kdim = kdim + 1;
                end
                kdim(kdim < 32) = 32;
                kdim(kdim > dimX) = dimX;
                calibSize = [kdim, kdim, 1];
                cSize = ['-d',num2str(calibSize(1)),':',num2str(calibSize(2)),':1'];
                app.TextMessage(strcat('Calibration size = ',{' '},num2str(kdim)));

                % Gradient delay vector from app
                dTotal(1) = app.GxDelayEditField.Value;
                dTotal(2) = app.GyDelayEditField.Value;
                dTotal(3) = app.GzDelayEditField.Value;
                app.DataOffsetRadialEditField.Value = 0;

                % Gradient delay calibration
                if app.GradDelayCalibrationCheckBox.Value

                    % Calculate k-space and trajectory of a typical frame
                    % and dynamic, remove all zero values

                    kSpacePicsSum = kSpacePics(:,:,:,:,:,:,:,:,:,:,floor(dimF/2)+1,floor(dimD/2)+1,:,floor(dimZ/2)+1);
                    trajPicsSum = trajPics(:,:,:,:,:,:,:,:,:,:,floor(dimF/2)+1,floor(dimD/2)+1,:,floor(dimZ/2)+1);

                    ze = squeeze(abs(kSpacePicsSum(1,floor(end/2),:))) > 0;
                    kSpacePicsSum = kSpacePicsSum(:,:,ze);
                    trajPicsSum = trajPicsSum(:,:,ze);

                    % Ring method
                    if app.RingMethodCheckBox.Value

                        % Sent gradient delay vector back to app
                        app.GxDelayEditField.Value = 0;
                        app.GyDelayEditField.Value = 0;
                        app.GzDelayEditField.Value = 0;

                        % Ring method using estdelay in Bart
                        try

                            dTotal = [];

                            delaysBart = bart(app,'estdelay -r4 ',trajPicsSum,kSpacePicsSum);

                            ff = strfind(delaysBart,"[0m");
                            if ~isempty(ff)
                                delaysBart = delaysBart(ff:end);
                                delaysBart = erase(delaysBart,"[0m");
                                delaysBart = erase(delaysBart,newline);
                            end

                            delaysBart = strrep(delaysBart,':',',');
                            dTotal = str2num(delaysBart);
                            dTotal(1) = -dTotal(1);

                        catch ME

                            app.TextMessage(ME.message);
                            app.TextMessage('Ring gradient delay estimation failed ...');
                            app.TextMessage('Trying iterative method ...');
                            app.SetStatus(1);
                            app.RingMethodCheckBox.Value = false;
                            dTotal = zeros(3,1);

                        end

                        % Sent gradient delay vector back to app
                        if ~isempty(dTotal)

                            app.GxDelayEditField.Value = double(round(dTotal(1),5));
                            app.GyDelayEditField.Value = double(round(dTotal(2),5));
                            app.GzDelayEditField.Value = double(round(dTotal(3),5));

                        else

                            app.TextMessage('Ring gradient delay estimation failed ...');
                            app.TextMessage('Trying iterative method ...');
                            app.SetStatus(1);
                            app.RingMethodCheckBox.Value = false;
                            dTotal = zeros(3,1);

                        end

                    end % Ring method

                    % Iterative method
                    if ~app.RingMethodCheckBox.Value

                        try

                            % Calibration size
                            kSize = [6,6];
                            kSkip = round(length(kSpacePicsSum)/3000);
                            kSkip(kSkip < 1) = 1;

                            % M1:M2 = indices in trajectory for which k-space value <= calibSize
                            [~,zm] = find(squeeze(trajPicsSum(1,:,:)) == max(trajPicsSum(:)),1,'first');
                            if isempty(zm)
                                [~,zm] = find(squeeze(trajPicsSum(2,:,:)) == max(trajPicsSum(:)),1,'first');
                            end
                            M = find(sqrt(trajPicsSum(1,:,zm).^2+trajPicsSum(2,:,zm).^2) <= calibSize(1)/2);
                            M1 = M(1);
                            M2 = M(end);

                            % Reduce size for gradient calibration
                            kTrajCalib = trajPicsSum(:,M1:M2,1:kSkip:end);
                            dataCalib = kSpacePicsSum(1,M1:M2,1:kSkip:end);
                            app.TextMessage(strcat('Calibration trajectory length = ',{' '},num2str(length(kTrajCalib))));

                            % Set values to zero first
                            dTotal = zeros(3,1);
                            kTraj = kTrajCalib;

                            % Initial image
                            imCalib = bart(app,['bart nufft -i -l 0.01 ',cSize,' -t'],kTrajCalib,dataCalib);

                            % Initialization
                            iteration = 0;
                            incre = 10;
                            kCalib = obj.fft2Dmri(imCalib);
                            wnRank = 0.4;
                            rank = floor(wnRank*prod(kSize));
                            app.TextMessage(strcat('Rank = ',{' '},num2str(rank)));

                            % Data consistency
                            y = squeeze(dataCalib);
                            xOld = kCalib;

                            % Prepare for manual stop
                            app.stopMovieFlag = false;
                            app.MovieStartButton.Enable = 'off';
                            app.MovieStopButton.Enable = 'on';

                            % Iterative method with Bart
                            while  (iteration<300) && (incre>0.001) && ~app.stopMovieFlag

                                % Iteration number
                                iteration = iteration + 1;
                                app.TextMessage(strcat('Iteration:',{' '},num2str(iteration)));

                                % Solve for X
                                rank(rank>prod(kSize)) = prod(kSize);
                                xNew = obj.lowRankThresh2D(xOld,kSize,rank);
                                rank = rank+0.2;

                                % NUFFT to get updated k-space data
                                kNew = obj.ifft2Dmri(xNew);
                                dataCalib = bart(app,'bart nufft',kTraj,kNew);
                                kNew  = reshape(dataCalib,[M2-M1+1 size(kTrajCalib,3) nrCoils]);

                                % Partial derivatives
                                [dydtx,dydty] = obj.partialDerivative2D(app,kTraj,xNew,calibSize);

                                % Direct solver
                                dydt = [real(obj.vec(dydtx)) real(obj.vec(dydty)) ; imag(obj.vec(dydtx)) imag(obj.vec(dydty))];
                                dStep = ((dydt)'*dydt)\(dydt' * [real(obj.vec(kNew - y)) ; imag(obj.vec(kNew - y))]);
                                dStep(isnan(dStep)) = 0;

                                % The accumalated delays
                                dTotal(1) = dTotal(1) + real(dStep(1));
                                dTotal(2) = dTotal(2) + real(dStep(2));
                                dTotal(dTotal > 10) = 10;
                                dTotal(dTotal < -10) = -10;

                                % Conversion criterium
                                incre = norm(real(dStep));

                                % Message
                                app.TextMessage(strcat('Estimated delays:',{' '},num2str(dTotal(1)),':',num2str(dTotal(2))));

                                % Sent gradient delay vector back to app
                                app.GxDelayEditField.Value = round(double(dTotal(1)),5);
                                app.GyDelayEditField.Value = round(double(dTotal(2)),5);
                                app.GzDelayEditField.Value = 0;

                                % Interpolation to update trajectory with new delays
                                kTraj = obj.trajInterpolation(kTrajCalib,dTotal);

                                % The new image with k-space updated for gradient delays
                                imCalib = bart(app,['bart nufft -i -l 0.01 ',cSize,' -t'],kTraj,reshape(y,[1 M2-M1+1 size(kTrajCalib,3) nrCoils]));

                                % Show image
                                im = squeeze(abs(imCalib(:,:)));
                                im = permute(im,[2 1]);
                                if app.r.PHASE_ORIENTATION
                                    im = rot90(im,-1);
                                    daspect(app.MovieFig,[1 app.r.aspectratio 1]);
                                else
                                    im = flip(im,1);
                                    daspect(app.MovieFig,[app.r.aspectratio 1 1]);
                                end
                                xlim(app.MovieFig, [0 size(im,2)+1]);
                                ylim(app.MovieFig, [0 size(im,1)+1]);
                                imshow(rot90(im),[],'Parent',app.MovieFig);

                                % Calculate k-space from new image
                                xOld = obj.fft2Dmri(squeeze(imCalib));

                            end

                        catch ME

                            app.TextMessage(ME.message);
                            app.TextMessage('Gradient delay estimation failed ...');
                            app.SetStatus(1);
                            dTotal = zeros(3,1);
                            app.GxDelayEditField.Value = double(dTotal(1));
                            app.GyDelayEditField.Value = double(dTotal(2));
                            app.GzDelayEditField.Value = double(dTotal(3));

                        end

                    end

                    % Hide stop button
                    app.MovieStartButton.Enable = 'off';
                    app.MovieStopButton.Enable = 'off';
                    app.GradDelayCalibrationCheckBox.Value = 0;

                end % Gradient calibration

                % Final gradient delay correction from optimization or values from app
                trajPics = permute(trajPics,[1 2 3 14 11 12 4 5 6 7 8 9 10 13]);
                trajPics = obj.trajInterpolation(trajPics,dTotal);
                trajPics = ipermute(trajPics,[1 2 3 14 11 12 4 5 6 7 8 9 10 13]);

                % Density correction
                if app.DensityCorrectCheckBox.Value
                    app.TextMessage('Calculating density correction ...');
                    denseOnes = ones(size(kSpacePics));
                    denseOnes = denseOnes.*avgPics; % Make sure denseOnes contains only 1's when data is available
                    denseOnes(denseOnes > 1) = 1;
                    tmpDense = bart(app,strcat('nufft -d',num2str(dimX),':',num2str(dimX),':1 -a'),trajPics,denseOnes);
                    densityPics = bart(app,'nufft ',trajPics,tmpDense);
                    densityPics = densityPics.^(-1/2);
                    densityPics(isnan(densityPics)) = 0;
                    densityPics(isinf(densityPics)) = 0;
                end

                % Sensitivity maps
                if nrCoils > 1 && ESPIRiT
                    kSpacePicsSum = sum(kSpacePics,[11,12,14]);
                    trajPicsSum = sum(trajPics,[11,12,14]);
                    ze = squeeze(abs(kSpacePicsSum(1,end,:))) > 0;
                    kSpacePicsSum = kSpacePicsSum(:,:,ze);
                    trajPicsSum = trajPicsSum(:,:,ze);
                    app.TextMessage('Calculating coil sensitivity maps ...');
                    lowResImage = bart(app,['nufft -i -l2 ',cSize,' -t'], trajPicsSum, kSpacePicsSum);
                    lowResKspace = bart(app,'fft -u 7', lowResImage);
                    kSpaceZeroFilled = bart(app,['resize -c 0 ',num2str(dimY),' 1 ',num2str(dimX)], lowResKspace);
                    sensitivities = bart(app,'ecalib -t0.002 -m1', kSpaceZeroFilled);
                else
                    sensitivities = ones(dimX,dimY,1,nrCoils,1,1,1,1,1,1,1,1,1,dimZ);
                end

                % Prepare the 2D radial PICS reconstruction
                app.TextMessage('PICS reconstruction ...');
                picsCommand = 'pics -S -i20 ';
                if Wavelet>0
                    picsCommand = [picsCommand,' -RW:6:0:',num2str(Wavelet)];
                end
                if TVxyz>0
                    picsCommand = [picsCommand,' -R',obj.totalVariation,':6:0:',num2str(TVxyz)];
                end
                if LR>0
                    % Locally low-rank in the spatial domain
                    blockSize = round(dimX/16);  % Block size
                    blockSize(blockSize<8) = 8;
                    picsCommand = [picsCommand,' -RL:6:6:',num2str(LR),' -b',num2str(blockSize),' -N '];
                end
                if TVt>0
                    picsCommand = [picsCommand,' -R',obj.totalVariation,':1024:0:',num2str(TVt)];
                end
                if TVd>0
                    picsCommand = [picsCommand,' -R',obj.totalVariation,':2048:0:',num2str(TVd)];
                end

                % Do the reco
                if app.DensityCorrectCheckBox.Value
                    igrid = bart(app,picsCommand,'-t',trajPics,'-p',densityPics,kSpacePics,sensitivities);
                else
                    igrid = bart(app,picsCommand,'-t',trajPics,kSpacePics,sensitivities);
                end

                % Root sum of squares over all coils
                recoImage = bart(app,'rss 8', igrid);

                % Rearrange to correct orientation: frames, x, y, z, dynamics
                imageReg = reshape(recoImage,[dimY,dimX,dimF,dimD,dimZ]);
                imageOut = permute(imageReg,[3,1,2,5,4]);

                % Flip for phase-orientation is vertical
                if obj.PHASE_ORIENTATION == 0
                    imageOut = flip(imageOut,3);
                end

                % Sense map orientations: x, y, slices, map1, map2
                senseMap1 = flip(permute(abs(sensitivities),[1,2,14,3,4,5,6,7,8,9,10,11,12,13]),2);

                % Normalize sense map to reasonable value range
                senseMapOut = senseMap1*4095/max(senseMap1(:));

                % Return the movie and sense-map objects
                obj.movieExp = imageOut;
                obj.senseMap = senseMapOut;

                app.validRecoFlag = true;
                app.validSenseMapFlag = true;
                app.ProgressGauge.Value = 100;
                drawnow;

            else

                % Reconstruction with NUFFT in Matlab
                app.TextMessage('BART toolbox not available ...');
                app.TextMessage('2D radial reconstruction using NUFFT ...');
                app.ProgressGauge.Value = 0;

                dimF = size(obj.kSpace{1},1);     % frames
                dimX = size(obj.kSpace{1},2);     % x
                dimY = dimX;                            % y
                dimZ = size(obj.kSpace{1},4);     % slices
                dimD = size(obj.kSpace{1},5);     % dynamics
                nrCoils = obj.nr_coils;             % coils
                loops = dimZ*dimF*dimD*nrCoils;         % number of NUFFT loops

                app.TextMessage('Slow reconstruction ...');
                app.TextMessage(strcat('Estimated reconstruction time >',{' '},num2str(loops/2),{' '},'min ...'));

                traj = obj.kSpaceTraj;                                        % trajectory
                image = zeros(dimF,dimZ,dimY,dimX,dimD,nrCoils);                    % image pre-allocation
                sensitivities = ones(dimX,dimY,1,nrCoils,1,1,1,1,1,1,1,1,1,dimZ);   % sensitivity maps = 1

                % Gradient delays from app
                dTotal(1) = app.GxDelayEditField.Value;
                dTotal(2) = app.GyDelayEditField.Value;
                dTotal(3) = app.GzDelayEditField.Value;
                app.DataOffsetRadialEditField.Value = 0;

                % Prepare the trajectory
                traj = permute(traj,[6,2,3,4,1,5]);
                traj = obj.trajInterpolation(traj,dTotal); % Shift with the gradient delay values
                traj = permute(traj,[5,3,2,4,6,1]);
                % frames, spokes, readout, slices, dynamics, 3 coordinates

                % Initialization
                maxit = 5;      % 0 or 1 for gridding, higher values for conjugate gradient
                damp = 0;       % Tikhonov penalty on ||x||
                weight = [];    % data weighting (optional)
                partial = 0.5;  % Tikhobov penalty on ||imag(x))||

                cnt = 1;

                for slice = 1:dimZ

                    for dynamic = 1:dimD

                        for frame = 1:dimF

                            % NOTE: coils can be incorporated in the NUFFT, need data first
                            for coil = 1:nrCoils

                                app.ProgressGauge.Value = round(100*cnt/loops);
                                drawnow;

                                om = permute(squeeze(traj(frame,:,:,slice,dynamic,:)),[3 2 1]);
                                obj = nufft_3d(om,dimX,app);

                                kSpaceData = squeeze(obj.kSpace{coil}(frame,:,:,slice,dynamic));
                                kSpaceData = kSpaceData(:);

                                reco = squeeze(obj.iNUFT(kSpaceData,maxit,damp,weight,'phase-constraint',partial,app));
                                reco = reshape(reco,[1,1,size(reco),1,1]);

                                image(frame,slice,:,:,dynamic,coil) = reco;

                            end

                            cnt = cnt + 1;

                        end

                    end

                end

                % Root sum of squares over coil dimension
                image = rssq(image,6);

                % ImageOut = frames, x, y, slices, dynamics
                imageOut = permute(image,[1,3,4,2,5]);

                % Flip for phase-orientation is vertical
                if obj.PHASE_ORIENTATION == 0
                    imageOut = flip(imageOut,3);
                end

                % sense map orientations: x, y, slices, map1, map2
                senseMap1 = flip(permute(abs(sensitivities),[2,1,14,3,4,5,6,7,8,9,10,11,12,13]),2);

                % normalize sense map to reasonable value range
                senseMapOut = senseMap1*4095/max(senseMap1(:));

                % Return the image and sense maps objects
                obj.movieExp = imageOut;
                obj.senseMap = senseMapOut;

                app.validRecoFlag = true;
                app.validSenseMapFlag = true;
                app.ProgressGauge.Value = 100;
                drawnow;

            end

        end % reco2Dradial




        % ---------------------------------------------------------------------------------
        % 3D UTE reconstruction with the Bart toolbox
        % ---------------------------------------------------------------------------------
        function obj = reco3Dute(obj, app)

            if app.bartDetectedFlag && ~app.MatlabRecoCheckBox.Value

                app.TextMessage('Reconstructing 3D UTE data with the BART toolbox ...');
                app.ProgressGauge.Value = 0;
                drawnow;

                % Wavelet = wavelet L1-norm regularization factor
                % TVxyz = total variation in xyz-dimension regularization
                % LR = low rank regularization
                % TVt = total variation in time regularization
                % TVd = total variation in dynamic dimension regularization
                % nrCoils = number of RF receiver coils

                Wavelet = app.WVxyzEditField.Value;
                TVxyz = app.TVxyzEditField.Value;
                LR = app.LLRxyzEditField.Value;
                TVt = app.TVcineEditField.Value;
                TVd = app.TVdynEditField.Value;
                nrCoils = obj.nr_coils;
                obj.totalVariation = 'T';

                % Retrieve the dimensions from the k-space object
                dimF = size(obj.kSpace{1},1);
                dimX = size(obj.kSpace{1},2);
                dimY = dimX;
                dimZ = dimX;
                dimd = size(obj.kSpace{1},5);

                % Retrieve k-space, trajectory, and averages
                kSpaceData = zeros(size(obj.kSpace{1}));
                for i = 1:nrCoils
                    kSpaceData(:,:,:,:,:,i) = obj.kSpace{i};
                end
                traj = obj.kSpaceTraj;
                averages = obj.kSpaceAvg;

                % Bart dimensions  Bart   Matlab
                % 	READ_DIM,       0       1   z
                % 	PHS1_DIM,       1       2   y
                % 	PHS2_DIM,       2       3   x
                % 	COIL_DIM,       3       4   coils
                % 	MAPS_DIM,       4       5   sense maps
                % 	TE_DIM,         5       6
                % 	COEFF_DIM,      6       7
                % 	COEFF2_DIM,     7       8
                % 	ITER_DIM,       8       9
                % 	CSHIFT_DIM,     9       10
                % 	TIME_DIM,       10      11  cardiac / respiratory CINE frames
                % 	TIME2_DIM,      11      12  dynamics
                % 	LEVEL_DIM,      12      13
                % 	SLICE_DIM,      13      14  slices
                % 	AVG_DIM,        14      15

                % Rearrange for BART         1  2  3  4  5  6  7  8  9 10 11 12 13 14
                kSpacePics = permute(kSpaceData,[7, 2, 3, 6,14, 8, 9,10,11,12,1, 5, 13, 4]);

                % Rearrange for BART        1  2  3  4  5  6  7  8  9 10 11 12 13 14
                avgPics = permute(averages,[7, 2, 3, 6,14, 8, 9,10,11,12,1, 5, 13, 4]);

                % Rearrange for BART     1  2  3  4  5  6  7  8  9 10 11 12 13 14
                trajPics = permute(traj,[6, 2, 3,14, 7, 8, 9,10,11,12, 1, 5,13, 4]);

                % Sum of k-space and trajectory over all frames and dynamics
                % This is the same as the originally acquired data with ommission of the resipratory windows
                kSpacePicsSum = sum(kSpacePics,[11,12]);
                trajPicsSum = sum(trajPics,[11,12]);

                % Gradient delay vector from app
                dTotal(1) = app.GxDelayEditField.Value;
                dTotal(2) = app.GyDelayEditField.Value;
                dTotal(3) = app.GzDelayEditField.Value;

                % Calibration and density correction size
                kdim = round(dimX/3);
                if mod(kdim,2) == 1
                    kdim = kdim + 1;
                end
                kdim(kdim < 32) = 32;
                kdim(kdim > dimX) = dimX;
                calibSize = [kdim, kdim, kdim];
                cSize = ['-d',num2str(calibSize(1)),':',num2str(calibSize(2)),':',num2str(calibSize(3))];
                app.TextMessage(strcat('Calibration size = ',{' '},num2str(kdim)));

                % Gradient delay calibration
                if app.GradDelayCalibrationCheckBox.Value

                    % Calibration size
                    kSize = [6,6,6];
                    kSkip = round(length(kSpacePicsSum)/3000);
                    kSkip(kSkip < 1) = 1;

                    % M = index in trajectory for which k-space value >= calibSize
                    [~,zm] = find(squeeze(trajPicsSum(1,:,:)) == max(trajPicsSum(:)),1,'first');
                    if isempty(zm)
                        [~,zm] = find(squeeze(trajPicsSum(2,:,:)) == max(trajPicsSum(:)),1,'first');
                    end
                    if isempty(zm)
                        [~,zm] = find(squeeze(trajPicsSum(3,:,:)) == max(trajPicsSum(:)),1,'first');
                    end
                    M = find(sqrt(trajPicsSum(1,:,zm).^2+trajPicsSum(2,:,zm).^2+trajPicsSum(3,:,zm).^2) >= calibSize(1,1)/2,1);

                    % Reduce size for gradient calibration
                    kTrajCalib = trajPicsSum(:,1:M,1:kSkip:end);
                    dataCalib = kSpacePicsSum(1,1:M,1:kSkip:end);
                    ze = squeeze(abs(dataCalib(1,1,:))) > 0;
                    kTrajCalib = kTrajCalib(:,:,ze);
                    app.TextMessage(strcat('Calibration trajectory length = ',{' '},num2str(length(kTrajCalib))));

                    % Interpolation to update trajectory with initial delays
                    kTrajCalib = obj.trajInterpolation(kTrajCalib,dTotal);
                    kTraj = kTrajCalib;

                    % Initial image
                    dataCalib = dataCalib(:,:,ze);
                    imCalib = bart(app,['bart nufft -i -l 0.01 ',cSize,' -t'],kTrajCalib,dataCalib);

                    % Initialization
                    iteration = 0;
                    incre = 10;
                    kCalib = obj.fft3Dmri(imCalib);
                    wnRank = 0.8;
                    rank = floor(wnRank*prod(kSize));
                    app.TextMessage(strcat('Rank = ',{' '},num2str(rank)));

                    % Data consistency
                    y = squeeze(dataCalib);
                    xOld = kCalib;

                    % Prepare for manual stop
                    app.stopMovieFlag = false;
                    app.MovieStartButton.Enable = 'off';
                    app.MovieStopButton.Enable = 'on';

                    % Calibration
                    while  (iteration<300)  && (incre>0.001) && ~app.stopMovieFlag

                        % Iteration number
                        iteration = iteration + 1;
                        app.TextMessage(strcat('Iteration:',{' '},num2str(iteration)));

                        % Solve for X
                        xNew = obj.lowRankThresh3D(xOld,kSize,rank);
                        rank = rank+0.2;
                        rank(rank>prod(kSize)) = prod(kSize);

                        % NUFFT to get updated k-space data
                        kNew = obj.ifft3Dmri(xNew);
                        dataCalib = bart(app,'bart nufft',kTraj,kNew);
                        kNew  = reshape(dataCalib,[M size(kTrajCalib,3) nrCoils]);

                        % Partial derivatives
                        [dydtx,dydty,dydtz] = obj.partialDerivative3D(app,kTraj,xNew,calibSize);

                        % Direct solver
                        dydt = [real(obj.vec(dydtx)) real(obj.vec(dydty)) real(obj.vec(dydtz)) ; imag(obj.vec(dydtx)) imag(obj.vec(dydty)) imag(obj.vec(dydtz))];
                        dStep = ((dydt)'*dydt)\(dydt' * [real(obj.vec(kNew - y)) ; imag(obj.vec(kNew - y))]);
                        dStep(isnan(dStep)) = 0;

                        % The accumalated delays
                        dTotal = dTotal + real(dStep);
                        dTotal(dTotal > 10) = 10;
                        dTotal(dTotal < -10) = -10;

                        % Conversion criterium
                        incre = norm(real(dStep));

                        % Message
                        app.TextMessage(strcat('Estimated delays:',{' '},num2str(dTotal(1)),':',num2str(dTotal(2)),':',num2str(dTotal(3))));

                        % Sent gradient delay vector back to app
                        app.GxDelayEditField.Value = double(dTotal(1));
                        app.GyDelayEditField.Value = double(dTotal(2));
                        app.GzDelayEditField.Value = double(dTotal(3));

                        % Interpolation to update trajectory with new delays
                        kTraj = obj.trajInterpolation(kTrajCalib,dTotal);

                        % The new image with k-space updated for gradient delays
                        imCalib = bart(app,['bart nufft -i -l 0.01 ',cSize,' -t'],kTraj,reshape(y,[1 M size(kTrajCalib,3) nrCoils]));

                        % Show image
                        im = squeeze(abs(imCalib(:,:,round(calibSize(3)/2),1)));
                        im = flip(flip(im,1),2);
                        if app.r.PHASE_ORIENTATION
                            im = rot90(im,-1);
                            daspect(app.MovieFig,[1 app.r.aspectratio 1]);
                        else
                            daspect(app.MovieFig,[app.r.aspectratio 1 1]);
                        end
                        xlim(app.MovieFig, [0 size(im,2)+1]);
                        ylim(app.MovieFig, [0 size(im,1)+1]);
                        imshow(rot90(im),[],'Parent',app.MovieFig);

                        % Calculate k-space from new image
                        xOld = obj.fft3Dmri(squeeze(imCalib));

                    end

                    % Hide stop button
                    app.MovieStartButton.Enable = 'off';
                    app.MovieStopButton.Enable = 'off';

                end % Gradient calibration

                % Final gradient delay correction from optimization or values from app
                trajPics = permute(trajPics,[1 2 3 4 11 12 5 6 7 8 9 10]);
                trajPics = obj.trajInterpolation(trajPics,dTotal);
                trajPics = ipermute(trajPics,[1 2 3 4 11 12 5 6 7 8 9 10]);

                % Coil sensitivities from sum of all frames and dynamics
                if nrCoils > 1
                    kSpacePicsSum = sum(kSpacePics,[11,12]);
                    trajPicsSum = sum(trajPics,[11,12]);
                    ze = squeeze(abs(kSpacePicsSum(1,end,:))) > 0;
                    kSpacePicsSum = kSpacePicsSum(:,:,ze);
                    trajPicsSum = trajPicsSum(:,:,ze);
                    app.TextMessage('Calculating coil sensitivity maps ...');
                    lowResImage = bart(app,['nufft -i ',cSize,' -t'], trajPicsSum, kSpacePicsSum);
                    lowResKspace = bart(app,'fft -u 7', lowResImage);
                    kSpaceZeroFilled = bart(app,['resize -c 0 ',num2str(dimZ),' 1 ',num2str(dimY),' 2 ',num2str(dimX)], lowResKspace);
                    sensitivities = bart(app,'ecalib -S -t0.0005 -m1', kSpaceZeroFilled);
                else
                    sensitivities = ones(dimZ,dimY,dimX,1,1,1,1,1,1,1,1,1,1,1);
                end

                % Density correction
                if app.DensityCorrectCheckBox.Value
                    app.TextMessage('Calculating density correction ...');
                    denseOnes = ones(size(kSpacePics));
                    denseOnes = denseOnes.*avgPics; % Make sure onesDense contains only 1's when data is available
                    denseOnes(denseOnes > 1) = 1;
                    denseTmp = bart(app,strcat('nufft -d',num2str(dimX),':',num2str(dimX),':',num2str(dimX),' -a'),trajPics,denseOnes);
                    densityPics = bart(app,'nufft ',trajPics,denseTmp);
                    densityPics = densityPics.^(-1/3);
                    densityPics(isnan(densityPics)) = 0;
                    densityPics(isinf(densityPics)) = 0;
                end

                % Prepare the PICS reconstruction
                app.TextMessage('PICS reconstruction ...');
                picsCommand = 'pics -i10 ';
                if Wavelet>0
                    picsCommand = [picsCommand,' -RW:7:0:',num2str(Wavelet)];
                end
                if TVxyz>0
                    picsCommand = [picsCommand,' -R',obj.totalVariation,':7:0:',num2str(TVxyz)];
                end
                if LR>0
                    % Locally low-rank in the spatial domain
                    blockSize = round(dimX/16);  % Block size
                    blockSize(blockSize<8) = 8;
                    picsCommand = [picsCommand,' -RL:7:7:',num2str(LR),' -b',num2str(blockSize)];
                end
                if TVt>0
                    picsCommand = [picsCommand,' -R',obj.totalVariation,':1024:0:',num2str(TVt)];
                end
                if TVd>0
                    picsCommand = [picsCommand,' -R',obj.totalVariation,':2048:0:',num2str(TVd)];
                end

                % Do the Bart reco
                if app.DensityCorrectCheckBox.Value
                    igrid = bart(app,picsCommand,'-t',trajPics,'-p',densityPics,kSpacePics,sensitivities);
                else
                    igrid = bart(app,picsCommand,'-t',trajPics,kSpacePics,sensitivities);
                end

                % Root sum of squares over all coils
                recoImage = bart(app,'rss 8', igrid);

                % Rearrange to correct orientation: frames, x, y, z, dynamics
                imageReg = reshape(recoImage,[dimZ,dimY,dimX,dimF,dimd]);
                imageOut = flip(flip(permute(imageReg,[4,1,2,3,5]),3),4);

                % Sense map orientations: x, y, z, map1, map2
                senseMap1 = flip(flip(permute(abs(sensitivities),[3,2,1,4,5,6,7,8,9,10,11,12,13,14]),2),3);

                % Normalize sense map to reasonable value range
                senseMap1 = senseMap1*4095/max(senseMap1(:));

                % Shift image in phase-encoding directions if needed
                obj.movieExp = circshift(imageOut,-obj.pixelshift1,3);
                obj.movieExp = circshift(imageOut,-obj.pixelshift2,4);
                obj.senseMap = circshift(senseMap1,-obj.pixelshift1,3);
                obj.senseMap = circshift(senseMap1,-obj.pixelshift2,4);

                app.validRecoFlag = true;
                app.validSenseMapFlag = true;
                app.ProgressGauge.Value = 100;
                drawnow;

            else

                app.TextMessage('BART toolbox not available ...');
                app.TextMessage('3D UTE reconstruction using 3D NUFFT ...');
                app.ProgressGauge.Value = 0;

                dimF = size(obj.kSpace{1},1);
                dimX = size(obj.kSpace{1},2);
                dimY = dimX;
                dimZ = dimX;
                dimd = size(obj.kSpace{1},5);
                nrCoils = obj.nr_coils;
                loops = dimF*dimd*nrCoils;
                app.TextMessage('Slow reconstruction ...');
                app.TextMessage(strcat('Estimated reconstruction time >',{' '},num2str(loops*2),{' '},'min ...'));

                traj = obj.kSpaceTraj;
                image = zeros(dimF,dimX,dimY,dimZ,dimd,nrCoils);
                sensitivities = 4095*ones(dimZ,dimY,dimX,1,1,1,1,1,1,1,1,1,1,1);

                maxit = 10;     % 0 or 1 for gridding, higher values for conjugate gradient
                damp = 0;       % Tikhonov penalty on ||x||
                weight = [];    % data weighting (optional)
                partial = 0.5;  % Tikhobov penalty on ||imag(x))||

                cnt = 0;

                for dynamic = 1:dimd

                    for frame = 1:dimF

                        % NOTE: coils can be incorporated in the NUFFT, need data first
                        for coil = 1:nrCoils

                            app.ProgressGauge.Value = round(100*cnt/loops);
                            drawnow;

                            om = permute(squeeze(traj(frame,:,:,1,dynamic,:)),[3,2,1]);
                            obj = nufft_3d(om,dimX,app);

                            kSpaceData = squeeze(obj.kSpace{coil}(frame,:,:,1,dynamic));
                            kSpaceData = permute(kSpaceData,[2 1]);
                            kSpaceData = kSpaceData(:);

                            image(frame,:,:,:,dynamic,coil) = obj.iNUFT(kSpaceData,maxit,damp,weight,'phase-constraint',partial,app);

                        end

                        cnt = cnt + 1;

                    end

                end

                % Root sum of squares over coil dimension and flip
                image = rssq(image,6);
                image = flip(image,3);

                % Shift image in phase-encoding directions if needed
                obj.movieExp = circshift(image,-obj.pixelshift1,3);
                obj.movieExp = circshift(image,-obj.pixelshift2,4);
                obj.senseMap = circshift(sensitivities,-obj.pixelshift1,3);
                obj.senseMap = circshift(sensitivities,-obj.pixelshift2,4);

                app.validRecoFlag = true;
                app.validSenseMapFlag = true;
                app.ProgressGauge.Value = 100;
                drawnow;

            end

        end % reco3Dute



        % ---------------------------------------------------------------------------------
        % PCA denoising
        % ---------------------------------------------------------------------------------
        function obj = PCAdenoise(obj,app)

            try

                app.TextMessage('PCA image denoising ...');

                % Image dimensions (x, y, z, frames, dynamics)
                imPCA = permute(obj.movieExp,[2,3,4,1,5]);

                % Image dimensions
                nSlices = size(imPCA,3);
                nFrames = size(imPCA,4);
                nDyn = size(imPCA,5);

                % Denoising window
                w = app.PCAwindowEditField.Value;
                window = [w w];
                if window(1) > size(imPCA,1)/2
                    window(1) = round(size(imPCA,1)/2);
                end
                if window(2) > size(imPCA,2)/2
                    window(2) = round(size(imPCA,2)/2);
                end

                % Loop over all slices, dynamics
                % Choose 2-dim image + extra dimension
                if nFrames > 1

                    for slice = 1:nSlices
                        for dyn = 1:nDyn
                            imPCA(:,:,slice,:,dyn) = denoise(double(squeeze(imPCA(:,:,slice,:,dyn))),window);
                        end
                    end

                elseif nDyn > 1

                    for slice = 1:nSlices
                        imPCA(:,:,slice,1,:) = denoise(double(squeeze(imPCA(:,:,slice,1,:))),window);
                    end

                elseif nSlices > 1

                    imPCA(:,:,:,1,1) = denoise(double(squeeze(imPCA(:,:,:,1,1))),window);

                else

                    app.TextMessage('Denoising requires multiple dimensions ...');
                    app.PCACheckBox.Value = 0;

                end

                % Return the images object
                obj.movieExp = permute(imPCA,[4,1,2,3,5]);


            catch ME

                app.TextMessage(ME.message);

            end

        end % PCAdenoise


        % ---------------------------------------------------------------------------------
        % Suppress Gibbs ringing
        % ---------------------------------------------------------------------------------
        function obj = unRing(obj,app)

            try

                app.TextMessage('Gibbs ringing suppression ...');

                im = obj.movieExp;

                % params - 3x1 array with [minW maxW nsh]
                % nsh discretization of subpixel spaceing (default 20)
                % minW  left border of window used for TV computation (default 1)
                % maxW  right border of window used for TV computation (default 3)
                params = [1 3 20];

                % image dimensions (frames, X, Y, Z, dynamics)
                nFrames = size(im,1);
                nSlices = size(im,4);
                nDyn = size(im,5);

                % unRing
                for dyn = 1:nDyn
                    for slice = 1:nSlices
                        for frame = 1:nFrames
                            for i = 1:3
                                im(frame,:,:,slice,dyn) = ringRm(double(squeeze(im(frame,:,:,slice,dyn))),params);
                            end
                        end
                    end
                end

                obj.movieExp = im;

            catch ME

                app.TextMessage(ME.message);

            end

        end % unRing




        % ---------------------------------------------------------------------------------
        % Normalize movie intensity
        % ---------------------------------------------------------------------------------
        function obj = normImages(obj)

            % Normalize the images to 2^15 range

            obj.rescaleIntercept = 0;                           % Dicom rescale intercept
            obj.rescaleSlope = max(obj.movieExp(:));        % Dicom rescale slope

            obj.movieExp = round(obj.maxImageValue*obj.movieExp/obj.rescaleSlope);
            obj.movieApp = obj.movieExp;

            if obj.PHASE_ORIENTATION == 0
                obj.movieApp = flip(obj.movieApp,2);
            end

        end % normImages



        % ---------------------------------------------------------------------------------
        % Apply sub-pixel 2D image shift (2D radial)
        % ---------------------------------------------------------------------------------
        function obj = shiftImages2D(obj, app)

            % Retrieve the in-plane image shifts
            obj = obj.get3DimageShift(obj.movieExp, app);

            % Apply the shift on sub-pixel level
            for frame = 1:size(obj.movieExp,1)
                for slice = 1:size(obj.movieExp,4)
                    for dynamic = 1:size(obj.movieExp,5)
                        obj.movieExp(frame,:,:,slice,dynamic) = retro.image2Dshift(squeeze(obj.movieExp(frame,:,:,slice,dynamic)),obj.yShift(slice),obj.xShift(slice));
                    end
                end
            end

            for slice = 1:size(obj.senseMap,3)
                for map1 = 1:size(obj.senseMap,4)
                    for map2 = 1:size(obj.senseMap,5)
                        obj.senseMap(:,:,slice,map1,map2) = retro.image2Dshift(squeeze(obj.senseMap(:,:,slice,map1,map2)),obj.yShift(slice),obj.xShift(slice));
                    end
                end
            end

        end


        % ---------------------------------------------------------------------------------
        % Apply sub-pixel 3D image shift (3D P2ROUD)
        % ---------------------------------------------------------------------------------
        function obj = shiftImages3D(obj, app)

            % Retrieve the in-plane image shifts
            obj = obj.get3DimageShift(obj.movieExp, app);

            % Apply the shift on sub-pixel level
            for frame = 1:size(obj.movieExp,1)
                for slice = 1:size(obj.movieExp,4)
                    for dynamic = 1:size(obj.movieExp,5)
                        obj.movieExp(frame,:,:,slice,dynamic) = retro.image2Dshift(squeeze(obj.movieExp(frame,:,:,slice,dynamic)),obj.yShift(1),obj.xShift(1));
                    end
                end
            end

            for frame = 1:size(obj.movieExp,1)
                for readout = 1:size(obj.movieExp,2)
                    for dynamic = 1:size(obj.movieExp,5)
                        obj.movieExp(frame,readout,:,:,dynamic) = retro.image2Dshift(squeeze(obj.movieExp(frame,readout,:,:,dynamic)),obj.zShift(1),0);
                    end
                end
            end

            for slice = 1:size(obj.senseMap,3)
                for map1 = 1:size(obj.senseMap,4)
                    for map2 = 1:size(obj.senseMap,5)
                        obj.senseMap(:,:,slice,map1,map2) = retro.image2Dshift(squeeze(obj.senseMap(:,:,slice,map1,map2)),obj.yShift(1),obj.xShift(1));
                    end
                end
            end

            for readout = 1:size(obj.senseMap,1)
                for map1 = 1:size(obj.senseMap,4)
                    for map2 = 1:size(obj.senseMap,5)
                        obj.senseMap(readout,:,:,map1,map2) = retro.image2Dshift(squeeze(obj.senseMap(readout,:,:,map1,map2)),obj.zShift(1),0);
                    end
                end
            end

        end


        % ---------------------------------------------------------------------------------
        % Determine whether movie has multiple slices and/or dynamics
        % ---------------------------------------------------------------------------------
        function obj = determineMultiDimensions(obj)

            % Multislice
            if size(obj.movieApp,4) > 1
                obj.multiSliceFlag = true;
            else
                obj.multiSliceFlag = false;
            end

            % Multidynamic
            if size(obj.movieApp,5) > 1
                obj.multiDynamicFlag = true;
            else
                obj.multiDynamicFlag = false;
            end

        end % determineMultiDimensions




        % ---------------------------------------------------------------------------------
        % Image reconstruction: SUR files
        % ---------------------------------------------------------------------------------
        function obj = recoSurFiles(obj, surpath, suffix, mrdfilename, rprfilename)

            % SUR file names
            surfiles = [surpath, suffix, '_00###.SUR'];

            % Link with the server
            m_Recon = actxserver('recon.Application');

            set(m_Recon,'Visible',1);
            set(m_Recon,'DisplayImages',1);

            % Filenames
            set(m_Recon,'DataFile',mrdfilename);
            set(m_Recon,'RPRFile',rprfilename);
            set(m_Recon,'ImageFile',surfiles);

            % Delete old SUR files
            scmd = ['del /Q ', surpath, '*.SUR'];
            system(scmd);

            % Do the reco
            invoke(m_Recon,'Run');

            % Wait for recon to complete
            while get(m_Recon,'StatusID')~=4
                if get(m_Recon,'StatusID')==5
                    break;
                end
                pause(0.1);
            end

            % Stop the link
            invoke(m_Recon,'Quit');

        end % recoSurFiles




        % ---------------------------------------------------------------------------------
        % L-curve calculation
        % ---------------------------------------------------------------------------------
        function obj = recoLcurve(obj, app)

            % Range of lamdas
            lambda = logspace(log10(app.Lambda1EditField.Value),log10(app.Lambda2EditField.Value),app.LambdaStepsEditField.Value);
            lambda = round(lambda,3);

            % Determine which parameters will be optimized
            par = [0 0 0 0 0];
            if app.LcurveWVxyzCheckBox.Value  == 1  par(1) = 1; end %#ok<*SEPEX>
            if app.LcurveTVxyzCheckBox.Value  == 1  par(2) = 1; end
            if app.LcurveLRxyzCheckBox.Value  == 1  par(3) = 1; end
            if app.LcurveTVcineCheckBox.Value == 1  par(4) = 1; end
            if app.LcurveTVdynCheckBox.Value  == 1  par(5) = 1; end

            % Intialize progress gauge
            totalNumberOfSteps = app.LambdaIterationsEditField.Value * sum(par) * length(lambda);
            app.ProgressGauge.Value = 0;
            cnt = 0;

            % Iterations
            iter = 0;
            while (iter<app.LambdaIterationsEditField.Value) && (app.stopLcurveFlag==false)
                iter = iter + 1;

                % Wavelet xyz
                if par(1)==1

                    % Reset graph
                    cla(app.LcurveFig);
                    cla(app.LcurveModelErrorFig);
                    app.LcurveFig.YLabel.String = "|| WVxyz x_{\lambda} ||_1";
                    drawnow;

                    x = [];
                    y = [];

                    i = 0;
                    while (i<length(lambda)) && (app.stopLcurveFlag==false)
                        i = i + 1;

                        % Wavelet
                        app.WVxyzEditField.Value = lambda(i);
                        app.TextMessage(strcat('L-curve, WVxyz, iteration =',num2str(iter),', lambda =',num2str(lambda(i))));

                        % CS reco
                        obj = obj.reco2D(obj, obj, app);
                        Lmovie = obj.movieExp;

                        % L1 norm wavelet
                        WVxyz = bart(app,'cdf97 6',Lmovie);
                        y(i) = norm(WVxyz(:),1);

                        % L2 norm
                        kSpaceIm = bart(app,'fft -u -i 6',Lmovie);
                        kSpaceIm = bart(app,['resize -c 1 ',num2str(size(obj.kSpace{1},2)),' 2 ',num2str(size(obj.kSpace{1},3))],kSpaceIm);

                        kSpaceSum = zeros(size(kSpaceIm));
                        for j = 1:obj.nr_coils
                            kSpaceSum = kSpaceSum + abs(obj.kSpace{j});
                        end
                        kSpaceSum = kSpaceSum/obj.nr_coils;

                        mask = kSpaceSum ~= 0;
                        kSpaceDiff = abs(kSpaceIm.*mask) - kSpaceSum;
                        x(i) = norm(kSpaceDiff(:),2);

                        % Plot the data
                        app.PlotLcurveFig(lambda,x,y);

                        % Update counter
                        cnt = cnt + 1;
                        app.ProgressGauge.Value = round(100*cnt/totalNumberOfSteps);

                    end

                end

                % Total Variation spatial domain
                if par(2)==1

                    % Reset graph
                    cla(app.LcurveFig);
                    cla(app.LcurveModelErrorFig);
                    app.LcurveFig.YLabel.String = "|| TVxyz x_{\lambda} ||_1";
                    drawnow;

                    x = [];
                    y = [];

                    i = 0;
                    while (i<length(lambda)) && (app.stopLcurveFlag==false)
                        i = i + 1;

                        % Total variation
                        app.TVxyzEditField.Value = lambda(i);
                        app.TextMessage(strcat('L-curve, TVxyz, iteration=',num2str(iter),', lambda=',num2str(lambda(i))));

                        % CS reco
                        obj = obj.reco2D(obj, obj, app);
                        Lmovie = obj.movieExp;

                        % L1 norm TV in spatial dimension
                        for dim1=1:size(Lmovie,1)
                            for dim4=1:size(Lmovie,4)
                                for dim5=1:size(Lmovie,5)
                                    for dim6=1:size(Lmovie,6)
                                        [Gx,Gy] = gradient(squeeze(Lmovie(dim1,:,:,dim4,dim5,dim6)));
                                        TVxyz(dim1,:,:,dim4,dim5,dim6) = sqrt(Gx.^2 + Gy.^2);
                                    end
                                end
                            end
                        end
                        y(i) = norm(TVxyz(:),1);

                        % L2 norm
                        kSpaceIm = bart(app,'fft -u -i 6',Lmovie);
                        kSpaceIm = bart(app,['resize -c 1 ',num2str(size(obj.kSpace{1},2)),' 2 ',num2str(size(obj.kSpace{1},3))],kSpaceIm);

                        kSpaceSum = zeros(size(kSpaceIm));
                        for j = 1:obj.nr_coils
                            kSpaceSum = kSpaceSum + abs(obj.kSpace{j});
                        end
                        kSpaceSum = kSpaceSum/obj.nr_coils;

                        mask = kSpaceSum ~= 0;
                        kSpaceDiff = abs(kSpaceIm.*mask) - kSpaceSum;
                        x(i) = norm(kSpaceDiff(:),2);

                        % Plot the data
                        app.PlotLcurveFig(lambda,x,y);

                        % Update counter
                        cnt = cnt + 1;
                        app.ProgressGauge.Value = round(100*cnt/totalNumberOfSteps);

                    end

                end

                % Low-Rank in spatial domain
                if par(3)==1

                    % Reset graph
                    cla(app.LcurveFig);
                    cla(app.LcurveModelErrorFig);
                    app.LcurveFig.YLabel.String = "|| LRxyz x_{\lambda} ||_1";
                    drawnow;

                    x = [];
                    y = [];

                    i = 0;
                    while (i<length(lambda)) && (app.stopLcurveFlag==false)
                        i = i + 1;

                        % Total variation
                        app.LLRxyzEditField.Value = lambda(i);
                        app.TextMessage(strcat('L-curve, LRxyz, iteration=',num2str(iter),', lambda=',num2str(lambda(i))));

                        % CS reco
                        obj = obj.reco2D(obj, obj, app);
                        Lmovie = obj.movieExp;

                        % L1 norm TV in spatial dimension
                        for dim1=1:size(Lmovie,1)
                            for dim4=1:size(Lmovie,4)
                                for dim5=1:size(Lmovie,5)
                                    for dim6=1:size(Lmovie,6)
                                        [Gx,Gy] = gradient(squeeze(Lmovie(dim1,:,:,dim4,dim5,dim6)));
                                        TVxyz(dim1,:,:,dim4,dim5,dim6) = sqrt(Gx.^2 + Gy.^2);
                                    end
                                end
                            end
                        end
                        y(i) = norm(TVxyz(:),1);

                        % L2 norm
                        kSpaceIm = bart(app,'fft -u -i 6',Lmovie);
                        kSpaceIm = bart(app,['resize -c 1 ',num2str(size(obj.kSpace{1},2)),' 2 ',num2str(size(obj.kSpace{1},3))],kSpaceIm);

                        kSpaceSum = zeros(size(kSpaceIm));
                        for j = 1:obj.nr_coils
                            kSpaceSum = kSpaceSum + abs(obj.kSpace{j});
                        end
                        kSpaceSum = kSpaceSum/obj.nr_coils;

                        mask = kSpaceSum ~= 0;
                        kSpaceDiff = abs(kSpaceIm.*mask) - kSpaceSum;
                        x(i) = norm(kSpaceDiff(:),2);

                        % Plot the data
                        app.PlotLcurveFig(lambda,x,y);

                        % Update counter
                        cnt = cnt + 1;
                        app.ProgressGauge.Value = round(100*cnt/totalNumberOfSteps);

                    end

                end

                % Total Variation CINE
                if par(4)==1

                    % Reset graph
                    cla(app.LcurveFig);
                    cla(app.LcurveModelErrorFig);
                    app.LcurveFig.YLabel.String = "|| TVcine x_{\lambda} ||_1";
                    drawnow;

                    x = [];
                    y = [];

                    i = 0;
                    while (i<length(lambda)) && (app.stopLcurveFlag==false)
                        i = i + 1;

                        % Total variation
                        app.TVcineEditField.Value = lambda(i);
                        app.TextMessage(strcat('L-curve, TVcine, iteration=',num2str(iter),', lambda=',num2str(lambda(i))));

                        % CS reco
                        obj = obj.reco2D(obj, obj, app);
                        Lmovie = obj.movieExp;

                        % L1 norm TV
                        TVcine = circshift(Lmovie,1) - Lmovie;
                        y(i) = norm(TVcine(:),1);

                        % L2 norm
                        kSpaceIm = bart(app,'fft -u -i 6',Lmovie);
                        kSpaceIm = bart(app,['resize -c 1 ',num2str(size(obj.kSpace{1},2)),' 2 ',num2str(size(obj.kSpace{1},3))],kSpaceIm);

                        kSpaceSum = zeros(size(kSpaceIm));
                        for j = 1:obj.nr_coils
                            kSpaceSum = kSpaceSum + abs(obj.kSpace{j});
                        end
                        kSpaceSum = kSpaceSum/obj.nr_coils;

                        mask = kSpaceSum ~= 0;
                        kSpaceDiff = abs(kSpaceIm.*mask) - kSpaceSum;
                        x(i) = norm(kSpaceDiff(:),2);

                        % Plot the data
                        app.PlotLcurveFig(lambda,x,y);

                        % Update counter
                        cnt = cnt + 1;
                        app.ProgressGauge.Value = round(100*cnt/totalNumberOfSteps);

                    end

                end

                % Total Variation dynamic dimension
                if par(5)==1

                    % Reset graph
                    cla(app.LcurveFig);
                    cla(app.LcurveModelErrorFig);
                    app.LcurveFig.YLabel.String = "|| TVdyn x_{\lambda} ||_1";
                    drawnow;

                    x = [];
                    y = [];

                    i = 0;
                    while (i<length(lambda)) && (app.stopLcurveFlag==false)
                        i = i + 1;

                        % Total variation
                        app.TVdynEditField.Value = lambda(i);
                        app.TextMessage(strcat('L-curve, TVdyn, iteration=',num2str(iter),', lambda=',num2str(lambda(i))));

                        % CS reco
                        obj = obj.reco2D(obj, obj, app);
                        Lmovie = obj.movieExp;

                        % L1 norm TV in dynamic dimension
                        TVtime = circshift(Lmovie,5) - Lmovie;
                        y(i) = norm(TVtime(:),1);

                        % L2 norm
                        kSpaceIm = bart(app,'fft -u -i 6',Lmovie);
                        kSpaceIm = bart(app,['resize -c 1 ',num2str(size(obj.kSpace{1},2)),' 2 ',num2str(size(obj.kSpace{1},3))],kSpaceIm);

                        kSpaceSum = zeros(size(kSpaceIm));
                        for j = 1:obj.nr_coils
                            kSpaceSum = kSpaceSum + abs(obj.kSpace{j});
                        end
                        kSpaceSum = kSpaceSum/obj.nr_coils;

                        mask = kSpaceSum ~= 0;
                        kSpaceDiff = abs(kSpaceIm.*mask) - kSpaceSum;
                        x(i) = norm(kSpaceDiff(:),2);

                        % Plot the data
                        app.PlotLcurveFig(lambda,x,y);

                        % Update counter
                        cnt = cnt + 1;
                        app.ProgressGauge.Value = round(100*cnt/totalNumberOfSteps);

                    end

                end


                % Analyze the L-curve and determine optimal value
                if ~app.stopLcurveFlag
                    LcurveAnalysisFnc(lambda,x,y,par);
                end


            end


            % Find the point of highest curvature
            function LcurveAnalysisFnc(lambda,x,y,par)

                % Calculate the curvature and find the maximum
                xq = log(x);
                yq = log(y);

                % Normalize curve
                xq = xq - min(xq);
                yq = yq - min(yq);
                xq = xq/max(xq);
                yq = yq/max(yq);

                curv = zeros(length(xq),1);

                for pnt = 2:length(xq)-1

                    P0 = [xq(pnt-1), yq(pnt-1), 0];
                    P1 = [xq(pnt)  , yq(pnt)  , 0];
                    P2 = [xq(pnt+1), yq(pnt+1), 0];
                    n1 = (P2 - P0) / norm(P2 - P0);
                    n2 = (P1 - P0) / norm(P1 - P0);

                    thetaInDegrees = (180/pi)*atan2(norm(cross(n1, n2)), dot(n1, n2));

                    curv(pnt) = thetaInDegrees;

                end

                % Set lambda to the value corresponding to maximum curvature
                [~,idx] = max(curv);

                [~,pidx] = max(par);

                switch pidx

                    case 1
                        app.WVxyzEditField.Value = lambda(idx);

                    case 2
                        app.TVxyzEditField.Value = lambda(idx);

                    case 3
                        app.LLRxyzEditField.Value = lambda(idx);

                    case 4
                        app.TVcineEditField.Value = lambda(idx);

                    case 5
                        app.TVdynEditField.Value = lambda(idx);

                end

                % Indicate the value in the L-curve and model-error curve
                hold(app.LcurveFig,"on");
                loglog(app.LcurveFig,x(idx),y(idx),'o','MarkerEdgeColor',app.orange,'MarkerFaceColor',app.brightred,'LineWidth',1.5);
                hold(app.LcurveFig,"off");

                hold(app.LcurveModelErrorFig,"on");
                loglog(app.LcurveModelErrorFig,[lambda(1) lambda(end)],[x(idx) x(idx)],'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',app.blue,'LineWidth',1.5);
                hold(app.LcurveModelErrorFig,"off");

                pause(1);

            end % LcurveAnalysisFnc

        end % recoLcurve





    end % Public methods










    % ---------------------------------------------------------------------------------
    % Static methods
    % ---------------------------------------------------------------------------------
    methods (Static)


        % ---------------------------------------------------------------------------------
        % 2D Tukey filter
        % ---------------------------------------------------------------------------------
        function output = circtukey2D(dimy, dimx, row, col, filterWidth)

            % Domain size is 256x256
            % Create a base matrix of all zeros.
            domain = 256;
            base = zeros(domain,domain);

            % generate a tukey window of size domain
            tukey1 = tukeywin(domain,filterWidth);

            % grab only the second half of the window
            tukey1 = tukey1(domain/2+1:domain);

            % compute the shift values for the row and column
            shifty = (row-dimy/2)*domain/dimy;
            shiftx = (col-dimx/2)*domain/dimx;

            % generate the x and y values for the domain
            y = linspace(-domain/2, domain/2, domain);
            x = linspace(-domain/2, domain/2, domain);

            % for each pixel in the domain, compute the radial distance from the center
            % and use that to index into the tukey window
            for i=1:domain
                for j=1:domain
                    rad = round(sqrt((shiftx-x(i))^2 + (shifty-y(j))^2));
                    if (rad <= domain/2) && (rad > 0)
                        base(j,i) = tukey1(rad);
                    end
                end
            end

            % resize the base window to the desired size
            output = imresize(base,[dimy dimx]);

        end % circtukey2D



        % ---------------------------------------------------------------------------------
        % 3D Tukey filter
        % ---------------------------------------------------------------------------------
        function output = circtukey3D(dimz,dimy,dimx,lev,row,col,filterWidth)

            % Domain size is 256x256x256.
            % Create a base matrix of all zeros.
            domain = 256;
            base = zeros(domain,domain,domain);

            % Generate a tukey window of size 256x1.
            tukey1 = tukeywin(domain,filterWidth);

            % Use only the second half of the window.
            tukey1 = tukey1(domain/2+1:domain);

            % Compute the shift values for the x, y, and z dimensions.
            shiftz = (lev-dimz/2)*domain/dimz;
            shifty = (row-dimy/2)*domain/dimy;
            shiftx = (col-dimx/2)*domain/dimx;

            % Generate the x, y, and z dimensions of the base matrix.
            z = linspace(-domain/2, domain/2, domain);
            y = linspace(-domain/2, domain/2, domain);
            x = linspace(-domain/2, domain/2, domain);

            % For each pixel in the base matrix, compute the distance from the center
            % pixel, and if the distance is within the filter width, set the pixel
            % value to the tukey window value.
            parfor i=1:domain
                for j=1:domain
                    for k = 1:domain
                        rad = round(sqrt((shiftx-x(i))^2 + (shifty-y(j))^2 + (shiftz-z(k))^2));
                        if (rad <= domain/2) && (rad > 0)
                            base(k,j,i) = tukey1(rad);
                        end
                    end
                end
            end

            % Resize the base matrix to the dimensions of the image.
            output = imresize3(base,[dimz dimy dimx]);

        end % circtukey3D



        % ---------------------------------------------------------------------------------
        % Centric K-space filling scheme
        % ---------------------------------------------------------------------------------
        function scheme = centricFilling(noViews2)

            ord2 = zeros(2*round(noViews2/2));

            for g= 1:round(noViews2/2) % Loop over the number of views divided by 2

                % Set the value of the gth element of ord2 to the number of views divided by 2 plus g
                ord2(2*g-1) = noViews2/2+g;
                % Set the value of the gth element of ord2 to the number of views divided by 2 minus g plus 1
                ord2(2*g) = noViews2/2-g+1;

            end

            scheme = round(ord2);

        end % centricFilling



        % ---------------------------------------------------------------------------------
        % Gauss function
        % ---------------------------------------------------------------------------------
        function y = gauss(x,s,m)

            % Gaussian function

            x = ((x-m).^2) ./ (s.^2);
            s = sqrt(2*pi) * s;
            y = exp(-x) ./ s;

        end % gauss



        % ---------------------------------------------------------------------------------
        % Fractional circshift
        % ---------------------------------------------------------------------------------
        function output = fracCircShift(input,shiftsize)
            % Shift an array on a sub-index level

            % Find the integer and fractional portions of the shift size for each dimension
            int = floor(shiftsize);     % Integer portions of shiftsize
            fra = shiftsize - int;      % Fractional portions of shiftsize

            % Get the number of dimensions of the input array
            dim = numel(shiftsize);

            % Initialize the output array as the input array
            output = input;

            % Loop through each dimension of the array
            for n = 1:numel(shiftsize)
                % Get the integer and fractional portions of the shift size for the current dimension
                intn = int(n);
                fran = fra(n);

                % Create two vectors for circular shift: one with the integer shift amount and one with the integer shift amount plus 1
                shift1 = zeros(dim,1);
                shift1(n) = intn;
                shift2 = zeros(dim,1);
                shift2(n) = intn+1;

                % Linear interpolation:
                % Perform a circular shift for both shift vectors and interpolate the results based on the fractional shift amount
                output = (1-fran)*circshift(output,shift1) + fran*circshift(output,shift2);
            end
        end % fracCircShift



        % ---------------------------------------------------------------------------------
        % peakFinder
        % ---------------------------------------------------------------------------------
        function varargout = peakFinder(x0, sel, thresh, extrema, includeEndpoints, interpolate)
            % peakFinder Find peaks and valleys in data
            %
            % This function is a simple implementation of a peak finder. It
            % finds peaks and valleys in the data, returning the indices of
            % those points. It can find local maxima and minima, or global
            % maxima and minima.
            %
            % This function is not particularly robust to noise. It is best
            % suited for finding peaks and valleys in relatively smooth data.
            %
            % Usage:
            %   ind = peakFinder(x0, sel, thresh, extrema, includeEndpoints, interpolate)
            %
            % Inputs:
            %   x0 - The data to search for peaks and valleys in.
            %   sel - The size of the window to use when searching for peaks.
            %       This value should be odd.
            %   thresh - The threshold to use when searching for peaks. This
            %       value is used to eliminate peaks that are too small to be
            %       significant. This value is a fraction of the maximum value
            %       in the data. A value of 0.1 is a good starting point.
            %   extrema - Set this to 'max' to find local maxima, or 'min' to
            %       find local minima. Set this to 'gmax' to find the global
            %       maximum, or 'gmin' to find the global minimum.
            %   includeEndpoints - Set this to true to include the two
            %       endpoints in the set of points to search for peaks. This
            %       can be useful when searching for peaks in data that has
            %       been smoothed.
            %   interpolate - Set this to true to use a parabolic fit to
            %       interpolate the location of each peak. This can be useful
            %       when the data is noisy, or when you need to find the peak
            %       location with greater accuracy than the sampling rate
            %       allows.
            %
            % Outputs:
            %   ind - The indices of the peaks and valleys in the data.

            narginchk(1, 6);
            nargoutchk(0, 2);

            s = size(x0);
            flipData =  s(1) < s(2);
            len0 = numel(x0);
            if len0 ~= s(1) && len0 ~= s(2)
            elseif isempty(x0)
                varargout = {[],[]};
                return;
            end
            if ~isreal(x0)
                x0 = abs(x0);
            end

            if nargin < 2 || isempty(sel)
                sel = (max(x0)-min(x0))/4;
            elseif ~isnumeric(sel) || ~isreal(sel)
                sel = (max(x0)-min(x0))/4;
            elseif numel(sel) > 1
                sel = sel(1);
            end

            if nargin < 3 || isempty(thresh)
                thresh = [];
            elseif ~isnumeric(thresh) || ~isreal(thresh)
                thresh = [];
            elseif numel(thresh) > 1
                thresh = thresh(1);
            end

            if nargin < 4 || isempty(extrema)
                extrema = 1;
            else
                extrema = sign(extrema(1)); % Should only be 1 or -1 but make sure
            end

            if nargin < 5 || isempty(includeEndpoints)
                includeEndpoints = true;
            end

            if nargin < 6 || isempty(interpolate)
                interpolate = false;
            end

            x0 = extrema*x0(:);                             % Make it so we are finding maxima regardless
            thresh = thresh*extrema;                        % Adjust threshold according to extrema.
            dx0 = diff(x0);                                 % Find derivative
            dx0(dx0 == 0) = -eps;                           % This is so we find the first of repeated values
            ind = find(dx0(1:end-1).*dx0(2:end) < 0)+1;     % Find where the derivative changes sign

            % Include endpoints in potential peaks and valleys as desired
            if includeEndpoints
                x = [x0(1);x0(ind);x0(end)];
                ind = [1;ind;len0];
                minMag = min(x);
                leftMin = minMag;
            else
                x = x0(ind);
                minMag = min(x);
                leftMin = min(x(1), x0(1));
            end

            % X only has the peaks, valleys, and possibly endpoints
            len = numel(x);

            if len > 2
                % Function with peaks and valleys
                % Set initial parameters for loop
                tempMag = minMag;
                foundPeak = false;

                if includeEndpoints
                    signDx = sign(diff(x(1:3)));
                    if signDx(1) <= 0
                        % The first point is larger or equal to the second
                        if signDx(1) == signDx(2) % Want alternating signs
                            x(2) = [];
                            ind(2) = [];
                            len = len-1;
                        end
                    else
                        % First point is smaller than the second
                        if signDx(1) == signDx(2) % Want alternating signs
                            x(1) = [];
                            ind(1) = [];
                            len = len-1;
                        end
                    end
                end

                % Skip the first point if it is smaller so we always start on a maxima
                if x(1) >= x(2)
                    ii = 0;
                else
                    ii = 1;
                end

                % Preallocate max number of maxima
                maxPeaks = ceil(len/2);
                peakLoc = zeros(maxPeaks,1);
                peakMag = zeros(maxPeaks,1);
                cInd = 1;

                % Loop through extrema which should be peaks and then valleys
                while ii < len

                    ii = ii+1;
                    % This is a peak
                    % Reset peak finding if we had a peak and the next peak is bigger
                    % than the last or the left min was small enough to reset.
                    if foundPeak
                        tempMag = minMag;
                        foundPeak = false;
                    end

                    % Found new peak that was lager than temp mag and selectivity larger
                    % than the minimum to its left.
                    if x(ii) > tempMag && x(ii) > leftMin + sel
                        tempLoc = ii;
                        tempMag = x(ii);
                    end

                    % Make sure we don't iterate past the length of our vector
                    if ii == len
                        break;
                        % We assign the last point differently out of the loop
                    end

                    ii = ii+1;
                    % Move onto the valley
                    % Come down at least sel from peak
                    if ~foundPeak && tempMag > sel + x(ii)
                        foundPeak = true; % We have found a peak
                        leftMin = x(ii);
                        peakLoc(cInd) = tempLoc; % Add peak to index
                        peakMag(cInd) = tempMag;
                        cInd = cInd+1;
                    elseif x(ii) < leftMin
                        % New left minima
                        leftMin = x(ii);
                    end

                end

                % Check end point
                if includeEndpoints

                    if x(end) > tempMag && x(end) > leftMin + sel
                        peakLoc(cInd) = len;
                        peakMag(cInd) = x(end);
                        cInd = cInd + 1;
                    elseif ~foundPeak && tempMag > minMag
                        % Check if we still need to add the last point
                        peakLoc(cInd) = tempLoc;
                        peakMag(cInd) = tempMag;
                        cInd = cInd + 1;
                    end

                elseif ~foundPeak

                    if x(end) > tempMag && x(end) > leftMin + sel
                        peakLoc(cInd) = len;
                        peakMag(cInd) = x(end);
                        cInd = cInd + 1;
                    elseif tempMag > min(x0(end), x(end)) + sel
                        peakLoc(cInd) = tempLoc;
                        peakMag(cInd) = tempMag;
                        cInd = cInd + 1;
                    end

                end

                % Create output
                if cInd > 1
                    peakInds = ind(peakLoc(1:cInd-1));
                    peakMags = peakMag(1:cInd-1);
                else
                    peakInds = [];
                    peakMags = [];
                end

            else

                % This is a monotone function where an endpoint is the only peak
                [peakMags,xInd] = max(x);
                if includeEndpoints && peakMags > minMag + sel
                    peakInds = ind(xInd);
                else
                    peakMags = [];
                    peakInds = [];
                end

            end

            % Apply threshold value.  Since always finding maxima it will always be
            % larger than the thresh.
            if ~isempty(thresh)
                m = peakMags>thresh;
                peakInds = peakInds(m);
                peakMags = peakMags(m);
            end

            if interpolate && ~isempty(peakMags)
                middleMask = (peakInds > 1) & (peakInds < len0);
                noEnds = peakInds(middleMask);

                magDiff = x0(noEnds + 1) - x0(noEnds - 1);
                magSum = x0(noEnds - 1) + x0(noEnds + 1)  - 2 * x0(noEnds);
                magRatio = magDiff ./ magSum;

                peakInds(middleMask) = peakInds(middleMask) - magRatio/2;
                peakMags(middleMask) = peakMags(middleMask) - magRatio .* magDiff/8;
            end

            % Rotate data if needed
            if flipData
                peakMags = peakMags.';
                peakInds = peakInds.';
            end

            % Change sign of data if was finding minima
            if extrema < 0
                peakMags = -peakMags;
                x0 = -x0;
            end

            % Plot if no output desired
            if nargout == 0
            else
                varargout = {peakInds,peakMags};
            end

        end % peakFinder



        % ---------------------------------------------------------------------------------
        % 3D partial derivatives
        % ---------------------------------------------------------------------------------
        function [dydtx,dydty,dydtz] = partialDerivative3D(app,kTraj,xNew,calibSize)

            % Get the number of coils
            nrCoils = size(xNew,4);

            % Extract k-space coordinates
            kx = squeeze(kTraj(1,:,:));
            ky = squeeze(kTraj(2,:,:));
            kz = squeeze(kTraj(3,:,:));

            % Calculate differences between adjacent k-space coordinates
            dkx = zeros(size(kx));
            dky = zeros(size(ky));
            dkz = zeros(size(kz));
            dkx(2:end,:) = kx(2:end,:)- kx(1:end-1,:);
            dky(2:end,:) = ky(2:end,:)- ky(1:end-1,:);
            dkz(2:end,:) = kz(2:end,:)- kz(1:end-1,:);

            % Perform 3D inverse FFT on input data
            xNewFFT = 1j*retro.ifft3Dmri(xNew);

            % Define arrays for NUFFT gridding
            repX = repmat([0:(calibSize(1)/2-1), 0, -calibSize(1)/2+1:-1]'/calibSize(1),[1 calibSize(1) calibSize(1) nrCoils]);
            repY = permute(repX,[2 1 3]);
            repZ = permute(repX,[3 1 2]);

            % Apply NUFFT gridding in x-direction
            tmp = xNewFFT.*repX;
            tmpCalib = bart(app,'bart nufft',kTraj,reshape(tmp,[calibSize nrCoils]));
            dydkx = reshape(tmpCalib,[size(kx) nrCoils]);

            % Apply NUFFT gridding in y-direction
            tmp = xNewFFT.*repY;
            tmpCalib = bart(app,'bart nufft',kTraj,reshape(tmp,[calibSize nrCoils]));
            dydky = reshape(tmpCalib,[size(kx) nrCoils]);

            % Apply NUFFT gridding in z-direction
            tmp = xNewFFT.*repZ;
            tmpCalib = bart(app,'bart nufft',kTraj,reshape(tmp,[calibSize nrCoils]));
            dydkz = reshape(tmpCalib,[size(kx) nrCoils]);

            % Calculate partial derivatives with respect to kx, ky, and kz
            dydtx = dydkx.*repmat(dkx,[1 1 1 nrCoils]);
            dydty = dydky.*repmat(dky,[1 1 1 nrCoils]);
            dydtz = -dydkz.*repmat(dkz,[1 1 1 nrCoils]); % negative sign used because positive doesn't converge

            % Set NaN values to 0
            dydtx(isnan(dydtx)) = 0;
            dydty(isnan(dydty)) = 0;
            dydtz(isnan(dydtz)) = 0;

        end % partialDerivative3D




        % ---------------------------------------------------------------------------------
        % 2D partial derivatives
        % ---------------------------------------------------------------------------------
        function [dydtx,dydty] = partialDerivative2D(app,kTraj,Xnew,calibSize)

            % Get the number of coils
            nrCoils = size(Xnew,3);

            % Extract k-space coordinates from the kTraj variable
            kx = squeeze(kTraj(1,:,:));
            ky = squeeze(kTraj(2,:,:));

            % Compute the difference between neighboring k-space points to get dkx and dky
            dkx = zeros(size(kx));
            dky = zeros(size(ky));
            dkx(2:end,:) = kx(2:end,:)- kx(1:end-1,:);
            dky(2:end,:) = ky(2:end,:)- ky(1:end-1,:);

            % Compute the derivative of the image along the x direction (dydkx)
            % by taking the inverse Fourier transform of the image Xnew and multiplying
            % it by a ramp filter (the frequency ramp filter is created using repmat)
            % and using the Bart toolbox to perform a non-uniform FFT (nufft)
            % on the resulting data, and reshaping the data to its original size
            tmp = (1j*retro.ifft2Dmri(Xnew).*repmat([0:(calibSize(1)/2-1), 0, -calibSize(1)/2+1:-1]'/calibSize(1),[1 calibSize(2) nrCoils]));
            tmpCalib = bart(app,'bart nufft',kTraj,reshape(tmp,[calibSize 1 nrCoils]));
            dydkx = reshape(tmpCalib,[size(kx) nrCoils]);

            % Compute the derivative of the image along the y direction (dydky)
            % in a similar way to the computation of dydkx, but using a different frequency ramp filter
            tmp = (1j*retro.ifft2Dmri(Xnew).*repmat([0:(calibSize(1)/2-1), 0, -calibSize(1)/2+1:-1]/calibSize(1),[calibSize(1) 1 nrCoils]));
            tmpCalib = bart(app,'bart nufft',kTraj,reshape(tmp,[calibSize 1 nrCoils]));
            dydky = reshape(tmpCalib,[size(kx) nrCoils]);

            % Compute the partial derivative of the image along the x direction (dydtx)
            % by multiplying dydkx by dkx, and replicating dkx along the third dimension
            % so that it can be multiplied element-wise with dydkx
            dydtx = dydkx.*repmat(dkx,[1 1 nrCoils]);

            % Compute the partial derivative of the image along the y direction (dydty)
            % by multiplying dydky by dky, and replicating dky along the third dimension
            % so that it can be multiplied element-wise with dydky
            dydty = dydky.*repmat(dky,[1 1 nrCoils]);

            % Set any NaN values in dydtx and dydty to 0
            dydtx(isnan(dydtx)) = 0;
            dydty(isnan(dydty)) = 0;

        end % partialDerivative2D



        % ---------------------------------------------------------------------------------
        % Low rank threshold 3D
        % ---------------------------------------------------------------------------------
        function Xnew = lowRankThresh3D(Xold,kSize,thresh)
            % This function performs thresholding on a 3D matrix using low-rank
            % approximation.
            %
            % Inputs:
            %   - Xold: the input 3D matrix
            %   - kSize: the size of the window used to partition the matrix into
            %     smaller blocks
            %   - thresh: the number of singular values to keep
            %
            % Output:
            %   - Xnew: the thresholded 3D matrix

            % round the threshold value to the nearest integer
            thresh = round(thresh);

            % create a vector of indices to keep the top singular values
            keep = 1:thresh;

            % Get the size of the input matrix
            [sx,sy,sz,~] = size(Xold);

            % Partition the matrix into smaller blocks using im2row3D
            tmp = retro.im2row3D(Xold,kSize);
            [tsx,tsy,Nc] = size(tmp);

            % Reshape the 3D matrix into a 2D matrix
            A = reshape(retro.im2row3D(Xold,kSize),tsx,tsy*Nc);

            % Perform SVD on the 2D matrix
            [U,S,V] = svd(A,'econ');

            % Keep only the largest singular values
            A = U(:,keep)*S(keep,keep)*V(:,keep)';

            % Reshape the 2D matrix back into a 3D matrix
            A = reshape(A,tsx,tsy*Nc);

            % Combine the smaller blocks back into the original matrix
            Xnew = retro.row2im3D(A,[sx,sy,sz,Nc],kSize);

        end % lowRankThresh3D


        % ---------------------------------------------------------------------------------
        % Low rank threshold 2D
        % ---------------------------------------------------------------------------------
        function Xnew = lowRankThresh2D(Xold,kSize,thresh)
            % This function applies a thresholding operation to a 2D image Xold to obtain
            % a low-rank approximation of Xold using singular value decomposition (SVD).
            % The output Xnew is the thresholded low-rank approximation of Xold.

            % Inputs:
            % Xold - 2D image to be thresholded
            % kSize - size of the square patch used for dividing the image into patches
            % thresh - threshold value for keeping the top singular values

            % Output:
            % Xnew - thresholded low-rank approximation of Xold

            % round the threshold value to the nearest integer
            thresh = round(thresh);

            % create a vector of indices to keep the top singular values
            keep = 1:thresh;

            % get the size and number of channels of the input image Xold
            [sx,sy,Nc] = size(Xold);

            % divide the image into patches of size kSize and stack them into a
            % matrix with each column representing a patch
            tmp = retro.im2row2D(Xold,kSize);

            % get the size of the temporary matrix tmp
            [tsx,tsy,tsz] = size(tmp);

            % reshape the matrix of patches into a matrix with each row
            % representing a patch
            A = reshape(retro.im2row2D(Xold,kSize),tsx,tsy*tsz);

            % perform SVD on matrix A and keep only the top singular values
            [U,S,V] = svd(A,'econ');
            A = U(:,keep)*S(keep,keep)*V(:,keep)';

            % reshape the thresholded matrix back into the shape of the matrix of
            % patches
            A = reshape(A,tsx,tsy,tsz);

            % combine the patches back into a 2D image
            Xnew = retro.row2im2D(A,[sx,sy,Nc],kSize);

        end % lowRankThresh2D



        % ---------------------------------------------------------------------------------
        % Trajectory interpolation
        % ---------------------------------------------------------------------------------
        function kSpaceUpdate = trajInterpolation(kSpaceOld,dShift)
            % Function that interpolates the k-space data using the specified shifts
            % in the readout direction. The resulting k-space data is returned as the
            % kSpaceUpdate variable.

            % Initialize the k-space update variable
            kSpaceUpdate = zeros(size(kSpaceOld));

            % Loop over all the dynamics, frames, slices, and spokes of the input k-space data
            for idx6 = 1:size(kSpaceOld,6) % dynamics
                for idx5 = 1:size(kSpaceOld,5) % frames
                    for idx4 = 1:size(kSpaceOld,4) % slices
                        for idx3 = 1:size(kSpaceOld,3) % spokes

                            % Perform linear interpolation of k-space data along the
                            % readout direction
                            kx = interp1((1:size(kSpaceOld,2))+dShift(1),kSpaceOld(1,:,idx3,idx4,idx5,idx6),1:size(kSpaceOld,2),'linear'); % Kx
                            ky = interp1((1:size(kSpaceOld,2))+dShift(2),kSpaceOld(2,:,idx3,idx4,idx5,idx6),1:size(kSpaceOld,2),'linear'); % Ky
                            kz = interp1((1:size(kSpaceOld,2))+dShift(3),kSpaceOld(3,:,idx3,idx4,idx5,idx6),1:size(kSpaceOld,2),'linear'); % Kz

                            % Replace NaN values with either 0 or the original k-space data
                            % value depending on the sign of the shift value
                            if dShift(1) > 0
                                kx(isnan(kx)) = 0;
                            else
                                kx(isnan(kx)) = kSpaceOld(1,isnan(kx),idx3,idx4,idx5,idx6);
                            end
                            if dShift(2) > 0
                                ky(isnan(ky)) = 0;
                            else
                                ky(isnan(ky)) = kSpaceOld(2,isnan(ky),idx3,idx4,idx5,idx6);
                            end
                            if dShift(3) > 0
                                kz(isnan(kz)) = 0;
                            else
                                kz(isnan(kz)) = kSpaceOld(3,isnan(kz),idx3,idx4,idx5,idx6);
                            end

                            % Update the k-space data with the interpolated values
                            kSpaceUpdate(1,:,idx3,idx4,idx5,idx6) = kx(:);
                            kSpaceUpdate(2,:,idx3,idx4,idx5,idx6) = ky(:);
                            kSpaceUpdate(3,:,idx3,idx4,idx5,idx6) = kz(:);

                        end % end loop over spokes
                    end % end loop over slices
                end % end loop over frames
            end % end loop over dynamics

        end % kSpaceInterpolation




        % ---------------------------------------------------------------------------------
        % image to rows 3D
        % ---------------------------------------------------------------------------------
        function res = im2row3D(im, winSize)
            % Function to convert 3D image into a matrix of patches of size winSize
            % Inputs:
            %   - im: 3D input image
            %   - winSize: patch size, [w1 w2 w3]
            % Outputs:
            %   - res: matrix of image patches

            % Get size of input image
            [sx,sy,sz,nc] = size(im);

            % Initialize output matrix with the right dimensions
            res = zeros((sx-winSize(1)+1)*(sy-winSize(2)+1)*(sz-winSize(3)+1),prod(winSize),nc);

            % Initialize count to keep track of patch number
            count=0;

            % Loop over all patches
            for z=1:winSize(3)
                for y=1:winSize(2)
                    for x=1:winSize(1)
                        count = count+1;
                        % Reshape each patch into a column vector and store it in the output matrix
                        res(:,count,:) = reshape(im(x:sx-winSize(1)+x,y:sy-winSize(2)+y,z:sz-winSize(3)+z),...
                            (sx-winSize(1)+1)*(sy-winSize(2)+1)*(sz-winSize(3)+1),nc);
                    end
                end
            end

        end % im2row3D



        % ---------------------------------------------------------------------------------
        % image to rows 2D
        % ---------------------------------------------------------------------------------
        function res = im2row2D(im, winSize)
            % This function takes a 2D image and a window size, and returns a matrix
            % where each row corresponds to a sliding window of pixels from the
            % original image.
            %
            % Inputs:
            % - im: the 2D image to be converted
            % - winSize: a 2-element vector specifying the size of the sliding window
            %
            % Outputs:
            % - res: a matrix with size (num_windows x window_size x num_channels)
            %   where num_windows is the number of sliding windows that can be
            %   extracted from the image, window_size is the product of the elements
            %   in winSize, and num_channels is the number of color channels in the
            %   image
            %

            % Get the size of the input image
            [sx,sy,sz] = size(im);

            % Allocate space for the output matrix
            res = zeros((sx-winSize(1)+1)*(sy-winSize(2)+1),prod(winSize),sz);

            % Initialize a counter for the rows of the output matrix
            count=0;

            % Loop over each pixel in the window
            for y=1:winSize(2)
                for x=1:winSize(1)

                    % Increment the row counter
                    count = count+1;

                    % Extract the pixels corresponding to the current
                    % window position and store them as a column in the
                    % output matrix
                    res(:,count,:) = reshape(im(x:sx-winSize(1)+x,y:sy-winSize(2)+y,:),...
                        (sx-winSize(1)+1)*(sy-winSize(2)+1),1,sz);

                end
            end

        end % im2row2D



        % ---------------------------------------------------------------------------------
        % rows to image 3D
        % ---------------------------------------------------------------------------------
        function [res,W] = row2im3D(mtx, imSize, winSize)

            % Extract the number of coils (channels) in the matrix
            nrCoils = size(mtx,4);

            % Extract the image size along the x, y, and z dimensions
            sx = imSize(1);
            sy = imSize(2);
            sz = imSize(3);

            % Initialize result and weight matrices
            res = zeros(imSize(1),imSize(2),imSize(3),nrCoils);
            W = res;

            % Loop over all windows of the specified size, moving one pixel at a time
            count=0;

            for z=1:winSize(3)
                for y=1:winSize(2)
                    for x=1:winSize(1)

                        % Increment the row counter
                        count = count+1;

                        % Compute the sum of matrix elements within the current window
                        % and add the resulting sum to the corresponding region in the result matrix
                        res(x:sx-winSize(1)+x,y:sy-winSize(2)+y,z:sz-winSize(3)+z,:) = res(x:sx-winSize(1)+x,y:sy-winSize(2)+y,z:sz-winSize(3)+z,:) + reshape(mtx(:,count,:,:),...
                            (sx-winSize(1)+1),(sy-winSize(2)+1),(sz-winSize(3)+1),nrCoils);

                        % Increment the weight matrix at the corresponding region
                        W(x:sx-winSize(1)+x,y:sy-winSize(2)+y,z:sz-winSize(3)+z,:) = W(x:sx-winSize(1)+x,y:sy-winSize(2)+y,z:sz-winSize(3)+z,:) + 1;

                    end
                end
            end

            % Compute the average of the result matrix using the weight matrix
            res = res./W;

        end % row2im3D


        % ---------------------------------------------------------------------------------
        % rows to image 2D
        % ---------------------------------------------------------------------------------

        function [res,W] = row2im2D(mtx,imSize, winSize)

            % Extract the number of channels (depth) in the matrix
            sz = size(mtx,3);

            % Extract the image size along the x and y dimensions
            sx = imSize(1);
            sy = imSize(2);

            % Initialize result and weight matrices
            res = zeros(imSize(1),imSize(2),sz);
            W = res;

            % Loop over all windows of the specified size, moving one pixel at a time
            count=0;
            for y=1:winSize(2)
                for x=1:winSize(1)
                    count = count+1;

                    % Compute the sum of matrix elements within the current window
                    % and add the resulting sum to the corresponding region in the result matrix
                    res(x:sx-winSize(1)+x,y:sy-winSize(2)+y,:) = res(x:sx-winSize(1)+x,y:sy-winSize(2)+y,:) + reshape(mtx(:,count,:),(sx-winSize(1)+1),(sy-winSize(2)+1),sz);

                    % Increment the weight matrix at the corresponding region
                    W(x:sx-winSize(1)+x,y:sy-winSize(2)+y,:) = W(x:sx-winSize(1)+x,y:sy-winSize(2)+y,:) + 1;
                end
            end

            % Compute the average of the result matrix using the weight matrix
            res = res./W;

        end % row2im




        % ---------------------------------------------------------------------------------
        % Vectorize
        % ---------------------------------------------------------------------------------
        function v = vec(x)

            % VEC   Vectorize.
            % VEC(X), where X is a vector, matrix, or N-D array, returns a column vector
            % Containing all of the elements of X; i.e., VEC(X)=X(:).

            v = reshape(x, numel(x), 1);

        end % vec



        % ---------------------------------------------------------------------------------
        % FFT 3D
        % ---------------------------------------------------------------------------------
        function X = fft3Dmri(x)

            X=fftshift(ifft(fftshift(x,1),[],1),1)*sqrt(size(x,1));
            X=fftshift(ifft(fftshift(X,2),[],2),2)*sqrt(size(x,2));
            X=fftshift(ifft(fftshift(X,3),[],3),3)*sqrt(size(x,3));

        end % FFT 3D



        % ---------------------------------------------------------------------------------
        % iFFT 3D
        % ---------------------------------------------------------------------------------
        function x=ifft3Dmri(X)

            x=fftshift(fft(fftshift(X,1),[],1),1)/sqrt(size(X,1));
            x=fftshift(fft(fftshift(x,2),[],2),2)/sqrt(size(X,2));
            x=fftshift(fft(fftshift(x,3),[],3),3)/sqrt(size(X,3));

        end



        % ---------------------------------------------------------------------------------
        % FFT 2D
        % ---------------------------------------------------------------------------------
        function X = fft2Dmri(x)

            X=fftshift(ifft(fftshift(x,1),[],1),1)*sqrt(size(x,1));
            X=fftshift(ifft(fftshift(X,2),[],2),2)*sqrt(size(x,2));

        end % FFT 2D



        % ---------------------------------------------------------------------------------
        % iFFT 2D
        % ---------------------------------------------------------------------------------
        function x=ifft2Dmri(X)

            x=fftshift(fft(fftshift(X,1),[],1),1)/sqrt(size(X,1));
            x=fftshift(fft(fftshift(x,2),[],2),2)/sqrt(size(X,2));

        end



        % ---------------------------------------------------------------------------------
        % Fractional 2D image shift
        % ---------------------------------------------------------------------------------
        function y = image2Dshift(im, xShift, yShift)

            % Shift in pixels, can be fractional
            H = retro.ifft2Dmri(im);

            % Create linear grid
            [xF,yF] = meshgrid(-size(im,2)/2:size(im,2)/2-1,-size(im,1)/2:size(im,1)/2-1);

            % Perform the shift
            H = H.*exp(-1i*2*pi.*(xF*xShift/size(im,2)+yF*yShift/size(im,1)));

            % Take the absolute value
            y = abs(retro.fft2Dmri(H));

        end % image2Dshift




    end % Static methods





end % retro

