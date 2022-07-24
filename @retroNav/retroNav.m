classdef retroNav
    
    % Navigator data class for retrospective app
    
    properties
      
        % Raw navigator
        navAmplitude                                            % navigator amplitude
        navPhase                                                % navigator phase (not used)
        upDown = 1                                              % flip the navigator up or down
        
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
        
    end
    
    
    % ---------------------------------------------------------------------------------
    % Public methods
    % -------------------------------------------------------------------------------
    
    methods (Access = public)
        
        
        
        % ---------------------------------------------------------------------------------
        % Object constructor
        % ---------------------------------------------------------------------------------
        function obj = retroNav
            
        end % retroNav
        
        
        
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
        % Extract the navigator signals
        % ---------------------------------------------------------------------------------
        function objNav = extractNavigator(objNav, objData)
            
            switch objData.dataType
                
                case '2D'
                    extractNavigator2D;
                    
                case {'3D','3Dp2roud'}
                    extractNavigator3D;
                    
                case '2Dms'
                    extractNavigator2Dms;
                    
                case '2Dradial'
                    extractNavigator2DRadial;

                case '3Dute'
                    extractNavigator3Dute;
                    
            end
            
            
            % ---------------------------------------------------------------------------------
            % ----- 2D single-slice data -------------------------------------------------------
            % ---------------------------------------------------------------------------------
            function amplitude = extractNavigator2D

                % Extracts the navigator data from the raw k-space data
                % Outputs a 1D array of doubles

                objNav.navAmplitude = cell(objData.nr_coils);

                for coilnr = 1:objData.nr_coils

                    % Size of the input data, 4th dimension is the readout direction which contains the navigator
                    [nrEperiments,dimz,dimy,dimx] = size(objData.data{coilnr});
                    
                    % Extract the navigator and put it in a long array
                    navdataAmplitude = reshape(permute(objData.data{coilnr},[3,2,1,4]),nrEperiments*dimy*dimz,dimx);
                    
                    if objData.nrNavPointsUsed > 1
                        
                        % Take the principal component of the data
                        data = navdataAmplitude(:,objData.primaryNavigatorPoint-objData.nrNavPointsUsed+1:objData.primaryNavigatorPoint);
                        [coeff,~,~] = pca(data);
                        dataPCA = data*coeff;
                        
                        % Take the principal component of the data
                        amplitude = abs(dataPCA(:,1))';
                        
                    else
                        
                        % Single nav point
                        amplitude = abs(navdataAmplitude(:,objData.primaryNavigatorPoint))';
                        
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
                    amplitude = amplitude * objNav.upDown;
                    
                    % Return the final nav amplitude array
                    objNav.navAmplitude{coilnr} = amplitude;
                    
                end
                
            end % extractNavigator2D
            
            
            
            
            % ---------------------------------------------------------------------------------
            % ----- 2D multi-slice data -------------------------------------------------------
            % ---------------------------------------------------------------------------------
            function amplitude = extractNavigator2Dms

                % Extracts the navigator data from the raw k-space data of multi-slice 2D data
                % Outputs a 1D array of doubles

                objNav.navAmplitude = cell(objData.nr_coils);

                for coilnr = 1:objData.nr_coils

                    % Size of the input data, 4th dimension is the readout direction which contains the navigator
                    [nrExperiments,dimz,dimy,dimx] = size(objData.data{coilnr});
                    
                    % Extract the navigator and put it in a long array
                    % Y-dimension, repetitions, slice, readout
                    navDataAmplitude = reshape(permute(objData.data{coilnr},[3,1,2,4]),nrExperiments*dimy*dimz,dimx);
                    
                    if objData.nrNavPointsUsed > 1
                        
                        % Take the principal component of the data
                        data = navDataAmplitude(:,objData.primaryNavigatorPoint-objData.nrNavPointsUsed+1:objData.primaryNavigatorPoint);
                        [coeff,~,~] = pca(data);
                        dataPCA = data*coeff;
                        
                        % Take the principal component of the data
                        amplitude = abs(dataPCA(:,1))';
                        
                    else
                        
                        % Single nav point
                        amplitude = abs(navDataAmplitude(:,objData.primaryNavigatorPoint))';
                        
                    end
                    
                    % Detrend
                    amplitude(1,:) = detrend(amplitude(1,:));
                    
                    % Make a guess whether the respiration peaks are positive or negative in the different slices
                    % This will not be needed with out-of-slice navigator
                    
                    nrElements = length(amplitude);
                    nrElementsPerSlice = round(nrElements/dimz);
                    
                    for i = 1:dimz
                        
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
                    amplitude = amplitude * objNav.upDown;
                    
                    % Return the final nav amplitude
                    objNav.navAmplitude{coilnr} = amplitude;
                    
                end
                
            end % extractNavigator2Dms
            
            

            
            % ---------------------------------------------------------------------------------
            % ----- 3D data -------------------------------------------------------------------
            % ---------------------------------------------------------------------------------
            function amplitude = extractNavigator3D

                % Extracts the navigator data from the raw k-space data
                % Outputs a 1D array of doubles

                objNav.navAmplitude = cell(objData.nr_coils);

                for coilnr = 1:objData.nr_coils

                    % Size of the input data, 4th dimension is the readout direction which contains the navigator
                    [nrRepetitions,dimz,dimy,dimx] = size(objData.data{coilnr});

                    % Extract the navigator and put it in a long array
                    navDataAmplitude = reshape(permute(objData.data{coilnr},[2,3,1,4]),nrRepetitions*dimy*dimz,dimx);
                    
                    if objData.nrNavPointsUsed > 1
                        
                        % Take the principal component of the data
                        data = navDataAmplitude(:,objData.primaryNavigatorPoint-objData.nrNavPointsUsed+1:objData.primaryNavigatorPoint);
                        [coeff,~,~] = pca(data);
                        dataPCA = data*coeff;
                        
                        % Take the principal component of the data
                        amplitude = abs(dataPCA(:,1))';
                        
                    else
                        
                        % Single nav point
                        amplitude = abs(navDataAmplitude(:,objData.primaryNavigatorPoint))';
                        
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
                    amplitude = amplitude * objNav.upDown;
                    
                    % Return the final nav amplitude
                    objNav.navAmplitude{coilnr} = amplitude;
                    
                end
                
            end % extractNavigator3D
            
            
            
            % ---------------------------------------------------------------------------------
            % ----- 2D radial data ------------------------------------------------------------
            % ---------------------------------------------------------------------------------
            function extractNavigator2DRadial

                % Extracts the navigator data from the raw k-space data
                % Outputs a 1D array of doubles

                objNav.navAmplitude = cell(objData.nr_coils);

                for coilnr = 1:objData.nr_coils

                    % Size of the input data, 4th dimension is the readout direction which contains the navigator
                    [nrRepetitions,dimz,dimy,dimx] = size(objData.data{coilnr});
                    
                    % Extract the navigator and put it in a long array
                    navDataAmplitude = reshape(permute(objData.data{coilnr},[3,2,1,4]),nrRepetitions*dimy*dimz,dimx);
                    
                    if objData.nrNavPointsUsed > 1
                        
                        range = round(objData.nrNavPointsUsed/2);
                        
                        % Take the principal component of the data
                        data = navDataAmplitude(:,objData.primaryNavigatorPoint-range:objData.primaryNavigatorPoint+range);
                        [coeff,~,~] = pca(data);
                        dataPCA = data*coeff;
                        
                        % Take the principal component of the data
                        amplitude = abs(dataPCA(:,1))';
                        phase = angle(dataPCA(:,1))';
                        
                    else
                        
                        % Single nav point
                        amplitude = abs(navDataAmplitude(:,objData.primaryNavigatorPoint))';
                        phase = angle(navDataAmplitude(:,objData.primaryNavigatorPoint))';
                        
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
                    amplitude = amplitude * objNav.upDown;
                    
                    % Return the final nav amplitude
                    objNav.navAmplitude{coilnr} = amplitude;
                    objNav.navPhase{coilnr} = phase;
                    
                end
                
            end % extractNavigatorRadial


            % ---------------------------------------------------------------------------------
            % ----- 3D UTE data ---------------------------------------------------------------
            % ---------------------------------------------------------------------------------
            function extractNavigator3Dute

                % Extracts the navigator data from the raw k-space data
                % Outputs a 1D array of doubles

                objNav.navAmplitude = cell(objData.nr_coils);

                for coilnr = 1:objData.nr_coils

                    % Size of the input data, 4th dimension is the readout direction which contains the navigator
                    [nrRepetitions,~,~,dimx] = size(objData.data{coilnr});
                    
                    % Extract the navigator and put it in a long array
                    navDataAmplitude = reshape(permute(objData.data{coilnr},[3,2,1,4]),nrRepetitions,dimx);
                    
                    if objData.nrNavPointsUsed > 1
                        
                        range = objData.nrNavPointsUsed;
                        
                        % Take the principal component of the data
                        data = navDataAmplitude(:,objData.primaryNavigatorPoint:objData.primaryNavigatorPoint+range);
                        [coeff,~,~] = pca(data);
                        dataPCA = data*coeff;
                        
                        % Take the principal component of the data
                        amplitude = abs(dataPCA(:,1))';
                        phase = angle(dataPCA(:,1))';
                        
                    else
                        
                        % Single nav point
                        amplitude = abs(navDataAmplitude(:,objData.primaryNavigatorPoint))';
                        phase = angle(navDataAmplitude(:,objData.primaryNavigatorPoint))';
                        
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
                    amplitude = amplitude * objNav.upDown;
                    
                    % Return the final nav amplitude
                    objNav.navAmplitude{coilnr} = amplitude;
                    objNav.navPhase{coilnr} = phase;
                    
                end
                
            end % extractNavigator3Dute


            
        end % extractNavigator
       
       
       
        % ---------------------------------------------------------------------------------
        % Determine the power-frequency spectrum of the navigator
        % ---------------------------------------------------------------------------------
        function  objNav = determinePowerSpectrumPCA(objNav, objData)
            
            if objData.nr_coils > 1
                
                data = zeros([length(objNav.navAmplitude{1}),objData.nr_coils]);
                for i = 1:objData.nr_coils
                    data(:,i) = objNav.navAmplitude{i};
                end
                
                % Take the principal component of the data
                [coeff,~,~] = pca(data);
                dataPCA = data*coeff;
                amplitude = dataPCA(:,1);
                
            else
                
                amplitude = objNav.navAmplitude{1}';
                
            end

            % Include only those navigators when includewindow == 1 and excludewindow == 1
            amplitude = amplitude.*objData.includeWindow.*objData.excludeWindow;

            % Determine the frequency power spectrum
            y = fft(amplitude);
            fs = 1000/objData.TR;                       % Sample frequency in Hz
            n = length(amplitude);                      % Number of samples
            objNav.frequency = (0:n-1)*(fs/n)*60;       % Frequency range in bpm
            power = abs(y).^2/n;                        % Power of the DFT

            % Determine frequency and harmonics of k-space trajectory and set those to zero
            if objData.NO_VIEWS > 1

                kfreq = 60/(0.001*objData.NO_VIEWS*objData.TR);
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
            objNav.powerSpectrum = power;
            
            % Detect heart rate
            minheartbpm = objNav.physioFilterSettings(1);
            maxheartbpm = objNav.physioFilterSettings(2);
            minidx = round(minheartbpm*n/(fs*60));
            maxidx = round(maxheartbpm*n/(fs*60));
            [~, idx] = max(power(minidx:maxidx));
            objNav.detectedHR = round(idx*fs*60/n + minheartbpm);
            
            % Detect respiratory rate
            minRRbpm = objNav.physioFilterSettings(3);
            maxRRbpm = objNav.physioFilterSettings(4);
            minidx = round(minRRbpm*n/(fs*60));
            maxidx = round(maxRRbpm*n/(fs*60));
            [~, idx] = max(power(minidx:maxidx));
            objNav.detectedRR = round(idx*fs*60/n + minRRbpm);
            
        end % determinePowerSpectrumPCA
        
        
        
        % ---------------------------------------------------------------------------------
        % Filter the navigator
        % ---------------------------------------------------------------------------------
        function objNav = filterNavPCA(objNav, objData)
            
            % Applies a bandwidth filter on the navigator data
            
            sf = 1000/objData.TR;                   % Sampling frequency in Hz = 1/TR[ms]
            respHarmonics = 2;                      % Number of higher order harmonics for respiratory frequency, 2 = main + 1 harmonic
            order = objNav.physioFilterSettings(5); % Butterworth filter order
            
            if objData.nr_coils > 1
                
                data = zeros([length(objNav.navAmplitude{1}),objData.nr_coils]);
                for i = 1:objData.nr_coils
                    data(:,i) = objNav.navAmplitude{i};
                end
                
                % Take the principal component of the data
                [coeff,~,~] = pca(data);
                dataPCA = data*coeff;
                
                % Take the principal component of the data
                amplitude = dataPCA(:,1);
                
            else
                
                amplitude = objNav.navAmplitude{1}';
                
            end
            
            % Filter for heart motion
            hrf = objNav.detectedHR/60;       % Expected heartrate in Hz = hr[bpm]/60
            bwh = objNav.bandwidthHR/60;      % Bandwidth heartrate in Hz = [bpm]/60
            [b,a] = butter(order,[hrf-0.5*bwh,hrf+0.5*bwh]/(sf/2),'bandpass');     % Butterworth bandpass filter
            heartOutputData = filtfilt(b,a,amplitude);                             % Apply zero-phase shift filtering
            
            % Detrend
            heartOutputData = detrend(heartOutputData);
            
            % Normalize envelope
            factor = round(100/hrf);      % Adapt the envelope setting to match the expected heart rate frequency
            [env,~] = envelope(heartOutputData,factor,'peak');
            objNav.heartNav = heartOutputData./abs(env);
            
            % Filter for respiration motion
            while true
                
                rrf = objNav.detectedRR/60;       % Expected resprate in Hz = rr[bpm]/60
                bwr = objNav.bandwidthRR/60;      % Bandwidth resprate in Hz = [bpm]/60
                
                respOutputData = zeros(size(amplitude));
                
                if objNav.detectedRR<45
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
                
                objNav.detectedRR = objNav.detectedRR+1;
                objNav.bandwidthRR = objNav.bandwidthRR+1;
                
            end
            
            % Detrend and normalize envelope
            respOutputData = detrend(respOutputData);
            
            % Normalize envelope
            factor = round(150/rrf);  % Adapt the envelope setting to match the expected respiration rate frequency
            [env,~] = envelope(respOutputData,factor,'peak');
            objNav.respNav = respOutputData./abs(env);
            
            % Return the principal component, or in case 1 coil the orignal navigator data, and detrend
            objNav.PCANav = detrend(amplitude);
            
        end % filterNavPCA
        
        
        
        % ---------------------------------------------------------------------------------
        % Heart filter shape for plotting
        % ---------------------------------------------------------------------------------
        function output = filterShape(objNav, objData, app, type)
            
            % Determines the filter shape for display purposes
            % Order = order of the filter
            % Type = 'heart' or 'resp'
            
            order = app.FilterOrderEditField.Value;
            sf = 1000/objData.TR;           % Sampling frequency in Hz = 1/TR[ms]
            
            if strcmp(type,'resp')
                hr = objNav.detectedRR;
                bw = objNav.bandwidthRR;
            else
                hr = objNav.detectedHR;
                bw = objNav.bandwidthHR;
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
        function objNav = trigPointsHeart(objNav, objData)
            
            % Extracts the ECG trigger points from the navigators
            
            % Find the peaks and locations = fast
            try
                objNav.heartTrigPoints = retroNav.peakFinder(objNav.heartNav',[],[],[],false,true);
            catch
                objNav.heartTrigPoints = []; % Failure
            end

            % Trigger points are in units of samples (actual time is thus heartTrigPoints*TR)
         
            % Backup plan in case peakFinder fails = slower
            if length(objNav.heartTrigPoints)<20
                
                % Minimal distance 50% of expected heart rate [in points]
                if(isnan(objNav.meanHeartRate))
                    objNav.meanHeartRate = 500;
                end
                dist = 0.50*(60/objNav.meanHeartRate)/(objData.TR/1000);
                interPolationFactor = objNav.splineFactor;
                
                % Cubic spline interpolation of the data
                nrl = length(objNav.heartNav);
                navi = interp1(1:nrl,objNav.heartNav(1:nrl),1:1/interPolationFactor:nrl,'spline');
                
                % Find the peaks and locations
                [~,locs]=findpeaks(navi,'MinPeakDistance',dist*interPolationFactor);
                locs = locs + interPolationFactor/2;
                
                % Recalculate orginal fractional time point peak positions
                objNav.heartTrigPoints = locs/interPolationFactor;
                
            end

            objNav.heartTrigTime = objNav.heartTrigPoints*objData.TR - objData.TR; % Time starts at zero
            
        end % trigPointsHeart
        
        
        
        % ---------------------------------------------------------------------------------
        % Determine the respiration trigger points
        % ---------------------------------------------------------------------------------
        function objNav = trigPointsResp(objNav, objData)
            
            % Extracts the ECG trigger points from the navigators
            
            % Find the peaks and locations
            try
                objNav.respTrigPoints = retroNav.peakFinder(objNav.respNav',[],[],[],false,true);
            catch
                objNav.respTrigPoints = []; % Failure
            end

            % Backup plan in case peakFinder fails
            if length(objNav.respTrigPoints)<10
                
                interPolationFactor = objNav.splineFactor;
                nrl = size(objNav.respNav,1);
                navi = interp1(1:nrl,objNav.respNav,1:1/interPolationFactor:nrl,'spline');
                
                % Find the peaks and locations
                [~,locs] = findpeaks(navi,'MinPeakProminence',0.1);
                locs = locs + interPolationFactor/2;
                
                % Recalculate orginal fractional time point peak positions
                objNav.respTrigPoints = locs/interPolationFactor;
                         
            end
         
            % Trigger points are in units of samples (actual time is thus heartTrigPoints*TR - TR)
            objNav.respTrigTime = objNav.respTrigPoints*objData.TR  - objData.TR; % Time starts at zero

        end % trigPointsResp
        
        
        
        % ---------------------------------------------------------------------------------
        % Determine average heart rate
        % ---------------------------------------------------------------------------------
        function objNav = calcHeartRate(objNav, objData)
            
            % Determine heart rate as function of time in bpm
            
            nrLocs = size(objNav.heartTrigPoints,2);
            
            rate = zeros(nrLocs-1,1);
            for i=1:nrLocs-1
                rate(i) = 60/((objNav.heartTrigPoints(i+1)-objNav.heartTrigPoints(i))*(objData.TR/1000));
            end
            
            hr = rate';
            hrf = movmedian(hr,32); % Smooth trendline
            
            % Determine the mean cardiac rate
            includeData = objData.includeWindow.*objData.excludeWindow;                         % Data-window which is included
            try
                includeData = round(resample(includeData,2*length(hr),length(includeData)));    
            catch
                includeData = round(resample(includeData,3*length(hr),length(includeData)));    % Resample the data-window to nr samples heartrate
            end
            includeData = round(resample(includeData,length(hr),length(includeData)));          % In 2 steps, to prevent an overflow error
            
            objNav.meanHeartRate = round(median(nonzeros(includeData'.*hr)));                   % Take the median of the heartrate
            objNav.heartRateTime = hr;
            objNav.heartRateTimeFiltered = hrf;
            
        end % calcHeartRate
        
        
        
        % ---------------------------------------------------------------------------------
        % Determine average respiration rate
        % ---------------------------------------------------------------------------------
        function objNav = calcRespRate(objNav, objData)
            
            % Determine respiration rate as function of time in bpm
            
            nrLocs = size(objNav.respTrigPoints,2);
            
            rate = zeros(nrLocs-1,1);
            for i=1:nrLocs-1
                rate(i) = 60/((objNav.respTrigPoints(i+1)-objNav.respTrigPoints(i))*(objData.TR/1000));
            end
            
            resp = rate';
            respf = movmedian(resp,32);     % Smooth trendline

            % Determine the mean respiration rate
            includeData = objData.includeWindow.*objData.excludeWindow;                         % Data-window which is included
            try
                includeData = round(resample(includeData,2*length(resp),length(includeData)));  % Resample the data-window to nr samples respiration rate
            catch
                includeData = round(resample(includeData,3*length(resp),length(includeData)));  % Resample the data-window to nr samples respiration rate
            end
            includeData = round(resample(includeData,length(resp),length(includeData)));        % In 2 steps, to prevent an overflow error

            objNav.meanRespRate = round(median(nonzeros(includeData'.*resp)));                  % Take the median of the respirationrate
            objNav.respRateTime = resp;
            objNav.respRateTimeFiltered = respf;
            
        end % calcRespRate
        
        
        
        % ---------------------------------------------------------------------------------
        % Determine respiration windows
        % ---------------------------------------------------------------------------------
        function objNav = makeRespWindow(objNav, app, objData)
            
            rLocs = objNav.respTrigPoints;
            meanResp = mean(objNav.respRateTimeFiltered);
            respPercent = objNav.respPercentage;
            nrKlines = objData.nrKlines;
            objNav.respWindowTime = [];
            
            % This function creates an array (time line) of rectangular boxes of 0's and 1's around the detected respiratory signals
            % Later on 1 means that there is a respiration, for which the k-lines will be discarded
            
            respWin = 0.5*(respPercent/100)*(60/meanResp)*1000/objData.TR;   % 1/2 window width around respiratory peak locations
            
            window = zeros(nrKlines,1);   % Fill array with zeros
            
            for i = 1 : size(rLocs,2)
                
                % Center of the respiration window
                center = rLocs(1,i);
                
                % Beginning and end of the respiration window in ms, time starts at zero
                objNav.respWindowTime(i,1) = (center - respWin)*objData.TR  - objData.TR;
                objNav.respWindowTime(i,2) = (center + respWin)*objData.TR  - objData.TR;

                % Set respiration window mask to 0 during respiration
                for j = round(center - respWin) : round(center + respWin)
                    if (j>0) && (j<=nrKlines)
                        window(j) = 1;
                    end
                end

            end
            
            % Set window to expiration (0) or inspiration (1)
            if app.RespirationToggleCheckBox.Value == 1
                window = 1 - window;
            end

            objNav.respWindow = window;
            
        end % makeRespWindow
        
        
        
    end % methods (public)
        
    
    
    
    % ---------------------------------------------------------------------------------
    % Static methods
    % ---------------------------------------------------------------------------------
    
    methods (Static)
        

        % ---------------------------------------------------------------------------------
        % peakFinder
        % ---------------------------------------------------------------------------------
    
        function varargout = peakFinder(x0, sel, thresh, extrema, includeEndpoints, interpolate)
            
            % This function finds peaks in a signal

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
                    if signDx(1) <= 0 % The first point is larger or equal to the second
                        if signDx(1) == signDx(2) % Want alternating signs
                            x(2) = [];
                            ind(2) = [];
                            len = len-1;
                        end
                    else % First point is smaller than the second
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
                    elseif x(ii) < leftMin % New left minima
                        leftMin = x(ii);
                    end

                end
                
                % Check end point
                if includeEndpoints

                    if x(end) > tempMag && x(end) > leftMin + sel
                        peakLoc(cInd) = len;
                        peakMag(cInd) = x(end);
                        cInd = cInd + 1;
                    elseif ~foundPeak && tempMag > minMag % Check if we still need to add the last point
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
                x0 = -x0; %#ok<NASGU> 
            end
            
            % Plot if no output desired
            if nargout == 0
            else
                varargout = {peakInds,peakMags};
            end
            
        end % peakFinder
        


        
    end % methods (static)
    
    
end % retroNav