classdef retroKspace
    
    % K-space data and sorting class for retrospective app
    
    properties
        
        % K-space data
        raw
        kSpace
        kSpaceMrd
        kSpaceAvg
        kSpaceTraj
        trajectory
        gradTrajectory
        
        % Cardiac and respiratory binning
        cardBins
        respBins
        cardBinNrs
        respBinNrs
        
    end
    
    
    
    methods
        
        
        % ---------------------------------------------------------------------------------
        % Object constructor
        % ---------------------------------------------------------------------------------
        function obj = retroKspace
            
        end % retroKspace
        
        
        
        % ---------------------------------------------------------------------------------
        % Extract unsorted k-data from raw data
        % ---------------------------------------------------------------------------------
        function objKspace = extractData(objKspace, objData)
            
            objKspace.raw = cell(objData.nr_coils);
            
            switch objData.dataType
                
                case {'2D','2Dms','3D','3Dp2roud'}
                    for i = 1:objData.nr_coils
                        % Cut off the navigator and spacer
                        objKspace.raw{i} = objData.data{i}(:,:,:,objData.primaryNavigatorPoint+objData.nrNavPointsDiscarded+1:end);
                    end
                    
                case {'2Dradial','3Dute'}
                    for i = 1:objData.nr_coils
                        objKspace.raw{i} = objData.data{i};
                    end
                    
            end
            
        end % extractData
        
        
        
        % ---------------------------------------------------------------------------------
        % Assign bin times, relative binning
        % ---------------------------------------------------------------------------------
        function objKspace = assignBinTimes(objKspace, objNav, app)
            
            cardLocations = objNav.heartTrigPoints;
            respLocations = objNav.respTrigPoints;
            nrCardFrames = app.nrCardFrames;
            nrRespFrames = app.nrRespFrames;
            
            nrCard = length(cardLocations); % Number of cardiac time points
            nrResp = length(respLocations); % Number of respiratory time points
            
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
            
            cBins = zeros(1,(nrCard-1)*nrCardFrames);
            
            cnt = 1;
            for i=1:nrCard-1
                for j=1:nrCardFrames
                    cBins(cnt)=cardLocations(i)+(j-1)*(cardLocations(i+1)-cardLocations(i))/nrCardFrames;
                    cnt = cnt + 1;
                end
            end
            
            % Respiratory binning
            
            rBins = zeros(1,(nrResp-1)*nrRespFrames);
            
            cnt = 1;
            for i=1:nrResp-1
                for j=1:nrRespFrames
                    rBins(cnt)=respLocations(i)+(j-1)*(respLocations(i+1)-respLocations(i))/nrRespFrames;
                    cnt = cnt + 1;
                end
            end
            
            objKspace.cardBins = cBins;
            objKspace.respBins = rBins;
            
        end % assignBinTimes
        
        
        
        % ---------------------------------------------------------------------------------
        % Assign bin times, absolute binning
        % ---------------------------------------------------------------------------------
        function objKspace = assignBinTimesAbs(objKspace, objNav, objData, app)
            
            cardLocations = objNav.heartTrigPoints;
            respLocations = objNav.respTrigPoints;
            heartRate = objNav.meanHeartRate;
            TR = objData.TR;
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
            frameDuration = (60/heartRate)*(1000/TR)/nrCardFrames;
            
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
            
            objKspace.cardBins = cBins;
            objKspace.respBins = rBins;
            
        end % assignBinTimesAbs
        
        
        
        % ---------------------------------------------------------------------------------
        % Assign bin frames, relative binning
        % ---------------------------------------------------------------------------------
        function objKspace = assignBinFrames(objKspace, objNav, objData, app)
            
            % Assigns all the measured k-lines to a specific cardiac phase and respiratory phase bin
            
            binTimesCard = objKspace.cardBins';
            binTimesResp = objKspace.respBins';
            respWindow = objNav.respWindow;
            nrKlines = objData.nrKlines;
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

            cardAssignments = zeros(nrKlines,1);                            % Zero = no assignment (breathing, begin/end of data)
            parfor i=startPoint:endPoint                                    % Start search to which heartbeat the measurement belongs
                j=locMaxCard;
                while j>1 && i<binTimesCard(j)  %#ok<*PFBNS>
                    j=j-1;
                end
                cardAssignments(i) = mod(j-1,nrCardFrames)+1;            % Assign to bin frame number = j modulus nr_frames
                if nrRespFrames==1 && respWindow(i)==1                   % If measurement is during respiration and only 1 resp state, put back to 0 to discard this k-line
                    cardAssignments(i) = 0;
                end
            end

            respAssignments = zeros(nrKlines,1);                        % Zero = no assignment (breathing, begin/end of data)
            parfor i=startPoint:endPoint                                % Start search to which respiration the measurement belongs
                j=locMaxResp;
                while j>1 && i<binTimesResp(j)
                    j=j-1;
                end
                respAssignments(i) = mod(j-1,nrRespFrames)+1;        % Assign to bin frame number = j modulus nr_frames
            end
            
            % Report back
            objKspace.cardBinNrs = cardAssignments;
            objKspace.respBinNrs = respAssignments;
            
        end % assignBinFrames
        
        
        
        % ---------------------------------------------------------------------------------
        % Assign bin frames, absolute binning
        % ---------------------------------------------------------------------------------
        function objKspace = assignBinFramesAbs(objKspace, objNav, objData, app)
            
            binTimesCard = objKspace.cardBins';
            binTimesResp = objKspace.respBins';
            respWindow = objNav.respWindow;
            nrKlines = objData.nrKlines;
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
           
            cardAssignments = zeros(nrKlines,1);                            % Zero = no assignment (breathing, begin/end of data)
            parfor i=startPoint:endPoint                                    % Start search to which heartbeat the measurement belongs
                j=locMaxCard;
                
                while j>1 && i<binTimesCard(j,1)
                    j=j-1;
                end
                
                cardAssignments(i) = binTimesCard(j,2);                     % Assign to bin frame number
                
                if nrRespFrames==1 && respWindow(i)==1
                    cardAssignments(i) = 0;                                 % If measurement is during respiration and only 1 resp state, put back to 0 to discard this k-line
                end
                
            end
            
            % Respiratory binning is still done using relative binning
            
            respAssignments = zeros(nrKlines,1);                            % Zero = no assignment (breathing, begin/end of data)
            parfor i=startPoint:endPoint                                    % Start search to which respiration the measurement belongs
                j=locMaxResp;
                while j>1 && i<binTimesResp(j)
                    j=j-1;
                end
                respAssignments(i) = mod(j-1,nrRespFrames)+1;               % Assign to bin frame number = j modulus nr_frames
            end
            
            objKspace.cardBinNrs = cardAssignments;
            objKspace.respBinNrs = respAssignments;
            
        end % assignBinFramesAbs
        
        
        
        % ---------------------------------------------------------------------------------
        % Fill K-space 2D
        % ---------------------------------------------------------------------------------
        function objKspace = fillKspace2D(objKspace, objData, app)

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
            objKspace.kSpace = cell(objData.nr_coils);
            objKspace.kSpaceAvg = [];
            includeWindow = objData.includeWindow.*objData.excludeWindow;
            cardBinAss = objKspace.cardBinNrs;
            respBinAss = objKspace.respBinNrs;
            dimx = objData.dimx;
            dimy = objData.dimy;
            dimz = objData.dimz;
            nrKsteps = objData.nrKsteps;
            traj = objKspace.trajectory;
            nrCardFrames = app.nrCardFrames;
            nrRespFrames = app.nrRespFrames;
            nrDynamics = app.nrDynamics;
 
            for coilnr = 1:objData.nr_coils
                
                rawdata = objKspace.raw{coilnr};
           
                nrReps = size(rawdata,1);                                       % Number of k-space repetitions
                unsortedKspace = reshape(rawdata,[1,size(rawdata),1]);
                
                % Dynamics assignment
                totalk = nrReps * nrKsteps * dimz;
                dynBinAss = round(linspace(0.5, nrDynamics+0.49, totalk));      % List of increasing integer number 1 .. nr_dynamics evenly spaced over the entire acquistion time

                % Sorting
                sortedKspace = complex(zeros(nrRespFrames,nrCardFrames,dimz,dimy,dimx,nrDynamics));   % Fill temp k-space with zeros
                sortedAverages = zeros(nrRespFrames,nrCardFrames,dimz,dimy,dimx,nrDynamics);          % Fill temp nr averages array with zeros
                cnt = 0;
                for slice=1:dimz                    % Loop over slices
                    for i=1:nrReps                  % Loop through all repetitions
                        for j=1:nrKsteps            % Loop through all the phase-encoding steps
                            cnt = cnt + 1;
                            if (cardBinAss(cnt) > 0) && (includeWindow(cnt) == 1)       % If assigment = 0, this acquisition is discarded
                                kline = traj(mod(cnt - 1,nrKsteps) + 1);                % The phase-encoding step using the trajectory info
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
                    Ry = round(dimy/share/2);
                    Rx = round(dimx/share/2);
                    [Y,X] = ndgrid(1:dimy,1:dimx);
                    L = zeros(share,dimy,dimx);
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
                            weights(i,j) = retroKspace.gauss(i+j-1,share,0);
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
                                    
                                    ROI = reshape(squeeze(C(j,:,:)),[1 1 1 dimy dimx 1])*weights(j,abs(i));
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
                tukeyFilter(1,1,1,:,:,1) = retroKspace.circtukey2D(dimy,dimx,row,col,filterWidth);
                tmpKspace = tmpKspace.*tukeyFilter;

                % Report back k-space per coil
                objKspace.kSpace{coilnr} = tmpKspace;

            end

            % Report back averages
            objKspace.kSpaceAvg = tmpAverages;

        end % fillKspace2D
        
        
        
        % ---------------------------------------------------------------------------------
        % Fill K-space 3D
        % ---------------------------------------------------------------------------------
        function objKspace = fillKspace3D(objKspace, objData, app)

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
            objKspace.kSpace = cell(objData.nr_coils);
            objKspace.kSpaceAvg = [];
            includeWindow = objData.includeWindow.*objData.excludeWindow;
            cardBinAss = objKspace.cardBinNrs;
            respBinAss = objKspace.respBinNrs;
            dimx = objData.dimx;
            dimy = objData.dimy;
            dimz = objData.dimz;
            nrKsteps = objData.nrKsteps;
            traj = objKspace.trajectory;
            nrCardFrames = app.nrCardFrames;
            nrRespFrames = app.nrRespFrames;
            nrDynamics = app.nrDynamics;

            for coilnr = 1:objData.nr_coils

                rawdata = objKspace.raw{coilnr};

                sortedKspace = complex(zeros(nrRespFrames,nrCardFrames,dimz,dimy,dimx,nrDynamics));     % Fill temp k-space with zeros
                sortedAverages = zeros(nrRespFrames,nrCardFrames,dimz,dimy,dimx,nrDynamics);            % Fill temp nr averages array with zeros
                nrReps = size(rawdata,1);                                                               % Number of k-space repetitions
                unsortedKspace = reshape(rawdata,[1,size(rawdata),1]);

                % Dynamics assignment
                totalk = nrReps * nrKsteps * dimz;
                dynBinAss = round(linspace(0.5, nrDynamics+0.49, totalk));       % List of increasing integer number 1 .. nr_dynamics evenly spaced over the entire acquistion time

                % Adapt trajectory for 3D acqusition

                % The y-dimension
                traj3Dy = zeros(nrReps * nrKsteps * dimz,1);
                cnt = 1;
                for i = 1:nrReps
                    for j = 1:nrKsteps
                        for k = 1:dimz
                            traj3Dy(cnt) = traj(j);
                            cnt = cnt + 1;
                        end
                    end
                end

                % The z-dimension
                traj3Dz = zeros(nrReps * nrKsteps * dimz,1);

                switch objData.pe2_centric_on

                    case 0

                        % Linear in the 3rd dimension
                        cnt = 1;
                        for i = 1:nrReps
                            for j = 1:nrKsteps
                                for k = 1:dimz
                                    traj3Dz(cnt) = k;  % Linear
                                    cnt = cnt + 1;
                                end
                            end
                        end

                    case 1

                        % Centric in the 3rd dimension
                        cnt = 1;
                        cf = retroKspace.centricFilling(dimz);
                        for i = 1:nrReps
                            for j = 1:nrKsteps
                                for k = 1:dimz
                                    traj3Dz(cnt) = cf(k);  % Centric
                                    cnt = cnt + 1;
                                end
                            end
                        end

                    case 2

                        % Special case, trajectory
                        cnt = 1;
                        for i = 1:nrReps
                            for j = 1:nrKsteps
                                for k = 1:dimz
                                    traj3Dz(cnt) = objData.pe2_traj(k) + round(dimz/2) + 1;
                                    cnt = cnt + 1;
                                end
                            end
                        end

                end

                % Do the filling of k-space
                cnt = 0;

                for i = 1:nrReps                    % Loop through all repetitions

                    for j = 1:nrKsteps              % Loop through all the phase-encoding steps

                        for k = 1:dimz              % Loop through phase-encoding 3rd dimension

                            cnt = cnt + 1;

                            if (cardBinAss(cnt) > 0) && (includeWindow(cnt) == 1)     % If assigment == 0, this acquisition is discarded

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
                    weights = retroKspace.gauss(1:share+1,share,0);
                    weights = weights/max(weights);

                    % Define ellipsoid regions
                    Rz = round(dimz/share/2);
                    Ry = round(dimy/share/2);
                    Rx = round(dimx/share/2);
                    [Z,Y,X] = ndgrid(1:dimz,1:dimy,1:dimx);
                    L = zeros(share,dimz,dimy,dimx);
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
                            weights(i,j) = retroKspace.gauss(i+j-1,share,0);
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

                                    ROI = reshape(squeeze(C(j,:,:,:)),[1 1 dimz dimy dimx 1])*weights(j,abs(i));
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
                tukeyFilter(1,1,:,:,:,1) = retroKspace.circtukey3D(dimz,dimy,dimx,lev,row,col,filterWidth);
                tmpKspace = tmpKspace.*tukeyFilter;

                % Report back
                objKspace.kSpace{coilnr} = tmpKspace;

            end

            objKspace.kSpaceAvg = tmpAverages;

        end % fillKspace3D
        
        
        
        % ---------------------------------------------------------------------------------
        % Fill K-space 2D real-time
        % ---------------------------------------------------------------------------------
        function objKspace = fillKspace2Drealtime(objKspace, objNav, objData, app)
            
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
            objKspace.kSpace = cell(objData.nr_coils);
            objKspace.kSpaceAvg = [];
            binTimesCard = objNav.heartTrigPoints;
            dimx = objData.dimx;
            dimy = objData.dimy;
            dimz = objData.dimz;
            nrKsteps = objData.nrKsteps;
            traj = objKspace.trajectory;
            nrCardFrames = app.nrCardFrames;
            nrRespFrames = app.nrRespFrames;
           
            for coilnr = 1:objData.nr_coils
                
                rawData = objKspace.raw{coilnr};
                
                nrDynamics = round(length(binTimesCard)/nrCardFrames);                                     % number of dynamics equals the number of heartbeats in the acquisition
                sortedKspace = complex(zeros(nrRespFrames,nrCardFrames,dimz,dimy,dimx,nrDynamics));       % fill temp k-space with zeros
                sortedAverages = zeros(nrRespFrames,nrCardFrames,dimz,dimy,dimx,nrDynamics);              % fill temp averages array with zeros
                nrReps = size(rawData,1);                                                                          % number of k-space repetitions
                nrKlines = nrReps * nrKsteps * dimz;                                                         % total number of k-lines
                
                % Unsorted k-lines
                rawData = permute(rawData,[3,1,2,4]);
                unsortedKlines = reshape(rawData,[1,1,1,nrKlines,dimx]);
                
                % Fill k-space
                fcnt = 1;
                
                for i = 1 : nrDynamics
                    
                    for j = 1 : nrCardFrames
                        
                        for k = 1:nrKlines
                            
                            if fcnt < length(binTimesCard)
                                
                                if k > binTimesCard(1,fcnt) && k < binTimesCard(1,fcnt+1)
                                    
                                    % The phase-encoding step using the trajectory info
                                    kline = traj(mod(k - 1,nrKsteps) + 1);
                                    
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
                    Ry = round(dimy/share/2);
                    Rx = round(dimx/share/2);
                    [Y,X] = ndgrid(1:dimy,1:dimx);
                    L = zeros(share,dimy,dimx);
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
                            weights(i,j) = retroKspace.gauss(i+j-1,share,0);
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
                                        
                                        ROI = reshape(squeeze(C(j,:,:)),[1 1 1 dimy dimx 1])*weights(j,abs(i));
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
                tukeyFilter(1,1,1,:,:,1) = retroKspace.circtukey2D(dimy,dimx,row,col,filterWidth);
                tmpKspace = tmpKspace.*tukeyFilter;

                % Report back
                objKspace.kSpace{coilnr} = tmpKspace;

            end

            objKspace.kSpaceAvg = tmpAverages;
            app.nrDynamics = nrDynamics;

        end % fillKspace2Drealtime
        
        
        
        % ---------------------------------------------------------------------------------
        % Fill K-space 3D P2ROUD
        % ---------------------------------------------------------------------------------
        function objKspace = fillKspace3Dp2roud(objKspace, objData, app)

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

            objKspace.kSpace = cell(objData.nr_coils);
            objKspace.kSpaceAvg = [];
            includeWindow = objData.includeWindow.*objData.excludeWindow;
            cardBinAss = objKspace.cardBinNrs;
            respBinAss = objKspace.respBinNrs;
            dimx = objData.dimx;
            dimy = objData.dimy;
            dimz = objData.dimz;
            nrKsteps = objData.nrKsteps;
            traj = objKspace.trajectory;
            nrCardFrames = app.nrCardFrames;
            nrRespFrames = app.nrRespFrames;
            nrDynamics = app.nrDynamics;
            
            for coilnr = 1:objData.nr_coils
                
                rawData = objKspace.raw{coilnr};
                
                sortedKspace = complex(zeros(nrRespFrames,nrCardFrames,dimz,dimy,dimx,nrDynamics));       % fill temp k-space with zeros
                sortedAverages = zeros(nrRespFrames,nrCardFrames,dimz,dimy,dimx,nrDynamics);              % fill temp nr averages array with zeros
                nrReps = size(rawData,1);                                                                      % number of k-space repetitions
                unsortedKspace = reshape(rawData,[1,size(rawData),1]);
                
                % Dynamics assignment
                totalk = nrReps * nrKsteps * dimz;
                dynBinAss = round(linspace(0.5, nrDynamics+0.49, totalk));       % List of increasing integer number 1 .. nr_dynamics evenly spaced over the entire acquistion time
                
                % Trajectory for 3D p2roud acquisition
                offset1 = floor(dimy/2 + 1);
                offset2 = floor(dimz/2 + 1);
                cnt1 = 1;
                for i = 1:nrReps
                    cnt2 = 1;
                    for j = 1:nrKsteps
                        for k = 1:dimz
                            traj3Dy(cnt1) = traj(cnt2)   + offset1; %#ok<*AGROW> 
                            traj3Dz(cnt1) = traj(cnt2+1) + offset2;
                            cnt1 = cnt1 + 1;
                            cnt2 = cnt2 + 2;
                        end
                    end
                end
                
                % Do the filling of k-space
                cnt = 0;
                
                for i = 1:nrReps                  % Loop through all repetitions
                    
                    for j = 1:nrKsteps            % Loop through all the phase-encoding steps
                        
                        for k = 1:dimz            % Loop through phase-encoding 3rd dimension
                            
                            cnt = cnt + 1;
                            
                            if (cardBinAss(cnt) > 0) && (includeWindow(cnt) == 1)     % If assigment == 0, this acquisition is discarded
                                
                                kLineY = traj3Dy(cnt);            % The phase-encoding step using the 3D trajectory info
                                kLineZ = traj3Dz(cnt);            % The 2nd phase-encoding
                                
                                sortedKspace(respBinAss(cnt),cardBinAss(cnt),kLineZ,kLineY,:,dynBinAss(cnt)) = sortedKspace(respBinAss(cnt),cardBinAss(cnt),kLineZ,kLineY,:,dynBinAss(cnt)) + unsortedKspace(1,i,k,j,:,1);     % add the data to the correct k-position
                                sortedAverages(respBinAss(cnt),cardBinAss(cnt),kLineZ,kLineY,:,dynBinAss(cnt)) = sortedAverages(respBinAss(cnt),cardBinAss(cnt),kLineZ,kLineY,:,dynBinAss(cnt)) + 1;                            % increase the number of averages with 1
                                
                            end
                            
                        end
                        
                    end
                    
                end
                
                % Normalize by number of averages
                sortedKspace = sortedKspace./sortedAverages;       % Normalize by number of averages
                sortedKspace(isnan(sortedKspace)) = complex(0);    % Correct for NaN because of division by zero in case of missing k-lines
                sortedKspace(isinf(sortedKspace)) = complex(0);
                
                % Apply a circular Tukey filter
                filterWidth = 0.1;
                tukeyFilter(1,1,:,:,:,1) = retroKspace.circtukey3D(dimz,dimy,dimx,lev,row,col,filterWidth);
                sortedKspace = sortedKspace.*tukeyFilter;
                
                % Report back
                objKspace.kSpace{coilnr} = sortedKspace;
                
            end

            % Report back
            objKspace.kSpaceAvg = sortedAverages;
            
        end % fillKspace3Dp2roud
        

        
        
        % ---------------------------------------------------------------------------------
        % Fill k-space 2D RADIAL
        % ---------------------------------------------------------------------------------
        function objKspace = fillKspace2Dradial(objKspace, objData, app)
            
            objKspace.kSpace = cell(objData.nr_coils);
            objKspace.kSpaceTraj = [];
            objKspace.kSpaceAvg = [];
         
            includeWindow = objData.includeWindow.*objData.excludeWindow;
            cardBinAss = objKspace.cardBinNrs;
            respBinAss = objKspace.respBinNrs;

            dimx = objData.dimx;
            dimz = objData.dimz;
 
            nrCardFrames = app.nrCardFrames;
            nrRespFrames = app.nrRespFrames;
            nrDynamics = app.nrDynamics;
            interpFactor = 16;

            % Create spokes from trajectory list of angles
            spoke = zeros(1,1,objData.nrKlines,1,dimx,1,3);
            for cnt = 1:objData.nrKlines
                spoke(1,1,cnt,1,:,1,1) = (-floor(dimx/2)+0.5:floor(dimx/2)-0.5)*cos(objKspace.trajectory(cnt)*pi/180);
                spoke(1,1,cnt,1,:,1,2) = (-floor(dimx/2)+0.5:floor(dimx/2)-0.5)*sin(objKspace.trajectory(cnt)*pi/180);
            end
       
            for coilnr = 1:objData.nr_coils
                
                rawData = permute(objKspace.raw{coilnr},[3 2 1 4]); % needed to be able to rearrange according to spokes, readout
                sortedKspace = complex(zeros(nrRespFrames,nrCardFrames,objData.nrKlines,dimz,dimx,nrDynamics));   % fill temp k-space with zeros
                sortedAverages = zeros(nrRespFrames,nrCardFrames,objData.nrKlines,dimz,dimx,nrDynamics);
                sortedTraj = zeros(nrRespFrames,nrCardFrames,objData.nrKlines,dimz,dimx,nrDynamics,3);
                unsortedKspace = reshape(rawData,[1,1,objData.nrKlines,1,dimx]);
      
                % Dynamics and slice assignment
                if dimz > 1
                    dynBinAss = ones(objData.nrKlines);
                    sliceAss = round(linspace(0.5, dimz+0.49, objData.nrKlines));

                else
                    dynBinAss = round(linspace(0.5, nrDynamics+0.49, objData.nrKlines));       % list of increasing integer number 1 .. nr_dynamics evenly spaced over the entire acquistion time
                    sliceAss = ones(objData.nrKlines);
                end

                % Loop over all acquired 2D spokes
                for cnt = 1:objData.nrKlines                

                    if (cardBinAss(cnt) > 0) && (includeWindow(cnt) == 1)

                        tmpKline1 = squeeze(unsortedKspace(1,1,cnt,1,:));

                        % Center echo
                        if app.CenterEchoCheckBox.Value
                            tmpKline2 = interp(tmpKline1,interpFactor);
                            [~,kCenter] = max(abs(tmpKline2));
                            kShift = floor(dimx/2)-kCenter/interpFactor;
                            tmpKline1 = retroKspace.fracCircShift(tmpKline1,kShift);
                        end

                        % Phase correction for k-space center
                        if app.PhaseCorrectCheckBox.Value
                            kCenterPhase = angle(tmpKline1(floor(dimx/2)+1));
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
                tukeyWindow(1,1,1,1,:,1) = tukeywin(dimx,filterWidth);
                sortedKspace = sortedKspace.*tukeyWindow;

                % Report back k-space
                objKspace.kSpace{coilnr} = sortedKspace;

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
            app.PlotTrajectoryMovieFrameFcn(xc,yc,[],dimx,true);

            % Report back averages and trajectory
            objKspace.kSpaceTraj = sortedTraj;
            objKspace.kSpaceAvg = sortedAverages;

        end % fillKspace2Dradial


        
        
        % ---------------------------------------------------------------------------------
        % Fill k-space 3D UTE
        % ---------------------------------------------------------------------------------
        function objKspace = fillKspace3Dute(objKspace, objData, app)

            % This function creates 2 arrays
            % (1) the kspace data sorted into the correct cardiac frames and phase-encoding positions
            % (2) the corresponding trajectories

            objKspace.kSpace = cell(objData.nr_coils);
            objKspace.kSpaceTraj = [];
            objKspace.kSpaceAvg = [];
            includeWindow = objData.includeWindow.*objData.excludeWindow;
            cardBinAss = objKspace.cardBinNrs;
            respBinAss = objKspace.respBinNrs;
            offset = app.DataOffsetRadialEditField.Value;
            dimx = length(objKspace.gradTrajectory);
          
            for coilnr = 1:objData.nr_coils

                % Unsorted data
                rawData = objKspace.raw{coilnr};

                % Dynamics assignment
                dynBinAss = round(linspace(0.5, app.nrDynamics+0.49, objData.nrKlines));       % list of increasing integer number 1 .. nr_dynamics evenly spaced over the entire acquistion time

                % K-space radial spokes
                for cnt = 1:dimx
                    spoke(1,1,:,1,cnt,1,1) = dimx*(objKspace.trajectory(1,:)/32767)*objKspace.gradTrajectory(cnt);
                    spoke(1,1,:,1,cnt,1,2) = dimx*(objKspace.trajectory(2,:)/32767)*objKspace.gradTrajectory(cnt);
                    spoke(1,1,:,1,cnt,1,3) = dimx*(objKspace.trajectory(3,:)/32767)*objKspace.gradTrajectory(cnt);
                end

                % Sort
                sortedKspace = complex(zeros(app.nrRespFrames,app.nrCardFrames,objData.nrKlines,1,dimx,app.nrDynamics));   % fill temp k-space with zeros
                sortedAverages = zeros(app.nrRespFrames,app.nrCardFrames,objData.nrKlines,1,dimx,app.nrDynamics);
                sortedTraj = zeros(app.nrRespFrames,app.nrCardFrames,objData.nrKlines,1,dimx,app.nrDynamics,3);            % fill temp trajectory with zeros
                unsortedKspace = reshape(rawData,[1,size(rawData),1]);

                % Check if offset does not lead to index larger than available
                offset(offset < 0) = 0;
                offset(offset > (size(unsortedKspace,5)-dimx)) = size(unsortedKspace,5)-dimx;
                app.DataOffsetRadialEditField.Value = offset;

                % Loop over acquired 3D spokes
                for cnt = 1:objData.nrKlines                

                    if (cardBinAss(cnt) > 0) && (includeWindow(cnt) == 1)

                        tmpKline = squeeze(unsortedKspace(1,cnt,1,1,1+offset:dimx+offset,1));

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
                tmpFilter = tukeywin(dimx*2,filterWidth);
                tmpFilter = tmpFilter(dimx+1:dimx*2);
                tukeyWindow(1,1,1,1,:,1) = tmpFilter;
                sortedKspace = sortedKspace.*tukeyWindow;

                % Report back k-space per coil
                objKspace.kSpace{coilnr} = sortedKspace;
           
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
            app.PlotTrajectoryMovieFrameFcn(xc,yc,zc,dimx,true);

            % Report back k-space trajectory and averages
            objKspace.kSpaceTraj = sortedTraj;
            objKspace.kSpaceAvg = sortedAverages;


        end % fillKspace3Dute


        
        % ---------------------------------------------------------------------------------
        % Reshape k-space
        % ---------------------------------------------------------------------------------
        function [objKspace,objReco] = reshapeKspace(objKspace,objReco,objData)
        
            switch objData.dataType

                case {'2D','2Dms','3D','3Dp2roud'}

                    % Reshape to either cardiac or respiratory CINE

                    [s1,s2,s3,s4,s5,s6] = size(objKspace.kSpaceAvg);
                    s = max([s1,s2]);
                    nrCoils = length(objKspace.kSpace);
                    for i = 1:nrCoils
                        objKspace.kSpace{i} = reshape(objKspace.kSpace{i},[s,s3,s4,s5,s6]);
                        objKspace.kSpace{i} = permute(objKspace.kSpace{i},[1,4,3,2,5,6]);
                    end
                    objKspace.kSpaceAvg = reshape(objKspace.kSpaceAvg,[s,s3,s4,s5,s6]);
                    objKspace.kSpaceAvg = permute(objKspace.kSpaceAvg,[1,4,3,2,5,6]);

                    % kSpace = frames, X, Y, Z(/slices), dynamics

                case {'2Dradial','3Dute'}

                    % Reshape to either cardiac or respiratory CINE
        
                    [s1,s2,s3,s4,s5,s6] = size(objKspace.kSpace{1});

                    s = max([s1,s2]);
                    nrCoils = length(objKspace.kSpace);
                    for i = 1:nrCoils
                        objKspace.kSpace{i} = reshape(objKspace.kSpace{i},[s,s3,s4,s5,s6]);
                        objKspace.kSpace{i} = permute(objKspace.kSpace{i},[1,4,2,3,5,6]);
                    end

                    objKspace.kSpaceTraj = reshape(objKspace.kSpaceTraj,[s,s3,s4,s5,s6,3]);
                    objKspace.kSpaceTraj = permute(objKspace.kSpaceTraj,[1,4,2,3,5,6]);
                    objKspace.kSpaceAvg = reshape(objKspace.kSpaceAvg,[s,s3,s4,s5,s6]);
                    objKspace.kSpaceAvg = permute(objKspace.kSpaceAvg,[1,4,2,3,5]);

                    % kSpace     = frames, X(readout), Y(spokes), Z(1 or slices), dynamics
                    % kSpaceTraj = frames, X(readout), Y(spokes), Z(1 or slices), dynamics, 3(trajectory coordinates)

                    % Set multi-slice and multi-dynamic flags for
                    % visualization of k-space trajectory
                    if size(objKspace.kSpace{1},4) > 1 
                        objReco.multiSliceFlag = true;
                    else
                        objReco.multiSliceFlag = false;
                    end

                    if size(objKspace.kSpace{1},5) > 1 
                        objReco.multiDynamicFlag = true;
                    else
                        objReco.multiDynamicFlag = false;
                    end                   

            end

        end % reshapeKspace
        

        
        % ---------------------------------------------------------------------------------
        % K-space statistics
        % ---------------------------------------------------------------------------------
        function [obj, app] = kSpaceStats(obj, objData, app)
            
            switch app.retroDataPars.dataType

                case {'2D','2Dms','3D','3Dp2roud'}

                    % Calculate effective number of averages
                    app.AveragesViewField.Value = mean2(obj.kSpaceAvg);
                    objData.NO_AVERAGES = round(app.AveragesViewField.Value);

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
                    objData.NO_AVERAGES = round(app.AveragesViewField.Value);

                    % Filling
                    filling = round(100*(nrNonZeros/nrFullNyquist)/(app.nrCardFrames*app.nrRespFrames*app.nrDynamics));
                    filling(filling>100) = 100;
                    app.FillingViewField.Value = filling;

                case '2Dradial'

                    % Kspace radial spokes
                    nrKpoints = objData.dimx;

                    % Kspace
                    sortedKspace = obj.kSpace{1};

                    % Number of averages
                    nrFullNyquist = (pi/2)*nrKpoints/2;
                    nrNonZeros = nnz(sortedKspace(:))/nrKpoints;
                    app.AveragesViewField.Value = (nrNonZeros/nrFullNyquist)/(app.nrCardFrames*app.nrRespFrames*app.nrDynamics);
                    objData.NO_AVERAGES = round(app.AveragesViewField.Value);

                    % Filling
                    filling = round(100*(nrNonZeros/nrFullNyquist)/(app.nrCardFrames*app.nrRespFrames*app.nrDynamics));
                    filling(filling>100) = 100;
                    app.FillingViewField.Value = filling;

            end

        end % kSpaceStats


        
        % ---------------------------------------------------------------------------------
        % K-space trajectory
        % ---------------------------------------------------------------------------------
        function [objKspace, objData] = getKspaceTrajectory(objKspace, objData, app)
            
            % Trajectory (linear, zigzag, user-defined, radial)
            switch objData.pe1_order
                
                case 0
                    objKspace.trajectory = linearTrajectory(objData.nrKsteps);
                    
                case 2
                    objKspace.trajectory = zigzagTrajectory(objData.nrKsteps);
                    
                case 3
                    objKspace.trajectory = arrayTrajectory(objData.NO_VIEWS,objData.gp_var_mul);
                    
                case 4 % reserved for 3D P2ROUD acquisition
                    flist = dir(fullfile(app.mrdImportPath,'*.txt'));
                    if ~isempty(flist)
                        objKspace.trajectory = load([flist(1).folder,filesep,flist(1).name]);
                        objData.dataType = '3Dp2roud';
                        app.TextMessage(strcat('P2ROUD trajectory: ',{' '},flist(1).name));
                    else
                        objData.dataType = '3D';
                        objKspace.trajectory = linearTrajectory(objData.nrKsteps);
                        app.TextMessage('WARNING: P2ROUD trajectory file not found, assuming linear k-space filling ...');
                        app.SetStatus(1);
                    end

                case 8 % Reserved for 2D radial
                    if app.HalfCircleButton.Value == 1              trajType = 1; end %#ok<SEPEX> 
                    if app.FullCircleButton.Value == 1              trajType = 2; end %#ok<SEPEX> 
                    if app.FullCircleInterleavedButton.Value == 1   trajType = 3; end %#ok<SEPEX> 
                    objKspace.trajectory = radialTrajectory(objData.nrKsteps,objData.nr_repetitions,trajType);

                case 9 % 3D UTE
                    flist = dir(fullfile(app.mrdImportPath,'lut*.txt'));
                    try
                        objKspace.trajectory = load([flist(1).folder,filesep,flist(1).name]);
                        objKspace.gradTrajectory = load('ktraj.txt');
                        objData.dataType = '3Dute';
                        app.TextMessage(strcat('3D UTE trajectory: ',{' '},flist(1).name));
                        objKspace.trajectory = reshape(objKspace.trajectory,[3,length(objKspace.trajectory)/3]);
                        objKspace.trajectory = objKspace.trajectory(:,1:objData.nr_repetitions);
                    catch ME
                        app.TextMessage(ME.message);
                        objKspace.trajectory = linearTrajectory(objData.nrKsteps);
                        app.TextMessage('ERROR: 3D UTE trajectory not found or invalid ...');
                        app.SetStatus(2);
                        app.objData.validDataFlag = false;
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
            function traject = radialTrajectory(nrPE,nrReps,trajType)

                % Radial k-space trajectory
                % The actual trajectory is later incorporated in the trajectory that is used by the Bart toolbox
                % Therefore the trajectory is simply linear here

                switch trajType

                    case 1

                        % Half circle linear
                        fullAngle = 180;

                        cnt = 1;
                        for i = 1:nrReps
                            for j = 1:nrPE
                                traject(cnt) = (j-1)*fullAngle/nrPE;
                                cnt = cnt + 1;
                            end
                        end

                    case 2

                        % Full circle linear
                        fullAngle = 360;

                        cnt = 1;
                        for i = 1:nrReps
                            for j = 1:nrPE
                                traject(cnt) = (j-1)*fullAngle/nrPE;
                                cnt = cnt + 1;
                            end
                        end

                    case 3
                        
                        % One of the trajectorie tried with Chayenne
                        cnt = 1;
                        for i = 1:nrReps
                            for j = 1:nrPE
                                if j < nrPE/2
                                     traject(cnt) = (j-1)*360/nrPE;
                                else
                                    traject(cnt) = (j-1)*360/nrPE + 180/nrPE;
                                end
                                cnt = cnt + 1;
                            end
                        end
                end

            end % radialTrajectory

        end % getKspaceTrajectory
        
        
        
        % ---------------------------------------------------------------------------------
        % Back to K-space
        % ---------------------------------------------------------------------------------
        function objKspace = backToKspace(objKspace, objData, objReco)

            objKspace.kSpaceMrd = [];
            [nf, ~, ~, dimz, nr] = size(objReco.movieExp);

            switch objData.dataType

                case {'2D','2Dms','2Dradial'}

                    for i = 1:nf
                        for j = 1:nr
                            for k = 1:dimz
                                objKspace.kSpaceMrd(i,:,:,k,j) = retroKspace.fft2Dmri(squeeze(objReco.movieExp(i,:,:,k,j)));
                            end
                        end
                    end

                    % Samples, views, views2, slices, echoes (frames), experiments
                    objKspace.kSpaceMrd = permute(objKspace.kSpaceMrd,[2,3,6,4,1,5]);
    
                case {'3D','3Dp2roud','3Dute'}

                    for i = 1:nf
                        for j = 1:nr
                            objKspace.kSpaceMrd(i,:,:,:,j) = retroKspace.fft3Dmri(squeeze(objReco.movieExp(i,:,:,:,j)));
                        end
                    end

                    % Samples, views, views2, slices, echoes (frames), experiments
                    objKspace.kSpaceMrd = permute(objKspace.kSpaceMrd,[2,3,4,6,1,5]);

            end

        end


    end % Public methods
    


    
    
    % ---------------------------------------------------------------------------------
    % Static methods
    % ---------------------------------------------------------------------------------
    methods (Static)
        
        
        % ---------------------------------------------------------------------------------
        % 2D Tukey filter
        % ---------------------------------------------------------------------------------
        function output = circtukey2D(dimy, dimx, row, col, filterWidth)
            
            domain = 256;
            base = zeros(domain,domain);
            
            tukey1 = tukeywin(domain,filterWidth);
            tukey1 = tukey1(domain/2+1:domain);
            
            shifty = (row-dimy/2)*domain/dimy;
            shiftx = (col-dimx/2)*domain/dimx;
            
            y = linspace(-domain/2, domain/2, domain);
            x = linspace(-domain/2, domain/2, domain);
            
            parfor i=1:domain
                for j=1:domain
                    rad = round(sqrt((shiftx-x(i))^2 + (shifty-y(j))^2));
                    if (rad <= domain/2) && (rad > 0)
                        base(j,i) = tukey1(rad);
                    end
                end
            end
            
            output = imresize(base,[dimy dimx]);
            
        end % circtukey2D
        
        
        
        % ---------------------------------------------------------------------------------
        % 3D Tukey filter
        % ---------------------------------------------------------------------------------
        function output = circtukey3D(dimz,dimy,dimx,lev,row,col,filterWidth)
            
            domain = 256;
            
            base = zeros(domain,domain,domain);
            
            tukey1 = tukeywin(domain,filterWidth);
            tukey1 = tukey1(domain/2+1:domain);
            
            shiftz = (lev-dimz/2)*domain/dimz;
            shifty = (row-dimy/2)*domain/dimy;
            shiftx = (col-dimx/2)*domain/dimx;
            
            z = linspace(-domain/2, domain/2, domain);
            y = linspace(-domain/2, domain/2, domain);
            x = linspace(-domain/2, domain/2, domain);
            
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
            
            output = imresize3(base,[dimz dimy dimx]);
            
        end % circtukey3D
        
        
        
        % ---------------------------------------------------------------------------------
        % Centric K-space filling scheme
        % ---------------------------------------------------------------------------------
        function scheme = centricFilling(noViews2)
            
            ord2 = zeros(2*round(noViews2/2));
            
            for g= 1:round(noViews2/2)
                
                ord2(2*g-1) = noViews2/2+g;
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
        % 2D FFT
        % ---------------------------------------------------------------------------------
        function Y = fft2Dmri(x)

            Y = fftshift(ifft(fftshift(x,1),[],1),1)*sqrt(size(x,1));
            Y = fftshift(ifft(fftshift(Y,2),[],2),2)*sqrt(size(x,2));
                        
        end % fft2Dmri



        % ---------------------------------------------------------------------------------
        % 3D FFT
        % ---------------------------------------------------------------------------------
        function Y = fft3Dmri(x)

            Y = fftshift(ifft(fftshift(x,1),[],1),1)*sqrt(size(x,1));
            Y = fftshift(ifft(fftshift(Y,2),[],2),2)*sqrt(size(x,2));
            Y = fftshift(ifft(fftshift(Y,3),[],3),3)*sqrt(size(x,3));

        end % fft3Dmri



        % ---------------------------------------------------------------------------------
        % Fractional circshift
        % ---------------------------------------------------------------------------------
        function output = fracCircShift(input,shiftsize)

            int = floor(shiftsize);     % Integer portions of shiftsize
            fra = shiftsize - int;      % Fractional portions of shiftsize
            dim = numel(shiftsize);
            output = input;
            for n = 1:numel(shiftsize)  % The dimensions are treated one after another.
                intn = int(n);
                fran = fra(n);
                shift1 = zeros(dim,1);
                shift1(n) = intn;
                shift2 = zeros(dim,1);
                shift2(n) = intn+1;
                % Linear intepolation:
                output = (1-fran)*circshift(output,shift1) + fran*circshift(output,shift2);
            end

        end % fracCircShift




    end % Static methods
    

end % retroKspace

