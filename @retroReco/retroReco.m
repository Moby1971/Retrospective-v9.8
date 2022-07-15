classdef retroReco
    
    % Reconstruction and movie class for retrospective app
    
    properties
        
        movieExp                        % movie for movie export
        movieApp                        % movie for viewing in the app
        senseMap                        % sense map
        sos = true                      % Sum of squares reco true or false
        rescaleSlope                    % Dicom info RescaleSlope for image scaling
        rescaleIntercept                % Dicom info RescaleIntercept for image scaling
        multiSliceFlag = false          % multi-slice true or false
        multiDynamicFlag = false        % mutli-dynamic true or false
        totalVariation = text           % total variation (T) or total generalized variation (G)
        
    end
    
    
    
    
    
    methods
        
        % ---------------------------------------------------------------------------------
        % Object constructor
        % ---------------------------------------------------------------------------------
        function obj = retroReco
            
        end % retroReco
        
        
        
        % ---------------------------------------------------------------------------------
        % 2D reconstruction
        % ---------------------------------------------------------------------------------
        function objReco = reco2D(objReco, objData, objKspace, app)
            
            if app.bartDetectedFlag
                
                % Perform CS reconstruction with the Bart toolbox, preferred option
                app.TextMessage('Reconstructing the data with the BART toolbox ...');
                app.ProgressGauge.Value = 0;
                drawnow;
                cs_reco2D_mc;
                app.ProgressGauge.Value = 100;
                app.validRecoFlag = true;
                app.validSenseMapFlag = true;
                drawnow;
                
            else
                
                % Perform CS reconstruction in matlab, slow but works
                app.TextMessage('WARNING: Reconstructing the data with MATLAB ...');
                app.TextMessage('Slow reconstruction, please be patient ...');
                app.ProgressGauge.Value = 0;
                app.SetStatus(1);
                drawnow;
                cs_reco2D_mat_mc;
                app.ProgressGauge.Value = 100;
                app.validRecoFlag = true;
                app.validSenseMapFlag = true;
                drawnow;
                
            end
            
            
            % ---------------------------------------------------------------------------------
            % 2D compressed sensing reconstruction with Matlab
            % ---------------------------------------------------------------------------------
            function cs_reco2D_mat_mc
                
                kSpaceIn = objKspace.kSpace;
                averages = objKspace.kSpaceAvg;
                lambdaTV = app.TVcineEditField.Value;
                
                % kspace = {coil}[frames, x, y, slice, dynamic]
                %                   1     2  3    4       5
                dimf = size(kSpaceIn{1},1);
                dimz = size(kSpaceIn{1},4);
                dimd = size(kSpaceIn{1},5);
                nrCoils = objData.nr_coils;
                
                % in case of 1 frame, duplicate that frame to facilitate reconstruction
                if dimf == 1
                    for i = 1:nrCoils
                        kSpaceIn{i}(2,:,:,:,:) = kSpaceIn{i}(1,:,:,:,:);
                    end
                end
                
                % kspace data x,y,frames,slices,dynamics
                for i = 1:nrCoils
                    kSpaceIn{i} = permute(kSpaceIn{i},[2,3,1,4,5]);
                end
                
                % kspace data x,y,frames,slices,dynamics,coils
                kSpace = zeros([size(kSpaceIn{1}),nrCoils]);
                for i = 1:nrCoils
                    kSpace(:,:,:,:,:,i) = kSpaceIn{i};
                end
                
                % averages data x,y,frames,slices,dynamics
                averages = permute(averages,[2,3,1,4,5]);
                
                % reset progress counter
                param.iteration = 0;
                
                % pad to next power of 2
                dimx = 2^nextpow2(size(kSpace,1));
                dimy = 2^nextpow2(size(kSpace,2));
                
                % for convenience make rectangular matrix size
                mdimxy = max([dimx,dimy]);
                dimx = mdimxy;
                dimy = mdimxy;
                dimf = size(kSpace,3);
                
                % pre-allocate memory for image_out
                imageOut = zeros(dimf,dimx,dimy,dimz,dimd);
                
                % slice and dynamic loop
                for slice = 1:dimz
                    
                    for dynamic = 1:dimd
                        
                        % kspace of slice and dynamic
                        kData = squeeze(kSpace(:,:,:,slice,dynamic,:));
                        
                        % zero padding
                        padSizex = round((dimx - size(kData,1))/2);
                        padSizey = round((dimy - size(kData,2))/2);
                        kData = padarray(kData,[padSizex,padSizey,0],'both');
                        kData = kData(1:dimx,1:dimy,:,:);
                        
                        % size of the data
                        [nx,ny,~,nrCoils]=size(kData);
                        
                        % normalize the data in the range of approx 0 - 1 for better numerical stability
                        kData = kData/max(abs(kData(:)));
                        
                        % kspace mask, 0 = nodata, 1 = data, zero-pad to same size as k-space
                        mask = squeeze(averages(:,:,:,slice,dynamic));
                        mask = padarray(mask,[padSizex,padSizey,0],'both');
                        mask = mask(1:dimx,1:dimy,:);
                        mask = mask./mask;
                        mask(isnan(mask)) = 1;
                        mask = logical(mask);
                        
                        % coil sensitivity map
                        b1 = ones(nx,ny,nrCoils);
                        
                        % data
                        param.y = kData;
                        
                        % reconstruction design matrix
                        param.E = Emat_yxt(mask,b1);
                        
                        % Total variation (TV) constraint in the temporal domain
                        % for 'consistency' with Bart reconstruction, TV seems to be scale empirically by a factor of 8
                        % TV only in the time domain
                        param.TV = TVOP;
                        param.TVWeight = lambdaTV/8;
                        
                        % number of iterations, 2 x 10 iterations
                        param.nite = 10;
                        param.nouter = 2;
                        param.totaliterations = dimz * dimd * param.nouter * param.nite;
                        
                        % linear reconstruction
                        kData1 = randn(size(kData))/2000 + kData;  % add a little bit of randomness, such that linear reco is not exactly right
                        reconDft = param.E'*kData1;
                        
                        % iterative reconstruction
                        reconCS = reconDft;
                        for n = 1:param.nouter
                            [reconCS,param.iteration] = CSL1NlCg(app,reconCS,param);
                        end
                        
                        % rearrange to correct orientation: frames, x, y
                        imageTmp = flip(permute(squeeze(abs(reconCS)),[3, 1, 2]),2);

                        % output reconstructed image
                        imageOut(:,:,:,slice,dynamic) = imageTmp;
                        
                    end
                    
                end
                
                % correct back to 1 frame reconstruction
                if dimf == 1
                    imageOut = imageOut(1,:,:,:,:);
                end
           
                % shift image in phase-encoding direction if needed
                objReco.movieExp = circshift(imageOut,-objData.pixelshift1,3);
                objReco.senseMap = ones(size(objReco.movieExp));
                
            end
            
            
            % ---------------------------------------------------------------------------------
            % 2D compressed sensing reconstruction with Bart toolbox
            % ---------------------------------------------------------------------------------
            function cs_reco2D_mc
                
                kSpaceIn = objKspace.kSpace;
                Wavelet = app.WVxyzEditField.Value;
                TVxy = app.TVxyzEditField.Value;
                LR = app.LLRxyzEditField.Value;
                TVt = app.TVcineEditField.Value;
                TVd = app.TVdynEditField.Value;
                ESPIRiT = app.ESPIRiTCheckBox.Value;
                SOS = objReco.sos;
                
                % app = matlab app
                % kspace_in = sorted k-space
                % nc = number of RF receiver coils
                % averages_in = number of averages per k-line
                % Wavelet = wavelet L1-norm regularization factor
                % TVt = total variation in CINE dimension
                % TVxy = total variation in xy-dimension regularization
                % TVd = total variation in dynamic dimension
                % ESPIRiT = reconstruction of multi-coil data with ESPIRiT true or false
                
                % kspace_in = {coil}[CINE, x, y, slice, dynamic]
                %                     1    2  3    4       5
                dimf = size(kSpaceIn{1},1);
                dimx = 2^nextpow2(size(kSpaceIn{1},2));
                dimy = 2^nextpow2(size(kSpaceIn{1},3));
                dimz = size(kSpaceIn{1},4);
                dimd = size(kSpaceIn{1},5);
                nrCoils = objData.nr_coils;
                
                % for convenience make rectangular matrix size
                mdimxy = max([dimx,dimy]);
                dimx = mdimxy;
                dimy = mdimxy;
                
                % resize k-space to next power of 2
                for i = 1:nrCoils
                    kSpaceIn{i} = bart(app,['resize -c 1 ',num2str(dimx),' 2 ',num2str(dimy)],kSpaceIn{i});
                end
                
                % put all data in a normal matrix
                kSpace = zeros(dimf,dimx,dimy,dimz,dimd);
                for l = 1:dimd
                    for k = 1:dimz
                        for j = 1:dimf
                            for i = 1:nrCoils
                                kSpace(j,:,:,k,l,i) = kSpaceIn{i}(j,:,:,k,l);
                            end
                        end
                    end
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
                
                kspacePics = permute(kSpace,[7,3,2,6,8,9,10,11,12,13,1,5,14,4]);
                
                % wavelet in spatial dimensions 2^1+2^2=6
                % total variation in spatial dimensions 2^1+2^2=6
                % total variation in cine dimension 2^10 = 1024
                % total variation in dynamic dimension 2^11 = 2048

                if ESPIRiT && nrCoils>1
                    
                    % ESPIRiT reconstruction
                    TextMessage(app,'ESPIRiT reconstruction ...');
                    
                    % Calculate coil sensitivity maps with ecalib bart function
                    kspace_pics_sum = sum(kspacePics,[11,12]);
                    sensitivities = bart(app,'ecalib -S -I -a', kspace_pics_sum);      % ecalib with softsense
                    
                    picsCommand = 'pics ';
                    if Wavelet>0
                        picsCommand = [picsCommand,' -RW:6:0:',num2str(Wavelet)];
                    end
                    if TVxy>0
                        picsCommand = [picsCommand,' -R',objReco.totalVariation,':6:0:',num2str(TVxy)];
                    end
                    if LR>0
                        % Locally low-rank in the spatial domain
                        blocksize = round(dimx/16);  % Block size
                        blocksize(blocksize < 8) = 8;
                        picsCommand = [picsCommand,' -RL:6:6:',num2str(LR),' -b',num2str(blocksize)];
                    end
                    if TVt>0
                        picsCommand = [picsCommand,' -R',objReco.totalVariation,':1024:0:',num2str(TVt)];
                    end
                    if TVd>0
                        picsCommand = [picsCommand,' -R',objReco.totalVariation,':2048:0:',num2str(TVd)];
                    end
                    imageReg = bart(app,picsCommand,kspacePics,sensitivities);
                    
                    % Sum of squares reconstruction
                    if SOS
                        imageReg = bart(app,'rss 16', imageReg);
                    end
                    imageReg = abs(imageReg);
                    
                else
                    
                    % Reconstruction without sensitivity correction, including coil scaling
                    sensitivities = ones(1,dimy,dimx,nrCoils,1,1,1,1,1,1,1,1,1,dimz);
                    
                    % regular reconstruction
                    picsCommand = 'pics ';
                    if Wavelet>0
                        picsCommand = [picsCommand,' -RW:6:0:',num2str(Wavelet)];
                    end
                    if TVxy>0
                        picsCommand = [picsCommand,' -R',objReco.totalVariation,':6:0:',num2str(TVxy)];
                    end
                    if LR>0
                        blocksize = round(dimx/16);  % Block size
                        blocksize(blocksize < 8) = 8;
                        picsCommand = [picsCommand,' -RL:6:6:',num2str(LR),' -b',num2str(blocksize)];
                    end
                    if TVt>0
                        picsCommand = [picsCommand,' -R',objReco.totalVariation,':1024:0:',num2str(TVt)];
                    end
                    if TVd>0
                        picsCommand = [picsCommand,' -R',objReco.totalVariation,':2048:0:',num2str(TVd)];
                    end
                    imageReg = bart(app,picsCommand,kspacePics,sensitivities);

                    % total-generalized-variation denoising. Seems not really useful
                    % imageReg2 = bart(app,'tgv 0.001 1030',imageReg);
                    % imageReg = imageReg2(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,1);
                    
                    % Take absolute values
                    imageReg = abs(imageReg);
                    
                end
                
                %                              y    x  frames dynamic slices
                imageReg = reshape(imageReg,[dimy,dimx,dimf,dimd,dimz]);

                % rearrange to correct orientation: frames, x, y, slices, dynamics
                imageOut = flip(permute(imageReg,[3,2,1,5,4]),3);

                % sense map orientations: x, y, slices, map1, map2
                senseMap1 = flip(permute(abs(sensitivities),[3,2,14,4,5,1,6,7,8,9,10,11,12,13]),2);
                
                % normalize sense map to reasonable value range
                senseMap1 = senseMap1*4095/max(senseMap1(:));

                % shift image in phase-encoding direction if needed
                objReco.movieExp = circshift(imageOut,-objData.pixelshift1,3);
                objReco.senseMap = circshift(senseMap1,-objData.pixelshift1,3);
                
            end
            
        end
        
                
        
        % ---------------------------------------------------------------------------------
        % 3D reconstruction
        % ---------------------------------------------------------------------------------
        function objReco = reco3D(objReco, objData, objKspace, app)
            
            if app.bartDetectedFlag
                
                % perform CS reconstruction with the Bart toolbox, preferred option
                app.TextMessage('Reconstructing the 3D data with the BART toolbox ...');
                app.ProgressGauge.Value = 0;
                drawnow;
                cs_reco3D_mc;
                app.ProgressGauge.Value = 100;
                app.validRecoFlag = true;
                app.validSenseMapFlag = true;
                drawnow;
                
            else
                
                % perform cs reconstruction in matlab, slow but works
                app.TextMessage('WARNING: Reconstructing the data with MATLAB ...');
                app.TextMessage('Slow reconstruction, please be patient ...');
                app.ProgressGauge.Value = 0;
                app.SetStatus(1);
                cs_reco3D_mat_mc;
                app.ProgressGauge.Value = 100;
                app.validRecoFlag = true;
                app.validSenseMapFlag = true;
                drawnow;
                
            end
            
            
            % ---------------------------------------------------------------------------------
            % 3D compressed sensing reconstruction with Matlab
            % ---------------------------------------------------------------------------------
            function cs_reco3D_mat_mc
                
                kSpaceIn = objKspace.kSpace;
                averages = objKspace.kSpaceAvg;
                lambdaTV = app.TVcineEditField.Value;
     
                % kdata_in = {coil}[frames, x, y, z, dynamics]
                %                      1    2  3  4     5
                dimf = size(kSpaceIn{1},1);
                dimd = size(kSpaceIn{1},5);
                nrCoils = objData.nr_coils;
                
                % in case of 1 frame, duplicate that frame to facilitate reconstruction
                if dimf == 1
                    for i = 1:nrCoils
                        kSpaceIn{i}(2,:,:,:,:) = kSpaceIn{i}(1,:,:,:,:);
                    end
                end
                
                % kspace data x,y,z,frames,dynamics
                for i = 1:nrCoils
                    kSpaceIn{i} = permute(kSpaceIn{i},[2,3,4,1,5]);
                end

                % kspace data x,y,z,frames,dynamics,coils
                kSpace = zeros([size(kSpaceIn{1}),nrCoils]);
                for i = 1:nrCoils
                    kSpace(:,:,:,:,:,i) = kSpaceIn{i};
                end
                
                % averages data x,y,z,frames,dynamics
                averages = permute(averages,[2,3,4,1,5]);
                
                % reset progress counter
                param.iteration = 0;
                
                % pad to next power of 2
                dimx = 2^nextpow2(size(kSpace,1));
                dimy = 2^nextpow2(size(kSpace,2));
                dimz = 2^nextpow2(size(kSpace,3));
                
                % for convenience make rectangular matrix size
                mdimxy = max([dimx,dimy]);
                dimx = mdimxy;
                dimy = mdimxy;
                dimf = size(kSpace,4);
                
                % pre-allocate memory for image_out
                imageOut = zeros(dimf,dimx,dimy,dimz,dimd);
                
                for dynamic = 1:dimd
                    
                    % kspace of slice and dynamic
                    kData = squeeze(kSpace(:,:,:,:,dynamic,:));
                    
                    % zero padding
                    padSizex = round((dimx - size(kData,1))/2);
                    padSizey = round((dimy - size(kData,2))/2);
                    padSizez = round((dimz - size(kData,3))/2);
                    kData = padarray(kData,[padSizex,padSizey,padSizez,0],'both');
                    kData = kData(1:dimx,1:dimy,1:dimz,:,:);
                    
                    % size of the data, z,y,x,frames,coils
                    [nx,ny,nz,~,nrCoils] = size(kData);
                    
                    % normalize the data in the range of approx 0 - 1 for better numerical stability
                    kData = kData/max(abs(kData(:)));
                    
                    % kspace mask, 0 = nodata, 1 = data, zero-pad to same size as k-space
                    mask = squeeze(averages(:,:,:,:,dynamic));
                    mask = padarray(mask,[padSizex,padSizey,padSizez,0],'both');
                    mask = mask(1:dimx,1:dimy,1:dimz,:);
                    mask = mask./mask;
                    mask(isnan(mask)) = 1;
                    mask = logical(mask);
                    
                    % coil sensitivity map
                    b1 = ones(nx,ny,nz,nrCoils);
                    
                    % data
                    param.y = kData;
                    
                    % reconstruction design matrix
                    param.E = Emat_zyxt(mask,b1);
                    
                    % Total variation (TV) constraint in the temporal domain
                    % for 'consistency' with Bart reconstruction, TV seems to be scale empirically by a factor of 8
                    % TV only in the time domain
                    param.TV = TVOP3D;
                    param.TVWeight = lambdaTV/8;
                    
                    % number of iterations, 2 x 10 iterations
                    param.nite = 10;
                    param.nouter = 2;
                    param.totaliterations = dimd * param.nouter * param.nite;
                    
                    % linear reconstruction
                    kData1 = randn(size(kData))/2000 + kData;  % add a little bit of randomness, such that linear reco is not exactly right
                    reconDft = param.E'*kData1;
                    
                    % iterative reconstruction
                    reconCs=reconDft;
                    for n=1:param.nouter
                        [reconCs,param.iteration] = CSL1NlCg(app,reconCs,param);
                    end
                    
                    % rearrange to correct orientation: frames, x, y, z
                    imageTmp = flip(permute(abs(reconCs),[4, 1, 2, 3]),2);

                    % output reconstructed image
                    imageOut(:,:,:,:,dynamic) = imageTmp;
                    
                end
                
                % correct back to 1 frame reconstruction
                if dimf == 1
                    imageOut = imageOut(1,:,:,:,:);
                end
                
                % shift image in phase-encoding directions if needed
                objReco.movieExp = circshift(imageOut,-objData.pixelshift1,3);
                objReco.movieExp = circshift(imageOut,-objData.pixelshift2,4);
                objReco.senseMap = ones(size(objReco.movieExp));
                
            end
            
            
            % ---------------------------------------------------------------------------------
            % 3D compressed sensing reconstruction with Bart toolbox
            % ---------------------------------------------------------------------------------
            function cs_reco3D_mc
                
                kSpaceIn = objKspace.kSpace;
                Wavelet = app.WVxyzEditField.Value;
                TVxyz = app.TVxyzEditField.Value;
                LR = app.LLRxyzEditField.Value;
                TVt = app.TVcineEditField.Value;
                TVd = app.TVdynEditField.Value;
                ESPIRiT = app.ESPIRiTCheckBox.Value;
                
                % app = matlab app
                % kspace_in = sorted k-space
                % nc = number of RF receiver coils
                % averages_in = number of averages per k-line
                % Wavelet = wavelet L1-norm regularization factor
                % TVt = total variation in time regularization
                % TVxyz = total variation in xyz-dimension regularization
                % VNSA = variable number of signal averages correction true or false
                % ESPIRiT = reconstruction of multi-coil data with ESPIRiT true or false
                
                % kspace_in = {coil}[frames, x, y, z, dynamics]
                %                     1      2  3  4     5
                dimf = size(kSpaceIn{1},1);
                dimx = 2^nextpow2(size(kSpaceIn{1},2));
                dimy = 2^nextpow2(size(kSpaceIn{1},3));
                dimz = 2^nextpow2(size(kSpaceIn{1},4));
                dimd = size(kSpaceIn{1},5);
                nrCoils = objData.nr_coils;
                
                % for convenience make rectangular matrix size
                mdimxy = max([dimx,dimy]);
                dimx = mdimxy;
                dimy = mdimxy;
                
                % resize to next power of 2
                for i = 1:nrCoils
                    kSpaceIn{i} = bart(app,['resize -c 1 ',num2str(dimx),' 2 ',num2str(dimy),' 3 ',num2str(dimz)],kSpaceIn{i});
                end
                
                % kspace suitable for bart
                kSpace = zeros(dimf,dimx,dimy,dimz,dimd,nrCoils);
                for l = 1:dimd
                    for k = 1:dimz
                        for j = 1:dimf
                            for i = 1:nrCoils
                                kSpace(j,:,:,k,l,i) = kSpaceIn{i}(j,:,:,k,l);
                            end
                        end
                    end
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
                
                kspace_pics = permute(kSpace,[4,3,2,6,7,8,9,10,11,12,1,5,13,14]);
                
                % wavelet in spatial dimensions 2^0+2^1+2^2=7
                % total variation in spatial dimensions 2^0+2^1+2^2=7
                % total variation in time 2^10 = 1024
                % total variation in dynamic dimension 2^11 = 2048
                
                if ESPIRiT && nrCoils>1
                    
                    TextMessage(app,'ESPIRiT reconstruction ...');
                    kspace_pics_sum = sum(kspace_pics,[11,12]);
                    sensitivities = bart(app,'ecalib -I -S -a', kspace_pics_sum);
                    
                    app.ProgressGauge.Value = 25;
                    drawnow;
                    
                    picsCommand = 'pics -S ';
                    if Wavelet>0
                        picsCommand = [picsCommand,' -RW:7:0:',num2str(Wavelet)];
                    end
                    if TVxyz>0
                        picsCommand = [picsCommand,' -R',objReco.totalVariation,':7:0:',num2str(TVxyz)];
                    end
                    if LR>0
                        % Locally low-rank in the spatial domain
                        blocksize = round(dimx/16);  % Block size
                        blocksize(blocksize < 8) = 8;
                        picsCommand = [picsCommand,' -RL:7:7:',num2str(LR),' -b',num2str(blocksize)];
                    end
                    if TVt>0
                        picsCommand = [picsCommand,' -R',objReco.totalVariation,':1024:0:',num2str(TVt)];
                    end
                    if TVd>0
                        picsCommand = [picsCommand,' -R',objReco.totalVariation,':2048:0:',num2str(TVd)];
                    end
                    imageReg = bart(app,picsCommand,kspace_pics,sensitivities);
                    
                    app.ProgressGauge.Value = 95;
                    drawnow;
                    imageReg = bart(app,'rss 16', imageReg);
                    imageReg = abs(imageReg);
                    
                else
                    
                    % Reconstruction without sensitivity correction
                    sensitivities = ones(dimz,dimy,dimx,nrCoils,1,1,1,1,1,1,1,1,1,1);
                    
                    picsCommand = 'pics -S ';
                    if Wavelet>0
                        picsCommand = [picsCommand,' -RW:7:0:',num2str(Wavelet)];
                    end
                    if TVxyz>0
                        picsCommand = [picsCommand,' -R',objReco.totalVariation,':7:0:',num2str(TVxyz)];
                    end
                    if LR>0
                        blocksize = round(dimx/16);  % Block size
                        blocksize(blocksize < 8) = 8;
                        picsCommand = [picsCommand,' -RL:7:7:',num2str(LR),' -b',num2str(blocksize)];
                    end
                    if TVt>0
                        picsCommand = [picsCommand,' -R',objReco.totalVariation,':1024:0:',num2str(TVt)];
                    end
                    if TVd>0
                        picsCommand = [picsCommand,' -R',objReco.totalVariation,':2048:0:',num2str(TVd)];
                    end
                    
                    imageReg = abs(bart(app,picsCommand,kspace_pics,sensitivities));
                    
                end
                
                % rearrange to correct orientation: frames, x, y, z, dynamics
                imageReg = reshape(imageReg,[dimz,dimy,dimx,dimf,dimd]);
                imageOut = flip(flip(permute(imageReg,[4,3,2,1,5]),3),4);
                imageOut = circshift(imageOut,1,4);
                
                % sense map orientations: x, y, z, map1, map2
                senseMap1 = flip(flip(permute(abs(sensitivities),[3,2,1,4,5,6,7,8,9,10,11,12,13,14]),2),3);
                senseMap1 = circshift(senseMap1,1,4);
                
                % normalize sense map to reasonable value range
                senseMap1 = senseMap1*4095/max(senseMap1(:));
                
                % shift image in phase-encoding directions if needed
                objReco.movieExp = circshift(imageOut,-objData.pixelshift1,3);
                objReco.movieExp = circshift(imageOut,-objData.pixelshift2,4);
                objReco.senseMap = circshift(senseMap1,-objData.pixelshift1,3);
                objReco.senseMap = circshift(senseMap1,-objData.pixelshift2,4);
                
            end
            
        end % reco3D
        
                
        
        % ---------------------------------------------------------------------------------
        % 2D radial reconstruction with the Bart toolbox or Matlab
        % ---------------------------------------------------------------------------------         
        function objReco = reco2Dradial(objReco, objData, objKspace, app)
      
            if app.bartDetectedFlag

                app.TextMessage('Reconstructing 2D radial data with the BART toolbox ...');
                app.ProgressGauge.Value = 0;
                drawnow;

                Wavelet = app.WVxyzEditField.Value;
                TVxyz = app.TVxyzEditField.Value;
                LR = app.LLRxyzEditField.Value;
                TVt = app.TVcineEditField.Value;
                TVd = app.TVdynEditField.Value;
                nrCoils = objData.nr_coils;
                ESPIRiT = app.ESPIRiTCheckBox.Value;

                % Wavelet = wavelet L1-norm regularization factor
                % TVxyz = total variation in xyz-dimension regularization
                % LR = low rank regularization
                % TVt = total variation in time regularization
                % TVd = total variation in dynamic dimension regularization
                % nCoils = number of RF receiver coils

                dimf = size(objKspace.kSpace{1},1);
                dimx = size(objKspace.kSpace{1},2);
                dimy = dimx;
                dimz = size(objKspace.kSpace{1},4);
                dimd = size(objKspace.kSpace{1},5);

                kSpace = zeros(size(objKspace.kSpace{1}));
                for i = 1:nrCoils
                    kSpace(:,:,:,:,:,i) = objKspace.kSpace{i};
                end
                traj = objKspace.kSpaceTraj;

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
                kSpacePics = permute(kSpace,[7, 2, 3, 6,14, 8, 9,10,11,12,1, 5, 13, 4]);

                % Rearrange for BART     1  2  3  4  5  6  7  8  9 10 11 12 13 14
                trajPics = permute(traj,[6, 2, 3,14, 7, 8, 9,10,11,12, 1, 5,13, 4]);

                % Calibration and density correction size
                kdim = round(dimx/2);
                if mod(kdim,2) == 1
                    kdim = kdim + 1;
                end
                kdim(kdim < 32) = 32;
                kdim(kdim > dimx) = dimx;
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

                    % Calculate k-space and trajectory sum of all frames
                    % and dynamics, remove all zero values
                    kSpacePicsSum = sum(kSpacePics(:,:,:,:,:,:,:,:,:,:,:,:,:,floor(dimz/2)+1),[11,12]);
                    trajPicsSum = sum(trajPics(:,:,:,:,:,:,:,:,:,:,:,:,:,floor(dimz/2)+1),[11,12]);

                    ze = squeeze(abs(kSpacePicsSum(1,floor(end/2),:))) > 0;
                    kSpacePicsSum = kSpacePicsSum(:,:,ze);
                    trajPicsSum = trajPicsSum(:,:,ze);

                    if app.RingMethodCheckBox.Value

                        % Sent gradient delay vector back to app
                        app.GxDelayEditField.Value = 0;
                        app.GyDelayEditField.Value = 0;
                        app.GzDelayEditField.Value = 0;

                        % Ring method using estdelay in Bart
                        try
                            delaysBart = bart(app,'estdelay -r4',trajPicsSum,kSpacePicsSum);
                        catch ME
                            app.TextMessage(ME.message);
                            app.TextMessage('Ring gradient delay estimation failed ...');
                            app.TextMessage('Trying iterative method ...');
                            app.SetStatus(1);
                            app.RingMethodCheckBox.Value = 0;
                            delaysBart = "0:0:0";
                        end

                        % Remove unknown warning
                        try
                            ff = strfind(delaysBart,"[0m");
                            delaysBart = delaysBart(ff:end);
                            delaysBart = erase(delaysBart,"[0m");
                            delaysBart = erase(delaysBart,newline);
                        catch
                        end

                        delaysBart = strrep(delaysBart,':',',');
                        dTotal = str2num(delaysBart); %#ok<ST2NM> 

                        % Sent gradient delay vector back to app
                        app.GxDelayEditField.Value = double(dTotal(1));
                        app.GyDelayEditField.Value = double(dTotal(2));
                        app.GzDelayEditField.Value = double(dTotal(3));

                    end

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

                            % Interpolation to update trajectory with initial delays
                            app.GzDelayEditField.Value = 0;
                            dTotal(3) = 0;
                            kTrajCalib = objReco.trajInterpolation(kTrajCalib,dTotal);
                            kTraj = kTrajCalib;

                            % Initial image
                            imCalib = bart(app,['bart nufft -i -l 0.01 ',cSize,' -t'],kTrajCalib,dataCalib);

                            % Initialization
                            iteration = 0;
                            incre = 10;
                            kCalib = objReco.fft2Dmri(imCalib);
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
                                xNew = objReco.lowRankThresh2D(xOld,kSize,rank);
                                rank = rank+0.2;

                                % NUFFT to get updated k-space data
                                kNew = objReco.ifft2Dmri(xNew);
                                dataCalib = bart(app,'bart nufft',kTraj,kNew);
                                kNew  = reshape(dataCalib,[M2-M1+1 size(kTrajCalib,3) nrCoils]);

                                % Partial derivatives
                                [dydtx,dydty] = objReco.partialDerivative2D(app,kTraj,xNew,calibSize);

                                % Direct solver
                                dydt = [real(objReco.vec(dydtx)) real(objReco.vec(dydty)) ; imag(objReco.vec(dydtx)) imag(objReco.vec(dydty))];
                                dStep = ((dydt)'*dydt)\(dydt' * [real(objReco.vec(kNew - y)) ; imag(objReco.vec(kNew - y))]);
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
                                app.GxDelayEditField.Value = double(dTotal(1));
                                app.GyDelayEditField.Value = double(dTotal(2));
                                app.GzDelayEditField.Value = 0;

                                % Interpolation to update trajectory with new delays
                                kTraj = objReco.trajInterpolation(kTrajCalib,dTotal);

                                % The new image with k-space updated for gradient delays
                                imCalib = bart(app,['bart nufft -i -l 0.01 ',cSize,' -t'],kTraj,reshape(y,[1 M2-M1+1 size(kTrajCalib,3) nrCoils]));

                                % Show image
                                im = squeeze(abs(imCalib(:,:)));
                                im = flip(permute(im,[2 1]),2);
                                if app.retroDataPars.PHASE_ORIENTATION
                                    im = rot90(im,-1);
                                    daspect(app.MovieFig,[1 app.retroDataPars.aspectratio 1]);
                                else
                                    daspect(app.MovieFig,[app.retroDataPars.aspectratio 1 1]);
                                end
                                xlim(app.MovieFig, [0 size(im,2)+1]);
                                ylim(app.MovieFig, [0 size(im,1)+1]);
                                imshow(rot90(im),[],'Parent',app.MovieFig);

                                % Calculate k-space from new image
                                xOld = objReco.fft3Dmri(squeeze(imCalib));

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
                trajPics = permute(trajPics,[1 2 3 4 11 12 5 6 7 8 9 10]);
                trajPics = objReco.trajInterpolation(trajPics,dTotal);
                trajPics = ipermute(trajPics,[1 2 3 4 11 12 5 6 7 8 9 10]);
            
                % Density correction
                % app.TextMessage('Calculating density correction ...');
                % onesDense = ones(size(kSpacePicsSum));
                % tmpDense = bart(app,strcat('nufft -d',num2str(dimx),':',num2str(dimx),':1 -a'),trajPicsSum,onesDense);
                % density = bart(app,'nufft ',trajPicsSum,tmpDense);
                % density = density.^(-1/2);
                % density(isnan(density)) = 0;
                % density(isinf(density)) = 0;

                % Sensitivity maps
                if nrCoils > 1 && ESPIRiT
                    kSpacePicsSum = sum(kSpacePics,[11,12]);
                    trajPicsSum = sum(trajPics,[11,12]);
                    ze = squeeze(abs(kSpacePicsSum(1,end,:))) > 0;
                    kSpacePicsSum = kSpacePicsSum(:,:,ze);
                    trajPicsSum = trajPicsSum(:,:,ze);
                    app.TextMessage('Calculating coil sensitivity maps ...');
                    lowResImage = bart(app,['nufft -i -l2 ',cSize,' -t'], trajPicsSum, kSpacePicsSum);
                    lowResKspace = bart(app,'fft -u 7', lowResImage);
                    kSpaceZeroFilled = bart(app,['resize -c 0 ',num2str(dimy),' 1 ',num2str(dimx)], lowResKspace);
                    sensitivities = bart(app,'ecalib -t0.002 -m1', kSpaceZeroFilled);
                else
                    sensitivities = ones(dimx,dimy,1,nrCoils,1,1,1,1,1,1,1,1,1,dimz);
                end

                % 2D radial PICS reconstruction
                app.TextMessage('PICS reconstruction ...');
                picsCommand = 'pics -e ';
                if Wavelet>0
                    picsCommand = [picsCommand,' -RW:6:0:',num2str(Wavelet)];
                end
                if TVxyz>0
                    picsCommand = [picsCommand,' -R',objReco.totalVariation,':6:0:',num2str(TVxyz)];
                end
                if LR>0
                    % Locally low-rank in the spatial domain
                    blockSize = round(dimx/16);  % Block size
                    blockSize(blockSize<8) = 8;
                    picsCommand = [picsCommand,' -RL:6:6:',num2str(LR),' -b',num2str(blockSize),' -N '];
                end
                if TVt>0
                    picsCommand = [picsCommand,' -R',objReco.totalVariation,':1024:0:',num2str(TVt)];
                end
                if TVd>0
                    picsCommand = [picsCommand,' -R',objReco.totalVariation,':2048:0:',num2str(TVd)];
                end

                %igrid = bart(app,picsCommand,'-t',trajPics,'-p',density,kSpacePics,sensitivities);
                igrid = bart(app,picsCommand,'-t',trajPics,kSpacePics,sensitivities);
            
                % Root sum of squares over all coils
                recoImage = bart(app,'rss 8', igrid);

                % Rearrange to correct orientation: frames, z, y, x, dynamics
                imageReg = reshape(recoImage,[dimz,dimy,dimx,dimf,dimd]);
                imageOut = flip(flip(permute(imageReg,[4,3,2,1,5]),3),4);

                % sense map orientations: x, y, slices, map1, map2
                senseMap1 = flip(permute(abs(sensitivities),[2,1,14,3,4,5,6,7,8,9,10,11,12,13]),2);

                % normalize sense map to reasonable value range
                senseMap1 = senseMap1*4095/max(senseMap1(:));

                % shift image in phase-encoding direction if needed
                objReco.movieExp = circshift(imageOut,-objData.pixelshift1,3);
                objReco.senseMap = circshift(senseMap1,-objData.pixelshift1,3);

                app.validRecoFlag = true;
                app.validSenseMapFlag = true;
                app.ProgressGauge.Value = 100;
                drawnow;

            else

                % Reconstruction with NUFFT in Matlab
                app.TextMessage('BART toolbox not available ...');
                app.TextMessage('2D radial reconstruction using NUFFT ...');
                app.ProgressGauge.Value = 0;
                
                dimf = size(objKspace.kSpace{1},1);
                dimx = size(objKspace.kSpace{1},2);
                dimy = dimx;
                dimz = size(objKspace.kSpace{1},4);
                dimd = size(objKspace.kSpace{1},5);
                nrCoils = objData.nr_coils;
                loops = dimz*dimf*dimd*nrCoils;
                app.TextMessage('Slow reconstruction ...');
                app.TextMessage(strcat('Estimated reconstruction time >',{' '},num2str(loops/2),{' '},'min ...'));
                
                traj = objKspace.kSpaceTraj;
                image = zeros(dimf,dimz,dimy,dimx,dimd,nrCoils);
                sensitivities = ones(1,dimy,dimx,1,1,1,1,1,1,1,1,1,14);

                % Gradient delays from app
                dTotal(1) = app.GxDelayEditField.Value;
                dTotal(2) = app.GyDelayEditField.Value;
                dTotal(3) = app.GzDelayEditField.Value;
                app.DataOffsetRadialEditField.Value = 0;

                % Prepare the trajectory
                traj = permute(traj,[6,2,3,4,1,5]);
                traj = objReco.trajInterpolation(traj,dTotal);
                traj = permute(traj,[5,3,2,4,6,1]);
                % frames, spokes, readout, slices, dynamics, 3 coordinates

                % Initialization
                maxit = 5;     % 0 or 1 for gridding, higher values for conjugate gradient
                damp = 0;       % Tikhonov penalty on ||x||
                weight = [];    % data weighting (optional)
                partial = 0.5;  % Tikhobov penalty on ||imag(x))||

                cnt = 1;

                for slice = 1:dimz

                    for dynamic = 1:dimd

                        for frame = 1:dimf

                            % NOTE: coils can be incorporated in the NUFFT, need data first
                            for coil = 1:nrCoils

                                app.ProgressGauge.Value = round(100*cnt/loops);
                                drawnow;

                                om = permute(squeeze(traj(frame,:,:,slice,dynamic,:)),[3 2 1]);
                                obj = nufft_3d(om,dimx);

                                data = squeeze(objKspace.kSpace{coil}(frame,:,:,slice,dynamic));
                                data = data(:);

                                reco = squeeze(obj.iNUFT(data,maxit,damp,weight,'phase-constraint',partial));
                                reco = reshape(reco,[1,1,size(reco),1,1]);

                                image(frame,slice,:,:,dynamic,coil) = reco;

                            end

                            cnt = cnt + 1;

                        end

                    end

                end

                % Root sum of squares over coil dimension and flip
                image = rssq(image,6);
                image = flip(permute(image,[1,4,3,2,5]),3);

                % sense map orientations: x, y, slices, map1, map2
                senseMap1 = permute(sensitivities,[3,2,14,4,5,1,6,7,8,9,10,11,12,13]);

                % normalize sense map to reasonable value range
                senseMap1 = senseMap1*4095/max(senseMap1(:));

                % shift image in phase-encoding direction if needed
                objReco.movieExp = circshift(image,-objData.pixelshift1,3);
                objReco.senseMap = circshift(senseMap1,-objData.pixelshift1,3);

                app.validRecoFlag = true;
                app.validSenseMapFlag = true;
                app.ProgressGauge.Value = 100;
                drawnow;

            end
            
        end % reco2Dradial
        
        
         
        
        % ---------------------------------------------------------------------------------
        % 3D UTE reconstruction with the Bart toolbox
        % ---------------------------------------------------------------------------------
        function objReco = reco3Dute(objReco, objData, objKspace, app)

            if app.bartDetectedFlag

                app.TextMessage('Reconstructing 3D UTE data with the BART toolbox ...');
                app.ProgressGauge.Value = 0;
                drawnow;

                Wavelet = app.WVxyzEditField.Value;
                TVxyz = app.TVxyzEditField.Value;
                LR = app.LLRxyzEditField.Value;
                TVt = app.TVcineEditField.Value;
                TVd = app.TVdynEditField.Value;
                nCoils = objData.nr_coils;

                % Wavelet = wavelet L1-norm regularization factor
                % TVxyz = total variation in xyz-dimension regularization
                % LR = low rank regularization
                % TVt = total variation in time regularization
                % TVd = total variation in dynamic dimension regularization
                % nCoils = number of RF receiver coils
                
                dimf = size(objKspace.kSpace{1},1);
                dimx = size(objKspace.kSpace{1},2);
                dimy = dimx;
                dimz = dimx;
                dimd = size(objKspace.kSpace{1},5);
               
                kSpace = zeros(size(objKspace.kSpace{1}));
                for i = 1:nCoils
                    kSpace(:,:,:,:,:,i) = objKspace.kSpace{i};
                end
                traj = objKspace.kSpaceTraj;
                
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
                kSpacePics = permute(kSpace,[7, 2, 3, 6,14, 8, 9,10,11,12,1, 5, 13, 4]);
                
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
                kdim = round(dimx/3);
                if mod(kdim,2) == 1
                    kdim = kdim + 1;
                end
                kdim(kdim < 32) = 32;
                kdim(kdim > dimx) = dimx;
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
                    kTrajCalib = objReco.trajInterpolation(kTrajCalib,dTotal);
                    kTraj = kTrajCalib;

                    % Initial image
                    dataCalib = dataCalib(:,:,ze);
                    imCalib = bart(app,['bart nufft -i -l 0.01 ',cSize,' -t'],kTrajCalib,dataCalib);
    
                    % Initialization
                    iteration = 0;
                    incre = 10;
                    kCalib = objReco.fft3Dmri(imCalib);
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
                        xNew = objReco.lowRankThresh3D(xOld,kSize,rank);
                        rank = rank+0.2;
                        rank(rank>prod(kSize)) = prod(kSize);

                        % NUFFT to get updated k-space data
                        kNew = objReco.ifft3Dmri(xNew);
                        dataCalib = bart(app,'bart nufft',kTraj,kNew);
                        kNew  = reshape(dataCalib,[M size(kTrajCalib,3) nCoils]);

                        % Partial derivatives
                        [dydtx,dydty,dydtz] = objReco.partialDerivative3D(app,kTraj,xNew,calibSize);

                        % Direct solver
                        dydt = [real(objReco.vec(dydtx)) real(objReco.vec(dydty)) real(objReco.vec(dydtz)) ; imag(objReco.vec(dydtx)) imag(objReco.vec(dydty)) imag(objReco.vec(dydtz))];
                        dStep = ((dydt)'*dydt)\(dydt' * [real(objReco.vec(kNew - y)) ; imag(objReco.vec(kNew - y))]); 
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
                        kTraj = objReco.trajInterpolation(kTrajCalib,dTotal);

                        % The new image with k-space updated for gradient delays
                        imCalib = bart(app,['bart nufft -i -l 0.01 ',cSize,' -t'],kTraj,reshape(y,[1 M size(kTrajCalib,3) nCoils]));

                        % Show image
                        im = squeeze(abs(imCalib(:,:,round(calibSize(3)/2),1)));
                        im = flip(flip(im,1),2);
                        if app.retroDataPars.PHASE_ORIENTATION
                            im = rot90(im,-1);
                            daspect(app.MovieFig,[1 app.retroDataPars.aspectratio 1]);
                        else
                            daspect(app.MovieFig,[app.retroDataPars.aspectratio 1 1]);
                        end
                        xlim(app.MovieFig, [0 size(im,2)+1]);
                        ylim(app.MovieFig, [0 size(im,1)+1]);
                        imshow(rot90(im),[],'Parent',app.MovieFig);
                
                        % Calculate k-space from new image
                        xOld = objReco.fft3Dmri(squeeze(imCalib));

                    end

                    % Hide stop button
                    app.MovieStartButton.Enable = 'off';
                    app.MovieStopButton.Enable = 'off';

                end % Gradient calibration

                % Final gradient delay correction from optimization or values from app
                trajPics = permute(trajPics,[1 2 3 4 11 12 5 6 7 8 9 10]);
                trajPics = objReco.trajInterpolation(trajPics,dTotal);
                trajPics = ipermute(trajPics,[1 2 3 4 11 12 5 6 7 8 9 10]);
            
                % Coil sensitivities from sum of all frames and dynamics
                if nCoils > 1
                    kSpacePicsSum = sum(kSpacePics,[11,12]);
                    trajPicsSum = sum(trajPics,[11,12]);
                    ze = squeeze(abs(kSpacePicsSum(1,end,:))) > 0;
                    kSpacePicsSum = kSpacePicsSum(:,:,ze);
                    trajPicsSum = trajPicsSum(:,:,ze);
                    app.TextMessage('Calculating coil sensitivity maps ...');
                    lowResImage = bart(app,['nufft -i ',cSize,' -t'], trajPicsSum, kSpacePicsSum);
                    lowResKspace = bart(app,'fft -u 7', lowResImage);
                    kSpaceZeroFilled = bart(app,['resize -c 0 ',num2str(dimz),' 1 ',num2str(dimy),' 2 ',num2str(dimx)], lowResKspace);
                    sensitivities = bart(app,'ecalib -S -t0.0005 -m1', kSpaceZeroFilled);
                else
                    sensitivities = ones(dimz,dimy,dimx,1,1,1,1,1,1,1,1,1,1,1);
                end

                % Density correction
                app.TextMessage('Calculating density correction ...');
                denseOnes = ones(size(kSpacePics));
                denseTmp = bart(app,strcat('nufft -d',num2str(dimx),':',num2str(dimx),':',num2str(dimx),' -a'),trajPics,denseOnes);
                density = bart(app,'nufft ',trajPics,denseTmp);
                density = density.^(-1/3);
                density(isnan(density)) = 0;
                density(isinf(density)) = 0;

                % PICS reconstruction
                app.TextMessage('PICS reconstruction ...');
                picsCommand = 'pics -i10 ';
                if Wavelet>0
                    picsCommand = [picsCommand,' -RW:7:0:',num2str(Wavelet)];
                end
                if TVxyz>0
                    picsCommand = [picsCommand,' -R',objReco.totalVariation,':7:0:',num2str(TVxyz)];
                end
                if LR>0
                    % Locally low-rank in the spatial domain
                    blockSize = round(dimx/16);  % Block size
                    blockSize(blockSize<8) = 8;
                    picsCommand = [picsCommand,' -RL:7:7:',num2str(LR),' -b',num2str(blockSize)];
                end
                if TVt>0
                    picsCommand = [picsCommand,' -R',objReco.totalVariation,':1024:0:',num2str(TVt)];
                end
                if TVd>0
                    picsCommand = [picsCommand,' -R',objReco.totalVariation,':2048:0:',num2str(TVd)];
                end
                igrid = bart(app,picsCommand,'-t',trajPics,'-p',density,kSpacePics,sensitivities);    
    
                % Root sum of squares over all coils
                recoImage = bart(app,'rss 8', igrid);

                % Rearrange to correct orientation: frames, x, y, z, dynamics
                imageReg = reshape(recoImage,[dimz,dimy,dimx,dimf,dimd]);
                imageOut = flip(flip(permute(imageReg,[4,1,2,3,5]),3),4);
                
                % Sense map orientations: x, y, z, map1, map2
                senseMap1 = flip(flip(permute(abs(sensitivities),[3,2,1,4,5,6,7,8,9,10,11,12,13,14]),2),3);
                
                % Normalize sense map to reasonable value range
                senseMap1 = senseMap1*4095/max(senseMap1(:));
                
                % Shift image in phase-encoding directions if needed
                objReco.movieExp = circshift(imageOut,-objData.pixelshift1,3);
                objReco.movieExp = circshift(imageOut,-objData.pixelshift2,4);
                objReco.senseMap = circshift(senseMap1,-objData.pixelshift1,3);
                objReco.senseMap = circshift(senseMap1,-objData.pixelshift2,4);

                app.validRecoFlag = true;
                app.validSenseMapFlag = true;
                app.ProgressGauge.Value = 100;
                drawnow;

            else

                app.TextMessage('BART toolbox not available ...');
                app.TextMessage('3D UTE reconstruction using 3D NUFFT ...');
                app.ProgressGauge.Value = 0;

                dimf = size(objKspace.kSpace{1},1);
                dimx = size(objKspace.kSpace{1},2);
                dimy = dimx;
                dimz = dimx;
                dimd = size(objKspace.kSpace{1},5);
                nCoils = objData.nr_coils;
                loops = dimf*dimd*nCoils;
                app.TextMessage('Slow reconstruction ...');
                app.TextMessage(strcat('Estimated reconstruction time >',{' '},num2str(loops*2),{' '},'min ...'));
                
                traj = objKspace.kSpaceTraj;
                image = zeros(dimf,dimx,dimy,dimz,dimd,nCoils);
                sensitivities = 4095*ones(dimz,dimy,dimx,1,1,1,1,1,1,1,1,1,1,1);
              
                maxit = 10;     % 0 or 1 for gridding, higher values for conjugate gradient
                damp = 0;       % Tikhonov penalty on ||x||
                weight = [];    % data weighting (optional)
                partial = 0.5;  % Tikhobov penalty on ||imag(x))||

                cnt = 0;

                for dynamic = 1:dimd

                    for frame = 1:dimf

                        % NOTE: coils can be incorporated in the NUFFT, need data first
                        for coil = 1:nCoils 

                            app.ProgressGauge.Value = round(100*cnt/loops);
                            drawnow;

                            om = permute(squeeze(traj(frame,:,:,1,dynamic,:)),[3,2,1]);
                            obj = nufft_3d(om,dimx);

                            data = squeeze(objKspace.kSpace{coil}(frame,:,:,1,dynamic));
                            data = permute(data,[2 1]);
                            data = data(:);
                       
                            image(frame,:,:,:,dynamic,coil) = obj.iNUFT(data,maxit,damp,weight,'phase-constraint',partial);
                    
                        end

                        cnt = cnt + 1;

                    end

                end

                % Root sum of squares over coil dimension and flip
                image = rssq(image,6);
                image = flip(image,3);
    
                % Shift image in phase-encoding directions if needed
                objReco.movieExp = circshift(image,-objData.pixelshift1,3);
                objReco.movieExp = circshift(image,-objData.pixelshift2,4);
                objReco.senseMap = circshift(sensitivities,-objData.pixelshift1,3);
                objReco.senseMap = circshift(sensitivities,-objData.pixelshift2,4);
                
                app.validRecoFlag = true;
                app.validSenseMapFlag = true;
                app.ProgressGauge.Value = 100;
                drawnow;

            end

        end % reco3Dute



        
        % ---------------------------------------------------------------------------------
        % Normalize movie intensity
        % ---------------------------------------------------------------------------------
        function objReco = normImages(objReco,objData)
            
            % Normalize the images to 2^15 range

            objReco.rescaleIntercept = 0;                           % Dicom rescale intercept
            objReco.rescaleSlope = max(objReco.movieExp(:));        % Dicom rescale slope
   
            objReco.movieExp = round(32765*objReco.movieExp/objReco.rescaleSlope);
            objReco.movieApp = objReco.movieExp;

            if objData.PHASE_ORIENTATION == 0
                objReco.movieApp = flip(objReco.movieApp,2);
            end
            
        end % normImages
        
        

        % ---------------------------------------------------------------------------------
        % Determine whether movie has multiple slices and/or dynamics
        % ---------------------------------------------------------------------------------
        function objReco = determineMultiDimensions(objReco)
            
            % Multislice
            if size(objReco.movieApp,4) > 1
                objReco.multiSliceFlag = true;
            else
                objReco.multiSliceFlag = false;
            end

            % Multidynamic
            if size(objReco.movieApp,5) > 1
                objReco.multiDynamicFlag = true;
            else
                objReco.multiDynamicFlag = false;
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
        function obj = recoLcurve(obj,objData,objKspace,app)

            % Range of lamda's
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
                        obj.sos = 0;
                        obj = obj.reco2D(objData, objKspace, app);
                        Lmovie = obj.movieExp;

                        % L1 norm wavelet
                        WVxyz = bart(app,'cdf97 6',Lmovie);
                        y(i) = norm(WVxyz(:),1);

                        % L2 norm
                        kSpaceIm = bart(app,'fft -u -i 6',Lmovie);
                        kSpaceIm = bart(app,['resize -c 1 ',num2str(size(objKspace.kSpace{1},2)),' 2 ',num2str(size(objKspace.kSpace{1},3))],kSpaceIm);
                        x(i) = 0;
                        for j = 1:objData.nr_coils
                            mask = objKspace.kSpace{j} ~= 0;
                            kSpaceDiff = abs(kSpaceIm(:,:,:,:,:,j).*mask) - abs(objKspace.kSpace{j});
                            x(i) = x(i) + norm(kSpaceDiff(:),2);
                        end

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
                        obj.sos = 0;
                        obj = obj.reco2D(objData, objKspace, app);
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
                        kSpaceIm = bart(app,['resize -c 1 ',num2str(size(objKspace.kSpace{1},2)),' 2 ',num2str(size(objKspace.kSpace{1},3))],kSpaceIm);
                        x(i) = 0;
                        for j = 1:objData.nr_coils
                            mask = objKspace.kSpace{j} ~= 0;
                            kSpaceDiff = abs(kSpaceIm(:,:,:,:,:,j).*mask) - abs(objKspace.kSpace{j});
                            x(i) = x(i) + norm(kSpaceDiff(:),2);
                        end

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
                        obj.sos = 0;
                        obj = obj.reco2D(objData, objKspace, app);
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
                        kSpaceIm = bart(app,['resize -c 1 ',num2str(size(objKspace.kSpace{1},2)),' 2 ',num2str(size(objKspace.kSpace{1},3))],kSpaceIm);
                        x(i) = 0;
                        for j = 1:objData.nr_coils
                            mask = objKspace.kSpace{j} ~= 0;
                            kSpaceDiff = abs(kSpaceIm(:,:,:,:,:,j).*mask) - abs(objKspace.kSpace{j});
                            x(i) = x(i) + norm(kSpaceDiff(:),2);
                        end

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
                        obj.sos = 0;
                        obj = obj.reco2D(objData, objKspace, app);
                        Lmovie = obj.movieExp;

                        % L1 norm TV
                        TVcine = circshift(Lmovie,1) - Lmovie;
                        y(i) = norm(TVcine(:),1);

                        % L2 norm
                        kSpaceIm = bart(app,'fft -u -i 6',Lmovie);
                        kSpaceIm = bart(app,['resize -c 1 ',num2str(size(objKspace.kSpace{1},2)),' 2 ',num2str(size(objKspace.kSpace{1},3))],kSpaceIm);
                        x(i) = 0;
                        for j = 1:objData.nr_coils
                            mask = objKspace.kSpace{j} ~= 0;
                            kSpaceDiff = abs(kSpaceIm(:,:,:,:,:,j).*mask) - abs(objKspace.kSpace{j});
                            x(i) = x(i) + norm(kSpaceDiff(:),2);
                        end

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
                        obj.sos = 0;
                        obj = obj.reco2D(objData, objKspace, app);
                        Lmovie = obj.movieExp;

                        % L1 norm TV in dynamic dimension
                        TVtime = circshift(Lmovie,5) - Lmovie;
                        y(i) = norm(TVtime(:),1);

                        % L2 norm
                        kSpaceIm = bart(app,'fft -u -i 6',Lmovie);
                        kSpaceIm = bart(app,['resize -c 1 ',num2str(size(objKspace.kSpace{1},2)),' 2 ',num2str(size(objKspace.kSpace{1},3))],kSpaceIm);
                        x(i) = 0; %#ok<*AGROW>
                        for j = 1:objData.nr_coils
                            mask = objKspace.kSpace{j} ~= 0;
                            kSpaceDiff = abs(kSpaceIm(:,:,:,:,:,j).*mask) - abs(objKspace.kSpace{j});
                            x(i) = x(i) + norm(kSpaceDiff(:),2);
                        end

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



    end % public methods




    % ---------------------------------------------------------------------------------
    % Static methods
    % ---------------------------------------------------------------------------------
    methods (Static)



        % ---------------------------------------------------------------------------------
        % 3D partial derivatives
        % ---------------------------------------------------------------------------------
        function [dydtx,dydty,dydtz] = partialDerivative3D(app,kTraj,xNew,calibSize)

            nCoils = size(xNew,4);

            kx = squeeze(kTraj(1,:,:));
            ky = squeeze(kTraj(2,:,:));
            kz = squeeze(kTraj(3,:,:));

            dkx = zeros(size(kx));
            dky = zeros(size(ky));
            dkz = zeros(size(kz));

            dkx(2:end,:) = kx(2:end,:)- kx(1:end-1,:);
            dky(2:end,:) = ky(2:end,:)- ky(1:end-1,:);
            dkz(2:end,:) = kz(2:end,:)- kz(1:end-1,:);
            xNewFFT = 1j*retroReco.ifft3Dmri(xNew);

            repX = repmat([0:(calibSize(1)/2-1), 0, -calibSize(1)/2+1:-1]'/calibSize(1),[1 calibSize(1) calibSize(1) nCoils]);
            repY = permute(repX,[2 1 3]);
            repZ = permute(repX,[3 1 2]);

            tmp = xNewFFT.*repX;
            tmpCalib = bart(app,'bart nufft',kTraj,reshape(tmp,[calibSize nCoils]));
            dydkx = reshape(tmpCalib,[size(kx) nCoils]);

            tmp = xNewFFT.*repY;
            tmpCalib = bart(app,'bart nufft',kTraj,reshape(tmp,[calibSize nCoils]));
            dydky = reshape(tmpCalib,[size(kx) nCoils]);

            tmp = xNewFFT.*repZ;
            tmpCalib = bart(app,'bart nufft',kTraj,reshape(tmp,[calibSize nCoils]));
            dydkz = reshape(tmpCalib,[size(kx) nCoils]);

            dydtx = dydkx.*repmat(dkx,[1 1 1 nCoils]);
            dydty = dydky.*repmat(dky,[1 1 1 nCoils]);
            dydtz = -dydkz.*repmat(dkz,[1 1 1 nCoils]); % positive does not converge

            dydtx(isnan(dydtx)) = 0;
            dydty(isnan(dydty)) = 0;
            dydtz(isnan(dydtz)) = 0;

        end % partialDerivative3D




        % ---------------------------------------------------------------------------------
        % 2D partial derivatives
        % ---------------------------------------------------------------------------------
        function [dydtx,dydty] = partialDerivative2D(app,kTraj,Xnew,calibSize)

            nCoils = size(Xnew,3);

            kx = squeeze(kTraj(1,:,:));
            ky = squeeze(kTraj(2,:,:));

            dkx = zeros(size(kx));
            dky = zeros(size(ky));

            dkx(2:end,:) = kx(2:end,:)- kx(1:end-1,:);
            dky(2:end,:) = ky(2:end,:)- ky(1:end-1,:);

            tmp =(1j*retroReco.ifft2Dmri(Xnew).*repmat([0:(calibSize(1)/2-1), 0, -calibSize(1)/2+1:-1]'/calibSize(1),[1 calibSize(2) nCoils]));
            tmpCalib = bart(app,'bart nufft',kTraj,reshape(tmp,[calibSize 1 nCoils]));
            dydkx = reshape(tmpCalib,[size(kx) nCoils]);

            tmp = (1j*retroReco.ifft2Dmri(Xnew).*repmat([0:(calibSize(1)/2-1), 0, -calibSize(1)/2+1:-1]/calibSize(1),[calibSize(1) 1 nCoils]));
            tmpCalib = bart(app,'bart nufft',kTraj,reshape(tmp,[calibSize 1 nCoils]));
            dydky = reshape(tmpCalib,[size(kx) nCoils]);

            dydtx = dydkx.*repmat(dkx,[1 1 nCoils]);
            dydty = dydky.*repmat(dky,[1 1 nCoils]);

            dydtx(isnan(dydtx)) = 0;
            dydty(isnan(dydty)) = 0;

        end % partialDerivative2D
        
        
        
        % ---------------------------------------------------------------------------------
        % Low rank threshold 3D
        % ---------------------------------------------------------------------------------
        function Xnew = lowRankThresh3D(Xold,kSize,thresh)

            thresh = round(thresh);

            keep = 1:thresh;

            [sx,sy,sz,~] = size(Xold);
            tmp = retroReco.im2row3D(Xold,kSize);
            [tsx,tsy,Nc] = size(tmp);
            A = reshape(retroReco.im2row3D(Xold,kSize),tsx,tsy*Nc);

            [U,S,V] = svd(A,'econ');
            A = U(:,keep)*S(keep,keep)*V(:,keep)';

            A = reshape(A,tsx,tsy*Nc);

            Xnew = retroReco.row2im3D(A,[sx,sy,sz,Nc],kSize);

        end % lowRankThresh3D



        % ---------------------------------------------------------------------------------
        % Low rank threshold 2D
        % ---------------------------------------------------------------------------------
        function Xnew = lowRankThresh2D(Xold,kSize,thresh)

            thresh = round(thresh);

            keep = 1:thresh;
       
            [sx,sy,Nc] = size(Xold);
            tmp = retroReco.im2row2D(Xold,kSize); 

            [tsx,tsy,tsz] = size(tmp);
            A = reshape(retroReco.im2row2D(Xold,kSize),tsx,tsy*tsz);

            [U,S,V] = svd(A,'econ');
            A = U(:,keep)*S(keep,keep)*V(:,keep)';

            A = reshape(A,tsx,tsy,tsz);       

            Xnew = retroReco.row2im2D(A,[sx,sy,Nc],kSize);

        end % lowRankThresh2D



        % ---------------------------------------------------------------------------------
        % Trajectory interpolation 
        % ---------------------------------------------------------------------------------
        function kSpaceUpdate = trajInterpolation(kSpaceOld,dShift)

            kSpaceUpdate = zeros(size(kSpaceOld));

            for idx6 = 1:size(kSpaceOld,6) % dynamics

                for idx5 = 1:size(kSpaceOld,5) % frames

                    for idx4 = 1:size(kSpaceOld,4) % slices

                        for idx3 = 1:size(kSpaceOld,3) % spokes

                            kx = interp1((1:size(kSpaceOld,2))+dShift(1),kSpaceOld(1,:,idx3,idx4,idx5,idx6),1:size(kSpaceOld,2),'linear'); %Kx
                            ky = interp1((1:size(kSpaceOld,2))+dShift(2),kSpaceOld(2,:,idx3,idx4,idx5,idx6),1:size(kSpaceOld,2),'linear'); %Ky
                            kz = interp1((1:size(kSpaceOld,2))+dShift(3),kSpaceOld(3,:,idx3,idx4,idx5,idx6),1:size(kSpaceOld,2),'linear'); %Kz

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

                            kSpaceUpdate(1,:,idx3,idx4,idx5,idx6) = kx(:);
                            kSpaceUpdate(2,:,idx3,idx4,idx5,idx6) = ky(:);
                            kSpaceUpdate(3,:,idx3,idx4,idx5,idx6) = kz(:);

                        end

                    end

                end

            end

        end % kSpaceInterpolation




        % ---------------------------------------------------------------------------------
        % image to rows 3D
        % ---------------------------------------------------------------------------------
        function res = im2row3D(im, winSize)

            [sx,sy,sz,nc] = size(im);
            res = zeros((sx-winSize(1)+1)*(sy-winSize(2)+1)*(sz-winSize(3)+1),prod(winSize),nc);

            count=0;
            for z=1:winSize(3)
                for y=1:winSize(2)
                    for x=1:winSize(1)
                        count = count+1;
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

            [sx,sy,sz] = size(im);
            res = zeros((sx-winSize(1)+1)*(sy-winSize(2)+1),prod(winSize),sz);

            count=0;
            for y=1:winSize(2)
                for x=1:winSize(1)
                    count = count+1;
                    res(:,count,:) = reshape(im(x:sx-winSize(1)+x,y:sy-winSize(2)+y,:),...
                        (sx-winSize(1)+1)*(sy-winSize(2)+1),1,sz);
                end
            end

        end % im2row2D



        % ---------------------------------------------------------------------------------
        % rows to image 3D
        % ---------------------------------------------------------------------------------
        function [res,W] = row2im3D(mtx, imSize, winSize)

            nCoils = size(mtx,4);
            sx = imSize(1);
            sy = imSize(2);
            sz = imSize(3);
            res = zeros(imSize(1),imSize(2),imSize(3),nCoils);
            W = res;

            count=0;
            for z=1:winSize(3)
                for y=1:winSize(2)
                    for x=1:winSize(1)
                        count = count+1;
                        res(x:sx-winSize(1)+x,y:sy-winSize(2)+y,z:sz-winSize(3)+z,:) = res(x:sx-winSize(1)+x,y:sy-winSize(2)+y,z:sz-winSize(3)+z,:) + reshape(mtx(:,count,:,:),(sx-winSize(1)+1),(sy-winSize(2)+1),(sz-winSize(3)+1),nCoils);
                        W(x:sx-winSize(1)+x,y:sy-winSize(2)+y,z:sz-winSize(3)+z,:) = W(x:sx-winSize(1)+x,y:sy-winSize(2)+y,z:sz-winSize(3)+z,:)+1;
                    end
                end
            end

            res = res./W;

        end % row2im3D


        % ---------------------------------------------------------------------------------
        % rows to image 2D
        % ---------------------------------------------------------------------------------
        function [res,W] = row2im2D(mtx,imSize, winSize)

            sz = size(mtx,3);
            sx = imSize(1); 
            sy = imSize(2);
            res = zeros(imSize(1),imSize(2),sz);
            W = res;

            count=0;
            for y=1:winSize(2)
                for x=1:winSize(1)
                    count = count+1;
                    res(x:sx-winSize(1)+x,y:sy-winSize(2)+y,:) = res(x:sx-winSize(1)+x,y:sy-winSize(2)+y,:) + reshape(mtx(:,count,:),(sx-winSize(1)+1),(sy-winSize(2)+1),sz);
                    W(x:sx-winSize(1)+x,y:sy-winSize(2)+y,:) = W(x:sx-winSize(1)+x,y:sy-winSize(2)+y,:)+1;
                end
            end

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






    end % Static methods


end % retroReco

