classdef retroData
    
    % Data and parameter class for retrospective app
    
    properties
        
        % Data
        data
        mrdFooter
        newMrdFooter
        rprFile
        newRprFile
        filename
        exportDir
        includeWindow                                                  
        excludeWindow 
        nrKlines
        
        % Sequence parameters mostly from MRD file header
        PPL
        NO_SAMPLES = 1
        NO_VIEWS = 1
        NO_VIEWS_2 = 1
        EXPERIMENT_ARRAY = 1
        nr_repetitions = 1
        NO_AVERAGES = 1
        NO_SLICES = 1
        SLICE_THICKNESS = 1
        SLICE_SEPARATION = 1
        SLICE_INTERLEAVE = 1
        r_angle_var
        p_angle_var
        s_angle_var
        nr_coils = 1
        FOV = 30
        PHASE_ORIENTATION = 0
        FOVf = 8
        aspectratio = 1
        alpha = 20
        te = 2
        te_us = 0
        TE
        tr = 10
        tr_extra_us = 0
        TR
        ti = 1000
        VFA_angles = []
        VFA_size = 0
        frame_loop_on
        radial_on = 0
        slice_nav = 0
        date
        pixelshift1 = 0
        pixelshift2 = 0
        coil_scaling = 1
        scanner = 'MRS'
        no_samples_nav = 10
        
        % K-space trajectory related
        pe1_order = 3
        pe2_centric_on = 1
        pe2_traj = 0
        gp_var_mul
        
        % Navigator related
        primaryNavigatorPoint = 10
        nrNavPointsDiscarded = 35
        nrNavPointsUsed = 5
        
        % Final dimensions
        dimx
        dimy
        dimz
        nrKsteps
        
        % Flags
        rprFlag = false;
        validDataFlag = false;
        multiCoilFlag = false
        multi2DFlag = false
        vfaDataFlag = false
        
        % Data and reconstruction type
        dataType = '2D'
        recoGuess = 'systolic function'
        
    end
    
    
    
    methods
        
        % ---------------------------------------------------------------------------------
        % Object constructor
        % ---------------------------------------------------------------------------------
        function obj = retroData(parameter, mridata)
            
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
                end
                
                if isfield(parameter,'pixelshift2')
                    obj.pixelshift2 = parameter.pixelshift2;
                end
                
                if isfield(parameter,'pe2_traj')
                    obj.pe2_traj = parameter.pe2_traj;
                end
                
                if isfield(parameter,'scanner')
                    obj.scanner = parameter.scanner;
                end
                
            end
            
        end
        
        
        
        
        % ---------------------------------------------------------------------------------
        % Check whether there are sufficient experiments to peform a valid reconstruction
        % ---------------------------------------------------------------------------------
        function obj = checkNumberOfExperiments(obj, app)
            
            if obj.EXPERIMENT_ARRAY * obj.NO_VIEWS * obj.NO_VIEWS_2 * obj.NO_SLICES < 100
                obj.validDataFlag = false;
                app.TextMessage('ERROR: Not enough k-lines for any reconstruction ...');
            end
            
        end
        
        
        

        % ---------------------------------------------------------------------------------
        % Check the number of averages
        % ---------------------------------------------------------------------------------
        function obj = checkNumberOfAverages(obj, app)
            
            if obj.NO_AVERAGES > 1
                obj.validDataFlag = false;
                app.TextMessage('ERROR: Number of averages should be 1 ...');
            end
            
        end
        
        

        
        % ---------------------------------------------------------------------------------
        % Check for multi-slab data
        % ---------------------------------------------------------------------------------
        function obj = checkForMultiSlab(obj, app)
            
            if strcmp(obj.dataType,'3D') && obj.NO_SLICES > 1
                obj.validDataFlag = false;
                app.TextMessage('ERROR: Only 3D single-slab data supported ...');
            end
            
        end
        


        % ---------------------------------------------------------------------------------
        % Check for TR > 0
        % ---------------------------------------------------------------------------------
        function obj = checkTR(obj, app)
            
            if obj.tr == 0
                obj.validDataFlag = false;
                app.TextMessage('ERROR: TR = 0 is not allowed ...');
            end
            
        end
        
        

        % ---------------------------------------------------------------------------------
        % Check for variable flip angle data
        % ---------------------------------------------------------------------------------
        function obj = checkForVFA(obj, app)
            
            if strcmp(obj.dataType,'3D') && obj.VFA_size > 0
                obj.vfaDataFlag = true;
                obj = setVariableFlipAngles(obj);
                app.TextMessage(strcat('INFO:',{' '},num2str(obj.VFA_size),{' '},'flip angles detected ...'));
            end
            
        end
        
        

        
        % ---------------------------------------------------------------------------------
        % Set the variable flip angles, sort the angles in groups
        % ---------------------------------------------------------------------------------
        function obj = setVariableFlipAngles(obj)
            
            a = obj.VFA_size;
            b = unique(obj.VFA_angles(1:obj.VFA_size),'Stable');
            c = length(b);
            d = obj.VFA_angles(1:obj.VFA_size);
            if a == c
                nrFlipAngles = a;     % FA1, FA2, FA3, FA4, .... = dym1, dyn2, dyn3, dyn4, ...
                lsFlipAngles = b;
            elseif mod(a,c) == 0
                nrFlipAngles = c;     % FA1, FA1, ..., FA2, FA2, ..., FA3, FA3, ... = dym1, dyn2, dyn3, dyn4, ...
                lsFlipAngles = b;
            else
                nrFlipAngles = a;     % each dynamic has its own flip-angle
                lsFlipAngles = d;
            end
            obj.VFA_size = nrFlipAngles;
            obj.VFA_angles = lsFlipAngles;
            
        end
        
        
        

        % ---------------------------------------------------------------------------------
        % Check if the data is a valid file with a motion navigator acquisition
        % ---------------------------------------------------------------------------------
        function obj = acquisitionType(obj, app)
            
            
            if obj.NO_VIEWS == 1 && obj.NO_VIEWS_2 == 1 && obj.EXPERIMENT_ARRAY > 1000

                obj.validDataFlag = true;

                obj.dataType = '3Dute';

                for i=1:obj.nr_coils
                    %                             spokes,1,1,X
                    obj.data{i} = permute(obj.data{i},[1,4,3,2]);
                end

                obj.primaryNavigatorPoint = 1;
                obj.nr_repetitions = size(obj.data{1},1);

            elseif obj.radial_on == 1
                
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
                    end
                end

                kSpaceSum = squeeze(sum(abs(obj.data{1}),[1,2,3]));
                [~,kCenter] = max(kSpaceSum);
                obj.primaryNavigatorPoint = kCenter;

                obj.nr_repetitions = size(obj.data{1},1);
                
            elseif (obj.slice_nav == 1) && (obj.no_samples_nav > 0)
                
                obj.validDataFlag = true;
                
                if obj.NO_VIEWS_2 > 1
                    
                    % 3D data
                    obj.dataType = '3D';                                        
                    for i=1:obj.nr_coils
                        if obj.EXPERIMENT_ARRAY > 1
                            %                                 NR Z Y X  
                            obj.data{i} = permute(obj.data{i},[1,3,2,4]);
                        else
                            %                                 NR Z Y X 
                            obj.data{i} = permute(obj.data{i},[4,2,1,3]);
                        end
                    end
                    
                else

                    % 2D single-slice data
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

                obj.primaryNavigatorPoint = obj.no_samples_nav;
                if obj.nrNavPointsUsed > obj.primaryNavigatorPoint
                    obj.nrNavPointsUsed = obj.primaryNavigatorPoint;
                end
                obj.nr_repetitions = size(obj.data{1},1);
                
            else
                
                obj.validDataFlag = false;

            end

            % Message the user on the type of data
            switch obj.dataType

                case '2D'
                    app.TextMessage('2D single-slice data ...');
                case '3D'
                    app.TextMessage('3D data ...');
                case '2Dms'
                    app.TextMessage('2D multi-slice data ...');
                case '2Dradial'
                    app.TextMessage('2D radial data ...');
                case '3Dute'
                    app.TextMessage('3D UTE data ...');
            
            end

        end
        
        

        
        % ---------------------------------------------------------------------------------
        % Guess which reconstruction type (systolic function, diastolic
        % function, scout, VFA
        % ---------------------------------------------------------------------------------
        function obj = guessRecoType(obj)

            switch obj.dataType

                case {'2D','2Dradial'}

                    if obj.nr_repetitions > 200
                        obj.recoGuess = 'diastolic function';
                    else
                        obj.recoGuess = 'systolic function';
                    end

                case '2Dms'

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

                case {'3D','3Dute'}

                    obj.recoGuess = 'systolic function';

                    if obj.vfaDataFlag
                        obj.recoGuess = 'variable flip-angle';
                    end

            end

        end
        



        % ---------------------------------------------------------------------------------
        % Read MRD file
        % ---------------------------------------------------------------------------------
        function [im,dim,par,unsortedkspace] = importMRD(obj, filename, reordering1, reordering2)

            % Description: Function to open multidimensional MRD/SUR files given a filename with PPR-parsing
            % Read in MRD and SUR file formats
            % Inputs: string filename, reordering1, reordering2
            % reordering1, 2 is 'seq' or 'cen'
            % reordering1 is for 2D (views)
            % reordering2 is for 3D (views2)
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
          
            num2read = no_expts*no_echoes*no_slices*no_views_2*no_views*no_samples*iscomplex; %*datasize;
            [m_total, count] = fread(fid,num2read,dataformat); % reading all the data at once
       
            if iscomplex == 2
                a=1:count/2;
                m_real = m_total(2*a-1);
                m_imag = m_total(2*a);
                clear m_total;
                m_C = m_real+m_imag*1i;
                clear m_real m_imag;
            else
                m_C = m_total;
                clear m_total;
            end

            unsortedkspace = m_C;

            n=0;
            % shaping the data manually:
            ord=1:no_views;
            if strcmp(reordering1,'cen')
                for g=1:no_views/2
                    ord(2*g-1)=no_views/2+g;
                    ord(2*g)=no_views/2-g+1;
                end
            end

            ord1 = 1:no_views_2;
            ord2 = ord1;
            if strcmp(reordering2,'cen')
                for g=1:no_views_2/2
                    ord2(2*g-1)=no_views_2/2+g;
                    ord2(2*g)=no_views_2/2-g+1;
                end
            end

            % pre-allocate the data matrix
            m_C_1=zeros(no_expts,no_echoes,no_slices,max(ord(:)),max(ord2(:)),no_samples);
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
                PPR_keywords = {'BUFFER_SIZE','DATA_TYPE','DECOUPLE_FREQUENCY','DISCARD','DSP_ROUTINE','EDITTEXT','EXPERIMENT_ARRAY','FOV','FOV_READ_OFF','FOV_PHASE_OFF','FOV_SLICE_OFF','GRADIENT_STRENGTH','MULTI_ORIENTATION','Multiple Receivers','NO_AVERAGES','NO_ECHOES','NO_RECEIVERS','NO_SAMPLES','NO_SLICES','NO_VIEWS','NO_VIEWS_2','OBLIQUE_ORIENTATION','OBSERVE_FREQUENCY','ORIENTATION','PHASE_CYCLE','READ/PHASE/SLICE_SELECTION','RECEIVER_FILTER','SAMPLE_PERIOD','SAMPLE_PERIOD_2','SCROLLBAR','SLICE_BLOCK','SLICE_FOV','SLICE_INTERLEAVE','SLICE_THICKNESS','SLICE_SEPARATION','SPECTRAL_WIDTH','SWEEP_WIDTH','SWEEP_WIDTH_2','VAR_ARRAY','VIEW_BLOCK','VIEWS_PER_SEGMENT','SMX','SMY','SWX','SWY','SMZ','SWZ','VAR','PHASE_ORIENTATION','X_ANGLE','Y_ANGLE','Z_ANGLE','PPL','IM_ORIENTATION','IM_OFFSETS'};
                %PPR_type_0 keywords have text fields only, e.g. ":PPL C:\ppl\smisim\1ge_tagging2_1.PPL"
                PPR_type_0 = [23 53];
                %PPR_type_1 keywords have single value, e.g. ":FOV 300"
                PPR_type_1 = [8 42:47];
                %PPR_type_2 keywords have single variable and single value, e.g. ":NO_SAMPLES no_samples, 16"
                PPR_type_2 = [4 7 9:11 15:21 25 31 33 41 49];
                PPR_type_3 = 48; % VAR keyword only (syntax same as above)
                PPR_type_4 = [28 29]; % :SAMPLE_PERIOD sample_period, 300, 19, "33.3 KHz  30 ?s" and SAMPLE_PERIOD_2 - read the first number=timeincrement in 100ns
                %PPR_type_5 keywords have single variable and two values, e.g. ":SLICE_THICKNESS gs_var, -799, 100"
                PPR_type_5 = [34 35];
                % KEYWORD [pre-prompt,] [post-prompt,] [min,] [max,] default, variable [,scale] [,further parameters ...];
                PPR_type_6 = [39 50:52]; % VAR_ARRAY and angles keywords
                PPR_type_7 = [54 55]; % IM_ORIENTATION and IM_OFFSETS (SUR only)

                par = struct('filename',filename);
                for j=1:size(cell_text{1},1)
                    char1 = char(cell_text{1}(j,:));
                    field_ = '';
                    if ~isempty(char1)
                        C = textscan(char1, '%*c%s %s', 1);
                        field_ = char(C{1});
                    end
                    % find matching number in PPR_keyword array:
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
                            numeric_field = [C{2}, C{3}, C{4}];
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
        function [rawData, parameters] = importB(obj, app) %#ok<*INUSL> 

            % Import path
            importPath = app.mrdImportPath;

            % Parameters
            info1 = jcampread(strcat(importPath,'acqp'));
            info2 = jcampread(strcat(importPath,'method'));

            % Scanner type
            parameters.scanner = 'B-type';

            % Slices
            parameters.NO_SLICES = str2num(info1.NSLICES);
            parameters.SLICE_THICKNESS = str2num(info2.pvm.slicethick) * parameters.NO_SLICES;

            % Matrix in readout direction
            parameters.NO_SAMPLES = info1.acq.size(1) / 2;
            if isfield(info2.pvm,"matrix")
                parameters.NO_VIEWS = info2.pvm.encmatrix(1);
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

            if isfield(info2.pvm,'spackarrreadorient')
                if strcmp(info2.pvm.spackarrreadorient,'L_R')
                    parameters.PHASE_ORIENTATION = 0;
                    flr =  0;
                    pm1 = +1;
                    pm2 = -1;
                end
                if strcmp(info2.pvm.spackarrreadorient,'A_P')
                    parameters.PHASE_ORIENTATION = 1;
                    flr =  0;
                    pm1 = -1;
                    pm2 = -1;
                end
                if strcmp(info2.pvm.spackarrreadorient,'H_F')
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
                if isfield(info2.pvm,'enczfaccel1') && isfield(info2.pvm,'encpftaccel1')
                    parameters.gp_var_mul = round(pm1 * info2.pvm.ppggradamparray1 * str2num(info2.pvm.enczfaccel1) * str2num(info2.pvm.encpftaccel1) * (parameters.NO_VIEWS / 2 - 0.5));
                else
                    parameters.gp_var_mul = round(pm1 * info2.pvm.ppggradamparray1 * (parameters.NO_VIEWS / 2 - 0.5));
                end
                parameters.pe1_order = 3;
            elseif isfield(info2.pvm,'encvalues1')
                if isfield(info2.pvm,'enczf') && isfield(info2.pvm,'encpft')
                    parameters.gp_var_mul = round(pm1 * info2.pvm.encvalues1 * info2.pvm.enczf(2) * info2.pvm.encpft(2) * (parameters.NO_VIEWS / 2 - 0.5));
                else
                    parameters.gp_var_mul = round(pm1 * info2.pvm.encvalues1 * (parameters.NO_VIEWS / 2 - 0.5));
                end
                parameters.pe1_order = 3;
            else
                % assume zigzag
                parameters.pe1_order = 2;
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
            kreal = dataRaw(1:2:end);
            kim = dataRaw(2:2:end);
            kspace = kreal + 1j*kim;

            % Read navigator
            if isfile(strcat(importPath,'fid.NavFid'))
                fileID = fopen(strcat(importPath,'fid.NavFid'));
            else
                fileID = fopen(strcat(importPath,'rawdata.job1'));
            end
            navdata = fread(fileID,datatype);
            fclose(fileID);
            kreal = navdata(1:2:end);
            kim = navdata(2:2:end);
            navkspace = kreal + 1j*kim;

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
                parameters.EXPERIMENT_ARRAY = floor(length(kspace)/singleRep);
                kspace = kspace(1:singleRep*parameters.EXPERIMENT_ARRAY,:);
                kspace = reshape(kspace,parameters.NO_SLICES,parameters.NO_SAMPLES,parameters.nr_coils,parameters.NO_VIEWS,parameters.EXPERIMENT_ARRAY);
                
                parameters.EXPERIMENT_ARRAY = size(kspace,5);
                kspace = permute(kspace,[3,5,1,4,2]); % nc, nr, ns, np, nf

                % Flip readout if needed
                if flr
                    kspace = flip(kspace,5);
                end

                % Coil intensity scaling
                if isfield(info2.pvm,'encchanscaling')
                    for i = 1:parameters.nr_coils
                        kspace(i,:) = kspace(i,:) * info2.pvm.encchanscaling(i);
                    end
                end

                % Navigator
                navkspace = navkspace(1:parameters.NO_SLICES*parameters.no_samples_nav*parameters.nr_coils*parameters.NO_VIEWS*parameters.EXPERIMENT_ARRAY);
                navkspace = reshape(navkspace,parameters.NO_SLICES,parameters.no_samples_nav,parameters.nr_coils,parameters.NO_VIEWS,parameters.EXPERIMENT_ARRAY);
                navkspace = permute(navkspace,[3,5,1,4,2]);

                % 34 point spacer
                kspacer = zeros(parameters.nr_coils,parameters.EXPERIMENT_ARRAY,parameters.NO_SLICES,parameters.NO_VIEWS,34);

                % Combine navigator + spacer + k-space
                raw = cat(5,navkspace,kspacer,kspace);
                rawData = cell(parameters.nr_coils);
                for i = 1:parameters.nr_coils
                    rawData{i} = squeeze(raw(i,:,:,:,:));
                    rawData{i} = reshape(rawData{i},parameters.EXPERIMENT_ARRAY,parameters.NO_SLICES,parameters.NO_VIEWS,parameters.NO_SAMPLES+34+parameters.no_samples_nav);
                end
     
            end

            % 3D data
            if strcmp(info2.pvm.spatdimenum,"3D") || strcmp(info2.pvm.spatdimenum,"<3D>")

                % 2nd phase encoding direction
                parameters.NO_VIEWS_2 = info1.acq.size(3);
                if isfield(info2.pvm,"matrix")
                    parameters.NO_VIEWS = info2.pvm.encmatrix(3);
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
                kspace = reshape(kspace,parameters.nr_coils,parameters.NO_SAMPLES,parameters.NO_VIEWS,parameters.NO_VIEWS_2,[]);
                parameters.EXPERIMENT_ARRAY = size(kspace,5);
                kspace = permute(kspace,[1,5,4,3,2]);

                % Flip readout if needed
                if flr
                    kspace = flip(kspace,5);
                end

                % Coil intesnity scaling
                if isfield(info2.pvm,'encchanscaling')
                    for i = 1:parameters.nr_coils
                        kspace(i,:) = kspace(i,:) * info2.pvm.encchanscaling(i);
                    end
                end

                % Navigator
                navkspace = reshape(navkspace,parameters.nr_coils,parameters.no_samples_nav,parameters.NO_VIEWS,parameters.NO_VIEWS_2,parameters.EXPERIMENT_ARRAY);
                navkspace = permute(navkspace,[1,5,4,3,2]);

                % 34 point spacer
                kspacer = zeros(parameters.nr_coils,parameters.EXPERIMENT_ARRAY,parameters.NO_VIEWS_2,parameters.NO_VIEWS,34);

                % Combine navigator + spacer + k-space
                raw = cat(5,navkspace,kspacer,kspace);
                for i = 1:parameters.nr_coils
                    rawData{i} = squeeze(raw(i,:,:,:,:));
                    rawData{i} = reshape(rawData{i},parameters.EXPERIMENT_ARRAY,parameters.NO_VIEWS_2,parameters.NO_VIEWS,parameters.NO_SAMPLES+34+parameters.no_samples_nav);
                end

            end


            % read reco files to a structure
            function struct = jcampread(filename) %#ok<STOUT> 

                % Open file read-only big-endian
                fid = fopen(filename,'r','b');
                skipline=0;

                % Loop through separate lines
                if fid~=-1
                    while 1
                        if skipline
                            line=nextline;
                            skipline=0;
                        else
                            line=fgetl(fid);
                        end
                        % Testing the text lines
                        while length(line)<2
                            line=fgetl(fid);
                        end
                        % Parameters and optional size of parameter are on lines starting with '##'
                        if line(1:2) == '##' %#ok<*BDSCA> 
                            % Parameter extracting and formatting
                            % Read parameter name
                            paramname = fliplr(strtok(fliplr(strtok(line,'=')),'#'));
                            % Check for illegal parameter names starting with '$' and correct (Matlab does not accepts variable names starting with $)
                            if paramname(1) == '$'
                                paramname = paramname(2:length(paramname));
                                % Check if EOF, if true return
                            elseif paramname(1:3) == 'END'
                                break
                            end
                            % Parameter value formatting
                            paramvalue = fliplr(strtok(fliplr(line),'='));

                            % Check if parameter values are in a matrix and read the next line
                            if paramvalue(1) == '('
                                paramvaluesize = str2num(fliplr(strtok(fliplr(strtok(paramvalue,')')),'(')));
                                % Create an empty matrix with size 'paramvaluesize' check if only one dimension
                                if ~isempty(paramvaluesize)
                                    if size(paramvaluesize,2) == 1
                                        paramvaluesize = [paramvaluesize,1];  
                                    end
                                    % Read the next line
                                    nextline = fgetl(fid);
                                    % See whether next line contains a character array
                                    if nextline(1) == '<'
                                        paramvalue = fliplr(strtok(fliplr(strtok(nextline,'>')),'<')); %#ok<*NASGU> 
                                    elseif strcmp(nextline(1),'L') || strcmp(nextline(1),'A') || strcmp(nextline(1),'H')
                                        paramvalue = nextline;
                                    else
                                        % Check if matrix has more then one dimension
                                        if paramvaluesize(2) ~= 1
                                            paramvaluelong = str2num(nextline);
                                            while (length(paramvaluelong)<(paramvaluesize(1)*paramvaluesize(2))) & (nextline(1:2) ~= '##') %#ok<*AND2> 
                                                nextline = fgetl(fid);
                                                paramvaluelong = [paramvaluelong str2num(nextline)];  
                                            end
                                            if (length(paramvaluelong)==(paramvaluesize(1)*paramvaluesize(2))) & (~isempty(paramvaluelong))
                                                paramvalue=reshape(paramvaluelong,paramvaluesize(1),paramvaluesize(2));
                                            else
                                                paramvalue=paramvaluelong;
                                            end
                                            if length(nextline)>1
                                                if (nextline(1:2) ~= '##')
                                                    skipline=1;
                                                end
                                            end
                                        else
                                            % If only 1 dimension just assign whole line to paramvalue
                                            paramvalue = str2num(nextline);
                                            if ~isempty(str2num(nextline))
                                                while length(paramvalue)<paramvaluesize(1)
                                                    line=fgetl(fid);
                                                    paramvalue = [paramvalue str2num(line)];  
                                                end
                                            end
                                        end
                                    end
                                else
                                    paramvalue='';
                                end
                            end

                            % Add paramvalue to structure.paramname
                            if isempty(findstr(paramname,'_'))
                                eval(['struct.' paramname '= paramvalue;']); %#ok<*EVLDOT>
                            else
                                try
                                    eval(['struct.' lower(paramname(1:findstr(paramname,'_')-1)) '.' lower(paramname(findstr(paramname,'_')+1:length(paramname))) '= paramvalue;']);
                                catch
                                    eval(['struct.' lower(paramname(1:findstr(paramname,'_')-1)) '.' datestr(str2num(paramname(findstr(paramname,'_')+1:findstr(paramname,'_')+2)),9) ...
                                        paramname(findstr(paramname,'_')+2:length(paramname)) '= paramvalue;']); %#ok<*FSTR>
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
                
                obj.validDataFlag = false;
                app.TextMessage(ME.message);
                
            end
            
        end
        
        
        
        
        % ---------------------------------------------------------------------------------
        % Write MRD file
        % ---------------------------------------------------------------------------------
        function obj = writeDataToMrd(obj, objKspace, filename, parameters)

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
            % footer - int8 data type footer as copied from an MRD containing a copy of
            % the PPR file, including preceeding 120-byte zeros
            % Output: 1 if write was successful, 0 if not, -1 if failed early (e.g. the dimension checks)

            % data: multidimensional, complex float, with dimensions arranged
            % dimensions: structure

            kSpaceMRDdata = objKspace.kSpaceMrd;
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
            
            if strcmp(obj.dataType,'3D') 
                kSpaceMRDdata = flip(permute(kSpaceMRDdata,[1,3,2,4,5,6,7]),1);
            else

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

            parameters = {':NO_SAMPLES no_samples, ',':NO_VIEWS no_views, ',':NO_VIEWS_2 no_views_2, ', ...
                ':NO_ECHOES no_echoes, ',':EXPERIMENT_ARRAY no_experiments, ',':NO_AVERAGES no_averages, ', ...
                ':VAR pe1_order, ',':VAR slice_nav, ',':VAR radial_on, ', ...
                ':VAR frame_loop_on, ',':VAR tr, ',':VAR te, ', ...
                ':BATCH_SLICES batch_slices, ',':NO_SLICES no_slices, ', ...
                ':VAR VFA_size, ',':VAR ti, ',':VAR pe2_centric_on, ', ...
                ':VAR tr_extra_us, '
                };

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

                txt = parameters{i};
                var = replacePars{i};
                pos = strfind(inputFooter,txt);

                if ~isempty(pos)
                    try
                        oldTextLength = strfind(inputFooter(pos+length(txt):pos+length(txt)+6),char(13))-1;
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
        function obj = readRprfile(obj, app, filename)
            
            try
                fid = fopen(filename,'r');
                obj.rprFile = char(fread(fid,Inf,'uchar')');
                fclose(fid);
                obj.rprFlag = true;
            catch
                obj.rprFile = '';
                obj.rprFlag = false;
                app.TextMessage('WARNING: rpr file not found ...');
            end
            
        end
        
        


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

            inputRpr = obj.rprFile;

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

            replacePars = {par.NoEchoes,par.NoEchoes, ...
                par.NoExperiments, par.NoExperiments, ...
                par.NoSamples, par.NoSamples, par.NoSamples, ...
                par.NoViews, par.NoViews, par.NoViews, ...
                par.NoViews2, par.NoViews2, par.NoViews2, ...
                par.NoSlices, par.NoSlices, ...
                par.NoSamples, par.NoViews, par.NoViews2, ...
                par.View1order, par.View2order
                };

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
        % Write physiological data to files
        % ---------------------------------------------------------------------------------
        function objData = writePhysLog(objData,objNav,exportPath)

            if ~isempty(objNav.heartTrigTime)
                writematrix(objNav.heartTrigTime',strcat(exportPath,filesep,'cardtrigpoints.txt'),'Delimiter','tab');
            end

            if ~isempty(objNav.respTrigTime)
                writematrix(objNav.respTrigTime',strcat(exportPath,filesep,'resptrigpoints.txt'),'Delimiter','tab');
            end

            if ~isempty(objNav.respWindowTime)
                writematrix(objNav.respWindowTime,strcat(exportPath,filesep,'respwindow.txt'),'Delimiter','tab');
            end

        end % writePhysLog



        % Export the reconstruction settings
        function objData = ExportRecoParametersFcn(objData,app)

            pars = strcat(...
                "------------------------- \n\n", ...
                "RETROSPECTIVE ", app.appVersion,"\n\n", ...
                "Gustav Strijkers\n", ...
                "Amsterdam UMC\n", ...
                "g.j.strijkers@amsterdamumc.nl\n\n", ...
                "------------------------- \n", ...
                "\nDATA \n\n", ...
                "file = ", app.DataFile.Value, "\n", ...
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
                "\nCINE\n\n", ...
                "type = ",app.RecoTypeDropDown.Value, "\n", ...
                "#frames = ",num2str(app.FramesEditField.Value), "\n", ...
                "#dynamics = ",num2str(app.DynamicsEditField.Value), "\n", ...
                "\nRECO \n\n", ...
                "WVxyz = ",num2str(app.WVxyzEditField.Value), "\n", ...
                "TVxyz = ",num2str(app.TVxyzEditField.Value), "\n", ...
                "LLRxyz = ",num2str(app.LLRxyzEditField.Value), "\n", ...
                "TVcine = ",num2str(app.TVcineEditField.Value), "\n", ...
                "TVdyn = ",num2str(app.TVdynEditField.Value), "\n", ...
                "ESPIRiT = ",num2str(app.ESPIRiTCheckBox.Value), "\n", ...
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

            if strcmp(app.retroDataPars.dataType,'3Dute')
                pars = strcat(pars,...
                    "\n3D UTE\n\n", ...
                    "Gx delay = ",num2str(app.GxDelayEditField.Value),"\n",...
                    "Gy delay = ",num2str(app.GyDelayEditField.Value),"\n",...
                    "Gz delay = ",num2str(app.GzDelayEditField.Value),"\n",...
                    "offset = ",num2str(app.DataOffsetRadialEditField.Value),"\n" ...
                    );
            end

            if strcmp(app.retroDataPars.dataType,'2Dradial')

                if app.HalfCircleButton.Value == 1              trajType = 1; end %#ok<SEPEX>
                if app.FullCircleButton.Value == 1              trajType = 2; end %#ok<SEPEX>
                if app.FullCircleInterleavedButton.Value == 1   trajType = 3; end %#ok<SEPEX>

                pars = strcat(pars,...
                    "\n2D radial\n\n", ...
                    "Gx delay = ",num2str(app.GxDelayEditField.Value),"\n",...
                    "Gy delay = ",num2str(app.GyDelayEditField.Value),"\n",...
                    "Gz delay = ",num2str(app.GzDelayEditField.Value),"\n",...
                    "trajectory = ",num2st(trajType)...
                    );
            end

            fid = fopen(strcat(objData.exportDir,filesep,'recoparmeters_',app.tag,'.txt'),'wt');
            fprintf(fid,pars);
            fclose(fid);

        end % ExportRecoParametersFcn





        
    end % methods
    
end % retroData

