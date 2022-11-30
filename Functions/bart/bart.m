% Call BART command from Matlab.

function [varargout] = bart(app, cmd, varargin)

%%%%% WINDOWS %%%%%%
if ispc

    if nargin==1 || all(cmd==0)
        app.TextMessage('Usage: bart <command> <arguments...>');
        return
    end

    bart_path = getenv('TOOLBOX_PATH'); %#ok<NASGU> 

    % clear the LD_LIBRARY_PATH environment variable (to work around
    % a bug in Matlab).

    setenv('LD_LIBRARY_PATH', '');
    name = strrep(tempname,' ','_');   % Windows user names with spaces give problems, replace with underscore

    in = cell(1, nargin-2);

    for i = 1:nargin-2
        if ischar(varargin{i})
            in{i} = varargin{i};
        else
            in{i} = strcat(name, 'in', num2str(i));
            writecfl(in{i}, varargin{i});
        end
    end

    in_str = sprintf(' %s', in{:});

    out = cell(1, nargout);

    for i=1:nargout
        out{i} = strcat(name, 'out', num2str(i));
    end

    out_str = sprintf(' %s', out{:});

    % For WSL and modify paths
    cmdWSL = WSLPathCorrection(cmd); %#ok<NASGU> 
    in_strWSL = WSLPathCorrection(in_str);
    out_strWSL =  WSLPathCorrection(out_str);

    % Execute the BART command
    [ERR,cmdout] = system(['wsl bart ',cmd,in_strWSL,out_strWSL]);

    % Return Version
    if contains(cmd,"version")
        app.TextMessage(cmdout);
        app.bartVersion = cmdout;
    end
   
    for i=1:nargin - 2
        if (exist(strcat(in{i}, '.cfl'),'file'))
            delete(strcat(in{i}, '.cfl'));
        end

        if (exist(strcat(in{i}, '.hdr'),'file'))
            delete(strcat(in{i}, '.hdr'));
        end
    end

    for i=1:nargout
        if ERR==0
            if contains(cmd,"estdelay") || contains(cmd,"-Rh") || contains(cmd,"version") 
                varargout{1} = cmdout;
            else
                varargout{i} = readcfl(out{i});
            end
        end
        if (exist(strcat(out{i}, '.cfl'),'file'))
            delete(strcat(out{i}, '.cfl'));
        end
        if (exist(strcat(out{i}, '.hdr'),'file'))
            delete(strcat(out{i}, '.hdr'));
        end
    end

    if ERR~=0
        app.bartVersion = 'none';
    end

end




%%%%% OSX / LINUX %%%%%%

if ismac

    if nargin==1 || all(cmd==0)
        app.TextMessage('Usage: bart <command> <arguments...>');
        return
    end

    bart_path = getenv('TOOLBOX_PATH');

    if isempty(bart_path)
        if exist('/usr/local/bin/bart', 'file')
            bart_path = '/usr/local/bin';
        elseif exist('/usr/bin/bart', 'file')
            bart_path = '/usr/bin';
        else
            app.TextMessage('Environment variable TOOLBOX_PATH is not set.');
        end
    end

    % clear the LD_LIBRARY_PATH environment variable (to work around
    % a bug in Matlab).

    if ismac==1
        setenv('DYLD_LIBRARY_PATH', '');
    else
        setenv('LD_LIBRARY_PATH', '');
    end

    name = tempname;

    in = cell(1, nargin - 2);

    for i = 1:nargin-2
        if ischar(varargin{i})
            in{i} = varargin{i};
        else
            in{i} = strcat(name, 'in', num2str(i));
            writecfl(in{i}, varargin{i});
        end
    end

    in_str = sprintf(' %s', in{:});

    out = cell(1, nargout);

    for i=1:nargout
        out{i} = strcat(name, 'out', num2str(i));
    end

    out_str = sprintf(' %s', out{:});

    % Execute the BART command
    [ERR,cmdout] = system([bart_path, '/bart ', cmd, ' ', in_str, ' ', out_str]);
    
    % Return Version
    if contains(cmd,"version")
        app.TextMessage(cmdout);
        app.bartVersion = cmdout;
    end

    for i=1:nargin - 2
        if (exist(strcat(in{i}, '.cfl'),'file'))
            delete(strcat(in{i}, '.cfl'));
        end

        if (exist(strcat(in{i}, '.hdr'),'file'))
            delete(strcat(in{i}, '.hdr'));
        end
    end

    % Output
    for i=1:nargout
        if ERR==0
            if contains(cmd,"estdelay") || contains(cmd,"-Rh") || contains(cmd,"version")
                varargout{1} = cmdout;
            else
                varargout{i} = readcfl(out{i}); %#ok<*AGROW> 
            end
        end
        if (exist(strcat(out{i}, '.cfl'),'file'))
            delete(strcat(out{i}, '.cfl'));
        end
        if (exist(strcat(out{i}, '.hdr'),'file'))
            delete(strcat(out{i}, '.hdr'));
        end
    end

    if ERR~=0
        app.bartVersion = 'none';
    end

end


