classdef NetCDFFile < NetCDFGroup
    % A class for reading and writing to NetCDF files
    %
    % NetCDF files are a standard file format for reading and writing data.
    % This class is designed to simplify the task of adding new dimensions,
    % variables, and attributes to a NetCDF file compared to using the
    % built-in `ncread` and `ncwrite` functions.
    %
    % ```matlab
    % ncfile = NetCDFFile('myfile.nc')
    %
    % % create two new dimensions and add them to the file
    % x = linspace(0,10,11);
    % y = linspace(-10,0,11);
    % ncfile.addDimension('x',x);
    % ncfile.addDimension('y',y);
    %
    % % Create new multi-dimensional variables, and add those to the file
    % [X,Y] = ncgrid(x,y);
    % ncfile.addVariable(X,{'x','y'});
    % ncfile.addVariable(Y,{'x','y'});
    % ```
    %
    % - Topic: Initializing
    % - Topic: Accessing file properties
    % - Topic: Working with dimensions
    % - Topic: Working with variables
    % - Topic: Working with groups
    % - Topic: Working with global attributes
    %
    % - Declaration: classdef NetCDFFile < NetCDFGroup
    properties
        % file path the NetCDF file
        % - Topic: Accessing file properties
        path

        % format
        % - Topic: Accessing file properties
        format = 'FORMAT_NETCDF4'
    end

    properties (Dependent)
        filename
    end
    methods
        function self = NetCDFFile(path,options)
            % NetCDFFile initialize an from existing or create new file
            %
            % Calling,
            %   ncfile = NetCDFFile(path)
            % will load an existing file (if one exists) or create a new
            % file (if none exists).
            %
            %   ncfile = NetCDFFile(path,shouldOverwriteExisting=1)
            % will delete any existing file and create a new file.
            %
            % - Topic: Initializing
            % - Declaration: ncfile = NetCDFFile(path,options)
            % - Parameter path: path to write file
            % - Parameter shouldOverwriteExisting: (optional) boolean indicating whether or not to overwrite an existing file at the path. Default 0.
            % - Returns: a new NetCDFFile instance
            arguments
                path char {mustBeNonempty}
                options.shouldOverwriteExisting logical = false
                options.shouldUseClassicNetCDF logical = false
                options.shouldReadOnly logical = false
            end

            if isfile(path) && options.shouldOverwriteExisting == 1
                delete(path);
                shouldCreateNewFile = 1;
            elseif ~isfile(path)
                shouldCreateNewFile = 1;
            else
                shouldCreateNewFile = 0;
            end

            if shouldCreateNewFile == 1
                if options.shouldUseClassicNetCDF == 1
                    ncid = netcdf.create(path, bitor(netcdf.getConstant('SHARE'),netcdf.getConstant('WRITE')));
                else  
                    ncid = netcdf.create(path, netcdf.getConstant('NETCDF4'));
                end
            else
                if options.shouldReadOnly
                    ncid = netcdf.open(path);
                else
                    ncid = netcdf.open(path, bitor(netcdf.getConstant('SHARE'),netcdf.getConstant('WRITE')));
                end
            end
            
            self@NetCDFGroup(id=ncid);
            self.path = path;
            self.format = netcdf.inqFormat(self.id);            
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Utilities
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function filename = get.filename(self)
            [~,name,ext] = fileparts(self.path);
            filename = name + ext;
        end


        function self = sync(self)
            % - Topic: Accessing file properties
            netcdf.sync(self.id);
        end

        % function self = open(self)
        %     % - Topic: Accessing file properties
        %     self.ncid = netcdf.open(self.path, bitor(netcdf.getConstant('SHARE'),netcdf.getConstant('WRITE')));
        %     self.format = netcdf.inqFormat(self.ncid);
        % end

        function self = close(self)
            % - Topic: Accessing file properties
            netcdf.close(self.id);
            self.id = [];
        end

        function delete(self)
            if ~isempty(self.id)
                netcdf.close(self.id);
                self.id = [];
                disp('NetCDFFile closed.');
            end
        end

        function ncfile = duplicate(self,path,options)
            arguments
                self 
                path 
                options.shouldOverwriteExisting logical = false
                options.indexRange = dictionary(string.empty,cell.empty)
            end
            if options.shouldOverwriteExisting == 1
                if isfile(path)
                    delete(path);
                end
            else
                if isfile(path)
                    error('A file already exists with that name.')
                end
            end
            ncfile = NetCDFFile(path);
            ncfile.addDuplicateGroup(self,shouldAddToSelf=true,indexRange=options.indexRange);
        end
    end
end