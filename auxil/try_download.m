function success = try_download(filename, dest_path)
% try_download  This function is part of the
% MATLAB toolbox for accessing Argo float data.
%
% USAGE:
%   success = try_download(filename, dest_path)
%
% DESCRIPTION:
%   This function attempts to download a file from any of the GDACs
%   specified in the Settings.hosts cell array.
%
% INPUTS:
%   filename  : name of the file at the GDAC
%   dest_path : full (relative or absolute) path to the local file
%
% OUTPUT:
%   success   : 1 for successul download; 2 for unsuccessful download,
%               but the file exists locally already; 0 for failure
%
% AUTHORS:
%   H. Frenzel, J. Sharp, A. Fassbender (NOAA-PMEL), N. Buzby (UW)
%
% CITATION:
%   H. Frenzel, J. Sharp, A. Fassbender, N. Buzby, 2022. OneArgo-Mat:
%   A MATLAB toolbox for accessing and visualizing Argo data.
%   Zenodo. https://doi.org/10.5281/zenodo.6588041
%
% LICENSE: oneargo_mat_license.m
%
% DATE: JUNE 1, 2022  (Version 1.0.1)
% Modified Aug 19, 2023 by Paul Lerner: added "options" to empty Certificate Filename becuase without this an error is thrown (see comments)
global Settings;

success = 0; % default: failure
for h = 1:length(Settings.hosts);
    options=weboptions; options.CertificateFilename=(''); %@PL added this because otherwise I get "error setting certificate verify locations:" with a path to a certificate bundle that does not exist!
    try
        if Settings.verbose
            fprintf('Attempting download of %s\nto %s ... ', ...
                [Settings.hosts{h}, filename], dest_path);
        end
        websave(dest_path, [Settings.hosts{h}, filename],options);
        [~,~,ext] = fileparts(dest_path);
        % this will throw an exception if it is not a valid nc file:
        if strcmp(ext, '.nc')
            ncinfo(dest_path);
        end
        if Settings.verbose
            fprintf('success!\n');
        end
        success = 1;
        break
    catch
        if Settings.verbose
            fprintf('failure!\n');
        end
        % delete bogus file if it was created during failed download
        % attempt
        if exist([dest_path, '.html'], 'file')
            delete([dest_path, '.html']);
        end
        % delete broken netcdf file if it was downloaded
        if exist(dest_path, 'file')
            delete(dest_path);
        end
        if h == length(Settings.hosts)
            if exist(dest_path, 'file') ~= 2
                return
            else
                % caller may distinguish between return values 1 and 2
                success = 2;
            end
        end
    end
end