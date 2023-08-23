% OneArgo_tutorial_underice.m
%
% Example WMO list from Daniel Koestner. This tutorial will show how to
% download and extract data if you already have a list of WMO ID's. 

% clean up workspace
clear; close all; clc; 

% ============== User Inputs ================================

% lat and lon bounds for North Atlantic. create example for polygon
lat_bound = [-50 -80];
lon_bound = [0 90];

% date range from 2020 to 2023
t1 = [2023, 1, 1]; 
t2 = [2023, 1, 31];

% list of variables to extract
vars2get = {'TEMP', 'PSAL', 'PRES', 'DOXY'};

% choose which QFs you want to extract
QF2use = [1 2 8]; % good and probably good data. Need to add 8 because for some Provor floats, salinity is interpolated onto DOXY so it gets a 8 flag. 

initialize_tf = false; % boolean to decide whether to initialzie or not. If it has been initialized recently, then make this false. 

csvname = 'underice_flts.csv'; 

% ============================================================
% main starts here

% initialize OneArgo toolbox. this may take a few minutes depending on internet speed
if(initialize_tf)
    fprintf('\n\nInitializing toolbox...\n\n');
    initialize_argo(); 
    fprintf('\nInitialize argo complete...\n');
end

% get float ID and profiles that have Oxygen, and is in Delayed mode
fprintf('\nSelecting Southern Ocean floats and profiles...\n\n');
[SO_flts, SO_profs] = select_profiles(lon_bound, lat_bound, t1, t2, ...
    'sensors', vars2get, ...
    'mode', 'AD', ...
    'outside', 'none'); % excludes profiles that are not in both space and time constraints. 
fprintf('\nFinished selecting DOXY floats and profiles...\n\n');

fprintf('\nStart plotting... \n\n');
show_trajectories(SO_flts, 'float_profs', SO_profs, ...
    'mark_estim', 'yes', 'interp_lonlat', 'yes', ...
    'legend', 'no', 'size', 50); 


% ============ Extract data into a table and write csv ====================
fprintf('\nCompiling DOXY data into table...\n');
% create structure 'data' with the variables vars2get
[data, mdata] = load_float_data(SO_flts, vars2get, SO_profs); 

fprintf('\nApplying QC filter...\n');
% apply qc filter based on QF flags chosen by QF2use
data_qc = qc_filter(data, vars2get, QF2use, 'adjusted', 'nostrict');

fprintf('\nCalculating MLD...\n');
% calculate mixed layer depth, based on temp and density thresholds
% data_qc = calc_mld(data_qc, 'calc_mld_temp', 1, 'calc_mld_dens',1);

fprintf('\nWriting table...\n');
% write data into a table
dataTable = data2table(data_qc);
% remove lines where all variables are NaNs. 
iremove = true(size(dataTable.DOXY));
for i = 1:length(vars2get)
    iremove = iremove & isnan(dataTable.(vars2get{i})); 
end
dataTable(iremove,:) = []; % delete all NaN rows

fprintf('\nWriting csv file...\n');
writetable(dataTable, csvname); 
fprintf('\nDone!\n');

























