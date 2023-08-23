% OneArgo_tutorial_NAtlantic_DOXY.m
%
%
% GO-BGC Float Data Workshop, 2023
%
% Created by: Yui Takeshita
% MBARI

% clean up workspace
clear; close all; clc; 

% ============== User Inputs ================================

% lat and lon bounds for North Atlantic. create example for polygon
lat_NA = [60,  30,  30, 60, 70];
lon_NA = [-70, -86, 0,  0,  -40];

% date range from 2020 to 2023
t1 = [2020, 1, 1]; 
t1_core = [2023, 1, 1]; % shorter time period for core to speed things up for tutorial
t2 = [2023, 12, 31];


% choose which QFs you want to extract
QF2use = [1 2 8]; % good and probably good data. Need to add 8 because for some Provor floats, salinity is interpolated onto DOXY so it gets a 8 flag. 

initialize_tf = false; % boolean to decide whether to initialzie or not. If it has been initialized recently, then make this false. 

% ============================================================
% main starts here

% initialize OneArgo toolbox. this may take a few minutes depending on internet speed
if(initialize_tf)
    initialize_argo(); 
    fprintf('\nInitialize argo complete...\n');
end

% get float ID and profiles that have Oxygen, and is in Delayed mode
fprintf('\nSelecting DOXY floats and profiles...\n\n');
[NA_DOXY_flts, NA_DOXY_profs] = select_profiles(lon_NA, lat_NA, t1, t2, ...
    'Sensor', 'DOXY', ...
    'mode', 'D', ... % can put 'RAD' for real-time, Adjusted, and Delayed mode (all data)
    'outside', 'none'); % excludes profiles that are not in both space and time constraints. 
fprintf('\nFinished selecting DOXY floats and profiles...\n\n');

% get core float in the same time period for mapping. This will take a
% while because 
fprintf('\nSelecting Core floats and profiles...\n\n');
[NA_core_flts, NA_core_profs] = select_profiles(lon_NA, lat_NA, t1_core, t2, ...
    'type', 'phys', ...
    'mode', 'D');
fprintf('\nFinished selecting Core floats and profiles...\n\n');



%% make plot of location
fprintf('\nPlotting all DOXY float locations...\n');
show_trajectories(NA_DOXY_flts,...
   'color','multiple', ... % this plots different floats in different colors
   'title', 'NA DOXY floats');

fprintf('\nPlotting last DOXY float locations...\n');
show_trajectories(NA_DOXY_flts,...
    'color','multiple', ... % this plots different floats in different colors
    'title', 'NA DOXY floats', 'position', 'last', 'size', 75);

fprintf('\nPlotting last Core float locations...\n');
%SLOWshow_trajectories(NA_core_flts, ...
%    'color', 'multiple', ...
%    'title', 'NA Core floats');
show_trajectories(NA_core_flts, ...
    'color', 'multiple', ...
    'title', 'NA Core floats', 'position', 'last', 'size', 75);

fprintf('\nPaused... hit any key to continue\n');
pause(); % pause to look at plots together

% extract data and put into a table and csv (slow; so will need to uncomment)
%% ================ DOXY Floats ====================
fprintf('\nCompiling DOXY data into table...\n');
[data_DOXY, meta_DOXY] = load_float_data(NA_DOXY_flts, {'TEMP', 'PSAL', 'DOXY', 'PRES'});
data_DOXYqc = qc_filter(data_DOXY, {'TEMP', 'PSAL', 'DOXY', 'PRES'}, QF2use, 'raw', 'no_strict');

data_DOXYqc = calc_mld(data_DOXYqc, 'calc_mld_temp', 1, 'calc_mld_dens',1);
dataTable_DOXY = data2table(data_DOXYqc);
% remove lines where DOXY, TEMP, PSAL, or PRES are NaNs. 
iremove = isnan(dataTable_DOXY.DOXY_ADJUSTED) | isnan(dataTable_DOXY.PRES_ADJUSTED) | ...
    isnan(dataTable_DOXY.PSAL_ADJUSTED) | isnan(dataTable_DOXY.TEMP_ADJUSTED);
dataTable_DOXY(iremove,:) = [];
%SLOW writetable(dataTable_DOXY, 'NAtlantic_DOXYdata.csv');
fprintf('\nDone!\n');

%% ================= Core Floats ===================
fprintf('\nCompiling Core data into table...\n');
data_core = load_float_data(NA_core_flts, {'TEMP', 'PSAL', 'PRES'});
data_coreqc = qc_filter(data_core, {'TEMP', 'PSAL', 'PRES'}, QF2use, 'raw', 'no_strict');
% Same error as above 
data_coreqc = calc_mld(data_coreqc, 'calc_mld_temp', 1, 'calc_mld_dens',1);
dataTable_core = data2table(data_coreqc);
iremove = isnan(dataTable_core.TEMP_ADJUSTED) | ...
    isnan(dataTable_core.PSAL_ADJUSTED) | ...
    isnan(dataTable_core.PRES_ADJUSTED);
dataTable_core(iremove,:) = [];
%SLOW writetable(dataTable_core, 'NAtlantic_coredata.csv');
fprintf('\nDone!\n');





