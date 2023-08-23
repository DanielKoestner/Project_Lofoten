% OneArgo_tutorial_WMOlist.m
%
% Example WMO list from Daniel Koestner. This tutorial will show how to
% download and extract data if you already have a list of WMO ID's. 
%
% GO-BGC Float Data Workshop, 2023
% Yui Takeshita
% MBARI

clear; close all; clc


% ============== User Inputs ================================

% example WMO list, provided by Daniel Koestner for the 2023 GO-BGC float
% data workshop. 
WMOlist = [6900799
6902545
6902547
6902548
6902549
6903549
6903550
6903553
6903568
6903570
6903577
6903578
6903590
7901028
]; 

vars2get = {'TEMP', 'PSAL', 'DOXY', 'PRES', 'CHLA', 'BBP700', ...
    'CHLA', 'DOWNWELLING_PAR', 'NITRATE'};
% for a full list of sensors/parameters that can be used here, use
% list_sensors();

% choose which QFs you want to extract
QF2use = [1 2 8]; % good and probably good data. Need to add 8 because for some Provor floats, salinity is interpolated onto DOXY so it gets a 8 flag. 

initialize_tf = false; % boolean to decide whether to initialzie or not. If it has been initialized recently, then make this false. 

csvname = 'NAtlantic_WMOlist_flts.csv'; % name of output csv file

% ============================================================
% main starts here

% initialize OneArgo toolbox. this may take a few minutes depending on internet speed
if(initialize_tf)
    fprintf('\n\nInitializing toolbox...\n\n');
    initialize_argo(); 
    fprintf('\nInitialization complete...\n');
end

% download multiple floats. For a single float, use download_float()
% instead
download_multi_floats(WMOlist);

% make plot of location
fprintf('\nPlotting float locations...\n');
show_trajectories(WMOlist,...
    'color','multiple', ... % this plots different floats in different colors
    'title', 'NA floats');

% ============ Extract data into a table and write csv ====================
fprintf('\nCompiling DOXY data into table...\n');
% create structure 'data' with the variables vars2get
data = load_float_data(WMOlist, vars2get); 

fprintf('\nApplying QC filter...\n');
% apply qc filter based on QF flags chosen by QF2use
data_qc = qc_filter(data, vars2get, QF2use, 'adjusted', 'nostrict');

fprintf('\nCalculating MLD...\n');
% calculate mixed layer depth, based on temp and density thresholds
data_qc = calc_mld(data_qc, 'calc_mld_temp', 1, 'calc_mld_dens',1);

show_sections(WMOlist(1), {'BBP700'},'isopyc',[27:0.1:28],...
    'mld', 0,'qc', [1 2],...   % tells the function to plot mixed layer depth using T
    'raw', 'no','time_label','s'); % tells the function to plot raw (unadjusted) data
caxis([0.0001 0.10]);set(gca,'ColorScale','log');


fprintf('\nWriting table...\n');
% write data into a table
dataTable = data2table(data_qc);
% remove lines where all variables are NaNs. 
iremove = true(size(dataTable.DOXY));
for i = 1:length(vars2get)
    iremove = iremove | isnan(dataTable.vars2get{i}); 
end
dataTable(iremove,:) = []; % delete all NaN rows

fprintf('\nWriting csv file...\n');
writetable(dataTable, csvname); 
fprintf('\nDone!\n');






































