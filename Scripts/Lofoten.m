%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lofoten.m %%
% Driver routine for the MATLAB toolbox for plotting BGC-ARGO data from the Lofoten Basin.

% Derived from examples in Tutorial.m from the OneArgo git repository.
% Aug 18, 2023
%
% Script for visualizing BGC-ARGO data in the Lofoten Basin, particularly near the Lofoten Vortex
%
% Author: Paul Lerner
%
% [Original AUTHORS of Tutorial.m:
%   J. Sharp, H. Frenzel, A. Fassbender (NOAA-PMEL), N. Buzby (UW)
%
% CITATION:
%   H. Frenzel, J. Sharp, A. Fassbender, N. Buzby, 2022. OneArgo-Mat:
%   A MATLAB toolbox for accessing and visualizing Argo data.
%   Zenodo. https://doi.org/10.5281/zenodo.6588042
%
% LICENSE: oneargo_mat_license.m
%
% DATE: JUNE 1, 2022  (Version 1.0.1)]

%% Close figures, clean up workspace, clear command window
close all; clear; clc

% These need to be changed to whatever the path to OneArgo and Copernicus data is
% load dataset
oneargo = '/discover/nobackup/plerner/GO-BGC-ARGO/BGC_ARGO_toolbox/OneArgo-Mat';
addpath(oneargo); % just in case, add to MATLAB search path

% path to sea level/absolute dynamic topography data
path_SLA = '/discover/nobackup/plerner/GO-BGC-ARGO/BGC_ARGO_toolbox/COPERNICUS_REMOTESENS/cmems_obs-sl_glo_phy-ssh_my_allsat-l4-duacs-0.25deg_P1D_1692531538213.nc'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function defines standard settings and paths and creates Index
% and Profiles folders in your current path. It also downloads the Sprof
% index file from the GDAC to your Index folder. The Sprof index is
% referenced when downloading and subsetting float data based on user
% specified criteria in other functions.

initialize_argo();
do_pause();

%Float_list = [7901028,6903590,6902548,6903578,6903577,6903568,6903553,6903550,6903570,6903549,6902549,6902547,6902545,6900799];
%for i=1:length(Float_list);
%  WMO = Float_list(i);
%  success = download_float(WMO);
%end

%% Examine global structures
% These global structures contain a variety of useful variables for
% downloading and manipulating float data. 'Sprof' contains fields
% with information for each profile, 'Float' contains fields with
% information for each float, 'Settings' contains settings to be used in
% the backgroud during plotting, etc. Variables in the global structures
% can be altered within the initialize_argo.m file.
global Sprof Float Settings;

% create a list of Float ids from the floats in "Profiles Directory"

%Float_dir = dir('Profiles/*.nc');

%for i=1:length(Float_dir);
%  Float_list(i) = str2num(Float_dir(i).name(1:end-9));
%end
%Float_list = Float_list';


%% set lat/lon limits
latlim = [68 72];
lonlim = [-5 10];
t1=[2010 1 1];
t2=[2023 12 31];


maplimits.lat_lim = latlim;
maplimits.lon_lim = lonlim;


%% downloads floats
[Float_list,Float_profs] = select_profiles(lonlim,latlim,t1,t2,...
    'sensor','BBP700',... % this selects only floats with nitrate sensors
    'outside','none'); % All floats that cross into the time/space limits


%% subselect floats within lat/lon limts %%%%%%%%%%%%%%%%%%%%%%%5

Settings.mapping = 'plain' % for adt data
% Domain of interest





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% load Copernicus SLA data

SLA_time = ncread(path_SLA,'adt'); % daily resolved absolute dyanmic topography
SLA_mean = nanmean(SLA_time,3); % time-average adt
lat_sla = double(ncread(path_SLA,'latitude'));
lon_sla = double(ncread(path_SLA,'longitude'));
time_sla = double(ncread(path_SLA,'time'));
dt_sla = datetime(1950,1,1)+calyears(fix(time_sla./365))+caldays(floor(mod(time_sla,365)))+hours((mod(time_sla,365)-floor(mod(time_sla,365)))*24); % converts time in netcdf to datetime array

% create structure GridProd to pass to showtrajectories

GridProd.var = SLA_mean;
GridProd.lat = lat_sla;
GridProd.lon = lon_sla;
GridProd.time = time_sla;


%show_trajectories(Float_list(1),'float_profs',Float_profs(1),...
%    'color','multiple','maplimits',maplimits,'showgridded',GridProd); % this plots different floats in different colors

show_trajectories(Float_list,'float_profs',Float_profs,...
    'color','multiple','showgridded',GridProd); % this plots different floats in different colors

saveas(figure(1),'Float_trajectory.png');

%% load selected data from profiles
vars2get = {'TEMP','PSAL','PRES','CHLA','DOXY','BBP700'};
vars2get_bbp = {'PRES','BBP700'};
vars2get_chla = {'TEMP','PSAL','PRES','CHLA','DOXY'};


fprintf('\nCompiling BBP700 data...\n');
% create structure 'data' with the variables vars2get
data = load_float_data(Float_list, vars2get, Float_profs);

% qc filter the data
qc_flags = [1 2 8];
Data_qc = qc_filter(data, vars2get, qc_flags)
Data_qc_BBP = qc_filter(data, vars2get_bbp, qc_flags,'raw','y');
Data_qc_Chl = qc_filter(data, vars2get_chla, qc_flags);

% save data table 

fprintf('\n saving data without NaNs...\n');

dataTable_BBP = data2table(Data_qc);
% remove lines where DOXY, TEMP, PSAL, or PRES are NaNs. 
iremove = isnan(dataTable_BBP.BBP700_ADJUSTED) | isnan(dataTable_BBP.PRES_ADJUSTED) | isnan(dataTable_BBP.CHLA_ADJUSTED);
dataTable_BBP(iremove,:) = [];

save([oneargo,'/matfiles/BBP_CHL_BGC-ARGO.mat'],'Data_qc','dataTable_BBP','Data_qc_BBP','Data_qc_Chl');

% get lat, lon, and date of each float

[lon, lat, date] = get_lon_lat_time(Float_list, Float_profs);

formatOut = 'yyyy-mm-dd'; % for some reason, mm is month for datestr but minute for datetime
formatOut_dt = 'yyyy-MM-dd'


%% show profiles by month, with number of proifles corresponding to the number of years of data for each month

% show profiles and sections wihtin domain of interest 
for i = length(Float_list)-7;%1:length(Float_list);

  %show_profiles(Float_list(i), {'PSAL'}, 'per_float', 0, ...
  %  'float_profs', Float_profs(i), ...
  %  'obs', 'on', ... % plot a marker at each observation
  %  'title_add', ' (LF basin)',...  % add this to the title
  %  'qc', [1 2]); % apply QC flags

  %saveas(figure(1),['PSAL_',num2str(Float_list(i)),'.png']);
  %close all
  %show_profiles(Float_list(i), {'TEMP'}, 'per_float', 0, ...
  %  'float_profs', Float_profs(i), ...
  %  'obs', 'on', ... % plot a marker at each observation
  %  'title_add', ' (LF basin)',...  % add this to the title
  %  'qc', [1 2]); % apply QC flags
  %saveas(figure(1),['TEMP_',num2str(Float_list(i)),'.png']);
  %close all

  show_profiles(Float_list(i), {'BBP700'}, 'per_float', 0, ...
    'float_profs', Float_profs(i), ...
    'obs', 'on', ... % plot a marker at each observation
    'title_add', ' (LF basin)',...  % add this to the title
    'qc', [1 2 8]); % apply QC flags
  saveas(figure(1),[oneargo,'/Figures/BBP700_',num2str(Float_list(i)),'.png']);
  close all

end

close all



% show trajectories of chosen floats


Settings.pad_lat = 3;
Settings.pad_lon = 3;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make videos of selected floats within domain of interest, drawing rectange around approximate location of LV
for i = 1:length(Float_list);
  dtimes = datetime(date{i},'ConvertFrom','datenum','Format',formatOut_dt);
  fprofs = Float_profs{i};
  lonsels = lon{i};
  latsels = lat{i};
  writerObj = VideoWriter([oneargo,'/Figures/',num2str(Float_list(i)),'_traj']);
  open(writerObj) ;
  for j=1:length(dtimes);
    dtime = dtimes(j);
    indt_sla = find(abs(dt_sla -dtime) == min(abs(dt_sla -dtime)));
    GridProd.var = SLA_time(:,:,indt_sla);
    fprof{1} = fprofs(j);
    lonsel = lonsels(j);
    latsel = latsels(j);      
    txt = datestr(dtime,formatOut)
    show_trajectories(Float_list(i),'float_profs',fprof,'size',60,'maplimits',maplimits,'txt',txt,'figshow','off','showgridded',GridProd); 
    hold on;

    F = getframe(gcf);         
    writeVideo(writerObj,F) 
    hold off;

    close all
  end
  close(writerObj);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% show sections

for i = length(Float_list);

% either seasonal or daily frequency for x-axis
if Float_list(i) == 7901028; % temporary solution, need to automate
 tlabel = 'm';
else;
 tlabel = 's';
end

show_sections(Float_list(i), {'PSAL'},'float_profs',Float_profs(i),'isopyc',[27:0.1:28],...
    'mld', 0,'qc', [1 2],...   % tells the function to plot mixed layer depth using T
    'raw', 'no','time_label',tlabel); % tells the function to plot raw (unadjusted) data
caxis([34.9 35.25])
saveas(figure(1),[oneargo,'/Figures/PSAL_',num2str(Float_list(i)),'_section.png']);
close all
show_sections(Float_list(i), {'TEMP'},'float_profs',Float_profs(i),'isopyc',[27:0.1:28],...
    'mld', 0,'qc', [1 2],...   % tells the function to plot mixed layer depth using T
    'raw', 'no','time_label',tlabel); % tells the function to plot raw (unadjusted) data
saveas(figure(1),[oneargo,'/Figures/TEMP',num2str(Float_list(i)),'_section.png']);
close all

show_sections(Float_list(i), {'BBP700'},'float_profs',Float_profs(i),'isopyc',[27:0.1:28],...
    'mld', 0,'qc', [1 2 8],...   % tells the function to plot mixed layer depth using T
    'raw', 'no','time_label',tlabel); % tells the function to plot raw (unadjusted) data
caxis([0.0001 0.10]);set(gca,'ColorScale','log');
saveas(figure(1),[oneargo,'/Figures/BBP700',num2str(Float_list(i)),'_section.png']);
close all
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clean up the workspace
clear data mdata S success WMO ans
clc; close all

