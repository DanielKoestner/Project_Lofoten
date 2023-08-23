%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to Process BGC-ARGO BBP700 and chla from a .mat file
% expects two Data structures (containing profiles, with each profiles contianing paramter names with 
% dimensions depth (presure) x cycle number):

% i.  a data structure that includes chlorophyll
% ii. a data structure that includes BBP700

% Currently, the script assumes that the pressure levels of chl and bbp are the same,
% and uses "PRES_ADJUSTED" from the chl data structure as the pressure variable.

% Follows Briggs et al., 2020, Science methodology
% Paul Lerner, Aug 22, 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

clear all
close all
clc

% load dataset
oneargo = '/discover/nobackup/plerner/GO-BGC-ARGO/BGC_ARGO_toolbox/OneArgo-Mat';
addpath(oneargo); % just in case, add to MATLAB search path

filename = 'BBP_CHL_BGC-ARGO';
load([oneargo,'/matfiles/',filename,'.mat']);


% Data structure names
Data_chl = 'Data_qc_Chl';
Data_bbp = 'Data_qc_BBP';


% percentile of bbp_bsr (small + "refractory" particles -> Briggs et al., 2020) between 850-900 m below which bbp is considerered to all be "refractory"
pctile_bsr = 0.05; %
pctile_chl = 0.05; 

% threshold above which to remove bp data 

thr = 0.01;

% loop over floats

fieldnms = fieldnames(Data_qc_BBP);
for i=1:length(fieldnames(Data_qc_BBP));
  fieldnm = fieldnms(i);
  fprintf(['\n Starting Float',char(fieldnm),'...\n']);
 
  Float_Struc_BBP = eval([Data_bbp,'.',char(fieldnm)]);
  Float_Struc_Chl = eval([Data_chl,'.',char(fieldnm)]);

% add to structure arrays for BBP_bs, BBP_bl

  Float_Struc_BBP.BBP700_bl = Float_Struc_BBP.BBP700;
  Float_Struc_BBP.BBP700_bs = Float_Struc_BBP.BBP700;

  Float_Struc_Chl.BBP700_bl = Float_Struc_Chl.CHLA_ADJUSTED;
  Float_Struc_Chl.chl_bs = Float_Struc_Chl.CHLA_ADJUSTED;
% loop over profiles, first to get instrument noise for each float

  BBP_resid_bin_tot = 0; % total vector of all binned BBP_resid profiles, needed to compute instrucment noise
  Pres_BBP_tmp_bin_tot = 0;

  chl_resid_bin_tot = 0; % total vector of all binned chl profiles, needed to compute instrucment noise
  Pres_chl_tmp_bin_tot = 0;


  Float_Struc_BBP.BBP700_despike = NaN.*ones(size(Float_Struc_BBP.BBP700));
  Float_Struc_Chl.CHLA_ADJUSTED_despike = NaN.*ones(size(Float_Struc_Chl.CHLA_ADJUSTED));

%%%%%% First for loop over is used to obtained background noise for the entire float (same noise applied to every prifle in the second loop)  %%%%%%%%%%

  for j=1:length(Float_Struc_Chl.PRES_ADJUSTED(1,:));


    
% filter BBP770

    BBP_tot_tmp = Float_Struc_BBP.BBP700(:,j);
    Pres_tmp = Float_Struc_Chl.PRES_ADJUSTED(:,j);
    chl_tot_tmp = Float_Struc_Chl.CHLA_ADJUSTED(:,j);
  
% zero out negative chl values

    chl_tot_tmp(find(chl_tot_tmp<0)) = 0;

% remove within 25 m of BBP_tot > thr 
    ind_hi = find(BBP_tot_tmp>thr);
    for k = 1:length(BBP_tot_tmp);
      if any(Pres_tmp(k) >Pres_tmp(ind_hi) -25) && any(Pres_tmp(k) <Pres_tmp(ind_hi) +25);
         BBP_tot_tmp(k) = NaN;
         chl_tot_tmp(k) = NaN;
         Pres_tmp(k) = NaN;
      end
    end

     
% remove NaN values
    indnonan = find(~isnan(BBP_tot_tmp) & ~isnan(chl_tot_tmp));
    
    Pres_tmp = Pres_tmp(indnonan);
    BBP_tot_tmp = BBP_tot_tmp(indnonan);
    chl_tot_tmp = chl_tot_tmp(indnonan);
   
 
% break out of loop if all NaNs
    if sum(~isnan(Pres_tmp)) == 0;
       continue
    end


% create arrays of BBP/chl after the initial despiking. These are stricly to compare with profiles after processing, to ensure only a small loss of data
% Note that despiking sets BBP of values in vcinity of thr to NaN, AS WELL AS Chl and Pressure

    Float_Struc_BBP.BBP700_despike(indnonan,j) = BBP_tot_tmp;
    Float_Struc_Chl.CHLA_ADJUSTED_despike(indnonan,j) = chl_tot_tmp;

% filter BBP.chl profiles    
    BBP_fil_tmp = movmin(BBP_tot_tmp,11); % 11 pt moving minimum filter
    BBP_fil_tmp = movmax(BBP_fil_tmp,11); % 11 pt moving minimum filter

    chl_fil_tmp = movmin(chl_tot_tmp,11); % 11 pt moving minimum filter
    chl_fil_tmp = movmax(chl_fil_tmp,11); % 11 pt moving minimum filter
% isolate BBP700 residual

    BBP_resid = BBP_tot_tmp - BBP_fil_tmp; % this is large particels + refractory particles, or residual spikes in Briggs et al., 2020
    chl_resid = chl_tot_tmp - chl_fil_tmp; % this is large particels + refractory particles, or residual spikes in Briggs et al., 2020

    if max(Pres_tmp) > 300;
      edges = (300:50:max(Pres_tmp));
      [~,~,loc]=histcounts(Pres_tmp,edges); % number of depths per bin
      Pres_b300 = Pres_tmp;
      BBP_resid_b300 = BBP_resid;
      BBP_resid_b300(find(loc==0)) = []; % remove any bins where there are no pressure values
      Pres_b300(find(loc==0)) = []; % remove any bins where there are no pressure values
      loc(find(loc==0)) = [];% remove any bins where there are no pressure values
   
      Pres_BBP_tmp_bin = accumarray(loc(:),Pres_b300(:),[],@mean);
      BBP_resid_bin = accumarray(loc(:),BBP_resid_b300(:),[],@mean);% depth-binned mean
      BBP_resid_bin_tot = vertcat(BBP_resid_bin_tot,BBP_resid_bin(:));
      Pres_BBP_tmp_bin_tot = vertcat(Pres_BBP_tmp_bin_tot,Pres_BBP_tmp_bin(:));

      [~,~,loc]=histcounts(Pres_tmp,edges); % number of depths per bin
      Pres_b300 = Pres_tmp;
      chl_resid_b300 = chl_resid;
      chl_resid_b300(find(loc==0)) = []; % remove any bins where there are no pressure values
      Pres_b300(find(loc==0)) = []; % remove any bins where there are no pressure values
      loc(find(loc==0)) = [];% remove any bins where there are no pressure values

      Pres_chl_tmp_bin = accumarray(loc(:),Pres_b300(:),[],@mean);
      chl_resid_bin = accumarray(loc(:),chl_resid_b300(:),[],@mean);% depth-binned mean
      chl_resid_bin_tot = vertcat(chl_resid_bin_tot,chl_resid_bin(:));
      Pres_chl_tmp_bin_tot = vertcat(Pres_chl_tmp_bin_tot,Pres_chl_tmp_bin(:));
    end

  end

% calcualte median of each composite vertical profile

  Pres_BBP_tmp_bin_tot = Pres_BBP_tmp_bin_tot(2:end);
  BBP_resid_bin_tot = BBP_resid_bin_tot(2:end);
  chl_resid_bin_tot = chl_resid_bin_tot(2:end);

  median_resid_bbp = median(BBP_resid_bin_tot);
  median_resid_chl = median(chl_resid_bin_tot);

% remove all bins 2x> than the median

  Pres_BBP_tmp_bin_tot = Pres_BBP_tmp_bin_tot(find(BBP_resid_bin_tot<= 2*median_resid_bbp));
  BBP_resid_bin_tot = BBP_resid_bin_tot(find(BBP_resid_bin_tot<= 2*median_resid_bbp));
  Pres_chl_tmp_bin_tot = Pres_chl_tmp_bin_tot(find(chl_resid_bin_tot<= 2*median_resid_chl));
  chl_resid_bin_tot = chl_resid_bin_tot(find(chl_resid_bin_tot<= 2*median_resid_chl));

  Float_Struc_BBP.BBP700_noise = median(BBP_resid_bin_tot); % Instrument noise for this float
  Float_Struc_Chl.CHLA_noise = median(chl_resid_bin_tot); % Instrument noise for this float

% now we filter profiles again and add the instrument noise
  BBP_bsr_tot = 0;  % an array of BBP for each float, used to find the 5th percentile of BBP_bsr between 850-900 m;
  chl_bsr_tot = 0;  % an array of chl for each float, used to find the 5th percentile of chl_bsr between 850-900 m;

% preallocated arrays in structures

  Float_Struc_BBP.BBP700_bl = NaN.*ones(size(Float_Struc_BBP.BBP700));
  Float_Struc_BBP.BBP700_bs = NaN.*ones(size(Float_Struc_BBP.BBP700));
  Float_Struc_Chl.chl_bl = NaN.*ones(size(Float_Struc_Chl.CHLA_ADJUSTED));
  Float_Struc_Chl.chl_bs = NaN.*ones(size(Float_Struc_Chl.CHLA_ADJUSTED));


%% Second for loop is used here after the background noise for the float is computed, to be applied to the filtered BBP_tot profiles %%%%%%%%%%
  for j=1:length(Float_Struc_Chl.PRES_ADJUSTED(1,:));
    BBP_tot_tmp = Float_Struc_BBP.BBP700(:,j);
    Pres_tmp = Float_Struc_Chl.PRES_ADJUSTED(:,j);
    chl_tot_tmp = Float_Struc_Chl.CHLA_ADJUSTED(:,j);
% zero out negative chl values

    chl_tot_tmp(find(chl_tot_tmp<0)) = 0;

% remove within 25 m of BBP_tot > thr 
    
    ind_hi = find(BBP_tot_tmp>thr);
    for k = 1:length(BBP_tot_tmp);
      if any(Pres_tmp(k) >Pres_tmp(ind_hi) -25) && any(Pres_tmp(k) <Pres_tmp(ind_hi) +25);
         BBP_tot_tmp(k) = NaN;
         Pres_tmp(k) = NaN;
         chl_tot_tmp(k) = NaN;
      end
    end

    indnonan = find(~isnan(BBP_tot_tmp) & ~isnan(chl_tot_tmp));

% remove NaN values
    Pres_tmp = Pres_tmp(indnonan);
    BBP_tot_tmp = BBP_tot_tmp(indnonan);
    chl_tot_tmp = chl_tot_tmp(indnonan);

% break out of loop if all NaNs
    if sum(~isnan(Pres_tmp)) == 0;
       continue
    end

    BBP_fil_tmp = movmin(BBP_tot_tmp,11); % 11 pt moving minimum filter
    BBP_fil_tmp = movmax(BBP_fil_tmp,11); % 11 pt moving minimum filter

    chl_fil_tmp = movmin(chl_tot_tmp,11); % 11 pt moving minimum filter
    chl_fil_tmp = movmax(chl_fil_tmp,11); % 11 pt moving minimum filter

    BBP_bsr = BBP_fil_tmp+Float_Struc_BBP.BBP700_noise; % refractory + small particles
    BBP_bl =  BBP_tot_tmp - BBP_bsr; % large particels
    BBP_bl(find(BBP_bl<0)) = 0;    

    chl_bsr = chl_fil_tmp+Float_Struc_Chl.CHLA_noise; % refractory + small particles
    chl_bl =  chl_tot_tmp - chl_bsr; % large particels
    chl_bl(find(chl_bl<0)) = 0;
    chl_bsr_850to900 = chl_bsr(find(Pres_tmp>850));
    chl_bsr_tot = vertcat(chl_bsr_tot,chl_bsr_850to900);

    BBP_bsr_850to900 = BBP_bsr(find(Pres_tmp>850));
    BBP_bsr_tot = vertcat(BBP_bsr_tot,BBP_bsr_850to900);

    Float_Struc_BBP.BBP700_bl(indnonan,j) = BBP_bl;
    Float_Struc_BBP.BBP700_bs(indnonan,j) = BBP_bsr; 
    Float_Struc_Chl.chl_bl(indnonan,j) = chl_bl;
    Float_Struc_Chl.chl_bs(indnonan,j) = chl_bsr;
  end

  BBP_bsr_tot = BBP_bsr_tot(2:end);
  chl_bsr_tot = chl_bsr_tot(2:end);


% calcualte pctile_bsr of BBP_bsr between 850-900 m. Use last float value if no deep data exist
  if sum(~isnan(BBP_bsr_tot)) ~= 0;
    P_bbp = prctile(BBP_bsr_tot,pctile_bsr);
    P_chl = prctile(chl_bsr_tot,pctile_chl);
  end
% apply correction to obtain BBP_bs;
  Float_Struc_BBP.BBP700_bs = Float_Struc_BBP.BBP700_bs - P_bbp;  
  Float_Struc_BBP.BBP700_bs(find(Float_Struc_BBP.BBP700_bs<0)) = 0; 
  Float_Struc_BBP.BBP700_br = P_bbp; % refractory material is a single value for the entire float: like instrument noise
  Float_Struc_Chl.chl_bs = Float_Struc_Chl.chl_bs - P_chl;
  Float_Struc_Chl.chl_br = P_chl; % refractory material is a single value for the entire float: like instrument noise
  Float_Struc_BBP.PRES_ADJUSTED = Float_Struc_Chl.PRES_ADJUSTED;
  Float_Struc_BBP.BBP_700_mean = nanmean((Float_Struc_BBP.BBP700_bs + Float_Struc_BBP.BBP700_bl),2);
  Float_Struc_Chl.chlmean = nanmean((Float_Struc_BBP.BBP700_bs + Float_Struc_BBP.BBP700_bl),2);

  eval([Data_bbp,'_proc.',char(fieldnm),'= Float_Struc_BBP;']);    
  eval([Data_chl,'_proc.',char(fieldnm),'= Float_Struc_Chl;']);

end
save([oneargo,'/matfiles/',filename,'_processed.mat'],[Data_bbp,'_proc'],[Data_chl,'_proc']); 



  
