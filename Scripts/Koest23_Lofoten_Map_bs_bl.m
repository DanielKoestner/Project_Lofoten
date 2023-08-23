% Get plot floats in Lofoten basin region from ~2010–2023
% Plotted variables include:
% bbp: total backscattering at 700 nm
% bbp_bs: backscattering of small particles at 700 nm
% bbp_bl: backscattering of large particles at 700 nm
% chla: chlorophyll-a from flourescence
% chla_bs: chlorophyll-a from flourescence in small particles
% chla_bl: chlorophyll-a from flourescence in large particles
% chla/bbp: rato of chlorophyll-a to total backscattering at 700 nm

% Requries OneArgo-Mat Package. Please change the path pointed to OneArgo
% directoires befor running
close all
clear all
curdir=cd;
oneargo = '/discover/nobackup/plerner/GO-BGC-ARGO/BGC_ARGO_toolbox/OneArgo-Mat';



addpath(oneargo); % just in case, add to MATLAB search path
addpath(strcat(oneargo,'/m_map'));

print_flag=0;
fs=11;
lw=1.5;
set(0, 'DefaultAxesFontName', 'Times');
set(0, 'DefaultTextFontName', 'Times');
alp=0.65;

hf1=figure();
set(hf1,'Units','inches','Position', [5 5 8 7], 'PaperPosition', [0 0 8 7], 'PaperSize', [8 7]);
ha1=iSubplot(1,1, 'Gap', [0.04 0], 'Min', [.05 0.05], 'Max', [0.98 0.96], 'XTickL', 'All', 'YTickL', 'All');

LATLIMS=[64 74];
LONLIMS=[-10 20];

m_proj('lambert','long',LONLIMS,'lat',LATLIMS);
% m_coast('patch',[1 .85 .7]);
%     m_gshhs_l('patch',[.5 .6 .5]);

axes(ha1(1)); hold on

[CS,CH]=m_elev('contourf',[-4000 -3000:100:0 ],'edgecolor','none');
m_gshhs_l('patch',[.8 .8 .8]);
  m_grid('linewi',2,'tickdir','out');
cols=flipud(bone);
% cols=cols(20:end,:);
% colormap(cols);
% colormap(m_colmap('blue',30));
colormap(gray(30));
  h=colorbar;
  set(get(h,'ylabel'),'String','bathymetry [m]');
  h.Location='southoutside';
    set(h,'tickdir','out');
% brighten(.5);

title(['Lofoten Basin, 2010 – 2023'])



%% add floats!
% valid as of 22 August 2023, used OneArgo_Lofoten code first

%load('Data/LBE_BGCDATA_2010_2023.mat');
load([oneargo,'/matfiles/BBP_CHL_BGC-ARGO_processed.mat']);
f=length(fieldnames(Data_qc_BBP_proc));
names=fieldnames(Data_qc_BBP_proc);
for i = 1:f
    temp=eval(['Data_qc_BBP_proc.',char(names{i})]);

    % inds=~isnan(temp.BBP700_ADJUSTED) | ~isnan(temp.CHLA_ADJUSTED);
    floats{i}.bbp700=temp.BBP700-nanmin(temp.BBP700);
    floats{i}.bbp700_bs = temp.BBP700_bs;
    floats{i}.bbp700_bl = temp.BBP700_bl; 
    temp=eval(['Data_qc_Chl_proc.',char(names{i})]);
    floats{i}.chla=temp.CHLA_ADJUSTED-nanmin(temp.CHLA_ADJUSTED);
    floats{i}.chla_bs = temp.chl_bs;
    floats{i}.chla_bl = temp.chl_bl;
    floats{i}.pres=temp.PRES_ADJUSTED;
    floats{i}.lat=temp.LATITUDE(1,:);
    floats{i}.lon=temp.LONGITUDE(1,:);
    floats{i}.juld=temp.TIME(1,:);
    floats{i}.wmo=str2num(names{i}(2:8));

end

colors=crameri('batlow',f);

cnt=1;
while cnt<=length(floats)
    floatcol=colors(cnt,:);

    m_track(floats{cnt}.lon,floats{cnt}.lat,...
        'ticks',259200,'times',0,'dates',0,'clip','off','color',...
        floatcol,'linew',1.1,'orien','upright');
    % ticks every 6 months, 

    text(0.97,1-cnt*.03,sprintf('%1.0f',str2num(names{cnt}(2:8))),'units','normalized','color',floatcol,'fontweight','bold');
    

    cnt=cnt+1;
end
% text(.06,1,'WMO ID','units','normalized','color','k','fontweight','bold');

cnt=1;
while cnt<=length(floats)

    floatcol=colors(cnt,:);
        m_plot(floats{cnt}.lon(1),floats{cnt}.lat(1),'marker','.','markerfacecolor',floatcol,'markeredgecolor',floatcol,'markersize',12);
        m_plot(floats{cnt}.lon(end),floats{cnt}.lat(end),'marker','.','markerfacecolor',floatcol,'markeredgecolor',floatcol,'markersize',12);

    cnt=cnt+1;
end

  m_grid('linewi',2,'tickdir','out');


  %% lims

  bbplims=([0.00001 0.03]);
  chlalims=([0.01 5]);
  complims=([2 2000]);


  %% Plot bbp_bs vs depth

hf5=figure();
set(hf5,'Units','inches','Position', [5 5 8 10], 'PaperPosition', [0 0 8 10], 'PaperSize', [8 10]);
ha5=iSubplot(7,4, 'Gap', [0 0.01], 'Min', [0.03 0.03], 'Max', [0.98 0.98], 'XTickL', 'All', 'YTickL', 'All');

c=crameri('romao',12);
colors2=[c(9:12,:); c(1:8,:)];

for i = 1:length(floats)
    floatcol=colors(i,:);
    axes(ha5(i));
    box on
    hold on
    ylim([-2000 -5]);

    if i == 1 || i ==5 || i == 9 || i ==13 || i ==17 || i ==21 || i==25
        ylabel('\itp \rm[dbar]')
        set(gca,'ytick',[-2000 -1000 -100 -5],'yticklabel',{'','1000','100','5'})
    else
        set(gca,'ytick',[-2000 -1000 -100 -5])
        set(gca,'yticklabel',{''})
    end

    if i == 25 || i == 26 || i == 27 || i ==28
        xlabel('b_{bpbs}(700) [m^{-1}]')
    end


    for ii = 1:length(floats{i}.juld)
        mn=month(datetime(floats{i}.juld(ii),'ConvertFrom','datenum','Format','yyyy-MM-dd'));
        s=scatter(floats{i}.bbp700_bs(:,ii),floats{i}.pres(:,ii)*-1,5,colors2(mn,:),'filled','MarkerFaceAlpha',alp);
        clear mc mn
    end

    title([num2str(floats{i}.wmo) '; ' datestr(floats{i}.juld(1),'yyyymmdd')] ,'color',floatcol)
hold off

set(gca,'yscale','log')
set(gca,'TickLength',[0.04 0.04])
set(gca,'fontsize',fs-2)

ha5(i).XLim=bbplims;
set(gca,'xscale','log')
set(gca,'xtick',[0.0001 0.001 0.01])

end

axes(ha5(27));
hold on
for i = 1:12
        s(i)=scatter(i,i,50,colors2(i,:),'filled');
end
set(gca,'yticklabel',{''})
ha5(15).XLim=([0.0001 0.03]);
ylim([-2000 0]);

lgd_text={'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
[hl,lines]=legendflex([s],lgd_text, 'ref', gca, ...
'anchor', {'nw', 'nw'}, 'buffer', [0 0], 'xscale', 0.5, 'ncol',4,'nrow',3, ...
'box', 'off', 'FontSize', fs);

ha5(27).Visible='off';
ha5(28).Visible='off';


  %% Plot bbp_bl vs depth

hf6=figure();
set(hf6,'Units','inches','Position', [5 5 8 10], 'PaperPosition', [0 0 8 10], 'PaperSize', [8 10]);
ha6=iSubplot(7,4, 'Gap', [0 0.01], 'Min', [0.03 0.03], 'Max', [0.98 0.98], 'XTickL', 'All', 'YTickL', 'All');

c=crameri('romao',12);
colors2=[c(9:12,:); c(1:8,:)];

for i = 1:length(floats)
    floatcol=colors(i,:);
    axes(ha6(i));
    box on
    hold on
    ylim([-2000 -5]);

    if i == 1 || i ==5 || i == 9 || i ==13 || i ==17 || i ==21 || i==25
        ylabel('\itp \rm[dbar]')
        set(gca,'ytick',[-2000 -1000 -100 -5],'yticklabel',{'','1000','100','5'})
    else
        set(gca,'ytick',[-2000 -1000 -100 -5])
        set(gca,'yticklabel',{''})
    end

    if i == 25 || i == 26 || i == 27 || i ==28
        xlabel('b_{bpbl}(700) [m^{-1}]')
    end


    for ii = 1:length(floats{i}.juld)
        mn=month(datetime(floats{i}.juld(ii),'ConvertFrom','datenum','Format','yyyy-MM-dd'));
        s=scatter(floats{i}.bbp700_bl(:,ii),floats{i}.pres(:,ii)*-1,5,colors2(mn,:),'filled','MarkerFaceAlpha',alp);
        clear mc mn
    end

    title([num2str(floats{i}.wmo) '; ' datestr(floats{i}.juld(1),'yyyymmdd')] ,'color',floatcol)
hold off

set(gca,'yscale','log')
set(gca,'TickLength',[0.04 0.04])
set(gca,'fontsize',fs-2)

ha6(i).XLim=bbplims;
set(gca,'xscale','log')
set(gca,'xtick',[0.0001 0.001 0.01])

end

axes(ha6(27));
hold on
for i = 1:12
        s(i)=scatter(i,i,50,colors2(i,:),'filled');
end
set(gca,'yticklabel',{''})
ha6(15).XLim=([0.0001 0.03]);
ylim([-2000 0]);

lgd_text={'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
[hl,lines]=legendflex([s],lgd_text, 'ref', gca, ...
'anchor', {'nw', 'nw'}, 'buffer', [0 0], 'xscale', 0.5, 'ncol',4,'nrow',3, ...
'box', 'off', 'FontSize', fs);

ha6(27).Visible='off';
ha6(28).Visible='off';






  %% Plot chla in small particles vs depth

hf7=figure();
set(hf7,'Units','inches','Position', [5 5 8 10], 'PaperPosition', [0 0 8 10], 'PaperSize', [8 10]);
ha7=iSubplot(7,4, 'Gap', [0 0.01], 'Min', [0.03 0.03], 'Max', [0.98 0.98], 'XTickL', 'All', 'YTickL', 'All');


for i = 1:length(floats)
    floatcol=colors(i,:);
    axes(ha7(i));
    box on
    hold on
    ylim([-2000 -5]);

    if i == 1 || i ==5 || i == 9 || i ==13 || i ==17 || i ==21 || i==25
        ylabel('\itp \rm[dbar]')
        set(gca,'ytick',[-2000 -1000 -100 -5],'yticklabel',{'','1000','100','5'})
    else
        set(gca,'ytick',[-2000 -1000 -100 -5])
        set(gca,'yticklabel',{''})
    end

    if i == 25 || i == 26 || i == 27 || i ==28
        xlabel('chla_{bs} [mg m^{-3}]')
    end


    for ii = 1:length(floats{i}.juld)
        mn=month(datetime(floats{i}.juld(ii),'ConvertFrom','datenum','Format','yyyy-MM-dd'));;
        s=scatter(floats{i}.chla_bs(:,ii),floats{i}.pres(:,ii)*-1,5,colors2(mn,:),'filled','MarkerFaceAlpha',alp);
        clear mc mn
    end

    title([num2str(floats{i}.wmo) '; ' datestr(floats{i}.juld(1),'yyyymmdd')] ,'color',floatcol)
hold off

set(gca,'yscale','log')
set(gca,'TickLength',[0.04 0.04])
set(gca,'fontsize',fs-2)

ha7(i).XLim=chlalims;
set(gca,'xscale','log')
set(gca,'xtick',[0.1 1 5])
end

axes(ha7(27));
hold on
for i = 1:12
        s(i)=scatter(i,i,50,colors2(i,:),'filled');
end
set(gca,'yticklabel',{''})
%ha7(15).XLim=([0.0001 0.03]);
ylim([-2000 0]);

lgd_text={'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
[hl,lines]=legendflex([s],lgd_text, 'ref', gca, ...
'anchor', {'nw', 'nw'}, 'buffer', [0 0], 'xscale', 0.5, 'ncol',4,'nrow',3, ...
'box', 'off', 'FontSize', fs);

ha7(27).Visible='off';
ha7(28).Visible='off';

  %% Plot chla_bl in small particles vs depth

hf8=figure();
set(hf8,'Units','inches','Position', [5 5 8 10], 'PaperPosition', [0 0 8 10], 'PaperSize', [8 10]);
ha8=iSubplot(7,4, 'Gap', [0 0.01], 'Min', [0.03 0.03], 'Max', [0.98 0.98], 'XTickL', 'All', 'YTickL', 'All');


for i = 1:length(floats)
    floatcol=colors(i,:);
    axes(ha8(i));
    box on
    hold on
    ylim([-2000 -5]);

    if i == 1 || i ==5 || i == 9 || i ==13 || i ==17 || i ==21 || i==25
        ylabel('\itp \rm[dbar]')
        set(gca,'ytick',[-2000 -1000 -100 -5],'yticklabel',{'','1000','100','5'})
    else
        set(gca,'ytick',[-2000 -1000 -100 -5])
        set(gca,'yticklabel',{''})
    end

    if i == 25 || i == 26 || i == 27 || i ==28
        xlabel('chla_{bl} [mg m^{-3}]')
    end


    for ii = 1:length(floats{i}.juld)
        mn=month(datetime(floats{i}.juld(ii),'ConvertFrom','datenum','Format','yyyy-MM-dd'));;
        s=scatter(floats{i}.chla_bl(:,ii),floats{i}.pres(:,ii)*-1,5,colors2(mn,:),'filled','MarkerFaceAlpha',alp);
        clear mc mn
    end

    title([num2str(floats{i}.wmo) '; ' datestr(floats{i}.juld(1),'yyyymmdd')] ,'color',floatcol)
hold off

set(gca,'yscale','log')
set(gca,'TickLength',[0.04 0.04])
set(gca,'fontsize',fs-2)

ha8(i).XLim=chlalims;
set(gca,'xscale','log')
set(gca,'xtick',[0.1 1 5])
end

axes(ha8(27));
hold on
for i = 1:12
        s(i)=scatter(i,i,50,colors2(i,:),'filled');
end
set(gca,'yticklabel',{''})
%ha8(15).XLim=([0.0001 0.03]);
ylim([-2000 0]);

lgd_text={'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
[hl,lines]=legendflex([s],lgd_text, 'ref', gca, ...
'anchor', {'nw', 'nw'}, 'buffer', [0 0], 'xscale', 0.5, 'ncol',4,'nrow',3, ...
'box', 'off', 'FontSize', fs);

ha8(27).Visible='off';
ha8(28).Visible='off';





%% depth integrate

%Datai = depth_interp(Data, qc_flags, varargin)
 %% Print
 cd(curdir)
if print_flag==1
    figure(hf5)
    saveas(figure(hf5),[oneargo,'/Figures/Lerner23_Lofoten_bbps_bs.pdf'])
    figure(hf6)
    saveas(figure(hf6),[oneargo,'/Figures/Lerner23_Lofoten_bbps_bl.pdf']);
    figure(hf7)
    saveas(figure(hf7),[oneargo,'/Figures/Lerner23_Lofoten_chlas_bs.pdf']);
    figure(hf8)
    saveas(figure(hf8),[onergo,'/Figures/Lerner23_Lofoten_chlas_bl.pdf']);

end
