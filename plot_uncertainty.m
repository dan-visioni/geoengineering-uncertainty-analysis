clear all

n_exp = {'G6sulfur','G6solar','ssp585','ssp245'};
% list models
n_mod = {'CNRM-ESM2-1','IPSL-CM6A-LR','CESM2-WACCM','UKESM1-0-LL','MPI-ESM1-2-LR','MPI-ESM1-2-HR','SOCOL-MPI'}; %Available models

tas = ncread('../GeoMIP_G6Sulfur/tas_timeseries_CMIP6_GeoMIP.nc','tas');
load ../GeoMIP_G6Sulfur/tas_timeseries_SOCOL

figure(1)
set(gcf, 'Position',  [200, 200, 1200, 800])



fc1 = brewermap(7,'Set2');
fc2 = brighten(fc1,.8);

time = 2020:2100;

SSP585 = squeeze(tas(:,3,:)); SSP585(7,:) = x_glob_m_save_SOCOL(3,5:85);
SSP245 = squeeze(tas(:,4,:)); SSP245(7,:) = x_glob_m_save_SOCOL(4,5:85);
G6sulf = squeeze(tas(:,1,:)); G6sulf(7,:) = x_glob_m_save_SOCOL(1,5:85);

subplot(2,2,1)
box on
hold on
for i=1:7
    plot(time,SSP585(i,:)-SSP245(i,:),'Linewidth',4,'Color',fc1(i,:));
end
set(gca,'Linewidth',2,'FontSize',16,'FontWeight','Bold')
set(gca,'XTick',[2020:10:2100],'YTick',[0:1:300],'FontSize',16,'FontWeight','Bold')
ylabel('\DeltaT(^{\circ}C)')
axis([2020 2100 -.5 3.5])
title('Changes in surface temperature')

subplot(2,2,2)
box on
hold on
for i=1:6
    plot(time,G6sulf(i,:)-SSP245(i,:),'Linewidth',4,'Color',fc1(i,:));
end
set(gca,'Linewidth',2,'FontSize',16,'FontWeight','Bold')
set(gca,'XTick',[2020:10:2100],'YTick',[0:1:300],'FontSize',16,'FontWeight','Bold')
ylabel('DeltaT(^{\circ}C)')
axis([2020 2100 -.5 3.5])
title('Changes in surface temperature')

load ../GeoMIP_G6Sulfur/CO2_concentrations.mat

subplot(2,2,3)
plot(2020:1:2100,CO2_SSP585(6:86),'b','Linewidth',2)
hold on
box on
plot(2020:1:2100,CO2_SSP245(6:86),'r','Linewidth',2)
plot(2020:1:2100,CO2_SSP370(6:86),'Color',[.7 .7 .7],'Linewidth',2)
plot(2020:1:2100,CO2_SSP126(6:86),'Color',[.7 .7 .7],'Linewidth',2)
plot(2020:1:2100,CO2_SSP534(6:86),'Color',[.7 .7 .7],'Linewidth',2)
plot(2020:1:2100,CO2_SSP585(6:86),'b','Linewidth',2)
plot(2020:1:2100,CO2_SSP245(6:86),'r','Linewidth',2)
set(gca,'Linewidth',2,'FontSize',16,'FontWeight','Bold')
set(gca,'XTick',[2020:10:2100],'YTick',[0:100:1200],'FontSize',18,'FontWeight','Bold')
axis([2020 2100 400 1200])
ylabel('ppm')
title('CO_2 concentrations in SSP scenarios')

load ../GeoMIP_G6Sulfur/global_AOD.mat
load ../GeoMIP_G6Sulfur/aod_SOCOL.mat

subplot(2,2,4)
box on
hold on
AOD(1,1) = AOD(2,1);
AOD(:,7) = NaN(1,81);
AOD(12:81,7) = aod_SOCOL;
for i=1:7
    plot(time,AOD(:,i),'Linewidth',4,'Color',fc1(i,:));
end
%plot(2031:2100,aod_SOCOL,'Linewidth',4,'Color',fc1(7,:));
ylabel('Global AOD')
set(gca,'Linewidth',2,'FontSize',16,'FontWeight','Bold')
set(gca,'XTick',[2020:10:2100],'YTick',[0:.05:1],'FontSize',18,'FontWeight','Bold')
axis([2020 2100 0 .4])
title('Changes in stratospheric AOD')

set(gcf,'renderer','painters')
print(gcf,'-depsc2',['uncertainties.eps'])

%%

figure(4)
set(gcf, 'Position',  [200, 200, 1200, 800])
time = 2030:2100;
subplot(2,2,1)
box on
hold on
for i=1:7
    p(i) = plot(time,movmean((SSP585(i,11:81)-G6sulf(i,11:81))./AOD(11:81,i)',3),'Linewidth',2,'Color',fc1(i,:));
end

tas = movmean(nanmean((SSP585(:,11:81)-G6sulf(:,11:81))./AOD(11:81,:)'),3);
tasS = movmean(nanstd((SSP585(:,11:81)-G6sulf(:,11:81))./AOD(11:81,:)',0,1),3)/sqrt(7);

mean(tas(51:71))
mean(tasS(51:71))

plot(time,tas,'Linewidth',4,'Color','k');

Y_std = [tas+tasS,fliplr(tas-tasS)];

X = [(time),fliplr(time)];

fill(X,(Y_std),'k','LineStyle','none','facealpha',.2)

hl=legend(p(1:7),n_mod,'Location','eastoutside');
hl.Position(1)=hl.Position(1)+0.15;

set(gca,'Linewidth',2,'FontSize',16,'FontWeight','Bold')
set(gca,'XTick',[2020:10:2100],'YTick',[0:1:300],'FontSize',16,'FontWeight','Bold')
ylabel('^{\circ}C/AOD')
axis([2040 2100 -1 10])
title('Cooling efficacy per unit AOD')

set(gcf,'renderer','painters')
print(gcf,'-depsc2',['cooling_eff.eps'])