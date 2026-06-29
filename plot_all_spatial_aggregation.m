clear all

load 1deg_tas.mat
sola_st1d = sola_st;
ssp_st1d = ssp_st;
sulf_st1d = sulf_st;

lat_int1 = -90:1:90;
lon_int1 = 0:1:360;

load 2deg_tas.mat
sola_st2d = sola_st;
ssp_st2d = ssp_st;
sulf_st2d = sulf_st;

lat_int2 = -90:1:90;
lon_int2 = 0:1:360;

load 5deg_tas.mat
sola_st5d = sola_st;
ssp_st5d = ssp_st;
sulf_st5d = sulf_st;

lat_int5 = -90:1:90;
lon_int5 = 0:1:360;

load 10deg_tas.mat
sola_st10d = sola_st;
ssp_st10d = ssp_st;
sulf_st10d = sulf_st;

lat_int10 = -90:1:90;
lon_int10 = 0:1:360;

load 15deg_tas.mat
sola_st15d = sola_st;
ssp_st15d = ssp_st;
sulf_st15d = sulf_st;

lat_int15 = -90:1:90;
lon_int15 = 0:1:360;

load 05deg_tas.mat
sola_st05d = sola_st;
ssp_st05d = ssp_st;
sulf_st05d = sulf_st;

lat_int05 = -90:1:90;
lon_int05 = 0:1:360;

figure(1)
set(gcf, 'Position',  [200, 200,1600,800])

l_colb = 4;
fc = brewermap(l_colb,'Spectral');

subplot(1,2,1)

hold on

ssp = abs(ssp_st1d)./sola_st1d*100;
sul = abs(sulf_st1d)./sola_st1d*100;

scatter(ssp,sul,15,'filled','MarkerFaceColor',fc(1,:))


ssp = abs(ssp_st2d)./sola_st2d*100;
sul = abs(sulf_st2d)./sola_st2d*100;

scatter(ssp,sul,15,'filled','MarkerFaceColor',fc(2,:))



ssp = abs(ssp_st5d)./sola_st5d*100;
sul = abs(sulf_st5d)./sola_st5d*100;

scatter(ssp,sul,15,'filled','MarkerFaceColor',fc(3,:))

ssp = abs(ssp_st10d)./sola_st10d*100;
sul = abs(sulf_st10d)./sola_st10d*100;

scatter(ssp,sul,15,'filled','MarkerFaceColor',fc(4,:))

%%

ssp = abs(ssp_st2d)./sola_st2d*100;
sul = abs(sulf_st2d)./sola_st2d*100;
median(ssp,'omitnan')
median(sul,'omitnan')

scatter(median(ssp,'omitnan'),median(sul,'omitnan'),150,'filled','k')
errorbar(median(ssp,'omitnan'),median(sul,'omitnan'),prctile(ssp,2),prctile(ssp,2),prctile(sul,2),prctile(sul,2),'k','linewidth',2)

set(gca,'XTick',[0:100:1000],'YTick',[0:100:1000],'FontSize',16,'FontWeight','Bold')
set(gca,'XTick',[50 100 150 200 500 1000],'YTick',[50 100 150 200 500 1000],'FontSize',16,'FontWeight','Bold')
plot([50 1000],[50 1000],'Color',[0 0 0],'Linewidth',2,'Linestyle','--')
plot([50 100],[100 100],'Color',[0 0 0],'Linewidth',2,'Linestyle','-')
plot([100 100],[50 100],'Color',[0 0 0],'Linewidth',2,'Linestyle','-')
axis([50 1000 50 1000])
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')

xlabel('Increase in variance due to GHG warming (%)');
ylabel('Increase in variance due to aerosol cooling (%)');
title({'a) Surface temperature uncertainties from GHG and SAI'},'Color','k');
set(gca,'Linewidth',2,'FontSize',18,'FontWeight','Bold')

box on

set(gcf, 'Color', 'w'); % Set figure background color to white
set(gca, 'Color', 'w')
%set(hl, 'Color', 'w','TextColor','k')
set(gca,'XColor','k','YColor','k');

subplot(1,2,2)

load 1deg_pr.mat
sola_st1d = sola_st;
ssp_st1d = ssp_st;
sulf_st1d = sulf_st;

lat_int1 = -90:1:90;
lon_int1 = 0:1:360;

load 2deg_pr.mat
sola_st2d = sola_st;
ssp_st2d = ssp_st;
sulf_st2d = sulf_st;

lat_int2 = -90:1:90;
lon_int2 = 0:1:360;

load 5deg_pr.mat
sola_st5d = sola_st;
ssp_st5d = ssp_st;
sulf_st5d = sulf_st;

lat_int5 = -90:1:90;
lon_int5 = 0:1:360;

load 10deg_pr.mat
sola_st10d = sola_st;
ssp_st10d = ssp_st;
sulf_st10d = sulf_st;

lat_int10 = -90:1:90;
lon_int10 = 0:1:360;

load 15deg_pr.mat
sola_st15d = sola_st;
sssp_st15d = ssp_st;
sulf_st15d = sulf_st;

lat_int15 = -90:1:90;
lon_int15 = 0:1:360;

load 05deg_pr.mat
sola_st05d = sola_st;
ssp_st05d = ssp_st;
sulf_st05d = sulf_st;

hold on

ssp = abs(ssp_st1d)./sola_st1d*100;
sul = abs(sulf_st1d)./sola_st1d*100;

p(1) = scatter(ssp,sul,15,'filled','MarkerFaceColor',fc(1,:));


ssp = abs(ssp_st2d)./sola_st2d*100;
sul = abs(sulf_st2d)./sola_st2d*100;

p(2) = scatter(ssp,sul,15,'filled','MarkerFaceColor',fc(2,:));


ssp = abs(ssp_st5d)./sola_st5d*100;
sul = abs(sulf_st5d)./sola_st5d*100;

p(3) = scatter(ssp,sul,15,'filled','MarkerFaceColor',fc(3,:));

ssp = abs(ssp_st10d)./sola_st10d*100;
sul = abs(sulf_st10d)./sola_st10d*100;

p(4) = scatter(ssp,sul,15,'filled','MarkerFaceColor',fc(4,:));


ssp_pr = abs(ssp_st2d)./sola_st2d*100;
sul_pr = abs(sulf_st2d)./sola_st2d*100;
median(ssp_pr,'omitnan')
median(sul_pr,'omitnan')

p(5) = scatter(median(ssp_pr,'omitnan'),median(sul_pr,'omitnan'),150,'filled','k');
errorbar(median(ssp_pr,'omitnan'),median(sul_pr,'omitnan'),prctile(ssp_pr,2),prctile(ssp_pr,2),prctile(sul_pr,2),prctile(sul_pr,2),'k','linewidth',2)

plot([50 1000],[50 1000],'Color',[0 0 0],'Linewidth',2,'Linestyle','--')
plot([50 100],[100 100],'Color',[0 0 0],'Linewidth',2,'Linestyle','-')
plot([100 100],[50 100],'Color',[0 0 0],'Linewidth',2,'Linestyle','-')

hl = legend({'1 degree resolution','2 degrees resolution','5 degrees resolution','10 degrees resolution','Median value (2 degrees)'},...
    'Location','NorthWest');


set(gca,'XTick',[0:100:1000],'YTick',[0:100:1000],'FontSize',16,'FontWeight','Bold')
set(gca,'XTick',[50 100 150 200 500 1000],'YTick',[50 100 150 200 500 1000],'FontSize',16,'FontWeight','Bold')

axis([50 1000 50 1000])
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')

xlabel('Increase in variance due to GHG warming (%)');
ylabel('Increase in variance due to aerosol cooling (%)');
title({'b) Precipitation uncertainties from GHG and SAI'},'Color','k');
set(gca,'Linewidth',2,'FontSize',18,'FontWeight','Bold')

box on

set(gcf, 'Color', 'w'); % Set figure background color to white
set(gca, 'Color', 'w')
set(hl, 'Color', 'w','TextColor','k')
set(gca,'XColor','k','YColor','k');

set(gcf,'renderer','painters')
print(gcf,'-depsc2',['GeoMIP_spatial_uncertainties.eps'])