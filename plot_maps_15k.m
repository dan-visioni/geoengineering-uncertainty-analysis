clear all

fold  = '/Volumes/D_Visioni2/GeoMIP/G6-1.5K-SAI/';

nam = 'b.e21.BW.f09_g17.SSP245-G6-1p5K-SAI.';
nam2 = 'b.e21.BWSSP245cmip6.f09_g17.CMIP6-SSP2-4.5-WACCM.';

ens = {'001','002','003','004','005','006','007','008','009','010'};

var = 'TREFHT';

for ie=1:3
    ta_CESM_g6(:,:,:,ie) = ncread([fold nam ens{ie} '.cam.h0.TREFHT.203501-208412.nc'],var);
end

lat_CESM = ncread([fold nam ens{ie} '.cam.h0.TREFHT.203501-208412.nc'],'lat');
lon_CESM = ncread([fold nam ens{ie} '.cam.h0.TREFHT.203501-208412.nc'],'lon');

for ie=1:3
    ta_CESM_ssp(:,:,:,ie) = ncread([fold nam2 ens{ie} '.cam.h0.TREFHT.201501-210012.nc'],var);
end

var = 'SurfT';

ens = {'r01','r01','r03','r04','r05','r06','r07','r08','r09','r10'};

for ie=1:10
    ta_MIROC_g6(:,:,:,ie) = ncread([fold 'SurfT_G6-1.5K-SAI_' ens{ie} '.nc'],var);
end

lat_MIROC = ncread([fold 'SurfT_G6-1.5K-SAI_' ens{ie} '.nc'],'lat');
lon_MIROC = ncread([fold 'SurfT_G6-1.5K-SAI_' ens{ie} '.nc'],'lon');

for ie=1:10
    ta_MIROC_ssp(:,:,:,ie) = ncread([fold 'SurfT_Baseline_' ens{ie} '.nc'],var);
end

var = 'temp';

ens = {'001','002','003','004','005','006','007','008','009','010'};

for ie=1:3
    ta_UKESM_g6(:,:,:,ie) = ncread([fold 'T1.5m_G6-1p5K-SAI.' ens{ie} '.nc'],var);
end

lat_UKESM = double(ncread([fold 'T1.5m_G6-1p5K-SAI.' ens{ie} '.nc'],'latitude'));
lon_UKESM = double(ncread([fold 'T1.5m_G6-1p5K-SAI.' ens{ie} '.nc'],'longitude'));

for ie=1:3
    ta_UKESM_ssp(:,:,:,ie) = ncread([fold 'T1.5m_SSP245.' ens{ie} '.nc'],var);
end



%%

ta_CESM_g6_em = squeeze(mean(ta_CESM_g6,4));
ta_CESM_ssp_em = squeeze(mean(ta_CESM_ssp,4));

ta_CESM_g6_zm = squeeze(mean(ta_CESM_g6_em(:,:,461:600),3));
ta_CESM_ssp_zm = squeeze(mean(ta_CESM_ssp_em(:,:,1:240),3));
ta_CESM_ssp_future_zm = squeeze(mean(ta_CESM_ssp_em(:,:,601:840),3));

ta_MIROC_g6_em = squeeze(mean(ta_MIROC_g6,4));
ta_MIROC_ssp_em = squeeze(mean(ta_MIROC_ssp,4));

ta_MIROC_g6_zm = squeeze(mean(ta_MIROC_g6_em(:,:,461:600),3));
ta_MIROC_ssp_zm = squeeze(mean(ta_MIROC_ssp_em(:,:,1:240),3));
ta_MIROC_ssp_future_zm = squeeze(mean(ta_MIROC_ssp_em(:,:,541:780),3));

ta_UKESM_g6_em = single(squeeze(mean(ta_UKESM_g6,4)));
ta_UKESM_ssp_em = single(squeeze(mean(ta_UKESM_ssp,4)));

ta_UKESM_g6_zm = squeeze(mean(ta_UKESM_g6_em(:,:,473:612),3));
ta_UKESM_ssp_zm = squeeze(mean(ta_UKESM_ssp_em(:,:,1:240),3));
ta_UKESM_ssp_future_zm = squeeze(mean(ta_UKESM_ssp_em(:,:,601:840),3));


l_colb = 22;
fc = brewermap(l_colb,'*RdBu');
mmin=-1.65*3;
mmax=1.65*3;
v = mmin:(mmax-mmin)/l_colb:mmax;
v2 = (mmin:(mmax-mmin)/(l_colb/2):mmax-(mmax-mmin)/(l_colb/2))+(mmax-mmin)/(l_colb);

l_colb = 10;
fcs = brewermap(l_colb,'Reds');
mmins=0;
mmaxs=5;
vs = mmins:(mmaxs-mmins)/l_colb:mmaxs;
v2s = mmins:(mmaxs-mmins)/(l_colb):mmaxs-(mmaxs-mmins)/(l_colb);

l_colb2 = 10;
fca = brewermap(10,'Oranges');
mmin2=0;
mmax2=.5;
v_2 = mmin2:(mmax2-mmin2)/l_colb2:mmax2;
v2_2 = mmin2:(mmax2-mmin2)/(l_colb2):mmax2-(mmax2-mmin2)/(l_colb2);

%%

figure(10)

colormap(fc)
set(gcf, 'Position',  [200, 200, 1600, 1600])

s1 = subplot(3,3,1);
colormap(gca,fc)
worldmap([-90 90],[-180 180])
box on
hold on
pcolorm(lat_CESM,lon_CESM,(ta_CESM_ssp_future_zm-ta_CESM_ssp_zm)','LineStyle','none')
load coastlines
geoshow(coastlat,coastlon)
caxis([mmin mmax])
set(gca,'Linewidth',2)
title({'a) SSP2-4.5';'CESM (2050-2069 compared to 1.5K)'},'FontSize',20,'FontWeight','Bold')
%s1.Position(1) = s1.Position(1)+.102;
hl = colorbar('FontSize',16,'YTick',(v2),'Linewidth',2,'Location','westoutside');
hl.Position(1) = hl.Position(1)-.01;
hl.Position(2) = hl.Position(2)-.2;
hl.Position(4) = hl.Position(4)*1.6;
ylabel(hl,'Surface air temperature (^{\circ}C)')

s3=subplot(3,3,4);

colormap(gca,fc)
worldmap([-90 90],[-180 180])
box on
hold on
pcolorm(lat_CESM',lon_CESM,(ta_CESM_g6_zm-ta_CESM_ssp_zm)','LineStyle','none')
load coastlines
geoshow(coastlat,coastlon)
caxis([mmin mmax])
set(gca,'Linewidth',2)
title({'a) G6-1.5K-SAI';'CESM (2050-2069 compared to 1.5K)'},'FontSize',20,'FontWeight','Bold')
%s1.Position(1) = s1.Position(1)+.102;
%hl = colorbar('FontSize',16,'YTick',(v2),'Linewidth',2,'Location','westoutside');
%hl.Position(1) = hl.Position(1)-.01;
%hl.Position(2) = hl.Position(2)-.2;
%hl.Position(4) = hl.Position(4)*1.6;
%ylabel(hl,'Surface air temperature (^{\circ}C)')

%set(gcf,'renderer','painters')
%print(gcf,'-depsc2',['maps_mean_15k.eps'])

s1 = subplot(3,3,2);
colormap(gca,fc)
worldmap([-90 90],[-180 180])
box on
hold on
pcolorm(lat_MIROC,lon_MIROC,(ta_MIROC_ssp_future_zm-ta_MIROC_ssp_zm)','LineStyle','none')
load coastlines
geoshow(coastlat,coastlon)
caxis([mmin mmax])
set(gca,'Linewidth',2)
title({'a) SSP2-4.5';'MIROC (2050-2069 compared to 1.5K)'},'FontSize',20,'FontWeight','Bold')
%s1.Position(1) = s1.Position(1)+.102;
hl = colorbar('FontSize',16,'YTick',(v2),'Linewidth',2,'Location','westoutside');
hl.Position(1) = hl.Position(1)-.01;
hl.Position(2) = hl.Position(2)-.2;
hl.Position(4) = hl.Position(4)*1.6;
ylabel(hl,'Surface air temperature (^{\circ}C)')

s3=subplot(3,3,5);

colormap(gca,fc)
worldmap([-90 90],[-180 180])
box on
hold on
pcolorm(lat_MIROC',lon_MIROC,(ta_MIROC_g6_zm-ta_MIROC_ssp_zm)','LineStyle','none')
load coastlines
geoshow(coastlat,coastlon)
caxis([mmin mmax])
set(gca,'Linewidth',2)
title({'a) G6-1.5K-SAI';'MIROC (2050-2069 compared to 1.5K)'},'FontSize',20,'FontWeight','Bold')



s3=subplot(3,3,3);
colormap(gca,fc)
worldmap([-90 90],[-180 180])
box on
hold on
pcolorm(lat_UKESM,lon_UKESM,(ta_UKESM_ssp_future_zm-ta_UKESM_ssp_zm)','LineStyle','none')
load coastlines
geoshow(coastlat,coastlon)
caxis([mmin mmax])
set(gca,'Linewidth',2)
title({'a) SSP2-4.5';'UKESM (2050-2069 compared to 1.5K)'},'FontSize',20,'FontWeight','Bold')
%s1.Position(1) = s1.Position(1)+.102;
hl = colorbar('FontSize',16,'YTick',(v2),'Linewidth',2,'Location','westoutside');
hl.Position(1) = hl.Position(1)-.01;
hl.Position(2) = hl.Position(2)-.2;
hl.Position(4) = hl.Position(4)*1.6;
ylabel(hl,'Surface air temperature (^{\circ}C)')

s3=subplot(3,3,6);

colormap(gca,fc)
worldmap([-90 90],[-180 180])
box on
hold on
pcolorm(lat_UKESM',lon_UKESM,(ta_UKESM_g6_zm-ta_UKESM_ssp_zm)','LineStyle','none')
load coastlines
geoshow(coastlat,coastlon)
caxis([mmin mmax])
set(gca,'Linewidth',2)
title({'a) G6-1.5K-SAI';'UKESM (2050-2069 compared to 1.5K)'},'FontSize',20,'FontWeight','Bold')

%%

lat_int1 = -90:1:90;
lon_int1 = 0:1:360;

nmodels_an = 3;

LAT = {lat_CESM,lat_MIROC,lat_UKESM};
LON = {lon_CESM,lon_MIROC,lon_UKESM};
ta_g6 = {ta_CESM_g6_zm,ta_MIROC_g6_zm,ta_UKESM_g6_zm};
ta_ba = {ta_CESM_ssp_zm,ta_MIROC_ssp_zm,ta_UKESM_ssp_zm};
ta_sp = {ta_CESM_ssp_future_zm,ta_MIROC_ssp_future_zm,ta_UKESM_ssp_future_zm};

for im=[1:nmodels_an]
    lat = LAT{im};
    lon = LON{im};
    lon(length(lon)+1) = 360;
     sulf_delta = ta_g6{im}-ta_ba{im};
     sulf_delta(length(lon),:) =  sulf_delta(1,:);
     ssp5_delta = ta_sp{im}-ta_ba{im};
     ssp5_delta(length(lon),:) =  ssp5_delta(1,:);
     sulf_delta_int(:,:,im) = interp2(lon',lat',sulf_delta',lon_int1',lat_int1,'linear');
     ssp5_delta_int(:,:,im) = interp2(lon',lat',ssp5_delta',lon_int1',lat_int1,'linear');
end

sulf_delta_all1 = nanmean(sulf_delta_int,3);
ssp5_delta_all1 = nanmean(ssp5_delta_int,3);
sulf_delta_std1 = nanstd(sulf_delta_int,0,3)/(nmodels_an);
ssp5_delta_std1 = nanstd(ssp5_delta_int,0,3)/(nmodels_an);

save('G615kSAI.mat','LAT','LON','ta_g6','ta_ba','ta_sp','sulf_delta_all1','ssp5_delta_all1')

%%

figure(15)

colormap(fc)
set(gcf, 'Position',  [200, 200, 1600, 1600])

s1 = subplot(3,2,1);
colormap(gca,fc)
worldmap([-90 90],[-180 180])
box on
hold on
pcolorm(lat_int1',lon_int1',ssp5_delta_all1,'LineStyle','none')
load coastlines
geoshow(coastlat,coastlon)
caxis([mmin mmax])
set(gca,'Linewidth',2)
title({'a) SSP2-4.5 - Baseline';'multi-model mean (2050-2069)'},'FontSize',20,'FontWeight','Bold')
s1.Position(1) = s1.Position(1)+.102;
hl = colorbar('FontSize',16,'YTick',(v2),'Linewidth',2,'Location','westoutside');
hl.Position(1) = hl.Position(1)-.01;
hl.Position(2) = hl.Position(2)-.2;
hl.Position(4) = hl.Position(4)*1.6;
ylabel(hl,'Surface air temperature (^{\circ}C)')

s2 = subplot(3,2,2);
colormap(gca,fca)
worldmap([-90 90],[-180 180])
box on
hold on
pcolorm(lat_int1',lon_int1',abs(ssp5_delta_std1),'LineStyle','none')
load coastlines
geoshow(coastlat,coastlon)
caxis([0 mmax2])
set(gca,'Linewidth',2)
title({'b) SSP2-4.5 - Baseline';'multi-model standard error'},'FontSize',20,'FontWeight','Bold')
s2.Position(1) = s2.Position(1)-.002;

s3=subplot(3,2,3);

colormap(gca,fc)
worldmap([-90 90],[-180 180])
box on
hold on
pcolorm(lat_int1',lon_int1',sulf_delta_all1,'LineStyle','none')
%scatterm(lat_sulf_all,lon_sulf_all,.8,'k','filled')
load coastlines
geoshow(coastlat,coastlon)
caxis([mmin mmax])
set(gca,'Linewidth',2)
title({'c) G6-1.5K-SAI - Baseline';'multi-model mean (2081-2100)'},'FontSize',20,'FontWeight','Bold')
s3.Position(1) = s3.Position(1)+.102;

s4=subplot(3,2,4);
colormap(gca,fca)
worldmap([-90 90],[-180 180])
box on
hold on
pcolorm(lat_int1',lon_int1',abs(sulf_delta_std1),'LineStyle','none')
load coastlines
geoshow(coastlat,coastlon)
caxis([0 mmax2])
set(gca,'Linewidth',2)
title({'d) G6-1.5K-SAI - Baseline';'multi-model standard error'},'FontSize',20,'FontWeight','Bold')
s4.Position(1) = s4.Position(1)-.002;
hl2 = colorbar('FontSize',16,'YTick',(v2_2),'Linewidth',2,'Location','eastoutside');
hl2.Position(2) = hl.Position(2);
hl2.Position(1) = hl2.Position(1)+0.03;
ylabel(hl2,'\sigma / n_{models} (^{\circ}C)')

set(gcf,'renderer','painters')
print(gcf,'-depsc2',['maps_mean_15k.eps'])