clear all

fold  = '/Volumes/D_Visioni2/GeoMIP/'; %Change with directory where you downloaded data
n_exp = {'G6sulfur','G6solar','ssp585','ssp245'}; %Experiments name

n_mod = {'SOCOL'};

ta_ssp2(:,:,:,1) = ncread([fold n_exp{1} '/' n_mod{1} '/socol4_fin1_2m_temp_mm_2015-2099.nc'],'var167');
ta_ssp2(:,:,:,2) = ncread([fold n_exp{1} '/' n_mod{1} '/socol4_fin2_2m_temp_mm_2015-2099.nc'],'var167');

ta_g6(:,:,:,1) = ncread([fold n_exp{1} '/' n_mod{1} '/socol4_g61_2m_temp_mm_2030-2099.nc'],'var167');
ta_g6(:,:,:,2) = ncread([fold n_exp{1} '/' n_mod{1} '/socol4_g62_2m_temp_mm_2030-2099.nc'],'var167');
ta_g6(:,:,:,3) = ncread([fold n_exp{1} '/' n_mod{1} '/socol4_g63_2m_temp_mm_2030-2099.nc'],'var167');

ta_ssp5(:,:,:,1) = ncread([fold n_exp{1} '/' n_mod{1} '/socol4_s851_2m_temp_mm_2015-2099.nc'],'var167');
ta_ssp5(:,:,:,2) = ncread([fold n_exp{1} '/' n_mod{1} '/socol4_s852_2m_temp_mm_2015-2099.nc'],'var167');
ta_ssp5(:,:,:,3) = ncread([fold n_exp{1} '/' n_mod{1} '/socol4_s853_2m_temp_mm_2015-2099.nc'],'var167');

lat = ncread([fold n_exp{1} '/' n_mod{1} '/socol4_s851_2m_temp_mm_2015-2099.nc'],'lat');
lat_SOCOL = lat;

lon = ncread([fold n_exp{1} '/' n_mod{1} '/socol4_s851_2m_temp_mm_2015-2099.nc'],'lon');
lon_SOCOL = lon;

%%

ww = cos(lat/180*pi);

ta_ssp2y = squeeze(sum(squeeze(mean(ta_ssp2,1)).*ww))/sum(ww); ta_ssp2y = squeeze(mean(reshape(ta_ssp2y,12,size(ta_ssp2y,1)/12,2),1));
ta_g6y = squeeze(sum(squeeze(mean(ta_g6,1)).*ww))/sum(ww); ta_g6y = squeeze(mean(reshape(ta_g6y,12,size(ta_g6y,1)/12,3),1));
ta_ssp5y = squeeze(sum(squeeze(mean(ta_ssp5,1)).*ww))/sum(ww); ta_ssp5y = squeeze(mean(reshape(ta_ssp5y,12,size(ta_ssp5y,1)/12,3),1));

time1 = 2015:2099;
time2 = 2030:2099;

figure(1)
box on
hold on

plot(time1,ta_ssp2y,'b')
plot(time1,ta_ssp5y,'r')
plot(time2,ta_g6y,'k')

x_glob_m_save_SOCOL = NaN(4,length(ta_ssp2y));

x_glob_m_save_SOCOL(4,:) = mean(ta_ssp2y,2);
x_glob_m_save_SOCOL(3,:) = mean(ta_ssp5y,2);
x_glob_m_save_SOCOL(1,16:85) = mean(ta_g6y,2);

save(['tas_timeseries_SOCOL.mat'],'x_glob_m_save_SOCOL')
figure(2)

l_colb = 11;
fc = brewermap(l_colb,'*RdBu');
mmin=-2.2;
mmax=2.2;
v = mmin:(mmax-mmin)/l_colb:mmax;
v2 = mmin:(mmax-mmin)/(l_colb):mmax;

lat_SOCOL = flip(lat_SOCOL);
lat = flip(lat);

colormap(fc)
worldmap([-90 90],[-180 180])
box on
hold on
contourfm(lat,lon,fliplr(double(squeeze(mean(mean(ta_g6(:,:,601:840,:),3),4))-squeeze(mean(mean(ta_ssp2(:,:,781:1020,:),3),4))))',v,'LineStyle','none')
load coastlines
geoshow(coastlat,coastlon)
caxis([mmin mmax])
colorbar('FontSize',16,'YTick',(v2),'Linewidth',2,'Location','southoutside');
set(gca,'Linewidth',2)

x_map_20_SOCOL{1} = fliplr(squeeze(mean(mean(ta_g6(:,:,601:840,:),3),4)));
x_map_20_SOCOL{4} = fliplr(squeeze(mean(mean(ta_ssp2(:,:,781:1020,:),3),4)));
x_map_20_SOCOL{3} = fliplr(squeeze(mean(mean(ta_ssp5(:,:,781:1020,:),3),4)));

save(['map_plots_GeoMIP_SOCOL.mat'],'x_map_20_SOCOL','lat_SOCOL','lon_SOCOL')

%%

aod_g6(:,:,:,1) = ncread([fold n_exp{1} '/' n_mod{1} '/SOCOL_AOD_440-625_g61_2030-2099.nc'],'od525str');
aod_g6(:,:,:,2) = ncread([fold n_exp{1} '/' n_mod{1} '/SOCOL_AOD_440-625_g62_2030-2099.nc'],'od525str');
aod_g6(:,:,:,3) = ncread([fold n_exp{1} '/' n_mod{1} '/SOCOL_AOD_440-625_g63_2030-2099.nc'],'od525str');

aod_ssp5(:,:,:,1) = ncread([fold n_exp{1} '/' n_mod{1} '/SOCOL_AOD_440-625_ssp51_2030-2099.nc'],'od525str');
aod_ssp5(:,:,:,2) = ncread([fold n_exp{1} '/' n_mod{1} '/SOCOL_AOD_440-625_ssp51_2030-2099.nc'],'od525str');
aod_ssp5(:,:,:,3) = ncread([fold n_exp{1} '/' n_mod{1} '/SOCOL_AOD_440-625_ssp51_2030-2099.nc'],'od525str');
%%

lat = ncread([fold n_exp{1} '/' n_mod{1} '/socol4_s851_2m_temp_mm_2015-2099.nc'],'lat'); ww = cos(lat/180*pi);
aod_g6s = squeeze(sum(squeeze(mean(aod_g6,1)).*ww))/sum(ww); aod_g6s = squeeze(mean(reshape(aod_g6s,12,size(aod_g6s,1)/12,3),1));
aod_ssp5s = squeeze(sum(squeeze(mean(aod_ssp5,1)).*ww))/sum(ww); aod_ssp5s = squeeze(mean(reshape(aod_ssp5s,12,size(aod_ssp5s,1)/12,3),1));

aod_SOCOL = mean(aod_g6s,2)-mean(aod_ssp5s,2);

plot(aod_SOCOL)

save(['aod_SOCOL.mat'],'aod_SOCOL')