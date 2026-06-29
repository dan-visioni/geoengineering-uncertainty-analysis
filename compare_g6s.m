clear all
load ../GeoMIP_G6Sulfur/map_plots_GeoMIP.mat %use saves from read_data_model
load map_plots_GeoMIP_SOCOL.mat %use saves from read_data_model
LAT{7} = lat_SOCOL;
LON{7} = lon_SOCOL;
for im=1:4
x_map_20{im,7} = x_map_20_SOCOL{im};
end
lat_int1 = -90:2:90;
lon_int1 = 0:2:360;

match(1,:) = [true,true];
match(2,:) = [true,true];
match(3,:) = [true,true];
match(4,:) = [true,true];
match(5,:) = [true,true];
match(6,:) = [true,true];
match(7,:) = [true,true];

solar_delta_all1 = nan(length(lat_int1),length(lon_int1),sum(match(:,1),1));
sulf_delta_all1 = nan(length(lat_int1),length(lon_int1),sum(match(:,2),1));
ssp5_delta_all1 = nan(length(lat_int1),length(lon_int1),sum(match(:,2),1));
solar_delta_int = nan(length(lat_int1),length(lon_int1),sum(match(:,1),1));
sulf_delta_int = nan(length(lat_int1),length(lon_int1),sum(match(:,2),1));
ssp5_delta_int = nan(length(lat_int1),length(lon_int1),sum(match(:,2),1));
% regridding everything on the same grid and averaging
c=1;

nmodels_an = 7;

for im=[1:nmodels_an]
    lat = LAT{im};
    lon = LON{im};
    lon(length(lon)+1) = 360;
    if match(im,1)==1
        sulf_delta = x_map_20{1,im} - x_map_20{4,im};
        sulf_delta(length(lon),:) =  sulf_delta(1,:);
        x = interp2(lon',lat,sulf_delta',lon_int1',lat_int1,'linear',NaN);
        sulf_delta_int(:,:,c) = x;
    end
    
    ssp5_delta = x_map_20{3,im} - x_map_20{4,im};
    ssp5_delta(length(lon),:) =  ssp5_delta(1,:);
    x = interp2(lon',lat',ssp5_delta',lon_int1',lat_int1,'linear');
    ssp5_delta_int(:,:,c) = x;
    
    c=c+1;
end


sulf_delta_all1 = nanmean(sulf_delta_int,3);
ssp5_delta_all1 = nanmean(ssp5_delta_int,3);
sulf_delta_std1 = nanstd(sulf_delta_int,0,3)/sqrt(nmodels_an);
ssp5_delta_std1 = nanstd(ssp5_delta_int,0,3)/sqrt(nmodels_an);

lat_int1 = -90:2:90;
lon_int1 = 0:2:360;

nmodels_an = 3;

load G615kSAI.mat
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

sulf_delta_15k = nanmean(sulf_delta_int,3);
ssp5_delta_15k = nanmean(ssp5_delta_int,3);

l_colb = 22;
fc = brewermap(l_colb,'*RdBu');
mmin=-1.65;
mmax=1.65;
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


figure(10)

colormap(fc)
set(gcf, 'Position',  [200, 200, 1600, 1600])

s3=subplot(3,2,1);

colormap(gca,fc)
worldmap([-90 90],[-180 180])
box on
hold on
pcolorm(lat_int1',lon_int1',sulf_delta_all1,'LineStyle','none')
load coastlines
geoshow(coastlat,coastlon)
caxis([mmin mmax])
set(gca,'Linewidth',2)
title({'c) G6sulfur-SSP2-4.5';'multi-model mean (2081-2100)'},'FontSize',20,'FontWeight','Bold')
s3.Position(1) = s3.Position(1)+.102;
hl = colorbar('FontSize',16,'YTick',(v2),'Linewidth',2,'Location','westoutside');
hl.Position(1) = hl.Position(1)-.01;
hl.Position(2) = hl.Position(2)-.2;
hl.Position(4) = hl.Position(4)*1.6;
ylabel(hl,'Surface air temperature (^{\circ}C)')


s3=subplot(3,2,3);

colormap(gca,fc)
worldmap([-90 90],[-180 180])
box on
hold on
pcolorm(lat_int1',lon_int1',sulf_delta_15k,'LineStyle','none')
%scatterm(lat_sulf_all,lon_sulf_all,.8,'k','filled')
load coastlines
geoshow(coastlat,coastlon)
caxis([mmin mmax])
set(gca,'Linewidth',2)
title({'c) G6-1.5K-SAI - Baseline';'multi-model mean (2066-2085)'},'FontSize',20,'FontWeight','Bold')
s3.Position(1) = s3.Position(1)+.102;

s3=subplot(3,2,2);

colormap(gca,fc)
worldmap([-90 90],[-180 180])
box on
hold on
pcolorm(lat_int1',lon_int1',sulf_delta_15k-sulf_delta_all1,'LineStyle','none')
%scatterm(lat_sulf_all,lon_sulf_all,.8,'k','filled')
load coastlines
geoshow(coastlat,coastlon)
caxis([mmin mmax])
set(gca,'Linewidth',2)
title({'c) G6-1.5K-SAI - G6sulfur';'changes'},'FontSize',20,'FontWeight','Bold')
s3.Position(2) = s3.Position(2)-.202;

set(gcf,'renderer','painters')
print(gcf,'-depsc2',['maps_mean_diff.eps'])
