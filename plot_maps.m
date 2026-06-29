clear all
load ../GeoMIP_G6Sulfur/map_plots_GeoMIP.mat %use saves from read_data_model
load map_plots_GeoMIP_SOCOL.mat %use saves from read_data_model

LAT{7} = lat_SOCOL;
LON{7} = lon_SOCOL;
for im=1:4
x_map_20{im,7} = x_map_20_SOCOL{im};
end

n_exp = {'G6sulfur','G6solar','ssp585','ssp245'};
n_mod = {'CNRM-ESM2-1','IPSL-CM6A-LR','CESM2-WACCM','UKESM1-0-LL','MPI-ESM1-2-LR','MPI-ESM1-2-HR','SOCOL'};

match = zeros(5,2);
match(1,:) = [true,true];
match(2,:) = [true,true];
match(3,:) = [true,true];
match(4,:) = [true,true];
match(5,:) = [true,true];
match(6,:) = [true,true];
match(7,:) = [true,true];

l_colb = 11;
fc = brewermap(l_colb,'*RdBu');
mmin=-2.2;
mmax=2.2;
v = mmin:(mmax-mmin)/l_colb:mmax;
v2 = mmin:(mmax-mmin)/(l_colb):mmax;

plot_single = false;

for im=[1:6]
    
    if plot_single
    
    figure(im)
    colormap(fc)
    set(gcf, 'Position',  [200, 200,1200,600])
    
    end
    
    if match(im,1)
        hs = subplot(1,2,1);
        lat = double(LAT{im});
        lon = double(LON{im});
        
        
        
        sulf_delta = x_map_20{1,im} - x_map_20{4,im};
        %         if im==7
        %             sulf_delta = fliplr(x_map_20{1,im}) - x_map_20{4,im};
        %         end
        %%%calculate T gradients%%%
        L0 = ones(length(lat),1);
        L1 = sind(lat);
        L2 = (1.5*L1.^2-.5);
        clat = cos(lat/180*pi)'; clatsum = sum(clat); slat = sin(lat/180*pi)';
        sulf_deltazm = mean(sulf_delta,1);
        T0_sulf{im} = sum(sulf_deltazm.*clat,2)/clatsum;
        T1_sulf{im} = sum(sulf_deltazm.*clat.*slat,2)/clatsum;
        T2_sulf{im} = sum(sulf_deltazm.*clat.*(slat.^2-1)/2,2)/clatsum;
        
        
        if lon(1)>=0
            lon2=lon; lon(1)=0; lon(2:length(lon2)+1)=lon2;
            sulf_delta2=sulf_delta; sulf_delta(2:size(sulf_delta,1)+1,:)=sulf_delta2;
        end
        lon(length(lon)+1) = 360;
        
        sulf_delta(length(lon),:) =  sulf_delta(1,:);
        
        if plot_single
            worldmap([-90 90],[-180 180])
            box on
            hold on
            contourfm(lat,lon,sulf_delta',v,'LineStyle','none')
            load coastlines
            geoshow(coastlat,coastlon)
            scatterm(lat_sulf{im},lon_sulf{im},.4,'k','filled')
            caxis([mmin mmax])
            %colorbar('FontSize',16,'YTick',(v2),'Linewidth',2,'Location','southoutside');
            set(gca,'Linewidth',2)
            title([n_mod{im} ' ' n_exp{1} '-' n_exp{4}],'FontSize',18,'FontWeight','Bold')
            hs.Position(2)=hs.Position(2)+.05;
            hs.Position(1)=hs.Position(1)+.05;
        end
    end
    
    if match(im,2)
        hs2 = subplot(1,2,2);
        lat = double(LAT{im});
        lon = double(LON{im});
        
        sulf_delta = x_map_20{2,im} - x_map_20{4,im};
        if im==7
           % sulf_delta = fliplr(x_map_20{2,im}) - x_map_20{4,im};
        end
        if lon(1)>0
            lon2=lon; lon(1)=0; lon(2:length(lon2)+1)=lon2;
            sulf_delta2=sulf_delta; sulf_delta(2:size(sulf_delta,1)+1,:)=sulf_delta2;
        end
        
        %%%calculate T gradients%%%
        L0 = ones(length(lat),1);
        L1 = sind(lat);
        L2 = (1.5*L1.^2-.5);
        clat = cos(lat/180*pi)'; clatsum = sum(clat); slat = sin(lat/180*pi)';
        sola_deltazm = mean(sulf_delta,1);
        T0_sola{im} = sum(sola_deltazm.*clat,2)/clatsum;
        T1_sola{im} = sum(sola_deltazm.*clat.*slat,2)/clatsum;
        T2_sola{im} = sum(sola_deltazm.*clat.*(slat.^2-1)/2,2)/clatsum;
        
        lon(length(lon)+1) = 360;
        
        sulf_delta(length(lon),:) =  sulf_delta(1,:);
        
        if plot_single
            worldmap([-90 90],[-180 180])
            box on
            hold on
            contourfm(lat,lon,sulf_delta',v,'LineStyle','none')
            load coastlines
            geoshow(coastlat,coastlon)
            scatterm(lat_sola{im},lon_sola{im},.4,'k','filled')
            caxis([mmin mmax])
            %colorbar('FontSize',16,'YTick',(v2),'Linewidth',2,'Location','southoutside');
            set(gca,'Linewidth',2)
            title([n_mod{im} ' ' n_exp{2} '-' n_exp{4}],'FontSize',18,'FontWeight','Bold')
            hs2.Position(2)=hs2.Position(2)+.05;
            hs2.Position(1)=hs2.Position(1)-.05;
        end
    end
    
    if plot_single
        set(gcf,'renderer','painters')
        print(gcf,'-depsc2',['figures/singlemodels' n_mod{im} '_maps.eps'])
    end
end


%%


lat_int1 = -90:2:90;
lon_int1 = 0:2:360;

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

sulf_agree = zeros(size(sulf_delta_all1)); 

for im=1:nmodels_an
    if match(im,1)==1
        field1 = squeeze(sulf_delta_int(:,:,im));
        for i=1:length(lat_int1)
            for j=1:length(lon_int1)
                if sign(field1(i,j))==sign(sulf_delta_all1(i,j))
                    sulf_agree(i,j) = sulf_agree(i,j)+1;
                end
            end
        end
    end
end

cl = 1;lat_sulf_all = 0; lon_sulf_all = 0;
cs = 1;lat_solar_all = 0; lon_solar_all = 0;
for i=1:length(lat_int1)
    for j=1:length(lon_int1)
        if sulf_agree(i,j)<ceil(nmodels_an*.66)
            cl=cl+1; lat_sulf_all(cl) = lat_int1(i); lon_sulf_all(cl) = lon_int1(j);
        end
    end
end

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
title({'a) SSP5-8.5-SSP2-4.5';'multi-model mean (2081-2100)'},'FontSize',20,'FontWeight','Bold')
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
title({'b) SSP5-8.5-SSP2-4.5';'multi-model standard error'},'FontSize',20,'FontWeight','Bold')
s2.Position(1) = s2.Position(1)-.002;

s3=subplot(3,2,3);

colormap(gca,fc)
worldmap([-90 90],[-180 180])
box on
hold on
pcolorm(lat_int1',lon_int1',sulf_delta_all1,'LineStyle','none')
scatterm(lat_sulf_all,lon_sulf_all,.8,'k','filled')
load coastlines
geoshow(coastlat,coastlon)
caxis([mmin mmax])
set(gca,'Linewidth',2)
title({'c) G6sulfur-SSP2-4.5';'multi-model mean (2081-2100)'},'FontSize',20,'FontWeight','Bold')
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
title({'d) G6sulfur-SSP2-4.5';'multi-model standard error'},'FontSize',20,'FontWeight','Bold')
s4.Position(1) = s4.Position(1)-.002;
hl2 = colorbar('FontSize',16,'YTick',(v2_2),'Linewidth',2,'Location','eastoutside');
hl2.Position(2) = hl.Position(2);
hl2.Position(1) = hl2.Position(1)+0.03;
ylabel(hl2,'\sigma / n_{models} (^{\circ}C)')

set(gcf,'renderer','painters')
print(gcf,'-depsc2',['maps_mean.eps'])
%%