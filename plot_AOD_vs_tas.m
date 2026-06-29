clear all

fc1 = brewermap(7,'Set2');
fc2 = brighten(fc1,.8);

load  ../GeoMIP_G6Sulfur/AOD_GeoMIP.mat

figure(1)
set(gcf, 'Position',  [200, 200, 600,700])

p(1) = plot(lat,x2_y_CNRM,'Color',fc1(1,:),'Linewidth',2);
hold on
p(2) = plot(LAT{2},x_lat{2,1}-x_lat{2,3}*0,'Color',fc1(2,:),'Linewidth',2);
p(3) = plot(LAT{3},x_lat{3,1}-x_lat{3,3},'Color',fc1(3,:),'Linewidth',2);
p(4) = plot(LAT{4},x_lat{4,1}-x_lat{4,3},'Color',fc1(4,:),'Linewidth',2);
p(5) = plot(lat2,x2_y_MPI,'Color',fc1(5,:),'Linewidth',2);
p(6) = plot(lat3,x2_y_MPIH,'Color',fc1(6,:),'Linewidth',2);
set(gca,'Linewidth',2,'FontSize',16,'FontWeight','Bold')
title('a) 2090-2100 zonal mean stratospheric AOD')
set(gca,'XTick',[-90:30:90],'YTick',[0:.05:1],'FontSize',18,'FontWeight','Bold')
axis([-90 90 0 .9])

lat_int1 = -90:2:90;

aod_all = zeros(size(lat_int1,2),6);

aod_all(:,1) = interp1(lat,x2_y_CNRM,lat_int1,'linear',NaN);
aod_all(:,2) = interp1(LAT{2},x_lat{2,1}-x_lat{2,3}*0,lat_int1,'linear',NaN);
aod_all(:,3) = interp1(LAT{3},x_lat{3,1}-x_lat{3,3},lat_int1,'linear',NaN);
aod_all(:,4) = interp1(LAT{4},x_lat{4,1}-x_lat{4,3},lat_int1,'linear',NaN);
aod_all(:,5) = interp1(lat2,x2_y_MPI,lat_int1,'linear',NaN);
aod_all(:,6) = interp1(lat3,x2_y_MPIH,lat_int1,'linear',NaN);

aod_all_std = std(aod_all,0,2);

load 2deg_tas.mat

figure(3)

fc1 = brewermap(81,'Spectral');
fc2 = brighten(fc1,.8);

hold on
box on

for i=1:91
    scatter(sulf_delta_std1(i,:),fc2(i,:))
end