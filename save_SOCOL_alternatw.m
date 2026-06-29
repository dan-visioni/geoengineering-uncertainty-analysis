clear all

fold  = '/Volumes/D_Visioni2/SOCOL/GEOMIP_SOCOL/';
ncdisp([fold 'GEOMIP_socol4_al1_surface.nc'])

ta_g6(:,:,:,1,2) = ncread([fold 'GEOMIP_socol4_al1_surface.nc'],'var167');
ta_g6(:,:,:,2,2) = ncread([fold 'GEOMIP_socol4_al2_surface.nc'],'var167');
ta_g6(:,:,:,3,2) = ncread([fold 'GEOMIP_socol4_al3_surface.nc'],'var167');

ta_g6(:,:,:,1,3) = ncread([fold 'GEOMIP_socol4_ca1_surface.nc'],'var167');
ta_g6(:,:,:,2,3) = ncread([fold 'GEOMIP_socol4_ca2_surface.nc'],'var167');
ta_g6(:,:,:,3,3) = ncread([fold 'GEOMIP_socol4_ca3_surface.nc'],'var167');

ta_g6(:,:,:,1,4) = ncread([fold 'GEOMIP_socol4_di1_surface.nc'],'var167');
ta_g6(:,:,:,2,4) = ncread([fold 'GEOMIP_socol4_di2_surface.nc'],'var167');
ta_g6(:,:,:,3,4) = ncread([fold 'GEOMIP_socol4_di3_surface.nc'],'var167');

ta_g6(:,:,:,1,1) = ncread([fold 'GEOMIP_socol4_g61_surface.nc'],'var167');
ta_g6(:,:,:,2,1) = ncread([fold 'GEOMIP_socol4_g62_surface.nc'],'var167');
ta_g6(:,:,:,3,1) = ncread([fold 'GEOMIP_socol4_g63_surface.nc'],'var167');

ta_ssp2(:,:,:,1) = ncread([fold 'GEOMIP_socol4_ssp21_surface.nc'],'var167');
ta_ssp2(:,:,:,2) = ncread([fold 'GEOMIP_socol4_ssp22_surface.nc'],'var167');
ta_ssp2(:,:,:,3) = ncread([fold 'GEOMIP_socol4_ssp23_surface.nc'],'var167');

ta_ssp5(:,:,:,1) = ncread([fold 'GEOMIP_socol4_ssp51_surface.nc'],'var167');
ta_ssp5(:,:,:,2) = ncread([fold 'GEOMIP_socol4_ssp52_surface.nc'],'var167');
ta_ssp5(:,:,:,3) = ncread([fold 'GEOMIP_socol4_ssp53_surface.nc'],'var167');

%%

pr_g6(:,:,:,1,2) = ncread([fold 'GEOMIP_socol4_al1_surface.nc'],'var142')+ncread([fold 'GEOMIP_socol4_al1_surface.nc'],'var143');
pr_g6(:,:,:,2,2) = ncread([fold 'GEOMIP_socol4_al2_surface.nc'],'var142')+ncread([fold 'GEOMIP_socol4_al1_surface.nc'],'var143');
pr_g6(:,:,:,3,2) = ncread([fold 'GEOMIP_socol4_al3_surface.nc'],'var142')+ncread([fold 'GEOMIP_socol4_al1_surface.nc'],'var143');

pr_g6(:,:,:,1,3) = ncread([fold 'GEOMIP_socol4_ca1_surface.nc'],'var142')+ncread([fold 'GEOMIP_socol4_ca1_surface.nc'],'var143');
pr_g6(:,:,:,2,3) = ncread([fold 'GEOMIP_socol4_ca2_surface.nc'],'var142')+ncread([fold 'GEOMIP_socol4_ca2_surface.nc'],'var143');
pr_g6(:,:,:,3,3) = ncread([fold 'GEOMIP_socol4_ca3_surface.nc'],'var142')+ncread([fold 'GEOMIP_socol4_ca3_surface.nc'],'var143');

pr_g6(:,:,:,1,4) = ncread([fold 'GEOMIP_socol4_di1_surface.nc'],'var142')+ncread([fold 'GEOMIP_socol4_di1_surface.nc'],'var143');
pr_g6(:,:,:,2,4) = ncread([fold 'GEOMIP_socol4_di2_surface.nc'],'var142')+ncread([fold 'GEOMIP_socol4_di2_surface.nc'],'var143');
pr_g6(:,:,:,3,4) = ncread([fold 'GEOMIP_socol4_di3_surface.nc'],'var142')+ncread([fold 'GEOMIP_socol4_di3_surface.nc'],'var143');

pr_g6(:,:,:,1,1) = ncread([fold 'GEOMIP_socol4_g61_surface.nc'],'var142')+ncread([fold 'GEOMIP_socol4_g61_surface.nc'],'var143');
pr_g6(:,:,:,2,1) = ncread([fold 'GEOMIP_socol4_g62_surface.nc'],'var142')+ncread([fold 'GEOMIP_socol4_g62_surface.nc'],'var143');
pr_g6(:,:,:,3,1) = ncread([fold 'GEOMIP_socol4_g63_surface.nc'],'var142')+ncread([fold 'GEOMIP_socol4_g63_surface.nc'],'var143');

pr_ssp2(:,:,:,1) = ncread([fold 'GEOMIP_socol4_ssp21_surface.nc'],'var142')+ncread([fold 'GEOMIP_socol4_ssp21_surface.nc'],'var143');
pr_ssp2(:,:,:,2) = ncread([fold 'GEOMIP_socol4_ssp22_surface.nc'],'var142')+ncread([fold 'GEOMIP_socol4_ssp22_surface.nc'],'var143');
pr_ssp2(:,:,:,3) = ncread([fold 'GEOMIP_socol4_ssp23_surface.nc'],'var142')+ncread([fold 'GEOMIP_socol4_ssp23_surface.nc'],'var143');

pr_ssp5(:,:,:,1) = ncread([fold 'GEOMIP_socol4_ssp51_surface.nc'],'var142')+ncread([fold 'GEOMIP_socol4_ssp51_surface.nc'],'var143');
pr_ssp5(:,:,:,2) = ncread([fold 'GEOMIP_socol4_ssp52_surface.nc'],'var142')+ncread([fold 'GEOMIP_socol4_ssp52_surface.nc'],'var143');
pr_ssp5(:,:,:,3) = ncread([fold 'GEOMIP_socol4_ssp53_surface.nc'],'var142')+ncread([fold 'GEOMIP_socol4_ssp53_surface.nc'],'var143');

%%

lat = ncread([fold 'GEOMIP_socol4_ssp21_surface.nc'],'lat');
lat_SOCOL = lat;

lon = ncread([fold 'GEOMIP_socol4_ssp21_surface.nc'],'lon');
lon_SOCOL = lon;

pr_g6_y = squeeze(mean(reshape(pr_g6,192,96,12,70,3,4),3));
ta_g6_y = squeeze(mean(reshape(ta_g6,192,96,12,70,3,4),3));

pr_ssp5_y = squeeze(mean(reshape(pr_ssp5,192,96,12,70,3),3));
ta_ssp5_y = squeeze(mean(reshape(ta_ssp5,192,96,12,70,3),3));

pr_ssp2_y = squeeze(mean(reshape(pr_ssp2,192,96,12,70,3),3));
ta_ssp2_y = squeeze(mean(reshape(ta_ssp2,192,96,12,70,3),3));

save('data_all_SOCOL.mat','pr_g6_y',"ta_g6_y","ta_ssp5_y","ta_ssp2_y","pr_ssp5_y","pr_ssp2_y","lat_SOCOL","lon_SOCOL")