clear all

load GeoMIP_tas.mat
load borderdata.mat

x_e_gm = NaN(6,4,3,81);

for im=1:6
    ww = cos(LAT{im}/180*pi);
    for ie=1:4
        for in=1:3
            if size(x_e{im,ie,in},1)>0
                x_e_zm{im,ie,in} = squeeze(mean(x_e{im,ie,in},1,'omitnan'));
                y = squeeze(sum(x_e_zm{im,ie,in}.*ww)/sum(ww));
                if size(y,2)==86
                    x_e_gm(im,ie,in,1:81) = y(6:86);
                else
                    x_e_gm(im,ie,in,:) = y;
                end
            end
        end
    end
end

%%
box on
hold on

%regions = shaperead('IPCC_WGII_continental_regions.shx');
regions = shaperead('IPCC-WGI-reference-regions-v4.shx');


for im=1:6
    ww = cos(LAT{im}/180*pi);
    
    [latM,lonM] = meshgrid(LAT{im},LON{im}-180);
    AA = meshgrid(ww,LON{im});
    
    ta_1 = (x_e{im,1,1}(:,:,1:81)-x_e{im,4,1}(:,:,1:81)).*AA;
    ta_2 = (x_e{im,2,1}(:,:,1:81)-x_e{im,4,1}(:,:,1:81)).*AA;
    
    for i=1:size(regions,1)
        fprintf([regions(i).Name '\n'])
        ind = inpolygon(latM,lonM,regions(i).Y,regions(i).X);
        wS = sum(AA(ind));
        
        for j=1:81
            fieldt = ta_1(:,:,j);
            ta_reg_sulf(im,j,i) = sum(fieldt(ind),'omitnan')/wS;
            fieldt = ta_2(:,:,j);
            ta_reg_solar(im,j,i) = sum(fieldt(ind),'omitnan')/wS;
        end
    end
    
end

for i=1:size(regions,1)

        x = mean(ta_reg_sulf(:,71:81,i),2,'omitnan');
        y = mean(ta_reg_solar(:,71:81,i),2,'omitnan');
        
        ta_var_sulf(i) = std((x),0).^2;
        ta_var_solar(i) = std((y),0).^2;

end

save('tas_regions1_std.mat','ta_var_sulf','ta_var_solar')
