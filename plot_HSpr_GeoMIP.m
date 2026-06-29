clear all

load GeoMIP_pr.mat

x_e_gm = NaN(6,4,3,81);

for im=1:6
    ww = cos(LAT{im}/180*pi);
    for ie=1:4
        for in=1:3
            if size(x_e{im,ie,in},1)>0
                x_e_zm{im,ie,in} = squeeze(mean(x_e{im,ie,in},1));
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

regions = shaperead('IPCC_WGII_continental_regions.shx');
regions = shaperead('IPCC-WGI-reference-regions-v4.shx');


for im=1:6
    ww = cos(LAT{im}/180*pi);
    
    [latM,lonM] = meshgrid(LAT{im},LON{im}-180);
    AA = meshgrid(ww,LON{im});
    
    ta_1 = (x_e{im,1,1}(:,:,1:81)-x_e{im,4,1}(:,:,1:81)).*AA;
    ta_2 = (x_e{im,2,1}(:,:,1:81)-x_e{im,4,1}(:,:,1:81)).*AA;
    ta_3 = (x_e{im,3,1}(:,:,1:81)-x_e{im,4,1}(:,:,1:81)).*AA;
    
    for i=1:size(regions,1)
        fprintf([regions([i]).Name '\n'])
        ind = inpolygon(latM,lonM,regions([i]).Y,regions([i]).X);
        wS = sum(AA(ind));
        
        for ex=1:7
            for j=1:81
                fieldt = ta_1(:,:,j);
                ta_reg_sulf(im,j,i) = sum(fieldt(ind),'omitnan')/wS;
                fieldt = ta_2(:,:,j);
                ta_reg_solar(im,j,i) = sum(fieldt(ind),'omitnan')/wS;
                fieldt = ta_3(:,:,j);
                ta_reg_ssp(im,j,i) = sum(fieldt(ind),'omitnan')/wS;
            end
        end
    end
    
end

%save('pr_regions2.mat','ta_reg_sulf','ta_reg_solar')

%%

% for i=1:size(regions,1)
%     time=2020:2100;
%     figure(i)
%     subplot(2,1,1)
%     for im=1:6
%         x = ta_reg_sulf(im,:,i);
%         mvm = 10;
%
%         y = ta_reg_solar(im,:,i);
%
%
%
%     hold on
%     box on
%
%     plot(time,movmean(x,mvm),'k')
%     plot(time,movmean(y,mvm),'r')
%
%     end
%     x = ta_reg_sulf(:,:,i);
%     y = ta_reg_solar(:,:,i);
%
%     plot(time,movmean(mean(x,1),mvm),'k','Linewidth',2)
%     plot(time,movmean(mean(y,1),mvm),'r','Linewidth',2)
%     set(gca,'Linewidth',2,'FontSize',18,'FontWeight','Bold')
%     title(regions([i]).Name,'FontSize',18,'FontWeight','Bold')
%
%     for t=1:81
%     x = ta_reg_sulf(:,t,i);
%     y = ta_reg_solar(:,t,i);
%
%     ta_var_sulf(t) = std((x),0).^2;
%     ta_var_solar(t) = std((y),0).^2;
%     end
%
%     subplot(2,1,2)
%     hold on
%     box on
%
%     tot_var = ta_var_sulf+ta_var_solar;
%
%     area(time,movmean((ta_var_solar+ta_var_sulf)./tot_var*100,10),'FaceColor','k')
%     area(time,movmean(ta_var_solar./tot_var*100,10),'FaceColor','r')
%     plot(time,ones(length(time)),'w')
%     set(gca,'Linewidth',2,'FontSize',18,'FontWeight','Bold')
%
% end

%%

for i=1:size(regions,1)
    time=2020:2100;
    %     figure(i)
    %     subplot(2,1,1)
    %     for im=1:6
    %         x = ta_reg_sulf(im,:,i);
    %         mvm = 10;
    %
    %         y = ta_reg_solar(im,:,i);
    %
    %
    %
    %         hold on
    %         box on
    %
    %         plot(time,movmean(x,mvm),'k')
    %         plot(time,movmean(y,mvm),'r')
    %
    %     end
    %     x = ta_reg_sulf(:,:,i);
    %     y = ta_reg_solar(:,:,i);
    %
    %     plot(time,movmean(mean(x,1),mvm),'k','Linewidth',2)
    %     plot(time,movmean(mean(y,1),mvm),'r','Linewidth',2)
    %     set(gca,'Linewidth',2,'FontSize',18,'FontWeight','Bold')
    %     title(regions([i]).Region,'FontSize',18,'FontWeight','Bold')
    
        x = mean(ta_reg_sulf(:,71:81,i),2);
        y = mean(ta_reg_solar(:,71:81,i),2);
        z = mean(ta_reg_ssp(:,71:81,i),2);
        
        ta_var_sulf(i) = std((x),0).^2;
        ta_var_solar(i) = std((y),0).^2;
        ta_var_ssp(i) = std((z),0).^2;

    
    %     subplot(2,1,2)
    %     hold on
    %     box on
    %
    %     tot_var = ta_var_sulf+ta_var_solar;
    %     tas_r2 = (ta_var_sulf./tot_var*100);
    %
    %     area(time,movmean((ta_var_solar+ta_var_sulf)./tot_var*100,1),'FaceColor','k')
    %     area(time,movmean(ta_var_solar./tot_var*100,1),'FaceColor','r')
    %     plot(time,ones(length(time)),'w')
    %     set(gca,'Linewidth',2,'FontSize',18,'FontWeight','Bold')
    %
end

save('pr_regions1_std.mat','ta_var_sulf','ta_var_solar','ta_var_ssp')