clear all

load GeoMIP_tas.mat
x_e_tas = x_e;
load GeoMIP_pr.mat
x_e_pr = x_e;

load data_all_SOCOL.mat

%%



x_e_gm_tas = NaN(6,4,81,3);
x_e_gm_pr = NaN(6,4,81,3);

%do global mean T and P

for im=1:6
    ww = cos(LAT{im}/180*pi);
    for ie=1:4
        for in=1:3
            if size(x_e{im,ie,in},1)>0
                x_e_zm_tas{im,ie,in} = squeeze(mean(x_e_tas{im,ie,in},1));
                x_e_zm_pr{im,ie,in} = squeeze(mean(x_e_pr{im,ie,in},1));
                y_tas_zm = squeeze(sum(x_e_zm_tas{im,ie,in}.*ww)/sum(ww));
                y_pr_zm = squeeze(sum(x_e_zm_pr{im,ie,in}.*ww)/sum(ww));

                if size(y_tas_zm,2)==86
                    x_e_gm_tas(im,ie,1:81,in) = y_tas_zm(6:86);
                    x_e_gm_pr(im,ie,1:81,in) = y_pr_zm(6:86);
                else
                    x_e_gm_tas(im,ie,:,in) = y_tas_zm;
                    x_e_gm_pr(im,ie,:,in) = y_pr_zm;
                end
            end
        end
    end
end

%%
%do regional means T and P
regions = shaperead('IPCC-WGI-reference-regions-v4.shx');

sizeR = size(regions,1);

for im=1:6
    ww = cos(LAT{im}/180*pi);

    [latM,lonM] = meshgrid(LAT{im},LON{im}-180);
    AA = meshgrid(ww,LON{im});

    ta_1 = (x_e_tas{im,1,1}(:,:,1:81)-x_e_tas{im,4,1}(:,:,1:81)).*AA;
    ta_2 = (x_e_tas{im,2,1}(:,:,1:81)-x_e_tas{im,4,1}(:,:,1:81)).*AA;
    ta_3 = (x_e_tas{im,3,1}(:,:,1:81)-x_e_tas{im,4,1}(:,:,1:81)).*AA;

    for i=1:sizeR
        fprintf([regions([i]).Name '\n'])
        ind = inpolygon(latM,lonM,regions([i]).Y,regions([i]).X);
        wS = sum(AA(ind));

        for ex=1:6
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



for im=1:6
    ww = cos(LAT{im}/180*pi);

    [latM,lonM] = meshgrid(LAT{im},LON{im}-180);
    AA = meshgrid(ww,LON{im});

    pr_1 = (x_e_pr{im,1,1}(:,:,1:81)-x_e_pr{im,4,1}(:,:,1:81)).*AA;
    pr_2 = (x_e_pr{im,2,1}(:,:,1:81)-x_e_pr{im,4,1}(:,:,1:81)).*AA;
    pr_3 = (x_e_pr{im,3,1}(:,:,1:81)-x_e_pr{im,4,1}(:,:,1:81)).*AA;

    for i=1:sizeR
        fprintf([regions([i]).Name '\n'])
        ind = inpolygon(latM,lonM,regions([i]).Y,regions([i]).X);
        wS = sum(AA(ind));

        for ex=1:6
            for j=1:81
                fieldt = pr_1(:,:,j);
                pr_reg_sulf(im,j,i) = sum(fieldt(ind),'omitnan')/wS;
                fieldt = pr_2(:,:,j);
                pr_reg_solar(im,j,i) = sum(fieldt(ind),'omitnan')/wS;
                fieldt = pr_3(:,:,j);
                pr_reg_ssp(im,j,i) = sum(fieldt(ind),'omitnan')/wS;
            end
        end
    end

end
%%
%do global mean T and P for SOCOL
ww = cos(lat_SOCOL/180*pi);

x_e_gm_tasS = NaN(4,4,70,3);
x_e_gm_prS = NaN(4,4,70,3);


ww = cos(lat_SOCOL/180*pi);

[latM,lonM] = meshgrid(lat_SOCOL,lon_SOCOL-180);
AA = meshgrid(ww,lon_SOCOL);

for ie=1:4
    for in=1:3
        x_e_zm_tasS{ie,1,in} = squeeze(mean(ta_g6_y(:,:,:,in,ie),1));
        x_e_zm_prS{ie,1,in} = squeeze(mean(pr_g6_y(:,:,:,in,ie),1))*86400;
        x_e_zm_tasS{ie,3,in} = squeeze(mean(ta_ssp5_y(:,:,:,in),1));
        x_e_zm_prS{ie,3,in} = squeeze(mean(pr_ssp5_y(:,:,:,in),1))*86400;
        x_e_zm_tasS{ie,4,in} = squeeze(mean(ta_ssp2_y(:,:,:,in),1));
        x_e_zm_prS{ie,4,in} = squeeze(mean(pr_ssp2_y(:,:,:,in),1))*86400;

        y_tas_zmS = squeeze(sum(x_e_zm_tasS{ie,1,in}.*ww)/sum(ww));
        y_pr_zmS = squeeze(sum(x_e_zm_prS{ie,1,in}.*ww)/sum(ww));

        x_e_gm_tasS(ie,1,:,in) = y_tas_zmS;
        x_e_gm_prS(ie,1,:,in) = y_pr_zmS;

        y_tas_zmS = squeeze(sum(x_e_zm_tasS{ie,3,in}.*ww)/sum(ww));
        y_pr_zmS = squeeze(sum(x_e_zm_prS{ie,3,in}.*ww)/sum(ww));

        x_e_gm_tasS(ie,3,:,in) = y_tas_zmS;
        x_e_gm_prS(ie,3,:,in) = y_pr_zmS;

        y_tas_zmS = squeeze(sum(x_e_zm_tasS{ie,4,in}.*ww)/sum(ww));
        y_pr_zmS = squeeze(sum(x_e_zm_prS{ie,4,in}.*ww)/sum(ww));

        x_e_gm_tasS(ie,4,:,in) = y_tas_zmS;
        x_e_gm_prS(ie,4,:,in) = y_pr_zmS;

    end
end
%%
%do regional means T and P for SOCOL

for ie=1:4
    for in=1:3

        ta_1S = (ta_g6_y(:,:,:,in,ie)-ta_ssp2_y(:,:,:,in)).*AA;
        ta_3S = (ta_ssp5_y(:,:,:,in)-ta_ssp2_y(:,:,:,in)).*AA;
        pr_1S = (pr_g6_y(:,:,:,in,ie)-pr_ssp2_y(:,:,:,in)).*AA;
        pr_3S = (pr_ssp5_y(:,:,:,in)-pr_ssp2_y(:,:,:,in)).*AA;        

        for i=1:sizeR
            fprintf([regions([i]).Name '\n'])
            ind = inpolygon(latM,lonM,regions([i]).Y,regions([i]).X);
            wS = sum(AA(ind));

            for j=1:70
                fieldt = ta_1S(:,:,j);
                ta_reg_sulfS(ie,j,in,i) = sum(fieldt(ind),'omitnan')/wS;
                fieldt = ta_3S(:,:,j);
                ta_reg_sspS(ie,j,in,i) = sum(fieldt(ind),'omitnan')/wS;
                fieldt = pr_1S(:,:,j)*86400;
                pr_reg_sulfS(ie,j,in,i) = sum(fieldt(ind),'omitnan')/wS;
                fieldt = pr_3S(:,:,j)*86400;
                pr_reg_sspS(ie,j,in,i) = sum(fieldt(ind),'omitnan')/wS;            
            end
        end

    end
end
%%
%plot global values
time = 2020:2100;
figure(1)
subplot(2,1,1)

l_colb = 5;
fc = brewermap(l_colb,'Pastel1');
fc2 = brighten(fc,-0.3);
fc3 = brighten(fc2,-0.8);

hold on
box on

for im=1:6
    for in=1:3
        mvm = 10;
        x = squeeze(x_e_gm_tas(im,1,:,in))-squeeze(x_e_gm_tas(im,4,:,in));
        y = squeeze(x_e_gm_tas(im,2,:,in))-squeeze(x_e_gm_tas(im,4,:,in));
        z = squeeze(x_e_gm_tas(im,3,:,in))-squeeze(x_e_gm_tas(im,4,:,in));


        plot(time,movmean(x,mvm),'Color',fc(1,:))
        plot(time,movmean(y,mvm),'Color',fc(2,:))
        plot(time,movmean(z,mvm),'Color',fc(3,:))

    end

    x = squeeze(mean(x_e_gm_tas(im,1,:,:),4,'omitnan'))-squeeze(mean(x_e_gm_tas(im,4,:,:),4,'omitnan'));
    y = squeeze(mean(x_e_gm_tas(im,2,:,:),4,'omitnan'))-squeeze(mean(x_e_gm_tas(im,4,:,:),4,'omitnan'));
    z = squeeze(mean(x_e_gm_tas(im,3,:,:),4,'omitnan'))-squeeze(mean(x_e_gm_tas(im,4,:,:),4,'omitnan'));


    plot(time,movmean(x,mvm),'k','LineWidth',2,'Color',fc2(1,:))
    plot(time,movmean(y,mvm),'b','LineWidth',2,'Color',fc2(2,:))
    plot(time,movmean(z,mvm),'r','LineWidth',2,'Color',fc2(3,:))

    tas_sulf(im,:) = x;
    tas_sola(im,:) = y;
    tas_ssp5(im,:) = z;

end

x = squeeze(mean(tas_sulf,1));
y = squeeze(mean(tas_sola,1));
z = squeeze(mean(tas_ssp5,1));

lnh = plot(time,movmean((x),mvm),'k','Linewidth',4,'Color',fc3(1,:));
lnh2 = plot(time,movmean((y),mvm),'b','Linewidth',4,'Color',fc3(2,:));
plot(time,movmean((z),mvm),'r','Linewidth',4,'Color',fc3(3,:))

time2 = 2031:2100;

for im=1:4
    for in=1:3
        mvm = 10;
        x = squeeze(x_e_gm_tasS(im,1,:,in))-squeeze(x_e_gm_tasS(im,4,:,in));
        z = squeeze(x_e_gm_tasS(im,3,:,in))-squeeze(x_e_gm_tasS(im,4,:,in));


        plot(time2,movmean(x,mvm),'Color',fc(4,:))
        if im==1
            plot(time2,movmean(z,mvm),'Color',fc(3,:))
        end

    end

    x = squeeze(mean(x_e_gm_tasS(im,1,:,:),4,'omitnan'))-squeeze(mean(x_e_gm_tasS(im,4,:,:),4,'omitnan'));
    z = squeeze(mean(x_e_gm_tasS(im,3,:,:),4,'omitnan'))-squeeze(mean(x_e_gm_tasS(im,4,:,:),4,'omitnan'));


    plot(time2,movmean(x,mvm),'k','LineWidth',2,'Color',fc2(4,:))
    if im==1
        plot(time2,movmean(z,mvm),'r','LineWidth',2,'Color',fc2(3,:))
    end

    tas_sulfS(im,:) = x;
    tas_ssp5S(im,:) = z;

end

x = squeeze(mean(tas_sulfS,1));
z = squeeze(mean(tas_ssp5S,1));

uistack(lnh,'top')
uistack(lnh2,'top')
plot(time2,movmean((x),mvm),'k','Linewidth',4,'Color',fc3(4,:))


set(gca,'Linewidth',2,'FontSize',18,'FontWeight','Bold')
title('Global mean surface temperature','FontSize',18,'FontWeight','Bold','Color','k')
ylabel('K')

set(gcf, 'Color', 'w'); % Set figure background color to white
set(gca, 'Color', 'w')
set(gca,'XColor','k','YColor','k');



subplot(2,1,2)

hold on
box on

for im=1:6
    for in=1:3
        mvm = 10;
        x = squeeze(x_e_gm_pr(im,1,:,in))-squeeze(x_e_gm_pr(im,4,:,in));
        y = squeeze(x_e_gm_pr(im,2,:,in))-squeeze(x_e_gm_pr(im,4,:,in));
        z = squeeze(x_e_gm_pr(im,3,:,in))-squeeze(x_e_gm_pr(im,4,:,in));


        plot(time,movmean(x,mvm),'Color',fc(1,:))
        plot(time,movmean(y,mvm),'Color',fc(2,:))
        plot(time,movmean(z,mvm),'Color',fc(3,:))

    end

    x = squeeze(mean(x_e_gm_pr(im,1,:,:),4,'omitnan'))-squeeze(mean(x_e_gm_pr(im,4,:,:),4,'omitnan'));
    y = squeeze(mean(x_e_gm_pr(im,2,:,:),4,'omitnan'))-squeeze(mean(x_e_gm_pr(im,4,:,:),4,'omitnan'));
    z = squeeze(mean(x_e_gm_pr(im,3,:,:),4,'omitnan'))-squeeze(mean(x_e_gm_pr(im,4,:,:),4,'omitnan'));


    plot(time,movmean(x,mvm),'k','LineWidth',2,'Color',fc2(1,:))
    plot(time,movmean(y,mvm),'b','LineWidth',2,'Color',fc2(2,:))
    plot(time,movmean(z,mvm),'r','LineWidth',2,'Color',fc2(3,:))
    pr_sulf(im,:) = x;
    pr_sola(im,:) = y;
    pr_ssp5(im,:) = z;

end

x = squeeze(mean(pr_sulf,1));
y = squeeze(mean(pr_sola,1));
z = squeeze(mean(pr_ssp5,1));

lnh = plot(time,movmean((x),mvm),'k','Linewidth',4,'Color',fc3(1,:));
lnh2 = plot(time,movmean((y),mvm),'b','Linewidth',4,'Color',fc3(2,:));
plot(time,movmean((z),mvm),'r','Linewidth',4,'Color',fc3(3,:))

time2 = 2031:2100;

for im=1:4
    for in=1:3
        mvm = 10;
        x = squeeze(x_e_gm_prS(im,1,:,in))-squeeze(x_e_gm_prS(im,4,:,in));
        z = squeeze(x_e_gm_prS(im,3,:,in))-squeeze(x_e_gm_prS(im,4,:,in));


        plot(time2,movmean(x,mvm),'Color',fc(4,:))
        if im==1
            plot(time2,movmean(z,mvm),'Color',fc(3,:))
        end

    end

    x = squeeze(mean(x_e_gm_prS(im,1,:,:),4,'omitnan'))-squeeze(mean(x_e_gm_prS(im,4,:,:),4,'omitnan'));
    z = squeeze(mean(x_e_gm_prS(im,3,:,:),4,'omitnan'))-squeeze(mean(x_e_gm_prS(im,4,:,:),4,'omitnan'));


    plot(time2,movmean(x,mvm),'k','LineWidth',2,'Color',fc2(4,:))
    if im==1
        plot(time2,movmean(z,mvm),'r','LineWidth',2,'Color',fc2(3,:))
    end

    pr_sulfS(im,:) = x;
    pr_ssp5S(im,:) = z;

end

x = squeeze(mean(pr_sulfS,1));
z = squeeze(mean(pr_ssp5S,1));

uistack(lnh,'top')
uistack(lnh2,'top')
plot(time2,movmean((x),mvm),'k','Linewidth',4,'Color',fc3(4,:))


set(gca,'Linewidth',2,'FontSize',18,'FontWeight','Bold')
title('Global mean precipitation','FontSize',18,'FontWeight','Bold','Color','k')
ylabel('mm/day')

set(gcf, 'Color', 'w'); % Set figure background color to white
set(gca, 'Color', 'w')
set(gca,'XColor','k','YColor','k');


%
box on
hold on

regions = shaperead('IPCC_WGII_continental_regions.shx');
regions = shaperead('IPCC-WGI-reference-regions-v4.shx');

%%
%
% %calculate regions for other GeoMIP
%
% for im=1:6
%     ww = cos(LAT{im}/180*pi);
%
%     [latM,lonM] = meshgrid(LAT{im},LON{im}-180);
%     AA = meshgrid(ww,LON{im});
%
%     ta_1 = (x_e{im,1,1}(:,:,1:81)-x_e{im,4,1}(:,:,1:81)).*AA;
%     ta_2 = (x_e{im,2,1}(:,:,1:81)-x_e{im,4,1}(:,:,1:81)).*AA;
%     ta_3 = (x_e{im,3,1}(:,:,1:81)-x_e{im,4,1}(:,:,1:81)).*AA;
%
%     for i=1:sizeR
%         fprintf([regions([i]).Name '\n'])
%         ind = inpolygon(latM,lonM,regions([i]).Y,regions([i]).X);
%         wS = sum(AA(ind));
%
%         for j=1:81
%             fieldt = ta_1(:,:,j);
%             ta_reg_sulf(im,j,i) = sum(fieldt(ind),'omitnan')/wS;
%             fieldt = ta_2(:,:,j);
%             ta_reg_solar(im,j,i) = sum(fieldt(ind),'omitnan')/wS;
%             fieldt = ta_3(:,:,j);
%             ta_reg_ssp(im,j,i) = sum(fieldt(ind),'omitnan')/wS;
%         end
%     end
%
% end


%%
fcs2 = brewermap(7,'Dark2');
fcs = brighten(fcs2,0.2);

figure(1)
set(gcf, 'Position',  [200, 200, 1600,800])
for i=[22:28]
    time=2020:2100;
    
    subplot(4,2,i-21)
    for im=1:6
        x = ta_reg_sulf(im,:,i);
        mvm = 10;

        y = ta_reg_solar(im,:,i);
        z = ta_reg_ssp(im,:,i);



        hold on
        box on

        plot(time,movmean(x,mvm),'Color',fcs(1,:))
        plot(time,movmean(y,mvm),'Color',fcs(2,:))
        plot(time,movmean(z,mvm),'Color',fcs(3,:))

    end
    x = ta_reg_sulf(:,:,i);
    y = ta_reg_solar(:,:,i);
    z = ta_reg_ssp(:,:,i);

    p(2)=plot(time,movmean(mean(x,1),mvm),'Color',fcs2(1,:),'Linewidth',2);
    p(3)=plot(time,movmean(mean(y,1),mvm),'Color',fcs2(2,:),'Linewidth',2);
    p(1)=plot(time,movmean(mean(z,1),mvm),'Color',fcs2(3,:),'Linewidth',2);
    set(gca,'Linewidth',2,'FontSize',14,'FontWeight','Bold')
    title(regions([i]).Name,'FontSize',14,'FontWeight','Bold','Color','k')

    x = mean(ta_reg_sulf(:,71:81,i),2);
    y = mean(ta_reg_solar(:,71:81,i),2);
    z = mean(ta_reg_ssp(:,71:81,i),2);

    ta_var_sulf(i) = std((x),0).^2;
    ta_var_solar(i) = std((y),0).^2;
    ta_var_ssp(i) = std((z),0).^2;

    errorbar(2101,mean(x),sqrt(ta_var_sulf(i)),'Color',fcs2(1,:),'Linewidth',2)
    errorbar(2102,mean(y),sqrt(ta_var_solar(i)),'Color',fcs2(2,:),'Linewidth',2)
    errorbar(2102,mean(z),sqrt(ta_var_ssp(i)),'Color',fcs2(3,:),'Linewidth',2)



    for im=1:4
        x = squeeze(mean(ta_reg_sulfS(im,:,:,i),3));
        y = squeeze(mean(ta_reg_sspS(im,:,:,i),3));
        mvm = 10;

        plot(time2,movmean(x,mvm),'Color',fcs(4,:))

    end
    x = squeeze(mean(mean(ta_reg_sulfS(:,:,:,i),3),1));
    y = squeeze(mean(ta_reg_sspS(:,:,:,i),3));


    p(4)=plot(time2,movmean(x,mvm),'Color',fcs2(4,:),'Linewidth',2);
    %plot(time2,movmean(mean(x,1),mvm),'b','Linewidth',2)

    ylabel('K')

    set(gcf, 'Color', 'w'); % Set figure background color to white
    set(gca, 'Color', 'w')
    set(gca,'XColor','k','YColor','k');

    k = mean(mean(ta_reg_sulfS(:,60:70,:,i),2),3);
    ta_var_sulfS(i) = std((k),0).^2;
    errorbar(2103,mean(k),sqrt(ta_var_sulfS(i)),'Color',fcs2(4,:),'Linewidth',2)
    axis([2030 2105 -1 2.5])

    if i==28
        hl=legend(p,{'SSP5-8.5','G6sulfur','G6solar','SOCOL'},'Location','westoutside','Color','w','TextColor','k','Location','west');
        hl.Position(1)=hl.Position(1)+.4;
    end

end

set(gcf,'renderer','painters')
print(gcf,'-depsc2',['GeoMIP_uncertainties_aerosols_SOCOL_tas_Africa.eps'])
%%
figure(2)
set(gcf, 'Position',  [200, 200, 1600,800])
for i=[22:28]
    time=2020:2100;
    
    subplot(4,2,i-21)
    for im=1:6
        x = pr_reg_sulf(im,:,i);
        mvm = 10;

        y = pr_reg_solar(im,:,i);
        z = pr_reg_ssp(im,:,i);



        hold on
        box on

        plot(time,movmean(x,mvm),'Color',fcs(1,:))
        plot(time,movmean(y,mvm),'Color',fcs(2,:))
        plot(time,movmean(z,mvm),'Color',fcs(3,:))

    end
    x = pr_reg_sulf(:,:,i);
    y = pr_reg_solar(:,:,i);
    z = pr_reg_ssp(:,:,i);

    p(2)=plot(time,movmean(mean(x,1),mvm),'Color',fcs2(1,:),'Linewidth',2);
    p(3)=plot(time,movmean(mean(y,1),mvm),'Color',fcs2(2,:),'Linewidth',2);
    p(1)=plot(time,movmean(mean(z,1),mvm),'Color',fcs2(3,:),'Linewidth',2);
    set(gca,'Linewidth',2,'FontSize',14,'FontWeight','Bold')
    title(regions([i]).Name,'FontSize',14,'FontWeight','Bold','Color','k')

    x = mean(pr_reg_sulf(:,71:81,i),2);
    y = mean(pr_reg_solar(:,71:81,i),2);
    z = mean(pr_reg_ssp(:,71:81,i),2);

    pr_var_sulf(i) = std((x),0).^2;
    pr_var_solar(i) = std((y),0).^2;
    pr_var_ssp(i) = std((z),0).^2;

    errorbar(2101,mean(x),sqrt(pr_var_sulf(i)),'Color',fcs2(1,:),'Linewidth',2)
    errorbar(2102,mean(y),sqrt(pr_var_solar(i)),'Color',fcs2(2,:),'Linewidth',2)
    errorbar(2104,mean(z),sqrt(pr_var_ssp(i)),'Color',fcs2(3,:),'Linewidth',2)



    for im=1:4
        x = squeeze(mean(pr_reg_sulfS(im,:,:,i),3));
        y = squeeze(mean(pr_reg_sspS(im,:,:,i),3));
        mvm = 10;

        plot(time2,movmean(x,mvm),'Color',fcs(4,:))

    end
    x = squeeze(mean(mean(pr_reg_sulfS(:,:,:,i),3),1));
    y = squeeze(mean(pr_reg_sspS(:,:,:,i),3));


    p(4)=plot(time2,movmean(x,mvm),'Color',fcs2(4,:),'Linewidth',2);
    %plot(time2,movmean(mean(x,1),mvm),'b','Linewidth',2)

    ylabel('mm/day')

    set(gcf, 'Color', 'w'); % Set figure background color to white
    set(gca, 'Color', 'w')
    set(gca,'XColor','k','YColor','k');

    k = mean(mean(pr_reg_sulfS(:,60:70,:,i),2),3);
    pr_var_sulfS(i) = std((k),0).^2;
    errorbar(2103,mean(k),sqrt(pr_var_sulfS(i)),'Color',fcs2(4,:),'Linewidth',2)
    axis([2030 2105 -2.5 2.5])

    if i==28
        hl=legend(p,{'SSP5-8.5','G6sulfur','G6solar','SOCOL'},'Location','westoutside','Color','w','TextColor','k','Location','west');
        hl.Position(1)=hl.Position(1)+.4;
    end

end

set(gcf,'renderer','painters')
print(gcf,'-depsc2',['GeoMIP_uncertainties_aerosols_SOCOL_pr_Africa.eps'])
%%
fcs = brewermap(l_colb,'Set1');
col_regs = {[1 45:46],[2:6],[7:9],[10:16],[17:20],[22:28],[29:39],[40:44],[47 58]};
spa_regs = {[45 46],[2 6],[7 9],[10 16],[17 20],[22 28],[29 39],[40 44],[47 58]};
spa_regs2 = {[44.5 46.5],[1.5 6.5],[6.5 9.5],[9.5 16.5],[16.5 20.5],[21.5 28.5],[28.5 39.5],[39.5 44.5],[46.5 58.5]};
mid_regs = [3,8,13,18,26,34,42,45,53];
mid_regN = {'NORTH AMERICA','CENTRAL AMERICA','SOUTH AMERICA','EUROPE','AFRICA','ASIA','OCEANIA','POLAR','OCEANS'};
figure(99)
set(gcf, 'Position',  [200, 200, 1600,800])
hold on
box on
tot_var = ta_var_sulf+ta_var_solar+ta_var_sulfS;

totsd_S = (ta_var_solar+ta_var_sulf)./tot_var*100;
%[totsd_S,I] = sort(totsd_S);

totsd_SS = (ta_var_solar)./tot_var*100;
%totsd_SS = (totsd_SS(I));

p(1)=area(1:sizeR,tot_var./tot_var*100,'FaceColor',fcs(5,:),'EdgeColor','none');
p(2)=area(1:sizeR,totsd_S,'FaceColor',fcs(3,:),'EdgeColor','none');
p(3)=area(1:sizeR,totsd_SS,'FaceColor',fcs(2,:),'EdgeColor','none');
for iR=1:size(col_regs,2)
plot([spa_regs2{iR}],[median(totsd_S(col_regs{iR})) median(totsd_S(col_regs{iR}))],'LineWidth',2,'Color','k')
plot([spa_regs2{iR}],[median(totsd_SS(col_regs{iR})) median(totsd_SS(col_regs{iR}))],'LineWidth',2,'Color','k')
end

axis([1 sizeR 0 100])
%area(time,(ta_var_solar./tot_var*100),'FaceColor','r')

set(gcf, 'Color', 'w'); % Set figure background color to white
set(gca, 'Color', 'w')
set(gca,'XColor','k','YColor','k');
set(gca,'Linewidth',2,'FontSize',18,'FontWeight','Bold')
title('Driver of temperature uncertainty in GeoMIP per IPCC region','FontSize',18,'FontWeight','Bold','Color','k')
ylabel('Fractional uncertainty')
xticks(mid_regs); xticklabels(mid_regN)
legend(p(:),{'Aerosol type (single model)','Aerosol impacts between models','Models climate response'},'Color','w','TextColor','k','Location','northwestoutside')

set(gcf,'renderer','painters')
print(gcf,'-depsc2',['GeoMIP_uncertainties_aerosols_SOCOL_tas.eps'])
%%

for i=1:sizeR
    time=2020:2100;
%    figure(i+1)
%    subplot(2,1,1)
    for im=1:6
        x = pr_reg_sulf(im,:,i);
        mvm = 10;

        y = pr_reg_solar(im,:,i);
        z = pr_reg_ssp(im,:,i);



        hold on
        box on

        plot(time,movmean(x,mvm),'k')
        plot(time,movmean(y,mvm),'r')
        %plot(time,movmean(z,mvm),'b')

    end
    x = pr_reg_sulf(:,:,i);
    y = pr_reg_solar(:,:,i);
    z = pr_reg_ssp(:,:,i);

    plot(time,movmean(mean(x,1),mvm),'k','Linewidth',2)
    plot(time,movmean(mean(y,1),mvm),'r','Linewidth',2)
    %plot(time,movmean(mean(z,1),mvm),'b','Linewidth',2)
    set(gca,'Linewidth',2,'FontSize',18,'FontWeight','Bold')
    title(regions([i]).Name,'FontSize',18,'FontWeight','Bold','Color','k')

    x = mean(pr_reg_sulf(:,71:81,i),2);
    y = mean(pr_reg_solar(:,71:81,i),2);
    z = mean(pr_reg_ssp(:,71:81,i),2);

    pr_var_sulf(i) = std((x),0).^2;
    pr_var_solar(i) = std((y),0).^2;
    pr_var_ssp(i) = std((z),0).^2;


    for im=1:4
        x = squeeze(mean(pr_reg_sulfS(im,:,:,i),3));
        y = squeeze(mean(pr_reg_sspS(im,:,:,i),3));
        mvm = 10;

        plot(time2,movmean(x,mvm),'b')

    end
    x = squeeze(mean(mean(pr_reg_sulfS(:,:,:,i),3),1));
    y = squeeze(mean(pr_reg_sspS(:,:,:,i),3));


    plot(time2,movmean(x,mvm),'b','Linewidth',2)
    %plot(time2,movmean(mean(x,1),mvm),'b','Linewidth',2)

    ylabel('K')

    set(gcf, 'Color', 'w'); % Set figure background color to white
    set(gca, 'Color', 'w')
    set(gca,'XColor','k','YColor','k');

    k = mean(mean(pr_reg_sulfS(:,60:70,:,i),2),3);
    pr_var_sulfS(i) = std((k),0).^2;
    %

end
%%
fcs = brewermap(l_colb,'Set1');

figure(99)

subplot(1,2,1)
hold on
box on
tot_var = ta_var_sulf+ta_var_solar+ta_var_sulfS;

totsd_S = (ta_var_solar+ta_var_sulf)./tot_var*100;
%[totsd_S,I] = sort(totsd_S);

totsd_SS = (ta_var_solar)./tot_var*100;
%totsd_SS = (totsd_SS(I));

p(1)=area(1:sizeR,tot_var./tot_var*100,'FaceColor',fcs(5,:),'EdgeColor','none');
p(2)=area(1:sizeR,totsd_S,'FaceColor',fcs(3,:),'EdgeColor','none');
p(3)=area(1:sizeR,totsd_SS,'FaceColor',fcs(2,:),'EdgeColor','none');
for iR=1:size(col_regs,2)
plot([spa_regs2{iR}],[median(totsd_S(col_regs{iR})) median(totsd_S(col_regs{iR}))],'LineWidth',2,'Color','k')
plot([spa_regs2{iR}],[median(totsd_SS(col_regs{iR})) median(totsd_SS(col_regs{iR}))],'LineWidth',2,'Color','k')
end

axis([1 sizeR 0 100])
%area(time,(ta_var_solar./tot_var*100),'FaceColor','r')

set(gcf, 'Color', 'w'); % Set figure background color to white
set(gca, 'Color', 'w')
set(gca,'XColor','k','YColor','k');
set(gca,'Linewidth',2,'FontSize',16,'FontWeight','Bold')
title('Driver of surface temperature uncertainty in GeoMIP simulations','FontSize',18,'FontWeight','Bold','Color','k')
ylabel('Fractional uncertainty')
xticks(mid_regs); xticklabels(mid_regN)
legend(p(:),{'Aerosol type (single model)','Aerosol impacts between models','Models climate response'},'Color','w','TextColor','k','Location','west')


set(gcf, 'Position',  [200, 200, 1600,800])
subplot(1,2,2)
hold on
box on
tot_var = pr_var_sulf+pr_var_solar+pr_var_sulfS;

totsd_S = (pr_var_solar+pr_var_sulf)./tot_var*100;
%[totsd_S,I] = sort(totsd_S);

totsd_SS = (pr_var_solar)./tot_var*100;
%totsd_SS = (totsd_SS(I));

p(1)=area(1:sizeR,tot_var./tot_var*100,'FaceColor',fcs(5,:),'EdgeColor','none');
p(2)=area(1:sizeR,totsd_S,'FaceColor',fcs(3,:),'EdgeColor','none');
p(3)=area(1:sizeR,totsd_SS,'FaceColor',fcs(2,:),'EdgeColor','none');
for iR=1:size(col_regs,2)
plot([spa_regs2{iR}],[median(totsd_S(col_regs{iR})) median(totsd_S(col_regs{iR}))],'LineWidth',2,'Color','k')
plot([spa_regs2{iR}],[median(totsd_SS(col_regs{iR})) median(totsd_SS(col_regs{iR}))],'LineWidth',2,'Color','k')
end

axis([1 sizeR 0 100])
%area(time,(ta_var_solar./tot_var*100),'FaceColor','r')

set(gcf, 'Color', 'w'); % Set figure background color to white
set(gca, 'Color', 'w')
set(gca,'XColor','k','YColor','k');
set(gca,'Linewidth',2,'FontSize',16,'FontWeight','Bold')
title('Driver of precipitation uncertainty in GeoMIP simulations','FontSize',18,'FontWeight','Bold','Color','k')
ylabel('Fractional uncertainty')
xticks(mid_regs); xticklabels(mid_regN)

%legend(p(:),{'Aerosol type (in one model)','Aerosol impacts between models','Internal models variability'},'Color','w','TextColor','k','Location','west')

set(gcf,'renderer','painters')
print(gcf,'-depsc2',['GeoMIP_uncertainties_aerosols_SOCOL_pr.eps'])