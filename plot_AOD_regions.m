exp = {'feedback_30N_30S','feedback_15N_15S','feedback_0N_1.0','feedback_60N_60S_1.0','SSP245-MA-GAUSS-DEFAULT','SSP245-MA-GAUSS-LOWER-0.5','SSP245-MA-GAUSS-LOWER-1.0'};
exp2 = {'1.0 above PI - 30-30 ','1.0 above PI - 15-15','1.0 above PI - EQ','1.0 above PI - 60-60','1.5 above PI - MS','1.0 above PI - MS','0.5 above PI - MS','SSP2-4.5'};

ww = cos(lat/180*pi);
ta_zm = squeeze(mean(ta,1)); 
for ex=1:7
    for en=1:nums(ex)
        ta_gm(:,ex,en) = squeeze(sum(ta_zm(:,:,ex,en).*ww)/sum(ww));
    end
end


fact = 1;
ta_gm_y = squeeze(mean(reshape(ta_gm,12,35,7,3),1))*fact;

time = 2035:2069;

l_colb = 7;
fc = brewermap(l_colb,'Set2');

figure(99)
set(gcf, 'Position',  [200, 200, 1200,800])
subplot(3,5,[1 2])

box on
hold on

for ex=1:7
    for en=1:nums(ex)
        plot(time,movmean(squeeze(ta_gm_y(:,ex,en)),10),'Color',fc(ex,:),'Linewidth',.1);
    end
    p(ex) = plot(time,mean(movmean(squeeze(ta_gm_y(:,ex,:)),10),2),'Color',fc(ex,:),'Linewidth',4);
end

plot([2035 2070],[0 0],'Color','k','Linewidth',1);
title({'Global mean stratospheric AOD compared'},'FontSize',18,'FontWeight','Bold')
hl=legend(p(:),exp2,'Location','northeastoutside');
hl.Position(2)=hl.Position(2)-0.04;
hl.Position(1)=hl.Position(1)+0.14;
set(gca,'Linewidth',2,'FontSize',18,'FontWeight','Bold')
axis([2035 2069 0 0.4])
ylabel('K')


ta_gm_tot = reshape(reshape(movmean(ta_gm_y,10),35,21),35,21);
ta_gm_int = squeeze(ta_gm_y(:,5,:));
ta_gm_stra = squeeze(ta_gm_y(:,[1:4 6],1));
ta_gm_sce = squeeze(ta_gm_y(:,5:7,1));

ta_gm_tot_var = std(ta_gm_tot,0,2).^2;

ta_gm_int_var = std(ta_gm_int,0,2).^2;

ta_gm_stra_var = std(ta_gm_stra,0,2).^2;

ta_gm_sce_var = std(ta_gm_sce,0,2).^2;

tot_var = ta_gm_int_var+ta_gm_stra_var+ta_gm_sce_var;

subplot(3,2,[3 5])

hold on
box on



area(time,movmean((ta_gm_stra_var+ta_gm_int_var+ta_gm_sce_var)./tot_var*100,10))
area(time,movmean((ta_gm_stra_var+ta_gm_int_var)./tot_var*100,10))
area(time,movmean(ta_gm_int_var./tot_var*100,10))

hl=legend('Scenario (how much cooling?)','Strategy (where do we inject to get that cooling?)','Internal (contribution of internal variability)','Location','north');
hl.Position(1)=hl.Position(1)+0.04;
hl.Position(2)=hl.Position(2)-0.004;
set(gca,'Linewidth',2,'FontSize',18,'FontWeight','Bold')
axis([2035 2069 0 100])

title({'Global mean AOD uncertainty partition'},'FontSize',18,'FontWeight','Bold')

set(gcf,'renderer','painters')
print(gcf,'-depsc2',['HS_plot_map_AOD.eps'])

%%

regions = shaperead('IPCC_WGII_continental_regions.shx');

[latM,lonM] = meshgrid(lat,lon-180);
AA = ncread('AREA.nc','AREA');
AA = squeeze(mean(AA,3));
ta_aa = ta.*AA/sum(sum(AA));

for i=1:8
    fprintf([regions([i]).Region '\n'])
    ind = inpolygon(latM,lonM,regions([i]).Y,regions([i]).X);
    %Area_nat(i) = nansum(AA(ind))*1e-6;
        for ex=1:7
            for en=1:nums(ex)
                for j=1:420
                    fieldt = ta(:,:,j,ex,en);
                    %ta_reg(j,ex,en,i) = sum(fieldt(ind).*AA(ind))/Area_nat(i);
                    ta_reg(j,ex,en,i) = nanmean(fieldt(ind));
                end
            end
        end
end



%%

l_colb = 7;
fc = brewermap(l_colb,'Set2');

for f=1:8

    ta_reg1 = ta_reg(:,:,:,f)*1;
    
ta_reg_y = reshape(ta_reg1,12,35,7,3);
ta_reg_y = squeeze(mean(ta_reg_y(1:12,:,:,:),1));

figure(f)
set(gcf, 'Position',  [200, 200, 1200,800])
subplot(3,5,[1 2])



ta_reg_tot = reshape(reshape(movmean(ta_reg_y,10),35,21),35,21);
ta_reg_int = squeeze(ta_reg_y(:,5,:));
ta_reg_stra = squeeze(ta_reg_y(:,[1:4 6],1));
ta_reg_sce = squeeze(ta_reg_y(:,5:7,1));

ta_reg_tot_var = std(ta_reg_tot,0,2).^2;

ta_reg_int_var = std(ta_reg_int,0,2).^2;

ta_reg_stra_var = std(ta_reg_stra,0,2).^2;

ta_reg_sce_var = std(ta_reg_sce,0,2).^2;

ta_reg_tot_var_all(:,f) = ta_reg_tot_var;
ta_reg_stra_var_all(:,f) = ta_reg_stra_var;
ta_reg_sce_var_all(:,f) = ta_reg_sce_var;
ta_reg_int_var_all(:,f) = ta_reg_int_var;

h1 = subplot(2,2,1);

box on
hold on

for ex=1:7
    for en=1:nums(ex)
        plot(time,movmean(squeeze(ta_reg_y(:,ex,en)),10),'Color',fc(ex,:),'Linewidth',.1);
    end
    p(ex) = plot(time,mean(movmean(squeeze(ta_reg_y(:,ex,:)),10),2),'Color',fc(ex,:),'Linewidth',4);
end

plot([2035 2070],[0 0],'Color','k','Linewidth',1);
title(regions([f]).Region,'FontSize',18,'FontWeight','Bold')
set(gca,'Linewidth',2,'FontSize',18,'FontWeight','Bold')
axis([2035 2069 0 0.4])
set(gca,'XTick',[])

tot_var_2 = ta_reg_int_var+ta_reg_stra_var+ta_reg_sce_var;

h2 = subplot(2,2,3);
h2.Position(2) = h2.Position(2) + 0.12;

hold on
box on



area(time,movmean((ta_reg_stra_var+ta_reg_int_var+ta_reg_sce_var)./tot_var_2*100,10))
area(time,movmean((ta_reg_stra_var+ta_reg_int_var)./tot_var_2*100,10))
area(time,movmean(ta_reg_int_var./tot_var_2*100,10))

save('AODVIS.mat','ta_reg_int_var_all','ta_reg_sce_var_all','ta_reg_stra_var_all','ta_reg_tot_var_all','tot_var','ta_gm_int_var','ta_gm_stra_var','ta_gm_sce_var')

set(gca,'Linewidth',2,'FontSize',18,'FontWeight','Bold')
axis([2035 2069 0 100])
%title({'EU JJA precipitation uncertainty partition'},'FontSize',18,'FontWeight','Bold')

set(gcf,'renderer','painters')
print(gcf,'-depsc2',['HS_plot_AODVIS_' regions([f]).Region '.eps'])

end
