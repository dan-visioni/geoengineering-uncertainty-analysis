clear all

fold  = '/Volumes/D_Visioni2/CESM2_MA/MA_GAUSS/';

nam = 'b.e21.BWSSP245.f09_g17.release-cesm2.1.3.WACCM-MA-1deg.';

exp = {'feedback_30N_30S','feedback_15N_15S','feedback_0N_1.0','feedback_60N_60S_1.0','SSP245-MA-GAUSS-DEFAULT','SSP245-MA-GAUSS-LOWER-0.5','SSP245-MA-GAUSS-LOWER-1.0'};
exp2 = {'1.0 above PI - 30-30 ','1.0 above PI - 15-15','1.0 above PI - EQ','1.0 above PI - 60-60','1.5 above PI - MS','1.0 above PI - MS','0.5 above PI - MS'};

id = '.cam.h0.';

var = 'AODVISstdn';

time = '.203501-206912.nc';

ens = {'001','002','003'};

nums = [3,3,3,3,3,3,3];

for ex=1:7
        for en=1:3
            list = dir([fold nam exp{ex} '.' ens{en} id var '*.nc']);
            t = ncread([list(1).folder '/' list(1).name],var);
            ta(:,:,:,ex,en) = t(:,:,1:420);
        end
end
%%

lat = ncread([list(1).folder '/' list(1).name],'lat');
lon = ncread([list(1).folder '/' list(1).name],'lon');

ww = cos(lat/180*pi);

Lm = 1; LM = 49;

ta_zm = squeeze(mean(ta,1)); 
ta_zm_reg = squeeze(mean(ta(Lm:LM,:,:,:,:),1)); 

lm = 129; lM = 171;

for ex=1:7
    for en=1:nums(ex)
        ta_gm(:,ex,en) = squeeze(sum(ta_zm(:,:,ex,en).*ww)/sum(ww));
        if ex==7
            tab_gm(:,en) = squeeze(sum(tab_zm(:,:,en).*ww)/sum(ww)); 
        end
    end
end

%%

ta_gm_y = squeeze(mean(reshape(ta_gm,12,35,7,3),1));
tab_gm_y = squeeze(mean(reshape(tab_gm,12,50,3),1));
ta_reg_y = reshape(ta_reg,12,35,7,3);
tab_reg_y = reshape(tab_reg,12,50,3);
ta_reg_y = squeeze(mean(ta_reg_y(1:3,:,:,:),1));
tab_reg_y = squeeze(mean(tab_reg_y(1:3,:,:),1));

time = 2035:2069;

l_colb = 7;
fc = brewermap(l_colb,'Set2');
figure(1)
set(gcf, 'Position',  [200, 200, 1200,800])
subplot(3,5,[1 2])

box on
hold on

for ex=1:7
    for en=1:nums(ex)
        plot(time,(movmean(squeeze(ta_gm_y(:,ex,en))-mean(tab_gm_y(6:25,en),1),10)),'Color',fc(ex,:),'Linewidth',.1);
    end
    p(ex) = plot(time,mean(movmean(squeeze(ta_gm_y(:,ex,:))-mean(tab_gm_y(6:25,:),1),10),2),'Color',fc(ex,:),'Linewidth',4);
end
plot([2035 2070],[0 0],'Color','k','Linewidth',1);
title({'Global mean cooling compared';'to reference period [2020-2039]'},'FontSize',18,'FontWeight','Bold')
hl=legend(p(:),exp2,'Location','northeastoutside');
hl.Position(2)=hl.Position(2)-0.04;
hl.Position(1)=hl.Position(1)+0.14;
set(gca,'Linewidth',2,'FontSize',18,'FontWeight','Bold')
axis([2035 2069 -2 1])

subplot(3,5,[4 5])

box on
hold on

for ex=1:7
    for en=1:nums(ex)
        plot(time,(movmean(squeeze(ta_reg_y(:,ex,en))-mean(tab_reg_y(6:25,en),1),10)),'Color',fc(ex,:),'Linewidth',.1);
    end
    p(ex) = plot(time,mean(movmean(squeeze(ta_reg_y(:,ex,:))-mean(tab_reg_y(6:25,:),1),10),2),'Color',fc(ex,:),'Linewidth',4);
end
plot([2035 2070],[0 0],'Color','k','Linewidth',1);
title({'EU DJF mean cooling compared';'to reference period [2020-2039]'},'FontSize',18,'FontWeight','Bold')
set(gca,'Linewidth',2,'FontSize',18,'FontWeight','Bold')
axis([2035 2069 -2 1])


ta_gm_tot = reshape(reshape(movmean(ta_gm_y,10),35,21),35,21);
ta_gm_int = squeeze(ta_gm_y(:,5,:));
ta_gm_stra = squeeze(ta_gm_y(:,[1:4 6],1));
ta_gm_sce = squeeze(ta_gm_y(:,5:7,1));

ta_gm_tot_var = std(ta_gm_tot,0,2).^2;

ta_gm_int_var = std(ta_gm_int,0,2).^2;

ta_gm_stra_var = std(ta_gm_stra,0,2).^2;

ta_gm_sce_var = std(ta_gm_sce,0,2).^2;

ta_reg_tot = reshape(reshape(movmean(ta_reg_y,10),35,21),35,21);
ta_reg_int = squeeze(ta_reg_y(:,5,:));
ta_reg_stra = squeeze(ta_reg_y(:,[1:4 6],1));
ta_reg_sce = squeeze(ta_reg_y(:,5:7,1));

ta_reg_tot_var = std(ta_reg_tot,0,2).^2;

ta_reg_int_var = std(ta_reg_int,0,2).^2;

ta_reg_stra_var = std(ta_reg_stra,0,2).^2;

ta_reg_sce_var = std(ta_reg_sce,0,2).^2;


hold on
box on

% plot(ta_gm_tot_var)
% plot(ta_gm_int_var)
% plot(ta_gm_stra_var)
% plot(ta_gm_sce_var)
% plot(ta_gm_int_var+ta_gm_stra_var+ta_gm_sce_var,'k')

tot_var_2 = ta_gm_int_var+ta_gm_stra_var+ta_gm_sce_var;

subplot(3,2,[3 5])

hold on
box on



area(time,movmean((ta_gm_stra_var+ta_gm_int_var+ta_gm_sce_var)./tot_var_2*100,10))
area(time,movmean((ta_gm_stra_var+ta_gm_int_var)./tot_var_2*100,10))
area(time,movmean(ta_gm_int_var./tot_var_2*100,10))

hl=legend('Scenario (how much cooling?)','Strategy (where do we inject to get that cooling?)','Internal (contribution of internal variability)','Location','north');
hl.Position(1)=hl.Position(1)+0.04;
hl.Position(2)=hl.Position(2)-0.004;
set(gca,'Linewidth',2,'FontSize',18,'FontWeight','Bold')
axis([2035 2069 0 100])

title({'Global mean temperature uncertainty partition'},'FontSize',18,'FontWeight','Bold')

tot_var_2 = ta_reg_int_var+ta_reg_stra_var+ta_reg_sce_var;

subplot(3,2,[4 6])

hold on
box on



area(time,movmean((ta_reg_stra_var+ta_reg_int_var+ta_reg_sce_var)./tot_var_2*100,10))
area(time,movmean((ta_reg_stra_var+ta_reg_int_var)./tot_var_2*100,10))
area(time,movmean(ta_reg_int_var./tot_var_2*100,10))

set(gca,'Linewidth',2,'FontSize',18,'FontWeight','Bold')
axis([2035 2069 0 100])
title({'EU DJF temperature uncertainty partition'},'FontSize',18,'FontWeight','Bold')

set(gcf,'renderer','painters')
print(gcf,'-depsc2',['HS_plot_TREFHT.eps'])