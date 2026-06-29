clear all
regions = shaperead('IPCC-WGI-reference-regions-v4.shx');

load tas_regions1_std.mat

tas_r1_sol = ta_var_solar;
tas_r1_sulf = ta_var_sulf;

load tas_regions2_std.mat

tas_r2_sol = ta_var_solar;
tas_r2_sulf = ta_var_sulf;

load tas_regions3_std.mat

tas_r3_sol = ta_var_solar;
tas_r3_sulf = ta_var_sulf;

load pr_regions1_std.mat

pr_r1_sol = ta_var_solar;
pr_r1_sulf = ta_var_sulf;

load pr_regions2_std.mat

pr_r2_sol = ta_var_solar;
pr_r2_sulf = ta_var_sulf;

load pr_regions3_std.mat

pr_r3_sol = ta_var_solar;
pr_r3_sulf = ta_var_sulf;

pr3 = (pr_r3_sulf./(pr_r3_sol+pr_r3_sulf)*100);
tas3 = (tas_r3_sulf./(tas_r3_sol+tas_r3_sulf)*100);

pr2 = (pr_r2_sulf./(pr_r2_sol+pr_r2_sulf)*100);
tas2 = (tas_r2_sulf./(tas_r2_sol+tas_r2_sulf)*100);

pr1 = (pr_r1_sulf./(pr_r1_sol+pr_r1_sulf)*100);
tas1 = (tas_r1_sulf./(tas_r1_sol+tas_r1_sulf)*100);
%%
figure(1)
set(gcf, 'Position',  [200, 200,700,600])

l_colb = 8;
fc = brewermap(l_colb,'Dark2');

col_regs = {[1 45:46],[2:6],[7:9],[10:16],[17:20],[22:28],[29:39],[40:44]};
%pr1(24)=NaN;
regions_n = {'Oceans','Polar','North America','Central America','South America','Europe','Africa','Asia','Oceania','Average'};
%scatter(tas2,pr2,250,'filled','MarkerFaceColor','k')
axis([0 100 0 100])
hold on
p(1) = plot(tas1,pr1,'o','MarkerFaceColor','b','MarkerEdgeColor','none','MarkerSize',10);
for i=1:8
    plot(tas1(col_regs{i}),pr1(col_regs{i}),'o','MarkerFaceColor',fc(i,:),'MarkerEdgeColor','none','MarkerSize',10);
    p(i+1) = plot(mean(tas1(col_regs{i}),'omitnan'),mean(pr1(col_regs{i}),'omitnan'),'o','MarkerFaceColor',fc(i,:),'MarkerEdgeColor','k','MarkerSize',15);
end

p(10) = plot(mean(tas1,'omitnan'),mean(pr1,'omitnan'),'o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',25);
%plot(mean(tas1),mean(pr1),'o','MarkerFaceColor','none','MarkerEdgeColor','k','MarkerSize',105);

%scatter(tas3,pr3,50,'filled','MarkerFaceColor','b')
box on
plot([50,50],[0,100],'k')
plot([0,100],[50,50],'k')


[hl]=legend(p(:),regions_n,'Location','southwest');
hl.FontSize = 15 ; hl.LineWidth = 2; hl.FontWeight = 'bold';

M = findobj(hl,'type','line'); % Find objects of type 'patch'
set(M,'MarkerSize',350) %Calculate marker size based on size of scatter points
%set(gca,'visible','off')
xlabel('Surface temperature fractional variance (%)','Color','k');
ylabel('Precipitation fractional variance (%)','Color','k');
title({'Fraction of modeling variance';'driven by aerosols uncertainties in GeoMIP'},'Color','k')
set(gca,'Linewidth',2,'FontSize',18,'FontWeight','Bold')

set(gcf, 'Color', 'w'); % Set figure background color to white
set(gca, 'Color', 'w')
set(hl, 'Color', 'w','TextColor','k')
set(gca,'XColor','k','YColor','k');

(mean(tas1,'omitnan')-50)/(100-mean(tas1,'omitnan'))*100
(mean(pr1,'omitnan')-50)/(100-mean(pr1,'omitnan'))*100

for i=2:8
   xx(i) =  mean(tas1(col_regs{i}));
t(i) = (mean(tas1(col_regs{i}))-50)/(100-mean(tas1(col_regs{i})))*100;
pr(i) = mean(pr1(col_regs{i})-50)/(100-mean(pr1(col_regs{i})))*100;
end

set(gcf,'renderer','painters')
print(gcf,'-depsc2',['GeoMIP_uncertainties_aerosols.eps'])