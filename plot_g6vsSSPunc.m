clear all
regions = shaperead('IPCC-WGI-reference-regions-v4.shx');

load tas_regions1_std.mat

tas_r1_sol = ta_var_ssp;
tas_r1_sulf = ta_var_sulf;


load pr_regions1_std.mat

pr_r1_sol = ta_var_ssp;
pr_r1_sulf = ta_var_sulf;


pr1 = (pr_r1_sol./(pr_r1_sol+pr_r1_sulf)*100);
tas1 = (tas_r1_sol./(tas_r1_sol+tas_r1_sulf)*100);
%%
figure(1)
set(gcf, 'Position',  [200, 200,700,600])

l_colb = 8;
fc = brewermap(l_colb,'Dark2');

col_regs = {[1 45:46],[2:6],[7:9],[10:16],[17:20],[22:28],[29:39],[40:44]};
regions_n = {'Oceans','Polar','North America','Central America','South America','Europe','Africa','Asia','Oceania','Average'};
%scatter(tas2,pr2,250,'filled','MarkerFaceColor','k')
%axis([0 100 0 100])
hold on
p(1) = plot(tas1,pr1,'o','MarkerFaceColor','b','MarkerEdgeColor','none','MarkerSize',12);
for i=1:8
    p(i+1) = plot(tas1(col_regs{i}),pr1(col_regs{i}),'o','MarkerFaceColor',fc(i,:),'MarkerEdgeColor','none','MarkerSize',12);
end

p(10) = plot(mean(tas1),mean(pr1),'o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',20);
%plot(mean(tas1),mean(pr1),'o','MarkerFaceColor','none','MarkerEdgeColor','k','MarkerSize',105);

%scatter(tas3,pr3,50,'filled','MarkerFaceColor','b')
box on
% plot([50,50],[0,100],'k')
% plot([0,100],[50,50],'k')


[hl]=legend(p(:),regions_n,'Location','southwest');
hl.FontSize = 15 ; hl.LineWidth = 2; hl.FontWeight = 'bold';

M = findobj(hl,'type','line'); % Find objects of type 'patch'
set(M,'MarkerSize',350) %Calculate marker size based on size of scatter points
%set(gca,'visible','off')
xlabel('tas fractional variance (%)');
ylabel('pr fractional variance (%)');
title({'Fraction of modeling variance';'driven by aerosols uncertainties in GeoMIP'})
set(gca,'Linewidth',2,'FontSize',18,'FontWeight','Bold')

(mean(tas1)-50)/(100-mean(tas1))*100
(mean(pr1)-50)/(100-mean(pr1))*100

for i=2:8
   xx(i) =  mean(tas1(col_regs{i}));
t(i) = (mean(tas1(col_regs{i}))-50)/(100-mean(tas1(col_regs{i})))*100;
pr(i) = mean(pr1(col_regs{i})-50)/(100-mean(pr1(col_regs{i})))*100;
end

set(gcf,'renderer','painters')
print(gcf,'-depsc2',['GeoMIP_uncertainties_ssp.eps'])