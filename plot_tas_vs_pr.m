clear all

load tas_regions1.mat

regions = shaperead('IPCC-WGI-reference-regions-v4.shx');

for i=1:size(regions,1)
    
    for t=1:81
        x = ta_reg_sulf(:,t,i);
        y = ta_reg_solar(:,t,i);
        
        tas_var1_sulf(i,t) = std((x),0).^2;
        tas_var1_solar(i,t) = std((y),0).^2;
    end
    tot_var_reg1_tas(i,:) = tas_var1_sulf(i,:)+tas_var1_solar(i,:);
end

regions = shaperead('IPCC_WGII_continental_regions.shx');

load tas_regions2.mat

for i=1:size(regions,1)
    
    x = ta_reg_sulf(:,:,i);
    y = ta_reg_solar(:,:,i);
    
    for t=1:81
        x = ta_reg_sulf(:,t,i);
        y = ta_reg_solar(:,t,i);
        
        tas_var2_sulf(i,t) = std((x),0).^2;
        tas_var2_solar(i,t) = std((y),0).^2;
    end
    tot_var_reg2_tas(i,:) = tas_var2_sulf(i,:)+tas_var2_solar(i,:);
    
end

load tas_regions3.mat

for i=1:238
    
    x = ta_reg_sulf(:,:,i);
    y = ta_reg_solar(:,:,i);
    
    for t=1:81
        x = ta_reg_sulf(:,t,i);
        y = ta_reg_solar(:,t,i);
        
        tas_var3_sulf(i,t) = std((x),0).^2;
        tas_var3_solar(i,t) = std((y),0).^2;
    end
    tot_var_reg3_tas(i,:) = tas_var3_sulf(i,:)+tas_var3_solar(i,:);
    
end

load pr_regions1.mat

regions = shaperead('IPCC-WGI-reference-regions-v4.shx');

for i=1:size(regions,1)
    
    x = ta_reg_sulf(:,:,i);
    y = ta_reg_solar(:,:,i);
    
    for t=1:81
        x = ta_reg_sulf(:,t,i);
        y = ta_reg_solar(:,t,i);
        
        pr_var1_sulf(i,t) = std((x),0).^2;
        pr_var1_solar(i,t) = std((y),0).^2;
    end
    tot_var_reg1_pr(i,:) = pr_var1_sulf(i,:)+pr_var1_solar(i,:);
end

load pr_regions2.mat

regions = shaperead('IPCC_WGII_continental_regions.shx');

for i=1:size(regions,1)
    
    x = ta_reg_sulf(:,:,i);
    y = ta_reg_solar(:,:,i);
    
    for t=1:81
        x = ta_reg_sulf(:,t,i);
        y = ta_reg_solar(:,t,i);
        
        pr_var2_sulf(i,t) = std((x),0).^2;
        pr_var2_solar(i,t) = std((y),0).^2;
    end
    tot_var_reg2_pr(i,:) = pr_var2_sulf(i,:)+pr_var2_solar(i,:);
end

load pr_regions3.mat


for i=1:238
    
    x = ta_reg_sulf(:,:,i);
    y = ta_reg_solar(:,:,i);
    
    for t=1:81
        x = ta_reg_sulf(:,t,i);
        y = ta_reg_solar(:,t,i);
        
        pr_var3_sulf(i,t) = std((x),0).^2;
        pr_var3_solar(i,t) = std((y),0).^2;
    end
    tot_var_reg3_pr(i,:) = pr_var3_sulf(i,:)+pr_var3_solar(i,:);
end

%%

pr3 = (pr_var3_sulf./tot_var_reg3_pr*100);
tas3 = (tas_var3_sulf./tot_var_reg3_tas*100);

pr2 = (pr_var2_sulf./tot_var_reg2_pr*100);
tas2 = (tas_var2_sulf./tot_var_reg2_tas*100);

pr1 = (pr_var1_sulf./tot_var_reg1_pr*100);
tas1 = (tas_var1_sulf./tot_var_reg1_tas*100);

scatter(mean(tas2(:,41:81),2),mean(pr2(:,41:81),2),250,'filled','MarkerFaceColor','k')
axis([45 70 45 70])
hold on
scatter(mean(tas1(:,41:81),2),mean(pr1(:,41:81),2),250,'filled','MarkerFaceColor','r')
scatter(mean(tas3(:,41:81),2),mean(pr3(:,41:81),2),50,'filled','MarkerFaceColor','b')
box on
plot([50,50],[0,100],'k')
plot([0,100],[50,50],'k')