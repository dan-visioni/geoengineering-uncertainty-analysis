clear all

fc1 = brewermap(5,'Set1');
fc2 = brighten(fc1,.8);
M_1 = fc1(3,:);
M_2 = fc1(1,:);
M_3 = fc1(2,:);
M_4 = fc1(4,:);
M_5 = fc1(5,:);
m_1 = fc2(1,:);
m_2 = fc2(2,:);
m_3 = fc2(3,:);

fold  = '/Volumes/D_Visioni2/GeoMIP/'; %Change with directory where you downloaded data
n_exp = {'G6sulfur','G6solar','ssp585','ssp245'}; %Experiments name
n_mod = {'CNRM-ESM2-1','IPSL-CM6A-LR','CESM2-WACCM','UKESM1-0-LL','MPI-ESM1-2-LR','MPI-ESM1-2-HR'}; %Available models

%%% Set up of available ensemble members %%%
n_ens{1,1}=[1,2,3];n_ens{1,2}=[1];n_ens{1,3}=n_ens{1,1};n_ens{1,4}=n_ens{1,1};
n_ens{2,1}=[1];n_ens{2,2}=n_ens{2,1};n_ens{2,3}=n_ens{2,1};n_ens{2,4}=n_ens{2,1};
n_ens{3,1}=[1,2];n_ens{3,2}=[1,2];n_ens{3,3}=[1];n_ens{3,4}=[1];
n_ens{4,1}=[1,4,8];n_ens{4,2}=[1,4,8];n_ens{4,3}=[1,4,8];n_ens{4,4}=[1,4,8];
n_ens{5,1}=[1,2,3];n_ens{5,2}=n_ens{5,1};n_ens{5,3}=n_ens{5,1};n_ens{5,4}=n_ens{5,1};
n_ens{6,1}=[1,2,3];n_ens{6,2}=n_ens{6,1};n_ens{6,3}=[1,2];n_ens{6,4}=n_ens{6,3};
n_ens{7,1}=[2,3];n_ens{7,2}=n_ens{7,1};n_ens{7,3}=n_ens{7,1};n_ens{7,4}=n_ens{7,1};

%if the beginning year is different than 2015 (some GeoMIP models, not all)
y_start = [0,5,5,5,0,5,5; 0,5,5,5,0,0,0; 0,0,0,0,0,0,0; 0,0,0,0,0,0,0];

%variable names
var1 = {'tas','pr'};
%variable explanation
var_l = {'TS','precipitation'};
%in case a change of unit is necessary
corr = [1,86400];

%select variable to analyse
iv = 2;

cc = 0;

fprintf('Model\tens size\t \n')

for im=1:6
    fprintf([n_mod{im} '\t(' ')\t'])
    
    for ie=1:4
        
        fprintf([n_exp{ie} ' '])
        nne = n_ens{im,ie};
        for i=1:length(nne)
            list = dir([fold n_exp{ie} '/' n_mod{im} '/' var1{iv} '_*' n_mod{im} ...
                '*' n_exp{ie} '_r' num2str(nne(i)) '*.nc']);
            if isempty(list)
                x = 0;
                if i==1, fprintf('(N) '); end
            else
                if i==1, fprintf('(Y) '); end
                for j=1:length(list)
                    x2 = ncread([list(j).folder '/' list(j).name],var1{iv})*corr(iv);
                    
                    if im==3 && j==1 && ie<3 %fixes a problem with CESM2 having weird length for the beginning
                        size_cesm = size(x2);
                        if ie==2 && i==1, x2 = x2(:,:,13:size_cesm(3)); end
                        if i==2, x2 = x2(:,:,61:size_cesm(3)); end
                    end
                    
                    if length(list)>1 && j>1
                        x = cat(3,x,x2);
                    else
                        x = x2;
                    end
                end
            end
            
            if im==5 && ie<3 %fixes a problem with MIP-LR sometimes missing the last year
                if size(x,3) == 1020
                    x(:,:,1021:1032) = x(:,:,1009:1020);
                end
            end
            if im==6 && ie<2
                if size(x,3) == 960
                    x(:,:,961:972) = x(:,:,949:960);
                end
            end
            if im==6 && ie>=2
                if size(x,3) == 1020
                    x(:,:,1021:1032) = x(:,:,1009:1020);
                end
            end 
            
        x = reshape(x,size(x,1),size(x,2),12,size(x,3)/12);
        x = squeeze(mean(x,3));
        x_e{im,ie,i} = x;
                
        end
        if isempty(list), else
            lat = ncread([list(1).folder '/' list(1).name],'lat'); LAT{im} = lat;
            lon = ncread([list(1).folder '/' list(1).name],'lon'); LON{im} = lon;
            time = ncread([list(1).folder '/' list(1).name],'time');
        end
        
        %%
        
        if isempty(list)
        else
            cc = cc+1;
            n_exp2{cc} = n_exp{ie}; %skip over missing experiments, if there are
            ww = cos(lat/180*pi);
            x_zon = squeeze(nanmean(x,1)); %average over longitude
            siz = size(x);
            
            
            
        end
        
        %clear x x_e
        
    end
    fprintf(['\n'])
end

save('GeoMIP_pr.mat','x_e','LAT','LON')

%%