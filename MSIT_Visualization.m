%% Figure 1d: Boxplot of log RT: conflict and group
% Get data from GLM
data=Mrt2.Variables;
data_ECLow = data(data.Conflict==categorical(0) & data.Group == categorical(0),:);
data_ADLow = data(data.Conflict==categorical(0) & data.Group == categorical(1),:);
data_ECHigh = data(data.Conflict==categorical(1) & data.Group == categorical(0),:);
data_ADHigh = data(data.Conflict==categorical(1) & data.Group == categorical(1),:);
grpData = {data_ECLow;data_ADLow;data_ECHigh;data_ADHigh};

%Get mean RTs of each participant
grpmeans=[];
for gp = 1:length(grpData)
    grpTable = grpData{gp};
    subList = unique(grpTable.Subject);
    submeans=[];
    for ss =1:length(subList)
        subData = grpTable(grpTable.Subject==subList(ss),:);
        subData = varfun(@(x) mean(x, "omitnan"), subData,"InputVariables",{'LogRT'},'GroupingVariables',{'Subject','Group','Conflict'});
        submeans(ss) = nanmean(subData.Fun_LogRT);
    end
    grpmeans{gp}=submeans;
end

%Make boxplots
X = 1:2:7;
figure();
bplot(data_ECLow.LogRT,X(1),'width',1,'linewidth',1,'color','k','facecolor','w','nomean');hold on
scatter(X(1).*ones(length(grpmeans{1}),1),grpmeans{1},30,'k','x','XJitter','rand');
bplot(data_ADLow.LogRT,X(2),'width',1,'linewidth',1,'color','k','facecolor',[0.93,0.83,0.29],'nomean');
scatter(X(2).*ones(length(grpmeans{2}),1),grpmeans{2},30,'k','x','XJitter','rand');
bplot(data_ECHigh.LogRT,X(3),'width',1,'linewidth',1,'color','k','facecolor','w','nomean');
scatter(X(3).*ones(length(grpmeans{3}),1),grpmeans{3},30,'k','x','XJitter','rand');
bplot(data_ADHigh.LogRT,X(4),'width',1,'linewidth',1,'color','k','facecolor',[0.93,0.83,0.29],'nomean');
scatter(X(4).*ones(length(grpmeans{4}),1),grpmeans{4},30,'k','x','XJitter','rand');

xlim([0 8])
xticks([2 6]);
xticklabels({'Low Conflict','High Conflict'})
ylim([-0.8 0.8]);
yticks([-0.8:0.4:0.8])
fontsize(gca,10,"points")
sigstar({[1,3],[1,5],[3,7]},[padj_posthoc(2:4)])
title('Response Time')
%% Figure 2a: Spectrograms for left dlPFC and right LTL
clear all
cd D:/Data

%Load participant labels
load('state_map_082924.mat', 'map_ad')
id_ec=map_ad((map_ad(:,2)==0),1);
id_ad=map_ad((map_ad(:,2)==1),1);

%Select regions to plot

%All regions included in analysis
%R=["L1","R1","L2","R2","L4","R4","L6","R6","L9","R9","L11","R11","L13","R13","L14","R14"];
%regionTitles=[{'L dlPFC'},{'R dlPFC'},{'L dmPFC'},{'R dmPFC'},{'L OFC'},{'R OFC'},{'L vlPFC'},{'R vlPFC'},{'L LTL'},{'R LTL'},{'L dACC'},{'R dACC'},{'L AMY'},{'R AMY'},{'L HC'},{'R HC'}];

%For paper: L1 = left dlPFC, R9 = right LTL
R=["L1","R9"];
regionTitles=[{'L dlPFC'},{'R LTL'}];

%Initialize power matrices
EC_psd_high_norm = cell(1,length(regionTitles));
EC_psd_low_norm = cell(1,length(regionTitles));
AD_psd_high_norm = cell(1,length(regionTitles));
AD_psd_low_norm = cell(1,length(regionTitles));

%Get log power TF values from EC participants
for ss=1:length(id_ec)
    file=['P',num2str(id_ec(ss)),'.mat'];
    load(file,'ft_freq_clean', 'Parcellation_Sided','TrialDet' );
    %Combine dlPFC, vlPFC, OFC
    for bb=1:length(Parcellation_Sided)
        switch Parcellation_Sided{bb}
            case "L5"
                Parcellation_Sided{bb} = 'L4'; %conbine left med/lat OFC
            case "R5"
                Parcellation_Sided{bb} = 'R4'; %combine right med/lat OFC
            case "L8"
                Parcellation_Sided{bb} = 'L6';%combine left vlFPC
            case "R8"
                Parcellation_Sided{bb} = 'R6';%combine right vlPFC
            case "L3"
                Parcellation_Sided{bb} = 'L1';%combine left dlPFC/dlPFCpost
            case "R3"
                Parcellation_Sided{bb} = 'R1';%combine right dlPFC/dlPFCpost
        end
    end
 
    freq=ft_freq_clean.freq; 
    time=ft_freq_clean.time;
    T0=find(time<=0 & time>=-0.5); 
    trial_low=find(~isnan(TrialDet(:,12)) & TrialDet(:,12)<2 & TrialDet(:,14)==min(TrialDet(:,14)) & TrialDet(:,27)==1); %Looks for response, conflict type, and Accuracy
    trial_high=find(~isnan(TrialDet(:,12)) & TrialDet(:,12)<2 & TrialDet(:,14)==max(TrialDet(:,14)) & TrialDet(:,27)==1);
   
    %Average over channels by participant
    for ii = 1:length(R)
        ch = find(Parcellation_Sided== R(ii));
        psd_high=ft_freq_clean.powspctrm(trial_high,ch,:,:); 
        psd0_high=nanmean(psd_high(:,:,:,T0),4);
        psd_high_norm=squeeze(nanmean(nanmean(psd_high./repmat(psd0_high,1,1,1,length(time)),1),2)); %Normalize power wrt power @ t=0 then average across trials then electrodes
        psd_low=ft_freq_clean.powspctrm(trial_low,ch,:,:); 
        psd0_low=nanmean(psd_low(:,:,:,T0),4);
        psd_low_norm=squeeze(nanmean(nanmean(psd_low./repmat(psd0_low,1,1,1,length(time)),1),2));
        EC_psd_high_norm{ii}(end+1,:,:) = psd_high_norm; %dims = {region} sub X freq X time
        EC_psd_low_norm{ii}(end+1,:,:) = psd_low_norm;
    end
end

%Get log power TF values from A/D participants
for ss=1:length(id_ad)
    file=['P',num2str(id_ad(ss)),'.mat'];
    load(file,'ft_freq_clean', 'Parcellation_Sided','TrialDet' );

    %Combine vlPFC, ACC, and OFC
    for bb=1:length(Parcellation_Sided)
        switch Parcellation_Sided{bb}
            case "L5"
                Parcellation_Sided{bb} = 'L4'; %conbine left med/lat OFC
            case "R5"
                Parcellation_Sided{bb} = 'R4'; %combine right med/lat OFC
            case "L8"
                Parcellation_Sided{bb} = 'L6';%combine left vlFPC
            case "R8"
                Parcellation_Sided{bb} = 'R6';%combine right vlPFC
            case "L3"
                Parcellation_Sided{bb} = 'L1';%combine left dlPFC/dlPFCpost
            case "R3"
                Parcellation_Sided{bb} = 'R1';%combine right dlPFC/dlPFCpost
        end
    end

    freq=ft_freq_clean.freq; 
    time=ft_freq_clean.time;
    T0=find(time<=0 & time>=-0.5); 
    trial_low=find(~isnan(TrialDet(:,12)) & TrialDet(:,12)<2 & TrialDet(:,14)==min(TrialDet(:,14)) & TrialDet(:,27)==1); %Looks for response, conflict type, and Accuracy
    trial_high=find(~isnan(TrialDet(:,12)) & TrialDet(:,12)<2 & TrialDet(:,14)==max(TrialDet(:,14)) & TrialDet(:,27)==1);
    
    for ii = 1:length(R)
        ch = find(Parcellation_Sided== R(ii));
        psd_high=ft_freq_clean.powspctrm(trial_high,ch,:,:); 
        psd0_high=nanmean(psd_high(:,:,:,T0),4);
        psd_high_norm=squeeze(nanmean(nanmean(psd_high./repmat(psd0_high,1,1,1,length(time)),1),2)); %Normalize power wrt power @ t=0 then average across trials then electrodes
        psd_low=ft_freq_clean.powspctrm(trial_low,ch,:,:); 
        psd0_low=nanmean(psd_low(:,:,:,T0),4);
        psd_low_norm=squeeze(nanmean(nanmean(psd_low./repmat(psd0_low,1,1,1,length(time)),1),2));
        AD_psd_high_norm{ii}(end+1,:,:)= psd_high_norm; %dims = {region} sub X freq X time
        AD_psd_low_norm{ii}(end+1,:,:) = psd_low_norm;
    end
end
%save(['D:/MSIT Analysis/Region_level/','Region_PSD_CorrectRespOnly_082924.mat'],"EC_psd_low_norm","EC_psd_high_norm","AD_psd_low_norm","AD_psd_high_norm",'-v7.3');

%Plot spectrograms
yvec = 1:length(freq);
ylab=[5:10:55,70:20:110];
fidx=find(ismember(freq,ylab));

for ii = 1:length(R)
    f = figure();
    subplot(3,3,1)
    imagesc(time,yvec,squeeze(nanmean(log(EC_psd_high_norm{ii}))),'Interpolation','bilinear');axis('xy');
    xlim([0 2]);ylim([1 fidx(end)]);yticks(fidx);yticklabels(ylab);
    cb = colorbar;set(cb,'YTick',[0:0.25:1.25]);clim([0 1.25]);fontsize(gca,10,'points');
    title('EC High Conflict')

    subplot(3,3,2)
    imagesc(time,yvec,squeeze(nanmean(log(EC_psd_low_norm{ii}))),'Interpolation','bilinear');axis('xy');
    xlim([0 2]);ylim([1 fidx(end)]);yticks(fidx);yticklabels(ylab);
    cb = colorbar;set(cb,'YTick',[0:0.25:1.25]);clim([0 1.25]);fontsize(gca,10,'points');
    title('EC Low Conflict')
    
    subplot(3,3,3)
    imagesc(time,yvec,squeeze(nanmean(log(EC_psd_high_norm{ii})))-squeeze(nanmean(log(EC_psd_low_norm{ii}))),'Interpolation','bilinear');axis('xy');
    xlim([0 2]);ylim([1 fidx(end)]);yticks(fidx);yticklabels(ylab);
    cb = colorbar;set(cb,'YTick',[-0.2:0.1:0.2]);clim([-0.2 0.2]);fontsize(gca,10,'points');
    title('EC High-Low Conflict')

    subplot(3,3,4)
    imagesc(time,yvec,squeeze(nanmean(log(AD_psd_high_norm{ii}))),'Interpolation','bilinear');axis('xy');
    xlim([0 2]);ylim([1 fidx(end)]);yticks(fidx);yticklabels(ylab);
    cb = colorbar;set(cb,'YTick',[0:0.25:1.25]);clim([0 1.25]);fontsize(gca,10,'points');
    title('AD High Conflict')
    
    subplot(3,3,5)
    imagesc(time,freq,squeeze(nanmean(log(AD_psd_low_norm{ii}))),'Interpolation','bilinear');axis('xy');
    xlim([0 2]);ylim([1 fidx(end)]);yticks(fidx);yticklabels(ylab);cb = colorbar;set(cb,'YTick',[0:0.25:1.25]);clim([0 1.25]);fontsize(gca,10,'points');
    title('AD Low Conflict')

    subplot(3,3,6)
    imagesc(time,yvec,squeeze(nanmean(log(AD_psd_high_norm{ii})))-squeeze(nanmean(log(AD_psd_low_norm{ii}))),'Interpolation','bilinear');axis('xy');
    xlim([0 2]);ylim([1 fidx(end)]);yticks(fidx);yticklabels(ylab);
    cb = colorbar;set(cb,'YTick',[-0.2:0.1:0.2]);clim([-0.2 0.2]);
    fontsize(gca,10,'points');
    title('AD High Conflict-Low Conflict')
    sgtitle(regionTitles(ii));    
end

%% Figure 2B-E: Interaction heatmap and Boxplots of power by conflict and group
%Define regions and frequency bands
F=[4,8;8,15;15,30;30,55;70,110];
FrequencyBands = [{'Theta'},{'Alpha'},{'Beta'},{'Gamma'},{'High Gamma'}];
R=["L1","R1","L2","R2","L4","R4","L6","R6","L9","R9","L11","R11","L13","R13","L14","R14"];
RegionNames=[{'L dlPFC'},{'R dlPFC'},{'L dmPFC'},{'R dmPFC'},{'L OFC'},{'R OFC'},{'L vlPFC'},{'R vlPFC'},{'L LTL'},{'R LTL'},{'L dACC'},{'R dACC'},{'L AMY'},{'R AMY'},{'L HC'},{'R HC'}];

%Sort regions into left then right
Regions_sorted = [{'L dlPFC'},{'L dmPFC'},{'L vlPFC'},{'L OFC'},{'L dACC'},{'L LTL'},{'L AMY'},{'L HC'},{'R dlPFC'},{'R dmPFC'},{'R vlPFC'},{'R OFC'},{'R dACC'},{'R LTL'},{'R AMY'},{'R HC'}];
label_sorted = [{'dlPFC'},{'dmPFC'},{'vlPFC'},{'OFC'},{'dACC'},{'LTL'},{'AMY'},{'HC'},{'dlPFC'},{'dmPFC'},{'vlPFC'},{'OFC'},{'dACC'},{'LTL'},{'AMY'},{'HC'}];
[~,ridx]=ismember(Regions_sorted,RegionNames);
M1s = M1(ridx,:);
Rs = R(ridx);
conflict_group_sorted = coef_conflict_group(ridx,:);
padj_conflict_group_sorted = padj_conflict_group(ridx,:);
sigsorted = padj_conflict_group_sorted<0.05;

%Get significant regions from power models
Rsig =[];
Fsig =[];
for fb = 1:length(F)
    intreg = sigInteractionModels{2,fb};
    for rr = 1:length(intreg)
        Rsig(end+1) = find(Rs==intreg(rr));
        Fsig(end+1) = fb;
    end
end


%Figure 2b: Interaction heatmap
imagesc([0.5 15.5],[0.5 4.5],(conflict_group_sorted)');
%xlabel('Band')
xticks(0.5:1:15.5)
yticks(0.5:1:4.5)
colormap(redblue);b=colorbar;
clim([-0.1 0.1])
%b.Label = 'Conflict x Group \beta';
yticklabels(FrequencyBands)
xticklabels(label_sorted)
%ylabel('Region')
set(gca,"XAxisLocation",'bottom','FontUnits','points','FontSize',10,'TickDir','none')
box on
        

%Figures 2c-e: Power Plots by conflict and group
labelloc=[];
for ff=1:5
    regs= sort(Rsig(Fsig==ff));
    if isempty(regs)
        continue;
    end
    labels = Regions_sorted(regs);
    x =[0,1.2,2.4,3.6];
    bandpvals = squeeze(padj_posthoc(ff,1:length(regs),:));
    if isscalar(regs)
        bandpvals=bandpvals';
    end
    figure()
    for rr=1:length(regs)
        locs = x+6*(rr);
        labelloc(rr)=median(locs);

        data=M1s{regs(rr),ff}.Variables; %Get original data table
        % Get data from each group and conflict level;
        data_ECLow = data(data.Conflict==categorical(0) & data.Group == categorical(0),:);
        data_ADLow = data(data.Conflict==categorical(0) & data.Group == categorical(1),:);
        data_ECHigh = data(data.Conflict==categorical(1) & data.Group == categorical(0),:);
        data_ADHigh = data(data.Conflict==categorical(1) & data.Group == categorical(1),:);
        grpData = {data_ECLow;data_ADLow;data_ECHigh;data_ADHigh};
        %Get individual participant data
        grpmeans=[];
        for gp = 1:length(grpData)
            grpTable = grpData{gp};
            subList = unique(grpTable.Subject);
            submeans=[];
            for ss =1:length(subList)
                subData = grpTable(grpTable.Subject==subList(ss),:);
                subData = varfun(@(x) mean(x, "omitnan"), subData,"InputVariables",{'LogPow'},'GroupingVariables',{'Subject','Group','Conflict','Electrode'});
                submeans(ss) = nanmean(subData.Fun_LogPow);
            end
            grpmeans{gp}=submeans;
        end

        %Make boxplots
        bplot(data_ECLow.LogPow,locs(1),'width',1,'linewidth',1,'color','k','facecolor','w','nomean');hold on
        scatter(locs(1).*ones(length(grpmeans{1}),1),grpmeans{1},30,'k','x','XJitter','rand');
        bplot(data_ADLow.LogPow,locs(2),'width',1,'linewidth',1,'color','k','facecolor',[0.93,0.83,0.29],'nomean'); hold on
        scatter(locs(2).*ones(length(grpmeans{2}),1),grpmeans{2},30,'k','x','XJitter','rand');
        bplot(data_ECHigh.LogPow,locs(3),'width',1,'linewidth',1,'color','r','facecolor','w','nomean'); hold on
        scatter(locs(3).*ones(length(grpmeans{3}),1),grpmeans{3},30,'k','x','XJitter','rand');
        bplot(data_ADHigh.LogPow,locs(4),'width',1,'linewidth',1,'color','r','facecolor',[0.93,0.83,0.29],'nomean');hold on
        scatter(locs(4).*ones(length(grpmeans{4}),1),grpmeans{4},30,'k','x','XJitter','rand');
        
        %Add significance markers
        comps = [locs(3),locs(4);locs(1),locs(2);locs(1),locs(3);locs(2),locs(4)];
        for k = 1:length(comps)
            if(bandpvals(rr,k)<0.05)
                sigstar(comps(k,:),bandpvals(rr,k))
            end
        end
    end
    set(gca,'xlim',[5 x(4)+6*rr+1],'xtick',labelloc,'xticklabel',labels,'TickDir','out');
    if ff<=2
        set(gca,'ytick',[-1.4:0.7:2.1],'ylim',[-1.4 2.1]);
    elseif ff>2
        set(gca,'ytick',[-1:0.5:1.5],'ylim',[-1 1.5]);
    end
    box off
    %title(FrequencyBands{ff});
end
%% Figure 3a: Coherence and PLV TF Plot for left dlPFC-right LTL
clear all
cd D:/Data
%Load participant list w/categorization
load('state_map.mat', 'map_ad')

id_ec=map_ad((map_ad(:,2)==0),1);
id_ad=map_ad((map_ad(:,2)==1),1);

%Select regions to plot
%All included regions
%R=["L1","R1","L2","R2","L4","R4","L6","R6","L9","R9","L10","R10","L13","R13","L14","R14"];
%regionTitles=[{'L dlPFC'},{'R dlPFC'},{'L dmPFC'},{'R dmPFC'},{'L OFC'},{'R OFC'},{'L vlPFC'},{'R vlPFC'},{'L LTL'},{'R LTL'},{'L ACC'},{'R ACC'},{'L AMY'},{'R AMY'},{'L HC'},{'R HC'}];

%Regions for paper: L1 = left dlPFC, R9 = right LTL
R=["L1","R9"];
regionTitles=[{'L dlPFC'},{'R LTL'}];

EC_coh_high = cell(length(R),length(R));
EC_coh_low = cell(length(R),length(R));
AD_coh_high = cell(length(R),length(R));
AD_coh_low = cell(length(R),length(R));

%Get coherence/PLV TF values for EC
for ss=1:length(id_ec)
    file=['P',num2str(id_ec(ss)),'.mat'];
    load(['D:\Data\Coherence Data\ChannelCoh_allbands_CorrectRespOnly\',file],'ft_coh_bands','regionIDs');

    % FOR PLOTTING COHERENCE
    %Combine theta-high gamma bands
    theta_low = ft_coh_bands{1,1}.cohspctrm;
    theta_high = ft_coh_bands{1,2}.cohspctrm;
    alpha_low = ft_coh_bands{2,1}.cohspctrm;
    alpha_high = ft_coh_bands{2,2}.cohspctrm;
    beta_low = ft_coh_bands{3,1}.cohspctrm;
    beta_high = ft_coh_bands{3,2}.cohspctrm;
    gamma_low = ft_coh_bands{4,1}.cohspctrm;
    gamma_high = ft_coh_bands{4,2}.cohspctrm;
    highgamma_low = ft_coh_bands{5,1}.cohspctrm;
    highgamma_high= ft_coh_bands{5,2}.cohspctrm;
    
    %FOR PLOTTING PLV  (coh below now refers to PLV)
    % load(['D:\Data\PLV Data\PLV_CorrectRespOnly\',file],'ft_plv_bands','regionIDs'); 
    %Combine theta-high gamma bands
    % theta_low = ft_plv_bands{1,1}.plvspctrm;
    % theta_high = ft_plv_bands{1,2}.plvspctrm;
    % alpha_low = ft_plv_bands{2,1}.plvspctrm;
    % alpha_high = ft_plv_bands{2,2}.plvspctrm;
    % beta_low = ft_plv_bands{3,1}.plvspctrm;
    % beta_high = ft_plv_bands{3,2}.plvspctrm;
    % gamma_low = ft_plv_bands{4,1}.plvspctrm;
    % gamma_high = ft_plv_bands{4,2}.plvspctrm;
    % highgamma_low = ft_plv_bands{5,1}.plvspctrm;
    % highgamma_high = ft_plv_bands{5,2}.plvspctrm;
    
    coh_low.cohspctrm = cat(3,theta_low,alpha_low(:,:,2:end,:),beta_low(:,:,2:end,:),gamma_low(:,:,2:end,:),highgamma_low(:,:,:,:));
    coh_high.cohspctrm = cat(3,theta_high,alpha_high(:,:,2:end,:),beta_high(:,:,2:end,:),gamma_high(:,:,2:end,:),highgamma_high(:,:,:,:));
    
    for dd = 1:length(R)
        R1 = find(regionIDs==R(dd));
        for jj = 1:length(R)
            R2 = find(regionIDs==R(jj));
            if ~isempty(R1) && ~isempty(R2)
            EC_coh_low{dd,jj}(end+1,:,:) = squeeze(nanmean(nanmean(abs(coh_low.cohspctrm(R1,R2,:,:)),1),2));
            EC_coh_high{dd,jj}(end+1,:,:) = squeeze(nanmean(nanmean(abs(coh_high.cohspctrm(R1,R2,:,:)),1),2));
            end
        end
    end
end

%Get coherence/PLV TF values for A/D
for ss=1:length(id_ad)
    file=['P',num2str(id_ad(ss)),'.mat'];
    load(['D:\Data\Coherence Data\ChannelCoh_allbands_CorrectRespOnly\',file],'ft_coh_bands','regionIDs');

    %Combine theta-high gamma bands
    %FOR PLOTTING COHERENCE
    theta_low = ft_coh_bands{1,1}.cohspctrm;
    theta_high = ft_coh_bands{1,2}.cohspctrm;
    alpha_low = ft_coh_bands{2,1}.cohspctrm;
    alpha_high = ft_coh_bands{2,2}.cohspctrm;
    beta_low = ft_coh_bands{3,1}.cohspctrm;
    beta_high = ft_coh_bands{3,2}.cohspctrm;
    gamma_low = ft_coh_bands{4,1}.cohspctrm;
    gamma_high = ft_coh_bands{4,2}.cohspctrm;
    highgamma_low = ft_coh_bands{5,1}.cohspctrm;
    highgamma_high= ft_coh_bands{5,2}.cohspctrm;
    
    %FOR PLOTTING PLV  (coh below now refers to PLV)
    % load(['D:\Data\PLV Data\PLV_CorrectRespOnly\',file],'ft_plv_bands','regionIDs');
    % theta_low = ft_plv_bands{1,1}.plvspctrm;
    % theta_high = ft_plv_bands{1,2}.plvspctrm;
    % alpha_low = ft_plv_bands{2,1}.plvspctrm;
    % alpha_high = ft_plv_bands{2,2}.plvspctrm;
    % beta_low = ft_plv_bands{3,1}.plvspctrm;
    % beta_high = ft_plv_bands{3,2}.plvspctrm;
    % gamma_low = ft_plv_bands{4,1}.plvspctrm;
    % gamma_high = ft_plv_bands{4,2}.plvspctrm;
    % highgamma_low = ft_plv_bands{5,1}.plvspctrm;
    % highgamma_high = ft_plv_bands{5,2}.plvspctrm;

    coh_low.cohspctrm = cat(3,theta_low,alpha_low(:,:,2:end,:),beta_low(:,:,2:end,:),gamma_low(:,:,2:end,:),highgamma_low(:,:,:,:));
    coh_high.cohspctrm = cat(3,theta_high,alpha_high(:,:,2:end,:),beta_high(:,:,2:end,:),gamma_high(:,:,2:end,:),highgamma_high(:,:,:,:));

    for dd = 1:length(R)
        R1 = find(regionIDs==R(dd));
        for jj = 1:length(R)
            R2 = find(regionIDs==R(jj));
            if ~isempty(R1) && ~isempty(R2)
            AD_coh_low{dd,jj}(end+1,:,:) = squeeze(nanmean(nanmean(abs(coh_low.cohspctrm(R1,R2,:,:)),1),2));
            AD_coh_high{dd,jj}(end+1,:,:) = squeeze(nanmean(nanmean(abs(coh_high.cohspctrm(R1,R2,:,:)),1),2));
            end
        end
    end
end

%save('Coherence_TF_082924.mat','AD_coh_low',"AD_coh_high","EC_coh_low","EC_coh_high","Regions_sorted",'R','median_low','median_high')

% Make TF plots of coherence or PLV
time = ft_coh_bands{1,1}.time;time(end)=2;%round up any timing differences
freq= [4:55,70:110];
yvec = 1:length(freq);
ylab=[5:5:25];
fidx=find(ismember(freq,ylab));

%Get median RTs from the
load('RTdata.mat','EC','AD','ECconflict','ADconflict');
AD1=AD(ADconflict==0);AD2=AD(ADconflict==1);
EC1=EC(ECconflict==0);EC2=EC(ECconflict==1);
median_low = [median(EC1),median(AD1)];
median_high = [median(EC2),median(AD2)];

for dd = 1:length(R)
    f = figure();
    for jj = 1:length(R)
    subplot(3,3,1)
    imagesc(time,yvec,squeeze(nanmean(EC_coh_high{dd,jj})),'Interpolation','bilinear');axis('xy');
    xlim([0 1.5]);ylim([1 fidx(end)]);yticks(fidx);yticklabels(ylab);
    colorbar;clim([0 .15]);
    title('EC High Conflict');
    
    subplot(3,3,2)
    imagesc(time,yvec,squeeze(nanmean(EC_coh_low{dd,jj})),'Interpolation','bilinear');axis('xy');
    xlim([0 2]);ylim([1 fidx(end)]);yticks(fidx);yticklabels(ylab);
    colorbar;clim([0 .035]);
    title('EC Low Conflict');

    subplot(3,3,3)
    imagesc(time,yvec,(squeeze(nanmean(EC_coh_high{dd,jj})))-(squeeze(nanmean(EC_coh_low{dd,jj}))),'Interpolation','bilinear');axis('xy');
    xlim([0 2]);ylim([1 fidx(end)]);yticks(fidx);yticklabels(ylab);
    cb = colorbar;set(cb,'YTick',[-0.01:0.005:0.01]);clim([-0.01 .01]);
    hold on
    xline(median_low(1),'--k');
    xline(median_high(1),'--r');
    fontsize(gca,28,'pixels');
    title('EC High-Low Conflict');

    subplot(3,3,4)
    imagesc(time,yvec,squeeze(nanmean(AD_coh_high{dd,jj})),'Interpolation','bilinear');axis('xy');
    xlim([0 2]);ylim([1 fidx(end)]);yticks(fidx);yticklabels(ylab);
    colorbar;clim([0 0.035])
    title('AD High Conflict')

    subplot(3,3,5)
    imagesc(time,yvec,squeeze(nanmean(AD_coh_low{dd,jj})),'Interpolation','bilinear');axis('xy');
    xlim([0 2]);ylim([1 fidx(end)]);yticks(fidx);yticklabels(ylab);
    colorbar;clim([0 .035])
    title('AD Low Conflict');

    subplot(3,3,6)
    imagesc(time,yvec,(squeeze(nanmean(AD_coh_high{dd,jj},1)))-(squeeze(nanmean(abs(AD_coh_low{dd,jj}),1))),'Interpolation','bilinear');axis('xy');
    xlim([0 2]);ylim([1 fidx(end)]);yticks(fidx);yticklabels(ylab);
    cb = colorbar;set(cb,'YTick',[-0.01:0.005:0.01]);clim([-0.01 .01]);
    hold on
    xline(median_low(2),'--k');
    xline(median_high(2),'--r');
    fontsize(gca,10,'points');
    title('AD High-Low Conflict');
    sgtitle([Regions_sorted{dd},'-',Regions_sorted{jj}])
    end
end

%% Figure 3b: Connectograms for significant conflict contrasts (high>low)
F=[4,8;8,15;15,30;30,55;70,110];
FrequencyBands = [{'Theta'},{'Alpha'},{'Beta'},{'Gamma'},{'High Gamma'}];
RegionNames_all=[{'L dlPFC'},{'R dlPFC'},{'L dmPFC'},{'R dmPFC'},{'L OFC'},{'R OFC'},{'L vlPFC'},{'R vlPFC'},{'L LTL'},{'R LTL'},{'L ACC'},{'R ACC'},{'L AMY'},{'R AMY'},{'L HC'},{'R HC'}];
GroupNames=[{'EC'},{'AD'}]; %Compare conflict types
num_nodes = length(RegionNames_all);
cmap = zeros(num_nodes,3);

FCmeasure = 'Coherence'; %Specified as 'Coherence' or 'PLV'

%Get significant connections from Coherence or PLV models
sigconns=cell(1,length(F));
for fb = 1:length(F)
    num_nodes=length(RegionNames{fb});
    sigconns{fb} = zeros(num_nodes,num_nodes);
    mdls = sigInteractionModels{1,fb};
    intregs = sigInteractionModels{2,fb};
    for rr = 1:length(intregs)
        regs = strsplit(intregs{rr},'-');
        idx = find(ismember(RegionNames{fb},regs));
        sigconns{fb}(idx(1),idx(2))=1;
    end
end

%Create connectivity graphs for follow-up PLV/AmpCorr contrasts
for fb = 1:length(F)
    regs = RegionNames{fb};
    num_nodes=length(regs);
    coef_ec = contrast_coef_posthoc(fb,:,3);
    pvals_ec = padj_posthoc(fb,:,3);
    coef_ad = contrast_coef_posthoc(fb,:,4);
    pvals_ad = padj_posthoc(fb,:,4);
    f=figure();
    fontsize(gca,20,'points');fontname(gca,'Arial')
    graph = circularGraph(ones(num_nodes,num_nodes),'Label',RegionNames{fb});
    graph.ShowButton.Visible=0;graph.HideButton.Visible=0;
    [graph.Node(:).Color]=deal([0 0 0]);
   
    idx=1;
    for jj = 1:length(sigconns{fb})
        for ii = 1:length(sigconns{fb})
            graph.Node(ii).Connection(jj).LineWidth=2;
            h = graph.Node(ii).Connection(jj);
            if sigconns{fb}(ii,jj)==0
                h.Color = [0.9 0.9 0.9];
            elseif sigconns{fb}(ii,jj)==1
                uistack(h,'top')
                h.LineWidth = 2;
                if pvals_ec(idx) < 0.05 && pvals_ad(idx) < 0.05
                    if (coef_ec(idx) > 0) && (coef_ad(idx) > 0) && (coef_ec(idx) > coef_ad(idx))
                        h.Color = [0 0 1];
                        h.LineStyle = '-';
                    elseif (coef_ec(idx) > 0) && (coef_ad(idx) > 0) && (coef_ec(idx) < coef_ad(idx))
                        h.Color = [1 0 0];
                        h.LineStyle = '-';
                    elseif (coef_ec(idx) > 0) && (coef_ad(idx) < 0) 
                        h.Color = [0.75 0 1];
                        h.LineStyle = '--';
                    elseif (coef_ec(idx) < 0) && (coef_ad(idx) > 0) 
                        h.Color = [0.75 0 1];
                        h.LineStyle = '-';
                    elseif (coef_ec(idx) < 0) && (coef_ad(idx) < 0) && (abs(coef_ec(idx)) > abs(coef_ad(idx)))
                        h.Color = [0 0 1];
                        h.LineStyle = "--";
                    elseif (coef_ec(idx) < 0) && (coef_ad(idx) < 0) && (abs(coef_ec(idx)) < abs(coef_ad(idx)))
                        h.Color = [1 0 0];
                        h.LineStyle = "--";
                    end
                elseif pvals_ec(idx) < 0.05 && pvals_ad(idx) > 0.05
                    if coef_ec(idx) > 0
                        h.Color = [0 0 0];
                        h.LineStyle = "-";
                    elseif coef_ec(idx) < 0
                        h.Color = [0 0 0];
                        h.LineStyle = "--";
                    end
                elseif pvals_ec(idx) > 0.05 && pvals_ad(idx) < 0.05
                        if coef_ad(idx) > 0
                        h.Color = [0.93,0.83,0.29];
                        h.LineStyle = "-";
                    elseif coef_ad(idx) < 0
                        h.Color = [0.93,0.83,0.29];
                        h.LineStyle = "--";
                        end
                end
                idx = idx+1;
            end
        end
    end
end