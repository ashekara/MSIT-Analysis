%% GLME analysis of Canonical Band Power
clear all
cd D:\Data
%Load participant list w/categorization
load('state_map.mat','map_ad')

%Median RT by subject
RT_High = [];
RT_Low = [];

for ii = 1:length(map_ad)
    filename=['P',num2str(map_ad(ii,1)),'.mat'];
    load(filename,'TrialDet');
    Trials_C=find(~isnan(TrialDet(:,12)) & TrialDet(:,12)<2 & TrialDet(:,14)==min(TrialDet(:,14)) & TrialDet(:,27)==1); %Looks for response, conflict type, and Accuracy
    Trials_I=find(~isnan(TrialDet(:,12)) & TrialDet(:,12)<2 & TrialDet(:,14)==max(TrialDet(:,14)) & TrialDet(:,27)==1);
    RT_Low=[RT_Low;median(TrialDet(Trials_C,12))];
    RT_High=[RT_High;median(TrialDet(Trials_I,12))];
end

%Specify frequency bands of interest
F=[4,8;8,15;15,30;30,55;70,110];
FrequencyBands = [{'Theta'},{'Alpha'},{'Beta'},{'Gamma'},{'High Gamma'}];

%Initialize data vectors
TwinLow=RT_Low;
TwinHigh=RT_High;
LogPow=[];
LogPowMed=[];
Conflict=[];
Subject=[];
Region=[];
Group=[];
FreqBand=[];
Electrode=[];
RT=[];
Trial=[];

%Build data table from Power avg
for ff=1:length(map_ad)
    filename=['P',num2str(map_ad(ff,1)),'.mat'];
    load(filename,'ft_freq_clean','Parcellation_Sided','TrialDet'); 
    respTime = TrialDet(:,12);
    channels=find(~contains(Parcellation_Sided,'NaN'));
    Trials_C=find(~isnan(TrialDet(:,12)) & TrialDet(:,12)<2 & TrialDet(:,14)==min(TrialDet(:,14)) & TrialDet(:,27)==1 & ismember(StimLocation,{'None','L Ventral','RVentral','LMiddle','RMiddle'})); %Looks for response, conflict type, and Accuracy
    Trials_I=find(~isnan(TrialDet(:,12)) & TrialDet(:,12)<2 & TrialDet(:,14)==max(TrialDet(:,14)) & TrialDet(:,27)==1);
    Id_time0=find(ft_freq_clean.time>-0.5 & ft_freq_clean.time<0);
    Id_time1_Low=find(ft_freq_clean.time>0.1 & ft_freq_clean.time<=TwinLow(ff));
    Id_time1_High=find(ft_freq_clean.time>0.1 & ft_freq_clean.time<=TwinHigh(ff));
    Mstatus=map_ad(ff,2);
    for fb=1:length(F)
        fIdx = find(ft_freq_clean.freq >= F(fb,1) & ft_freq_clean.freq <= F(fb,2));
        %Single trial RT
        for rr=1:length(channels)
            R=Parcellation_Sided{channels(rr)};
            PowC0=[];PowC=[];PowI0=[];PowI=[];
            % Average log power between 0.1s to pre-response RT by channel
            for tt = 1:length(Trials_C)
              Id_time1 = find(ft_freq_clean.time>0.1 & ft_freq_clean.time<=respTime(Trials_C(tt)));
              PowC0(tt)=nanmean(nanmean(ft_freq_clean.powspctrm(Trials_C(tt),channels(rr),fIdx,Id_time0),4),3); %average across timepoints then across frequency band for specific channel
              PowC(tt)=nanmean(nanmean(ft_freq_clean.powspctrm(Trials_C(tt),channels(rr),fIdx,Id_time1),4),3);
            end
            for tt = 1:length(Trials_I)
              Id_time1 = find(ft_freq_clean.time>0.1 & ft_freq_clean.time<=respTime(Trials_I(tt)));
              PowI0(tt)=nanmean(nanmean(ft_freq_clean.powspctrm(Trials_I(tt),channels(rr),fIdx,Id_time0),4),3);
              PowI(tt)=nanmean(nanmean(ft_freq_clean.powspctrm(Trials_I(tt),channels(rr),fIdx,Id_time1),4),3);
            end
            Pow1=log(PowC./PowC0)'; %log ratio of power normalized to baseline
            Pow2=log(PowI./PowI0)';
            
            %Add data to matrices
            LogPow=[LogPow;Pow1;Pow2];
            RT = [RT;respTime(Trials_C);respTime(Trials_I)];
            Trial = [Trial;Trials_C;Trials_I];
            Conflict=[Conflict;zeros(length(Pow1),1);ones(length(Pow2),1)];
            sub=strings(size(Pow1,1)+size(Pow2,1),1);
            sub(:)=['P',num2str(map_ad(ff,1))];
            Subject=[Subject;sub];
            Region=[Region;repmat(string(R),size(Pow1,1)+size(Pow2,1),1)];
            Group=[Group;Mstatus.*ones(size(Pow1,1)+size(Pow2,1),1)];
            FreqBand=[FreqBand;fb.*ones(size(Pow1,1)+size(Pow2,1),1)];
        end
    end
end

% Adjust region labels
Region(Region=='L3')='L1'; %combine left dlPFC/dlPFCpost
Region(Region=='R3')='R1'; %combine right dlPFC/dlPFCpost
Region(Region=='L8')='L6'; %combine left vlPFCpt/vlPFCpo
Region(Region=='R8')='R6'; %combine right vlPFCpt/vlPFCpo
Region(Region=='L5')='L4'; %combine left med/lat OFC
Region(Region=='R5')='R4'; %combine right med/lat OFC

%::::::::::::::::::::::::::GLME Code::::::::::::::::::::::::::
%All regions
%R=["L1","R1","L2","R2","L3","R3","L4","R4","L5","R5","L6","R6","L9","R9","L10","R10","L11","R11","L13","R13","L14","R14","L16","R16","L17","R17"];
%RegionNames=[{'dlPFC'},{'dmPFC'},{'dlPFCpost'},{'MedOFC'},{'LatOFC'},{'vlPFC'},{'Acc'},{'LTL'},{'rACC'},{'dACC'},{'PPC'},{'Amyg'},{'Hipp'},{'insula'},{'caudate'}];

%Set Regions to model
R=["L1","R1","L2","R2","L4","R4","L6","R6","L9","R9","L11","R11","L13","R13","L14","R14"];
RegionNames=[{'L dlPFC'},{'R dlPFC'},{'L dmPFC'},{'R dmPFC'},{'L OFC'},{'R OFC'},{'L vlPFC'},{'R vlPFC'},{'L LTL'},{'R LTL'},{'L dACC'},{'R dACC'},{'L AMY'},{'R AMY'},{'L HC'},{'R HC'}];

%Initialize model variables
M1=cell(length(R),length(F));
pval_conflict1=nan(length(R),length(F));
pval_group1=nan(length(R),length(F));
pval_conflict_group=nan(length(R),length(F));
coef_conflict1=nan(length(R),length(F));
coef_group1=nan(length(R),length(F));
coef_conflict_group=nan(length(R),length(F));
tstat_conflict1=nan(length(R),length(F));
tstat_group1=nan(length(R),length(F));
tstat_conflict_group=nan(length(R),length(F));

%Run models for each frequency band and region
for fb=1:length(F)
    for rr=1:length(R)
        rid=find(Region==R(rr) & FreqBand==fb);
        if ~isempty(rid)
        dataTable = table(LogPow(rid),log(RT(rid)),categorical(Conflict(rid)),categorical(Group(rid)),categorical(Subject(rid)),categorical(Electrode(rid)),...
         'VariableNames',{'LogPow','RT','Conflict','Group','Subject','Electrode'});
        M1{rr,fb}=fitglme(dataTable,'LogPow ~ Conflict*Group  + (1|Subject)','Link','identity');
        coef_conflict1(rr,fb)=M1{rr,fb}.Coefficients.Estimate(2);
        coef_group1(rr,fb)=M1{rr,fb}.Coefficients.Estimate(3);
        coef_conflict_group(rr,fb)=M1{rr,fb}.Coefficients.Estimate(4);
        tstat_conflict1(rr,fb)=M1{rr,fb}.Coefficients.tStat(2);
        tstat_group1(rr,fb)=M1{rr,fb}.Coefficients.tStat(3);
        tstat_conflict_group(rr,fb)=M1{rr,fb}.Coefficients.tStat(4);
        pval_conflict1(rr,fb)=M1{rr,fb}.Coefficients.pValue(2);
        pval_group1(rr,fb)=M1{rr,fb}.Coefficients.pValue(3);
        pval_conflict_group(rr,fb)=M1{rr,fb}.Coefficients.pValue(4);  
        end
    end
end

%FDR Correction
[~,critp_conflict1,~,padj_conflict1]=fdr_bh(pval_conflict1);
[~,critp_group1,~,padj_group1]=fdr_bh(pval_group1);
[~,critp_conflict_group,~,padj_conflict_group]=fdr_bh(pval_conflict_group);

%Get regions with and without sigificant interactions from GLMEs
sigRegionInteractions={};
sigRegionConflict={};
sigRegionGroup={};
conflictGLMRegions={};
for dd = 1:length(F)
sigRegionInteractions{dd} = R(padj_conflict_group(:,dd)<0.05); %Region w/sig interaction
sigRegionConflict{dd} = R(padj_conflict1(:,dd)<0.05);%Regions w/significant effect of conflict
sigRegionGroup{dd} = R(padj_group1(:,dd)<0.05);%Regions w/significant effect of group
conflictGLMRegions{dd} = R(~ismember(R,sigRegionInteractions{dd})); %Region w/o significant interaction
end

%Initialize variables for reduced GLMEs
maxRegNum = max(cell2mat(cellfun(@length,conflictGLMRegions,'uni',false))); %get # regions in band with most nonsignificant regions
M2=cell(maxRegNum,length(F));
coef_conflict2=nan(maxRegNum,length(F));
tstat_conflict2=nan(maxRegNum,length(F));
pval_conflict2=nan(maxRegNum,length(F));
coef_group2=nan(maxRegNum,length(F));
tstat_group2=nan(maxRegNum,length(F));
pval_group2=nan(maxRegNum,length(F));

%Run reduced GLMEs for main effect of conflict and group on region/bands without sig. interaction
for fb=1:length(F)
    R_conflict = conflictGLMRegions{fb};
    for rr=1:length(R_conflict)
        rid=find(Region==R_conflict(rr) & FreqBand==fb);
        dataTable = table(LogPow(rid),categorical(Conflict(rid)),categorical(Group(rid)),categorical(Subject(rid)),categorical(Electrode(rid)),...
                'VariableNames',{'LogPow','Conflict','Group','Subject','Electrode'});
        M2{rr,fb}=fitglme(dataTable,'LogPow ~ Conflict + Group + (1|Subject)','Link','identity');
        coef_conflict2(rr,fb)=M2{rr,fb}.Coefficients.Estimate(2);
        coef_group2(rr,fb)=M2{rr,fb}.Coefficients.Estimate(3);
        tstat_conflict2(rr,fb)=M2{rr,fb}.Coefficients.tStat(2);
        tstat_group2(rr,fb)=M2{rr,fb}.Coefficients.tStat(3);
        pval_conflict2(rr,fb)=M2{rr,fb}.Coefficients.pValue(2); 
        pval_group2(rr,fb)=M2{rr,fb}.Coefficients.pValue(3); 
   
    end
end
[~,critp_conflict2,~,padj_conflict2]=fdr_bh(pval_conflict2);
[~,critp_group2,~,padj_group2]=fdr_bh(pval_group2);


%Get regions with significant conflict or group effects from GLMEs
sigRegionConflict2={};
sigRegionGroup2={};
for dd = 1:length(F)
    R_conflict = conflictGLMRegions{dd};
    sigRegionConflict2{dd} = R_conflict(padj_conflict2(:,dd)<0.05);%Regions w/significant effect of conflict
    sigRegionGroup2{dd} = R_conflict(padj_group2(:,dd)<0.05);%Regions w/significant effect of group
end

%Initialize results tables
%Full model table
coeffTable1 = table('Size',[length(R).*length(F).*3 7],'VariableTypes',{'string','string','string','double','double','double','double'},'VariableNames', {'Region','Frequency Band','Predictor','Coefficient','tStat','pValue','FDRpValue'});
%Reduced model table
coeffTable2 = table('Size',[sum(cellfun(@(x) ~isempty(x),M2),'all').*2 7],'VariableTypes',{'string','string','string','double','double','double','double'},'VariableNames', {'Region','Frequency Band','Predictor','Coefficient','tStat','pValue','FDRpValue'});

% Get regions included in reduced models
RegionNames_Conflict={};
for dd = 1:length(F)
RegionNames_Conflict{dd} = RegionNames(ismember(R,conflictGLMRegions{dd}));
end

%Build results tables for full and reduced models
idx = 1;
for fb = 1:length(F)
    for rr = 1:length(R)
        coeffTable1{idx,'Region'}= string(RegionNames{rr});
        coeffTable1{idx+1,'Region'}= string(RegionNames{rr});
        coeffTable1{idx+2,'Region'}= string(RegionNames{rr});
        coeffTable1{idx,'Frequency Band'}=string(FrequencyBands{fb});
        coeffTable1{idx+1,'Frequency Band'}=string(FrequencyBands{fb});
        coeffTable1{idx+2,'Frequency Band'}=string(FrequencyBands{fb});
        coeffTable1{idx,'Predictor'}="Conflict";
        coeffTable1{idx+1,'Predictor'}="Group";
        coeffTable1{idx+2,'Predictor'}="Conflict x Group";
        coeffTable1{idx,'Coefficient'}=coef_conflict1(rr,fb);
        coeffTable1{idx+1,'Coefficient'}=coef_group1(rr,fb);
        coeffTable1{idx+2,'Coefficient'}=coef_conflict_group(rr,fb);
        coeffTable1{idx,'tStat'}=tstat_conflict1(rr,fb);
        coeffTable1{idx+1,'tStat'}=tstat_group1(rr,fb);
        coeffTable1{idx+2,'tStat'}=tstat_conflict_group(rr,fb);
        coeffTable1{idx,'pValue'}=pval_conflict1(rr,fb);
        coeffTable1{idx+1,'pValue'}=pval_group1(rr,fb);
        coeffTable1{idx+2,'pValue'}=pval_conflict_group(rr,fb);
        coeffTable1{idx,'FDRpValue'}=padj_conflict1(rr,fb);
        coeffTable1{idx+1,'FDRpValue'}=padj_group1(rr,fb);
        coeffTable1{idx+2,'FDRpValue'}=padj_conflict_group(rr,fb);
        idx = idx+3;
    end
end

idx = 1;
for fb = 1:length(F)
    R_conflict = conflictGLMRegions{fb};
    for rr = 1:length(R_conflict)
        coeffTable2{idx,'Region'}= string(RegionNames_Conflict{fb}{rr});
        coeffTable2{idx+1,'Region'}= string(RegionNames_Conflict{fb}{rr});
        coeffTable2{idx,'Frequency Band'}=string(FrequencyBands{fb});
        coeffTable2{idx+1,'Frequency Band'}=string(FrequencyBands{fb});
        coeffTable2{idx,'Predictor'}="Conflict";
        coeffTable2{idx+1,'Predictor'}="Group";
        coeffTable2{idx,'Coefficient'}=coef_conflict2(rr,fb);
        coeffTable2{idx+1,'Coefficient'}=coef_group2(rr,fb);
        coeffTable2{idx,'tStat'}=tstat_conflict2(rr,fb);
        coeffTable2{idx+1,'tStat'}=tstat_group2(rr,fb);
        coeffTable2{idx,'pValue'}=pval_conflict2(rr,fb);       
        coeffTable2{idx+1,'pValue'}=pval_group2(rr,fb);
        coeffTable2{idx,'FDRpValue'}=padj_conflict2(rr,fb);       
        coeffTable2{idx+1,'FDRpValue'}=padj_group2(rr,fb);
        idx = idx+2;
    end
end

% Make tables for significant effects
sigInteractionTable = coeffTable1(coeffTable1.Predictor=='Conflict x Group' & coeffTable1.FDRpValue<0.05,:);
sigConflict1Table = coeffTable1(coeffTable1.Predictor=='Conflict' & coeffTable1.FDRpValue<0.05,:);
sigGroup1Table = coeffTable1(coeffTable1.Predictor=='Group' & coeffTable1.FDRpValue<0.05,:);
sigConflict2Table = coeffTable2(coeffTable2.Predictor=='Conflict' & coeffTable2.FDRpValue<0.05,:);
sigGroup2Table = coeffTable2(coeffTable2.Predictor=='Group' & coeffTable2.FDRpValue<0.05,:);

%::::::::::::::::::::::::::Posthoc Code::::::::::::::::::::::::::
% Posthoc contrast tests for significant interactions

%Get models with significant interactions
sigInteractionModels={}; %Row 1 = Models, Row 2 = Regions IDs
for dd = 1:length(F)
sigInteractionModels{1,dd} = M1(padj_conflict_group(:,dd)<0.05,dd);
sigInteractionModels{2,dd} = sigRegionInteractions{dd};
end

%Define Contrasts
EC_Low_Conflict = [1 0 0 0]; 
EC_High_Conflict = [1 1 0 0];  
AD_Low_Conflict = [1 0 1 0];   
AD_High_Conflict = [1 1 1 1];
C1 = AD_High_Conflict-EC_High_Conflict;
C2 = AD_Low_Conflict-EC_Low_Conflict;
C3 = EC_High_Conflict-EC_Low_Conflict;
C4 = AD_High_Conflict-AD_Low_Conflict;

maxRegNum2 = max(cell2mat(cellfun(@length,sigInteractionModels(2,:),'uni',false))); %get # regions in band with most significant regions
Fval_posthoc = nan(length(F),maxRegNum2,4);
pval_posthoc=nan(length(F),maxRegNum2,4);
contrast_coef_posthoc=nan(length(F),maxRegNum2,4);
DF1_posthoc = nan(length(F),maxRegNum2,4);
DF2_posthoc = nan(length(F),maxRegNum2,4);

%Run posthoc tests of coefficient contrasts
for fb = 1:length(F)
    mdls = sigInteractionModels{1,fb};
    R_int = sigInteractionModels{2,fb};
    for rr = 1:length(R_int)
        [pval_posthoc(fb,rr,1),Fval_posthoc(fb,rr,1),DF1_posthoc(fb,rr,1),DF2_posthoc(fb,rr,1)]=coefTest(mdls{rr},C1);
        [pval_posthoc(fb,rr,2),Fval_posthoc(fb,rr,2),DF1_posthoc(fb,rr,2),DF2_posthoc(fb,rr,2)]=coefTest(mdls{rr},C2);
        [pval_posthoc(fb,rr,3),Fval_posthoc(fb,rr,3),DF1_posthoc(fb,rr,3),DF2_posthoc(fb,rr,3)]=coefTest(mdls{rr},C3);
        [pval_posthoc(fb,rr,4),Fval_posthoc(fb,rr,4),DF1_posthoc(fb,rr,4),DF2_posthoc(fb,rr,4)]=coefTest(mdls{rr},C4);
        contrast_coef_posthoc(fb,rr,1) = C1*fixedEffects(mdls{rr});
        contrast_coef_posthoc(fb,rr,2) = C2*fixedEffects(mdls{rr});
        contrast_coef_posthoc(fb,rr,3) = C3*fixedEffects(mdls{rr});
        contrast_coef_posthoc(fb,rr,4) = C4*fixedEffects(mdls{rr});
    end
end
[~,critp_posthoc,~,padj_posthoc]=fdr_bh(pval_posthoc);

%Build table of posthoc comparisons:
posthocTable = table('Size',[sum(cellfun(@(x) size(x,2),sigInteractionModels(2,:)))*4 9],'VariableTypes',{'string','string','string','double','double','double','double','double','double'},'VariableNames', {'Region','Frequency Band','Contrast','Coefficient','Fstat','DF1','DF2','pValue','FDRpValue'});
idx = 1;
for fb=1:length(F)
    mdls = sigInteractionModels{1,fb};
    R_int = sigInteractionModels{2,fb};
    if ~isempty(mdls)
        for rr = 1:length(R_int)
            posthocTable{idx:idx+3,'Region'}= RegionNames(R==string(R_int{rr}));
            posthocTable{idx:idx+3,'Frequency Band'}=string(FrequencyBands{fb});
            posthocTable{idx,'Contrast'}="AD-EC High Conflict";
            posthocTable{idx+1,'Contrast'}="AD-EC Low Conflict";
            posthocTable{idx+2,'Contrast'}="EC High-Low";
            posthocTable{idx+3,'Contrast'}="AD High-Low";
            posthocTable{idx:idx+3,'Coefficient'}=squeeze(contrast_coef_posthoc(fb,rr,:));
            posthocTable{idx:idx+3,'Fstat'}=squeeze(Fval_posthoc(fb,rr,:));
            posthocTable{idx:idx+3,'DF1'}=squeeze(DF1_posthoc(fb,rr,:));
            posthocTable{idx:idx+3,'DF2'}=squeeze(DF2_posthoc(fb,rr,:));
            posthocTable{idx:idx+3,'pValue'}=squeeze(pval_posthoc(fb,rr,:));
            posthocTable{idx:idx+3,'FDRpValue'}=squeeze(padj_posthoc(fb,rr,:));
            idx = idx+4;
         end
    end
end

%Tables for signficant contrasts
sigGroupContrast = posthocTable((posthocTable.Contrast=="AD-EC High Conflict"|posthocTable.Contrast=="AD-EC Low Conflict") & posthocTable.FDRpValue<0.05,:);
sigConflictContrast = posthocTable((posthocTable.Contrast=="AD High-Low"|posthocTable.Contrast=="EC High-Low") & posthocTable.FDRpValue<0.05,:);

%Define regions to be included in following Coherence analyses
coh_R = cell(1,length(F));
coh_RegionNames = cell(1,length(F));
for ii = 1:length(F)
    regIDs=string.empty;
    regNames={};
    for jj = 1:length(R)
        if ismember(R(jj),sigInteractionModels{2,ii}) || ismember(R(jj),sigRegionConflict2{ii}) || ismember(R(jj),sigRegionGroup2{ii}) 
            regIDs(end+1)=R(jj);
            regNames(end+1)=RegionNames(jj);
        end
    end
    coh_R{ii}=regIDs;
    coh_RegionNames{ii}=regNames;
end

save('INSERT SAVE PATH HERE','LogPow','Conflict','Group','RT','Subject','Electrode','FreqBand','Region','Trial','M1','M2','coeffTable1','coeffTable2','conflictGLMRegions','padj_conflict1','padj_conflict2','padj_conflict_group','padj_group1','padj_group2','coef_conflict1','coef_group1','coef_conflict_group','coef_conflict2','coef_group2','tstat_conflict1','tstat_group1','tstat_conflict_group','tstat_conflict2','tstat_group2','sigInteractionModels','contrast_coef_posthoc','Fval_posthoc','padj_posthoc','coh_R','coh_RegionNames','sigConflict1Table','sigGroup1Table','sigConflict2Table','sigGroup2Table','sigInteractionTable','posthocTable','sigGroupContrast','sigConflictContrast','-v7.3');
%% Check number of trials/Electrodes included in GLMEs
eleccounts= zeros(length(map_ad),length(R),3);
trialcounts=zeros(length(map_ad),length(R),3);
for ff=1:length(map_ad)
    for fb=1:3
        for rr = 1:length(R)
            rid = find(Subject==['P',num2str(map_ad(ff,1))] & FreqBand == fb & Region == R(rr) & ~isnan(LogPow));
            if ~isempty(rid)
                [elecs, ia, ic]= unique(Electrode(rid));
                a_counts = accumarray(ic,1);
                numtrials = max(a_counts);
                numelec = numel(elecs);
                eleccounts(ff,rr,fb) = numelec;
                trialcounts(ff,rr,fb)=numtrials;
            end
        end
    end
end
trialsum = squeeze(sum(trialcounts,1));
trialtotal=max(trialsum,[],'all')
elecsum = squeeze(sum(eleccounts,1));
elecstotal=sum(elecsum(:,1))