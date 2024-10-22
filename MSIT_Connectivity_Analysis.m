%% Coherence/PLV Estimation: All channel-channel connections
clear;
addpath C:\Users\shekaraa\Documents\MATLAB\fieldtrip-20230613
cd D:/Data

%Load participant list w/categorization
load('state_map.mat','map_ad')

%Select regions for coherence/PLV estimation
R=["L1","R1","L2","R2","L3","R3","L4","R4","L6","R6","L9","R9","L11","R11","L13","R13","L14","R14"];
FrequencyBands = [{'Theta'},{'Alpha'},{'Beta'},{'Gamma'},{'High Gamma'}];
F=[4,8;8,15;15,30;30,55;70,110];

for ss = 1:length(map_ad)
    filename=['P',num2str(map_ad(ss,1)),'.mat'];
    subject =['P',num2str(map_ad(ss,1))];
    load(['D:\Data\Coherence Data\ChannelCoh_allbands_CorrectRespOnly\',filename],'ft_data3_clean','TrialDet','Parcellation_Sided') %Load fourier spectrum created from coherence analysis

    trial_low=find(~isnan(TrialDet(:,12)) & TrialDet(:,12)<2 & TrialDet(:,14)==min(TrialDet(:,14)) & TrialDet(:,27)==1); %Looks for response, conflict type, and Accuracy, exlude RT>2
    trial_high=find(~isnan(TrialDet(:,12)) & TrialDet(:,12)<2 & TrialDet(:,14)==max(TrialDet(:,14)) & TrialDet(:,27)==1);

    %Combine vlPFC, OFC, and dlPFC
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
                Parcellation_Sided{bb} = 'L1';%combine left dlFPC
            case "R3"
                Parcellation_Sided{bb} = 'R1';%combine right dlPFC
        end
    end

    %Include ROIs only
    chans = ft_data3_clean.label(find(ismember(Parcellation_Sided,R)));
    regionIDs = (Parcellation_Sided(ismember(Parcellation_Sided,R)));

    %Calculate Fourier transform
    cfg = [];
    cfg.channel = chans;
    cfg.output = 'fourier';
    cfg.keeptrials = 'yes';
    cfg.method = 'wavelet';      
    cfg.pad='nextpow2';
    cfg.foi = [4:110];
    cfg.toi = [0:(1/512*24):2];
    ft_freq_fourier = ft_freqanalysis(cfg,ft_data3_clean);

    %Set artifactual trials marked with NaNs to zero
    ft_freq_fourier.fourierspctrm(isnan(ft_freq_fourier.fourierspctrm))=0;

    %Select band to calculate coherence
    ft_coh_bands = cell(5,2);
    ft_plv_bands=cell(5,2);
    for ii = 1:length(F)
        cfg=[];
        cfg.frequency = F(ii,:);
        freqBand_fourier = ft_selectdata(cfg,ft_freq_fourier);
        
        %Calculate coherence for low and high conflict
        cfg = [];
        cfg.method = 'coh';
        cfg.complex = 'complex';
    
        cfg.trials=trial_low;
        coh_low = ft_connectivityanalysis(cfg, freqBand_fourier);
    
        cfg.trials = trial_high;
        coh_high = ft_connectivityanalysis(cfg, freqBand_fourier);
        
        ft_coh_bands{ii,1}=coh_low;
        ft_coh_bands{ii,2}=coh_high;
    
        %Calculate PLV for low and high conflict
        cfg = [];
        cfg.method = 'plv';
        cfg.complex = 'complex';
    
        cfg.trials=trial_low;
        plv_low = ft_connectivityanalysis(cfg, freqBand_fourier);
    
        cfg.trials = trial_high;
        plv_high = ft_connectivityanalysis(cfg, freqBand_fourier);
    
        ft_plv_bands{ii,1}=plv_low;
        ft_plv_bands{ii,2}=plv_high;
    end
    save(['D:\Data\Coherence Data\ChannelCohBL_allbands_CorrectRespOnly\',filename],'ft_freq_fourier','ft_coh_bands','regionIDs','-v7.3');
    clear ft_coh_bands 
    save(['D:\Data\PLV Data\PLVBL_CorrectRespOnly\',filename],'ft_freq_fourier','ft_plv_bands','regionIDs','-v7.3');
    clear ft_plv_bands ft_freq_fourier ft_data3_clean
end
%% Coherence GLMEs for conflict-encoding regions
clear all
%Load participant list w/categorization
load('state_map.mat','map_ad')
%Load conflict-encoding regions/bands from power GLMEs
load('D:\Data\GLME Results\Combined dlPFC\dACC Only\Power_ConflictxGroup_Canonical_GLME_CorrectRespOnly_dACCOnly_ExactRT_082924', 'coh_R','coh_RegionNames')
R = coh_R;RegionNames=coh_RegionNames;
F=[4,8;8,15;15,30;30,55;70,110];

for ii = 1:length(R)
regcmb{ii} = triu(ones(length(R{ii}),length(R{ii})),1); %Index of connections to consider excluding reciprocal ones
end

%Get median RTs
RT_High = [];
RT_Low = [];
for dd = 1:length(map_ad)
    filename=['P',num2str(map_ad(dd,1)),'.mat'];
    load(filename,'TrialDet');
    Trials_C=find(~isnan(TrialDet(:,12)) & TrialDet(:,12)<2 & TrialDet(:,14)==min(TrialDet(:,14)) & TrialDet(:,27)==1); %Looks for response, conflict type, and Accuracy
    Trials_I=find(~isnan(TrialDet(:,12)) & TrialDet(:,12)<2 & TrialDet(:,14)==max(TrialDet(:,14)) & TrialDet(:,27)==1);
    RT_Low=[RT_Low;median(TrialDet(Trials_C,12))];
    RT_High=[RT_High;median(TrialDet(Trials_I,12))];
end

TwinLow=RT_Low;
TwinHigh=RT_High;

Conflict = [];
Subject=[];
Group = [];
Coh=[];
FreqBand = [];
Conn={};

% Average Coherence for each channel pair per region-region connection and build data table
for ff=1:length(map_ad)
    filename=['P',num2str(map_ad(ff,1)),'.mat'];
    load(['D:\Data\Coherence Data\ChannelCoh_allbands_CorrectRespOnly\',filename],'ft_coh_bands','regionIDs');
    Mstatus=map_ad(ff,2);
        
    for fb=1:length(F)
        coh_low = ft_coh_bands{fb,1};
        coh_high = ft_coh_bands{fb,2};
        Id_time1_Low=find(coh_low.time>0.1 & coh_low.time<=TwinLow(ff));
        Id_time1_High=find(coh_high.time>0.1 & coh_high.time<=TwinHigh(ff));
        for mm=1:length(R{fb})
            for nn = 1:length(R{fb})
                if regcmb{fb}(mm,nn) == 1
                    chans1 = find(regionIDs==R{fb}(mm));
                    chans2 = find(regionIDs == R{fb}(nn));
                    cohlow=nanmean(nanmean(abs(coh_low.cohspctrm(chans1,chans2,:,Id_time1_Low)),4),3); % select chans, average across freq and time
                    cohlow = reshape(cohlow,[],1);
                    cohhigh=nanmean(nanmean(abs(coh_high.cohspctrm(chans1,chans2,:,Id_time1_High)),4),3);
                    cohhigh = reshape(cohhigh,[],1);
                    
                    Coh=[Coh;cohlow;cohhigh];
                    Conflict = [Conflict;zeros(size(cohlow,1),1);ones(size(cohhigh,1),1)];
                    regs = cell(size(cohlow,1)+size(cohhigh,1),1);
                    regs(:) ={[RegionNames{fb}{mm},'-',RegionNames{fb}{nn}]};
                    Conn = [Conn;regs];
                    Group = [Group;repmat(map_ad(ff,2),size(cohlow,1)+size(cohhigh,1),1)];
                    sub=strings(size(cohlow,1)+size(cohhigh,1),1);
                    sub(:)=['P',num2str(map_ad(ff,1))];
                    Subject=[Subject;sub];
                    FreqBand=[FreqBand;fb.*ones(size(cohlow,1)+size(cohhigh,1),1)];
                end
            end
        end
    end
end

%Log transform coherence values
Coh = log(Coh);

% :::::::::::::::GLME CODE::::::::::::::::

for ii = 1:length(F)
    M1{ii}= cell(size(regcmb{ii}));
end
maxReg = max(cell2mat(cellfun(@length,R,'uni',false)));
coef_conflict1=nan(length(F),maxReg,maxReg);
coef_group1=nan(length(F),maxReg,maxReg);
coef_conflict_group=nan(length(F),maxReg,maxReg);
tstat_conflict1=nan(length(F),maxReg,maxReg);
tstat_group1=nan(length(F),maxReg,maxReg);
tstat_conflict_group=nan(length(F),maxReg,maxReg);
pval_conflict1=nan(length(F),maxReg,maxReg);
pval_group1=nan(length(F),maxReg,maxReg);
pval_conflict_group=nan(length(F),maxReg,maxReg);

%Run full GLMEs
for fb=1:length(F)
    for mm = 1:length(R{fb})
        for nn = 1:length(R{fb})
            rid = find(ismember(Conn,[RegionNames{fb}{mm},'-',RegionNames{fb}{nn}]) & FreqBand ==fb);
            if ~isempty(rid)
                dataTable = table(Coh(rid),categorical(Conflict(rid)),categorical(Group(rid)),Subject(rid),'VariableNames',{'Coh','Conflict','Group','Subject'});
                M1{fb}{mm,nn}=fitglme(dataTable,'Coh ~ Conflict*Group + (1|Subject)','Link','identity');
                coef_conflict1(fb,mm,nn)=M1{fb}{mm,nn}.Coefficients.Estimate(2);
                coef_group1(fb,mm,nn)=M1{fb}{mm,nn}.Coefficients.Estimate(3);
                coef_conflict_group(fb,mm,nn)=M1{fb}{mm,nn}.Coefficients.Estimate(4);
                tstat_conflict1(fb,mm,nn)=M1{fb}{mm,nn}.Coefficients.tStat(2);
                tstat_group1(fb,mm,nn)=M1{fb}{mm,nn}.Coefficients.tStat(3);
                tstat_conflict_group(fb,mm,nn)=M1{fb}{mm,nn}.Coefficients.tStat(4);
                pval_conflict1(fb,mm,nn)=M1{fb}{mm,nn}.Coefficients.pValue(2);
                pval_group1(fb,mm,nn)=M1{fb}{mm,nn}.Coefficients.pValue(3);
                pval_conflict_group(fb,mm,nn)=M1{fb}{mm,nn}.Coefficients.pValue(4);
            else
                M1{fb}{mm,nn}=NaN;
            end        
        end
    end
end

%FDR-correct p values
[~,~,~,padj_conflict1]=fdr_bh(pval_conflict1);
[~,~,~,padj_group1]=fdr_bh(pval_group1);
[~,~,~,padj_conflict_group]=fdr_bh(pval_conflict_group);

%Get connection list
for fb=1:length(F)
    connections{fb} =cell(length(R{fb}),length(R{fb}));
    for mm = 1:length(R{fb})
        for nn = 1:length(R{fb})
            if regcmb{fb}(mm,nn) == 1
                connections{fb}{mm,nn} = [RegionNames{fb}{mm},'-',RegionNames{fb}{nn}];
            end
        end
    end
end

%Reduced GLMEs for regions without interactions
for fb = 1:length(F)
sigConnectionInteractions{fb} = connections{fb}(squeeze(padj_conflict_group(fb,1:size(M1{fb},1),1:size(M1{fb},2))<0.05)); %Region w/sig interaction
sigConnectionConflict1{fb} = connections{fb}(squeeze(padj_conflict1(fb,1:size(M1{fb},1),1:size(M1{fb},2))<0.05));%Regions w/significant effect of conflict
sigConnectionGroup1{fb} = connections{fb}(squeeze(padj_group1(fb,1:size(M1{fb},1),1:size(M1{fb},2))<0.05));%Regions w/significant effect of group

%Select only significant regions w/o interactions
regcmb_conflict{fb}=  regcmb{fb}-(squeeze(padj_conflict_group(fb,1:size(regcmb{fb},1),1:size(regcmb{fb},2)))<0.05);
M2{fb}= cell(size(regcmb{fb}));
end
coef_conflict2=nan(length(F),maxReg,maxReg);
pval_conflict2=nan(length(F),maxReg,maxReg);
tstat_conflict2=nan(length(F),maxReg,maxReg);
coef_group2=nan(length(F),maxReg,maxReg);
pval_group2=nan(length(F),maxReg,maxReg);
tstat_group2=nan(length(F),maxReg,maxReg);

%Run reduced GLMEs
for fb=1:length(F)
    for mm = 1:length(R{fb})
        for nn = 1:length(R{fb})
            rid = find(ismember(Conn,[RegionNames{fb}{mm},'-',RegionNames{fb}{nn}]) & FreqBand==fb);
            if ~isempty(rid) && regcmb_conflict{fb}(mm,nn)==1
                dataTable = table(Coh(rid),categorical(Conflict(rid)),categorical(Group(rid)),Subject(rid),'VariableNames',{'Coh','Conflict','Group','Subject'});
                M2{fb}{mm,nn}=fitglme(dataTable,'Coh ~ Conflict + Group + (1|Subject)','Link','identity');
                coef_conflict2(fb,mm,nn) = M2{fb}{mm,nn}.Coefficients.Estimate(2);
                tstat_conflict2(fb,mm,nn)=M2{fb}{mm,nn}.Coefficients.tStat(2);
                pval_conflict2(fb,mm,nn) = M2{fb}{mm,nn}.Coefficients.pValue(2);
                coef_group2(fb,mm,nn) = M2{fb}{mm,nn}.Coefficients.Estimate(3);
                tstat_group2(fb,mm,nn)=M2{fb}{mm,nn}.Coefficients.tStat(3);
                pval_group2(fb,mm,nn) = M2{fb}{mm,nn}.Coefficients.pValue(3);
            else
                M2{fb}{mm,nn}=NaN;       
                coef_conflict2(fb,mm,nn) = NaN;
                tstat_conflict2(fb,mm,nn)= NaN;
                pval_conflict2(fb,mm,nn) = NaN;
                coef_group2(fb,mm,nn) = NaN;
                tstat_group2(fb,mm,nn)= NaN;
                pval_group2(fb,mm,nn) = NaN;
            end        
        end
    end
end

[~,~,~,padj_conflict2]=fdr_bh(pval_conflict2);
[~,~,~,padj_group2]=fdr_bh(pval_group2);

%Create Data tables
%Full GLMEs
coeffTable1 = table('Size',[sum(cellfun(@(x) sum(x==1,'all'),regcmb))*3 7],'VariableTypes',{'string','string','string','double','double','double','double'},'VariableNames', {'Connection','Frequency Band','Predictor','Coefficient','tstat','pValue','FDRpValue'});
%Reduced GLMEs
coeffTable2 = table('Size',[sum(cellfun(@(x) sum(x==1,'all'),regcmb_conflict))*2 7],'VariableTypes',{'string','string','string','double','double','double','double'},'VariableNames', {'Connection','Frequency Band','Predictor','Coefficient','tstat','pValue','FDRpValue'});
FrequencyBands = [{'Theta'},{'Alpha'},{'Beta'},{'Gamma'},{'High Gamma'}];
idx = 1;
for fb=1:length(F)
    for mm = 1:length(M1{fb})
        for nn = 1:length(M1{fb})
            if isobject(M1{fb}{mm,nn})
            coeffTable1{idx,'Connection'}= string(connections{fb}{mm,nn});
            coeffTable1{idx+1,'Connection'}= string(connections{fb}{mm,nn});
            coeffTable1{idx+2,'Connection'}= string(connections{fb}{mm,nn});
            coeffTable1{idx,'Frequency Band'}=string(FrequencyBands{fb});
            coeffTable1{idx+1,'Frequency Band'}=string(FrequencyBands{fb});
            coeffTable1{idx+2,'Frequency Band'}=string(FrequencyBands{fb});
            coeffTable1{idx,'Predictor'}="Conflict";
            coeffTable1{idx+1,'Predictor'}="Group";
            coeffTable1{idx+2,'Predictor'}="Conflict x Group";
            coeffTable1{idx,'Coefficient'}=coef_conflict1(fb,mm,nn);
            coeffTable1{idx+1,'Coefficient'}=coef_group1(fb,mm,nn);
            coeffTable1{idx+2,'Coefficient'}=coef_conflict_group(fb,mm,nn);
            coeffTable1{idx,'tstat'}=tstat_conflict1(fb,mm,nn);
            coeffTable1{idx+1,'tstat'}=tstat_group1(fb,mm,nn);
            coeffTable1{idx+2,'tstat'}=tstat_conflict_group(fb,mm,nn);
            coeffTable1{idx,'pValue'}=pval_conflict1(fb,mm,nn);
            coeffTable1{idx+1,'pValue'}=pval_group1(fb,mm,nn);
            coeffTable1{idx+2,'pValue'}=pval_conflict_group(fb,mm,nn);
            coeffTable1{idx,'FDRpValue'}=padj_conflict1(fb,mm,nn);
            coeffTable1{idx+1,'FDRpValue'}=padj_group1(fb,mm,nn);
            coeffTable1{idx+2,'FDRpValue'}=padj_conflict_group(fb,mm,nn);
            idx = idx+3;
            end
        end
    end
end

idx = 1;
for fb=1:length(F)
    for mm = 1:length(M2{fb})
        for nn = 1:length(M2{fb})
            if isobject(M2{fb}{mm,nn})
                coeffTable2{idx,'Connection'}= string(connections{fb}{mm,nn});
                coeffTable2{idx+1,'Connection'}= string(connections{fb}{mm,nn});
                coeffTable2{idx,'Frequency Band'}=string(FrequencyBands{fb});
                coeffTable2{idx+1,'Frequency Band'}=string(FrequencyBands{fb});
                coeffTable2{idx,'Predictor'}="Conflict";
                coeffTable2{idx+1,'Predictor'}="Group";
                coeffTable2{idx,'Coefficient'}=coef_conflict2(fb,mm,nn);
                coeffTable2{idx+1,'Coefficient'}=coef_group2(fb,mm,nn);
                coeffTable2{idx,'tstat'}=tstat_conflict2(fb,mm,nn);
                coeffTable2{idx+1,'tstat'}=tstat_group2(fb,mm,nn);
                coeffTable2{idx,'pValue'}=pval_conflict2(fb,mm,nn);
                coeffTable2{idx+1,'pValue'}=pval_group2(fb,mm,nn);
                coeffTable2{idx,'FDRpValue'}=padj_conflict2(fb,mm,nn);
                coeffTable2{idx+1,'FDRpValue'}=padj_group2(fb,mm,nn);
                idx = idx + 2;
            end
        end
    end
end

%:::::::::::::::POSTHOC ANALYSIS CODE::::::::::::::::

% Posthoc Coefficient Test for significant interactions
sigInteractionModels={}; %Row 1 = Models, Row 2 = Connections
for fb = 1:length(F)
    sigInteractionModels{1,fb} = M1{fb}(squeeze(padj_conflict_group(fb,1:size(M1{fb},1),1:size(M1{fb},2)))<0.05);
    sigInteractionModels{2,fb} = sigConnectionInteractions{fb};
end

%Contrasts
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
for fb = 1:length(F)
    mdls = sigInteractionModels{1,fb};
    R_int = sigInteractionModels{2,fb};
    if ~isempty(mdls)
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
end
[~,~,~,padj_posthoc]=fdr_bh(pval_posthoc);

sigInteractionTable = coeffTable1(coeffTable1.Predictor=='Conflict x Group' & coeffTable1.FDRpValue<0.05,:);
sigConflict1Table = coeffTable1(coeffTable1.Predictor=='Conflict' & coeffTable1.FDRpValue<0.05,:);
sigGroup1Table = coeffTable1(coeffTable1.Predictor=='Group' & coeffTable1.FDRpValue<0.05,:);
sigConflict2Table = coeffTable2(coeffTable2.Predictor=='Conflict' & coeffTable2.FDRpValue<0.05,:);
sigGroup2Table = coeffTable2(coeffTable2.Predictor=='Group' & coeffTable2.FDRpValue<0.05,:);


%Build table of posthoc connections:
posthocTable = table('Size',[sum(cellfun(@(x) size(x,1),sigInteractionModels(2,:)))*4 9],'VariableTypes',{'string','string','string','double','double','double','double','double','double'},'VariableNames', {'Connection','Frequency Band','Contrast','Coefficient','Fstat','DF1','DF2','pValue','FDRpValue'});
idx = 1;
for fb=1:length(F)
    mdls = sigInteractionModels{1,fb};
    R_int = sigInteractionModels{2,fb};
    if ~isempty(mdls)
        for rr = 1:length(R_int)
            posthocTable{idx:idx+3,'Connection'}= string(R_int{rr});
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

%Make tables for significant contrasts
sigGroupContrast = posthocTable((posthocTable.Contrast=="AD-EC High Conflict"|posthocTable.Contrast=="AD-EC Low Conflict") & posthocTable.FDRpValue<0.05,:);
sigConflictContrast = posthocTable((posthocTable.Contrast=="AD High-Low"|posthocTable.Contrast=="EC High-Low") & posthocTable.FDRpValue<0.05,:);

cohDataTable = table(Coh,Conflict,Group,Subject,Conn,FreqBand,'VariableNames',{'Coh','Conflict','Group','Subject','Connection','FreqBand'});
%save('D:\Data\GLME Results\Combined dlPFC\dACC Only\Coherence_ConflictxGroup_Canonical_GLME_CorrectRespOnly_dACCOnly_ExactRTRegs_082924','Coh','Conflict','Conn','Group','Subject','FreqBand','cohDataTable','F','R','RegionNames','M1','M2','coeffTable1','coeffTable2','sigConnectionConflict1','sigConnectionGroup1','coef_conflict1','coef_group1','coef_conflict_group','coef_conflict2','coef_group2','tstat_conflict1','tstat_group1','tstat_conflict_group','tstat_conflict2','tstat_group2','padj_conflict1','padj_conflict2','padj_group1','padj_group2','padj_conflict_group','Fval_posthoc','padj_posthoc','sigInteractionModels','contrast_coef_posthoc','posthocTable','sigInteractionTable','sigConflict1Table','sigConflict2Table','sigGroup1Table','sigGroup2Table','sigConflictContrast','sigGroupContrast','-v7.3')

%% PLV GLMEs for significant coherence interaction effects
clear all
cd D:/Data
load('state_map_082924.mat','map_ad')

%Get significant coherence interaction models
load('D:\Data\GLME Results\Combined dlPFC\dACC Only\Coherence_ConflictxGroup_Canonical_GLME_CorrectRespOnly_dACCOnly_ExactRTRegs_082924.mat', 'sigInteractionModels')
R=["L1","R1","L2","R2","L4","R4","L6","R6","L9","R9","L11","R11","L13","R13","L14","R14"];
RegionNames=[{'L dlPFC'},{'R dlPFC'},{'L dmPFC'},{'R dmPFC'},{'L OFC'},{'R OFC'},{'L vlPFC'},{'R vlPFC'},{'L LTL'},{'R LTL'},{'L dACC'},{'R dACC'},{'L AMY'},{'R AMY'},{'L HC'},{'R HC'}];
sigIntConns = sigInteractionModels(2,:);
F=[4,8;8,15;15,30;30,55;70,110];

%Get median RTs
RT_High = [];
RT_Low = [];
for dd = 1:length(map_ad)
    filename=['P',num2str(map_ad(dd,1)),'.mat'];
    load(filename,'TrialDet');
    Trials_C=find(~isnan(TrialDet(:,12)) & TrialDet(:,12)<2 & TrialDet(:,14)==min(TrialDet(:,14)) & TrialDet(:,27)==1); %Looks for response, conflict type, and Accuracy
    Trials_I=find(~isnan(TrialDet(:,12)) & TrialDet(:,12)<2 & TrialDet(:,14)==max(TrialDet(:,14)) & TrialDet(:,27)==1);
    RT_Low=[RT_Low;median(TrialDet(Trials_C,12))];
    RT_High=[RT_High;median(TrialDet(Trials_I,12))];
end

TwinLow=RT_Low;
TwinHigh=RT_High;

Conflict = [];
Subject=[];
Group = [];
Plv=[];
Amp=[];
FreqBand = [];
Conn={};

% Average PLV and ampcorr for each channel pair per region-region connection and build data table
for ff=1:length(map_ad)
    filename=['P',num2str(map_ad(ff,1)),'.mat'];
    load(['D:\Data\PLV Data\PLV_CorrectRespOnly\',filename],'ft_plv_bands','regionIDs');
    Mstatus=map_ad(ff,2);
        
    for fb=1:length(F)
        plv_low = ft_plv_bands{fb,1};
        plv_high = ft_plv_bands{fb,2};
        Id_time1_Low=find(plv_low.time>0.1 & plv_low.time<=TwinLow(ff));
        Id_time1_High=find(plv_high.time>0.1 & plv_high.time<=TwinHigh(ff));
        
        for rr = 1:size(sigIntConns{fb},1)
            sigregs = split(sigIntConns{fb}{rr},'-');
            R1 = R(ismember(RegionNames,sigregs{1}));
            R2 = R(ismember(RegionNames,sigregs{2}));
            chans1 = find(regionIDs == R1);
            chans2 = find(regionIDs == R2);
            
            % Get log-transformed PLV
            plvlow=log(nanmean(nanmean(abs(plv_low.plvspctrm(chans1,chans2,:,Id_time1_Low)),4),3)); % select chans, average across freq and time
            plvlow = reshape(plvlow,[],1);
            plvhigh=log(nanmean(nanmean(abs(plv_high.plvspctrm(chans1,chans2,:,Id_time1_High)),4),3));
            plvhigh = reshape(plvhigh,[],1);
            Plv=[Plv;plvlow;plvhigh];

            Conflict = [Conflict;zeros(size(plvlow,1),1);ones(size(plvhigh,1),1)];
            regs = cell(size(plvlow,1)+size(plvhigh,1),1);
            regs(:) =sigIntConns{fb}(rr);
            Conn = [Conn;regs];
            Group = [Group;repmat(map_ad(ff,2),size(plvlow,1)+size(plvhigh,1),1)];
            sub=strings(size(plvlow,1)+size(plvhigh,1),1);
            sub(:)=['P',num2str(map_ad(ff,1))];
            Subject=[Subject;sub];
            FreqBand=[FreqBand;fb.*ones(size(plvlow,1)+size(plvhigh,1),1)];
        end
    end
end

%:::::::::::::::::::GLME CODE:::::::::::::::::::

for ii = 1:length(F)
    M1{ii}= cell(size(sigIntConns{ii}));
end
maxReg = max(cell2mat(cellfun(@length,sigIntConns,'uni',false)));
coef_conflict1=nan(maxReg,length(F));
coef_group1=nan(maxReg,length(F));
coef_conflict_group=nan(maxReg,length(F));
tstat_conflict1=nan(maxReg,length(F));
tstat_group1=nan(maxReg,length(F));
tstat_conflict_group=nan(maxReg,length(F));
pval_conflict1=nan(maxReg,length(F));
pval_group1=nan(maxReg,length(F));
pval_conflict_group=nan(maxReg,length(F));

%Run full GLMEs for PLV
for fb=1:length(F)
    for rr = 1:size(sigIntConns{fb},1)
        rid = find(ismember(Conn,sigIntConns{fb}{rr}) & FreqBand==fb);
        if ~isempty(rid)
            dataTable = table(Plv(rid),categorical(Conflict(rid)),categorical(Group(rid)),Subject(rid),'VariableNames',{'Plv','Conflict','Group','Subject'});
            M1{fb}{rr}=fitglme(dataTable,'Plv ~ Conflict*Group + (1|Subject)','Link','identity');
            coef_conflict1(rr,fb)=M1{fb}{rr}.Coefficients.Estimate(2);
            coef_group1(rr,fb)=M1{fb}{rr}.Coefficients.Estimate(3);
            coef_conflict_group(rr,fb)=M1{fb}{rr}.Coefficients.Estimate(4);
            tstat_conflict1(rr,fb)=M1{fb}{rr}.Coefficients.tStat(2);
            tstat_group1(rr,fb)=M1{fb}{rr}.Coefficients.tStat(3);
            tstat_conflict_group(rr,fb)=M1{fb}{rr}.Coefficients.tStat(4);
            pval_conflict1(rr,fb)=M1{fb}{rr}.Coefficients.pValue(2);
            pval_group1(rr,fb)=M1{fb}{rr}.Coefficients.pValue(3);
            pval_conflict_group(rr,fb)=M1{fb}{rr}.Coefficients.pValue(4);
        else
            M1{fb}{rr}=NaN;
            coef_conflict1(rr,fb)=NaN;
            coef_group1(rr,fb)=NaN;
            coef_conflict_group(rr,fb)=NaN;
            tstat_conflict1(rr,fb)=NaN;
            tstat_group1(rr,fb)=NaN;
            tstat_conflict_group(rr,fb)=NaN;
            pval_conflict1(rr,fb)=NaN;
            pval_group1(rr,fb)=NaN;
            pval_conflict_group(rr,fb)=NaN;
        end        
    end
end

%FDR-correct p values
[~,critp_conflict1,~,padj_conflict1]=fdr_bh(pval_conflict1);
[~,critp_group1,~,padj_group1]=fdr_bh(pval_group1);
[~,critp_conflict_group,~,padj_conflict_group]=fdr_bh(pval_conflict_group);

%Get regions without interactions for reduced GLMEs
conflictGLMconns=cell(1,length(F));
for fb = 1:length(F)
sigConnectionInteractions{fb} = sigIntConns{fb}(padj_conflict_group(1:size(M1{fb},1),fb)<0.05); %Region w/sig interaction
sigConnectionConflict1{fb} = sigIntConns{fb}(padj_conflict1(1:size(M1{fb},1),fb)<0.05);%Regions w/significant effect of conflict
sigConnectionGroup1{fb} = sigIntConns{fb}(padj_group1(1:size(M1{fb},1),fb)<0.05);%Regions w/significant effect of group
%Select only significant regions w/o interactions
if ~isequal(sigIntConns{fb},sigConnectionInteractions{fb})
    conflictGLMconns{fb}= sigIntConns{fb}(~ismember(sigIntConns{fb},sigConnectionInteractions{fb}));
end
M2{fb}=  cell(size(conflictGLMconns{fb},1),1);
end

coef_conflict2=nan(maxReg,length(F));
coef_group2=nan(maxReg,length(F));
tstat_conflict2=nan(maxReg,length(F));
tstat_group2=nan(maxReg,length(F));
pval_conflict2=nan(maxReg,length(F));
pval_group2=nan(maxReg,length(F));

%Run reduced models
for fb=1:length(F)
    if ~isempty(conflictGLMconns{fb})
        for rr = 1:size(conflictGLMconns{fb},1)
            rid = find(ismember(Conn,conflictGLMconns{fb}{rr}) & FreqBand==fb);
            if ~isempty(rid)
                dataTable = table(Plv(rid),categorical(Conflict(rid)),categorical(Group(rid)),Subject(rid),'VariableNames',{'Plv','Conflict','Group','Subject'});
                M2{fb}{rr}=fitglme(dataTable,'Plv ~ Conflict + Group + (1|Subject)','Link','identity');
                coef_conflict2(rr,fb) = M2{fb}{rr}.Coefficients.Estimate(2);
                coef_group2(rr,fb) = M2{fb}{rr}.Coefficients.Estimate(3);
                tstat_conflict2(rr,fb) = M2{fb}{rr}.Coefficients.tStat(2);
                tstat_group2(rr,fb) = M2{fb}{rr}.Coefficients.tStat(3);       
                pval_conflict2(rr,fb) = M2{fb}{rr}.Coefficients.pValue(2);
                pval_group2(rr,fb) = M2{fb}{rr}.Coefficients.pValue(3);
            else
                M2{fb}{rr}=NaN;
                coef_conflict2(rr,fb) = NaN;
                coef_group2(rr,fb) = NaN;
                tstat_conflict2(rr,fb) = NaN;
                tstat_group2(rr,fb) = NaN;    
                pval_conflict2(rr,fb) = NaN;
                pval_group2(rr,fb) = NaN;
            end        
        end
    end
end

[~,critp_conflict2,~,padj_conflict2]=fdr_bh(pval_conflict2);
[~,critp_group2,~,padj_group2]=fdr_bh(pval_group2);

%Create Data tables
%Full models
coeffTable1 = table('Size',[sum(cell2mat(cellfun(@length,sigIntConns,'uni',false)))*3 7],'VariableTypes',{'string','string','string','double','double','double','double'},'VariableNames', {'Connection','Frequency Band','Predictor','Coefficient','tstat','pValue','FDRpValue'});
%Reduced models
coeffTable2 = table('Size',[sum(cell2mat(cellfun(@length,conflictGLMconns,'uni',false)))*2 7],'VariableTypes',{'string','string','string','double','double','double','double'},'VariableNames', {'Connection','Frequency Band','Predictor','Coefficient','tstat','pValue','FDRpValue'});
FrequencyBands = [{'Theta'},{'Alpha'},{'Beta'},{'Gamma'},{'High Gamma'}];
idx = 1;
for fb=1:length(F)
    for rr = 1:size(M1{fb},1)
        if isobject(M1{fb}{rr})
        coeffTable1{idx,'Connection'}= string(sigIntConns{fb}{rr});
        coeffTable1{idx+1,'Connection'}= string(sigIntConns{fb}{rr});
        coeffTable1{idx+2,'Connection'}= string(sigIntConns{fb}{rr});
        coeffTable1{idx,'Frequency Band'}=string(FrequencyBands{fb});
        coeffTable1{idx+1,'Frequency Band'}=string(FrequencyBands{fb});
        coeffTable1{idx+2,'Frequency Band'}=string(FrequencyBands{fb});
        coeffTable1{idx,'Predictor'}="Conflict";
        coeffTable1{idx+1,'Predictor'}="Group";
        coeffTable1{idx+2,'Predictor'}="Conflict x Group";
        coeffTable1{idx,'Coefficient'}=coef_conflict1(rr,fb);
        coeffTable1{idx+1,'Coefficient'}=coef_group1(rr,fb);
        coeffTable1{idx+2,'Coefficient'}=coef_conflict_group(rr,fb);
        coeffTable1{idx,'tstat'}=tstat_conflict1(rr,fb);
        coeffTable1{idx+1,'tstat'}=tstat_group1(rr,fb);
        coeffTable1{idx+2,'tstat'}=tstat_conflict_group(rr,fb);
        coeffTable1{idx,'pValue'}=pval_conflict1(rr,fb);
        coeffTable1{idx+1,'pValue'}=pval_group1(rr,fb);
        coeffTable1{idx+2,'pValue'}=pval_conflict_group(rr,fb);
        coeffTable1{idx,'FDRpValue'}=padj_conflict1(rr,fb);
        coeffTable1{idx+1,'FDRpValue'}=padj_group1(rr,fb);
        coeffTable1{idx+2,'FDRpValue'}=padj_conflict_group(rr,fb);
        idx = idx+3;
        end
    end
end

idx = 1;
for fb=1:length(F)
    for rr = 1:size(M2{fb},1)
        if isobject(M2{fb}{rr})
            coeffTable2{idx,'Connection'}= string(conflictGLMconns{fb}{rr});
            coeffTable2{idx+1,'Connection'}= string(conflictGLMconns{fb}{rr});
            coeffTable2{idx,'Frequency Band'}=string(FrequencyBands{fb});
            coeffTable2{idx+1,'Frequency Band'}=string(FrequencyBands{fb});
            coeffTable2{idx,'Predictor'}="Conflict";
            coeffTable2{idx+1,'Predictor'}="Group";
            coeffTable2{idx,'Coefficient'}=coef_conflict2(rr,fb);
            coeffTable2{idx+1,'Coefficient'}=coef_group2(rr,fb);
            coeffTable2{idx,'tstat'}=tstat_conflict2(rr,fb);
            coeffTable2{idx+1,'tstat'}=tstat_group2(rr,fb);
            coeffTable2{idx,'pValue'}=pval_conflict2(rr,fb);
            coeffTable2{idx+1,'pValue'}=pval_group2(rr,fb);
            coeffTable2{idx,'FDRpValue'}=padj_conflict2(rr,fb);
            coeffTable2{idx+1,'FDRpValue'}=padj_group2(rr,fb);
            idx = idx + 2;
        end
    end
end

%:::::::::::::::::::POSTHOC ANALYSIS CODE:::::::::::::::::::

% Posthoc Coefficient Test for significant interactions
sigInteractionModels={}; %Row 1 = Models, Row 2 = Connections
for fb = 1:length(F)
    sigInteractionModels{1,fb} = M1{fb}(padj_conflict_group(1:size(M1{fb},1),fb)<0.05);
    sigInteractionModels{2,fb} = sigConnectionInteractions{fb};
end

%Contrasts
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
for fb = 1:length(F)
    mdls = sigInteractionModels{1,fb};
    R_int = sigInteractionModels{2,fb};
    if ~isempty(mdls)
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
end
[~,critp_posthoc,~,padj_posthoc]=fdr_bh(pval_posthoc);

sigInteractionTable = coeffTable1(coeffTable1.Predictor=='Conflict x Group' & coeffTable1.FDRpValue<0.05,:);
sigConflict1Table = coeffTable1(coeffTable1.Predictor=='Conflict' & coeffTable1.FDRpValue<0.05,:);
sigConflict2Table = coeffTable2(coeffTable2.Predictor=='Conflict' & coeffTable2.FDRpValue<0.05,:);
sigGroup1Table = coeffTable1(coeffTable1.Predictor=='Group' & coeffTable1.FDRpValue<0.05,:);
sigGroup2Table = coeffTable2(coeffTable2.Predictor=='Group' & coeffTable2.FDRpValue<0.05,:);

%Build table of posthoc connections:
posthocTable = table('Size',[sum(cellfun(@(x) size(x,1),sigInteractionModels(2,:)))*4 9],'VariableTypes',{'string','string','string','double','double','double','double','double','double'},'VariableNames', {'Connection','Frequency Band','Contrast','Coefficient','Fstat','DF1','DF2','pValue','FDRpValue'});
idx = 1;
for fb=1:length(F)
    mdls = sigInteractionModels{1,fb};
    R_int = sigInteractionModels{2,fb};
    if ~isempty(mdls)
        for rr = 1:length(R_int)
            posthocTable{idx:idx+3,'Connection'}= string(R_int{rr});
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

%Make tables of significant contrasts
sigGroupContrast = posthocTable((posthocTable.Contrast=="AD-EC High Conflict"|posthocTable.Contrast=="AD-EC Low Conflict") & posthocTable.FDRpValue<0.05,:);
sigConflictContrast = posthocTable((posthocTable.Contrast=="AD High-Low"|posthocTable.Contrast=="EC High-Low") & posthocTable.FDRpValue<0.05,:);

PLVTable = table(Plv,categorical(Conflict),categorical(Group),Subject,Conn,FreqBand,'VariableNames',{'Plv','Conflict','Group','Subject','Conn','Freqband'});
%save('PLVFollowup_ConflictxGroup_Canonical_GLME_CorrectRespOnly_dACCOnly_ExactRTRegs_082924','sigConflictContrast','sigGroupContrast','PLVTable','F','sigIntConns','conflictGLMconns','M1','M2','coeffTable1','coeffTable2','sigConnectionConflict1','sigConnectionGroup1','coef_conflict1','coef_group1','coef_conflict_group','coef_conflict2','coef_group2','tstat_conflict1','tstat_group1','tstat_conflict_group','tstat_conflict2','tstat_group2','padj_conflict1','padj_conflict2','padj_group1','padj_group2','padj_conflict_group','Fval_posthoc','padj_posthoc','sigInteractionModels','contrast_coef_posthoc','posthocTable','sigInteractionTable','sigConflict1Table','sigConflict2Table','sigGroup1Table','sigGroup2Table','-v7.3')

%% Q-Q plots to assess residual normality
res={};
figure();
for ii = 1:length(M1)
    for jj = 1:numel(M1{ii})
    if isobject(M1{ii}{jj})
       res{ii,jj}=plotResiduals(M1{ii}{jj},'probability');
        pause;
    else
        res{ii,jj}=NaN;
    end
    end
end
%Power models
figure();
for ii = 1:size(M1,1)
    for jj = 1:size(M1,2)
    if isobject(M1{ii,jj})
       res{ii,jj}=plotResiduals(M1{ii,jj},'probability');
        pause;
    else
        res{ii,jj}=NaN;
    end
    end
end