%% Behavior analysis: Accuracy
clear all
%Load participant list w/categorization
load('state_map.mat','map_ad')

% Build behavior stats table
behaviorTable={};
for ii = 1:length(map_ad)
    filename=['P',num2str(map_ad(ii,1)),'.mat'];
    load(filename,'TrialDet');
    low =find(TrialDet(:,14)==min(TrialDet(:,14)));
    high= find(TrialDet(:,14)==max(TrialDet(:,14)));
    behaviorTable{ii,1}=map_ad(ii,1);
    behaviorTable{ii,2}=map_ad(ii,2);
    behaviorTable{ii,3}=sum(~isnan(TrialDet(:,12)) & TrialDet(:,27)==1); % # Correct Trials
    behaviorTable{ii,4}=sum(~isnan(TrialDet(low,12)) & TrialDet(low,27)==1); % # Correct Low conflict Trials
    behaviorTable{ii,5}=sum(~isnan(TrialDet(high,12)) & TrialDet(high,27)==1); % # Correct High conflict Trials
    behaviorTable{ii,6}=sum(~isnan(TrialDet(:,12)) & TrialDet(:,27)==0); % #Error trials
    behaviorTable{ii,7}=sum(~isnan(TrialDet(low,12)) & TrialDet(low,27)==0); % #Error low conflict trials
    behaviorTable{ii,8}=sum(~isnan(TrialDet(high,12)) & TrialDet(high,27)==0); % #Error high conflict trials
    behaviorTable{ii,9}=sum(isnan(TrialDet(:,12))|TrialDet(:,27)==-1); % #Missing/Omitted trials
    behaviorTable{ii,10}=sum(isnan(TrialDet(low,12))|TrialDet(low,27)==-1); % #Missing/Omistted low conflict trials
    behaviorTable{ii,11}=sum(isnan(TrialDet(high,12))|TrialDet(high,27)==-1); % #Missing/Omistted high conflict trials
    behaviorTable{ii,12}=behaviorTable{ii,3}./size(TrialDet,1); %Total Correct %
    behaviorTable{ii,13}=behaviorTable{ii,4}./length(low); %Low conflict correct %
    behaviorTable{ii,14}=behaviorTable{ii,5}./length(high); %High conflict correct %
    behaviorTable{ii,15}=mean(log(TrialDet(~isnan(TrialDet(:,12)),12))); %mean logRT on correct or error trials
    behaviorTable{ii,16}=mean(log(TrialDet(~isnan(TrialDet(:,12)) & TrialDet(:,14)==min(TrialDet(:,14)),12))); %mean logRT on correct or error low conflict trials
    behaviorTable{ii,17}=mean(log(TrialDet(~isnan(TrialDet(:,12)) & TrialDet(:,14)==max(TrialDet(:,14)),12))); %mean logRT on correct or error high conflict trials
    behaviorTable{ii,18}=median(TrialDet(~isnan(TrialDet(:,12)),12)); % median RT on correct or error trials
    behaviorTable{ii,19}=median(TrialDet(~isnan(TrialDet(:,12)) & TrialDet(:,14)==min(TrialDet(:,14)),12)); % median RT on correct or error low conflict trials
    behaviorTable{ii,20}=median(TrialDet(~isnan(TrialDet(:,12)) & TrialDet(:,14)==max(TrialDet(:,14)),12)); % median RT on correct or error high conflict trials
    behaviorTable{ii,21}=behaviorTable{ii,3}+behaviorTable{ii,6}+behaviorTable{ii,9}; %Total trial #
    behaviorTable{ii,22}=behaviorTable{ii,4}+behaviorTable{ii,7}+behaviorTable{ii,10}; %Total trial #
    behaviorTable{ii,23}=behaviorTable{ii,5}+behaviorTable{ii,8}+behaviorTable{ii,11}; %Total trial #
end
behaviorTable = cell2table(behaviorTable,'VariableNames',{'SubjectID','Group','CorrectTotal','CorrectLow','CorrectHigh','ErrorTotal','ErrorLow','ErrorHigh','MissTotal','MissLow','MissHigh','AccuracyTotal','AccuracyLow','AccuracyHigh','AvgLogRTTotal','AvgLogRTLow','AvgLogRTHigh','MedianRT','MedianRTLow','MedianRTHigh','TotalTrialNum','LowTrialNum','HighTrialNum'});

% Get accuracy values
Accuracy = behaviorTable.AccuracyTotal;
Accuracy_Low = behaviorTable.AccuracyLow;
Accuracy_High = behaviorTable.AccuracyHigh;
Accuracy_EC=behaviorTable.AccuracyTotal(behaviorTable.Group==0);
Accuracy_AD= behaviorTable.AccuracyTotal(behaviorTable.Group==1);
Accuracy_ECLow=behaviorTable.AccuracyLow(behaviorTable.Group==0);
Accuracy_ADLow= behaviorTable.AccuracyLow(behaviorTable.Group==1);
Accuracy_ECHigh=behaviorTable.AccuracyHigh(behaviorTable.Group==0);
Accuracy_ADHigh= behaviorTable.AccuracyHigh(behaviorTable.Group==1);

%Perform Wilcoxon rank-sum tests comparing accuracy between groups
[p,h,stats]=ranksum(Accuracy_EC,Accuracy_AD)
[p,h,stats]=ranksum(Accuracy_ECLow,Accuracy_ADLow)
[p,h,stats]=ranksum(Accuracy_ECHigh,Accuracy_ADHigh)
%% Behavior analysis: Reaction time GLME
clear all
%Load participant list w/categorization
load('state_map.mat')
id_ec=map_ad((map_ad(:,2)==0),1);
id_ad=map_ad((map_ad(:,2)==1),1);

subject=[];
EC=[];
ECconflict=[];

AD=[];
ADconflict=[];

for ss=1:length(id_ec)
    filename=['P',num2str(id_ec(ss)),'.mat'];
    load(filename,'TrialDet');
    Trials_C=find(~isnan(TrialDet(:,12)) & TrialDet(:,12)<2 & TrialDet(:,14)==min(TrialDet(:,14)) & TrialDet(:,27)==1); %Looks for response, conflict type, and Accuracy
    Trials_I=find(~isnan(TrialDet(:,12)) & TrialDet(:,12)<2 & TrialDet(:,14)==max(TrialDet(:,14)) & TrialDet(:,27)==1);
    EC=[EC;TrialDet(Trials_C,12);TrialDet(Trials_I,12)];
    ECconflict=[ECconflict;zeros(length(Trials_C),1);ones(length(Trials_I),1)];
    sub=strings(length(Trials_C)+length(Trials_I),1);
    sub(:)=['P',num2str(id_ec(ss))];
    subject=[subject;sub];
end

for ss=1:length(id_ad)
    filename=['P',num2str(id_ad(ss)),'.mat'];
    load(filename,'TrialDet');
    Trials_C=find(~isnan(TrialDet(:,12)) & TrialDet(:,12)<2 & TrialDet(:,14)==min(TrialDet(:,14)) & TrialDet(:,27)==1);
    Trials_I=find(~isnan(TrialDet(:,12)) & TrialDet(:,12)<2 & TrialDet(:,14)==max(TrialDet(:,14)) & TrialDet(:,27)==1);
    AD=[AD;TrialDet(Trials_C,12);TrialDet(Trials_I,12)];
    ADconflict=[ADconflict;zeros(length(Trials_C),1);ones(length(Trials_I),1)];
    sub=strings(length(Trials_C)+length(Trials_I),1);
    sub(:)=['P',num2str(id_ad(ss))];
    subject=[subject;sub];
end

RT = [EC;AD];
LogRT=[log(EC);log(AD)];
Conflict=[ECconflict;ADconflict];
Group=[zeros(length(EC),1);ones(length(AD),1)];

save('RTdata.mat','EC','AD','ECconflict','ADconflict');

%GLME for RT
dataTable = table(LogRT,categorical(Conflict),categorical(Group),categorical(subject),...
            'VariableNames',{'LogRT','Conflict','Group','Subject'});
Mrt=fitglme(dataTable,'LogRT ~ Conflict + Group + (1|Subject)','Link','identity');
Mrt2 = fitglme(dataTable,'LogRT ~ Conflict*Group + (1|Subject)','Link','identity');

%Test contrasts of conflict and group
EC_Low_Conflict = [1 0 0 0]; 
EC_High_Conflict = [1 1 0 0];  
AD_Low_Conflict = [1 0 1 0];   
AD_High_Conflict = [1 1 1 1];
C1 = AD_High_Conflict-EC_High_Conflict;
C2 = AD_Low_Conflict-EC_Low_Conflict;
C3 = EC_High_Conflict-EC_Low_Conflict;
C4 = AD_High_Conflict-AD_Low_Conflict;
%Columns = [pval,Fval,DF1,DF2]
[posthocStats(1,1),posthocStats(1,2),posthocStats(1,3),posthocStats(1,4)]=coefTest(Mrt2,C1);
[posthocStats(2,1),posthocStats(2,2),posthocStats(2,3),posthocStats(2,4)]=coefTest(Mrt2,C2);
[posthocStats(3,1),posthocStats(3,2),posthocStats(3,3),posthocStats(3,4)]=coefTest(Mrt2,C3);
[posthocStats(4,1),posthocStats(4,2),posthocStats(4,3),posthocStats(4,4)]=coefTest(Mrt2,C4);
[~,critp,~,padj_posthoc]=fdr_bh(posthocStats(:,1));
coef_posthoc(1) = C1*fixedEffects(Mrt2);
coef_posthoc(2) = C2*fixedEffects(Mrt2);
coef_posthoc(3) = C3*fixedEffects(Mrt2);
coef_posthoc(4) = C4*fixedEffects(Mrt2);