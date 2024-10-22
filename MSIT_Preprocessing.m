%% Plot Trials/channels and exclude artifacts
clear all
addpath C:\Users\shekaraa\Documents\MATLAB\fieldtrip-20230613
%Load raw data frpm participant
load('P1.mat','ft_data3')

%Remove trials with NaNs 
cfg=[];
cfg.trials = ~cell2mat(cellfun(@anynan, ft_data3.trial, 'UniformOutput', false)); %finds trials that have nans, selects all others to be included
ft_data3 = ft_selectdata(cfg,ft_data3);

%Check if high pass filtering is needed
cfg=[];
ft_databrowser(cfg,ft_data3);

%Resample data if needed and bandstop line noise + 0.5Hz highpass
cfg = [];
cfg.hpfilter = 'yes';
cfg.hpfreq = 0.5;
cfg.hpfiltord = 5;
cfg.bsfilter = 'yes';
cfg.bsfreq = [55 65];
cfg.bsfiltord = 4;
ft_data3_rs = ft_preprocessing(cfg,ft_data3);

cfg=[];
ft_databrowser(cfg,ft_data3_rs);

if ft_data3_rs.fsample ~=512
cfg=[];
cfg.resamplefs = 512;
cfg.detrend = 'no';
ft_data3_rs = ft_resampledata(cfg,ft_data3_rs);
end
save('P29.mat','ft_data3_rs','-append');

%Plot Power spectra for all channels
cfg = [];
cfg.output = 'pow';
cfg.keeptrials = 'yes';
cfg.method = 'wavelet';
cfg.pad='nextpow2';
cfg.foi = [2:50 55:5:200];
cfg.toi = [-1:0.01:2];
ft_freq_rs = ft_freqanalysis(cfg,ft_data3_rs);

freq = ft_freq_rs.freq;
chanPSD = squeeze(nanmean(ft_freq_rs.powspctrm(:,:,:,:)));
plot(freq,nanmean(chanPSD,3)');xlim([0,200]); set(gca,'Yscale','log');legend(ft_freq_rs.label)
title('Power Spectra of all channels')

% Select threshold for artifact detection
cfg = [];
cfg.artfctdef.zvalue.channel = 'all';
cfg.artfctdef.zvalue.interactive = 'yes';
cfg.artfctdef.zvalue.cutoff = 7;
cfg.artfctdef.zvalue.artpadding = 0.1;
cfg.artfctdef.zvalue.rectify = 'yes';
ft_artifact_zvalue(cfg, ft_data3_rs);

%Now run artifact detection
load('P29.mat', 'ParcellationValues');
include_chans = find(~isnan(ParcellationValues(:,8))); %Select chans that aren't white matter
data={};
for ch = 1:length(ft_data3_rs.label)
    %Select channel data
    cfg = [];
    cfg.channel = ft_data3_rs.label{ch};
    data{ch}= ft_selectdata(cfg,ft_data3_rs);
    disp(ch)
    if any(include_chans == ch) 
        %Define artifact conditions
        cfg = [];
        cfg.artfctdef.zvalue.channel = 'all';
        %cfg.artfctdef.zvalue.interactive = 'yes';
        cfg.artfctdef.zvalue.cutoff = 7; %input zscore threshold here
        cfg.artfctdef.zvalue.artpadding = 0.1;
        cfg.artfctdef.zvalue.rectify = 'yes';
        [~, artifact.zvalue] = ft_artifact_zvalue(cfg, data{ch});
        %Scrub through data and highlight remaining artifacts
        cfg = [];
        cfg.continuous = 'no';
        cfg.ylim = [-1650 1650]; %adjust as needed
        cfg.artfctdef.zvalue.artifact = artifact.zvalue;
        cfg.artfctdef.visual.artifact = [];
        f = figure();
        f.Position = [900 200 900 800];
        cfg.figure = f;
        artf = ft_databrowser(cfg,data{ch});
        % replace artifacts with NaNs
        cfg = [];
        cfg.artfctdef.zvalue.artifact = artifact.zvalue;
        cfg.artfctdef.visual.artifact = artf.artfctdef.visual.artifact;
        cfg.artfctdef.reject = 'nan';
        data{ch} = ft_rejectartifact(cfg,data{ch});
    end
end

% Rejoin all the channels
cfg = [];
data_post = ft_appenddata([],data{:});

% Browse data again to make sure channels look good
cfg = [];
cfg.channel = ft_data3_rs.label(include_chans);
ft_databrowser(cfg,data_post);

%Save if everything looks good
ft_data3_clean = data_post;
save('P1.mat','ft_data3_clean','-append');

%:::::::::::::::::::::Optional steps below:::::::::::::::::::::

%Remove any artifacts that were missed
redo_chan = []; %idx of chans you need to redo
for ch = redo_chan
    cfg = [];
    cfg.channel = ft_data3_rs.label{ch};
    data{ch}= ft_selectdata(cfg,ft_data3_rs);
    %Define artifact conditions
    cfg = [];
    cfg.artfctdef.zvalue.channel = 'all';
    %cfg.artfctdef.zvalue.interactive = 'yes';
    cfg.artfctdef.zvalue.cutoff = 13;
    cfg.artfctdef.zvalue.artpadding = 0.1;
    cfg.artfctdef.zvalue.rectify = 'yes';
    [~, artifact.zvalue] = ft_artifact_zvalue(cfg, data{ch});
    %Scrub through data and highlight remaining artifacts
    cfg = [];
    cfg.continuous = 'no';
    cfg.ylim = [-1650 1650]; %adjust as needed
    cfg.artfctdef.zvalue.artifact = artifact.zvalue;
    cfg.artfctdef.visual.artifact = [];
    f = figure();
    f.Position = [900 200 900 800];
    cfg.figure = f;
    artf = ft_databrowser(cfg,data{ch});
    % replace artifacts with NaNs
    cfg = [];
    cfg.artfctdef.zvalue.artifact = artifact.zvalue;
    cfg.artfctdef.visual.artifact = artf.artfctdef.visual.artifact;
    cfg.artfctdef.reject = 'nan';
    data{ch} = ft_rejectartifact(cfg,data{ch});
end
cfg = [];
data_post = ft_appenddata([],data{:});

% Insert re-cleaned data back into ft_data3_clean
for ch = 1:length(redo_chan)
    for ii = 1:length(ft_data3_clean.trial)
        ft_data3_clean.trial{ii}(redo_chan(ch),:) = data_post.trial{ii}(ch,:);
    end
    ft_data3_clean.cfg.previous{redo_chan(ch)} = data_post.cfg.previous{ch};
end

%Save data if everything looks good
save('P1.mat','ft_data3_clean','-append');


% Exclude any bad channels found during artifact cleanup
preproc_ch_excluded={}; %Enter bipolar channel labels to exlcude
cfg = [];
cfg.channel = ['all', preproc_ch_excluded]; %input channels you want to exclude with a "-" before each label
ft_data3_clean = ft_selectdata(cfg,data_post);

%Save bad channels to mat file
cfg = [];
cfg.channel = {}; %input channels you excluded above
ft_data3_bad = ft_selectdata(cfg,data_post);
save ('P1.mat', 'ft_data3_bad','-append');

%Exclude bad channels from ParcellationValues
badidx = find(ismember(ft_data3_rs.label,cfg.channel));
ParcellationValues(badidx,:)=[];
save('P1.mat','ParcellationValues','-append')

%Save data
save('P1.mat','ft_data3_clean','-append');

%Re-run preprocessing and remove outliers
data={};
for ch = 1:length(ft_data3_rs.label)
    %Select channel data
    cfg = [];
    cfg.channel = ft_data3_rs.label{ch};
    data{ch}= ft_selectdata(cfg,ft_data3_rs);
    disp(ch)
    if any(include_chans == ch) 
        % replace artifacts with NaNs
        cfg = [];
        cfg.artfctdef.zvalue.artifact = ceil(data_post.cfg.previous{ch}.artfctdef.zvalue.artifact./4);
        cfg.artfctdef.visual.artifact = ceil(data_post.cfg.previous{ch}.artfctdef.visual.artifact./4);
        cfg.artfctdef.reject = 'nan';
        data{ch} = ft_rejectartifact(cfg,data{ch});
    end
end

% Rejoin all the channels and double check
cfg = [];
data_post_new = ft_appenddata([],data{:});
cfg = [];
cfg.channel = ft_data3_rs.label(include_chans);
ft_databrowser(cfg,data_post_new);
%% Spectral power after cleanup
clear all
%Load participant list w/categorization
load('state_map.mat','map_ad')

%Get TF data before and after artifact cleaning
for ii = 1:length(map_ad)
    filename=['P',num2str(map_ad(ii,1)),'.mat'];
    load(filename,'ft_data3', 'ft_data3_clean')
    cfg = [];
    cfg.output = 'pow';
    cfg.keeptrials = 'yes';
    cfg.method = 'wavelet';
    cfg.pad='nextpow2';
    cfg.foi = [2:50 55:5:200];
    cfg.toi = [-1:(1/512*12):2];
    ft_freq_clean = ft_freqanalysis(cfg,ft_data3_clean);
    ft_freq = ft_freqanalysis(cfg,ft_data3);

save(filenames(kk).name,'ft_freq','ft_freq_clean','-append')
end

%% Add Left/Right hemisphere Information to Data
clear all
%Load participant list w/categorization
load('state_map.mat','map_ad')

for ii = 1:length(map_ad)
    filename=['P',num2str(map_ad(ii,1)),'.mat'];
    load(filename, 'ft_freq_clean','ParcellationValues');
    chan = ft_freq_clean.label;
    Parcellation_Sided={};
    for jj = 1:length(chan)
        if startsWith(chan{jj},'L') == 1 || contains(chan{jj},'''') == 1
            Parcellation_Sided{jj} = ['L' num2str(ParcellationValues(jj,8))];
        else
            Parcellation_Sided{jj} = ['R' num2str(ParcellationValues(jj,8))];
        end
    end
    save(filename,'Parcellation_Sided','-append');
end
