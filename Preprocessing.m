clear all; close all; clc; 
%Open eeglab functions - must be included in the path
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
%Load .set from the path
EEG = pop_loadset( 'Gating_Jens_64_channels.set', 'D:\7.semester\Projekt\Jens\Gat');
%Import channels
EEG.chanlocs=pop_chanedit(EEG.chanlocs, 'load',{'D:\7.semester\Projekt\10-20_64_channels.locs','filetype','autodetect'});
eeglab redraw;
%% Import events 
[EEG, eventnumbers] = pop_importevent(EEG, 'event', 'D:\7.semester\Projekt\Jens\Gat\event_gat_jens.txt', 'fields', {'latency','type','position'},'append', 'no', 'skipline', 1, 'timeunit',1)

%Store the dataset in the eeglab
[ALLEEG EEG CURRENTSET ] = eeg_store(ALLEEG, EEG);
%Update the interface of the EEGLab to follow whether all processing looks
%fine
eeglab redraw;
%% Downsampling
EEG = pop_select(EEG, 'nochannel', {'EEG', 'VEO', 'HEO'});

% Downsample to 250Hz
EEG = pop_resample( EEG, 250);

% Remove DC offset
EEG = pop_rmbase(EEG, 'mean',0);
pop_eegplot(EEG,1)

%% Filtering

% High pass filter at 0.23Hz

%EEG = pop_eegfiltnew(EEG, [], 0.3, [], true, [], 0);
EEG = pop_eegfiltnew(EEG, [], 0.235, [], true, [], 0);
%EEG = pop_eegfiltnew(EEG, [], 0.1, [], true, [], 0);

% Low pass filter at 45Hz
EEG = pop_eegfiltnew(EEG, [], 45, 66, 0, [], 0);



% Create a new dataset with the name Filtered Continuous EEG Data (Gating). Now CURRENTSET = 2
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'setname', 'Filtered Continuous EEG Data (Gating)');
% Re-refrence the new dataset such that filtered dataset is used for the
% upcoming procedures
EEG = pop_reref( EEG, [], 'refstate',0);
% Adding comments to the dataset (See comments by typing >> EEG.comments)
EEG.comments = pop_comments(EEG.comments,'','Dataset was bandpass filtered at 0.23-45 Hz and rereferenced.',1);

%Store the dataset in the eeglab
[ALLEEG EEG CURRENTSET ] = eeg_store(ALLEEG, EEG);
[ALLEEG EEG] = pop_saveset(EEG, 'filename', 'Gating_Jens_Continous_Filtered.set', 'filepath', 'D:\7.semester\Projekt\Jens\Gat' )
EEG = eeg_checkset(EEG)

eeglab redraw
% Keep original EEG as will need channels from here to interpolate
% after run clean_rawdata
originalEEG = EEG;

%% Clean data
rawEEG = EEG;
cleanLinedEEG = pop_cleanline(rawEEG) %Enter [50 100] to remove electrical line noise at 50 Hz in Europe and its harmonics at 100 Hz

% Run clean_rawdata which will get rid of flat or noisy channels and
% get rid of blinks, etc in continuous data
%EEG = clean_rawdata(EEG, 10, [0.25 0.75], 0.7, 4, 5, 0.5);

%artifact_rm_EEG = clean_artifacts(cleanLinedEEG,'ChannelCriterionMaxBadTime', 0.7, 4, 5, 0.5, [0.25 0.75])

artifact_rm_EEG = clean_artifacts(cleanLinedEEG,'WindowCriterion',0.25)


% Now interpolate the channels that were deleted in the previous step:
artifact_rm_EEG = pop_interp(artifact_rm_EEG, originalEEG.chanlocs, 'spherical');

% Reference electrodes to average:
% EEG = pop_reref( EEG, []);
% EEG = pop_saveset(EEG, 'filename', 'Gat_Jens_clean.set', 'filepath', 'C:\Users\ahaan\Desktop\7.semester\Projekt\Jens\Gat' )
% EEG = eeg_checkset(EEG)
%Update the interface of the EEGLab to follow whether all processing looks
%fine
eeglab redraw;  
%     % Finally, save the new file that has been cleaned:
% 	EEG = pop_saveset( EEG, 'filename',[filename, '_worked.set'],'filepath',[directory_name '\Worked']);

%% ICA 
EEG = pop_runica( EEG,'runica' ); % pops-up a data entry window
EEG=pop_selectcomps(EEG, [1:55])
EEG=pop_subcomp(EEG)
%pop_timtopo(EEG, [], [NaN], 'ERP data and scalp maps of Continuous EEG Data epochs');
SecondEEG=EEG;
%% Epoch division 

%% Epochs condition C_c
% Extract epochs time locked to the event - 'C_c' and 'T_c', from 100 ms before to 150 ms after those events.
EEG = pop_epoch( EEG, { 'C_c' }, [0 0.7], 'newname', 'Continuous EEG Data epochs Condition Click', 'epochinfo', 'yes');
%Create new dataset such that the filtered continous set can be utilized if
%necessary, at a later point
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'setname', 'Continuous EEG Data epochs C_c');

% Remove baseline based on the 100 ms before the event
% EEG = pop_rmbase( EEG, [-100 0]);
% Add a description of the epoch extraction to EEG.comments fro later use
EEG.comments = pop_comments(EEG.comments,'','Extracted ''C_c'' epochs [-0.100 0.150] ms, and corrected baseline based on 100 ms before event.',1); 
%Store the dataset in the eeglab
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET); 
%Update the interface of the EEGLab to follow whether all processing looks
%fine
[ALLEEG EEG] = pop_saveset(EEG, 'filename', 'Gat_Jens_Epochs_C_c.set', 'filepath', 'D:\7.semester\Projekt\Jens\Gat' )
EEG = eeg_checkset(EEG)
eeglab redraw;
%% Go into EEGlab, select the Continous EEG Data epochs before continueing with the next step
filename_for_T_c_epoching = 'Gating_Jens_Continous_Filtered.set'
filepath_for_T_c_epoching = 'C:\Users\ahaan\Desktop\7.semester\Projekt\Jens'
%EEG = pop_epoch( EEG, { 'C_c' }, [-0.100 0.150], 'newname', 'Continuous EEG Data epochs Condition Click', 'epochinfo', 'yes');
SecondEEG = pop_epoch( SecondEEG, { 'T_c' }, [0 0.7], 'newname', 'Continuous EEG Data epochs Test Click', 'epochinfo', 'yes');
%Create new dataset such that the filtered continous set can be utilized if
%necessary, at a later point
[ALLEEG SecondEEG CURRENTSET] = pop_newset(ALLEEG, SecondEEG, CURRENTSET, 'setname', 'Continuous EEG Data epochs T_c');

% Remove baseline based on the 100 ms before the event
% EEG = pop_rmbase( EEG, [-100 0]);
% Add a description of the epoch extraction to EEG.comments fro later use
EEG.comments = pop_comments(EEG.comments,'','Extracted ''T_c'' epochs [-0.100 0.150] ms, and corrected baseline based on 100 ms before event.',1); 
%Store the dataset in the eeglab
[ALLEEG SecondEEG] = eeg_store(ALLEEG, SecondEEG, CURRENTSET); 

SecondEEG = pop_saveset(SecondEEG, 'filename', 'Gat_Jens_Epochs_T_c.set', 'filepath', 'D:\7.semester\Projekt\Jens\Gat' )
SecondEEG = eeg_checkset(SecondEEG)
%Update the interface of the EEGLab to follow whether all processing looks
%fine
eeglab redraw;


%% Plot the ERPs based on the Epoching to make sure it looks fine before ICA for C_c

filename_C_c_epoched = 'Gat_Jens_Epochs_C_c.set'
filepath_C_c_epoched = 'D:\7.semester\Projekt\Jens\Gat'
EEG = pop_loadset('filename', filename_C_c_epoched ,'filepath',filepath_C_c_epoched);
[ALLEEG EEG CURRENTSET ] = eeg_store(ALLEEG, EEG);

%figure of scal activations in time interval -100ms to 100ms 
figure;
pop_timtopo(EEG)
%, [30 70], [NaN], 'ERP data and scalp maps of Continuous EEG Data epochs C_c');

eeglab redraw

%% %% Plot the ERPs based on the Epoching to make sure it looks fine before ICA for T_c
filename_T_c_epoched = 'Gat_Jens_Epochs_T_c.set'
filepath_T_c_epoched = 'D:\7.semester\Projekt\Jens\Gat'
SecondEEG = pop_loadset('filename', filename_T_c_epoched ,'filepath',filepath_T_c_epoched);
[ALLEEG SecondEEG CURRENTSET ] = eeg_store(ALLEEG, SecondEEG);

%figure of scal activations in time interval -100ms to 100ms 
figure;
pop_timtopo(SecondEEG, [0 70], [NaN], 'ERP data and scalp maps of Continuous EEG Data epochs T_c');

eeglab redraw


%% Averaging across epochs 
figure,
pop_plottopo( EEG, [], 'Average C_c', 0)
% AVG=mean(mean(SecondEEG(:,:,:),3))
% figure,
%  plot(EEG.times, mean(mean(EEG([1 3 5 7 9],:,:),3),1)); %Averages epoched channels, here 1, 3, 5, 9, 10. 
% %  
% Peaks=findpeaks(EEG.times,mean(mean(SecondEEG([1 3 5 7 9 11],:,:),3)));
% Amplitude=max(Peaks)
figure,
pop_plottopo(SecondEEG, [], 'Average T_c', 0)

pop_comperp( ALLEEG, 1,3,[],'addavg','on','addstd','off','addall','on','diffavg','off','diffstd','off','chans',[],'tplotopt',{'ydir' -1});
[EEGOUT,event_indices] = pop_selectevent( EEG, 'position', 1);% For C_c 
%'key2', value2, ... )
[EEGOUT1,event_indices1] = pop_selectevent( SecondEEG, 'position', 2);% For T_c 
pop_comperp(EEGOUT1,1)

%%

AVG_channel_20 = mean(EEG.data([20],:,:),3)
AVG_channel_21 = mean(EEG.data([21],:,:),3)
AVG_channel_22 = mean(EEG.data([22],:,:),3)
AVG_channel_23 = mean(EEG.data([23],:,:),3)
AVG_channel_24 = mean(EEG.data([24],:,:),3)
AVG_channel_29 = mean(EEG.data([29],:,:),3)
AVG_channel_30 = mean(EEG.data([30],:,:),3)
AVG_channel_31 = mean(EEG.data([31],:,:),3)
AVG_channel_32 = mean(EEG.data([32],:,:),3)
AVG_channel_33 = mean(EEG.data([33],:,:),3)
AVG_channel_38 = mean(EEG.data([38],:,:),3)
AVG_channel_39 = mean(EEG.data([39],:,:),3)
AVG_channel_40 = mean(EEG.data([40],:,:),3)
AVG_channel_41 = mean(EEG.data([41],:,:),3)
AVG_channel_42 = mean(EEG.data([42],:,:),3)

figure,
plot(EEG.times, AVG_channel_20)
hold on
plot(EEG.times, AVG_channel_21)
hold on
plot(EEG.times, AVG_channel_22)
hold on
plot(EEG.times, AVG_channel_23)
hold on
plot(EEG.times, AVG_channel_24)
hold on
plot(EEG.times, AVG_channel_29)
hold on
plot(EEG.times, AVG_channel_30)
hold on
plot(EEG.times, AVG_channel_31)
hold on
plot(EEG.times, AVG_channel_32)
hold on
plot(EEG.times, AVG_channel_33)
hold on
plot(EEG.times, AVG_channel_38)
hold on
plot(EEG.times, AVG_channel_39)
hold on
plot(EEG.times, AVG_channel_40)
hold on
plot(EEG.times, AVG_channel_41)
hold on
plot(EEG.times, AVG_channel_42)
legend
axis([0 80 0 2])
% 
t_int = [40, 70];                                     % Time Interval
idx = find((EEG.times >= t_int(1)) & (EEG.times <= t_int(2)));          % Indices Correspoinding To Time Interval
%[pks,locs] = findpeaks(AVG_channel_20(idx)) % Find Peaks In Time Imterval
%adjlocs = locs + idx(1)-1;   % Adjust ‘locs’ To Correct For Offset
Amplitude_ch_20=max(AVG_channel_20(idx))
%[pks1,locs1] = findpeaks(AVG_channel_21(idx))
Amplitude_ch_21=max(AVG_channel_21(idx))
%[pks2,locs2] = findpeaks(AVG_channel_22(idx)) 
Amplitude_ch_22=max(pks2)
%[pks3,locs3] = findpeaks(AVG_channel_23(idx)) 
Amplitude_ch_23=max(AVG_channel_23(idx))
%[pks4,locs4] = findpeaks(AVG_channel_24(idx)) 
Amplitude_ch_24=max(AVG_channel_24(idx))
%[pks5,locs5] = findpeaks(AVG_channel_30(idx)) 
Amplitude_ch_30=max(pks5)
%[pks6,locs6] = findpeaks(AVG_channel_31(idx)) 
Amplitude_ch_31=max(AVG_channel_31(idx))
%[pks7,locs7] = findpeaks(AVG_channel_32(idx)) 
Amplitude_ch_32=max(AVG_channel_32(idx))
%[pks8,locs8] = findpeaks(AVG_channel_33(idx)) 
Amplitude_ch_33=max(AVG_channel_33(idx))
%[pks9,locs9] = findpeaks(AVG_channel_38(idx)) 
Amplitude_ch_38=max(AVG_channel_38(idx))
%[pks10,locs10] = findpeaks(AVG_channel_39(idx)) 
Amplitude_ch_39=max(AVG_channel_39(idx))
%[pks11,locs11] = findpeaks(AVG_channel_40(idx)) 
Amplitude_ch_40=max(AVG_channel_40(idx))
%[pks12,locs12] = findpeaks(AVG_channel_41(idx))
Amplitude_ch_41=max(AVG_channel_41(idx))
%[pks13,locs13] = findpeaks(AVG_channel_42(idx)) 
Amplitude_ch_42=max(AVG_channel_42(idx))

%% Latencies 

 
pop_erpimage(EEG,1) 
 % [~,~,~,~,~,erp] = pop_erpimage(EEG,1, [31],[[]],'C_z','smooth_factor',1,{[31] EEG.chanlocs EEG.chaninfo})
%[], [], 1, 1, {}, [], '', 'erp','on','noshow','on');

%Only use the part of the ERP that represents the component you're analyzing, e.g.
samples = 100:200;
erp = erp(samples);

%Get the index in the ERP, as above.
[~,latency_index] = min(abs(cumsum(erp)-sum(erp)/2))

%Get the index in milliseconds from stimulus onset.
latency_ms = (EEG.xmin + (samples(1)+latency_index)/EEG.srate)*1000;


%% Extra 
% find_peaks_ch20=findpeaks(AVG_channel_20)
% Amplitude_ch_20=max(find_peaks_ch20)
% find_peaks_ch21=findpeaks(AVG_channel_21)
% Amplitude_ch_21=max(find_peaks_ch21)
% find_peaks_ch22=findpeaks(AVG_channel_22)
% Amplitude_ch_22=max(find_peaks_ch22)
% find_peaks_ch23=findpeaks(AVG_channel_23)
% Amplitude_ch_23=max(find_peaks_ch23)
% find_peaks_ch24=findpeaks(AVG_channel_24)
% Amplitude_ch_24=max(find_peaks_ch24)
% find_peaks_ch29=findpeaks(AVG_channel_29)
% Amplitude_ch_29=max(find_peaks_ch29)
% find_peaks_ch30=findpeaks(AVG_channel_30)
% Amplitude_ch_30=max(find_peaks_ch30)
% find_peaks_ch31=findpeaks(AVG_channel_31)
% Amplitude_ch_31=max(find_peaks_ch31)
% find_peaks_ch32=findpeaks(AVG_channel_32)
% Amplitude_ch_32=max(find_peaks_ch32)
% find_peaks_ch33=findpeaks(AVG_channel_33)
% Amplitude_ch_33=max(find_peaks_ch33)
% find_peaks_ch38=findpeaks(AVG_channel_38)
% Amplitude_ch_38=max(find_peaks_ch38)
% find_peaks_ch39=findpeaks(AVG_channel_39)
% Amplitude_ch_39=max(find_peaks_ch39)
% find_peaks_ch40=findpeaks(AVG_channel_40)
% Amplitude_ch_40=max(find_peaks_ch40)
% find_peaks_ch41=findpeaks(AVG_channel_41)
% Amplitude_ch_41=max(find_peaks_ch41)
% find_peaks_ch42=findpeaks(AVG_channel_42)
% Amplitude_ch_42=max(find_peaks_ch42)
% 
% 
% %% Amplitude 
% % 
% %  subjav1 = mean(EEG.data(31,:),3); %average across epochs ->  subave = chan x time points
% % subjav2 = mean(EEG.data(32,:),3); %average across epochs ->  subave = chan x time points
% % subjav2 = mean(EEG.data(38,:),3); %average across epochs ->  subave = chan x time points
% % subjav4 = mean(EEG.data(40,:),3); %average across epochs ->  subave = chan x time points
% % subjav5 = mean(EEG.data(39,:),3); %average across epochs ->  subave = chan x time points
% % subjav6 = mean(EEG.data(38,:),3); %average across epochs ->  subave = chan x time points
% % subjav7 = mean(EEG.data(41,:),3); %average across epochs ->  subave = chan x time points
% % subjav8 = mean(EEG.data(48,:),3); %average across epochs ->  subave = chan x time points
% 
% % AVG=mean(mean(EEG.data([20 21 22 23 30 31 32 38 39 40 41 48 49 50],:,:),3))
% % figure,
 plot(EEG.times, mean(mean(EEG.data([20 21 22 23 24 30 31 32 33 38 39 40 41 42],:,:),3),1)); %Averages epoched channels, here 1, 3, 5, 9, 10.
 axis([0 80 0 1])
 % % % 
% % mean_1 = mean(EEG.data,3);
% % mean_2 = mean(mean_1); 
% % mean_2 = mean_2(169:193);
% % 
% % % plot(EEG.times(169:223),mean_2(169:223))
% % % plot(EEG.times(169:193),mean_2)
% % 
% % % [pks, locs, widths, ~] = findpeaks(mean_2);
% % % peak = max(pks); 
% % 
% % % plot(EEG.times(169:193), mean_2); 
% % % hold on; 
% % % plot(locs, peak, 'or')
% % % plot(EEG.times(169:193),mean_2); hold on; plot(peak,'or')
% % 
% % [pks, locs, widths, ~] = findpeaks(mean_2);
% % peak = max(pks);
% % Number_peak = find(pks==peak); % finds placement in matrix 
% % Loc_peak = locs(1,Number_peak); % finds location of max peak 
% % Widths = widths(1,Number_peak); % finds width of max peak 
% % 
% % %Loc_peak=find(locs(peak));
% %  
% % % plot(EEG.times(169:193),mean_2); hold on; plot(locs, peak, 'or')
% % % plot(EEG.times(169:193),mean_2); hold on; plot(peak,'or')
% % 
% % lat=round(Loc_peak-(Widths/2),1); 
% % 
% % %Lat_half=find(round(peak/2)==mean_2);
% % Lat_half=find(round(peak/2),1)==mean_2;
% % 
% % %find en værdi som er peak/2 
% % %mellem lat og Loc_peak
% % %%
% % %Lat_half=find(peak/2); 
% % % 
% % while Lat_half >lat && Lat_half < Loc_peak
% % Lat_half=find(peak/2);        % finds location of half the latency time 
% % end 
% % 
% % 
% % Number_peak=find(pks==peak);
% % 
% % % 
% % % erp = mean( EEG.data, 3 );
% % % window = [0 70];
% % 
% % baseline = -50;
% % mean_amplitude = mean(  erp( :, (window(1)-baseline):(window(2)-baseline)), 2 )