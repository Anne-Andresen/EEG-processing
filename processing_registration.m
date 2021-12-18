%% Steps in data processing according to Introduction to the Event Related potential technique
%by Steven J luck (2005) - find it on aub

%High-pass to continous EEG (our case butterworth)
%Epoch extraction
%Baseline correction
%   -offset
%Artifact Rejection
%   - eyeblinks
%   - movement etc
%Averaging - across signals see page 259 
%   - Signal to Noise Ratio
%   - Variability in latency
%   - Variability in amplitude
% Time analysis using Wavelets 
%Quantifying ERP amplitudes and Latencies
%   - Peak amplitude , peak latency, mean amplitude
%% Filtering continous EEG registration
%If the data we get can be converted into .set, then we can use following: 
registration_signal = EEG.data;

%https://irenevigueguix.wordpress.com/2016/07/06/eeglab-tutorial-import-data/
%To import CSV data
%If the data we get is on csv format, then we can use following for loading
%data
%registration_data = importdata('data/data/0i4o2.csv');  %remember to
%change path
%timestamp = registration_data.data(:,1);
%all_waves= registration_data.data(:,2:22); %adjust to our channels.
%Possibly transpose if necessary


x_reg=registration_signal;


%% Epoching registration
% loop for each channel and input single channel into buffer

%% Filtering the signal 
EEG = pop_eegfilt(EEG,10,50,4);% When performing the real change 50 to 200. 

%% Ephochs of the signal 
EEG = pop_epoch(EEG, { 'square' }, [-1 8], 'newname', 'Continuous EEG Data epochs','epochinfo' ,  'yes');
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'setname','Continuous EEG Data epochs', 'overwrite', 'on');
%% Baseline correction in accordance with the 
EEG = pop_rmbase( EEG, [0 100]);

%% ICA 
EEG = pop_runica(EEG, 'runica')
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
eeglab redraw;

pop_timtopo(EEG,EEG.chanlocs)


%% Averaging 
AVG=mean(mean(EEG.data([1 3 5 9 10],:,:),3))
figure;
plot(EEG.times, mean(mean(EEG.data([1 3 5 9 10],:,:),3),1))





%% Time analysis using Wavelets 
%% Quantifying ERP amplitudes and Latencies
%   - Peak amplitude , peak latency, mean amplitude