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

%% Loading appropriate data, 1 participant at the time
%If the data we get can be converted into .set, then we can use following: 
gating_signal = EEG.data;

%EEG = pop_chanedit(EEG, 'convert')
%https://irenevigueguix.wordpress.com/2016/07/06/eeglab-tutorial-import-data/
%To import CSV data
%If the data we get is on csv format, then we can use following for loading
%data
%registration_data = importdata('data/data/0i4o2.csv');  %remember to
%change path
%timestamp = registration_data.data(:,1);
%all_waves= registration_data.data(:,2:22); %adjust to our channels.
%Possibly transpose if necessary

%% Filtering with bandpass filter 
EEG = pop_eegfilt(EEG,0.23,30,4); %0.23Hz-30Hz, 4th order
%pop_eegplot(EEG,1,1,1); Use this function to see the raw EEG data

EEG = pop_epoch(EEG, { 'square' }, [0 2], 'newname', 'Continuous EEG Data epochs','epochinfo' ,  'yes'); % function for epoching the data. Should
%only be run once, as it divides into epochs every time the script is run
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'setname','Continuous EEG Data epochs', 'overwrite', 'on'); % Overwrites EEG data that is now epoched by creating new set
%% Baseline correction
EEG = pop_rmbase( EEG, [0 161]); %removes channel baseline from the epoched data within limits

%% Artifact Rejection by the use of ICA
    
EEG = pop_runica(EEG, 'runica') %Runs ICA algorithm on epoched data. Run
%only once, else runs ICA on ICA etc. 
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG); %Creates new set based on the calculated ICA weights
eeglab redraw; %updates the EEGlab interface where the ICA weights are changed to available. Important for later, as otherwise, ERPs from the signals cannot be accessed
figure; 
topoplot([],EEG.chanlocs, 'style', 'blank', 'electrodes', 'labelpoint'); %Creates plot with the electrodes from which the 

pop_topoplot(EEG,1,[0:100:500],'ERP image',[2:3] ,0, 'electrodes', 'on'); %Plots scalp ERP activation maps
times = [0:100:500]; % Chooses the start and end ms of scalp ERP plots
pos = eeg_lat2point(times/1000, 1, EEG.srate, [EEG.xmin EEG.xmax]); %converts latencies in time untils to the latencies relative to data in epochs
mean_data = mean(EEG.data,3); %calculates the mean of the data used for plotting
maxlim = max(mean_data(:)); %sets limit for plots
minlim = min(mean_data(:)); %sets limit for plots
maplimits = [ -max(maxlim, -minlim) max(maxlim, -minlim)]; %places the limits onto the map
%scalp map series
figure
for k = 1:6
sbplot(2,3,k);% flexible version of subplot
topoplot( mean_data(:,k), EEG.chanlocs, 'maplimits', maplimits, 'electrodes','on', 'style', 'both'); %plots topographic map of EEG in 2D circular view
title([ num2str(times(k)) 'ms']);
end
figure, 
pop_timtopo(EEG,EEG.chanlocs) %Plots epoch mean for each channel on the axis and the scalp map for latencies 

figure;
plot(EEG.times, mean(mean(EEG.data([1 3 5 9 10],:,:),3),1)); %Averages epoched channels, here 1, 3, 5, 9, 10. 

%% Time analysis using Wavelets 
%% Quantifying ERP amplitudes and Latencies
%   - Peak amplitude , peak latency, mean amplitude
