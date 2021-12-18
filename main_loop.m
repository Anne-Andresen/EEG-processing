% For loop iterating through protocol 1+2 (step 5,6,7,8 in docs)

%for loop condition for when to end data acquistion 
    % when gating iteration = 120 & registration iteration = 400
%create time/trigger vector  
% while trigger(activated) & (trigger_gating >= 120 OR trigger_gating = 120 & trigger registration >=400)
    %map trigger value to the time vector
%end 

% This is a demo script for the use of g.HIamp in the g.NEEDaccess MATLAB
% API. It records data for 10 seconds from all analog channels available
% using the internal test signal.




%% Create connection + config
% create gtecDeviceInterface object
gds_interface = gtecDeviceInterface();

% define connection settings (loopback)
gds_interface.IPAddressHost = '127.0.0.1';
gds_interface.IPAddressLocal = '127.0.0.1';
gds_interface.LocalPort = 50224; %% check the local port and host port in settings
gds_interface.HostPort = 50223;

% get connected devices
connected_devices = gds_interface.GetConnectedDevices();

% create g.HIamp configuration object
ghiamp_config = gHIampDeviceConfiguration();
% set serial number in g.HIamp device configuration
ghiamp_config.Name = connected_devices(1,1).Name;

% set configuration to use functions in gds interface which require device
% connection
gds_interface.DeviceConfigurations = ghiamp_config;

% get available channels
available_channels = gds_interface.GetAvailableChannels()

% edit configuration to have a sampling rate of 256Hz, 4 scans,all
% available analog channels as well as the Counter (recorded on first
% analog channel).
% Acquire the internal test signal of g.HIamp
ghiamp_config.SamplingRate = 256;
ghiamp_config.NumberOfScans = 8;
ghiamp_config.CounterEnabled = true;
ghiamp_config.TriggerLinesEnabled = true;
%ghiamp_config.CounterEnabled = true;
ghiamp_config.HoldEnabled = false;
%ghiamp_siggen = gHIampInternalSignalGenerator();
%ghiamp_siggen.Enabled = true;
%ghiamp_siggen.Frequency = 10;
%ghiamp_config.InternalSignalGenerator = ghiamp_siggen;
for i=1:size(ghiamp_config.Channels,2)
    if (available_channels(1,i))
    	ghiamp_config.Channels(1,i).Available = true;
        if (i <= 80)
            ghiamp_config.Channels(1,i).Acquire = true;
        else
            ghiamp_config.Channels(1,i).Acquire = false;
        end
        % do not use filters
        ghiamp_config.Channels(1,i).BandpassFilterIndex = -1;
        ghiamp_config.Channels(1,i).NotchFilterIndex = -1;
        % do not use a bipolar channel
        ghiamp_config.Channels(1,i).ReferenceChannel = 0;
    end
end



%to showcase a signal
%scope_handle =dsp.TimeScope(2,256, 'BufferLength', 25600, 'TimeAxisLabels', 'Bottom', 'YLimits', [0 2560], 'TimeSpan',10,'LayoutDimensions', [2,1], 'ReduceUpdates', true, 'YLabel', 'Counter');
global T2Iterations TriggerPoint
TriggerPoint=[0,clock];
T2Iterations=0;

%% 
trigger_gating = 0;
trigger_registration = 0;
gating_iteration = 0; 
registration_iteration = 0;
samples_acquired = 0;
Vector=samples_acquired[];
trigger_activated = 'false'; % boolean operator 

while(gating_iteration <= 120 && registration_iteration <= 400)
    %time/trigger vector
    while trigger_activated == true & (gating_trigger <= 120 || registration_trigger <= 400) % conditional such that it can only be activated if trigger is true and holds either gating or registration condition iteration
        samples_acquired = samples_acquired/samples_acquired % overwrites the EEG data creating a peak at 1
    end
    
    
%% Gating
 for z = 1:2 
     if z = 1 
        u1 = size(samples_acquired)+1;
        %samples_acquired = 0;
data_received = single(zeros(size(u1), 80));
while (samples_acquired < size(u1)) % no of itterations of loop for stimuli
    try
        [scans_received, data] = gds_interface.GetData(8);
    catch ME
        disp(ME.message);
        break;
    end

    for gating_iteration <= 1:120; gating_iteration =120;  % Runs 60 times, thus playing 120 gating scenarios 
        gating_trigger= gating_trigger+1;
        trigger_activated ='true';
        soundsc(a,fs,slim);
        trigger_activated='false';
        pause(0.5); % waiting 500ms between clicks 
        trigger_activated='true';
        soundsc(a,fs,slim);
        trigger_activated = 'false';
     gating_iteration= gating_iteration+1;
     pause (10); % waiting 10s between paired clicks
    end  
% record data for 10 seconds, allocate memory for data acquired, data type
% is single
    data_received(((samples_acquired + 1) : (samples_acquired + scans_received)),:) = data;
  %  step(scope_handle, data(:,2), data(:,1));
    samples_acquired = samples_acquired + scans_received;  
end

     else
%% Registration 

u1 = size(samples_acquired)+1;
data_received = single(zeros(size(u1), 80));
while (samples_acquired < size(u1)) % no of itterations of loop for stimuli
    try
        [scans_received, data] = gds_interface.GetData(8);
    catch ME
        disp(ME.message);
        break;
    end

    for registration_iteration <= 1:4; registration_iteration =4;  % Runs 4 times, thus playing 400 registration scenarios 
        registration_trigger= registration_trigger+1;
        s=1;
        v=1;
        for s=1:100
        load('Random_matrix1.mat', 'Rand_Matrix1'); 
        trigger_activated ='true';
        for v=1:size(Random_matrix1)
       soundsc(Rand_Matrix1(v,:)); % play first matix sound 
        end
        trigger_activated='false';
         v= v+1;
        pause(2); % waiting 2 ssec between each sound 
        end
        s=s+1;
     registration_iteration= registration_iteration+1;
     pause (30); % waiting 30 sec between blocks 
    end
    end  
% record data for 10 seconds, allocate memory for data acquired, data type
% is single
    data_received(((samples_acquired + 1) : (samples_acquired + scans_received)),:) = data;
  %  step(scope_handle, data(:,2), data(:,1));
    samples_acquired = samples_acquired + scans_received;  
     end
 end
end

