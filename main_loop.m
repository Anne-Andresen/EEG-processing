%% For loop iterating through protocol 1+2 (step 5,6,7,8 in docs)

%for loop condition for when to end data acquistion 
    % when gating iteration = 120 & registration iteration = 400
%create time/trigger vector  
% while trigger(activated) & (trigger_gating >= 120 OR trigger_gating = 120 & trigger registration >=400)
    %map trigger value to the time vector
%end 


trigger_gating = 0;
trigger_registration = 0;
gating_iteration = 0; 
registration_iteration = 0;
samples_acquired = 0;

%trigger = (120*2)+400
trigger_activated = 'false'; % boolean operator 

while(gating_iteration <= 120 & registration_iteration <= 400)
    %time/trigger vector
    while trigger_activated == true & (gating_trigger <= 120 || registration_trigger <= 400) % conditional such that it can only be activated if trigger is true and holds either gating or registration condition iteration
        samples_acquired = samples_acquired/samples_acquired % overwrites the EEG data creating a peak at 1
    end 
    
    
 %% Gating
 
u1 = size(samples_acquired)+1;
samples_acquired = 0;
data_received = single(zeros(size(u1), 80));
while (samples_acquired < size(u1)) % no of itterations of loop for stimuli
    try
        [scans_received, data] = gds_interface.GetData(8);
    catch ME
        disp(ME.message);
        break;
    end

    for gating_iteration =< 1:120; gating_iteration =120;  % Runs 60 times, thus playing 120 gating scenarios 
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
end
