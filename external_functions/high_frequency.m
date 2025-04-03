%% High frequency

eeglab;

% Define files and bad trials
name = 'TMSEEG_X35488_S27_B70_BiphasicPPC_110RMT_APPA.vhdr';

badTrials =  [32 33 58];
badChannels = [53];

TMSremoval_range =[-1.2 2.0];

EEG = pop_loadbv('/mnt/projects/iTEPs/AP-PA,PA-AP - MEP-latency/nobackup/EEG_Raw/X35488/', name);

    % 2) Epoch
    EEG = pop_epoch( EEG, {  'R  8'  }, [-0.5         0.5]);
    % 8) Remove bad trials 
    EEG = pop_select(EEG,'rmtrial',badTrials);
    % 3) Remove bad channels
    EEG = pop_select( EEG, 'rmchannel',badChannels);
    % 4) Remove TMS arifact
    EEG = pop_tesa_removedata( EEG, TMSremoval_range );
    % 5) Interpolate
    EEG = pop_tesa_interpdata( EEG, ['cubic'], [abs(TMSremoval_range(1)) TMSremoval_range(2)] );
    EEG2 = EEG;
    % 6) Bandpass filter - conservative
    %EEGfast = mybutter(EEG2, 500, 1000, 2, 'bandpass' );
    EEGslow = mybutter(EEG2, 0.1, 500, 2, 'bandpass');
    EEGall = mybutter(EEG2, 0.1, 1000, 2, 'bandpass');
    % Re-referencing
    %EEGfast = pop_reref( EEGfast, []);
    EEGslow = pop_reref( EEGslow, []);
    EEGall = pop_reref( EEGall, []);
    % 7) Baseline correction
    %EEGfast = pop_rmbase( EEGfast, [-110 -10] ,[]);
    EEGslow = pop_rmbase( EEGslow, [-110 -10] ,[]);
    EEGall = pop_rmbase( EEGall, [-110 -10] ,[]);
    % Store dataset
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);

   
%% Define 
X66650_M1_APPA_cubic_01_500 = EEGslow;   
X66650_M1_APPA_cubic_01_1000 = EEGall;   

save('X66650_M1_APPA_cubic_01_500', 'EEGslow');
save('X66650_M1_APPA_cubic_01_1000', 'EEGall');

%%
EEG_vis_channel = 16;
X66650_M1_APPA_16_rmrange_12_2_linear = mean(EEG.data(EEG_vis_channel, :, :),3);
%%time_vector = EEG.times;

%% Find peaks and troughs

start_index = find(time_vector >= 2, 1); % 2 ms in timevector
end_index = find(time_vector >= 6, 1); % 6 ms in timevector

[pks, locs] = findpeaks(EEG_average(start_index:end_index),time_vector(start_index:end_index),'MinPeakDistance',0.9, 'MinPeakProminence',0.1);
trough_signal = -EEG_average;
[trough, idx] = findpeaks(trough_signal(start_index:end_index),time_vector(start_index:end_index),'MinPeakDistance',0.9,'MinPeakProminence',0.1);
trough = -trough;

% Define latencies for further calculations
peak_latencies = locs;   % Latencies of peaks
trough_latencies = idx;  % Latencies of troughs

%% Plot
plot(time_vector,EEG_average, locs, pks, 'X', idx, trough, 'O')
hold on 

% Annotate peaks with latency
for i = 1:length(peak_latencies)
    text(peak_latencies(i), pks(i), sprintf('%.1f ms', peak_latencies(i)), ...
         'FontSize', 6, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
end

% Annotate troughs with  latency
for i = 1:length(trough_latencies)
    text(trough_latencies(i), trough(i), sprintf('%.1f ms', trough_latencies(i)), ...
         'FontSize', 6, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center');
end



xlim([-5 10])
xlabel('time [ms]')
ylim([-50 70])
ylabel('Amplitude [µV]')
title(['Peaks and throughs of i-TEPs in channel ' EEG.chanlocs(EEG_vis_channel).labels])
legend({'EEG Signal', 'Peaks', 'Troughs'}, 'Location', 'best');

%hold off;

%% Plot without latencies
plot(time_vector,EEG_average, locs, pks, 'X', idx, trough, 'O')
hold on 
xlim([-5 10])
xlabel('time [ms]')
ylim([-50 70])
ylabel('Amplitude [µV]')
title(['Peaks and throughs of i-TEPs in channel ' EEG.chanlocs(EEG_vis_channel).labels])
legend({'EEG Signal', 'Peaks', 'Troughs'}, 'Location', 'best');

%%

% Peak to trough

p1_t1 = pks(1)-trough(1);
t1_p2 = pks(2)-trough(1);
p2_t2 = pks(2)-trough(2);
t2_p3 = pks(3)-trough(2);