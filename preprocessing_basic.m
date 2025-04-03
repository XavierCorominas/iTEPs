%% Preprocessing pipeline for early iTEPs.

% DRCMR 2025 - XCT - BrainStim Methods group

% Preprocessing pipeline OPTION 1 - BASIC


% --> GENERAL INFO:
% The current preprocessing pieplpile is intended to provide a simple and fast eeg
% preprocessing focusing on the essential steps. The preprocessing might be
% sufficient to explore iTEPS in datasets non contaminated with large
% muscular (or other type of) artefacts.


% ---> STEPS:
% 0) Load data
% 1) Epoch
% 2) Remove bad trials 
% 3) Remove TMS arifact
% 4) Remove bad channels
% 5) Interpolate TMS removed data
% 6) Bandpass filter - conservative
% 7) Re-referencing
% 8) Baseline correction
% 9) Plot iTEP data. Segment fas oscillations and slow oscillations
% 10) Save preprocessed dataset

% --> EXTERNAL TOOLS TO INSTALL
% EEGlab and Fieltrip are NOT distributed with this code. Install them
% on your system and add the paths to matlab.

% EEGLAB install:  https://eeglab.org/tutorials/01_Install/Install.html
% -- Once EEglab installed, you will need to download the TESA plugin. From
% the eeglab interface go to : FILE-->MANAGE EEGLAB EXTENSIONS--> SEARCH TESA TOOLBOX AND INSTALL.

% FIELDTRIP install: https://www.fieldtriptoolbox.org/download/
% Fieltrip in only used for future analyses and not for the present
% preprocessing. IThere is no need to install if you are only runing this
% preprocessing code.

%% ADDPATHS 
clear all;
close all;
clc

%path to external functions.
addpath(genpath('/home/xavi/Documents/PROJECTS/iTEPS/eeg_analyses_tool/TMS_EEG_preprocessing/'))

%path eto eeglab
addpath('/home/xavi/Documents/PROJECTS/iTEPS/eeg_analyses_tool/TMS_EEG_preprocessing/external_functions/eeglab2025.0.0/')

%path to fieltrip
addpath('/home/xavi/Documents/PROJECTS/iTEPS/eeg_analyses_tool/TMS_EEG_preprocessing/external_functions/fieldtrip-20240110/')


%% 0. LOAD DATA

% Open EEGlab 
eeglab;

% Define file to load 
name_dataset = 'TMSEEG_X35488_S27_B70_BiphasicPPC_95RMT_APPA.vhdr';
path_dataset = '/home/xavi/Documents/PROJECTS/iTEPS/eeg_analyses_tool/data/X35488/';

%load file
EEG = pop_loadbv(path_dataset, name_dataset);

eeglab redraw

%% 1. EPOCH

%Epoch
EEG = pop_epoch( EEG, {  'R  8'  }, [-0.5         0.5]); % TMS pulses are stamped on the EEG as 'R  8' markers. Modify your marker if necessary for epoching.
eeglab redraw



% Figures
figure; pop_timtopo(EEG, [-5  10], [2.3         4.5         4.8], 'ERP data and scalp maps');
%
figure;
plot(EEG.times, mean(EEG.data, 3)', 'b'); 
xlim([-5 10]);
% ylim([-100 100])
set(gca, 'FontSize', 24, 'FontName', 'Arial');
xlabel('Time (ms)', 'FontSize', 24, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 24, 'FontName', 'Arial');
%title('EEG Signal', 'FontSize', 24, 'FontName', 'Arial');

%% 2. Remove bad trials
%% Trials rejection - 2.1 (identify bad trials) --> Select the bad trials
TMPREJ=[];
eegplot(EEG.data,'winlength',5,'command','pippo','srate',EEG.srate,'limits',[EEG.times(1) EEG.times(end)]);
'data have been displayed for first round of trials rejection'

 
%% Trials rejection - 2.2 (remove bad trials)
if ~isempty(TMPREJ)
    [trialrej elecrej]=eegplot2trial(TMPREJ,size(EEG.data,2),size(EEG.data,3));
else
    trialrej=[];
end

EEG.BadTr =find(trialrej==1);
EEG = pop_rejepoch( EEG, EEG.BadTr ,0); 

eeglab redraw



% Figures
figure; pop_timtopo(EEG, [-5  10], [2.3         4.5         4.8], 'ERP data and scalp maps');
%
figure;
plot(EEG.times, mean(EEG.data, 3)', 'b'); 
xlim([-5 10]);
% ylim([-100 100])
set(gca, 'FontSize', 24, 'FontName', 'Arial');
xlabel('Time (ms)', 'FontSize', 24, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 24, 'FontName', 'Arial');
%title('EEG Signal', 'FontSize', 24, 'FontName', 'Arial');

%% Trials rejection - 2.3 (de mean the data - substract the mean value of the data from each data point to remove DC offset and prepare data for filtering)

EEG = pop_rmbase( EEG, [-500 499.9] ,[]); % De meaning assuming a epoch from -500 499.9

%% 3. REMOVE TMS ARTIFACT

% TMS artefact lengths in ms
TMSremoval_range =[-1.2 2.0]; % in ms

% Remoe TMS pulse
EEG = pop_tesa_removedata( EEG, TMSremoval_range );
eeglab redraw


% Figures
figure; pop_timtopo(EEG, [-5  10], [2.3         4.5         4.8], 'ERP data and scalp maps');
%
figure;
plot(EEG.times, mean(EEG.data, 3)', 'b'); 
xlim([-5 10]);
% ylim([-100 100])
set(gca, 'FontSize', 24, 'FontName', 'Arial');
xlabel('Time (ms)', 'FontSize', 24, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 24, 'FontName', 'Arial');
%title('EEG Signal', 'FontSize', 24, 'FontName', 'Arial');

%% 4. REMOVE BAD CHANNELS


% To identify bad channels we run the following section applying a Wiener filter to automatically detect bad channels:
EEG_evo = mean(EEG.data,3);
[~, sigmas] = DDWiener(EEG_evo);
figure;
plot(sigmas,'*')

% labeling the very worst channels to not affect the ICA run
badC = find(sigmas > (median(sigmas) + 5*std(sigmas)));
goodC = setdiff(1:length(sigmas),badC);

EEG = pop_select( EEG, 'nochannel', badC);
badC
eeglab redraw

% Figures
figure; pop_timtopo(EEG, [-5  10], [2.3         4.5         4.8], 'ERP data and scalp maps');
%
figure;
plot(EEG.times, mean(EEG.data, 3)', 'b'); 
xlim([-5 10]);
% ylim([-100 100])
set(gca, 'FontSize', 24, 'FontName', 'Arial');
xlabel('Time (ms)', 'FontSize', 24, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 24, 'FontName', 'Arial');
%title('EEG Signal', 'FontSize', 24, 'FontName', 'Arial');

% Run wiener again to verify the bad channel has been removed:
EEG_evo = mean(EEG.data,3);
[~, sigmas] = DDWiener(EEG_evo);
figure;
plot(sigmas,'*')



%{
% Alternatively do it manually

% Indentify bad channel numbers
badChannels = [53];

% Remove channels from the data
EEG = pop_select( EEG, 'rmchannel',badChannels);
eeglab redraw

%}



%% 5. INTERPOLATE REMOVED WINDOW

%
EEG = pop_tesa_interpdata( EEG, ['cubic'], [abs(TMSremoval_range(1)) TMSremoval_range(2)] ); % Cubic interpolation. Alternatively do linear interpolation, or substitute tms segment with baseline data window
eeglab redraw


% Figures
figure; pop_timtopo(EEG, [-5  10], [2.3         4.5         4.8], 'ERP data and scalp maps');
%
figure;
plot(EEG.times, mean(EEG.data, 3)', 'b'); 
xlim([-5 10]);
% ylim([-100 100])
set(gca, 'FontSize', 24, 'FontName', 'Arial');
xlabel('Time (ms)', 'FontSize', 24, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 24, 'FontName', 'Arial');
%title('EEG Signal', 'FontSize', 24, 'FontName', 'Arial');

%% 6.BANDPASS FILTER - segment fast and total fqs. to not afect each other on the baseline correction

EEGslow = mybutter(EEG, 2, 500, 2, 'bandpass'); 
EEGall = mybutter(EEG, 2, 1000, 2, 'bandpass');

%% 7.REREFERENCING

EEGslow = pop_reref( EEGslow, []);
EEGall = pop_reref( EEGall, []);

%% 8.BASELINE CORRECTION

EEGslow = pop_rmbase( EEGslow, [-110 -10] ,[]);
EEGall = pop_rmbase( EEGall, [-110 -10] ,[]);
EEGfast = EEGall;
EEGfast.data = EEGall.data - EEGslow.data;

%% 9. PLOT FINAL DATA. SEGMENT FAST OSCILLATION around iTEP and plot

%Find the indices within this range
time_range = [-5, 10];
idx_range = find(EEG.times >= time_range(1) & EEG.times <= time_range(2));
new_time_vector = EEG.times(idx_range);

% Correct data to store only ITEPS segment
EEG_slow = EEGslow;
EEG_slow.data = EEG_slow.data(:,idx_range,:); 
EEG_slow.times = new_time_vector;


% Correct data to store only ITEPS segment
EEG_all = EEGall;
EEG_all.data = EEG_all.data(:,idx_range,:); 
EEG_all.times = new_time_vector;

% Get fast oscillations
EEG_fast = EEG_all;
EEG_fast.data = EEG_all.data - EEG_slow.data;

%
% Plot time series of channels under the stimulation location; mean of 3 channels on the left motor cortex: 17, 16 , 6) ;

figure;
hold on
set(gcf, 'Position', [100, 100, 717, 492]); % [x y width height]

% Define time vector (assuming 250 points and 15 ms total)
num_points = size(EEG_fast.data, 2); % Number of time points
time_axis = linspace(-5, 10, num_points); % Adjusted time axis from -5 to 10 ms

% Data processing
% Data fast
data = mean(EEG_fast.data(:,:,:), 3);
mean_channel_data(1,:) = data(6,:);
mean_channel_data(2,:) = data(16,:);
mean_channel_data(3,:) = data(17,:);
mean_channel_data_fast = mean(mean_channel_data, 1);

% Data slow
data = mean(EEG_slow.data(:,:,:), 3);
mean_channel_data(1,:) = data(6,:);
mean_channel_data(2,:) = data(16,:);
mean_channel_data(3,:) = data(17,:);
mean_channel_data_slow = mean(mean_channel_data, 1);

% Data all
data = mean(EEG_all.data(:,:,:), 3);
mean_channel_data(1,:) = data(6,:);
mean_channel_data(2,:) = data(16,:);
mean_channel_data(3,:) = data(17,:);
mean_channel_data_all = mean(mean_channel_data, 1);


% Plot vertical grey bar from 0 ms to 1.5 ms
x_patch = [0 2 2 0];  % X-coordinates (start and end of the bar)
y_patch = [-5 -5 5 5];     % Y-coordinates (full height of the plot)
patch(x_patch, y_patch, [0.7 0.7 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.8,'HandleVisibility', 'off'); % Grey with transparency

% Plot with correct time axis
plot(time_axis, mean_channel_data_fast, 'b', 'LineWidth', 4);
plot(time_axis, mean_channel_data_slow, 'r', 'LineWidth', 4);
plot(time_axis, mean_channel_data_all, 'k', 'LineWidth', 4);

% Formatting
xlim([-5 10]); % Set x-axis limits
xticks(-5:5:10); % Set x-axis ticks at -5, 0, 5, 10
ylim([-6 6]); % Set y-axis limits (if needed)
set(gca, 'FontSize', 20, 'FontName', 'Arial');

% Labels and Legend
xlabel('Time (ms)', 'FontSize', 20, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 20, 'FontName', 'Arial');
legend({'500-1000Hz', '2-500Hz', '2-1000Hz'}, 'FontSize', 20, 'Location', 'northeast');
title('Time series EEG Signal Comparison', 'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Arial');


% Add "TMS Pulse" text at (0, maxY) where maxY is a bit above the plot
maxY = 3.25; % Adjust this value if needed
text(-1.25, maxY, 'TMS', 'FontSize', 20, 'FontWeight', 'bold', 'Color', [0, 0, 0], 'HorizontalAlignment', 'center');

% Add red arrow pointing to 0 ms
annotation('arrow', ...
    [0.35 0.385], ...  % X-coordinates in normalized figure units (0.5 centers it at 0ms)
    [0.7 0.65], ... % Y-coordinates in normalized figure units (adjust as needed)
    'Color', [0, 0, 0], 'LineWidth', 4, 'HeadStyle', 'vback2'); % Customize arrow

hold off;


%%

%% SAVE DATASETS

EEG = pop_saveset( EEGall, 'filename',[name_dataset,'_all_oscill_cleaned_pipeline_1.set'],'filepath',[path_dataset]);
EEG = pop_saveset( EEGslow, 'filename',[name_dataset,'_slow_oscill_cleaned_pipeline_1.set'],'filepath',[path_dataset]);
EEG = pop_saveset( EEGfast, 'filename',[name_dataset,'_fast_oscill_cleaned_pipeline_1.set'],'filepath',[path_dataset]);


%% END

