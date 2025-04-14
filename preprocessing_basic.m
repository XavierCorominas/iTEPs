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
name_dataset = 'AC_STATE_6.1.90CO25.vhdr';
path_dataset = '/mnt/projects/iTEPs/iTEPs_State/Pilots/RAW/EEG_Raw/AC/';

%load file
EEG = pop_loadbv(path_dataset, name_dataset);

eeglab redraw

%% 1. EPOCH

%Epoch
EEG = pop_epoch( EEG, {  'R  8'  }, [-0.5         0.5]); % TMS pulses are stamped on the EEG as 'R  8' markers. Modify your marker if necessary for epoching.
eeglab redraw



% Figures
figure; pop_timtopo(EEG, [-5  10], [2.3         3         4.8], 'ERP data and scalp maps');
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
figure; pop_timtopo(EEG, [-5  10], [2.3         3         4.8], 'ERP data and scalp maps');
%
figure;
plot(EEG.times, mean(EEG.data, 3)', 'b'); 
xlim([-5 10]);
% ylim([-100 100])
set(gca, 'FontSize', 24, 'FontName', 'Arial');
xlabel('Time (ms)', 'FontSize', 24, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 24, 'FontName', 'Arial');
%title('EEG Signal', 'FontSize', 24, 'FontName', 'Arial');


%% Plot raw data

% Data
EEGr = EEG;

% Demean
EEGr = pop_rmbase( EEGr, [-500 499.9] ,[]); % De meaning assuming a epoch from -500 499.9
EEGr = pop_rmbase( EEGr, [-110 -10] ,[]);

% Mean of Channels of itnerest
data_to_plot (1,:,:) = EEGr.data(6,:,:);
data_to_plot (2,:,:) = EEGr.data(16,:,:);
data_to_plot (3,:,:) = EEGr.data(17,:,:);
mean_roi_raw = mean(data_to_plot,1);

% Figure
figure;
set(gcf, 'Position', [100, 100, 2000, 600]); % Manually adjust figure size [x y width height]

subplot(1,2,1)

% Plot individual sensors (mean over trials), in grey
h1 = plot(EEGr.times, mean(EEGr.data, 3)', 'Color', [0.5 0.5 0.5]); 
hold on;

% Plot mean of selected ROI channels, in red
h2 = plot(EEGr.times, mean(mean_roi_raw, 3)', 'r'); 

% Axes settings
xlim([-5 10]);
ylim([-80 80]);
set(gca, 'FontSize', 24, 'FontName', 'Arial');
xlabel('Time (ms)', 'FontSize', 24, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 24, 'FontName', 'Arial');
title('Raw EEG Signal 2-50000Hz', 'FontSize', 24, 'FontName', 'Arial');

% Legend
% Create dummy lines for legend
dummy1 = plot(NaN, NaN, 'Color', [0.5 0.5 0.5], 'LineWidth', 2.5);
dummy2 = plot(NaN, NaN, 'r', 'LineWidth', 2.5);

% Add legend using dummies
legend([dummy1, dummy2], {'All sensors', 'Mean left M1 sensors'}, 'FontSize', 18);
set(gca,'FontName','Arial','fontsize',20,'FontWeight','bold','LineWidth',1.5)

hold off;



% Plot 2: raw data 2-2000Hz
% Filter raw data
EEGrf = pop_tesa_filtbutter( EEG, 0.2, 2000, 4, 'bandpass' );

% Demean
EEGrf = pop_rmbase( EEGrf, [-500 499.9] ,[]); % De meaning assuming a epoch from -500 499.9
EEGrf = pop_rmbase( EEGrf, [-110 -10] ,[]);


% Mean of Channels of itnerest
data_to_plot2 (1,:,:) = EEGrf.data(6,:,:);
data_to_plot2 (2,:,:) = EEGrf.data(16,:,:);
data_to_plot2 (3,:,:) = EEGrf.data(17,:,:);
mean_roi_raw_fil = mean(data_to_plot2,1);

% Figure
subplot(1,2,2)
h1 = plot(EEGrf.times, mean(EEGrf.data, 3)', 'Color', [0.5 0.5 0.5]); 
xlim([-5 10]);
ylim([-80 80])
set(gca, 'FontSize', 24, 'FontName', 'Arial');
xlabel('Time (ms)', 'FontSize', 24, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 24, 'FontName', 'Arial');
title('Raw EEG Signal 2-2000Hz', 'FontSize', 24, 'FontName', 'Arial');
hold on
h2 = plot(EEGr.times, mean(mean_roi_raw_fil, 3)', 'r'); 

% Legend
% Create dummy lines for legend
dummy1 = plot(NaN, NaN, 'Color', [0.5 0.5 0.5], 'LineWidth', 2.5);
dummy2 = plot(NaN, NaN, 'r', 'LineWidth', 2.5);

% Add legend using dummies
legend([dummy1, dummy2], {'All sensors', 'Mean left M1 sensors'}, 'FontSize', 18);
set(gca,'FontName','Arial','fontsize',20,'FontWeight','bold','LineWidth',1.5)

hold off;

%{
% In case you want topographies: modify data file to plot accordingly:
% Define the time points you want to plot (in ms)
time_points = [2.3, 3.5, 4.8]; % Example time points

% Number of time points
num_plots = length(time_points);

% Create a figure
figure;
set(gcf, 'Position', [100, 100, 2000, 900]); % Adjust figure size to fit both plots

% Create a subplot for topoplots (this will take up the top half of the figure)
subplot(2, 1, 1); % 2 rows, 1 column (top row)
% Loop over each time point for topography plots
for i = 1:num_plots
    time_point = time_points(i); % Get the current time point in ms

    % Find the nearest time index
    [~, time_idx] = min(abs(EEG.times - time_point)); % Find the closest index
    
    % Create a subplot for the current time point (topographies on the left)
    % We're using `subplot(1, num_plots, i)` for all topoplots together in one big subplot.
    subplot(1, num_plots, i); % In the 1 row, 'num_plots' columns (placing topoplots in one row)
    
    % Check if the selected index is valid and there is data at that time
    if ~isempty(time_idx) && ~all(isnan(EEG.data(:, time_idx))) 
        % Plot topography
        topoplot(EEG.data(:, time_idx), EEG.chanlocs, 'electrodes', 'on');
        title([num2str(EEG.times(time_idx)) ' ms'],'FontSize', 15); % Title with the actual time point
        % Add colorbar with custom settings
        colorbar; % Add colorbar
        caxis([-20 20]); % Set color limits (adjust based on your data range)
    else
        % If no valid data, plot a message
        topoplot([], EEG.chanlocs, 'electrodes', 'off'); % Empty topography
        title('No Data'); % Indicate no data for this time point
    end
end



%}
%% Trials rejection - 2.3 (de mean the data - substract the mean value of the data from each data point to remove DC offset and prepare data for filtering)

EEG = pop_rmbase( EEG, [-500 499.9] ,[]); % De meaning assuming a epoch from -500 499.9


% Figures
figure; pop_timtopo(EEG, [-5  10], [2.3         3         4.8], 'ERP data and scalp maps');
%
figure;
plot(EEG.times, mean(EEG.data, 3)', 'b'); 
xlim([-5 10]);
% ylim([-100 100])
set(gca, 'FontSize', 24, 'FontName', 'Arial');
xlabel('Time (ms)', 'FontSize', 24, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 24, 'FontName', 'Arial');
%title('EEG Signal', 'FontSize', 24, 'FontName', 'Arial');

%% 3. REMOVE TMS ARTIFACT

% TMS artefact lengths in ms
TMSremoval_range =[-1.2 2.0]; % in ms

% Remoe TMS pulse
EEG = pop_tesa_removedata( EEG, TMSremoval_range );
eeglab redraw


% Figures
figure; pop_timtopo(EEG, [-5  10], [2.3         3         4.8], 'ERP data and scalp maps');
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
figure; pop_timtopo(EEG, [-5  10], [2.3         3         4.8], 'ERP data and scalp maps');
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
figure; pop_timtopo(EEG, [-5  10], [2.3         3         4.8], 'ERP data and scalp maps');
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
set(gcf, 'Position', [100, 100, 900, 492]); % Manually adjust figure size [x y width height]

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

% Plot with correct time axis
plot(time_axis, mean_channel_data_fast, 'b', 'LineWidth', 4);
plot(time_axis, mean_channel_data_slow, 'r', 'LineWidth', 4);
plot(time_axis, mean_channel_data_all, 'k', 'LineWidth', 4);

% Formatting
xlim([-5 10]); % Set x-axis limits
xticks(-5:5:10); % Set x-axis ticks at -5, 0, 5, 10
set(gca, 'FontSize', 20, 'FontName', 'Arial');

% Get current y-axis limits and adjust them (add +1 to the upper limit)
y_limits = ylim; 
ylim([y_limits(1), y_limits(2) + 1]);  % Increase the upper limit by 1

% Plot vertical grey bar from 0 ms to 1.5 ms
x_patch = [0 2 2 0];  % X-coordinates (start and end of the bar)
y_patch = [y_limits(1) y_limits(1) y_limits(2)+1 y_limits(2)+1]; % Y-coordinates based on updated y-limits
patch(x_patch, y_patch, [0.7 0.7 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.8, 'HandleVisibility', 'off'); % Grey with transparency

%
% Plot a black dotted line at y = 0
plot([-5 10], [0 0], 'k--', 'LineWidth', 2); % Dotted black line at y = 0

% Labels and Legend
xlabel('Time (ms)', 'FontSize', 20, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 20, 'FontName', 'Arial');

% Adjust the legend location to be outside the plot area
legend({'500-1000Hz', '2-500Hz', '2-1000Hz'}, 'FontSize', 20, 'Location', 'northeastoutside', 'Orientation', 'vertical');

% Title
title('Time series EEG Signal Comparison', 'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Arial');

hold off


%{
% Add "TMS Pulse" text at (0, maxY) where maxY is a bit above the plot
maxY = 3.25; % Adjust this value if needed
text(-1.25, maxY, 'TMS', 'FontSize', 20, 'FontWeight', 'bold', 'Color', [0, 0, 0], 'HorizontalAlignment', 'center');

% Add red arrow pointing to 0 ms
annotation('arrow', ...
    [0.35 0.385], ...  % X-coordinates in normalized figure units (0.5 centers it at 0ms)
    [0.7 0.65], ... % Y-coordinates in normalized figure units (adjust as needed)
    'Color', [0, 0, 0], 'LineWidth', 4, 'HeadStyle', 'vback2'); % Customize arrow

hold off;
%}

%%

%% SAVE DATASETS

EEG = pop_saveset( EEGall, 'filename',[name_dataset,'_all_oscill_cleaned_pipeline_1.set'],'filepath',[path_dataset]);
EEG = pop_saveset( EEGslow, 'filename',[name_dataset,'_slow_oscill_cleaned_pipeline_1.set'],'filepath',[path_dataset]);
EEG = pop_saveset( EEGfast, 'filename',[name_dataset,'_fast_oscill_cleaned_pipeline_1.set'],'filepath',[path_dataset]);


%% END

path_dataset = '/home/xavi/Documents/PROJECTS/iTEPS/eeg_analyses_tool/TMS_EEG_preprocessing/'

EEG = pop_saveset( EEGfast, 'filename',[name_dataset,'_fast_oscill_cleaned_pipeline_1.set'],'filepath',[path_dataset]);


