%% Preprocessing pipeline for early iTEPs conserving the TMS pulse

% DRCMR 2025 - YS - BrainStim Methods group

% --> GENERAL INFO:
% The current preprocessing pieplpile is intended to provide a simple and fast eeg
% preprocessing focusing on the essential steps. The preprocessing might be
% sufficient to explore iTEPS in datasets non contaminated with large
% muscular (or other type of) artefacts.

% ---> STEPS:
% 0) Load data
% 1) Epoch
% 2) Demean
% 3) Remove bad trials or channels
% 4) Baseline correction
% 5) Bandpass filtering
% 6) Save data 


% --> EXTERNAL TOOLS TO INSTALL
% EEGlab and Fieltrip are NOT distributed with this code. Install them
% on your system and add the paths to matlab.

% EEGLAB install:  https://eeglab.org/tutorials/01_Install/Install.html
% -- Once EEglab installed, you will need to download the TESA plugin. From
% the eeglab interface go to : FILE-->MANAGE EEGLAB EXTENSIONS--> SEARCH TESA TOOLBOX AND INSTALL.


%% ADDPATHS 
clear all;
close all;
clc


%path eto eeglab
addpath('.../eeglab2025.0.0/')



%% 0. LOAD DATA

% Open EEGlab 
eeglab;

% Define file to load 
name_dataset = 'AC_STATE_6.1.90CO25.vhdr';
save_name = 'AC_STATE_6.1.90CO25';
path_dataset = 'path_to_data_folder/';

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

%% 2. Demean

EEG = pop_rmbase( EEG, [-500 499.9] ,[]); % De meaning assuming a epoch from -500 499.9
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

%% 3. Remove bad trials or channels
%% trials rejection - 3.1 (identify bad trials) --> Select the bad trials
TMPREJ=[];
eegplot(EEG.data,'winlength',5,'command','pippo','srate',EEG.srate,'limits',[EEG.times(1) EEG.times(end)]);
'data have been displayed for first round of trials rejection'

 
%% trials rejection - 3.2 (remove bad trials)
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

%% 3.3 manual rejection and interpolation of channels

% Indentify bad channel numbers
badChannels = [53]; % example: channel 53

% Save original channel locations before removal
originalChanLocs = EEG.chanlocs;

% Remove bad channels
EEG = pop_select(EEG, 'nochannel', badChannels);
eeglab redraw

% Interpolate the removed channels
EEG = pop_interp(EEG, originalChanLocs, 'spherical'); % or use 'invdist' or 'nearest'
eeglab redraw

%% 4.Baseline correction

EEG = pop_rmbase( EEG, [-210 -10] ,[]);
eeglab redraw

% Plot data with scalp maps
figure; pop_timtopo(EEG, [-500 499], [-200 0 100 200 300 400], 'After Re-Baseline Correction');

% Figure for paper
figure;
plot(EEG.times, mean(EEG.data, 3)', 'b'); 
xlim([-500 499]); ylim([-100 100])
% Set font size and font name for all axis labels and ticks
set(gca, 'FontSize', 24, 'FontName', 'Arial');
% Add axis labels and title
xlabel('Time (ms)', 'FontSize', 24, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 24, 'FontName', 'Arial');
title('EEG Signal After Re-Baseline Correction', 'FontSize', 24, 'FontName', 'Arial');

%% 5. Bandpass filtering.

EEG = pop_tesa_filtbutter(EEG, 0.1, [], 2, 'highpass');
EEG = tesa_filtbutter(EEG, [], 2000, 2, 'lowpass');

figure; pop_timtopo(EEG, [-500 499], [-200 0 100 200 300 400], 'After Re-Baseline Correction');

% Figure for paper
figure;
plot(EEG.times, mean(EEG.data, 3)', 'b'); 
xlim([-500 499]); ylim([-100 100])
% Set font size and font name for all axis labels and ticks
set(gca, 'FontSize', 24, 'FontName', 'Arial');
% Add axis labels and title
xlabel('Time (ms)', 'FontSize', 24, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 24, 'FontName', 'Arial');
title('EEG Signal After Re-Baseline Correction', 'FontSize', 24, 'FontName', 'Arial');

%% Save Dataset.

EEG = pop_saveset( EEG, 'filename',[save_name,'_clean.set'],'filepath',[path_dataset]);

%% End
