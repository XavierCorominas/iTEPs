%% Preprocessing pipeline for early iTEPs.

% DRCMR 2025 - XCT - BrainStim Methods group

% Preprocessing pipeline OPTION 2 - ADVANCED

% The preprocessing is inspired in the following references:
%https://www.sciencedirect.com/science/article/pii/S2666166725000280?via%3Dihub
%https://www.sciencedirect.com/science/article/pii/S1053811921005486

% --> GENERAL INFO:
% The current preprocessing pieplpile is intended to provide a elaborated and meticulous eeg
% preprocessing focusing on removing backgrond noise, TMS induced artifacts, muscular artifacts and eye movements. The preprocessing might be
% sufficient to explore iTEPS in datasets contaminated with muscular (or other type of) artifacts.

% ---> STEPS:
% 0) Load data
% 1) Epoch
% 2) De-mean data
% 3) Remove 50Hz noise
% 4) Remove bad trials
% 5) Remove TMS artifact
% 6) Interpolate removed TMS signal
% 7) Remove recharging artefact
% 8) First roud ICA - remove only eyes blinks and movements
% 9) High pass filter from 2Hz
% 10) Baseline correct the data
% 11) Remove recording noise and interpolate bad electrodes with SOUND
% 12) De-mean data
% 13) Supres time-locked muscular artifacts with SSP-SIR
% 14) Second round of ICA - remove persistent residual artifacts
% 15) Remove noisiest trials
% 16) Bandpass filter to segment fast and slow oscillation
% 17) Rereference to the average
% 18) Final baseline correction
% 19) Plot fast and slow oscillations
% 20) Save datasets


%% ADDPATHS 
clear all;
close all;
clc

% Path to external functions
addpath(genpath('/home/xavi/Documents/PROJECTS/iTEPS/eeg_analyses_tool/TMS_EEG_preprocessing/'))

% Path to eeglab
addpath('/home/xavi/Documents/PROJECTS/iTEPS/eeg_analyses_tool/TMS_EEG_preprocessing/external_functions/eeglab2025.0.0/')

% Path to fieltrip
addpath('/home/xavi/Documents/PROJECTS/iTEPS/eeg_analyses_tool/TMS_EEG_preprocessing/external_functions/fieldtrip-20240110/')

%
%% 0. LOAD DATA

% Open EEGlab 
eeglab;

% Define file to load 
name_dataset = 'AC_STATE_6.1.90CO25.vhdr';
path_dataset = '/mnt/projects/iTEPs/iTEPs_State/Pilots/RAW/EEG_Raw/AC/';

% Load file
EEG = pop_loadbv(path_dataset, name_dataset);

eeglab redraw


%% 1. EPOCH

% Epoch
EEG = pop_epoch( EEG, {  'R  8'  }, [-0.5         0.5]); % TMS pulses are stamped on the EEG as 'R  8' markers. Modify your marker if necessary for epoching.
eeglab redraw


% Plot data
figure; pop_timtopo(EEG, [-5  10], [2.3         4.5         4.8], 'ERP data and scalp maps');

% Figure for paper
figure;
plot(EEG.times, mean(EEG.data, 3)', 'b'); 
xlim([-5 10]); % ylim([-100 100])
% Set font size and font name for all axis labels and ticks
set(gca, 'FontSize', 24, 'FontName', 'Arial');
% Add axis labels and title (optional)
xlabel('Time (ms)', 'FontSize', 24, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 24, 'FontName', 'Arial');
%title('EEG Signal', 'FontSize', 24, 'FontName', 'Arial');

%% 2. DE MEAN

EEG = pop_rmbase( EEG, [-500  499.9]);

%Plot data
figure; pop_timtopo(EEG, [-5  10], [2.3         4.5         4.8], 'ERP data and scalp maps');

% Figure for paper
figure;
plot(EEG.times, mean(EEG.data, 3)', 'b'); 
xlim([-5 10]); % ylim([-100 100])
% Set font size and font name for all axis labels and ticks
set(gca, 'FontSize', 24, 'FontName', 'Arial');
% Add axis labels and title (optional)
xlabel('Time (ms)', 'FontSize', 24, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 24, 'FontName', 'Arial');
%title('EEG Signal', 'FontSize', 24, 'FontName', 'Arial');

%% 3. REMOVE LINE NOISE 50hZ FITTING A 50HZ OSCILLATION

EEG.data = fit_and_remove_line_noise_from_trials(EEG.data, EEG.srate, 50);

%Plot data
figure; pop_timtopo(EEG, [-5  10], [2.3         4.5         4.8], 'ERP data and scalp maps');

% Figure for paper
figure;
plot(EEG.times, mean(EEG.data, 3)', 'b'); 
xlim([-5 10]); % ylim([-100 100])
% Set font size and font name for all axis labels and ticks
set(gca, 'FontSize', 24, 'FontName', 'Arial');
% Add axis labels and title (optional)
xlabel('Time (ms)', 'FontSize', 24, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 24, 'FontName', 'Arial');
%title('EEG Signal', 'FontSize', 24, 'FontName', 'Arial');

%% 4. Remove bad trials
%% Trials rejection - 4.1 (identify bad trials) --> Select the bad trials

TMPREJ=[];
eegplot(EEG.data,'winlength',5,'command','pippo','srate',EEG.srate,'limits',[EEG.times(1) EEG.times(end)]);
'data have been displayed for first round of trials rejection'

 
%% Trials rejection - 4.2 (remove bad trials)

if ~isempty(TMPREJ)
    [trialrej elecrej]=eegplot2trial(TMPREJ,size(EEG.data,2),size(EEG.data,3));
else
    trialrej=[];
end

EEG.BadTr =find(trialrej==1);
EEG = pop_rejepoch( EEG, EEG.BadTr ,0); 

eeglab redraw


%Plot data
figure; pop_timtopo(EEG, [-5  10], [2.3         4.5         4.8], 'ERP data and scalp maps');

% Figure for paper
figure;
plot(EEG.times, mean(EEG.data, 3)', 'b'); 
xlim([-5 10]); % ylim([-100 100])
% Set font size and font name for all axis labels and ticks
set(gca, 'FontSize', 24, 'FontName', 'Arial');
% Add axis labels and title (optional)
xlabel('Time (ms)', 'FontSize', 24, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 24, 'FontName', 'Arial');
%title('EEG Signal', 'FontSize', 24, 'FontName', 'Arial');


%% 5. REMOVE TMS ARTIFACT
EEG = pop_tesa_removedata( EEG, [-1,2], [], {'R  8'} ); %remove from -1 to 2 ms
eeglab redraw 


%Plot data
figure; pop_timtopo(EEG, [-5  10], [2.3         4.5         4.8], 'ERP data and scalp maps');

% Figure for paper
figure;
plot(EEG.times, mean(EEG.data, 3)', 'b'); 
xlim([-5 10]); % ylim([-100 100])
% Set font size and font name for all axis labels and ticks
set(gca, 'FontSize', 24, 'FontName', 'Arial');
% Add axis labels and title (optional)
xlabel('Time (ms)', 'FontSize', 24, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 24, 'FontName', 'Arial');
%title('EEG Signal', 'FontSize', 24, 'FontName', 'Arial');

%% 6. TMS Interpolation

EEG = pop_tesa_interpdata( EEG, 'cubic', [1,2] );
eeglab redraw 

%Plot data
figure; pop_timtopo(EEG, [-5  10], [2.3         4.5         4.8], 'ERP data and scalp maps');

% Figure for paper
figure;
plot(EEG.times, mean(EEG.data, 3)', 'b'); 
xlim([-5 10]); % ylim([-100 100])
% Set font size and font name for all axis labels and ticks
set(gca, 'FontSize', 24, 'FontName', 'Arial');
% Add axis labels and title (optional)
xlabel('Time (ms)', 'FontSize', 24, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 24, 'FontName', 'Arial');
%title('EEG Signal', 'FontSize', 24, 'FontName', 'Arial');

%% 7. Recharging artefact

EEG = EEGLAB_remove_and_interpolate_recharging_artifact(EEG, 1); % 1 + single pulse, 2 = repetitive pulse

%% 8.*******First round of ICA --> identify oscular movements


% Resample for ica if necessary:
EEG = pop_resample( EEG, 25000);
eeglab redraw


%% Prepare the data for ICA: First remove for a while the very worst channels!

EEG_evo = mean(EEG.data,3);
[~, sigmas] = DDWiener(EEG_evo);
figure;
plot(sigmas,'*')


%% Labeling the very worst channels to not affect the ICA run

badC = find(sigmas > (median(sigmas) + 5*std(sigmas)));
goodC = setdiff(1:length(sigmas),badC);

EEG2ICA = pop_select( EEG, 'nochannel', badC);
badC


%% Before applying ICA to reduce artifacts, check which is the maximum number of independent components via singular value decomposition


tmpdata=reshape(EEG2ICA.data,[size(EEG2ICA.data,1) size(EEG2ICA.data,2)*size(EEG2ICA.data,3)]); 
tmpdata=tmpdata-repmat(mean(tmpdata,2),[1,size(EEG2ICA.data,2)*size(EEG2ICA.data,3)]);
tempCpca=svd(tmpdata);
th=0.01;% this threshold is arbitrary but it works in most cases
figure,semilogy(tempCpca,'.-')% plot of the eigenvalues of the data matrix
a=ylim;
hold on
Cpca=max(find(tempCpca>th));
line([Cpca Cpca],[a(1) a(2)],'Color','k','LineStyle','--');
hold off
'check whether the vertical dashed line correctly separates the eigenvalues close to zero from the ones much higher than zero, otherwise change the threshold th'

% if the vertical dotted line correctly separates the eigenvalues close to
% zero from the ones much higher than zero, Cpca correctly corresponds to
% the total number of independent components that can be estimated from the
% data matrix. if this is not the case, you have to change the threshold above or to manually correct Cpca


%% ICA - First round

% Set the number of dimentions accordingly
% You might need to go one below the previous value
% Here we only remove eye-blinks and horizontal movements: to make your
% choise, if unsure check: 
% https://link.springer.com/article/10.1007/s11517-011-0748-9
% https://www.sciencedirect.com/science/article/abs/pii/S105381191400620X?via%3Dihub
% https://www.sciencedirect.com/science/article/pii/S1053811916305845?via%3Dihub

EEG2ICA = pop_tesa_pcacompress( EEG2ICA, 'compVal', 60, 'plot','on' ); % Set the number of dimentions here
EEG2ICA = pop_tesa_fastica( EEG2ICA, 'approach', 'symm', 'g', 'tanh', 'stabilization', 'off' );
 

%% Plot the independent components following ICA. Here the user can verify the results of ICA
% as well as manually classify components. For details on classifying components, see
% https://nigelrogasch.gitbook.io/tesa-user-manual/remove_minimise_tms_muscle_activity/auto_comp_select
EEG2ICA = pop_tesa_compplot( EEG2ICA,'figSize','medium','plotTimeX',...
    [-500 499],'plotFreqX',[1 100], 'freqScale','log', 'saveWeights','off' );

EEG = EEG2ICA;


%Plot data
figure; pop_timtopo(EEG, [-5  10], [2.3         4.5         4.8], 'ERP data and scalp maps');

% Figure for paper
figure;
plot(EEG.times, mean(EEG.data, 3)', 'b'); 
xlim([-5 10]); % ylim([-100 100])
% Set font size and font name for all axis labels and ticks
set(gca, 'FontSize', 24, 'FontName', 'Arial');
% Add axis labels and title (optional)
xlabel('Time (ms)', 'FontSize', 24, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 24, 'FontName', 'Arial');
%title('EEG Signal', 'FontSize', 24, 'FontName', 'Arial');

%% 9. 2Hz High pass filter

EEG = pop_tesa_filtbutter( EEG, 2, [], 4, 'highpass' );

%Plot data
figure; pop_timtopo(EEG, [-5  10], [2.3         4.5         4.8], 'ERP data and scalp maps');

% Figure for paper
figure;
plot(EEG.times, mean(EEG.data, 3)', 'b'); 
xlim([-5 10]); % ylim([-100 100])
% Set font size and font name for all axis labels and ticks
set(gca, 'FontSize', 24, 'FontName', 'Arial');
% Add axis labels and title (optional)
xlabel('Time (ms)', 'FontSize', 24, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 24, 'FontName', 'Arial');
%title('EEG Signal', 'FontSize', 24, 'FontName', 'Arial');

%% 10. Baseline correction

EEG = pop_rmbase( EEG, [-210 -10] ,[]);

%Plot data
figure; pop_timtopo(EEG, [-5  10], [2.3         4.5         4.8], 'ERP data and scalp maps');

% Figure for paper
figure;
plot(EEG.times, mean(EEG.data, 3)', 'b'); 
xlim([-5 10]); % ylim([-100 100])
% Set font size and font name for all axis labels and ticks
set(gca, 'FontSize', 24, 'FontName', 'Arial');
% Add axis labels and title (optional)
xlabel('Time (ms)', 'FontSize', 24, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 24, 'FontName', 'Arial');
%title('EEG Signal', 'FontSize', 24, 'FontName', 'Arial');


%% 11. Deal with channel noise if necessary with SOUND
%make a copy
data_before_sound = EEG;

% Modify channel type name for runing sound
for i = 1:length(EEG.chanlocs)
EEG.chanlocs(i).type = 'EEG';
end 

% OPTION 1: USE SOUND
EEG = pop_tesa_sound(EEG, 'lambdaValue', 0.1, 'iter', 10 ); % spherical 3-layer reference model


%Plot data
figure; pop_timtopo(EEG, [-5  10], [2.3         4.5         4.8], 'ERP data and scalp maps');

% Figure for paper
figure;
plot(EEG.times, mean(EEG.data, 3)', 'b'); 
xlim([-5 10]); % ylim([-100 100])
% Set font size and font name for all axis labels and ticks
set(gca, 'FontSize', 24, 'FontName', 'Arial');
% Add axis labels and title (optional)
xlabel('Time (ms)', 'FontSize', 24, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 24, 'FontName', 'Arial');
%title('EEG Signal', 'FontSize', 24, 'FontName', 'Arial');

% Figure for paper
figure;
set(gcf, 'Position', [100, 100, 900, 492]); % Manually adjust figure size [x y width height]

plot(data_before_sound.times, mean(data_before_sound.data, 3)', 'r'); 
hold on
plot(EEG.times, mean(EEG.data, 3)', 'b'); 
xlim([-5 10]); % ylim([-100 100])
% Set font size and font name for all axis labels and ticks
set(gca, 'FontSize', 24, 'FontName', 'Arial');
% Add axis labels and title (optional)
xlabel('Time (ms)', 'FontSize', 24, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 24, 'FontName', 'Arial');
%title('EEG Signal', 'FontSize', 24, 'FontName', 'Arial');
%legend({'Before SOUND', 'After SOUND'}, 'FontSize', 20, 'Location', 'northeastoutside', 'Orientation', 'vertical');
hold off


% Figure for paper
figure;
set(gcf, 'Position', [100, 100, 900, 492]); % Manually adjust figure size [x y width height]

plot(data_before_sound.times, mean(data_before_sound.data, 3)', 'r'); 
hold on
plot(EEG.times, mean(EEG.data, 3)', 'b'); 
%xlim([-5 10]); % ylim([-100 100])
% Set font size and font name for all axis labels and ticks
set(gca, 'FontSize', 24, 'FontName', 'Arial');
% Add axis labels and title (optional)
xlabel('Time (ms)', 'FontSize', 24, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 24, 'FontName', 'Arial');
%title('EEG Signal', 'FontSize', 24, 'FontName', 'Arial');
%legend({'Before SOUND', 'After SOUND'}, 'FontSize', 20, 'Location', 'northeastoutside', 'Orientation', 'vertical');
hold off


%
%
%  OPTION 2: INTERPOLATE IT WITH EEGLAB. In the GUI go to TOOLS-->INTERPOLATE
% Interpolate electrode : check manually previous outsanding electrode and change it
%{
EEG = pop_interp(EEG, [7], 'spherical'); % Example to interpolate electrode numb 7

clear data

%Plot data
figure; pop_timtopo(EEG, [-5  10], [2.3         4.5         4.8], 'ERP data and scalp maps');

% Figure for paper
figure;
plot(EEG.times, mean(EEG.data, 3)', 'b'); 
xlim([-5 10]); % ylim([-100 100])
% Set font size and font name for all axis labels and ticks
set(gca, 'FontSize', 24, 'FontName', 'Arial');
% Add axis labels and title (optional)
xlabel('Time (ms)', 'FontSize', 24, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 24, 'FontName', 'Arial');
%title('EEG Signal', 'FontSize', 24, 'FontName', 'Arial');

%}

%% 12. De-mean the data

EEG = pop_rmbase( EEG, [-500 499] ,[]);

%% 13. SSP-SiR for muscular evoked activity

%make a copy
data_before_sspsir = EEG;

%
[EEG] = pop_tesa_sspsir(EEG, 'artScale', 'automatic');

%Plot data
figure; pop_timtopo(EEG, [-5  10], [2.3         4.5         4.8], 'ERP data and scalp maps');

% Figure for paper
figure;
plot(EEG.times, mean(EEG.data, 3)', 'b'); 
xlim([-5 10]); % ylim([-100 100])
% Set font size and font name for all axis labels and ticks
set(gca, 'FontSize', 24, 'FontName', 'Arial');
% Add axis labels and title (optional)
xlabel('Time (ms)', 'FontSize', 24, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 24, 'FontName', 'Arial');
%title('EEG Signal', 'FontSize', 24, 'FontName', 'Arial');


% Figure for paper
figure;
set(gcf, 'Position', [100, 100, 900, 492]); % Manually adjust figure size [x y width height]

plot(data_before_sspsir.times, mean(data_before_sspsir.data, 3)', 'k'); 
hold on
plot(EEG.times, mean(EEG.data, 3)', 'b'); 
xlim([-5 10]); % ylim([-100 100])
% Set font size and font name for all axis labels and ticks
set(gca, 'FontSize', 24, 'FontName', 'Arial');
% Add axis labels and title (optional)
xlabel('Time (ms)', 'FontSize', 24, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 24, 'FontName', 'Arial');
%title('EEG Signal', 'FontSize', 24, 'FontName', 'Arial');
%legend({'Before SSP-SIR', 'After SSP-SIR'}, 'FontSize', 20, 'Location', 'northeastoutside', 'Orientation', 'vertical');
hold off


% Figure for paper
figure;
set(gcf, 'Position', [100, 100, 900, 492]); % Manually adjust figure size [x y width height]

plot(data_before_sspsir.times, mean(data_before_sspsir.data, 3)', 'k'); 
hold on
plot(EEG.times, mean(EEG.data, 3)', 'b'); 
%xlim([-5 10]); % ylim([-100 100])
% Set font size and font name for all axis labels and ticks
set(gca, 'FontSize', 24, 'FontName', 'Arial');
% Add axis labels and title (optional)
xlabel('Time (ms)', 'FontSize', 24, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 24, 'FontName', 'Arial');
%title('EEG Signal', 'FontSize', 24, 'FontName', 'Arial');
%legend({'Before SSP-SIR', 'After SSP-SIR'}, 'FontSize', 20, 'Location', 'northeastoutside', 'Orientation', 'vertical');
hold off



% Figure for paper
figure;
set(gcf, 'Position', [100, 100, 900, 492]); % Manually adjust figure size [x y width height]

plot(data_before_sound.times, mean(data_before_sound.data, 3)', 'r'); 
hold on
plot(data_before_sspsir.times, mean(data_before_sspsir.data, 3)', 'k'); 

plot(EEG.times, mean(EEG.data, 3)', 'b'); 
%xlim([-5 10]); % ylim([-100 100])
% Set font size and font name for all axis labels and ticks
set(gca, 'FontSize', 24, 'FontName', 'Arial');
% Add axis labels and title (optional)
xlabel('Time (ms)', 'FontSize', 24, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 24, 'FontName', 'Arial');
%title('EEG Signal', 'FontSize', 24, 'FontName', 'Arial');
%legend({'Before SSP-SIR','Before SSP-SIR', 'After SSP-SIR'}, 'FontSize', 20, 'Location', 'northeastoutside', 'Orientation', 'vertical');
hold off


%% 14.*******Second round of ICA --> indentify residual muscular and persisting artifacts
%% Prepare the data for ICA: First remove for a while the very worst channels!

EEG_evo = mean(EEG.data,3);
[~, sigmas] = DDWiener(EEG_evo);
figure;
plot(sigmas,'*')


%% Labeling the very worst channels to not affect the ICA run

badC = find(sigmas > (median(sigmas) + 5*std(sigmas)));
goodC = setdiff(1:length(sigmas),badC);

EEG2ICA = pop_select( EEG, 'nochannel', badC);
badC


%% Before applying ICA to reduce artifacts, check which is the maximum number of independent components via singular value decomposition

tmpdata=reshape(EEG2ICA.data,[size(EEG2ICA.data,1) size(EEG2ICA.data,2)*size(EEG2ICA.data,3)]); 
tmpdata=tmpdata-repmat(mean(tmpdata,2),[1,size(EEG2ICA.data,2)*size(EEG2ICA.data,3)]);
tempCpca=svd(tmpdata);
th=0.01;% this threshold is arbitrary but it works in most cases
figure,semilogy(tempCpca,'.-')% plot of the eigenvalues of the data matrix
a=ylim;
hold on
Cpca=max(find(tempCpca>th));
line([Cpca Cpca],[a(1) a(2)],'Color','k','LineStyle','--');
hold off
'check whether the vertical dashed line correctly separates the eigenvalues close to zero from the ones much higher than zero, otherwise change the threshold th'

% if the vertical dotted line correctly separates the eigenvalues close to
% zero from the ones much higher than zero, Cpca correctly corresponds to
% the total number of independent components that can be estimated from the
% data matrix. if this is not the case, you have to change the threshold above or to manually correct Cpca


%% ICA - second round

% Set the number of dimentions accordingly
% You might need to go one below the previous value
% Here we  remove persistent muscular artifacts and noise
% choise, if unsure check: 
% https://link.springer.com/article/10.1007/s11517-011-0748-9
% https://www.sciencedirect.com/science/article/abs/pii/S105381191400620X?via%3Dihub
% https://www.sciencedirect.com/science/article/pii/S1053811916305845?via%3Dihub

EEG = pop_tesa_pcacompress( EEG, 'compVal', 55, 'plot','on' ); % Set the number of dimentions here
EEG = pop_tesa_fastica( EEG, 'approach', 'symm', 'g', 'tanh', 'stabilization', 'off' );
 

%% Plot the independent components following ICA. Here the user can verify the results of ICA
% as well as manually classify components. For details on classifying components, see
% https://nigelrogasch.gitbook.io/tesa-user-manual/remove_minimise_tms_muscle_activity/auto_comp_select
EEG = pop_tesa_compplot( EEG,'figSize','medium','plotTimeX',...
    [-500 499],'plotFreqX',[1 100], 'freqScale','log', 'saveWeights','off' );



%Plot data
figure; pop_timtopo(EEG, [-5  10], [2.3         4.5         4.8], 'ERP data and scalp maps');

% Figure for paper
figure;
plot(EEG.times, mean(EEG.data, 3)', 'b'); 
xlim([-5 10]); % ylim([-100 100])
% Set font size and font name for all axis labels and ticks
set(gca, 'FontSize', 24, 'FontName', 'Arial');
% Add axis labels and title (optional)
xlabel('Time (ms)', 'FontSize', 24, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 24, 'FontName', 'Arial');
%title('EEG Signal', 'FontSize', 24, 'FontName', 'Arial');



%% 15. Remove noisiest trials
%% Trials rejection - 15.1 (identify bad trials) --> Select the bad trials
TMPREJ=[];
eegplot(EEG.data,'winlength',5,'command','pippo','srate',EEG.srate,'limits',[EEG.times(1) EEG.times(end)]);
'data have been displayed for first round of trials rejection'

 
%% Trials rejection - 15.2 (remove bad trials)
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

%% 16.BANDPASS FILTER - segment fast and total fqs. to not afect each other on the baseline correction

EEGslow = mybutter(EEG, 2, 500, 2, 'bandpass'); 
EEGall = mybutter(EEG, 2, 1000, 2, 'bandpass');

%{
% Alternatively use TESA bandpass

EEGslow = pop_tesa_filtbutter( EEG, 2, 500, 2, 'bandpass' );;
EEGall = pop_tesa_filtbutter( EEG, 2, 1000, 2, 'bandpass' );;

%}
%% 17.REREFERENCING

EEGslow = pop_reref( EEGslow, []);
EEGall = pop_reref( EEGall, []);

%% 18. BASELINE CORRECTION

EEGslow = pop_rmbase( EEGslow, [-210 -10] ,[]);
EEGall = pop_rmbase( EEGall, [-210 -10] ,[]);
EEGfast = EEGall;
EEGfast.data = EEGall.data - EEGslow.data;


% Plot data
figure; pop_timtopo(EEG, [-5  10], [2.3         4.5         4.8], 'ERP data and scalp maps');

% Figure for paper
figure;
plot(EEG.times, mean(EEG.data, 3)', 'b'); 
xlim([-5 10]); % ylim([-100 100])
% Set font size and font name for all axis labels and ticks
set(gca, 'FontSize', 24, 'FontName', 'Arial');
% Add axis labels and title (optional)
xlabel('Time (ms)', 'FontSize', 24, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 24, 'FontName', 'Arial');
%title('EEG Signal', 'FontSize', 24, 'FontName', 'Arial');



%% 19. PLOT FINAL DATA. SEGMENT FAST OSCILLATION around iTEP and plot

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
title('Time series EEG Mean Signal Comparison', 'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Arial');

hold off

%% PLOT MEAN AND STD
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


% Plot time series of channels under the stimulation location; mean of 3 channels on the left motor cortex: 17, 16 , 6) ;
% Define time vector (assuming 250 points and 15 ms total)
num_points = size(EEG_fast.data, 2); % Number of time points
time_axis = linspace(-5, 10, num_points); % Adjusted time axis from -5 to 10 ms


% Data fast
data2 = EEG_fast.data;
mean_channel_data_roi(1,:,:) = data2(6,:,:);
mean_channel_data_roi(2,:,:) = data2(16,:,:);
mean_channel_data_roi(3,:,:) = data2(17,:,:);
mean_channel_data_roi = mean(mean_channel_data_roi, 1);
mean_eeg_fast = squeeze(mean_channel_data_roi)';

% Data slow
data3 = EEG_slow.data;
mean_channel_data_roi(1,:,:) = data3(6,:,:);
mean_channel_data_roi(2,:,:) = data3(16,:,:);
mean_channel_data_roi(3,:,:) = data3(17,:,:);
mean_channel_data_roi = mean(mean_channel_data_roi, 1);
mean_eeg_slow = squeeze(mean_channel_data_roi)';

% Data all
data4 = EEG_all.data;
mean_channel_data_roi(1,:,:) = data4(6,:,:);
mean_channel_data_roi(2,:,:) = data4(16,:,:);
mean_channel_data_roi(3,:,:) = data4(17,:,:);
mean_channel_data_roi = mean(mean_channel_data_roi, 1);
mean_eeg_all = squeeze(mean_channel_data_roi)';

% PLOT rhythmic conditions
x = (1:376);
figure() 

h1 = shadedErrorBar(x, mean(mean_eeg_fast,1), std(mean_eeg_fast), 'lineProps', 'b');
h1.mainLine.LineWidth = 2;
h1.mainLine.Color = [0 0 1]; % Blue
h1.patch.FaceColor = h1.mainLine.Color;

h2 = shadedErrorBar(x, mean(mean_eeg_slow,1), std(mean_eeg_slow), 'lineProps', 'r');
h2.mainLine.LineWidth = 2;
h2.mainLine.Color = [1 0 0]; % Red
h2.patch.FaceColor = h2.mainLine.Color;

h3 = shadedErrorBar(x, mean(mean_eeg_all,1), std(mean_eeg_all), 'lineProps', 'k');
h3.mainLine.LineWidth = 2;
h3.mainLine.Color = [0 0 0]; % Black
h3.patch.FaceColor = h3.mainLine.Color;

hold on
set(gcf, 'Position', [100, 100, 900, 492]); % Manually adjust figure size [x y width height]

% Manually set the positions of the ticks and labels
tick_positions = [0 125 250 376]; % Positions where the ticks will appear
tick_labels = {'-5', '0', '5', '10'}; % Labels for each tick

% Apply the tick positions and labels
xticks(tick_positions);      % Set tick positions
xticklabels(tick_labels);    % Set custom labels for each tick position

set(gca, 'FontSize', 20, 'FontName', 'Arial'); % Set font size and font for axes

% Optional: labels
xlabel('Time (ms)', 'FontSize', 20, 'FontName', 'Arial');
ylabel('Amplitude (µV)', 'FontSize', 20, 'FontName', 'Arial');

%
% Get current y-axis limits and adjust them (add +1 to the upper limit)
y_limits = ylim; 
ylim([y_limits(1), y_limits(2) + 1]);  % Increase the upper limit by 1

% Plot vertical grey bar from 0 ms to 1.5 ms
x_patch = [125 175 175 125];  % X-coordinates (start and end of the bar)
y_patch = [y_limits(1) y_limits(1) y_limits(2)+1 y_limits(2)+1]; % Y-coordinates based on updated y-limits
patch(x_patch, y_patch, [0.7 0.7 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.8, 'HandleVisibility', 'off'); % Grey with transparency

%
% Plot a black dotted line at y = 0
plot([0 376], [0 0], 'k--', 'LineWidth', 2); % Dotted black line at y = 0


% Adjust the legend location to be outside the plot area
legend({'500-1000Hz', '2-500Hz', '2-1000Hz'}, 'FontSize', 20, 'Location', 'northeastoutside', 'Orientation', 'vertical');

% Title
title('Time series EEG Mean+STD Signal Comparison', 'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Arial');

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



%% 20. SAVE DATASETS

EEG = pop_saveset( EEGall, 'filename',[name_dataset,'_all_oscill_cleaned_pipeline_2.set'],'filepath',[path_dataset]);
EEG = pop_saveset( EEGslow, 'filename',[name_dataset,'_slow_oscill_cleaned_pipeline_2.set'],'filepath',[path_dataset]);
EEG = pop_saveset( EEGfast, 'filename',[name_dataset,'_fast_oscill_cleaned_pipeline_2.set'],'filepath',[path_dataset]);


%% END
