%% -----Remove High frequency noise (5kHz in iTEPs)----
%%%% Step1: design the gauss filter
% Design Gaussian Kernel
fs = EEG.srate;
fwhm_ms   = 0.16;   % requested FWHM in ms 
k         = 12;     % window half-length in samples around center
% time vector for kernel (ms)
gtime     = ( -k:k )/fs * 1e3;  % in ms 

% create Gaussian window
gausswin  = exp(-4*log(2)*(gtime/fwhm_ms).^2);
% compute empirical FWHM ()
gstHalf   = gtime(k+1:end);                  % post-peak segment
preHalfIdx = dsearchn(gausswin(1:k)', 0.5);
pstHalfIdx = k + dsearchn(gausswin(k+1:end)', 0.5);
% Full width at half maximum (FWHM)
empFWHM   = gtime(pstHalfIdx) - gtime(preHalfIdx);
% Normalize Gaussian to unit energy, area=1
gausswinN = gausswin / sum(gausswin);

%%%% Step2: apply to all trials
clearvars originalEEG gaulpEEG
originalEEG = double(EEG.data);
[channels,timepoints, trials] = size(originalEEG);
gaulpEEG = NaN(channels, timepoints, trials);
for ch = 1:channels
    for tr = 1:trials
        clear datain
        datain = squeeze(originalEEG(ch, :, tr));
        gaulpEEG(ch, :, tr) = filtfilt(gausswinN,1,datain);
    end
end

fprintf('\n=====================================\n');
fprintf('\nGaussian Low Pass, FWHM = %.2f ms \n\n',  empFWHM);
fprintf('=====================================\n\n');

EEG.data = gaulpEEG;
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'gui','off',...
   'overwrite','on', 'setname', 'LP_filtered');
eeglab redraw 