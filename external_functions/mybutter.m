function EEG = mybutter(EEG, low_cutoff, high_cutoff, order, filter_type)

% Function to filter EEG data channel by channel and trial by trial using
% zero-phase butterworth filter 

    N_ch = EEG.nbchan; % Defining number of channels
    Fs = EEG.srate; % Defining sampling frequency
    N_epochs = size(EEG.data,3); % Defining number of epochs

    
 if strcmp(filter_type, 'lowpass')
        [b,a] = butter(order, high_cutoff / (Fs/2), 'low');
    elseif strcmp(filter_type, 'highpass')
        [b,a] = butter(order, low_cutoff / (Fs/2), 'high');
    elseif strcmp(filter_type, 'bandpass')
        [b,a] = butter(order, [low_cutoff high_cutoff] / (Fs/2));
    else
        error('Invalid filter type. Choose lowpass, highpass, or bandpass.');
    end

    % Looping over each channel and each epoch to filter the data
  
        for ch = 1:N_ch

            for epoch = 1:N_epochs

            EEG.data(ch,:,epoch) = filtfilt(b,a,double(EEG.data(ch,:,epoch)));

            end
        end

end