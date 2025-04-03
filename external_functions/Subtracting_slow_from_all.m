% Define your time range
time_range = [-5, 10];

% Find the indices within this range
idx_range = find(EEGfast.times >= time_range(1) & EEGfast.times <= time_range(2));

% Extract the EEG data within this time range
EEG_slow = X35488_M1_PAAP_cubic_01_500.data(:, idx_range,:);
EEG_all = X35488_M1_PAAP_cubic_01_1000.data(:, idx_range,:);

new_time_vector = EEGfast.times(idx_range);

%%
fast = EEGall.data - EEGslow.data; 


%%
EEG_fast = EEG_all - EEG_slow;

%%
figure;
hold on
plot(EEG.times, mean(EEGslow.data(17,:,:),3))
plot(EEG.times, mean(EEGall.data(17,:,:),3))
xlim([-5 10])
hold off 

figure;
plot(EEG.times, mean(fast(17,:,:),3))
xlim([-5 10])