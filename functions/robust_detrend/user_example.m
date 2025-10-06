% Robust detrend to remove slow drift in TMS-EEG data
%
% Functions developed by Johanna Metsomaa, BNP, University of Tuebingen 
% ref: Hernandez-Pavon, J. C., Kugiumtzis, D., Zrenner, C., Kimiskidis, V. K., & Metsomaa, J. (2022). 
% Removing artifacts from TMS-evoked EEG: A methods review and a unifying theoretical framework. 
% 
% require functions: 
% robust_detrend_EEG.m, 
% detrendAllData.m, 
% removePolyTrendlineTEP_robust
% 
% Note:
% 1- epoched data in eeglab structure
% epoch length is prefered to be long, if 2sec long, polynomial_order may
% need to be small
% 2- try out polynomial_order, eg., 2 or 3
% 3- [-20 600] ms is set to be the outlier segment. This can be changed in
% the robust_detrend_EEG.m
% 4- re-epoch your data after robust detrend is done to cut off edge
% artifacts, e.g., 500 ms at both the start and end of the epoch
% 5- Computationally heavy, better go through cluster

tic
polynomial_order = 3;
warning('off') % disable warnings 

EEG.data = double(EEG.data);
[EEG] =  robust_detrend_EEG(EEG, polynomial_order);
warning('on')

elapsedTime = toc; 
fprintf('Robust Detrend Done in : %.6f min\n', elapsedTime./60); %(in min)

