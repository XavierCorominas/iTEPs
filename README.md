# **TMS-EEG Preprocessing Pipelines for iTEPS**

Yufei Song, Lasse Christiansen, Xavier Corominas
**Repository under construction**.




Danish Research Center for Magnetic Resonance (DRCMR) - 2025 - https://www.drcmr.dk/

The present repository contains **example preprocessing TMS-EEG scripts** to recover immediate transcranial evoked potentials (iTEP) from raw EEG recordings.

In the case the TMS pulse artifact wants to be conserved in the dataset, a **simple preprocessing** (e.g., **preprocessing_simple.m**) is necessary and sufficient when TMS-EEG data is not contaminated with large physiological and non-physiological artifacts. "Preprocessing simple.m" employs basic processing features including bandpass filtering. Given that TMS pulse is not removed, the filtering produces ringing signal artifacts.If that represents a serious concern for your analyses and conclusions, a simple processing but employing detrending instead of filtering is recomended (**preprocessing_simple_detrending.m**).


In the case the TMS pulse artifact wants to be removed from the dataset, a **basic preprocessing** (e.g., **preprocessing_basic.m**) is necessary and sufficient when TMS-EEG data is not contaminated with large physiological and non-physiological artifacts, and iTEPS are visible during evoked responses observed online during TMS experiments. A **advanced preprocessing** (e.g., **preprocessing_advanced.m**) is necessary and sufficient when TMS-EEG data is contaminated with large muscular or other types of artifacts overimposing to the iTEPS.

![eeg](https://github.com/user-attachments/assets/088f3c75-3389-4539-8fe7-463d663d705a)

For more information about iTEPS, please see:

https://www.brainstimjrnl.com/article/S1935-861X(24)00114-1/fulltext

https://www.sciencedirect.com/science/article/pii/S1053811925004495

No warranty or guarantee of any kind is provided. All rights reserved.

