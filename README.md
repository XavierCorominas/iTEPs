# **TMS-EEG Preprocessing Pipelines for iTEPS**

**Xavier Corominas-Teruel, Mikkel Malling Beck, Lasse Christiansen**

**Repository under construction**.




Danish Research Center for Magnetic Resonance (DRCMR) - 2025 - https://www.drcmr.dk/

The present repository contains **example preprocessing TMS-EEG scripts** to recover immediate transcranial evoked potentials (iTEP) from raw EEG recordings.

In the case the TMS pulse artifact wants to be conserved in the dataset, a **simple preprocessing** (e.g., preprocessing_simple.m) is necessary and sufficient when TMS-EEG data is not contaminated with large physiological and non-physiological artifacts. In the case the TMS pulse artifact wants to be removed from the dataset, a **basic preprocessing** (e.g., preprocessing_basic.m) is necessary and sufficient when TMS-EEG data is not contaminated with large physiological and non-physiological artifacts, and iTEPS are visible during evoked responses observed online during TMS experiments. A **advanced preprocessing** (e.g., preprocessing_advanced.m) is necessary and sufficient when TMS-EEG data is contaminated with large muscular or other types of artifacts overimposing to the iTEPS.

![Picture2](https://github.com/user-attachments/assets/ae51f46c-015b-4cea-9692-9f2786246e71)

For more information about iTEPS, please see:

https://www.brainstimjrnl.com/article/S1935-861X(24)00114-1/fulltext

https://www.biorxiv.org/content/10.1101/2025.02.14.638272v1.full


The scripts are currently under testing and may change in the future.

