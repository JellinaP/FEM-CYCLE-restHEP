clear; clc
format shortG

addpath /data/u_prinsen_software/eeglab2024.2.1
eeglab nogui; 

basePath = '/data/pt_02938/rsHEP_analysis/';
cd(basePath)

method = 'icCFA';
files = dir(fullfile(basePath,'**', ['*postICA2_' method '.set']));

% BrainBeats HEP Feature extraction from continuous EEG recording
for i = 65:numel(files)
    fileCheck = dir(fullfile(basePath, files(i).folder, ['*_HEP2_' method '.set']));

    if isempty(fileCheck)
        EEG = pop_loadset(files(i).name, files(i).folder);
        EEG = eeg_checkset(EEG);

        % HEP MODE: Preprocess ECG signal and compute the HEP
        EEG = brainbeats_process(EEG, ...
            'analysis','hep','heart_signal','ecg','heart_channels',{ 'ECG' }, ...
            'clean_heart',1, 'rr_correct','pchip','clean_eeg',0, ...
            'parpool',0, 'gpu',0, 'vis_cleaning',0,'vis_outputs',0, 'save',0);

        % Save the file
        newFile = strsplit(files(i).name, '_');
        newFile = [strjoin(newFile(1:3),'_'), '_HEP2.set'];
        EEG = pop_saveset(EEG, 'filename', newFile, 'filepath', files(i).folder);
    end
end