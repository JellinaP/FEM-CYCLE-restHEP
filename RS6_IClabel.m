clc; clear
basePath = '/data/pt_02938/rsHEP_analysis/';
cd(basePath)

addpath /data/u_prinsen_software/LIFE_scripts_Paul/hep_rest/functions_settings/
addpath /data/pt_02938/rsHEP_analysis/1_scripts/helpers/
addpath /data/pt_02938/rsHEP_analysis/1_scripts/MINT_RS_EEG_scripts/helperFunctions/
addpath /data/u_prinsen_software/eeglab2024.2.1
eeglab('nogui');

% === User choice: which ICA rejection method to apply ===
% 'iclabel'   - use ICLabel only (eye + heart IC rejection via ICflag)
%               + Alpha power protection (preserves ICs with >50% alpha power)
% 'icCFA'     - use ICLabel + ECG template correlation for heart ICs and CFA
%               + Alpha power protection (preserves ICs with >50% alpha power)
%
% ALPHA PROTECTION: Components classified as cardiac BUT with >50% power in
% alpha band (7-15 Hz) are PRESERVED, as they may represent true brain-HEP
% activity that is heartbeat-locked rather than pure cardiac artifacts.
%
method = 'icCFA';
QC_dashboard = "no";
PSD_plot = "no";

% IClabel rejection thresholds
thr = struct();
thr.Brain        = [0 0];
thr.Muscle       = [0.9 1];
thr.Eye          = [0.7 1];
thr.Heart        = [0.9 1];
thr.Line         = [0.9 1];
thr.Channel      = [0.9 1];
thr.Other        = [0.95 1];

% time-window for ECG template building
time_window = [-0.2 0.2]; %in s, cfr Buot et al., 2021
% ========================================================

% Start group iteration
files = dir(fullfile(basePath, '**', '*_clean2.set'));

for si = 1:length(files)
    fileCheck = dir(fullfile(files(si).folder, ['*_postICA2_' method '.set']));

    if isempty(fileCheck)
        EEG_orig = pop_loadset(files(si).name, files(si).folder);
        EEG = pop_select(EEG_orig, 'channel', (1:EEG_orig.nbchan-1));  % extract EEG+vEOG only
        
        ECG = pop_select(EEG_orig, 'channel', EEG_orig.nbchan);    % extract and store ECG
        ECG.event = []; % start clean to avoid event duplication

        amicaPath = fullfile(files(si).folder, 'AMICA');
        EEG.etc.amica  = loadmodout15(amicaPath);

        % --- Assign AMICA weights
        EEG.etc.amica.S = EEG.etc.amica.S(1:EEG.etc.amica.num_pcs, :); % Weirdly, I saw size(S,1) be larger than rank. This process does not hurt anyway.
        EEG.icaweights = EEG.etc.amica.W;
        EEG.icasphere  = EEG.etc.amica.S;
        EEG.icawinv = pinv(EEG.icaweights * EEG.icasphere);
        EEG = eeg_checkset(EEG, 'ica');

        % --- Label ICA components
        EEG = iclabel(EEG);
        EEG = eeg_checkset(EEG);

        % --- Run selected rejection method
        switch method
            case 'iclabel'
                EEG = pop_icflag(EEG, [thr.Brain; thr.Muscle; thr.Eye; thr.Heart; ...
                    thr.Line; thr.Channel; thr.Other]);

                heartICs_iclabel = find(EEG.etc.ic_classification.ICLabel.classifications(:,4) >= thr.Heart(1));

                % Alpha power protection for ICLabel-only method
                alphaProtected = [];
                if ~isempty(heartICs_iclabel)
                    fprintf('  Checking alpha power for %d cardiac ICs (ICLabel method)...\n', length(heartICs_iclabel));

                    for ic_idx = 1:length(heartICs_iclabel)
                        comp = heartICs_iclabel(ic_idx);

                        % Get continuous IC activity for this component
                        IC_continuous = EEG.icaact(comp, :);

                        % Compute relative power in alpha band (7-15 Hz)
                        alpha_rel_pow = compute_RP(IC_continuous', EEG_orig.srate, [1 45], [7 15]);

                        % If >50% power in alpha band, likely brain activity (potential HEP)
                        if alpha_rel_pow > 0.5
                            alphaProtected(end+1) = comp;
                            fprintf('    IC %d: alpha_RP = %.2f (PROTECTED - potential HEP)\n', comp, alpha_rel_pow);
                        end
                    end

                    % Remove alpha-protected components from cardiac removal list
                    if ~isempty(alphaProtected)
                        heartICs_iclabel = setdiff(heartICs_iclabel, alphaProtected);
                        EEG.etc.ic_protected.alpha = alphaProtected;
                        EEG.etc.ic_protected.alpha_reason = 'High alpha power (>0.5 RP), likely HEP';
                        fprintf('  Protected %d ICs with high alpha power from cardiac removal\n', length(alphaProtected));
                    end
                end

                EEG.etc.ic_removed.Heart = heartICs_iclabel;

            case 'icCFA'
                % Step 1: standard rejection based on IClabel
                EEG = pop_icflag(EEG, [thr.Brain; thr.Muscle; thr.Eye; thr.Heart; ...
                    thr.Line; thr.Channel; thr.Other]);

                % Step 2: ECG template correlation
                % get RR times using a fast R-peak detector
                RR_samples = local_detect_rpeaks(ECG.data, EEG_orig.srate);  % returns sample

                % add R events to ECG
                ECG.event = struct('type',{},'latency',{},'duration',{});
                for b = 1:numel(RR_samples)
                    ECG.event(end+1).type = 'R';
                    ECG.event(end).latency = RR_samples(b); % in samples
                    ECG.event(end).duration = 0;
                end
                ECG = eeg_checkset(ECG, 'eventconsistency');

                % create ECG template
                ECG_template = pop_epoch(ECG, 'R' , time_window);
                ECG_template = mean(ECG_template.data,3); % [1 x samples] template

                % add R events to ICA-labeled EEG
                ICA = EEG;
                ICA.data = EEG.icaact;
                for b = 1:numel(RR_samples)
                    ICA.event(end+1).type = 'R';
                    ICA.event(end).latency = RR_samples(b); % in samples
                    ICA.event(end).duration = 0;
                end

                % pad ICA data with dummy channels for easier epoching
                if size(EEG.data,1) ~= size(ICA.data,1)
                    add_chans =  size(EEG.data,1) - size(ICA.data,1);
                    ICA.data = [ICA.data; zeros(add_chans, length(ICA.data))];
                end
                ICA = eeg_checkset(ICA, 'eventconsistency');

                % epoch ICA data
                ICA = pop_epoch(ICA, 'R' , time_window);
                ICA.data(end-add_chans:end,:,:) = [];
                ICA = mean(ICA.data,3);

                % correlation ICs with ECG template
                SDcrit = 1.5;
                corrVals = zeros(size(ICA,1),1);
                for ic = 1:size(ICA,1)
                    corrVals(ic) = abs(corr(ICA(ic,:)', ECG_template'));
                end
                corthreshV = mean(corrVals) + SDcrit * std(corrVals);
                EEG.etc.ic_ecg_corr = corrVals;

                % identify ICs over adaptive threshold and with template correlation >.8
                heartICs_template = find((corrVals > corthreshV) + (corrVals>0.8) == 2);

                % exclude brain-like ICs
                heartICs_template = heartICs_template(~(EEG.etc.ic_classification.ICLabel.classifications(heartICs_template,1)>0.5));

                % Step 3: Alpha power protection (preserve HEP-related components)
                % Check if cardiac-correlated ICs have high alpha power (may be HEP)
                alphaProtected = [];
                if ~isempty(heartICs_template)
                    fprintf('  Checking alpha power for %d cardiac ICs...\n', length(heartICs_template));

                    for ic_idx = 1:length(heartICs_template)
                        comp = heartICs_template(ic_idx);

                        % Get IC activity for this component (already epoched and averaged)
                        IC_data = ICA(comp, :);

                        % Compute relative power in alpha band (7-15 Hz)
                        % Parameters: signal, fs, frequency_range, bands_of_interest
                        alpha_rel_pow = compute_RP(IC_data', EEG_orig.srate, [1 45], [7 15]);

                        % If >50% power in alpha band, likely brain activity (potential HEP)
                        if alpha_rel_pow > 0.5
                            alphaProtected(end+1) = comp;
                            fprintf('    IC %d: alpha_RP = %.2f (PROTECTED - potential HEP)\n', comp, alpha_rel_pow);
                        end
                    end

                    % Remove alpha-protected components from cardiac removal list
                    if ~isempty(alphaProtected)
                        heartICs_template = setdiff(heartICs_template, alphaProtected);
                        EEG.etc.ic_protected.alpha = alphaProtected;
                        EEG.etc.ic_protected.alpha_reason = 'High alpha power (>0.5 RP), likely HEP';
                        fprintf('  Protected %d ICs with high alpha power from cardiac removal\n', length(alphaProtected));
                    end
                end

                % merge ICLabel + template
                heartICs_iclabel = find(EEG.etc.ic_classification.ICLabel.classifications(:,4) >= thr.Heart(1));
                EEG.etc.ic_removed.Heart = unique([heartICs_iclabel; heartICs_template]);

                % mark for rejection
                EEG.reject.gcompreject(EEG.etc.ic_removed.Heart) = 1;
        end

        % --- Log all other categories (same for all methods)
        EEG.etc.ic_n = numel(find(EEG.reject.gcompreject > 0));
        EEG.etc.ic_method = method;
        
        EEG.etc.ic_classification.ICLabel.thresholds.Eye     = thr.Eye(1);
        EEG.etc.ic_classification.ICLabel.thresholds.Heart   = thr.Heart(1);
        EEG.etc.ic_classification.ICLabel.thresholds.Muscle  = thr.Muscle(1);
        EEG.etc.ic_classification.ICLabel.thresholds.Line    = thr.Line(1);
        EEG.etc.ic_classification.ICLabel.thresholds.Channel = thr.Channel(1);
        EEG.etc.ic_classification.ICLabel.thresholds.Other   = thr.Other(1);

        EEG.etc.ic_removed.Muscle     = find(EEG.etc.ic_classification.ICLabel.classifications(:,2) >= thr.Muscle(1));
        EEG.etc.ic_removed.Eye        = find(EEG.etc.ic_classification.ICLabel.classifications(:,3) >= thr.Eye(1));
        EEG.etc.ic_removed.Line       = find(EEG.etc.ic_classification.ICLabel.classifications(:,5) >= thr.Line(1));
        EEG.etc.ic_removed.Channel    = find(EEG.etc.ic_classification.ICLabel.classifications(:,6) >= thr.Channel(1));
        EEG.etc.ic_removed.Other      = find(EEG.etc.ic_classification.ICLabel.classifications(:,7) >= thr.Other(1));
        
        if strcmpi(string(QC_dashboard), "yes")
            % --- Make a copy with removed IC still included for dasboard
            EEG_qc = EEG;
            % Add ECG again to QC copy
            EEG_qc.data(EEG_qc.nbchan+1,:) = ECG;
            EEG_qc.nbchan = EEG_qc.nbchan+1;
            EEG_qc.chanlocs(end+1).labels = 'ECG';EEG_qc.chanlocs(end).type   = 'ECG';

            % -- Create IC QC html dashboard
            generate_ica_qc_html_dashboard(EEG_qc, files(si).folder, [EEG.subject '_ses-' EEG.session '_ICA_dashboard']);
        end

        % --- Remove marked ICs, then remove vEOG channel
        EEG = pop_subcomp(EEG, find(EEG.reject.gcompreject), 0);
        EEG = pop_select(EEG, 'rmchannel', EEG.nbchan);  % remove vEOG
        EEG = eeg_checkset(EEG);

        % --- Plot PSD after removal of flagged components
        if strcmpi(string(PSD_plot), "yes")
            fig = figure('Visible','off');
            plot_spec(EEG.data(1:EEG.nbchan,:)', EEG.srate, 'f_max', 70);
            saveas(fig, fullfile(files(si).folder, [files(si).name(1:5), '_PSD_post_ICA.png']), 'png');
            close(fig);
        end

        % --- Add back ECG channel
        if size(EEG.data, 2) == size(ECG.data, 2)
            EEG.data(EEG.nbchan+1,:) = ECG.data(1,:);
            EEG.nbchan = EEG.nbchan+1;
            EEG.chanlocs(end+1).labels = 'ECG'; EEG.chanlocs(end).type = 'ECG';
        end
        EEG = eeg_checkset(EEG);

        % --- Save
        pop_saveset(EEG, [files(si).name(1:9), '_postICA2_' method '.set'], files(si).folder);
    end
end