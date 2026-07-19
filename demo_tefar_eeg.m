% DEMO_TEFAR_EEG  Simulate ordinary EEG with known artefacts and run TEFAR_eeg.
%
% Same idea as demo_tefar_tms, but for non-TMS EEG. The 'eeg' profile injects:
% line, blink, lateral eye-movement, muscle, cardiac (no decay/recharge).
% We simulate labelled sources -> mix to channels (x = A*s) -> hand the raw
% data to TEFAR_eeg, which runs FastICA and flags artefact components.
%
% Requires FieldTrip on the path (for ft_componentanalysis / ft_freqanalysis).

clc; clear;

%% 0. FieldTrip on the path
addpath /Users/christinadelta/software/neuroscience/fieldtrip-20250414
ft_defaults

%% 1. Simulate EEG with ground-truth artefacts
scfg          = [];
scfg.profile  = 'eeg';   % injects: line, blink, eyemove, muscle, cardiac
scfg.ntrials  = 40;      % number of epochs
scfg.fs       = 1000;    % sampling rate (Hz)
scfg.seed     = 7;       % reproducible
[comp_true, truth] = simulate_tefar_data(scfg);

fprintf('\nInjected sources (ground truth):\n');
for i = 1:numel(truth.label)
    fprintf('  source %2d : %s\n', i, truth.label{i});
end

%% 2. Mix sources to channels  ->  a normal FieldTrip raw data structure
data          = [];
data.label    = comp_true.topolabel;      % channel names
data.fsample  = comp_true.fsample;
data.time     = comp_true.time;
data.trial    = cell(1, numel(comp_true.trial));
for k = 1:numel(comp_true.trial)
    data.trial{k} = comp_true.topo * comp_true.trial{k};   % x = A * s
end
data.trialinfo = ones(numel(data.trial), 1);

%% 3. Run the tool  (runs FastICA internally, then classifies)
% Note the simpler EEG signature: [comp, artifacts].
[comp_eeg, artifacts] = TEFAR_eeg(data);

%% 4. Inspect what it flagged
fprintf('\n--- TEFAR_eeg suggestions ---\n');
fprintf('line        : %s\n', mat2str(artifacts.line));
fprintf('blink       : %s\n', mat2str(artifacts.blink));
fprintf('eyemove     : %s\n', mat2str(artifacts.eyemove));
fprintf('muscle(HF)  : %s\n', mat2str(artifacts.muscle));
fprintf('muscle(topo): %s\n', mat2str(artifacts.muscle_topo));
fprintf('cardiac     : %s\n', mat2str(artifacts.cardiac));
fprintf('=> suggested reject (union): %s\n', mat2str(artifacts.reject));

%% 5. Prepare a layout for topoplots
% The simulated montage uses standard 10-20 names, so we build a layout from
% the same standard_1005 template your pipeline already loads.
elec_file    = fullfile('/Users/christinadelta/software/neuroscience/fieldtrip-20250414', ...
                        'template','electrode','standard_1005.elc');
cfg_lay      = [];
cfg_lay.elec = ft_read_sens(elec_file);
lay          = ft_prepare_layout(cfg_lay);

%% 6. Plot ALL components as topographies (overview)
cfg           = [];
cfg.component = 1:numel(comp_eeg.label);
cfg.layout    = lay;
cfg.comment   = 'no';
figure('Name','All components');
ft_topoplotIC(cfg, comp_eeg);

%% 7. Plot only the FLAGGED components, titled by artefact category
cats   = {'line','muscle_topo','muscle','blink','eyemove','cardiac'};
reject = artifacts.reject;
nrej   = numel(reject);
if nrej > 0
    ncols = ceil(sqrt(nrej)); nrows = ceil(nrej/ncols);
    figure('Name','Flagged components');
    for i = 1:nrej
        c = reject(i);
        tags = {};
        for j = 1:numel(cats)
            if isfield(artifacts, cats{j}) && any(artifacts.(cats{j}) == c)
                tags{end+1} = cats{j}; %#ok<SAGROW>
            end
        end
        subplot(nrows, ncols, i);
        cfg           = [];
        cfg.component = c;
        cfg.layout    = lay;
        cfg.comment   = 'no';
        cfg.figure    = 'gca';    % draw into this subplot, no new window
        cfg.interactive = 'no';
        ft_topoplotIC(cfg, comp_eeg);
        title(sprintf('IC %d: %s', c, strjoin(tags, '+')), 'Interpreter','none');
    end
    colormap(jet);
else
    fprintf('No components were flagged.\n');
end

%% 8. (optional) browse the components as time courses
cfg          = [];
cfg.viewmode = 'component';
cfg.layout   = lay;
ft_databrowser(cfg, comp_eeg);

%% 9. (optional) reject and reconstruct
cfg           = [];
cfg.component = artifacts.reject;
cfg.demean    = 'no';
data_clean    = ft_rejectcomponent(cfg, comp_eeg);
fprintf('\nRejected %d components; data_clean is ready.\n', numel(artifacts.reject));
