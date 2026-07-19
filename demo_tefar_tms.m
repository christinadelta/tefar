% DEMO_TEFAR_TMS  Simulate TMS-EEG with known artefacts and run TEFAR_tms.
%
% Shows the full realistic path: simulate labelled sources -> mix to channels
% (x = A*s) -> hand the raw data to TEFAR_tms, which runs FastICA and flags
% artefact components. Because the ground truth is known, we can also check
% which flagged ICs correspond to which injected artefact.
%
% Requires FieldTrip on the path (for ft_componentanalysis / ft_freqanalysis).

clc; clear;

%% 0. FieldTrip on the path
addpath /Users/christinadelta/software/neuroscience/fieldtrip-20250414
ft_defaults

%% 1. simulate TMS-EEG with ground-truth artefacts
% profile 'tms' injects: line, blink, muscle, cardiac, decay, recharge
scfg          = [];
scfg.profile  = 'tms';
scfg.ntrials  = 40;      % number of epochs
scfg.fs       = 1000;    % sampling rate (Hz)
scfg.seed     = 7;       % reproducible
[comp_true, truth] = simulate_tefar_data(scfg);

fprintf('\nInjected sources (ground truth):\n');
for i = 1:numel(truth.label)
    fprintf('  source %2d : %s\n', i, truth.label{i});
end

%% 2. mix sources to channels: a normal FieldTrip raw data structure
% This is what your real recording looks like: channels x time, no labels.
data          = [];
data.label    = comp_true.topolabel;      % channel names
data.fsample  = comp_true.fsample;
data.time     = comp_true.time;
data.trial    = cell(1, numel(comp_true.trial));

for k = 1:numel(comp_true.trial)
    data.trial{k}   = comp_true.topo * comp_true.trial{k};   % x = A * s
end
data.trialinfo      = ones(numel(data.trial), 1);

%% 3. run the tool  (this runs FastICA internally, then classifies)
% trl is not needed by the core; pass [] for the demo.
[comp_tms, line_c, muscle_c, decay_c, addmuscle_c, recharge_c, blink_c, artifacts] = ...
    TEFAR_tms(data, []);

%% 4. inspect what it flagged
fprintf('\n--- TEFAR_tms suggestions ---\n');
fprintf('line       : %s\n', mat2str(line_c));
fprintf('muscle(topo): %s\n', mat2str(muscle_c));
fprintf('muscle(HF) : %s\n', mat2str(addmuscle_c));
fprintf('decay      : %s\n', mat2str(decay_c));
fprintf('recharge   : %s\n', mat2str(recharge_c));
fprintf('blink      : %s\n', mat2str(blink_c));
fprintf('cardiac    : %s\n', mat2str(artifacts.cardiac));
fprintf('=> suggested reject (union): %s\n', mat2str(artifacts.reject));

%% 5. prepare a layout for topoplots
% The simulated montage uses standard 10-20 names, so we build a layout from
% the same standard_1005 template my actual pipelines already load.
elec_file    = fullfile('/Users/christinadelta/software/neuroscience/fieldtrip-20250414', ...
                        'template','electrode','standard_1005.elc');
cfg_lay      = [];
cfg_lay.elec = ft_read_sens(elec_file);
lay          = ft_prepare_layout(cfg_lay);   % topoplot keeps only matching channels

%% 6. plot ALL components as topographies (like your original overview)
cfg           = [];
cfg.component = 1:numel(comp_tms.label);
cfg.layout    = lay;
cfg.colormap    = 'jet';
cfg.comment   = 'no';
figure('Name','All components');
ft_topoplotIC(cfg, comp_tms);

%% 7. plot only the FLAGGED components, titled by artefact category
cats   = {'line','muscle_topo','muscle','decay','recharge','blink','cardiac','eyemove'};
reject = artifacts.reject;
nrej   = numel(reject);
if nrej > 0
    ncols = ceil(sqrt(nrej)); nrows = ceil(nrej/ncols);
    figure('Name','Flagged components');
    for i = 1:nrej
        c = reject(i);
        % collect every category that flagged this component
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
        ft_topoplotIC(cfg, comp_tms);
        title(sprintf('IC %d: %s', c, strjoin(tags, '+')), 'Interpreter','none');
    end
    colormap(jet);
else
    fprintf('No components were flagged.\n');
end

%% 8. (optional) visualise the flagged components as time courses
cfg          = [];
comp_avg     = ft_timelockanalysis(cfg, comp_tms);
figure; cfg = []; cfg.viewmode = 'butterfly';
ft_databrowser(cfg, comp_avg);

%% 9. (optional) reject components and reconstruct the epochs
cfg              = [];
cfg.component    = artifacts.reject;
cfg.demean       = 'no';
data_clean       = ft_rejectcomponent(cfg, comp_tms);
fprintf('\nRejected %d components; data_clean is ready.\n', numel(artifacts.reject));
