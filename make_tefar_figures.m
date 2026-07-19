function make_tefar_figures(profile, nseed, savedir)
% MAKE_TEFAR_FIGURES  Generate the validation figures for the TEFAR paper.
%
%   make_tefar_figures(profile, nseed, savedir)
%     profile : 'tms' (default) | 'eeg'
%     nseed   : number of random seeds for the performance panel (default 30)
%     savedir : if given, figures are saved there as PNG + PDF
%
% Produces three figures from the ground-truth simulation:
%   FIG 1  Artefact signatures : topography + power spectrum + time course,
%          one row per injected artefact class.
%   FIG 2  Feature separation  : the discriminating metric for each detector,
%          with the target artefact highlighted against the brain population
%          and the decision threshold drawn in.
%   FIG 3  Detection performance across seeds : sensitivity / specificity /
%          accuracy (boxplot) over nseed independent simulations.
%
% Requires FieldTrip on the path (topographies + ft_freqanalysis).

if nargin < 1 || isempty(profile), profile = 'tms'; end
if nargin < 2 || isempty(nseed),   nseed   = 30;    end
if nargin < 3, savedir = ''; end

FTPATH = '/Users/christinadelta/software/neuroscience/fieldtrip-20250414';
addpath(FTPATH); ft_defaults;

%% simulate (fixed seed for the illustrative panels)
[comp, truth] = simulate_tefar_data(struct('profile',profile,'ntrials',40,'fs',1000,'seed',7));

% layout from the standard template
lay_cfg      = []; lay_cfg.elec = ft_read_sens(fullfile(FTPATH,'template','electrode','standard_1005.elc'));
lay          = ft_prepare_layout(lay_cfg);

% spectra of the true sources
fcfg = []; fcfg.method='mtmfft'; fcfg.taper='hanning'; fcfg.output='pow'; fcfg.foilim=[1 100];
freq = ft_freqanalysis(fcfg, comp);

% one representative component per artefact class (first occurrence)
classes = artefact_order(truth.label, profile);
repidx  = zeros(1,numel(classes));
for i = 1:numel(classes)
    repidx(i) = find(strcmp(truth.label, classes{i}), 1, 'first');
end

% diagnostic band to shade per class (empty = time-domain feature)
band = containers.Map();
band('line')=[48 52]; band('muscle')=[20 90]; band('blink')=[1 4];
band('eyemove')=[1 4]; band('cardiac')=[]; band('decay')=[]; band('recharge')=[];

%% ---------------- FIG 1: artefact signatures ----------------
f1 = figure('Name','TEFAR artefact signatures','Color','w', ...
            'Position',[80 80 760 140*numel(classes)]);
for i = 1:numel(classes)
    c = repidx(i);

    % col 1: topography
    subplot(numel(classes),3,(i-1)*3+1);
    tcfg=[]; tcfg.component=c; tcfg.layout=lay; tcfg.comment='no';
    tcfg.figure='gca'; tcfg.interactive='no';   % draw into this subplot, no new window
    ft_topoplotIC(tcfg, comp); colormap(jet);
    % label the row by artefact NAME (overrides ft_topoplotIC's "component N";
    % ylabel is hidden because topoplots use axis off, so we use title)
    title(sprintf('%s (IC %d)', classes{i}, c), 'FontWeight','bold','Interpreter','none');

    % col 2: power spectrum (log-log) with diagnostic band shaded
    subplot(numel(classes),3,(i-1)*3+2);
    loglog(freq.freq, freq.powspctrm(c,:), 'k', 'LineWidth', 1.2); hold on;
    yl = ylim;
    if isKey(band, classes{i}) && ~isempty(band(classes{i}))
        b = band(classes{i});
        patch([b(1) b(2) b(2) b(1)],[yl(1) yl(1) yl(2) yl(2)], ...
              [1 .85 .3],'FaceAlpha',.35,'EdgeColor','none');
    end
    xlim([1 100]); grid on;
    if i==1, title('power spectrum'); end
    if i==numel(classes), xlabel('Hz'); end

    % col 3: example time course (trial 1)
    subplot(numel(classes),3,(i-1)*3+3);
    plot(comp.time{1}, comp.trial{1}(c,:), 'k'); hold on;
    if strcmpi(profile,'tms'), xline(0,'r--'); end   % TMS pulse
    axis tight; grid on;
    if i==1, title('time course'); end
    if i==numel(classes), xlabel('time (s)'); end
end
sgtitle(sprintf('Artefact signatures (%s profile)', upper(profile)));
save_fig(f1, savedir, sprintf('fig1_signatures_%s', profile));

%% ---------------- FIG 2: feature separation ----------------
res = verify_tefar_logic(profile);           % metrics + truth on seed 7
m   = res.metrics; lbl = res.truth.label; nc = numel(lbl);

% (metric vector, threshold, is-robust-z, target CLASS, y-label, panel title)
% Note: target class is a ground-truth label. 'muscle_topo' is a detector,
% not a class -- the muscle artefact is caught by both HF fraction and
% focality, so both panels target the 'muscle' class.
panels = {
    rz(m.line), 5,  true,  'line',   'line prominence (rz)', 'line'
    rz(m.emg),  5,  true,  'muscle', 'HF fraction (rz)',     'muscle (HF power)'
    rz(m.foc),  5,  true,  'muscle', 'focality (rz)',        'muscle (topography)'
    m.kurt,     4,  false, 'blink',  'kurtosis',             'blink'
    rz(m.per),  5,  true,  'cardiac','periodicity (rz)',     'cardiac'
};
if strcmpi(profile,'tms')
    % decay & recharge are temporal: window RMS relative to baseline SD
    panels(end+1,:) = { m.decay ./ max(m.bsd,eps), 2, false, 'decay',    'decay RMS / baseline SD',    'decay' };
    panels(end+1,:) = { m.rech  ./ max(m.bsd,eps), 4, false, 'recharge', 'recharge RMS / baseline SD', 'recharge' };
else
    panels(end+1,:) = { m.asym, 0.15, false, 'eyemove', 'front. asymmetry', 'eyemove' };
end

f2 = figure('Name','TEFAR feature separation','Color','w','Position',[100 100 900 520]);
np = size(panels,1); ncol = 3; nrow = ceil(np/ncol);
for p = 1:np
    val = panels{p,1}; thr = panels{p,2}; tgt = panels{p,4};
    ylab = panels{p,5}; ptitle = panels{p,6};
    isTgt  = strcmp(lbl, tgt);
    subplot(nrow,ncol,p); hold on;
    plot(find(~isTgt), val(~isTgt), 'o', 'Color',[.6 .6 .6], 'MarkerFaceColor',[.85 .85 .85]);
    plot(find(isTgt),  val(isTgt),  'o', 'Color',[.8 0 0],   'MarkerFaceColor',[1 .3 .3],'MarkerSize',8);
    yline(thr,'r--','threshold');
    xlim([0 nc+1]); grid on; ylabel(ylab,'Interpreter','none');
    title(ptitle,'Interpreter','none'); if p> (nrow-1)*ncol, xlabel('component #'); end
end
sgtitle(sprintf('Feature separation: artefact (red) vs brain (grey) — %s', upper(profile)));
save_fig(f2, savedir, sprintf('fig2_separation_%s', profile));

%% ---------------- FIG 3: performance across seeds ----------------
S = zeros(nseed,3);
for s = 1:nseed
    evalc('rr = verify_tefar_logic(profile, s);');     % suppress printout
    S(s,:) = [rr.sens rr.spec rr.acc];
end
% toolbox-free "boxplot": jittered per-seed points + mean +/- SD
f3 = figure('Name','TEFAR performance','Color','w','Position',[120 120 520 420]);
hold on;
names = {'sensitivity','specificity','accuracy'};
for j = 1:3
    xj = j + (rand(nseed,1)-0.5)*0.18;                 % jitter
    plot(xj, S(:,j), 'o', 'Color',[.6 .6 .6], 'MarkerFaceColor',[.85 .85 .85], 'MarkerSize',4);
    mu = mean(S(:,j)); sd = std(S(:,j));
    plot([j-0.25 j+0.25], [mu mu], 'k', 'LineWidth', 2);          % mean
    plot([j j], [mu-sd mu+sd], 'k', 'LineWidth', 1);              % +/- SD
    text(j+0.30, mu, sprintf('%.3f', mu), 'FontSize',9);
end
set(gca,'XTick',1:3,'XTickLabel',names); xlim([0.5 3.5]); ylim([0 1.03]);
grid on; ylabel('score'); yline(1,'k:');
title(sprintf('Detection performance over %d seeds (%s)', nseed, upper(profile)));
fprintf('\n%s: mean sens=%.3f spec=%.3f acc=%.3f (n=%d seeds)\n', ...
        upper(profile), mean(S(:,1)), mean(S(:,2)), mean(S(:,3)), nseed);
save_fig(f3, savedir, sprintf('fig3_performance_%s', profile));
end

% ============================ helpers ==================================
function order = artefact_order(labels, profile)
% stable, sensible ordering of the artefact classes present
if strcmpi(profile,'tms')
    pref = {'line','blink','muscle','cardiac','decay','recharge'};
else
    pref = {'line','blink','eyemove','muscle','cardiac'};
end
order = pref(ismember(pref, unique(labels)));
end

function y = rz(x)
med = median(x); madv = median(abs(x-med));
y = (x-med) ./ (1.4826*madv + eps);
end

function save_fig(h, savedir, name)
if isempty(savedir), return; end
if ~exist(savedir,'dir'), mkdir(savedir); end
try
    exportgraphics(h, fullfile(savedir,[name '.png']), 'Resolution',200);
    exportgraphics(h, fullfile(savedir,[name '.pdf']));
catch
    print(h, fullfile(savedir,[name '.png']), '-dpng','-r200');  % older MATLAB
    print(h, fullfile(savedir,[name '.pdf']), '-dpdf');
end
end