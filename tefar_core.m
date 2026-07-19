function [comp, artifacts] = tefar_core(cfg, data)
% TEFAR_CORE  Configurable ICA-based artefact classification for (TMS-)EEG.
%
%   [comp, artifacts] = tefar_core(cfg, data)
%
% Runs ICA (default FastICA) on FieldTrip-formatted data and scores every
% independent component against a configurable set of artefact detectors.
% It is the shared engine behind TEFAR_tms (TMS-EEG) and TEFAR_eeg (EEG).
%
% Design goals
%   * Generic: no assumption about trial/segment structure. Every temporal
%     detector locates its window from comp.time, per trial.
%   * FieldTrip-only: no MATLAB toolboxes required (kurtosis/skewness are
%     implemented locally; see subfunctions).
%   * Data-driven: spectral/topographic thresholds are expressed as z-scores
%     across the component distribution, with sensible absolute fallbacks.
%
% INPUT
%   data   FieldTrip raw/epoched data (output of ft_preprocessing), OR omit
%          and pass precomputed components via cfg.comp.
%
% CONFIGURATION (all optional; defaults in brackets)
%   cfg.comp          precomputed components (skips ICA)                 []
%   cfg.method        ICA method                                 ['fastica']
%   cfg.fastica       struct of FastICA options   [approach='symm', g='gauss']
%   cfg.demean        demean before ICA                             ['yes']
%   cfg.numcomponent  number of components to estimate                  []
%   cfg.detectors     cell array of detectors to enable. Any of:
%                     {'line','muscle','muscle_topo','blink','eyemove',
%                      'cardiac','decay','recharge'}
%                     Default: {'line','muscle','muscle_topo','blink'}
%   cfg.reject_rule   'union' (any detector flags) | 'score'      ['union']
%   cfg.reject_threshold  score needed if reject_rule=='score'          [2]
%   cfg.verbose       print a per-component report                  [true]
%
%   Spectral / mains
%   cfg.line_freq     mains frequency (Hz)                             [50]
%   cfg.line.bw       half-bandwidth around mains & harmonics (Hz)      [2]
%   cfg.line.z        z-threshold on line-band power                     [4]
%   cfg.muscle.band   EMG band [lo hi] (hi clipped to Nyquist)     [20 90]
%   cfg.muscle.z      z-threshold on EMG-band power                      [4]
%   cfg.blink.band    low-freq blink band                            [1 4]
%   cfg.blink.kurt    kurtosis threshold                                 [5]
%   cfg.blink.frontal_ratio  frontal topo concentration threshold     [0.5]
%   cfg.cardiac.band  cardiac fundamental band (Hz)                 [0.8 2]
%   cfg.cardiac.z     z-threshold on rhythmicity                         [3]
%
%   Topography
%   cfg.topo.z        z-threshold on peak topo weight (focal/muscle)     [3]
%   cfg.frontal_labels    frontal channel names (see defaults below)
%
%   Temporal (TMS)
%   cfg.baseline      pre-event window for baseline std (s)     [-Inf 0]
%   cfg.decay.window     early post-pulse window (s)         [0.010 0.050]
%   cfg.decay.k          multiples of baseline std                       [2]
%   cfg.recharge.window  recharge window (s)                 [0.100 0.500]
%   cfg.recharge.k       multiples of baseline std                       [4]
%
% OUTPUT
%   comp        the FieldTrip component structure (from ft_componentanalysis)
%   artifacts   struct with fields:
%     .line .muscle .muscle_topo .blink .eyemove .cardiac .decay .recharge
%                 -> vectors of component indices flagged by each detector
%     .reject     -> final suggested rejection list (per cfg.reject_rule)
%     .score      -> [1 x ncomp] number of detectors flagging each component
%     .metrics    -> struct of per-component metric vectors
%     .cfg        -> the fully-populated cfg actually used
%
% @christinadelta - 2026 
% ---------------------------------------------------------------- defaults
if nargin < 1 || isempty(cfg), cfg = struct(); end

cfg           = ft_setdefault(cfg, 'method',      'fastica');
cfg           = ft_setdefault(cfg, 'demean',      'yes');
cfg           = ft_setdefault(cfg, 'comp',        []);
cfg           = ft_setdefault(cfg, 'numcomponent',[]);
cfg           = ft_setdefault(cfg, 'detectors',   {'line','muscle','muscle_topo','blink'});
cfg           = ft_setdefault(cfg, 'reject_rule', 'union');
cfg           = ft_setdefault(cfg, 'reject_threshold', 2);
cfg           = ft_setdefault(cfg, 'verbose',     true);
cfg           = ft_setdefault(cfg, 'line_freq',   50);
cfg           = ft_setdefault(cfg, 'baseline',    [-Inf 0]);

if ~isfield(cfg,'fastica'),  cfg.fastica  = struct('approach','symm','g','gauss'); end
if ~isfield(cfg,'line'),     cfg.line     = struct();  end
if ~isfield(cfg,'muscle'),   cfg.muscle   = struct();  end
if ~isfield(cfg,'blink'),    cfg.blink    = struct();  end
if ~isfield(cfg,'cardiac'),  cfg.cardiac  = struct();  end
if ~isfield(cfg,'topo'),     cfg.topo     = struct();  end
if ~isfield(cfg,'decay'),    cfg.decay    = struct();  end
if ~isfield(cfg,'recharge'), cfg.recharge = struct();  end

% Thresholds marked (robust-z) are in median/MAD units; see local_rz().
cfg.line     = ft_setdefault(cfg.line,     'bw', 2);
cfg.line     = ft_setdefault(cfg.line,     'z',  5);   % robust-z on mains prominence
cfg.muscle   = ft_setdefault(cfg.muscle,   'band', [20 90]);
cfg.muscle   = ft_setdefault(cfg.muscle,   'z',  5);   % robust-z on HF power fraction
cfg.muscle   = ft_setdefault(cfg.muscle,   'hf_min', 0.5); % absolute HF-fraction floor
cfg.blink    = ft_setdefault(cfg.blink,    'kurt', 4);
cfg.blink    = ft_setdefault(cfg.blink,    'frontal_ratio', 0.5);
cfg.cardiac  = ft_setdefault(cfg.cardiac,  'lag', [0.5 1.2]); % periodicity lag window (s)
cfg.cardiac  = ft_setdefault(cfg.cardiac,  'z', 5);    % robust-z on periodicity
cfg.cardiac  = ft_setdefault(cfg.cardiac,  'minac', 0.15); % min autocorr peak
cfg.cardiac  = ft_setdefault(cfg.cardiac,  'kurt', 3.5); % QRS spikiness (excludes sinusoids)
cfg.topo     = ft_setdefault(cfg.topo,     'z', 5);    % robust-z on spatial focality
cfg.decay    = ft_setdefault(cfg.decay,    'window', [0.010 0.050]);
cfg.decay    = ft_setdefault(cfg.decay,    'k', 2);
cfg.recharge = ft_setdefault(cfg.recharge, 'window', [0.100 0.500]);
cfg.recharge = ft_setdefault(cfg.recharge, 'k', 4);
cfg          = ft_setdefault(cfg, 'frontal_labels', ...
    {'Fp1','Fp2','Fpz','AF7','AF8','AF3','AF4','F7','F8','F5','F6',...
     'F3','F4','Fz','F1','F2'});

% ------------------------------------------------------------------- ICA
if ~isempty(cfg.comp)
    comp = cfg.comp;
else
    cfg_ica              = [];
    cfg_ica.demean       = cfg.demean;
    cfg_ica.method       = cfg.method;
    if ~isempty(cfg.numcomponent), cfg_ica.numcomponent = cfg.numcomponent; end
    if strcmpi(cfg.method,'fastica'), cfg_ica.fastica = cfg.fastica; end
    comp                 = ft_componentanalysis(cfg_ica, data);
end

ncomp   = numel(comp.label);
nchan   = numel(comp.topolabel);

% --------------------------------------------------- spectral analysis
fs        = 1 / mean(diff(comp.time{1}));
nyq       = fs/2;
cfg_freq          = [];
cfg_freq.method   = 'mtmfft';
cfg_freq.taper    = 'hanning';
cfg_freq.output   = 'pow';
cfg_freq.foilim   = [1 min(100, floor(nyq)-1)];
freq              = ft_freqanalysis(cfg_freq, comp);
P                 = freq.powspctrm;      % [ncomp x nfreq]
f                 = freq.freq;

% --------------------------------------------------- per-component metrics
frontal_idx = ismember(comp.topolabel, cfg.frontal_labels);
if ~any(frontal_idx)
    warning('tefar_core:noFrontal', ...
        ['None of cfg.frontal_labels match comp.topolabel; blink/eye ' ...
         'detection will be disabled. Check your montage naming.']);
end

m.line_prom     = zeros(1,ncomp);   % mains fundamental / sideband power (peak prominence)
m.hf_frac       = zeros(1,ncomp);   % high-frequency power fraction (EMG)
m.periodicity   = zeros(1,ncomp);   % autocorrelation peak in cardiac lag window
m.topo_peak     = zeros(1,ncomp);   % raw peak |topo| (diagnostic)
m.topo_focal    = zeros(1,ncomp);   % peak / sum(|topo|): spatial focality
m.frontal_ratio = zeros(1,ncomp);   % frontal |topo| energy / total
m.frontal_asym  = zeros(1,ncomp);   % L-R anti-symmetry (lateral eye moves)
m.kurt          = zeros(1,ncomp);
m.decay_amp     = zeros(1,ncomp);
m.recharge_amp  = zeros(1,ncomp);
m.baseline_sd   = zeros(1,ncomp);

% spectral index sets
f0        = cfg.line_freq; bw = cfg.line.bw;
lineband  = local_line_band(f0, bw, f);               % mains + harmonics (to exclude from HF)
lfund     = abs(f - f0) <= bw;                        % mains fundamental
lbg       = (f>=f0-6 & f<=f0-2) | (f>=f0+2 & f<=f0+6);% mains sidebands (local background)
hiband    = f >= cfg.muscle.band(1) & f <= min(cfg.muscle.band(2),nyq-1) & ~lineband;
totband   = f >= 1;

% frontal L/R split for eye-movement asymmetry
[isL,isR] = local_lr_frontal(comp.topolabel, frontal_idx);

for c = 1:ncomp
    totpow            = mean(P(c,totband)) + eps;
    m.line_prom(c)    = local_bandpow(P(c,:), f, lfund) / (local_bandpow(P(c,:), f, lbg) + eps);
    m.hf_frac(c)      = local_bandpow(P(c,:), f, hiband) / totpow;
    m.periodicity(c)  = local_periodicity(comp, c, fs, cfg.cardiac.lag);

    topo             = comp.topo(:,c);
    m.topo_peak(c)   = max(abs(topo));
    m.topo_focal(c)  = max(abs(topo)) / (sum(abs(topo)) + eps);
    if any(frontal_idx)
        m.frontal_ratio(c)      = sum(abs(topo(frontal_idx))) / (sum(abs(topo)) + eps);
        if any(isL) && any(isR)
            m.frontal_asym(c)   = -(mean(topo(isL)) * mean(topo(isR))) / (m.topo_peak(c)^2 + eps);
        end
    end

    % concatenate this component's time course across all trials
    tc          = local_concat(comp, c);
    m.kurt(c)   = local_kurtosis(tc);

    % temporal windows (generic, per-trial)
    m.baseline_sd(c)  = local_window_std(comp, c, cfg.baseline);
    m.decay_amp(c)    = local_window_rms(comp, c, cfg.decay.window);
    m.recharge_amp(c) = local_window_rms(comp, c, cfg.recharge.window);
end

% data-driven robust-z (median / MAD) — resistant to the artefact's own outlier
zline    = local_rz(m.line_prom);
zemg     = local_rz(m.hf_frac);
ztopo    = local_rz(m.topo_focal);
zcardiac = local_rz(m.periodicity);

% ------------------------------------------------------------- detectors
A = struct('line',[],'muscle',[],'muscle_topo',[],'blink',[],...
           'eyemove',[],'cardiac',[],'decay',[],'recharge',[]);
score = zeros(1,ncomp);
det   = cfg.detectors;

for c = 1:ncomp
    hit = false;

    if want(det,'line') && zline(c) > cfg.line.z
        A.line(end+1) = c; hit = true;
    end
    if want(det,'muscle') && zemg(c) > cfg.muscle.z ...
            && m.hf_frac(c) > cfg.muscle.hf_min
        A.muscle(end+1) = c; hit = true;
    end
    if want(det,'muscle_topo') && ztopo(c) > cfg.topo.z
        A.muscle_topo(end+1) = c; hit = true;
    end
    if want(det,'blink') && any(frontal_idx) ...
            && m.kurt(c) > cfg.blink.kurt ...
            && m.frontal_ratio(c) > cfg.blink.frontal_ratio
        A.blink(end+1) = c; hit = true;
    end
    if want(det,'eyemove') && any(frontal_idx) ...
            && m.frontal_asym(c) > 0.15 ...
            && m.frontal_ratio(c) > cfg.blink.frontal_ratio
        A.eyemove(end+1) = c; hit = true;
    end
    if want(det,'cardiac') && zcardiac(c) > cfg.cardiac.z ...
            && m.periodicity(c) > cfg.cardiac.minac ...
            && m.kurt(c) > cfg.cardiac.kurt
        A.cardiac(end+1) = c; hit = true;
    end
    if want(det,'decay') && m.baseline_sd(c) > 0 ...
            && m.decay_amp(c) > cfg.decay.k * m.baseline_sd(c)
        A.decay(end+1) = c; hit = true;
    end
    if want(det,'recharge') && m.baseline_sd(c) > 0 ...
            && m.recharge_amp(c) > cfg.recharge.k * m.baseline_sd(c)
        A.recharge(end+1) = c; hit = true;
    end

    if hit, score(c) = score(c) + ...
        ismember(c,A.line)+ismember(c,A.muscle)+ismember(c,A.muscle_topo)+...
        ismember(c,A.blink)+ismember(c,A.eyemove)+ismember(c,A.cardiac)+...
        ismember(c,A.decay)+ismember(c,A.recharge);
    end
end

% ------------------------------------------------------- final rejection
switch lower(cfg.reject_rule)
    case 'union'
        reject = unique([A.line A.muscle A.muscle_topo A.blink ...
                         A.eyemove A.cardiac A.decay A.recharge]);
    case 'score'
        reject = find(score >= cfg.reject_threshold);
    otherwise
        error('tefar_core:rejectRule','Unknown cfg.reject_rule "%s".',cfg.reject_rule);
end

% ------------------------------------------------------------- assemble out
artifacts          = A;
artifacts.reject   = reject;
artifacts.score    = score;
artifacts.metrics  = m;
artifacts.cfg      = cfg;

if cfg.verbose
    local_report(A, reject, score, cfg);
end
end 

% ======================================================================
% Subfunctions (all toolbox-free)
% ======================================================================
function s = ft_setdefault(s, field, val)
if ~isfield(s, field) || isempty(s.(field)), s.(field) = val; end
end

function tf = want(list, name)
tf = any(strcmp(list, name));
end

function band = local_line_band(f0, bw, f)
% logical mask of mains fundamental + harmonics up to available spectrum
fmax  = max(f);
harm  = f0:f0:fmax;
band  = false(size(f));
for h = harm
    band = band | (f >= h-bw & f <= h+bw);
end
end

function p = local_bandpow(Prow, f, band)
% BAND may be a logical mask over f, or a [lo hi] frequency range.
if islogical(band)
    sel = band;
else
    sel = f >= band(1) & f <= band(2);
end
if ~any(sel), p = 0; else, p = mean(Prow(sel)); end
end

function rzx = local_rz(x)
% robust z-score using median and MAD (median absolute deviation)
med = median(x);
madv = median(abs(x - med));
rzx = (x - med) ./ (1.4826*madv + eps);
end

function val = local_periodicity(comp, c, fs, lagwin)
% mean normalised autocorrelation peak in a physiological lag window
peaks = [];
l0 = round(lagwin(1)*fs);
for k = 1:numel(comp.trial)
    s = comp.trial{k}(c,:); s = s - mean(s); n = numel(s);
    if n < 4, continue; end
    ac = xcorr_local(s);           % 0..n-1 lags, normalised
    l1 = min(round(lagwin(2)*fs), n-1);
    if l1 > l0 && l0 >= 1
        peaks(end+1) = max(ac(l0+1:l1+1)); %#ok<AGROW>  (+1: lag 0 at index 1)
    end
end
if isempty(peaks), val = 0; else, val = mean(peaks); end
end

function ac = xcorr_local(s)
% one-sided normalised autocorrelation via FFT (no Signal toolbox)
n  = numel(s);
nf = 2^nextpow2(2*n-1);
S  = fft(s, nf);
r  = real(ifft(S .* conj(S)));
ac = r(1:n) / (r(1) + eps);
end

function tc = local_concat(comp, c)
tc = [];
for k = 1:numel(comp.trial)
    tc = [tc comp.trial{k}(c,:)]; %#ok<AGROW>
end
end

function v = local_window_rms(comp, c, win)
acc = [];
for k = 1:numel(comp.trial)
    t   = comp.time{k};
    sel = t >= win(1) & t <= win(2);
    if any(sel)
        seg = comp.trial{k}(c, sel);
        acc(end+1) = sqrt(mean(seg.^2)); %#ok<AGROW>   RMS (unsigned)
    end
end
if isempty(acc), v = 0; else, v = mean(acc); end
end

function s = local_window_std(comp, c, win)
acc = [];
for k = 1:numel(comp.trial)
    t   = comp.time{k};
    sel = t >= win(1) & t <= win(2);
    if any(sel), acc = [acc comp.trial{k}(c, sel)]; end %#ok<AGROW>
end
if isempty(acc), s = 0; else, s = std(acc); end
end

function k = local_kurtosis(x)
% excess-free (Fisher=false) kurtosis, no Statistics toolbox
x  = x(~isnan(x));
mu = mean(x);
s2 = mean((x-mu).^2);
if s2 == 0, k = 0; return; end
k  = mean((x-mu).^4) / (s2^2);
end

function [isL,isR] = local_lr_frontal(labels, frontal_idx)
% classify frontal channels as left/right by trailing odd/even digit
isL = false(size(labels)); isR = false(size(labels));
idx = find(frontal_idx);
for i = idx(:)'
    d = regexp(labels{i}, '\d+$', 'match', 'once');
    if isempty(d), continue; end
    if mod(str2double(d),2)==1, isL(i) = true; else, isR(i) = true; end
end
end

function local_report(A, reject, score, cfg)
fprintf('\n==== TEFAR core: flagged components (profile detectors: %s) ====\n', ...
        strjoin(cfg.detectors, ', '));
flds = {'line','muscle','muscle_topo','blink','eyemove','cardiac','decay','recharge'};
for i = 1:numel(flds)
    if want(cfg.detectors, flds{i})
        fprintf('  %-12s : %s\n', flds{i}, mat2str(A.(flds{i})));
    end
end
fprintf('  ---------------------------------------------\n');
fprintf('  suggested reject (%s): %s\n', cfg.reject_rule, mat2str(reject));
hi = find(score >= max(1,cfg.reject_threshold));
if ~isempty(hi)
    fprintf('  high-score comps (>=%d): %s\n', cfg.reject_threshold, mat2str(hi));
end
fprintf('=================================================================\n\n');
end
