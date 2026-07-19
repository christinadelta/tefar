function [comp, truth] = simulate_tefar_data(cfg)
% SIMULATE_TEFAR_DATA  Ground-truth simulation for validating TEFAR detectors.
%
%   [comp, truth] = simulate_tefar_data(cfg)
%
% Builds data directly under the ICA generative model:  x = A * s.
% Each latent source s has a KNOWN artefact label and a KNOWN topography A,
% so the returned structure can be fed straight into tefar_core via
% cfg.comp, and every flag can be scored against ground truth.
%
% This function uses base MATLAB only (no FieldTrip, no toolboxes) so it can
% be run and unit-tested anywhere, including Octave.
%
% cfg fields (all optional)
%   .profile   'tms' (adds decay+recharge) | 'eeg'            ['tms']
%   .ntrials   number of trials                                   [40]
%   .fs        sampling rate (Hz)                                [1000]
%   .nbrain    number of brain sources                            [12]
%   .line_freq mains frequency (Hz)                                [50]
%   .seed      RNG seed                                             [1]
%   .noise     per-source additive noise SD (fraction of source)  [0.05]
%
% OUTPUT
%   comp   struct with fields topo, topolabel, label, trial, time, unmixing,
%          topodimord, unmixingdimord — compatible with tefar_core(cfg.comp).
%   truth  struct with:
%          .label   {1 x nsource} artefact class per component
%                   (one of 'brain','line','blink','eyemove','muscle',
%                    'cardiac','decay','recharge')
%          .isartefact  logical [1 x nsource]

if nargin < 1, cfg = struct(); end
cfg = d(cfg,'profile','tms');
cfg = d(cfg,'ntrials',40);
cfg = d(cfg,'fs',1000);
cfg = d(cfg,'nbrain',12);
cfg = d(cfg,'line_freq',50);
cfg = d(cfg,'seed',1);
cfg = d(cfg,'noise',0.05);

rng_local(cfg.seed);

% ---- montage (30ch subset of 10-20) -----------------------------------
chan = {'Fp1','Fp2','AF3','AF4','F7','F3','Fz','F4','F8','FC5','FC1',...
        'FC2','FC6','T7','C3','Cz','C4','T8','CP5','CP1','CP2','CP6',...
        'P7','P3','Pz','P4','P8','O1','Oz','O2'};
nchan = numel(chan);

% ---- time base ---------------------------------------------------------
if strcmpi(cfg.profile,'tms')
    t = (-1 : 1/cfg.fs : 1 - 1/cfg.fs);      % event-locked at 0
else
    t = (0 : 1/cfg.fs : 2 - 1/cfg.fs);       % 2 s epochs
end
nt = numel(t);

% ---- assemble sources --------------------------------------------------
S      = {};   % per-source time course generator handles produce [1 x nt]
A      = [];   % topographies [nchan x nsource]
labels = {};

% brain sources: 1/f + an oscillation, smooth-ish random topography
for b = 1:cfg.nbrain
    fpk  = 3 + 11*rand;                        % theta/alpha/low-beta (3-14 Hz)
    topo = 0.3*randn(nchan,1);                % distributed, modest
    S{end+1}    = @() brain_source(nt,cfg.fs,fpk);
    A(:,end+1)  = topo;
    labels{end+1} = 'brain';
end

% line: mains sinusoid (+2nd harmonic), broad topography
A(:,end+1)  = 0.6 + 0.2*randn(nchan,1);
S{end+1}    = @() 3.0*(sin(2*pi*cfg.line_freq*t) + 0.3*sin(2*pi*2*cfg.line_freq*t));
labels{end+1} = 'line';

% blink: sparse large low-frequency deflections, frontal same-sign topo
A(:,end+1)  = topo_on(chan,{'Fp1','Fp2','AF3','AF4','Fz'},[1 1 .8 .8 .5],nchan);
S{end+1}    = @() blink_source(t);
labels{end+1} = 'blink';

% muscle/EMG: broadband high-frequency, focal at a temporal edge
A(:,end+1)  = topo_on(chan,{'T7','FC5'},[1 .6],nchan);
S{end+1}    = @() emg_source(nt,cfg.fs);
labels{end+1} = 'muscle';

% cardiac: quasi-periodic QRS-like spikes ~1.1 Hz, broad posterior-ish topo
A(:,end+1)  = 0.4 + topo_on(chan,{'T7','T8','P7','P8'},[.5 .5 .4 .4],nchan);
S{end+1}    = @() cardiac_source(t);
labels{end+1} = 'cardiac';

if strcmpi(cfg.profile,'eeg')
    % lateral eye movement: step-like, L(+)/R(-) antisymmetric frontal topo
    A(:,end+1)  = topo_on(chan,{'Fp1','F7','AF3','Fp2','F8','AF4'},...
                          [1 .8 .6 -1 -.8 -.6],nchan);
    S{end+1}    = @() eyemove_source(t);
    labels{end+1} = 'eyemove';
end

if strcmpi(cfg.profile,'tms')
    % decay: large exponential right after t=0, focal near "coil" (C3/FC5)
    A(:,end+1)  = topo_on(chan,{'C3','FC5','FC1'},[1 .7 .5],nchan);
    S{end+1}    = @() decay_source(t,0.020);
    labels{end+1} = 'decay';

    % recharge: slower bump in a later window, focal
    A(:,end+1)  = topo_on(chan,{'C3','FC1'},[.8 .5],nchan);
    S{end+1}    = @() recharge_source(t,[0.15 0.45]);
    labels{end+1} = 'recharge';
end

nsource = numel(S);

% ---- build per-trial component activations ----------------------------
trial = cell(1,cfg.ntrials);
time  = cell(1,cfg.ntrials);
for k = 1:cfg.ntrials
    s = zeros(nsource, nt);
    for j = 1:nsource
        sj = S{j}();
        sj = sj + cfg.noise * std(sj) * randn(1,nt);   % measurement jitter
        s(j,:) = sj;
    end
    trial{k} = s;
    time{k}  = t;
end

% ---- pack as a FieldTrip component structure --------------------------
comp = struct();
comp.topo         = A;
comp.topolabel    = chan(:);
comp.unmixing     = pinv(A);
comp.label        = arrayfun(@(i) sprintf('runica%03d',i), 1:nsource, ...
                             'UniformOutput', false)';
comp.trial        = trial;
comp.time         = time;
comp.fsample      = cfg.fs;
comp.topodimord   = 'chan_comp';
comp.unmixingdimord = 'comp_chan';
comp.cfg          = struct('simulated', true);

truth = struct();
truth.label       = labels;
truth.isartefact  = ~strcmp(labels,'brain');
end % ===================================================================

% ---------------------------- source generators -----------------------
function y = brain_source(nt,fs,fpk)
y = pinknoise(nt, 1.0);                        % steeper 1/f (cortical-like)
% add a narrowband oscillation around fpk
osc = sin(2*pi*fpk*(0:nt-1)/fs + 2*pi*rand);
y = y + 0.6*osc;
y = y / std(y);
end

function y = blink_source(t)
nt = numel(t); y = zeros(1,nt);
nb = max(1, round(2*(t(end)-t(1))));          % ~2 blinks/s of epoch span
w  = 0.08;                                    % ~80 ms blink width
for i = 1:nb
    tc = t(1) + (t(end)-t(1))*rand;
    y  = y + 6*exp(-((t-tc).^2)/(2*w^2));
end
y = y + 0.05*pinknoise(nt);
end

function y = emg_source(nt,fs)
y = randn(1,nt);
y = highpass_fft(y, fs, 20);                  % broadband > 20 Hz
y = 1.5 * y / std(y);
end

function y = cardiac_source(t)
nt = numel(t); y = zeros(1,nt);
hr = 1.1;                                      % ~66 bpm
tc = t(1);
while tc < t(end)
    y = y + 3*exp(-((t-tc).^2)/(2*0.010^2));   % sharp QRS
    tc = tc + 1/hr;
end
end

function y = eyemove_source(t)
nt = numel(t); y = zeros(1,nt);
ns = max(1, round(1.5*(t(end)-t(1))));
lvl = 0;
edges = sort(t(1) + (t(end)-t(1))*rand(1,ns));
for i = 1:nt
    if ~isempty(edges) && t(i) >= edges(1)
        lvl = 3*(2*rand-1); edges(1) = [];
    end
    y(i) = lvl;
end
y = y + 0.05*pinknoise(nt);
end

function y = decay_source(t,tau)
y = zeros(size(t));
sel = t >= 0;
y(sel) = 20*exp(-t(sel)/tau);                  % huge early transient
end

function y = recharge_source(t,win)
y = zeros(size(t));
c  = mean(win); w = (win(2)-win(1))/2;
y = 5*exp(-((t-c).^2)/(2*(w/2)^2));
end

% ------------------------------- helpers -------------------------------
function v = topo_on(chan, names, weights, nchan)
v = 0.05*randn(nchan,1);
for i = 1:numel(names)
    idx = find(strcmp(chan, names{i}));
    if ~isempty(idx), v(idx) = weights(i); end
end
end

function y = pinknoise(nt, expo)
% 1/f^expo amplitude-spectrum noise, base MATLAB/Octave. expo=0.5 -> pink.
if nargin < 2, expo = 0.5; end
X = fft(randn(1,nt));
f = (0:nt-1); f(1) = 1;
X = X ./ (f.^expo);
y = real(ifft(X));
y = y / std(y);
end

function y = highpass_fft(x, fs, fc)
nt = numel(x);
X  = fft(x);
f  = (0:nt-1)*(fs/nt);
f(f>fs/2) = fs - f(f>fs/2);
X(f < fc) = 0;
y = real(ifft(X));
end

function s = d(s, field, val)
if ~isfield(s, field) || isempty(s.(field)), s.(field) = val; end
end

function rng_local(seed)
% Octave/MATLAB-safe RNG seeding
try
    rng(seed);
catch
    rand('state', seed); randn('state', seed); %#ok<RAND>
end
end
