function [comp, artifacts] = TEFAR_eeg(data, cfg)
% TEFAR_EEG  Artefact classification for continuous / epoched EEG (no TMS).
%
%   [comp, artifacts] = TEFAR_eeg(data, cfg)
%
% Sister function to TEFAR_tms for ordinary (resting / task) EEG. Same core
% engine, EEG-appropriate detector set. Detects:
%
%   line         mains (50/60 Hz) + harmonics
%   blink        eye blinks (frontal topography + high kurtosis + low-freq)
%   eyemove      lateral/saccadic eye movement (fronto-lateral anti-symmetry)
%   muscle       EMG (broadband high-frequency power)
%   muscle_topo  focal/edge topography
%   cardiac      cardiac rhythm (quasi-periodic ~1 Hz band power)
%
% INPUT
%   data   FieldTrip data (continuous or epoched; ft_preprocessing output)
%   cfg    optional overrides (see tefar_core for the full list). Common:
%            cfg.line_freq   (default 50; set 60 for US mains)
%            cfg.detectors   (override the default set below)
%            cfg.reject_rule 'union' (default) | 'score'
%
% OUTPUT
%   comp        FieldTrip component structure
%   artifacts   full core result: per-detector index lists, .reject, .score,
%               .metrics, .cfg
%
% @christinadelta 2026

if nargin < 2 || isempty(cfg), cfg = struct(); end

cfg = setdef(cfg, 'detectors', ...
        {'line','blink','eyemove','muscle','muscle_topo','cardiac'});
cfg = setdef(cfg, 'line_freq',   50);
cfg = setdef(cfg, 'reject_rule', 'union');

% For non-epoched EEG there is no event-locked baseline; the temporal
% decay/recharge detectors are off by default (not in the detector list),
% so baseline is irrelevant here. Spectral/topographic/statistical
% detectors need no time-zero.

[comp, artifacts] = tefar_core(cfg, data);
end

function s = setdef(s, field, val)
if ~isfield(s, field) || isempty(s.(field)), s.(field) = val; end
end
