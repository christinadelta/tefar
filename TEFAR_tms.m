function [comp_tms, line_components, muscle_components, decay_components, ...
          additional_muscle_comp, recharge_components, eye_blinks_components, ...
          artifacts] = TEFAR_tms(data_filtered, trl, cfg)
% TEFAR_TMS  TMS-EEG artefact classification (TEFAR, TMS profile).
%
%   [comp_tms, line_components, muscle_components, decay_components, ...
%    additional_muscle_comp, recharge_components, eye_blinks_components, ...
%    artifacts] = TEFAR_tms(data_filtered, trl, cfg)
%
% Backward-compatible drop-in replacement for TEFAR_v3: the first seven
% outputs match the original signature exactly, so existing pipelines keep
% working. The eighth output (artifacts) exposes the full core result
% (scores, metrics, populated cfg) for anyone who wants it.
%
% This is a thin wrapper: it sets TMS-appropriate defaults and calls
% tefar_core. Override any default by passing a cfg struct (third arg).
%
% Output mapping (old name -> core detector)
%   line_components        <- line        (mains + harmonics power)
%   muscle_components      <- muscle_topo (focal topography near coil)
%   additional_muscle_comp <- muscle      (high-frequency EMG power)
%   decay_components       <- decay       (early post-pulse amplitude)
%   recharge_components    <- recharge    (later recharge-window amplitude)
%   eye_blinks_components  <- blink       (frontal + high kurtosis)
%
% trl is accepted for signature compatibility; the core is trial-structure
% agnostic and does not require it.
%
% @christinadelta - 2026

if nargin < 3 || isempty(cfg), cfg = struct(); end

% TMS-profile defaults (only set what the user did not) 
cfg = setdef(cfg, 'detectors', ...
        {'line','muscle_topo','muscle','decay','recharge','blink','cardiac'});
cfg = setdef(cfg, 'line_freq',   50);
cfg = setdef(cfg, 'reject_rule', 'union');
cfg = setdef(cfg, 'baseline',    [-0.5 -0.01]);

if ~isfield(cfg,'decay'),    cfg.decay    = struct(); end
if ~isfield(cfg,'recharge'), cfg.recharge = struct(); end
cfg.decay    = setdef(cfg.decay,    'window', [0.010 0.050]);
cfg.decay    = setdef(cfg.decay,    'k', 2);
cfg.recharge = setdef(cfg.recharge, 'window', [0.100 0.500]);
cfg.recharge = setdef(cfg.recharge, 'k', 4);

% run core 
[comp_tms, artifacts] = tefar_core(cfg, data_filtered);

% map to legacy outputs 
line_components        = artifacts.line;
muscle_components      = artifacts.muscle_topo;
additional_muscle_comp = artifacts.muscle;
decay_components       = artifacts.decay;
recharge_components    = artifacts.recharge;
eye_blinks_components  = artifacts.blink;
end

function s = setdef(s, field, val)
if ~isfield(s, field) || isempty(s.(field)), s.(field) = val; end
end
