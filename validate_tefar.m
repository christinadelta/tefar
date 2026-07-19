function results = validate_tefar(profile, mode)
% VALIDATE_TEFAR  Validate TEFAR detectors against ground-truth simulation.
%
%   results = validate_tefar(profile, mode)
%     profile : 'tms' (default) | 'eeg'
%     mode    : 'components' (default) | 'endtoend'
%
% Requires FieldTrip on the path (ft_defaults). Two modes:
%
%   'components'  Feed the simulated ground-truth sources straight into
%                 tefar_core via cfg.comp. This isolates DETECTOR accuracy
%                 from ICA separation quality -- the cleanest validation of
%                 the classification rules themselves.
%
%   'endtoend'    Mix the sources to channels (x = A*s), hand raw data to
%                 TEFAR_tms / TEFAR_eeg (which run real FastICA), then match
%                 recovered ICs back to the true sources by topography
%                 correlation before scoring. Tests the full pipeline.
%
% Reports an artefact-vs-brain confusion matrix and sensitivity /
% specificity / accuracy / F1, plus a per-detector breakdown.
%
% @christinadelta 2026

if nargin < 1 || isempty(profile), profile = 'tms'; end
if nargin < 2 || isempty(mode),    mode    = 'components'; end

scfg = struct('profile',profile,'ntrials',40,'fs',1000,'seed',7);
[comp, truth] = simulate_tefar_data(scfg);

switch lower(mode)
    case 'components'
        cfg = struct('comp', comp, 'verbose', true);
        if strcmpi(profile,'tms')
            [~,~,~,~,~,~,~,art] = TEFAR_tms([], [], cfg);
        else
            [~,art] = TEFAR_eeg([], cfg);
        end
        pred = false(1,numel(comp.label)); pred(art.reject) = true;
        A = art;

    case 'endtoend'
        data = comp_to_raw(comp);                 % x = A*s as FieldTrip raw
        if strcmpi(profile,'tms')
            [rc,~,~,~,~,~,~,A] = TEFAR_tms(data, []);
        else
            [rc,A] = TEFAR_eeg(data);
        end
        % match recovered ICs to true sources by |topo correlation|
        map = match_sources(rc.topo, comp.topo);  % map(recovered)=truthIdx
        pred_truth = false(1,numel(truth.label));
        pred_truth(map(A.reject)) = true;
        pred = pred_truth;

    otherwise
        error('Unknown mode "%s".', mode);
end

pos = truth.isartefact;
TP=sum(pred&pos); FP=sum(pred&~pos); FN=sum(~pred&pos); TN=sum(~pred&~pos);
sens=TP/max(TP+FN,1); spec=TN/max(TN+FP,1);
acc=(TP+TN)/numel(pos); prec=TP/max(TP+FP,1);
f1=2*prec*sens/max(prec+sens,eps);

fprintf('\n=========== validate_tefar (%s / %s) ===========\n',profile,mode);
fprintf('TP=%d FP=%d FN=%d TN=%d\n',TP,FP,FN,TN);
fprintf('sensitivity=%.3f specificity=%.3f accuracy=%.3f F1=%.3f\n',sens,spec,acc,f1);
fprintf('================================================\n\n');

results = struct('sens',sens,'spec',spec,'acc',acc,'f1',f1,...
                 'artifacts',A,'truth',truth);
end

% ---------------------------------------------------------------- helpers
function data = comp_to_raw(comp)
% reconstruct channel data x = A*s for each trial
data = struct();
data.label   = comp.topolabel;
data.fsample = comp.fsample;
data.trial   = cell(1,numel(comp.trial));
data.time    = comp.time;
for k = 1:numel(comp.trial)
    data.trial{k} = comp.topo * comp.trial{k};
end
data.trialinfo = ones(numel(comp.trial),1);
end

function map = match_sources(recovered_topo, true_topo)
% for each recovered IC, index of best-correlated true source
nr = size(recovered_topo,2);
map = zeros(1,nr);
for i = 1:nr
    r = zeros(1,size(true_topo,2));
    for j = 1:size(true_topo,2)
        r(j) = abs(corr_local(recovered_topo(:,i), true_topo(:,j)));
    end
    [~,map(i)] = max(r);
end
end

function r = corr_local(a,b)
a=a-mean(a); b=b-mean(b);
r = (a'*b)/(sqrt(a'*a)*sqrt(b'*b)+eps);
end
