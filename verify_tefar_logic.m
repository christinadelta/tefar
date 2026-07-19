function results = verify_tefar_logic(profile, seed)
% VERIFY_TEFAR_LOGIC  Toolbox-free reproduction of the TEFAR detector math.
%
%   results = verify_tefar_logic(profile, seed)   % profile = 'tms' | 'eeg'
%   seed is optional (default 7); vary it to test robustness across draws.
%
% Independent re-implementation of the tefar_core detector logic WITHOUT
% FieldTrip (spectra via a local Hann-tapered periodogram) so the detector
% parameterisation can be validated against the ground-truth simulation in
% plain MATLAB or Octave. It is a cross-check / test harness, not part of
% the pipeline. It mirrors tefar_core exactly:
%   line     -> mains fundamental / sideband prominence   (robust-z)
%   muscle   -> high-frequency power fraction (excl mains) (robust-z + floor)
%   muscle_topo -> spatial focality peak/sum(|topo|)       (robust-z)
%   blink    -> kurtosis + frontal ratio
%   eyemove  -> fronto-lateral anti-symmetry + frontal ratio
%   cardiac  -> autocorrelation periodicity + spikiness    (robust-z)
%   decay/recharge -> windowed RMS vs baseline SD
%
% Prints an artefact-vs-brain confusion matrix and sensitivity/specificity/
% accuracy/F1, plus a per-detector breakdown.

if nargin < 1, profile = 'tms'; end
if nargin < 2 || isempty(seed), seed = 7; end

[comp, truth] = simulate_tefar_data(struct('profile',profile,'ntrials',40,'fs',1000,'seed',seed));
cfg = default_cfg(profile);
frontal = {'Fp1','Fp2','Fpz','AF7','AF8','AF3','AF4','F7','F8','F5','F6',...
           'F3','F4','Fz','F1','F2'};

ncomp = numel(comp.label);
fs    = comp.fsample; nyq = fs/2;
[P,f] = avg_periodogram(comp, [1 min(100,floor(nyq)-1)]);

frontal_idx = ismember(comp.topolabel, frontal);
[isL,isR]   = lr_frontal(comp.topolabel, frontal_idx);

f0=cfg.line_freq; bw=cfg.line.bw;
lineband = line_band(f0,bw,f);
lfund = abs(f-f0)<=bw;
lbg   = (f>=f0-6 & f<=f0-2) | (f>=f0+2 & f<=f0+6);
hi    = f>=cfg.muscle.band(1) & f<=min(cfg.muscle.band(2),nyq-1) & ~lineband;
tot   = f>=1;

m.line=zeros(1,ncomp); m.emg=zeros(1,ncomp); m.kurt=zeros(1,ncomp);
m.foc=zeros(1,ncomp);  m.fr=zeros(1,ncomp);  m.asym=zeros(1,ncomp);
m.per=zeros(1,ncomp);  m.decay=zeros(1,ncomp); m.rech=zeros(1,ncomp); m.bsd=zeros(1,ncomp);

for c = 1:ncomp
    totpow  = mean(P(c,tot))+eps;
    m.line(c)= (mean(P(c,lfund)))/(mean(P(c,lbg))+eps);
    m.emg(c) = (mean(P(c,hi)))/totpow;
    m.per(c) = periodicity(comp,c,fs,cfg.cardiac.lag);
    topo = comp.topo(:,c);
    m.foc(c) = max(abs(topo))/(sum(abs(topo))+eps);
    if any(frontal_idx)
        m.fr(c) = sum(abs(topo(frontal_idx)))/(sum(abs(topo))+eps);
        if any(isL)&&any(isR)
            m.asym(c) = -(mean(topo(isL))*mean(topo(isR)))/(max(abs(topo))^2+eps);
        end
    end
    m.kurt(c) = kurt_local(concat(comp,c));
    m.bsd(c)  = win_std(comp,c,cfg.baseline);
    m.decay(c)= win_rms(comp,c,cfg.decay.window);
    m.rech(c) = win_rms(comp,c,cfg.recharge.window);
end

zline=rz(m.line); zemg=rz(m.emg); ztopo=rz(m.foc); zcard=rz(m.per);

A = struct('line',[],'muscle',[],'muscle_topo',[],'blink',[],...
           'eyemove',[],'cardiac',[],'decay',[],'recharge',[]);
for c = 1:ncomp
    if want(cfg.detectors,'line')        && zline(c)>cfg.line.z,  A.line(end+1)=c; end
    if want(cfg.detectors,'muscle')      && zemg(c)>cfg.muscle.z && m.emg(c)>cfg.muscle.hf_min, A.muscle(end+1)=c; end
    if want(cfg.detectors,'muscle_topo') && ztopo(c)>cfg.topo.z,  A.muscle_topo(end+1)=c; end
    if want(cfg.detectors,'blink') && any(frontal_idx) && ...
            m.kurt(c)>cfg.blink.kurt && m.fr(c)>cfg.blink.frontal_ratio, A.blink(end+1)=c; end
    if want(cfg.detectors,'eyemove') && any(frontal_idx) && ...
            m.asym(c)>0.15 && m.fr(c)>cfg.blink.frontal_ratio, A.eyemove(end+1)=c; end
    if want(cfg.detectors,'cardiac') && zcard(c)>cfg.cardiac.z && ...
            m.per(c)>cfg.cardiac.minac && m.kurt(c)>cfg.cardiac.kurt, A.cardiac(end+1)=c; end
    if want(cfg.detectors,'decay') && m.bsd(c)>0 && m.decay(c)>cfg.decay.k*m.bsd(c), A.decay(end+1)=c; end
    if want(cfg.detectors,'recharge') && m.bsd(c)>0 && m.rech(c)>cfg.recharge.k*m.bsd(c), A.recharge(end+1)=c; end
end
reject = unique([A.line A.muscle A.muscle_topo A.blink A.eyemove A.cardiac A.decay A.recharge]);

pred=false(1,ncomp); pred(reject)=true; pos=truth.isartefact;
TP=sum(pred&pos); FP=sum(pred&~pos); FN=sum(~pred&pos); TN=sum(~pred&~pos);
sens=TP/max(TP+FN,1); spec=TN/max(TN+FP,1); acc=(TP+TN)/ncomp;
prec=TP/max(TP+FP,1); f1=2*prec*sens/max(prec+sens,eps);

fprintf('\n========== VERIFY TEFAR LOGIC (profile: %s) ==========\n',profile);
fprintf('components: %d (artefact %d, brain %d)\n',ncomp,sum(pos),sum(~pos));
fprintf('TP=%d FP=%d FN=%d TN=%d\n',TP,FP,FN,TN);
fprintf('sensitivity=%.3f specificity=%.3f accuracy=%.3f F1=%.3f\n',sens,spec,acc,f1);
fprintf('------------------------------------------------------\n');
flds={'line','muscle','muscle_topo','blink','eyemove','cardiac','decay','recharge'};
for i=1:numel(flds)
    if want(cfg.detectors,flds{i})
        idx=A.(flds{i});
        fprintf('  %-12s flagged %-12s true: [%s]\n',flds{i},mat2str(idx),strjoin(truth.label(idx),','));
    end
end
fprintf('  missed: %s   false alarms: %s\n', ...
    lbls(truth,find(pos&~pred)), lbls(truth,find(~pos&pred)));
fprintf('======================================================\n\n');

results=struct('sens',sens,'spec',spec,'acc',acc,'f1',f1,'A',A,'reject',reject,'metrics',m,'truth',truth);
end

% ============================ helpers ==================================
function cfg = default_cfg(profile)
cfg.line_freq=50; cfg.baseline=[-0.5 -0.01];
cfg.line.bw=2; cfg.line.z=5;
cfg.muscle.band=[20 90]; cfg.muscle.z=5; cfg.muscle.hf_min=0.5;
cfg.blink.kurt=4; cfg.blink.frontal_ratio=0.5;
cfg.cardiac.lag=[0.5 1.2]; cfg.cardiac.z=5; cfg.cardiac.minac=0.15; cfg.cardiac.kurt=3.5;
cfg.topo.z=5;
cfg.decay.window=[0.010 0.050]; cfg.decay.k=2;
cfg.recharge.window=[0.100 0.500]; cfg.recharge.k=4;
if strcmpi(profile,'tms')
    cfg.detectors={'line','muscle_topo','muscle','decay','recharge','blink','cardiac'};
else
    cfg.detectors={'line','blink','eyemove','muscle','muscle_topo','cardiac'};
end
end

function tf=want(list,name), tf=any(strcmp(list,name)); end

function band=line_band(f0,bw,f)
fmax=max(f); harm=f0:f0:fmax; band=false(size(f));
for h=harm, band=band|(f>=h-bw & f<=h+bw); end
end

function [P,f]=avg_periodogram(comp,foilim)
K=numel(comp.trial); nt=size(comp.trial{1},2); nc=size(comp.trial{1},1);
w=hann_local(nt); U=sum(w.^2); fs=comp.fsample;
f=(0:floor(nt/2))*(fs/nt); sel=f>=foilim(1)&f<=foilim(2); f=f(sel);
P=zeros(nc,sum(sel));
for k=1:K
    X=fft(bsxfun(@times,comp.trial{k},w'),[],2);
    Pk=(abs(X(:,1:floor(nt/2)+1)).^2)/(fs*U);
    P=P+Pk(:,sel);
end
P=P/K;
end

function w=hann_local(n), w=0.5-0.5*cos(2*pi*(0:n-1)'/(n-1)); end

function rzx=rz(x)
med=median(x); madv=median(abs(x-med)); rzx=(x-med)./(1.4826*madv+eps);
end

function val=periodicity(comp,c,fs,lagwin)
peaks=[]; l0=round(lagwin(1)*fs);
for k=1:numel(comp.trial)
    s=comp.trial{k}(c,:); s=s-mean(s); n=numel(s);
    nf=2^nextpow2(2*n-1); S=fft(s,nf); r=real(ifft(S.*conj(S)));
    ac=r(1:n)/(r(1)+eps); l1=min(round(lagwin(2)*fs),n-1);
    if l1>l0 && l0>=1, peaks(end+1)=max(ac(l0+1:l1+1)); end
end
if isempty(peaks), val=0; else val=mean(peaks); end
end

function tc=concat(comp,c)
tc=[]; for k=1:numel(comp.trial), tc=[tc comp.trial{k}(c,:)]; end
end
function v=win_rms(comp,c,win)
acc=[]; for k=1:numel(comp.trial)
    t=comp.time{k}; sel=t>=win(1)&t<=win(2);
    if any(sel), acc(end+1)=sqrt(mean(comp.trial{k}(c,sel).^2)); end
end
if isempty(acc), v=0; else v=mean(acc); end
end
function s=win_std(comp,c,win)
acc=[]; for k=1:numel(comp.trial)
    t=comp.time{k}; sel=t>=win(1)&t<=win(2);
    if any(sel), acc=[acc comp.trial{k}(c,sel)]; end
end
if isempty(acc), s=0; else s=std(acc); end
end
function k=kurt_local(x)
x=x(~isnan(x)); mu=mean(x); s2=mean((x-mu).^2);
if s2==0, k=0; else k=mean((x-mu).^4)/(s2^2); end
end
function [isL,isR]=lr_frontal(labels,frontal_idx)
isL=false(size(labels)); isR=false(size(labels));
for i=find(frontal_idx)'
    dd=regexp(labels{i},'\d+$','match','once');
    if isempty(dd), continue; end
    if mod(str2double(dd),2)==1, isL(i)=true; else isR(i)=true; end
end
end
function s=lbls(truth,idx)
if isempty(idx), s='[]'; else s=['[' strjoin(truth.label(idx),',') ']']; end
end