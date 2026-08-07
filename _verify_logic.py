"""Faithful Python port of simulate_tefar_data + verify_tefar_logic.
Purpose: confirm the TEFAR detector math/thresholds actually separate
artefact components from brain components on the ground-truth simulation.
Mirrors the MATLAB formulas exactly (bandpower=mean power in band,
z-score across components, same thresholds)."""
import numpy as np

CHAN = ['Fp1','Fp2','AF3','AF4','F7','F3','Fz','F4','F8','FC5','FC1',
        'FC2','FC6','T7','C3','Cz','C4','T8','CP5','CP1','CP2','CP6',
        'P7','P3','Pz','P4','P8','O1','Oz','O2']
FRONTAL = {'Fp1','Fp2','Fpz','AF7','AF8','AF3','AF4','F7','F8','F5','F6',
           'F3','F4','Fz','F1','F2'}

def pink(nt, rng, expo=0.5):
    X = np.fft.fft(rng.standard_normal(nt))
    f = np.arange(nt).astype(float); f[0]=1
    y = np.real(np.fft.ifft(X/(f**expo)))
    return y/np.std(y)

def highpass(x, fs, fc):
    nt=len(x); X=np.fft.fft(x); f=np.arange(nt)*(fs/nt)
    f=np.where(f>fs/2, fs-f, f); X[f<fc]=0
    return np.real(np.fft.ifft(X))

def topo_on(names, weights, rng):
    v=0.05*rng.standard_normal(len(CHAN))
    for n,w in zip(names,weights):
        if n in CHAN: v[CHAN.index(n)]=w
    return v

def simulate(profile='tms', ntrials=40, fs=1000, nbrain=12, line_freq=50, seed=7, noise=0.05):
    rng=np.random.default_rng(seed)
    if profile=='tms': t=np.arange(-1,1,1/fs)
    else:              t=np.arange(0,2,1/fs)
    nt=len(t)
    gens=[]; A=[]; labels=[]
    for _ in range(nbrain):
        fpk=3+11*rng.random()                 # theta/alpha/low-beta (3-14 Hz)
        A.append(0.3*rng.standard_normal(len(CHAN)))
        ph=2*np.pi*rng.random()
        gens.append(lambda fpk=fpk,ph=ph: _brain(nt,fs,fpk,ph,rng))
        labels.append('brain')
    # line
    A.append(0.6+0.2*rng.standard_normal(len(CHAN)))
    gens.append(lambda: 3.0*(np.sin(2*np.pi*line_freq*t)+0.3*np.sin(2*np.pi*2*line_freq*t)))
    labels.append('line')
    # blink
    A.append(topo_on(['Fp1','Fp2','AF3','AF4','Fz'],[1,1,.8,.8,.5],rng))
    gens.append(lambda: _blink(t,rng)); labels.append('blink')
    # muscle
    A.append(topo_on(['T7','FC5'],[1,.6],rng))
    gens.append(lambda: _emg(nt,fs,rng)); labels.append('muscle')
    # cardiac
    A.append(0.4+topo_on(['T7','T8','P7','P8'],[.5,.5,.4,.4],rng))
    gens.append(lambda: _cardiac(t)); labels.append('cardiac')
    if profile=='eeg':
        A.append(topo_on(['Fp1','F7','AF3','Fp2','F8','AF4'],[1,.8,.6,-1,-.8,-.6],rng))
        gens.append(lambda: _eyemove(t,rng)); labels.append('eyemove')
    if profile=='tms':
        A.append(topo_on(['C3','FC5','FC1'],[1,.7,.5],rng))
        gens.append(lambda: _decay(t,0.020)); labels.append('decay')
        A.append(topo_on(['C3','FC1'],[.8,.5],rng))
        gens.append(lambda: _recharge(t,[0.15,0.45])); labels.append('recharge')
    A=np.array(A).T  # nchan x nsource
    nsrc=A.shape[1]
    trials=[]
    for k in range(ntrials):
        s=np.zeros((nsrc,nt))
        for j,g in enumerate(gens):
            sj=g(); sj=sj+noise*np.std(sj)*rng.standard_normal(nt); s[j]=sj
        trials.append(s)
    comp=dict(topo=A, trial=trials, time=[t]*ntrials, fs=fs, label=list(range(nsrc)))
    truth=dict(label=labels, isart=np.array([l!='brain' for l in labels]))
    return comp, truth

def _brain(nt,fs,fpk,ph,rng):
    y=pink(nt,rng,expo=1.0)+0.6*np.sin(2*np.pi*fpk*np.arange(nt)/fs+ph); return y/np.std(y)
def _blink(t,rng):
    nt=len(t); y=np.zeros(nt); nb=max(1,round(2*(t[-1]-t[0]))); w=0.08
    for _ in range(nb):
        tc=t[0]+(t[-1]-t[0])*rng.random(); y+=6*np.exp(-((t-tc)**2)/(2*w*w))
    return y+0.05*pink(nt,rng)
def _emg(nt,fs,rng):
    y=highpass(rng.standard_normal(nt),fs,20); return 1.5*y/np.std(y)
def _cardiac(t):
    nt=len(t); y=np.zeros(nt); hr=1.1; tc=t[0]
    while tc<t[-1]: y+=3*np.exp(-((t-tc)**2)/(2*0.010**2)); tc+=1/hr
    return y
def _eyemove(t,rng):
    nt=len(t); y=np.zeros(nt); ns=max(1,round(1.5*(t[-1]-t[0])))
    edges=np.sort(t[0]+(t[-1]-t[0])*rng.random(ns)).tolist(); lvl=0
    for i in range(nt):
        if edges and t[i]>=edges[0]: lvl=3*(2*rng.random()-1); edges.pop(0)
        y[i]=lvl
    return y+0.05*pink(nt,rng)
def _decay(t,tau):
    y=np.zeros_like(t); sel=t>=0; y[sel]=20*np.exp(-t[sel]/tau); return y
def _recharge(t,win):
    c=np.mean(win); w=(win[1]-win[0])/2
    return 5*np.exp(-((t-c)**2)/(2*(w/2)**2))

# ---- detector math (mirrors tefar_core) ----
def avg_periodogram(comp, foilim):
    trials=comp['trial']; nt=trials[0].shape[1]; nc=trials[0].shape[0]; fs=comp['fs']
    w=0.5-0.5*np.cos(2*np.pi*np.arange(nt)/(nt-1)); U=np.sum(w**2)
    f=np.arange(nt//2+1)*(fs/nt); sel=(f>=foilim[0])&(f<=foilim[1])
    P=np.zeros((nc,sel.sum()))
    for x in trials:
        X=np.fft.fft(x*w,axis=1)
        Pk=(np.abs(X[:,:nt//2+1])**2)/(fs*U); P+=Pk[:,sel]
    return P/len(trials), f[sel]

def bandpow(Prow,f,b):
    sel=(f>=b[0])&(f<=b[1]); return Prow[sel].mean() if sel.any() else 0.0
def line_band(f0,bw,f):
    band=np.zeros(len(f),bool); h=f0
    while h<=f.max(): band|=(f>=h-bw)&(f<=h+bw); h+=f0
    return band
def kurt(x):
    x=x[~np.isnan(x)]; mu=x.mean(); s2=((x-mu)**2).mean()
    return 0 if s2==0 else ((x-mu)**4).mean()/s2**2
def z(x): return (x-x.mean())/(x.std()+1e-12)
def rz(x):
    med=np.median(x); mad=np.median(np.abs(x-med))
    return (x-med)/(1.4826*mad+1e-12)          # robust z (median/MAD)
def periodicity(comp,c,fs,lagwin=(0.5,1.2)):
    # mean normalised autocorrelation peak in a physiological lag window
    peaks=[]
    for x in comp['trial']:
        s=x[c]-np.mean(x[c]); n=len(s)
        ac=np.correlate(s,s,'full')[n-1:]; ac=ac/(ac[0]+1e-12)
        l0=int(lagwin[0]*fs); l1=min(int(lagwin[1]*fs),n-1)
        if l1>l0: peaks.append(ac[l0:l1].max())
    return float(np.mean(peaks)) if peaks else 0.0

def default_cfg(profile):
    # thresholds are in ROBUST-z (median/MAD) units unless noted
    c=dict(line_freq=50,baseline=(-0.5,-0.01),
           line=dict(bw=2,zt=5),                 # narrowband mains power fraction
           muscle=dict(band=(20,90),zt=5,hf_min=0.5),  # high-frequency power fraction
           blink=dict(kurt=4.0,fr=0.5),          # kurtosis(abs) + frontal ratio
           cardiac=dict(lag=(0.5,1.2),zt=5,minac=0.15),  # temporal periodicity
           topo=dict(zt=5),                       # spatial focality
           decay=dict(window=(0.010,0.050),k=2),
           recharge=dict(window=(0.100,0.500),k=4))
    c['det']=(['line','muscle_topo','muscle','decay','recharge','blink','cardiac'] if profile=='tms'
              else ['line','blink','eyemove','muscle','muscle_topo','cardiac'])
    return c

def win_rms(comp,c,win):
    acc=[]
    for x,t in zip(comp['trial'],comp['time']):
        sel=(t>=win[0])&(t<=win[1])
        if sel.any(): acc.append(np.sqrt((x[c,sel]**2).mean()))
    return np.mean(acc) if acc else 0.0
def win_std(comp,c,win):
    acc=[]
    for x,t in zip(comp['trial'],comp['time']):
        sel=(t>=win[0])&(t<=win[1])
        if sel.any(): acc.append(x[c,sel])
    return np.std(np.concatenate(acc)) if acc else 0.0

def verify(profile):
    comp,truth=simulate(profile)
    cfg=default_cfg(profile); nc=comp['topo'].shape[1]; fs=comp['fs']; nyq=fs/2
    P,f=avg_periodogram(comp,[1,min(100,int(nyq)-1)])
    fidx=np.array([ch in FRONTAL for ch in CHAN])
    def lr():
        isL=np.zeros(len(CHAN),bool); isR=np.zeros(len(CHAN),bool)
        import re
        for i,ch in enumerate(CHAN):
            if not fidx[i]: continue
            m=re.search(r'\d+$',ch)
            if not m: continue
            (isL if int(m.group())%2 else isR)[i]=True
        return isL,isR
    isL,isR=lr()
    lb=line_band(cfg['line_freq'],cfg['line']['bw'],f)
    f0=cfg['line_freq']; bw=cfg['line']['bw']
    lfund=(np.abs(f-f0)<=bw)                                  # mains fundamental
    lbg=(((f>=f0-6)&(f<=f0-2))|((f>=f0+2)&(f<=f0+6)))         # sidebands (background)
    hi=(f>=cfg['muscle']['band'][0])&(f<=min(cfg['muscle']['band'][1],nyq-1))&(~lb)  # HF, exclude mains
    tot=(f>=1)
    mL=np.zeros(nc);mE=np.zeros(nc);mK=np.zeros(nc);mT=np.zeros(nc)
    mF=np.zeros(nc);mAsym=np.zeros(nc);mC=np.zeros(nc);mD=np.zeros(nc);mR=np.zeros(nc);mB=np.zeros(nc)
    for c in range(nc):
        totpow=P[c,tot].mean()+1e-12
        mL[c]=(P[c,lfund].mean() if lfund.any() else 0)/(P[c,lbg].mean()+1e-12)  # mains PROMINENCE (peak/sidebands)
        mE[c]=(P[c,hi].mean() if hi.any() else 0)/totpow      # HF power FRACTION
        mC[c]=periodicity(comp,c,fs,cfg['cardiac']['lag'])    # temporal periodicity
        topo=comp['topo'][:,c]; mT[c]=np.abs(topo).max()/(np.abs(topo).sum()+1e-12)  # focality
        if fidx.any():
            mF[c]=np.abs(topo[fidx]).sum()/(np.abs(topo).sum()+1e-12)
            if isL.any() and isR.any():
                pk=np.abs(topo).max()
                mAsym[c]=-(topo[isL].mean()*topo[isR].mean())/(pk**2+1e-12)
        tc=np.concatenate([x[c] for x in comp['trial']]); mK[c]=kurt(tc)
        mB[c]=win_std(comp,c,cfg['baseline']); mD[c]=win_rms(comp,c,cfg['decay']['window'])
        mR[c]=win_rms(comp,c,cfg['recharge']['window'])
    zL,zE,zT,zC=rz(mL),rz(mE),rz(mT),rz(mC)   # ROBUST z (median/MAD)
    A={k:[] for k in ['line','muscle','muscle_topo','blink','eyemove','cardiac','decay','recharge']}
    det=cfg['det']
    for c in range(nc):
        if 'line' in det and zL[c]>cfg['line']['zt']: A['line'].append(c)
        if 'muscle' in det and zE[c]>cfg['muscle']['zt'] and mE[c]>cfg['muscle']['hf_min']: A['muscle'].append(c)
        if 'muscle_topo' in det and zT[c]>cfg['topo']['zt']: A['muscle_topo'].append(c)
        if 'blink' in det and fidx.any() and mK[c]>cfg['blink']['kurt'] and mF[c]>cfg['blink']['fr']: A['blink'].append(c)
        if 'eyemove' in det and fidx.any() and mAsym[c]>0.15 and mF[c]>cfg['blink']['fr']: A['eyemove'].append(c)
        if 'cardiac' in det and zC[c]>cfg['cardiac']['zt'] and mC[c]>cfg['cardiac']['minac'] and mK[c]>3.5: A['cardiac'].append(c)
        if 'decay' in det and mB[c]>0 and mD[c]>cfg['decay']['k']*mB[c]: A['decay'].append(c)
        if 'recharge' in det and mB[c]>0 and mR[c]>cfg['recharge']['k']*mB[c]: A['recharge'].append(c)
    reject=sorted(set(sum(A.values(),[])))
    pred=np.zeros(nc,bool); pred[reject]=True; pos=truth['isart']
    TP=int((pred&pos).sum());FP=int((pred&~pos).sum());FN=int((~pred&pos).sum());TN=int((~pred&~pos).sum())
    sens=TP/max(TP+FN,1);spec=TN/max(TN+FP,1);acc=(TP+TN)/nc;prec=TP/max(TP+FP,1)
    f1=2*prec*sens/max(prec+sens,1e-12)
    print(f"\n===== profile: {profile} =====")
    print(f"components={nc} artefact={int(pos.sum())} brain={int((~pos).sum())}")
    print(f"TP={TP} FP={FP} FN={FN} TN={TN}")
    print(f"sensitivity={sens:.3f} specificity={spec:.3f} accuracy={acc:.3f} F1={f1:.3f}")
    for k in det:
        idx=A[k]; print(f"  {k:12s} flagged {idx}  labels={[truth['label'][i] for i in idx]}")
    miss=[i for i in range(nc) if pos[i] and not pred[i]]
    fa=[i for i in range(nc) if not pos[i] and pred[i]]
    print(f"  MISSED: {[(i,truth['label'][i]) for i in miss]}")
    print(f"  FALSE ALARMS: {[(i,truth['label'][i]) for i in fa]}")
    return sens,spec,acc

if __name__=='__main__':
    for p in ('tms','eeg'):
        verify(p)
