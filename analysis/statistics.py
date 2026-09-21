"""Explicit numerical summaries; block uncertainty is not a convergence test."""
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from .geometry import PAIRS

# Archived literature-direction convention: +1 means distance increases on activation.
# These signs describe the reference expectation, not fitted mutant outcomes.
ORIENTATION=dict(zip(PAIRS,[1,-1,1,1,-1,1,1,-1,-1,1,-1,1,-1,-1,-1,-1,1,1,-1,-1,-1]))

def classify(value,threshold):
    return int(value>=threshold)-int(value<=-threshold)

def block_size(times,block_ns):
    times=np.asarray(times,dtype=float)
    if len(times)<2 or np.any(np.diff(times)<=0):raise ValueError('Increasing times required')
    dt=np.diff(times)
    if not np.allclose(dt,dt[0],atol=1e-5):raise ValueError('Uniform frame spacing required')
    n=round(block_ns/dt[0])
    if n<1 or not np.isclose(n*dt[0],block_ns,atol=1e-5):raise ValueError('Invalid block duration')
    return n

def median_block_ci(active,inactive,active_size,inactive_size,rng,resamples=10000,factor=1):
    """CI for DIFFERENCE OF MEAN BLOCK MEDIANS, not the pooled median difference.

    Follows archived convention: retain a last partial block only if at least
    half a block. States are resampled independently. At least two blocks each.
    """
    blocks=[]
    for values,size in [(active,active_size),(inactive,inactive_size)]:
        if size<1:raise ValueError('Positive block size required')
        v=np.asarray(values,dtype=float)
        b=np.array([np.median(v[i:i+size]) for i in range(0,len(v),size)
                    if len(v[i:i+size])>=max(1,size/2)])
        if len(b)<2 or not np.isfinite(b).all():raise ValueError('At least two finite blocks per state required')
        blocks.append(b)
    a,b=blocks
    differences=factor*(rng.choice(a,(resamples,len(a))).mean(axis=1)-rng.choice(b,(resamples,len(b))).mean(axis=1))
    lo,hi=np.percentile(differences,[2.5,97.5])
    return {'mean_block_median_difference_A':float(factor*(a.mean()-b.mean())),
            'block_median_difference_ci_low_A':float(lo),'block_median_difference_ci_high_A':float(hi),
            'Ag_blocks':len(a),'IAg_blocks':len(b)}

def bh_adjust(pvalues):
    """Benjamini-Hochberg adjustment among finite tests; undefined tests stay NaN."""
    p=np.asarray(pvalues,dtype=float);out=np.full(p.shape,np.nan);ids=np.flatnonzero(np.isfinite(p))
    order=ids[np.argsort(p[ids])];n=len(order)
    if n:out[order]=np.minimum(1,np.minimum.accumulate((p[order]*n/np.arange(1,n+1))[::-1])[::-1])
    return out

def oriented_switches(frame,scores,response_threshold=.1):
    """Orient median/peak shifts by reference expectation, then subtract same-dataset WT."""
    frame=frame.copy();rows=[]
    for (dataset,a,b),group in frame.groupby(['dataset','resid1','resid2'],sort=True):
        wt=group[group.system=='WT']
        if len(wt)!=1:raise ValueError('One WT baseline per microswitch and dataset required')
        factor=ORIENTATION[a,b]
        for r in group.to_dict('records'):
            r['orientation_factor']=factor;r['S_A']=scores[dataset][r['system']]
            r['S_class']=classify(r['S_A'],1.)
            for metric in ['median','peak']:
                value=factor*r[f'{metric}_shift_A'];baseline=factor*float(wt.iloc[0][f'{metric}_shift_A'])
                response=value-baseline
                r[f'oriented_{metric}_shift_A']=value;r[f'R_{metric}_A']=response
                r[f'R_{metric}_class']=classify(response,response_threshold)
            rows.append(r)
    values=pd.DataFrame(rows);correlations=[]
    for (d,a,b),group in values[values.system!='WT'].groupby(['dataset','resid1','resid2'],sort=True):
        x=group.S_A.to_numpy();y=group.R_median_A.to_numpy()
        valid=len(x)>=3 and np.ptp(x)>0 and np.ptp(y)>0
        sr,sp=spearmanr(x,y) if valid else (np.nan,np.nan)
        pr,pp=pearsonr(x,y) if valid else (np.nan,np.nan)
        correlations.append({'dataset':d,'resid1':a,'resid2':b,'n_mutants':len(x),
                             'spearman_rho':sr,'spearman_p':sp,'pearson_r':pr,'pearson_p':pp})
    corr=pd.DataFrame(correlations)
    if not corr.empty:
        for d,group in corr.groupby('dataset'):
            for method in ['spearman','pearson']:
                corr.loc[group.index,method+'_q']=bh_adjust(group[method+'_p'])
    return values,corr

def half_summaries(traces):
    """Full, first-half and second-half means/SDs, with independently WT-referenced S."""
    rows=[];stats={}
    for (d,s,st),values in traces.items():
        v=np.asarray(values,dtype=float)
        if len(v)<4 or len(v)%2:raise ValueError('Even frame count >=4 required for equal halves')
        for window,x in [('full',v),('first_half',v[:len(v)//2]),('second_half',v[len(v)//2:])]:
            stats[d,s,st,window]=(float(x.mean()),float(x.std(ddof=1)),len(x))
    for d,s in sorted({(d,s) for d,s,st in traces}):
        for window in ['full','first_half','second_half']:
            ag,ia=[stats[d,s,st,window] for st in ['Ag','IAg']]
            wt=stats[d,'WT','Ag',window][0]-stats[d,'WT','IAg',window][0]
            delta=ag[0]-ia[0]
            rows.append({'dataset':d,'system':s,'window':window,'Ag_mean_A':ag[0],'Ag_sd_A':ag[1],
                         'IAg_mean_A':ia[0],'IAg_sd_A':ia[1],'Ag_frames':ag[2],'IAg_frames':ia[2],
                         'delta_A':delta,'WT_delta_A':wt,'S_A':delta-wt,'S_class':classify(delta-wt,1.)})
    return pd.DataFrame(rows)
