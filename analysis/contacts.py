"""GetContacts events -> per-type CP -> Ag-minus-IAg -> majority geometric mean -> F."""
from collections import Counter
import math
from pathlib import Path
import re
import numpy as np
import pandas as pd

TYPES=('wb','hbss','hbsb','sb','ts','ps','pc','hbbb','vdw','hp')
KEYS=['res1_num','res2_num','contact_type']
FREQUENCY_COLUMNS=KEYS+['res1_name','res2_name','frequency','count','frames']
MEMBRANE=((1,36),(39,68),(77,108),(119,142),(169,197),(209,238),(245,268),(270,280))
DEFAULTS={'sb_cutoff_dist':4.,'pc_cutoff_dist':6.,'pc_cutoff_ang':60.,'ps_cutoff_dist':7.,
          'ps_cutoff_ang':30.,'ps_psi_ang':45.,'ts_cutoff_dist':5.,'ts_cutoff_ang':30.,
          'ts_psi_ang':45.,'hbond_cutoff_dist':3.5,'hbond_cutoff_ang':70.,'vdw_epsilon':.5}

def residue(atom):
    tokens=atom.split(':')
    if len(tokens)<4:raise ValueError(f'Invalid atom identifier: {atom}')
    return ':'.join(tokens[:-3]),tokens[-3],int(tokens[-2])

def event_frequencies(path,frame_indices,chain=None):
    """At most one event per residue pair/type/frame. Empty frames remain in denominator.

    GetContacts water-bridge endpoints are its first two atom fields; later
    fields describe mediating water atoms, not additional receptor contacts.
    Numeric ordering removes residue-name-dependent orientation in mutants.
    """
    expected=set(map(int,frame_indices))
    if not expected:raise ValueError('No selected frames')
    counts=Counter();names={};seen=set();last=None;declared=None;chain_seen=set()
    with Path(path).open(encoding='utf-8') as stream:
        for line in stream:
            if not line.strip():continue
            if line.startswith('#'):
                match=re.search(r'total_frames\s*:\s*(\d+)',line)
                if match:declared=int(match[1])
                continue
            fields=line.split()
            if len(fields)<4:raise ValueError('Malformed GetContacts event')
            frame=int(fields[0]);kind=fields[1]
            if frame not in expected:raise ValueError(f'Unexpected frame {frame}')
            if last is not None and frame<last:raise ValueError('Events must be frame-sorted')
            if frame!=last:counts.update(seen);seen=set();last=frame
            if kind not in TYPES:continue
            first,second=residue(fields[2]),residue(fields[3])
            if chain is not None and (first[0]!=chain or second[0]!=chain):continue
            chain_seen.update([first[0],second[0]])
            if len(chain_seen)>1:raise ValueError('Multiple protein chains: specify a single receptor chain')
            if not (1<=first[2]<=282 and 1<=second[2]<=282):continue
            if first[2]==second[2]:continue
            a,b=sorted([first,second],key=lambda r:r[2])
            key=(a[2],b[2],kind)
            seen.add(key);names[key]=(a[1],b[1])
    counts.update(seen)  # includes final frame, even if it has just one event
    if declared!=len(expected):raise ValueError(f'Header frames {declared} != requested {len(expected)}')
    return pd.DataFrame([dict(zip(KEYS,k),res1_name=names[k][0],res2_name=names[k][1],
                              frequency=n/len(expected),count=n,frames=len(expected)) for k,n in sorted(counts.items())],columns=FREQUENCY_COLUMNS)

def majority_signed(values,epsilon=1e-6):
    """Archived majority-signed geometric mean; exact epsilon convention preserved."""
    v=np.asarray(values,dtype=float)
    if not np.isfinite(v).all():raise ValueError('Non-finite contact deltas')
    v=v[np.abs(v)>epsilon];pos=int((v>0).sum());neg=int((v<0).sum());total=len(v)
    sign=int(pos>neg)-int(neg>pos);agree=max(pos,neg)
    value=sign*math.exp(float(np.log(abs(v[np.sign(v)==sign])+epsilon).mean())) if sign else 0.
    return {'msgm':value,'maj_sign':sign,'n_agree':agree,'n_total':total,'agree_frac':agree/total if total else 0.}

def score(tables,scores,criterion='default',classes=None,min_agree=(3,3)):
    """tables maps (dataset,system,state) to CP tables; missing states are errors.

    Scores use complete Ag/IAg pairs. WT-like observations are excluded by S.
    Explicit classes permit reproduction of a documented historical grouping.
    """
    observations=sorted({(d,s) for d,s,state in tables})
    if classes is None:
        classes={}
        for d,s in observations:
            value=0 if s=='WT' else scores[d][s]
            classes[d,s]='AL' if value>=1 else ('IL' if value<=-1 else 'excluded')
    for d,s in observations:
        if any((d,s,st) not in tables for st in ['Ag','IAg']):raise ValueError(f'Missing state: {d}/{s}')
        for st in ['Ag','IAg']:
            table=tables[d,s,st]
            if table.duplicated(KEYS).any():raise ValueError('Duplicate contact keys')
            if not table.frequency.between(0,1).all():raise ValueError('CP outside [0,1]')
    if any(sum(classes.get(k)==c for k in observations)<n for c,n in zip(['AL','IL'],min_agree)):
        raise ValueError('Insufficient class observations for the requested consensus')
    all_keys=sorted(set().union(*(set(map(tuple,t[KEYS].to_numpy())) for t in tables.values())))
    lookup={k:t.set_index(KEYS).frequency.to_dict() for k,t in tables.items()}
    names={}
    for t in tables.values():
        for r in t.itertuples():names.setdefault((r.res1_num,r.res2_num,r.contact_type),(r.res1_name,r.res2_name))
    result=[];deltas=[]
    for key in all_keys:
        byclass={'AL':[],'IL':[]}
        for d,s in observations:
            ag=lookup[d,s,'Ag'].get(key,0.);ia=lookup[d,s,'IAg'].get(key,0.);delta=ag-ia
            cls=classes.get((d,s),'excluded')
            deltas.append(dict(zip(KEYS,key),dataset=d,system=s,CP_Ag=ag,CP_IAg=ia,delta_CP=delta,functional_class=cls))
            if cls in byclass:byclass[cls].append(delta)
        row=dict(zip(KEYS,key),res1_name=names[key][0],res2_name=names[key][1],criterion=criterion)
        for cls,suffix in [('AL','act'),('IL','inact')]:
            row.update({f'{k}_{suffix}':v for k,v in majority_signed(byclass[cls]).items()})
        f=row['msgm_act']-row['msgm_inact'];a,b,_=key
        row.update(functional_shift_msgm=f,abs_functional_shift_msgm=abs(f))
        row['eligible_edge']=row['n_agree_act']>=min_agree[0] and row['n_agree_inact']>=min_agree[1] and abs(a-b)>1 and all(any(lo<=r<=hi for lo,hi in MEMBRANE) for r in [a,b])
        row['display_edge_f03']=row['eligible_edge'] and abs(f)>=.3
        result.append(row)
    columns=KEYS+['res1_name','res2_name','criterion']+[f'{k}_{s}' for s in ['act','inact'] for k in ['msgm','maj_sign','n_agree','n_total','agree_frac']]+['functional_shift_msgm','abs_functional_shift_msgm','eligible_edge','display_edge_f03']
    return pd.DataFrame(result,columns=columns),pd.DataFrame(deltas)
