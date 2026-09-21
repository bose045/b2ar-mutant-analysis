"""Numerical geometry; Å, ns and degrees. Input coordinates are never modified."""
import numpy as np
from scipy.signal import find_peaks

REGIONS={'TM1':(2,36),'TM2':(39,68),'TM3':(77,108),'TM4':(119,142),
         'TM5':(169,197),'TM6':(209,238),'TM7':(245,268)}
# Ligand alignment follows the original selection (TM1 begins at residue 1).
FIT_REGIONS={**REGIONS,'TM1':(1,36)}
PAIRS=[(26,265),(26,266),(44,265),(22,263),(51,92),(96,218),(96,219),
       (96,262),(96,265),(99,215),(99,265),(103,215),(103,265),(184,222),
       (188,219),(195,215),(218,262),(222,258),(226,258),(226,93),(266,273)]

def windows():
    result={}
    for plane in ['extracellular','middle','intracellular']:
        result[plane]={}
        for helix,(a,b) in REGIONS.items():
            if plane=='middle':start=a+(b-a)//2-1
            else:
                first=(int(helix[-1])%2==1)==(plane=='extracellular')
                start=a if first else b-2
            result[plane][helix]=np.arange(start,start+3)-1
    return result

def regional_distances(ca):
    """Arithmetic three-CA centroids, using the archived region windows."""
    out={}
    for plane,groups in windows().items():
        centers={h:np.asarray(ca)[indices].mean(axis=0) for h,indices in groups.items()}
        for i in range(1,8):
            for j in range(i+1,8):
                out[f'{plane}_TM{i}_TM{j}']=float(np.linalg.norm(centers[f'TM{i}']-centers[f'TM{j}']))
    return out

def kink(ca,site):
    axes=[]
    for points in [ca[site-7:site-1],ca[site:site+6]]:
        if points.shape!=(6,3):raise ValueError('Kink requires six CA atoms per side')
        centered=np.asarray(points,dtype=float)-np.mean(points,axis=0)
        _,s,v=np.linalg.svd(centered,full_matrices=False)
        if s[0]<1e-10:raise ValueError('Degenerate helix axis')
        axes.append(v[0])
    return float(np.degrees(np.arccos(np.clip(abs(np.dot(*axes)),0,1))))

def minimum_distance(first,second):
    if not len(first) or not len(second):raise ValueError('Empty distance selection')
    return float(np.linalg.norm(np.asarray(first)[:,None,:]-np.asarray(second)[None,:,:],axis=-1).min())

def rigid_fit(mobile,reference):
    """Row-vector Kabsch fit, reflection corrected; no mass weighting."""
    mobile=np.asarray(mobile,dtype=float);reference=np.asarray(reference,dtype=float)
    if mobile.shape!=reference.shape or mobile.ndim!=2 or mobile.shape[1]!=3 or len(mobile)<3:
        raise ValueError('Matching fit atoms required')
    mc=mobile.mean(axis=0);rc=reference.mean(axis=0)
    u,_,vt=np.linalg.svd((mobile-mc).T@(reference-rc))
    correction=np.eye(3);correction[-1,-1]=np.linalg.det(u@vt)
    rotation=u@correction@vt
    return rotation,rc-mc@rotation

def ligand_metrics(fit,ref_fit,ligand,ref_ligand,cation,asp,cell=None):
    """Whole-ligand nearest-image translation, then receptor-fit ligand RMSD.

    cell uses row lattice vectors. For rotated/fitted coordinates, do not supply
    the original unrotated box. Use cell=None only for a validated whole complex.
    """
    ligand=np.array(ligand,dtype=float,copy=True);cation=np.asarray(cation,dtype=float)
    asp=np.asarray(asp,dtype=float)
    if asp.shape!=(2,3) or ligand.shape!=np.shape(ref_ligand):raise ValueError('Ligand/Asp atom mismatch')
    shift=np.zeros(3)
    if cell is not None:
        cell=np.asarray(cell,dtype=float)
        delta=cation-asp.mean(axis=0)
        shift=-np.round(np.linalg.solve(cell.T,delta))@cell
    distance=minimum_distance(asp,[cation+shift])
    rotation,translation=rigid_fit(fit,ref_fit)
    aligned=(ligand+shift)@rotation+translation
    rmsd=float(np.sqrt(np.mean(np.sum((aligned-ref_ligand)**2,axis=1))))
    return rmsd,distance,int(distance<=4)

def summarize(values,times,block_ns=100):
    values=np.asarray(values,dtype=float);times=np.asarray(times,dtype=float)
    if len(values)<2 or not np.isfinite(values).all() or not np.isfinite(times).all() or np.any(np.diff(times)<=0):
        raise ValueError('Finite observations and increasing times required')
    dt=np.diff(times)
    if not np.allclose(dt,dt[0],atol=1e-5):raise ValueError('Nonuniform frame spacing')
    size=round(block_ns/dt[0])
    if size<1 or not np.isclose(size*dt[0],block_ns,atol=1e-5) or len(values)%size:
        raise ValueError('Block length must divide the selected trajectory exactly')
    means=values.reshape(-1,size).mean(axis=1)
    if len(means)<2:raise ValueError('At least two blocks required for block SD')
    return {'mean':float(values.mean()),'sd':float(values.std(ddof=1)),
            'block_means':means.tolist(),'block_sd':float(means.std(ddof=1))}

def dominant_peak(values,edges):
    density,_=np.histogram(values,bins=edges,density=True)
    peaks,_=find_peaks(density,prominence=.01)
    index=int(peaks[np.argmax(density[peaks])]) if len(peaks) else int(np.argmax(density))
    return float((edges[index]+edges[index+1])/2)

def pair_summary(active,inactive,threshold=.25):
    """60 shared histogram bins per pair/system, matching the archived method."""
    lo=min(min(active),min(inactive));hi=max(max(active),max(inactive))
    edges=np.linspace(lo if hi>lo else lo-.5,hi if hi>lo else hi+.5,61)
    delta=dominant_peak(active,edges)-dominant_peak(inactive,edges)
    return {'median_shift_A':float(np.median(active)-np.median(inactive)),
            'peak_shift_A':delta,'direction':int(delta>=threshold)-int(delta<=-threshold),
            'threshold_A':threshold}
