"""Calculate coordinate observables, run GetContacts and build contact scores.

Run: python calculate.py study.json --out inputs/recalculated --stages coordinates contacts score
The output is a plotting input root. Study files use explicit paths and frame windows.
"""
from pathlib import Path
import argparse
import hashlib
import json
import subprocess
import sys
import numpy as np
import pandas as pd
from analysis.geometry import (PAIRS,FIT_REGIONS,kink,regional_distances,minimum_distance,
                               ligand_metrics,summarize,pair_summary)
from analysis.contacts import DEFAULTS,TYPES,event_frequencies,score
from analysis.statistics import half_summaries,oriented_switches,median_block_ci,block_size,ORIENTATION

def dump(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2),encoding='utf-8')

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()

def path(base,value):
    p=Path(value);return p.resolve() if p.is_absolute() else (base/p).resolve()

def jobs_checked(config):
    jobs=config['systems'];labels=[]
    for job in jobs:
        d,s,st=job['dataset'],job['system'],job['state']
        label=f'{s}_{d}_{st}'
        if d not in ['A','B'] or st not in ['Ag','IAg'] or not s.isalnum():raise ValueError('Invalid system label')
        if label in labels:raise ValueError(f'Duplicate system: {label}')
        labels.append(label)
        if d=='A' and s=='P183H' and job.get('source_revision')!='replacement_2':
            raise ValueError('P183H A requires explicit source_revision: replacement_2')
        if job.get('first',0)<0 or job['last']<job.get('first',0) or job.get('stride',1)<1:raise ValueError('Invalid frame window')
    return jobs

def atom_ids(group):
    return [(int(a.resid),str(a.name)) for a in group]

def selections(universe,job):
    protein=universe.select_atoms(job.get('protein_selection','segid PROA'))
    ca=protein.select_atoms('name CA')
    if len(ca)!=282 or not np.array_equal(ca.resids,np.arange(1,283)) or ca.resnames[46]!='LEU' or ca.resnames[214]!='LEU':
        raise ValueError('Expected ordered receptor CA residues 1..282 and marker leucines')
    heavy={r:protein.select_atoms(f'resid {r} and not name H*') for pair in PAIRS for r in pair}
    if any(not len(g) for g in heavy.values()):raise ValueError('Empty microswitch residue selection')
    return ca,heavy

def coordinates(config,base,out):
    import MDAnalysis as mda
    current=out/'data/current';current.mkdir(parents=True,exist_ok=True)
    if (current/'statistics_and_provenance.json').exists():raise FileExistsError('Coordinate output exists; select a new output root')
    records=[];regional={'A':{},'B':{}};pairs={};provenance=[];ligand_expected=[];traces={};time_arrays={}
    for job in jobs_checked(config):
        d,s,st=job['dataset'],job['system'],job['state'];label=f'{s}_{d}_{st}'
        top=path(base,job['topology']);files=[path(base,p) for p in job['trajectories']]
        if not files or any(not p.is_file() for p in [top,*files]):raise FileNotFoundError(label)
        # File list order is explicit; no silent skipping of incomplete numbered segments.
        u=mda.Universe(str(top),*[str(p) for p in files]);ca,heavy=selections(u,job)
        first,last,stride=job.get('first',0),job['last'],job.get('stride',1)
        if last>=len(u.trajectory):raise ValueError(f'{label}: requested frame {last} exceeds trajectory')
        selected=range(first,last+1,stride);times=[];dist=[];angles={site:[] for site in [60,183,228]}
        region_values={};pair_values={pair:[] for pair in PAIRS};ligand_rows=[]
        ligand=job.get('ligand');reference=None
        if config.get('require_ligand',False) and not ligand:raise ValueError(f'{label}: ligand settings required')
        if ligand:
            reference=mda.Universe(str(path(base,ligand.get('reference_topology',job['topology']))),str(path(base,ligand['reference'])))
            refca,_=selections(reference,job)
            fit_ids=np.array([r-1 for a,b in FIT_REGIONS.values() for r in range(a,b+1)])
            ref_fit=refca.positions[fit_ids].copy()
            lig=u.select_atoms(ligand['selection']);rl=reference.select_atoms(ligand['selection'])
            if not len(lig) or atom_ids(lig)!=atom_ids(rl):raise ValueError(f'{label}: ligand atom identity/order mismatch')
            ref_lig=rl.positions.copy()
            cation=u.select_atoms(ligand['cation_selection'])
            asp=u.select_atoms(job.get('protein_selection','segid PROA')+' and resid 85 and name OD1 OD2')
            if len(cation)!=1 or len(asp)!=2:raise ValueError('Expected ligand cation and two Asp85 oxygens')
            if not set(cation.indices)<=set(lig.indices):raise ValueError('Cation must belong to ligand selection')
            if ligand.get('periodic_mode') not in ['whole','minimum_image']:raise ValueError('Specify ligand periodic_mode')
            if ligand['periodic_mode']=='minimum_image' and job.get('coordinates_rotated',False):
                raise ValueError('Rotated coordinates cannot use an unrotated periodic box')
        print('CALCULATE',label,flush=True)
        try:
            for frame in selected:
                ts=u.trajectory[frame];xyz=ca.positions.astype(float)
                if not np.isfinite(xyz).all() or np.linalg.norm(np.diff(xyz,axis=0),axis=1).max()>5:
                    raise ValueError(f'{label}/{frame}: receptor is not whole or coordinates invalid')
                times.append(float(ts.time)/1000);dist.append(float(np.linalg.norm(xyz[46]-xyz[214])))
                for site in angles:angles[site].append(kink(xyz,site))
                for key,value in regional_distances(xyz).items():region_values.setdefault(key,[]).append(value)
                for pair in PAIRS:pair_values[pair].append(minimum_distance(heavy[pair[0]].positions,heavy[pair[1]].positions))
                if ligand:
                    cell=ts.triclinic_dimensions if ligand['periodic_mode']=='minimum_image' else None
                    if ligand['periodic_mode']=='minimum_image' and cell is None:raise ValueError('Missing periodic cell')
                    rmsd,ionic,occupied=ligand_metrics(xyz[fit_ids],ref_fit,lig.positions,ref_lig,cation.positions[0],asp.positions,cell)
                    ligand_rows.append({'system':label,'frame':frame,'time_ns':times[-1],
                                        'ligand_rmsd_A':rmsd,'asp113_ionic_distance_A':ionic,'salt_bridge_le_4A':occupied})
        finally:
            u.trajectory.close()
            if reference is not None:reference.trajectory.close()
        times=np.array(times);dist=np.array(dist);block=config.get('block_ns',100.)
        traces[d,s,st]=dist;time_arrays[d,s,st]=times
        stats=summarize(dist,times,block)
        angle_stats={str(site):summarize(v,times,block) for site,v in angles.items()}
        values={'time_ns':times,'distance_A':dist,**{f'kink_{site}_deg':v for site,v in angles.items()}}
        (current/'cache').mkdir(exist_ok=True);np.savez_compressed(current/'cache'/f'{label}.npz',**values)
        record={'label':label,'system':s,'dataset':d,'state':st,'frames':len(times),
                'distance_mean_A':stats['mean'],'distance_sd_A':stats['sd'],'angles':angle_stats}
        records.append(record)
        legacy='Active' if st=='Ag' else 'Inactive'
        for key,v in region_values.items():regional[d][f'{s}_{legacy}_{key}']=np.array(v,dtype=np.float32)
        pairs[d,s,st]=pair_values
        if ligand:
            folder=out/'data/historical/ligand'/f'dataset_{d}';folder.mkdir(parents=True,exist_ok=True)
            pd.DataFrame(ligand_rows).to_csv(folder/f'{s.replace("P183H","P183HSP")}_{legacy}_ligand_pose.csv',index=False)
            ligand_expected.append(label)
        provenance.append({'label':label,'topology_sha256':sha(top),'trajectory_order':[str(p) for p in files],
                           'files':[{'path':str(p),'bytes':p.stat().st_size,'mtime_ns':p.stat().st_mtime_ns} for p in files],
                           'frame_indices':[first,last,stride],'ligand':ligand,'source_revision':job.get('source_revision')})
    expected=[f'{s}_{d}_{st}' for d in ['A','B'] for s in ['WT','P60F','P60M','P183A','P183F','P183H','P228A','P228M'] for st in ['Ag','IAg']]
    dump(current/'statistics_and_provenance.json',{'results':records,'missing':sorted(set(expected)-{r['label'] for r in records})})
    scores,rows=structural_scores(records)
    dump(out/'regional_scores.json',scores);pd.DataFrame(rows).to_csv(out/'structural_shift_scores.csv',index=False)
    half_summaries(traces).to_csv(out/'full_and_half_trajectory_scores.csv',index=False)
    folder=out/'data/historical/regions';folder.mkdir(parents=True,exist_ok=True)
    for d,values in regional.items():
        if values:np.savez_compressed(folder/f'dataset_{d}_framewise_interhelix_regions.npz',**values)
    motif_rows=[];folder=out/'data/microswitches';folder.mkdir(parents=True,exist_ok=True)
    for d in ['A','B']:
        values={f'{s}_{st}_{a}_{b}':np.array(v) for (ds,s,st),pv in pairs.items() if ds==d for (a,b),v in pv.items()}
        if values:np.savez_compressed(folder/f'dataset_{d}.npz',**values)
    rng=np.random.default_rng(config.get('bootstrap_seed',20260903))
    for d,s in sorted({(d,s) for d,s,st in pairs}):
        if any((d,s,st) not in pairs for st in ['Ag','IAg']):continue
        for pair in PAIRS:
            active,inactive=[pairs[d,s,st][pair] for st in ['Ag','IAg']]
            sizes=[block_size(time_arrays[d,s,st],config.get('microswitch_block_ns',block)) for st in ['Ag','IAg']]
            motif_rows.append({'dataset':d,'system':s,'resid1':pair[0],'resid2':pair[1],
                               **pair_summary(active,inactive,config.get('pair_threshold_A',.25)),
                               **median_block_ci(active,inactive,*sizes,rng,
                                                 config.get('bootstrap_resamples',10000),ORIENTATION[pair])})
    motif_table=pd.DataFrame(motif_rows)
    motif_table.to_csv(folder/'peak_and_median_shifts.csv',index=False)
    oriented,correlations=oriented_switches(motif_table,scores,config.get('response_threshold_A',.1))
    oriented.to_csv(folder/'WT_referenced_oriented_shifts.csv',index=False)
    correlations.to_csv(folder/'S_correlations.csv',index=False)
    dump(out/'coordinate_provenance.json',{'systems':provenance,'ligand_calculated':ligand_expected,
         'method':'whole receptor; Euclidean distances; minimum heavy-atom pairs; six-CA SVD; three-CA centroids'})

def structural_scores(records):
    means={(r['dataset'],r['system'],r['state']):r['distance_mean_A'] for r in records}
    scores={};rows=[]
    for d,s in sorted({(d,s) for d,s,st in means}):
        required=[(d,x,st) for x in [s,'WT'] for st in ['Ag','IAg']]
        if not all(k in means for k in required):raise ValueError(f'Both WT and system states required for S: {d}/{s}')
        delta=means[d,s,'Ag']-means[d,s,'IAg'];wt=means[d,'WT','Ag']-means[d,'WT','IAg'];value=delta-wt
        scores.setdefault(d,{})[s]=value
        rows.append({'dataset':d,'system':s,'delta_A':delta,'WT_delta_A':wt,'S_A':value})
    return scores,rows

def contacts(config,base,out,execute=False):
    settings=config.get('getcontacts',{});engine=path(base,settings['directory'])
    if not (engine/'get_dynamic_contacts.py').is_file():raise FileNotFoundError('GetContacts entry point not found')
    executable=settings.get('python',sys.executable)
    if settings.get('expected_commit'):
        revision=subprocess.check_output(['git','-C',str(engine),'rev-parse','HEAD'],text=True).strip()
        if revision!=settings['expected_commit']:raise ValueError('GetContacts revision mismatch')
        if subprocess.check_output(['git','-C',str(engine),'diff','--name-only'],text=True).strip():
            raise ValueError('GetContacts tracked files modified; review engine before use')
    dump(out/'getcontacts_engine.json',{'python':executable,'expected_commit':settings.get('expected_commit'),
         'python_file_hashes':{str(p.relative_to(engine)):sha(p) for p in sorted(engine.rglob('*.py'))}})
    plan=[]
    for job in jobs_checked(config):
        label=f'{job["system"]}_{job["dataset"]}_{job["state"]}'
        detector=job.get('contacts',{})
        first,last,stride=contact_window(job)
        if last==0:raise ValueError('The archived detector treats --end 0 as unlimited; choose at least two frames')
        trajectory=detector.get('trajectory')
        if trajectory is None:
            if len(job['trajectories'])!=1:raise ValueError('GetContacts requires one explicitly selected trajectory; do not silently concatenate')
            trajectory=job['trajectories'][0]
        for criterion,factor in [('strict',.95),('default',1),('lenient',1.05)]:
            dest=out/'contacts'/criterion/label;dest.mkdir(parents=True,exist_ok=True)
            events=dest/'events.tsv'
            args=[executable,str(engine/'get_dynamic_contacts.py'),'--topology',str(path(base,detector.get('topology',job['topology']))),
                  '--trajectory',str(path(base,trajectory)),'--beg',str(first),'--end',str(last),
                  '--stride',str(stride),'--cores',str(settings.get('cores',1)),'--itypes','all',
                  '--sele',detector.get('selection','protein'),'--solv',detector.get('solvent','resname WAT'),
                  '--hbond_res_diff','1','--vdw_res_diff','2','--output',str(events)]
            for key,value in DEFAULTS.items():args.extend(['--'+key,str(value*factor)])
            plan.append({'label':label,'criterion':criterion,'argv':args})
            dump(out/'getcontacts_plan.json',plan)
            if execute:
                if events.exists():raise FileExistsError(events)
                with (dest/'detector.log').open('w') as log:subprocess.run(args,stdout=log,stderr=subprocess.STDOUT,check=True)
    dump(out/'getcontacts_plan.json',plan)
    print('GetContacts executed' if execute else 'GetContacts commands written only; pass --execute-contacts to run')

def contact_window(job):
    override=job.get('contacts',{})
    first=override.get('first',job.get('first',0));last=override.get('last',job['last']);stride=override.get('stride',job.get('stride',1))
    if first<0 or last<first or stride<1:raise ValueError('Invalid contact frame window')
    return first,last,stride

def contact_scores(config,base,out):
    scores=json.loads((out/'regional_scores.json').read_text());jobs=jobs_checked(config)
    folder=out/'data/historical/contacts';folder.mkdir(parents=True,exist_ok=True)
    for criterion in ['strict','default','lenient']:
        tables={}
        for job in jobs:
            label=f'{job["system"]}_{job["dataset"]}_{job["state"]}'
            provided=job.get('contact_events',{}).get(criterion)
            events=path(base,provided) if provided else out/'contacts'/criterion/label/'events.tsv'
            first,last,stride=contact_window(job)
            selected=range(first,last+1,stride)
            frame_mode=job.get('contact_frame_labels','absolute')
            if frame_mode=='ordinal':selected=range(len(selected))
            elif frame_mode!='absolute':raise ValueError('contact_frame_labels must be absolute or ordinal')
            t=event_frequencies(events,selected,job.get('contact_chain'))
            (out/'contact_probabilities'/criterion).mkdir(parents=True,exist_ok=True)
            t.to_csv(out/'contact_probabilities'/criterion/f'{label}.csv',index=False)
            tables[job['dataset'],job['system'],job['state']]=t
        profile=config.get('contact_profile','pooled_S')
        classes=None;agreement=(3,3)
        if profile=='legacy_A':
            tables={k:v for k,v in tables.items() if k[0]=='A'}
            activating={'WT','P183H','P60F','P60M','P228A','P228M'}
            classes={(d,s):'AL' if s in activating else 'IL' for d,s,st in tables};agreement=(4,2)
        elif profile!='pooled_S':raise ValueError('Unknown contact profile')
        result,deltas=score(tables,scores,criterion,classes,agreement)
        result.to_csv(folder/f'functional_shift_scores_{criterion}.csv',index=False)
        deltas.to_csv(folder/f'observation_deltas_{criterion}.csv',index=False)
        display=result[result.display_edge_f03.astype(bool)]
        display.to_csv(folder/f'display_edges_F0.3_{criterion}.csv',index=False)
        if criterion=='default':
            display[display.contact_type!='vdw'].to_csv(folder/'combined_AB_nonVDW_edges.csv',index=False)
    dump(folder/'analysis_profile.json',{'classification':profile,'min_agree':agreement,
         'frequency_counter':'corrected frame-set counter v1; not the archived res_contacts_xl',
         'consensus':'majority sign excluding |delta|<=1e-6; geometric mean of matching |delta|+1e-6',
         'sensitivity':'12 geometric float criteria scaled +/-5%; sequence-separation exclusions unchanged'})

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('config',type=Path);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--stages',nargs='+',choices=['coordinates','contacts','score'],default=['coordinates'])
    p.add_argument('--execute-contacts',action='store_true')
    args=p.parse_args();config=json.loads(args.config.read_text());base=args.config.resolve().parent
    jobs_checked(config);out=args.out.resolve();out.mkdir(parents=True,exist_ok=True)
    signature=sha(args.config)
    record=out/'calculation_config.json'
    if record.exists() and json.loads(record.read_text())['sha256']!=signature:raise ValueError('Output belongs to a different config')
    dump(record,{'sha256':signature,'config':config})
    if config.get('synthetic',False):dump(out/'EXAMPLE_DATA.json',{'synthetic':True})
    for stage in args.stages:
        if stage=='coordinates':coordinates(config,base,out)
        elif stage=='contacts':contacts(config,base,out,args.execute_contacts)
        else:contact_scores(config,base,out)
    print('Calculation outputs:',out)

if __name__=='__main__':main()
