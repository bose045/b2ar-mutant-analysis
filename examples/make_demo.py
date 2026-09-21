"""Generate deterministic, fictional inputs; no manuscript data are used."""
from pathlib import Path
import argparse
import itertools
import json
import hashlib
import numpy as np
import pandas as pd

SYSTEMS=['WT','P183A','P183F','P183H','P60F','P60M','P228M','P228A']
SITES={228:['WT','P228A','P228M'],183:['WT','P183A','P183F','P183H'],60:['WT','P60F','P60M']}
SCORES={'A':{'P60F':2.,'P60M':2.5,'P228M':1.5,'P228A':1.2,'P183H':.3,'P183A':-2.,'P183F':-1.5},
        'B':{'P183H':1.7,'P60M':-1.4,'P183F':-2.,'P228M':-1.8,'P60F':.1,'P228A':.2,'P183A':.3}}

def write_json(path,obj):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(obj,indent=2),encoding='utf-8')

def build(root):
    if root.exists() and any(root.iterdir()):
        raise FileExistsError(f'Choose an empty output directory: {root}')
    rng=np.random.default_rng(42);root.mkdir(parents=True,exist_ok=True)
    write_json(root/'EXAMPLE_DATA.json',{'synthetic':True,'seed':42,'description':'Fictional demonstration data; no scientific inference permitted.'})
    write_json(root/'regional_scores.json',SCORES)
    records=[];cache=root/'data/current/cache';cache.mkdir(parents=True)
    times=np.arange(10000)*.1
    for dataset in ['A','B']:
        regions={}
        for index,system in enumerate(SYSTEMS):
            score=SCORES[dataset].get(system,0)
            for state,legacy in [('Ag','Active'),('IAg','Inactive')]:
                label=f'{system}_{dataset}_{state}'
                distance=12+(4+score if state=='Ag' else 0)+.25*np.sin(times/60)+rng.normal(0,.2,len(times))
                values={'time_ns':times,'distance_A':distance};angles={}
                for site,systems in SITES.items():
                    if system not in systems:continue
                    angle=30+index+.5*score+3*np.sin(times/150)+rng.normal(0,1,len(times))
                    means=angle.reshape(10,1000).mean(axis=1)
                    values[f'kink_{site}_deg']=angle
                    angles[str(site)]={'mean':float(angle.mean()),'sd':float(angle.std(ddof=1)),
                                       'block_sd':float(means.std(ddof=1)),'block_means':means.tolist()}
                np.savez_compressed(cache/f'{label}.npz',**values)
                records.append({'label':label,'system':system,'dataset':dataset,'state':state,
                                'frames':len(times),'distance_mean_A':float(distance.mean()),
                                'distance_sd_A':float(distance.std(ddof=1)),'angles':angles,'synthetic':True})
                for plane in ['extracellular','middle','intracellular']:
                    for a,b in itertools.combinations(range(1,8),2):
                        shift=.25*score*(1 if (a+b)%2 else -1)
                        regions[f'{system}_{legacy}_{plane}_TM{a}_TM{b}']=rng.normal(15+b-a+(shift if state=='Ag' else 0),.3,200)
        path=root/'data/historical/regions';path.mkdir(parents=True,exist_ok=True)
        np.savez_compressed(path/f'dataset_{dataset}_framewise_interhelix_regions.npz',**regions)
    write_json(root/'data/current/statistics_and_provenance.json',{'results':records,'missing':[],'synthetic':True})
    # Synthetic contact tables exercise graph rendering and criterion comparisons.
    nodes=[(44,'ILE'),(96,'LEU'),(99,'ILE'),(102,'ASP'),(218,'ILE'),(222,'PHE'),(226,'TRP'),(262,'ASN'),(265,'ILE'),(272,'PHE')]
    rows=[]
    for i in range(len(nodes)-1):
        a,an=nodes[i];b,bn=nodes[i+1];f=(.4+.02*i)*(1 if i%2 else -1)
        rows.append({'res1_num':a,'res2_num':b,'res1_name':an,'res2_name':bn,'contact_type':'hp',
                     'functional_shift_msgm':f,'abs_functional_shift_msgm':abs(f),
                     'msgm_act':f/2,'msgm_inact':-f/2,'maj_sign_act':int(np.sign(f)),
                     'maj_sign_inact':int(-np.sign(f)),'n_agree_act':4,'n_agree_inact':4,
                     'eligible_edge':True,'display_edge_f03':True})
    contacts=root/'data/historical/contacts';contacts.mkdir(parents=True)
    pd.DataFrame(rows).to_csv(contacts/'combined_AB_nonVDW_edges.csv',index=False)
    for criterion,factor in [('strict',.92),('default',1),('lenient',1.08)]:
        table=pd.DataFrame(rows)
        for col in ['functional_shift_msgm','abs_functional_shift_msgm','msgm_act','msgm_inact']:table[col]*=factor
        for prefix in ['functional_shift_scores','display_edges_F0.3']:table.to_csv(contacts/f'{prefix}_{criterion}.csv',index=False)
    # Illustrative assay-shaped matrix, not Jones experimental measurements.
    rows=[{'AA':aa,'Pos':pos,'Norm':float(rng.uniform(.2,2.5))} for aa in 'GAVLIMSTCNQFWYDERKH' for pos in [88,168,211,288,323]]
    write_json(root/'data/jones_EC100.json',{'synthetic':True,'rows':rows})
    ligand=root/'data/historical/ligand';ligand.mkdir(parents=True)
    for system in SYSTEMS:
        for state in ['Active','Inactive']:
            table=pd.DataFrame({'ligand_rmsd_A':rng.uniform(1,3,1000),
                                'asp113_ionic_distance_A':rng.uniform(2.5,4.5,1000)})
            table['salt_bridge_le_4A']=(table.asp113_ionic_distance_A<=4).astype(int)
            table.to_csv(ligand/f'{system.replace("P183H","P183HSP")}_{state}_ligand_pose.csv',index=False)
    manifest=[]
    for path in sorted(root.rglob('*')):
        if path.is_file():manifest.append({'destination':path.relative_to(root).as_posix(),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
    write_json(root/'source_manifest.json',manifest)
    return records

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,default=Path(__file__).resolve().parent/'demo_inputs')
    args=p.parse_args();build(args.out.resolve());print('Synthetic inputs:',args.out.resolve())
