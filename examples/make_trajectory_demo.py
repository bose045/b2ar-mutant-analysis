"""Create artificial receptor-like coordinates and contact EVENT fixtures (not detector output).

Requires MDAnalysis. These are mathematical fixtures, not physically valid simulations.
"""
from pathlib import Path
import argparse
import json
import warnings
import numpy as np
import MDAnalysis as mda

SCORES={'A':{'WT':0,'P60F':2,'P60M':2.5,'P228M':1.5,'P228A':1.2,'P183H':.3,'P183A':-2,'P183F':-1.5},
        'B':{'WT':0,'P183H':1.7,'P60M':-1.4,'P183F':-2,'P228M':-1.8,'P60F':.1,'P228A':.2,'P183A':.3}}

def build(root):
    root=Path(root)
    if root.exists() and any(root.iterdir()):raise FileExistsError(root)
    root.mkdir(parents=True,exist_ok=True)
    # 282 CA plus two Asp85 O atoms and two ligand atoms; protein residue numbers preserved.
    resindices=np.r_[np.arange(282),84,84,282,282]
    u=mda.Universe.empty(286,n_residues=283,atom_resindex=resindices,trajectory=True)
    u.add_TopologyAttr('names',['CA']*282+['OD1','OD2','N1','C1'])
    u.add_TopologyAttr('resids',np.arange(1,284))
    names=['ALA']*283
    names[46]=names[214]='LEU';names[84]='ASP';names[282]='LIG'
    u.add_TopologyAttr('resnames',names);u.add_TopologyAttr('segids',['PROA'])
    u.add_TopologyAttr('chainIDs',['A']*286)
    u.add_TopologyAttr('elements',['C']*282+['O','O','N','C'])
    angle=np.arange(282)*2*np.pi/282
    ca=np.column_stack([9*np.cos(angle),9*np.sin(angle),np.zeros(282)])
    xyz=np.vstack([ca,ca[84]+[0,0,1],ca[84]+[0,0,-1],ca[84]+[2.8,0,0],ca[84]+[3.8,0,0]])
    u.atoms.positions=xyz
    with warnings.catch_warnings():
        warnings.simplefilter('ignore');u.atoms.write(str(root/'receptor.pdb'))
    direction=(ca[214]-ca[46]);direction/=np.linalg.norm(direction)
    weights=np.maximum(0,1-np.abs(np.arange(282)-214)/15)
    systems=[]
    for dataset in ['A','B']:
        for system,value in SCORES[dataset].items():
            for state in ['Ag','IAg']:
                label=f'{system}_{dataset}_{state}';traj=root/f'{label}.dcd'
                with mda.Writer(str(traj),n_atoms=286,dt=100) as writer:
                    for frame in range(20):
                        positions=xyz.copy()
                        positions[:282]+=weights[:,None]*direction*((2+value if state=='Ag' else 0)+.01*np.sin(frame))
                        u.atoms.positions=positions;writer.write(u.atoms)
                job={'dataset':dataset,'system':system,'state':state,'topology':'receptor.pdb',
                     'trajectories':[traj.name],'first':0,'last':19,'protein_selection':'segid PROA and resid 1:282',
                     'source_revision':'replacement_2' if system=='P183H' and dataset=='A' else 'synthetic',
                     'ligand':{'reference':'receptor.pdb','selection':'resname LIG and not name H*',
                               'cation_selection':'resname LIG and name N1','periodic_mode':'whole'},
                     'contact_events':{},'contact_chain':'A'}
                for criterion in ['strict','default','lenient']:
                    events=root/f'{label}_{criterion}.tsv'
                    # At one edge, AL gains contact; IL loses contact; WT-like has unchanged CP.
                    count=(16 if state=='Ag' else 4) if value>=1 else ((4 if state=='Ag' else 16) if value<=-1 else 10)
                    lines=['# total_frames:20\n']+[f'{f}\thp\tA:LEU:96:CA\tA:ILE:218:CA\n' for f in range(count)]
                    events.write_text(''.join(lines),encoding='utf-8');job['contact_events'][criterion]=events.name
                systems.append(job)
    config={'synthetic':True,'block_ns':1,'pair_threshold_A':.25,'require_ligand':True,'contact_profile':'pooled_S','systems':systems}
    (root/'study.json').write_text(json.dumps(config,indent=2),encoding='utf-8')
    return root/'study.json'

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,default=Path('generated/trajectory_demo'))
    args=p.parse_args();print(build(args.out))
