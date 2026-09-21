"""Read completed GPCRmd packages without changing any coordinates or files in them."""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np

SYSTEMS=['WT','P183A','P183F','P183H','P60F','P60M','P228M','P228A']
SITES={228:['WT','P228A','P228M'],183:['WT','P183A','P183F','P183H'],60:['WT','P60F','P60M']}

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def validate(package):
    manifest=json.loads((package/'manifest.json').read_text())
    completion=json.loads((package/'COMPLETE.json').read_text())
    assert completion.get('checksums'), f'{package}: no completion checksums'
    assert manifest['label']==package.name
    assert manifest['frame_count']==10000, f'{package}: expected full 10,000-frame trajectory'
    sources=manifest['source_files']['trajectory']
    if isinstance(sources,str):sources=[sources]
    if package.name.startswith('P183H_A_'):
        expected='P183HSP_Active_2' if package.name.endswith('_Ag') else 'P183HSP_Inactive_2'
        assert sources and all(expected in p for p in sources), 'Dataset A P183H must use replacement _2 sources'
    top=package/'upload'/manifest['gpcrmd_files']['topology']
    traj=package/'upload'/manifest['gpcrmd_files']['trajectory']
    assert top.is_file() and traj.is_file()
    assert sha(top)==completion['checksums'][top.name], 'Topology checksum mismatch'
    return manifest,top,traj

def observables(xyz,times,system):
    assert xyz.shape==(10000,282,3) and np.isfinite(xyz).all()
    assert np.allclose(np.diff(times),.1,atol=1e-5)
    adjacent=np.linalg.norm(np.diff(xyz,axis=1),axis=2)
    assert adjacent.max()<5, 'Receptor may be split across periodic images'
    data={'time_ns':times,'distance_A':np.linalg.norm(xyz[:,46].astype(float)-xyz[:,214],axis=1)}
    angles={}
    block=np.floor((times-times[0])/100+1e-8).astype(int)
    assert np.array_equal(np.bincount(block),np.full(10,1000))
    for site,systems in SITES.items():
        if system not in systems:continue
        axes=[]
        for selection in [slice(site-7,site-1),slice(site,site+6)]:
            positions=xyz[:,selection].astype(float)
            positions-=positions.mean(axis=1,keepdims=True)
            axes.append(np.linalg.svd(positions,full_matrices=False)[2][:,0])
        angle=np.degrees(np.arccos(np.clip(np.abs(np.einsum('ij,ij->i',*axes)),0,1)))
        means=[float(angle[block==i].mean()) for i in range(10)]
        data[f'kink_{site}_deg']=angle
        angles[str(site)]={'mean':float(angle.mean()),'sd':float(angle.std(ddof=1)),
                           'block_means':means,'block_sd':float(np.std(means,ddof=1))}
    d=data['distance_A']
    stats={'frames':len(times),'first_time_ns':float(times[0]),'last_time_ns':float(times[-1]),
           'distance_mean_A':float(d.mean()),'distance_sd_A':float(d.std(ddof=1)),
           'distance_min_A':float(d.min()),'distance_max_A':float(d.max()),'angles':angles,
           'maximum_adjacent_CA_distance_A':float(adjacent.max())}
    return data,stats

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('packages',type=Path,help='Directory containing dataset_A and dataset_B')
    p.add_argument('--out',type=Path,default=Path(__file__).resolve().parent/'data/refreshed')
    p.add_argument('--labels',nargs='+',help='Optional subset; missing systems remain explicitly missing')
    p.add_argument('--check-only',action='store_true',help='Validate manifests/topologies, without reading trajectories')
    args=p.parse_args();results=[];missing=[]
    if not args.check_only:
        import MDAnalysis as mda
        if (args.out/'statistics_and_provenance.json').exists():
            raise FileExistsError('Use a new --out directory to preserve earlier results')
        (args.out/'cache').mkdir(parents=True,exist_ok=True)
    for dataset in ['A','B']:
        for system in SYSTEMS:
            for state in ['Ag','IAg']:
                label=f'{system}_{dataset}_{state}'
                package=args.packages/f'dataset_{dataset}'/label
                if (args.labels and label not in args.labels) or not (package/'COMPLETE.json').exists():
                    missing.append(label);continue
                manifest,top,traj=validate(package)
                print(('CHECK ' if args.check_only else 'READ ')+label,flush=True)
                if args.check_only:continue
                u=mda.Universe(str(top),str(traj));ca=u.select_atoms('segid PROA and name CA')
                assert len(ca)==282 and np.array_equal(ca.resids,np.arange(1,283))
                assert ca.resnames[46]==ca.resnames[214]=='LEU'
                xyz=np.empty((len(u.trajectory),282,3),dtype=np.float32);times=np.empty(len(u.trajectory))
                for i,ts in enumerate(u.trajectory):xyz[i]=ca.positions;times[i]=ts.time/1000
                u.trajectory.close()
                data,stats=observables(xyz,times,system)
                stats.update(label=label,system=system,dataset=dataset,state=state,
                             original_sources=manifest['source_files'],
                             fingerprint={'manifest_sha256':sha(package/'manifest.json'),
                                          'topology_sha256':sha(top),'trajectory_path':str(traj),
                                          'trajectory_size':traj.stat().st_size,
                                          'trajectory_mtime_ns':traj.stat().st_mtime_ns})
                np.savez_compressed(args.out/'cache'/f'{label}.npz',**data)
                results.append(stats)
    if not args.check_only:
        (args.out/'statistics_and_provenance.json').write_text(json.dumps({'results':results,'missing':missing},indent=2),encoding='utf-8')
    print('Missing/not selected:',missing)

if __name__=='__main__':main()
