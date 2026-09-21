"""Small, offline figure runner. See README.md for current versus historical inputs."""
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import shutil
import sys

CODE_ROOT=Path(__file__).resolve().parent
ROOT=CODE_ROOT
from example_label import annotate_example, set_example
sys.path.insert(0,str(CODE_ROOT/'vendor'))

def module(name):
    spec=importlib.util.spec_from_file_location(name,CODE_ROOT/'vendor'/f'{name}.py')
    value=importlib.util.module_from_spec(spec);spec.loader.exec_module(value)
    return value

def save(fig,path):
    import matplotlib.pyplot as plt
    annotate_example(fig)
    for ext in ['png','svg','pdf']:fig.savefig(path.with_suffix('.'+ext),dpi=220,bbox_inches='tight',facecolor='white')
    plt.close(fig)

def snapshot(out):
    inventory=json.loads((ROOT/'figure_inventory.json').read_text())
    for doc in inventory.values():
        for row in doc['figures']:
            for i,asset in enumerate(row['assets']):
                src=ROOT/asset
                dest=out/f'Figure_{row["figure"]}{"_"+str(i+1) if len(row["assets"])>1 else ""}{src.suffix}'
                shutil.copy2(src,dest)
            (out/f'Figure_{row["figure"]}_caption.txt').write_text(row['caption'],encoding='utf-8')
    print('Exported supplied embedded figures unchanged. These are manuscript snapshots, not recalculations.')

def current(out,data_dir):
    import numpy as np
    m=module('current_replot')
    m.ROOT=ROOT;m.OUT=out
    data=json.loads((data_dir/'statistics_and_provenance.json').read_text())
    m.CACHE_ROOT=data_dir/'cache'
    differences=m.plots(data['results'],data['missing'])
    lookup={(r['dataset'],r['system'],r['state']):r['distance_mean_A'] for r in data['results']}
    scores=[]
    for r in differences:
        if any((r['dataset'],'WT',state) not in lookup for state in ['Ag','IAg']):continue
        wt=lookup[r['dataset'],'WT','Ag']-lookup[r['dataset'],'WT','IAg']
        scores.append({**r,'S_A':r['Ag_minus_IAg_A']-wt})
    (out/'structural_shift_scores.json').write_text(json.dumps(scores,indent=2),encoding='utf-8')
    print('Current Figure 2c and Figure 3: missing observations:',data['missing'])

def regions(out,figures):
    import numpy as np
    import pandas as pd
    a=module('analyze_interhelix_regions');h=module('plot_regional_histogram_matrices')
    a.RESULTS=out
    scores=json.loads((ROOT/'regional_scores.json').read_text())
    a.SCORES=scores;h.SCORES=scores
    caches={d:np.load(ROOT/f'data/historical/regions/dataset_{d}_framewise_interhelix_regions.npz') for d in ['A','B']}
    if '4' in figures:
        summary=pd.concat([a.summaries(caches[d],d) for d in ['A','B']],ignore_index=True)
        assoc=a.associations(summary)
        a.paired_matrix_plot_with_sign_disagreement_x(assoc)
        assoc.to_json(out/'Figure_4_correlations.json',orient='records',indent=2)
    consensus,records=h.build_class_consensus(caches)
    (out/'regional_consensus.json').write_text(json.dumps(records,indent=2),encoding='utf-8')
    for figure,region in [('S3','extracellular'),('S4','middle'),('S5','intracellular')]:
        if figure in figures:save(h.make_region_comparison_figure(caches,region,consensus),out/f'Figure_{figure}')
    for c in caches.values():c.close()

def ec100(out):
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm
    data=json.loads((ROOT/'data/jones_EC100.json').read_text())
    residues=list('GAVLIMSTCNQFWYDERKH');positions=[88,168,211,288,323]
    x=np.full((19,5),np.nan)
    for r in data['rows']:x[residues.index(r['AA']),positions.index(int(r['Pos']))]=float(r['Norm'])
    # Values are supplied by the input table, not enforced against a study result.
    fig,ax=plt.subplots(figsize=(8,8.5))
    cmap=plt.get_cmap('RdBu_r').copy();cmap.set_bad('white')
    im=ax.imshow(x,cmap=cmap,norm=TwoSlopeNorm(vmin=0,vcenter=1,vmax=4),aspect='auto')
    ax.set(xticks=range(5),xticklabels=[f'P{p}\n(TM{tm})' for p,tm in zip(positions,[2,4,5,6,7])],yticks=range(19),yticklabels=residues,title='EC₁₀₀ (0.625 μM)')
    for i in range(19):
        for j in range(5):
            if np.isfinite(x[i,j]):ax.text(j,i,f'{x[i,j]:.2f}',ha='center',va='center',fontsize=10,color='white' if x[i,j]<.4 else 'black')
    for end in [1,5,8,10,13,15]:ax.axhline(end+.5,color='#444444',lw=1)
    fig.colorbar(im,ax=ax,shrink=.8,ticks=[0,.5,1,1.5,2,2.5,3,3.5,4],label='Activity relative to mean frameshift')
    save(fig,out/'Figure_S1')

def contacts(out,figures):
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from scipy.stats import spearmanr
    n=module('network_plotting_common');c=module('make_datasetB_figure5')
    folder=ROOT/'data/historical/contacts'
    if 'S6' in figures:
        edges=pd.read_csv(folder/'combined_AB_nonVDW_edges.csv')
        if edges.empty or (edges.contact_type=='vdw').any() or (edges.functional_shift_msgm.abs()<.3).any():
            raise ValueError('S6 requires nonempty non-VDW edges with |F| >= 0.3')
        fig=plt.figure(figsize=(12,13));grid=fig.add_gridspec(3,1,height_ratios=[7,1.2,2.5])
        ax=fig.add_subplot(grid[0]);_,tms=n.draw_network(ax,edges,'',edges.abs_functional_shift_msgm.min(),edges.abs_functional_shift_msgm.max(),panel_label='(a)',label_scores=True,tm_label_radius=15,axis_limit=16.3,score_font_scale=.8)
        leg=fig.add_subplot(grid[1]);leg.axis('off')
        handles=n.legend_handles(edges,tms)
        leg.legend(handles=handles,loc='center',ncol=5,frameon=False,fontsize=9)
        c.draw_clusters(fig.add_subplot(grid[2]),edges,font_scale=.95,node_scale=.8,layout_scale=1.1)
        save(fig,out/'Figure_S6')
    if 'S9' in figures:
        tables={k:pd.read_csv(folder/f'functional_shift_scores_{k}.csv') for k in ['strict','default','lenient']}
        fig,axes=plt.subplots(1,2,figsize=(11,5));stats=[]
        for ax,k in zip(axes,['strict','lenient']):
            columns=n.KEYS+['functional_shift_msgm']
            joined=tables['default'][columns].merge(tables[k][columns],on=n.KEYS,how='outer',suffixes=('_default','_other')).fillna(0)
            x=joined.functional_shift_msgm_default.to_numpy();y=joined.functional_shift_msgm_other.to_numpy()
            nonzero=(x!=0)|(y!=0);both=(x!=0)&(y!=0)
            rho=float(spearmanr(x[nonzero],y[nonzero]).statistic)
            ax.scatter(x,y,s=3,alpha=.2,rasterized=True);ax.plot([-1,1],[-1,1],'k--',lw=.8)
            ax.set(xlabel='Default F',ylabel=f'{k.capitalize()} F',title=f'Spearman ρ = {rho:.3f}')
            stats.append({'criterion':k,'n_union':len(x),'n_nonzero':int(nonzero.sum()),'rho':rho,'sign_agreement_nonzero_both':float(np.mean(np.sign(x[both])==np.sign(y[both])))})
        save(fig,out/'Figure_S9');(out/'S9_statistics.json').write_text(json.dumps(stats,indent=2))
    if 'S10' in figures:
        edges=[pd.read_csv(folder/'display_edges_F0.3_default.csv'),pd.read_csv(folder/'display_edges_F0.3_lenient.csv'),n.make_robust_table(folder)]
        fig,axes=plt.subplots(1,3,figsize=(24,9))
        low=min(v.abs_functional_shift_msgm.min() for v in edges);high=max(v.abs_functional_shift_msgm.max() for v in edges)
        for ax,e,title in zip(axes,edges,['Default','Lenient','Robust across all three']):n.draw_network(ax,e,title,low,high,label_scores=False,tm_label_radius=15,axis_limit=16.3)
        union=pd.concat(edges,ignore_index=True)
        present={n.assign_tm(int(v)) for v in np.r_[union.res1_num,union.res2_num]}
        tms=[tm for tm in n.TM_ORDER if tm in present]
        fig.legend(handles=n.legend_handles(union,tms),loc='lower center',ncol=7,frameon=False,fontsize=10)
        fig.subplots_adjust(bottom=.18)
        save(fig,out/'Figure_S10')

def ligand(out,dataset='A'):
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    order=['WT','P60F','P60M','P183A','P183F','P183HSP','P228A','P228M']
    fig,axes=plt.subplots(1,3,figsize=(15,5));rng=np.random.default_rng(20260820)
    metrics=[('ligand_rmsd_A','Ligand RMSD (Å)',1),('asp113_ionic_distance_A','Ligand N–Asp85 distance (Å)',1),('salt_bridge_le_4A','Salt-bridge occupancy (%)',100)]
    for i,s in enumerate(order):
        for state,color,offset in [('Active','#ef5a63',-.14),('Inactive','#4d8bea',.14)]:
            folder=ROOT/'data/historical/ligand'
            if (folder/f'dataset_{dataset}').is_dir():folder=folder/f'dataset_{dataset}'
            d=pd.read_csv(folder/f'{s}_{state}_ligand_pose.csv')
            for ax,(key,label,scale) in zip(axes,metrics):
                v=d[key].to_numpy()*scale
                blocks=np.array([v[j:j+500].mean() for j in range(0,len(v),500)])
                lo,hi=np.percentile(rng.choice(blocks,(20000,len(blocks)),replace=True).mean(axis=1),[2.5,97.5])
                ax.errorbar(i+offset,v.mean(),yerr=[[max(0,v.mean()-lo)],[max(0,hi-v.mean())]],fmt='o',color=color,capsize=3)
                ax.set(ylabel=label,xticks=range(8),xticklabels=[s.replace('HSP','H') for s in order])
                ax.tick_params(axis='x',rotation=50);ax.grid(axis='y',alpha=.2)
    axes[1].axhline(4,color='gray',ls='--')
    axes[2].set_ylim(0,105)
    axes[2].ticklabel_format(axis='y',style='plain',useOffset=False)
    axes[0].legend(handles=[Line2D([],[],marker='o',ls='',color='#ef5a63',label='Ag-bound'),Line2D([],[],marker='o',ls='',color='#4d8bea',label='IAg-bound')],frameon=False)
    fig.suptitle('Ligand metrics from supplied tables; verify trajectory provenance',color='#9b2525',fontsize=12)
    fig.tight_layout();save(fig,out/f'Figure_S7_dataset_{dataset}')

def verify():
    records=json.loads((ROOT/'source_manifest.json').read_text())
    failures=[]
    for r in records:
        p=ROOT/r['destination']
        if not p.exists() or hashlib.sha256(p.read_bytes()).hexdigest()!=r['sha256']:failures.append(r['destination'])
    if failures:raise RuntimeError('Changed/missing inputs: '+str(failures))
    print(f'Verified {len(records)} copied inputs.')

def main():
    global ROOT
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=['list','verify','snapshot','current','historical','structural'])
    parser.add_argument('--figures',nargs='+',default=['4','S1','S3','S4','S5','S6','S7','S9','S10'])
    parser.add_argument('--data-dir',type=Path,default=None,help='Current observables directory, e.g. data/refreshed')
    parser.add_argument('--input-root',type=Path,default=CODE_ROOT/'inputs')
    parser.add_argument('--output-root',type=Path,default=CODE_ROOT/'generated')
    parser.add_argument('--ligand-dataset',choices=['A','B'],default='A')
    args=parser.parse_args()
    ROOT=args.input_root.resolve()
    set_example((ROOT/'EXAMPLE_DATA.json').exists())
    if args.data_dir is None:args.data_dir=ROOT/'data/current'
    if args.mode=='list':
        print((CODE_ROOT/'docs/FIGURE_MAP.md').read_text(encoding='utf-8'));return
    if args.mode=='verify':verify();return
    import matplotlib
    matplotlib.use('Agg');matplotlib.rcParams['svg.fonttype']='none'
    out=args.output_root.resolve()/args.mode;out.mkdir(parents=True,exist_ok=True)
    if args.mode=='snapshot':snapshot(out)
    elif args.mode=='current':current(out,args.data_dir.resolve())
    elif args.mode=='historical':
        print('HISTORICAL REPLAY: not a corrected source-consistent manuscript. See FIGURE_MAP.md.')
        f=set(args.figures)
        if f&{'4','S3','S4','S5'}:regions(out,f)
        if 'S1' in f:ec100(out)
        if f&{'S6','S9','S10'}:contacts(out,f)
        if 'S7' in f:ligand(out,args.ligand_dataset)
        (out/'INPUT_WARNING.txt').write_text('Plots use supplied precomputed tables. Verify source trajectories, frame windows and class definitions before scientific interpretation. Synthetic examples are not simulation results.')
    elif args.mode=='structural':
        import subprocess
        work=out/'S12';work.mkdir(exist_ok=True)
        for p in (ROOT/'structural').iterdir():
            if p.is_file():shutil.copy2(p,work/p.name)
        shutil.copy2(CODE_ROOT/'structural/compose_publication_map.py',work/'compose_publication_map.py')
        subprocess.run([sys.executable,str(work/'compose_publication_map.py')],check=True)
    print('Outputs:',out)

if __name__=='__main__':main()
