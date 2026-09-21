"""Plot study-specific cached marker distances and kink angles; no trajectory IO."""
from example_label import annotate_example
from pathlib import Path
import concurrent.futures as futures
import hashlib
import json
import sys
import time
import argparse
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
CACHE_ROOT = ROOT/'data/current/cache'
ORDER = ['WT', 'P183A', 'P183F', 'P183H', 'P60F', 'P60M', 'P228M', 'P228A']
COLORS = dict(zip(ORDER, ['#1f77b4','#ff7f0e','#2ca02c','#7f7f7f','#e377c2','#d62728','#9467bd','#8c564b']))
SITES = {228:['WT','P228A','P228M'],183:['WT','P183A','P183F','P183H'],60:['WT','P60F','P60M']}



def plots(results,missing):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':12,'axes.labelsize':13,
                         'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none',
                         'axes.titleweight':'bold','axes.labelweight':'bold','axes.linewidth':.8})
    records={r['label']:r for r in results}
    def get(dataset,system,state): return records.get(f'{system}_{dataset}_{state}')
    fig,axs=plt.subplots(1,3,figsize=(16.5,4.9),sharey=True)
    for ax,(site,systems) in zip(axs,SITES.items()):
        for i,system in enumerate(systems):
            for dataset,offset,marker in [('A',-.045,'o'),('B',.045,'D')]:
                for state,shift,color in [('Ag',-.16,'#ef425c'),('IAg',.16,'#397de9')]:
                    row=get(dataset,system,state)
                    if row is None: continue
                    v=row['angles'][str(site)]
                    ax.errorbar(i+shift+offset,v['mean'],yerr=v['block_sd'],fmt=marker,color=color,
                                mfc=color if dataset=='A' else 'white',mec=color,mew=1.5,
                                markersize=8,capsize=3.5,elinewidth=1.4,zorder=3)
        ax.set(title=f'Residue {site}',xticks=range(len(systems)),xticklabels=systems,
               xlim=(-.55,len(systems)-.45),ylim=(0,75),yticks=np.arange(0,71,10))
        ax.grid(axis='y',color='#dddddd',lw=.75)
        ax.set_axisbelow(True)
        ax.tick_params(labelsize=13)
    axs[0].set_ylabel('Kink angle (degrees)')
    fig.text(.013,.97,'(c)',fontsize=24,va='top',fontfamily='DejaVu Serif')
    handles=[Line2D([],[],color='#ef425c',lw=3,label='Ag-bound'),Line2D([],[],color='#397de9',lw=3,label='IAg-bound'),
             Line2D([],[],color='#444444',marker='o',linestyle='',markersize=9,label='Dataset A'),
             Line2D([],[],color='#444444',marker='D',mfc='white',linestyle='',markersize=9,label='Dataset B')]
    fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.065),ncol=4,frameon=False,fontsize=14,columnspacing=2.4)
    note='Error bars: SD of ten non-overlapping 100-ns block means.'
    if missing: note+='   Missing inputs: '+', '.join(missing)
    fig.text(.5,.025,note,ha='center',fontsize=10,color='#666666')
    fig.subplots_adjust(left=.065,right=.99,top=.82,bottom=.25,wspace=.14)
    for ext in ['png','svg','pdf']: annotate_example(fig);fig.savefig(OUT/f'kink_angles_processed.{ext}',dpi=260,facecolor='white')
    plt.close(fig)
    traces={r['label']:np.load(CACHE_ROOT/f'{r["label"]}.npz') for r in results}
    # One structural schematic; A and B each occupy two rows with common scales.
    fig=plt.figure(figsize=(24,11.8),facecolor='white')
    leg=fig.add_axes([.025,.64,.24,.29]);leg.axis('off')
    leg.legend(handles=[Line2D([],[],color=COLORS[s],lw=2.3,label=s) for s in ORDER],
               ncol=2,loc='center',frameon=False,fontsize=17,labelspacing=1.7,columnspacing=2.8,handlelength=2.7)
    struct=fig.add_axes([.004,.10,.26,.52]);struct.axis('off')
    image_path=ROOT/'assets/marker_structure.png'
    if image_path.exists():
        struct.imshow(plt.imread(image_path))
    else:
        struct.text(.5,.5,'TM2–TM6 opening marker\nCα(Leu47) – Cα(Leu215)\n\nStructural image not supplied',ha='center',va='center',fontsize=15)
    bins=np.arange(8,25.0001,.1)
    histmax=max(np.histogram(v['distance_A'],bins=bins,density=True)[0].max() for v in traces.values())
    histlim=max(1.05,np.ceil(histmax*1.05*10)/10)
    differences=[]
    for dataset,ys,heading,letters in [('A',[.76,.545],.955,'bcdefgh'),('B',[.29,.075],.485,'ijklmno')]:
        fig.text(.30,heading,f'Dataset {dataset}',fontsize=20,fontweight='bold')
        for rowidx,state in enumerate(['Ag','IAg']):
            y=ys[rowidx]
            ax=fig.add_axes([.30,y,.22,.15])
            meanax=fig.add_axes([.566,y,.118,.15])
            histax=fig.add_axes([.735,y,.118,.15])
            for a,l in zip([ax,meanax,histax],letters[rowidx*3:rowidx*3+3]):
                a.text(-.1,1.09,f'({l})',transform=a.transAxes,fontsize=15,fontweight='bold')
            for idx,s in enumerate(ORDER):
                r=get(dataset,s,state)
                if r is None:
                    meanax.text(idx,10,'pending',ha='center',rotation=90,fontsize=8,color='#666666')
                    continue
                t=traces[r['label']]
                d=t['distance_A']
                ax.plot(t['time_ns'],d,color=COLORS[s],lw=.40,rasterized=True)
                meanax.bar(idx,d.mean(),yerr=d.std(ddof=1),color=COLORS[s],width=.73,
                           error_kw={'elinewidth':.8,'capsize':2,'ecolor':'#444444'})
                counts,_=np.histogram(d,bins=bins,density=True)
                histax.stairs(counts,bins,color=COLORS[s],lw=.85)
            ax.set(xlim=(0,1000),ylim=(8,25),xticks=[0,250,500,750,1000],yticks=[10,15,20,25],
                   xlabel='Time (ns)',ylabel='Distance (Å)')
            ax.grid(axis='y',alpha=.13)
            ax.text(.015,.97,f'{state}-bound',transform=ax.transAxes,va='top',fontsize=10,
                    bbox={'facecolor':'white','edgecolor':'none','alpha':.8,'pad':1})
            meanax.set(ylim=(8,25),yticks=[10,15,20,25],ylabel='Distance (Å)',xticks=range(8),xticklabels=ORDER)
            meanax.tick_params(axis='x',rotation=60,labelsize=9)
            for t in meanax.get_xticklabels():t.set_ha('right')
            histax.set(xlim=(8,25),ylim=(0,histlim),xticks=[10,15,20,25],xlabel='Distance (Å)',ylabel='Density (Å⁻¹)')
        da=fig.add_axes([.925,ys[1]+.012,.061,.35])
        valid=[]
        for s in ORDER:
            ag,ia=get(dataset,s,'Ag'),get(dataset,s,'IAg')
            if ag and ia:
                delta=ag['distance_mean_A']-ia['distance_mean_A']
                valid.append((s,delta))
                differences.append({'dataset':dataset,'system':s,'Ag_minus_IAg_A':delta})
        valid.sort(key=lambda v:v[1])
        da.barh(np.arange(len(valid)),[v for s,v in valid],color=[COLORS[s] for s,v in valid],height=.67)
        da.axvline(0,color='#444444',lw=.8)
        da.set(yticks=range(len(valid)),yticklabels=[s for s,v in valid],xlim=(-1.2,10),xticks=[0,3,6,9],
               xlabel='⟨d$_{Ag}$⟩−⟨d$_{IAg}$⟩ (Å)')
        da.tick_params(axis='y',labelsize=10)
        da.xaxis.label.set_size(11)
        da.text(-.12,1.05,f'({letters[-1]})',transform=da.transAxes,fontsize=15,fontweight='bold')
        for i,(s,v) in enumerate(valid): da.text(max(v+.15,.18),i,f'{v:+.2f}',va='center',fontsize=9)
    foot='Supplied framewise observations; raw traces, no smoothing. Error bars: framewise SD. Common 0.10 Å histogram bins.'
    if missing:foot+='   Missing inputs: '+', '.join(missing)
    fig.text(.5,.011,foot,ha='center',fontsize=10,color='#555555')
    for ext in ['png','svg','pdf']:annotate_example(fig);fig.savefig(OUT/f'activation_marker_processed.{ext}',dpi=230,facecolor='white')
    plt.close(fig)
    for v in traces.values():v.close()
    return differences
