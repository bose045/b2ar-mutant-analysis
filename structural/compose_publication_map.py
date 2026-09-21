"""Lay out a genuine VMD ray trace with coordinate-derived vector annotations."""
from pathlib import Path
import csv
import re
import tkinter
import numpy as np
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parent
tcl = tkinter.Tcl()
matrices = {}
for line in (OUT/'publication_matrices.txt').read_text().splitlines():
    name, value = line.split(' ', 1)
    rows = tcl.splitlist(tcl.splitlist(value)[0])
    matrices[name] = np.array([list(map(float, tcl.splitlist(row))) for row in rows])
transform = matrices['global_matrix'] @ matrices['scale_matrix'] @ matrices['rotate_matrix'] @ matrices['center_matrix']
coords = {}
names = {}
with (OUT/'publication_coordinates.tsv').open() as fp:
    for row in csv.DictReader(fp, delimiter='\t'):
        r = int(row['resid'])
        coords[r] = (transform @ np.array([float(row[k]) for k in ('x','y','z')] + [1]))[:2]
        names[r] = row['resname']
with (OUT/'mapped_pairs.csv').open() as fp:
    pairs = list(csv.DictReader(fp))
regions = {'TM1':(2,36),'TM2':(39,68),'TM3':(77,108),'TM4':(119,142),
           'TM5':(169,197),'TM6':(209,238),'TM7':(245,268),'H8':(271,282)}
colors = dict(zip(regions, ['#2475b4','#ef8b23','#cf3e40','#987469','#c95fad','#a3ab24','#19a9b6','#737373']))
groups = [('TM3–TM6: hydrophobic lock, DRY, PIF/CWxP', [6,7,10,12,20], '#383838'),
          ('TM3–TM7 / NPxxY region', [8,9,11,13], '#78409e'),
          ('TM5–TM6 and TM6–TM7', [14,15,16,17,18,19], '#08787d'),
          ('TM1/TM2–TM7, Na-pocket and H8', [1,2,3,4,5,21], '#af6420')]
groupcolor = {i:color for _,ids,color in groups for i in ids}
used = sorted({int(re.search(r'\d+',row[k]).group()) for row in pairs for k in ('residue_1','residue_2')})
img = Image.open(OUT/'publication_switch_map_highres.tga').convert('RGB')
img.save(OUT/'publication_switch_map_highres.png')
w,h = img.size
scene = (OUT/'publication_switch_map_highres.dat').read_text()[:2000]
zoom = float(re.search(r'Zoom\s+([\d.]+)', scene).group(1))
# Tachyon orthographic camera: full vertical span = 1/Zoom, aspect = width/height.
halfheight = 0.5/zoom
halfwidth = halfheight*w/h
plt.rcParams.update({'font.family':'DejaVu Sans','svg.fonttype':'none','font.size':10})
fig = plt.figure(figsize=(12.2,9.5), facecolor='white')
ax = fig.add_axes([.035,.105,.57,.79])
ax.imshow(img, extent=(-halfwidth,halfwidth,-halfheight,halfheight), zorder=0)
ax.set_xlim(-.96,.96); ax.set_ylim(-.90,.94); ax.set_aspect('equal'); ax.axis('off')
# Overlay all pair guides so contacts behind the cartoon remain identifiable.
for row in pairs:
    i = int(row['table_row'])
    a,b = [int(re.search(r'\d+', row[k]).group()) for k in ('residue_1','residue_2')]
    xy = np.array([coords[a],coords[b]])
    ax.plot(xy[:,0],xy[:,1],color=groupcolor[i],lw=.75,alpha=.65,zorder=2)
for r in used:
    helix = next(k for k,(a,b) in regions.items() if a<=r<=b)
    ax.scatter(*coords[r],s=19,c=colors[helix],edgecolors='white',linewidths=.55,zorder=3)
# Names are tied to CA coordinates, with separate left/right label columns.
left = sorted(used, key=lambda r:coords[r][0])[:len(used)//2]
right = [r for r in used if r not in left]
for side,resids in [(-1,left),(1,right)]:
    resids.sort(key=lambda r:coords[r][1])
    raw = np.array([coords[r][1] for r in resids])
    positions = raw.copy()
    for j in range(1,len(positions)):
        positions[j] = max(positions[j],positions[j-1]+.107)
    positions -= (positions.mean()-raw.mean())
    if positions[-1]>.77: positions-=positions[-1]-.77
    if positions[0]<-.74: positions+=-.74-positions[0]
    for r,y in zip(resids,positions):
        x = side*.76
        p = coords[r]
        elbow = side*.59
        ax.plot([p[0],elbow,x-side*.025],[p[1],y,y],color='#9a9a9a',lw=.5,zorder=1)
        helix = next(k for k,(a,b) in regions.items() if a<=r<=b)
        ax.text(x,y,f'{names[r]}{r}',ha='right' if side<0 else 'left',va='center',fontsize=9.5,
                color=colors[helix],weight='bold',bbox=dict(facecolor='white',edgecolor='none',pad=.6),zorder=4)
ax.text(0,.88,'Intracellular side',ha='center',fontsize=10,color='#555555')
ax.text(0,-.85,'Extracellular side',ha='center',fontsize=10,color='#555555')
fig.text(.035,.932,'(a) Structural locations',fontsize=14,weight='bold')
fig.text(.64,.932,'(b) Analyzed residue pairs',fontsize=14,weight='bold')
key = fig.add_axes([.64,.12,.35,.765]); key.axis('off')
y = .99
for title,ids,color in groups:
    key.text(0,y,title,fontsize=11,weight='bold',color=color,va='top')
    y -= .044
    for i in ids:
        row = pairs[i-1]
        key.text(.01,y,f'{i:02d}',fontsize=9.4,color=color,va='top',weight='bold')
        key.text(.11,y,f"{row['residue_1']}–{row['residue_2']}",fontsize=10,va='top')
        key.text(.98,y,f"{row['generic_1']}–{row['generic_2']}",fontsize=9,color='#666666',va='top',ha='right')
        y -= .0315
    y -= .022
for j,(helix,color) in enumerate(colors.items()):
    x = .055+j*.071
    fig.add_artist(plt.Line2D([x,x+.019],[.083,.083],color=color,lw=5,transform=fig.transFigure))
    fig.text(x+.024,.083,helix,fontsize=9,va='center')
fig.text(.035,.04,'Prepared WT Ag-bound β₂AR • simulation residue numbering • all 21 switch-table pairs',fontsize=10,color='#444444')
fig.text(.035,.018,'Connections indicate pair locations, not contact occupancy or switching direction.',fontsize=9,color='#555555')
for suffix in ('png','svg','pdf'):
    fig.savefig(OUT/f'microswitch_structural_map_publication.{suffix}',dpi=350)
plt.close(fig)
print(f'Composed publication map from {w} × {h} VMD ray trace.')
print(f'All {len(pairs)} pairs and {len(used)} unique residue endpoints included.')
