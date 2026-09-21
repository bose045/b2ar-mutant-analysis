"""Sequential calculation-to-figure launcher; fails immediately on an unsuccessful stage."""
from pathlib import Path
import argparse
import subprocess
import sys

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('config',type=Path)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--contact-mode',choices=['existing-events','run-detector'],required=True)
    p.add_argument('--figures',nargs='+',default=['4','S6','S7','S9','S10'])
    p.add_argument('--ligand-dataset',choices=['A','B'],default='A')
    args=p.parse_args();root=Path(__file__).resolve().parent
    inputs=args.out.resolve()/'calculated_inputs';plots=args.out.resolve()/'plots'
    stages=['coordinates','score'] if args.contact_mode=='existing-events' else ['coordinates','contacts','score']
    command=[sys.executable,str(root/'calculate.py'),str(args.config.resolve()),'--out',str(inputs),'--stages',*stages]
    if args.contact_mode=='run-detector':command.append('--execute-contacts')
    subprocess.run(command,check=True)
    common=['--input-root',str(inputs),'--output-root',str(plots)]
    subprocess.run([sys.executable,str(root/'reproduce.py'),'current',*common],check=True)
    subprocess.run([sys.executable,str(root/'reproduce.py'),'historical',*common,'--figures',*args.figures,
                    '--ligand-dataset',args.ligand_dataset],check=True)
    print('Completed selected calculation/plot stages. Structural assets and external assay inputs are separate.')

if __name__=='__main__':main()
