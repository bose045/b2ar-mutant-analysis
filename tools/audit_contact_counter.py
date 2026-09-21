"""Read-only diagnostic of an explicitly supplied GetContacts transformations.py.

The supplied module is Python code and must be trusted. No external files are changed.
"""
from pathlib import Path
import argparse
import importlib.util
import json
import sys
import tempfile
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from analysis.contacts import event_frequencies

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('transformations',type=Path)
    args=p.parse_args()
    spec=importlib.util.spec_from_file_location('archived_transformations',args.transformations)
    old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
    lines=['# total_frames:3\n']+[f'{i}\thp\tA:ALA:1:CA\tA:LEU:9:CA\n' for i in range(3)]
    records,total=old.res_contacts_xl(iter(lines),{'hp'})
    with tempfile.TemporaryDirectory() as tmp:
        events=Path(tmp)/'events.tsv';events.write_text(''.join(lines))
        corrected=event_frequencies(events,range(3),'A').iloc[0]
    print(json.dumps({'expected_contact_frames':3,'archived_contact_frames':len(records),
                      'archived_CP':len(records)/total,'corrected_contact_frames':int(corrected['count']),
                      'corrected_CP':float(corrected.frequency)},indent=2))

if __name__=='__main__':main()
