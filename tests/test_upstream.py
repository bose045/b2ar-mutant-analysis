from pathlib import Path
import tempfile
import unittest
import numpy as np
import pandas as pd
from analysis.contacts import event_frequencies,majority_signed,score,FREQUENCY_COLUMNS
from analysis.geometry import regional_distances,windows,kink,rigid_fit,ligand_metrics,pair_summary
from analysis.statistics import half_summaries,median_block_ci,bh_adjust,oriented_switches
from calculate import contact_window

class UpstreamTests(unittest.TestCase):
    def test_first_and_final_contact_and_duplicate_atoms(self):
        with tempfile.TemporaryDirectory() as tmp:
            f=Path(tmp)/'events.tsv'
            f.write_text('# total_frames:4\n0\thp\tA:ALA:1:CA\tA:LEU:9:CA\n'
                         '0\thp\tA:LEU:9:CB\tA:ALA:1:CB\n'
                         '2\thp\tA:ALA:1:CA\tA:LEU:9:CA\n'
                         '3\twb\tA:ALA:1:O\tA:LEU:9:N\tW:WAT:900:O\n')
            table=event_frequencies(f,range(4),'A').set_index('contact_type')
            self.assertEqual(table.loc['hp','frequency'],.5)
            self.assertEqual(table.loc['wb','frequency'],.25)
            with self.assertRaises(ValueError):event_frequencies(f,range(5),'A')

    def test_zero_frames_in_denominator(self):
        with tempfile.TemporaryDirectory() as tmp:
            f=Path(tmp)/'events.tsv';f.write_text('# total_frames:5\n')
            self.assertTrue(event_frequencies(f,range(5)).empty)

    def test_geometric_majority_and_tie(self):
        result=majority_signed([.2,.8,-.3,0,0])
        self.assertAlmostEqual(result['msgm'],np.sqrt((.2+1e-6)*(.8+1e-6)))
        self.assertEqual(result['n_agree'],2)
        self.assertEqual(majority_signed([.2,-.8])['msgm'],0)

    def test_score_complete_pairs_and_F(self):
        tables={};scores={'A':{}}
        for i in range(6):
            s=f'M{i}';scores['A'][s]=2 if i<3 else -2
            for st in ['Ag','IAg']:
                cp=.8 if (i<3)==(st=='Ag') else .2
                tables['A',s,st]=pd.DataFrame([[96,218,'hp','LEU','ILE',cp,8,10]],columns=FREQUENCY_COLUMNS)
        result,deltas=score(tables,scores)
        self.assertAlmostEqual(result.iloc[0].functional_shift_msgm,1.200002)
        self.assertTrue(result.iloc[0].display_edge_f03)
        del tables['A','M0','IAg']
        with self.assertRaises(ValueError):score(tables,scores)

    def test_region_windows(self):
        self.assertEqual(windows()['extracellular']['TM1'].tolist(),[1,2,3])
        self.assertEqual(windows()['intracellular']['TM2'].tolist(),[38,39,40])
        ca=np.column_stack([np.arange(282),np.zeros((282,2))])
        self.assertEqual(len(regional_distances(ca)),63)
        self.assertEqual(kink(ca,60),0)

    def test_rotation_translation_and_periodic_ligand(self):
        ref=np.array([[0.,0,0],[1,0,0],[0,1,0],[0,0,1]])
        rotation=np.array([[0,1,0],[-1,0,0],[0,0,1.]])
        moving=ref@rotation+3
        r,t=rigid_fit(moving,ref)
        np.testing.assert_allclose(moving@r+t,ref,atol=1e-12)
        lig=np.array([[2.8,0,0],[3.8,0,0]])
        rmsd,distance,occ=ligand_metrics(ref,ref,lig+np.array([10,0,0]),lig,lig[0]+[10,0,0],np.array([[0,0,0],[0,0,1]]),np.eye(3)*10)
        self.assertAlmostEqual(rmsd,0);self.assertAlmostEqual(distance,2.8);self.assertEqual(occ,1)

    def test_peak_shift(self):
        result=pair_summary(np.full(20,5.),np.full(20,3.))
        self.assertEqual(result['median_shift_A'],2)
        self.assertEqual(result['direction'],1)

    def test_halves_use_corresponding_WT(self):
        traces={('A','WT','Ag'):[4,4,6,6],('A','WT','IAg'):[1,1,1,1],
                ('A','M1','Ag'):[7,7,8,8],('A','M1','IAg'):[2,2,2,2]}
        values=half_summaries(traces)
        result=values[values.system=='M1'].set_index('window')
        self.assertEqual(result.loc['full','S_A'],1.5)
        self.assertEqual(result.loc['first_half','S_A'],2)
        self.assertEqual(result.loc['second_half','S_A'],1)

    def test_block_CI_and_BH(self):
        result=median_block_ci(np.full(20,4),np.ones(20),10,10,np.random.default_rng(1),100,factor=-1)
        self.assertEqual(result['mean_block_median_difference_A'],-3)
        self.assertEqual(result['block_median_difference_ci_low_A'],-3)
        self.assertEqual(result['block_median_difference_ci_high_A'],-3)
        np.testing.assert_allclose(bh_adjust([.01,.04,.03,np.nan]),[.03,.04,.04,np.nan],equal_nan=True)

    def test_oriented_response_subtracts_WT(self):
        records=[{'dataset':'A','system':s,'resid1':96,'resid2':265,
                  'median_shift_A':v,'peak_shift_A':v} for s,v in [('WT',1),('M1',2),('M2',3),('M3',4)]]
        values,corr=oriented_switches(pd.DataFrame(records),{'A':{'WT':0,'M1':1,'M2':2,'M3':3}})
        self.assertEqual(values[values.system=='M1'].iloc[0].R_median_A,-1)
        self.assertAlmostEqual(corr.iloc[0].spearman_rho,-1)

    def test_contact_window_override(self):
        self.assertEqual(contact_window({'last':9999,'contacts':{'last':8999,'stride':2}}),(0,8999,2))

if __name__=='__main__':unittest.main()
