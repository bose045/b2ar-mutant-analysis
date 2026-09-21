"""Deterministic numerical tests without study data or MDAnalysis."""
from pathlib import Path
import json
import tempfile
import unittest
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import reproduce
from examples.make_demo import build, SCORES
from extract_observables import observables
from example_label import annotate_example,set_example

class AnalysisTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp=tempfile.TemporaryDirectory()
        cls.root=Path(cls.temp.name)/'inputs'
        cls.records=build(cls.root)

    @classmethod
    def tearDownClass(cls):cls.temp.cleanup()

    def test_demo_is_explicit_and_complete(self):
        self.assertTrue(json.loads((self.root/'EXAMPLE_DATA.json').read_text())['synthetic'])
        self.assertEqual(len(self.records),32)
        with self.assertRaises(FileExistsError):build(self.root)

    def test_summary_matches_arrays(self):
        for record in self.records:
            with np.load(self.root/'data/current/cache'/f'{record["label"]}.npz') as arrays:
                self.assertAlmostEqual(arrays['distance_A'].mean(),record['distance_mean_A'],places=12)
                for site,stats in record['angles'].items():
                    angle=arrays[f'kink_{site}_deg']
                    self.assertAlmostEqual(angle.mean(),stats['mean'],places=12)
                    self.assertAlmostEqual(angle.reshape(10,1000).mean(axis=1).std(ddof=1),stats['block_sd'],places=12)

    def test_straight_chain_has_zero_bend(self):
        xyz=np.zeros((10000,282,3),dtype=np.float32)
        xyz[:,:,0]=np.arange(282)*3.8
        arrays,stats=observables(xyz,np.arange(10000)*.1,'WT')
        for site in ['60','183','228']:self.assertAlmostEqual(stats['angles'][site]['mean'],0,places=7)
        np.testing.assert_allclose(arrays['distance_A'],168*3.8,atol=1e-4)

    def test_rejects_wrong_frame_count(self):
        with self.assertRaises(AssertionError):observables(np.zeros((2,282,3)),np.array([0,.1]),'WT')

    def test_regional_consensus(self):
        helper=reproduce.module('plot_regional_histogram_matrices');helper.SCORES=SCORES
        caches={d:np.load(self.root/f'data/historical/regions/dataset_{d}_framewise_interhelix_regions.npz') for d in ['A','B']}
        try:
            _,records=helper.build_class_consensus(caches)
            self.assertEqual(len(records),126)
            self.assertTrue(all(r['n_observations']==5 for r in records))
        finally:
            for cache in caches.values():cache.close()

    def test_robust_contact_intersection(self):
        helper=reproduce.module('network_plotting_common')
        table=helper.make_robust_table(self.root/'data/historical/contacts')
        self.assertEqual(len(table),9)
        self.assertTrue((table.abs_functional_shift_msgm>=.3).all())

    def test_all_plot_modules_import(self):
        for name in ['current_replot','analyze_interhelix_regions','plot_regional_histogram_matrices',
                     'network_plotting_common','make_datasetB_figure5']:
            self.assertIsNotNone(reproduce.module(name))

    def test_synthetic_label_only_once(self):
        set_example(True);fig=plt.figure()
        annotate_example(fig);annotate_example(fig)
        self.assertEqual(len(fig.texts),1)
        self.assertIn('SYNTHETIC',fig.texts[0].get_text())
        plt.close(fig);set_example(False)

if __name__=='__main__':unittest.main()
