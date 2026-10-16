"""Run with SciPy/NumPy and python_reference.mat exported by MATLAB."""
import importlib.util
import pathlib
import sys
import unittest
from unittest.mock import patch

import numpy as np
from numpy.testing import assert_allclose
from scipy.io import loadmat

ROOT = pathlib.Path(__file__).resolve().parent
PACKAGE = pathlib.Path(sys.argv.pop(1)) if len(sys.argv) > 1 and not sys.argv[1].startswith('-') else (ROOT.parent if ROOT.name == 'tests' else ROOT/'python_updated')
spec = importlib.util.spec_from_file_location('audit_package', PACKAGE/'__init__.py', submodule_search_locations=[str(PACKAGE)])
pkg = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = pkg
spec.loader.exec_module(pkg)
REF = loadmat(ROOT/'python_reference.mat', simplify_cells=True)


class Equivalence(unittest.TestCase):
    def test_full_pipelines(self):
        for case in REF['cases']:
            for corrected in (True, False):
                with self.subTest(method=case['method'], hop=case['window']-case['overlapLength'], voices=case['VoicesPerOctave'], corrected=corrected):
                    out = pkg.periodicity_analysis(REF['signal'], REF['Fs'], case['window_BF'], case['method'], case['base'],
                        window=case['window'], overlapLength=case['overlapLength'], FFTLength=case['FFTLength'],
                        VoicesPerOctave=case['VoicesPerOctave'], CI=corrected, BF_corr=corrected)
                    expected = case['result'] if corrected else case['uncorrected']
                    self.assertEqual(set(out), set(expected))
                    for key in out:
                        assert_allclose(out[key], expected[key], rtol=2e-9, atol=2e-10, err_msg=key)

    def test_direct_bayes_factors(self):
        for tail, ref in zip(('right', 'left', 'both'), REF['direct']):
            out = pkg.ttest_bf('T', REF['Tgrid'], 'N', REF['Ngrid'], 'tail', tail, return_diagnostics=True)
            for actual, key in zip(out[:5], ('b', 'p', 't', 'n', 'l')):
                assert_allclose(actual, ref[key], rtol=2e-9, atol=2e-10, err_msg=tail+'/'+key)
            self.assertTrue(np.all(out[5]['converged']))
            self.assertTrue(np.all(out[5]['relativeErrorEstimate'] <= 1e-10))

    def test_data_missing_and_ess(self):
        X, Y = REF['X'], REF['Y']
        options = [dict(rho1=[.5,-.5,.99,0], tail='right'), dict(rho1=[.5,-.5,.99,0], tail='left'), dict(paired=False)]
        for k in range(3):
            args = (X,) if k == 0 else (X, Y if k == 1 else Y[:8])
            out = pkg.ttest_bf(*args, return_log=True, **options[k])
            for actual, key in zip(out, ('b', 'p', 't', 'n', 'l')):
                assert_allclose(actual, REF['data'][k][key], rtol=2e-9, atol=2e-10)
        for rho in (-.5, 0, .5, .99):
            ne = max(2, min(20, 20*(1-rho)/(1+rho)))
            a = pkg.ttest_bf(T=2, N=20, rho1=rho, return_log=True)
            b = pkg.ttest_bf(T=2*np.sqrt(ne/20), N=ne, return_log=True)
            for x,y in zip(a,b): assert_allclose(x,y,atol=1e-12)

    def test_ci_ar1_and_orientation(self):
        X = REF['CI_X']
        for data, dim in ((X,1), (X.T,2)):
            out = pkg.credible_interval_individual_spectrum(data, np.arange(4), SampleDimension=dim, return_rho=True)
            for a, key in zip(out, ('mu','ci','df','rho')):
                assert_allclose(a, REF['CI_ref'][key], rtol=1e-11, atol=1e-12)
        for data in (np.ones((5,2)), np.zeros((5,2))):
            with self.assertRaisesRegex(ValueError, 'DegenerateData'):
                pkg.credible_interval_individual_spectrum(data, [1,2])
        with self.assertRaisesRegex(ValueError, 'AmbiguousOrientation'):
            pkg.credible_interval_individual_spectrum(np.eye(4), np.arange(4))
        pkg.credible_interval_individual_spectrum(np.eye(4), np.arange(4), SampleDimension=2)

    def test_postprocessing_and_comparison(self):
        assert_allclose(pkg.posterior_odds_smoothing(REF['oddsX'],REF['oddsF'],KernelSigma=2),REF['smooth'],atol=1e-12)
        assert_allclose(pkg.posterior_odds_temperature(REF['oddsX'],Temp=.7),REF['temperature'],atol=1e-12)
        assert_allclose(pkg.spectra_comparison(REF['X'].T,REF['Y'].T,'right'),REF['comparison'],rtol=1e-10)

    def test_simulation(self):
        for method, hop in (('STFT',1), ('STFT',3), ('CWT',1)):
            args = dict(signal_duration=601, Fs=20, window_BF=50, method=method, hop=hop, N_simulations=10, random_state=9)
            ref = pkg.non_oscillatory_periodicity_spectrum(**args)
            repeat = pkg.non_oscillatory_periodicity_spectrum(**args)
            for key in ref: assert_allclose(ref[key],repeat[key],rtol=0,atol=0)
            self.assertTrue(np.all(np.isfinite(ref['mean_psi_non_osc_fitted'])))
            if method == 'STFT':
                assert_allclose(ref['mean_psi_non_osc_fitted'],np.mean(ref['mean_psi_non_osc'][1:-1]))
            else:
                assert_allclose(ref['mean_psi_non_osc_fitted'],np.exp(np.polyval(ref['fit_mean'],np.log(ref['freqs']))))
            result = pkg.periodicity_analysis(np.random.default_rng(1).normal(size=601),20,50,method,ref,
                                               overlapLength=20-hop, CI=True)
            self.assertTrue(np.all(np.isfinite(result['bf_psi'])))

    def test_simulation_matlab_identical_noise(self):
        for k, method in enumerate(('STFT','CWT')):
            class Replay:
                column=0
                def standard_normal(self, length):
                    self_data=REF['noise'][:,self.column].copy()
                    self.column+=1
                    return self_data
            with patch('audit_package.simulation.np.random.default_rng', return_value=Replay()):
                result=pkg.non_oscillatory_periodicity_spectrum(601,20,50,method,hop=3 if k==0 else 1,N_simulations=10)
            for key in result:
                assert_allclose(result[key],REF['simulation_ref'][k][key],rtol=2e-8,atol=2e-9,err_msg=method+'/'+key)

    def test_validation_and_legacy_api(self):
        legacy = pkg.ttest_bf_directional([1,2,4,3])
        modern = pkg.ttest_bf([1,2,4,3],tail='right')
        self.assertEqual(len(legacy),4)
        for a,b in zip(legacy,modern): assert_allclose(a,b)
        self.assertEqual(len(pkg.credible_interval_individual_spectrum(REF['CI_X'],range(4))),3)
        for kwargs in (dict(T=1,N=20,tail='bad'),dict(T=1,N=20,rho1=1),dict(T=1,N=[10,20],rho1=.1),dict(T=np.inf,N=20)):
            with self.assertRaises(ValueError): pkg.ttest_bf(**kwargs)
        with self.assertRaises(ValueError): pkg.ttest_bf([1,1,1])
        with self.assertRaises(ValueError):
            pkg.periodicity_analysis(np.ones(600),20,1,'STFT',{'mean_psi_non_osc_fitted':np.ones(9)})


if __name__ == '__main__': unittest.main(verbosity=2)
