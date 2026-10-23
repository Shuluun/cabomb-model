"""Mathematical tests retained from the retrospective methods study."""
import inspect
import json
import math
import unittest
from unittest.mock import patch
import numpy as np
from cabomb_release import methods
from cabomb_release.bridge import cdf, cohort_response

class MathematicalChecks(unittest.TestCase):
    def test_shapes_bounds_and_structural_zero(self):
        x=np.r_[0,np.geomspace(1e-8,1e5,200)]
        rng=np.random.default_rng(3)
        for family in methods.FITTED:
            p=methods.decode(family,rng.uniform(0,1,(20,len(methods.CONFIG['parameter_order'][family]))))
            y=methods.predict(family,x,p)
            self.assertTrue(np.isfinite(y).all(),family)
            self.assertTrue((y>=0).all() and (y<=100+1e-12).all(),family)
            self.assertTrue(np.all(np.diff(y,axis=1)>=-1e-12),family)
            np.testing.assert_array_equal(y[:,0],0)

    def test_naive_anchoring(self):
        x=np.array([.1,.2,.4,.7,1.])
        y=np.array([2.,5.,8.,11.,16.])
        q=methods.naive('last3OLS_slope_anchored_last',x,y,np.array([1.,2.,100.]))
        self.assertEqual(q[0],16)
        self.assertEqual(q[-1],100)
        decreasing=methods.naive('last3OLS_slope_anchored_last',x,np.array([10,9,8,7,6]),np.array([2,3]))
        np.testing.assert_array_equal(decreasing,[6,6])

    def test_synthetic_exp_recovery_and_seed(self):
        x=np.array([.05,.1,.2,.4,.6,1.])
        y=methods.predict('exp100',x,np.array([.7]))[0]
        a=methods.fit('exp100',x,y,'synthetic_exp')
        b=methods.fit('exp100',x,y,'synthetic_exp')
        self.assertTrue(a['usable'])
        self.assertLess(a['rmse_pp'],1e-8)
        self.assertEqual(a,b)

    def test_rolling_common_folds_and_fallback(self):
        fits={f:{'usable':True} for f in methods.FITTED}
        e={f:[1.,2.] for f in methods.FITTED+methods.NAIVE}
        e['exp100']=[0.,None]
        result=methods.select_rolling(methods.OLD,e,fits,'exp100',2)
        self.assertEqual(result['selected_family'],'expA')
        self.assertIn('exp100',result['excluded'])
        self.assertFalse(result['fallback'])
        low=methods.select_rolling(methods.OLD,{f:[1.] for f in methods.OLD},fits,'exp100',1)
        self.assertTrue(low['fallback'])
        self.assertEqual(low['selected_family'],'exp100')
        all_failed=methods.select_rolling(methods.OLD,{f:[None,None] for f in methods.OLD},fits,None,2)
        self.assertTrue(all_failed['unresolved'])

    def test_forecast_early_information_only(self):
        seen=[]
        def fake(family,x,y,curve_id,fold_suffix=''):
            seen.append((fold_suffix,len(y),float(x[-1])))
            return {'family':family,'usable':True,'unresolved':False,'rmse_pp':1.,
                    'parameters_scaled':methods.decode(family,np.full(len(methods.CONFIG['parameter_order'][family]),.5)).tolist()}
        t=np.array([0,1,2,4,6,8,12.])
        y=np.array([0,1,2,4,6,8,12.])
        with patch.object(methods,'fit',side_effect=fake):
            a=methods.forecast('fake',t,y,np.array([24.,48.]))
            nofold=methods.forecast('fake',t[:5],y[:5],np.array([24.,48.]))
            json.dumps(nofold,allow_nan=False)
            self.assertEqual(nofold['selectors']['old4Rolling']['cv_scores_pp'],{})
        self.assertEqual(a['selectors']['expandedRolling']['n_folds'],2)
        self.assertEqual(sorted(set(n for suffix,n,last in seen if suffix)),[4,5])
        self.assertTrue(all(last<=1 for suffix,n,last in seen))
        self.assertNotIn('late',inspect.signature(methods.forecast).parameters)
        # Held-out values are available only to scoring after the prefix fit.
        self.assertEqual(sorted(set(r['heldout_time_h'] for r in a['fold_rows'])),[8.,12.])

    def test_exponential_convolution_analytic(self):
        t=np.linspace(0,24,1921)
        q,p,total=cohort_response(t,np.ones_like(t),1)
        exact=t-4*(-np.expm1(-t/4))
        self.assertLess(float(np.max(np.abs(q-exact))),1e-5)
        np.testing.assert_allclose(q+p,total,rtol=0,atol=1e-12)
        self.assertTrue((p>=0).all())

    def test_zero_input_and_shifted_cohort(self):
        t=np.linspace(0,20,201)
        for beta in (.5,1,2):
            q,p,total=cohort_response(t,np.zeros_like(t),beta)
            np.testing.assert_array_equal(q,0)
            np.testing.assert_array_equal(p,0)
            # A unit cohort arriving at s=5 has zero response before arrival.
            arrival=5.
            age=np.maximum(0,t-arrival)
            response=cdf(age,beta)
            np.testing.assert_array_equal(response[t<=arrival],0)
            self.assertGreater(response[-1],0)

    def test_common_mean_wait(self):
        for beta in (.5,1,2):
            tau=4/math.gamma(1+1/beta)
            self.assertAlmostEqual(tau*math.gamma(1+1/beta),4.)
        self.assertAlmostEqual(float(cdf(np.array([4*math.log(2)]),1)[0]),.5)
