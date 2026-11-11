import argparse
import numpy as np
import pytest
from swmtrack import cli
from swmtrack.preprocessing import PreprocessConfig
from swmtrack.open_preprocessing import smoothing_mass_matrix,mesh_operators,robust_unit_costs


def test_recommended_profile_preserves_explicit_baseline_and_overrides():
    defaults=PreprocessConfig()
    assert (defaults.curvature_method,defaults.smoothing_mass,defaults.segmentation_energy)==('cotan','lumped','area-weighted')
    candidate=PreprocessConfig.agreement_candidate()
    candidate.validate()
    assert (candidate.smoothing_mass,candidate.segmentation_energy)==('consistent','robust-unit')
    assert candidate.graph_regularization==1
    assert (candidate.implicit_grid_spacing_mm,candidate.implicit_derivative_sigma_mm,candidate.implicit_band_mm)==(1,1.5,6)
    assert (candidate.skeleton_method,candidate.skeleton_relabeling)==('hj-flux','pruned-branches')
    assert candidate.hj_endpoint_flux_percentile==defaults.hj_endpoint_flux_percentile==5
    assert candidate.gyral_pruning_mm==defaults.gyral_pruning_mm
    assert candidate.sulcal_pruning_mm==defaults.sulcal_pruning_mm
    assert PreprocessConfig.agreement_candidate(implicit_derivative_sigma_mm=1.25).implicit_derivative_sigma_mm==1.25


@pytest.mark.parametrize('option,value',[('smoothing_mass','invalid'),('segmentation_energy','invalid'),('backend','legacy-native')])
def test_candidate_rejects_invalid_or_native_only_combinations(option,value):
    with pytest.raises(ValueError):PreprocessConfig.agreement_candidate(**{option:value}).validate()


def test_consistent_mass_matches_analytic_triangle_and_integrates_constants():
    p=np.array([[0.,0.,0.],[2.,0.,0.],[0.,1.,0.]])
    f=np.array([[0,1,2]])
    _,area,_=mesh_operators(p,f)
    actual=smoothing_mass_matrix(p,f,area,'consistent')
    np.testing.assert_allclose(actual@np.ones(3),area,atol=1e-14)
    np.testing.assert_allclose(actual.toarray(),np.array([[2,1,1],[1,2,1],[1,1,2]])/12.)
    assert np.linalg.eigvalsh(actual.toarray()).min()>0


def test_candidate_energy_matches_gaussian_likelihood_and_ignores_unknown_fit():
    rng=np.random.default_rng(82)
    feature=np.r_[rng.normal(-.3,.1,80),rng.normal(.2,.08,80),[-.8,.9,0.]]
    area=rng.uniform(.2,2,len(feature));valid=np.arange(len(feature))<160
    actual,fit=robust_unit_costs(feature,area,valid)
    x=np.clip(feature,*fit['winsorization_bounds'])
    means=np.asarray(fit['means']);variances=np.asarray(fit['variances'])
    priors=np.asarray(fit['weights'])
    expected=.5*np.log(2*np.pi*variances)-np.log(priors)+.5*(x[:,None]-means)**2/variances
    expected=np.minimum(expected-expected.min(axis=1,keepdims=True),20)
    np.testing.assert_allclose(actual,expected,rtol=1e-12,atol=1e-12)
    moved=feature.copy();moved[~valid]=[-30,20,2]
    changed,new_fit=robust_unit_costs(moved,area,valid)
    for key in ('means','variances','weights'):assert fit[key]==new_fit[key]
    np.testing.assert_array_equal(actual[valid],changed[valid])
    assert np.any(actual[~valid]>0)
    # Exact cortical weighted-cost / area recovery can exceed 20 by roundoff.
    assert np.isfinite(actual).all() and np.all((actual>=0)&(actual<=20+1e-12))


@pytest.mark.parametrize('ordering',[0,1])
def test_cli_profile_options_override_candidate_in_either_order(tmp_path,monkeypatch,ordering):
    captured={}
    def preprocess(*args,**kwargs):captured.update(kwargs);return {'status':'test'}
    monkeypatch.setattr(cli,'preprocess',preprocess)
    profile=['--preprocess-profile','agreement-candidate']
    override=['--implicit-derivative-sigma-mm','1.25','--graph-regularization','.7']
    cli.main(['preprocess','--surface','unused','--output-dir',str(tmp_path)]+(profile+override if ordering==0 else override+profile))
    config=captured['config'];config.validate()
    assert config.implicit_derivative_sigma_mm==1.25 and config.graph_regularization==.7
    assert config.smoothing_mass=='consistent' and config.skeleton_relabeling=='pruned-branches'


def test_run_command_uses_candidate_profile(monkeypatch,tmp_path):
    captured={}
    def run(*args,**kwargs):captured.update(kwargs);return {'status':'test','hemispheres':{}}
    monkeypatch.setattr(cli,'run_pipeline',run)
    cli.main(['run','--subject-dir','unused','--fod','unused','--output-dir',str(tmp_path),
              '--sh-basis','legacy-real-even','--assume-scanner-aligned','--preprocess-profile','agreement-candidate'])
    assert captured['preprocess_config']==PreprocessConfig.agreement_candidate()


def test_cli_defaults_to_recommended_and_previous_profile_remains_explicit():
    parser=argparse.ArgumentParser();cli.add_preprocessing_options(parser)
    assert cli.preprocessing_config(parser.parse_args([]))==PreprocessConfig.recommended()
    assert cli.preprocessing_config(parser.parse_args(['--preprocess-profile','previous-defaults']))==PreprocessConfig()
