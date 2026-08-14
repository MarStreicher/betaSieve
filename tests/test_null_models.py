import numpy as np
import pandas as pd
import pytest
from scipy import stats

from betasieve.analysis import Col, _add_statistics
from betasieve.cg_probe_table import DesignGroup
from betasieve.null_models import (
    BetaBinomialNull,
    BinomialNull,
    EmpiricalNull,
    NullModels,
)


def test_binomial_sf_is_inclusive_upper_tail() -> None:
    model = BinomialNull(n=10, p0=0.2)

    np.testing.assert_allclose(model.sf([0, 3]), stats.binom.sf([-1, 2], 10, 0.2))
    np.testing.assert_allclose(model.pmf(np.arange(11)).sum(), 1.0)


def test_binomial_critical_p_hat_is_first_attainable_value_beyond_alpha() -> None:
    model = BinomialNull(n=10, p0=0.2)

    critical = model.critical_p_hat(0.05)

    assert critical == pytest.approx(0.5)
    assert model.sf(critical * model.n) < 0.05
    assert model.sf(critical * model.n - 1) >= 0.05


def test_critical_p_hat_is_none_when_alpha_is_unreachable() -> None:
    assert BinomialNull(n=4, p0=0.5).critical_p_hat(1e-6) is None


def test_beta_binomial_moments_reproduce_the_p_hat_mean_and_variance() -> None:
    reference = np.array([0.0, 0.05, 0.1, 0.2, 0.35])
    p_mean, variance = reference.mean(), reference.var()

    model = BetaBinomialNull.fit_alpha_beta(20, reference)

    assert model.p0 == pytest.approx(p_mean)
    assert model.a + model.b == pytest.approx(p_mean * (1 - p_mean) / variance - 1)
    assert model.a == pytest.approx(p_mean * (model.a + model.b))
    assert model.b == pytest.approx((1 - p_mean) * (model.a + model.b))


def test_empirical_sf_and_cutoff_agree() -> None:
    model = EmpiricalNull.fit(n=4, p_hat_reference=[0.0, 0.0, 0.25, 0.5])

    np.testing.assert_allclose(model.sf([0.0, 0.25, 0.75]), [1.0, 0.6, 0.2])

    critical = model.critical_p_hat(0.25)
    assert critical == pytest.approx(0.75)
    assert model.sf(critical) < 0.25


def test_null_models_round_trip_to_dict() -> None:
    models = NullModels.fit(
        n=8,
        p0=0.25,
        p_hat_reference=[0.0, 0.125, 0.25, 0.5],
        threshold=0.1,
        confidence=0.95,
    )

    payload = models.to_dict()

    assert models.alpha == pytest.approx(0.05)
    assert models.n == 8
    assert payload["binom"] == {"n": 8, "p0": 0.25}
    assert payload["empir"] == {"n": 8, "m": 4}
    assert set(payload["bb"]) == {"n", "a", "b", "p0"}


def test_add_statistics_columns_match_the_returned_models(
    diff_frame: pd.DataFrame,
) -> None:
    frame, models = _add_statistics(diff_frame, 0.1, "fdr_bh", 0.95)

    reference = frame.loc[
        frame[Col.GROUP] == DesignGroup.EXACT_REPLICATES.value, Col.P_HAT
    ]
    assert models.threshold == 0.1
    assert models.binomial.p0 == pytest.approx(frame[Col.P0].iloc[0])
    assert models.binomial.p0 == pytest.approx(reference.mean())

    np.testing.assert_allclose(
        frame[Col.P_EMPIR].to_numpy(),
        models.empirical.sf(frame[Col.P_HAT].to_numpy()),
    )
    np.testing.assert_allclose(
        frame[Col.P_BETA].to_numpy(),
        models.beta_binomial.sf(frame[Col.ABOVE].to_numpy()),
    )
