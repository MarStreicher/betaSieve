from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Optional

import numpy as np
from scipy import stats


def _critical_p_hat(ppf_value: float, n: int) -> Optional[float]:
    k = int(ppf_value) + 1
    return k / n if k <= n else None


@dataclass(frozen=True)
class BinomialNull:
    n: int
    p0: float

    def pmf(self, k: Any) -> np.ndarray:
        return stats.binom.pmf(k, self.n, self.p0)

    def sf(self, k: Any) -> np.ndarray:
        # P(K>=k) (continuous)
        return stats.binom.sf(np.asarray(k) - 1, self.n, self.p0)

    def critical_p_hat(self, alpha: float) -> Optional[float]:
        return _critical_p_hat(stats.binom.ppf(1.0 - alpha, self.n, self.p0), self.n)

    def to_dict(self) -> Dict[str, float]:
        return {"n": self.n, "p0": self.p0}


@dataclass(frozen=True)
class BetaBinomialNull:
    n: int
    a: float
    b: float

    @classmethod
    def fit_alpha_beta(cls, n: int, p_hat_reference: Any) -> "BetaBinomialNull":
        values = np.asarray(p_hat_reference, dtype=float)
        p_mean = float(values.mean()) if values.size else 0.0
        p_variance = np.sum((p_hat_reference - p_mean) ** 2) / p_hat_reference.size()

        p_t = (p_mean * (1 - p_mean) / p_variance) - 1
        p_alpha = p_mean * p_t
        p_beta = (1 - p_mean) * p_t
        # expected_successes = n_samples * p_alpha / (p_alpha + p_beta)

        return cls(
            n=int(n),
            a=p_alpha,
            b=p_beta,
        )

    @property
    def p0(self) -> float:
        return self.a / (self.a + self.b)

    def _dist(self):
        return stats.betabinom(self.n, self.a, self.b)

    def pmf(self, k: Any) -> np.ndarray:
        # P(K=k)
        return self._dist().pmf(k)

    def sf(self, k: Any) -> np.ndarray:
        # P(K>k) bzw. P(K>=k) == P(K>(k-1)) (discrete)
        return self._dist().sf(np.asarray(k) - 1)

    def critical_p_hat(self, alpha: float) -> Optional[float]:
        return _critical_p_hat(self._dist().ppf(1.0 - alpha), self.n)

    def to_dict(self) -> Dict[str, float]:
        return {"n": self.n, "a": self.a, "b": self.b, "p0": self.p0}


@dataclass(frozen=True, eq=False)
class EmpiricalNull:
    n: int
    reference_p_hat: np.ndarray

    @classmethod
    def fit(cls, n: int, p_hat_reference: Any) -> "EmpiricalNull":
        return cls(
            n=int(n),
            reference_p_hat=np.sort(np.asarray(p_hat_reference, dtype=float)),
        )

    @property
    def m(self) -> int:
        return int(self.reference_p_hat.size)

    def sf(self, p_hat: Any) -> np.ndarray:
        counts_ge = self.m - np.searchsorted(
            self.reference_p_hat, np.asarray(p_hat), side="left"
        )
        return (1 + counts_ge) / (1 + self.m)

    def critical_p_hat(self, alpha: float) -> Optional[float]:
        attainable = np.arange(self.n + 1) / self.n
        beyond = np.flatnonzero(self.sf(attainable) < alpha)
        return float(attainable[beyond[0]]) if beyond.size else None

    def to_dict(self) -> Dict[str, float]:
        return {"n": self.n, "m": self.m}


@dataclass(frozen=True, eq=False)
class NullModels:
    threshold: float
    confidence: float
    binomial: BinomialNull
    beta_binomial: BetaBinomialNull
    empirical: EmpiricalNull

    @classmethod
    def fit(
        cls,
        *,
        n: int,
        p0: float,
        p_hat_reference: Any,
        threshold: float,
        confidence: float,
    ) -> "NullModels":
        return cls(
            threshold=float(threshold),
            confidence=float(confidence),
            binomial=BinomialNull(n=int(n), p0=float(p0)),
            beta_binomial=BetaBinomialNull.fit_alpha_beta(n, p_hat_reference),
            empirical=EmpiricalNull.fit(n, p_hat_reference),
        )

    @property
    def alpha(self) -> float:
        return 1.0 - self.confidence

    @property
    def n(self) -> int:
        return self.binomial.n

    def to_dict(self) -> Dict[str, Any]:
        return {
            "threshold": self.threshold,
            "confidence": self.confidence,
            "binom": self.binomial.to_dict(),
            "bb": self.beta_binomial.to_dict(),
            "empir": self.empirical.to_dict(),
        }


__all__ = [
    "BetaBinomialNull",
    "BinomialNull",
    "EmpiricalNull",
    "NullModels",
]
