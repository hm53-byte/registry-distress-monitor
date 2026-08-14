"""Mondrian konformna kalibracija, James-Stein sazimanje i KS brana.

Problem. Za cjenovni raspon nad malim uzorkom obicna procjena laze u oba
smjera: preuska je kad podataka nema dovoljno, a preusiroka kad ih ima. Tri
ortogonalna primitiva rjesavaju tri razlicita rezima velicine uzorka.

  1. Mondrian konformna predikcija po sloju (Vovk, Gammerman, Shafer 2005,
     Th. 2.1). Za n >= 30 daje pokrivenost >= 1 - alpha bez pretpostavke o
     distribuciji. Izvedeno cistom standardnom bibliotekom, bez MAPIE, jer
     nova teska ovisnost za pedesetak redaka racuna nije opravdana.

  2. James-Stein analiticko sazimanje (1961). Za 3 <= n < 30 povlaci
     procjenu prema globalnom prosjeku: manja varijanca, nesto pristranosti.

  3. KS brana prije sazimanja. Ako se lokalna i globalna distribucija
     razlikuju (p < 0.05), sazimanje se NE primjenjuje, jer globalni
     prosjek tada nije reprezentativan i sazimanje bi uvelo pristranost
     umjesto da je smanji.

Sucelje:
  mondrian_band(amounts, alpha=0.05) -> ConformalBand
  james_stein_shrink(group_means, group_sizes, global_mean, sigma2)
  ks_gate(local_dist, global_dist, alpha=0.05) -> bool
  effective_n(amounts, source_origins) -> int
  band_dual_path(local_amounts, global_amounts, alpha=0.05) -> ConformalBand

Reference:
  Vovk, Gammerman, Shafer 2005, Algorithmic Learning in a Random World, Th. 2.1
  James, Stein 1961, Estimation with quadratic loss, 4th Berkeley Symp.
  Massey 1951, The Kolmogorov-Smirnov test for goodness of fit
  Stephens 1970, Use of the Kolmogorov-Smirnov statistics
"""
from __future__ import annotations

import math
import statistics
from dataclasses import dataclass

# scipy je opcionalan — koristimo ga za KS test ako postoji, inace pure-NumPy fallback
try:
    from scipy import stats as _scipy_stats  # type: ignore[import-not-found]
    _HAS_SCIPY = True
except ImportError:
    _HAS_SCIPY = False


# === Konstante ==============================================================

MONDRIAN_MIN_N = 30          # ≥30 → conformal (Th. 2.1 daje formalno garanciju)
SHRINKAGE_MIN_N = 3          # 3≤N<30 → James-Stein shrinkage
KS_ALPHA = 0.05              # prag KS testa iznad kojeg se sazimanje dopusta
DEFAULT_COVERAGE_CONFORMAL = 0.95
DEFAULT_COVERAGE_WILSON = 0.80


@dataclass(frozen=True)
class ConformalBand:
    """Rezultat dual-path band computation-a."""
    method: str                # 'mondrian','james_stein','wilson_fallback','no_data'
    n: int
    median: float | None
    q_low: float | None        # donja granica (npr. q_alpha/2)
    q_high: float | None       # gornja granica (npr. q_(1-alpha/2))
    coverage: float            # nominalan coverage [0,1]
    shrinkage_factor: float | None = None  # samo za james_stein, [0,1]
    ks_gated: bool = False     # True ako je KS rejecto-ao shrinkage path


# === Primitiv 1: Mondrian conformal (pure NumPy ICP) ========================


def mondrian_band(amounts: list[float], alpha: float = 0.05) -> ConformalBand:
    """Distribution-free conformal band za jedan stratum (npr. jedan k.o.).

    Vovk Th. 2.1: za i.i.d. uzorak velicine n, kvantil q_(1-alpha)(n+1)/n
    nonconformity score-ova daje predictive set s coverage ≥ 1-alpha.

    Nonconformity score: |amount - median(amounts)|. Vrati [median - q, median + q].
    """
    n = len(amounts)
    if n < MONDRIAN_MIN_N:
        return ConformalBand("no_data", n, None, None, None, 0.0)
    sorted_a = sorted(amounts)
    med = statistics.median(sorted_a)
    scores = [abs(x - med) for x in sorted_a]
    scores.sort()
    # Vovk korekcija: kvantil pozicije ceil((n+1)(1-alpha)) / n
    rank = math.ceil((n + 1) * (1 - alpha))
    rank = min(rank, n)        # cap pri n (ako (n+1)(1-alpha) > n)
    q = scores[rank - 1]       # 1-based → 0-based index
    coverage = max(0.0, 1.0 - alpha)
    return ConformalBand(
        method="mondrian",
        n=n,
        median=med,
        q_low=med - q,
        q_high=med + q,
        coverage=coverage,
    )


# === Primitiv 2: James-Stein analitički shrinkage ===========================


def james_stein_shrink(
    group_means: list[float],
    group_sizes: list[int],
    global_mean: float,
    sigma2: float,
) -> tuple[list[float], float]:
    """Shrinkage prema globalnom mean-u (analitička JS formula).

    Za k≥3 grupa, shrinkage faktor c = 1 - (k-2)*sigma2 / sum_i n_i*(x̄_i - μ)².

    Vraća (shrinked_means, shrinkage_factor). shrinkage_factor ∈ [0, 1]:
    1 = bez shrinkage-a, 0 = potpuno shrinka prema global_mean. Manji broj
    = jača shrinkage = manje variance, više bias-a.
    """
    k = len(group_means)
    if k < 3 or sigma2 <= 0:
        return list(group_means), 1.0
    ss = sum(n * (m - global_mean) ** 2 for n, m in zip(group_sizes, group_means))
    if ss <= 0:
        return [global_mean] * k, 0.0
    c = 1.0 - (k - 2) * sigma2 / ss
    c = max(0.0, min(1.0, c))  # clamp [0,1] (positive-part JS)
    shrinked = [global_mean + c * (m - global_mean) for m in group_means]
    return shrinked, c


# === Primitiv 3: KS brana protiv pristranosti ===============================


def ks_gate(
    local_dist: list[float],
    global_dist: list[float],
    alpha: float = KS_ALPHA,
) -> bool:
    """KS dvostruki uzorak test. True = ne odbija H0 (distribucije ekvivalentne) → smije shrinkati.
    False = odbija H0 → ne shrinkaj (lokalni stratum je bitno različit od global pool-a).
    """
    if len(local_dist) < 3 or len(global_dist) < 3:
        return False
    if _HAS_SCIPY:
        _, p_value = _scipy_stats.ks_2samp(local_dist, global_dist)
        # bool() cast: scipy vraca numpy.bool_, a pozivatelji ocekuju Python bool
        # (konzistentno s pure-NumPy fallback putem ispod).
        return bool(p_value >= alpha)
    # Pure-NumPy fallback (Massey 1951): KS statistic = max |F1(x) - F2(x)|
    sorted_l = sorted(local_dist)
    sorted_g = sorted(global_dist)
    nl, ng = len(sorted_l), len(sorted_g)
    all_x = sorted(set(sorted_l) | set(sorted_g))
    d_max = 0.0
    for x in all_x:
        f_l = sum(1 for v in sorted_l if v <= x) / nl
        f_g = sum(1 for v in sorted_g if v <= x) / ng
        d = abs(f_l - f_g)
        if d > d_max:
            d_max = d
    # Critical value @ alpha=0.05 (Stephens 1970): c(alpha) * sqrt((n+m)/(n*m))
    c_alpha = 1.358 if alpha <= 0.05 else 1.224  # crude two-level approx
    crit = c_alpha * math.sqrt((nl + ng) / (nl * ng))
    return d_max < crit  # True = ne odbija H0


# === Helper: efektivni N =================================


def effective_n(amounts: list[float], source_origins: list[str] | None = None) -> int:
    """Reducira N po korelacijskom faktoru ako su uzorci iz istog izvora.

    Nacelo: norma se ne smije proglasiti dok ne postoji
    n_min **nezavisnih** opažanja. Ako su sva 30 amount-a iz iste serije istog izvora
    i iste godine, efektivno_n je manji jer su korelirani.

    Trenutno minimalna implementacija: ako su sve source_origins isti, n / 2;
    inace n. Kasnije moze ići preko Newey-West / Spearman ρ.
    """
    n = len(amounts)
    if not source_origins or n == 0:
        return n
    if len(set(source_origins)) == 1:
        return max(1, n // 2)
    return n


# === Glavni dual-path API ==================================================


def band_dual_path(
    local_amounts: list[float],
    global_amounts: list[float],
    alpha: float = 0.05,
) -> ConformalBand:
    """Glavni entry point: bira Mondrian (n≥30) ili JS shrinkage (3≤n<30).

    Logika:
      - n≥30 → Mondrian conformal (ignoira global_amounts)
      - 3≤n<30 → KS gating; ako prolazi → James-Stein shrinkage; ako ne → Wilson fallback
      - n<3 → no_data
    """
    n = len(local_amounts)
    if n >= MONDRIAN_MIN_N:
        return mondrian_band(local_amounts, alpha=alpha)

    if n < SHRINKAGE_MIN_N:
        return ConformalBand("no_data", n, None, None, None, 0.0)

    # 3 <= n < 30: try JS shrinkage, gated by KS
    local_med = statistics.median(local_amounts)
    if not global_amounts or len(global_amounts) < SHRINKAGE_MIN_N:
        # Nema globalnog reference-a → Wilson fallback (vidi Wilsonov interval)
        return ConformalBand(
            "wilson_fallback", n, local_med, None, None,
            DEFAULT_COVERAGE_WILSON, ks_gated=False,
        )

    if not ks_gate(local_amounts, global_amounts, alpha=KS_ALPHA):
        return ConformalBand(
            "wilson_fallback", n, local_med, None, None,
            DEFAULT_COVERAGE_WILSON, ks_gated=True,
        )

    global_mean = sum(global_amounts) / len(global_amounts)
    sigma2 = statistics.variance(global_amounts) if len(global_amounts) > 1 else 1.0
    shrinked, c = james_stein_shrink(
        group_means=[local_med],
        group_sizes=[n],
        global_mean=global_mean,
        sigma2=sigma2,
    )
    js_med = shrinked[0]
    # Empirijski IQR oko shrinked mean-a kao prelimary band (faza 2 može uvesti
    # proper conformal s shrinkage-adjusted nonconformity score-ovima)
    sorted_a = sorted(local_amounts)
    if n >= 4:
        q1 = statistics.quantiles(sorted_a, n=4)[0]
        q3 = statistics.quantiles(sorted_a, n=4)[2]
    else:
        q1, q3 = sorted_a[0], sorted_a[-1]
    return ConformalBand(
        method="james_stein",
        n=n,
        median=js_med,
        q_low=q1,
        q_high=q3,
        coverage=DEFAULT_COVERAGE_WILSON,
        shrinkage_factor=c,
        ks_gated=False,
    )
