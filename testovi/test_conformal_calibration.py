"""Testovi konformne kalibracije.

Dio testova provjerava matematicko svojstvo (pokrivenost), a dio pribija
zateceno ponasanje ukljucujuci jedno mjesto na kojem sucelje obecava vise
nego sto izvedba radi. Vidi test_sazimanje_u_dual_path_ne_saziva.
"""

from __future__ import annotations

import random
import statistics
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from primjer.conformal_calibration import (  # noqa: E402
    MONDRIAN_MIN_N,
    band_dual_path,
    effective_n,
    james_stein_shrink,
    ks_gate,
    mondrian_band,
)


# --- Mondrian ---------------------------------------------------------------


def test_ispod_praga_nema_pojasa():
    """Ispod 30 opazanja teorem ne daje jamstvo, pa se pojas ne izmislja."""
    b = mondrian_band([1.0] * (MONDRIAN_MIN_N - 1))
    assert b.method == "no_data"
    assert b.q_low is None and b.q_high is None
    assert b.coverage == 0.0


def test_pojas_je_simetrican_oko_medijana():
    amounts = [float(i) for i in range(MONDRIAN_MIN_N)]
    b = mondrian_band(amounts)
    assert b.method == "mondrian"
    assert b.median == statistics.median(amounts)
    assert b.q_high - b.median == pytest.approx(b.median - b.q_low)


def test_pokrivenost_je_barem_nominalna():
    """Glavno svojstvo. Za i.i.d. uzorak udio opazanja unutar pojasa mora
    biti barem 1 - alpha. Provjerava se na sto uzoraka."""
    rnd = random.Random(20260815)
    ispod = 0
    for _ in range(100):
        uzorak = [rnd.gauss(1000, 250) for _ in range(60)]
        b = mondrian_band(uzorak, alpha=0.10)
        unutra = sum(1 for x in uzorak if b.q_low <= x <= b.q_high) / len(uzorak)
        if unutra < 0.90:
            ispod += 1
    # Jamstvo je nad ponovljenim uzorkovanjem, pa poneki uzorak smije pasti
    # ispod. Sustavno probijanje ne smije.
    assert ispod <= 10, f"{ispod} od 100 uzoraka ispod nominalne pokrivenosti"


def test_uzi_alpha_daje_siri_pojas():
    uzorak = [float(i) for i in range(100)]
    uski = mondrian_band(uzorak, alpha=0.50)
    siroki = mondrian_band(uzorak, alpha=0.01)
    assert (siroki.q_high - siroki.q_low) > (uski.q_high - uski.q_low)


def test_konstantan_uzorak_daje_pojas_nula():
    b = mondrian_band([500.0] * 40)
    assert b.q_low == b.q_high == 500.0


def test_rang_ne_izlazi_izvan_uzorka():
    """Vovkova korekcija trazi kvantil na poziciji ceil((n+1)(1-alpha)),
    sto za male alpha prelazi n. Bez ogranicenja bi indeks izasao iz polja."""
    b = mondrian_band([float(i) for i in range(30)], alpha=0.001)
    assert b.q_high is not None


# --- James-Stein ------------------------------------------------------------


def test_manje_od_tri_grupe_ne_saziva():
    """Formula trazi barem tri grupe. Ispod toga vraca ulaz nepromijenjen."""
    for k in (1, 2):
        ulaz = [100.0] * k
        izlaz, c = james_stein_shrink(ulaz, [10] * k, 500.0, 25.0)
        assert izlaz == ulaz
        assert c == 1.0


def test_sazimanje_povlaci_prema_globalnom_prosjeku():
    prosjeci = [10.0, 20.0, 30.0]
    izlaz, c = james_stein_shrink(prosjeci, [5, 5, 5], 20.0, 500.0)
    assert 0.0 <= c < 1.0
    # Svaka grupa mora se pomaknuti prema 20, nijedna preko njega.
    assert izlaz[0] > 10.0 and izlaz[0] < 20.0
    assert izlaz[2] < 30.0 and izlaz[2] > 20.0
    assert izlaz[1] == pytest.approx(20.0)


def test_identicne_grupe_potpuno_sazimaju():
    izlaz, c = james_stein_shrink([7.0, 7.0, 7.0], [3, 3, 3], 7.0, 1.0)
    assert c == 0.0
    assert izlaz == [7.0, 7.0, 7.0]


def test_faktor_je_ogranicen_na_jedinicni_raspon():
    """Pozitivni dio James-Steina: negativan faktor bi prebacio procjenu na
    drugu stranu globalnog prosjeka, sto nema smisla."""
    _, c = james_stein_shrink([1.0, 2.0, 3.0], [1, 1, 1], 2.0, 10_000.0)
    assert c == 0.0
    _, c2 = james_stein_shrink([1.0, 500.0, 1000.0], [50, 50, 50], 500.0, 0.001)
    assert c2 == pytest.approx(1.0)


# --- KS brana ---------------------------------------------------------------


def test_iste_distribucije_prolaze_branu():
    rnd = random.Random(11)
    a = [rnd.gauss(100, 10) for _ in range(60)]
    b = [rnd.gauss(100, 10) for _ in range(60)]
    assert ks_gate(a, b) is True


def test_razlicite_distribucije_padaju_na_brani():
    a = [float(x) for x in range(0, 50)]
    b = [float(x) for x in range(1000, 1050)]
    assert ks_gate(a, b) is False


def test_premali_uzorak_ne_prolazi():
    """Brana pusta sazimanje samo kad ima cime obraniti odluku."""
    assert ks_gate([1.0, 2.0], [1.0, 2.0, 3.0, 4.0]) is False
    assert ks_gate([1.0, 2.0, 3.0], []) is False


def test_brana_vraca_python_bool():
    """Kad je scipy prisutan, vraca numpy.bool_. Pozivatelji ocekuju bool."""
    a = [float(x) for x in range(10)]
    assert type(ks_gate(a, a)) is bool


# --- efektivni N ------------------------------------------------------------


def test_isti_izvor_polovi_efektivni_n():
    amounts = [1.0] * 30
    assert effective_n(amounts, ["fina"] * 30) == 15
    assert effective_n(amounts, ["fina"] * 15 + ["sudreg"] * 15) == 30
    assert effective_n(amounts, None) == 30


def test_efektivni_n_nikad_nije_nula_za_neprazan_uzorak():
    assert effective_n([1.0], ["a"]) == 1


# --- dual path --------------------------------------------------------------


def test_velik_uzorak_ide_na_mondriana():
    b = band_dual_path([float(i) for i in range(40)], [1.0, 2.0, 3.0])
    assert b.method == "mondrian"


def test_premalen_uzorak_nema_podatka():
    assert band_dual_path([1.0, 2.0], [1.0] * 50).method == "no_data"


def test_bez_globalnog_referenca_pada_na_wilsona():
    b = band_dual_path([1.0, 2.0, 3.0, 4.0], [])
    assert b.method == "wilson_fallback"
    assert b.ks_gated is False


def test_odbijena_brana_se_biljezi():
    """Kad KS odbije, razlog mora ostati vidljiv u rezultatu, inace se
    poslije ne moze razlikovati 'nije bilo globalnog uzorka' od
    'globalni uzorak nije bio usporediv'."""
    lokalno = [100.0, 105.0, 110.0, 115.0, 120.0]
    globalno = [10_000.0 + i for i in range(50)]
    b = band_dual_path(lokalno, globalno)
    assert b.method == "wilson_fallback"
    assert b.ks_gated is True


def test_sazimanje_u_dual_path_ne_saziva():
    """Zateceno ponasanje, ne zeljeno.

    band_dual_path predaje James-Steinu jednu jedinu grupu. Formula trazi
    barem tri, pa vraca ulaz nepromijenjen i faktor 1,0. Rezultat se oznaci
    kao 'james_stein', ali sazimanja nema: vraceni medijan jednak je
    lokalnom medijanu, a globalni prosjek ne utjece ni na sto.

    Test postoji da promjena tog ponasanja bude vidljiva, a ne tiha.
    """
    rnd = random.Random(7)
    globalno = [rnd.gauss(100, 15) for _ in range(200)]
    lokalno = [rnd.gauss(100, 15) for _ in range(10)]
    b = band_dual_path(lokalno, globalno)
    assert b.method == "james_stein"
    assert b.shrinkage_factor == 1.0
    assert b.median == pytest.approx(statistics.median(lokalno))
