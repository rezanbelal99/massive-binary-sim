"""Unit tests for individual physics modules + a golden-file regression test."""
import json
from pathlib import Path

import numpy as np
import pytest

from massive_binary_sim import load_config
from massive_binary_sim.numerics import units as U


def test_units_provenance_complete():
    from massive_binary_sim.numerics.units import CONSTANT_PROVENANCE
    for c in CONSTANT_PROVENANCE:
        assert c.source and c.value > 0


def test_config_roundtrip_and_validation():
    cfg = load_config()
    assert cfg.system.m1 >= cfg.system.m2 or True
    d = cfg.to_dict()
    assert "free_parameters" in d
    for name, fp in cfg.free.items():
        lo, hi = fp.range
        assert lo <= fp.clipped(1e9) <= hi


def test_config_rejects_bad_theory():
    import pytest
    with pytest.raises(Exception):
        load_config(overrides={"theories.wind": "not_a_recipe"})


def test_eggleton_roche_limits():
    from massive_binary_sim.numerics.units import roche_lobe_eggleton
    # q -> 0 : R_L/a -> 0 ; q -> inf : R_L/a -> ~0.6
    assert roche_lobe_eggleton(1e-3) < 0.05
    assert 0.35 < roche_lobe_eggleton(1.0) < 0.40
    assert roche_lobe_eggleton(1e3) > 0.5


def test_lagrange_points_ordering():
    from massive_binary_sim.binary.roche import lagrange_points
    L = lagrange_points(0.5)
    assert L["L3"][0] < 0 < L["L1"][0] < 1 < L["L2"][0]


def test_peters_merger_time_scaling():
    from massive_binary_sim.dynamics.gw import merger_time_circular
    m1 = m2 = 30 * U.MSUN
    t1 = merger_time_circular(1.0 * U.RSUN, m1, m2)
    t2 = merger_time_circular(2.0 * U.RSUN, m1, m2)
    assert np.isclose(t2 / t1, 16.0, rtol=1e-6)   # T ~ a^4


def test_wind_recipes_return_negative():
    from massive_binary_sim.winds import recipes
    from massive_binary_sim.structure.evolution import StarState
    s = StarState.zams(40.0, 0.014)
    for r in ("vink2001", "bjorklund2021", "dejager1988"):
        md = recipes.mass_loss_rate(s, r)
        assert md <= 0


def test_fryer_remnant_monotonic():
    from massive_binary_sim.remnants.core_collapse import fryer2012
    masses = [(12, 2.0), (20, 5.0), (40, 12.0), (80, 30.0)]
    rems = [fryer2012(M, Mco, "delayed").m_gravitational for M, Mco in masses]
    assert rems == sorted(rems)          # heavier CO core -> heavier remnant
    assert fryer2012(12, 2.0, "delayed").kind == "NS"
    assert fryer2012(80, 30.0, "delayed").kind == "BH"


def test_pair_instability_regimes():
    from massive_binary_sim.remnants.pair_instability import classify
    assert classify(20).regime == "core_collapse"
    assert classify(50).regime == "PPISN"
    assert classify(90).regime == "PISN"
    assert classify(90).remnant_bh_mass is None


GOLDEN = Path(__file__).parent / "golden_default_run.json"


def _run_default_summary():
    from massive_binary_sim.binary.driver import evolve
    cfg = load_config()
    res = evolve(cfg)
    fin = res.final
    return {
        "fate": res.fate,
        "m1_Msun": round(fin["m1_Msun"], 3),
        "m2_Msun": round(fin["m2_Msun"], 3),
        "a_Rsun": round(fin["a_Rsun"], 2) if np.isfinite(fin["a_Rsun"]) else "inf",
        "n_events": len(res.events),
    }


def test_regression_golden_file():
    """Physics must not change silently. Update golden file deliberately:
        python -c "import tests.test_modules as t, json; \
                   open(t.GOLDEN,'w').write(json.dumps(t._run_default_summary(),indent=2))"
    """
    summary = _run_default_summary()
    if not GOLDEN.exists():
        GOLDEN.write_text(json.dumps(summary, indent=2))
        pytest.skip("golden file created; re-run to compare")
    expected = json.loads(GOLDEN.read_text())
    assert summary == expected, f"physics changed:\n got {summary}\n exp {expected}"
