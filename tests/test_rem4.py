import importlib.util
import sys
from pathlib import Path
import numpy as np

MODULE_PATH = Path(__file__).parents[1] / "src" / "rem4_numerical.py"
spec = importlib.util.spec_from_file_location("rem4_numerical", MODULE_PATH)
rem4 = importlib.util.module_from_spec(spec)
assert spec.loader is not None
sys.modules[spec.name] = rem4
spec.loader.exec_module(rem4)


def test_hamiltonian_is_hermitian():
    h, bonds = rem4.xy_chain_hamiltonian(3, [1.5, 0.6], 0.2)
    assert np.allclose(h, h.conj().T)
    assert len(bonds) == 2


def test_ground_state_is_normalized():
    h, _ = rem4.xy_chain_hamiltonian(3, [1.5, 0.6], 0.2)
    _, psi = rem4.ground_state(h)
    assert np.isclose(np.linalg.norm(psi), 1.0)


def test_mutual_information_is_nonnegative():
    h, _ = rem4.xy_chain_hamiltonian(3, [1.5, 0.6], 0.2)
    _, psi = rem4.ground_state(h)
    for cut in (1, 2):
        assert rem4.mutual_information_for_cut(psi, 3, cut) >= -1e-10


def test_phase_scan_writes_output(tmp_path):
    result = rem4.make_three_qubit_phase_scan(
        outdir=str(tmp_path),
        j23=0.6,
        h=0.2,
        lambda_max=1.0,
        grid_r=8,
        grid_l=8,
        alpha=0.85,
        beta=0.35,
        state_mode="ground",
        evolve_time=0.8,
    )
    assert Path(result["phase_png"]).exists()
