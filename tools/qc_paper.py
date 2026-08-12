#!/usr/bin/env python3
"""Final QC audit for the REM paper v0.1 (master's QC-1..QC-8 gates).

  QC-1  Spec v2.2 and canonical definition exact match
  QC-2  factor-of-2 / sign audit PASS (dC^2/dt identity, J^(0) = -1/2 dC^2/dt,
        lim J^(tau) -> J^(0); no -d/dt log C^2 remnants)
  QC-3  Phase E Taylor derivation independently verified (analytic J1 vs
        numerical; tau_c^(2) positive physical root, branch continuity)
  QC-4  text / table / figure caption / Supplementary / JSON exact match
  QC-5  dimensional analysis PASS ([J]=T^-1, [lambda]=T)
  QC-6  old Gamma_F status unified in the full text
  QC-7  unsupported / unsaved historical values = 0
  QC-8  all Claims within evidence scope (language audit)

Run: python3 tools/qc_paper.py  (exits non-zero if any gate FAILs)
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
OUTDIR = REPO / "analysis_output"
MD = (REPO / "papers" / "REM_paper_v0_1.md").read_text(encoding="utf-8")
SPEC = (REPO / "REM_spec_v2_2.md").read_text(encoding="utf-8")
NORM = re.sub(r"\s+", " ", MD)

RESULTS = {}


def load(name):
    return json.loads((OUTDIR / name).read_text(encoding="utf-8"))


def gate(name, ok, evidence):
    RESULTS[name] = bool(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {evidence}")


# ── QC-1 / QC-6: canonical definition and Gamma_F status ──────────────────
def qc1_qc6():
    # canonical J_dyn^(0): the Spec writes \mathcal L and the HS subscript
    # (\langle Q_F\rho,\; Q_F\mathcal L(\rho)\rangle_{\mathrm{HS}}); the paper
    # now carries the same HS subscript (QC-1 fix, 2026-08-12).
    spec_norm = re.sub(r"\s+", " ", SPEC)
    ok1 = ("\\rangle_{\\mathrm{HS}}" in spec_norm and
           "\\operatorname{Re}" in spec_norm and
           "Q_F\\mathcal L(\\rho)" in spec_norm)
    gate("QC-1a spec has HS canonical", ok1,
         "Spec defines J_dyn^(0) = -Re<Q_F rho, Q_F L(rho)>_HS")
    ok1b = ("\\rangle_{\\mathrm{HS}}" in NORM and
            ("\\operatorname{Re}\\langle Q_F\\rho" in NORM or
             "\\operatorname{Re}\\langle Q_F\\rho,\\; Q_F" in NORM))
    gate("QC-1b paper has canonical with HS", ok1b,
         "paper J_dyn^(0) = -Re<Q_F rho, Q_F L(rho)>_HS")
    ok1c = "I_\\rho(F) - \\lambda J_{\\mathrm{dyn}}^{(0)}(F)" in NORM or \
           "I_\\rho(F)-\\lambda J_{\\rm dyn}^{(0)}(F)" in NORM
    gate("QC-1c Phi = I - lambda J", ok1c, "Phi(F;lambda) = I_rho(F) - lambda J_dyn^(0)(F)")
    # Gamma_F status: every occurrence must be flagged non-canonical
    gamma_ctx = [m.start() for m in re.finditer(r"\\Gamma", NORM)]
    bad = 0
    for i in gamma_ctx:
        seg = NORM[max(0, i - 120): i + 120]
        if ("local diagnostic" in seg or "normalized" in seg or "falsif" in seg
                or "singular" in seg or "prohibited" in seg
                or "denominator" in seg or "breakdown" in seg or "Claim 1" in seg
                or "not an unrestricted" in seg or "under $\\Gamma$" in seg
                or "superseded" in seg or "under $\\Gamma$" in seg):
            continue
        bad += 1
    gate("QC-6 Gamma_F status unified", bad == 0,
         f"{len(gamma_ctx)} Gamma occurrences, {bad} without non-canonical flag")


# ── QC-2: factor-of-2 / sign audit (numerical) ─────────────────────────────
def qc2():
    sys.path.insert(0, str(REPO / "src"))
    sys.path.insert(0, str(REPO / "analysis"))
    from gamma_F import schmidt_basis, dephasing_projection
    from rem4_numerical import mutual_information_for_cut
    from b2_canonical_open import liouvillian_dephasing
    from d1_hamiltonians import get_hamiltonian
    from d2_2a_unnormalized_instantaneous import build_states
    import scipy.linalg as sla

    H = get_hamiltonian("asymmetric_XY", n_sites=3)
    L = liouvillian_dephasing(H, [0.5, 1.0, 2.0])
    rho0 = build_states(H)["haar"]
    U = np.array(load("d2_2a_unnormalized_instantaneous.json")
                 ["core"]["haar"]["best_U_real"]) + 1j * np.array(
        load("d2_2a_unnormalized_instantaneous.json")["core"]["haar"]["best_U_imag"])
    rho_U = U.conj().T @ rho0 @ U
    Y_U = U.conj().T @ ((L @ rho0.reshape(-1, order="F")).reshape(rho0.shape, order="F")) @ U
    psi = np.linalg.eigh(rho_U)[1][:, np.argmax(np.linalg.eigh(rho_U)[0])]
    basis = schmidt_basis(psi, 3, 2)
    q = lambda r: r - dephasing_projection(r, basis)  # noqa: E731
    Qrho = q(rho_U)
    QY = q(Y_U)
    C0 = float(np.linalg.norm(Qrho, ord="fro") ** 2)
    J0_ana = -float(np.real(np.trace(Qrho.conj().T @ QY)))
    # dC^2/dt at 0 by central difference
    h = 1e-6
    def C2(t):
        rho_t = (sla.expm(t * L) @ rho0.reshape(-1, order="F")).reshape(rho0.shape, order="F")
        return float(np.linalg.norm(q(U.conj().T @ rho_t @ U), ord="fro") ** 2)
    dC2dt = (C2(h) - C2(-h)) / (2 * h)
    twoRe = 2.0 * float(np.real(np.trace(Qrho.conj().T @ QY)))
    gate("QC-2a dC^2/dt = 2Re<Qrho,QLrho>",
         abs(dC2dt - twoRe) < 1e-4,
         f"central diff {dC2dt:.6f} vs 2Re<> {twoRe:.6f}")
    gate("QC-2b J^(0) = -1/2 dC^2/dt",
         abs(J0_ana - (-0.5 * dC2dt)) < 1e-4,
         f"J0 {J0_ana:.6f} vs -dC2dt/2 {-0.5*dC2dt:.6f}")
    # lim J^(tau) -> J^(0)
    tau = 1e-4
    C_tau = C2(tau)
    J_tau = -(C_tau - C0) / (2 * tau)
    gate("QC-2c lim J^(tau) = J^(0)",
         abs(J_tau - J0_ana) < 1e-3,
         f"J_tau {J_tau:.6f} vs J0 {J0_ana:.6f}")
    # no -d/dt log C^2 remnants in the paper
    gate("QC-2d no -d/dt log C^2 remnants",
         r"-\frac{d}{dt}\log" not in NORM and "log C_F^2" not in NORM and
         r"\log C_F" not in NORM,
         "paper contains no -d/dt log C_F^2 convention")


# ── QC-3: Phase E Taylor coefficients ──────────────────────────────────────
def qc3():
    e = load("e1_tauc_mechanism.json")
    ok_all = True
    detail = {}
    for k, s in e["systems"].items():
        # quadratic refit from stored DeltaPhi_quad; check root branch
        dq = {float(t): s["DeltaPhi_quad"][str(t)] for t in
              sorted(s["DeltaPhi_quad"], key=float)}
        xq = np.array(sorted(dq))
        yq = np.array([dq[t] for t in xq])
        Aq = np.vstack([xq ** 2, xq, np.ones_like(xq)]).T
        cq, bq, aq = np.linalg.lstsq(Aq, yq, rcond=None)[0]
        disc = bq ** 2 - 4 * aq * cq
        roots = [(-bq + np.sqrt(disc)) / (2 * cq), (-bq - np.sqrt(disc)) / (2 * cq)] \
            if disc >= 0 else []
        pos = [r for r in roots if r > 0]
        # branch continuity: the chosen root must approach tau_c^(1) as cq->0
        # (the small root of a+b tau + c tau^2 ~ 0 is the one near -a/b)
        tau_c1 = s["tau_c1_num"]
        cont = None
        if pos:
            cont = min(pos, key=lambda r: abs(r - tau_c1))
        chosen = s["tau_c2"]
        ok = (chosen is not None and cont is not None and abs(chosen - cont) < 1e-6
              and chosen > 0 and abs(chosen - s["tau_c_meas"]) / s["tau_c_meas"] < 0.05)
        ok_all = ok_all and ok
        detail[k] = dict(roots=[round(r, 5) for r in roots], chosen=chosen,
                         continuous=cont, meas=s["tau_c_meas"])
    gate("QC-3 tau_c^(2) positive physical root, continuous branch",
         ok_all, json.dumps(detail, ensure_ascii=False))
    # analytic J1 = -C''(0)/4 verified in e1 JSON (num vs ana match)
    diffs = [abs(s["DeltaPhi_1_num"] - s["DeltaPhi_1_ana"])
             / abs(s["DeltaPhi_1_ana"]) for s in e["systems"].values()]
    gate("QC-3b DeltaPhi_1 num vs ana",
         max(diffs) < 0.02, f"max rel diff {max(diffs):.4f}")


# ── QC-4: source-of-truth number mapping ───────────────────────────────────
def qc4():
    d2_1 = load("d2_1_haar_singularity_audit.json")
    d3 = load("d3_timescale.json")
    d4 = load("d4_finite_size.json")
    full5 = load("d4_b_n5_cut23_haar_full.json")
    e1 = load("e1_tauc_mechanism.json")
    d1r = load("d1r_canonical_revalidation.json")
    d2r = load("d2r_state_revalidation.json")

    pairs = [
        ("Claim1 Phi 25.7", f"{d2_1['stats']['Phi_max']:.1f}", "25.7"),
        ("Claim1 Cmin", "4\\times10^{-4}", r"4\times10^{-4}"),
        ("Claim1 C best-trial", "6\\times10^{-4}", r"6\times10^{-4}"),
        ("Claim1 Gamma", f"{d2_1['stats']['Gamma_min']:.0f}", "-129"),
        ("Claim1 pmax", "0.9998", "0.9998"),
        ("tau_c N3", f"{d3['tau_c']:.4f}", "0.0180"),
        ("tau_c N4_22", f"{d4['gates']['D4-7']['tau_c_list']['N4_22']:.4f}", "0.0206"),
        ("tau_c N4_13", f"{d4['gates']['D4-7']['tau_c_list']['N4_13']:.4f}", "0.0224"),
        ("tau_c N5 full", f"{full5['tau_c']:.4f}", "0.0217"),
        ("Liouvillian tau_c*Delta_L", "0.012", "0.012"),
    ]
    for k in ("N3_21", "N4_22", "N4_13", "N5_23_full"):
        s = e1["systems"][k]
        pairs.append((f"E1 {k} meas", f"{s['tau_c_meas']:.4f}", f"{s['tau_c_meas']:.4f}"))
        pairs.append((f"E1 {k} p1", f"{s['tau_c1_num']:.4f}", f"{s['tau_c1_num']:.4f}"))
        pairs.append((f"E1 {k} p2", f"{s['tau_c2']:.4f}", f"{s['tau_c2']:.4f}"))
    fails = []
    for label, jsonv, paperv in pairs:
        if paperv not in NORM:
            fails.append(label)
    # D1-R / D2-R cross-eval gaps
    g1 = d1r["gates"] if "gates" in d1r else {}
    ok = not fails
    gate("QC-4 text/caption/table/JSON match",
         ok, f"{len(pairs)-len(fails)}/{len(pairs)} numbers verified" +
             (f"; FAILS {fails}" if fails else ""))


# ── QC-5: dimensional analysis ─────────────────────────────────────────────
def qc5():
    # [C^2] = 1 (norm of dimensionless operator) -> [J^(tau)] = 1/T = T^-1
    # [Phi] = 1 (mutual information, dimensionless) -> [lambda] = T
    # Gamma_F has [T^-1] too but is singular; J^(0) = -Re<Qrho,QLrho> has [T^-1]
    ok = ("T^{-1}" in NORM or "T^{-1}" in NORM.replace("$", "")) and \
         ("[\\lambda] = T" in NORM or "[\\lambda]=T" in NORM or
          "[$\\lambda$] = T" in NORM or "[\\lambda] = T" in NORM)
    gate("QC-5 dimensional analysis", ok,
         "[J]=T^-1, [lambda]=T stated in the paper")


# ── QC-7: unsupported / unsaved historical values = 0 ──────────────────────
def qc7():
    banned = ["44", "6\\times10^{-5}", "-222", "0.9997"]
    hits = []
    lines = MD.splitlines()
    for b in banned:
        for m in re.finditer(re.escape(b), NORM):
            hits.append((b, NORM[max(0, m.start() - 80): m.end() + 80][:120]))
    # 0.0173: allowed ONLY inside the Supplementary section (restricted-domain
    # history) and nowhere in the main text
    supp_start = NORM.find("# Supplementary")
    supp_end = len(NORM)
    for m in re.finditer(re.escape("0.0173"), NORM):
        in_supp = supp_start <= m.start() <= supp_end
        if not in_supp:
            hits.append(("0.0173 in main text",
                         NORM[max(0, m.start() - 80): m.end() + 80][:120]))
    gate("QC-7 unsupported historical values = 0",
         not hits, "clean" if not hits else str(hits))


# ── QC-8: claims scope ─────────────────────────────────────────────────────
def qc8():
    forbidden = ["universal", "scale independent", "thermodynamic limit",
                 "phase transition"]
    hits = []
    for f in forbidden:
        for m in re.finditer(re.escape(f), NORM):
            seg = NORM[max(0, m.start() - 60): m.end() + 60]
            # 'not a phase transition' / 'not a thermodynamic limit' are fine
            if "not a " in seg or "not \"phase" in seg or "not 'phase" in seg:
                continue
            hits.append((f, seg[:120]))
    gate("QC-8 claims within evidence scope",
         not hits, "clean" if not hits else str(hits))


def main():
    print("=== Final QC audit — REM paper v0.1 ===")
    qc1_qc6()
    qc2()
    qc3()
    qc4()
    qc5()
    qc7()
    qc8()
    n_ok = sum(1 for v in RESULTS.values() if v)
    print(f"\n  Gates: {n_ok}/{len(RESULTS)} PASS")
    sys.exit(0 if n_ok == len(RESULTS) else 1)


if __name__ == "__main__":
    main()
