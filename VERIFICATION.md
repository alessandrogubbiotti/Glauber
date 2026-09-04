# Verification status — scaling-limit claims of `glauber_paper.tex`

**Motivation.** The scaling-limit section of the paper makes a sequence of
quantitative claims (Theorem `thm:pfaffian` with the complex hierarchy kernel
`e:Kpf`, the boxed density `e:fmain`, the finite-N Pfaffian `e:multipf` and its
centered form `c:centered`, the regime family of `ss:numregimes`). Several are
backed by numerics that live in *this* repository. This document makes the
correspondence "paper claim ⇄ script ⇄ measured number" **explicit**: for every
claim, what verifies it, what the measured outcome was, and — for the open
items — *why* the check still matters and what concretely is missing.

The paper itself (`../Interfacce/glauber_paper.tex`) is a sibling directory and
is **not** tracked by this repo. Where the paper already carries an inline
`%NOTE`/`%TODO`, it is cross-referenced below.

Run the fast gates (~10 s on an unloaded machine):

```
python python_scripts/theory/run_checks.py
```

Fourteen gates, exit code 0 iff all pass: finite-N Pfaffian identity (operator
route, free-fermion engine, 300-configuration randomised fuzz), negative
controls, centered-vs-connected at k=4, continuum k=2 anchor, continuum k=3
(validated kernel vs discrete limit, wrong-sign control on the superseded
reference, guard trip), and the box-averaged density `box_fk`.
`pfaffian_verify.py` and `ff_engine.py` are also directly runnable and exit
nonzero on regression. The annihilating-walk dual cross-check is
`dual_multitime_check.py` (~2 min; `--k4` adds the k=4 rung, ~15 s more). The
heavyweight Monte-Carlo collapse and the persistence run are separate long
campaigns, noted as such.

*Last full audit: 2026-08-12 — every fast check below re-run, then the whole
cleanup adversarially reviewed (5 independent lenses, 36 verified findings —
all addressed, see R5); MC campaigns not re-run.*

> **Reproducibility.** The whole import closure of the gate scripts is tracked
> (this was not true before 2026-09-04, when only `run_checks.py` and
> `scaling_check_3time.py` were). Verified by cloning the branch into an empty
> directory and running `run_checks.py` there: **14/14 gates pass, exit 0**.
> Dependencies are in `requirements.txt` (`mpmath` is needed only by
> `droplet_lifetime_check.py`).
>
> Deliberately **not** tracked, and so reproducible only on the machine that
> produced them: the Monte-Carlo output directories (`campaign_eq*/`,
> `results_*/`, `grid_T05/`, `probe_T08/`, `persistence_test/`, …). Only their
> small `.log` files are kept, recording the parameters each grid was run with.
> Re-deriving a published χ² therefore means re-running the campaign, not just
> re-reading a file — see ⚠️ 3.

---

## ✅ Verified

| Paper claim | Statement | Script(s) | Outcome |
|---|---|---|---|
| `thm:majorana`, `thm:bond`, reversibility | Clifford algebra, `μ = i a a`, `L~` symmetric, stationarity | `pfaffian_verify.py` (gate [1]) | all ≤ **3e-16** |
| `e:multipf` / `cor:multitime` (finite N, **pre-scaling**) | `⟨μ…μ⟩ = i^k Pf[Φ]`, k=2,3,4 | `pfaffian_verify.py` (exact 2^M operator route), `ff_engine.py` (closed-form C, G_sp vs brute force at 4 (M,η)) | ≤ **2e-15** everywhere (worst: C-builder at M=10, η=0.9); plus **randomised fuzz** (2026-08-12, gate [1]): **300 configs** incl. 192 with coinciding times, 122 with a repeated bond, 55 with an identical (bond,time) pair, over 4 (M,η), k=2,3,4 — worst **1.7e-15** on both routes |
| `c:centered` + `e:pfhier` reading | zero-diag Pf = **centered** moment (every k); cumulant = single-2k-cycle sum; they **differ** from k=4 (offset = the 2+2 / 2+3 κ-products) | `k4_zerodiag_check.py` (k=3,4,5 at η=0.6, 0.9), gate [3] | identities ≤ **5e-15**; k=4 offset = pair products to 2e-17 (97% of the centered value at η=0.6, 99% at η=0.9). The paper states `e:pfhier` as centered — correct; **do not revert** to "cumulant at every k" (the pre-2026-07-05 draft error) |
| Cor. `thm:multitime` via the **dual** | multi-time bond Pfaffian vs the annihilating-random-walk dual (nested `e:stat`), fully independent of the Jordan–Wigner/Majorana pipeline | `dual_multitime_check.py` (re-run 2026-08-12) | ladder: k=1 statics, σσ two-time vs Bessel `D`, k=2 vs Glauber determinant, k=2,3 vs `ff_engine` `iᵏPf[Φ]` at **both** γ=0.6 and γ=0.9 — all ≤ **4e-16**; k=4 rung (`--k4`, γ=0.6 only by design, window L=24) **3.8e-17**. Window truncation of ℤ, exact subset-space `expm` — no MC, no fermions on the dual side |
| negative controls | Pf[D] (scalar spin kernel) must FAIL in sign **and** magnitude — the anomalous matrix contraction is genuinely needed; multi-time spin Wick must FAIL while equal-time (static Ising Wick) holds — why the theory is stated for **bonds** | `pfaffian_verify.py` (`D_universality_gaps`, `spin_wick_gaps`), gate [2] | Pf[D] off by 2.6e-3 / 1.5e-2 / 9.9e-3 at k=2/3/4 (magnitudes too); spin Wick: equal-time gap **0**, multi-time gaps 1.8e-5, 1.7e-4 |
| `e:Hlimit` | two-point spin → killed heat kernel `H(x,τ)` | `finiteN_validation.py`, `theory_glauber.py`; box collapse `plot_box_vary_h.py` (`campaign_eq/plot_spin.txt`, N=200, 50 000/pt) | χ²/dof **0.67–1.42** over the (h,τ) grid |
| `e:Htorus` | length-`L` torus spin correlator (mode sum) | `chi2_torus_compare.py`, `theory_continuum.H_torus` | collapses exact ring to O(1/N); χ²≈1 |
| `e:fmain` / `e:f` | two-point **interface density** `f = ¼[(∂ₓH)² − H∂ₓ²H]` (quartered Wronskian) | `cont3_final.py` (k=2 anchor), `validate_iface.py`, `theory_continuum.f`; gate [4]; `campaign_eq/plot_iface.txt` (N=200, 50 000/pt) | **sharp prefactor test**: data sit on `f`, a clean **factor 4 below** the un-quartered draft form; χ²/dof **0.49–1.53** over the (h,τ) grid. Triangulated by box MC + exact finite-N `H_interface` + free-fermion engine. Gate [4] is the *algebraic* anchor only (same H,P primitives); the closed forms are tested by the MC rows |
| `thm:pfaffian`(i)+(ii), k=3 **continuum** | connected 3-point density = `(−i/2)³ Pf` of the complex kernel `e:Kpf` — **the kernel now printed in the theorem** | `cont3_final.py`, `paper_kernel_check.py`; gate [5] | matches the discrete N³ Richardson limit within 8% at configs A, B, C (N=20,40; ratios 0.977 / 0.978 / 1.002), → 1.00 at N=80 (`ff_scale.py`). **Negative control** (gate [5b]): the superseded real-kernel reference has the **wrong sign** at A and B (ratios −1.67, −0.94) and is 10% off at C — config-dependent, which is why it went unnoticed. *(Not a gate: `fk(real kernel).real = 0` at odd k is structural — imaginary prefactor × real Pfaffian, Remark `r:Kstruct` — any real kernel satisfies it.)* |
| `box_fk` (box-averaged `e:pfhier` density) | replacement for the superseded `cont_box_cumulant` | `cont3_final.box_fk`; gate [6] | equals the old function at k=2 (where both are valid) to **0.0**; converges to pointwise `fk` as h→0 at k=3 (rel. −7.8e-2 at h=0.1, −3.2e-3 at h=0.02). Line geometry only — the torus `L=` option is deliberately not carried over |
| App. `s:numerics` | finite-N exact forms `C_eq`, `G_two_time`, `H_interface` | `theory_glauber.py`, `finiteN_validation.py` | χ²/dof ≈ **1** vs engine MC |
| equilibrium IC | even-conditioned Bernoulli(p) sampler (Floyd) | `mc_box_correlation.py`, `spin_init.c` | reproduces `C_eq_ring`; box residual ∼1/√nsims |
| k=0 colour bug | cross-time correlator sign loss | fixed, commit `8a8cc21` | all post-fix cross-time data trustworthy |

---

## 🔧 Resolved since the last revision of this file

### R1. The k ≥ 3 kernel "blocker" is closed (paper fixed 2026-06-27)
The previous revision carried a `[BLOCKER]`: Theorem `thm:pfaffian` printed the
real block `[[H,−∂ₓH],[∂ₓH,−∂ₓ²H]]`, identically zero for connected k=3. The
paper has since been corrected: the theorem now prints the **complex** kernel
`e:Kpf` (`P` on the diagonal), with a proof of part (i), and `e:Kblock`
survives only as the k=2 covariance ("observable shadow") in Cor. `c:k2` —
which is legitimate, since at k=2 only `det K` is observable. `cont3_final.py`
keeps the real block under its historical name `Kpaper` (alias `Kblock`) as a
negative control. `run_checks.py`'s former "theorem is false as printed"
message was stale and is removed (2026-08-12), as are the matching stale
headers of `scaling_check_3time.py` and `k4_zerodiag_check.py`.

### R2. Contaminated continuum references repointed (2026-08-12)
`ff_scale.py` and `freefermion_largeN.py` took their k=3 continuum reference
from the *superseded* `scaling_check_3time.py` (real kernel with the real
prefactor — wrong for k≥3). Measured against the validated
`cont3_final.fk(·, Kcont)` at the standard configs: wrong **sign** at A and B
(ratios **−1.63**, **−0.92**), 10% low at C. Both scripts are repointed (`fk` /
new `cont3_final.box_fk`); `ff_scale` now converges cleanly (k=3 ratio
**1.003** at N=80, Richardson 0.93–1.02). The discrete sides were never
affected — only the comparison target. The superseded functions now **raise**
for k≥3 unless called with `allow_wrong_k3=True` (and `ValueError` for k∉{2,3},
where their hard-coded loops would silently drop points); gate [5c] checks the
guard, gate [5b] the sign error. No other caller exists in the repo
(`ff_scale_big.py` already used `cont3_final`).

### R3. `freefermion_largeN.py` marked SUPERSEDED (2026-08-12)
Independent of R2 it is broken: (a) `scaling()` returns **nan** at every N —
its propagator exponentiates the growing modes, overflowing at t = N²τ
(`ff_engine.FF._RT` is the fix: drop the λ>0 modes, which carry zero weight);
(b) its own builder validation fails away from η=0.5 (‖G1−extract‖ = **5.4e-1**
at η=0.7) and it only prints this, never asserts. The working path to the same
physics is `ff_engine.py` + `ff_scale.py`. Header updated accordingly.

### R4. The fast checks are now enforced, not just printable (2026-08-12)
Previously `run_checks.py` gated only the 3 continuum checks; the finite-N
validations printed numbers with no PASS/FAIL and no exit code. Now:
`run_checks.py` = 14 gates, exit 0 iff all pass; `pfaffian_verify.py` and
`ff_engine.py` exit nonzero on regression; the formerly dead `spin_wick_test`
is integrated as a documented control.

### R5. Adversarial review of the cleanup (2026-08-12)
Five independent review lenses (refactor integrity, mathematical conventions,
gate strength, documentation accuracy, guard airtightness), every finding then
re-verified by a separate agent: 36 findings confirmed, 1 refuted (a timing
measured under machine load). All confirmed items are fixed, notably:
- the first version of gate [5b] ("`e:Kblock` gives zero at k=3") was a
  **tautology** — any real kernel passes it — and is replaced by the wrong-sign
  control on the superseded reference plus the guard-trip gate;
- the Pf[D] control now gates the magnitude gap as well as the signed gap;
- the fuzz no longer skips identical (bond,time) pairs — the identity holds
  there too (d ≤ 7e-17), so coverage is 300 configs, not 245;
- `box_fk` is gated; `fk`/`box_fk`/`cont3_final` docstrings say "centered"
  where they said "connected" (they coincide for k≤3 only);
- the persistence row (below) was mis-attributed and rests on 4 trajectories;
- the "~750k samples/point" figure in the old open item had no source on disk
  (the grids are 50 000/pt, see ⚠️ 3).

---

## ⚠️ Open — still to verify (with motivation)

### 1. **k ≥ 4 in the continuum** is unconfirmed
The finite-N side is settled to k=5 (`e:multipf` moments, `c:centered`, cycle
sums — `k4_zerodiag_check.py`), and the continuum hierarchy is checked at k=2
and k=3. What is missing is the **limit** at k≥4: that N⁴×(centered / cycle
sum) converges to the corresponding `e:Kpf`-block Pfaffian objects.
- **Motivation:** the paper asserts an all-orders Pfaffian point field; only
  the first two non-trivial orders of the *limit* are tested.
- **Plan:** `NUMERICS_PROGRAM.md` carries a sharpened, costed plan (centered
  route avoids the catastrophic cancellation; geometry with a dominant 2+2
  product term).

### 2. **Regime-family probes (`ss:numregimes`) are PLANNED, not run**
The paper text itself says *"The planned probes are…"*. Only single
representative kymographs (`f:regimes`) exist. Missing:
- (a) the interface density `f` as an **order parameter for nucleation
  survival** across the deformation `aₙ ∼ Nᵃ`, `cₙ ∼ N^{a−2}`, `a ∈ ℝ`;
- (b) **tagged-colour persistence** across the regime family. Even the `a=0` /
  zero-T anchor (`3/8`, `DHP95`) is only **indicative**: the numbers on disk
  come from `persistence_test/persistence.txt` (`Persistence` binary, N=100,
  a=0, **nsims=4**, `configuration.json`) — late slope **−0.33** (τ≥0.1) vs
  −0.375, amplitude √(2N²τ)ρ ≈ **0.26–0.29** vs 1/√(4π)=0.282. Four
  trajectories give no usable error bar; the previous revision of this file
  credited these numbers to `interface_survival.py`, which computes neither.
- **Motivation:** the regime section currently shows qualitative pictures only;
  no quantitative collapse backs the conjectured transition at `a=−1`.

### 3. **Precision of the interface collapse**
The published χ² grids (`campaign_eq/`, N=200) rest on **50 000 samples/point**.
The one high-statistics run on disk, `campaign_eq2/` (**800 000/pt**, h=0.25
only), gives interface χ²/dof **0.85 / 1.97 / 1.63** at τ=0.05/0.1/0.2 (spin
1.04 / 0.84 / 0.85 / 0.47 at τ=0.05…0.4): at 16× the statistics a finite-N
residual becomes resolvable at the sharper time corners of the interface
observable.
- **Motivation:** decide whether the residual is the expected O(1/N) correction
  (then quote it) or a systematic; extend `campaign_eq2` to the other box
  widths (`combine_passes.py`).

### 4. Out of scope here
Process-level convergence (uniform control of higher cumulants) is deferred to
the companion paper `[BGL]`; it is not a numerical task for this manuscript.
Whether the paper's *proof* of `thm:pfaffian` is fully rigorous is likewise a
referee's question, not a numerical one — the numerics certify the statement.
