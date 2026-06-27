# Verification status — scaling-limit claims of `glauber_paper.tex`

**Motivation.** The scaling-limit section of the paper makes a sequence of
quantitative claims (Theorem `thm:pfaffian`, the boxed density `e:fmain`, the
matrix kernel `e:Kblock`, the regime family of `ss:numregimes`). Several are
backed by numerics that live in *this* repository, but the correspondence
"paper claim ⇄ script ⇄ measured number" is currently scattered across script
docstrings, campaign folders, and `%NOTE` comments inside the `.tex`. This
document makes that correspondence **explicit**: for every claim, what verifies
it, what the measured outcome was, and — for the open items — *why* the check
still matters and what concretely is missing.

The paper itself (`../Interfacce/glauber_paper.tex`) is a sibling directory and
is **not** tracked by this repo; only the verification code is. Where the paper
already carries an inline `%NOTE`/`%TODO`, it is cross-referenced below.

Run the fast analytic gates:

```
python python_scripts/theory/run_checks.py
```

(seconds; free-fermion + closed-form only). The heavyweight Monte-Carlo
collapse and the persistence run are separate long campaigns, noted as such.

---

## ✅ Verified

| Paper claim | Statement | Script(s) | Outcome |
|---|---|---|---|
| `e:Hlimit` | two-point spin → killed heat kernel `H(x,τ)` | `python_scripts/theory/finiteN_validation.py`, `theory_glauber.py`; box collapse `plot_box_vary_h.py` | χ²/dof **0.7–1.4** over the (h,τ) grid |
| `e:Htorus` | length-`L` torus spin correlator (mode sum) | `chi2_torus_compare.py`, `theory_continuum.H_torus` | collapses exact ring to O(1/N); χ²≈1 |
| `e:fmain` / `e:f` | two-point **interface density** `f = ¼[(∂ₓH)² − H∂ₓ²H]` (quartered Wronskian) | `cont3_final.py` (k=2 anchor), `validate_iface.py`, `theory_continuum.f` | **sharp prefactor test**: data sit on `f`, a clean **factor 4 below** the un-quartered draft form; χ²/dof **0.5–1.5**. Triangulated by box MC + exact finite-N `H_interface` + free-fermion engine |
| `thm:pfaffian`(ii), k=3 **discrete** | `⟨μμμ⟩ = i³ Pf` (anomalous, P on diagonal) | `ff_engine.py`, `pfaffian_verify.py` | agreement to **1e-17**; `Pf[D]` (number-conserving) fails ⇒ genuinely anomalous |
| `thm:pfaffian`(ii), k=3 **continuum** | connected 3-point density = continuum Pfaffian | `cont3_final.py` | corrected kernel matches discrete N³ limit, Richardson ratio **0.90–1.00** (see ⚠️ #1 — uses the *corrected* kernel, not the one printed in the theorem) |
| `e:pfhier`, k=4 **discrete** | Pfaffian structure persists | `ff_engine.py` | to **1e-17** (discrete only) |
| App. `s:numerics` | finite-N exact forms `C_eq`, `G_two_time`, `H_interface` | `theory_glauber.py`, `finiteN_validation.py` | χ²/dof ≈ **1** vs engine MC |
| equilibrium IC | even-conditioned Bernoulli(p) sampler (Floyd) | `mc_box_correlation.py`, `spin_init.c` | reproduces `C_eq_ring`; box residual ∼1/√nsims |
| `3/8` persistence (`DHP95`, a=0 / zero-T) | annihilating-walk persistence exponent | `Persistence` binary, `interface_survival.py` | late slope **−0.35** vs −0.375; amplitude √(2N²τ)ρ ≈ **0.27** vs 1/√(4π)=0.282 |
| k=0 colour bug | cross-time correlator sign loss | fixed, commit `8a8cc21` | all post-fix cross-time data trustworthy |

---

## ⚠️ Open — still to verify (with motivation)

### 1. **[BLOCKER] Theorem `thm:pfaffian` prints the WRONG kernel for k ≥ 3**
The theorem states the real block
`K = [[H, −∂ₓH], [∂ₓH, −∂ₓ²H]]` (`e:Kblock`). For the connected **k = 3**
correlator this is **identically zero**: `fₖ = (−i/2)ᵏ Pf[…]`, and at k=3 the
prefactor `(−i/2)³` is imaginary while the Pfaffian of a real antisymmetric
matrix is real, so the physical (real) part vanishes. Gate `[2]` of
`run_checks.py` reproduces this: `fk(Kpaper).real == 0` exactly, at every test
configuration.

The validated replacement is the **complex** kernel (P on the diagonal):

```
K(z,dt) = [[ P,              i(2H − P − H') ],
           [ i(P − 2H − H'),  P             ]]      P = −¼∂_τH,  H' = ∂_zH
```

which reduces to the same k=2 `f` and reproduces the discrete N³ limit.

- **Motivation:** as printed, the theorem is *literally false* for k ≥ 3.
- **Status:** validated numerically (Richardson), **not yet a rigorous proof**.
- **Action:** promote the corrected kernel from the `.tex` `%NOTE` (just after
  `e:Kblock`) into the *statement* of `thm:pfaffian`; ideally prove it.
- **Also:** `python_scripts/theory/scaling_check_3time.py` still hard-codes the
  wrong kernel and is **superseded by `cont3_final.py`** (flagged in its header).

### 2. **k ≥ 4 in the continuum** is unconfirmed
Discrete Pfaffian structure is checked to k=4; the *continuum* hierarchy
`e:pfhier` is checked only at k=2 and k=3.
- **Motivation:** the paper asserts an all-orders Pfaffian point field; only the
  first two non-trivial orders of the continuum limit are tested.

### 3. **Regime-family probes (`ss:numregimes`) are PLANNED, not run**
The paper text itself says *"The planned probes are…"*. Only single
representative kymographs (`f:regimes`) exist. Missing:
- (a) the interface density `f` as an **order parameter for nucleation
  survival** across the deformation `aₙ ∼ Nᵃ`, `cₙ ∼ N^{a−2}`, `a ∈ ℝ`;
- (b) **tagged-colour persistence** across the regime family (only `a=0` /
  zero-T is done — item ✅ above).
- **Motivation:** the regime section currently shows qualitative pictures only;
  no quantitative collapse backs the conjectured transition at `a=−1`.

### 4. **Precision of the interface collapse**
The published interface χ²≈1 rests on ~750k samples/point; the high-statistics
campaign (target 1e7/pt, `campaign_eq/` + `combine_passes.py`) was still
accumulating.
- **Motivation:** tighten error bars at the widest box (h=2) / sharpest time
  corner, where finite-N residual is largest.

### 5. Out of scope here
Process-level convergence (uniform control of higher cumulants) is deferred to
the companion paper `[BGL]`; it is not a numerical task for this manuscript.
