---
layout: default
title: The Programs
---

[Back to the index](index.html)

# The Programs

There are three programs, built on [the Gillespie engine](engine.html): `Gillespie` for the trajectories and the statistics of an ensemble, `Corr2pt` for the correlations in time on the lattice, `BoxCorr` for the correlations of the box averages. The Python scripts draw the outputs and compare them with the exact formulas.

---

## 🔹 Compiling

The programs need a C compiler and the **GSL**. The scripts need Python with the packages of `requirements.txt`.

```bash
bash compile_gillespie.sh      # ./Gillespie
bash compile_corr2pt.sh        # ./Corr2pt
bash compile_box_corr.sh       # ./BoxCorr
```

The compile scripts look for the GSL in `/usr/include` and `/usr/lib64`. Another place is given with `GSL_INC` and `GSL_LIB`:

```bash
GSL_INC=/usr/local/include GSL_LIB=/usr/local/lib bash compile_gillespie.sh
```

They compile with `-fopenmp`, and the simulations of an ensemble run on several threads (`OMP_NUM_THREADS`). With a compiler without OpenMP, remove `-fopenmp` from the script: the programs run on one thread.

The sticky Brownian motion is a submodule, see [its page](sticky.html):

```bash
git clone --recurse-submodules https://github.com/alessandrogubbiotti/Glauber.git
```

---

## 🔹 The file config.txt

A line is `key=value`, and a line that starts with `#` is a comment. The keys that a program does not know are ignored.

- **`N`**: the scaling parameter, the number of sites in a unit of macroscopic length.
- **`L`**: the macroscopic length of the ring, an integer. The ring has `N*L` sites.
- **`T`**: the final macroscopic time.
- **`nu`**: the intensity $\nu$ of the invariant measure. The default is `1.0`.
- **`rate_mode`**: `glauber`, `scaled` or `explicit`, see [the model](model.html). The default is `glauber`.
- **`slowdown_exponent`**: the exponent $a$ of the `scaled` mode. The default is `0.0`.
- **`annihilation`**, **`creation`**: the two rates, read in the `explicit` mode only.
- **`resolution`**: the number of saved times in a unit of macroscopic time.
- **`N_simulations`**: the number of independent simulations.
- **`init_function`**: the initial configuration, from `1` to `6`.
  1. The first spin is $+1$ or $-1$ with probability $1/2$, and each spin is different from the previous one with probability `init_param`.
  2. $+1$ on the first half of the ring and $-1$ on the second: two interfaces.
  3. All the spins $-1$: no interfaces.
  4. Independent spins, $+1$ with probability `init_param`.
  5. Alternating spins: an interface on every bond.
  6. The invariant measure: an interface on each bond independently with probability `init_param`, conditioned to an even number of interfaces. Put `init_param` equal to $\nu/(N+\nu)$.
- **`init_param`**: the parameter of the initial configuration.

`Gillespie` reads all the keys. `Corr2pt` reads `N`, `L`, `nu` and the rates. `BoxCorr` does not read the file: its parameters are on the command line.

---

## 🔹 Gillespie

```bash
./Gillespie
```

It reads `config.txt` in the current directory and runs `N_simulations` independent simulations up to the time `T`, saving the configuration `resolution` times in a unit of time. The seed is taken from the system. With `GILLESPIE_SEED=<n>` the run is repeated exactly.

The output is in the directory `results_N<N>_T<T>_L<L>_Ann<ann>_create<create>_gillespie`:

- **`configuration.json`**: the parameters and the derived quantities, $\nu$, $\beta$, $\gamma$, $\eta$, $p$, the two rates and the factor $2N^2$ of the time dictionary.
- **`trajectory/spins_output.bin`**: the trajectory of the first simulation. At every saved time the `N*L` spins as 32-bit integers, then the marker `-9999`.
- **`mean_spin_t.txt`**: at every saved time, the mean and the standard deviation over the simulations of the spin at the site `N*L/2`.
- **`corr_pair_t.txt`**: at every saved time, the mean and the standard deviation of the product of the spins at the sites `N*L/4` and `3*N*L/4`.
- **`time_autocorr.txt`**: the correlation in time of the spin at the site `N*L/2`, averaged over the time and the simulations.
- **`spin_corr_sel_t.txt`**: the spin correlation $C(r) = \frac{1}{\ell}\sum_x \sigma_x\sigma_{x+r}$ at 20 saved times, averaged over the simulations.
- **`spin_corr.txt`**, **`iface_corr.txt`**: the spin correlation and the number of pairs of interfaces at distance $r$ per site, averaged over all the saved times and the simulations.
- **`timing_stats.txt`**: the wall time and the number of events of a simulation.
- **`Statistics.txt`**: only its header is written.

At the end it calls `python_scripts/plot_gillespie.py`, that draws the pictures of these files in the same directory. The picture of the trajectory is drawn by

```bash
python3 python_scripts/Plot_Glauber_tidy_Binary.py <results>/trajectory
```

---

## 🔹 Corr2pt

```bash
./Corr2pt [config.txt] [n_sims] [equilibrium|disordered]
```

The defaults are `config.txt`, `2000` and `equilibrium`. The times are $T = 0.25, 0.5, 1, 2$ and are fixed in the source. The correlations are averaged over the translations and over the `n_sims` simulations. With `CORR2PT_SEED=<n>` the run is repeated exactly.

- **`equilibrium`**: the start is the invariant measure. The observables are the time-delayed correlations $G_n(T) = \langle \sigma_0(0)\,\sigma_n(T)\rangle$ and $H_n(T) = \langle \eta_0(0)\,\eta_n(T)\rangle$, where $\eta_b = (1 - \sigma_b\sigma_{b+1})/2$ is the interface at the bond $b$.
- **`disordered`**: the start is independent uniform spins. The observable is the correlation at a single time after the quench, $C_n(T) = \langle \sigma_0(T)\,\sigma_n(T)\rangle$.

The output is in the directory `results_N<N>_nu<nu>_L<L>_<start>_Ann<ann>_create<create>_corr2pt`:

- **`corr2pt_multi.txt`**: $G_n(T)$, with its standard error, for each $T$.
- **`iface_corr2pt_multi.txt`**: $H_n(T)$.
- **`corr_eqtime_multi.txt`**: the correlation at a single time, at the times $0$ and $T$.
- **`corr2pt_pairs.txt`**: the correlations between two positive times.
- **`corr2pt_cov.bin`**: the second moments of $G$, for a $\chi^2$ with the covariance.
- **`configuration.json`**.

At the end it calls `python_scripts/plot_theory_vs_mc.py`, that overlays the exact formulas and prints a $\chi^2/\text{dof}$ per curve. The formulas are exact for the `glauber` rates.

---

## 🔹 BoxCorr

```bash
./BoxCorr [--N 12] [--L 4] [--nu 1.0] [--h 0.5] [--M 15]
          [--T 0.05,0.1,0.2,0.4] [--nsims 1000] [--seed 0]
          [--observable spin|iface] [--checkpoint-interval 5] [--no-resume]
          [--out box_corr_c]
```

The box magnetization $m_{\text{box}}(x,t)$ is the average of the spins over the box of macroscopic width $h$ centred at $x$. The program estimates

$$
C(x, T) = \langle\, m_{\text{box}}(L/2, 0)\; m_{\text{box}}(x, T)\,\rangle
$$

at `M` points $x$ and at the times of `--T`, starting from the invariant measure, with the `glauber` rates. Each point $(x, T)$ has its own ensemble of `nsims` simulations, so the estimates are independent. With `--observable iface` the spins are replaced by the interfaces minus their density $p$, and the correlation is multiplied by $N^2$.

`--h` accepts a list, `--h 0.1,0.25,0.5`: all the widths are measured on the same trajectories.

The run is saved in `<out>/checkpoint.bin` every `--checkpoint-interval` seconds, and a run that has been interrupted restarts from there. `--no-resume` starts again. The seed is `--seed`; without it one is drawn and kept in the checkpoint.

The output, for each width, is `<out>/h<h>/box_corr.txt`, with $x$ and, for each $T$, the estimate and its standard error, and `<out>/h<h>/configuration.json`.

---

## 🔹 The Scripts

They are in `python_scripts`, and each one starts with a docstring that says what it does and how it is called.

- **The trajectory:** `Plot_Glauber_tidy_Binary.py`.
- **The output of `Gillespie`:** `plot_gillespie.py`.
- **The output of `Corr2pt`:** `plot_theory_vs_mc.py`, `plot_time_delayed_corr.py`.
- **The three regimes:** `plot_dynamics.py`, `plot_corr_comparison.py`, see [the three regimes](regimes.html).
- **The output of `BoxCorr`:** `mc_box_correlation.py --plot-only <dir>` overlays the exact formula; `validate_iface.py` for the interfaces; `plot_box_vary_h.py` draws one picture for each width; `combine_passes.py` and `combine_shards.py` put together runs with different seeds; `analyze_nsim_h.py` and `chi2_torus_compare.py` compute the $\chi^2$ against the formulas.
- **The formulas:** `theory_glauber.py` on the lattice, `theory_continuum.py` in the limit.
- **The checks of the paper:** the directory `theory`.

`plot_corr2pt.py` reads a file of an earlier version of `Corr2pt`.

---

## 🔹 The Checks

```bash
python python_scripts/theory/run_checks.py
```

These are the fourteen checks of the formulas of the paper, and the exit code is $0$ if all of them pass. The file [`VERIFICATION.md`](https://github.com/alessandrogubbiotti/Glauber/blob/main/VERIFICATION.md) says, for each claim of the paper, which script verifies it and with which result, and what is still open.

```bash
ruff check
```

This is the automatic checker for docstrings, with the rules of `ruff.toml`. GitHub runs it at every push on `main`.
