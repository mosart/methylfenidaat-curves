# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A single-file [marimo](https://marimo.io) reactive notebook that models methylphenidate (MPH) plasma concentration curves. Users configure body parameters and a dosing schedule; the notebook computes and visualises personalised PK curves, peak times, and a hypothetical rebound moment.

## Running the notebook

Dependencies are declared inline via PEP 723 (the `# /// script` block at the top of the file), so no separate `requirements.txt` or `pyproject.toml` is needed.

```bash
# Preferred: run with uv (resolves deps automatically)
uvx marimo run methylfenidaat-curves.py

# Or edit interactively
uvx marimo edit methylfenidaat-curves.py

# Or if marimo is already installed
marimo run methylfenidaat-curves.py
marimo edit methylfenidaat-curves.py

# Or running in a sandbox and watch changes done by agents
marimo edit methylfenidaat-curves.py --sandbox --watch
```

## Checking the notebook
```bash
# Running marimo check on this notebook reveals several issues
marimo check methylfenidaat-curves.py

# marimo check includes a --fix flag that automatically resolves simple issues
marimo check methylfenidaat-curves.py --fix 

# For more complex scenarios where the fix might change code behavior, we provide --unsafe-fixes
marimo check methylfenidaat-curves.py --unsafe-fixes

```

## Notebook architecture

Marimo notebooks are Python files where each `@app.cell` function is a **reactive cell**. Cells re-execute automatically when their inputs change. Return values become variables available to downstream cells; cells that return nothing are pure UI side-effects.

The execution flow in this notebook:

1. **`imports`** — marimo, numpy, altair, pandas
2. **`lichaamskenmerken_defaults`** — UI sliders/radio for weight, height, sex, metabolism
3. **`pk_parameters`** — computes derived PK scalars (ke_scale, scale_factor, half-lives) from body params
4. **`pk_functies`** — defines `ir_curve(t, dose, t0)` and `er_curve(t, dose, t0)` as closures over PK params
5. **`y_schaal`** — pre-computes Y_MAX across all schemas to keep the axis stable when switching tabs
6. **`schema_sliders`** + **`ontwaaktijd_invoer`** — UI for wake time and per-schema timing offsets
7. **`grafiek`** — main chart cell: builds dose list from active tab, runs curves, detects rebound, renders Altair chart
8. **`grafiek_beschrijving`** — text summary accordion
9. **`schema_tabs_cel`** — tab UI wiring dosing schemas to sliders and alarm callouts
10. **`lichaamskenmerken_ui`** — layout cell for body-param widgets + PK parameter accordion
11. **`rebound_info`** + **`achtergrond`** + **`disclaimer`** — informational UI cells

## PK model

- **IR curve**: AUC-normalised Bateman function, ka=3.5 h⁻¹, ke₀=0.277 h⁻¹ → Tmax ~1.3 h, t½ ~2.5 h
- **ER curve (Concerta 27 mg)**: 22% IR coating (ka=3.5) + 78% osmotic pump (ka=0.38, ke₀=0.198)
- **Body scaling**: ke scales with `(W/70)^0.75` (allometric law); Vd scales with LBM via Janmahasatian formula; CES1 metabolism variability via `exp(slider)`
- **Rebound detection**: first time-point after the peak where total curve drops below 25% of peak value

## Key invariants

- `Y_MAX` is computed once over all four schemas so the Y-axis never rescales when the user switches tabs — do not recompute it inside `grafiek`.
- `ir_curve` and `er_curve` are zero outside their active windows (14 h and 16 h respectively) via `np.where` guards.
- Time is represented as a float (hours since midnight, e.g. 7.5 = 07:30); `fmt(h)` converts to HH:MM strings.
