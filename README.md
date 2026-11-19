# Laguerre-Based MOR Numerical Examples

This directory contains the MATLAB files supplied for the numerical examples.

## Files

- `run_numerical_examples.m` — main script containing the three numerical
  scenarios: interpolation/model-based Laguerre MOR, data-driven realization,
  and minimum-phase preservation.
- `laguerre_SL.m` — constructs the Laguerre generator matrix and the
  corresponding Laguerre input matrix.

## Running the main script

Place all MATLAB source files required by the main script in this same
directory, then run:

```matlab
run_numerical_examples
```

The main script adds its own directory to the MATLAB search path, so the
project can be kept flat without a separate `src` folder.

## Additional dependencies

The uploaded main script calls the following project-specific functions,
which were not included among the two files supplied for this organization
task and therefore must also be placed in this directory before the script
can run:

- `RLC_Ladder.m`
- `rk_icop.m`
- `Zlag_laguerre.m`
- `relative_H2_error.m`

The script also uses MATLAB functionality for state-space models, impulse
responses, Bode plots, balanced truncation, Lyapunov equations, and
transmission-zero calculations. The relevant MATLAB toolboxes must be
available in the MATLAB installation.

## Notes on preservation of the source code

The numerical statements and computational sequence of the supplied main
program have been retained. Comments and section descriptions have been
added for readability. The path setup has been adapted to the requested
flat, single-directory upload layout.

The two `norm(...,'inf')` calls in Scenario C are retained unchanged.
They are included for comparison purposes only; no improvements or
optimizations have been applied. For large-scale systems, more appropriate
algorithms may be selected to enhance computational efficiency.

The source files have been annotated but have not been executed or
numerically validated in MATLAB as part of this packaging step.
