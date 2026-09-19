# Cycloidal Drive Optimization

This project performs a reproducible mixed-integer, multiobjective screening
optimization for a one-tooth-difference cycloidal reducer intended for robotic
actuator studies.

## Requirements

- MATLAB R2025b or a compatible release
- Global Optimization Toolbox (gamultiobj)
- Java enabled is recommended for SHA-256 source hashes; runs still save if
  hashing is unavailable

Parallel execution is disabled by default because parallel genetic-algorithm
population generation is not reproducible. Enable it only when repeatability
is not required.

## Run

    cd('C:\Users\User\Desktop\CycloidalDriveOpt')
    run_tests
    results = main_optimize();

Results are saved under results/ using a timestamped filename. Each result
contains the parameters, solver state, source-file hashes, RNG state, robust
validation statistics, all returned candidates, the globally nondominated
robust-feasible set, and the ranked best candidate.

## Live generation feedback

Live feedback is enabled in defaultParams. At generation zero and each
completed generation, the command window and a figure show the selected
design, its profile and holes, dimensions, mass, backlash, contact pressure,
elastic stiffness, safety factors, and constraint-violation history.

The selected design minimizes the existing weighted, normalized score on
the current generation's nominal feasible Pareto front. This is a preference
choice among competing objectives, not an absolute best. Normalization is
recomputed for each front, so scores are not comparable between generations.
If none is feasible, the display explicitly shows the least-violating design.
Robustness is still checked after optimization.

Set p.livePlot=false for console-only feedback, or p.liveFeedback=false to
disable the callback entirely. Closing the plot disables its updates without
stopping optimization. Generation snapshots are stored in results.liveHistory.
The callback re-evaluates unique population members to check all constraints;
this adds nominal evaluation time but never runs Monte Carlo.

## Design variables

The optimizer uses:

    [Rp, e, dRing, dOut, tDisc, ringClearance, rOutPitch, outputClearance]

All lengths are metres. Legacy seven-value calls remain accepted and use
p.outputClear_nom as the eighth value.

The two discrete variables are indices into p.Zr_set and p.nOut_set.
Zc is always enforced as Zr-1.

## Objectives

1. disc mass in kg (minimize);
2. no-load lost-motion surrogate in degrees (minimize);
3. Hertz line-contact pressure in MPa (minimize);
4. negative elastic torsional stiffness in N*m/rad (minimize, which maximizes
   positive stiffness).

All feasibility conditions are dimensionless nonlinear constraints. Penalties
are not mixed into objective units.

## Validation strategy

- The nominal evaluator is deterministic and sweeps multiple input phases.
- The optimization uses one mixed-integer Pareto run, not disconnected fronts.
- Robustness is evaluated only after optimization using fixed common random
  samples and a Wilson 95% confidence interval.
- A candidate passes robust screening only if the upper confidence bound on
  failure probability is no greater than p.maxFailureProb.
- Final candidates use NthetaValidation; inner optimization and uncertainty
  runs use lower configurable resolutions.

## Engineering limitations

This is a screening model, not certification. Root bending and clearance lost
motion remain analytical surrogates and are explicitly labeled in outputs.
Before releasing hardware, correlate them with:

- independent CAD motion/contact analysis;
- nonlinear contact FEA over input phase;
- printed-material coupon data at operating temperature and humidity;
- pin, shaft, bearing, carrier, and housing checks;
- fatigue, creep, wear, lubrication, efficiency, and thermal tests;
- prototype transmission-error, reversal, stiffness, shock, and life tests.

The optimized mass covers the cycloidal disc only. Actuator-level mass,
reflected inertia, motor integration, bearings, housing, carrier, and a
possible dual-disc architecture must be assessed separately.
