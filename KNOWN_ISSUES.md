# Known issues

Defects that are known and not yet fixed, one entry each, with its kind and its evidence.

## KI-1 · docs · the JET testset labels do not name every equilibrium

`test/quality/jet.jl:41` labels each testset `"$(nameof(typeof(equ)))"`. `SolovevEquilibriumFRC`,
`SolovevEquilibriumITER` and `SolovevEquilibriumNSTX` all construct a `SolovevEquilibrium`
(`src/analytic/solovev.jl:225`), and the three X-point constructors all construct a
`SolovevXpointEquilibrium`. So a report on one of them does not say which equilibrium it is. The
fix: carry `(name, equ)` pairs, as `test/analytic/equilibria.jl` does.

## KI-2 · docs · the JET file's comment says too much

`test/quality/jet.jl:10-11` says that the analysis reports "a runtime dispatch or a captured
variable". `report_opt` reports a runtime dispatch; it reports a captured variable only where that
variable causes one. The fix: say "a runtime dispatch".

## KI-3 · upstream · Revise prints EMFILE errors in the test log

JET 0.12 loads Revise. In a `Pkg.test()` run, Revise's file watcher prints 8
`UNHANDLED TASK ERROR: IOError: FolderMonitor: too many open files (EMFILE)` stack traces into the
log. The test totals do not change (4922 pass on the branch that adds `test/quality/jet.jl`).
