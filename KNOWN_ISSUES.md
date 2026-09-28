# Known issues

## Upstream

### K1 · Revise prints EMFILE errors in the test log

- **location:** `test/quality/jet.jl`
- **evidence:** JET 0.12 loads Revise, and its file watcher runs out of file handles.
  `grep -c 'UNHANDLED TASK ERROR'` on a `run-tests.jl full` log of the branch that adds
  `test/quality/jet.jl` counts 8 blocks, each an
  `IOError: FolderMonitor: too many open files (EMFILE)` stack trace. The same count on the log of
  its base gives 0. The test totals do not change.
- **kind:** upstream
- **found:** 2026-09-28

### K2 · The fatou suppression covers the whole `plot_equilibrium` docstring

- **location:** `src/plots.jl:2`
- **evidence:** `# fatou-ignore invalid-docstring-code` goes above the docstring opener, so it
  applies to every code fence in that docstring, not only to the `size = ...` placeholder at
  `src/plots.jl:7`. A new invalid fence in the same docstring gets no report. fatou has no
  suppression for one line inside a docstring.
- **kind:** upstream
- **found:** 2026-09-28
