A short document on frameworks for workflows

## preferences/desiderata

- hash-based rather than date-based (and ability to manually mark as "outdated"; this is via `touch` for date-based frameworks)
- file-level granularity for better version control logging (`targets` deprecates this)
- convenience of a 'side channel' for large/slow outputs
- easy to do intermediate debugging (load in dependencies for a given step, run code manually/interactively)

## targets

### pro

- new hotness
- active, relatively well documented/widely used
- supports dynamic branching
- supports HPC
- R only (other machinery could be called via `system()` etc. but?
- `targetopia` ???

## shellpipes

- nice automatic rules
- poorly documented

## plain old make

## targets questions

- correct way to modify earlier part of pipeline but not re-run later stuff? (`tar_invalidate()` ? )
- literate programming/docstrings for targets etc.
- is there a `tar_make(ncores= ...)`
