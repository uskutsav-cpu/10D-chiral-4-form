# P12 projection: resource analysis, and a forecast I got wrong

## 1. What I predicted, and why it was wrong

Partway through the degree-12 projection I extrapolated from three timing
points:

    cols  1-12:   76 s   ( 6.3 s/col)
    cols 13-24:  609 s   (50.8 s/col)
    cols 25-36: 1086 s   (90.5 s/col)

and concluded the deceleration was roughly linear, projecting ~1.7 h for the
first prime, ~4 h overall, and an OOM kill before completion given RSS at
1566 MB against a 1.5 GB ceiling with ~60 MB free.

The full curve:

| columns | seconds | s/col |
|---|---:|---:|
| 1-12 | 76 | 6.3 |
| 13-24 | 609 | 50.8 |
| 25-36 | 1086 | 90.5 |
| 37-48 | 807 | 67.3 |
| 49-60 | 172 | 14.3 |
| 61-72 | 168 | 14.0 |
| **total** | **2925** | **40.6** |

The cost per column **peaked in the middle and then fell by a factor of six**.
Prime 32749 finished in 49 minutes, not 1.7 hours. RSS *decreased* from
1566 MB to 962 MB rather than growing without bound, and the process was never
close to being killed.

**Why the extrapolation failed.** Column cost tracks the structural complexity
of each atlas graph, and the atlas is not ordered by complexity. Three
consecutive points from a rising stretch say nothing about the shape of the
whole curve. Fitting a trend to a monotone prefix of a non-monotone series is
a straightforward error and I made it.

**What I should have done**: sampled cost across the *whole* column range
before forecasting, or simply reported "unknown remaining time, monitoring"
rather than converting three points into an hours-scale prediction and an
OOM claim.

## 2. What was nonetheless real

- Free memory genuinely reached ~54-62 MB with two jobs running, and three
  pytest runs were silently killed at that level earlier in this project.
- Terminating the cheap, rerunnable degree-10 job to protect the expensive one
  was the correct triage regardless of the forecast, and it recovered
  62 MB -> 1030 MB.
- The projection script genuinely had no checkpointing, so *had* it died at
  column 60, roughly 45 minutes would have been unrecoverable. That gap was
  real; the probability I attached to it was not.

## 3. Status of the checkpoint work

`src/sdinv/projection_checkpoint.py` and its 8 tests are worth keeping: the
degree-12 run takes ~50 min/prime, and a future run over five fit primes plus
two holdouts is ~6 hours, where an interruption is a genuine risk.

But it is **not yet wired into the projection scripts**, and the run it was
built to protect completed without it. It should be integrated before the next
multi-prime run, not treated as urgent mid-flight work.

## 4. Rules adopted

1. Do not extrapolate a completion time from fewer than ~5 points spread
   across the full work range.
2. State "unknown, monitoring" rather than producing a number that will be
   read as a forecast.
3. Distinguish *"this job has no checkpointing"* (a fact) from *"this job will
   die"* (a prediction). The first justifies building checkpoints; only the
   second justifies killing a healthy run.
4. Memory pressure on this machine is real; forecasts about it still need
   evidence.
