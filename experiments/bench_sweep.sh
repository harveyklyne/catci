#!/bin/sh
# Size + power, MLP vs XGBoost vs oracle, on identical data.
#
# run.py seeds from SeedSequence(seed) and the config does not enter the seed
# stream, so a given (config, seed) draws the *same* datasets whichever learner
# is named -- the comparisons below are paired, not two independent Monte Carlo
# runs, which is what makes small power gaps readable at these rep counts.
set -e

WT=/Users/harveyklyne/Documents/Github/catci/.claude/worktrees/mlp-learner
PY=/Users/harveyklyne/miniforge3/envs/catci/bin/python
export PYTHONPATH="$WT/src:$WT/experiments"
export OMP_NUM_THREADS=1

SIZE_REPS=${SIZE_REPS:-1000}
POWER_REPS=${POWER_REPS:-200}
WORKERS=${WORKERS:-8}
STRENGTHS=${STRENGTHS:-0.2,0.6,1.0,1.4,1.8}

echo "########## SIZE (reps=$SIZE_REPS) ##########"
for learner in oracle xgb mlp; do
  for setting in lin_lin sin_sin; do
    echo ">>> size $setting learner=$learner"
    $PY "$WT/experiments/run_size.py" "$setting" \
        --reps "$SIZE_REPS" --workers "$WORKERS" --learner "$learner"
  done
done

echo "########## POWER (reps=$POWER_REPS, strengths=$STRENGTHS) ##########"
for learner in xgb mlp; do
  for cfg in lin_lin_step sin_sin_binary_tree; do
    echo ">>> power $cfg learner=$learner"
    $PY "$WT/experiments/run.py" "$cfg" \
        --reps "$POWER_REPS" --workers "$WORKERS" \
        --strengths "$STRENGTHS" --learner "$learner"
  done
done

echo "########## SWEEP DONE ##########"
