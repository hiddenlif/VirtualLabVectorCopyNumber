# Virtual Lab Task: Digital PCR Vector Copy Number

Allocate digital PCR wells across three synthetic DNA lots, choose dilution (DNA input), controls and culture time, and use the returned partition counts to estimate integrated vector copies per cell-equivalent and decide whether each lot supports an integration claim. Choices change the simulated counts and use up a well budget, DNA stock and virtual days.

**Task type:** Interactive simulated experimentation. **Version:** 1.0.0.

All lots, measurements and plate conditions are synthetic. The model is a simplified droplet-counting model with differential inhibition and dilution of non-integrated DNA by culture; it is not a validated assay and a vector copy count here is not proof of integration in any real sample.

## Contents

- `TASK.md`: task, deliverable and acceptance rules.
- `app/SPEC.md`: protocol and model. `app/lab.py`: simulator, ledger and scoring. `app/harness.py`, `app/labclient.py`: the JSON-lines protocol.
- `app/solver.py`: starter to replace. `app/practice.py`, `app/practice_case.json`: practice case with known truth.
- `SOURCES.md`: source-to-mechanism table, invented fixtures and modeling limitations.
- `Dockerfile`: dependency environment.

## Use

```sh
docker build -t virtual-lab-vector-copy-number .
docker run --rm virtual-lab-vector-copy-number python3 /app/practice.py
```

The starter raises `NotImplementedError`; the practice runner reports metrics once `solver.py` submits a result.

## Execution evidence

In a local Docker run of the grader, a reference solution was accepted on 8 of 8 held-out cases (reward 1) and the unchanged starter scored 0. These are local simulator results. No blind solver trial or external review has been run, and nothing here is validated against laboratory data. A fixed acquisition plan is a valid solution, so the task does not require adaptivity.

## Reference and verifier

The reference solution, held-out cases, grader, comparison routes and correctness checks are kept in a separate private repository. They are not part of this repository.
