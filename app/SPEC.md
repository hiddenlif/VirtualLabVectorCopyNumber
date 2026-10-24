# Digital PCR virtual lab: specification

All data are synthetic. The graded cases are other draws from `lab.make_case`.

## Protocol

Same JSON-lines exchange as in `labclient.py`: a `hello` line carries the config, the solver sends
`{"action": ..., ...}` lines and reads `{"ok": true, "result": ...}` or `{"ok": false, "error": ...}`,
and ends with `{"action": "submit", "result": {...}}`.

The grader runs `solver.py` as an unprivileged user with `/app` read-only and no network access; each case has a 900-second wall-clock limit and at most 5000 messages. Standard output carries only protocol lines, so write diagnostics to standard error. A message that is not a JSON object with an `action` field ends the session.

### Actions

- `status`: wells left, culture days used, DNA stock left per lot and passage, late-passage availability. Changes nothing.
- `culture` with `lot`: prepares late-passage DNA for that lot. Uses 7 of 21 virtual days and cannot be repeated for a lot. Late wells are refused before it.
- `run_wells` with `lot`, `passage` (`early` or `late`), `assay`, `input_ng` and `wells` (1 to 24):
  - `duplex` returns, per well, accepted droplets, `vector_positive` and `reference_positive`.
  - `spike` returns accepted droplets, `spike_positive` and `reference_positive`. The spike assay detects 4000 synthetic template copies added to the well and is affected by inhibition exactly like the vector assay.
  - `ntc` (no template; `lot`, `passage`, `input_ng` ignored) returns accepted droplets, `vector_positive` and `reference_positive`.
  - `input_ng` is one of 1.25, 2.5, 5, 10, 20, 40. A well uses `input_ng + 2` ng from the lot's stock for that passage (300 ng each) and one of 96 wells.

## Model

Per well, 20000 droplets are generated and n_acc of them are accepted, n_acc = Normal(17500, 1200) rounded down and clipped to [13000, 19800]. With GE = input_ng x 303 x load (load log-normal, 4% spread; 303 cell-equivalents per ng is a fixture of this model, see SOURCES.md):

- reference copies 2 x GE; vector copies GE x (v + n x d); d = 1 for early DNA, 0.5^4 = 0.0625 for late DNA
- efficiency eta = 1 / (1 + k x input_ng / 10) for the vector and spike assays; 1 for the reference
- lambda per droplet = copies x efficiency / 20000, plus the plate background rate on the vector and spike channels
- `ntc` wells: lambda is the plate background rate on the vector channel and 1e-6 on the reference channel
- positives ~ Binomial(accepted, 1 - exp(-lambda))

v is the integrated VCN (the estimand), n the non-integrated vector copies per cell-equivalent in early DNA, k the lot's inhibitor strength.
With lambda = -ln(1 - positive/accepted), a well estimates ((lambda_vector - background) / efficiency) / (lambda_reference / 2),
which is v + n for early DNA and v + 0.0625 n for late DNA.

## Submission

```
{"lots": {"L1": {"integrated_vcn": x, "ci": [lo, hi], "integration_supported": bool}, ...},
 "plate": {"ntc_valid": bool}}
```

## Scoring quantities (`lab.evaluate`, `lab.check_rules`)

Estimate error against the true integrated VCN, interval coverage and width, the integration call and the plate call, with the acceptance rules in `TASK.md`.
