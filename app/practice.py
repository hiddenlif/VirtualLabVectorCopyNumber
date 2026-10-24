"""Run solver.py against the disclosed practice case and print metrics with known truth."""
import json
import sys
from pathlib import Path

import harness
import lab

HERE = Path(__file__).parent


def main():
    case = json.loads((HERE / "practice_case.json").read_text())
    L = lab.Lab(case)
    out = harness.run_session(L, str(HERE / "solver.py"), str(HERE))
    if out["error"]:
        print("session failed:", out["error"])
        return 1
    m = lab.evaluate(L, out["result"])
    for lid, r in m["lots"].items():
        if r.get("estimate") is None:
            print(lid, "no estimate")
            continue
        print("%s true=%.3f est=%.3f rel_err=%.3f ci_covers=%s integration_call_correct=%s (non-integrated early %.2f, inhibitor k %.2f)"
              % (lid, r["true_integrated_vcn"], r["estimate"], r["rel_error"], r["ci_covers"],
                 r["integration_call_correct"], r["true_nonintegrated_early"], r["k_inhibitor"]))
    ok, why = lab.check_rules(m)
    print("plate control call correct:", m["ntc_call_correct"])
    print("meets the published rules on this case:", ok, why)
    print("wells used:", m["wells_used"], "culture days used:", m["days_used"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
