"""Synthetic digital PCR: vector copies per cell-equivalent and integration evidence.

Units and definitions
    cell-equivalent  1/303 ng of DNA, carrying 2 reference copies (a simplified fixture: 3.3 pg is
                     the haploid human genome mass, so real DNA has about 152 diploid genomes per ng;
                     the copies-per-cell-equivalent ratio does not depend on this constant)
    reference target RPP30-like locus, 2 copies per cell-equivalent, amplifies fully
    vector target    amplifies with efficiency eta (inhibition), see below
    VCN              vector copies per cell-equivalent. The estimand is the INTEGRATED VCN v:
                     integrated copies only. Early-passage DNA also carries non-integrated
                     copies n_ni; culture dilutes them to n_ni * 0.0625 but does not remove them,
                     so late-passage DNA measures v + 0.0625 * n_ni, not v.

One well
    20000 droplets are generated; n_acc of them (about 17500, varying per well) are accepted.
    Copies loaded: 2*GE reference, GE*(v + n_ni*decay) vector, with GE = input_ng*303*load.
    Per-droplet mean lambda = copies * efficiency / 20000 (+ background for the vector channel);
    a well returns positive counts ~ Binomial(n_acc, 1 - exp(-lambda)).
    efficiency (vector and spike) = 1 / (1 + k_lot * input_ng / 10); reference efficiency = 1.
    Early DNA: decay = 1. Late DNA (after `culture`): decay = 0.5**4 = 0.0625.
"""
import numpy as np

from harness import LabError

DROPLETS = 20000
GE_PER_NG = 303.0
SPIKE_COPIES = 4000.0
DECAY_LATE = 0.5 ** 4
INPUT_MENU = [1.25, 2.5, 5.0, 10.0, 20.0, 40.0]
ASSAYS = ["duplex", "spike", "ntc"]
WELL_BUDGET = 96
DAY_BUDGET = 21
CULTURE_DAYS = 7
STOCK_NG = 300.0
OVERHEAD_NG = 2.0
MAX_WELLS_PER_CALL = 24
INTEGRATION_THRESHOLD = 0.5


def make_case(seed):
    r = np.random.default_rng(seed)
    bg = 1e-5 if r.random() < 0.65 else float(r.uniform(2e-4, 6e-4))
    lots = []
    # lot 1 is always inhibited, lot 2 always has high non-integrated carryover and low integrated VCN
    for i in range(3):
        inhibited = True if i == 0 else bool(r.random() < 0.5)
        carry = True if i == 1 else bool(r.random() < 0.4)
        low = True if i == 1 else bool(r.random() < 0.2)
        v = float(r.uniform(0.1, 0.3)) if low else float(np.exp(r.uniform(np.log(0.8), np.log(4.0))))
        lots.append({"id": "L%d" % (i + 1), "v_int": v,
                     "n_ni": float(r.uniform(0.8, 2.5)) if carry else float(r.uniform(0.0, 0.15)),
                     "k_inhibitor": float(r.uniform(1.5, 4.0)) if inhibited else float(r.uniform(0.0, 0.2))})
    order = r.permutation(3)
    return {"seed": int(seed), "background_rate": bg, "lots": [lots[i] | {"id": "L%d" % (j + 1)} for j, i in enumerate(order)]}


class Lab:
    def __init__(self, case):
        self.case = case
        self.lots = {l["id"]: l for l in case["lots"]}
        self.wells_left = WELL_BUDGET
        self.days_used = 0
        self.stock = {(l, p): STOCK_NG for l in self.lots for p in ("early",)}
        self.late_ready = set()
        self.log = []
        self.well_counter = 0

    def public_config(self):
        return {"lots": sorted(self.lots), "assays": ASSAYS, "input_ng_menu": INPUT_MENU,
                "droplets_per_well": DROPLETS, "ge_per_ng": GE_PER_NG, "reference_copies_per_ge": 2,
                "spike_copies_per_well": SPIKE_COPIES, "late_passage_decay": DECAY_LATE,
                "well_budget": self.wells_left, "day_budget": DAY_BUDGET, "culture_days": CULTURE_DAYS,
                "stock_ng_per_passage": STOCK_NG, "overhead_ng_per_well": OVERHEAD_NG,
                "max_wells_per_call": MAX_WELLS_PER_CALL, "integration_threshold": INTEGRATION_THRESHOLD}

    def status(self):
        return {"wells_left": self.wells_left, "days_used": self.days_used,
                "stock_ng": {"%s:%s" % k: round(v, 2) for k, v in self.stock.items()},
                "late_available": sorted(self.late_ready), "wells_run": len(self.log)}

    def handle(self, action, msg):
        if action == "status":  # passive read: no state change
            return self.status()
        if action == "culture":
            lot = msg.get("lot")
            if lot not in self.lots:
                raise LabError("unknown lot %r" % (lot,))
            if lot in self.late_ready:
                raise LabError("late-passage DNA already prepared for %s" % lot)
            if self.days_used + CULTURE_DAYS > DAY_BUDGET:
                raise LabError("day budget exhausted")
            self.days_used += CULTURE_DAYS
            self.late_ready.add(lot)
            self.stock[(lot, "late")] = STOCK_NG
            return self.status()
        if action != "run_wells":
            raise LabError("unknown action %r" % action)
        lot, passage, assay = msg.get("lot"), msg.get("passage", "early"), msg.get("assay")
        n_w, ng = msg.get("wells", 1), msg.get("input_ng")
        if assay not in ASSAYS:
            raise LabError("unknown assay %r" % (assay,))
        if not isinstance(n_w, int) or isinstance(n_w, bool) or not 1 <= n_w <= MAX_WELLS_PER_CALL:
            raise LabError("wells must be an integer in [1, %d]" % MAX_WELLS_PER_CALL)
        if n_w > self.wells_left:
            raise LabError("well budget exhausted")
        if assay != "ntc":
            if lot not in self.lots:
                raise LabError("unknown lot %r" % (lot,))
            if passage not in ("early", "late"):
                raise LabError("passage must be early or late")
            if passage == "late" and lot not in self.late_ready:
                raise LabError("late-passage DNA for %s is not available; run culture first" % lot)
            if ng not in INPUT_MENU:
                raise LabError("input_ng must be one of %s" % INPUT_MENU)
            need = n_w * (ng + OVERHEAD_NG)
            if need > self.stock[(lot, passage)] + 1e-9:
                raise LabError("not enough %s DNA for %s" % (passage, lot))
            self.stock[(lot, passage)] -= need
        self.wells_left -= n_w
        wells = []
        for _ in range(n_w):
            rng = np.random.default_rng([self.case["seed"], 3, self.well_counter])
            self.well_counter += 1
            wells.append(self._well(rng, lot, passage, assay, ng))
        self.log.append({"lot": lot if assay != "ntc" else None, "passage": passage if assay != "ntc" else None,
                         "assay": assay, "input_ng": ng if assay != "ntc" else None, "wells": n_w})
        return {"wells": wells, "wells_left": self.wells_left}

    def _well(self, rng, lot, passage, assay, ng):
        bg = self.case["background_rate"]
        n_acc = int(np.clip(rng.normal(17500, 1200), 13000, 19800))
        if assay == "ntc":
            lam_vec, lam_ref, lam_spk = bg, 1e-6, None
        else:
            L = self.lots[lot]
            ge = ng * GE_PER_NG * float(np.exp(0.04 * rng.standard_normal()))
            eta = 1.0 / (1.0 + L["k_inhibitor"] * ng / 10.0)
            lam_ref = 2.0 * ge / DROPLETS
            vec_copies = ge * (L["v_int"] + L["n_ni"] * (DECAY_LATE if passage == "late" else 1.0))
            lam_vec = vec_copies * eta / DROPLETS + bg
            lam_spk = SPIKE_COPIES * eta / DROPLETS + bg
        def pos(lam):
            return int(rng.binomial(n_acc, 1.0 - np.exp(-lam)))
        if assay == "spike":
            return {"accepted": n_acc, "spike_positive": pos(lam_spk), "reference_positive": pos(lam_ref)}
        return {"accepted": n_acc, "vector_positive": pos(lam_vec), "reference_positive": pos(lam_ref)}


RULES = {"max_rel_error": 0.10, "min_abs_tolerance": 0.03, "max_ci_rel_width": 0.5, "min_ci_abs_width": 0.1,
         "ntc_valid_background_rate": 1e-4}


def evaluate(lab, result):
    out = {"lots": {}, "wells_used": WELL_BUDGET - lab.wells_left, "days_used": lab.days_used,
           "truth_ntc_valid": lab.case["background_rate"] <= RULES["ntc_valid_background_rate"]}
    lots = result.get("lots") if isinstance(result, dict) and isinstance(result.get("lots"), dict) else {}
    for lid, L in lab.lots.items():
        rec = {"true_integrated_vcn": L["v_int"], "true_nonintegrated_early": L["n_ni"], "k_inhibitor": L["k_inhibitor"]}
        try:
            r = lots[lid]
            est, lo, hi = float(r["integrated_vcn"]), float(r["ci"][0]), float(r["ci"][1])
            if not (np.all(np.isfinite([est, lo, hi])) and lo <= hi):
                raise ValueError("estimate and interval must be finite with lo <= hi")
            rec.update(estimate=est, ci=[lo, hi], rel_error=abs(est - L["v_int"]) / L["v_int"],
                       abs_error=abs(est - L["v_int"]), ci_covers=bool(lo <= L["v_int"] <= hi),
                       ci_width=hi - lo, integration_call_correct=bool(r["integration_supported"]) == (L["v_int"] >= INTEGRATION_THRESHOLD))
        except Exception:
            rec["estimate"] = None
        out["lots"][lid] = rec
    try:
        out["ntc_call_correct"] = bool(result["plate"]["ntc_valid"]) == out["truth_ntc_valid"]
    except Exception:
        out["ntc_call_correct"] = False
    return out


def check_rules(m):
    why = []
    for lid, r in m["lots"].items():
        if r.get("estimate") is None:
            why.append("%s: no valid estimate and interval" % lid)
            continue
        if r["abs_error"] > max(RULES["max_rel_error"] * r["true_integrated_vcn"], RULES["min_abs_tolerance"]):
            why.append("%s: estimate %.3f outside tolerance" % (lid, r["estimate"]))
        if not r["ci_covers"]:
            why.append("%s: interval misses the true value" % lid)
        if r["ci_width"] > max(RULES["max_ci_rel_width"] * r["true_integrated_vcn"], RULES["min_ci_abs_width"]):
            why.append("%s: interval too wide" % lid)
        if not r["integration_call_correct"]:
            why.append("%s: integration call wrong" % lid)
    if not m["ntc_call_correct"]:
        why.append("plate control validity call wrong")
    return not why, why
