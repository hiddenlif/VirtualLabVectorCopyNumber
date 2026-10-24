# Sources, mechanisms and limitations

| Modeled mechanism | Source | Status in this simulator |
|---|---|---|
| Partition occupancy follows a Poisson model; concentration from -ln(1 - positive fraction) | Hindson BJ et al. High-throughput droplet digital PCR system for absolute quantitation of DNA copy number. Anal Chem 83:8604-8610 (2011) | Sourced mechanism |
| Copy number per genome from a duplex target/reference measurement with a two-copy reference | Hindson et al. 2011 (above) | Sourced practice |
| Reporting of partitions, controls and no-template controls | dMIQE Group. The digital MIQE guidelines update. Clin Chem 66:1012-1029 (2020) | Used for the record and control set |
| Inhibitors reduce amplification; digital PCR is more tolerant than qPCR but not immune | Dingle TC et al. Tolerance of droplet-digital PCR versus real-time quantitative PCR to inhibitory substances. Clin Chem 59:1670-1672 (2013) | Sourced qualitatively; the efficiency curve 1/(1 + k x ng/10) and its target-specificity are invented |
| Non-integrated vector DNA is lost as cells divide, while integrated copies persist | General mechanism | Modeling assumption; the dilution factor 0.5^4 after culture is invented |

Estimand: integrated vector copies per cell-equivalent (one diploid genome, 303 per ng). This is not copies per transduced cell (that would need the transduced fraction), and the integration call here is a threshold on the late-passage estimate within this model, not a molecular proof of integration.

Limitations: no droplet-size variation, no rain or threshold-setting step (counts are given directly), a fixed spike-in, no locus-specific or copy-number variation in the reference, and lots are single specimens.
