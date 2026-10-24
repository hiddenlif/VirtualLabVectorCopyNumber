# Sources, mechanisms and limitations

| Modeled mechanism | Source | Status in this simulator |
|---|---|---|
| Partition occupancy follows a Poisson model; concentration from -ln(1 - positive fraction) | Hindson BJ et al. High-throughput droplet digital PCR system for absolute quantitation of DNA copy number. Anal Chem 83:8604-8610 (2011) | Sourced mechanism |
| Copy number per genome from a duplex target/reference measurement with a two-copy reference | Hindson et al. 2011 (above) | Sourced practice |
| Reporting of partitions, controls and no-template controls | dMIQE Group. The digital MIQE guidelines update. Clin Chem 66:1012-1029 (2020) | Used for the record and control set |
| Inhibitors reduce amplification; digital PCR is more tolerant than qPCR but not immune | Dingle TC et al. Tolerance of droplet-digital PCR versus real-time quantitative PCR to inhibitory substances. Clin Chem 59:1670-1672 (2013) | Sourced qualitatively; the efficiency curve 1/(1 + k x ng/10) and its target-specificity are invented |
| Non-integrated vector DNA is lost as cells divide, while integrated copies persist | General mechanism | Modeling assumption; the dilution factor 0.5^4 after culture is invented |

Estimand: integrated vector copies per cell-equivalent. A cell-equivalent here is 1/303 ng of DNA carrying two reference copies; this is a simplified fixture (3.3 pg is the haploid human genome mass, so real human DNA holds about 152 diploid genomes per ng), and the copies-per-cell-equivalent ratio does not depend on it. This is not copies per transduced cell (that would need the transduced fraction). Late-passage DNA still carries 0.0625 of the early non-integrated copies, so the integrated VCN is not the late-passage VCN. The integration call is a threshold (0.5) on the integrated VCN within this model, not a molecular proof of integration.

Limitations: no droplet-size variation, no rain or threshold-setting step (counts are given directly), a fixed spike-in, no locus-specific or copy-number variation in the reference, and lots are single specimens.
