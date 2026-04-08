**Table 1. GWAS Datasets Used in Mendelian Randomization Analyses**

| **Trait** | **Category** | **Role** | **GWAS ID** | **Sample size** | **No. of SNPs** |
|:--|:--|:--|:--|--:|--:|
| Apolipoprotein A1 levels | Metabolic | Exposure | ebi-a-GCST90025955 | 398,508 | 284 |
| Apolipoprotein B levels | Metabolic | Exposure | ebi-a-GCST90025952 | 435,744 | 189 |
| Fasting glucose | Metabolic | Exposure | ebi-a-GCST90002232 | 200,622 | 69 |
| Fasting insulin | Metabolic | Exposure | ebi-a-GCST90002238 | 151,013 | 38 |
| HDL cholesterol | Metabolic | Exposure | ieu-b-109 | 403,943 | 356 |
| Insulin-like growth factor 1 levels | Metabolic | Exposure | ebi-a-GCST90025989 | 435,516 | 395 |
| LDL cholesterol | Metabolic | Exposure | ieu-b-110 | 440,546 | 179 |
| Ratio of apolipoprotein B to apolipoprotein A1 | Metabolic | Exposure | ebi-a-GCST90092810 | 115,082 | 72 |
| Remnant cholesterol | Metabolic | Exposure | ebi-a-GCST90092943 | 115,082 | 53 |
| Body mass index (BMI) | Metabolic | Exposure | ieu-b-40 | 681,275 | 501 |
| Triglycerides | Metabolic | Exposure | ieu-b-111 | 441,016 | 313 |
| C-reactive protein levels | Inflammatory | Mediator | ebi-a-GCST90029070 | 575,531 | 264 |
| Interleukin-6 levels | Inflammatory | Exposure | ebi-a-GCST90012005 | 21,758 | 2 |
| Interleukin-6 receptor subunit alpha levels | Inflammatory | Exposure | ebi-a-GCST90012025 | 21,758 | 5 |
| Diastolic blood pressure | Blood pressure | Exposure | ieu-b-39 | 757,601 | 460 |
| Systolic blood pressure | Blood pressure | Exposure | ieu-b-38 | 757,601 | 461 |
| Smoking initiation | Lifestyle | Exposure | ieu-b-4877 | 607,291 | 93 |
| LUAD | Outcome | Outcome | ieu-a-984 | 65,864 | 14 |
| Overall lung cancer | Outcome | Outcome | ebi-a-GCST004748 | 85,716 | 15 |
| LUSC | Outcome | Outcome | ieu-a-989 | 62,467 | 9 |

*European-ancestry GWAS sources; columns retained for MR interpretation only.*

**Table 2. Causal Associations Between Exposures and Lung Cancer Outcomes (IVW MR)**

| **Exposure** | **Outcome** | **OR (95% CI)** | **P-value** | **Effect direction** | **Significance** | **Subtype specificity (LUSC > LUAD)** |
|:--|:--|:--|:--:|:--:|:--:|:--:|
| Smoking initiation | Overall lung cancer | 1.639 (1.421–1.890) | **\<0.001** | Risk-increasing | *** | NA |
| Smoking initiation | LUAD | 1.447 (1.227–1.706) | **\<0.001** | Risk-increasing | *** | No |
| Smoking initiation | LUSC | 1.863 (1.497–2.317) | **\<0.001** | Risk-increasing | *** | Yes |
| Body mass index (BMI) | Overall lung cancer | 1.211 (1.102–1.330) | **\<0.001** | Risk-increasing | *** | NA |
| Body mass index (BMI) | LUSC | 1.428 (1.245–1.637) | **\<0.001** | Risk-increasing | *** | Yes |
| C-reactive protein (CRP) | LUSC | 1.188 (1.051–1.343) | 0.006 | Risk-increasing | ** | Yes |
| Alcohol drinks/week | Overall lung cancer | 1.560 (1.174–2.072) | 0.002 | Risk-increasing | ** | NA |
| Alcohol drinks/week | LUAD | 1.673 (1.145–2.445) | 0.008 | Risk-increasing | ** | No |
| Diastolic blood pressure | LUSC | 0.985 (0.972–0.998) | 0.021 | Protective |  | Yes |
| Systolic blood pressure | LUSC | 0.992 (0.984–0.999) | 0.026 | Protective |  | Yes |
| Fasting insulin | LUAD | 1.567 (1.052–2.334) | 0.027 | Risk-increasing |  | No |

**Table 3. Integrated Multi-Omics Evidence for Prioritized Genes**

| **Gene** | **Genetic (MR)** | **Epigenetic (mQTL)** | **Coloc (PP.H4)** | **Transcript** | **Protein** | **Immune** | **Clinical** | **Biological role** | **Evidence level** | **Key driver** |
|:--|:--:|:--:|--:|:--:|:--:|:--:|:--:|:--|:--:|:--:|
| MFAP2 | ✔ | ✔ | 0.305 | Decrease | Increase | ✔ | ✔ | ECM remodeling and immune interaction | ★★★★★ | Yes |
| CDK11A | ✔ | ✔ | 0.612 | Not significant | Increase | ✔ | ✔ | Cell-cycle regulation | ★★★☆☆ | No |
| WRAP73 | ✔ | ✔ | 0.489 | Decrease | Decrease | ✔ | ✖ | DNA repair | ★★★☆☆ | No |

*Mechanism order: Genetic → Epigenetic → Transcript → Protein → Immune → Clinical. Evidence definition: ★★★★★ = ≥5 concordant layers; ★★★☆☆ = 3–4 concordant layers. All evidence layers were weighted equally without hierarchical prioritization.*

**Table 4. Final M12 CpG Prognostic Signature Definition**

| **CpG ID** | **Gene annotation** | **β** | **HR (95% CI)** | **P-value** | **Direction** |
|:--|:--|--:|:--|:--:|:--:|
| cg06238570 | BRWD1 | 0.210 | 1.23 (1.04–1.46) | 0.015 | Risk-increasing |
| cg07318658 | MTMR4 | 0.140 | 1.15 (0.99–1.34) | 0.077 | Risk-increasing |
| cg13144823 | CXXC5 | 0.104 | 1.11 (0.90–1.36) | 0.319 | Risk-increasing |
| cg14526939 | TOR1AIP1 | -0.222 | 0.80 (0.66–0.97) | 0.024 | Protective |
| cg07151966 | intergenic | -0.294 | 0.75 (0.63–0.88) | 4.58 × 10⁻⁴ | Protective |
| cg12672785 | EED | 0.041 | 1.04 (0.87–1.25) | 0.656 | Risk-increasing |
| cg05801080 | LNP1; TOMM70A | 0.075 | 1.08 (0.90–1.30) | 0.428 | Risk-increasing |
| cg07154958 | intergenic | -0.205 | 0.81 (0.62–1.07) | 0.144 | Protective |
| cg16098170 | JSRP1 | 0.252 | 1.29 (1.05–1.57) | 0.013 | Risk-increasing |
| cg09010499 | BDKRB2 | 0.260 | 1.30 (1.06–1.58) | 0.010 | Risk-increasing |
| cg05984948 | intergenic | -0.229 | 0.80 (0.66–0.96) | 0.017 | Protective |
| cg00646216 | intergenic | -0.152 | 0.86 (0.74–0.99) | 0.041 | Protective |

*Model: M12_Top12 derived from penalized Cox selection and finalized by parsimony-performance balance. HR is reported from multivariable Cox with z-scored CpG features.*

**Table 5. Performance of the M12 Model Across Cohorts**

| **Dataset** | **Endpoint (OS / RFS)** | **N** | **HR (high vs low risk, Cox model, 95% CI)** | **P-value** | **C-index** | **AUC (3-year)** | **Calibration slope (3Y)** | **Calibration R² (3Y)** | **Cutoff method** |
|:--|:--:|--:|:--|:--:|:--:|:--:|:--:|:--:|:--|
| TCGA (training) | Overall survival (OS) | 364 | 2.88 (2.05–4.05) | P = 2.32 × 10⁻¹⁰ | 0.664 | 0.685 | 1.227 | 0.843 | Median cutoff derived from the training cohort |
| TCGA (internal validation) | Overall survival (OS) | 110 | 2.66 (1.43–4.94) | P = 0.001 | Not available | Not available | Not available | Not available | Median cutoff derived from the training cohort |
| GSE39279 (external) | Recurrence-free survival (RFS) | 122 | 2.10 (1.10–4.03) | P = 0.022 | Not available | Not available | Not available | Not available | Maxstat-derived optimal cutoff (external cohorts) |
| GSE30219 (external) | Overall survival (OS) | 274 | 1.57 (1.17–2.10) | P = 0.002 | Not available | Not available | Not available | Not available | Maxstat-derived optimal cutoff (external cohorts) |

*Figure/Table consistency note: all metrics in this table are aligned with the final M12 model outputs used in Figure 5–7 and Supplementary Figure S6. HR in this table represents high versus low risk groups; HR per SD is reported in Supplementary Table S8.*
