# Figure Names and Legends for Submission

This document is generated from the current finalized figure outputs in `E:/GWAS/Final_Submission_Figures` and the corresponding plotting scripts (including `scripts/rebuild_supplement_s1.py`, `scripts/rebuild_supplement_s2.py`, `scripts/rebuild_supplement_s3.py`, `scripts/regen_s4_pub.py`, `scripts/rebuild_supplement_s5.py`, `scripts/regen_s6_enhanced.py`, `scripts/generate_s7_bridge.py`, `scripts/rebuild_supplement_s8.py`, `scripts/optimize_figure2_panels.py`, `create_figure4_final_multilevel.py`, and `scripts/Figure5_Main_UNIFIED.R`).

---

## Main Figures

**Figure 1. Study design and analytical workflow integrating genetic, epigenetic, transcriptomic, and clinical validation layers.**  
Schematic overview of the full study pipeline from GWAS instrument selection and Mendelian randomization to multi-omics prioritization and prognostic model development. The workflow defines data sources, statistical modules, candidate-gene prioritization strategy, and validation stages used for downstream translational interpretation.

**Figure 2. Causal inference and mediation architecture for lung cancer risk.**  
Panel A shows a forest plot of key univariable MR associations (odds ratios with 95% confidence intervals) for representative exposures and lung cancer outcomes. Panel B compares effect estimates before and after multivariable adjustment, highlighting direction and magnitude shifts under confounder control. Panel C presents a mediation path diagram centered on CRP, quantifying exposure-to-mediator and mediator-to-outcome effects and mediation proportion for core risk pathways.

**Figure 3. Colocalization and functional enrichment landscape of prioritized risk loci.**  
This figure integrates locus-level colocalization evidence and pathway-level enrichment results to summarize shared genetic regulation between molecular traits and disease risk. Panels jointly show evidence strength, locus prioritization, and biological-context interpretation for candidate mechanisms.

**Figure 4. Multi-layer biological validation across transcriptomic, immune, and proteomic dimensions.**  
Panel A summarizes differential expression patterns in prioritized genes. Panel B presents inflammatory/immune-association signals. Panel C integrates protein-level and pathway-level validation. Together, these panels provide convergent orthogonal support linking upstream genetic-epigenetic findings to biologically plausible downstream programs.

**Figure 5. Model-construction workflow and training-cohort survival stratification.**  
Panel A shows LASSO coefficient profiles across log(λ), illustrating variable shrinkage and selection dynamics. Panel B shows cross-validation error behavior with λmin/λ1se reference, supporting penalty stability. Panel C summarizes final risk-score coefficients (β) for selected CpG features and the score-construction form (Risk score = Σ[β × CpG methylation]). Panel D shows Kaplan-Meier survival stratification in the TCGA training cohort only (log-rank test), emphasizing in-cohort discrimination at the model-construction stage.

**Figure 6. Biological mechanisms linking risk score to tumor microenvironment remodeling.**  
Panel A (Immune features associated with risk score) displays selected representative immune populations (one algorithm per feature) ranked by effect direction and magnitude. Panel B (Functional pathway activation associated with risk score) summarizes pathways grouped into platelet-related, ECM/cytoskeleton, and inflammation modules, reflecting functional activation at the tumor level. This pathway-activity view is distinct from gene-level enrichment (Figure 3). Significance marks indicate pathway-level association strength (* P < 0.05, ** P < 0.01, *** P < 0.001). Together, the two panels support a coherent mechanism axis: immune alteration → platelet/ECM/inflammation programs → risk-associated tumor microenvironment remodeling.

**Figure 7. Robust prognostic validation of the M12 model across internal and external cohorts.**  
Panel A shows continuous Cox effect size (hazard ratio per SD increase in risk score). Panel B shows the dose-response relationship from restricted cubic spline (RCS), with low/intermediate/high-risk background zones and non-linearity statistics. Panel C shows the continuous risk-score distribution with median cutoff for orientation (not dichotomization inference). Panel D maps event distribution across the same risk ordering, linking risk continuum to survival outcomes. Together, these panels emphasize that risk operates continuously rather than as a purely dichotomous construct.

---

## Supplementary Figures

**Supplementary Figure S1. Sensitivity analyses for Mendelian randomization robustness.**  
Panel A is a funnel plot for the primary MR association (smoking initiation to lung cancer), showing SNP-level estimates against precision with IVW and MR-Egger reference lines. Panel B is a leave-one-out analysis (top influential instruments), where each row reports OR (95% CI) after removing one SNP; the global estimate line indicates that no single instrument materially changes inference.

**Supplementary Figure S2. Multivariable MR heatmap across exposures and lung cancer subtypes.**  
Heatmap of adjusted MVMR effect sizes (β) for exposure-outcome combinations after mutual adjustment. Color encodes effect direction/magnitude (risk-increasing versus protective), and symbols indicate significance thresholds (** P < 0.05, * P < 0.10). This figure highlights independent exposure signals under multivariable control.

**Supplementary Figure S3. Representative regional association plots for colocalization-prioritized loci.**  
Three panels display locus-level GWAS/eQTL overlap for MFAP2, MMEL1, and CDK11A, with PP.H4 labels and standardized visual annotation. Symbols distinguish GWAS versus eQTL variants, and gene regions are explicitly marked, enabling direct visual appraisal of shared-signal plausibility.

**Supplementary Figure S4. Full forest plot of MR estimates across exposure categories.**  
Two grouped sections summarize inflammatory and metabolic/lifestyle exposures with OR (95% CI), using significance-aware color coding (FDR-significant, nominal-significant, and non-significant). The figure provides complete transparency for effect direction, confidence bounds, and category-level patterning.

**Supplementary Figure S5. Diagnostic performance of the M12 model in the TCGA training cohort.**  
Panel A shows recurrence-free survival (RFS) stratification in GSE39279 under predefined risk grouping. Panel B shows overall-survival (OS) stratification in GSE30219. Together, these panels demonstrate external transportability of risk-group separation across cohorts with different endpoint definitions.

**Supplementary Figure S6. Extended model diagnostics beyond main-figure ROC and KM.**  
Panel A: 3-year calibration in TCGA with binomial confidence intervals and calibration slope statistics. Panel B: decision curve analysis in TCGA, comparing net benefit of model versus treat-all and treat-none across threshold probabilities. Panel C: TCGA risk-score distribution around the median cutoff. Panel D: external 3-year calibration in GSE30219. This supplementary figure focuses on calibration, clinical utility, and score-distribution diagnostics without duplicating main ROC/KM content.

**Supplementary Figure S7. Multi-layer bridge from genetic variation to clinical prognosis.**  
Conceptual bridge diagram linking five layers: Genetic (GWAS/mQTL), Epigenetic (regulatory methylation), Transcriptome (RNA-seq/eQTL with MFAP2 as key bridge), Biological (protein and immune context), and Clinical (survival validation). The diagram encodes directed translational flow from causal inference to functional validation and clinical translation.

**Supplementary Figure S8. Risk stratification and survival-pattern visualization in TCGA.**  
Panel 1 shows ordered risk-score ranking with median cutoff and low/high-risk partitions. Panel 2 shows survival-time scatter aligned to risk ordering with event-status encoding. Panel 3 shows standard Kaplan-Meier curves on time axis for low- and high-risk groups, demonstrating reduced survival probability in the high-risk stratum. This figure links score distribution, event distribution, and survival probability in one coherent diagnostic narrative.

---

## Figure File Mapping (Final Submission Folder)

### Main figures
- Figure 1: `Main_Figures/Figure1.tif`
- Figure 2: `Main_Figures/Figure2.tif` (components: `Figure2A_Forest_Plot.png`, `Figure2B_MVMR.png`, `Figure2C_Mediation.png`)
- Figure 3: `Main_Figures/Figure3.tif`
- Figure 4: `Main_Figures/Figure4.tif` (component: `Figure4_RNA_Immune.png`)
- Figure 5: `Main_Figures/Figure5.tif` (components: `Figure5_Clinical_Model.pdf/.png/.tif`)
- Figure 6: `Main_Figures/Figure6.tif` (components: `Figure6.pdf/.png/.tif`; panels: immune infiltration + pathway activity)
- Figure 7: `Main_Figures/Figure7.tif`

### Supplementary figures
- S1: `Supplementary_Figures/S1.tif` (`S1_MR_Sensitivity_Analysis.pdf/.png`)
- S2: `Supplementary_Figures/S2.tif` (`S2_MVMR_Heatmap.pdf/.png`)
- S3: `Supplementary_Figures/S3.tif` (`S3_Colocalization_Prioritized.pdf/.png`)
- S4: `Supplementary_Figures/S4.tif` (`S4_MR_Full_Forest_Plot.pdf/.png`)
- S5: `Supplementary_Figures/S5.tif` (`S5_External_Validation_KM.pdf/.png`)
- S6: `Supplementary_Figures/S6.tif` (`S6_ROC_Calibration.pdf/.png`)
- S7: `Supplementary_Figures/S7.tif` (`S7_Bridge_Diagram.pdf/.png`)
- S8: `Supplementary_Figures/S8.tif` (`S8_Risk_Score_Distribution.pdf/.png`)
