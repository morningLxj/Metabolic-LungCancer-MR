**2. Methods**

**2.1 Study design and analytical framework**

We established a causal-to-translational multi-omics analytical
framework to identify biologically informed epigenetic biomarkers and to
assess their clinical relevance in lung squamous cell carcinoma (LUSC).
The workflow integrated Mendelian randomization (MR), methylation
quantitative trait locus (mQTL) mapping, Bayesian colocalization, and
multi-layer omics validation, thereby enabling stepwise prioritization
from upstream genetic liability to downstream prognostic application. A
schematic overview of the study design is provided in Figure 1.

**2.2 Genetic liability analysis using Mendelian randomization**

Two-sample MR analyses were conducted to estimate the putative causal
effects of metabolic, inflammatory, and lifestyle-related exposures on
lung cancer outcomes. Instrumental variables were selected at
genome-wide significance (P \< 5 × 10⁻⁸) and pruned for linkage
disequilibrium (r² \< 0.001) to ensure independence among retained
variants.

The primary causal estimate was derived using the inverse-variance
weighted (IVW) approach:

> *β̂\_IVW = ∑ \[ w_i · (β_Y,i / β_X,i) \] / ∑ w_i*

where w_i denotes the inverse variance of the SNP-outcome association,
defined as 1 / (SE_Y,i)².

To evaluate the robustness of the IVW estimates and to assess potential
horizontal pleiotropy, complementary sensitivity analyses were performed
using MR-Egger regression, the weighted median estimator, and MR-PRESSO.
The MR-Egger model was specified as:

> *β_Y,i = α + β_MR · β_X,i + ε_i*

**2.3 Multivariable and mediation MR**

Multivariable MR (MVMR) was performed to estimate the independent
effects of correlated exposures within a joint modeling framework. To
further explore biologically plausible intermediate pathways, two-step
MR analyses were undertaken by separately modeling exposure-mediator and
mediator-outcome associations, with particular emphasis on the potential
mediating role of inflammatory markers in metabolic risk transmission.
Indirect effects were interpreted according to the
product-of-coefficients framework.

**2.4 Epigenetic triangulation framework**

To prioritize CpG sites that potentially translate inherited liability
into tumor-level molecular regulation, we implemented an epigenetic
triangulation framework integrating mQTL mapping, MR directionality, and
Bayesian colocalization analysis. mQTL datasets were first used to
identify genetically regulated CpG loci. Colocalization analysis was
then applied to evaluate whether the exposure-associated and
methylation-associated signals were likely attributable to a shared
causal variant; posterior probability for hypothesis 4 (PP.H4) \> 0.30
was considered suggestive evidence of colocalization. CpG sites were
retained only when the inferred directions of effect were concordant
across the MR and mQTL layers.

**2.5 Multi-omics integration and functional annotation**

Prioritized CpG sites and their annotated genes were subsequently
evaluated across transcriptomic, proteomic, and immune-related datasets
to strengthen biological interpretation and translational relevance.
Functional enrichment analyses were performed to identify pathways
associated with the candidate genes, with particular emphasis on
coagulation, platelet activation, extracellular matrix organization, and
immune-related processes implicated in tumor microenvironment
remodeling.

**2.6 Prognostic model construction (M12 Signature)**

DNA methylation profiles and corresponding clinical follow-up data were
obtained from the TCGA-LUSC cohort for prognostic model development.
Feature selection was performed using a penalized Cox proportional
hazards model with LASSO regularization and 10-fold cross-validation.
Candidate signatures were subsequently compared using a composite
selection score, and the final 12-CpG signature (M12_Top12) was selected
on the basis of overall parsimony-performance balance. The Cox
proportional hazards model was defined as:

> *h(t\|X) = h_0(t) exp(β_1X_1 + \... + β_pX_p)*

For each patient, the prognostic risk score was calculated from the
multivariable Cox model using z-scored CpG features, where Z(CpG)
represents the standardized methylation value of a given CpG site. The
general form of the risk score was:

> *Risk Score = ∑ (β_j × Z(CpG_j))*

The final M12 signature was specified as:

> *Risk Score = (0.210 × Z\[cg06238570\]) + (0.140 × Z\[cg07318658\]) +
> (0.104 × Z\[cg13144823\]) - (0.222 × Z\[cg14526939\]) - (0.294 ×
> Z\[cg07151966\]) + (0.041 × Z\[cg12672785\]) + (0.075 ×
> Z\[cg05801080\]) - (0.205 × Z\[cg07154958\]) + (0.253 ×
> Z\[cg16098170\]) + (0.260 × Z\[cg09010499\]) - (0.230 ×
> Z\[cg05984948\]) - (0.152 × Z\[cg00646216\])*

**2.7 Model evaluation and validation**

Model performance was comprehensively evaluated using Kaplan-Meier
survival analysis, the concordance index (C-index), time-dependent
receiver operating characteristic (ROC) analysis, calibration curves,
and decision curve analysis (DCA). In the TCGA-LUSC training cohort,
patients were dichotomized into high- and low-risk groups according to
the median training-cohort risk score. External validation was performed
in two independent microarray cohorts (GSE39279 and GSE30219), in which
optimal cutoffs were determined using the maxstat (maximally selected
rank statistics) procedure. Given the heterogeneity in endpoint
definition across validation cohorts, namely recurrence-free survival in
GSE39279 and overall survival in GSE30219, external validation primarily
focused on discriminative performance, whereas calibration and DCA were
principally assessed in the TCGA training cohort.

**2.8 Statistical analysis**

All statistical analyses were performed using R software (version
4.5.2). Unless otherwise specified, all statistical tests were
two-sided, and P-values \< 0.05 were considered statistically
significant. Where appropriate, results were interpreted in conjunction
with multiple-testing considerations and consistency across analytical
layers.
