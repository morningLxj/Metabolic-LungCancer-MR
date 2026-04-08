**From genetic liability to epigenetic risk: a multi-omics framework for
prognostic stratification in lung squamous cell carcinoma**

Xiongjie Li¹†, Fengyue Zhang²†, Zhenyao Wu³, Xuan Xu³, Xiaoyan Zhang⁴,
Xianghui Wang¹\*

¹ Department of Cardiothoracic and Vascular Surgery, Jingzhou Hospital
Affiliated to Yangtze University, No. 26 Chuyuan Avenue, Jingzhou
District, Jingzhou 434020, Hubei, China

² Department of Gastroenterology, Qianjiang Jianghan Oilfield General
Hospital, Qianjiang 433100, Hubei, China

³ Department of Cardiothoracic, Thyroid and Breast Surgery, Qianjiang
Jianghan Oilfield General Hospital, Qianjiang 433100, Hubei, China

⁴ Department of Clinical Laboratory, Qianjiang Jianghan Oilfield General
Hospital, Qianjiang 433100, Hubei, China

*† These authors contributed equally to this work.*

\*Corresponding author: Xianghui Wang (wangxianghui8104@163.com)

Xiongjie Li (ORCID: 0009-0005-9507-6829; Email: surobert518@163.com)

**Abstract**

**Background:** Metabolic and inflammatory liabilities have been
implicated in lung cancer risk, yet the molecular mechanisms linking
upstream exposures to tumor behavior and clinical outcomes remain
unclear. Integrative frameworks bridging causal inference and
translational biomarkers are needed.

**Methods:** We developed a causal-to-translational multi-omics
framework integrating Mendelian randomization (MR), methylation
quantitative trait loci (mQTL), Bayesian colocalization, and multi-layer
omics validation. Genetically informed CpG sites were prioritized and
incorporated into a methylation-based prognostic model using penalized
Cox regression. Internal validation was performed in TCGA-LUSC, with
external validation in GSE39279 and GSE30219.

**Results:** MR analyses revealed a subtype-specific
metabolic--inflammatory liability architecture, with stronger effects
observed in lung squamous cell carcinoma (LUSC) than lung
adenocarcinoma. Two-step MR suggested C-reactive protein (CRP) as a
potential mediator linking metabolic traits to LUSC risk. Epigenetic
triangulation prioritized MFAP2-associated CpG sites supported by
concordant mQTL, MR, and colocalization evidence. Multi-omics analyses
converged on a coherent tumor program characterized by platelet
activation, coagulation, and extracellular matrix remodeling. A 12-CpG
methylation signature (M12) demonstrated stable prognostic performance
in the TCGA cohort (C-index = 0.664; 3-year AUC = 0.685), with good
calibration (slope = 1.23; R² = 0.84). External validation confirmed
significant survival stratification across independent cohorts.

**Conclusions:** This study establishes a causal-to-translational
framework linking genetic liability to biologically grounded epigenetic
markers. The M12 CpG-based signature enables clinically relevant risk
stratification in LUSC and provides a scalable strategy for
mechanistically informed biomarker development.

**Keywords:** Mendelian randomization; multi-omics integration; DNA
methylation; epigenetics; MFAP2; lung squamous cell carcinoma; risk
stratification; prognostic model

**1. Introduction**

Lung cancer remains a leading cause of cancer-related mortality
worldwide \[33,36\]. Recurrence after curative-intent resection is a
major determinant of outcome in non-small cell lung cancer, including
LUSC \[27,28\]. Current clinicopathological risk stratification,
including TNM staging, does not fully capture inter-individual
variability in recurrence risk \[14\]. Despite advances in treatment,
the identification of patients at high risk of relapse remains limited
in LUSC \[23\], particularly in the context of the heterogeneous
molecular landscape of the disease \[5,37\].

Accumulating evidence suggests that systemic metabolic and inflammatory
liabilities contribute to lung cancer development and progression
\[15,25\]. However, such associations are often difficult to interpret
because of residual confounding, reverse causation, and population
heterogeneity. Mendelian randomization (MR) has emerged as a rigorous
framework for strengthening causal inference by leveraging genetic
variants as instrumental variables \[4,18\]. Nevertheless, etiological
signals derived from MR alone are insufficient for clinical translation
because they do not directly resolve tumor-level molecular mechanisms or
prognostic biomarkers.

A key missing component is a molecular "translation layer" linking
upstream liabilities to stable and clinically measurable tumor features.
DNA methylation is well suited to this role because of its biological
relevance, relative stability, and mechanistic coupling to cellular
metabolism \[9,22\]. Its translational appeal is further supported by
established links between epigenetic remodeling, cancer biology, and
biomarker development \[20,42\]. Partial genetic regulation through
methylation quantitative trait loci (mQTL) further enables causal
anchoring \[11,26\]. Integrating MR with mQTL mapping and Bayesian
colocalization therefore provides an opportunity to prioritize
epigenetic signals that are both genetically anchored and functionally
plausible \[13,29\], while improving the specificity of candidate
biomarker discovery through shared-signal evaluation \[30,43\].

In this study, we propose a causal-to-translational multi-omics
framework that integrates Mendelian randomization, mQTL mapping,
Bayesian colocalization, and multi-layer omics validation to identify
biologically grounded epigenetic signals \[17\]. We further translate
these genetically informed CpG sites into a parsimonious and externally
validated methylation-based prognostic model (M12) for risk
stratification in lung squamous cell carcinoma. By bridging causal
inference with tumor-level molecular features and clinical outcomes,
this approach aims to establish a mechanistically informed pathway from
genetic liability to prognostic biomarkers in precision oncology.

**2. Methods**

**2.1 Study design and analytical framework**

We established a causal-to-translational multi-omics analytical
framework to identify biologically informed epigenetic biomarkers and to
assess their clinical relevance in lung squamous cell carcinoma (LUSC).
The workflow integrated Mendelian randomization (MR) with methylation
quantitative trait locus (mQTL) mapping to anchor upstream liabilities
to CpG-level regulation \[4,13\]. Bayesian colocalization and
multi-layer omics validation were then used to refine candidate signals
and support downstream prognostic application \[17,30\]. A schematic
overview of the study design is provided in Figure 1.

***\[Insert Figure 1 here\]***

**2.2 Genetic liability analysis using Mendelian randomization**

Two-sample MR analyses were conducted to estimate the putative causal
effects of metabolic, inflammatory, and lifestyle-related exposures on
lung cancer outcomes. Instrumental variable selection followed
established MR guidance, with variants retained at genome-wide
significance (P \< 5 × 10⁻⁸) and pruned for linkage disequilibrium (r²
\< 0.001) to ensure independence among retained variants \[4\]. Summary
statistics were obtained from large-scale public GWAS consortia and
harmonized through established MR-oriented resources \[1,18,39\].

The primary causal estimate was derived using the inverse-variance
weighted (IVW) approach:

*β̂\_IVW = ∑ \[ w_i · (β_Y,i / β_X,i) \] / ∑ w_i*

where w_i denotes the inverse variance of the SNP-outcome association,
defined as 1 / (SE_Y,i)².

To evaluate the robustness of the IVW estimates and to assess potential
horizontal pleiotropy, complementary sensitivity analyses were performed
using MR-Egger regression, the weighted median estimator, and MR-PRESSO
\[3,41\]. The MR-Egger model was specified as:

*β_Y,i = α + β_MR · β_X,i + ε_i*

Detailed MR sensitivity analyses and supporting results are provided in
Tables S1--S2.

**2.3 Multivariable and mediation MR**

Multivariable MR (MVMR) was performed to estimate the independent
effects of correlated exposures within a joint modeling framework
\[32\]. To further explore biologically plausible intermediate pathways,
two-step MR analyses were undertaken by separately modeling
exposure-mediator and mediator-outcome associations, with particular
emphasis on the potential mediating role of inflammatory markers in
metabolic risk transmission \[6,29\]. Indirect effects were interpreted
according to the product-of-coefficients framework. Extended MVMR and
mediation analyses are shown in Tables S2--S3.

**2.4 Epigenetic triangulation framework**

To prioritize CpG sites that potentially translate inherited liability
into tumor-level molecular regulation, we implemented an epigenetic
triangulation framework integrating mQTL mapping, MR directionality, and
Bayesian colocalization analysis. mQTL resources were first used to
identify genetically regulated CpG loci \[11,26\]. These candidates were
then filtered according to concordant MR directionality within a
triangulation framework \[13,29\]. Colocalization analysis was
subsequently applied to evaluate whether the exposure-associated and
methylation-associated signals were likely attributable to a shared
causal variant \[30,43\]; posterior probability for hypothesis 4 (PP.H4)
\> 0.30 was considered suggestive evidence of colocalization. CpG sites
were retained only when the inferred directions of effect were
concordant across the MR and mQTL layers. Full colocalization results
are provided in Table S4, with regional plots shown in Figure S3.

**2.5 Multi-omics integration and functional annotation**

Prioritized CpG sites and their annotated genes were subsequently
evaluated across transcriptomic, proteomic, and immune-related datasets
to strengthen biological interpretation and translational relevance
\[17\]. Transcriptomic processing and differential analyses were aligned
with established RNA-seq and microarray workflows \[24,31\], with batch
correction performed using empirical Bayes adjustment where required
\[19\]. Functional enrichment analyses were performed using established
pathway resources, including Reactome and GSVA \[10,16\], with
tissue-level contextualization from GTEx \[38\] and particular emphasis
on tumor immune microenvironment remodeling \[2\]. Enrichment
visualization and annotation were implemented with clusterProfiler
\[45\]. CpG-to-gene mapping details are summarized in Table S5.

**2.6 Prognostic model construction (M12 Signature)**

DNA methylation profiles and corresponding clinical follow-up data were
obtained from the TCGA-LUSC cohort for prognostic model development
\[37\]. Feature selection was performed using a penalized Cox
proportional hazards model with LASSO regularization and 10-fold
cross-validation \[34,44\]. Model development and reporting were aligned
with contemporary principles for biomedical prediction modeling
\[7,35\]. Candidate signatures were subsequently compared using a
composite selection score, and the final 12-CpG signature (M12_Top12)
was selected on the basis of overall parsimony-performance balance. The
Cox proportional hazards model was defined as:

*h(t\|X) = h_0(t) exp(β_1X_1 + \... + β_pX_p)*

For each patient, the prognostic risk score was calculated from the
multivariable Cox model using z-scored CpG features, where Z(CpG)
represents the standardized methylation value of a given CpG site. The
general form of the risk score was:

*Risk Score = ∑ (β_j × Z(CpG_j))*

The final M12 signature was specified as:

*Risk Score = (0.210 × Z\[cg06238570\]) + (0.140 × Z\[cg07318658\]) +
(0.104 × Z\[cg13144823\]) - (0.222 × Z\[cg14526939\]) - (0.294 ×
Z\[cg07151966\]) + (0.041 × Z\[cg12672785\]) + (0.075 ×
Z\[cg05801080\]) - (0.205 × Z\[cg07154958\]) + (0.253 ×
Z\[cg16098170\]) + (0.260 × Z\[cg09010499\]) - (0.230 ×
Z\[cg05984948\]) - (0.152 × Z\[cg00646216\])*

Full model coefficients and feature-selection procedures are provided in
Tables S6--S7.

**2.7 Model evaluation and validation**

Model performance was comprehensively evaluated using Kaplan-Meier
survival analysis, nonparametric time-to-event estimation under
censoring, the concordance index (C-index), time-dependent receiver
operating characteristic (ROC) analysis, calibration curves, and
decision curve analysis (DCA) \[40\]. In the TCGA-LUSC training cohort,
patients were dichotomized into high- and low-risk groups according to
the median training-cohort risk score. The TCGA cohort was randomly
split into training and internal validation sets using stratified
sampling based on survival status. External validation was performed in
two independent microarray cohorts (GSE39279 and GSE30219), in which
optimal cutoffs were determined using the maxstat (maximally selected
rank statistics) procedure. Given the heterogeneity in endpoint
definition across validation cohorts, namely recurrence-free survival in
GSE39279 and overall survival in GSE30219, external validation primarily
focused on discriminative performance, whereas calibration and DCA were
principally assessed in the TCGA training cohort. Model development and
reporting were aligned with contemporary prediction-model guidance
\[7,35,44\]. Additional diagnostic analyses are shown in Figure S5 and
Tables S8--S9.

**2.8 Statistical analysis**

All statistical analyses were performed using R software (version
4.5.2). Unless otherwise specified, all statistical tests were
two-sided, and P-values \< 0.05 were considered statistically
significant. Where appropriate, results were interpreted in conjunction
with multiple-testing considerations and consistency across analytical
layers.

**3. Results**

**3.1 Genetic liability reveals subtype-specific risk architecture**

We first assessed the causal effects of metabolic, inflammatory, and
lifestyle exposures on lung cancer outcomes using Mendelian
randomization. Smoking initiation exhibited the strongest and most
consistent associations across outcomes, with larger effect sizes
observed in lung squamous cell carcinoma (LUSC) compared with lung
adenocarcinoma (LUAD). Body mass index (BMI) showed significant
associations with LUSC but not LUAD, suggesting subtype-specific
susceptibility patterns (Figure 2; Table 2).

Detailed effect estimates, including odds ratios and confidence
intervals, are provided in Table 2.

Multivariable MR analyses attenuated these associations, consistent with
shared genetic architecture among cardiometabolic traits. Collectively,
these findings support a networked metabolic--inflammatory liability
framework rather than a single dominant causal factor (Figure S2; Tables
S1--S2).

***\[Insert Figure 2 and Table 2 here\]***

**3.2 Inflammatory pathways may mediate metabolic risk in LUSC**

To explore potential mechanisms linking metabolic traits to LUSC, we
performed two-step Mendelian randomization analyses. In particular,
C-reactive protein (CRP) emerged as a potential mediator, with the
BMI→CRP→LUSC pathway showing consistent exposure--mediator and
mediator--outcome effects.

Although these associations did not survive multiple testing correction,
the recurrent involvement of CRP across multiple pathways suggests a
biologically plausible inflammatory axis underlying LUSC risk,
consistent with prior evidence linking systemic inflammation to lung
carcinogenesis \[15\] (Figure 2C; Tables S2--S3).

**3.3 Epigenetic triangulation prioritizes MFAP2-associated CpG sites**

To bridge genetic liability to tumor-level regulation, we implemented an
epigenetic triangulation framework integrating mQTL mapping, Mendelian
randomization, and colocalization analysis (see Methods). This approach
enabled the identification of CpG-specific regulatory signals.

At the MFAP2 locus, one CpG site demonstrated concordant evidence across
all layers, whereas adjacent probes did not, indicating probe-specific
functional relevance. Colocalization analysis further supported a shared
genetic signal at this locus.

These findings prioritized MFAP2 as an upstream candidate node linking
inherited liabilities to epigenetic regulation and downstream tumor
biology (Figure 3; Table 3; Figure S3; Table S4).

***\[Insert Figure 3 and Table 3 here\]***

**3.4 Multi-omics integration reveals a coagulation--ECM tumor program**

To further characterize the biological relevance of prioritized loci, we
performed multi-omics integration across transcriptomic, proteomic, and
immune-related datasets. These analyses converged on a coherent tumor
program characterized by platelet activation, coagulation processes, and
extracellular matrix remodeling, with Figure 6 further illustrating
CpG--pathway associations and immune-related signatures (Figures 4 and
6).

Notably, transcriptomic and proteomic analyses showed partial
discordance at the gene expression level but consistent activation at
the pathway level, suggesting post-transcriptional regulation.
Functional enrichment further supported a platelet--ECM axis and
microtubule-related transport processes, in line with prior evidence
implicating platelet activation, extracellular matrix remodeling, and
metastatic support functions in tumor progression \[12,21\].

Together, these results provide a mechanistic bridge linking systemic
metabolic--inflammatory liabilities to tumor microenvironment remodeling
(Table S5).

***\[Insert Figures 4 and 6 here\]***

**3.5 A 12-CpG methylation signature enables clinically relevant risk
stratification**

We translated the prioritized epigenetic signals into a
methylation-based prognostic model using penalized LASSO-Cox regression.
The final model comprised 12 CpG sites (M12) and demonstrated stable
prognostic performance in the TCGA-LUSC training cohort, consistent with
contemporary principles for parsimonious prognostic modeling \[35\]
(Figure 5A--C; Table 4; Tables S6--S7).

***\[Insert Figure 5A--C and Table 4 here\]***

Using the median training-cohort cutoff, Kaplan--Meier analysis showed
significantly worse overall survival in the high-risk group (Figure 5D).
Continuous risk modeling further supported the prognostic relevance of
the M12 signature. Model coefficients and feature contributions are
summarized in Table 4.

***\[Insert Figure 5D here\]***

**3.6 External validation confirms model portability across independent
cohorts**

We next evaluated model performance across the TCGA internal validation
cohort and two independent external cohorts (GSE39279 and GSE30219).
Consistent survival stratification was observed in the internal
validation subset and across external datasets using maxstat-derived
optimal cutoffs, supporting the portability of the M12 model (Figure
7B--D).

In the TCGA training cohort, the model achieved a C-index of 0.664 and a
3-year AUC of 0.685, with good calibration performance (slope = 1.227;
R² = 0.843). Detailed cohort-level performance metrics are summarized in
Table 5 and Table S8, whereas calibration analysis demonstrated good
agreement between predicted and observed survival in the training cohort
(Figure S5A).

***\[Insert Figure 7 and Table 5 here\]***

Taken together, these findings indicate that the M12 signature captures
clinically meaningful prognostic information and supports both
group-based stratification and continuous risk assessment; additional
diagnostic analyses further supported model robustness (Figure S5; Table
S9).

**4. Discussion**

In this study, we developed a causal-to-translational multi-omics
framework that links genetic liability to clinically relevant epigenetic
markers in lung squamous cell carcinoma (LUSC). By integrating Mendelian
randomization with epigenetic triangulation and multi-omics
interpretation, our approach moves beyond conventional association-based
biomarker discovery and aligns with broader efforts to connect causal
inference with molecular profiling in complex disease biology \[17\].
These signals were subsequently translated into a 12-CpG methylation
signature (M12), which demonstrated consistent prognostic stratification
in both training and independent validation cohorts, with acceptable
discrimination and good calibration. Together, these findings support a
unified framework in which upstream metabolic--inflammatory liabilities
are connected to tumor microenvironment remodeling and ultimately to
clinically actionable risk prediction.

A key insight from this work is that metabolic--inflammatory risk in
LUSC is best conceptualized as a networked liability rather than a
single-factor effect \[15,25\]. The recurrent involvement of CRP across
multiple pathways supports an inflammatory axis linking systemic
metabolic perturbations to tumor risk. Although these findings should be
interpreted as hypothesis-generating rather than definitive mediation,
they provide a coherent biological framework for understanding
subtype-specific susceptibility.

The epigenetic triangulation strategy represents an important refinement
step by integrating genetic and molecular evidence at the CpG level. The
probe-specific signal observed at the MFAP2 locus highlights the
importance of fine-grained epigenetic resolution, as adjacent CpG sites
within the same locus may exhibit distinct functional roles. This
finding underscores the value of combining mQTL, MR, and colocalization
to improve the specificity of biomarker discovery \[11,30\].

Multi-omics integration further demonstrated that these genetically
anchored signals converge on a coherent tumor program characterized by
coagulation, platelet activation, extracellular matrix remodeling, and
immune-contexture shifts within the tumor microenvironment \[2\]. This
program provides a plausible mechanistic bridge linking systemic
liabilities to tumor microenvironment dynamics, consistent with
foundational links between inflammatory microenvironments and tumor
promotion \[8\]. It is also biologically compatible with prior evidence
implicating platelet-mediated interactions in metastatic progression
\[12,21\].

From a clinical perspective, the M12 methylation signature illustrates
that biologically informed multi-feature models can capture meaningful
prognostic information. Although its discriminative performance was
moderate, calibration was good in the training cohort and significant
survival stratification was preserved across independent datasets. The
stronger performance of continuous risk assessment further suggests that
the signature may be most informative when interpreted as a quantitative
prognostic measure, consistent with established recommendations for
clinically usable prediction models and epigenetic biomarker translation
\[35,42\].

Importantly, beyond the specific findings, this study provides a
generalizable framework for biomarker discovery that integrates causal
inference with multi-omics validation. This design enables a systematic
transition from genetic association to functional prioritization and
clinical translation, in line with broader efforts to integrate
molecular profiling across complex diseases \[17\], and may be
applicable to other complex diseases.

Several limitations should be acknowledged. First, MR is designed for
etiological inference rather than prognostic modeling, and the
transition from genetic liability to recurrence outcomes should be
interpreted with caution \[6\]. Second, colocalization analyses may be
affected by locus complexity and the presence of multiple causal
variants \[43\]. Third, external validation was conducted in relatively
modest cohorts, which may limit generalizability. Future studies should
focus on prospective validation and experimental characterization of
prioritized CpG sites.

In conclusion, our findings support a causal-to-translational paradigm
for biomarker discovery, demonstrating how genetically informed
epigenetic signals can be leveraged to bridge upstream risk factors and
clinical outcomes in LUSC, thereby providing a generalizable blueprint
for translating causal inference into clinically actionable biomarkers
in complex diseases.

**5. Conclusion**

This study establishes a causal-to-translational paradigm linking
genetic liability to tumor biology and clinical outcomes in lung
squamous cell carcinoma. Through integrative analyses, we identify
MFAP2-linked epigenetic regulation as a key node connecting systemic
metabolic--inflammatory risk to tumor progression. The derived M12
CpG-based signature demonstrates that biologically informed models can
capture clinically relevant prognostic signals, particularly when
interpreted as continuous measures. Beyond the specific findings, this
work provides a generalizable framework for translating causal inference
into mechanistically grounded biomarkers for precision oncology.

**Acknowledgements**

We acknowledge the investigators and participants of the International
Lung Cancer Consortium (ILCCO), The Cancer Genome Atlas (TCGA), and the
Clinical Proteomic Tumor Analysis Consortium (CPTAC) for generating and
making these data publicly available. We also thank the developers of
the analytical tools and open-source resources that made this study
possible.

**Author contributions**

Xiongjie Li (XL): Conceived the study, designed the methodological
framework, performed data acquisition and curation, conducted Mendelian
randomization analyses and prognostic modeling, and drafted the original
manuscript.

Fengyue Zhang (FZ): Participated in data acquisition and curation,
Mendelian randomization analyses, and prognostic modeling.

Xuan Xu (XX): Participated in multi-omics data integration and data
visualization.

Zhenyao Wu (ZW): Participated in multi-omics data integration and data
visualization.

Xiaoyan Zhang (XZ): Provided clinical interpretation and assisted with
the external validation analyses.

Xianghui Wang (XW): Supervised the study, provided administrative
support, and critically revised the manuscript.

All authors read and approved the final manuscript.

Xiongjie Li and Fengyue Zhang contributed equally to this work.

**Funding**

This research did not receive any specific grant from funding agencies
in the public, commercial, or not-for-profit sectors.

**Declarations**

**Competing interests**

The authors declare no competing interests.

**Ethics approval and consent to participate**

Ethical review and approval were waived for this study because all
analyses were based on publicly available, de-identified summary
statistics and open-access molecular datasets (e.g., GWAS, TCGA, CPTAC).
All original studies obtained appropriate institutional ethical
approvals and informed consent from participants, in accordance with the
Declaration of Helsinki.

**Availability of data and materials**

Publicly available datasets were analyzed in this study. GWAS summary
statistics were obtained from the IEU OpenGWAS platform
(https://gwas.mrcieu.ac.uk/). Transcriptomic and methylation data were
accessed through The Cancer Genome Atlas data portal
(https://portal.gdc.cancer.gov/). Proteomic data were obtained from the
Clinical Proteomic Tumor Analysis Consortium data portal. Code for
Mendelian randomization analyses, prognostic modeling, and visualization
is available at: https://github.com/morningLxj/Metabolic-LungCancer-MR.
All data sources used in this study are publicly accessible, and
detailed accession information is provided in the Methods section.

**References**

1\. 23andMe Research Team, HUNT All-In Psychiatry, Liu M, et al.
Association studies of up to 1.2 million individuals yield new insights
into the genetic etiology of tobacco and alcohol use. Nat Genet.
2019;51(2):237-244. doi:10.1038/s41588-018-0307-5

2\. Binnewies M, Roberts EW, Kersten K, et al. Understanding the tumor
immune microenvironment (TIME) for effective therapy. Nat Med.
2018;24(5):541-550. doi:10.1038/s41591-018-0014-x

3\. Bowden J, Davey Smith G, Burgess S. Mendelian randomization with
invalid instruments: Effect estimation and bias detection through egger
regression. International Journal of Epidemiology. 2015;44(2):512-525.
doi:10.1093/ije/dyv080

4\. Burgess S, Davey Smith G, Davies NM, et al. Guidelines for
performing mendelian randomization investigations: Update for summer
2023. Wellcome Open Res. 2023;4:186.
doi:10.12688/wellcomeopenres.15555.3

5\. Campbell JD, Alexandrov A, Kim J, et al. Distinct patterns of
somatic genome alterations in lung adenocarcinomas and squamous cell
carcinomas. Nat Genet. 2016;48(6):607-616. doi:10.1038/ng.3564

6\. Carter AR, Sanderson E, Hammerton G, et al. Mendelian randomisation
for mediation analysis: Current methods and challenges for
implementation. Eur J Epidemiol. 2021;36(5):465-478.
doi:10.1007/s10654-021-00757-1

7\. Collins GS, Reitsma JB, Altman DG, Moons KGM. Transparent reporting
of a multivariable prediction model for individual prognosis or
diagnosis (TRIPOD): The TRIPOD statement. BMJ. 2015;350:g7594.
doi:10.1136/bmj.g7594

8\. Coussens LM, Werb Z. Inflammation and cancer. Nature.
2002;420(6917):860-867. doi:10.1038/nature01322

9\. Dor Y, Cedar H. Principles of DNA methylation and their implications
for biology and medicine. Lancet. 2018;392(10149):777-786.
doi:10.1016/S0140-6736(18)31268-6

10\. Fabregat A, Sidiropoulos K, Garapati P, et al. The reactome pathway
knowledgebase. Nucleic Acids Res. 2016;44(D1):D481-D487.
doi:10.1093/nar/gkv1351

11\. Gaunt TR, Shihab HA, Hemani G, et al. Systematic identification of
genetic influences on methylation across the human life course. Genome
Biol. 2016;17(1):61. doi:10.1186/s13059-016-0926-z

12\. Gay LJ, Felding-Habermann B. Contribution of platelets to tumour
metastasis. Nat Rev Cancer. 2011;11(2):123-134. doi:10.1038/nrc3004

13\. Giambartolomei C, Vukcevic D, Schadt EE, et al. Bayesian test for
colocalisation between pairs of genetic association studies using
summary statistics. Williams SM, ed. PLoS Genet. 2014;10(5):e1004383.
doi:10.1371/journal.pgen.1004383

14\. Goldstraw P, Chansky K, Crowley J, et al. The IASLC lung cancer
staging project: Proposals for revision of the TNM stage groupings in
the forthcoming (eighth) edition of the TNM classification for lung
cancer. Journal of Thoracic Oncology. 2016;11(1):39-51.
doi:10.1016/j.jtho.2015.09.009

15\. Greten FR, Grivennikov SI. Inflammation and cancer: Triggers,
mechanisms, and consequences. Immunity. 2019;51(1):27-41.
doi:10.1016/j.immuni.2019.06.025

16\. Hänzelmann S, Castelo R, Guinney J. GSVA: Gene set variation
analysis for microarray and RNA-seq data. BMC Bioinformatics.
2013;14(1):7. doi:10.1186/1471-2105-14-7

17\. Hasin Y, Seldin M, Lusis A. Multi-omics approaches to disease.
Genome Biol. 2017;18(1):83. doi:10.1186/s13059-017-1215-1

18\. Hemani G, Zheng J, Elsworth B, et al. The MR-base platform supports
systematic causal inference across the human phenome. eLife.
2018;7:e34408. doi:10.7554/eLife.34408

19\. Johnson WE, Li C, Rabinovic A. Adjusting batch effects in
microarray expression data using empirical bayes methods. Biostatistics.
2007;8(1):118-127. doi:10.1093/biostatistics/kxj037

20\. Kaelin WG, McKnight SL. Influence of metabolism on epigenetics and
disease. Cell. 2013;153(1):56-69. doi:10.1016/j.cell.2013.03.004

21\. Labelle M, Begum S, Hynes RO. Direct signaling between platelets
and cancer cells induces an epithelial-mesenchymal-like transition and
promotes metastasis. Cancer Cell. 2011;20(5):576-590.
doi:10.1016/j.ccr.2011.09.009

22\. Laird PW. The power and the promise of DNA methylation markers. Nat
Rev Cancer. 2003;3(4):253-266. doi:10.1038/nrc1045

23\. Lau SCM, Pan Y, Velcheti V, Wong KK. Squamous cell lung cancer:
Current landscape and future therapeutic options. Cancer Cell.
2022;40(11):1279-1293. doi:10.1016/j.ccell.2022.09.018

24\. Love MI, Huber W, Anders S. Moderated estimation of fold change and
dispersion for RNA-seq data with DESeq2. Genome Biol. 2014;15(12):550.
doi:10.1186/s13059-014-0550-8

25\. Mantovani A, Allavena P, Sica A, Balkwill F. Cancer-related
inflammation. Nature. 2008;454(7203):436-444. doi:10.1038/nature07205

26\. Min JL, Hemani G, Hannon E, et al. Genomic and phenotypic insights
from an atlas of genetic effects on DNA methylation. Nat Genet.
2021;53(9):1311-1321. doi:10.1038/s41588-021-00923-x

27\. Potter AL, Costantino CL, Suliman RA, et al. Recurrence after
complete resection for non-small cell lung cancer in the national lung
screening trial. The Annals of Thoracic Surgery. 2023;116(4):684-692.
doi:10.1016/j.athoracsur.2023.06.004

28\. Rajaram R, Huang Q, Li RZ, et al. Recurrence-free survival in
patients with surgically resected non-small cell lung cancer. CHEST.
2024;165(5):1260-1270. doi:10.1016/j.chest.2023.11.042

29\. Relton CL, Davey Smith G. Two-step epigenetic mendelian
randomization: A strategy for establishing the causal role of epigenetic
processes in pathways to disease. Int J Epidemiol. 2012;41(1):161-176.
doi:10.1093/ije/dyr233

30\. Richardson TG, Haycock PC, Zheng J, et al. Systematic mendelian
randomization framework elucidates hundreds of CpG sites which may
mediate the influence of genetic variants on disease. Human Molecular
Genetics. 2018;27(18):3293-3304. doi:10.1093/hmg/ddy210

31\. Robinson MD, McCarthy DJ, Smyth GK. edgeR : A bioconductor package
for differential expression analysis of digital gene expression data.
Bioinformatics. 2010;26(1):139-140. doi:10.1093/bioinformatics/btp616

32\. Sanderson E, Davey Smith G, Windmeijer F, Bowden J. An examination
of multivariable mendelian randomization in the single-sample and
two-sample summary data settings. International Journal of Epidemiology.
2019;48(3):713-727. doi:10.1093/ije/dyy262

33\. Siegel RL, Giaquinto AN, Jemal A. Cancer statistics, 2024. CA
Cancer J Clin. 2024;74(1):12-49. doi:10.3322/caac.21820

34\. Simon N, Friedman J, Hastie T, Tibshirani R. Regularization paths
for cox's proportional hazards model via coordinate descent. J Stat
Softw. 2011;39(5):1-13. doi:10.18637/jss.v039.i05

35\. Steyerberg EW. Clinical Prediction Models: A Practical Approach to
Development, Validation, and Updating. Springer; 2019.

36\. Sung H, Ferlay J, Siegel RL, et al. Global cancer statistics 2020:
GLOBOCAN estimates of incidence and mortality worldwide for 36 cancers
in 185 countries. CA A Cancer J Clinicians. 2021;71(3):209-249.
doi:10.3322/caac.21660

37\. The Cancer Genome Atlas Research Network. Comprehensive genomic
characterization of squamous cell lung cancers. Nature.
2012;489(7417):519-525. doi:10.1038/nature11404

38\. The GTEx Consortium, Aguet F, Anand S, et al. The GTEx consortium
atlas of genetic regulatory effects across human tissues. Science.
2020;369(6509):1318-1330. doi:10.1126/science.aaz1776

39\. The LifeLines Cohort Study, The ADIPOGen Consortium, The AGEN-BMI
Working Group, et al. Genetic studies of body mass index yield new
insights for obesity biology. Nature. 2015;518(7538):197-206.
doi:10.1038/nature14177

40\. Uno H, Cai T, Pencina MJ, D'Agostino RB, Wei LJ. On the
C-statistics for evaluating overall adequacy of risk prediction
procedures with censored survival data. Stat Med. 2011;30(10):1105-1117.
doi:10.1002/sim.4154

41\. Verbanck M, Chen CY, Neale B, Do R. Publisher correction: Detection
of widespread horizontal pleiotropy in causal relationships inferred
from mendelian randomization between complex traits and diseases. Nat
Genet. 2018;50(8):1196-1196. doi:10.1038/s41588-018-0164-2

42\. Wagner W. How to translate DNA methylation biomarkers into clinical
practice. Front Cell Dev Biol. 2022;10:854797.
doi:10.3389/fcell.2022.854797

43\. Wallace C. A more accurate method for colocalisation analysis
allowing for multiple causal variants. Cordell HJ, ed. PLoS Genet.
2021;17(9):e1009440. doi:10.1371/journal.pgen.1009440

44\. Wolff RF, Moons KGM, Riley RD, et al. PROBAST: A tool to assess the
risk of bias and applicability of prediction model studies. Ann Intern
Med. 2019;170(1):51-58. doi:10.7326/M18-1376

45\. Yu G, Wang LG, Han Y, He QY. clusterProfiler: An R package for
comparing biological themes among gene clusters. OMICS: A Journal of
Integrative Biology. 2012;16(5):284-287. doi:10.1089/omi.2011.0118
