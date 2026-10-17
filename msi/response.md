# Response to Reviewers

**Manuscript:** OSSpRe: An Open Source Spatial Reactomics Workflow for Trapped Ion Mobility Spectrometry Coupled with High Resolution Mass Spectrometry Imaging

**Authors:** Yu et al.

We thank both reviewers for careful and constructive comments. The reviews prompted us to (i) add a concrete, biologically interpreted PMD example with quantitative spatial evidence, (ii) add empirical performance metrics (CCS accuracy, spatial statistics), (iii) correct several technical statements about TIMS/CCS and isomer resolution, (iv) document the HPC deployment concretely, and (v) revise overclaims ("standalone", "ensures annotation", Table S1 comparisons). Below, reviewer comments are in *italics*; our responses follow. New or revised manuscript text is shown in quotation blocks. All analyses described below use only the public demo data (Zenodo) and can be reproduced with a script we have added to the GitHub repository.

---

## Reviewer 1

**General comment:** *"...validations and convincing demonstration of the most important scientific steps (feature detection, molecular annotation, and reaction inference) are very limited or absent... The scientific claims currently go considerably beyond what the workflow indeed delivers."*

**Response:** We agree the original draft overclaimed in places and underreported quantitative evidence. We have (1) added empirical validation using the public demo dataset wherever ground truth or reference values exist (CCS accuracy, annotation discrimination, spatial reproducibility), (2) added a fully worked biological example of spatial reactomics (new subsection + Figure X), (3) reworded all overclaims, and (4) reframed the remaining unbenchmarked aspects explicitly as limitations and future work. Details per point below.

### R1-1. Performance and computational claims

*"There is no validation of the performance and computational claims... CCS accuracy? ... What is the accuracy of spatial segmentation? ... runtime/memory scaling..."*

**Response:** We have added the following quantitative evidence to the revised manuscript and Supplementary Information:

- **CCS annotation residuals (renamed from "CCS accuracy").** Across the 227 unique annotated reference features in the demo dataset, the |CCS deviation| from database values has a median of 0.82% (mean 0.86%; 90th percentile 1.70%; maximum 2.25% ≈ 5 Å²). We now report these as *annotation residuals* rather than "accuracy": by construction, matched features lie within the ±5 Å² annotation window (see response to Comment 7), so the distribution is window-truncated and cannot be read as an independent accuracy estimate. We additionally report: (i) a systematic offset — 82% of signed CCS deviations are negative (median −0.84%) and m/z deviations show a common negative bias (e.g., −8.9 ppm for both ions of the PC 36:4/PC 36:3 pair), consistent with a calibration offset rather than random error; and (ii) a decoy (database-shift) check: of 7,442 features, 683 match the target lipid CCS database within ±20 ppm/±5 Å², whereas shifting the database CCS by +15 Å² or −15 Å² yields 103 and 438 matches respectively, and m/z-only matching yields 1,072 — i.e., the CCS filter removes ~36% of m/z-only candidate matches, with a non-negligible residual false-match rate at this window width. [CHECK: report the same residual statistics restricted to the high-confidence subset of ions whose annotations agree with the vendor/software annotations shown in Fig. S6.]
- **Super-pixel vs pixel-wise extraction.** We now state the design intent precisely: pixel-wise peak picking is not a gold standard but a slower, noisier alternative on low-abundance pixels — the super pixel sums all frames so that low-abundance analytes accumulate SNR before peak detection. We report the existing timing comparison (436 s pixel-wise vs 30 s super-pixel for 2,000 pixels) plus a new feature-recall check: [RUN: extract peaks pixel-wise on a 2,000-pixel subset and report the fraction of super-pixel reference peaks recovered, and the number of pixel-wise peaks not present in the reference list (expected: mostly noise/low-abundance)].
- **Runtime/memory scaling.** [RUN: add a Supplementary table of runtime and peak RAM vs number of pixels (e.g., 10k / 50k / 128k pixels of the demo dataset) for the reference-peak and quantification steps, single core and multi-core.] We agree a single example is insufficient; the revised Supplementary includes this scaling table.
- **Segmentation accuracy and isomer resolution.** We respectfully note that no ground truth exists for "segmentation accuracy" of unannotated MSI data beyond histological overlay, and validated isomer resolution requires standards and MS/MS, which is beyond the scope of a software Application Note. We have softened these claims and now present segmentation as exploratory (with the UMAP/PCA descriptions moved to the Supplementary), and we explicitly state that isomer-level assignments require orthogonal confirmation (see also Comment 5).

### R1-2. Biological demonstration

*"Including at least one concrete example with biological interpretation would substantially improve the manuscript..."* and *"be careful about statistical language referring to 'conditions', 'samples' and ANOVA."*

**Response:** We have added a new subsection "Example: spatial reactomics of a desaturation pair" with a new figure (Figure X) built entirely from the public demo dataset. In brief:

- **The pair.** PC 36:4 (m/z 804.5447, CCS 285.19, [M+Na]⁺) and PC 36:3 (m/z 806.5603, CCS 285.18, [M+Na]⁺), connected by PMD 2.0157 (H₂; C₄₄H₈₀NO₈P ↔ C₄₄H₈₂NO₈P), consistent with a (de)saturation/hydrogenation-type transformation. We have unified the terminology for PMD 2.0157 throughout the manuscript (previously described as "oxidation" in the Introduction and "double bond breaking" in the Results; now consistently "H₂ (de)saturation").
- **Annotation confidence.** At m/z 806.56 the competing database annotation is C12 lactosyl ceramide (d18:1/12:0) [M+H]⁺ (LacCer/Hex₂Cer; CCS deviation −1.36%, m/z deviation −2.6 ppm). Its CCS deviation lies within the 1–2% typical of TIMS measurements, so CCS alone does **not** disfavor it. Our confidence in the PC 36:3 assignment instead rests on converging, orthogonal evidence: (1) the PC 36:4 and PC 36:3 [M+Na]⁺ ions show *identical* m/z deviations (both −8.9 ppm), consistent with a shared systematic calibration offset and same-adduct chemistry; (2) their CCS values are mutually consistent (285.19 vs 285.18 Å²); and (3) they are connected by the PMD 2.0157 relationship expected for a (de)saturation pair within the same lipid class and adduct. We present this explicitly as multiple weak pieces of evidence combining into a candidate assignment, not as CCS-based exclusion. Notably, at m/z 804.54 the same molecule appears as two database records ("18:2 PC" and "PC 36:04") whose reference CCS values differ by 1.7% (290.0 vs 285.1 Å²) — direct evidence that inter-database CCS variability is of the same order as the deviations discussed here (see also Comment 7).
- **Spatial evidence.** Across 127,577 pixels where both ions co-occur, the per-pixel Pearson correlation is r = 0.81 (independently recomputed; consistent with the network export, 0.808). For context, this ranks in the top 0.013% of all 15,133 exported network edges, whereas PMD 2.0157 edges have a median correlation of only 0.04. We note that strong co-localization of two same-class, same-adduct lipids is expected from shared tissue density and ionization conditions alone, so we phrase the inference as *consistent with* local interconversion/candidate transformation — not as evidence that the reaction is occurring.
- **ROI contrast.** Using the H&E-validated islet ROI from the paper's own UMAP segmentation (segment 5; 4,289 pixels), the PC 36:4/PC 36:3 intensity-ratio median is 0.88 inside islets vs 1.48 outside — a stronger and threshold-free contrast. We have replaced the previous SM 34:1-abundance threshold definition (which we found sensitive to the arbitrary quantile cut: high-abundance-region ratios 1.43/1.34/1.20/1.16 at the 50/75/90/95% quantiles) with this segmentation-based ROI, and now cite the established islet-marker literature where marker-based definitions are discussed. [ADD ref for SM 34:1 as islet marker, or drop the marker-based definition.]
- **Inference level.** Consistent with the PMD framework, this is presented as a *candidate* biochemical transformation supported by accurate mass, CCS consistency, same-adduct pairing, and spatial covariance; acyl composition and reaction direction require MS/MS.
- **Statistical language.** We have removed "conditions/samples/ANOVA" phrasing associated with the single demo sample; group-level tests are now described as capabilities of the workflow for multi-sample designs, and the demo analysis uses descriptive ROI-level contrasts only.

We also provide an interactive, self-contained Google Colab notebook (`msi/OSSpRe_Reviewer_Response_Reproducibility.ipynb`) in the repository so reviewers can execute and reproduce the biological example, statistical tables, and distribution figures directly in the cloud in under 3 minutes without manual data staging.

### R1-3. PMD and reactomics background

*"PMD and reactomics are central concepts in this work, but they are currently introduced too briefly..."*

**Response:** Agreed. The revised Introduction now: (1) defines a PMD operationally (mass difference between co-detected features that matches a chemical transformation at instrument mass accuracy); (2) summarizes prior PMD/reactomics methodology and applications in LC-MS based untargeted analysis with citations (Yu et al. 2019; Yu & Petrick 2020; and further methodological/application studies to be added); (3) explains how specific mass differences map to putative transformations (H₂ = desaturation, O = hydroxylation/oxidation, C₂H₄ = chain elongation) and how reaction databases (KEGG reaction; HMDB) are used; and (4) explicitly distinguishes candidate/putative transformations from experimentally validated reactions, including a statement of the inference level supported by PMD relationships alone. [ADD: 5–8 new references on PMD/reactomics applications.]

### R1-4. TIMS measures inverse mobility, not CCS directly

*"This sentence is incorrect... Is OSSpRe truly calculating CCS values or reading it/using functions from a package?"*

**Response:** The reviewer is correct, and we thank them for this correction. We have also re-checked our description against the code and corrected it accordingly. The actual implementation: 1/K₀ values are read from the raw data through the OpenTIMS library (Łącki et al. 2021); CCS values are then obtained by calling the Bruker TDF SDK function `tims_oneoverk0_to_ccs_for_mz` (which applies the Mason–Schamp-based conversion using the acquisition's calibration), loaded at runtime from the proprietary `libtimsdata` library under EULA acceptance. OSSpRe does not derive CCS from first principles: the C++ file in the repository (`one_over_k0_to_ccs.cpp`) is a thin Rcpp/OpenMP wrapper around that Bruker routine, and the conversion assumes z = 1 for all ions (reasonable for MALDI, and now stated explicitly). Revised text:

> "Raw data are accessed through the OpenTIMS library (Łącki et al. 2021), which returns, per pixel (frame), the intensity, m/z, and inverse ion mobility (1/K₀) of each detected ion. 1/K₀ values are converted to collisional cross sections (CCS) by the Bruker TDF SDK conversion routine (based on the Mason–Schamp equation and the instrument calibration exported with the acquisition), which OSSpRe calls via a wrapper supplied with the repository; all ions are treated as singly charged. OSSpRe does not implement the physical conversion itself."

Two consequences are now stated in the Implementation/dependency sections (see also responses to Reviewer 2, Comments 1 and 8): (a) the CCS conversion depends on the proprietary Bruker library (accepted under EULA as an explicit workflow parameter); (b) `libtimsdata` is provided for Windows and Linux, so raw-data conversion (Stage 1) is restricted to those platforms, while all downstream analysis stages run anywhere R runs.

The corresponding Introduction sentence now reads "...m/z values are augmented with ion mobility separation, from which collisional cross section (CCS) values are derived via calibration."

### R1-5. Isomer differentiation overstatement

*"Using IM does not ensure proper annotation, it can improve annotation confidence..."*

**Response:** Correct; revised. The sentence now reads:

> "Ion mobility provides an orthogonal molecular descriptor (CCS) that reduces candidate annotations and, where mobility separation and reference CCS differences are sufficient, can discriminate isomeric and isobaric species, thereby increasing — but not guaranteeing — annotation confidence."

The new PC 36:4/PC 36:3 example (Comment 2) — where CCS is one of several converging pieces of evidence rather than a discriminator by itself — demonstrates this framing concretely.

### R1-6. How are .d files accessed?

**Response:** Via OpenTIMS/TimsR (Łącki et al. 2021, currently ref. 16), including Bruker library loading under EULA acceptance (an explicit workflow parameter); the 1/K₀→CCS conversion itself is performed by the Bruker TDF SDK routine as described in our response to Comment 4. We have made both citations prominent at first mention of raw-data access.

### R1-7. CCS tolerance

*"Clarify CCS tolerance, an example mentions 5, 5%? Also, a 5% CCS window is broad... typical CCS error reported for TIMS is 2%."*

**Response:** We thank the reviewer for prompting us to re-check this, and we must correct the units in our own description. The annotation step does **not** use a 5% window: the code (`enviGCMS::getalign(..., ppm = 20, deltart = 5)`) uses an **absolute ±5 Å² CCS window** (plus ±20 ppm m/z). Over the CCS range of the annotated features (192–324 Å²), ±5 Å² corresponds to a *relative* 2.6%–1.5% (≈1.8% at 280 Å²) — i.e., the effective annotation window is already close to the 2% the reviewer recommends, and the largest observed |ΔCCS| among matches (≈5 Å²) simply reflects this window. We have removed the "permissive 5% screening window" wording (the unit was omitted in the original draft), and corrected the Supplementary example that read "5% of the current CCS" to "±5 Å²". The three distinct roles are now stated precisely:

> "Annotation against databases uses an absolute CCS window of ±5 Å² together with ±20 ppm m/z (equivalent to ≈1.5–2.6% relative over the 190–330 Å² range covered by the databases); peak picking and peak alignment use a relative CCS tolerance of 3% (`ccs_tolerance = 0.03`); and the observed annotation residuals on this dataset are tighter still (median 0.82%; see response to Comment 1). For TIMS-to-TIMS comparisons of the same species we recommend a 2% relative window; when screening against databases compiled across platforms, an absolute window avoids the scale dependence of a percentage window at low CCS."

We also note (cf. our response to Comment 2) that two database records of the *same* molecule (PC 18:2/18:2, listed as "18:2 PC" and "PC 36:04") differ by 1.7% in reference CCS (290.0 vs 285.1 Å²) — inter-database CCS variability is itself of this order, which is why a window much tighter than ±5 Å² would discard correct matches.

---

## Reviewer 2

**General comment:** *"...more like a collection of R-based scripts... HPC claims... not adequately supported..."*

**Response:** We have revised the manuscript to (1) accurately describe OSSpRe as an R-based, open-source computational workflow (terminology corrected throughout, including the title/abstract if needed); (2) document the parallelization strategy and ship a concrete SLURM deployment script (see Comment 3); (3) add the quantitative evidence described to Reviewer 1 (CCS accuracy, scaling table, worked biological example); and (4) restructure the software comparison (Comment 5). We respectfully note that a full software-engineering program (GUI integration, GPU utilization reporting, scheduler-agnostic deployment frameworks) is beyond the scope of an Application Note; we believe the revised manuscript meets the standard of prior Application Notes for R-based MSI software (e.g., Cardinal).

### R2-1. "Standalone workflow" and usability

**Response:** We have replaced "standalone workflow" with "R-based open-source computational workflow" throughout (see also Comment 8). Our intended users are computational MS/metabolomics researchers working in R — the same audience as Cardinal and MALDIquant — and the revised manuscript states this explicitly, together with the usage procedure (renv-restorable environment, stepwise scripts, video tutorials). Dependencies are now listed explicitly in the Implementation section: R (≥4.x), OpenTIMS for raw-data access, the Bruker `libtimsdata` library for 1/K₀→CCS conversion (proprietary, accepted under EULA as an explicit workflow parameter, and available for Windows/Linux only — Stage-1 raw processing is therefore restricted to those platforms, while all downstream stages are platform-independent), and the m/z–CCS reference databases. For non-programmers, the two Shiny applications cover annotation (m/z–CCS search) and molecular-network exploration. We do not claim, and no longer imply, turn-key GUI operation of the full workflow.

### R2-2. Algorithmic detail

**Response:** The revised Supplementary now contains: the complete 2D local-maximum algorithm (dynamic relative m/z-ppm × CCS window, SNR donut filter) with pseudocode and parameter defaults; the super-pixel aggregation and batching scheme; the TIC/RMS normalization and missing-value filtering rules; and the PMD quantification definition (sum of normalized intensities of paired features, low/high/total per PMD). Complexity is O(N·w) per candidate peak with window size w, and memory is bounded by the user-defined frame batch size. We believe main-text brevity plus this Supplementary level is consistent with Application Note conventions.

### R2-3. HPC claims and implementation

**Response:** The parallelization model is embarrassingly parallel at two levels: (1) *within sample*, the super-pixel construction processes frames in user-defined batches (default 100,000 frames) and the 1/K₀→CCS conversion runs as a parallel apply; (2) *across samples*, each sample/dataset is an independent job with no inter-job communication. The repository now includes (and the Supplementary documents) a concrete SLURM submission script specifying the environment used for the demo dataset: 1 node, 70 cores, 700 GB RAM, 48 h wall time, R restored reproducibly via renv, followed by the two top-level scripts (reference peak picking; quantification). There is no GPU dependency; the workflow is pure R/C++ (CPU only) and runs on any Linux HPC with a SLURM-compatible scheduler (direct execution on workstations is also supported). We have correspondingly revised Table S1 from "Designed for HPC scalability" to "Parallel/HPC batch processing supported", which reflects the actual capability without overclaiming.

### R2-4. Computational performance evaluation

**Response:** We agree the single 436 s vs 30 s example is insufficient. The revised Supplementary includes a benchmark table reporting runtime and peak memory for both major steps across [RUN: three dataset sizes, e.g., ~10k, ~50k, and ~128k pixels, single- and multi-core], on the same hardware, with the exact session information. We have clarified what the 436 s vs 30 s comparison measures (per-pixel peak picking vs reference-peak extraction on 2,000 pixels). CPU/GPU utilization traces and scaling-efficiency curves would support claims we no longer make; we describe the workflow's parallelization factually (Comment 3) instead.

### R2-5. Comparison with other MSI software

**Response:** We have carefully re-verified each tool's current capabilities against vendor documentation and the recent IM-MSI software review of Liu et al. (BioTech 2025, 14:98), and restructured Supplementary Table S1 accordingly:

- **SCiLS Lab.** Our original entry ("no ion mobility support") was incorrect. SCiLS Lab provides drift-time-resolved ion images, 4D peak picking, segmentation, and statistical analysis for timsTOF fleX datasets (Bruker documentation; Liu et al. 2025). This entry has been corrected.
- **Cardinal v3.** Cardinal can import imzML files containing ion mobility information, but treats them as standard MSI data: it "lacks IM-specific functionality and does not support drift time/CCS feature extraction" (Liu et al. 2025). We have therefore retained our entry with a clarifying footnote ("no native IM/CCS-aware processing; IM-containing imzML handled as conventional MSI") rather than an unqualified "No".
- **METASPACE.** We found no evidence of native ion-mobility/CCS-aware processing in METASPACE or METASPACE-ML; the entry is retained with a footnote and we would welcome correction with documentation.
- **Pyxis** (Anal. Bioanal. Chem., Dec 2025), a vendor-neutral open-source tool for IM-MSI analysis, appeared after our submission; we cite it as concurrent work in the revised manuscript, which further confirms the timeliness and demand for open IM-resolved MSI analysis.

Notably, Liu et al. 2025 independently survey the IM-MSI software landscape and already cite OSSpRe as an example of open-source IM-MSI workflows (their ref. 108), while concluding that "current IM-MSI workflows remain fragmented" with "absence of an integrated processing pipeline" — supporting the motivation for OSSpRe. The revised Table S1 differentiates tools along dimensions reflecting genuine design differences (open source; native IM/CCS-aware processing; end-to-end reactomics/PMD analysis; direct vendor .d ingestion), and we no longer claim uniqueness in ion-mobility support per se. A full quantitative head-to-head benchmark (accuracy, runtime, memory, usability) across tools with different input formats, scopes, and output semantics would itself constitute a separate benchmarking study; we now scope the comparison accordingly and note complementarity (e.g., OSSpRe ingests vendor .d data directly and adds reactomics; Cardinal/METASPACE provide mature MS/MS-based annotation). Where capability overlaps, we provide the runtime comparison of Comment R2-4 as an honest but limited data point.

### R2-6. Validation across tissues/conditions

**Response:** The demo study is intentionally a single-tissue demonstration of the workflow, consistent with Application Note conventions; we have revised the wording to avoid implying broad multi-tissue validation ("demonstration on one mouse pancreas dataset"). The workflow itself is tissue-agnostic: it operates on any TIMS-MSI acquisition. We acknowledge in the Discussion that robustness across tissues, sample sizes, and acquisition conditions remains to be established systematically, and we frame multi-dataset evaluation as future work.

### R2-7. Statistical justification

*"multiple-testing correction, assumptions... spatial autocorrelation... inflated significance."*

**Response:** This is a fair and important point. Revisions: (1) the workflow documentation and Supplementary now note that ANOVA/t-tests on PMD intensities are exploratory and, in multi-sample designs, should be accompanied by multiple-testing correction (e.g., Benjamini–Hochberg FDR across PMDs and features); (2) we now report the spatial autocorrelation of PMD maps on the demo dataset explicitly (lag-1 pixel Pearson r = 0.82; lag-20 r = 0.41) and state that treating individual pixels as independent inflates effective significance; we therefore recommend and demonstrate ROI/segment-level (aggregated) inference, which reduces the dependence problem; (3) assumptions of the tests are stated. Revised text:

> "Because MSI intensities exhibit strong spatial autocorrelation (e.g., lag-1 r ≈ 0.8 for PMD maps in the demo dataset), pixel-level tests should be regarded as descriptive; we recommend aggregating intensities within histologically or segmentation-defined ROIs before inferential testing, with multiple-testing correction across PMDs."

### R2-8. "Standalone" terminology

**Response:** Agreed; all occurrences changed to "R-based open-source computational workflow" (title, abstract, main text, Table S1 header), with dependencies stated explicitly in the Implementation section: R, OpenTIMS (raw-data access), the Bruker `libtimsdata` proprietary library for CCS conversion (EULA-gated; Windows/Linux only), and the m/z–CCS annotation databases (see also our response to Reviewer 1, Comment 4, and Reviewer 2, Comment 1).

---

## Summary of manuscript changes

1. New subsection + figure: worked spatial-reactomics example (PC 36:4/PC 36:3, PMD 2.0157) with annotation table, converging-evidence annotation rationale (identical −8.9 ppm m/z offset, mutual CCS consistency, PMD pairing; LacCer competing annotation acknowledged, not claimed excluded), spatial correlation placed in network-wide context (top 0.013% of 15,133 edges; PMD 2.0157 edge median r = 0.04), segmentation-based islet ROI ratio contrast (0.88 inside vs 1.48 outside, H&E-validated segment 5); reproducible script added to the repository.
2. New empirical metrics: CCS *annotation residuals* (median 0.82%) with window-truncation caveat, systematic-offset report (82% negative, median −0.84%), decoy-based false-match assessment (683 target vs 103/+15 Å² and 438/−15 Å² decoy matches; 1,072 m/z-only), spatial autocorrelation statistics, [RUN: runtime/memory scaling table; super-pixel vs pixel-wise recall].
3. Technical corrections: TIMS measures 1/K₀; CCS converted by the Bruker TDF SDK routine (z = 1 assumed) called through a repository wrapper — not by OpenTIMS and not derived from first principles; proprietary-library and Windows/Linux-only constraints stated; isomer claims reworded; CCS tolerance units corrected (annotation = ±5 Å² absolute ≈ 1.5–2.6% relative; peak picking/alignment = 3% relative); PMD 2.0157 terminology unified as H₂ (de)saturation.
4. Terminology: "R-based open-source computational workflow" throughout; Table S1 corrected (Cardinal IM support; SCiLS timsTOF support) and reframed.
5. Statistics: FDR guidance; ROI-level inference recommendation; exploratory framing of demo analyses; removal of "conditions/ANOVA" language for the single demo sample.
6. HPC: SLURM script (msi.job) documented in the Supplementary; parallelization model described factually; overclaims removed.
7. Introduction expanded on PMD/reactomics with additional references. [ADD refs]

## Remaining tasks before submission

- [RUN] Pixel-wise vs super-pixel feature-recall comparison on a 2,000-pixel subset.
- [RUN] Runtime/memory scaling table (≥3 dataset sizes).
- [ADD] PMD/reactomics application references to the Introduction.
- [ADD] Overlay the PMD-ratio panel with the H&E/optical image for the new example figure.
- [FIX] SI text: replace the "5% of the current CCS" example with "±5 Å² absolute (≈1.5–2.6% relative)"; check every other mention of CCS windows for the same unit error.
- [FIX] Manuscript terminology: PMD 2.0157 described as "oxidation" (Introduction) and "double bond breaking" (Results) → unify as "H₂ (de)saturation"; formula arrow for the PC 36:4/PC 36:3 pair should be ↔ (both directions possible).
- [FIX] Rename "CCS accuracy" phrasing to "annotation residuals" everywhere it appears (main text + SI), and add the window-truncation caveat.
- [CHECK] High-confidence residual subset: restrict CCS-residual statistics to ions whose annotations agree with the vendor/software annotations (Fig. S6) and report separately.
- [ADD] Citation for SM 34:1 as an established islet marker (or drop the marker-based ROI definition entirely in favor of the segmentation-based ROI).
- [CHECK] ~~Verify final Table S1 claims against current Cardinal/SCiLS/METASPACE versions.~~ Done (2026-10): SCiLS = supports IM (correct entry); Cardinal = imzML import only, no IM-specific functionality (footnote); METASPACE = no native IM found (footnote). Add Pyxis (ABC 2025) and Liu et al. 2025 (BioTech 14:98, cites OSSpRe) to the manuscript.
- [DECIDE] Whether to also mention the SM 34:1 → SM 36:1 elongation pair (compartment-crossing example) as a second, contrast example.
