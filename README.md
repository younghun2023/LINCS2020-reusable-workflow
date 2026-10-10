# Reusable LINCS2020 screening workflow

Use your own human upregulated/downregulated gene signature. No predefined
160-gene module, JDM dataset, muscle/PBMC comparison, or fixed gene count is required.
The method requires at least one mapped gene in each direction. Non-human
identifiers need an independently documented ortholog-mapping step upstream.

## Inputs

Copy `gene_signature_template.csv` to `inputs/gene_signature.csv` and populate:
- `symbol`: unique human gene symbol.
- `direction`: Up or Down relative to the phenotype you wish to reverse.
- `gene_id`: optional user-verified Entrez ID, particularly for legacy symbols.

Choose your signature and significance/effect-size criteria before screening.
This workflow does not perform differential expression or require a particular
method for choosing genes. No module-size restriction is imposed. Both directions
are required; a one-direction-only query is not supported by this implementation.

## Reference files and run order

Set the R working directory to this package folder. All paths in config.R are
portable. Place already downloaded reference files in data/LINCS2020, or set
download_data <- TRUE. Exact filenames and public URLs are in reference_downloads.csv.
No GEO dataset is required: users bring their own upstream results.
Allow at least 45 GB of disk space for downloads. Existing downloads are reused.
Run 00_install_packages.R if necessary, then source("01_screen.R").
The full analysis can take hours; checkpoints permit resuming matching runs.

## Optional network results

Leave network_file empty for transcriptomic screening alone. To integrate network
evidence, provide a CSV matching network_template.csv and set its path in config.R.
Scores must be calculated independently for your own disease/module. This package
does not construct a PPI network, calculate distances, randomize protein sets,
or adjust P values. Document those methods and testing families separately.
Only unique normalized names are matched; ambiguous entries are preserved.

## Scores and defaults

All matrix-supported query genes are used, including inferred features; mapping
audits report membership. Compound signatures require qc_pass==1 by default;
is_hiq is optional. Codes 2/3 are recognized but excluded without interpreting them.
Review these QC choices against release documentation and your scientific question.
Unweighted enrichment is computed over the full matrix gene universe, then averaged
across retained signatures per perturbagen ID. At least five signatures are required
for ranking. Tau-like = 200 * rank / eligible_count - 100, with average ranks for ties.
These are local rank statistics, not official CMap Tau/NCS or P values.
Selection requires negative mean connectivity and bottom-5% rank by default.
Optional joint selection additionally requires Z <= -0.15; this is exploratory.
Network significance requires negative Z and the supplied adjusted P < 0.05.
Thresholds are configurable; justify changes before inspecting selected compounds.

## Outputs

The run folder contains query mapping, exclusions, signature-level scores,
compound and cell summaries, counts, and software session information.
Its report subfolder contains Compound_results.csv and Transcriptomic_screen.xlsx.
Optional network inputs additionally produce integrated/network tables and a PNG.
Missing scores remain NA; eligible transcriptomic-only results are not discarded.
Labels at tau-like <= -50 and the network threshold affect display only.
Network-only records remain in tables and are not assigned fabricated plot positions.

To regenerate reports only, run 02_report.R. It uses the latest run under results,
or set completed_run to a specific run directory first. It does not rescore the matrix.
These reporting outputs are general-purpose and do not replace manuscript figures
or supplements without review. Existing study-specific files are left unchanged.

## Reproducibility notes

Keep upstream gene-selection code, contrasts, annotation versions and network code.
Mapping, query composition, reference coverage, filtering and aggregation all affect
results. Uneven cell/dose/time representation influences compound means. Inspect
individual signature and cell summaries. Computational selection is not evidence
of therapeutic efficacy. These scripts are not a complete GWAS or expression pipeline.
Source scripts contain no personal paths. Run-produced provenance may contain your
local paths and should be reviewed separately before uploading. Packages are not pinned;
retain session information. No full matrix rerun was performed for this generalization.
