# DiabeticWoundHealing
Code to reproduce the analyses and figures in the diabetic wound healing study of Wietecha, et al. 2026 <BR>

Code is split up into R scripts for mouse and human data analyses, which are annotated according to figure numbers in the manuscript: <BR>
https://www.biorxiv.org/content/10.64898/2026.04.26.720829v1.full <BR>

List and description of relevant files: <BR>
- mouseDataset_allFigures.R : R script file with code to reproduce all figures related to the mouse wound healing datasets <BR>
- humanDataset_allFigures.R : R script file with code to reproduce all figures related to the human wound healing datasets <BR>
- mmFB_markers.txt : Tab-delimited text file with gene signature markers of wound fibroblast sub-types identified in Almet, et al. J Invest Dermatol (2025) <BR>
- wFB_pheno_Markers.txt : Tab-delimited text file with gene signature markers of early vs late wound fibroblast sub-types identified in Wietecha, et al. Matrix Biol (2023), and genes from relevant REACTOME pathways <BR>
- CommonMoMarkers.txt : Tab-delimited text file with gene signature markers of in vitro M1, M2a and M2c polarized macrophages, and of monocyte-derived macrophage phenotypes identified in Sanin, et al. Sci Immunol (2022) <BR>
- CommonNeutMarkers.txt : Tab-delimited text file with gene signature markers in vitro N1 and N2 polarized neutrophils identified in Mihaila, et al. Front Immunol (2021), and gene signature markers of Neutrotime phases of neutrophil activation identified in Grieshaber-Bouyer, et al. Nat Commun (2021) <BR>

RDS files of Seurat and CellChat objects needed to reproduce the results can be found on Zenodo: <BR>
https://zenodo.org/records/20549537 (DOI: 10.5281/zenodo.20549537) <BR>

Single-cell sequencing raw data of mouse wounds are uploaded to NCBI Sequence Read Archive (SRA): <BR>
https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1467292 (BioProject number PRJNA1467292) <BR>
Single-cell sequencing processed datasets of mouse wounds are uploaded to NCBI Gene Expression Omnibus (GEO): <BR>
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE333591 (Accession number GSE333591; Reviewer token: kdafmyoknzappof)
