Before running this, we need to identify marker genes for each cluster (cell type) in the ovary based on single-cell ovary RNAseq studies. This was done by Melody Parker following scripts here - https://github.com/melparker101/p50_Infertility/tree/main/cluster_marker_genes.

Scripts in this directory assess partitioned heritability for ovary cell-type-specific enrichment. Full instructions to do this are in the tutorial from Bulik-Sullivan here: https://github.com/bulik/ldsc/wiki/Cell-type-specific-analyses.

First run **0_create_ldscores_baseline.sh** to create custom baseline LD scores.

1. **1_create_custom_annot_files.sh** - Create custom annot files using marker gene-sets identified from single-cell ovary RNAseq differential gene expression studies. Batch submit these files using **batch_submit_create_custom_annot_files.R**. The LDSC-provided **make_annot.py** script has to be tweaked to create annotations, which is documented in **custom_make_annot.py**.
2. **2_partitioned_heritability_custom.sh** - Submit partitioned heritability analyses for custom annotations created above.
3. **3_plot_partitioned_heritability_gtex.R** - Scripts to plot the results as presented in the manuscript.