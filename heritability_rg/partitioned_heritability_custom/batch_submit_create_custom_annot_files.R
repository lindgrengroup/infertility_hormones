# Author: Samvida S. Venkatesh
# Date: 08/02/2022

mainpath <- "/well/lindgren/samvida/hormones_infertility/ovary_celltype_markers_melody"

list_of_genesets <- read.table(paste0(mainpath, "/all_filenames.txt"),
                               header = F, sep = "\t", stringsAsFactors = F)
colnames(list_of_genesets) <- c("filename", "cluster_set")

submission_script <- paste0(mainpath, "/scripts/1_create_custom_annot_files.sh")

for (cs in list_of_genesets$cluster_set) {
  job_options <- paste0(
    "--export=",
    paste0(
      "CLUSTER_GENESET=\"", cs, "\""
    )
  )
  job_submission <- paste("sbatch", job_options, submission_script)
  system(job_submission)
  print(job_submission)
}
