# Subworkflows and modules

Here we store the subworkflows and modules used by the pipeline described in main.nf.

Workflows for the processing of proteins:
- prot_down_workflow inside getProteins.nf
- build_protein_db subworkflow inside build_dmnd_db.nf
- alignGenome_Proteins subworkflow inside runDMND_BLASTx.nf

Workflow for autotraining procedures:
- matchAssessment workflow inside getTrainingSeq.nf

  Also several subworkflows inside the following files.
  - introns_estimates.nf
  - modifyParamFile.nf
  - CDS_estimates.nf

Workflow for gene prediction:
- geneid_WORKFLOW inside geneid.nf

Workflows for concatenation:
- prep_concat subworkflow inside prepare_concatenation.nf
- concatenate_Outputs_once subworkflow inside geneid_concatenate.nf

Several tools for different functionalities inside the tools.nf file:
- unzipFASTA
- compress_n_indexFASTA
- gff34portal
