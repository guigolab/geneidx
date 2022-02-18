# Subworkflows and modules

Here we store the subworkflows and modules used by the pipeline described in main.nf:

- Module list_files_to_download
- Module DownloadFASTA_fromID
- Sub-workflow geneid_WORKFLOW_single including:
   - Module UncompressFASTA
   - Module runGeneid_single
