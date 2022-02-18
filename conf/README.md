# Conf folder

Here we have a list of Nextflow configuration files, each one for a different profile.
You can switch among them by using the Nextflow parameter -profile

- standard: run in the local computer [default]
- hpc_sge: submit jobs to a HPC using SGE as scheduler
- hpc_slurm: submit jobs to HPC using SLURM as scheduler
- retry: submit a job and retry if fails
- cloud: submit jobs to Amazon AWS batch system 
