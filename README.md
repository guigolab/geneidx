# Geneid+BLASTx in Nextflow

Using the [Nexflow tutorial for the Elixir workflow workshop](https://nextflow-io.github.io/elixir-workshop-21/docs/) taking place on November 2021, as an initial template, we are integrating Geneid+BLASTx into Nextflow so that it can be run easily with almost no requirements.


## Schema
1. **Make sure ( Docker or Singularity ) and Nextflow are installed in your computer.**
  - [Docker](https://docs.docker.com/engine/install/)
  - [Singularity](https://sylabs.io/guides/3.0/user-guide/installation.html#)
  - [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation)

2. **Receive the FASTA file, the parameter file and the set of proteins to use as a reference as input in the params.config file.**
3. **Uncompress the FASTA file (keeping the compressed file).**

4. MISSING: RepeatMasking!! I need to work on this as it is important for producing a good annotation.

5. **Create the protein database for DIAMOND to use it as a source.**
6. **Align the provided genome against the database created using DIAMOND BLASTx flavour.**
7. **Run Geneid using the matches provided by the previous alignment as an additional evidence.**
This is done in parallel for each independent sequence inside the genome FASTA file. It consists of several steps:
  - Pre-process the matches to join those that overlap and re-score them appropriately.
  - Run Geneid using the evidence and obtain the GFF3 file.
  - Remove the files from the internal steps.
8. **Concatenate all the outputs into a single GFF3 file that is sorted by coordinates and without the lines starting with #.**

DETAILS:
- **The name of the sequences in the FASTA file cannot contain unusual characters.**
- **The input genome file must be a gzip-compressed FASTA file.**


DOING:
-  Optimizing the memory and CPU requirements for it to run smoothly on the cluster.


PENDING:
- Define accurate memory limits (dynamic based on genome sizes?).
- Tune DIAMOND parameters to make the most of the resources available and to adapt to the capacity of each computer.
- **Choose the parameter file appropriately:**
	- **Have a dictionary of taxid2parameter file.**
	- **Using the taxid of the species to annotate, identify the closest parameter file.**
- **Auto-train the parameter file.**
