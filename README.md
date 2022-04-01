# Geneid+BLASTx in Nextflow

Using the [Nexflow tutorial for the Elixir workflow workshop](https://nextflow-io.github.io/elixir-workshop-21/docs/) taking place on November 2021, as an initial template, we are integrating Geneid+BLASTx into Nextflow so that it can be run easily with almost no requirements.


## Schema
1. **Make sure ( Docker or Singularity ) and Nextflow are installed in your computer.**
  - [Docker](https://docs.docker.com/engine/install/)
  - [Singularity](https://sylabs.io/guides/3.0/user-guide/installation.html#)
  - [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation)

2. Receive the FASTA file, the parameter file and the set of proteins to use as a reference as input in the params.config file.
3. Uncompress and index the FASTA file.

4. *MISSING: RepeatMasking!!*

5. Create the protein database for DIAMOND to use it as a source.
6. Align the provided genome against the database created using DIAMOND BLASTx flavour.
7. Run the auto-training process:
  - Use matches to estimate the coding sections, look for open reading frames between stop codons.
  - Use matches from the same protein to predict the potential introns.
  - From the sequences of both previous steps compute the initial and transition probability matrices required for the computation of the coding potential of the genome that will be annotated.
8. Build the parameter file with the parameters indicated in the *params.config* file and also the matrices automatically generated from the matches.
9. Run Geneid with the new parameter file and the matches from the previous steps as additional evidence.
This is done in parallel for each independent sequence inside the genome FASTA file. It consists of several steps:
  - Pre-process the matches to join those that overlap and re-score them appropriately.
  - Run Geneid using the evidence and obtain the GFF3 file.
  - Remove the files from the internal steps.
10. Concatenate all the outputs into a single GFF3 file that is sorted by coordinates and without the lines starting with #.

## DETAILS:
- **The name of the sequences in the FASTA file cannot contain unusual characters.**
- **The input genome file must be a gzip-compressed FASTA file. (.fa.gz)**
- **Auto-train the parameter file always.**


## DOING:
-  Optimizing the memory and CPU requirements for it to run smoothly on the cluster.
-  Optimization of the auto-training additional parameters missing


## PENDING:
- Define accurate memory limits (dynamic based on genome sizes?).
- Tune DIAMOND parameters to make the most of the resources available and to adapt to the capacity of each computer.
- DIAMOND now uses a lot of RAM memory, we have to adjust the execution to reduce the amount of resources used. This may cause an increase in the execution time.
