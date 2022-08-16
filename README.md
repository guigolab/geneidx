# Geneid+BLASTx
We provide a fast and accurate prediction of the protein-coding genes in a genome, taking as input the genome assembly and a set of proteins from closely related species.

This repository contains the implementation of the Geneid+BLASTx *ab initio* gene prediction method.

In the description here, you can find our preliminary results, a schema of our method and a description of the minimal requirements and commands required for running it.

Stay tuned for an article with detailed descriptions and feel free to [contact us](mailto:ferriol.calvet@crg.eu) if you are trying it and find any problem.


## Preliminary results
 The results of an initial benchmarking using vertebrates genomes annotated in Ensembl show that our method is **as accurate as the top *ab initio* gene predictors**, but it is **between 10 and 100 times faster**.
![Summary of vertebrate benchmark](images/Benchmarking4GitHub.svg)



## Running Geneid+BLASTx
Having defined the parameters and ensuring that Nextflow and a container technology.

`nextflow run main.nf -with-(docker|singularity)`


## Before running Geneid+BLASTx
1. **Make sure ( Docker or Singularity ) and Nextflow are installed in your computer.**
  - [Docker](https://docs.docker.com/engine/install/)
  - [Singularity](https://sylabs.io/guides/3.0/user-guide/installation.html#)
  - [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation)

2. Define the input parameters:
  - Compressed FASTA file with `.fa.gz` termination.
  - Taxid of the species to annotate or from the closest relative described.
  - Proteins from closely related species: (choose one of these options)
      - Provide a compressed FASTA file with the proteins selected.
      - Let the pipeline download the proteins automatically from UniRef90 (nothing should be provided in this case.)
  - Tune parameters for the gene prediction process.
      - Indicate the PWMs defining the start, end, donor and acceptor sites. (this can be obtained from the parameter files listed here: [GENEID parameter files](https://genome.crg.es/software/geneid/index.html#parameters))



## Schema
Which steps are taking place as part of this pipeline?
This is a graphical summary, and the specific steps are outlined below.
![Summary of vertebrate benchmark](images/SchemaWhite.png)
1. Uncompress and index the FASTA file.
2. Get the set of proteins to be used for the protein-to-genome alignments.
3. Create the protein database for DIAMOND to use it as a source.
4. Align the provided genome against the database created using DIAMOND BLASTx flavour.
5. Run the auto-training process:
  - Use matches to estimate the coding sections, look for open reading frames between stop codons.
  - Use matches from the same protein to predict the potential introns.
  - From the sequences of both previous steps compute the initial and transition probability matrices required for the computation of the coding potential of the genome that will be annotated.
6. Update the parameter file with the parameters indicated in the *params.config* file and also the matrices automatically generated from the protein-DNA matches.
7. Run Geneid with the new parameter file and the protein-DNA matches from the previous steps as additional evidence.
This is done in parallel for each independent sequence inside the genome FASTA file, and it consists of the following steps:
  - Pre-process the matches to join those that overlap and re-score them appropriately.
  - Run Geneid using the evidence and obtain the GFF3 file.
  - Remove the files from the internal steps.
8. Concatenate all the outputs into a single GFF3 file that is sorted by coordinates.


## DETAILS:
- **The name of the sequences in the FASTA file cannot contain unusual characters.**
- **The input genome file must be a gzip-compressed FASTA file. (.fa.gz)**
- **Auto-train the parameter file always.**


## DOING:
-  Optimizing the memory and CPU requirements for it to run smoothly on the cluster.
-  Optimization of the auto-training additional parameters missing.


## PENDING:
- Define accurate memory limits (dynamic based on genome sizes?).
- Tune DIAMOND parameters to make the most of the resources available and to adapt to the capacity of each computer.
- DIAMOND now uses a lot of RAM memory, we have to adjust the execution to reduce the amount of resources used. This may cause an increase in the execution time.


Follow us on Twitter ([@GuigoLab](https://twitter.com/GuigoLab)) for updates in the article describing our method.
