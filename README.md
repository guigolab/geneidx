# Geneid in Nextflow

Using the [Nexflow tutorial for the Elixir workflow workshop](https://nextflow-io.github.io/elixir-workshop-21/docs/) taking place on November 2021, as an initial template, we are integrating Geneid into Nextflow so that it can be run easily in parallel.


## Schema
1. Make sure geneid is installed, and fastafetch and fastaindex are also available. (all this will go in a container)
2. **Receive the FASTA file and the parameter file as input.**
3. **Uncompress the FASTA file (keeping the compressed file).**

4. MISSING: RepeatMasking!! I need to work on this as it is important for producing a good annotation.

5. **Run the geneidWorkflow using the FASTA file and the parameter file selected.**
6. **Remove the uncompressed FASTA file.**
7. **Provide a single GFF3 file as output of the run.** (it could be sorted and without the lines starting with "#")

This is the schema of the process, but we should also take care of:
- Making it suitable for running in the cluster.
- Choosing the number of CPUs and the memory.
- Explore if it would be possible to use Nextflow for splitting the FASTA file in each independent sequence so that more instances of the geneid predictions can go in parallel. If it is, we should change some steps of the pipeline. We have an older implementation that tried to take advantage of this but we did not manage to make it work.
- https://hub.docker.com/r/guigolab/geneid/tags
  It looks like there is a container of Geneid here, I don't know if we could take it from here or if it would be good to put it into Biocontainers

DETAILS:
- **Adjust the name of the sequences in the FASTA file we download from the ENA, so that it contains the sequence name (chromosome number) and not the identifier of the sequence.**
- Explore the possibility of speeding up the download and uncompression steps in order to make the most of the CPUs that will be requested by the process.
- Define accurate memory limits (dynamic or static based on genome sizes?)


PENDING:
- Discuss the criteria used for selecting a genome to annotate. By now we only select chromosome level assemblies and we only download the sequences belonging to a chromosome.
- Automatically fetch new species and decide which ones should be annotated
- **Choose the parameter file appropriately:**
	- **Have a dictionary of taxid2parameter file.**
	- **Using the taxid of the species to annotate, identify the closest parameter file.**
- Keep track somewhere of which species&assemblies were run, with which parameter files and also the statistics of the run.

(the sections that have been implemented are in **bold**)
