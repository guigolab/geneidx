//
// This file holds several functions specific to the workflow/genomeannotator.nf in the nf-core/genomeannotator pipeline
//

class WorkflowGenomeannotator {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log) {

        //genomeExistsError(params, log)
        
        //  TODO revise if we need to add some other parameter here

        if (!params.tsv) {
            log.error "TSV input file not specified with e.g. '--tsv input_genomes.tsv'"
            System.exit(1)
        }

        // if (params.assembly.contains('*')) {
        //     log.error "This pipeline is not currently designed to annotate multiple assemblies in parallel. Please start separate pipeline runs instead."
        //     System.exit(1)
        // }

        // if (!params.taxid) {
        //     log.error "Taxid of the genome assembly not specified with e.g. '--assembly 9606'\nFind the taxid required by searching for your species here: https://www.ncbi.nlm.nih.gov/taxonomy/"
        //     System.exit(1)
        // }
    }

    //
    // Exit pipeline if incorrect --genome key provided
    //
    private static void genomeExistsError(params, log) {
        if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
            log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
                "  Currently, the available genome keys are:\n" +
                "  ${params.genomes.keySet().join(", ")}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            System.exit(1)
        }
    }
}
