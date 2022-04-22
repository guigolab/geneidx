// Parameter definitions
params.CONTAINER = "ferriolcalvet/filter-genome"


// it would be good to keep track of the restriction that we applied
process filter_Fasta_by_length {

    // where to store the results and in which way
    // publishDir(params.OUTPUT, mode: "copy", pattern : "*.fa")

    // indicates to use as a label the value indicated in the parameter
    // label (params.LABEL)

    // indicates to use as container the value indicated in the parameter
    container params.CONTAINER

    // show in the log which input file is analysed
    tag "filtering ${ref_genome}"


    input:
    path (ref_genome)
    val (min_seq_length)

    output:
    path ("${ref_genome_name}.filtered.fa")

    script:
    ref_genome_name = ref_genome.BaseName
    """
    infoseq -only -name -length -noheading ${ref_genome} > seqs_length;
    seq_threshold=\$(awk '{s+=\$2}END{print s*${min_seq_length}}' seqs_length)
    awk -v var=\$seq_threshold '\$2<var' seqs_length > remove_seqs;
    FastaToTbl ${ref_genome} | grep -Fwvf remove_seqs | TblToFasta > ${ref_genome_name}.filtered.fa;
    rm remove_seqs;
    """
    // infoseq -only -name -length -noheading ${ref_genome} | awk '\$2<${min_seq_length} {print \$1}' > remove_seqs;
    // FastaToTbl ${ref_genome} | grep -Fwvf remove_seqs | TblToFasta > ${ref_genome_name}.${min_seq_length}.fa;
    // rm remove_seqs;

}
