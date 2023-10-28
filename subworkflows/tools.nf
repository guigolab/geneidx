/*
 * Use gffread to get the sequence of the introns
 */
 process getFASTA {

     // indicates to use as a container the value indicated in the parameter
     container "quay.io/biocontainers/gffread:0.12.7--hd03093a_1"

     // show in the log which input file is analysed
     tag "${name}"

    //  // indicates to use as a label the value indicated in the parameter
    //  label (params.LABEL)

     input:
     tuple val(id), path(gff3), path(genome)

     output:
     tuple val(id), path(name)

     script:
     name = "${id}_gffread.fa"
     """
     gffread -x ${name} -g ${genome} ${gff3};
     """
}


process parseFastaHeader {

    tag "${id}"

    input:
    tuple val(id), path(fasta)

    output:
    tuple val(id), path(out_fasta)

    script:
    out_fasta = "parsed_${fasta.BaseName}.fa"
    
    """
    #!/usr/bin/env python3

    with open("${fasta}", 'r') as infile, open("${out_fasta}", 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                header = line.strip()
                if 'ENA' in header:
                    parts = header.split('|')
                    if len(parts) > 2:
                        new_header = f'>{parts[-1].split()[0]}\\n'
                        outfile.write(new_header)
                    else:
                        outfile.write(header + '\\n')
                else:
                    outfile.write(header + '\\n')
            else:
                outfile.write(line)
    """
}