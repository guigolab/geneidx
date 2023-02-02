
process createDB {

    publishDir(params.OUTPUT, mode : 'copy', pattern : '*.dmnd')

    label 'diamond'

    tag "building ${file_name} database"

    input:
    val(file_name)
    path(file_path)

    output:
    path ("${file_name}.dmnd")

    script:
    """
        echo "Building ${file_name}.dmnd database"
        diamond makedb --in ${file_path} -d ${file_name};
    """
}


//get hsp in gff3
process runDiamond {

    // where to store the results and in which way
    publishDir(params.OUTPUT, mode : 'copy', pattern : '*.gff')

    label "diamond"

    // show in the log which input file is analysed
    tag "${reference_genome_name} against ${dmnd_database_name}"

    input:
    tuple val(dmnd_database_name), path(dmnd_database_file)
    tuple val(reference_genome_name), path(reference_genome_file)

    output:
    path ("${reference_genome_name}.${dmnd_database_name}.hsp.gff")

    script:
    """
    if [ ! -s ${params.OUTPUT}/${reference_genome_name}.${dmnd_database_name}.hsp.gff ]; then
        echo "Running matches ${reference_genome_name}.${dmnd_database_name}"
        fmt6_custom='6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qstrand qframe'
        diamond blastx --db ${dmnd_database_file} \
                       --query ${reference_genome_file} \
                       --max-target-seqs 0 \
                       --max-hsps 0  \
                       --outfmt \$fmt6_custom \
                       --evalue 0.0001 \
                       --block-size 0.8 \
                       --threads ${task.cpus} \
                       --out ${reference_genome_name}.${dmnd_database_name}.hsp.out

        awk 'BEGIN{OFS="\t"}{if (\$13=="-"){frame=-(\$14);print \$1,"blastx","hsp",\$8,\$7,\$12,\$13,frame,\$2}else if(\$13=="+"){print \$1,"blastx","hsp",\$7,\$8,\$12,\$13,\$14,\$2}}' \
                    ${reference_genome_name}.${dmnd_database_name}.hsp.out | sort -k1,1 -k4,5n -k9,9 \
                    | egrep -vw '^MT' | cut -d '|' -f1 > ${reference_genome_name}.${dmnd_database_name}.hsp.gff

    else
        echo "${reference_genome_name}.${dmnd_database_name}.hsp.gff already built"
        ln -s ${params.OUTPUT}/${reference_genome_name}.${dmnd_database_name}.hsp.gff ${reference_genome_name}.${dmnd_database_name}.hsp.gff;
    fi
    rm ${reference_genome_file}
    """
}

process getProteinFileName {


    label "geneidx"

    tag "${taxon}"

    input:
    val lower_lim_proteins
    val upper_lim_proteins

    output:
    stdout emit: filename

    script:
    taxon = params.taxid
    """
    #!/usr/bin/env python3
    # coding: utf-8

    import pandas as pd
    import requests
    from lxml import etree
    import time    

    sp_taxon_id = ${taxon}

    lower_lim = ${lower_lim_proteins}
    upper_lim = ${upper_lim_proteins}

    identity = ${params.uniref_identity}
    ini_clu_size = 24
    clu_size = ini_clu_size

    max_iterations = 30

    def parse_taxon_mod(xml):
        root = etree.fromstring(xml)
        species = (root[0].attrib["rank"], root[0].attrib["taxId"])
        lineage = []
        for taxon in root[0]:
            if taxon.tag == 'lineage':
                lin_pos = 0
                for node in taxon:
                    if "rank" in node.attrib.keys():
                        lineage.append( (node.attrib["rank"], node.attrib["taxId"]) )
                    else:
                        lineage.append((lin_pos, node.attrib["taxId"]))
                    lin_pos += 1

        lineage.insert(0,species)
        return lineage

    response = requests.get(f"https://www.ebi.ac.uk/ena/browser/api/xml/{sp_taxon_id}?download=false")

    lineage_l = parse_taxon_mod(response.content)
    rank_l = []
    taxid_l = []
    for a, b in lineage_l:
        rank_l.append(a)
        taxid_l.append(b)

    if rank_l.index("class"):
        pos_class = rank_l.index("class")
    else:
        pos_class = len(rank_l)//2

    taxon = taxid_l[pos_class]
    # print(taxon)


    # build initial query
    search_type = "search"
    format_type = "list"
    initial_query = "https://rest.uniprot.org/uniref/{}?format={}".format(search_type, format_type)


    # define parameters controlling the loop
    i=0
    in_range = False

    counter = 0
    while not in_range and i < max_iterations:
        query = initial_query + "&query=(taxonomy_id:{})%20AND%20(identity:{})%20AND%20(count:[{}%20TO%20*])".format(taxon, identity, clu_size)
        counter = counter + 1
        if counter > 3:
            time.sleep(2)
            counter= 0
        r = requests.get(query)
        # if we have more proteins than the minimum required
        if lower_lim <= int(r.headers["X-Total-Results"]):

            # check if we are inside the target range of proteins
            if int(r.headers["X-Total-Results"]) <= upper_lim:
                in_range = True

            # if we are above the target range, increase the cluster size required, this will reduce the number of proteins
            elif clu_size <= 50:
                clu_size += 3

            # if we already have a very big cluster size, go one taxonomic rank lower to reduce the diversity of proteins available
            else :
                # get one taxonomic rank lower
                pos_class = pos_class - 1
                taxon = taxid_l[pos_class]
                clu_size = ini_clu_size

        # as we don't reach the minimum number of proteins, we need to get less stringent in terms of cluster size
        elif clu_size >= 6:
            clu_size -= 3

        # if even when reducing cluster size a lot, we don't get enough proteins, let's get one taxonomic rank higher.
        else:
            # get one taxonomic rank higher
            pos_class = pos_class + 1
            taxon = taxid_l[pos_class]
            clu_size = ini_clu_size

        # update i as we want to have a maximum number of queries
        # do we want this maximum number??
        i += 1

    # this is included to have a solution for the worst case scenario
    if not in_range:
        taxon = "2759" # eukaryotes taxid
        clu_size = 60  # change appropriately

    taxon = str(taxon)
    clu_size = str(clu_size).rstrip()
    print(clu_size)
    """

}

process getUniRefQuery {

    label "geneidx"

    tag "${file_name}"

    input:
    val file_name

    output:
    stdout emit: query

    script:
    """
    #!/usr/bin/env python3
    # coding: utf-8

    taxon, clu_size = "${file_name}".split(",")
    identity = ${params.uniref_identity}
    search_type = "stream"
    format_type = "fasta"
    initial_query = "https://rest.uniprot.org/uniref/{}?format={}".format(search_type, format_type)
    query = initial_query + "&compressed=true&query=(taxonomy_id:{})%20AND%20(identity:{})%20AND%20(count:[{}%20TO%20*])".format(taxon, identity, clu_size)

    print(query)    
    """
}



workflow diamond_blastx {
    
    take:
    genome

    main:

    name = getProteinFileName(params.proteins_lower_lim, params.proteins_upper_lim)

    name.view()

    path = getUniRefQuery(name) 
    
    // | splitFasta( file: params.OUTPUT )

    path.view()

    protein_db = createDB( name, path )

    hsp_found = runDiamond( protein_db.map { file -> tuple(file.baseName, file)}, genome.map { file -> tuple(file.baseName, file)} )
    
    emit:
    hsp_found
}
