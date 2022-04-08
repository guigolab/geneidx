
/*
 * Get specific rank taxid from species taxid
 */

process get_tax_rank_ID {

    container 'ferriolcalvet/python-modules' // this one has requests but not lxml

    // where to store the results and in which way
    // publishDir(params.OUTPUT)

    // // indicates to use as a container the value indicated in the parameter
    // container params.CONTAINER

    // show in the log which input file is analysed
    tag "Searching ${tax_rank} of ${taxid}"

    input:
    val taxid
    val tax_rank

    output:
    // stdout emit: taxid_of_interest
    stdout

    script:
    """
    #!/usr/bin/env python

    from lxml import etree
    import requests

    def get_rank_taxid(taxon_id):
        response = requests.get(f"https://www.ebi.ac.uk/ena/browser/api/xml/{taxon_id}?download=false")
        if response.status_code == 200:
            root = etree.fromstring(response.content)
            species = root[0].attrib
            lineage = []
            for taxon in root[0]:
                if taxon.tag == 'lineage':
                    for node in taxon:
                        if "rank" in node.attrib.keys():
                            if node.attrib["rank"] == "${tax_rank}":
                                return node.attrib["taxId"]
                        lineage.append(node.attrib["taxId"])
        return lineage[0]

    taxid_of_interest = get_rank_taxid(${taxid})
    print(taxid_of_interest, end = "")
    """
}






/*
*  Download Protein File
*/
process getProteinFile {

    // container 'ferriolcalvet/wget-download'
    container 'ferriolcalvet/wget-download'

    // where to store the results and in which way
    publishDir(params.OUTPUT, mode : 'copy', pattern : '*.fa.gz')

    // show in the log which input file is analysed
    tag "Download UniRef90.${taxid}.150+.fa.gz"

    input:
    val taxid

    output:
    path ("UniRef90.${taxid}.150+.fa.gz")

    script:
    """
    wget 'https://www.uniprot.org/uniref/?query=taxonomy:${taxid}+AND+count:[150 TO *]+AND+identity:0.9&format=fasta&compress=yes' \
                  -O UniRef90.${taxid}.150+.fa.gz
    """
                   // --max-hsps 0  \
    // mkdir -p tmp;
    // --tmpdir ./tmp/ \
    // https://www.metagenomics.wiki/tools/blast/evalue
}
