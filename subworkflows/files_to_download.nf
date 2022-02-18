import groovy.json.JsonSlurper

/*
*  Geneid module.
*/

// Parameter definitions
// params.CONTAINER = "geneid_path"
// params.OUTPUT = "geneid_output"
// params.LABEL = ""

/*
 * Get assemblies from ENA
 */

process list_files_to_download {

    // container 'python:3.6-stretch'

    // where to store the results and in which way
    publishDir(params.OUTPUT)

    // // indicates to use as a container the value indicated in the parameter
    // container params.CONTAINER

    // show in the log which input file is analysed
    tag "Looking for assemblies in ENA"

    input:
    val accession

    output:
    // val("${link_to_ref}"), emit: fasta_name
    stdout emit: taxons

    script:
    """
    #!/usr/bin/env python
    import os
    import json

    import requests
    from lxml import etree
    import pandas as pd

    LIMIT=2

    # PRJEB33226 25 genomes
    # PRJEB40665 DTOL
    # PRJNA489243 VGP
    # PRJNA312960 200 Mammals
    # PRJNA512907 DNA Zoo

    def parse_related_projects(root_xml):
        secondary_projects = list()
        related_projects = root_xml.findall('.//RELATED_PROJECTS')
        if related_projects:
            for project in related_projects[0]:
                if project[0].tag == 'CHILD_PROJECT':
                    secondary_projects.append(project[0].attrib['accession'])
        return secondary_projects

    def get_xml(accession):
        resp = requests.get(f'https://www.ebi.ac.uk/ena/browser/api/xml/{accession}?download=true')
        return etree.fromstring(resp.content)

    def get_chromosomes_and_taxon(project_accession, taxons):
        resp = requests.get(f'https://www.ebi.ac.uk/ena/portal/api/links/study?accession={project_accession}&format=JSON&result=assembly')
        if resp.status_code != 200:
            return
        else:
            for assembly in resp.json():
                ass_xml = get_xml(assembly['accession'])
                for el in ass_xml[0]:
                    if el.tag == 'CHROMOSOMES':
                        taxon = dict()
                        taxon['assembly'] = assembly['accession']
                        taxon['chromosomes'] = get_chromosomes(el)[-2:]
                        taxon_tag = ass_xml[0].find('TAXON')
                        taxon['taxid'] = taxon_tag.find('TAXON_ID').text
                        taxon['name'] = taxon_tag.find('SCIENTIFIC_NAME').text
                        if not any(tax['taxid'] == taxon['taxid'] for tax in taxons):
                            taxon['param'] = get_param_file(taxon['taxid'])
                            taxons.append(taxon)
            return

    def get_chromosomes(element):
        chromosomes = list()
        for t in element:
            chromosomes.append(t.attrib['accession'])
        return chromosomes

    def find_organism_infos(element):
        organism = element.findall('.//TAXON')
        if organism:
            return organism
        else:
            if len(element) > 0:
                for el in element:
                    find_organism_infos(el)
            else:
                return

    def crawl_assemblies(project_accession, taxons):
        project = get_xml(project_accession)
        related_projects = parse_related_projects(project)
        if related_projects:
            for pr in related_projects[:LIMIT]:
                crawl_assemblies(pr,taxons)
        else:
            get_chromosomes_and_taxon(project_accession, taxons)

    def get_organism(taxon_id):
        response = requests.get(f"https://www.ebi.ac.uk/ena/browser/api/xml/{taxon_id}?download=false")
        if response.status_code == 200:
            root = etree.fromstring(response.content)
            species = root[0].attrib
            lineage = []
            for taxon in root[0]:
                if taxon.tag == 'lineage':
                    for node in taxon:
                        lineage.append(node.attrib["taxId"])
        return lineage

    def get_param_file(taxid):
        table_with_param_data = "$projectDir/data/select-param-files/taxid_to_links.tsv"
        scored_links = pd.read_csv(table_with_param_data, header = 0, sep = "\t")
        query = pd.DataFrame(get_organism(taxid))
        query.columns = ["taxid"]
        query.loc[:,"taxid"] = query.loc[:,"taxid"].astype(int)
        intersected_params = query.merge(scored_links, on = "taxid").sort_values(by = "rank_pos", ascending = False)
        if intersected_params.shape[0] > 0:
            selected_param = str(intersected_params.loc[0,"param_taxid"])
        else:
            selected_param = str(9606)
        #return selected_param
        return f"$projectDir/data/param-files/{selected_param}.param"


    taxons = list()
    crawl_assemblies("$accession", taxons)
    print(json.dumps(taxons))
    """
}


/*
 * Uncompressing if needed
 */

process parse_json {

    // where to store the results and in which way
    publishDir(params.OUTPUT)

    // show in the log which input file is analysed
    tag "Parsing JSON output.."

    input:
    val taxons

    output:
    val result

    exec:
    result = new JsonSlurper().parseText(taxons)
}
