import requests
from lxml import etree

def parse_related_projects(root_xml):
    related_projects = list()
    for project in root_xml[0][5]:
        if project[0].tag == 'CHILD_PROJECT':
            related_projects.append(project[0].attrib['accession'])
    return related_projects

def get_xml(accession):
    resp = requests.get(f'https://www.ebi.ac.uk/ena/browser/api/xml/{accession}?download=true')
    return etree.fromstring(resp.content)

def get_assembly(project_accession):
    resp = requests.get(f'https://www.ebi.ac.uk/ena/portal/api/links/study?accession={project_accession}&format=JSON&result=assembly')
    if resp.status_code != 200:
        return
    return resp.json()[0]['accession']

def get_chromosomes(element):
    chromosomes = list()
    for t in element:
        chromosomes.append(t.attrib['accession'])
    return chromosomes
# TODO create assemblies crawler (caveats -> link taxid to assemblies)
# PRJEB33226 25 genomes
# PRJEB40665 DTOL
root_xml = get_xml('PRJEB40665')
LIMIT = 3
primary_related_projects = parse_related_projects(root_xml)
secondary_related_projects = list()
taxons = list()
for project in primary_related_projects[:LIMIT]:
    taxon = dict()
    secondary_root_xml = get_xml(project)
    taxon['taxid'] = secondary_root_xml[0][4][0][0].text
    taxon['name'] = secondary_root_xml[0][4][0][1].text
    taxon['assemblies'] = list()
    for p in parse_related_projects(secondary_root_xml):
        assembly = dict()
        resp = get_assembly(p)
        if resp:
            assembly['accession'] = resp
            assembly_root = get_xml(assembly['accession'])
            for tag in assembly_root[0]:
                if tag.tag == 'CHROMOSOMES':
                    assembly['chromosomes'] = get_chromosomes(tag)
                    taxon['assemblies'].append(assembly)
    if len(taxon['assemblies']) > 0:
        taxons.append(taxon)
print(taxons)
