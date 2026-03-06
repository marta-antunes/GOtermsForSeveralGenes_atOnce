from Bio import Entrez  #required to access ncbi API
import requests  #required to access the Gene Ontology (GO) API
import csv  #required to read csv files


def eSearch(terms):
    """search for gene Ids given a query
    Requires: term is string like "Escherichia coli[orgn]+rpoS[Gene Name]"
    Ensures: gene ID"""

    GeneIDs = 0
    
    Entrez.email = "Your.Name.Here@example.org"
    handle = Entrez.esearch(db="gene", retmax= 1, term = terms)

    import xml.etree.ElementTree as ET
    handle_as_string = handle.read()
    root = ET.fromstring(handle_as_string)
    
    for neighbor in root.iter('Id'):
        GeneIDs=neighbor.text

    handle.close()
    return (GeneIDs)


def get_go_terms_by_gene_id(gene_id):
    """search for GO terms given a gene Id
    Requires: gene Id is int
    Ensures: GO terms are strings"""
    url = f"https://api.geneontology.org/api/bioentity/gene/{gene_id}/function"
    headers = {
        "Accept": "application/json"
    }

    try:
        response = requests.get(url, headers=headers)
        response.raise_for_status()  # Raise an error for non-2xx status codes
        data = response.json()
        return data
    except requests.exceptions.RequestException as e:
        print("Error:", e)
        return None

def get_flybase_id(gene_id):
    """search for flybase Ids given a gene id
    Requires: gene id is int
    Ensures: flybase id"""
    Entrez.email = "Your.Name.Here@example.org"

    try:
        handle = Entrez.esearch(db="gene", term=gene_id)
        record = Entrez.read(handle)
    except Exception as e:
        print(f"Error fetching data for gene_id {gene_id}: {e}")
        return None

    if record['Count'] == '0':
        return None

    gene_id_list = record['IdList']
    gene_id = gene_id_list[0]  # Assuming only one gene found

    try:
        handle = Entrez.efetch(db="gene", id=gene_id, retmode="xml")
        record = Entrez.read(handle)
    except Exception as e:
        print(f"Error fetching data for gene_id {gene_id}: {e}")
        return None
    #print(record)

    for data in record:
        if 'Entrezgene_gene' in data and 'Gene-ref' in data['Entrezgene_gene']:
            gene_ref = data['Entrezgene_gene']['Gene-ref']
            if 'Gene-ref_db' in gene_ref:
                for db in gene_ref['Gene-ref_db']:
                    if 'Dbtag_db' in db and db['Dbtag_db'] == 'FLYBASE' and 'Dbtag_tag' in db and 'Object-id' in db[
                        'Dbtag_tag']:
                        flybase_id = db['Dbtag_tag']['Object-id']['Object-id_str']
                        #print("flybase:" flybase_id)
                        return(flybase_id)
    return None


def get_go_terms_by_flybase_id(go_domain, fbgn_id):
    """search for go terms given a flybase id and a GO category e.g biological_process
    Requires: flybase_id and fbgn_id are strings
    Ensures: GO terms"""
    url = f"https://api.flybase.org/api/v1.0/ribbon/go/{go_domain}/{fbgn_id}"
    headers = {
        "Accept": "application/json"
    }

    try:
        response = requests.get(url, headers=headers)
        response.raise_for_status()  # Raise an error for non-2xx status codes
        data = response.json()
        return data
    except requests.exceptions.RequestException as e:
        print("Error fetching data:", e)
        return None
    except ValueError as e:
        print("Error decoding JSON:", e)
        return None


def get_SGD_id(gene_id):
    # Entrez email
    Entrez.email = "your@email.com"

    # Entrez API query to fetch data
    try:
        handle = Entrez.efetch(db="gene", id=gene_id, retmode="xml")
        # Parse the XML data
        record = Entrez.read(handle)
    except Exception as e:
        print(f"Error fetching data for gene_id {gene_id}: {e}")
        return None


    # Extract SGD ID
    sgd_id = None
    for entry in record:
        for field in entry.get('Entrezgene_unique-keys', []):
            if field.get('Dbtag_db') == 'SGD':
                sgd_id = field.get('Dbtag_tag', {}).get('Object-id', {}).get('Object-id_str')
                break

    # Output the SGD ID
    return (sgd_id)

def fetch_gene_go_terms(sgd_id):
    url = f"https://stage.alliancegenome.org/alliancemine/service/template/results?name=Gene_GOTerms&constraint1=Gene&op1=LOOKUP&value1={sgd_id}&extra1=&format=tab&size=10"
    response = requests.get(url)

    if response.status_code == 200:
        # Split the response text into lines and return as list of lists
        data=[line.split('\t') for line in response.text.strip().split('\n')]
        if data:
            return(data)

    else:
        print("Error:", response.status_code)
        return None

def main():
    """gets a file name corresponding to a txt file with one gene name per line"""
    fileName = input("Enter the file name: ")
    output_file = fileName.strip(".csv") + "_output.csv"  # Output file name
    print(f"Processing file: {fileName}")

    # Open the tab-delimited file
    with open(fileName, 'r', encoding='ISO-8859-1') as file, open(output_file, 'w', encoding='ISO-8859-1') as out_file:
        # Create a CSV reader object with comma delimiter
        reader = csv.reader(file, delimiter=',')
        #reader = csv.reader(file, delimiter='')

        # Skip the first four rows
        next(reader)
        next(reader)
        next(reader)
        next(reader)

        # Iterate through each row in the file
        for row in reader:
            # Access the 3rd column of the current row that contains the species name
            column_3 = row[2]  # Python is zero-indexed

            # Access the 20th column of the current row that contains several gene names
            column_20 = row[19]  # Python is zero-indexed

            # Access the 24th column of the current row that contains the reference
            column_24 = row[23]  # Python is zero-indexed

            # Split the cell content by comma
            elements_in_column_20 = column_20.split('+')

            # Iterate over the elements
            for element in elements_in_column_20:
                print(f"Processing element: {element}")
                listOfDescendentTerm=[]
                species_and_gene = column_3 + "[orgn]+" + element.strip("\n") + "[Gene Name]"  #spaces.strip("\n")
                print(f"Searching for gene id for: {species_and_gene}")
                geneID = eSearch(species_and_gene)  # search for gene id
                if geneID == 0:  # situations when GeneID was not found
                    print(f"No gene id found for: {species_and_gene}")
                    out_file.write(species_and_gene + "," + str(geneID) + "," + "NA" + "," + "NA" + "," + "NA" + "," + column_24 + "\n")
                else:
                    print(f"Found gene id: {geneID} for: {species_and_gene}")
                    if column_3 == "Drosophila melanogaster":
                        FB = get_flybase_id(geneID)
                        if FB is None:
                            print(f"Skipping geneID {geneID} due to error or no data found")
                            continue
                        if FB == 0:
                            out_file.write(species_and_gene + "," + str(geneID) + "," + "NA" + "," + "NA" + "," + "NA" + "," + column_24 + "\n")
                        else:
                            go_terms = get_go_terms_by_flybase_id("biological_process", FB)
                            if go_terms is not None:
                                for result in go_terms['resultset']['result']:
                                    for slim_id, slim_info in result['ribbon'].items():
                                        for descendant_term in slim_info['descendant_terms']:
                                            listOfDescendentTerm.append(descendant_term)
                                            go_term = descendant_term['name']
                                            go_id = descendant_term['id']
                                            if go_term[8:].lstrip().startswith("annotation(s) using"):
                                                pass
                                            else:
                                                out_file.write(species_and_gene + "," + str(geneID) + "," + go_id + "," + go_term[8:] + "," + "biological_process" + "," + column_24 + "\n")
                            else:
                                pass
                            go_terms_MF = get_go_terms_by_flybase_id("molecular_function", FB)
                            if go_terms_MF is not None:
                                for result in go_terms_MF['resultset']['result']:
                                    for slim_id, slim_info in result['ribbon'].items():
                                        for descendant_term in slim_info['descendant_terms']:
                                            listOfDescendentTerm.append(descendant_term)
                                            go_term = descendant_term['name']
                                            go_id = descendant_term['id']
                                            if go_term[8:].lstrip().startswith("annotation(s) using"):
                                                pass
                                            else:
                                                out_file.write(species_and_gene + "," + str(geneID) + "," + go_id + "," + go_term[8:] + "," + "molecular_function"+ "," + column_24 + "\n")
                            else:
                                pass
                            if listOfDescendentTerm == []:  # when all descendent terms are empty
                                out_file.write(species_and_gene + "," + str(geneID) + "," + "NA" + "," + "NA" + "," + "NA"+ "," + column_24 + "\n")
                    elif column_3 == "Saccharomyces cerevisiae":
                        SGD = get_SGD_id(geneID)
                        if SGD == 0:
                            out_file.write(species_and_gene + "," + str(geneID) + "," + "NA" + "," + "NA" + "," + "NA" + "," + column_24 + "\n")
                        else:
                            data_set = fetch_gene_go_terms(SGD)
                            for go_term in data_set:
                                if go_term[4]=="cellular_component":
                                    pass
                                else:
                                    out_file.write(species_and_gene + "," + str(geneID) + "," + go_term[2] + "," + go_term[3] + "," + go_term[4]+ "," + column_24 + "\n")
                    else:
                        go_terms = get_go_terms_by_gene_id(geneID)  # search go terms
                        if all(not value for value in go_terms.values()):  # Check if there are no GO terms
                            out_file.write(species_and_gene + "," + str(geneID) + "," + "NA" + "," + "NA" + "," + "NA" + "," + column_24 + "\n")
                        elif go_terms:
                            for annotation in go_terms['associations']:
                                go_id = annotation['object']['id']
                                go_term = annotation['object']['label']
                                go_category = annotation['object']['category']
                                if 'cellular_component' in go_category:  # ignore the category "cellular component"
                                    pass
                                else:
                                    out_file.write(species_and_gene + "," + str(geneID) + "," + go_id + "," + go_term + "," + str(go_category) + "," + column_24 + "\n")



if __name__ == "__main__":
    main()