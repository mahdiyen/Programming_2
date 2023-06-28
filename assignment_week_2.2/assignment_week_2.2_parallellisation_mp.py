import multiprocessing as mp
from Bio import Entrez
import ssl
import time

# the next two lines are needed to create an environment in which the 
# ssl doesn't complain about non-existing public keys...
ssl._create_default_https_context = ssl._create_unverified_context

#enter your email here; the one you used to create an api key in step 0
Entrez.email = 'm.nabilou@st.hanze.nl' 

#using the Entrez module from the Biopython library to perform a database link query.
file = Entrez.elink(dbfrom="pubmed",
                   db="pmc",
                   LinkName="pubmed_pmc_refs",
                   id="30049270",
                   api_key='ddd9478cf718e4f5eb5150e7de20f7679d08')
results = Entrez.read(file)
print ('results= ', results)
print ('###############################')

references = [f'{link["Id"]}' for link in results[0]["LinkSetDb"][0]["Link"]]
print ('references= ', references)
print ('###############################')

def article_details(id = '30049270'):
    """
    This function fetches the XML data of a PubMed article using its ID and saves it as a file.
    input: id (str): The PubMed ID of the article to fetch.
    """
    handle = Entrez.efetch(db="pubmed",
                    id= id,
                    retmode="xml",
                    api_key='ddd9478cf718e4f5eb5150e7de20f7679d08')
    xml_data = handle.read()
    print(id)
    
    #create the .xml file about the given id information
    file_name = str(id) +'.xml'
    with open(file_name, "wb") as xml_file:
        xml_file.write(xml_data)

def article_details_with_mp(reference = references[0:9]):
    """
    Fetches the XML data of multiple PubMed articles in parallel using multiprocessing.
    input: reference (list): A list of IDs of the articles to fetch.
    """
    with mp.Pool() as p:
        res = p.map(article_details, reference)


if __name__ == '__main__':
    start_time = time.time()
    article_details_with_mp(reference=references[0:9])
    end_time = time.time()
    execution_time = end_time - start_time
    print("Execution time with multiprocessing:", execution_time, "seconds")