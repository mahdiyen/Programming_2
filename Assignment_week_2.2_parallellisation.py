import multiprocessing as mp
from Bio import Entrez
import ssl

# the next two lines are needed to create an environment in which the 
# ssl doesn't complain about non-existing public keys...
ssl._create_default_https_context = ssl._create_unverified_context

#enter your email here; the one you used to create an api key in step 0
Entrez.email = 'm.nabilou@st.hanze.nl' 

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

    handle = Entrez.efetch(db="pubmed",
                    id= id,
                    retmode="xml",
                    api_key='ddd9478cf718e4f5eb5150e7de20f7679d08')
    xml_data = handle.read()
    print(id)
    
    file_name = str(id) +'.xml'
    with open(file_name, "w") as xml_file:
        xml_file.write(xml_data)

def article_details_with_mp(reference = references[0:9]):
    with mp.Pool() as p:
        res = p.map(article_details, reference)

if __name__ == '__main__':
    article_details_with_mp(reference = references[0:9])

