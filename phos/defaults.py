from phos.util import read_json
from os.path import dirname, join, isfile
ROOT_DIR = dirname(dirname(__file__))
WORK_DIR = join(ROOT_DIR, 'work')
EXAMPLE_DIR = join(ROOT_DIR, 'example')
_DB_DIR = join(ROOT_DIR, 'db')
parameters = read_json(join(ROOT_DIR, 'parameters.json'))
descriptions = {
        'go' : 'Gene Ontology enrichment',
        'max_category_size' : 'Maximal category size',
        'min_category_size' : 'Minimal category size',
        'max_category_depth' : 'Maximal ontolgy-tree depth',
        'anat' : 'Network reconstruction',
        'anchor' : 'Signaling source (anchor)',
        'alpha' : 'Global-local optimization tradeoff ($\\alpha$)',
        'activity' : 'Signaling functionality scoring',
        'permutations' : 'Permutations ($r$)',
        'min_size' : 'Required number of observations ($n_0$)',
        'ppi-network' : 'Protein-protein interaction network',
        'degree_threshold' : 'Maximal protein degree',
        'confidence' : 'Edge-confidence threshold'}

steinprt = join(ROOT_DIR, 'bin', 'steinprt')
db = {
        "geneinfo" : join(_DB_DIR, 'geneinfo', 'Homo_sapiens.gene_info.gz'),
        "uniprot_mapping" : join(_DB_DIR, 'uniprot',
            'HUMAN_9606_idmapping.dat.gz'),
        'uniprot_humsavar' : join(_DB_DIR, 'uniprot', 'humsavar.txt'),
        'ppi-network' : join(_DB_DIR, 'anat', 'H_sapiens.net'),
        'go' : {
            'ontology' : join(_DB_DIR, 'go', 'go-basic.obo'),
            'annotations' : join(_DB_DIR, 'go', 'gene2go.gz')
            },
        'phosphosite' : {
            'all_sites' : join(_DB_DIR, 'phosphosite',
                'Phosphorylation_site_dataset.gz'),
            'disease_associated_sites' : join(_DB_DIR, 'phosphosite',
                'Disease-associated_sites.csv'),
            'regulatory_sites' : join(_DB_DIR, 'phosphosite',
                'Regulatory_sites.csv'),
            'kinase_substrate_interactions' : join(_DB_DIR, 'phosphosite',
                'Kinase_Substrate_Dataset.gz')
            }
        }
