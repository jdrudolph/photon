from phos.util import read_json
from os.path import dirname, join, isfile

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
        'confidence' : 'Edge-confidence threshold',
        'side' : 'Test sidedness (greater, twosided, lesser)'
        }

def make_defaults(root_dir, work='work', example='example', db='db'):
    work_dir = join(root_dir, work)
    example_dir = join(root_dir, example)
    db_dir = join(root_dir, db)
    steinprt = join(root_dir, 'bin', 'steinprt')
    db = { "geneinfo" : join(db_dir, 'geneinfo', 'Homo_sapiens.gene_info.gz'),
            "uniprot_mapping" : join(db_dir, 'uniprot',
                'HUMAN_9606_idmapping.dat.gz'),
            'uniprot_humsavar' : join(db_dir, 'uniprot', 'humsavar.txt'),
            'ppi-network' : join(db_dir, 'anat', 'H_sapiens.net'),
            'go' : {
                'ontology' : join(db_dir, 'go', 'go-basic.obo'),
                'annotations' : join(db_dir, 'go', 'gene2go.gz')
                },
            'phosphosite' : {
                'all_sites' : join(db_dir, 'phosphosite',
                    'Phosphorylation_site_dataset.gz'),
                'disease_associated_sites' : join(db_dir, 'phosphosite',
                    'Disease-associated_sites.csv'),
                'regulatory_sites' : join(db_dir, 'phosphosite',
                    'Regulatory_sites.csv'),
                'kinase_substrate_interactions' : join(db_dir, 'phosphosite',
                    'Kinase_Substrate_Dataset.gz')
                }
            }
    return {'root': root_dir, 'work': work_dir, 'example': example_dir, 'db': db, 'steinprt': steinprt}
