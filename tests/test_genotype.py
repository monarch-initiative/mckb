from dipper import curie_map
from mckb.sources.CGD import CGD
from rdflib.namespace import URIRef
from dipper.utils.CurieUtil import CurieUtil
import unittest
import logging
import datetime

logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)


class DiseaseDrugGenotypeTestCase(unittest.TestCase):
    """
    Test triples created from genotype modelling functions
    """

    def setUp(self):

        self.curie_map = curie_map.get()
        cu = CurieUtil(self.curie_map)
        # Fake credentials as these tests do not require a database connection
        database = 'foo'
        user = 'bar'
        password = 'baz'

        self.cgd = CGD(database, user, password)

        mapping_file = '../../resources/mappings/gene.tsv'
        self.cgd.gene_map = self.cgd.set_gene_map(mapping_file)

        # Sample output from _get_genotype_protein_info() where genotype
        # is a missense mutation
        self.test_set_1 = ((2, 'CSF3R Q741X  missense mutation', 'p.Q741X ',
                            None, 'CCDS413.1', 'Primary', None,
                            'gain-of-function', None, 'CSF3R', None),)

        # Sample output from _get_genotype_cdna_info()
        self.test_set_2 = ((19, 'ABL1 T315I missense mutation', 'p.T315I',
                            315, 'CCDS35166.1', 'Secondary',
                            'nonsynonymous - missense', 'gain-of-function',
                            None, 'ABL1', None, 'ABL1', 944, 'c.944C>T',
                            'COSM12560', 'rs121913459', 133748283, 133748283,
                            'C', 'T', 'Ex6', 'nonsynonymous - missense',
                            'Substitution', 'chr9', 'hg19', 'hg19',
                            datetime.datetime(2009, 2, 1, 0, 0)),)

        return

    def tearDown(self):
        self.cgd.graph = None
        self.cgd = None
        return

    def test_missense_variant_protein_model(self):
        """
        Test missense variant with only protein information
        """
        from dipper.utils.TestUtils import TestUtils

        self.cgd.add_genotype_info_to_graph(self.test_set_1)

        # Make testutils object and load bindings
        test_env = TestUtils(self.cgd.graph)
        cu = CurieUtil(self.curie_map)
        self.cgd.load_bindings()

        (genotype_key, genotype_label, amino_acid_variant, amino_acid_position,
         transcript_id, transcript_priority, protein_variant_type,
         functional_impact, stop_gain_loss, transcript_gene,
         protein_variant_source) = self.test_set_1[0][0:11]

        gene_id = self.cgd.gene_map[transcript_gene]

        genotype_id = self.cgd.make_id('cgd-genotype{0}'.format(genotype_key))
        transcript = self.cgd.make_id('cgd-transcript{0}'.format(transcript_id))
        genotype_uri = URIRef(cu.get_uri(genotype_id))
        transcript_uri = URIRef(cu.get_uri(transcript))
        gene_uri = URIRef(cu.get_uri(gene_id))

        sparql_query = """
                       SELECT ?genotype ?gene ?transcript
                       WHERE {{
                           ?genotype a OBO:SO_0001059;
                               a OBO:SO_0001583 ;
                               rdfs:label "{0}" ;
                               OBO:GENO_0000408 ?gene ;
                               OBO:SO_transcribed_to ?transcript .
                           ?transcript a OBO:SO_0000185 ;
                               rdfs:label "{1}" .
                       }}
                       """.format(genotype_label, transcript_id)

        # Expected Results
        expected_results = [[genotype_uri, gene_uri, transcript_uri]]
        # Query graph
        sparql_output = test_env.query_graph(sparql_query)

        self.assertEqual(expected_results, sparql_output)

    def test_missense_variant_cdna_model(self):
        """
        Test missense variant with only protein information
        """
        from dipper.utils.TestUtils import TestUtils

        self.cgd.add_genotype_info_to_graph(self.test_set_2)

        # Make testutils object and load bindings
        test_env = TestUtils(self.cgd.graph)
        cu = CurieUtil(self.curie_map)
        self.cgd.load_bindings()

        (genotype_key, genotype_label, amino_acid_variant, amino_acid_position,
         transcript_id, transcript_priority, protein_variant_type,
         functional_impact, stop_gain_loss, transcript_gene,
         protein_variant_source) = self.test_set_2[0][0:11]

        gene_id = self.cgd.gene_map[transcript_gene]

        genotype_id = self.cgd.make_id('cgd-genotype{0}'.format(genotype_key))
        transcript = self.cgd.make_id('cgd-transcript{0}'.format(transcript_id))
        genotype_uri = URIRef(cu.get_uri(genotype_id))
        transcript_uri = URIRef(cu.get_uri(transcript))
        gene_uri = URIRef(cu.get_uri(gene_id))

        sparql_query = """
                       SELECT ?genotype ?gene ?transcript
                       WHERE {{
                           ?genotype a OBO:SO_0001059;
                               a OBO:SO_0001583 ;
                               rdfs:label "{0}" ;
                               OBO:GENO_0000408 ?gene ;
                               OBO:SO_transcribed_to ?transcript .
                           ?transcript a OBO:SO_0001596 ;
                               rdfs:label "{1}" .
                       }}
                       """.format(genotype_label, transcript_id)

        # Expected Results
        expected_results = [[genotype_uri, gene_uri, transcript_uri]]
        # Query graph
        sparql_output = test_env.query_graph(sparql_query)

        self.assertEqual(expected_results, sparql_output)

if __name__ == '__main__':
    unittest.main()