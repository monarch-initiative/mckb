from dipper import curie_map
from mckb.sources.CGD import CGD
from rdflib.namespace import URIRef
from dipper.utils.CurieUtil import CurieUtil
import unittest
import logging
import datetime
import re

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
        Using test data set 1, and the function add_genotype_info_to_graph()
        We want to test the following triples:

        MONARCH:GenotypeID is an instance of OBO:SO_0001059
        MONARCH:GenotypeID is an instance of OBO:SO_0001583
        MONARCH:GenotypeID has the label "CSF3R Q741X  missense mutation"
        MONARCH:GenotypeID is_sequence_variant_instance_of (OBO:GENO_0000408) NCBIGene:1441
        MONARCH:GenotypeID has location (faldo:location) MONARCH:PositionID
        MONARCH:GenotypeID OBO:SO_transcribed_to MONARCH:TranscriptID

        MONARCH:TranscriptID is an instance of OBO:GENO_primary
        MONARCH:TranscriptID has the label "CCDS413.1"
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
        aa_position_id = self.cgd.make_id('cgd-aa-pos{0}{1}'.format(genotype_key, amino_acid_variant))
        genotype_uri = URIRef(cu.get_uri(genotype_id))
        transcript_uri = URIRef(cu.get_uri(transcript))
        gene_uri = URIRef(cu.get_uri(gene_id))
        aa_position_uri = URIRef(cu.get_uri(aa_position_id))

        sparql_query = """
                       SELECT ?genotype ?gene ?position ?transcript
                       WHERE {{
                           ?genotype a OBO:SO_0001059;
                               a OBO:SO_0001583 ;
                               rdfs:label "{0}" ;
                               OBO:GENO_0000408 ?gene ;
                               faldo:location ?position ;
                               OBO:SO_transcribed_to ?transcript .
                           ?transcript a OBO:GENO_primary ;
                               rdfs:label "{1}" .
                       }}
                       """.format(genotype_label, transcript_id)

        # Expected Results
        expected_results = [[genotype_uri, gene_uri, aa_position_uri, transcript_uri]]
        # Query graph
        sparql_output = test_env.query_graph(sparql_query)

        self.assertEqual(expected_results, sparql_output)

    def test_missense_variant_cdna_model(self):
        """
        Test missense variant with cdna information
        Using test data set 2, and the function add_genotype_info_to_graph()
        We want to test the following triples:

        MONARCH:GenotypeID is an instance of OBO:SO_0001059
        MONARCH:GenotypeID is an instance of OBO:SO_0001583
        MONARCH:GenotypeID has the label "ABL1 T315I missense mutation"
        MONARCH:GenotypeID is_sequence_variant_instance_of (OBO:GENO_0000408) NCBIGene:25
        MONARCH:GenotypeID has location (faldo:location) MONARCH:PositionID
        MONARCH:GenotypeID OBO:SO_transcribed_to MONARCH:TranscriptID

        MONARCH:TranscriptID is an instance of OBO:GENO_secondary
        MONARCH:TranscriptID has the label "CCDS35166.1"
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
        aa_position_id = self.cgd.make_id('cgd-aa-pos{0}{1}'.format(genotype_key, amino_acid_variant))
        genotype_uri = URIRef(cu.get_uri(genotype_id))
        transcript_uri = URIRef(cu.get_uri(transcript))
        gene_uri = URIRef(cu.get_uri(gene_id))
        aa_position_uri = URIRef(cu.get_uri(aa_position_id))

        sparql_query = """
                       SELECT ?genotype ?gene ?position ?transcript
                       WHERE {{
                           ?genotype a OBO:SO_0001059;
                               a OBO:SO_0001583 ;
                               rdfs:label "{0}" ;
                               OBO:GENO_0000408 ?gene ;
                               faldo:location ?position ;
                               OBO:SO_transcribed_to ?transcript .
                           ?transcript a OBO:GENO_secondary ;
                               rdfs:label "{1}" .
                       }}
                       """.format(genotype_label, transcript_id)

        # Expected Results
        expected_results = [[genotype_uri, gene_uri, aa_position_uri, transcript_uri]]
        # Query graph
        sparql_output = test_env.query_graph(sparql_query)

        self.assertEqual(expected_results, sparql_output)

    def test_amino_acid_position_region_model(self):
        """
        Test modelling of amino acid positions


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

        transcript = self.cgd.make_id('cgd-transcript{0}'.format(transcript_id))
        aa_position_id = self.cgd.make_id('cgd-aa-pos{0}{1}'.format(genotype_key, amino_acid_variant))
        region_id = ":_{0}Region".format(aa_position_id)

        transcript_uri = URIRef(cu.get_uri(transcript))
        aa_position_uri = URIRef(cu.get_uri(aa_position_id))
        region_uri = URIRef(cu.get_uri(region_id))

        # Get position
        amino_acid_regex = re.compile(r'^p\.([A-Za-z]{1,3})(\d+)([A-Za-z]{1,3})$')
        match = re.match(amino_acid_regex, amino_acid_variant.rstrip())
        position = match.group(2)


        sparql_query = """
                       SELECT ?position ?region ?bsPosition ?transcript
                       WHERE {{
                           ?position a faldo:Position ;
                               rdfs:label "{0}" ;
                               faldo:location ?region .
                           ?region a faldo:Region ;
                               faldo:begin ?bsPosition ;
                               faldo:end ?bsPosition .
                           ?region a faldo:BothStrandPosition ;
                               a faldo:Position ;
                               faldo:position {1} ;
                               faldo:reference ?transcript .
                       }}
                       """.format(amino_acid_variant, position)

        # Expected Results
        expected_results = [[aa_position_uri, region_uri, transcript_uri]]
        print(region_uri)
        # Query graph
        sparql_output = test_env.query_graph(sparql_query)

        self.assertEqual(expected_results, sparql_output)

if __name__ == '__main__':
    unittest.main()