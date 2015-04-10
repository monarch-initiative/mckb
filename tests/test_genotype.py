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
        Using test data set 1, and the function add_genotype_info_to_graph()
        We want to test the following triples:

        MONARCH:GenotypeID is an instance of OBO:SO_0001059
        MONARCH:GenotypeID is an instance of OBO:SO_0001583
        MONARCH:GenotypeID has the label "CSF3R Q741X  missense mutation"
        MONARCH:GenotypeID is_sequence_variant_instance_of (OBO:GENO_0000408) NCBIGene:1441
        MONARCH:GenotypeID has location (faldo:location) MONARCH:PositionID
        MONARCH:GenotypeID OBO:GENO_reference_amino_acid "Q"
        MONARCH:GenotypeID OBO:GENO_results_in_amino_acid_change "X"
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
        ref_amino_acid = "Q"
        altered_amino_acid = "X"

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
                               OBO:GENO_reference_amino_acid "{1}" ;
                               OBO:GENO_results_in_amino_acid_change "{2}" ;
                               OBO:SO_transcribed_to ?transcript .

                           ?transcript a OBO:GENO_primary ;
                               rdfs:label "{3}" .
                       }}
                       """.format(genotype_label, ref_amino_acid,
                                  altered_amino_acid, transcript_id)

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
        MONARCH:GenotypeID has location (faldo:location) MONARCH:PositionID1 (amino acid location)
        MONARCH:GenotypeID has location (faldo:location) MONARCH:PositionID2 (location of chromosome)
        MONARCH:GenotypeID OBO:GENO_reference_amino_acid "T"
        MONARCH:GenotypeID OBO:GENO_results_in_amino_acid_change "I"
        MONARCH:GenotypeID OBO:SO_transcribed_to MONARCH:TranscriptID

        MONARCH:TranscriptID is an instance of OBO:GENO_secondary
        MONARCH:TranscriptID has the label "CCDS35166.1"

        MONARCH:PositionID1 (amino acid location) has the label "p.T315I"
        MONARCH:PositionID2 (chromosome location) hase the label "ABL1 genomic location"
        """
        from dipper.utils.TestUtils import TestUtils

        self.cgd.add_genotype_info_to_graph(self.test_set_2)

        # Make testutils object ande  load bindings
        test_env = TestUtils(self.cgd.graph)
        cu = CurieUtil(self.curie_map)
        self.cgd.load_bindings()

        (genotype_key, genotype_label, amino_acid_variant, amino_acid_position,
         transcript_id, transcript_priority, protein_variant_type,
         functional_impact, stop_gain_loss, transcript_gene,
         protein_variant_source, variant_gene, bp_pos, genotype_cdna,
         cosmic_id, db_snp_id, genome_pos_start, genome_pos_end, ref_base,
         variant_base, primary_transcript_exons,
         primary_transcript_variant_sub_types, variant_type, chromosome,
         genome_build, build_version, build_date) = self.test_set_2[0]

        gene_id = self.cgd.gene_map[transcript_gene]
        ref_amino_acid = "T"
        altered_amino_acid = "I"
        variant_position_label = '{0} genomic location'.format(variant_gene)

        genotype_id = self.cgd.make_id('cgd-genotype{0}'.format(genotype_key))
        transcript = self.cgd.make_id('cgd-transcript{0}'.format(transcript_id))
        aa_position_id = self.cgd.make_id('cgd-aa-pos{0}{1}'.format(genotype_key, amino_acid_variant))
        variant_position_id = self.cgd.make_id(
            'cgd-var-pos{0}{1}'.format(genotype_key, variant_gene))
        genotype_uri = URIRef(cu.get_uri(genotype_id))
        transcript_uri = URIRef(cu.get_uri(transcript))
        gene_uri = URIRef(cu.get_uri(gene_id))
        aa_position_uri = URIRef(cu.get_uri(aa_position_id))
        chr_position_uri = URIRef(cu.get_uri(variant_position_id))

        sparql_query = """
                       SELECT ?genotype ?gene ?aaPosition ?chrPosition ?transcript
                       WHERE {{
                           ?genotype a OBO:SO_0001059;
                               a OBO:SO_0001583 ;
                               rdfs:label "{0}" ;
                               OBO:GENO_0000408 ?gene ;
                               faldo:location ?aaPosition ;
                               faldo:location ?chrPosition ;
                               OBO:GENO_reference_amino_acid "{1}" ;
                               OBO:GENO_results_in_amino_acid_change "{2}" ;
                               OBO:SO_transcribed_to ?transcript .

                           ?transcript a OBO:GENO_secondary ;
                               rdfs:label "{3}" .

                           ?aaPosition rdfs:label "{4}" .
                           ?chrPosition rdfs:label "{5}" .
                       }}
                       """.format(genotype_label, ref_amino_acid,
                                  altered_amino_acid, transcript_id,
                                  amino_acid_variant, variant_position_label)

        # Expected Results
        expected_results = [[genotype_uri, gene_uri, aa_position_uri, chr_position_uri, transcript_uri]]
        # Query graph
        sparql_output = test_env.query_graph(sparql_query)

        self.assertEqual(expected_results, sparql_output)

    def test_amino_acid_position_region_model(self):
        """
        Test modelling of amino acid positions
        Using test data set 1, and the function add_genotype_info_to_graph()
        We want to test the following triples:

        MONARCH:PositionID is an instance of faldo:Position
        MONARCH:PositionID rdfs:label "p.Q741X"
        MONARCH:PositionID faldo:location MONARCH:RegionID

        MONARCH:RegionID is an instance of faldo:Region
        MONARCH:RegionID faldo:begin MONARCH:BothStrandPositionID
        MONARCH:RegionID faldo:end MONARCH:BothStrandPositionID

        MONARCH:BothStrandPositionID is an instance of faldo:BothStrandPosition
        MONARCH:BothStrandPositionID is an instance of faldo:Position
        MONARCH:BothStrandPositionID faldo:position 741
        MONARCH:BothStrandPositionID faldo:reference MONARCH:TranscriptID
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

        position = 741

        transcript = self.cgd.make_id('cgd-transcript{0}'.format(transcript_id))
        aa_position_id = self.cgd.make_id('cgd-aa-pos{0}{1}'.format(genotype_key, amino_acid_variant))
        region_id = ":_{0}Region".format(aa_position_id)
        both_strand_id = ":_{0}-{1}".format(transcript, position)

        transcript_uri = URIRef(cu.get_uri(transcript))
        aa_position_uri = URIRef(cu.get_uri(aa_position_id))
        region_uri = URIRef(cu.get_uri(region_id))
        both_strand_uri = URIRef(cu.get_uri(both_strand_id))

        sparql_query = """
                       SELECT ?position ?region ?bsPosition ?transcript
                       WHERE {{
                           ?position a faldo:Position ;
                               rdfs:label "{0}" ;
                               faldo:location ?region .

                           ?region a faldo:Region ;
                               faldo:begin ?bsPosition ;
                               faldo:end ?bsPosition .

                           ?bsPosition a faldo:BothStrandPosition ;
                               a faldo:Position ;
                               faldo:position {1} ;
                               faldo:reference ?transcript .
                       }}
                       """.format(amino_acid_variant, position)

        # Expected Results
        expected_results = [[aa_position_uri, region_uri, both_strand_uri, transcript_uri]]

        # Query graph
        sparql_output = test_env.query_graph(sparql_query)

        self.assertEqual(expected_results, sparql_output)

    def test_genome_build_chromosome_model(self):
        """
        Test modelling of genome, builds, and chromosomes
        Using test data set 2, and the function add_genotype_info_to_graph()
        """
        from dipper.utils.TestUtils import TestUtils
        self.cgd.add_genotype_info_to_graph(self.test_set_2)

        # Make testutils object and load bindings
        test_env = TestUtils(self.cgd.graph)
        cu = CurieUtil(self.curie_map)
        self.cgd.load_bindings()

        genome = ":9606genome"
        genome_label = "Human genome"
        chromosome = ":9606chr9"
        chromosome_label = "chr9 (Human)"
        build_curie = "UCSC:hg19"
        build_label = "hg19"
        chrom_on_build = ":hg19chr9"
        chrom_build_label = "chr9 (hg19)"

        genome_uri = URIRef(cu.get_uri(genome))
        chromosome_uri = URIRef(cu.get_uri(chromosome))
        build_uri = URIRef(cu.get_uri(build_curie))
        chrom_on_build_uri = URIRef(cu.get_uri(chrom_on_build))

        sparql_query = """
                       SELECT ?genome ?chromosome ?build ?chromOnBuild
                       WHERE {{
                           ?genome a owl:Class ;
                               rdfs:label "{0}" ;
                               OBO:RO_0002162 OBO:NCBITaxon_9606 ;
                               OBO:RO_0002351 ?chromosome ;
                               rdfs:subClassOf OBO:SO_0001026 .

                           ?chromosome a owl:Class ;
                               rdfs:label "{1}" ;
                               OBO:RO_0002350 ?genome ;
                               rdfs:subClassOf OBO:SO_0000340 .

                           ?build a OBO:SO_0001505 ;
                               a ?genome ;
                               rdfs:label "{2}" ;
                               OBO:RO_0002351 ?chromOnBuild ;
                               rdfs:subClassOf ?genome .

                           ?chromOnBuild a ?chromosome ;
                               rdfs:label "{3}" ;
                               OBO:RO_0002350 ?build .
                       }}
                       """.format(genome_label, chromosome_label,
                                  build_label, chrom_build_label)

        # Expected Results
        expected_results = [[genome_uri, chromosome_uri,
                             build_uri, chrom_on_build_uri]]

        # Query graph
        sparql_output = test_env.query_graph(sparql_query)

        self.assertEqual(expected_results, sparql_output)

    def test_genomic_position_model(self):
        """
        Test modelling of genomic positions
        Using test data set 2, and the function add_genotype_info_to_graph()
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
         protein_variant_source, variant_gene, bp_pos, genotype_cdna,
         cosmic_id, db_snp_id, genome_pos_start, genome_pos_end, ref_base,
         variant_base, primary_transcript_exons,
         primary_transcript_variant_sub_types, variant_type, chromosome,
         genome_build, build_version, build_date) = self.test_set_2[0]

        chromosome = ":hg19chr9"
        variant_position_id = self.cgd.make_id(
            'cgd-var-pos{0}{1}'.format(genotype_key, variant_gene))
        variant_position_label = '{0} genomic location'.format(variant_gene)
        region_id = ":_{0}Region".format(variant_position_id)
        start_id = ":_{0}-{1}".format(chromosome, genome_pos_start)
        end_id = ":_{0}-{1}".format(chromosome, genome_pos_end)

        position_uri = URIRef(cu.get_uri(variant_position_id))
        region_uri = URIRef(cu.get_uri(region_id))
        start_uri = URIRef(cu.get_uri(start_id))
        end_uri = URIRef(cu.get_uri(end_id))
        chromosome_uri = URIRef(cu.get_uri(chromosome))

        sparql_query = """
                       SELECT ?position ?region ?startPosition ?endPosition ?chromosome
                       WHERE {{
                           ?position a faldo:Position ;
                               rdfs:label "{0}" ;
                               faldo:location ?region .

                           ?region a faldo:Region ;
                               faldo:begin ?startPosition ;
                               faldo:end ?endPosition .

                           ?startPosition a faldo:BothStrandPosition ;
                               a faldo:Position ;
                               faldo:position {1} ;
                               faldo:reference ?chromosome .

                           ?endPosition a faldo:BothStrandPosition ;
                               a faldo:Position ;
                               faldo:position {2} ;
                               faldo:reference ?chromosome .
                       }}
                       """.format(variant_position_label, genome_pos_start,
                                  genome_pos_end,)

        # Expected Results
        expected_results = [[position_uri, region_uri, start_uri, end_uri, chromosome_uri]]

        # Query graph
        sparql_output = test_env.query_graph(sparql_query)

        self.assertEqual(expected_results, sparql_output)

if __name__ == '__main__':
    unittest.main()