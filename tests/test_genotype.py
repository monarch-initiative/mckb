from dipper import curie_map
from mckb.sources.CGD import CGD
from mckb.sources.CGDOntologyMap import CGDOntologyMap
from rdflib.namespace import URIRef
from dipper.utils.CurieUtil import CurieUtil
import unittest
import logging
import datetime

logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)


class DiseaseDrugVariantTestCase(unittest.TestCase):
    """
    Test triples created from variant modelling functions

    Here we define a series of functional tests where we import RDFLib,
    create a test data set, run the data through a single function, and
    test the in memory RDF graph with a sparql query
    """

    def setUp(self):

        self.curie_map = curie_map.get()

        # Fake credentials as these tests do not require a database connection
        database = 'foo'
        user = 'bar'
        password = 'baz'

        self.cgd = CGD(database, user, password)
        ontology_map = CGDOntologyMap('cgd-ontology-mappings')
        ontology_map.parse()

        self.cgd.gene_map = ontology_map.gene_map

        # Sample output from _get_variant_protein_info() where variant
        # is a missense mutation
        self.test_set_1 = ((2, 'CSF3R Q741X  missense mutation', 'p.Q741X ',
                            None, 'CCDS413.1', 'Primary', None,
                            'gain-of-function', None, 'CSF3R', None),)

        # Sample output from _get_variant_cdna_info()
        self.test_set_2 = ((19, 'ABL1 T315I missense mutation', 'p.T315I',
                            315, 'CCDS35166.1', 'Secondary',
                            'nonsynonymous - missense', 'gain-of-function',
                            None, 'ABL1', None, 'ABL1', 944, 'c.944C>T',
                            'COSM12560', 'rs121913459', 133748283, 133748283,
                            'C', 'T', 'Ex6', 'nonsynonymous - missense',
                            'Substitution', 'chr9', 'hg19', 'hg19',
                            datetime.datetime(2009, 2, 1, 0, 0)),)

        self.cgd.transcript_xrefs = {
            'RefSeq':  {'CCDS35166.1': 'NP_005148.2',
                        'CCDS413.1': 'NP_000751.1'},
            'UniProt': {'CCDS35166.1': 'P00519-1',
                        'CCDS413.1': 'Q99062-1'}
        }

        return

    def tearDown(self):
        self.cgd.graph = None
        self.cgd = None
        return

    def test_missense_variant_protein_model(self):
        """
        Test missense variant with only protein information
        Using test data set 1, and the function add_variant_info_to_graph()
        We want to test the following triples:

        CGD:VariantID is an instance of OBO:SO_0001059
        CGD:VariantID is an instance of OBO:SO_0001583
        CGD:VariantID has the label "CSF3R Q741X  missense mutation"
        CGD:VariantID is_sequence_variant_instance_of (OBO:GENO_0000408) NCBIGene:1441
        CGD:VariantID has location (faldo:location) CGD:RegionID
        CGD:VariantID OBO:GENO_reference_amino_acid "Q"
        CGD:VariantID OBO:GENO_results_in_amino_acid_change "X"
        CGD:VariantID RO:0002205 CCDS:413.1

        CCDS:413.1 is an instance of OBO:GENO_primary
        CCDS:413.1 has the label "CCDS413.1"
        """
        from dipper.utils.TestUtils import TestUtils

        self.cgd.add_variant_info_to_graph(self.test_set_1)

        # Make testutils object and load bindings
        test_env = TestUtils(self.cgd.graph)
        cu = CurieUtil(self.curie_map)
        self.cgd.load_bindings()

        (variant_key, variant_label, amino_acid_variant, amino_acid_position,
         transcript_id, transcript_priority, protein_variant_type,
         functional_impact, stop_gain_loss, transcript_gene,
         protein_variant_source) = self.test_set_1[0][0:11]

        gene_id = self.cgd.gene_map[transcript_gene]
        ref_amino_acid = "Q"
        altered_amino_acid = "X"
        position = 741
        uniprot_curie = "UniProtKB:Q99062#Q99062-1"

        variant_id = self.cgd.make_cgd_id('variant{0}'.format(variant_key))
        transcript = "CCDS:413.1"
        region_id = ":_{0}{1}{2}Region".format(position, position, uniprot_curie)
        variant_uri = URIRef(cu.get_uri(variant_id))
        transcript_uri = URIRef(cu.get_uri(transcript))
        gene_uri = URIRef(cu.get_uri(gene_id))
        region_uri = URIRef(cu.get_uri(region_id))

        sparql_query = """
                       SELECT ?variant ?gene ?region ?transcript
                       WHERE {{
                           ?variant a OBO:SO_0001059;
                               a OBO:SO_0001583 ;
                               rdfs:label "{0}" ;
                               OBO:GENO_0000408 ?gene ;
                               faldo:location ?region ;
                               OBO:GENO_reference_amino_acid "{1}" ;
                               OBO:GENO_results_in_amino_acid_change "{2}" ;
                               RO:0002205 ?transcript .

                           ?transcript a OBO:SO_0000233 ;
                               rdfs:label "{3}" .
                       }}
                       """.format(variant_label, ref_amino_acid,
                                  altered_amino_acid, transcript_id)

        # Expected Results
        expected_results = [[variant_uri, gene_uri, region_uri, transcript_uri]]
        # Query graph
        sparql_output = test_env.query_graph(sparql_query)

        self.assertEqual(expected_results, sparql_output)

    def test_missense_variant_cdna_model(self):
        """
        Test missense variant with cdna information
        Using test data set 2, and the function add_variant_info_to_graph()
        We want to test the following triples:

        CGD:VariantID is an instance of OBO:SO_0001059
        CGD:VariantID is an instance of OBO:SO_0001583
        CGD:VariantID has the label "ABL1 T315I missense mutation"
        CGD:VariantID is_sequence_variant_instance_of (OBO:GENO_0000408) NCBIGene:25
        CGD:VariantID has location (faldo:location) AminoAcidRegionID
        CGD:VariantID has location (faldo:location) CDNARegionID
        CGD:VariantID has location (faldo:location) ChromosomalRegionID
        CGD:VariantID OBO:GENO_reference_amino_acid "T"
        CGD:VariantID OBO:GENO_results_in_amino_acid_change "I"
        CGD:VariantID owl:sameAs dbSNP:rs121913459
        CGD:VariantID owl:sameAs COSMIC:12560
        CGD:VariantID RO:0002205 (transcribed_to) CCDS:35166.1

        CCDS:35166.1 is an instance of OBO:SO_0000233
        CCDS:35166.1 has the label "CCDS35166.1"
        CCDS:35166.1 OBO:RO_0002513 (translates_to) UniProtKB:P00519#P00519-1
        CCDS:35166.1 OBO:RO_0002513 (translates_to) NCBIProtein:NP_005148.2

        UniProtKB:P00519#P00519-1 owl:sameAs NCBIProtein:NP_005148.2

        UniProtKB:P00519#P00519-1 is an instance of OBO:SO_0000104 (polypeptide)
        UniProtKB:P00519#P00519-1 has the label "P00519#P00519-1"

        NCBIProtein:NP_005148.2 is an instance of OBO:SO_0000104 (polypeptide)
        NCBIProtein:NP_005148.2 has the label "NP_005148.2"
        """
        from dipper.utils.TestUtils import TestUtils

        self.cgd.add_variant_info_to_graph(self.test_set_2)

        # Make testutils object and load bindings
        test_env = TestUtils(self.cgd.graph)
        cu = CurieUtil(self.curie_map)
        self.cgd.load_bindings()

        (variant_key, variant_label, amino_acid_variant, amino_acid_position,
         transcript_id, transcript_priority, protein_variant_type,
         functional_impact, stop_gain_loss, transcript_gene,
         protein_variant_source, variant_gene, bp_pos, variant_cdna,
         cosmic_id, db_snp_id, genome_pos_start, genome_pos_end, ref_base,
         variant_base, primary_transcript_exons,
         primary_transcript_variant_sub_types, variant_type, chromosome,
         genome_build, build_version, build_date) = self.test_set_2[0]

        gene_id = self.cgd.gene_map[transcript_gene]
        ref_amino_acid = "T"
        altered_amino_acid = "I"
        db_snp_curie = "dbSNP:121913459"
        cosmic_curie = "COSMIC:12560"
        uniprot_curie = "UniProtKB:P00519#P00519-1"
        uniprot_id = "P00519#P00519-1"
        refseq_curie = "NCBIProtein:NP_005148.2"
        transcript_curie = "CCDS:35166.1"
        ccds_id = "35166.1"
        position = 315
        chromosome_curie = "hg19chr9"

        variant_id = self.cgd.make_cgd_id('variant{0}'.format(variant_key))
        aa_region_id = ":_{0}{1}{2}Region".format(position, position, uniprot_curie)
        cdna_region_id = ":_{0}Region".format(transcript_curie)
        chr_region_id = ":_{0}{1}Region-{2}-{3}".format(genome_build,
                                                        chromosome,
                                                        genome_pos_start,
                                                        genome_pos_end)
        aa_coord_id = ":_{0}-{1}".format(uniprot_id, position)
        cdna_coord_id = ":_{0}-{1}".format(ccds_id, bp_pos)
        # chr_coord_id = "CHR:{0}-{1}".format(chromosome_curie, genome_pos_start)
        chr_coord_id = ":_{0}-{1}".format(chromosome_curie, genome_pos_start)

        variant_uri = URIRef(cu.get_uri(variant_id))
        transcript_uri = URIRef(cu.get_uri(transcript_curie))
        gene_uri = URIRef(cu.get_uri(gene_id))
        db_snp_uri = URIRef(cu.get_uri(db_snp_curie))
        cosmic_uri = URIRef(cu.get_uri(cosmic_curie))
        uniprot_uri = URIRef(cu.get_uri(uniprot_curie))
        refseq_uri = URIRef(cu.get_uri(refseq_curie))
        aa_region_uri = URIRef(cu.get_uri(aa_region_id))
        cdna_region_uri = URIRef(cu.get_uri(cdna_region_id))
        chr_region_uri = URIRef(cu.get_uri(chr_region_id))
        aa_coord_uri = URIRef(cu.get_uri(aa_coord_id))
        cdna_coord_uri = URIRef(cu.get_uri(cdna_coord_id))
        chr_coord_uri = URIRef(cu.get_uri(chr_coord_id))


        sparql_query = """
                       SELECT ?cosmic ?gene ?aaRegion ?cdnaRegion ?chrRegion
                              ?dbSNP ?transcript ?uniprot ?refseq
                              ?aaCoord ?cdnaCoord ?chrCoord
                       WHERE {{
                           ?cosmic a OBO:SO_0001059;
                               a OBO:SO_0001583 ;
                               OBO:GENO_0000408 ?gene ;
                               faldo:location ?aaRegion ;
                               faldo:location ?cdnaRegion ;
                               faldo:location ?chrRegion ;
                               OBO:GENO_reference_amino_acid "{0}" ;
                               OBO:GENO_reference_nucleotide "{1}" ;
                               OBO:GENO_altered_nucleotide "{2}" ;
                               OBO:GENO_results_in_amino_acid_change "{3}" ;
                               owl:sameAs ?dbSNP ;
                               RO:0002205 ?transcript .

                           ?cosmic owl:sameAs ?dbSNP .

                           ?transcript a OBO:SO_0000233 ;
                               rdfs:label "{4}" ;
                               OBO:RO_0002513 ?uniprot ;
                               OBO:RO_0002513 ?refseq .

                           ?uniprot a OBO:SO_0000104 ;
                               rdfs:label "P00519-1" .

                           ?refseq a OBO:SO_0000104 ;
                               rdfs:label "NP_005148.2" .

                           ?refseq owl:sameAs ?uniprot .

                           ?aaRegion faldo:begin ?aaCoord .
                           ?cdnaRegion faldo:begin ?cdnaCoord .
                           ?chrRegion faldo:begin ?chrCoord .

                           ?aaCoord faldo:position {5} .
                           ?cdnaCoord faldo:position {6} .
                           ?chrCoord faldo:position {7} .

                           ?dbSNP rdfs:label "{8}" .
                       }}
                       """.format(ref_amino_acid, ref_base,
                                  variant_base, altered_amino_acid,
                                  transcript_id, position, bp_pos,
                                  genome_pos_start, db_snp_id)

        # Expected Results
        expected_results = [[cosmic_uri, gene_uri, aa_region_uri,
                             cdna_region_uri, chr_region_uri, db_snp_uri,
                             transcript_uri, uniprot_uri,
                             refseq_uri, aa_coord_uri, cdna_coord_uri,
                             chr_coord_uri]]
        # Query graph
        sparql_output = test_env.query_graph(sparql_query)

        self.assertEqual(expected_results, sparql_output)

    def test_amino_acid_position_region_model(self):
        """
        Test modelling of amino acid positions
        Using test data set 1, and the function add_variant_info_to_graph()
        We want to test the following triples:

        CGD:RegionID is an instance of faldo:Region
        CGD:RegionID faldo:begin BothStrandPositionID
        CGD:RegionID faldo:end BothStrandPositionID

        CGD:BothStrandPositionID is an instance of faldo:BothStrandPosition
        CGD:BothStrandPositionID is an instance of faldo:Position
        CGD:BothStrandPositionID faldo:position 741
        CGD:BothStrandPositionID faldo:reference UniProtID
        """
        from dipper.utils.TestUtils import TestUtils
        self.cgd.add_variant_info_to_graph(self.test_set_1)

        # Make testutils object and load bindings
        test_env = TestUtils(self.cgd.graph)
        cu = CurieUtil(self.curie_map)
        self.cgd.load_bindings()

        (variant_key, variant_label, amino_acid_variant, amino_acid_position,
         transcript_id, transcript_priority, protein_variant_type,
         functional_impact, stop_gain_loss, transcript_gene,
         protein_variant_source) = self.test_set_1[0][0:11]

        position = 741
        variant_id = self.cgd.make_cgd_id('variant{0}'.format(variant_key))

        uniprot_curie = "UniProtKB:Q99062#Q99062-1"
        uniprot_id = "Q99062#Q99062-1"
        region_id = ":_{0}{1}{2}Region".format(position, position, uniprot_curie)
        both_strand_id = ":_{0}-{1}".format(uniprot_id, position)

        region_uri = URIRef(cu.get_uri(region_id))
        both_strand_uri = URIRef(cu.get_uri(both_strand_id))
        uniprot_uri = URIRef(cu.get_uri(uniprot_curie))

        sparql_query = """
                       SELECT ?region ?bsPosition ?protein
                       WHERE {{
                           ?region a faldo:Region ;
                               faldo:begin ?bsPosition ;
                               faldo:end ?bsPosition .

                           ?bsPosition a faldo:Position ;
                               faldo:position {0} ;
                               faldo:reference ?protein .
                       }}
                       """.format(position)

        # Expected Results
        expected_results = [[region_uri, both_strand_uri, uniprot_uri]]

        # Query graph
        sparql_output = test_env.query_graph(sparql_query)

        self.assertEqual(expected_results, sparql_output)

    def test_variant_position_region_model(self):
        """
        Test modelling of variant positions on a transcript
        Using test data set 2, and the function add_variant_info_to_graph()
        We want to test the following triples:

        CGD:RegionID is an instance of faldo:Region
        CGD:RegionID faldo:begin BothStrandPositionID
        CGD:RegionID faldo:end BothStrandPositionID

        CGD:BothStrandPositionID is an instance of faldo:BothStrandPosition
        CGD:BothStrandPositionID is an instance of faldo:Position
        CGD:BothStrandPositionID faldo:position 944
        CGD:BothStrandPositionID faldo:reference CGD:TranscriptID
        """
        from dipper.utils.TestUtils import TestUtils
        self.cgd.add_variant_info_to_graph(self.test_set_2)

        # Make testutils object and load bindings
        test_env = TestUtils(self.cgd.graph)
        cu = CurieUtil(self.curie_map)
        self.cgd.load_bindings()

        (variant_key, variant_label, amino_acid_variant, amino_acid_position,
         transcript_id, transcript_priority, protein_variant_type,
         functional_impact, stop_gain_loss, transcript_gene,
         protein_variant_source, variant_gene, bp_pos, variant_cdna,
         cosmic_id, db_snp_id, genome_pos_start, genome_pos_end, ref_base,
         variant_base, primary_transcript_exons,
         primary_transcript_variant_sub_types, variant_type, chromosome,
         genome_build, build_version, build_date) = self.test_set_2[0]

        transcript_curie = self.cgd._make_transcript_curie(transcript_id)
        ccds_id = "35166.1"
        variant_id = self.cgd.make_cgd_id('variant{0}'.format(variant_key))

        region_id = ":_{0}Region".format(transcript_curie)
        both_strand_id = ":_{0}-{1}".format(ccds_id, bp_pos)

        region_uri = URIRef(cu.get_uri(region_id))
        both_strand_uri = URIRef(cu.get_uri(both_strand_id))
        ccds_uri = URIRef(cu.get_uri(transcript_curie))

        sparql_query = """
                       SELECT ?region ?bsPosition ?transcript
                       WHERE {{
                           ?region a faldo:Region ;
                               faldo:begin ?bsPosition ;
                               faldo:end ?bsPosition .

                           ?bsPosition a faldo:Position ;
                               faldo:position {0} ;
                               faldo:reference ?transcript .
                       }}
                       """.format(bp_pos)

        # Expected Results
        expected_results = [[region_uri, both_strand_uri, ccds_uri]]

        # Query graph
        sparql_output = test_env.query_graph(sparql_query)

        self.assertEqual(expected_results, sparql_output)

    def test_genome_build_chromosome_model(self):
        """
        Test modelling of genome, builds, and chromosomes
        Using test data set 2, and the function add_variant_info_to_graph()
        """
        from dipper.utils.TestUtils import TestUtils
        self.cgd.add_variant_info_to_graph(self.test_set_2)

        # Make testutils object and load bindings
        test_env = TestUtils(self.cgd.graph)
        cu = CurieUtil(self.curie_map)
        self.cgd.load_bindings()

        genome = ":9606genome"
        genome_label = "Human genome"
        chromosome = "CHR:9606chr9"
        chromosome_label = "chr9 (Human)"
        build_curie = "UCSC:hg19"
        build_label = "hg19"
        chrom_on_build = ":MONARCH_hg19chr9"
        chrom_build_label = "chr9 (hg19)"

        genome_uri = URIRef(cu.get_uri(genome))
        chromosome_uri = URIRef(cu.get_uri(chromosome))
        build_uri = URIRef(cu.get_uri(build_curie))
        chrom_on_build_uri = URIRef(cu.get_uri(chrom_on_build))
        '''
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
        '''
        sparql_query = """
                       SELECT ?genome ?chromosome ?build ?chromOnBuild
                       WHERE {{
                           ?genome a owl:Class ;
                               rdfs:label "{0}" ;
                               rdfs:subClassOf OBO:SO_0001026 .

                           ?chromosome a owl:Class ;
                               rdfs:label "{1}" ;
                               rdfs:subClassOf OBO:SO_0000340 .

                           ?build a OBO:SO_0001505 ;
                               a ?genome ;
                               rdfs:label "{2}" ;
                               OBO:RO_0002162 OBO:NCBITaxon_9606 ;
                               OBO:RO_0002351 ?chromOnBuild .

                           ?chromOnBuild a ?chromosome ;
                               a OBO:SO_0000340 ;
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

    def test_chromosome_position_model(self):
        """
        Test modelling of genomic positions
        Using test data set 2, and the function add_variant_info_to_graph()
        """
        from dipper.utils.TestUtils import TestUtils
        self.cgd.add_variant_info_to_graph(self.test_set_2)

        # Make testutils object and load bindings
        test_env = TestUtils(self.cgd.graph)
        cu = CurieUtil(self.curie_map)
        self.cgd.load_bindings()

        (variant_key, variant_label, amino_acid_variant, amino_acid_position,
         transcript_id, transcript_priority, protein_variant_type,
         functional_impact, stop_gain_loss, transcript_gene,
         protein_variant_source, variant_gene, bp_pos, variant_cdna,
         cosmic_id, db_snp_id, genome_pos_start, genome_pos_end, ref_base,
         variant_base, primary_transcript_exons,
         primary_transcript_variant_sub_types, variant_type, chromosome,
         genome_build, build_version, build_date) = self.test_set_2[0]

        variant_id = self.cgd.make_cgd_id('variant{0}'.format(variant_key))

        chromosome_curie = ":MONARCH_hg19chr9"
        region_id = ":_{0}{1}Region-{2}-{3}".format(genome_build,
                                                    chromosome,
                                                    genome_pos_start,
                                                    genome_pos_end)
        start_id = ":_hg19chr9-{0}".format(genome_pos_start)
        end_id = ":_hg19chr9-{0}".format(genome_pos_end)

        region_uri = URIRef(cu.get_uri(region_id))
        start_uri = URIRef(cu.get_uri(start_id))
        end_uri = URIRef(cu.get_uri(end_id))
        chromosome_uri = URIRef(cu.get_uri(chromosome_curie))

        sparql_query = """
                       SELECT ?region ?startPosition ?endPosition ?chromosome
                       WHERE {{
                           ?region a faldo:Region ;
                               faldo:begin ?startPosition ;
                               faldo:end ?endPosition .

                           ?startPosition a faldo:Position ;
                               faldo:position {0} ;
                               faldo:reference ?chromosome .

                           ?endPosition a faldo:Position ;
                               faldo:position {1} ;
                               faldo:reference ?chromosome .
                       }}
                       """.format(genome_pos_start, genome_pos_end,)

        # Expected Results
        expected_results = [[region_uri, start_uri, end_uri, chromosome_uri]]

        # Query graph
        sparql_output = test_env.query_graph(sparql_query)

        self.assertEqual(expected_results, sparql_output)

if __name__ == '__main__':
    unittest.main()