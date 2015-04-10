from dipper import curie_map
from mckb.sources.CGD import CGD
from rdflib.namespace import URIRef
from dipper.utils.CurieUtil import CurieUtil
import unittest
import logging

logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)


class DiseaseDrugGenotypeTestCase(unittest.TestCase):
    """
    Test triples created from add_disease_drug_genotype_to_graph()
    Sample data for testing should resemble output from
    _get_disease_drug_genotype_relationship()
    """

    def setUp(self):

        self.curie_map = curie_map.get()
        cu = CurieUtil(self.curie_map)
        # Fake credentials as these tests do not require a database connection
        database = 'foo'
        user = 'bar'
        password = 'baz'

        self.cgd = CGD(database, user, password)
        test_data = ((387, 'MLH1 any mutation', 13, 'Adenocarcinoma',
                     None, 'Colon', 'detrimental effect', 1,
                     '5FU-based adjuvant therapy', 'late trials', '20498393'),)
        self.cgd.add_disease_drug_genotype_to_graph(test_data)

        (genotype_key, genotype_label, diagnoses_key, diagnoses,
             specific_diagnosis, organ, relationship,
             drug_key, drug, therapy_status, pubmed_id) = test_data[0]

        source_id = "PMID:{0}".format(pubmed_id)
        population_id = self.cgd.make_id('cgd{0}{1}'.format(genotype_key,
                                                            genotype_label))
        genotype_id = self.cgd.make_id('cgd-genotype{0}'.format(genotype_key))
        phenotype_id = self.cgd.make_id('cgd-phenotype{0}'.format(diagnoses_key))
        relationship_id = ("MONARCH:{0}".format(relationship)).replace(" ", "_")
        drug_id = self.cgd.make_id('cgd-drug{0}'.format(drug_key))
        genotype_annot = self.cgd.make_id("{0}{1}".format(population_id, genotype_label))
        phenotype_annot = self.cgd.make_id("{0}{1}".format(population_id, diagnoses))
        drug_annot = self.cgd.make_id("{0}{1}".format(population_id, drug))

        # Set up URIs
        self.source_uri = URIRef(cu.get_uri(source_id))
        self.population_uri = URIRef(cu.get_uri(population_id))
        self.genotype_uri = URIRef(cu.get_uri(genotype_id))
        self.phenotype_uri = URIRef(cu.get_uri(phenotype_id))
        self.relationship_uri = URIRef(cu.get_uri(relationship_id))
        self.drug_uri = URIRef(cu.get_uri(drug_id))
        self.genotype_annot_uri = URIRef(cu.get_uri(genotype_annot))
        self.phenotype_annot_uri = URIRef(cu.get_uri(phenotype_annot))
        self.drug_annot_uri = URIRef(cu.get_uri(drug_annot))

        self.genotype_label = genotype_label
        self.population_label = "Patient population diagnosed with {0} with" \
                                " genotype {1}".format(diagnoses, genotype_label)
        self.disease_label = diagnoses
        self.drug_label = drug

        return

    def tearDown(self):
        self.cgd.graph = None
        self.cgd = None
        return

    def test_classes_indiv_properties(self):
        """
        Given the above sample input, produce the following:
        A Monarch:PopulationID is an individual of type OBO:GENO_0000110
        A Monarch:PopulationID rdfs:label "Patient population diagnosed with Adenocarcinoma with genotype MLH1 any mutation"
        A Monarch:DiseaseID is an OWL Class
        A Monarch:Disease rdfs:label "Adenocarcinoma"
        A Monarch:DrugID is an OWL Class
        A Monarch:DrugID rdfs:label "5FU-based adjuvant therapy"
        A Monarch:RelationID is an object property

        Testing Note: Testing associations (assoc a Annotation) will
        occur elsewhere, but could also be grouped into this
        test if needed.
        """
        from dipper.utils.TestUtils import TestUtils

        # Make testutils object and load bindings
        test_env = TestUtils(self.cgd.graph)
        self.cgd.load_bindings()

        sparql_query = """
                       SELECT ?pop ?phenotype ?drug ?source
                       WHERE {{
                           ?pop a OBO:GENO_0000110 ;
                               rdfs:label "{0}" .
                           ?phenotype a owl:Class ;
                               rdfs:label "{1}" .
                           ?drug a owl:Class ;
                               rdfs:label "{2}" .
                           <{3}> a owl:ObjectProperty .
                           ?source a owl:NamedIndividual .
                       }}
                       """.format(self.population_label, self.disease_label,
                                  self.drug_label, self.relationship_uri)

        # Expected Results
        expected_results = [[self.population_uri, self.phenotype_uri,
                             self.drug_uri, self.source_uri]]
        # Query graph
        sparql_output = test_env.query_graph(sparql_query)

        self.assertEqual(expected_results, sparql_output)

    def test_population_triples(self):
        """
        Given the above sample input, produce the following:
        A population has_phenotype (OBO:RO_0002200) a disease
        A population has_genotype(OBO:RO_0002200) a genotype
        A population MONARCH:has_detrimental_effect to a drug/chemical
        """
        from dipper.utils.TestUtils import TestUtils

        # Make testutils object and load bindings
        test_env = TestUtils(self.cgd.graph)
        self.cgd.load_bindings()

        sparql_query = """
                       SELECT ?pop ?phenotype ?genotype ?drug
                       WHERE {{
                           ?pop OBO:RO_0002200 ?phenotype .
                           ?pop OBO:GENO_0000222 ?genotype .
                           ?pop <{0}> ?drug .
                       }}
                       """.format(self.relationship_uri)

        # Expected Results
        expected_results = [[self.population_uri, self.phenotype_uri,
                             self.genotype_uri, self.drug_uri]]
        # Query graph
        sparql_output = test_env.query_graph(sparql_query)

        self.assertEqual(expected_results, sparql_output)

    def test_associations(self):
        """
        Given the above sample input, produce the following:
        A Monarch:AssociationID dc:evidence Traceable Author Statement (ECO:0000033)
        A Monarch:AssociationID dc:source PMID:20498393
        A Monarch:AssociationID :hasSubject Monarch:PublicationID
        A Monarch:AssociationID :hasPredicate Genotype (GENO:0000222)
        A Monarch:AssociationID :hasObject Monarch:GenotypeID

        And two additional associations with the same evidence and source
        with hasSubject, hasPredicate, hasObject links documented in
        test_population_triples
        """
        from dipper.utils.TestUtils import TestUtils

        # Make testutils object and load bindings
        cu = CurieUtil(self.curie_map)
        test_env = TestUtils(self.cgd.graph)
        evidence = 'OBO:ECO_0000033'
        evidence_uri = URIRef(cu.get_uri(evidence))
        self.cgd.load_bindings()

        sparql_query = """
                       SELECT ?genotype_annot ?phenotype_annot
                              ?drug_annot ?evidence ?source
                              ?pop ?genotype ?phenotype ?drug
                       WHERE {{
                           ?genotype_annot a Annotation: ;
                               dc:evidence ?evidence ;
                               dc:source ?source ;
                               :hasObject ?genotype ;
                               :hasPredicate OBO:GENO_0000222 ;
                               :hasSubject ?pop .
                           ?phenotype_annot a Annotation: ;
                               dc:evidence ?evidence ;
                               dc:source ?source ;
                               :hasObject ?phenotype ;
                               :hasPredicate OBO:RO_0002200 ;
                               :hasSubject ?pop .
                           ?drug_annot a Annotation: ;
                               dc:evidence ?evidence ;
                               dc:source ?source ;
                               :hasObject ?drug ;
                               :hasPredicate <{0}> ;
                               :hasSubject ?pop .
                       }}
                       """.format(self.relationship_uri)

        # Expected Results
        expected_results = [[self.genotype_annot_uri, self.phenotype_annot_uri,
                             self.drug_annot_uri, evidence_uri,
                             self.source_uri, self.population_uri,
                             self.genotype_uri, self.phenotype_uri,
                             self.drug_uri]]
        # Query graph
        sparql_output = test_env.query_graph(sparql_query)

        self.assertEqual(expected_results, sparql_output)

if __name__ == '__main__':
    unittest.main()