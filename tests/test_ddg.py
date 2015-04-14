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
        genotype_id = self.cgd.make_id('cgd-genotype{0}'.format(genotype_key))
        disease_id = self.cgd.make_id('cgd-disease{0}{1}'.format(diagnoses_key,
                                                                 diagnoses))
        relationship_id = ("MONARCH:{0}".format(relationship)).replace(" ", "_")
        drug_id = self.cgd.make_id('cgd-drug{0}'.format(drug_key))
        disease_instance_id = self.cgd.make_id('cgd-disease{0}{1}'.format(
            diagnoses, genotype_key))
        disease_genotype_annot = self.cgd.make_id("assoc{0}{1}".format(
            disease_instance_id, genotype_key))

        # Set up URIs
        self.source_uri = URIRef(cu.get_uri(source_id))
        self.genotype_uri = URIRef(cu.get_uri(genotype_id))
        self.disease_uri = URIRef(cu.get_uri(disease_id))
        self.disease_ind_uri = URIRef(cu.get_uri(disease_instance_id))
        self.relationship_uri = URIRef(cu.get_uri(relationship_id))
        self.drug_uri = URIRef(cu.get_uri(drug_id))
        self.dg_annot_uri = URIRef(cu.get_uri(disease_genotype_annot))

        self.genotype_label = genotype_label
        self.population_label = "Patient population diagnosed with {0} with" \
                                " genotype {1}".format(diagnoses, genotype_label)
        self.disease_label = diagnoses
        self.disease_instance_label = "{0} caused by variant {1}".format(
            diagnoses, genotype_label)
        self.drug_label = drug

        return

    def tearDown(self):
        self.cgd.graph = None
        self.cgd = None
        return

    def test_classes_indiv_properties(self):
        """
        Given the above sample input, produce the following:
        A Monarch:DiseaseID is an OWL Class
        A Monarch:Disease rdfs:label "Adenocarcinoma"
        A Monarch:DiseaseInstance is an individual of Monarch:DiseaseID
        A Monarch:DiseaseInstance rdfs:label "Adenocarcinoma caused by variant MLH1 any mutation"
        A Monarch:DrugID is an OWL Class
        A Monarch:DrugID rdfs:label "5FU-based adjuvant therapy"
        A Monarch:RelationID is an object property
        PMID:12345 is a named individual
        """
        from dipper.utils.TestUtils import TestUtils

        # Make testutils object and load bindings
        test_env = TestUtils(self.cgd.graph)
        self.cgd.load_bindings()

        sparql_query = """
                       SELECT ?disease ?diseaseInd ?drug ?source
                       WHERE {{
                           ?disease a owl:Class ;
                               rdfs:label "{0}" .
                           ?diseaseInd a ?disease ;
                               rdfs:label "{1}" .
                           ?drug a owl:Class ;
                               rdfs:label "{2}" .
                           <{3}> a owl:ObjectProperty .
                           ?source a owl:NamedIndividual .
                       }}
                       """.format(self.disease_label, self.disease_instance_label,
                                  self.drug_label, self.relationship_uri)

        # Expected Results
        expected_results = [[self.disease_uri, self.disease_ind_uri,
                             self.drug_uri, self.source_uri]]
        # Query graph
        sparql_output = test_env.query_graph(sparql_query)

        self.assertEqual(expected_results, sparql_output)

    def test_associations(self):
        """
        Given the above sample input, produce the following:
        A Monarch:DiseaseInstance RO:caused_by Monarch:GenotypeID

        A Monarch:DrugID has_relationship_to Monarch:AssociationID

        A Monarch:AssociationID dc:evidence Traceable Author Statement (ECO:0000033)
        A Monarch:AssociationID dc:source PMID:20498393
        A Monarch:AssociationID :hasSubject A Monarch:DiseaseInstance
        A Monarch:AssociationID :hasPredicate RO:caused_by
        A Monarch:AssociationID :hasObject Monarch:GenotypeID
        """
        from dipper.utils.TestUtils import TestUtils

        # Make testutils object and load bindings
        cu = CurieUtil(self.curie_map)
        test_env = TestUtils(self.cgd.graph)
        self.cgd.load_bindings()
        evidence = 'OBO:ECO_0000033'
        evidence_uri = URIRef(cu.get_uri(evidence))

        sparql_query = """
                       SELECT ?diseaseInd ?genotype ?dgannot ?drug ?source ?evidence
                       WHERE {{
                           ?diseaseInd OBO:RO_caused_by ?genotype .

                           ?drug <{0}> ?dgannot .

                           ?dgannot a Annotation: ;
                               dc:evidence ?evidence ;
                               dc:source ?source ;
                               :hasObject ?genotype ;
                               :hasPredicate OBO:RO_caused_by ;
                               :hasSubject ?diseaseInd .
                       }}
                       """.format(self.relationship_uri)

        # Expected Results
        expected_results = [[self.disease_ind_uri, self.genotype_uri,
                             self.dg_annot_uri, self.drug_uri,
                             self.source_uri, evidence_uri]]
        # Query graph
        sparql_output = test_env.query_graph(sparql_query)

        self.assertEqual(expected_results, sparql_output)

if __name__ == '__main__':
    unittest.main()