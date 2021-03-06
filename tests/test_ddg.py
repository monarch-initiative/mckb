from dipper import curie_map
from mckb.sources.CGD import CGD
from rdflib.namespace import URIRef
from dipper.utils.CurieUtil import CurieUtil
import unittest
import logging

logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)


class DiseaseDrugVariantTestCase(unittest.TestCase):
    """
    Test triples created from add_disease_drug_variant_to_graph()

    Here we define a series of functional tests where we import RDFLib,
    create a test data set, run the data through a single
    function (add_disease_drug_variant_to_graph()), and
    test the in memory RDF graph with a sparql query
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
                     None, 'Colon', 'no response', 1,
                     '5FU-based adjuvant therapy', 'late trials', '20498393'),)
        self.cgd.add_disease_drug_variant_to_graph(test_data)

        (variant_key, variant_label, diagnoses_key, diagnoses,
         specific_diagnosis, organ, relationship,
         drug_key, drug, therapy_status, pubmed_id) = test_data[0]

        source_id = "PMID:{0}".format(pubmed_id)
        variant_id = self.cgd.make_cgd_id('variant{0}'.format(variant_key))
        disease_id = self.cgd.make_cgd_id('disease{0}{1}'.format(diagnoses_key,
                                                                 diagnoses))
        relationship_id = "RO:has_environment"
        disease_quality = ("CGD:{0}".format(relationship)).replace(" ", "_")
        has_quality_property = "BFO:0000159"
        drug_id = self.cgd.make_cgd_id('drug{0}'.format(drug_key))
        disease_instance_id = self.cgd.make_cgd_id('phenotype{0}{1}{2}'.format(
            diagnoses, variant_key, relationship))

        variant_disease_annot = self.cgd.make_cgd_id("assoc{0}{1}".format(variant_key, diagnoses))

        # Set up URIs
        self.source_uri = URIRef(cu.get_uri(source_id))
        self.variant_uri = URIRef(cu.get_uri(variant_id))
        self.disease_uri = URIRef(cu.get_uri(disease_id))
        self.disease_ind_uri = URIRef(cu.get_uri(disease_instance_id))
        self.relationship_uri = URIRef(cu.get_uri(relationship_id))
        self.drug_uri = URIRef(cu.get_uri(drug_id))
        self.vd_annot_uri = URIRef(cu.get_uri(variant_disease_annot))
        self.disease_quality_uri = URIRef(cu.get_uri(disease_quality))

        self.variant_label = variant_label
        self.disease_label = diagnoses
        self.disease_instance_label = "{0} with {1} to therapy".format(diagnoses, relationship)
        self.drug_label = drug

        return

    def tearDown(self):
        self.cgd.graph = None
        self.cgd = None
        return

    def test_classes_indiv_properties(self):
        """
        Given the above sample input, produce the following:
        A CGD:DiseaseID is an OWL Class
        A CGD:DiseaseID is a subclass of DOID:4
        A CGD:Disease rdfs:label "Adenocarcinoma"
        A CGD:DiseaseInstance is an individual of CGD:DiseaseID
        A CGD:DiseaseInstance rdfs:label "Adenocarcinoma with response {1} to therapy"
        A CGD:DrugID is an OWL Class
        A CGD:DrugID is a subclass of CHEBI:23888
        A CGD:DrugID rdfs:label "5FU-based adjuvant therapy"
        A CGD:RelationID is an object property
        PMID:12345 is a IAO:0000013 (journal article)
        """
        from dipper.utils.TestUtils import TestUtils

        # Make testutils object and load bindings
        test_env = TestUtils(self.cgd.graph)
        self.cgd.load_bindings()

        sparql_query = """
                       SELECT ?disease ?diseaseInd ?diseaseQual ?drug ?source
                       WHERE {{
                           ?disease a owl:Class ;
                               rdfs:subClassOf DOID:4 ;
                               rdfs:label "{0}" .
                           ?diseaseInd a ?disease ;
                               rdfs:label "{1}" ;
                               BFO:0000159 ?diseaseQual .
                           ?drug a owl:Class ;
                               rdfs:subClassOf CHEBI:23888 ;
                               rdfs:label "{2}" .
                           <{3}> a owl:ObjectProperty .
                           ?source a IAO:0000013 .
                       }}
                       """.format(self.disease_label, self.disease_instance_label,
                                  self.drug_label, self.relationship_uri)

        # Expected Results
        expected_results = [[self.disease_uri, self.disease_ind_uri,
                             self.disease_quality_uri, self.drug_uri,
                             self.source_uri]]
        # Query graph
        sparql_output = test_env.query_graph(sparql_query)

        self.assertEqual(expected_results, sparql_output)

    def test_associations(self):
        """
        Given the above sample input, produce the following:
        CGD:VariantID has_phenotype(RO:0002200) CGD:DiseaseInstance

        A CGD:AssociationID OBO:RO_0002558 Traceable Author Statement (ECO:0000033)
        A CGD:AssociationID dc:source PMID:20498393
        A CGD:AssociationID has_environment CGD:DrugID
        A CGD:AssociationID OBAN:association_has_subject CGD:VariantID
        A CGD:AssociationID OBAN:association_has_object_property has_phenotype
        A CGD:AssociationID OBAN:association_has_object CGD:DiseaseInstance
        """
        from dipper.utils.TestUtils import TestUtils

        # Make testutils object and load bindings
        cu = CurieUtil(self.curie_map)
        test_env = TestUtils(self.cgd.graph)
        self.cgd.load_bindings()
        evidence = 'OBO:ECO_0000033'
        evidence_uri = URIRef(cu.get_uri(evidence))

        sparql_query = """
                       SELECT ?diseaseInd ?variant ?drug ?vdannot ?source ?evidence
                       WHERE {{
                           ?variant OBO:RO_0002200 ?diseaseInd .

                           ?vdannot a OBAN:association ;
                               OBO:RO_0002558 ?evidence ;
                               dc:source ?source ;
                               <{0}> ?drug ;
                               OBAN:association_has_object ?diseaseInd ;
                               OBAN:association_has_object_property OBO:RO_0002200 ;
                               OBAN:association_has_subject ?variant .
                       }}
                       """.format(self.relationship_uri)

        # Expected Results
        expected_results = [[self.disease_ind_uri, self.variant_uri, self.drug_uri,
                             self.vd_annot_uri,
                             self.source_uri, evidence_uri]]
        # Query graph
        sparql_output = test_env.query_graph(sparql_query)

        self.assertEqual(expected_results, sparql_output)

if __name__ == '__main__':
    unittest.main()