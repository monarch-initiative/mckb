from mckb.sources.CGD import CGD
import unittest
import logging

logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)


class CGDTestCase(unittest.TestCase):
    """
    Test connection, loading, and querying of CGD snapshot
    """

    def setUp(self):
        database = 'cgd_test'
        user = 'travis'
        self.cgd_test = CGD(database, user)
        self.connection, self.cursor = self.cgd_test._connect_to_database()
        return

    def tearDown(self):
        self.cgd_test._disconnect_from_database(self.cursor, self.connection)
        self.cgd_test = None
        return

    def test_queries(self):
        """
        Just checking that these run without errors, probably needs
        do so some actual checking of things
        :return:
        """
        self.cgd_test.check_if_db_is_empty(self.cursor)
        self.cgd_test._get_disease_drug_variant_relationship(self.cursor)
        self.cgd_test._get_variant_protein_info(self.cursor)
        self.cgd_test._get_variant_cdna_info(self.cursor)
        self.cgd_test._get_fusion_copy_any_mutation_genotypes(self.cursor)
        self.cgd_test._get_genotypes_with_no_gene_protein_cdna_mapping(self.cursor)
        return

    def test_fetch(self):
        """
        Just checking that we can fetch sources without errors
        :return:
        """
        self.cgd_test.fetch()
        return

    def test_parse(self):
        """
        Just checking that we can parse sources without errors
        :return:
        """
        self.cgd_test.parse()
        return



if __name__ == '__main__':
    unittest.main()