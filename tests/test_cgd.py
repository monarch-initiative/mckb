from mckb.sources.CGD import CGD
import unittest
import logging

logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)


class MySQLSourceTestCase(unittest.TestCase):
    """
    Test triples created from add_disease_drug_variant_to_graph()
    Sample data for testing should resemble output from
    _get_disease_drug_variant_relationship()
    """

    def setUp(self):
        database = 'cgd_test'
        user = 'travis'
        self.cgd_test = CGD(database, user)
        self.cursor, self.connection = self.cgd_test._connect_to_database()
        return

    def tearDown(self):
        self.cgd_test._disconnect_from_database(self.cursor, self.connection)
        self.cgd_test = None
        return

    def test_connection(self):
        self.assertTrue(True)

if __name__ == '__main__':
    unittest.main()