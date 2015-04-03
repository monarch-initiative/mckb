from dipper import curie_map
from rdflib import Graph
from mckb.sources.CGD import CGD
import json
import os.path
import unittest
import logging

logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)


class DiseaseDrugGenotypeTestCase(unittest.TestCase):

    def setUp(self):

        self.graph = Graph()
        self.curie_map = curie_map.get()
        # There's probably a better way to do this
        if os.path.exists(os.path.join(os.path.dirname(__file__), 'conf/conf.json')):
            credentials = json.load(open(os.path.join(os.path.dirname(__file__),
                                         'conf/conf.json'), 'r'))
            host = credentials['dbauth']['g2p']['host']
            database = credentials['dbauth']['g2p']['database']
            user = credentials['dbauth']['g2p']['user']
            password = credentials['dbauth']['g2p']['password']

        self.cgd = CGD(database, user, password, host)
        self.cgd.graph = Graph()

        return

    def tearDown(self):
        self.cgd.graph = None
        self.cgd = None
        return

    def test_foo(self):
        self.assertEqual(1, 1)

if __name__ == '__main__':
    unittest.main()