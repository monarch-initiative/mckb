from mckb.sources.MySQLSource import MySQLSource
from dipper.models.Dataset import Dataset
import tempfile
import gzip
import logging
import os


logger = logging.getLogger(__name__)


class G2P(MySQLSource):
    """
    Test data source for cancer knowledge base
    """
    static_files = {
        'test_data': {'file': 'g2p.sql.gz'}
    }

    def __init__(self, database, username, password, host=None):
        super().__init__('g2p', database, username, password, host)
        self.dataset = Dataset('g2p', 'G2P', 'http://ga4gh.org')
        self.rawdir = 'resources/'

    def parse(self):
        """
        Override Source.parse()
        Args:
            :param
        Returns:
            :return None
        """
        file = '/'.join((self.rawdir,
                              self.static_files['test_data']['file']))
        self.load_data_from_dump_file(file)
        return

    def load_data_from_dump_file(self, file):
        """
        Assumes dump file is gzipped
        :param file:
        :return: None
        """
        logger.debug("Loading data into database from file {0}".format(file))
        gz_file = gzip.open(file, 'rb')
        with tempfile.NamedTemporaryFile(mode='w+b') as f:
            f.write(gz_file.read())
            gz_file.close()
            os.system("mysql -h {0} -p{1} -u {2} -D {3} < {4}".format(self.host, self.password,
                      self.username, self.database, f.name))
        return

