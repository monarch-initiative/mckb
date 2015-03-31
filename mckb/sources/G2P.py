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
        self.rawdir = 'resources'

    def parse(self):
        """
        Override Source.parse()
        Args:
            :param
        Returns:
            :return None
        """
        logger.debug("Checking if database is empty")
        is_db_empty = self.check_if_db_is_empty()
        if is_db_empty:
            file = '/'.join((self.rawdir,
                                  self.static_files['test_data']['file']))
            logger.debug("Loading data into database from file {0}".format(file))
            self.load_data_from_dump_file(file)
        else:
            logger.debug("Database contains tables, "
                         "skipping load from dump file")
        return

    def load_data_from_dump_file(self, file):
        """
        Assumes dump file is gzipped
        Note: Here we load the database via a system command as
        cursor.execute(source file.sql) does not work and there is
        no robust way to parse and send the entire sql file via
        the execute function.

        If security might be an issue, it is probably best to
        pre-load the database before running to avoid having the
        password displayed from the command line
        :param file:
        :return: None
        """
        gz_file = gzip.open(file, 'rb')
        with tempfile.NamedTemporaryFile(mode='w+b') as f:
            f.write(gz_file.read())
            gz_file.close()
            os.system("mysql -h {0} -p{1} -u {2} -D {3} < {4}".format(self.host, self.password,
                      self.username, self.database, f.name))
        return

