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
            self._load_data_from_dump_file(file)
        else:
            logger.debug("Database contains tables, "
                         "skipping load from dump file")

        print(self._get_disease_drug_genotype_relationship())
        return

    def _get_disease_drug_genotype_relationship(self):
        """
        Query database to get disease-drug-genotype associations
        :return: tuple of query results
        """

        sql = """
            SELECT distinct
              tg.id as genotype_id,
              tg.comment as genotype_label,
              diagnoses.id as diagnoses_id,
              diagnoses.description as diagnoses,
              specific_diagnosis.id as specific_diagnosis_id,
              specific_diagnosis.description as specific_diagnosis,
              organs.id as organ_id,
              organs.description as organ,
              ta.id as relationship_id,
              ta.description as relationship,
              tc.id as drug_id,
              tc.description as drug,
              tgp.pub_med_id as pubmed_id
            FROM therapy_genotype tg
            JOIN diagnoses ON tg.diagnosis = diagnoses.id
            LEFT OUTER JOIN specific_diagnosis
            ON tg.specific_diagnosis = specific_diagnosis.id
            LEFT OUTER JOIN therapy_genotype_publication as tgp
            ON tg.id = tgp.therapy_genotype
            LEFT OUTER JOIN organs
            ON tg.organ = organs.id
            JOIN therapeutic_association as ta
            ON tg.therapeutic_association = ta.id
            JOIN therapeutic_context tc
            ON tg.therapeutic_context = tc.id;
        """
        self.cursor.execute(sql)
        results = self.cursor.fetchall()
        return results

    def _load_data_from_dump_file(self, file):
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

