from dipper.sources.Source import Source
import pymysql
import logging
import os

logger = logging.getLogger(__name__)


class MySQLSource(Source):
    """
    Abstract class for interacting with remote or local MySQL databases
    """

    def __init__(self, source_name, database, username, password=None, host=None):
        super().__init__(source_name)
        self.database = database
        self.username = username
        self.password = password
        self.host = host

        if self.host is None:
            self.host = "localhost"

    def _connect_to_database(self):
        logger.info("Connecting to database %s on %s", self.database, self.host)

        connection = pymysql.connect(host=self.host, user=self.username,
                                     passwd=self.password, db=self.database)
        cursor = connection.cursor()
        logger.info("Connected to %s:%s", self.host, self.database)
        return connection, cursor

    def _disconnect_from_database(self, cursor, connection):
        logger.info("Disconnecting from to database %s", self.database)
        cursor.close()
        connection.close()
        return

    def check_if_db_is_empty(self, cursor):
        is_db_empty = True
        query = "SELECT count(*) FROM information_schema.tables " \
                "WHERE table_type = 'BASE TABLE' " \
                "AND table_schema = '{0}'".format(self.database)

        cursor.execute(query)
        results = cursor.fetchone()
        if results[0] > 0:
            is_db_empty = False

        return is_db_empty

    @staticmethod
    def execute_query(cursor, file):
        """
        Query database from query stored in external file
        :param: PyMySQL Cursor object generated
                from self._connect_to_database()
        :param: file path of file relative to this class file
        :return: results of query as a String
        """

        try:
            query_file = open(os.path.join(os.path.dirname(__file__), file), 'rb')
            query = query_file.read()
        except IOError:
            print("Unable to open sql file")

        cursor.execute(query)
        results = cursor.fetchall()
        return results
