from dipper.sources.Source import Source
import pymysql
import logging

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

        (self.connection, self.cursor) = self._connect_to_database()

    def _connect_to_database(self):
        logger.debug("Connecting to database %s on %s", self.database, self.host)

        connection = pymysql.connect(host=self.host, user=self.username,
                                     passwd=self.password, db=self.database)
        cursor = connection.cursor()
        logger.debug("Connected to %s:%s", self.host, self.database)
        return connection, cursor

    def disconnect_from_database(self):
        logger.debug("Disconnecting from to database %s", self.database)
        self.cursor.close()
        self.connection.close()
        return

    def check_if_db_is_empty(self):
        is_db_empty = True
        query = "SELECT count(*) FROM information_schema.tables " \
                "WHERE table_type = 'BASE TABLE' " \
                "AND table_schema = '{0}'".format(self.database)

        self.cursor.execute(query)
        results = self.cursor.fetchone()
        if results[0] > 0:
            is_db_empty = False

        return is_db_empty
