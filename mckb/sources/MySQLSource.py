from dipper.sources.Source import Source
import logging

logger = logging.getLogger(__name__)


class MySQLSource(Source):
    """
    Abstract class for interacting with remote or local MySQL databases
    """

    def __init__(self):
        super().__init__()

    def load_database_from_dump_file(self, dump, host, database, user, password):
        pass

