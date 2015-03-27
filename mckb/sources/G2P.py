from mckb.sources import MySQLSource
from dipper.models.Dataset import Dataset
import logging

logger = logging.getLogger(__name__)


class G2P(MySQLSource):
    """
    Test data source for cancer knowledge base
    """

    def __init__(self):
        super().__init__(self, 'g2p')
        self.dataset = Dataset('g2p', 'G2P', 'http://ga4gh.org')

    def fetch(self, ):
        """
        Override Source.fetch()
        Fetches data from mysql dump
        Args:
            :param
        Returns:
            :return None
        """
        return

    def parse(self):
        """
        Override Source.parse()
        Args:
            :param
        Returns:
            :return None
        """
        return