from dipper.sources.Source import Source
import logging

logger = logging.getLogger(__name__)


class CuratedSource(Source):
    """
    Abstract class for interacting with curated data, for example, where some parts
    of the data are from one source and other parts are curated (manually or automatically)

    Right now this class serves as a placeholder.  The below functions could be implemented
    here or in a child class.  The idea here would be that these wrap a Dipper model class that
    generates a provenance chain in the graph
    """

    def __init__(self, source_name, reference_dataset):
        super().__init__(source_name)
        self.reference_dataset = reference_dataset

    def addCurator(self):
        pass

    def addCurationDateTime(self):
        pass

    def addReferenceSource(self):
        pass
