from mckb.sources.CuratedSource import CuratedSource
import logging
import os
import csv

logger = logging.getLogger(__name__)


class CgdOntologyMap(CuratedSource):
    """
    Ontology mappings for CGD diseases and drugs, manually curated
    Entrez gene ID mappings from HGNC symbols, automatically/manually curated
    """

    static_files = {
        'disease_map': {
            'file': '../../resources/mappings/disease-map.tsv'
        },
        'gene_map': {
            'file': '../../resources/mappings/gene.tsv'
        },
        'drug_map': {
            'file': '../../resources/mappings/drug-map.tsv'
        }
    }

    def __init__(self, source_name):
        super().__init__(source_name)
        self.disease_map = {}
        self.gene_map = {}
        self.drug_map = {}

    def parse(self):
        self.disease_map = self._parse_mapping_file(
            self.static_files.disease_map.file)

        self.gene_map = self._parse_mapping_file(
            self.static_files.gene_map.file)

        self.drug_map = self._parse_mapping_file(
            self.static_files.drug_map.file)

    def _parse_mapping_file(self, file):
        """
        :param file: String, path to file containing label-id mappings in
                             the first two columns of each row
        :return: dict where keys are labels and values are ids
        """
        id_map = {}
        if os.path.exists(os.path.join(os.path.dirname(__file__), file)):
            with open(os.path.join(os.path.dirname(__file__), file)) as tsvfile:
                reader = csv.reader(tsvfile, delimiter="\t")
                for row in reader:
                    label = row[1]
                    id = row[2]
                    id_map[label] = id

        return id_map