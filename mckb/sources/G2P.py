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
            JOIN diagnoses
            ON tg.diagnosis = diagnoses.id

            JOIN therapeutic_association as ta
            ON tg.therapeutic_association = ta.id

            JOIN therapeutic_context tc
            ON tg.therapeutic_context = tc.id

            LEFT OUTER JOIN specific_diagnosis
            ON tg.specific_diagnosis = specific_diagnosis.id

            LEFT OUTER JOIN therapy_genotype_publication as tgp
            ON tg.id = tgp.therapy_genotype

            LEFT OUTER JOIN organs
            ON tg.organ = organs.id;
        """
        self.cursor.execute(sql)
        results = self.cursor.fetchall()
        return results

    def _get_genotype_protein_info(self):
        """
        STATUS: incomplete
        Query database to therapy genotypes that have been mapped to
        a protein variant (rather than a cdna variant)
        :return: tuple of query results
        """

        sql = """
            SELECT distinct
              tg.id as therapy_genotype_id,
              tg.comment as genotype_label,
              pv.genotype_amino_acid_onel as aa_var,
              transcript.description as transcript_id,
              transcript_priority.description as transcript_priority,
              protein_variant_type.description as protein_variant_type,
              functional_impact.description as functional_impact,
              stop_gain_loss.description as stop_gain_loss,
              trg.description as transcript_gene,
              pv.pub_med_ids as pubmed_ids

            FROM therapy_genotype tg
            JOIN therapy_variant tv
            ON tg.id = tv.therapy_genotype

            JOIN protein_variant pv
            ON tv.protein_variant = pv.id

            LEFT OUTER JOIN transcript
            ON pv.transcript = transcript.id

            LEFT OUTER JOIN transcript_priority
            ON transcript.transcript_priority = transcript_priority.id

            LEFT OUTER JOIN protein_variant_type
            ON pv.protein_variant_type = protein_variant_type.id

            LEFT OUTER JOIN functional_impact
            ON pv.functional_impact = functional_impact.id

            LEFT OUTER JOIN stop_gain_loss
            ON pv.stop_gain_loss = stop_gain_loss.id

            LEFT OUTER JOIN gene trg
            ON transcript.gene = trg.id;
        """
        self.cursor.execute(sql)
        results = self.cursor.fetchall()
        return results

    def _get_genotype_cdna_info(self):
        """
        STATUS: incomplete
        Query database to therapy genotypes that have been mapped to
        a cdna variant
        :return: tuple of query results
        """

        sql = """
            SELECT distinct
              tg.id as therapy_genotype_id,
              tg.comment as genotype_label

            FROM therapy_genotype tg
            JOIN therapy_variant tv
            ON tg.id = tv.therapy_genotype

            JOIN protein_variant pv
            ON tv.protein_variant = pv.id

            JOIN cdna_variant cdna
            ON pv.id = cdna.protein_variant;
        """
        self.cursor.execute(sql)
        results = self.cursor.fetchall()
        return results

    def _get_fusion_copy_any_mutations(self):
        """
        STATUS: incomplete
        Query database to therapy genotypes that have been mapped to
        a protein variant (rather than a cdna variant)
        :return: tuple of query results
        """

        sql = """
            SELECT distinct
              tg.id as therapy_genotype_id,
              tg.comment as genotype_label,
              gene.description as ref_gene_for_fusion_or_copy,
              gf.description as gene_fusion,
              cg.description as copy_gene

            FROM therapy_genotype tg
            JOIN therapy_variant tv
            ON tg.id = tv.therapy_genotype

            LEFT OUTER JOIN gene
            ON tv.gene = gene.id

            LEFT OUTER JOIN gene gf
            ON tv.gene_fusion = gf.id

            LEFT OUTER JOIN gene cg
            ON tv.copy_gene = cg.id

            WHERE gene.description IS NOT NULL;
        """
        self.cursor.execute(sql)
        results = self.cursor.fetchall()
        return results

    def _get_genotypes_with_no_protein_cdna_mapping(self):
        """
        STATUS: incomplete
        Query database to therapy genotypes that have been mapped to
        a protein variant (rather than a cdna variant)
        :return: tuple of query results
        """

        sql = """
            SELECT distinct
              tg.id as therapy_genotype_id,
              tg.comment as genotype_label,
              tv.amino_acid_start,
              tv.amino_acid_end,
              transcript.description as transcript_id,
              protein_variant_type.description as protein_variant_type,
              gene.description as gene_fusion


            FROM therapy_genotype tg
            JOIN therapy_variant tv
            ON tg.id = tv.therapy_genotype

            LEFT OUTER JOIN transcript
            ON tv.transcript = transcript.id

            LEFT OUTER JOIN protein_variant_type
            ON tv.protein_variant_type = protein_variant_type.id

            LEFT OUTER JOIN gene
            ON tv.gene_fusion = gene.id

            WHERE tv.protein_variant IS NULL
            AND tv.gene IS NULL;
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

