from mckb.sources.MySQLSource import MySQLSource
from dipper.models.Dataset import Dataset
from dipper.utils.GraphUtils import GraphUtils
from dipper.models.Genotype import Genotype
from dipper.models.InteractionAssoc import InteractionAssoc
from dipper import curie_map
import tempfile
import gzip
import logging
import os

logger = logging.getLogger(__name__)


class CGD(MySQLSource):
    """
    Test data source for cancer knowledge base
    """
    static_files = {
        'test_data': {'file': 'g2p.sql.gz'}
    }

    def __init__(self, database, username, password, host=None):
        super().__init__('cgd', database, username, password, host)
        self.dataset = Dataset('cgd', 'cgd', 'http://ga4gh.org')
        self.rawdir = 'resources'

    def parse(self):
        """
        Override Source.parse()
        Args:
            :param
        Returns:
            :return None
        """
        (connection, cursor) = self._connect_to_database()
        logger.debug("Checking if database is empty")
        is_db_empty = self.check_if_db_is_empty(cursor)
        if is_db_empty:
            file = '/'.join((self.rawdir,
                                  self.static_files['test_data']['file']))
            logger.debug("Loading data into database from file {0}".format(file))
            self._load_data_from_dump_file(file)
        else:
            logger.debug("Database contains tables, "
                         "skipping load from dump file")

        disease_drug_geno_list = self._get_disease_drug_genotype_relationship(cursor)
        self.add_disease_drug_genotype_to_graph(disease_drug_geno_list)
        self.load_bindings()
        self._disconnect_from_database(cursor, connection)
        return

    def add_genotype_info_to_graph(self, table):
        """
        :param table: iterable of iterables
        :return: None
        """
        gu = GraphUtils(curie_map.get())
        geno = Genotype(self.graph)

        for row in table:
            pass
        return

    def add_disease_drug_genotype_to_graph(self, table):
        """
        :param table: iterable of iterables, for example, a tuple of tuples
                      from _get_disease_drug_genotype_relationship
        :return: None
        """
        gu = GraphUtils(curie_map.get())
        geno = Genotype(self.graph)

        for row in table:
            (genotype_key, genotype_label, diagnoses_key, diagnoses,
             specific_diagnosis_id, specific_diagnosis, organ, relationship,
             drug_key, drug, therapy_status, pubmed_id) = row

            # Arbitrary IDs to be replaced by ontology mappings
            population_id = self.make_id('cgd{0}{1}'.format(genotype_key,
                                                            genotype_label))
            genotype_id = self.make_id('cgd-genotype{0}'.format(genotype_key))
            phenotype_id = self.make_id('cgd-phenotype{0}'.format(diagnoses_key))
            relationship_id = ("MONARCH:{0}".format(relationship)).replace(" ", "_")
            drug_id = self.make_id('cgd-drug{0}'.format(drug_key))

            # Add individuals/classes
            gu.addIndividualToGraph(self.graph, population_id, None,
                                    geno.genoparts['population'])
            gu.addClassToGraph(self.graph, phenotype_id, diagnoses)
            gu.addClassToGraph(self.graph, drug_id, drug)
            gu.loadObjectProperties(self.graph, {relationship:relationship_id})
            geno.addGenotype(genotype_id, genotype_label)

            # Add triples
            gu.addTriple(self.graph, population_id,
                         geno.properties['has_genotype'], genotype_id)
            gu.addTriple(self.graph, population_id,
                         geno.properties['has_phenotype'], phenotype_id)
            gu.addTriple(self.graph, population_id, relationship_id, drug_id)

            # Add 1 association per above triple,
            # see https://github.com/monarch-initiative/mckb/issues/1
            # refactor using generic associations,
            # see https://github.com/monarch-initiative/dipper/issues/96
            if pubmed_id is not None:
                source_id = "PMID:{0}".format(pubmed_id)
                evidence = 'ECO:0000033'

                genotype_annot = self.make_id("{0}{1}".format(population_id, genotype_label))
                phenotype_annot = self.make_id("{0}{1}".format(population_id, diagnoses))
                drug_annot = self.make_id("{0}{1}".format(population_id, drug))

                genotype_assoc = InteractionAssoc(genotype_annot, population_id,
                                                  genotype_id, source_id, evidence)
                phenotype_assoc = InteractionAssoc(phenotype_annot, population_id,
                                                   phenotype_id, source_id, evidence)
                drug_assoc = InteractionAssoc(drug_annot, population_id,
                                              drug_id, source_id, evidence)

                genotype_assoc.rel = geno.properties['has_genotype']
                phenotype_assoc.rel = geno.properties['has_phenotype']
                drug_assoc.rel = relationship_id

                genotype_assoc.addAssociationToGraph(self.graph)
                phenotype_assoc.addAssociationToGraph(self.graph)
                drug_assoc.addAssociationToGraph(self.graph)

        return

    def _get_disease_drug_genotype_relationship(self, cursor):
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
              organs.description as organ,
              ta.description as relationship,
              tc.id as drug_id,
              tc.description as drug,
              therapy_status.description as therapy_status,
              tgp.pub_med_id as pubmed_id
            FROM therapy_genotype tg
            JOIN diagnoses
            ON tg.diagnosis = diagnoses.id

            JOIN therapeutic_association as ta
            ON tg.therapeutic_association = ta.id

            JOIN therapeutic_context tc
            ON tg.therapeutic_context = tc.id

            LEFT OUTER JOIN therapy_status
            ON tg.therapy_status = therapy_status.id

            LEFT OUTER JOIN specific_diagnosis
            ON tg.specific_diagnosis = specific_diagnosis.id

            LEFT OUTER JOIN therapy_genotype_publication as tgp
            ON tg.id = tgp.therapy_genotype

            LEFT OUTER JOIN organs
            ON tg.organ = organs.id;
        """
        cursor.execute(sql)
        results = cursor.fetchall()
        return results

    def _get_genotype_protein_info(self, cursor):
        """
        Select out genotypes that have protein variant information but no
        cdna variant information
        :return: tuple of query results
        """

        sql = """
            SELECT distinct
              tg.id as therapy_genotype_id,
              tg.comment as genotype_label,
              pv.genotype_amino_acid_onel as aa_var,
              tv.amino_acid_start,
              tv.amino_acid_end,
              tv.genomic_start,
              tv.genomic_end,
              pv.amino_acid_position,
              transcript.description as transcript_id,
              transcript_priority.description as transcript_priority,
              protein_variant_type.description as protein_variant_type,
              functional_impact.description as functional_impact,
              stop_gain_loss.description as stop_gain_loss,
              trg.description as transcript_gene,
              pv.pub_med_ids as protein_variant_pubmed_ids

            FROM therapy_genotype tg
            JOIN therapy_variant tv
            ON tg.id = tv.therapy_genotype

            JOIN protein_variant pv
            ON tv.protein_variant = pv.id

            LEFT OUTER JOIN cdna_variant cdna
            ON pv.id = cdna.protein_variant

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
            ON transcript.gene = trg.id

            WHERE cdna.protein_variant IS NULL;
        """
        cursor.execute(sql)
        results = cursor.fetchall()
        return results

    def _get_genotype_cdna_info(self, cursor):
        """
        Query database to therapy genotypes that have been mapped to
        a cdna variant
        :return: tuple of query results
        """

        sql = """
            SELECT distinct
              tg.id as therapy_genotype_id,
              tg.comment as genotype_label,
              pv.genotype_amino_acid_onel as aa_var,
              tv.amino_acid_start,
              tv.amino_acid_end,
              tv.genomic_start,
              tv.genomic_end,
              pv.amino_acid_position,
              transcript.description as transcript_id,
              transcript_priority.description as transcript_priority,
              protein_variant_type.description as protein_variant_type,
              functional_impact.description as functional_impact,
              stop_gain_loss.description as stop_gain_loss,
              trg.description as transcript_gene,
              pv.pub_med_ids as protein_variant_pubmed_ids
              gene.description as variant_gene,
              cdna.base_pair_position,
              cdna.genotype_cdna,
              genomic_variant.cosmic_id,
              genomic_variant.db_snp_id,
              genomic_variant.position_start,
              genomic_variant.position_end,
              genomic_variant.reference_base,
              genomic_variant.variant_base,
              genomic_variant.primary_transcript_exons,
              genomic_variant.primary_transcript_variant_sub_types,
              variant_type.description as variant_type,
              chromosome.description as chromosome,
              genome_build.description as genome_build,
              genome_build.build_version as build_version,
              genome_build.build_date as build_date,

            FROM therapy_genotype tg
            JOIN therapy_variant tv
            ON tg.id = tv.therapy_genotype

            JOIN protein_variant pv
            ON tv.protein_variant = pv.id

            JOIN cdna_variant cdna
            ON pv.id = cdna.protein_variant

            LEFT OUTER JOIN transcript
            ON cdna.transcript = transcript.id

            LEFT OUTER JOIN genomic_variant
            on cdna.genomic_variant = genomic_variant.id

            LEFT OUTER JOIN transcript_priority
            ON transcript.transcript_priority = transcript_priority.id

            LEFT OUTER JOIN protein_variant_type
            ON pv.protein_variant_type = protein_variant_type.id

            LEFT OUTER JOIN functional_impact
            ON pv.functional_impact = functional_impact.id

            LEFT OUTER JOIN stop_gain_loss
            ON pv.stop_gain_loss = stop_gain_loss.id

            LEFT OUTER JOIN variant_type
            ON genomic_variant.variant_type = variant_type.id

            LEFT OUTER JOIN chromosome
            ON genomic_variant.chromosome = chromosome.id

            LEFT OUTER JOIN genome_build
            ON genomic_variant.genome_build = genome_build.id

            LEFT OUTER JOIN gene trg
            ON transcript.gene = trg.id

            LEFT OUTER JOIN gene
            ON genomic_variant.gene = gene.id;
        """
        cursor.execute(sql)
        results = cursor.fetchall()
        return results

    def _get_fusion_copy_any_mutation_genotypes(self, cursor):
        """
        Get genotypes with a gene mapping but no protein variant mapping
        Typically, this will capture fusion genes and copy/gain loss
        along with mutations labelled "any" (any mutation of gene X)
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

            JOIN gene
            ON tv.gene = gene.id

            LEFT OUTER JOIN gene gf
            ON tv.gene_fusion = gf.id

            LEFT OUTER JOIN gene cg
            ON tv.copy_gene = cg.id

            WHERE tv.protein_variant IS NULL;
        """
        cursor.execute(sql)
        results = cursor.fetchall()
        return results

    def _get_genotypes_with_no_gene_protein_cdna_mapping(self, cursor):
        """
        Get genotypes with no protein_variant mapping or gene
        mapping.  This will capture rearrangement fusion genes,
        missense mutations where there are multiple mutations or
        the specific mutation is unknown, amplification copy number genes,
        and indels

        :return: tuple of query results
        """

        sql = """
            SELECT distinct
              tg.id as therapy_genotype_id,
              tg.comment as genotype_label,
              tv.amino_acid_start,
              tv.amino_acid_end,
              variant_type.description,
              transcript.description as transcript_id,
              protein_variant_type.description as protein_variant_type,
              gene.description as gene_fusion,
              g.description as copy_gene,
              cns.description as copy_number_result


            FROM therapy_genotype tg
            JOIN therapy_variant tv
            ON tg.id = tv.therapy_genotype

            LEFT OUTER JOIN transcript
            ON tv.transcript = transcript.id

            LEFT OUTER JOIN protein_variant_type
            ON tv.protein_variant_type = protein_variant_type.id

            LEFT OUTER JOIN variant_type
            ON tv.variant_type_aa_coords = variant_type.id

            LEFT OUTER JOIN gene
            ON tv.gene_fusion = gene.id

            LEFT OUTER JOIN gene g
            ON tv.copy_gene = g.id

            LEFT OUTER JOIN copy_number_result cns
            ON tv.copy_number_result = cns.id

            WHERE tv.protein_variant IS NULL
            AND tv.gene IS NULL;
        """
        cursor.execute(sql)
        results = cursor.fetchall()
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

