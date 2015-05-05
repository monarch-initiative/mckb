from mckb.sources.MySQLSource import MySQLSource
from mckb.sources.CGDOntologyMap import CGDOntologyMap
from dipper.models.Dataset import Dataset
from dipper.utils.GraphUtils import GraphUtils
from dipper.models.Genotype import Genotype
from dipper.models.G2PAssoc import G2PAssoc
from dipper.models.GenomicFeature import Feature, makeChromID
from dipper import curie_map
import tempfile
import gzip
import logging
import os
import csv
import re
import hashlib

logger = logging.getLogger(__name__)


class CGD(MySQLSource):
    """
    Test data source for cancer knowledge base
    """
    files = {
        'transcript_xrefs': {
            'file': 'CCDS2UniProtKB.current.txt',
            'url': 'ftp://ftp.ncbi.nlm.nih.gov/pub/CCDS/current_human/CCDS2UniProtKB.current.txt'
        }
    }
    static_files = {
        'cgd': {
            'file': '../../resources/g2p.sql.gz'
        },
        'ncbi_gene_mappings': {
            'file': '../../resources/mappings/gene.tsv'
        }
    }

    def __init__(self, database, username, password, host=None):
        super().__init__('cgd', database, username, password, host)
        self.dataset = Dataset('cgd', 'cgd', 'http://ga4gh.org')
        self.gene_map = {}
        self.disease_map = {}
        self.drug_map = {}
        self.transcript_xrefs = {'RefSeq': {}, 'UniProt': {}}

    def fetch(self, is_dl_forced):
        """
        Override Source.fetch()
        Fetches resources from CTD using the CTD.files dictionary
        Args:
            :param is_dl_forced (bool): Force download
        Returns:
            :return None
        """
        self.get_files(is_dl_forced)
        return

    def parse(self):
        """
        Override Source.parse()
        :param
        :return None
        """
        (connection, cursor) = self._connect_to_database()

        logger.info("Checking if database is empty")
        is_db_empty = self.check_if_db_is_empty(cursor)
        if is_db_empty:
            file = self.static_files['cgd']['file']
            logger.info("Loading data into database from file {0}".format(file))
            self._load_data_from_dump_file(file)
        else:
            logger.info("Database contains tables, "
                        "skipping load from dump file")

        ontology_map = CGDOntologyMap('cgd-ontology-mappings')
        ontology_map.parse()
        self.gene_map = ontology_map.gene_map
        self.disease_map = ontology_map.disease_map
        self.drug_map = ontology_map.drug_map

        ccds_xref_file = '/'.join((self.rawdir,
                                   self.files['transcript_xrefs']['file']))

        self.set_transcript_xrefs(ccds_xref_file)

        variant_protein_assocs = self. _get_variant_protein_info(cursor)
        self.add_variant_info_to_graph(variant_protein_assocs)
        variant_cdna_assocs = self._get_variant_cdna_info(cursor)
        self.add_variant_info_to_graph(variant_cdna_assocs)

        disease_drug_geno_list = self._get_disease_drug_variant_relationship(cursor)
        self.add_disease_drug_variant_to_graph(disease_drug_geno_list)

        self.load_bindings()
        self._disconnect_from_database(cursor, connection)
        return

    def add_variant_info_to_graph(self, table):
        """
        Takes an iterable or iterables as input with one of two structures,
        Structure 1: Only protein variant information
        optional indices can be null:
        [[variant_key, variant_label, amino_acid_variant,
          amino_acid_position (optional), transcript_id, transcript_priority,
          protein_variant_type (optional), functional_impact,
          stop_gain_loss (optional), transcript_gene,
          protein_variant_source (optional)]]

        Structure 2: Protein and cDNA information
        optional indices can be null:
        [[variant_key variant_label, amino_acid_variant,
          amino_acid_position (optional), transcript_id, transcript_priority,
          protein_variant_type (optional), functional_impact,
          stop_gain_loss (optional), transcript_gene,
          protein_variant_source (optional), variant_gene, bp_pos,
          variant_cdna, cosmic_id (optional), db_snp_id (optional),
          genome_pos_start, genome_pos_end, ref_base, variant_base,
          primary_transcript_exons, primary_transcript_variant_sub_types,
          variant_type, chromosome, genome_build, build_version, build_date]]

        :param table: iterable of iterables
        :return: None
        """
        gu = GraphUtils(curie_map.get())
        geno = Genotype(self.graph)

        for row in table:
            self._add_variant_protein_variant_assoc_to_graph(row)
            if len(row) > 11:
                self._add_variant_cdna_variant_assoc_to_graph(row)

        return

    def _add_variant_protein_variant_assoc_to_graph(self, row):
        """
        Generates relationships between variants and protein variants
        given a row of data
        :param iterable: row of data, see add_variant_info_to_graph()
                                      docstring for expected structure
        :return None
        """
        gu = GraphUtils(curie_map.get())
        geno = Genotype(self.graph)
        is_missense = False
        is_literal = True

        (variant_key, variant_label, amino_acid_variant, amino_acid_position,
         transcript_id, transcript_priority, protein_variant_type,
         functional_impact, stop_gain_loss, transcript_gene,
         protein_variant_source) = row[0:11]

        variant_id = self.make_cgd_id('variant{0}'.format(variant_key))

        transcript_curie = self._make_transcript_curie(transcript_id)
        uniprot_curie = self._make_uniprot_polypeptide_curie(transcript_id)
        ncbi_protein_curie = self._make_ncbi_polypeptide_curie(transcript_id)

        geno.addGenotype(variant_id, variant_label,
                         geno.genoparts['sequence_alteration'])

        # Make fake amino acid sequence in case we
        # can't get a CCDS to Uniprot and/or NCBI Protein mapping
        aa_seq_id = self.make_cgd_id('transcript{0}'.format(amino_acid_variant))

        # Add Transcript:
        geno.addTranscript(variant_id, transcript_curie, transcript_id,
                           geno.genoparts['transcript'])

        # Add polypeptide
        if ncbi_protein_curie is not None:
            geno.addPolypeptide(ncbi_protein_curie,
                                self.transcript_xrefs['RefSeq'][transcript_id],
                                transcript_curie)
            aa_seq_id = ncbi_protein_curie
        if uniprot_curie is not None:
            geno.addPolypeptide(uniprot_curie,
                                self.transcript_xrefs['UniProt'][transcript_id],
                                transcript_curie)
            # Overrides ncbi_protein_curie,
            # but we set them as equal individuals below
            aa_seq_id = uniprot_curie

        if ncbi_protein_curie is not None and uniprot_curie is not None:
            gu.addSameIndividual(self.graph, ncbi_protein_curie, uniprot_curie)
        else:
            aa_seq_id = self.make_cgd_id('transcript{0}'.format(amino_acid_variant))

        if protein_variant_type == 'nonsynonymous - missense' \
                or re.search(r'missense', variant_label):
            is_missense = True
            geno.addGenotype(variant_id, variant_label,
                             geno.genoparts['missense_variant'])

        # Get gene ID from gene map
        self._add_variant_gene_relationship(variant_id, transcript_gene)

        amino_acid_regex = re.compile(r'^p\.([A-Za-z]{1,3})(\d+)([A-Za-z]{1,3})$')

        if is_missense:
            match = re.match(amino_acid_regex, amino_acid_variant.rstrip())
        else:
            match = None

        if match is not None:
            ref_amino_acid = match.group(1)
            position = match.group(2)
            altered_amino_acid = match.group(3)
        else:
            logger.debug("Could not parse amino acid information"
                         " from {0} variant:"
                         " {1} type: {2}".format(amino_acid_variant,
                                                 variant_label,
                                                 protein_variant_type))

        # Add amino acid change to model
        if is_missense is True and match is not None:
            gu.addTriple(self.graph, variant_id,
                         geno.properties['reference_amino_acid'],
                         ref_amino_acid, is_literal)
            gu.addTriple(self.graph, variant_id,
                         geno.properties['results_in_amino_acid_change'],
                         altered_amino_acid, is_literal)

            aa_region_id = ":_{0}{1}Region".format(variant_id, aa_seq_id)
            self._add_feature_with_coords(variant_id, position,
                                          position, aa_seq_id, aa_region_id)

        return

    def _add_variant_cdna_variant_assoc_to_graph(self, row):
        """
        Generates relationships between variants and cDNA variants
        given a row of data
        :param iterable: row of data, see add_variant_info_to_graph()
                                      docstring for expected structure.
                                      Only applicable for structure 2.
        :return None
        """
        gu = GraphUtils(curie_map.get())
        geno = Genotype(self.graph)
        is_literal = True

        (variant_key, variant_label, amino_acid_variant, amino_acid_position,
         transcript_id, transcript_priority, protein_variant_type,
         functional_impact, stop_gain_loss, transcript_gene,
         protein_variant_source, variant_gene, bp_pos, variant_cdna,
         cosmic_id, db_snp_id, genome_pos_start, genome_pos_end, ref_base,
         variant_base, primary_transcript_exons,
         primary_transcript_variant_sub_types, variant_type, chromosome,
         genome_build, build_version, build_date) = row

        variant_id = self.make_cgd_id('variant{0}'.format(variant_key))

        # Add gene
        self._add_variant_gene_relationship(variant_id, variant_gene)

        # Transcript reference for nucleotide position
        transcript_curie = self._make_transcript_curie(transcript_id)

        # Make region IDs
        cdna_region_id = ":_{0}{1}Region".format(variant_id, transcript_curie)
        chrom_region_id = ":_{0}{1}{2}Region".format(variant_id, genome_build,
                                                     chromosome)

        # Add the genome build
        genome_label = "Human"
        build_id = "UCSC:{0}".format(genome_build)
        taxon_id = 'NCBITaxon:9606'
        geno.addGenome(taxon_id, genome_label)
        geno.addReferenceGenome(build_id, genome_build, taxon_id)

        # Add chromosome
        chromosome_id = makeChromID(chromosome, genome_build)
        geno.addChromosome(chromosome, taxon_id, genome_label,
                           build_id, genome_build)

        # Add variant coordinates in reference to chromosome
        self._add_feature_with_coords(variant_id,genome_pos_start,
                                      genome_pos_end, chromosome_id, chrom_region_id)

        # Add mutation coordinates in reference to gene
        self._add_feature_with_coords(variant_id, bp_pos,
                                      bp_pos, transcript_curie, cdna_region_id)

        # Add nucleotide mutation
        gu.addTriple(self.graph, variant_id,
                     geno.properties['reference_nucleotide'],
                     ref_base, is_literal)
        gu.addTriple(self.graph, variant_id,
                     geno.properties['altered_nucleotide'],
                     variant_base, is_literal)

        # Add SNP xrefs
        if cosmic_id is not None:
            cosmic_id_list = cosmic_id.split(', ')
            for c_id in cosmic_id_list:
                cosmic_curie = re.sub(r'COSM(\d+)', r'COSMIC:\1', c_id)
                gu.addIndividualToGraph(self.graph, cosmic_curie, c_id,
                                        geno.genoparts['missense_variant'])
                gu.addSameIndividual(self.graph, variant_id, cosmic_curie)
        if db_snp_id is not None:
            db_snp_curie = re.sub(r'rs(\d+)', r'dbSNP:\1', db_snp_id)
            gu.addIndividualToGraph(self.graph, db_snp_curie, db_snp_id,
                                    geno.genoparts['missense_variant'])
            gu.addSameIndividual(self.graph, variant_id, db_snp_curie)

        if (db_snp_id is not None) and (cosmic_id is not None):
            cosmic_id_list = cosmic_id.split(', ')
            for c_id in cosmic_id_list:
                cosmic_curie = re.sub(r'COSM(\d+)', r'COSMIC:\1', c_id)
                gu.addSameIndividual(self.graph, cosmic_curie, db_snp_curie)

        return

    def _add_feature_with_coords(self, feature_id, start_pos, end_pos, reference, region_id):
        """
        :param feature_id: URIRef or Curie - instance of faldo:Position
        :param feature_label: String
        :param feature_type: Object Property
        :param start_pos: int, starting coordinate
        :param end_pos: int, ending coordinate
        :param reference: URIRef or Curie - reference Node (gene, transcript, genome)
        :return: None
        """
        add_region = True
        feature = Feature(feature_id, None, None)
        feature.addFeatureStartLocation(start_pos, reference)
        feature.addFeatureEndLocation(end_pos, reference)
        feature.addFeatureToGraph(self.graph, add_region, region_id)
        return

    def _add_variant_gene_relationship(self, variant_id, hgnc_symbol):
        """
        :param variant_id
        :param hgnc_symbol
        :return: None
        """
        gu = GraphUtils(curie_map.get())
        geno = Genotype(self.graph)
        if hgnc_symbol in self.gene_map:
            gene_id = self.gene_map[hgnc_symbol]
        else:
            gene_id = self.make_cgd_id("{0}{1}".format(variant_id, hgnc_symbol))
            logger.warn("Can't map gene symbol {0} "
                        "to entrez ID".format(hgnc_symbol))
        gu.addClassToGraph(self.graph, gene_id, hgnc_symbol)
        geno.addAlleleOfGene(variant_id, gene_id)
        return

    def add_disease_drug_variant_to_graph(self, table):
        """
        Takes an iterable of iterables as input with the following structure,
        optional indices can be Null:
        [[variant_key, variant_label, diagnoses_key, diagnoses,
          specific_diagnosis, organ, relationship,
          drug_key, drug, therapy_status (optional), pubmed_id(optional)]]

        See ongoing discussion of how to best model here:
        https://github.com/monarch-initiative/mckb/issues/9

        :param table: iterable of iterables, for example, a tuple of tuples
                      from _get_disease_drug_variant_relationship
        :return: None
        """
        gu = GraphUtils(curie_map.get())
        geno = Genotype(self.graph)

        for row in table:
            (variant_key, variant_label, diagnoses_key, diagnoses,
             specific_diagnosis, organ, relationship,
             drug_key, drug_label, therapy_status, pubmed_id) = row

            if specific_diagnosis is not None:
                diagnoses_label = specific_diagnosis
            else:
                diagnoses_label = diagnoses

            # Arbitrary IDs to be replaced by ontology mappings
            variant_id = self.make_cgd_id('variant{0}'.format(variant_key))
            disease_id = self._get_disease_id(diagnoses_key, diagnoses_label)
            relationship_id = ("RO:{0}".format(relationship)).replace(" ", "_")
            drug_id = self._get_drug_id(drug_key, drug_label)

            """
            Here we can either create an instance of the disease class
            or link the disease class directly.  For now we will link
            the class directly as this is easier to query out with the
            current scigraph services, but down the line this should be
            added back in along with updating the cypher query in the
            scigraph layer
            """

            # disease_instance_id = self.make_cgd_id('disease{0}{1}'.format(
            #    diagnoses_label, variant_key))
            # disease_instance_label = "{0} linked to variant {1}".format(diagnoses_label, variant_label)

            # Reified association for disease caused_by genotype
            variant_disease_annot = self.make_cgd_id("assoc{0}{1}".format(variant_key, diagnoses_label))

            # Add individuals/classes
            gu.addClassToGraph(self.graph, disease_id, diagnoses_label, 'DOID:4')

            gu.addClassToGraph(self.graph, drug_id, drug_label, 'CHEBI:23888')
            #gu.addIndividualToGraph(self.graph, disease_instance_id, disease_instance_label,
            #                        disease_id)
            gu.loadObjectProperties(self.graph, {relationship: relationship_id})

            if pubmed_id is not None:
                source_id = "PMID:{0}".format(pubmed_id)
                evidence = 'ECO:0000033'
            else:
                source_id = None
                evidence = None

            variant_phenotype_assoc = G2PAssoc(variant_disease_annot,
                                               variant_id,
                                               disease_id,
                                               source_id, evidence)
            variant_phenotype_assoc.addAssociationToGraph(self.graph)
            gu.addTriple(self.graph, variant_disease_annot, relationship_id, drug_id)

        return

    def _make_transcript_curie(self, transcript_id):
        """
        :param transript_id:
        :return: transcript_curie
        """
        if re.match(r'^CCDS', transcript_id):
            transcript_curie = re.sub(r'(CCDS)(\d+)', r'\1:\2', transcript_id)
        elif re.match(r'^NM', transcript_id):
            transcript_curie = re.sub(r'(NM_)(\d+)', r'GenBank:\1\2', transcript_id)

        return transcript_curie

    def _make_uniprot_polypeptide_curie(self, transcript_id):
        """
        :param transript_id:
        :return: uniprot_curie
        """
        uniprot_curie = None
        if transcript_id in self.transcript_xrefs['UniProt']:
            uniprot_id = self.transcript_xrefs['UniProt'][transcript_id]
            if re.search(r'\-', uniprot_id):
                uniprot_curie = re.sub(r'([\w\d]+)(\-\d+)', r'UniProtKB:\1#\1\2',
                                       uniprot_id)
            else:
                uniprot_curie = "UniProtKB:{0}".format(uniprot_id)
        else:
            # raise Exception("Could not find {0} in CCDS cross reference"
            #                " dictionary".format(transcript_id))
            pass
        return uniprot_curie

    def _make_ncbi_polypeptide_curie(self, transcript_id):
        """
        :param transript_id:
        :return: ncbi_protein_curie
        """
        ncbi_protein_curie = None
        if transcript_id in self.transcript_xrefs['RefSeq']:
            protein_id = self.transcript_xrefs['RefSeq'][transcript_id]
            ncbi_protein_curie = "NCBIProtein:{0}".format(protein_id)
        else:
            # raise Exception("Could not find {0} in CCDS cross reference"
            #                " dictionary".format(transcript_id))
            pass
        return ncbi_protein_curie

    def _get_disease_drug_variant_relationship(self, cursor):
        """
        Query database to get disease-drug-variant associations
        :param: PyMySQL Cursor object generated
                from self._connect_to_database()
        :return: tuple of query results
        """

        sql = """
            SELECT distinct
              tg.id as genotype_id,
              tg.comment as genotype_label,
              diagnoses.id as diagnoses_id,
              diagnoses.description as diagnoses,
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

    def _get_variant_protein_info(self, cursor):
        """
        Select out genotypes that have protein variant information but no
        cdna variant information
        :param: PyMySQL Cursor object generated
                from self._connect_to_database()
        :return: tuple of query results
        """

        sql = """
            SELECT distinct
              tg.id as therapy_genotype_id,
              tg.comment as genotype_label,
              pv.genotype_amino_acid_onel as aa_var,
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

    def _get_variant_cdna_info(self, cursor):
        """
        Query database to therapy variants that have been mapped to
        a cdna variant
        :param: PyMySQL Cursor object generated
                from self._connect_to_database()
        :return: tuple of query results
        """

        sql = """
            SELECT distinct
              tg.id as therapy_genotype_id,
              tg.comment as genotype_label,
              pv.genotype_amino_acid_onel as aa_var,
              pv.amino_acid_position,
              transcript.description as transcript_id,
              transcript_priority.description as transcript_priority,
              protein_variant_type.description as protein_variant_type,
              functional_impact.description as functional_impact,
              stop_gain_loss.description as stop_gain_loss,
              trg.description as transcript_gene,
              pv.pub_med_ids as protein_variant_pubmed_ids,
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
              genome_build.build_date as build_date

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
        :param: PyMySQL Cursor object generated
                from self._connect_to_database()
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

        :param: PyMySQL Cursor object generated
                from self._connect_to_database()
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
        :param file: file path containing sql dump file
        :return: None
        """
        gz_file = gzip.open(file, 'rb')
        with tempfile.NamedTemporaryFile(mode='w+b') as f:
            f.write(gz_file.read())
            gz_file.close()
            os.system("mysql -h {0} -p{1} -u {2} -D {3} < {4}".format(self.host, self.password,
                      self.username, self.database, f.name))
        return

    def set_transcript_xrefs(self, mapping_file):
        """
        Sets transcript_xref instance variable in this structure:
                {
                     'RefSeq': {
                         'CCDS1234.5': 'NP_12345',
                         ...
                     }
                     'UniProt': {
                         'CCDS1234.5': 'Q12345'.
                         ...
                     }
                 }
        :param mapping_file: String, local path to file containing
                                     CCDS xref mappings
        :return: None
        """
        with open(mapping_file, 'rt') as tsvfile:
            reader = csv.reader(tsvfile, delimiter="\t")
            for row in reader:
                (ccds_id, ncbi_protein_id, uniprot_id) = row[0:3]
                if re.match('^#', ccds_id):
                    next
                else:
                    self.transcript_xrefs['RefSeq'][ccds_id] = ncbi_protein_id
                    self.transcript_xrefs['UniProt'][ccds_id] = uniprot_id

        return

    def make_cgd_id(self, long_string):
        """
        Make identifier to CGD prefix
        :param long_string:
        :return: Curie with CGD prefix
        """
        byte_string = long_string.encode("utf-8")
        md5sum = ':'.join(('CGD', hashlib.md5(byte_string).hexdigest()))
        md5sum = md5sum[0:12]
        return md5sum

    def _get_disease_id(self, diagnoses_key, diagnoses_label):
        if (diagnoses_label in self.disease_map) \
                and (self.disease_map[diagnoses_label] != ''):
            disease_id = self.disease_map[diagnoses_label]
        else:
            logger.debug("Can't map disease {0}"
                         " to ontology".format(diagnoses_label))
            disease_id = self.make_cgd_id('disease{0}{1}'.format(
                diagnoses_key, diagnoses_label))

        return disease_id

    def _get_drug_id(self, drug_key, drug_label):
        if (drug_label in self.drug_map) \
                and (self.drug_map[drug_label] != ''):
            drug_id = self.drug_map[drug_label]
        else:
            logger.debug("Can't map drug {0}"
                         " to ontology".format(drug_label))
            drug_id = self.make_cgd_id('drug{0}'.format(drug_key))

        return drug_id


