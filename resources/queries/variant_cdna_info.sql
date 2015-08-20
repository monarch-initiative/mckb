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
ON cdna.genomic_variant = genomic_variant.id

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