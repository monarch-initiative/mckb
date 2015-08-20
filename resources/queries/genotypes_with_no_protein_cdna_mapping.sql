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