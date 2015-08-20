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