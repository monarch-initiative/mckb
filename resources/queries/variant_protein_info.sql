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