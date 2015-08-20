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