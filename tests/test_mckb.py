from dipper import curie_map
from rdflib import Graph

import unittest
import logging

logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)

# Example 1
"""
indiv1 rdf:type glioma_of_the_brain
  rdfs:label “glioma of the brain with MGMT any promoter methylation”
  GENO:has_genotype gtindiv1
  monarch:responds_with_increased_benefit DrugBank:DB00853

DrugBank:DB00853 a owl:Class  #drugs will have to be curated.
DrugBank:DB00853 rdfs:label “temozolomide”

gtindiv1 rdf:type GENO:methylation  #we might need to curate these kinds of genotypes
  RO:contained_in SO:promoter
  RO:targets NCBIGene:4255    #need to curate gene ids/fetch from NCBI
  refs:label “MGMT any promoter methylation”
"""

# Example 2
"""
indiv2 rdf:type ductal_carcinoma
indiv2 rdfs:label “ductal carcinoma of the breast with ERBB2 G309A missense mutation”
indiv2 GENO:has_genotype gtindiv2

indiv2 monarch:reduced_sensitivity_toward CHEBI:49603      #lapatinib

gtindiv2 rdf:type GENO:sequence_alteration,SO:missense_variant
  rdfs:label “ERBB2 G309A missense mutation”
  GENO:is_sequence_variant_of NCBIGene:ERBB2
  faldo:location [
    a faldo:position,
    faldo:begin 390,  #preferable to use proper genomic location
    faldo:end 390
    reference: protein_id]   #if we don’t map to the genomic level, could use the protein id here for relative position
  GENO:reference_nucleotide “X”  #not sure we have these relationships?
  GENO:altered_nucleotide “Y”
  GENO:results_in_amino_acid_change “A”

annotation180 a annotation
  hasSubject indiv2
  hasObject gtindiv2
  hasPredicate monarch:reduced_sensitivity_toward
  dc:source PMID:23220880
"""


class VariantDiseaseDrugTestCase(unittest.TestCase):

    def setUp(self):
        return

    def tearDown(self):
        return

    def test_foo(self):
        self.assertEqual(1, 1)

if __name__ == '__main__':
    unittest.main()