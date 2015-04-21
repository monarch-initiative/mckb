NCIT = http://ncicb.nci.nih.gov/xml/owl/EVS/Thesaurus.owl#

neoplasm.owl: Thesaurus.owl
	owltools $< --make-subset-by-properties //  --remove-annotation-assertions -l --extract-module -d Neoplasm -o -f omn --prefix ncit $(NCIT) $@


drug.obo:
	blip -debug "ontol_db:subclassRT(+,+)" ontol-subset -r chebi -query "class(Drug,drug),subclassRT(ID,X),parent(X,has_role,Role),subclassRT(Role,Drug)" -rel subclass -to obo > $@

drug.owl: drug.obo
	owltools $< --remove-axiom-annotations --remove-dangling -o -f omn --prefix obo http://purl.obolibrary.org/obo/ --prefix oio http://www.geneontology.org/formats/oboInOwl#  $@

drug-lite.owl: drug.owl
	owltools $< --remove-annotation-assertions -l  -o -f omn --prefix obo http://purl.obolibrary.org/obo/ --prefix oio http://www.geneontology.org/formats/oboInOwl#  $@
