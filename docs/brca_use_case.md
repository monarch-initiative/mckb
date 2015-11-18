# BRCA Genotype-Phenotype Association and Evidence Modeling

### Table of Contents
1. [Background](#i-background)
2. [Requirements Analysis for Modeling Evidence and Provenance](#ii-requirements-analysis-for-modeling-evidence-and-provenance)
3. [Candidate Model for Describing Evidence and Provenance](#iii-candidate-model-for-describing-evidence-and-provenance)
4. [Collaborators](#iv-collaborators)
5. [References](#v-references)


## I. Background

The scientific and medical communities have a tremendous and urgent need for a comprehensive data store of variation in the *BRCA1* and *2* genes.  Variation in these genes can indicate genetic predisposition to breast and ovarian cancer, leading causes of death that claim greater than fifty thousand lives in the United States annually [1].  There is enormous public interest in *BRCA* testing, both because of growing public interest in genetic testing and following Angelina Jolie’s double mastectomy and the ensuing "Angelina Effect".  Yet not all *BRCA* variants are pathogenic, and some cause diseases other than the commonly attributed breast cancer.  There is therefore a significant scientific need to catalog the *BRCA* variants and their pathogenicity, and an unprecedented opportunity to do so in the wake of the Supreme Court’s decision to strike down the Myriad Genetics patent.  However, partly due to the patent litigation history, there is currently no comprehensive data source on variation within the *BRCA* genes.  As a result clinicians are frequently working with incomplete knowledge when determining a patient’s risk.  In addition, datasets containing patient level phenotyping data, such as TCGA, do not leverage this pathogenicity information (figure 1).

![Data Landscape](https://raw.githubusercontent.com/monarch-initiative/mckb/master/docs/image/data_landscape.png)

**Figure 1**. Catalog of public and controlled use datasets containing information about the BRCA variant

Another challenge faced by researchers and clinicians in understanding and acting on BRCA data is the diversity and complexity of evidence that is used to support characterization of a given variant’s pathogenicity. Agents must be able to evaluate all lines of evidence to determine their confidence in the pathogenicity of a given variant. In support of this, curation efforts around variant classification typically follow detailed and specific protocols in documenting and classifying evidence for pathogenicity. This level of rigor is critical due to the direct clinical implications for pathogenicity calls in the treatment of patients, and the fact that different lines of evidence for the pathogenicity of a given variant is often conflicting. Protocols for variant classification clearly define evidence by which variants are assessed, and rules for combining multiple lines of evidence to arrive at a pathogenicity call along a spectrum of confidence from decisively 'pathogenic' to definitively 'benign'. 

While many labs and initiatives have defined their own protocols for variant classification, the guidelines set by the American College of Medical Genetics and Genomics (ACMG) represent a domain standard that helps clinicians and curators define variant pathogenicity in a more informed and consistent manner [2]. The guidelines define 16 criteria that support a variant being pathogenic, which are categorized as representing Very Strong (PVS), Strong (PS), Moderate, or Supporting (PP) evidence. Likewise, they define 12 criteria that support a variant being benign, which are categorized as representing Stand-alone (BA), Strong (BS), or Supporting (BP) evidence. In a classification pipeline, variants are annotated with the criteria they meet. For example a variant predicted to cause a loss-of-function in its gene product get a PVS1 annotation as this is very strong evidence that the variant will be pathogenic), and these scores are then combined according to a set of specific rules to arrive at a final classification along a five point scale from 'Pathogenic', 'Likely Pathogenic', 'Uncertain Significance', 'Likely Benign', and 'Benign'. This approach is followed by most all efforts to classify cancer variant pathogenicity, using the ACMG or similar set of criteria and heuristics.

## II. Requirements Analysis for Modeling Evidence and Provenance

A model for evidence and provenance around cancer variant classification should be capable of capturing the detail reflected in the complex classification protocols for cancer variant classification - including the distinct lines of evidence that support a pathogenicity call, and the techniques and agents involved in generating them. This is especially important in efforts to integrate pathogenicity evidence across systems, where diverse and conflicting evidence may exist. There are few existing standards and schema to draw from in describing provenance and evidence for scientific claims. Those that do exist are narrowly scoped and inconsistent in their treatment of provenance and evidence - often confusing and conflating these concepts. A critical first step towards defining a broadly applicable and rigorous community standard for describing the evidence and provenance of scientific claims is to clearly define and disentangle these concepts. 

### A. Defining and Distinguishing ‘Provenance’ and ‘Evidence’

Below we present a coherent characterization of these evidence and provenance as related to scientific claims, which will inform our modeling efforts and facilitate its adoption and understanding by the broader community. 

**Evidence** is a collection of information that is used to support a scientific claim or association.  Evidence as ‘information’ can include primary data produced through experimentation or computational analyses, derived data such as statistical calculations and confidence scores that are about this primary data, summaries or reports based on such data, or even intrinsic knowledge or opinion shared by a domain expert.

**Provenance** is commonly defined as a history (or a record of a history) of where an ‘object’ came from and who has owned or modified it. The provenance of an object can be viewed as having two stages: 

1. A history leading to its creation - going back in time from its point of creation to understand how it came to be (i.e. its ancestral lineage)

2. A history since its creation - going forward in time from its point of creation to understand where it has been and how it has evolved (i.e. its life story).  

When the object in question is a scientific claim, the provenance we concern ourselves with is primarily the history leading to its creation. That is, we want to know what processes led to the claim being made, and what entities participated in these processes. Treating provenance as a history (a process), rather than a record of this history (an information artifact), helps us to distinguish it more clearly from evidence (an information artifact). The provenance of a scientific claim then, is the set of processes leading to (i.e. producing evidence for) its assertion.

## B. Use Cases and Competency Questions 

With respect to the evidence and provenance for claims around BRCA variant data, we have defined the following use cases that will inform the scope and structure of our model.

**1. Representation of all key types of evidence and provenance information described in BRCA variant datasets**

Review all BRCA variant datasets will reveal the diversity of evidence and provenance information each source provides. Examples we have encountered in our efforts to date include:

**Evidence (information that supports a claim, which is output from a line of provenance)**

* raw data items (measured values, image data, sequence data, etc)

* evidence codes

* publications

* statistical confidence measures (p-values, z-scores, standard errors, etc)

* conclusions from interpreting data (e.g. pathogenicity criteria from ACMG guidelines)

* tacit knowledge of a domain expert

**Provenance (the process history leading to the claim and key participants)**

* the *types *of assay/technique/study that generated evidence (or actual *instances *of these types of processes)

* the agent(s) who produced evidence for the claim (performed studies/assays, or analyses/calculations)

* the agent(s) who asserted the claim itself (reviewed evidence and deduced the assertion)

* time and place evidence was produced, or the claim was made

* materials used to produce evidence (e.g. model systems, reagents, instruments)

**2. Distinguish separate lines of provenance/evidence for a given association/claim**

An association can have more than one provenance line, each producing a set of evidence. We need a way to group evidence coming separate lines of provenance.  This requires structuring evidence information around organizational nodes representing a single line of evidence.

**3. Capturing conflicting/contradictory evidence**

It is not uncommon for there to be lines of evidence that support and lines that refute a given claim/association.  We need to be able to capture this using our model.  For example, variant-phenotype associations that are supported using one assay but not another.

**4. Distinguish attribution for the experimental evidence from attribution for the claim**

Often the agent performing supporting experiments or writing the supporting paper is not the one who generated the association. The model needs to support attributing these separately - to  track agents responsible for different processes that create evidence for the association, including the annotation process itself, studies creating primary data, statistical analyses creating derived data, publications, assignment of evidence codes.

**5. Transitive provenance**

Tracking the lineage of evidence  supporting an association? (e.g. raw dataset1 -> derived dataset2 -> derived dataset3 -> association). We will need to evaluate how precisely we want to track contributors across previous artifacts (datasets, publications) that were used directly or indirectly in making a new claim (e.g.  datasets re-analyzed in a new context, or that informed experiments that created evidence for an association).

**6. Address core competency questions**

Competency questions represent diverse and valuable queries that users might ask of the data.  Our model should support all such queries in as efficient and effective a way as possible. Below, we provide some example competency questions for BRCA variant pathogenicity calls that will inform our modeling.

**A. Questions about evidence and provenance**

**1. Given a variant in the BRCA gene:**

A. What is this variant’s pathogenicity?

B. Who/what organization classified this variant?

C. How was the variant classified? What are the lines of evidence?

D. Given the lines of evidence, what is the level of review of this evidence?

**2. Given a variant and classification**

A. What sources agree/disagree with this classification?

**3. Return all variants with more than one source where two or more sources have a different classification of pathogenicity**

A. Would require mapping classifications to each other  (ACMG vs. ENIGMA)

**4. Given a variant, classification, and evidence type**

A. Is this evidence type compatible with other evidence models?  In other words, can we combine this new evidence with other models?  If so, does this change the classification of the variant?

**B. Questions about clinical phenotypes**

**1. Find all variants found in patients with ‘infiltrating carcinoma of breast’**

**2. Find overlapping phenotypes in patients with brca1 variant of type x, in a specific exon.**

A. Overlapping Phenotypes - phenotypes common in every patient with the input brca   variation

B. Variant type - missense variant, frameshift variant

C. Exon number, reference transcript

**3. Find variants in studied cell lines and/or patients with reactivity to drug X**

## III. Candidate Model for Describing Evidence and Provenance

In defining a first iteration model to capture evidence and provenance around BRCA variant-disease associations, we extend the model we have previously defined around genotype-phenotype associations for the Monarch Initiative (Appendix 1)  to accommodate evidence and provenance information inherent in the ACMG variant classification guidelines [2].  

We view each ACMG criteria as representing a different class of ‘evidence’ from the Evidence Code Ontology [3], and capture any additional information about how the evidence was produced as its provenance. Most variant associations have evidence based on more than one criteria, resulting in a set of evidence instances that when combined represents the total set of evidence for a given pathogenicity call.  This set of different types of evidence is what the Evidence Code Ontology terms 'combinatorial evidence'.  Notably, another component of this 'combinatorial evidence' set is the necessary association between the gene containing the variant, and the disease being investigated. This is because a causative role for the gene in the disease is a prerequisite for other criteria being evidence for pathogenicity of the variant  - for example, a prediction of the variant causing loss of gene function (criteria PVS1) is evidence for its pathogenicity only if there is also evidence that the gene plays a mechanistic role in the disease.

To represent this state of affairs, the model we define organizes all of these evidence instances around a central evidence node that represents an ECO:'combinatorial evidence' (Figure 2). This combinatorial evidence has as parts each atomic piece of evidence corresponding to a pathogenicity criteria from the classification guidelines, which can be linked to any techniques or data items related to the production of this evidence.  Attribution can be made for entities at any level of this graph, including the individual evidences, the combinatorial evidence, and the association representing the pathogenicity call.

![BRCA Model](https://raw.githubusercontent.com/monarch-initiative/mckb/master/docs/image/brca_evidence_model.jpg)
![cmap key](https://raw.githubusercontent.com/monarch-initiative/mckb/master/docs/image/cmap_key.png)

**Figure 2**. Model to describe evidence and provenance of pathogenicity classifications for cancer variants. Nodes in green represent evidence information.  Nodes in darker blue represent represent provenance-related    information, with nodes in lighter blue representing attribution components of provenance (agents responsible for a given process or evidence). The relations/properties used are not final specifications. The aim initially is to define the structure of the model.  The exact properties to be used will be defined once the structure of the model is settled.

## IV. Collaborators

The following collaborators have contributed to this use case document:

**Matthew Brush, PhD**: Ontologist, OHSU

**Melissa Haendel, PhD**: Director of the Ontology Development Group, OHSU

**Kent Shefchek, MS**: Bioinformatics Software Engineer, OHSU

**Mark Diekhans, PhD**: Technical Project Manager, UCSC

**Mary Goldman, BS**: Design and Usability Engineer, UCSC

**Benedict Paten, PhD**: Assistant Research Scientist, Comparative Genomics, UCSC

**Melissa Cline, PhD**:** **Associate Project Scientist, UCSC

**Pascale Gaudet, PhD**: neXtProt Scientific Manager, SIB

**Rebecca Jacobson MD, MS**: Professor of Biomedical Informatics, UPITT

**Harry Hochheiser, PhD**: Assistant Professor, UPITT

**Chris Mungall, PhD**: Bioinformatics Scientist, LBNL

**Nicole Washington, PhD**:  Computational Biologist Research Scientist, LBNL

## V. References

1. Statistics from cdc.gov/cancer/overian/statistics and cdc.gove.cancer/breast/statistics, numbers for 2011 (latest reported available)

2. Richards, Sue, et al. "Standards and guidelines for the interpretation of sequence variants: a joint consensus recommendation of the American College of Medical Genetics and Genomics and the Association for Molecular Pathology." *Genetics in Medicine* (2015).

3. http://www.evidenceontology.org/

**Appendix I: Candidate Model for Monarch Evidence and Provenance**

![Monarch Model](https://github.com/monarch-initiative/mckb/blob/master/docs/image/monarch_evidence_model.jpg)

