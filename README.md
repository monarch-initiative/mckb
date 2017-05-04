## Monarch Cancer Knowledge Base
MCKB is a pure python application to extract and transform clinically actionable cancer linked variants and metadata into a directed graph.

[![Build Status](https://travis-ci.org/monarch-initiative/mckb.svg?branch=master)](https://travis-ci.org/monarch-initiative/mckb)
[![Coverage Status](https://coveralls.io/repos/monarch-initiative/mckb/badge.svg)](https://coveralls.io/r/monarch-initiative/mckb)

### Motivation
Represent cancer data using controlled vocabularies (ontologies) and output as a directed graph serialized as RDF triples.  MCKB is currently a proof of concept to explore the benefits and challenges of mapping cancer data to available ontologies and storing both the output data and ontologies in a single datastore.  As a test set, MCKB is using a subset of a [physician curated dataset](https://www.synapse.org/#!Synapse:syn2370773/wiki/) by Dr. Rodrigo Dienstmann which was curated and transformed into a RDMS by the OHSU Clinical Genomics Database team.

#### Note: MCKB is in the pre-alpha stage of development.  Data models and any output should be considered unstable.  We will indicate here when models have stabilized.

### Applications

##### [SciGraph Implementation](https://github.com/monarch-initiative/mckb/tree/master/conf/SciGraph)
While the output data files can be stored in various databases, we also provide [configuration files] (https://github.com/monarch-initiative/mckb/tree/master/conf/SciGraph) to insert data into a Neo4J graph database using the [SciGraph](https://github.com/SciGraph/SciGraph) application.  These files also contain extensions to the default REST services that include queries specific to cancer use cases using the Cypher query language and SciGraph query expansion.

### Requirements
* MCKB requires Python 3 and the following packages:
    * [Dipper](https://github.com/monarch-initiative/dipper)
    * [PyMySQL](https://github.com/PyMySQL/PyMySQL)
    * [RDFLib](https://github.com/RDFLib/rdflib)
    
### Quickstart

1. Load a MySQL database using the dump file in the resources directory
2. Create a configuration file in the conf directory using the example_conf.json as a template
3. Run:

        ./GraphGenerator.py --config conf/conf.json

4. This will create a directory called out in your working directory containing the output turtle files

### Output
Example output can be found here:
https://raw.githubusercontent.com/monarch-initiative/mckb/master/ttl/cgd.ttl

