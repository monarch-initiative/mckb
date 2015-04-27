[![Build Status](https://travis-ci.org/monarch-initiative/mckb.svg?branch=master)](https://travis-ci.org/monarch-initiative/mckb)

* Required python packages:
    * [Dipper](https://github.com/monarch-initiative/dipper)
    * [PyMySQL](https://github.com/PyMySQL/PyMySQL)
    * [RDFLib](https://github.com/RDFLib/rdflib)
    
##Quickstart

1. Load a MySQL database using the dump file in the resources directory
2. Create a configuration file in the conf directory using the example_conf.json as a template
3. Run:

        ./GraphGenerator.py --config conf/conf.json

4. This will create a directory called out in your working directory containing the output turtle files
    

