#!/usr/bin/env python3

from setuptools import setup, find_packages

setup(name='mckb',
      version='0.0.1',
      packages=find_packages(),
      install_requires=['dipper', 'pymysql', 'rdflib'],
      include_package_data=True,
      dependency_links=['https://github.com/monarch-initiative/dipper/tarball/master#egg=Dipper-0.0.1']
      )