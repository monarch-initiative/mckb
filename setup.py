#!/usr/bin/env python3

from setuptools import setup, find_packages

setup(name='mckb',
      version='0.0.1',
      packages=find_packages(),
      install_requires=['dipper', 'pymysql'],
      include_package_data=True
      )