#!/usr/bin/env python3
from mckb.sources.CGD import CGD
import argparse
import logging
import getpass
import sys
import json


def main():
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)

    parser = argparse.ArgumentParser(description='Cancer Knowledge Base Graph'
                                                 ' Generator',
                                     formatter_class=
                                     argparse.RawTextHelpFormatter)
    parser.add_argument('--host', '-H', type=str, default="localhost",
                        help='Location of MySQL Server')
    parser.add_argument('--database', '-D', type=str, help='Name of database')
    parser.add_argument('--user', '-u', help='Username')
    parser.add_argument('--password', '-p', help='Password')
    parser.add_argument('--config', '-c', help='Config file, see example '
                                               'formatting in conf directory')

    args = parser.parse_args()

    # Config file overrides command line credentials
    # We need to refactor the Dipper config.py so it is reusable here
    if args.config is not None:
        credentials = json.load(open(args.config, 'r'))
        args.host = credentials['dbauth']['cgd']['host']
        args.database = credentials['dbauth']['cgd']['database']
        args.user = credentials['dbauth']['cgd']['user']
        args.password = credentials['dbauth']['cgd']['password']

    if args.password is None:
        if sys.stdin.isatty():
            args.password = getpass.getpass(prompt="Enter your password: ",
                                            stream=None)
        else:
            args.password = input("Enter your password: ")

    # Parse test source
    cgd = CGD(args.database, args.user, args.password, args.host)
    #cgd.fetch(False)
    cgd.parse()
    cgd.write(format='turtle')
    return

if __name__ == "__main__":
    main()
