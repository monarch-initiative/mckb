#!/usr/bin/env python3
from mckb.sources.G2P import G2P
import argparse
import logging


def main():
    logging.basicConfig(level=logging.DEBUG)
    logger = logging.getLogger(__name__)

    parser = argparse.ArgumentParser(description='Cancer Knowledge Base Graph'
                                                 ' Generator',
                                     formatter_class=
                                     argparse.RawTextHelpFormatter)
    parser.add_argument('--host', '-H', type=str, default="localhost",
                        help='Location of MySQL Server')
    parser.add_argument('--database', '-D', type=str, required=True,
                        help='Name of database')
    parser.add_argument('--user', '-u', help='Username')
    parser.add_argument('--password', '-p', help='Password')

    args = parser.parse_args()

    test_source = G2P(args.database, args.user, args.password, args.host)
    test_source.parse()
    test_source.disconnect_from_database()
    return

if __name__ == "__main__":
    main()