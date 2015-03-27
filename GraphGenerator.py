#!/usr/bin/env python3

import argparse
import logging


def main():
    logger = logging.getLogger(__name__)

    parser = argparse.ArgumentParser(description='Cancer Knowledge Base Graph'
                                                 ' Generator',
                                     formatter_class=
                                     argparse.RawTextHelpFormatter)
    parser.add_argument('--host', '-H', type=str, default="localhost:8080",
                        help='Location of MySQL Server')
    parser.add_argument('--database', '-D', type=str, required=True,
                        help='Name of database')
    parser.add_argument('--user', '-u', help='Username')
    parser.add_argument('--password', '-p', help='Password')

    args = parser.parse_args()
    return

if __name__ == "__main__":
    main()
