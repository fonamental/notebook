#!/usr/bin/env python

from __future__ import print_function

import argparse


def main():
    parser = argparse.ArgumentParser(
        description='Convert multiple .seq files into a multi-sequence '
                    ' .fas file.')

    parser.add_argument('seqs', metavar="SEQ", nargs='+',
                        help='.seq file containing a single line')

    args = parser.parse_args()
    print(seq_to_fas(args.seqs), end='')


def seq_to_fas(seq_paths):
    fas_text = ''

    for seq_path in seq_paths:
        fas_text += '>' + seq_path + '\n'

        with open(seq_path, "r") as seq_file:
            fas_text += seq_file.read().strip() + '\n'

    return fas_text


if __name__ == "__main__":
    main()
