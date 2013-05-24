#!/usr/bin/env python
# -*- coding: utf-8 *-*
from query import *
from range_type import *
import argparse

MAX_RECORDS = 10000;

def parse_cmd():
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--database", help="Database filename", required = True)
    parser.add_argument("-m", "--max-records", help="maximum number of records to process", default=MAX_RECORDS, type=int )
    parser.add_argument("-s", "--system-range", help="Range of systems to display",type=RangeType(int),default=Range.universal());
    parser.add_argument("-t", "--time-range", help="Range of time to display", type=RangeType(float), default=Range.universal())
    parser.add_argument("-b", "--body-range", help="Range of bodies to display", type=RangeType(int), default=Range.universal())
    parser.add_argument("-e", "--evt-id", help="The type of event to display (provide codes)",type=RangeType(int), default=Range.universal());
    parser.add_argument("-k", "--keplerian", help="Keplerian output. The argument defines the type of coordinate system used", choices = Keplerian.choices)
    parser.add_argument("--initial-conditions",  action="store_true",default=False)
    parser.add_argument("--final-conditions", action="store_true", default=False)

    args = parser.parse_args()
    return args


##################### MAIN #################################
if __name__ == '__main__':
    args = parse_cmd()
    print args
    run_with_args(args)