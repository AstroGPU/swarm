# -*- coding: utf-8 *-*
from logdb import IndexedLogDB
from bsddb3.db import *
from log import LogRecord
from struct import pack, unpack


def dump_a_file(fileName, start, end):
    d = IndexedLogDB(fileName)

    c = d.Cursor(d.primary.cursor())

    for i in range(0,end):
        l = c.next()
        if(l == None):
            break;
        elif(i >= start):
            print l,"\n"


from sys import argv
if(len(argv) < 4):
    print("Usage dump.py baseDbName start end")
else:
    dump_a_file(argv[1], int(argv[2]), int(argv[3]))

