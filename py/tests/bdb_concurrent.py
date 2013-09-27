#!/usr/bin/env python2
# -*- coding: utf8 -*-

## @file bdb_concurrent.py Testing that concurrent reads from a BDB log file works
# 
# swarm command is run to generate a log file and it is scheduled to remove it
# immedaitely after the integration is done.
#
# this test only passes if the query is quick enough to open the file and get the result.
#
#
# 

import swarmng
import sys
import os
import time
import psutil
import swarmng.logdb
from swarmng.range_type import Range
from swarmng.query import truncate
from argparse import Namespace
import unittest
from bsddb3.db import DBError

class BDBConcurrencyTest(unittest.TestCase):
  def runTest(self):
    ## Setting up the output log file
    output_file_name='Testing/testing_log.db'
    final_time = 100
    
    try:
      os.remove(output_file_name)
    except OSError:
      pass

    def do_query():
      d = swarmng.logdb.IndexedLogDB(output_file_name)
      q = list(truncate(1,d.final_conditions(Range.single(0))))
      self.assertEqual(len(q), 1)
      key, record = q[0]
      return record.time
      

    def do_integration():
      cfg = swarmng.config(
	nsys=64,
	nbod=3,
	log_writer='bdb',
	log_output_db=output_file_name,
	log_interval=1,
	integrator='hermite_cpu_log',
	time_step=0.01,
	nogpu=1
      );
      ref = swarmng.generate_ensemble( cfg )
      ens = ref.clone()

      swarmng.init(cfg)

      integ = swarmng.Integrator.create( cfg )

      integ.ensemble = ens
    
      swarmng.sync
      for i in range(1,final_time+1):
	integ.destination_time = i
	integ.integrate()
	

    pid = os.fork()
    if pid == 0:
      do_integration()
      os.abort()
    else:
      child = psutil.Process(pid)
      
      time.sleep(.05)
      counter = 0
      while child.status != psutil.STATUS_ZOMBIE:
	try:
	    # Get the last_time for system 0, it
	    # should be greater than 0
	    last_time = do_query()
	    self.assertGreater(last_time,0)
	    
	    # Only increment the counter if it is not at
	    # the final time
	    if(last_time < final_time):
	      counter += 1
	except DBError as e:
	    print("Querying failed", e)
	time.sleep(.05)
	
      # The query must have run before the integration
      # over at least once
      assert(counter > 0)
      
      # Wait for the child process to finish
      os.waitpid(pid, 0)
