from query import *
from range_type import *
from argparse import Namespace

args = Namespace()
args.system_range = Range.universal()
args.time_range = Range.universal()
args.body_range = Range.single(1)
args.evt_id = Range.single(1)
args.max_records = 1
args.database = "/scratch/hpc/salehqt/hg/swarm-build/d23.db"
args.keplerian = None
args.initial_conditions = False
args.final_conditions = False

run_with_args(args)



