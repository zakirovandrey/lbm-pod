import shutil
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--drop_dir', dest='drop_into', default="./drop/", help='drop directory')
parser.add_argument('--nocomp', dest='iscompile', action='store_false', help='not recompile')
parser.add_argument('--norun', dest='norun', action='store_true', help='do not make clean')
parser.add_argument('--redefine', dest='newNxyz', default="", help='redefine Nx Ny Nz dr and dt')
arguments, unknown_args = parser.parse_known_args()
drop_into = arguments.drop_into
python_site_packages=['/s/ls4/users/zakirov/.local/lib/python2.6/site-packages',]
for sp in python_site_packages: sys.path.append(sp)
if not os.path.exists(drop_into): os.makedirs(drop_into)
shutil.copy2(__file__, drop_into)

constdefs=[]
if FloatPrecision==1: constdefs.append("USE_FLOAT=1")
if FloatPrecision==2: constdefs.append("USE_DOUBLE=1")

constdefs=" ".join(constdefs)

if arguments.iscompile:
  retval = os.system("SET_FROM_PYTHON=1 NX=%d NY=%d NZ=%d %s make -j4"%(Nx,Ny,Nz,constdefs))
  if retval: exit(retval)
if arguments.norun: exit(0)
if MPIon==1: comm.Barrier()

