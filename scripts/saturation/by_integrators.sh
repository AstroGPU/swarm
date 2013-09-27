DTIME=2.0 # Destination time
NBOD="3 4 6 8 10" # a list of number of bodies
INTEGRATORS="Hermite" # a list of integrators
SOURCE=`dirname $0`/../../
NSYS=`python $SOURCE/scripts/saturation/ranges.py` # auto generated list of nsys range

[ -z "$BUILD" ] && BUILD=.

# Guard against missing arguments
if [ $# -lt 1 ]
then
    echo "Syntax: $0 <output directory>"
    exit -1
fi

OUTPUTDIR="$1"
mkdir -p $OUTPUTDIR

for I in $INTEGRATORS
do
for N in $NBOD
do
    echo "Starting $N bodies, $I --------------------------------"
    mkdir -p $OUTPUTDIR/$I
    O=$OUTPUTDIR/$I/$N.csv
    $BUILD/bin/swarm benchmark -c $SOURCE/samples/$I.cfg --range nsys=$NSYS verbose=1 \
        nsys=2000 nbod=$N -d $DTIME verbose=1 -v 100 2> $OUTPUTDIR/$I/$N-log.txt | tail -n +10 > $O 
done
done
