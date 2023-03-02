# Copy information to correct location
dst=../../wrk_swe
mkdir $dst 
mkdir $dst/bathy 

# Copy compute etc to wrk directory
cp compute.py $dst/
cp misc.py $dst/
cp values.py $dst/
cp plot.py $dst/

# Copy equations and MPI modified file
cp dnami_mpi.py ../../src/src_py
cp equations.py ../../src/generate
cp genRhs.py ../../src/generate
cp functions.py ../../src/generate
