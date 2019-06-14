export OMP_NUM_THREADS=4
export MPICH_MAX_THREAD_SAFETY=multiple
nblock=$1
mem_size=18
sites=8
verb=4
#VALGRIND=valgrind
LOG=mpirun.log."$nblock"
NP=16
GEOM='2 2 2 2'
echo mpirun -n $NP $VALGRIND ./BlockScramble.x $sites $mem_size $nblock $verb $GEOM 2>&1 |tee  $LOG
mpirun -n $NP $VALGRIND ./BlockScramble.x $sites $mem_size $nblock $verb $GEOM 2>&1 |tee  $LOG
#mpirun -n 16 $VALGRIND ./BlockScramble.x $sites $mem_size $nblock --mpi 2.2.2.2 2>&1 |tee  $LOG
grep '\-1' $LOG
#grep send_before mpirun.log |sort
#grep send_after mpirun.log |sort
#grep recv_buf mpirun.log |sort
