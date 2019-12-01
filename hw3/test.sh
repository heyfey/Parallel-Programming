case=$1

make

time srun -N1 -n1 -c12 ./hw3 cases/${case} ${case}.out

diff cases/${case}.out ${case}.out

# ./print_out cases/${case}.out
