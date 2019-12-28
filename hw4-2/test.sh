case=$1

make

srun -n 1 -p pp --gres=gpu:2  ./hw4-2 cases/${case} ${case}.out
#time ./seq cases/${case} ${case}.out
#time ./hw4-2 cases/${case} ${case}.out
#time ./hw4-1 cases/${case} ${case}.out

diff ${case}.out cases/${case}.out
