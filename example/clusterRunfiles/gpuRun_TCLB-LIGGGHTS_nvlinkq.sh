#!/bin/bash -l

TCLB=/group/pawsey0267/uqtmitc3/upstream_dev/CLB
MODEL=${1}
RUNFILE=${2}
SIMPLEFILE=${3}

sbatch <<EOT
#!/bin/bash -l
#SBATCH --account=director2188
#SBATCH --job-name=${6}	
#SBATCH --partition=nvlinkq
#SBATCH --nodes=1
#SBATCH --gres=gpu:${4}
#SBATCH --ntasks-per-node=28

#SBATCH --time=24:00:00

export LD_LIBRARY_PATH=/pawsey/centos7.6/devel/cascadelake/gcc/8.3.0/openmpi-ucx/4.0.2/hdf5-parallel/1.10.5/lib:/pawsey/centos7.6/devel/cascadelake/gcc/8.3.0/openmpi-ucx/4.0.2/lib:/pawsey/centos7.6/devel/gcc/4.8.5/ucx/1.6.0/lib:/pawsey/centos7.6/devel/binary/cuda/10.2/extras/CUPTI/lib64:/pawsey/centos7.6/devel/binary/cuda/10.2/lib64:/pawsey/centos7.6/devel/gcc/4.8.5/gcc/8.3.0/lib64:/home/uqtmitc3/bin/lib:/pawsey/centos7.6/apps/cascadelake/gcc/8.3.0/r/3.6.3/lib64/R/lib

export PMI_NO_FORK=1
module load cuda/10.2 r/3.6.3 openmpi-ucx-gpu/4.0.2
nvidia-smi

MPMDF=mpmd_conf/mpmd_\${SLURM_JOB_ID}.conf
cat > \${MPMDF} <<EOF
0-$(expr ${4} - 1) ${TCLB}/${1}/main ${RUNFILE}
${4}-$(expr ${4} + ${5} - 1) ${TCLB}/${1}/simplepart ${SIMPLEFILE}
EOF

srun -u -N 1 -n $(expr ${4} + ${5}) --mem=0 --gres=gpu:${4} --exclusive -l --cpu_bind=verbose,cores --multi-prog \${MPMDF} >tclb.\${SLURM_JOB_ID}.out 2>&1 &

wait

EOT
