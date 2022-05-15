#!/bin/bash -l

TCLB=/group/pawsey0267/uqtmitc3/upstream_dev/CLB
MODEL=${1}
RUNFILE=${2}
SIMPLEFILE=${3}

sbatch <<EOT
#!/bin/bash -l
#SBATCH --account=pawsey0267
#SBATCH --job-name=${6}	
#SBATCH --partition=gpuq
#SBATCH --nodes=1
#SBATCH --gres=gpu:${4}
#SBATCH --ntasks-per-node=12

#SBATCH --output=logs/out.%6j
#SBATCH --error=logs/out.%6j

#SBATCH --time=4:00:00

export LD_LIBRARY_PATH=/pawsey/centos7.6/devel/cascadelake/gcc/8.3.0/openmpi-ucx/4.0.2/hdf5-parallel/1.10.5/lib:/pawsey/centos7.6/devel/cascadelake/gcc/8.3.0/openmpi-ucx/4.0.2/lib:/pawsey/centos7.6/devel/gcc/4.8.5/ucx/1.6.0/lib:/pawsey/centos7.6/devel/binary/cuda/10.2/extras/CUPTI/lib64:/pawsey/centos7.6/devel/binary/cuda/10.2/lib64:/pawsey/centos7.6/devel/gcc/4.8.5/gcc/8.3.0/lib64:/home/uqtmitc3/bin/lib:/pawsey/centos7.6/apps/cascadelake/gcc/8.3.0/r/3.6.3/lib64/R/lib
export PMI_NO_FORK=1

module load cuda/10.1 r/3.6.3 openmpi-ucx-gpu/4.0.2

## Create MPMD configuration file
cat > mpmd_conf/mpmd_\${SLURM_JOB_ID}.conf <<EOF
0-$(expr ${4} - 1) ${TCLB}/${1}/main ${RUNFILE}
${4}-$(expr ${5} - 1) ${TCLB}/${1}/simplepart ${SIMPLEFILE}
EOF

srun -l --cpu_bind=verbose,cores --multi-prog mpmd_conf/mpmd_\${SLURM_JOB_ID}.conf



#srun -n ${4} ${TCLB}/${1}/main ${RUNFILE} : -n $(expr ${5} - ${4}) ${TCLB}/${1}/simplepart ${SIMPLEFILE}

wait

EOT
