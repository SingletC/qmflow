import os

from ase.atoms import Atoms
from ase.io import read


class OrcaNEB:
    batch = """#!/bin/bash
#SBATCH -t 80:00:00
#SBATCH --mem 150G
#SBATCH -p exclusive
#SBATCH -n 48
#SBATCH -N 1
ulimit -s unlimited
module purge
module load orca
export job=orca
export RSH_COMMAND="/usr/bin/ssh -x"
export scratchlocation=/local/
export ORCA_SCRDIR=$scratchlocation
tdir=$(mktemp -d $scratchlocation//orcajob__$SLURM_JOB_ID-XXXX)
cp  $SLURM_SUBMIT_DIR/*.inp $tdir/

cp  $SLURM_SUBMIT_DIR/*.gbw $tdir/

cp  $SLURM_SUBMIT_DIR/*.xyz $tdir/
cd $tdir
command_path=$(which orca)
$command_path $job.inp > $job.out
cp $tdir/*.xyz $SLURM_SUBMIT_DIR
cp $tdir/*.out $SLURM_SUBMIT_DIR
rm $tdir/ -rf
    """

    def __init__(self, method, label='orca_temp', srun_command=None):
        self.inp = f"""!{method} NEB-TS
%PAL NPROCS 48 END
%maxcore 3000
%NEB NEB_END_XYZFILE "init.xyz"
preopt false
NImages 8
END
* XYZfile 0 1 final.xyz
"""
        self.label = label
        self.srun = srun_command or 'sbatch --wait orca.sh'

    def run_neb(self, r: Atoms, p: Atoms):
        # make directory if not exist
        if not os.path.exists(self.label):
            os.makedirs(self.label)
        # write input file
        with open(f'{self.label}/orca.inp', 'w') as f:
            f.write(self.inp)
        # write bash file
        with open(f'{self.label}/orca.sh', 'w') as f:
            f.write(self.batch)
        # write initial and final xyz files
        r.write(f'{self.label}/init.xyz')
        p.write(f'{self.label}/final.xyz')
        # run orca
        os.system(f'cd {self.label}; {self.srun}')

    def get_ts(self):
        return read(f'{self.label}/orca_NEB-TS_converged.xyz')

    def get_reactant(self):
        return read(f'{self.label}/init.xyz')

    def get_product(self):
        return read(f'{self.label}/final.xyz')