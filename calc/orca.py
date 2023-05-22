import os

from ase.atoms import Atoms
from ase.io import read


class OrcaNEB:
    batch = """#!/bin/bash
ulimit -s unlimited
module purge
module load orca
export job=orca
export RSH_COMMAND="/usr/bin/ssh -x"
export scratchlocation=/scratch/
tdir=$(mktemp -d $scratchlocation/$USER/orcajob__$SLURM_JOB_ID-XXXX)
cp  $SLURM_SUBMIT_DIR/*.inp $tdir/

cp  $SLURM_SUBMIT_DIR/*.gbw $tdir/

cp  $SLURM_SUBMIT_DIR/*.xyz $tdir/
cd $tdir
/usr/public/orca/orca_5_0_1_linux_x86-64_shared_openmpi411///orca $job.inp >  $SLURM_SUBMIT_DIR/$job.out
cp $tdir/* $SLURM_SUBMIT_DIR
"""

    def __init__(self, method, label='orca_temp', srun_command=None):
        self.inp = f"""!{method} NEB-TS
%PAL NPROCS 64 END
%NEB NEB_END_XYZFILE "init.xyz"
preopt true
NImages 8
END
* XYZfile 0 1 final.xyz
"""
        self.label = label
        self.srun = srun_command or os.getenv('ORCA_CMD') or 'srun -n 64 ./orca.sh'

    def run_neb(self, r: Atoms, p: Atoms):
        # make directory if not exist
        if not os.path.exists(self.label):
            os.makedirs(self.label)
        # write input file
        with open(f'{self.label}/orca.inp', 'w') as f:
            f.write(self.inp)
        # write bash file and make executable
        with open(f'{self.label}/orca.sh', 'w') as f:
            f.write(self.batch)
        os.system(f'chmod +x {self.label}/orca.sh')
        # write initial and final xyz files
        r.write(f'{self.label}/init.xyz')
        p.write(f'{self.label}/final.xyz')
        # run orca
        os.system(f'cd {self.label}; {self.srun}')

    def get_ts(self):
        return read(f'{self.label}/NEB-TS.xyz')

    def get_reactant(self):
        return read(f'{self.label}/orca_reactant.xyz')

    def get_product(self):
        return read(f'{self.label}/orca_product.xyz')