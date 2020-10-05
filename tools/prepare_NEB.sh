#!/usr/bin/env bash
read_E() {
  grep -i entropy OUTCAR >tmp1
  tail -1 tmp1 >tmp2
  E=$(awk '{print $7}' tmp2)
  rm tmp1 tmp2
}
main() {
  max=$(ls | sed 's/e_//' | sed 's/s//' | sort -n | tail -1)
  for i in $(seq 0 "${max}"); do
    mkdir -p NEB_$i/00 NEB_$i/06
    cp s/CONTCAR NEB_$i/POSCAR0
    cp e_$i/CONTCAR NEB_$i/POSCAR1
    cp s/OUTCAR NEB_$i/00/OUTCAR
    cp s/OUTCAR NEB_$i/OUTCAR0
    cp e_$i/OUTCAR NEB_$i/06/OUTCAR
    cp e_$i/OUTCAR NEB_$i/OUTCAR1
    cp s/KPOINTS s/POTCAR NEB_$i/
    cd NEB_$i
    nebmake.pl POSCAR0 POSCAR1 5
    cd 00
    read_E
    E1=$E
    cd ..
    cd 06
    read_E
    E2=$E
    cd ..
    cat <<EOF >INCAR
NWRITE = 2

Electronic Relaxation 1:
PREC   = Acc
ISYM   = 2
NELM =  240
NELMIN = 4

Ionic Relaxation:
NSW    = 10000
NBLOCK = 1
KBLOCK = 1
IBRION = 3
POTIM  = 0
IOPT = 3
ISIF   = 2

DOS Related Values:
ISMEAR = 1
SIGMA  = 0.4

Electronic Relaxation 2:
IALGO  = 48
LREAL  = AUTO
ENCUT  = 450
ENAUG  = 600.0

LCLIMB = .TRUE.
ICHAIN = 0
IMAGES = 5
SPRING = -5

LWAVE = .FALSE.
EFIRST =  lala
ELAST =   lblb

EDIFFG = -0.05
EOF

    sed -i s/lala/${E1}/g INCAR
    sed -i s/lblb/${E2}/g INCAR
    cat <<EOF >>va_neb.stampede2
#!/bin/bash
#SBATCH -J vasp
#SBATCH -o vasp_neb.%j.out
#SBATCH -e vasp_neb.%j.err
#SBATCH -n 60
#SBATCH -N 2
#SBATCH -p normal
#SBATCH -t 36:00:00
#SBATCH -A  TG-DMR190035

#TG-MSS160003 TG-DMR190035

module load vasp/5.4.4

run_mpirun()
{
    ibrun -n 60 vasp_std_vtst > vasp_test.out
}

main()
{
run_mpirun
clear_vasp_out
perl ~/goali/vtstscripts/nebbarrier.pl
perl ~/goali/vtstscripts/nebspline.pl
gnuplot ~/goali/vtstscripts/nebplot.gnu
perl ~/goali/vtstscripts/vfin.pl stored
}
time main;
EOF
    sbatch va_neb.stampede2
    cd ..
  done
}
time main