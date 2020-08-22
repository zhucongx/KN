getEBarrier() {
  Estart=$(sed '1,1!d' neb.dat | awk '{print $3}')
  Esaddle=$(awk 'BEGIN{a=   0}{if ($3>0+a) a=$3} END{print a}' neb.dat)
  Fsaddle=$(awk 'BEGIN{a=0; b=0}{if ($3>0+a) {a=$3; b=$4} } END{print b}' neb.dat)
  Eend=$(sed '7,7!d' neb.dat | awk '{print $3}')
  Eforward=$(echo "($Esaddle)-($Estart)" | bc)
  Ebackward=$(echo "($Esaddle)-($Eend)" | bc)
  Estartshow=$(grep EFIRST INCAR | awk '{print $3}')
  Eendshow=$(grep ELAST INCAR | awk '{print $3}')
  Dist=$(dist.pl POSCAR0 POSCAR1)
  DistAlong=$(sed '7,7!d' neb.dat | awk '{print $2}')
  Position=$(awk 'BEGIN{a=0; b=0}{if ($3>0+a) {a=$3; b=$2} } END{print b}' neb.dat)
}

main() {
  rm -rf ./E_neb_all_start_end.dat
  for j in $(seq 0 62); do
    cd config$j/
    max=$(ls | sed 's/NEB_//' | sort -n | tail -1)
    if [ ${max} = 's' ]; then
      echo "#" "neb was not run"
      cd ..
      continue
    fi
    for k in $(seq 0 "${max}"); do
      cd NEB_$k
      if [ ! -f "neb.dat" ]; then
        echo "#" $i $j $k "neb.dat not found"
        cd ..
        continue
      fi
      getEBarrier
      if (($(echo $Dist | awk '{if ($1 < 3.5 && $1 > 0.5) print 1;}'))); then
        read Dis0 Dis1 Dis2 Dis3 Dis4 Dis5 Dis6 Dis7 Dis8 Dis9 Dis10 <<<$(/Users/zhucongx/Program/goali/KN/bin/kn.exe)
        cd ../..
        echo $i $j $k $Eforward $Ebackward $Estartshow $Eendshow $Fsaddle $Dist $DistAlong $Position $Dis0 $Dis1 $Dis2 $Dis3 $Dis4 $Dis5 $Dis6 $Dis7 $Dis8 $Dis9 $Dis10
        echo $i $j $k $Eforward $Ebackward $Estartshow $Eendshow $Fsaddle $Dist $DistAlong $Position $Dis0 $Dis1 $Dis2 $Dis3 $Dis4 $Dis5 $Dis6 $Dis7 $Dis8 $Dis9 $Dis10 >>./E_neb_all_start_end.dat
      else
        cd ../..
        echo "#" $i $j $k "Dist is" $Dist
      fi
      cd config$j/
    done
    cd ..
  done
}

main
