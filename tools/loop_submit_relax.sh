for i in `seq 0 0`
do
    cd config$i/
    for j in $(ls -d */)
    do
        cd ${j}
        sbatch submit_stampede2.sh
        cd ..
    done
    cd ..
done
