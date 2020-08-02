for i in `seq 0 0`
do
	cp prepare_NEB.sh config$i/
	cd config$i
	sh prepare_NEB.sh
	cd ..
done
