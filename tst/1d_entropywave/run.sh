#Create grid and run calculation

# - Run dnami with actual nozzle case
if [[ "$1" == "recomp"  ]]; then  
	./copy.sh
	cd ../../src
	./install_clean.sh
	cd ../wrk/
else
	./copy.sh
	cd ../../wrk/
fi
mpirun -np 1  python3 ./compute.py
