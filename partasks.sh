#!/bin/bash
# PBS job script for launching the 'critical_bilayer' simulation.
#
# Written by Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
# Based on an example script by Roy Dragseth and Egil Holsvik.

##################################################################
# Configure the PBS queue system
##################################################################
#PBS -M jabirali@switzerlandmail.ch
#PBS -A acc-fysikk
#PBS -q default	
#PBS -lnodes=1:ppn=12
#PBS -lwalltime=600:00:00
#PBS -lpmem=2000MB

maxtasks=12



##################################################################
# Set the simulation parameters
##################################################################
exchange=(  0.0 1.0 3.0 6.0 )
spinorbit=( 0.0 0.5 1.0 2.0 )
width=( 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.75 1.00 )



##################################################################
# Perform the simulation
##################################################################

# Launch the simulation script for each parameter defined above
module load matlab
for h in ${exchange[@]}; do
	for a in ${spinorbit[@]}; do
		for d in ${width[@]}; do
			# Make a work directory for Matlab
			workdir="/work/$LOGNAME/matlab/bilayer_${h}_${a}/${d}"
			mkdir -p "$workdir"
			cd "$workdir"

			# Link the required scripts into the work directory
			ln -st . "$PBS_O_WORKDIR/BVP" "$PBS_O_WORKDIR/Classes" "$PBS_O_WORKDIR/initialize.m" "$PBS_O_WORKDIR/critical_bilayer.m"

			# Execute the launch script
			echo $(uname -n) "Task starting: h=${h}, a=${a}, d=${d}."
			echo matlab -nodisplay -r "critical_bilayer(1,${d},0.2,[${h},0,0],${a},pi/4)" 1>output.txt 2>error.txt &
			echo $(uname -n) "Task complete: h=${h}, a=${a}, d=${d}."

			# Wait until activetasks < maxtasks before continuing the loop
  			activetasks=$(jobs | wc -l)
			while [ $activetasks -ge $maxtasks ]; do
				sleep 5
				activetasks=$(jobs | wc -l)
			done
		done
	done
done
module unload matlab

# Wait for all tasks to finish
echo "Waiting for tasks to complete"
wait
echo "All jobs complete!"
