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
#PBS -lwalltime=150:00:00
#PBS -lpmem=2000MB

maxtasks=12

##################################################################
# Set the simulation parameters
##################################################################
commands=( 'critical_bilayer(1, 0.15, 0.2, [ 10 0 0 ], 2, pi/4)'
           'critical_bilayer(1, 0.17, 0.2, [ 10 0 0 ], 2, pi/4)'
           'critical_bilayer(1, 0.19, 0.2, [ 10 0 0 ], 2, pi/4)'
           'critical_bilayer(1, 0.22, 0.2, [ 10 0 0 ], 2, pi/4)'
           'critical_bilayer(1, 0.25, 0.2, [ 10 0 0 ], 2, pi/4)'
           'critical_bilayer(1, 0.27, 0.2, [ 10 0 0 ], 2, pi/4)'
           'critical_bilayer(1, 0.32, 0.2, [ 10 0 0 ], 2, pi/4)'
           'critical_bilayer(1, 0.35, 0.2, [ 10 0 0 ], 2, pi/4)'
           'critical_bilayer(1, 0.37, 0.2, [ 10 0 0 ], 2, pi/4)'
           'critical_bilayer(1, 0.42, 0.2, [ 10 0 0 ], 2, pi/4)'
           'critical_bilayer(1, 0.45, 0.2, [ 10 0 0 ], 2, pi/4)'
           'critical_bilayer(1, 0.47, 0.2, [ 10 0 0 ], 2, pi/4)' )

##################################################################
# Perform the simulation
##################################################################

# Launch the simulation script once for each parameter defined above
module load matlab
for n in $(seq 0 11); do
	# Make a work directory for Matlab
    	workdir="/work/$LOGNAME/matlab/df_hx_a2_part2/${n}"
    	mkdir -p "$workdir"
    	cd "$workdir"

    	# Link the required scripts into the work directory
    	ln -st . "$PBS_O_WORKDIR/BVP" "$PBS_O_WORKDIR/Classes" "$PBS_O_WORKDIR/initialize.m" "$PBS_O_WORKDIR/critical_bilayer.m"

    	# Execute the simulation script in the background
    	echo $(uname -n) $(date +'%H:%M:%S') "Starting simulation:"
	echo "::  ${commands[n]}"
    	{
		echo "${commands[n]}";
    		matlab -nodisplay -r "${commands[n]}";
    	} 1>output.log 2>error.log &

    	# Wait until activetasks < maxtasks before continuing the loop
    	activetasks=$(jobs | wc -l)
    	while [ $activetasks -ge $maxtasks ]; do
    		sleep 5
    		activetasks=$(jobs | wc -l)
    	done
done
module unload matlab

# Wait for all tasks to finish
echo "Waiting for tasks to complete"
wait
echo "All jobs complete!"
