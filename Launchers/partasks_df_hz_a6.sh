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
#PBS -lwalltime=800:00:00
#PBS -lpmem=2000MB

maxtasks=10

##################################################################
# Set the simulation parameters
##################################################################
commands=( 'critical_bilayer(1, 0.1, 0.2, [ 0 0 10 ], 6, pi/4)'
           'critical_bilayer(1, 0.2, 0.2, [ 0 0 10 ], 6, pi/4)'
           'critical_bilayer(1, 0.3, 0.2, [ 0 0 10 ], 6, pi/4)'
           'critical_bilayer(1, 0.4, 0.2, [ 0 0 10 ], 6, pi/4)'
           'critical_bilayer(1, 0.5, 0.2, [ 0 0 10 ], 6, pi/4)'
           'critical_bilayer(1, 0.6, 0.2, [ 0 0 10 ], 6, pi/4)'
           'critical_bilayer(1, 0.7, 0.2, [ 0 0 10 ], 6, pi/4)'
           'critical_bilayer(1, 0.8, 0.2, [ 0 0 10 ], 6, pi/4)'
           'critical_bilayer(1, 0.9, 0.2, [ 0 0 10 ], 6, pi/4)'
           'critical_bilayer(1, 1.0, 0.2, [ 0 0 10 ], 6, pi/4)' )

##################################################################
# Perform the simulation
##################################################################

# Launch the simulation script once for each parameter defined above
module load matlab
for n in $(seq 0 9); do
	# Make a work directory for Matlab
    	workdir="/work/$LOGNAME/matlab/df_hz_a6/${n}"
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
