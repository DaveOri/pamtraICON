#!/bin/bash

source /home/dori/.bashrc

ICON_PATH='/data/inscape/icon/experiments/juelich/testbed/testbed_'
DATA_PATH='/data/optimice/pamtra_runs/tripex-pol/data/'
CODE_PATH='/home/dori/develop/pamtraICON/tripex-pol/'

FIRST_DAY='20181101'
TODAY=`date +%Y%m%d`

declare -a hydro_combo=("all_hydro" "no_snow" "only_snow" "only_liquid" "only_ice")
declare -a radar_names=("Joyrad10" "Joyrad35" "Grarad94")

newdata=0

DAY=$FIRST_DAY
echo $DAY $TODAY
until [[ ${DAY} > ${TODAY} ]]; do
	# Check if there is ICON output
	if [ -f ${ICON_PATH}${DAY}/METEOGRAM_patch001_${DAY}_joyce.nc ]; then
		echo ${DAY}
		for hydro in "${hydro_combo[@]}"; do
			for radar in "${radar_names[@]}"; do
				if [ -f ${DATA_PATH}${DAY}${hydro}_mom_${radar}.nc  ]; then
					echo "Already processed " ${DAY} ${hydro} ${radar}
				else
					echo "Running "${DAY} ${hydro} ${radar}
					python ${CODE_PATH}pamtra_radar.py --date ${DAY} --hydroset ${hydro} --radarset ${radar} > ${CODE_PATH}pamtra${DAY}_${hydro}_${radar}.out
					newdata=1
				fi
			done
			if [ "$newdata" -eq "1" ]; then
				python ${CODE_PATH}plot_results_separated.py --date ${DAY} --hydroset ${hydro} > ${CODE_PATH}plot${DAY}_${hydro}.out
				newdata=0
			fi
		done
	else
		echo "no ICON data for "${DAY}
	fi

	DAY=$(date -d "$DAY + 1 day" +%Y%m%d)
done
