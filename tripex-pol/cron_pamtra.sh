#!/bin/bash

source /home/dori/.bashrc

ICON_PATH='/data/inscape/icon/experiments/juelich/testbed/testbed_'
DATA_PATH='/data/optimice/pamtra_runs/tripex-pol/data/'
PLOT_PATH='/data/optimice/pamtra_runs/tripex-pol/plots/'
CODE_PATH='/net/ora/develop/pamtraICON/tripex-pol/'

FIRST_DAY='20181101'
TODAY=`date +%Y%m%d`

declare -a hydro_combo=("all_hydro" "no_snow" "only_snow" "only_liquid" "only_ice" "only_graupel_hail")
declare -a radar_names=("Joyrad10" "Joyrad35" "Grarad94")

newactive=0
newpassive=1

DAY=$FIRST_DAY
echo $DAY $TODAY
until [[ ${DAY} > ${TODAY} ]]; do
	# Check if there is ICON output
	if [ -f ${ICON_PATH}${DAY}/METEOGRAM_patch001_${DAY}_joyce.nc ]; then
		echo ${DAY}
		if [ -f ${DATA_PATH}${DAY}hatpro.nc ]; then
			echo "passive "${DAY}" already done"
		else
			python run_pamtra_hatpro.py --date ${DAY} > ${CODE_PATH}pamtra${DAY}_hatpro.out
			newpassive=1
		fi
		if [ "$newpassive" -eq "1" ]; then
			echo "New passive data ... plotting"
			python ${CODE_PATH}plot_pamtra_hatpro.py --date ${DAY} >> ${CODE_PATH}plot${DAY}_hatpro.out
			newpassive=1
		fi
		if [ -f ${PLOT_PATH}${DAY}hatpro.png ]; then
			echo "Already plotted passive " ${DAY}
		else
			echo "Found no passive plot ... plotting"
			python ${CODE_PATH}plot_pamtra_hatpro.py --date ${DAY} >> ${CODE_PATH}plot${DAY}_hatpro.out
		fi
		for hydro in "${hydro_combo[@]}"; do
			for radar in "${radar_names[@]}"; do
				if [ -f ${DATA_PATH}${DAY}${hydro}_mom_${radar}.nc  ]; then
					echo "Already processed " ${DAY} ${hydro} ${radar}
				else
					echo "Running "${DAY} ${hydro} ${radar}
					python ${CODE_PATH}pamtra_radar.py --date ${DAY} --hydroset ${hydro} --radarset ${radar} > ${CODE_PATH}pamtra${DAY}_${hydro}_${radar}.out
					newactive=1
				fi
			done
			if [ "$newactive" -eq "1" ]; then
				echo "Newdata active ... plotting"
				python ${CODE_PATH}plot_results_separated.py --date ${DAY} --hydroset ${hydro} >> ${CODE_PATH}plot${DAY}_${hydro}.out
				newactive=0
			fi
			if [ -f ${PLOT_PATH}${DAY}${hydro}_Ze.png ]; then
				echo "Already plotted " ${DAY} ${hydro}
			else
				echo "found no plot ... plotting"
				python ${CODE_PATH}plot_results_separated.py --date ${DAY} --hydroset ${hydro} >> ${CODE_PATH}plot${DAY}_${hydro}.out
			fi
		done
	else
		echo "no ICON data for "${DAY}
	fi

	DAY=$(date -d "$DAY + 1 day" +%Y%m%d)
done
