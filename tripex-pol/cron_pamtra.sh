#!/bin/bash

source /home/dori/.bashrc

ICON_PATH=/data/inscape/icon/experiments/juelich/testbed/testbed_
ROOT_PATH=/data/optimice/pamtra_runs/tripex-pol/
DATA_PATH=${ROOT_PATH}data/
PLOT_PATH=${ROOT_PATH}plots/
CODE_PATH=/net/ora/develop/pamtraICON/tripex-pol/

FIRST_DAY='20181101'
TODAY=`date +%Y%m%d`

declare -a hydro_combo=("all_hydro" "no_snow" "only_snow" "only_liquid" "only_ice" "only_graupel_hail")
declare -a radar_names=("Joyrad10" "Joyrad35" "Grarad94")

newactive=0
newpassive=0

DAY=$FIRST_DAY
echo $DAY $TODAY
until [[ ${DAY} > ${TODAY} ]]; do
	# Check if there is ICON output
	if [ -f ${ICON_PATH}${DAY}/METEOGRAM_patch001_${DAY}_joyce.nc ]; then
		echo ${DAY}
		if [ -f ${DATA_PATH}${DAY}hatpro.nc ]; then
			echo "passive "${DAY}" already done"
		else
			python2 ${CODE_PATH}pamtra_radar+hatpro.py --date ${DAY} --datapath ${DATA_PATH} --hydroset all_hydro --radarset hatpro > ${CODE_PATH}pamtra${DAY}_hatpro.out
			newpassive=1
		fi
		if [ "$newpassive" -eq "1" ]; then
			echo "New passive data ... plotting"
			python2 ${CODE_PATH}plot_pamtra_hatpro.py --date ${DAY} --rootpath ${ROOT_PATH} >> ${CODE_PATH}plot${DAY}_hatpro.out
			newpassive=0
		fi
		if [ -f ${PLOT_PATH}${DAY}hatpro.png ]; then
			echo "Already plotted passive " ${DAY}
		else
			echo "Found no passive plot ... plotting"
			python2 ${CODE_PATH}plot_pamtra_hatpro.py --date ${DAY} --rootpath ${ROOT_PATH} >> ${CODE_PATH}plot${DAY}_hatpro.out
		fi
		for hydro in "${hydro_combo[@]}"; do
			for radar in "${radar_names[@]}"; do
				if [ -f ${DATA_PATH}${hydro}/${DAY}${hydro}_mom_${radar}.nc  ]; then
					echo "Already processed " ${DAY} ${hydro} ${radar}
				else
					echo "Running "${DAY} ${hydro} ${radar}
					python2 ${CODE_PATH}pamtra_radar+hatpro.py --date ${DAY} --hydroset ${hydro} --radarset ${radar} --datapath ${DATA_PATH} > ${CODE_PATH}pamtra${DAY}_${hydro}_${radar}.out
					newactive=1
				fi
			done
			if [ "$newactive" -eq "1" ]; then
				echo "Newdata active ... plotting"
				python2 ${CODE_PATH}plot_results_separated.py --date ${DAY} --hydroset ${hydro} --rootpath ${ROOT_PATH} >> ${CODE_PATH}plot${DAY}_${hydro}.out
				newactive=0
			fi
			if [ -f ${PLOT_PATH}${hydro}/${DAY}${hydro}_Ze.png ]; then
				echo "Already plotted " ${DAY} ${hydro}
			else
				echo "found no plot ... plotting"
				python2 ${CODE_PATH}plot_results_separated.py --date ${DAY} --hydroset ${hydro} --rootpath ${ROOT_PATH} >> ${CODE_PATH}plot${DAY}_${hydro}.out
			fi
		done
	else
		echo "no ICON data for "${DAY}
	fi

	DAY=$(date -d "$DAY + 1 day" +%Y%m%d)
done

# NOW DO THE SAME FOR OLD TRIPEX

ICON_PATH=/data/inscape/icon/experiments/juelich/testbed/testbed_
ROOT_PATH=/data/optimice/pamtra_runs/tripex
DATA_PATH=${ROOT_PATH}/data
PLOT_PATH=${ROOT_PATH}/plots
CODE_PATH=/net/ora/develop/pamtraICON/tripex-pol

FIRST_DAY='20151111'
LAST_DAY='20160105'
#LAST_DAY='20151215'

declare -a hydro_combo=("all_hydro" "no_snow" "only_snow" "only_liquid" "only_ice" "only_graupel_hail")
#declare -a hydro_combo=("all_hydro")
declare -a radar_names=("KiXPol" "Joyrad35" "Joyrad94")

newactive=0
newpassive=0

DAY=$FIRST_DAY
echo $DAY $LAST_DAY
until [[ ${DAY} > ${LAST_DAY} ]]; do
	# Check if there is ICON output
	if [ -f ${ICON_PATH}${DAY}/METEOGRAM_patch001_${DAY}_joyce.nc ]; then
		echo ${DAY} ${DATA_PATH}/${DAY}hatpro.nc
		if [ -f ${DATA_PATH}/${DAY}hatpro.nc ]; then
			echo "passive "${DAY}" already done"
		else
			python2 ${CODE_PATH}/pamtra_radar+hatpro.py --date ${DAY} --datapath ${DATA_PATH}/ --hydroset all_hydro --radarset hatpro > ${CODE_PATH}/pamtra${DAY}_hatpro.out
			newpassive=1
		fi
		if [ "$newpassive" -eq "1" ]; then
			echo "New passive data ... plotting"
			python2 ${CODE_PATH}/plot_pamtra_hatpro.py --date ${DAY} --rootpath ${ROOT_PATH}/ >> ${CODE_PATH}/plot${DAY}_hatpro.out
			newpassive=0
		fi
		if [ -f ${PLOT_PATH}/${DAY}hatpro.png ]; then
			echo "Already plotted passive " ${DAY}
		else
			echo "Found no passive plot ... plotting"
			python2 ${CODE_PATH}/plot_pamtra_hatpro.py --date ${DAY} --rootpath ${ROOT_PATH}/ >> ${CODE_PATH}/plot${DAY}_hatpro.out
		fi
		for hydro in "${hydro_combo[@]}"; do
			for radar in "${radar_names[@]}"; do
				if [ -f ${DATA_PATH}/${hydro}/${DAY}${hydro}_mom_${radar}.nc  ]; then
					echo "Already processed " ${DAY} ${hydro} ${radar}
				else
					echo "Running "${DAY} ${hydro} ${radar}
					python2 ${CODE_PATH}/pamtra_radar+hatpro.py --date ${DAY} --hydroset ${hydro} --radarset ${radar} --datapath ${DATA_PATH}/ > ${CODE_PATH}/pamtra${DAY}_${hydro}_${radar}.out
					newactive=1
				fi
			done
			if [ "$newactive" -eq "1" ]; then
				echo "Newdata active ... plotting"
				python2 ${CODE_PATH}/plot_results_separated.py --date ${DAY} --hydroset ${hydro} --rootpath ${ROOT_PATH}/ >> ${CODE_PATH}/plot${DAY}_${hydro}.out
				newactive=0
			fi
			if [ -f ${PLOT_PATH}/${hydro}/${DAY}${hydro}_Ze.png ]; then
				echo "Already plotted " ${DAY} ${hydro}
			else
				echo "found no plot ... plotting"
				python2 ${CODE_PATH}/plot_results_separated.py --date ${DAY} --hydroset ${hydro} --rootpath ${ROOT_PATH}/ >> ${CODE_PATH}/plot${DAY}_${hydro}.out
			fi
		done
	else
		echo "no ICON data for "${DAY}
	fi
	#echo ${ICON_PATH}${DAY}_1mom/METEOGRAM_patch001_${DAY}_joyce.nc 
	if [ -f ${ICON_PATH}${DAY}_1mom/METEOGRAM_patch001_${DAY}_joyce.nc ]; then
		echo ${DAY}_1mom
		#if [ -f ${DATA_PATH}/${DAY}hatpro.nc ]; then
		#	echo "passive "${DAY}" already done"
		#else
		#	python2 ${CODE_PATH}/pamtra_radar+hatpro.py --date ${DAY} --datapath ${DATA_PATH}/ --hydroset all_hydro --radarset hatpro > ${CODE_PATH}/pamtra${DAY}_hatpro.out
		#	newpassive=1
		#fi
		#if [ "$newpassive" -eq "1" ]; then
		#	echo "New passive data ... plotting"
		#	python2 ${CODE_PATH}/plot_pamtra_hatpro.py --date ${DAY} --rootpath ${ROOT_PATH}/ >> ${CODE_PATH}/plot${DAY}_hatpro.out
		#	newpassive=0
		#fi
		#if [ -f ${PLOT_PATH}/${DAY}hatpro.png ]; then
		#	echo "Already plotted passive " ${DAY}
		#else
		#	echo "Found no passive plot ... plotting"
		#	python2 ${CODE_PATH}/plot_pamtra_hatpro.py --date ${DAY} --rootpath ${ROOT_PATH}/ >> ${CODE_PATH}/plot${DAY}_hatpro.out
		#fi
		for hydro in "${hydro_combo[@]}"; do
			for radar in "${radar_names[@]}"; do
				if [ -f ${DATA_PATH}_1mom/${hydro}/${DAY}${hydro}_mom_${radar}.nc  ]; then
					echo "Already processed " ${DAY} ${hydro} ${radar}
				else
					echo "Running "${DAY}_1mom ${hydro} ${radar}
					python2 ${CODE_PATH}/pamtra_radar+hatpro.py --date ${DAY}_1mom --hydroset ${hydro} --radarset ${radar} --datapath ${DATA_PATH}_1mom/ > ${CODE_PATH}/pamtra${DAY}_${hydro}_${radar}.out
					newactive=1
				fi
			done
			if [ "$newactive" -eq "1" ]; then
				echo "Newdata active ... plotting"
				python2 ${CODE_PATH}/plot_results_separated.py --date ${DAY} --hydroset ${hydro} --rootpath ${ROOT_PATH}/ -m1 _1mom >> ${CODE_PATH}/plot${DAY}_${hydro}.out
				newactive=0
			fi
			if [ -f ${PLOT_PATH}${hydro}/${DAY}${hydro}_Ze.png ]; then
				echo "Already plotted " ${DAY} ${hydro}
			else
				echo "found no plot ... plotting"
				python2 ${CODE_PATH}plot_results_separated.py --date ${DAY} --hydroset ${hydro} --rootpath ${ROOT_PATH}/ -m1 _1mom >> ${CODE_PATH}/plot${DAY}_${hydro}.out
			fi
		done
	else
		echo "no ICON 1 moment data for "${DAY}
	fi
	DAY=$(date -d "$DAY + 1 day" +%Y%m%d)
done

