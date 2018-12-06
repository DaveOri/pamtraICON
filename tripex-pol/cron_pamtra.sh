#!/bin/bash

ICON_PATH='/data/inscape/icon/experiments/juelich/testbed/testbed_'
DATA_PATH='/data/optimice/pamtra_runs/tripex-pol/data/'
CODE_PATH='/home/dori/develop/pamtraICON/tripex-pol/'

FIRST_DAY='20181101'
TODAY=`date +%Y%m%d`

DAY=$FIRST_DAY
echo $DAY $TODAY
until [[ ${DAY} > ${TODAY} ]]; do
	# Check if there is ICON output
	if [ -f ${ICON_PATH}${DAY}/METEOGRAM_patch001_${DAY}_joyce.nc ]; then
		echo ${DAY}
		if [ -f ${DATA_PATH}${DAY}all_hydro_mom.nc  ]; then
			echo "Day already processed "${DAY}
		else
			echo "Running "${DAY}
			python ${CODE_PATH}run_pamtra_group_hydro.py --date ${DAY} > pamtra${DAY}.out
			python ${CODE_PATH}plot_results.py --date ${DAY} > plot${DAY}.out
		fi
	else
		echo "no ICON data for "${DAY}
	fi
	DAY=$(date -d "$DAY + 1 day" +%Y%m%d)
done