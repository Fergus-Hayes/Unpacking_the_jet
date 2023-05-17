#!/bin/bash

INPATH=$1
NAME=$2
OUTPATH=$3

if [ "${INPATH: -1}" != '/' ]; then
INPATH="${INPATH}/"
fi

declare -a JSS=("GJ" "TH" "PL" "DG")

for STRUCT in "${JSS[@]}" ; do

LOC3="${INPATH}rates_${STRUCT}/`ls ${INPATH}rates_${STRUCT}/`/temp/${STRUCT}/1/trace/rates_${STRUCT}_result.json"
LOC1="${INPATH}GW170817_${STRUCT}/`ls ${INPATH}GW170817_${STRUCT}/`/temp/${STRUCT}/1/trace/GW170817_${STRUCT}_result.json"
LOC2="${INPATH}GW190425_${STRUCT}/`ls ${INPATH}GW190425_${STRUCT}/`/temp/${STRUCT}/1/trace/GW190425_${STRUCT}_result.json"
LOC4="${INPATH}rates_GW170817_${STRUCT}/`ls ${INPATH}rates_GW170817_${STRUCT}/`/temp/${STRUCT}/1/trace/rates_GW170817_${STRUCT}_result.json"
LOC5="${INPATH}rates_GW170817_GW190425_${STRUCT}/`ls ${INPATH}rates_GW170817_GW190425_${STRUCT}/`/temp/${STRUCT}/1/trace/rates_GW170817_GW190425_${STRUCT}_result.json"
LAB3="sGRB_rate"
LAB1="GW170817"
LAB2="GW190425"
LAB4="sGRB_rate_+_GW170817"
LAB5="sGRB_rate_+_GW170817_+_GW190425"
python3 plot_corner_multiplot.py --loc ${LOC1} ${LOC2} ${LOC3} ${LOC4} ${LOC5} --labels ${LAB1} ${LAB2} ${LAB3} ${LAB4} ${LAB5} --outpath ${OUTPATH}${STRUCT}_${NAME}_ --JS ${STRUCT}

done
