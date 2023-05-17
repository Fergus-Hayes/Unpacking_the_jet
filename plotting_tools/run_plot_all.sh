#!/bin/bash

IP=${1}
OP=${2}

HLP="${IP}hierarchical_luminosity/"
SLP="${IP}schechter_luminosity/"

bash plot_all_joint_corners.sh ${HLP} hiecrachical_luminosity ${OP}Hierarchical_
bash plot_all_joint_corners.sh ${SLP} schechter_luminosity ${OP}Schecter_

python3 plot_rate.py --loc ${HLP}rates_GW170817_GW190425_DG/temp/DG/1/trace/rates_GW170817_GW190425_DG_result.json ${HLP}rates_GW170817_GW190425_PL/temp/PL/1/trace/rates_GW170817_GW190425_PL_result.json ${HLP}rates_GW170817_GW190425_GJ/temp/GJ/1/trace/rates_GW170817_GW190425_GJ_result.json ${HLP}rates_GW170817_GW190425_TH/temp/TH/1/trace/rates_GW170817_GW190425_TH_result.json --labels DG PL GJ TH --outpath ${OP}Hierarchical_

python3 plot_lum.py --loc ${HLP}rates_GW170817_GW190425_DG/temp/DG/1/trace/rates_GW170817_GW190425_DG_result.json ${HLP}rates_GW170817_GW190425_PL/temp/PL/1/trace/rates_GW170817_GW190425_PL_result.json ${HLP}rates_GW170817_GW190425_GJ/temp/GJ/1/trace/rates_GW170817_GW190425_GJ_result.json ${HLP}rates_GW170817_GW190425_TH/temp/TH/1/trace/rates_GW170817_GW190425_TH_result.json --labels DG PL GJ TH --outpath ${OP}Hierarchical_

python3 plot_rate.py --loc ${SLP}rates_GW170817_GW190425_DG/temp/DG/1/trace/rates_GW170817_GW190425_DG_result.json ${SLP}rates_GW170817_GW190425_PL/temp/PL/1/trace/rates_GW170817_GW190425_PL_result.json ${SLP}rates_GW170817_GW190425_GJ/temp/GJ/1/trace/rates_GW170817_GW190425_GJ_result.json ${SLP}rates_GW170817_GW190425_TH/temp/TH/1/trace/rates_GW170817_GW190425_TH_result.json --labels DG PL GJ TH  --outpath ${OP}Schecter_
