#!/bin/sh

## STEP4 collect bohrium jobs 

tool_file=/home/jywang/tool_file/DP_conf_gen/denovo_gen_workflow
LBG_UTIL=/home/jywang/tool_file/lbg_util/g16

cp ${LBG_UTIL}/collect.sh ./
bash collect.sh

#mv 2*/input_out.sdf ./
#rm -rf 2023*

## STEP5 processing it for final results
cp ${tool_file}/run4ShermoCalc.py ./
python run4ShermoCalc.py 

## FINAL CLEAN
rm -f *.py *.sh *.yaml *.txt