#!/bin/sh


sdf_tag=${1}
tool_file=/home/jywang/tool_file/DP_conf_gen/denovo_gen_workflow
LBG_UTIL=/home/jywang/tool_file/lbg_util/g16

## STEP1 prepare and run auto3D
#conda init
cp ${tool_file}/run4auto3D.py ./
python run4auto3D.py --sdf_tag ${sdf_tag}
cp ${tool_file}/auto3D_param/auto3D.* ./
./auto3D.sh
cp 2*/input_out.sdf ./
#conda deactivate

rm -f ${sdf_tag}*.sdf

## STEP2 processing for g16 run
cp ${tool_file}/run4g16.py ./
python run4g16.py --input_sdf input_out.sdf
rm -f input_out.sdf
rm -f *.xyz *.sdf

## STEP3 run in bohrium
cp ${LBG_UTIL}/submit.sh ./
bash submit.sh

##
