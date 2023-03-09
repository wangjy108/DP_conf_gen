#!/bin/sh

#cd ${1}
#runfile=`ls *.sdf`
folder=`ls -l | grep "^d" | awk '{print $9}'`

for ff in ${folder}
do
	cd $ff
	conformer=`ls *.log`
	touch list.txt
	for each in ${conformer}
	do
		EE=`grep 'E(RB2PLYPD3)' ${each} | python -c "n=input();print(n.split()[4])"`
		echo ${each}';'${EE} >> list.txt
		Shermo list.txt | grep 'Boltzmann weight' | awk '{print $9}' > Dis.dat
	done
	cd ../
done
