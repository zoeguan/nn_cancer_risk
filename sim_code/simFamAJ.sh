#!/bin/bash

for i in `seq 1 50`;
do
	let num=i;
	argString="--args "$num
	R CMD BATCH "$argString" simFamAJ.R;
	sleep 10
done