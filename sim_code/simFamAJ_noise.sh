#!/bin/bash

for i in `seq 1 50`;
do
	let num=i;
	argString="--args "$num
	R CMD BATCH "$argString" simFamAJ_noise.R;
	sleep 1
done