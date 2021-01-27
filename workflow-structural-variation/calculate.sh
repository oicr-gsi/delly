#!/bin/bash
cd $1
for z in *gz;do echo $z;zcat $z | grep -v ^# | wc -l;done | sed 's!.*/!!'
