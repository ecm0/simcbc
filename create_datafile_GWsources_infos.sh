#!/usr/bin/env bash

for d in `find $1 -maxdepth 1 -type d -name 'mdc*[0-9]'`; do
    simdir=`basename $d`
    python -W ignore create_datafile_GWsources_infos.py $d/${simdir}.xml $d/coinc.xml 
done
