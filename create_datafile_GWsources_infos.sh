#!/usr/bin/env bash

for d in `find . -type d -name 'mdc*[0-9]'`; do
    python -W ignore create_datafile_GWsources_infos.py $d/$d.xml $d/coinc.xml 
done
