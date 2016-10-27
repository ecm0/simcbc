#!/usr/bin/env bash

for d in `find $1 -maxdepth 1 -type d -name 'mdc*[0-9]'`; do
    for f in $d/*.fits.gz; do
	./compute_search_area $f
    done
done
