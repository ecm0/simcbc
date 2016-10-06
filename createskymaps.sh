#!/usr/bin/env bash

for d in `find $1 -maxdepth 1 -type d -name 'mdc*[0-9]'`; do
echo "Entering $d"
for f in $d/*.fits.gz; do
    echo "Converting $f"
    bayestar_plot_allsky $f --contour 90 --radec 0.0 0.0 -o ${f%.*}.png
done
done
