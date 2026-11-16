#!/bin/bash
# After prep_barrier_cont.sh:  continue the unconverged images (list ads/list_bar_cont.txt) in one mediumq job (3 nodes, 20 tasks).
sbatch --partition=mediumq --nodes=3 --exclude=compute28 --mem=60G --time=80:00:00 --export=ALL,LISTS="ads/list_bar_cont.txt" submit_ads.sh
