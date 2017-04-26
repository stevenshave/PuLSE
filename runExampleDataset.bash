#!/usr/bin/env bash

if [ ! -f pulse ]; then
    echo -e "pulse not found, please compile it first with \"./configure && make\".\nWith the pulse executable in this directory, this script may be used to run the included example."
    exit
fi

./pulse sample-pulse-5merCyclic-CGTTGCXXXXXXXXXXXXXXXTGTGCT.fastq.gz CGTTGCXXXXXXXXXXXXXXXTGTGCT UAG Q

echo -e "Pulse run completed, you should now be able to open the newly generated HTML \nand txt reports - sample-pulse-5merCyclic-CGTTGCXXXXXXXXXXXXXXXTGTGCT.html\n & sample-pulse-5merCyclic-CGTTGCXXXXXXXXXXXXXXXTGTGCT.txt"
