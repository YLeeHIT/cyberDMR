#!/bin/bash
bash ../cyberDMR.sh -i ./input -o ./output -chr chr21 -g1 lethal -g2 normal -d 0.1 -q 0.05
bash ../simulate_data.sh -o ./simulate_data -t 100
