#!/bin/bash
./run_SC2300771.sh
mkdir -p /mnt/pikachu/prod-12-27/SC2300771
mv /storage/production/SC2300771 /mnt/pikachu/prod-12-27/SC2300771
./run_SC2300772.sh
