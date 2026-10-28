#!/bin/bash
# Build GORILLA and run the MATLAB and Python examples with their plotting
# scripts and tutorials. Needs MATLAB (matlab-batch) and the Python packages
# installed by the nightly workflow.
set -euo pipefail
make
export LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libcurl.so.4
matlab-batch "cd MATLAB; example_1; example_2; example_3; example_4; example_5; example_6; example_7; example_8"
matlab-batch "cd MATLAB; plotting_tutorial"
unset LD_PRELOAD
cd PYTHON
for i in 1 2 3 4 5 6 7; do python3 "example_$i.py"; done
for i in 1 2 3 4 5 6 7; do python3 "plot_example_$i.py"; done
jupyter nbconvert --to script plotting_tutorial.ipynb
python3 plotting_tutorial.py
