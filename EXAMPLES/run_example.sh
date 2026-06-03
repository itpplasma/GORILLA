#!/bin/bash

# Script to run a GORILLA example
# Usage from parent directory: .EXAMPLES/run_example.sh <example_number>

if [ -z "$1" ]; then
    echo "Usage: $0 <example_number>"
    echo "Example: $0 1"
    exit 1
fi

EXAMPLE_NUM=$1
EXAMPLE_DIR="EXAMPLES/example_$EXAMPLE_NUM"

# Check if example directory exists
if [ ! -d "$EXAMPLE_DIR" ]; then
    echo "Error: $EXAMPLE_DIR does not exist"
    exit 1
fi

# Check if executable exists
if [ ! -f "BUILD/test_gorilla_main.x" ]; then
    echo "Error: Executable BUILD/test_gorilla_main.x not found. Please run 'make build' first."
    exit 1
fi

echo "Running example $EXAMPLE_NUM from $EXAMPLE_DIR"
cd "$EXAMPLE_DIR"
./test_gorilla_main.x

if [ $? -eq 0 ]; then
    echo "Example $EXAMPLE_NUM completed successfully"
    
    # Validate expected output files exist based on example number
    case $EXAMPLE_NUM in
        1)
            EXPECTED_FILES=(
                "poincare_plot_vpar_0_rphiz.dat"
                "poincare_plot_vpar_0_sthetaphi.dat"
                "J_par.dat"
            )
            ;;
        2)
            EXPECTED_FILES=(
                "poincare_plot_phi_0_rphiz.dat"
                "poincare_plot_phi_0_sthetaphi.dat"
                "e_tot.dat"
            )
            ;;
        3)
            EXPECTED_FILES=(
                "poincare_plot_phi_0_rphiz.dat"
                "poincare_plot_phi_0_sthetaphi.dat"
                "p_phi.dat"
            )
            ;;
        4)
            EXPECTED_FILES=(
                "poincare_plot_phi_0_rphiz.dat"
                "poincare_plot_phi_0_sthetaphi.dat"
            )
            ;;
        5)
            EXPECTED_FILES=(
                "full_orbit_plot_rphiz_trapped.dat"
                "poincare_plot_phi_0_rphiz_trapped.dat"
            )
            ;;
        6)
            EXPECTED_FILES=(
                "full_orbit_plot_rphiz_passing.dat"
                "orbit_start_rphizlambda_passing.dat"
            )
            ;;
        7)
            EXPECTED_FILES=(
                "poincare_plot_vpar_0_rphiz_adaptive.dat"
                "poincare_plot_vpar_0_sthetaphi_adaptive.dat"
                "J_par_adaptive.dat"
                "e_tot_adaptive.dat"
            )
            ;;
        8)
            EXPECTED_FILES=(
                "poincare_plot_phi_0_rphiz.dat"
            )
            ;;
        *)
            echo "Warning: Unknown example number, skipping output validation"
            EXPECTED_FILES=()
            ;;
    esac
    
    for file in "${EXPECTED_FILES[@]}"; do
        if [ ! -f "$file" ]; then
            echo "Error: Expected output file $file not found"
            exit 1
        fi
    done
    echo "Output validation passed: all expected files present"
    
    echo "Plotting results for example $EXAMPLE_NUM from $EXAMPLE_DIR"
    cd ../../PYTHON
    python3 plot_example_$EXAMPLE_NUM.py
    
    # Validate plot was created
    if [ ! -f "example_$EXAMPLE_NUM.png" ]; then
        echo "Error: Plot file example_$EXAMPLE_NUM.png not created"
        exit 1
    fi
    echo "Plot validation passed: example_$EXAMPLE_NUM.png created"
else
    echo "Error: Example $EXAMPLE_NUM failed"
    exit 1
fi
    