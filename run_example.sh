#!/bin/bash

# Script to run a GORILLA example
# Usage: ./run_example.sh <example_number>

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
    echo "Error: Executable BUILD/SRC/test_gorilla_main.x not found. Please run 'make build' first."
    exit 1
fi

echo "Running example $EXAMPLE_NUM from $EXAMPLE_DIR"
cd "$EXAMPLE_DIR"
./test_gorilla_main.x

if [ $? -eq 0 ]; then
    echo "Example $EXAMPLE_NUM completed successfully"
else
    echo "Error: Example $EXAMPLE_NUM failed"
    exit 1
fi
    
echo "Plotting results for example $EXAMPLE_NUM from $EXAMPLE_DIR"
cd ../../PYTHON
python3 plot_example_$EXAMPLE_NUM.py