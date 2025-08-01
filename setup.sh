#!/bin/bash

echo "Loading required modules..."

source /project/ptgroup/spinquest/this-e1039.sh

# Load dependencies for Python
module load gcc/11.4.0
module load openmpi/4.1.4

# Load Python 3.11.4
module load python/3.11.4

# Load Apptainer for TensorFlow
module load apptainer/1.3.4

# Load TensorFlow (inside Apptainer)
module load tensorflow/2.10.0

# Load ROOT
module load root/6.32.06

# Get the directory of this script (project root)
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Add the project root to PYTHONPATH
export PYTHONPATH="$PROJECT_ROOT:$PYTHONPATH"

echo "PYTHONPATH set to: $PYTHONPATH"