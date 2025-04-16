#!/bin/bash

echo "Loading required modules..."

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