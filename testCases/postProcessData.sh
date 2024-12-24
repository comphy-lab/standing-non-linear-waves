#!/bin/bash

# Post-Processing Script for Standing Non-Linear Waves Analysis
# Author: Vatsal Sanjay
# Last Updated: Dec 24, 2024
#
# Description:
#   This script processes simulation data for standing non-linear waves by:
#   1. Copying required processing scripts from postProcessScripts directory
#   2. Running video generation and data analysis
#
# Usage: 
#   ./postProcessData.sh <folderToProcess>
#
# Arguments:
#   folderToProcess: Directory containing the simulation data to be processed
#
# Required Files (from postProcessScripts/):
#   - video.py: Python script for visualization and video generation
#   - getData: Executable for extracting field data
#   - getFacets: Executable for extracting interface geometry

# Validate input
if [ -z "$1" ]; then
    echo "Error: Please provide the folder to process"
    echo "Usage: ./postProcessData.sh <folderToProcess>"
    exit 1
fi

folderToProcess=$1

# Copy required processing scripts
scp ../postProcessScripts/video.py .
scp ../postProcessScripts/getData .
scp ../postProcessScripts/getFacets .

# Run video processing script
# Additional options can be configured in video.py
python video.py --caseToProcess $folderToProcess

rm video.py
rm getData
rm getFacets