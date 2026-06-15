#!/bin/bash
#### works on uni pcs and remote access (emailed it support for login)
# Basic HPC access / environment check for C. auris annotation project

echo "=== HPC ACCESS CHECK ==="
echo "User:"
whoami

echo ""
echo "Host:"
hostname

echo ""
echo "Current directory:"
pwd

echo ""
echo "Date:"
date

echo ""
echo "Available modules:"
module avail 2>&1 | head -50

echo ""
echo "Python version:"
python --version 2>/dev/null || echo "Python not found"

echo ""
echo "Conda version:"
conda --version 2>/dev/null || echo "Conda not found"

echo ""
echo "GPU check:"
nvidia-smi 2>/dev/null || echo "No GPU visible or nvidia-smi unavailable"

echo ""
echo "Disk space:"
df -h . | tail -1

echo ""
echo "Memory:"
free -h 2>/dev/null || echo "free command unavailable"

echo ""
echo "Suggested next step:"
echo "Create separate conda environments for ProtNote/GO2Sum if dependencies conflict."