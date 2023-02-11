#!/bin/bash
# @author: sharib
# install miniconda or anaconda from here:
# creates the conda environment
create_env() {
    conda create -n "$1" python=3.7
    conda activate "$1"
    pip install -r requirements.txt
}

# checks if the conda environment exists
env_exists() {
conda info --envs | awk '{print $1}' | tail -n +3 | grep -w "$1" > /dev/null
}

create_if () {
    env_exists "$1" && echo "environment exists" || create_env "$1"
    }
    
create_if myenv

#conda env remove --name my_env
