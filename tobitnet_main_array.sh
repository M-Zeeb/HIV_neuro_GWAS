#!/bin/bash

echo "$1"

Rscript --vanilla ./tobit_crossval_cluster.R $1

