#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Illegal number of parameters."
    echo "Usage: runExperiments.sh pathToDatasetsFolder pathToCppBinary"
    exit 1
fi

pathToDatasetsFolder=$1
pathToCppBinary=$2

hostname
strings $pathToCppBinary | grep " -m"

mappersExceptLinear="--segmented --succinctPgm --leMonVl"
allMappers="--linear $mappersExceptLinear"

$pathToCppBinary $allMappers --numQueries 5M --type int32 --filename "$pathToDatasetsFolder/5GRAM_1"
$pathToCppBinary $mappersExceptLinear --numQueries 5M --filename "$pathToDatasetsFolder/fb_200M_uint64"
$pathToCppBinary --linear --numQueries 10k --filename "$pathToDatasetsFolder/fb_200M_uint64"
$pathToCppBinary $allMappers --numQueries 5M --filename "$pathToDatasetsFolder/osm_cellids_800M_uint64"

$pathToCppBinary $allMappers --numQueries 5M --filename "$pathToDatasetsFolder/uniform_uint64"
$pathToCppBinary $allMappers --numQueries 5M --filename "$pathToDatasetsFolder/exponential_uint64"
$pathToCppBinary $mappersExceptLinear --numQueries 5M --filename "$pathToDatasetsFolder/normal_uint64"
$pathToCppBinary --linear --numQueries 10k --filename "$pathToDatasetsFolder/normal_uint64"

