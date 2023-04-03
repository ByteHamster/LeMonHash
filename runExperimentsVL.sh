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

$pathToCppBinary --indexed --withoutAlphabetMaps --fixedEpsilon --numQueries 5M --filename "$pathToDatasetsFolder/dna-31-mer.txt"
$pathToCppBinary --indexed --withoutAlphabetMaps --fixedEpsilon --numQueries 5M --filename "$pathToDatasetsFolder/trec-text.terms"
$pathToCppBinary --indexed --withoutAlphabetMaps --fixedEpsilon --numQueries 5M --filename "$pathToDatasetsFolder/uk-2007-05.urls"

