#!/usr/bin/env bash

if [ "$(uname)" == "Darwin" ]; then
    curl -o blast.tar.gz 'https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.7.1/ncbi-blast-2.7.1+-x64-macosx.tar.gz'
    gunzip -c blast.tar.gz | tar xopf -
    mv ncbi-blast-2.7.1+ ncbi-blast
elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
    wget -O blast.tar.gz https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.7.1/ncbi-blast-2.7.1+-x64-linux.tar.gz
    tar -xzf blast.tar.gz
    mv ncbi-blast-2.7.1+ ncbi-blast
fi
