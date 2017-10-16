#!/bin/bash

mkdir resources 

#wget -P resources/ "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/*"

gunzip resources/hs37d5.fa.gz

#bwa index resources/hs37d5.fa
