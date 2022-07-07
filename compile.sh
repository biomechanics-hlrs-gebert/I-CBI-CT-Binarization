#!/bin/bash
# ----------------------------------------------------------------------------------------
# Wrapper to always acces a makefile
#
# Author: Johannes Gebert - HLRS - NUM - gebert@hlrs.de
# Created on: 05.07.2022
# ----------------------------------------------------------------------------------------
#
SEP="-------------------------------------------------------------------------------"
echo $SEP
# ----------------------------------------------------------------------------------------
#
if [ $# -ne 2 ]; then
    echo "Usage: ./compile <optional makefile> <optional make target>"
fi
# ----------------------------------------------------------------------------------------
#
if [ -f "$1" ]; then
    make -f "$1" "$2"
else
    # ----------------------------------------------------------------------------------------
    #
    if [ -f makefile ] || [ -f Makefile ] || [ -f GNUmakefile ]; then
        make
    else
        echo $SEP
        #
        for file in make.*
        do 
            echo "Makefile found: $file"
        done
        #
        echo $SEP
        echo "You may compile the file with 'make -f <make.*> <target>"
        echo ""
        echo "Examples: make -f make.CIF cleanall"
        echo "          make -f make.CIF all"
    fi
fi
#
echo $SEP