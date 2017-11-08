#!/bin/bash

#Convert all pdb files of a directory
#in mae files in the parent directory

var=1
for f in *.pdb; do
     /opt/schrodinger2016-4/utilities/pdbconvert -ipdb $f -omae ../${var}.mae
     var=$((var+1))
 done