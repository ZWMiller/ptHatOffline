#!/bin/bash
if [ "$1" == "-h" ]; then
    printf "\nUsage, source scpPythia.sh with 3 options: \nTemplate type with ptbin: C0_1\nDate: Aug27\nFile Number: 1 (first file on date)\n\n"
    return
fi

if [ "$#" -eq 3 ]; then
scp zamiller@rftpexp.rhic.bnl.gov:/star/u/zamiller/simu/ptHatTemplates/output/npe$1/pythia_tree_$2_$3_$1.root outputs/.
else
    echo "Wrong number of arguments. Need 'C0_1', 'Aug23', '2'."
    return
fi

if [ $? -eq 0 ]; then
    echo "pythia_tree_$2_$3.root moved to local directory."
else
    echo "---SCP FAILED---"
fi

