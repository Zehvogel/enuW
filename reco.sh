#! /usr/bin/bash
set -e

COMPACT_FILE=$K4GEO/FCCee/CLD/compact/CLD_o2_v05/CLD_o2_v05.xml
FILE=../../E250-SetA.P4f_sw_sl.Gwhizard-2_8_5.eL.pR.I500106.0
NEVT=10

pushd CLDConfig/CLDConfig

k4run CLDReconstruction.py -n $NEVT --inputFiles $FILE.10000.SIM.edm4hep.root --outputBasename $FILE.$NEVT.REC

popd
