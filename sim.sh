#! /usr/bin/bash
set -e

COMPACT_FILE=$K4GEO/FCCee/CLD/compact/CLD_o2_v05/CLD_o2_v05.xml
FILE=E250-SetA.P4f_sw_sl.Gwhizard-2_8_5.eL.pR.I500106.0
NEVT=10

ddsim --compactFile $COMPACT_FILE \
      --outputFile $FILE.$NEVT.SIM.edm4hep.root \
      --inputFile $FILE.slcio \
      --steeringFile CLDConfig/CLDConfig/cld_steer.py \
      --random.seed 0123456789 \
      --numberOfEvents $NEVT
