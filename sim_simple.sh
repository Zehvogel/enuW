#! /usr/bin/bash
set -e

COMPACT_FILE=$K4GEO/FCCee/CLD/compact/CLD_o2_v06/CLD_o2_v06.xml
FILE=data/gen/simple_whizard/enuqq
NEVT=1000

ddsim --compactFile $COMPACT_FILE \
      --outputFile $FILE.$NEVT.SIM.edm4hep.root \
      --inputFile $FILE.slcio \
      --steeringFile CLDConfig/CLDConfig/cld_steer.py \
      --random.seed 0123456789 \
      --numberOfEvents $NEVT
