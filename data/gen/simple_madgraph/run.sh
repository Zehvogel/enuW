#!/usr/bin/env bash
mg5_aMC lnuqq.dat
gunzip Output/Events/run_01/unweighted_events.lhe.gz
ln -sf Output/Events/run_01/unweighted_events.lhe unweighted_events.lhe
sed -i '/<header>/,/<\/header>/{//!d}' unweighted_events.lhe
$CONVERTHEPMC2EDM4HEP/convertHepMC2EDM4HEP -i lhe -o hepmc3 unweighted_events.lhe enuqq.hepmc
$CONVERTHEPMC2EDM4HEP/convertHepMC2EDM4HEP -i hepmc3 -o edm4hep enuqq.hepmc enuqq.edm4hep.root
