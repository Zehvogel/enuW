#
# Copyright (c) 2014-2024 Key4hep-Project.
#
# This file is part of Key4hep.
# See https://key4hep.github.io/key4hep-doc/ for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
import os
from Gaudi.Configuration import INFO, WARNING, DEBUG

from Configurables import k4DataSvc, MarlinProcessorWrapper
from k4MarlinWrapper.inputReader import create_reader, attach_edm4hep2lcio_conversion
from k4FWCore.parseArgs import parser


parser.add_argument("--inputFiles", action="extend", nargs="+", metavar=("file1", "file2"), help="One or multiple input files")
parser.add_argument("--outputBasename", help="Basename of the output file(s)", default="output")
reco_args = parser.parse_known_args()[0]

algList = []
svcList = []

evtsvc = k4DataSvc("EventDataSvc")
svcList.append(evtsvc)

if reco_args.inputFiles:
    read = create_reader(reco_args.inputFiles, evtsvc)
    read.OutputLevel = INFO
    algList.append(read)
else:
    read = None

MyAIDAProcessor = MarlinProcessorWrapper("MyAIDAProcessor")
MyAIDAProcessor.OutputLevel = WARNING
MyAIDAProcessor.ProcessorType = "AIDAProcessor"
MyAIDAProcessor.Parameters = {
                              "Compress": ["1"],
                              "FileName": [f"{reco_args.outputBasename}_aida"],
                              "FileType": ["root"]
                              }

# setup AIDA histogramming and add eventual background overlay
# algList.append(MyAIDAProcessor)

EventNumber = MarlinProcessorWrapper("EventNumber")
EventNumber.OutputLevel = INFO
EventNumber.ProcessorType = "Statusmonitor"
EventNumber.Parameters = {
                          "HowOften": ["1"]
                          }
algList.append(EventNumber)

CheatedOverlayRemoval = MarlinProcessorWrapper("CheatedOverlayRemoval")
CheatedOverlayRemoval.OutputLevel = WARNING
CheatedOverlayRemoval.ProcessorType = "CheatedMCOverlayRemoval"
CheatedOverlayRemoval.Parameters = {
    "RecoParticleCollection": ["PandoraPFOs"],
    "MCParticleCollection": ["MCParticlesSkimmed"],
    "RecoMCTruthLink": ["RecoMCTruthLink"],
    "MCTruthRecoLink": ["MCTruthRecoLink"],
    "OutputPfoCollection": ["PFOsWithoutMCOverlay"],
    "OutputOverlayCollection": ["PFOsFromOverlay"],
}
algList.append(CheatedOverlayRemoval)

truejet = MarlinProcessorWrapper("truejet")
truejet.OutputLevel = WARNING
truejet.ProcessorType = "TrueJet"
truejet.Parameters = {
    # "MCParticleCollection": ["MCParticle"]
    "MCParticleCollection": ["MCParticlesSkimmed"]
}
algList.append(truejet)

from Configurables import Lcio2EDM4hepTool
lcioConvTool = Lcio2EDM4hepTool("lcio2EDM4hep")
lcioConvTool.convertAll = True
lcioConvTool.collNameMapping = {
    # "MCParticle": "MCParticles"
}
lcioConvTool.OutputLevel = WARNING
# attach to the last non output processor
truejet.Lcio2EDM4hepTool = lcioConvTool

from Configurables import PodioOutput
out = PodioOutput("PodioOutput", filename = f"{reco_args.outputBasename}.truejet.edm4hep.root")
out.outputCommands = ["keep *"]
algList.append(out)

# We need to convert the inputs in case we have EDM4hep input
attach_edm4hep2lcio_conversion(algList, read)

from Configurables import ApplicationMgr
ApplicationMgr( TopAlg = algList,
                EvtSel = 'NONE',
                EvtMax = 3, # Overridden by the --num-events switch to k4run
                ExtSvc = svcList,
                OutputLevel=WARNING
              )
