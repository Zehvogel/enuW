################################################
#
#  Simple configuration using the MarlinWrapper
#  to convert an LCIO file to EDM4hep
#  F.Gaede, DESY
#  Nov 2022
#################################################

from Gaudi.Configuration import *

from Configurables import LcioEvent, EventDataSvc, MarlinProcessorWrapper
from Configurables import PodioOutput, ToolSvc, Lcio2EDM4hepTool
from Configurables import k4DataSvc

podioEvtSvc = k4DataSvc('EventDataSvc')

#====================================================================

#fileNameBase = "E250-SetA.P4f_sw_sl.Gwhizard-2_8_5.eL.pR.I500106.0"
fileNameBase = "E250-SetA.P4f_sw_sl.Gwhizard-2_8_5.eR.pL.I500108.0"

inFileName = fileNameBase +".slcio"
outFileName = fileNameBase +"_edm4hep.root"

maxEvt = 100000

#====================================================================
# --- for now we need to get the number of events in the file first (should be fixed soon)
import subprocess
result = subprocess.run( ['lcio_event_counter', inFileName] , stdout=subprocess.PIPE ).stdout.decode('utf-8')
nEvt=int(result)
if nEvt > maxEvt:
    nEvt = maxEvt
#====================================================================


# --- conversion tool:  configure the collections to convert
edmConvTool = Lcio2EDM4hepTool('LCIO2EDM4hep')
edmConvTool.convertAll = False
edmConvTool.collNameMapping = {
    "MCParticle" : "MCParticles",
    "MarlinTrkTracks" : "MarlinTrkTracks",
    "MarlinTrkTracksMCTruthLink" : "MarlinTrkTracksMCTruthLink",
    "RecoMCTruthLink" : "RecoMCTruthLink",
    "PandoraClusters" : "PandoraClusters",
    "PandoraPFOs" : "PandoraPFOs"
}
#edmConvTool.OutputLevel = DEBUG

#====================================================================


# --- reader and writer :
read = LcioEvent()
read.OutputLevel = INFO
read.Files = [ inFileName ]

write = PodioOutput('EDM4hepOutput', filename=outFileName )
write.outputCommands = ["keep *"]


# --- need to attach the conversion tool to one Marlin processor - use StatusMonitor here
status = MarlinProcessorWrapper("MyStatusmonitor")
status.OutputLevel = INFO
status.ProcessorType = "Statusmonitor"
status.Parameters = {"HowOften": ["1"]}

status.Lcio2EDM4hepTool = edmConvTool


algList = []
algList.append(read)
algList.append(status)
algList.append(write)


from Configurables import ApplicationMgr
ApplicationMgr( TopAlg = algList,
                EvtSel = 'NONE',
                EvtMax   = nEvt,
                ExtSvc = [podioEvtSvc],
                OutputLevel=INFO
              )
# -----------------------------------------------------------------------------------
