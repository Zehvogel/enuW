from DIRAC.Core.Base import Script
Script.parseCommandLine()

from ILCDIRAC.Interfaces.API.DiracILC import DiracILC
from ILCDIRAC.Interfaces.API.NewInterface.UserJob import UserJob
from ILCDIRAC.Interfaces.API.NewInterface.Applications import GaudiApp

from glob import glob

dIlc = DiracILC()

full = True

files = []
if full:
    paths = glob("/eos/experiment/fcc/prod/fcc/user/L/LReichenbac/merged*/*.slcio")
    files = [path.removeprefix("/eos/experiment/fcc/prod") for path in paths]
else:
    files = ["/fcc/user/L/LReichenbac/merged1/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500105.P4f_sw_sl.eL.pL.n000.d_dstm_15064_0.slcio"]

inputFiles = ["LFN:" + file for file in files]

# print(inputFiles)

file_basenames = [f.split("/")[-1].strip(".slcio") for f in files]
outputFiles = [[f"{name}.miniDST.slcio"] for name in file_basenames]

job = UserJob()
job.setSplitDoNotAlterOutputFilename()
job.setName('paris_miniDST_%n')
job.setSplitInputData(files)
job.setSplitParameter('outputFile', outputFiles)
job.setSplitOutputData(outputFiles, f'paris/miniDSTs{"-test2" if not full else ""}', 'CERN-DST-EOS')
job.setOutputSandbox('*.log')
# too big
# job.setInputSandbox(["../../ILDConfig/StandardConfig/production"])
job.setInputSandbox(["LFN:/fcc/user/L/LReichenbac/paris/input-env/ild-std-prod.tar.gz", "Ref_singleP_16bins_conservative.txt"])

gaudi = GaudiApp()
gaudi.setExecutableName("k4run")
gaudi.setVersion("key4hep_nightly")
gaudi.setSteeringFile("MiniDST.py")
gaudi.setInputFileFlag("--inputFiles")
gaudi.setOutputFileFlag("--outputFile")
gaudi.setOutputFile("%(outputFile)s")
gaudi.setExtraCLIArguments(
    "-n -1"
    )


job.append(gaudi)
job.submit(dIlc, mode="wms")