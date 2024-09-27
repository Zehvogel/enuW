from Gaudi.Configuration import *

from Configurables import LcioEvent, EventDataSvc, MarlinProcessorWrapper
from k4MarlinWrapper.parseConstants import *
import os

from k4FWCore.parseArgs import parser

parser.add_argument("--inputFiles", action="extend", nargs="+", metavar=("file1", "file2"), help="One or multiple input files")
parser.add_argument("--outputFile", help="Name of the output file", default="output")

reco_args = parser.parse_known_args()[0]

algList = []
evtsvc = EventDataSvc()


CONSTANTS = {
             'OutputFile': reco_args.outputFile,
             'ProductionDir': "./ILDConfig/StandardConfig/production",
             'RundEdxCorrections': "true",
             'lcgeo_DIR': "/path/to/lcgeo_DIR",
             'LCFIPlusConfig_DIR': "./ILDConfig/StandardConfig/production/LCFIPlusConfig",
             'LCFIPlusVertexPrefix': "ildl5_4q250_ZZ",
             'LCFIPlusD0ProbFile': "%(LCFIPlusConfig_DIR)s/vtxprob/d0probv2_%(LCFIPlusVertexPrefix)s.root",
             'LCFIPlusZ0ProbFile': "%(LCFIPlusConfig_DIR)s/vtxprob/z0probv2_%(LCFIPlusVertexPrefix)s.root",
             'LCFIPlusWeightsPrefix': "4q250_v04_p00_ildl5",
             'LCFIPlusWeightsDir': "%(LCFIPlusConfig_DIR)s/lcfiweights/4q250_ZZ_v4_p00_ildl5",
             'IsolatedLeptonTagging_DIR': "./ILDConfig/StandardConfig/production/IsolatedLeptonTagging",
             'ElectronIsolationWeightsDir': "%(IsolatedLeptonTagging_DIR)s/weights/e1e1h_gg_qqqq_250",
             'MuonIsolationWeightsDir': "%(IsolatedLeptonTagging_DIR)s/weights/e2e2h_gg_qqqq_250",
}

parseConstants(CONSTANTS)

read = LcioEvent()
read.OutputLevel = ERROR
read.Files = reco_args.inputFiles
algList.append(read)

from Configurables import GeoSvc
geoservice = GeoSvc("GeoSvc")
geoservice.detectors = [os.environ["K4GEO"]+"/ILD/compact/ILD_l5_v02/ILD_l5_v02.xml"]
geoservice.OutputLevel = INFO
geoservice.EnableGeant4Geo = False

Statusmonitor = MarlinProcessorWrapper("Statusmonitor")
Statusmonitor.OutputLevel = ERROR
Statusmonitor.ProcessorType = "Statusmonitor"
Statusmonitor.Parameters = {
                            "HowOften": ["1000"]
                            }

MyTrueJet = MarlinProcessorWrapper("MyTrueJet")
MyTrueJet.OutputLevel = ERROR
MyTrueJet.ProcessorType = "TrueJet"
MyTrueJet.Parameters = {
                        "MCParticleCollection": ["MCParticlesSkimmed"]
                        }

MyCheatedMCOverlayRemoval = MarlinProcessorWrapper("MyCheatedMCOverlayRemoval")
MyCheatedMCOverlayRemoval.OutputLevel = ERROR
MyCheatedMCOverlayRemoval.ProcessorType = "CheatedMCOverlayRemoval"
MyCheatedMCOverlayRemoval.Parameters = {
                                        "MCParticleCollection": ["MCParticlesSkimmed"],
                                        "MCTruthRecoLink": ["MCTruthRecoLink"],
                                        "OutputOverlayCollection": ["PFOsFromOverlay"],
                                        "OutputPfoCollection": ["PFOsWithoutMCOverlay"],
                                        "RecoMCTruthLink": ["RecoMCTruthLink"],
                                        "RecoParticleCollection": ["PandoraPFOs"]
                                        }

FastJetOverlay = MarlinProcessorWrapper("FastJetOverlay")
FastJetOverlay.OutputLevel = ERROR
FastJetOverlay.ProcessorType = "FastJetProcessor"
FastJetOverlay.Parameters = {
                             "algorithm": ["ee_genkt_algorithm", "3.0", "1.0"],
                             "clusteringMode": ["ExclusiveNJets", "2"],
                             "jetOut": ["PFOsminusoverlayJets"],
                             "recParticleIn": ["PandoraPFOs"],
                             "recombinationScheme": ["E_scheme"]
                             }

ExpandJet = MarlinProcessorWrapper("ExpandJet")
ExpandJet.OutputLevel = ERROR
ExpandJet.ProcessorType = "ExpandJetProcessor"
ExpandJet.Parameters = {
                        "InputCollection": ["PFOsminusoverlayJets"],
                        "OutputCollection": ["PFOsminusoverlay"]
                        }

Thrust = MarlinProcessorWrapper("Thrust")
Thrust.OutputLevel = ERROR
Thrust.ProcessorType = "ThrustReconstruction"
Thrust.Parameters = {
                     "inputCollectionName": ["PandoraPFOs"],
                     "typeOfThrustFinder": ["2"]
                     }

Sphere = MarlinProcessorWrapper("Sphere")
Sphere.OutputLevel = ERROR
Sphere.ProcessorType = "Sphere"
Sphere.Parameters = {
                     "CollectionName": ["PandoraPFOs"],
                     "r_value": ["2.0"]
                     }

Fox = MarlinProcessorWrapper("Fox")
Fox.OutputLevel = ERROR
Fox.ProcessorType = "Fox"
Fox.Parameters = {
                  "NameOfReconstructedParticlesCollection": ["PandoraPFOs"]
                  }

IsolatedMuonTagging = MarlinProcessorWrapper("IsolatedMuonTagging")
IsolatedMuonTagging.OutputLevel = ERROR
IsolatedMuonTagging.ProcessorType = "IsolatedLeptonTaggingProcessor"
IsolatedMuonTagging.Parameters = {
                                  "CosConeLarge": ["0.95"],
                                  "CosConeSmall": ["0.98"],
                                  "CutOnTheISOElectronMVA": ["2"],
                                  "CutOnTheISOMuonMVA": ["0.7"],
                                  "DirOfISOElectronWeights": ["%(ElectronIsolationWeightsDir)s" % CONSTANTS],
                                  "DirOfISOMuonWeights": ["%(MuonIsolationWeightsDir)s" % CONSTANTS],
                                  "InputPandoraPFOsCollection": ["PandoraPFOs"],
                                  "InputPrimaryVertexCollection": ["PrimaryVertex"],
                                  "IsSelectingOneIsoLep": ["false"],
                                  "MaxD0SigForElectron": ["50"],
                                  "MaxD0SigForMuon": ["20"],
                                  "MaxEOverPForElectron": ["1.3"],
                                  "MaxEOverPForMuon": ["0.3"],
                                  "MaxZ0SigForElectron": ["50"],
                                  "MaxZ0SigForMuon": ["20"],
                                  "MinEOverPForElectron": ["0.5"],
                                  "MinEecalOverTotEForElectron": ["0.9"],
                                  "MinEyokeForMuon": ["1.2"],
                                  "MinPForElectron": ["5"],
                                  "MinPForMuon": ["5"],
                                  "OutputIsoLeptonsCollection": ["IsolatedMuons"],
                                  "OutputPFOsWithoutIsoLepCollection": ["PFOsminusmu"],
                                  "UseIP": ["true"],
                                  "UseYokeForMuonID": ["true"]
                                  }

IsolatedElectronTagging = MarlinProcessorWrapper("IsolatedElectronTagging")
IsolatedElectronTagging.OutputLevel = ERROR
IsolatedElectronTagging.ProcessorType = "IsolatedLeptonTaggingProcessor"
IsolatedElectronTagging.Parameters = {
                                      "CosConeLarge": ["0.95"],
                                      "CosConeSmall": ["0.98"],
                                      "CutOnTheISOElectronMVA": ["0.5"],
                                      "CutOnTheISOMuonMVA": ["2"],
                                      "DirOfISOElectronWeights": ["%(ElectronIsolationWeightsDir)s" % CONSTANTS],
                                      "DirOfISOMuonWeights": ["%(MuonIsolationWeightsDir)s" % CONSTANTS],
                                      "InputPandoraPFOsCollection": ["PFOsminusmu"],
                                      "InputPrimaryVertexCollection": ["PrimaryVertex"],
                                      "IsSelectingOneIsoLep": ["false"],
                                      "MaxD0SigForElectron": ["10"],
                                      "MaxD0SigForMuon": ["10"],
                                      "MaxEOverPForElectron": ["1.3"],
                                      "MaxEOverPForMuon": ["0.3"],
                                      "MaxZ0SigForElectron": ["10"],
                                      "MaxZ0SigForMuon": ["10"],
                                      "MinEOverPForElectron": ["0.5"],
                                      "MinEecalOverTotEForElectron": ["0.9"],
                                      "MinEyokeForMuon": ["1.2"],
                                      "MinPForElectron": ["5"],
                                      "MinPForMuon": ["5"],
                                      "OutputIsoLeptonsCollection": ["IsolatedElectrons"],
                                      "OutputPFOsWithoutIsoLepCollection": ["PFOsminuse"],
                                      "UseIP": ["true"],
                                      "UseYokeForMuonID": ["true"]
                                      }

IsolatedTauTagging = MarlinProcessorWrapper("IsolatedTauTagging")
IsolatedTauTagging.OutputLevel = ERROR
IsolatedTauTagging.ProcessorType = "TaJetClustering"
IsolatedTauTagging.Parameters = {
                                 "AcceptFlexibleLowEnergyTrack": ["1"],
                                 "ConeMaxCosAngle": ["1"],
                                 "ConeMaxEnergyFrac": ["0.1"],
                                 "ConeMinCosAngle": ["0.9"],
                                 "MinimumJetEnergy": ["3"],
                                 "MinimumTrackEnergy": ["2"],
                                 "MinimumTrackEnergyAssoc": ["2"],
                                 "NoSelection": ["0"],
                                 "OutputTauCollection": ["IsolatedTaus"],
                                 "PFOCollection": ["PFOsminuse"],
                                 "RemainPFOCollection": ["PFOsminustau"],
                                 "TauCosAngle": ["0.98"],
                                 "TauMass": ["2"]
                                 }

IsolatedPhotonTagging = MarlinProcessorWrapper("IsolatedPhotonTagging")
IsolatedPhotonTagging.OutputLevel = ERROR
IsolatedPhotonTagging.ProcessorType = "IsolatedPhotonTaggingProcessor"
IsolatedPhotonTagging.Parameters = {
                                    "CosConeLarge": ["0.95"],
                                    "CosConeSmall": ["0.98"],
                                    "CutOnTheISOPhotonMVA": ["0.5"],
                                    "DirOfIsoPhotonWeights": ["isolated_photon_weights"],
                                    "InputPandoraPFOsCollection": ["PFOsminustau"],
                                    "IsSelectingOneIsoPhoton": ["false"],
                                    "MinEForPhoton": ["5."],
                                    "OutputIsoPhotonsCollection": ["IsolatedPhotons"],
                                    "OutputPFOsWithoutIsoPhotonCollection": ["PFOsminusphoton"]
                                    }

JC2FT = MarlinProcessorWrapper("JC2FT")
JC2FT.OutputLevel = ERROR
JC2FT.ProcessorType = "LcfiplusProcessor"
JC2FT.Parameters = {
                    "Algorithms": ["JetClustering", "JetVertexRefiner", "FlavorTag", "ReadMVA"],
                    "FlavorTag.BookName": ["bdt"],
                    "FlavorTag.CategoryDefinition1": ["nvtx==0"],
                    "FlavorTag.CategoryDefinition2": ["nvtx==1&&nvtxall==1"],
                    "FlavorTag.CategoryDefinition3": ["nvtx==1&&nvtxall==2"],
                    "FlavorTag.CategoryDefinition4": ["nvtx>=2"],
                    "FlavorTag.CategoryPreselection1": ["trk1d0sig!=0"],
                    "FlavorTag.CategoryPreselection2": ["trk1d0sig!=0"],
                    "FlavorTag.CategoryPreselection3": ["trk1d0sig!=0"],
                    "FlavorTag.CategoryPreselection4": ["trk1d0sig!=0"],
                    "FlavorTag.CategorySpectators1": ["aux", "nvtx"],
                    "FlavorTag.CategorySpectators2": ["aux", "nvtx"],
                    "FlavorTag.CategorySpectators3": ["aux", "nvtx"],
                    "FlavorTag.CategorySpectators4": ["aux", "nvtx"],
                    "FlavorTag.CategoryVariables1": ["trk1d0sig", "trk2d0sig", "trk1z0sig", "trk2z0sig", "trk1pt_jete", "trk2pt_jete", "jprobr25sigma", "jprobz25sigma", "d0bprob2", "d0cprob2", "d0qprob2", "z0bprob2", "z0cprob2", "z0qprob2", "nmuon", "nelectron", "trkmass"],
                    "FlavorTag.CategoryVariables2": ["trk1d0sig", "trk2d0sig", "trk1z0sig", "trk2z0sig", "trk1pt_jete", "trk2pt_jete", "jprobr2", "jprobz2", "vtxlen1_jete", "vtxsig1_jete", "vtxdirang1_jete", "vtxmom1_jete", "vtxmass1", "vtxmult1", "vtxmasspc", "vtxprob", "d0bprob2", "d0cprob2", "d0qprob2", "z0bprob2", "z0cprob2", "z0qprob2", "trkmass", "nelectron", "nmuon"],
                    "FlavorTag.CategoryVariables3": ["trk1d0sig", "trk2d0sig", "trk1z0sig", "trk2z0sig", "trk1pt_jete", "trk2pt_jete", "jprobr2", "jprobz2", "vtxlen1_jete", "vtxsig1_jete", "vtxdirang1_jete", "vtxmom1_jete", "vtxmass1", "vtxmult1", "vtxmasspc", "vtxprob", "1vtxprob", "vtxlen12all_jete", "vtxmassall"],
                    "FlavorTag.CategoryVariables4": ["trk1d0sig", "trk2d0sig", "trk1z0sig", "trk2z0sig", "trk1pt_jete", "trk2pt_jete", "jprobr2", "jprobz2", "vtxlen1_jete", "vtxsig1_jete", "vtxdirang1_jete", "vtxmom1_jete", "vtxmass1", "vtxmult1", "vtxmasspc", "vtxprob", "vtxlen2_jete", "vtxsig2_jete", "vtxdirang2_jete", "vtxmom2_jete", "vtxmass2", "vtxmult2", "vtxlen12_jete", "vtxsig12_jete", "vtxdirang12_jete", "vtxmom_jete", "vtxmass", "vtxmult", "1vtxprob"],
                    "FlavorTag.D0ProbFileName": ["%(LCFIPlusD0ProbFile)s" % CONSTANTS],
                    "FlavorTag.JetCollectionName": ["Refined2Jets"],
                    "FlavorTag.PIDAlgo": ["lcfiplus"],
                    "FlavorTag.WeightsDirectory": ["%(LCFIPlusWeightsDir)s" % CONSTANTS],
                    "FlavorTag.WeightsPrefix": ["%(LCFIPlusWeightsPrefix)s" % CONSTANTS],
                    "FlavorTag.Z0ProbFileName": ["%(LCFIPlusZ0ProbFile)s" % CONSTANTS],
                    "JetClustering.InputVertexCollectionName": ["BuildUpVertex"],
                    "JetClustering.JetAlgorithm": ["DurhamVertex"],
                    "JetClustering.MuonIDExternal": ["0"],
                    "JetClustering.MuonIDMaximum3DImpactParameter": ["5."],
                    "JetClustering.MuonIDMinimumD0Significance": ["5."],
                    "JetClustering.MuonIDMinimumProbability": ["0.5"],
                    "JetClustering.MuonIDMinimumZ0Significance": ["5."],
                    "JetClustering.NJetsRequested": ["2"],
                    "JetClustering.OutputJetCollectionName": ["Vertex2Jets"],
                    "JetClustering.PrimaryVertexCollectionName": ["PrimaryVertex"],
                    "JetClustering.UseBeamJets": ["0"],
                    "JetClustering.UseMuonID": ["1"],
                    "JetClustering.VertexSelectionK0MassWidth": ["0.02"],
                    "JetClustering.VertexSelectionMaximumDistance": ["30."],
                    "JetClustering.VertexSelectionMinimumDistance": ["0.3"],
                    "JetClustering.YCut": ["0."],
                    "JetVertexRefiner.InputJetCollectionName": ["Vertex2Jets"],
                    "JetVertexRefiner.InputVertexCollectionName": ["BuildUpVertex"],
                    "JetVertexRefiner.MaxAngleSingle": ["0.5"],
                    "JetVertexRefiner.MaxCharmFlightLengthPerJetEnergy": ["0.1"],
                    "JetVertexRefiner.MaxPosSingle": ["30."],
                    "JetVertexRefiner.MaxSeparationPerPosSingle": ["0.1"],
                    "JetVertexRefiner.MinEnergySingle": ["1."],
                    "JetVertexRefiner.MinPosSingle": ["0.3"],
                    "JetVertexRefiner.OneVertexProbThreshold": ["0.001"],
                    "JetVertexRefiner.OutputJetCollectionName": ["Refined2Jets"],
                    "JetVertexRefiner.OutputVertexCollectionName": ["RefinedVertex2Jets"],
                    "JetVertexRefiner.PrimaryVertexCollectionName": ["PrimaryVertex"],
                    "JetVertexRefiner.V0VertexCollectionName": ["BuildUpVertex_V0"],
                    "JetVertexRefiner.mind0sigSingle": ["5."],
                    "JetVertexRefiner.minz0sigSingle": ["5."],
                    "MCPCollection": [],
                    "MCPFORelation": [],
                    "MakeNtuple.AuxiliaryInfo": ["-1"],
                    "PFOCollection": ["PFOsminusphoton"],
                    "PrimaryVertexCollectionName": ["PrimaryVertex"],
                    "PrintEventNumber": ["0"],
                    "ReadSubdetectorEnergies": ["1"],
                    "TrackHitOrdering": ["1"],
                    "UpdateVertexRPDaughters": ["0"],
                    "UseMCP": ["0"]
                    }

EF2 = MarlinProcessorWrapper("EF2")
EF2.OutputLevel = ERROR
EF2.ProcessorType = "ErrorFlow"
EF2.Parameters = {
                  "InputMCTruthLinkCollection": ["RecoMCTruthLink"],
                  "InputPFOCollection": ["Refined2Jets"],
                  "OutputPFOCollection": ["Refined2JetsEF"]
                  }

JC3FT = MarlinProcessorWrapper("JC3FT")
JC3FT.OutputLevel = ERROR
JC3FT.ProcessorType = "LcfiplusProcessor"
JC3FT.Parameters = {
                    "Algorithms": ["JetClustering", "JetVertexRefiner", "FlavorTag", "ReadMVA"],
                    "FlavorTag.BookName": ["bdt"],
                    "FlavorTag.CategoryDefinition1": ["nvtx==0"],
                    "FlavorTag.CategoryDefinition2": ["nvtx==1&&nvtxall==1"],
                    "FlavorTag.CategoryDefinition3": ["nvtx==1&&nvtxall==2"],
                    "FlavorTag.CategoryDefinition4": ["nvtx>=2"],
                    "FlavorTag.CategoryPreselection1": ["trk1d0sig!=0"],
                    "FlavorTag.CategoryPreselection2": ["trk1d0sig!=0"],
                    "FlavorTag.CategoryPreselection3": ["trk1d0sig!=0"],
                    "FlavorTag.CategoryPreselection4": ["trk1d0sig!=0"],
                    "FlavorTag.CategorySpectators1": ["aux", "nvtx"],
                    "FlavorTag.CategorySpectators2": ["aux", "nvtx"],
                    "FlavorTag.CategorySpectators3": ["aux", "nvtx"],
                    "FlavorTag.CategorySpectators4": ["aux", "nvtx"],
                    "FlavorTag.CategoryVariables1": ["trk1d0sig", "trk2d0sig", "trk1z0sig", "trk2z0sig", "trk1pt_jete", "trk2pt_jete", "jprobr25sigma", "jprobz25sigma", "d0bprob2", "d0cprob2", "d0qprob2", "z0bprob2", "z0cprob2", "z0qprob2", "nmuon", "nelectron", "trkmass"],
                    "FlavorTag.CategoryVariables2": ["trk1d0sig", "trk2d0sig", "trk1z0sig", "trk2z0sig", "trk1pt_jete", "trk2pt_jete", "jprobr2", "jprobz2", "vtxlen1_jete", "vtxsig1_jete", "vtxdirang1_jete", "vtxmom1_jete", "vtxmass1", "vtxmult1", "vtxmasspc", "vtxprob", "d0bprob2", "d0cprob2", "d0qprob2", "z0bprob2", "z0cprob2", "z0qprob2", "trkmass", "nelectron", "nmuon"],
                    "FlavorTag.CategoryVariables3": ["trk1d0sig", "trk2d0sig", "trk1z0sig", "trk2z0sig", "trk1pt_jete", "trk2pt_jete", "jprobr2", "jprobz2", "vtxlen1_jete", "vtxsig1_jete", "vtxdirang1_jete", "vtxmom1_jete", "vtxmass1", "vtxmult1", "vtxmasspc", "vtxprob", "1vtxprob", "vtxlen12all_jete", "vtxmassall"],
                    "FlavorTag.CategoryVariables4": ["trk1d0sig", "trk2d0sig", "trk1z0sig", "trk2z0sig", "trk1pt_jete", "trk2pt_jete", "jprobr2", "jprobz2", "vtxlen1_jete", "vtxsig1_jete", "vtxdirang1_jete", "vtxmom1_jete", "vtxmass1", "vtxmult1", "vtxmasspc", "vtxprob", "vtxlen2_jete", "vtxsig2_jete", "vtxdirang2_jete", "vtxmom2_jete", "vtxmass2", "vtxmult2", "vtxlen12_jete", "vtxsig12_jete", "vtxdirang12_jete", "vtxmom_jete", "vtxmass", "vtxmult", "1vtxprob"],
                    "FlavorTag.D0ProbFileName": ["%(LCFIPlusD0ProbFile)s" % CONSTANTS],
                    "FlavorTag.JetCollectionName": ["Refined3Jets"],
                    "FlavorTag.PIDAlgo": ["lcfiplus"],
                    "FlavorTag.WeightsDirectory": ["%(LCFIPlusWeightsDir)s" % CONSTANTS],
                    "FlavorTag.WeightsPrefix": ["%(LCFIPlusWeightsPrefix)s" % CONSTANTS],
                    "FlavorTag.Z0ProbFileName": ["%(LCFIPlusZ0ProbFile)s" % CONSTANTS],
                    "JetClustering.InputVertexCollectionName": ["BuildUpVertex"],
                    "JetClustering.JetAlgorithm": ["DurhamVertex"],
                    "JetClustering.MuonIDExternal": ["0"],
                    "JetClustering.MuonIDMaximum3DImpactParameter": ["5."],
                    "JetClustering.MuonIDMinimumD0Significance": ["5."],
                    "JetClustering.MuonIDMinimumProbability": ["0.5"],
                    "JetClustering.MuonIDMinimumZ0Significance": ["5."],
                    "JetClustering.NJetsRequested": ["3"],
                    "JetClustering.OutputJetCollectionName": ["Vertex3Jets"],
                    "JetClustering.PrimaryVertexCollectionName": ["PrimaryVertex"],
                    "JetClustering.UseBeamJets": ["0"],
                    "JetClustering.UseMuonID": ["1"],
                    "JetClustering.VertexSelectionK0MassWidth": ["0.02"],
                    "JetClustering.VertexSelectionMaximumDistance": ["30."],
                    "JetClustering.VertexSelectionMinimumDistance": ["0.3"],
                    "JetClustering.YAddedForJetLeptonLepton": ["100"],
                    "JetClustering.YAddedForJetLeptonVertex": ["100"],
                    "JetClustering.YAddedForJetVertexVertex": ["100"],
                    "JetClustering.YCut": ["0."],
                    "JetVertexRefiner.InputJetCollectionName": ["Vertex3Jets"],
                    "JetVertexRefiner.InputVertexCollectionName": ["BuildUpVertex"],
                    "JetVertexRefiner.MaxAngleSingle": ["0.5"],
                    "JetVertexRefiner.MaxCharmFlightLengthPerJetEnergy": ["0.1"],
                    "JetVertexRefiner.MaxPosSingle": ["30."],
                    "JetVertexRefiner.MaxSeparationPerPosSingle": ["0.1"],
                    "JetVertexRefiner.MinEnergySingle": ["1."],
                    "JetVertexRefiner.MinPosSingle": ["0.3"],
                    "JetVertexRefiner.OneVertexProbThreshold": ["0.001"],
                    "JetVertexRefiner.OutputJetCollectionName": ["Refined3Jets"],
                    "JetVertexRefiner.OutputVertexCollectionName": ["RefinedVertex3Jets"],
                    "JetVertexRefiner.PrimaryVertexCollectionName": ["PrimaryVertex"],
                    "JetVertexRefiner.V0VertexCollectionName": ["BuildUpVertex_V0"],
                    "JetVertexRefiner.mind0sigSingle": ["5."],
                    "JetVertexRefiner.minz0sigSingle": ["5."],
                    "MCPCollection": [],
                    "MCPFORelation": [],
                    "MakeNtuple.AuxiliaryInfo": ["-1"],
                    "PFOCollection": ["PFOsminusphoton"],
                    "PrimaryVertexCollectionName": ["PrimaryVertex"],
                    "PrintEventNumber": ["0"],
                    "ReadSubdetectorEnergies": ["1"],
                    "TrackHitOrdering": ["1"],
                    "UpdateVertexRPDaughters": ["0"],
                    "UseMCP": ["0"]
                    }

EF3 = MarlinProcessorWrapper("EF3")
EF3.OutputLevel = ERROR
EF3.ProcessorType = "ErrorFlow"
EF3.Parameters = {
                  "InputMCTruthLinkCollection": ["RecoMCTruthLink"],
                  "InputPFOCollection": ["Refined3Jets"],
                  "OutputPFOCollection": ["Refined3JetsEF"]
                  }

JC4FT = MarlinProcessorWrapper("JC4FT")
JC4FT.OutputLevel = ERROR
JC4FT.ProcessorType = "LcfiplusProcessor"
JC4FT.Parameters = {
                    "Algorithms": ["JetClustering", "JetVertexRefiner", "FlavorTag", "ReadMVA"],
                    "FlavorTag.BookName": ["bdt"],
                    "FlavorTag.CategoryDefinition1": ["nvtx==0"],
                    "FlavorTag.CategoryDefinition2": ["nvtx==1&&nvtxall==1"],
                    "FlavorTag.CategoryDefinition3": ["nvtx==1&&nvtxall==2"],
                    "FlavorTag.CategoryDefinition4": ["nvtx>=2"],
                    "FlavorTag.CategoryPreselection1": ["trk1d0sig!=0"],
                    "FlavorTag.CategoryPreselection2": ["trk1d0sig!=0"],
                    "FlavorTag.CategoryPreselection3": ["trk1d0sig!=0"],
                    "FlavorTag.CategoryPreselection4": ["trk1d0sig!=0"],
                    "FlavorTag.CategorySpectators1": ["aux", "nvtx"],
                    "FlavorTag.CategorySpectators2": ["aux", "nvtx"],
                    "FlavorTag.CategorySpectators3": ["aux", "nvtx"],
                    "FlavorTag.CategorySpectators4": ["aux", "nvtx"],
                    "FlavorTag.CategoryVariables1": ["trk1d0sig", "trk2d0sig", "trk1z0sig", "trk2z0sig", "trk1pt_jete", "trk2pt_jete", "jprobr25sigma", "jprobz25sigma", "d0bprob2", "d0cprob2", "d0qprob2", "z0bprob2", "z0cprob2", "z0qprob2", "nmuon", "nelectron", "trkmass"],
                    "FlavorTag.CategoryVariables2": ["trk1d0sig", "trk2d0sig", "trk1z0sig", "trk2z0sig", "trk1pt_jete", "trk2pt_jete", "jprobr2", "jprobz2", "vtxlen1_jete", "vtxsig1_jete", "vtxdirang1_jete", "vtxmom1_jete", "vtxmass1", "vtxmult1", "vtxmasspc", "vtxprob", "d0bprob2", "d0cprob2", "d0qprob2", "z0bprob2", "z0cprob2", "z0qprob2", "trkmass", "nelectron", "nmuon"],
                    "FlavorTag.CategoryVariables3": ["trk1d0sig", "trk2d0sig", "trk1z0sig", "trk2z0sig", "trk1pt_jete", "trk2pt_jete", "jprobr2", "jprobz2", "vtxlen1_jete", "vtxsig1_jete", "vtxdirang1_jete", "vtxmom1_jete", "vtxmass1", "vtxmult1", "vtxmasspc", "vtxprob", "1vtxprob", "vtxlen12all_jete", "vtxmassall"],
                    "FlavorTag.CategoryVariables4": ["trk1d0sig", "trk2d0sig", "trk1z0sig", "trk2z0sig", "trk1pt_jete", "trk2pt_jete", "jprobr2", "jprobz2", "vtxlen1_jete", "vtxsig1_jete", "vtxdirang1_jete", "vtxmom1_jete", "vtxmass1", "vtxmult1", "vtxmasspc", "vtxprob", "vtxlen2_jete", "vtxsig2_jete", "vtxdirang2_jete", "vtxmom2_jete", "vtxmass2", "vtxmult2", "vtxlen12_jete", "vtxsig12_jete", "vtxdirang12_jete", "vtxmom_jete", "vtxmass", "vtxmult", "1vtxprob"],
                    "FlavorTag.D0ProbFileName": ["%(LCFIPlusD0ProbFile)s" % CONSTANTS],
                    "FlavorTag.JetCollectionName": ["Refined4Jets"],
                    "FlavorTag.PIDAlgo": ["lcfiplus"],
                    "FlavorTag.WeightsDirectory": ["%(LCFIPlusWeightsDir)s" % CONSTANTS],
                    "FlavorTag.WeightsPrefix": ["%(LCFIPlusWeightsPrefix)s" % CONSTANTS],
                    "FlavorTag.Z0ProbFileName": ["%(LCFIPlusZ0ProbFile)s" % CONSTANTS],
                    "JetClustering.InputVertexCollectionName": ["BuildUpVertex"],
                    "JetClustering.JetAlgorithm": ["DurhamVertex"],
                    "JetClustering.MuonIDExternal": ["0"],
                    "JetClustering.MuonIDMaximum3DImpactParameter": ["5."],
                    "JetClustering.MuonIDMinimumD0Significance": ["5."],
                    "JetClustering.MuonIDMinimumProbability": ["0.5"],
                    "JetClustering.MuonIDMinimumZ0Significance": ["5."],
                    "JetClustering.NJetsRequested": ["4"],
                    "JetClustering.OutputJetCollectionName": ["Vertex4Jets"],
                    "JetClustering.PrimaryVertexCollectionName": ["PrimaryVertex"],
                    "JetClustering.UseBeamJets": ["0"],
                    "JetClustering.UseMuonID": ["1"],
                    "JetClustering.VertexSelectionK0MassWidth": ["0.02"],
                    "JetClustering.VertexSelectionMaximumDistance": ["30."],
                    "JetClustering.VertexSelectionMinimumDistance": ["0.3"],
                    "JetClustering.YAddedForJetLeptonLepton": ["100"],
                    "JetClustering.YAddedForJetLeptonVertex": ["100"],
                    "JetClustering.YAddedForJetVertexVertex": ["100"],
                    "JetClustering.YCut": ["0."],
                    "JetVertexRefiner.InputJetCollectionName": ["Vertex4Jets"],
                    "JetVertexRefiner.InputVertexCollectionName": ["BuildUpVertex"],
                    "JetVertexRefiner.MaxAngleSingle": ["0.5"],
                    "JetVertexRefiner.MaxCharmFlightLengthPerJetEnergy": ["0.1"],
                    "JetVertexRefiner.MaxPosSingle": ["30."],
                    "JetVertexRefiner.MaxSeparationPerPosSingle": ["0.1"],
                    "JetVertexRefiner.MinEnergySingle": ["1."],
                    "JetVertexRefiner.MinPosSingle": ["0.3"],
                    "JetVertexRefiner.OneVertexProbThreshold": ["0.001"],
                    "JetVertexRefiner.OutputJetCollectionName": ["Refined4Jets"],
                    "JetVertexRefiner.OutputVertexCollectionName": ["RefinedVertex4Jets"],
                    "JetVertexRefiner.PrimaryVertexCollectionName": ["PrimaryVertex"],
                    "JetVertexRefiner.V0VertexCollectionName": ["BuildUpVertex_V0"],
                    "JetVertexRefiner.mind0sigSingle": ["5."],
                    "JetVertexRefiner.minz0sigSingle": ["5."],
                    "MCPCollection": [],
                    "MCPFORelation": [],
                    "MakeNtuple.AuxiliaryInfo": ["-1"],
                    "PFOCollection": ["PFOsminusphoton"],
                    "PrimaryVertexCollectionName": ["PrimaryVertex"],
                    "PrintEventNumber": ["0"],
                    "ReadSubdetectorEnergies": ["1"],
                    "TrackHitOrdering": ["1"],
                    "UpdateVertexRPDaughters": ["0"],
                    "UseMCP": ["0"]
                    }

EF4 = MarlinProcessorWrapper("EF4")
EF4.OutputLevel = ERROR
EF4.ProcessorType = "ErrorFlow"
EF4.Parameters = {
                  "InputMCTruthLinkCollection": ["RecoMCTruthLink"],
                  "InputPFOCollection": ["Refined4Jets"],
                  "OutputPFOCollection": ["Refined4JetsEF"]
                  }

JC5FT = MarlinProcessorWrapper("JC5FT")
JC5FT.OutputLevel = ERROR
JC5FT.ProcessorType = "LcfiplusProcessor"
JC5FT.Parameters = {
                    "Algorithms": ["JetClustering", "JetVertexRefiner", "FlavorTag", "ReadMVA"],
                    "FlavorTag.BookName": ["bdt"],
                    "FlavorTag.CategoryDefinition1": ["nvtx==0"],
                    "FlavorTag.CategoryDefinition2": ["nvtx==1&&nvtxall==1"],
                    "FlavorTag.CategoryDefinition3": ["nvtx==1&&nvtxall==2"],
                    "FlavorTag.CategoryDefinition4": ["nvtx>=2"],
                    "FlavorTag.CategoryPreselection1": ["trk1d0sig!=0"],
                    "FlavorTag.CategoryPreselection2": ["trk1d0sig!=0"],
                    "FlavorTag.CategoryPreselection3": ["trk1d0sig!=0"],
                    "FlavorTag.CategoryPreselection4": ["trk1d0sig!=0"],
                    "FlavorTag.CategorySpectators1": ["aux", "nvtx"],
                    "FlavorTag.CategorySpectators2": ["aux", "nvtx"],
                    "FlavorTag.CategorySpectators3": ["aux", "nvtx"],
                    "FlavorTag.CategorySpectators4": ["aux", "nvtx"],
                    "FlavorTag.CategoryVariables1": ["trk1d0sig", "trk2d0sig", "trk1z0sig", "trk2z0sig", "trk1pt_jete", "trk2pt_jete", "jprobr25sigma", "jprobz25sigma", "d0bprob2", "d0cprob2", "d0qprob2", "z0bprob2", "z0cprob2", "z0qprob2", "nmuon", "nelectron", "trkmass"],
                    "FlavorTag.CategoryVariables2": ["trk1d0sig", "trk2d0sig", "trk1z0sig", "trk2z0sig", "trk1pt_jete", "trk2pt_jete", "jprobr2", "jprobz2", "vtxlen1_jete", "vtxsig1_jete", "vtxdirang1_jete", "vtxmom1_jete", "vtxmass1", "vtxmult1", "vtxmasspc", "vtxprob", "d0bprob2", "d0cprob2", "d0qprob2", "z0bprob2", "z0cprob2", "z0qprob2", "trkmass", "nelectron", "nmuon"],
                    "FlavorTag.CategoryVariables3": ["trk1d0sig", "trk2d0sig", "trk1z0sig", "trk2z0sig", "trk1pt_jete", "trk2pt_jete", "jprobr2", "jprobz2", "vtxlen1_jete", "vtxsig1_jete", "vtxdirang1_jete", "vtxmom1_jete", "vtxmass1", "vtxmult1", "vtxmasspc", "vtxprob", "1vtxprob", "vtxlen12all_jete", "vtxmassall"],
                    "FlavorTag.CategoryVariables4": ["trk1d0sig", "trk2d0sig", "trk1z0sig", "trk2z0sig", "trk1pt_jete", "trk2pt_jete", "jprobr2", "jprobz2", "vtxlen1_jete", "vtxsig1_jete", "vtxdirang1_jete", "vtxmom1_jete", "vtxmass1", "vtxmult1", "vtxmasspc", "vtxprob", "vtxlen2_jete", "vtxsig2_jete", "vtxdirang2_jete", "vtxmom2_jete", "vtxmass2", "vtxmult2", "vtxlen12_jete", "vtxsig12_jete", "vtxdirang12_jete", "vtxmom_jete", "vtxmass", "vtxmult", "1vtxprob"],
                    "FlavorTag.D0ProbFileName": ["%(LCFIPlusD0ProbFile)s" % CONSTANTS],
                    "FlavorTag.JetCollectionName": ["Refined5Jets"],
                    "FlavorTag.PIDAlgo": ["lcfiplus"],
                    "FlavorTag.WeightsDirectory": ["%(LCFIPlusWeightsDir)s" % CONSTANTS],
                    "FlavorTag.WeightsPrefix": ["%(LCFIPlusWeightsPrefix)s" % CONSTANTS],
                    "FlavorTag.Z0ProbFileName": ["%(LCFIPlusZ0ProbFile)s" % CONSTANTS],
                    "JetClustering.InputVertexCollectionName": ["BuildUpVertex"],
                    "JetClustering.JetAlgorithm": ["DurhamVertex"],
                    "JetClustering.MuonIDExternal": ["0"],
                    "JetClustering.MuonIDMaximum3DImpactParameter": ["5."],
                    "JetClustering.MuonIDMinimumD0Significance": ["5."],
                    "JetClustering.MuonIDMinimumProbability": ["0.5"],
                    "JetClustering.MuonIDMinimumZ0Significance": ["5."],
                    "JetClustering.NJetsRequested": ["5"],
                    "JetClustering.OutputJetCollectionName": ["Vertex5Jets"],
                    "JetClustering.PrimaryVertexCollectionName": ["PrimaryVertex"],
                    "JetClustering.UseBeamJets": ["0"],
                    "JetClustering.UseMuonID": ["1"],
                    "JetClustering.VertexSelectionK0MassWidth": ["0.02"],
                    "JetClustering.VertexSelectionMaximumDistance": ["30."],
                    "JetClustering.VertexSelectionMinimumDistance": ["0.3"],
                    "JetClustering.YAddedForJetLeptonLepton": ["100"],
                    "JetClustering.YAddedForJetLeptonVertex": ["100"],
                    "JetClustering.YAddedForJetVertexVertex": ["100"],
                    "JetClustering.YCut": ["0."],
                    "JetVertexRefiner.InputJetCollectionName": ["Vertex5Jets"],
                    "JetVertexRefiner.InputVertexCollectionName": ["BuildUpVertex"],
                    "JetVertexRefiner.MaxAngleSingle": ["0.5"],
                    "JetVertexRefiner.MaxCharmFlightLengthPerJetEnergy": ["0.1"],
                    "JetVertexRefiner.MaxPosSingle": ["30."],
                    "JetVertexRefiner.MaxSeparationPerPosSingle": ["0.1"],
                    "JetVertexRefiner.MinEnergySingle": ["1."],
                    "JetVertexRefiner.MinPosSingle": ["0.3"],
                    "JetVertexRefiner.OneVertexProbThreshold": ["0.001"],
                    "JetVertexRefiner.OutputJetCollectionName": ["Refined5Jets"],
                    "JetVertexRefiner.OutputVertexCollectionName": ["RefinedVertex5Jets"],
                    "JetVertexRefiner.PrimaryVertexCollectionName": ["PrimaryVertex"],
                    "JetVertexRefiner.V0VertexCollectionName": ["BuildUpVertex_V0"],
                    "JetVertexRefiner.mind0sigSingle": ["5."],
                    "JetVertexRefiner.minz0sigSingle": ["5."],
                    "MCPCollection": [],
                    "MCPFORelation": [],
                    "MakeNtuple.AuxiliaryInfo": ["-1"],
                    "PFOCollection": ["PFOsminusphoton"],
                    "PrimaryVertexCollectionName": ["PrimaryVertex"],
                    "PrintEventNumber": ["0"],
                    "ReadSubdetectorEnergies": ["1"],
                    "TrackHitOrdering": ["1"],
                    "UpdateVertexRPDaughters": ["0"],
                    "UseMCP": ["0"]
                    }

EF5 = MarlinProcessorWrapper("EF5")
EF5.OutputLevel = ERROR
EF5.ProcessorType = "ErrorFlow"
EF5.Parameters = {
                  "InputMCTruthLinkCollection": ["RecoMCTruthLink"],
                  "InputPFOCollection": ["Refined5Jets"],
                  "OutputPFOCollection": ["Refined5JetsEF"]
                  }

JC6FT = MarlinProcessorWrapper("JC6FT")
JC6FT.OutputLevel = ERROR
JC6FT.ProcessorType = "LcfiplusProcessor"
JC6FT.Parameters = {
                    "Algorithms": ["JetClustering", "JetVertexRefiner", "FlavorTag", "ReadMVA"],
                    "FlavorTag.BookName": ["bdt"],
                    "FlavorTag.CategoryDefinition1": ["nvtx==0"],
                    "FlavorTag.CategoryDefinition2": ["nvtx==1&&nvtxall==1"],
                    "FlavorTag.CategoryDefinition3": ["nvtx==1&&nvtxall==2"],
                    "FlavorTag.CategoryDefinition4": ["nvtx>=2"],
                    "FlavorTag.CategoryPreselection1": ["trk1d0sig!=0"],
                    "FlavorTag.CategoryPreselection2": ["trk1d0sig!=0"],
                    "FlavorTag.CategoryPreselection3": ["trk1d0sig!=0"],
                    "FlavorTag.CategoryPreselection4": ["trk1d0sig!=0"],
                    "FlavorTag.CategorySpectators1": ["aux", "nvtx"],
                    "FlavorTag.CategorySpectators2": ["aux", "nvtx"],
                    "FlavorTag.CategorySpectators3": ["aux", "nvtx"],
                    "FlavorTag.CategorySpectators4": ["aux", "nvtx"],
                    "FlavorTag.CategoryVariables1": ["trk1d0sig", "trk2d0sig", "trk1z0sig", "trk2z0sig", "trk1pt_jete", "trk2pt_jete", "jprobr25sigma", "jprobz25sigma", "d0bprob2", "d0cprob2", "d0qprob2", "z0bprob2", "z0cprob2", "z0qprob2", "nmuon", "nelectron", "trkmass"],
                    "FlavorTag.CategoryVariables2": ["trk1d0sig", "trk2d0sig", "trk1z0sig", "trk2z0sig", "trk1pt_jete", "trk2pt_jete", "jprobr2", "jprobz2", "vtxlen1_jete", "vtxsig1_jete", "vtxdirang1_jete", "vtxmom1_jete", "vtxmass1", "vtxmult1", "vtxmasspc", "vtxprob", "d0bprob2", "d0cprob2", "d0qprob2", "z0bprob2", "z0cprob2", "z0qprob2", "trkmass", "nelectron", "nmuon"],
                    "FlavorTag.CategoryVariables3": ["trk1d0sig", "trk2d0sig", "trk1z0sig", "trk2z0sig", "trk1pt_jete", "trk2pt_jete", "jprobr2", "jprobz2", "vtxlen1_jete", "vtxsig1_jete", "vtxdirang1_jete", "vtxmom1_jete", "vtxmass1", "vtxmult1", "vtxmasspc", "vtxprob", "1vtxprob", "vtxlen12all_jete", "vtxmassall"],
                    "FlavorTag.CategoryVariables4": ["trk1d0sig", "trk2d0sig", "trk1z0sig", "trk2z0sig", "trk1pt_jete", "trk2pt_jete", "jprobr2", "jprobz2", "vtxlen1_jete", "vtxsig1_jete", "vtxdirang1_jete", "vtxmom1_jete", "vtxmass1", "vtxmult1", "vtxmasspc", "vtxprob", "vtxlen2_jete", "vtxsig2_jete", "vtxdirang2_jete", "vtxmom2_jete", "vtxmass2", "vtxmult2", "vtxlen12_jete", "vtxsig12_jete", "vtxdirang12_jete", "vtxmom_jete", "vtxmass", "vtxmult", "1vtxprob"],
                    "FlavorTag.D0ProbFileName": ["%(LCFIPlusD0ProbFile)s" % CONSTANTS],
                    "FlavorTag.JetCollectionName": ["Refined6Jets"],
                    "FlavorTag.PIDAlgo": ["lcfiplus"],
                    "FlavorTag.WeightsDirectory": ["%(LCFIPlusWeightsDir)s" % CONSTANTS],
                    "FlavorTag.WeightsPrefix": ["%(LCFIPlusWeightsPrefix)s" % CONSTANTS],
                    "FlavorTag.Z0ProbFileName": ["%(LCFIPlusZ0ProbFile)s" % CONSTANTS],
                    "JetClustering.InputVertexCollectionName": ["BuildUpVertex"],
                    "JetClustering.JetAlgorithm": ["DurhamVertex"],
                    "JetClustering.MuonIDExternal": ["0"],
                    "JetClustering.MuonIDMaximum3DImpactParameter": ["5."],
                    "JetClustering.MuonIDMinimumD0Significance": ["5."],
                    "JetClustering.MuonIDMinimumProbability": ["0.5"],
                    "JetClustering.MuonIDMinimumZ0Significance": ["5."],
                    "JetClustering.NJetsRequested": ["6"],
                    "JetClustering.OutputJetCollectionName": ["Vertex6Jets"],
                    "JetClustering.PrimaryVertexCollectionName": ["PrimaryVertex"],
                    "JetClustering.UseBeamJets": ["0"],
                    "JetClustering.UseMuonID": ["1"],
                    "JetClustering.VertexSelectionK0MassWidth": ["0.02"],
                    "JetClustering.VertexSelectionMaximumDistance": ["30."],
                    "JetClustering.VertexSelectionMinimumDistance": ["0.3"],
                    "JetClustering.YAddedForJetLeptonLepton": ["100"],
                    "JetClustering.YAddedForJetLeptonVertex": ["100"],
                    "JetClustering.YAddedForJetVertexVertex": ["100"],
                    "JetClustering.YCut": ["0."],
                    "JetVertexRefiner.InputJetCollectionName": ["Vertex6Jets"],
                    "JetVertexRefiner.InputVertexCollectionName": ["BuildUpVertex"],
                    "JetVertexRefiner.MaxAngleSingle": ["0.5"],
                    "JetVertexRefiner.MaxCharmFlightLengthPerJetEnergy": ["0.1"],
                    "JetVertexRefiner.MaxPosSingle": ["30."],
                    "JetVertexRefiner.MaxSeparationPerPosSingle": ["0.1"],
                    "JetVertexRefiner.MinEnergySingle": ["1."],
                    "JetVertexRefiner.MinPosSingle": ["0.3"],
                    "JetVertexRefiner.OneVertexProbThreshold": ["0.001"],
                    "JetVertexRefiner.OutputJetCollectionName": ["Refined6Jets"],
                    "JetVertexRefiner.OutputVertexCollectionName": ["RefinedVertex6Jets"],
                    "JetVertexRefiner.PrimaryVertexCollectionName": ["PrimaryVertex"],
                    "JetVertexRefiner.V0VertexCollectionName": ["BuildUpVertex_V0"],
                    "JetVertexRefiner.mind0sigSingle": ["5."],
                    "JetVertexRefiner.minz0sigSingle": ["5."],
                    "MCPCollection": [],
                    "MCPFORelation": [],
                    "MakeNtuple.AuxiliaryInfo": ["-1"],
                    "PFOCollection": ["PFOsminusphoton"],
                    "PrimaryVertexCollectionName": ["PrimaryVertex"],
                    "PrintEventNumber": ["0"],
                    "ReadSubdetectorEnergies": ["1"],
                    "TrackHitOrdering": ["1"],
                    "UpdateVertexRPDaughters": ["0"],
                    "UseMCP": ["0"]
                    }

EF6 = MarlinProcessorWrapper("EF6")
EF6.OutputLevel = ERROR
EF6.ProcessorType = "ErrorFlow"
EF6.Parameters = {
                  "InputMCTruthLinkCollection": ["RecoMCTruthLink"],
                  "InputPFOCollection": ["Refined6Jets"],
                  "OutputPFOCollection": ["Refined6JetsEF"]
                  }

ComputeCorrectAngulardEdX = MarlinProcessorWrapper("ComputeCorrectAngulardEdX")
ComputeCorrectAngulardEdX.OutputLevel = ERROR
ComputeCorrectAngulardEdX.ProcessorType = "AngularCorrection_dEdxProcessor"
ComputeCorrectAngulardEdX.Parameters = {
                                        "AngularCorrectionParameters": ["0.970205", "0.0007506", "4.41781e-8", "5.8222e-8"],
                                        "LDCTrackCollection": ["MarlinTrkTracks"]
                                        }

LikelihoodPID = MarlinProcessorWrapper("LikelihoodPID")
LikelihoodPID.OutputLevel = ERROR
LikelihoodPID.ProcessorType = "LikelihoodPIDProcessor"
LikelihoodPID.Parameters = {
                            "CostMatrix": ["1.0e-50", "1.0", "1.5", "1.0", "1.5", "1.0", "1.0e-50", "3.0", "1.0", "1.0", "1.0", "1.0", "1.0e-50", "1.0", "3.0", "1.0", "1.0", "4.0", "1.0e-50", "2.0", "1.0", "1.0", "5.0", "1.0", "1.0e-50"],
                            "Debug": ["0"],
                            "EnergyBoundaries": ["0", "1.0e+07"],
                            "FilePDFName": ["%(ProductionDir)s/HighLevelReco/PIDFiles/LikelihoodPID_Standard_l5_v01.root" % CONSTANTS],
                            "FileWeightFormupiSeparationName": ["%(ProductionDir)s/HighLevelReco/PIDFiles/LowMomMuPiSeparation/TMVAClassification_BDTG_l5_02GeVP_clusterinfo.weights.xml" % CONSTANTS, "%(ProductionDir)s/HighLevelReco/PIDFiles/LowMomMuPiSeparation/TMVAClassification_BDTG_l5_03GeVP_clusterinfo.weights.xml" % CONSTANTS, "%(ProductionDir)s/HighLevelReco/PIDFiles/LowMomMuPiSeparation/TMVAClassification_BDTG_l5_04GeVP_clusterinfo.weights.xml" % CONSTANTS, "%(ProductionDir)s/HighLevelReco/PIDFiles/LowMomMuPiSeparation/TMVAClassification_BDTG_l5_05GeVP_clusterinfo.weights.xml" % CONSTANTS, "%(ProductionDir)s/HighLevelReco/PIDFiles/LowMomMuPiSeparation/TMVAClassification_BDTG_l5_06GeVP_clusterinfo.weights.xml" % CONSTANTS, "%(ProductionDir)s/HighLevelReco/PIDFiles/LowMomMuPiSeparation/TMVAClassification_BDTG_l5_07GeVP_clusterinfo.weights.xml" % CONSTANTS, "%(ProductionDir)s/HighLevelReco/PIDFiles/LowMomMuPiSeparation/TMVAClassification_BDTG_l5_08GeVP_clusterinfo.weights.xml" % CONSTANTS, "%(ProductionDir)s/HighLevelReco/PIDFiles/LowMomMuPiSeparation/TMVAClassification_BDTG_l5_09GeVP_clusterinfo.weights.xml" % CONSTANTS, "%(ProductionDir)s/HighLevelReco/PIDFiles/LowMomMuPiSeparation/TMVAClassification_BDTG_l5_10GeVP_clusterinfo.weights.xml" % CONSTANTS, "%(ProductionDir)s/HighLevelReco/PIDFiles/LowMomMuPiSeparation/TMVAClassification_BDTG_l5_11GeVP_clusterinfo.weights.xml" % CONSTANTS, "%(ProductionDir)s/HighLevelReco/PIDFiles/LowMomMuPiSeparation/TMVAClassification_BDTG_l5_12GeVP_clusterinfo.weights.xml" % CONSTANTS, "%(ProductionDir)s/HighLevelReco/PIDFiles/LowMomMuPiSeparation/TMVAClassification_BDTG_l5_13GeVP_clusterinfo.weights.xml" % CONSTANTS, "%(ProductionDir)s/HighLevelReco/PIDFiles/LowMomMuPiSeparation/TMVAClassification_BDTG_l5_14GeVP_clusterinfo.weights.xml" % CONSTANTS, "%(ProductionDir)s/HighLevelReco/PIDFiles/LowMomMuPiSeparation/TMVAClassification_BDTG_l5_15GeVP_clusterinfo.weights.xml" % CONSTANTS, "%(ProductionDir)s/HighLevelReco/PIDFiles/LowMomMuPiSeparation/TMVAClassification_BDTG_l5_16GeVP_clusterinfo.weights.xml" % CONSTANTS, "%(ProductionDir)s/HighLevelReco/PIDFiles/LowMomMuPiSeparation/TMVAClassification_BDTG_l5_17GeVP_clusterinfo.weights.xml" % CONSTANTS, "%(ProductionDir)s/HighLevelReco/PIDFiles/LowMomMuPiSeparation/TMVAClassification_BDTG_l5_18GeVP_clusterinfo.weights.xml" % CONSTANTS, "%(ProductionDir)s/HighLevelReco/PIDFiles/LowMomMuPiSeparation/TMVAClassification_BDTG_l5_19GeVP_clusterinfo.weights.xml" % CONSTANTS, "%(ProductionDir)s/HighLevelReco/PIDFiles/LowMomMuPiSeparation/TMVAClassification_BDTG_l5_20GeVP_clusterinfo.weights.xml" % CONSTANTS],
                            "PIDMethodsToRun_version": ["v2"],
                            "RecoParticleCollection": ["PandoraPFOs"],
                            "UseBayesian": ["2"],
                            "UseLowMomentumMuPiSeparation": ["true"],
                            "dEdxErrorFactor": ["7.55"],
                            "dEdxNormalization": ["1.350e-7"],
                            "dEdxParameter_electron": ["-0.00232937", "-3.88424e+13", "-37881.1", "-1.56837", "0"],
                            "dEdxParameter_kaon": ["0.0792784", "3798.12", "4.06952e+07", "0.450671", "0.00050169"],
                            "dEdxParameter_muon": ["0.0717375", "-16596.5", "-4.84028e+07", "0.356728", "0.000371431"],
                            "dEdxParameter_pion": ["0.0733683", "51678.4", "8.19644e+07", "0.453505", "0.000404984"],
                            "dEdxParameter_proton": ["0.0770318", "1053.24", "4.95076e+06", "0.281489", "0.000168616"]
                            }

LeptonID = MarlinProcessorWrapper("LeptonID")
LeptonID.OutputLevel = ERROR
LeptonID.ProcessorType = "LeptonIDProcessor"
LeptonID.Parameters = {
                       "BuildTree": ["false"],
                       "EvalMVA": ["true"],
                       "weightfile": ["%(ProductionDir)s/HighLevelReco/PIDFiles/CPID/LeptonID_multi_jet_dEdx_800t_3d_cm50r50_BDTG.weights.xml" % CONSTANTS]
                       }

ComprehensivePID = MarlinProcessorWrapper("ComprehensivePID")
ComprehensivePID.OutputLevel = ERROR
ComprehensivePID.ProcessorType = "ComprehensivePIDProcessor"
ComprehensivePID.Parameters = {
                               "PFOCollection": ["PandoraPFOs"],
                               "RecoMCTruthLink": ["RecoMCTruthLink"],
                               "TMVA_BDT_MC_12bins_singleP.S": ["!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=multiclass", "SplitMode=Random:NormMode=NumEvents:!V", "!V:!H:NTrees=1000:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.50:nCuts=30:MaxDepth=5", "dEdx_RCD_piDis>-900&&dEdx_RCD_kaDis>-900"],
                               "TOF100.S": ["TOFEstimators100ps"],
                               "TTreeFileName": [],
                               "backgroundPDGs": [],
                               "cutD0": ["0"],
                               "cutLamMax": ["0"],
                               "cutLamMin": ["0"],
                               "cutNTracksMax": ["-1"],
                               "cutNTracksMin": ["1"],
                               "cutZ0": ["0"],
                               "dEdx_RCD.F": ["-1.28883368e-02", "2.72959919e+01", "1.10560871e+01", "-1.74534200e+00", "-9.84887586e-07", "6.49143971e-02", "1.55775592e+03", "9.31848047e+08", "2.32201725e-01", "2.50492066e-04", "6.54955215e-02", "8.26239081e+04", "1.92933904e+07", "2.52743206e-01", "2.26657525e-04", "7.52235689e-02", "1.59710415e+04", "1.79625604e+06", "3.15315795e-01", "2.30414997e-04", "7.92251260e-02", "6.38129720e+04", "3.82995071e+04", "2.80793601e-01", "7.14371743e-04", "1"],
                               "fileFormat": [".png"],
                               "inputAlgoSpecs": ["dEdx_RCD:dEdx_RCD", "TOF:TOF100", "Pandora", "LeptonID"],
                               "modeExtract": ["true"],
                               "modeInfer": ["true"],
                               "modeTrain": ["false"],
                               "momLog": ["true"],
                               "momMax": ["100"],
                               "momMin": ["1"],
                               "momNBins": ["12"],
                               "plotFolder": ["CPID_Plots"],
                               "reffile": ["Ref_singleP_16bins_conservative.txt"],
                               "signalPDGs": ["11", "13", "211", "321", "2212"],
                               "trainModelSpecs": ["TMVA_BDT_MC:TMVA_BDT_MC_12bins_singleP"],
                               "trainingObservables": []
                               }

WWCategorisation = MarlinProcessorWrapper("WWCategorisation")
WWCategorisation.OutputLevel = ERROR
WWCategorisation.ProcessorType = "WWCategorisationProcessor"
WWCategorisation.Parameters = {
                               "IsolatedElectrons": ["IsolatedElectrons"],
                               "IsolatedMuons": ["IsolatedMuons"],
                               "IsolatedPhotons": ["IsolatedPhotons"],
                               "IsolatedTaus": ["IsolatedTaus"],
                               "PFOsminusphoton": ["PFOsminusphoton"],
                               }

LCIOOutputProcessor = MarlinProcessorWrapper("LCIOOutputProcessor")
LCIOOutputProcessor.OutputLevel = ERROR
LCIOOutputProcessor.ProcessorType = "LCIOOutputProcessor"
LCIOOutputProcessor.Parameters = {
                                  "CompressionLevel": ["6"],
                                  "DropCollectionNames": ["ClusterMCTruthLink", "MCTruthClusterLink", "MCTruthMarlinTrkTracksLink", "MarlinTrkTracksMCTruthLink", "DistilledPFOs", "GammaGammaCandidateEtaPrimes", "GammaGammaCandidateEtas", "GammaGammaCandidatePi0s", "GammaGammaParticles", "V0RecoParticles", "V0Vertices", "ProngRecoParticles", "ProngVertices", "KinkRecoParticles", "KinkVertices", "SplitRecoParticles", "SplitVertices", "BuildUpVertex_RP", "BuildUpVertex", "BuildUpVertex_V0_RP", "BuildUpVertex_V0", "PFOsminusmu", "PFOsminuse", "PFOsminustau", "PFOsminusphoton", "Vertex2Jets", "Vertex3Jets", "Vertex4Jets", "Vertex5Jets", "Vertex6Jets", "Refined2Jets_rel", "Refined3Jets_rel", "Refined4Jets_rel", "Refined5Jets_rel", "Refined6Jets_rel", "Refined2Jets_vtx", "Refined3Jets_vtx", "Refined4Jets_vtx", "Refined5Jets_vtx", "Refined6Jets_vtx", "Refined2Jets_vtx_RP", "Refined3Jets_vtx_RP", "Refined4Jets_vtx_RP", "Refined5Jets_vtx_RP", "Refined6Jets_vtx_RP", "RefinedVertex2Jets", "RefinedVertex3Jets", "RefinedVertex4Jets", "RefinedVertex5Jets", "RefinedVertex6Jets", "RefinedVertex2Jets_RP", "RefinedVertex3Jets_RP", "RefinedVertex4Jets_RP", "RefinedVertex5Jets_RP", "RefinedVertex6Jets_RP", "Refined2JetsEF", "Refined3JetsEF", "Refined4JetsEF", "Refined5JetsEF", "Refined6JetsEF", "BCALMCTruthLink", "MCTruthBcalLink", "MCTruthTrackLink", "TrackMCTruthLink"],
                                  "DropCollectionTypes": ["Track", "Cluster"],
                                  "LCIOOutputFile": ["%(OutputFile)s" % CONSTANTS],
                                  "LCIOWriteMode": ["WRITE_NEW"]
                                  }

algList.append(Statusmonitor)
algList.append(MyTrueJet)
algList.append(MyCheatedMCOverlayRemoval)
algList.append(Thrust)
algList.append(Sphere)
algList.append(Fox)
algList.append(IsolatedMuonTagging)
algList.append(IsolatedElectronTagging)
algList.append(IsolatedTauTagging)
algList.append(IsolatedPhotonTagging)
algList.append(JC2FT)
algList.append(EF2)
algList.append(JC3FT)
algList.append(EF3)
algList.append(JC4FT)
algList.append(EF4)
algList.append(JC5FT)
algList.append(EF5)
algList.append(JC6FT)
algList.append(EF6)
algList.append(ComputeCorrectAngulardEdX)
algList.append(LikelihoodPID)
algList.append(LeptonID)
algList.append(ComprehensivePID)
algList.append(WWCategorisation)
algList.append(LCIOOutputProcessor)

from Configurables import ApplicationMgr
ApplicationMgr( TopAlg = algList,
                EvtSel = 'NONE',
                EvtMax   = 10,
                ExtSvc = [evtsvc, geoservice],
                OutputLevel=ERROR
              )
