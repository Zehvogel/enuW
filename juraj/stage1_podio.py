#!/usr/bin/env python
# coding: utf-8

# In[1]:


import ROOT
import glob


# In[2]:


# get_ipython().run_line_magic('jsroot', 'on')


# In[3]:


ROOT.EnableImplicitMT(6)
ROOT.TH1.SetDefaultSumw2()


# In[ ]:


ROOT.gInterpreter.ProcessLine("using namespace ROOT::VecOps;")


# In[5]:


ROOT.gInterpreter.Declare("#include <podio/DataSource.h>")
ROOT.gInterpreter.Declare("#include <edm4hep/utils/kinematics.h>")
ROOT.gSystem.Load("libedm4hep")
_edm  = ROOT.edm4hep.ReconstructedParticleData()


# In[6]:


plots = False
# plots = True
# all = False
all = True
df = None
if all:
    paths = glob.glob("/eos/user/l/lreichen/paris/*.edm4hep.root")
    files = [f"root://eosuser.cern.ch/{path}" for path in paths]
    df = ROOT.podio.CreateDataFrame(files)
else:
    df = ROOT.podio.CreateDataFrame("root://eosuser.cern.ch//eos/user/l/lreichen/paris/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500105.P4f_sw_sl.eL.pL.n000.d_dstm_15064_0.miniDST.edm4hep.root")


# In[7]:


ROOT.gInterpreter.ProcessLine("""
// returns the first particle whose absolute type and pdg values are fulfilling type and pdg
edm4hep::ReconstructedParticle getParticleWithPID(const edm4hep::ParticleIDCollection& pid_col, std::int32_t type, std::int32_t pdg)
{
    for (const edm4hep::ParticleID& pid : pid_col) {
        if (abs(pid.getType()) == type && abs(pid.getPDG()) == pdg) {
            return pid.getParticle();
        }
    }
    return edm4hep::ReconstructedParticle();
}
""")


# In[8]:


ROOT.gInterpreter.ProcessLine("""
std::array<edm4hep::ReconstructedParticle, 5> getTrueJets(const edm4hep::ParticleIDCollection& tj_pid_col)
{
    edm4hep::ReconstructedParticle e = getParticleWithPID(tj_pid_col, 2, 11);
    edm4hep::ReconstructedParticle nu = getParticleWithPID(tj_pid_col, 2, 12);

    std::vector<edm4hep::ReconstructedParticle> isr;
    for (const auto& pid : tj_pid_col) {
        if (abs(pid.getType()) == 4) {
            isr.push_back(pid.getParticle());
        }
    }

    edm4hep::ReconstructedParticle isr1 = isr[0];
    edm4hep::ReconstructedParticle isr2 = isr[1];

    edm4hep::ReconstructedParticle overlay = getParticleWithPID(tj_pid_col, 5, 0);

    return std::array<edm4hep::ReconstructedParticle, 5>{e, nu, isr1, isr2, overlay};
}
""")


# In[9]:


ROOT.gInterpreter.ProcessLine("""
ROOT::VecOps::RVec<edm4hep::ReconstructedParticle> getQuarkTrueJets(const edm4hep::ParticleIDCollection& icn_pid_col)
{
    RVec<edm4hep::ReconstructedParticle> res;
    for (const auto& pid : icn_pid_col) {
        if (pid.getType() == 1 || pid.getType() == 3) {
            // should only happen once anyway
            // fingers crossed
            auto particles = pid.getParticle().getParticles();
            return RVec<edm4hep::ReconstructedParticle>(particles.begin(), particles.end());
        }
    }
    return res;
}
""")


# In[10]:


ROOT.gInterpreter.ProcessLine("""
ROOT::VecOps::RVec<edm4hep::MCParticle> getEnergySafeMCParticles(const edm4hep::RecoMCParticleLinkCollection& links,const edm4hep::ReconstructedParticle& reco)
{
    RVec<edm4hep::MCParticle> res;
    for (const auto& link : links) {
        if (link.getFrom() == reco && link.getWeight() == 1) {
            res.push_back(link.getTo());
        }
    }
    return res;
}
""")


# In[11]:


df = df.Define("tj_pfos", "getTrueJets(TrueJets_PID_TrueJetPID)")
df = df.Define("tj_quarks", "getQuarkTrueJets(InitialColourNeutrals_PID_TrueJet_fafpi)")

simple_tjs = ["e", "nu", "isr1", "isr2", "ovl"]

for i, tj in enumerate(simple_tjs):
    df = df.Define(tj, f"tj_pfos[{i}]")
    df = df.Define(f"{tj}_lvec", f"edm4hep::utils::p4({tj}, edm4hep::utils::UseEnergy)")

# for the quark truejets everything needs one more loop because their number can vary
df = df.Define("had_lvecs", "return Map(tj_quarks, [] (const auto& el) {return edm4hep::utils::p4(el, edm4hep::utils::UseEnergy);});")
df = df.Define("had_sum_lvec", "Sum(had_lvecs, ROOT::Math::PxPyPzEVector())")

df = df.Define("had_lvecs_px", "return Map(had_lvecs, [] (const auto& el) {return el.Px();})")
df = df.Define("had_lvecs_py", "return Map(had_lvecs, [] (const auto& el) {return el.Py();})")
df = df.Define("had_lvecs_pz", "return Map(had_lvecs, [] (const auto& el) {return el.Pz();})")
df = df.Define("had_lvecs_e",  "return Map(had_lvecs, [] (const auto& el) {return el.E();})")

df = df.Define("e_charge", "e.getCharge()")


# In[12]:


df = df.Define("e_energy", "e_lvec.energy()")
df = df.Define("nu_energy", "nu_lvec.energy()")
df = df.Define("isr1_energy", "isr1_lvec.energy()")
df = df.Define("isr2_energy", "isr2_lvec.energy()")
df = df.Define("ovl_energy", "ovl_lvec.energy()")
df = df.Define("had_energy", "had_sum_lvec.energy()")


# In[13]:


h_e = df.Histo1D(("", "; reco energy", 150, 0., 150.), "e_energy")
h_nu = df.Histo1D(("", "; reco energy", 150, 0., 150.), "nu_energy")
h_isr1 = df.Histo1D(("", "; reco energy", 150, 0., 150.), "isr1_energy")
h_isr2 = df.Histo1D(("", "; reco energy", 150, 0., 150.), "isr2_energy")
h_ovl = df.Histo1D(("", "; reco energy", 150, 0., 150.), "ovl_energy")
h_had = df.Histo1D(("", "; reco energy", 150, 0., 150.), "had_energy")


# In[14]:


if plots:
    c_e = ROOT.TCanvas()
    h_e.Draw()
    c_e.Draw()

    c_nu = ROOT.TCanvas()
    h_nu.Draw()
    c_nu.Draw()

    c_isr1 = ROOT.TCanvas()
    h_isr1.Draw()
    c_isr1.Draw()

    c_isr2 = ROOT.TCanvas()
    h_isr2.Draw()
    c_isr2.Draw()

    c_ovl = ROOT.TCanvas()
    h_ovl.Draw()
    c_ovl.Draw()

    c_had = ROOT.TCanvas()
    h_had.Draw()
    c_had.Draw()


# In[15]:


comps = ["px", "py", "pz", "e"]
branches = [f"{tj}_lvec" for tj in simple_tjs] + [f"had_lvecs_{c}" for c in comps] + ["e_charge"]
if all:
    df.Snapshot("events", "podio_all.nano.root", branches)
else:
    df.Snapshot("events", "podio_test.nano.root", branches)


# In[ ]:




