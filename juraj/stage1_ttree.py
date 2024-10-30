#!/usr/bin/env python
# coding: utf-8

# In[1]:


import ROOT


# In[2]:


# get_ipython().run_line_magic('jsroot', 'on')


# In[3]:


ROOT.EnableImplicitMT(6)
ROOT.TH1.SetDefaultSumw2()


# In[4]:


# wasn't necessary in the file that I copied all of this from to make a reproducer :shrugs:
ROOT.gInterpreter.Declare("#include <edm4hep/ReconstructedParticleData.h>")


# In[5]:


ROOT.gInterpreter.ProcessLine("using namespace ROOT::VecOps;")


# In[6]:


plots = False
# plots = True
# all = False
all = True
df = None
if all:
    df = ROOT.RDataFrame("events", "root://eosuser.cern.ch//eos/user/l/lreichen/paris/*.edm4hep.root")
else:
    df = ROOT.RDataFrame("events", "root://eosuser.cern.ch//eos/user/l/lreichen/paris/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500105.P4f_sw_sl.eL.pL.n000.d_dstm_15064_0.miniDST.edm4hep.root")


# In[7]:


ROOT.gInterpreter.ProcessLine("""
template <typename coll_t>
ROOT::VecOps::RVec<coll_t> index2entry(RVec<coll_t> collection, RVec<int> indices)
{
   RVec<coll_t> res;
   res.reserve(indices.size());

   for (auto i : indices)
   {
      res.push_back(collection[i]);
   }

   return res;
}
""")


# In[ ]:


simple_tjs = ["e", "nu", "isr1", "isr2", "ovl"]

# ugly [0] at the right place is needed to obtain index... as the mask returns an RVec with one element
df = df.Define("e_idx", "_TrueJets_PID_TrueJetPID_particle.index[abs(TrueJets_PID_TrueJetPID.PDG) == 11 && abs(TrueJets_PID_TrueJetPID.type) == 2][0]")
df = df.Define("e", "TrueJets[_TrueJets_PID_TrueJetPID_particle.index[abs(TrueJets_PID_TrueJetPID.PDG) == 11 && abs(TrueJets_PID_TrueJetPID.type) == 2][0]]")
df = df.Define("nu_idx", "_TrueJets_PID_TrueJetPID_particle.index[abs(TrueJets_PID_TrueJetPID.PDG) == 12 && abs(TrueJets_PID_TrueJetPID.type) == 2][0]")
df = df.Define("nu", "TrueJets[_TrueJets_PID_TrueJetPID_particle.index[abs(TrueJets_PID_TrueJetPID.PDG) == 12 && abs(TrueJets_PID_TrueJetPID.type) == 2][0]]")
df = df.Define("isr", "index2entry(TrueJets, _TrueJets_PID_TrueJetPID_particle.index[abs(TrueJets_PID_TrueJetPID.PDG) == 12 && abs(TrueJets_PID_TrueJetPID.type) == 2])")
df = df.Define("isr1", "isr[0]")
df = df.Define("isr2", "isr[1]")
df = df.Define("ovl", "TrueJets[_TrueJets_PID_TrueJetPID_particle.index[abs(TrueJets_PID_TrueJetPID.PDG) == 0 && abs(TrueJets_PID_TrueJetPID.type) == 5][0]]")

# quarks:
df = df.Define("quark_icn_pid_mask", "InitialColourNeutrals_PID_TrueJet_fafpi.type == 1 || InitialColourNeutrals_PID_TrueJet_fafpi.type == 3")
df = df.Define("quark_icn_idx", "_InitialColourNeutrals_PID_TrueJet_fafpi_particle.index[quark_icn_pid_mask][0]")
df = df.Define("tj_quarks", "index2entry(TrueJets, Range(InitialColourNeutrals.particles_begin[quark_icn_idx], InitialColourNeutrals.particles_end[quark_icn_idx]))")
df = df.Define("tj_quarks_px", "index2entry(TrueJets.momentum.x, Range(InitialColourNeutrals.particles_begin[quark_icn_idx], InitialColourNeutrals.particles_end[quark_icn_idx]))")
df = df.Define("tj_quarks_py", "index2entry(TrueJets.momentum.y, Range(InitialColourNeutrals.particles_begin[quark_icn_idx], InitialColourNeutrals.particles_end[quark_icn_idx]))")
df = df.Define("tj_quarks_pz", "index2entry(TrueJets.momentum.z, Range(InitialColourNeutrals.particles_begin[quark_icn_idx], InitialColourNeutrals.particles_end[quark_icn_idx]))")
df = df.Define("tj_quarks_e", "index2entry(TrueJets.energy, Range(InitialColourNeutrals.particles_begin[quark_icn_idx], InitialColourNeutrals.particles_end[quark_icn_idx]))")
df = df.Define("had_lvecs", "Construct<ROOT::Math::PxPyPzEVector>(tj_quarks_px, tj_quarks_py, tj_quarks_pz, tj_quarks_e)")
df = df.Define("had_sum_lvec", "Sum(had_lvecs, ROOT::Math::PxPyPzEVector())")

for tj in simple_tjs:
    df = df.Define(f"{tj}_lvec", f"ROOT::Math::PxPyPzEVector({tj}.momentum.x, {tj}.momentum.y, {tj}.momentum.z, {tj}.energy)")

df = df.Alias("had_lvecs_px", "tj_quarks_px")
df = df.Alias("had_lvecs_py", "tj_quarks_py")
df = df.Alias("had_lvecs_pz", "tj_quarks_pz")
df = df.Alias("had_lvecs_e", "tj_quarks_e")

df = df.Define("e_charge", "e.charge")


# In[ ]:


df = df.Define("e_energy", "e_lvec.energy()")
df = df.Define("nu_energy", "nu_lvec.energy()")
df = df.Define("isr1_energy", "isr1_lvec.energy()")
df = df.Define("isr2_energy", "isr2_lvec.energy()")
df = df.Define("ovl_energy", "ovl_lvec.energy()")
df = df.Define("had_energy", "had_sum_lvec.energy()")


# In[ ]:


h_e = df.Histo1D(("", "; reco energy", 150, 0., 150.), "e_energy")
h_nu = df.Histo1D(("", "; reco energy", 150, 0., 150.), "nu_energy")
h_isr1 = df.Histo1D(("", "; reco energy", 150, 0., 150.), "isr1_energy")
h_isr2 = df.Histo1D(("", "; reco energy", 150, 0., 150.), "isr2_energy")
h_ovl = df.Histo1D(("", "; reco energy", 150, 0., 150.), "ovl_energy")
h_had = df.Histo1D(("", "; reco energy", 150, 0., 150.), "had_energy")


# In[ ]:


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


# In[ ]:


comps = ["px", "py", "pz", "e"]
branches = [f"{tj}_lvec" for tj in simple_tjs] + [f"had_lvecs_{c}" for c in comps] + ["e_charge"]
if all:
    df.Snapshot("events", "ttree_all.nano.root", branches)
else:
    df.Snapshot("events", "ttree_test.nano.root", branches)


# In[ ]:




