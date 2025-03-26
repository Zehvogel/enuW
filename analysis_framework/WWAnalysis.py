from .Analysis import Analysis
import ROOT


class WWAnalysis(Analysis):
    def __init__(self, dataset):
        super().__init__(dataset)


    def define_reco_objects(self, x_angle: float):
        # select isolated lepton and the two jets and build Ws and neutrino
        self.Define("iso_lep_idx", "IsolatedElectrons_objIdx.index[0]")
        self.Define("iso_lep_charge", "PandoraPFOs.charge[iso_lep_idx]")
        # TODO: figure out: are isolated electrons also always Pandora electrons? :( (taking a look into the IsolatedLeptonTagging shows that this is not required!)
        self.Define("iso_lep_lvec", "ROOT::Math::PxPyPzEVector(PandoraPFOs.momentum.x[iso_lep_idx], PandoraPFOs.momentum.y[iso_lep_idx], PandoraPFOs.momentum.z[iso_lep_idx], PandoraPFOs.energy[iso_lep_idx])")
        self.Define("jet1_lvec", "R2jet_lvecs[0]")
        self.Define("jet2_lvec", "R2jet_lvecs[1]")
        self.Define("sqrt_s_E", "params_Energy")
        # XXX: ignoring the electron mass here
        self.Define("sqrt_s_lvec", f"ROOT::Math::PxPyPzEVector(sqrt_s_E * sin({x_angle}/2.), 0, 0, sqrt_s_E)")
        self.Define("hadronic_W_lvec", "jet1_lvec + jet2_lvec")
        # TODO: overlay contamination etc are still in this...
        self.Define("leptonic_W_lvec", "sqrt_s_lvec - hadronic_W_lvec")
        self.Define("nu_lvec", "leptonic_W_lvec - iso_lep_lvec")


    def remove_x_angle(self, x_angle: float):
        # now we would want to remove the crossing angle
        #small experiment using an actual dd4hep crossing angle boosted e+e- pair:
        # e_lvec = ROOT.Math.PxPyPzMVector(+8.750143e-01, 0., +1.250000e+02, +5.109968e-04)
        # p_lvec = ROOT.Math.PxPyPzMVector(+8.750143e-01, 0., -1.250000e+02, +5.109968e-04)
        # s_lvec = e_lvec + p_lvec
        # print(s_lvec)
        # beta = s_lvec.BoostToCM()
        # print(beta.x(), beta.y(), beta.z())
        # boost = ROOT.Math.BoostX(beta.x())
        # print(boost(s_lvec))
        # print(boost(e_lvec))
        # print(ROOT.Math.sin(0.007))
        # boost2 = ROOT.Math.BoostX(-ROOT.Math.sin(0.007))
        # print(boost2(s_lvec))
        # print(boost2(e_lvec))
        ROOT.gInterpreter.Declare(f"ROOT::Math::BoostX unboost_xangle(-std::sin({x_angle/2}));")
        lvec_list = ["iso_lep_lvec", "jet1_lvec", "jet2_lvec", "hadronic_W_lvec", "leptonic_W_lvec", "nu_lvec"]
        for lvec in lvec_list:
            self.Define(f"ub_{lvec}", f"unboost_xangle({lvec})")


    def define_hagiwara_angles(self):
        # now we want to calculate the angles as in the Hagiwara parametrisation
        # XXX: there is this one decision to make, do we use both measured daughter values or just one and flip it? :(
        # TODO: is this well defined for the e+nuqq case?
        # I need to be a bit careful, there can be some replacements that are only valid in zero-width CC03
        # It all depends on the systems of the measurement actually being the rest-frame of the parent i.e. that the decay are back-to-back
        # Exactly in this condition ang(W-, e-) == ang(W+, e+) == Theta_W
        # In the W- rest-frame f has theta and phi and fbar has -theta and phi, what are these angles actually
        # theta ang(f*, z') with f* in a system rotated such that z' direction of W- in the CMS and y' = P_e- x P_W-
        # boosted into the rest-frame of the parent of f
        # for now we just act as if we are in the double pole case as in Hagiwara and proceed.
        # TODO: guaranteed to not get 0 charge?
        self.Define("Wm_lvec", "iso_lep_charge < 0. ? ub_leptonic_W_lvec : ub_hadronic_W_lvec")
        self.Define("Wp_lvec", "iso_lep_charge > 0. ? ub_leptonic_W_lvec : ub_hadronic_W_lvec")
        # they literally did not put a cheaper cosTheta accessor into Genvectors...
        self.Define("Wm_cosTheta", "cos(Wm_lvec.Theta())")

        # need to distinguish: f == fermion child of W-, fbar == anti-fermion child of W+
        # i.e. Wm == W_lep -> f = iso_lep, fbar = not determinable, need to fold quarks
        # i.e. Wm == W_had -> fbar = nu, f = not determinable, need to fold quarks
        # so I could make the quark stuff into vectors so that the histograms work but I would also have to divide by two
        # or I make both histograms, add them and then divide by two

        ROOT.gInterpreter.Declare("#include <WWTools.h>")
        # let's act like we are in W_lep == Wm case to proceed
        self.Define("iso_lep_star_lvec", "WWTools::starVectorHagiwara(ub_leptonic_W_lvec, ub_iso_lep_lvec, {0, 0, 1}, Wm_lvec)")
        self.Define("nu_star_lvec", "WWTools::starVectorHagiwara(ub_leptonic_W_lvec, ub_nu_lvec, {0, 0, 1}, Wm_lvec)")

        self.Define("jet1_star_lvec", "WWTools::starVectorHagiwara(ub_hadronic_W_lvec, ub_jet1_lvec, {0, 0, 1}, Wm_lvec)")
        self.Define("jet2_star_lvec", "WWTools::starVectorHagiwara(ub_hadronic_W_lvec, ub_jet2_lvec, {0, 0, 1}, Wm_lvec)")

        for lvec in ["iso_lep", "nu", "jet1", "jet2"]:
            self.Define(f"{lvec}_co", f"cos({lvec}_star_lvec.Theta())")
            self.Define(f"{lvec}_ph", f"{lvec}_star_lvec.Phi()")


    def define_OO(self):
        self.Define("co", "Wm_cosTheta")
        ROOT.gInterpreter.Declare("""
        double proba_fold_Wm(double co, double co11, double co12, double co2, double ph11, double ph12, double ph2, int cpl, int order)
        {
                                  double res = 0.;
                                  res += OOTools::proba(co, co11, co2, ph11, ph2, cpl, false, order);
                                  res += OOTools::proba(co, co12, co2, ph12, ph2, cpl, false, order);
                                  return res;
        }
        double proba_fold_Wp(double co, double co1, double co21, double co22, double ph1, double ph21, double ph22, int cpl, int order)
        {
                                  double res = 0.;
                                  res += OOTools::proba(co, co1, co21, ph1, ph21, cpl, false, order);
                                  res += OOTools::proba(co, co1, co22, ph1, ph22, cpl, false, order);
                                  return res;
        }
        """)

        self.Define("S_0",   f"iso_lep_charge < 0. ? proba_fold_Wp(co, iso_lep_co, jet1_co, jet2_co, iso_lep_ph, jet1_ph, jet2_ph, 0, 0) : proba_fold_Wm(co, jet1_co, jet2_co, nu_co, jet1_ph, jet2_ph, nu_ph, 0, 0)")
        self.Define("S_1_1",   f"iso_lep_charge < 0. ? proba_fold_Wp(co, iso_lep_co, jet1_co, jet2_co, iso_lep_ph, jet1_ph, jet2_ph, 1, 1) : proba_fold_Wm(co, jet1_co, jet2_co, nu_co, jet1_ph, jet2_ph, nu_ph, 1, 1)")
        self.Define("S_1_2",   f"iso_lep_charge < 0. ? proba_fold_Wp(co, iso_lep_co, jet1_co, jet2_co, iso_lep_ph, jet1_ph, jet2_ph, 2, 1) : proba_fold_Wm(co, jet1_co, jet2_co, nu_co, jet1_ph, jet2_ph, nu_ph, 2, 1)")
        self.Define("S_1_3",   f"iso_lep_charge < 0. ? proba_fold_Wp(co, iso_lep_co, jet1_co, jet2_co, iso_lep_ph, jet1_ph, jet2_ph, 3, 1) : proba_fold_Wm(co, jet1_co, jet2_co, nu_co, jet1_ph, jet2_ph, nu_ph, 3, 1)")

        self.Define("O_1", "S_1_1 / S_0")
        self.Define("O_2", "S_1_2 / S_0")
        self.Define("O_3", "S_1_3 / S_0")

        # for debug
        self.Define("co1", "iso_lep_charge < 0. ? iso_lep_co : -iso_lep_co")
        self.Define("ph1", "iso_lep_charge < 0. ? iso_lep_ph : iso_lep_ph + ROOT::Math::Pi() <= ROOT::Math::Pi() ? iso_lep_ph + ROOT::Math::Pi() : iso_lep_ph - ROOT::Math::Pi()")
