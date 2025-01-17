#pragma once
#include <cmath>

extern "C"
{
    double observ_(double& co, double& co1, double& co2, double& ph1, double& ph2);
    double proba_(double& si, double& co, double& si1, double& co1, double& si2, double& co2, double& ph1, double& ph2);
    void init_();

    extern struct {double Pi;} pi_;
    extern struct {double mw; double mz; double z; double be; double ga; double xi;} energy_;
    extern struct {int order; int formf; int CPodd; int absorptive; int photon; int Zbasis;} mode_;
    // extern struct __attribute__((packed)) {int order; int formf; bool CPodd; bool absorptive; bool photon; bool Zbasis;} mode_;
    // extern struct __attribute__((packed)) {int order; int formf; unsigned int flags;} mode_;
}

namespace OOTools {
    // for debug checks
    double mw() {return energy_.mw;}
    double mz() {return energy_.mz;}

    double proba(double co, double co1, double co2, double ph1, double ph2, int formf, bool photon, int order)
    {
        mode_.formf = formf;
        mode_.photon = photon;
        mode_.order = order;
        double si = std::sqrt(1. - co*co);
        double si1 = std::sqrt(1. - co1*co1);
        double si2 = std::sqrt(1. - co2*co2);

        return proba_(si, co, si1, co1, si2, co2, ph1, ph2);
    }

    double proba_sm(double co, double co1, double co2, double ph1, double ph2)
    {
        return proba(co, co1, co2, ph1, ph2, 0, true, 0);
    }

    double sw2()
    {
        return energy_.z;
    }

    double cw2()
    {
        return 1. - energy_.z;
    }

    // TODO: check? the couplings are linear combinations of each other so this should work??
    // FIXME: this just cannot be correct!!!!
    // technically the results of proba still need to be multiplied by the input value
    // for the coupling f_i calculated backwards from an input value for the coupling a...
    // etc...
    // so I need the opposite of (19) in convert.pdf which should be recoverable from (18)
    // but I am a bit confused with the (delta, x, y)_gamma
    double proba_abph(double co, double co1, double co2, double ph1, double ph2)
    {
        double fz3 = proba(co, co1, co2, ph1, ph2, 3, false, 1);
        double fg3 = proba(co, co1, co2, ph1, ph2, 3, true, 1);

        return -sw2() / cw2() * fz3 + fg3;
    }

    double proba_awph(double co, double co1, double co2, double ph1, double ph2)
    {
        double fz1 = proba(co, co1, co2, ph1, ph2, 1, false, 1);
        double fz3 = proba(co, co1, co2, ph1, ph2, 3, false, 1);
        double fg3 = proba(co, co1, co2, ph1, ph2, 3, true, 1);

        return 1 / cw2() * fz1 + (2 - sw2()) / cw2() * fz3 + fg3;
    }

    double proba_aw(double co, double co1, double co2, double ph1, double ph2)
    {
        double fz1 = proba(co, co1, co2, ph1, ph2, 1, false, 1);
        double fz2 = proba(co, co1, co2, ph1, ph2, 2, false, 1);
        double fz3 = proba(co, co1, co2, ph1, ph2, 3, false, 1);
        double fg1 = proba(co, co1, co2, ph1, ph2, 1, true, 1);
        double fg2 = proba(co, co1, co2, ph1, ph2, 2, true, 1);
        double fg3 = proba(co, co1, co2, ph1, ph2, 3, true, 1);

        return 2. * energy_.ga * energy_.ga * (fz1 + fg1) + fz2 + fg2 + fz3 + fg3;
    }

    double proba_dg1z(double co, double co1, double co2, double ph1, double ph2)
    {
        return proba_awph(co, co1, co2, ph1, ph2) / cw2();
    }

    double proba_dkg(double co, double co1, double co2, double ph1, double ph2)
    {
        return proba_awph(co, co1, co2, ph1, ph2) + proba_abph(co, co1, co2, ph1, ph2);
    }

    double proba_lz(double co, double co1, double co2, double ph1, double ph2)
    {
        return proba_aw(co, co1, co2, ph1, ph2);
    }

    // FIXME: only for the case where W- is leptonic
    // and for tagged quark charge
    double observ(double co, double co1, double co2, double ph1, double ph2, int formf, bool photon)
    {
        mode_.formf = formf;
        // FIXME
        mode_.photon = photon;
        // TODO: test effect and figure out if false here is also false in fortran...
        mode_.absorptive = false;
        double si = std::sqrt(1. - co*co);
        double si1 = std::sqrt(1. - co1*co1);
        double si2 = std::sqrt(1. - co2*co2);

        mode_.order = 0;
        // double res = proba_(si, co, si1, co1, si2, co2, ph1, ph2);
        // necessary because the values are passed by reference to fortran...
        double nco2 = -co2;
        double ph2_pi = ph2 + pi_.Pi;
        double res = proba_(si, co, si1, co1, si2, co2, ph1, ph2) + proba_(si, co, si1, co1, si2, nco2, ph1, ph2_pi);

        mode_.order = 1;
        // res = proba_(si, co, si1, co1, si2, co2, ph1, ph2) / res;
        res = (proba_(si, co, si1, co1, si2, co2, ph1, ph2) + proba_(si, co, si1, co1, si2, nco2, ph1, ph2_pi)) / res;
        return res;
    }

    double fg1(double co, double co1, double co2, double ph1, double ph2) { return observ(co, co1, co2, ph1, ph2, 1, true); }
    double fg2(double co, double co1, double co2, double ph1, double ph2) { return observ(co, co1, co2, ph1, ph2, 2, true); }
    double fg3(double co, double co1, double co2, double ph1, double ph2) { return observ(co, co1, co2, ph1, ph2, 3, true); }
    double fg4(double co, double co1, double co2, double ph1, double ph2) { return observ(co, co1, co2, ph1, ph2, 4, true); }
    double fg5(double co, double co1, double co2, double ph1, double ph2) { return observ(co, co1, co2, ph1, ph2, 5, true); }
    double fg6(double co, double co1, double co2, double ph1, double ph2) { return observ(co, co1, co2, ph1, ph2, 6, true); }
    double fg7(double co, double co1, double co2, double ph1, double ph2) { return observ(co, co1, co2, ph1, ph2, 7, true); }

    double fz1(double co, double co1, double co2, double ph1, double ph2) { return observ(co, co1, co2, ph1, ph2, 1, false); }
    double fz2(double co, double co1, double co2, double ph1, double ph2) { return observ(co, co1, co2, ph1, ph2, 2, false); }
    double fz3(double co, double co1, double co2, double ph1, double ph2) { return observ(co, co1, co2, ph1, ph2, 3, false); }
    double fz4(double co, double co1, double co2, double ph1, double ph2) { return observ(co, co1, co2, ph1, ph2, 4, false); }
    double fz5(double co, double co1, double co2, double ph1, double ph2) { return observ(co, co1, co2, ph1, ph2, 5, false); }
    double fz6(double co, double co1, double co2, double ph1, double ph2) { return observ(co, co1, co2, ph1, ph2, 6, false); }
    double fz7(double co, double co1, double co2, double ph1, double ph2) { return observ(co, co1, co2, ph1, ph2, 7, false); }

    // conversions that I calculated myself from the Hagiwara 87 paper
    double dg1z_new(double co, double co1, double co2, double ph1, double ph2)
    {
        // FIXME: check if -1 is needed...
        // return fz1(co, co1, co2, ph1, ph2) - 2. * energy_.ga * energy_.ga * fz2(co, co1, co2, ph1, ph2) - 1.;
        return fz1(co, co1, co2, ph1, ph2) - 2. * energy_.ga * energy_.ga * fz2(co, co1, co2, ph1, ph2);
    }

    double dkg_new(double co, double co1, double co2, double ph1, double ph2)
    {
        // FIXME: check if -1 is needed...
        // return -fg1(co, co1, co2, ph1, ph2) * (2. * energy_.ga * energy_.ga - 1.) * fg2(co, co1, co2, ph1, ph2) + fg3(co, co1, co2, ph1, ph2) - 1.;
        return -fg1(co, co1, co2, ph1, ph2) * (2. * energy_.ga * energy_.ga - 1.) * fg2(co, co1, co2, ph1, ph2) + fg3(co, co1, co2, ph1, ph2);
    }

    double lz_new(double co, double co1, double co2, double ph1, double ph2)
    {
        return fz2(co, co1, co2, ph1, ph2);
    }

    // see convert.tex/.pdf
    double abph(double co, double co1, double co2, double ph1, double ph2)
    {
        return -sw2() / cw2() * fz3(co, co1, co2, ph1, ph2) + fg3(co, co1, co2, ph1, ph2);
    }

    double awph(double co, double co1, double co2, double ph1, double ph2)
    {
        return 1. / cw2() * fz1(co, co1, co2, ph1, ph2) + (2.-sw2())/cw2() * fz3(co, co1, co2, ph1, ph2) + fg3(co, co1, co2, ph1, ph2);
    }

    double aw(double co, double co1, double co2, double ph1, double ph2)
    {
        return 2. * energy_.ga * energy_.ga * (fz1(co, co1, co2, ph1, ph2) + fg1(co, co1, co2, ph1, ph2))
        + fz2(co, co1, co2, ph1, ph2) + fg2(co, co1, co2, ph1, ph2)
        + fz3(co, co1, co2, ph1, ph2) + fg3(co, co1, co2, ph1, ph2);
    }

    // TODO: link where this is from...
    // lol all of these were converted wrongly by me
    double dg1z(double co, double co1, double co2, double ph1, double ph2)
    {
        // return awph(co, co1, co2, ph1, ph2) / cw2();
        return cw2() * awph(co, co1, co2, ph1, ph2) - cw2() * abph(co, co1, co2, ph1, ph2);
    }

    double dkg(double co, double co1, double co2, double ph1, double ph2)
    {
        // return awph(co, co1, co2, ph1, ph2) + abph(co, co1, co2, ph1, ph2);
        return abph(co, co1, co2, ph1, ph2);
    }

    double lz(double co, double co1, double co2, double ph1, double ph2)
    {
        return aw(co, co1, co2, ph1, ph2);
    }

    double test()
    {
        init_();
        mode_.CPodd = true;
        // mode_.CPodd = false; //debug
        mode_.absorptive = false;
        // mode_.photon = false; //debug
        mode_.photon = true;
        mode_.Zbasis = true;
        // mode_.Zbasis = false; //debug
        mode_.formf = 4;
        // mode_.flags = -1;
        double co = 5.54e-1;
        double co1 = -2.35e-1;
        double co2 = 9.97e-1;
        double ph1 = 3.04;
        double ph2 = -6.95;

        return observ_(co, co1, co2, ph1, ph2);
    }

    // bool get_photon(){ return mode_.photon;}

    void init(double sqrt_s, double mw = 80.22, double mz = 91.17)
    {
        energy_.mw = mw;
        energy_.mz = mz;

        double s = sqrt_s * sqrt_s;
        energy_.z = 1. - (mw/mz) * (mw/mz); // squared sine of Weinberg angle

        // relativistic factors
        energy_.ga = sqrt_s / (2. * mw);
        energy_.be = std::sqrt(1. - 4. * (mw * mw)/s);
        energy_.xi = s / (s - mz * mz);

        mode_.absorptive = false; // don't worry about complex for now
        mode_.Zbasis = true;
        mode_.CPodd = true; // does nothing
    }
}
