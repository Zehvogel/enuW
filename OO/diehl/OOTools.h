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

    // FIXME: only for the case where W- is leptonic
    double observ(double co, double co1, double co2, double ph1, double ph2, int formf, bool photon)
    {
        mode_.formf = formf;
        // FIXME
        mode_.photon = photon;
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

    double O_fg1(double co, double co1, double co2, double ph1, double ph2) { return observ(co, co1, co2, ph1, ph2, 1, true); }
    double O_fg2(double co, double co1, double co2, double ph1, double ph2) { return observ(co, co1, co2, ph1, ph2, 2, true); }
    double O_fg3(double co, double co1, double co2, double ph1, double ph2) { return observ(co, co1, co2, ph1, ph2, 3, true); }
    double O_fg4(double co, double co1, double co2, double ph1, double ph2) { return observ(co, co1, co2, ph1, ph2, 4, true); }
    double O_fg5(double co, double co1, double co2, double ph1, double ph2) { return observ(co, co1, co2, ph1, ph2, 5, true); }
    double O_fg6(double co, double co1, double co2, double ph1, double ph2) { return observ(co, co1, co2, ph1, ph2, 6, true); }
    double O_fg7(double co, double co1, double co2, double ph1, double ph2) { return observ(co, co1, co2, ph1, ph2, 7, true); }

    double O_fz1(double co, double co1, double co2, double ph1, double ph2) { return observ(co, co1, co2, ph1, ph2, 1, false); }
    double O_fz2(double co, double co1, double co2, double ph1, double ph2) { return observ(co, co1, co2, ph1, ph2, 2, false); }
    double O_fz3(double co, double co1, double co2, double ph1, double ph2) { return observ(co, co1, co2, ph1, ph2, 3, false); }
    double O_fz4(double co, double co1, double co2, double ph1, double ph2) { return observ(co, co1, co2, ph1, ph2, 4, false); }
    double O_fz5(double co, double co1, double co2, double ph1, double ph2) { return observ(co, co1, co2, ph1, ph2, 5, false); }
    double O_fz6(double co, double co1, double co2, double ph1, double ph2) { return observ(co, co1, co2, ph1, ph2, 6, false); }
    double O_fz7(double co, double co1, double co2, double ph1, double ph2) { return observ(co, co1, co2, ph1, ph2, 7, false); }

    // see convert.tex/.pdf
    double O_abph(double co, double co1, double co2, double ph1, double ph2)
    {
        return -sw2() / cw2() * O_fz3(co, co1, co2, ph1, ph2) + O_fg3(co, co1, co2, ph1, ph2);
    }

    double O_awph(double co, double co1, double co2, double ph1, double ph2)
    {
        return 1. / cw2() * O_fz1(co, co1, co2, ph1, ph2) + (2.-sw2())/cw2() * O_fz3(co, co1, co2, ph1, ph2) + O_fg3(co, co1, co2, ph1, ph2);
    }

    double O_aw(double co, double co1, double co2, double ph1, double ph2)
    {
        return 2. * energy_.ga * energy_.ga * (O_fz1(co, co1, co2, ph1, ph2) + O_fg1(co, co1, co2, ph1, ph2))
        + O_fz2(co, co1, co2, ph1, ph2) + O_fg2(co, co1, co2, ph1, ph2)
        + O_fz3(co, co1, co2, ph1, ph2) + O_fg3(co, co1, co2, ph1, ph2);
    }

    // TODO: link where this is from...
    double O_dg1z(double co, double co1, double co2, double ph1, double ph2)
    {
        return cw2() * O_awph(co, co1, co2, ph1, ph2) - cw2() * O_abph(co, co1, co2, ph1, ph2);
    }

    double O_dkg(double co, double co1, double co2, double ph1, double ph2)
    {
        return O_abph(co, co1, co2, ph1, ph2);
    }

    double O_lz(double co, double co1, double co2, double ph1, double ph2)
    {
        return O_aw(co, co1, co2, ph1, ph2);
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
