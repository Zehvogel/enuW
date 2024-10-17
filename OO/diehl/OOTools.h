#pragma once
#include <cmath>

extern "C"
{
    double observ_(double& co, double& co1, double& co2, double& ph1, double& ph2);
    void init_();

    extern struct {double Pi;} pi_;
    extern struct {double mw; double mz; double z; double be; double ga; double xi;} energy_;
    extern struct {int order; int formf; bool CPodd; bool absorptive; bool photon; bool Zbasis;} mode_;
}

namespace OOTools {
    // for debug checks
    double mw() {return energy_.mw;}
    double mz() {return energy_.mz;}

    // FIXME: only for the case where W- is leptonic
    double observ(double co, double co1, double co2, double ph1, double ph2, int formf, bool photon)
    {
        mode_.formf = formf;
        mode_.photon = photon;
        return observ_(co, co1, co2, ph1, ph2);
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

    double sw2()
    {
        return energy_.z;
    }

    double cw2()
    {
        return 1. - energy_.z;
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
    double dg1z(double co, double co1, double co2, double ph1, double ph2)
    {
        return awph(co, co1, co2, ph1, ph2) / cw2();
    }

    double dkg(double co, double co1, double co2, double ph1, double ph2)
    {
        return awph(co, co1, co2, ph1, ph2) + abph(co, co1, co2, ph1, ph2);
    }

    double lz(double co, double co1, double co2, double ph1, double ph2)
    {
        return aw(co, co1, co2, ph1, ph2);
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
    }
}
