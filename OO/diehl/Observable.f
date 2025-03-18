      PROGRAM test
C     ----- more comments below

      IMPLICIT NONE

      INTEGER  order, formf
      LOGICAL  CPodd, absorptive, photon, Zbasis
      COMMON/MODE/    order, formf, CPodd,absorptive, photon,Zbasis

      DOUBLE PRECISION  observ
      EXTERNAL          observ

      DOUBLE PRECISION  co,co1,co2,ph1,ph2
      DOUBLE PRECISION  res

C     ----- initialise masses, energy and related quantities
      CALL init

C     ----- set control parameters
      CPodd = .true.
      absorptive = .false.
      photon = .true.
      Zbasis = .true.
      formf = 4

C     ----- example set of phase space variables
      co = 5.54D-1                  ! calling program should ensure that
      co1 = -2.35D-1                ! absolute values of cosines do not
      co2 = 9.97D-1                 ! exceed 1.
      ph1 = 3.04D0
      ph2 = -6.95D0

      res = observ(co,co1,co2,ph1,ph2)
      write(*,*)   res

      END
C       of PROGRM test

C     ******************************************************************

      SUBROUTINE init
C     ---- initialise masses, c.m. energy	and Pi

      IMPLICIT NONE

      DOUBLE PRECISION  Pi, mw,mz, z, ga,be,xi, s
      COMMON/PI/        Pi
      COMMON/ENERGY/    mw,mz, z, be,ga,xi

      DATA  Pi/3.1415926535898D0/
      DATA  mw/80.22D0/, mz/91.17D0/      ! in GeV
C                              from Particle Data PR D45 (1992) No 11II

      s = 190D0**2                  ! CM energy squared in GeV**2

      z = 1.- (mw/mz)**2            ! squared sine of Weinberg angle

      ga = sqrt(s)/(2.*mw)          ! relativistic factors
      be = sqrt(1.- 4.*(mw**2)/s)
      xi = s/(s - mz**2)

      RETURN
      END
C       of SUBROUTINE Init

C     ******************************************************************

C        1         2         3         4         5         6         7
C    6789012345678901234567890123456789012345678901234567890123456789012

C     ------------------------------------------------------------------
C     -----	subroutines containing the optimal observables as
C           defined in preprint HD-THEP-93-37
C              "Optimal observables for the measurement of
C               three gauge boson couplings in e+e- -> W+W-"
C           published in Z. Phys. C62 (1994) 397
C
C     ---------   Markus Diehl, Inst. f. theoret. Physik
C                 Philosophenweg 16, D-69120 Heidelberg
C     ---------   September 1993
C     ------------------------------------------------------------------


      DOUBLE PRECISION FUNCTION observ(co,co1,co2,ph1,ph2)
C     --- optimised observable, with summation over jet charge ambiguity
C         (for hadronic decay of the W+)
C
C     ----- control parameters in COMMON/MODE:
C              order:   selects specified part in function proba
C                       0 = standard model part of squared amplitude
C                       1 = linear part multiplying specified
C                           form factor
C                    needs no initialisation, but is set within
C                    this function
C                    the following parameters have to be initialised:
C              formf    1 ... 7      index of form factor
C		   absorptive
C                       true  = imaginary part of form factor
C                       false = real part
C             Zbasis    true  = basis of gamma/Z form factors
C                       false = left/right form factors fL and fR
C                               as defined in preprint
C             photon    true  = photon form factor for Zbasis=true,
C                               left handed factor for Zbasis=false
C                       false = Z form factor for Zbasis=true,
C                               right handed factor for for Zbasis=false
C     ----- parameters in COMMON/ENERGY and COMMON/PI need
C           to be initialised, done here by subroutine init

      IMPLICIT NONE

      DOUBLE PRECISION  si,co,si1,co1,si2,co2,ph1,ph2
C     --- phase space variables as defined in Hagiwara et al.
C         Nucl. Phys. B282 (1987) on p. 270
C           si,co       sine and cosine of angle (large) theta
C                       between e- and W-
C          si1,co1      sin and cos of polar angle (small) theta and
C              ph1         azimutal angle phi for negative decay lepton
C                          or d type quark
C          si2,co2,ph2        same for W+ decay,
C                       i.e. for angles theta bar, phi bar
C     --- the sines are calculated from the cosines in this function

      DOUBLE PRECISION  proba
      EXTERNAL          proba

      DOUBLE PRECISION  Pi
      COMMON/PI/        Pi

      INTEGER  order, formf
      LOGICAL  CPodd, absorptive, photon, Zbasis
      COMMON/MODE/    order, formf, CPodd,absorptive, photon,Zbasis


C     ----- calculate sines
      si = sqrt(1.- co**2)
      si1 = sqrt(1.- co1**2)
      si2 = sqrt(1.- co2**2)


C     ----- calculate observable
      order = 0                                 ! for denominator
      observ =    proba(si,co,si1,co1,si2,co2,ph1,ph2)
     |          + proba(si,co,si1,co1,si2,-co2,ph1,ph2+Pi)
C     ----- summation over jet ambiguity,
C     ----- for hadronic decay of the W- replace second term with
C                 proba(si,co,si1,-co1,si2,co2,ph1+Pi,ph2)
C     ----- with charge tagging drop second term

      IF (observ .EQ. 0.D+00) THEN
            WRITE(*,*)  'warning from function observ: ',
     |                  'denominator was zero.'
            observ = 1.D-10                     ! has never occured
      END IF

      order = 1                                 ! for numerator
      observ =    (proba(si,co,si1,co1,si2,co2,ph1,ph2)
     |           + proba(si,co,si1,co1,si2,-co2,ph1,ph2+Pi)
     |            )
     |           / observ
C     ----- symmetrisation as above

C     ----- debug output
C     WRITE(*,*)  'CPodd: ', CPodd
C     WRITE(*,*)  'absorptive: ', absorptive
C     WRITE(*,*)  'photon: ', photon
C     WRITE(*,*)  'Zbasis: ', Zbasis

      RETURN
      END
C      of FUNCTION observ

C     ==================================================================

      DOUBLE PRECISION FUNCTION proba(si,co,si1,co1,si2,co2,ph1,ph2)
C     ----- combines amplitudes to squared amplitudes
C           normalisation: the differential cross section for
C           e+e- -> W+W-is given by
C       d(sigma) =  (alpha/z)**2 *(be/s)* (3/64)**2 *(1/ (2*Pi))* proba*
C            (appropriate form factor if not zeroth order)*
C            d(phase space)
C     with
C       d(phase space) = d(co) * d(co1)*d(phi1) * d(co2)*d(phi2) ,
C     the fine structure constant alpha and z,be,s defined as in
C     suboutine init.

      IMPLICIT NONE

      DOUBLE PRECISION  si,co,si1,co1,si2,co2,ph1,ph2
C     ---- sines are also passed as arguments here

      DOUBLE PRECISION
     |      ReStL,ImStL,ReStR,ImStR,
     |      ReAnL,ImAnL,ReAnR,ImAnR,
     |      ReAnL1,ImAnL1,ReAnR1,ImAnR1, ReAnL2,ImAnL2,ReAnR2,ImAnR2,
     |      ReAnL3,ImAnL3,ReAnR3,ImAnR3, ReAnL4,ImAnL4,ReAnR4,ImAnR4,
     |      ReAnL5,ImAnL5,ReAnR5,ImAnR5, ReAnL6,ImAnL6,ReAnR6,ImAnR6,
     |      ReAnL7,ImAnL7,ReAnR7,ImAnR7
      DOUBLE PRECISION  aStL00,aStL01,aStL10,aStL11

      INTEGER  order, formf
      LOGICAL  CPodd, absorptive, photon, Zbasis
      COMMON/MODE/    order, formf, CPodd,absorptive, photon,Zbasis

      DOUBLE PRECISION  Pi, mw,mz, z, ga,be,xi
      COMMON/PI/        Pi
      COMMON/ENERGY/    mw,mz, z, be,ga,xi

C     ------------------------------------------------------------------

C     ----- define amplitude functions

C     ----- The product of production and decay amplitudes for W+W- is
C           (mw**2) * (e**4) / (8* z**2) times the standard model terms
C           plus the same factor times the anomalous terms times the
C           left or right handed form factors as defined in the preprint.
C           e is the positron charge. (For nonleptonic W decay, this has
C           to be multiplied with the appropriate KMM element.)

C     --- Re(al part of the) St(andard model term for)
C         L(eft handed electrons):

      ReStL(si,co,si1,co1,si2,co2,ph1,ph2) =
     |  (aStL00*si*si1*si2 + (-aStL10*(1. - co)*(1. - co1) +
     |       aStL01*(1. + co)*(1. + co1))*si2*Cos(ph1) +
     |    2.*aStL11*(1. + co1*co2)*si*Cos(ph1 - ph2) +
     |    (-aStL01*(1. + co)*(1. - co2) +
     |       aStL10*(1. - co)*(1. + co2))*si1*Cos(ph2) +
     |    4.*(-co1 + co2 + co*(-1. + co1*co2))*si*
     |     Cos(ph1 + ph2))/(1. + be**2 - 2.*be*co)

      ImStL(si,co,si1,co1,si2,co2,ph1,ph2) =
     |  ((-aStL10*(1. - co)*(1. - co1) - aStL01*(1. + co)*
     |  (1. + co1))*si2*
     |     Sin(ph1) - 2.*aStL11*(co1 + co2)*si*Sin(ph1 - ph2) +
     |    (aStL01*(1. + co)*(1. - co2) + aStL10*(1. - co)*(1. + co2))*
     |     si1*Sin(ph2) + 4.*
     |     (1. + co*(co1 - co2) - co1*co2)*si*Sin(ph1 + ph2))/
     |  (1. + be**2 - 2.*be*co)

      ReStR(si,co,si1,co1,si2,co2,ph1,ph2) =
     |  4.*be*(1. + 2.*ga**2)*si*si1*si2*(-1. + xi)*z -
     |  8.*be*(-1. + co*co1)*ga*si2*(-1. + xi)*z*Cos(ph1) -
     |  4.*be*(1. + co1*co2)*si*(-1. + xi)*z*Cos(ph1 - ph2) -
     |  8.*be*(1. + co*co2)*ga*si1*(-1. + xi)*z*Cos(ph2)

      ImStR(si,co,si1,co1,si2,co2,ph1,ph2) =
     |  -8.*be*(-co + co1)*ga*si2*(-1. + xi)*z*Sin(ph1) +
     |  4.*be*(co1 + co2)*si*(-1. + xi)*z*Sin(ph1 - ph2) -
     |  8.*be*(co + co2)*ga*si1*(-1. + xi)*z*Sin(ph2)


C     --- Re(al part of the) An(omalous term multiplying the)
C         L(eft handed form factor) fL1 in the expression of the
C         amplitude (note that fL1 may be complex itself):

      ReAnL1(si,co,si1,co1,si2,co2,ph1,ph2) =
     |  be*(-1. + 2.*ga**2)*si*si1*si2 +
     |  be*(1. + co1*co2)*si*Cos(ph1 - ph2)

      ImAnL1(si,co,si1,co1,si2,co2,ph1,ph2) =
     |  -be*(co1 + co2)*si*Sin(ph1 - ph2)

      ReAnR1(si,co,si1,co1,si2,co2,ph1,ph2) =
     |  be*(-1. + 2.*ga**2)*si*si1*si2 +
     |  be*(1. + co1*co2)*si*Cos(ph1 - ph2)

      ImAnR1(si,co,si1,co1,si2,co2,ph1,ph2) =
     |  -be*(co1 + co2)*si*Sin(ph1 - ph2)

      ReAnL2(si,co,si1,co1,si2,co2,ph1,ph2) =
     |  -4.*be**3*ga**4*si*si1*si2

      ImAnL2(si,co,si1,co1,si2,co2,ph1,ph2) =
     |  0

      ReAnR2(si,co,si1,co1,si2,co2,ph1,ph2) =
     |  -4.*be**3*ga**4*si*si1*si2

      ImAnR2(si,co,si1,co1,si2,co2,ph1,ph2) =
     |  0

      ReAnL3(si,co,si1,co1,si2,co2,ph1,ph2) =
     |  -2.*be*ga**2*si*si1*si2 + be*(1. + co*co1)*ga*si2*Cos(ph1) -
     |  be*(1. - co*co2)*ga*si1*Cos(ph2)

      ImAnL3(si,co,si1,co1,si2,co2,ph1,ph2) =
     |  -be*(co + co1)*ga*si2*Sin(ph1) - be*(-co + co2)*ga*si1*Sin(ph2)

      ReAnR3(si,co,si1,co1,si2,co2,ph1,ph2) =
     |  -2.*be*ga**2*si*si1*si2 - be*(1. - co*co1)*ga*si2*Cos(ph1) +
     |  be*(1. + co*co2)*ga*si1*Cos(ph2)

      ImAnR3(si,co,si1,co1,si2,co2,ph1,ph2) =
     |  -be*(co - co1)*ga*si2*Sin(ph1) + be*(co + co2)*ga*si1*Sin(ph2)

      ReAnL4(si,co,si1,co1,si2,co2,ph1,ph2) =
     |  -be*(co + co1)*ga*si2*Sin(ph1) + be*(-co + co2)*ga*si1*Sin(ph2)

      ImAnL4(si,co,si1,co1,si2,co2,ph1,ph2) =
     |  -be*(1. + co*co1)*ga*si2*Cos(ph1) +
     |  be*(-1. + co*co2)*ga*si1*Cos(ph2)

      ReAnR4(si,co,si1,co1,si2,co2,ph1,ph2) =
     |  be*(-co + co1)*ga*si2*Sin(ph1) - be*(co + co2)*ga*si1*Sin(ph2)

      ImAnR4(si,co,si1,co1,si2,co2,ph1,ph2) =
     |  be*(1. - co*co1)*ga*si2*Cos(ph1) +
     |  be*(1. + co*co2)*ga*si1*Cos(ph2)

      ReAnL5(si,co,si1,co1,si2,co2,ph1,ph2) =
     |  -be**2*(co + co1)*ga*si2*Cos(ph1) -
     |  be**2*(-co + co2)*ga*si1*Cos(ph2)

      ImAnL5(si,co,si1,co1,si2,co2,ph1,ph2) =
     |  be**2*(1. + co*co1)*ga*si2*Sin(ph1) -
     |  be**2*(1. - co*co2)*ga*si1*Sin(ph2)

      ReAnR5(si,co,si1,co1,si2,co2,ph1,ph2) =
     |  -be**2*(co - co1)*ga*si2*Cos(ph1) +
     |  be**2*(co + co2)*ga*si1*Cos(ph2)

      ImAnR5(si,co,si1,co1,si2,co2,ph1,ph2) =
     |  -be**2*(1. - co*co1)*ga*si2*Sin(ph1) +
     |  be**2*(1. + co*co2)*ga*si1*Sin(ph2)

      ReAnL6(si,co,si1,co1,si2,co2,ph1,ph2) =
     |  -(1. + co*co1)*ga*si2*Sin(ph1) -
     |  (1. + co1*co2)*si*Sin(ph1 - ph2) -
     |  (1. - co*co2)*ga*si1*Sin(ph2)

      ImAnL6(si,co,si1,co1,si2,co2,ph1,ph2) =
     |  -(co + co1)*ga*si2*Cos(ph1) - (co1 + co2)*si*Cos(ph1 - ph2) -
     |  (co - co2)*ga*si1*Cos(ph2)

      ReAnR6(si,co,si1,co1,si2,co2,ph1,ph2) =
     |  -(-1. + co*co1)*ga*si2*Sin(ph1) -
     |  (1. + co1*co2)*si*Sin(ph1 - ph2) +
     |  (1. + co*co2)*ga*si1*Sin(ph2)

      ImAnR6(si,co,si1,co1,si2,co2,ph1,ph2) =
     |  -(co - co1)*ga*si2*Cos(ph1) - (co1 + co2)*si*Cos(ph1 - ph2) -
     |  (co + co2)*ga*si1*Cos(ph2)

      ReAnL7(si,co,si1,co1,si2,co2,ph1,ph2) =
     |  -4.*be**2*(1. + co1*co2)*ga**2*si*Sin(ph1 - ph2)

      ImAnL7(si,co,si1,co1,si2,co2,ph1,ph2) =
     |  -4.*be**2*(co1 + co2)*ga**2*si*Cos(ph1 - ph2)

      ReAnR7(si,co,si1,co1,si2,co2,ph1,ph2) =
     |  -4.*be**2*(1. + co1*co2)*ga**2*si*Sin(ph1 - ph2)

      ImAnR7(si,co,si1,co1,si2,co2,ph1,ph2) =
     |  -4.*be**2*(co1 + co2)*ga**2*si*Cos(ph1 - ph2)


C     ---- auxillary functions for left handed amplitudes
      aStL00 =
     |  2/be + 2*be - 2/(be*ga**2) + 2*be/ga**2 +
     |  co*(-4 - 4/ga**2 + 2*be**2*(-1 + xi)*(2 - 4*z) +
     |  4*be**2*ga**2*(-1 + xi)*(2 - 4*z)) - be*(-1 + xi)*(2 - 4*z) -
     |  be**3*(-1 + xi)*(2 - 4*z) - 2*be*ga**2*(-1 + xi)*(2 - 4*z) -
     |  2*be**3*ga**2*(-1 + xi)*(2 - 4*z)

      aStL01 =
     |  -2/(be*ga) + 2*(1 - be)/(be*ga) - 2*be/ga +
     |  co*(4/ga - 2*be**2*ga*(-1 + xi)*(2 - 4*z)) +
     |  be*ga*(-1 + xi)*(2 - 4*z) + be**3*ga*(-1 + xi)*(2 - 4*z)

      aStL10 =
     |  2/(be*ga) + 2*be/ga - 2*(1 + be)/(be*ga) +
     |  co*(-4/ga + 2*be**2*ga*(-1 + xi)*(2 - 4*z)) -
     |  be*ga*(-1 + xi)*(2 - 4*z) - be**3*ga*(-1 + xi)*(2 - 4*z)

      aStL11 =
     |  -(be/ga**2) + co*(2/ga**2 - be**2*(-1 + xi)*(2 - 4*z)) +
     |  be*(-1 + xi)*(2 - 4*z)/2 + be**3*(-1 + xi)*(2 - 4*z)/2

C     at last ...

C     ------------------------------------------------------------------

C     -------- squared T-matrix to zeroth order

      IF ( order .EQ. 0 ) THEN
            proba = ReStL(si,co,si1,co1,si2,co2,ph1,ph2)**2 +
     |              ImStL(si,co,si1,co1,si2,co2,ph1,ph2)**2 +
     |              ReStR(si,co,si1,co1,si2,co2,ph1,ph2)**2 +
     |              ImStR(si,co,si1,co1,si2,co2,ph1,ph2)**2
            GO TO 900
      END IF


C     ----- anomalous amplitudes as function of form factor

 100  IF (formf .EQ. 1) THEN
            ReAnL = ReAnL1(si,co,si1,co1,si2,co2,ph1,ph2)
      ELSE IF (formf .EQ. 2) THEN
            ReAnL = ReAnL2(si,co,si1,co1,si2,co2,ph1,ph2)
      ELSE IF (formf .EQ. 3) THEN
            ReAnL = ReAnL3(si,co,si1,co1,si2,co2,ph1,ph2)
      ELSE IF (formf .EQ. 4) THEN
            ReAnL = ReAnL4(si,co,si1,co1,si2,co2,ph1,ph2)
      ELSE IF (formf .EQ. 5) THEN
            ReAnL = ReAnL5(si,co,si1,co1,si2,co2,ph1,ph2)
      ELSE IF (formf .EQ. 6) THEN
            ReAnL = ReAnL6(si,co,si1,co1,si2,co2,ph1,ph2)
      ELSE IF (formf .EQ. 7) THEN
            ReAnL = ReAnL7(si,co,si1,co1,si2,co2,ph1,ph2)
      END IF

 200  IF (formf .EQ. 1) THEN
            ImAnL = ImAnL1(si,co,si1,co1,si2,co2,ph1,ph2)
      ELSE IF (formf .EQ. 2) THEN
            ImAnL = ImAnL2(si,co,si1,co1,si2,co2,ph1,ph2)
      ELSE IF (formf .EQ. 3) THEN
            ImAnL = ImAnL3(si,co,si1,co1,si2,co2,ph1,ph2)
      ELSE IF (formf .EQ. 4) THEN
            ImAnL = ImAnL4(si,co,si1,co1,si2,co2,ph1,ph2)
      ELSE IF (formf .EQ. 5) THEN
            ImAnL = ImAnL5(si,co,si1,co1,si2,co2,ph1,ph2)
      ELSE IF (formf .EQ. 6) THEN
            ImAnL = ImAnL6(si,co,si1,co1,si2,co2,ph1,ph2)
      ELSE IF (formf .EQ. 7) THEN
            ImAnL = ImAnL7(si,co,si1,co1,si2,co2,ph1,ph2)
      END IF

 300  IF (formf .EQ. 1) THEN
            ReAnR = ReAnR1(si,co,si1,co1,si2,co2,ph1,ph2)
      ELSE IF (formf .EQ. 2) THEN
            ReAnR = ReAnR2(si,co,si1,co1,si2,co2,ph1,ph2)
      ELSE IF (formf .EQ. 3) THEN
            ReAnR = ReAnR3(si,co,si1,co1,si2,co2,ph1,ph2)
      ELSE IF (formf .EQ. 4) THEN
            ReAnR = ReAnR4(si,co,si1,co1,si2,co2,ph1,ph2)
      ELSE IF (formf .EQ. 5) THEN
            ReAnR = ReAnR5(si,co,si1,co1,si2,co2,ph1,ph2)
      ELSE IF (formf .EQ. 6) THEN
            ReAnR = ReAnR6(si,co,si1,co1,si2,co2,ph1,ph2)
      ELSE IF (formf .EQ. 7) THEN
            ReAnR = ReAnR7(si,co,si1,co1,si2,co2,ph1,ph2)
      END IF

 400  IF (formf .EQ. 1) THEN
            ImAnR = ImAnR1(si,co,si1,co1,si2,co2,ph1,ph2)
      ELSE IF (formf .EQ. 2) THEN
            ImAnR = ImAnR2(si,co,si1,co1,si2,co2,ph1,ph2)
      ELSE IF (formf .EQ. 3) THEN
            ImAnR = ImAnR3(si,co,si1,co1,si2,co2,ph1,ph2)
      ELSE IF (formf .EQ. 4) THEN
            ImAnR = ImAnR4(si,co,si1,co1,si2,co2,ph1,ph2)
      ELSE IF (formf .EQ. 5) THEN
            ImAnR = ImAnR5(si,co,si1,co1,si2,co2,ph1,ph2)
      ELSE IF (formf .EQ. 6) THEN
            ImAnR = ImAnR6(si,co,si1,co1,si2,co2,ph1,ph2)
      ELSE IF (formf .EQ. 7) THEN
            ImAnR = ImAnR7(si,co,si1,co1,si2,co2,ph1,ph2)
      END IF


      IF (Zbasis) THEN
C     -------- amplitudes proportional to boson form factors (photon or Z)
            IF (photon) THEN                          ! photon
                  ReAnL = 4.*z*ReAnL
                  ImAnL = 4.*z*ImAnL
                  ReAnR = 4.*z*ReAnR
                  ImAnR = 4.*z*ImAnR
            ELSE                                      ! Z boson
                  ReAnL = (2.-4.*z)*xi*ReAnL
                  ImAnL = (2.-4.*z)*xi*ImAnL
                  ReAnR = -4.*z*xi*ReAnR
                  ImAnR = -4.*z*xi*ImAnR
            END IF
      ELSE
C     ------- left or right handed form factors
            IF (photon) THEN              ! left handed amplitude only
                  ReAnR = 0.
                  ImAnR = 0.
            ELSE                          ! right handed only
                  ReAnL = 0.
                  ImAnL = 0.
            END IF
	END IF


C     -------- squared T-matrix: first order part

      IF (absorptive) THEN                ! imaginary form factor
            proba = 2.*(ReAnL*ImStL(si,co,si1,co1,si2,co2,ph1,ph2)
     |               -  ImAnL*ReStL(si,co,si1,co1,si2,co2,ph1,ph2)
     |               +  ReAnR*ImStR(si,co,si1,co1,si2,co2,ph1,ph2)
     |               -  ImAnR*ReStR(si,co,si1,co1,si2,co2,ph1,ph2))

      ELSE                                ! real form factor
            proba = 2.*(ReAnL*ReStL(si,co,si1,co1,si2,co2,ph1,ph2)
     |               +  ImAnL*ImStL(si,co,si1,co1,si2,co2,ph1,ph2)
     |               +  ReAnR*ReStR(si,co,si1,co1,si2,co2,ph1,ph2)
     |               +  ImAnR*ImStR(si,co,si1,co1,si2,co2,ph1,ph2))
      END IF


 900  CONTINUE

      RETURN
      END
C       of FUNCTION proba
