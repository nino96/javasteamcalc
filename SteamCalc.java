/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package com.spayc;

/**
 *
 * @author Kehinde
 */
public class SteamCalc {
    
    /** Molar gas constant for water
     *  R = 0.461526 kJ/kg K
     */
    final static double R = 0.461526;

    // Critical Parameters
    /** Critical Temperature
     *  Tc = 647.096 K
     */
    final static double Tc = 647.096;
    /** Critical Pressure
     *  Pc = 22.064 MPa
     */
    final static double Pc = 22.064;
    /** Density at critical point
     *  ρc = 322 kg/m3
     */
    final static double RHOc = 322;

    // Saturation parameters
    final static double nsat[] = {1.1670521452767E3,-7.2421316703206E5,
    -1.7073846940092E1,1.2020824702470E4,-3.2325550322333E6,1.4915108613530E1,
    -4.8232657361591E3,4.0511340542057E5,-0.23855557567849,6.5017534844798E2};

    /** Fundamental equation: specific Gibbs Free Energy
     *  g(p,T)
     *
     * @param p pressure
     * @param t absolute temperature
     * @return double
     */
    public static double g(double p, double t) {
        switch (getRegion(p,t)) {
            case 1:
                return Region1.g(p,t);
            case 2:
                return Region2.g(p, t);
            case 5:
                return Region5.g(p, t);
            default:
                return Double.NaN;
        }
    }

    /** Gibbs Free Energy for Wet Steam
     *
     * @param p pressure
     * @param x dryness fraction
     * @return double
     */

    public static double gsat(double p, float x) {
        // Check validity of x and p
        if (x < 0 || x > 1.0f) return Double.NaN;
        double t = tsat(p);
        if (t == Double.NaN) return t;
        return (1-x) * Region1.g(p, t) + x * Region2.gsat(p, t);
    }

    /** Gibbs Free Energy for Wet Steam
     *
     * @param x dryness fraction
     * @param t absolute temperature
     * @return double
     */

    public static double gsat(float x, double t) {
        // Check validity of x and t
        if (x < 0 || x > 1.0f) return Double.NaN;
        double p = psat(t);
        if (p == Double.NaN) return p;
        return (1-x) * Region1.g(p, t) + x * Region2.gsat(p, t);
    }

    /** Specific volume v(p,T)
     *
     * @param p pressure
     * @param t absolute temperature
     * @return double
     */
    public static double v(double p, double t) {
        switch (getRegion(p,t)) {
            case 1:
                return Region1.v(p,t);
            case 2:
                return Region2.v(p,t);
            case 3:
                return Region3.v(p,t);
            case 5:
                return Region5.v(p, t);
            default:
                return Double.NaN;
        }
    }

    /** Specific volume for Wet Steam
     *
     * @param p pressure
     * @param x dryness fraction
     * @return double
     */

    public static double vsat(double p, float x) {
        // Check validity of x and p
        if (x < 0 || x > 1.0f) return Double.NaN;
        double t = tsat(p);
        if (t == Double.NaN) return t;
        return (1-x) * Region1.v(p, t) + x * Region2.vsat(p, t);
    }

    /** Specific volume for Wet Steam
     *
     * @param x dryness fraction
     * @param t absolute temperature
     * @return double
     */

    public static double vsat(float x, double t) {
        // Check validity of x and t
        if (x < 0 || x > 1.0f) return Double.NaN;
        double p = psat(t);
        if (p == Double.NaN) return p;
        return (1-x) * Region1.v(p, t) + x * Region2.vsat(p, t);
    }

    /** Specific enthalpy h(p,T)
     *
     * @param p pressure
     * @param t absolute temperature
     * @return double
     */
    public static double h(double p, double t) {
        switch (getRegion(p,t)) {
            case 1:
                return Region1.h(p,t);
            case 2:
                return Region2.h(p,t);
            case 3:
                return Region3.h(p,t);
            case 5:
                return Region5.h(p, t);
            default:
                return Double.NaN;
        }
    }

    /** Specific enthalpy for Wet Steam
     *
     * @param p pressure
     * @param x dryness fraction
     * @return double
     */

    public static double hsat(double p, float x) {
        // Check validity of x and p
        if (x < 0 || x > 1.0f) return Double.NaN;
        double t = tsat(p);
        if (t == Double.NaN) return t;
        return (1-x) * Region1.h(p, t) + x * Region2.hsat(p, t);
    }

    /** Specific enthalpy for Wet Steam
     *
     * @param x dryness fraction
     * @param t absolute temperature
     * @return double
     */

    public static double hsat(float x, double t) {
        // Check validity of x and t
        if (x < 0 || x > 1.0f) return Double.NaN;
        double p = psat(t);
        if (p == Double.NaN) return p;
        return (1-x) * Region1.h(p, t) + x * Region2.hsat(p, t);
    }

    /** Specific isochoric heat capacity cv(p,T)
     *
     * @param p pressure
     * @param t absolute temperature
     * @return double
     */
    public static double cv(double p, double t) {
        switch (getRegion(p,t)) {
            case 1:
                return Region1.cv(p,t);
            case 2:
                return Region2.cv(p,t);
            case 3:
                return Region3.cv(p,t);
            case 5:
                return Region5.cv(p, t);
            default:
                return Double.NaN;
        }
    }

    /** Specific isochoric heat capacity for Wet Steam
     *
     * @param p pressure
     * @param x dryness fraction
     * @return double
     */

    public static double cvsat(double p, float x) {
        // Check validity of x and p
        if (x < 0 || x > 1.0f) return Double.NaN;
        double t = tsat(p);
        if (t == Double.NaN) return t;
        return (1-x) * Region1.cv(p, t) + x * Region2.cvsat(p, t);
    }

    /** Specific isochoric heat capacity for Wet Steam
     *
     * @param x dryness fraction
     * @param t absolute temperature
     * @return double
     */

    public static double cvsat(float x, double t) {
        // Check validity of x and t
        if (x < 0 || x > 1.0f) return Double.NaN;
        double p = psat(t);
        if (p == Double.NaN) return p;
        return (1-x) * Region1.cv(p, t) + x * Region2.cvsat(p, t);
    }

    /** Specific isobaric heat capacity cp(p,T)
     *
     * @param p pressure
     * @param t absolute temperature
     * @return double
     */
    public static double cp(double p, double t) {
        switch (getRegion(p,t)) {
            case 1:
                return Region1.cp(p,t);
            case 2:
                return Region2.cp(p,t);
            case 3:
                return Region3.cp(p,t);
            case 5:
                return Region5.cp(p, t);
            default:
                return Double.NaN;
        }
    }

    /** Specific isobaric heat capacity for Wet Steam
     *
     * @param p pressure
     * @param x dryness fraction
     * @return double
     */

    public static double cpsat(double p, float x) {
        // Check validity of x and p
        if (x < 0 || x > 1.0f) return Double.NaN;
        double t = tsat(p);
        if (t == Double.NaN) return t;
        return (1-x) * Region1.cp(p, t) + x * Region2.cpsat(p, t);
    }

    /** Specific isobaric heat capacity for Wet Steam
     *
     * @param x dryness fraction
     * @param t absolute temperature
     * @return double
     */

    public static double cpsat(float x, double t) {
        // Check validity of x and t
        if (x < 0 || x > 1.0f) return Double.NaN;
        double p = psat(t);
        if (p == Double.NaN) return p;
        return (1-x) * Region1.cp(p, t) + x * Region2.cpsat(p, t);
    }

    /** Speed of sound w(p,T)
     *
     * @param p pressure
     * @param t absolute temperature
     * @return double
     */
    public static double w(double p, double t) {
        switch (getRegion(p,t)) {
            case 1:
                return Region1.w(p, t);
            case 2:
                return Region2.w(p, t);
            case 3:
                return Region3.w(p, t);
            case 5:
                return Region5.w(p, t);
            default:
                return Double.NaN;
        }
    }

    /** Speed of sound for Wet Steam
     *
     * @param p pressure
     * @param x dryness fraction
     * @return double
     */

    public static double wsat(double p, float x) {
        // Check validity of x and p
        if (x < 0 || x > 1.0f) return Double.NaN;
        double t = tsat(p);
        if (t == Double.NaN) return t;
        return (1-x) * Region1.w(p, t) + x * Region2.wsat(p, t);
    }

    /** Speed of sound for Wet Steam
     *
     * @param x dryness fraction
     * @param t absolute temperature
     * @return double
     */

    public static double wsat(float x, double t) {
        // Check validity of x and t
        if (x < 0 || x > 1.0f) return Double.NaN;
        double p = psat(t);
        if (p == Double.NaN) return p;
        return (1-x) * Region1.w(p, t) + x * Region2.wsat(p, t);
    }

    /** Specific internal energy u(p,T)
     *
     * @param p pressure
     * @param t absolute temperature
     * @return double
     */
    public static double u(double p, double t) {
        switch (getRegion(p,t)) {
            case 1:
                return Region1.u(p,t);
            case 2:
                return Region2.u(p,t);
            case 3:
                return Region3.u(p,t);
            case 5:
                return Region5.u(p, t);
            default:
                return Double.NaN;
        }
    }

    /** Specific internal energy for Wet Steam
     *
     * @param p pressure
     * @param x dryness fraction
     * @return double
     */

    public static double usat(double p, float x) {
        // Check validity of x and p
        if (x < 0 || x > 1.0f) return Double.NaN;
        double t = tsat(p);
        if (t == Double.NaN) return t;
        return (1-x) * Region1.u(p, t) + x * Region2.usat(p, t);
    }

    /** Specific Internal Energy for Wet Steam
     *
     * @param x dryness fraction
     * @param t absolute temperature
     * @return double
     */

    public static double usat(float x, double t) {
        // Check validity of x and t
        if (x < 0 || x > 1.0f) return Double.NaN;
        double p = psat(t);
        if (p == Double.NaN) return p;
        return (1-x) * Region1.u(p, t) + x * Region2.usat(p, t);
    }

    /** Specific entropy s(p,T)
     *
     * @param p pressure
     * @param t absolute temperature
     * @return double
     */
    public static double s(double p, double t) {
        switch (getRegion(p,t)) {
            case 1:
                return Region1.s(p,t);
            case 2:
                return Region2.s(p,t);
            case 3:
                return Region3.s(p,t);
            case 5:
                return Region5.s(p, t);
            default:
                return Double.NaN;
        }
    }

    /** Specific entropy for Wet Steam
     *
     * @param p pressure
     * @param x dryness fraction
     * @return double
     */

    public static double ssat(double p, float x) {
        // Check validity of x and p
        if (x < 0 || x > 1.0f) return Double.NaN;
        double t = tsat(p);
        if (t == Double.NaN) return t;
        return (1-x) * Region1.s(p, t) + x * Region2.ssat(p, t);
    }

    /** Specific entropy for Wet Steam
     *
     * @param x dryness fraction
     * @param t absolute temperature
     * @return double
     */

    public static double ssat(float x, double t) {
        // Check validity of x and t
        if (x < 0 || x > 1.0f) return Double.NaN;
        double p = psat(t);
        if (p == Double.NaN) return p;
        return (1-x) * Region1.s(p, t) + x * Region2.ssat(p, t);
    }

    /** Fundamental equation: specific Helmholtz Free Energy
     *  f(ρ,T)
     *
     * @param rho density
     * @param t absolute temperature
     * @return double
     */
    public static double f(double rho, double t) {
        return Region3.f(rho, t);
    }

    /** Fundamental equation: Saturation Pressure Ps(T)
     *
     * @param t absolute temperature
     * @return double
     */
    public static double psat(double t) {
        double p__ = 1; // MPa
        double t__ = 1; // K
        // test for validity of temperature:
        if (t < 273.15 || t > Tc) return Double.NaN;

        double theta = t/t__ + nsat[9-1]/(t/t__ - nsat[10-1]);
        double a = theta*theta + nsat[1-1]*theta + nsat[2-1];
        double b = nsat[3-1]*theta*theta + nsat[4-1]*theta + nsat[5-1];
        double c = nsat[6-1]*theta*theta + nsat[7-1]*theta + nsat[8-1];
        double pp = 2*c/(-b + Math.sqrt(b*b - 4*a*c));
        return p__ * pp*pp*pp*pp;
    }

    /** Backward equation: Saturation Temperature Ts(p)
     * 
     * @param p pressure
     * @return
     */
    public static double tsat(double p) {
        double p__ = 1; // MPa
        double t__ = 1; // K
        // test for validity of pressure:
        if (p < 611.213E-6 || p > Pc) return Double.NaN;

        double beta = Math.pow(p/p__, .25);
        double e = beta*beta + nsat[3-1]*beta + nsat[6-1];
        double f = nsat[1-1]*beta*beta + nsat[4-1]*beta + nsat[7-1];
        double g = nsat[2-1]*beta*beta + nsat[5-1]*beta + nsat[8-1];
        double d = 2*g/(-f - Math.sqrt(f*f - 4*e*g));

        return (nsat[10-1] + d - Math.sqrt((nsat[10-1]+d)*(nsat[10-1]+d) - 4*(nsat[9-1] + nsat[10-1]*d)))*t__/2;
    }

    /** Backward equation: Temperature T(p,h)
     *
     * @param p pressure
     * @param h specific enthalpy
     * @return double
     */
    public static double t_ph(double p, double h) {
        // Need to determine in which region the (p,h) coordinates lie
        if (p <= psat(623.15)) {
            double h1 = hsat(p, 0f);
            double h11 = hsat(p, 1.0f);
            if (h < h1) return Region1.t_ph(p, h);
            if (h > h11) return Region2.t_ph(p, h);
            return tsat(p); // 2-phase region
        }

        double h13 = Region1.h(p, 623.15);
        double t = b23_(p);
        double h23 = Region2.h(p, t);

        if (h <= h13) return Region1.t_ph(p, h);
        if (h >= h23) return Region2.t_ph(p, h);
        if (p > Pc) return Region3.t_ph(p, h);

        // Test for 2-phase region
        double h1 = hsat(p, 0f);
        double h11 = hsat(p, 1.0f);
        if (h < h1 || h > h11) return Region3.t_ph(p, h);
        return tsat(p);
    }

    /** Backward equation: Temperature T(p,s)
     *
     * @param p pressure
     * @param s specific entropy
     * @return double
     */
    public static double t_ps(double p, double s) {
        // Need to determine in which region the (p,h) coordinates lie
        if (p <= psat(623.15)) {
            double s1 = ssat(p, 0f);
            double s11 = ssat(p, 1.0f);
            if (s < s1) return Region1.t_ps(p, s);
            if (s > s11) return Region2.t_ps(p, s);
            return tsat(p); // 2-phase region
        }

        double s13 = Region1.s(p, 623.15);
        double t = b23_(p);
        double s23 = Region2.s(p, t);

        if (s <= s13) return Region1.t_ps(p, s);
        if (s >= s23) return Region2.t_ps(p, s);
        if (p > Pc) return Region3.t_ps(p, s);

        // Test for 2-phase region
        double s1 = ssat(p, 0f);
        double s11 = ssat(p, 1.0f);
        if (s < s1 || s > s11) return Region3.t_ps(p, s);
        return tsat(p);
    }

    /**
     * Tsat(h, s) for the 2-phase region
     * @param h specific enthalpy
     * @param s specific entropy
     * @return temperature
     */
    public static double tsat_hs(double h, double s) {
        int i_[] = {0,0,0,1,1,1,1,2,2,2,3,3,3,3,4,4,5,5,5,5,6,6,6,8,10,10,12,14,
        14,16,16,18,18,18,20,28};
        int j_[] = {0,3,12,0,1,2,5,0,5,8,0,2,3,4,0,1,1,2,4,16,6,8,22,1,20,36,24,
        1,28,12,32,14,22,36,24,36};
        double n_[] = {0.179882673606601,-0.267507455199603,1.16276722612600,
        0.147545428713616,-0.512871635973248,0.421333567697984,0.563749522189870,
        0.429274443819153,-3.35704552142140,10.8890916499278,-0.248483390456012,
        0.304153221906390,-0.494819763939905,1.07551674933261,7.33888415457688E-2,
        1.40170545411085E-2,-0.106110975998808,1.68324361811875E-2,1.25028363714877,
        1013.16840309509,-1.51791558000712,52.4277865990866,2.30495545563912E4,
        2.49459806365456E-2,2.10796467412137E6,3.66836848613065E8,-1.44814105365163E8,
        -1.79276373003590E-3,4.89955602100459E9,471.262212070518,-8.29294390198652E10,
        -1715.45662263191,3.55777682973575E6,5.86062760258436E11,-1.29887635078195E7,
        3.17247449371057E10};
        double t__ = 550; // K
        double h__ = 2800; // kJ/kg
        double s__ = 9.2; // kJ/(kg K)
        double eta = h/h__;
        double sigma = s/s__;
        double theta = 0;

        for (int k=0; k<36; k++) {
            theta += n_[k] * Math.pow(eta-0.119, i_[k]) * Math.pow(sigma-1.07, j_[k]);
        }

        return theta * t__;
    }

    /**
     * Pressure in 2 phase region p(h,s)
     * @param h
     * @param s
     * @return
     */
    public static double psat_hs(double h, double s) {
        return psat(tsat_hs(h, s));
    }

    /**
     * Dryness fraction x(h,s)
     * @param h specific enthalpy
     * @param s specific entropy
     * @return
     */
    public static double x_hs(double h, double s) {
        double tsat = tsat_hs(h,s);
//        double psat = psat(tsat);
//        double h1 = Region1.h(psat, tsat);
//        double h2 = Region2.h(h, s);
        double h1 = hsat(0f, tsat);
        double h2 = hsat(1f, tsat);
        return (h - h1)/(h2 - h1);
    }

    /**
     * Backward equation: p(h, s)
     * @param h specific enthalpy
     * @param s specific entropy
     * @return pressure
     */
    public static double p_hs(double h, double s) {
        // First identify which region the (h,s) coordinates lie in
        switch (getRegion_hs(h, s)) {
            case 1:
                return Region1.p_hs(h, s);
            case 2:
                return Region2.p_hs(h, s);
            case 3:
                return Region3.p_hs(h, s);
            case 4:
                return psat_hs(h, s);
            default:
                return Double.NaN;
        }
    }

    /** Determines which region of the p-T diagram
     *  the given (p,T) coordinates lie.
     * Required in order to select the right set of equations
     * for calculation.
     * @param p pressure (MPa)
     * @param t absolute temperature (K)
     * @return
     */
    protected static int getRegion(double p, double t) { // TODO: how to determine metastable region?
        if (t < 273.15) return -1;
        if (p > 100 || p <= 0) return -1;
        if (t > 623.15) { // test for regions 2,3,5
            if (t > 1073.15) { // test for region 5
                if (t <= 2273.15 && p <= 50) return 5;
                else return -1;
            }
            else { // test for regions 2 or 3
                if (t > 863.15) return 2;
                double p23 = b23(t);
                if (p23 < p) return 3;
                else return 2;
            }
        }
        else { // test for region 1
            double ps = psat(t);
            if (p > ps) return 1;
            if (p < ps) return 2;
            return 4; // saturation line
        }
    }

    /**
     * Determines which region of the h-s diagram any given (h,s) coordinates
     * lie.
     * Required in order to select the right set of backward equations for
     * calculation
     * @param h specific enthalpy
     * @param s specific entropy
     * @return
     */
    protected static int getRegion_hs(double h, double s) {
        if (s > 9.155759395) return 2; // because h"2ab(s) doesn't apply beyond
                                        // this value of s

        if (s >= 5.85) {
            if (h > hpprime2ab_s(s)) return 2;
            else return 4;
        }

        if (s >= 4.41202148223476) {
            if (h > hpprime2c3b_s(s)) {
                if (s <= 5.048096828 || h <= 2.563592004E3) return 3;
                if (s >= 5.260578707 || h >= 2.812942061E3) return 2;
                
                double t = tb23_hs(h, s);
                double p = Region2.p_hs2c(h, s);
                if (p > b23(t)) return 3;
                else return 2;
            }
            else return 4;
        }

        if (s > 3.778281340) {
            if (h > hprime3a_s(s)) return 3;
            else return 4;
        }

        if (h > hprime1_s(s)) {
            if (s > 3.397782955) {
                if (h > hb13_s(s)) return 3;
            }
            return 1;
        }
        else return 4;
    }

    // Boundary Equations

    /**
     * Boundary between region 2 and region 3
     * @param t temperature
     * @return
     */
    protected static double b23(double t) {
        double n1 = 348.05185628969;
        double n2 = -1.1671859879975;
        double n3 = 1.0192970039326E-3;

        return n1 + n2*t + n3*t*t;
    }

    /**
     * Boundary between region 2 and region 3
     * @param p pressure
     * @return
     */
    protected static double b23_(double p) {
        double n5 = 13.918839778870;
        double n4 = 572.54459862746;
        double n3 = 1.0192970039326E-3;

        return n4 + Math.sqrt((p-n5)/n3);
    }

    /**
     * Saturated liquid line between Region 1 and 2-phase Region 4
     * @param s specific entropy
     * @return specific enthalpy h'
     */
    protected static double hprime1_s(double s) {
        int i_[] = {0,0,1,1,2,2,3,3,4,4,4,5,5,7,8,12,12,14,14,16,20,20,22,24,28,
        32,32};
        int j_[] = {14,36,3,16,0,5,4,36,4,16,24,18,24,1,4,2,4,1,22,10,12,28,8,3,
        0,6,8};
        double n_[] = {0.332171191705237,6.11217706323496E-4,-8.82092478906822,
        -0.455628192543250,-2.63483840850452E-5,-22.3949661148062,-4.28398660164013,
        -0.616679338856916,-14.6823031104040,284.523138727299,-113.398503195444,
        1156.71380760859,395.551267359325,-1.54891257229285,19.4486637751291,
        -3.57915139457043,-3.35369414148819,-0.664426796332460,3.23321885383934E4,
        3317.66744667084,-2.23501257931087E4,5.73953875852936E6,173.226193407919,
        -3.63968822121321E-2,8.34596332878346E-7,5.03611916682674,65.5444787064505};
        double h__ = 1700; // kJ/kg
        double s__ = 3.8; // kJ/(kg K)
        double sigma = s/s__;
        double eta = 0;

        for (int k=0; k<27; k++) {
            eta += n_[k] * Math.pow(sigma-1.09, i_[k]) * Math.pow(sigma+3.66E-5, j_[k]);
        }

        return eta * h__;
    }

    /**
     * Saturated liquid line between Region 3 and 2-phase region 4
     * @param s
     * @return specific enthalpy h'
     */
    protected static double hprime3a_s(double s) {
        int i_[] = {0,0,0,0,2,3,4,4,5,5,6,7,7,7,10,10,10,32,32};
        int j_[] = {1,4,10,16,1,36,3,16,20,36,4,2,28,32,14,32,36,0,6};
        double n_[] = {0.822673364673336,0.181977213534479,-1.12000260313624E-2,
        -7.46778287048033E-4,-0.179046263257381,4.24220110836657E-2,-0.341355823438768,
        -2.09881740853565,-8.22477343323596,-4.99684082076008,0.191413958471069,
        5.81062241093136E-2,-1655.05498701029,1588.70443421201,-85.0623535172818,
        -3.17714386511207E4,-9.45890406632871E4,-1.39273847088690E-6,0.631052532240980};
        double h__ = 1700; // kJ/kg
        double s__ = 3.8; // kJ/(kg K)
        double sigma = s/s__;
        double eta = 0;

        for (int k=0; k<19; k++) {
            eta += n_[k] * Math.pow(sigma-1.09, i_[k]) * Math.pow(sigma+3.66E-5, j_[k]);
        }

        return eta * h__;
    }

    /**
     * Saturated vapour line adjacent to regions 2a and 2b
     * @param s
     * @return specific enthalpy h"
     */
    protected static double hpprime2ab_s(double s) {
        int i_[] = {1,1,2,2,4,4,7,8,8,10,12,12,18,20,24,28,28,28,28,28,32,32,32,
        32,32,36,36,36,36,36};
        int j_[] = {8,24,4,32,1,2,7,5,12,1,0,7,10,12,32,8,12,20,22,24,2,7,12,14,
        24,10,12,20,22,28};
        double n_[] = {-524.581170928788,-9.26947218142218E6,-237.385107491666,
        2.10770155812776E10,-23.9494562010986,221.802480294197,-5.10472533393438E6,
        1.24981396109147E6,2.00008436996201E9,-815.158509791035,-157.612685637523,
        -1.14200422332791E10,6.62364680776872E15,-2.27622818296144E18,-1.71048081348406E31,
        6.60788766938091E15,1.66320055886021E22,-2.18003784381501E29,-7.87276140295618E29,
        1.51062329700346E31,7.95732170300541E6,1.31957647355347E15,-3.25097068299140E23,
        -4.18600611419248E25,2.97478906557467E34,-9.53588761745473E19,1.66957699620939E24,
        -1.75407764869978E32,3.47581490626396E34,-7.10971318427851E38};
        double h__ = 2800; // kJ/kg
        double s_1 = 5.21; // kJ/(kg K)
        double s_2 = 9.2; // kJ/(kg K)
        double sigma1 = s/s_1;
        double sigma2 = s/s_2;
        double eta = 0;

        for (int k=0; k<30; k++) {
            eta += n_[k] * Math.pow(1/sigma1 -0.513, i_[k]) * Math.pow(sigma2 -0.524, j_[k]);
        }

        return h__ * Math.exp(eta);
    }

    /**
     * Saturated vapor line between regions 2c and 3b
     * @param s
     * @return specific enthalpy h"
     */
    protected static double hpprime2c3b_s(double s) {
        int i_[] = {0,0,0,1,1,5,6,7,8,8,12,16,22,22,24,36};
        int j_[] = {0,3,4,0,12,36,12,16,2,20,32,36,2,32,7,20};
        double n_[] = {1.04351280732769,-2.27807912708513,1.80535256723202,
        0.420440834792042,-1.05721244834660E5,4.36911607493884E24,-3.28032702839753E11,
        -6.78686760804270E15,7439.57464645363,-3.56896445355761E19,1.67590585186801E31,
        -3.55028625419105E37,3.96611982166538E11,-4.14716268484468E40,3.59080103867382E18,
        -1.16994334851995E40};
        double h__ = 2800; // kJ/kg
        double s__ = 5.9; // kJ/(kg K)
        double sigma = s/s__;
        double eta = 0;

        for (int k=0; k<16; k++) {
            eta += n_[k] * Math.pow(sigma-1.02, i_[k]) * Math.pow(sigma-0.726, j_[k]);
        }

        return h__ * eta*eta*eta*eta;
    }

    /**
     * 623.15 K Isotherm between Region 1 and Region 3
     * @param s
     * @return specific enthalpy
     */
    protected static double hb13_s(double s) {
        int i_[] = {0,1,1,3,5,6};
        int j_[] = {0,-2,2,-12,-4,-3};
        double n_[] = {0.913965547600543,-4.30944856041991E-5,60.3235694765419,
        1.17518273082168E-18,0.220000904781292,-69.0815545851641};
        double h__ = 1700; // kJ/kg
        double s__ = 3.8;  // kJ/(kg K)
        double sigma = s/s__;
        double eta = 0;

        for (int k=0; k<6; k++) {
            eta += n_[k] * Math.pow(sigma-0.884, i_[k]) * Math.pow(sigma-0.864, j_[k]);
        }

        return eta * h__;
    }

    protected static double tb23_hs(double h, double s) {
        int i_[] = {-12,-10,-8,-4,-3,-2,-2,-2,-2,0,1,1,1,3,3,5,6,6,8,8,8,12,12,
        14,14};
        int j_[] = {10,8,3,4,3,-6,2,3,4,0,-3,-2,10,-2,-1,-5,-6,-3,-8,-2,-1,-12,
        -1,-12,1};
        double n_[] = {6.29096260829810E-4,-8.23453502583165E-4,5.15446951519474E-8,
        -1.17565945784945,3.48519684726192,-5.07837382408313E-12,-2.84637670005479,
        -2.36092263939673,6.01492324973779,1.48039650824546,3.60075182221907E-4,
        -1.26700045009952E-2,-1.22184332521413E6,0.149276502463272,0.698733471798484,
        -2.52207040114321E-2,1.47151930985213E-2,-1.08618917681849,-9.36875039816322E-4,
        81.9877897570217,-182.041861521835,2.61907376402688E-6,-2.91626417025961E4,
        1.40660774926165E-5,7.83237062349385E6};
        double t__ = 900; // K
        double h__ = 3000; // kJ/kg
        double s__ = 5.3; // kJ/(kg K)
        double eta = h/h__;
        double sigma = s/s__;
        double theta = 0;

        for (int k=0; k<25; k++) {
            theta += n_[k] * Math.pow(eta-0.727, i_[k]) * Math.pow(sigma-0.864, j_[k]);
        }

        return theta * t__;
    }
}
