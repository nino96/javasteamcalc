/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package com.spayc;

/**
 *
 * @author Kehinde
 */
class Region1 {
    // Shared coefficients for Region I
    protected static int i[] = {0,0,0,0,0,0,0,0,1,1,1,1,1,1,2,2,2,2,2,3,3,3,4,
    4,4,5,8,8,21,23,29,30,31,32};
    protected static int j[] = {-2,-1,0,1,2,3,4,5,-9,-7,-1,0,1,3,-3,0,1,3,17,
    -4,0,6,-5,-2,10,-8,-11,-6,-29,-31,-38,-39,-40,-41};
    protected static double n[] = {0.14632971213167,-0.84548187169114,
    -3.7563603672040,3.3855169168385,-0.95791963387872,0.15772038513228,
    -1.6616417199501E-2,8.1214629983568E-4,2.8319080123804E-4,-6.0706301565874E-4,
    -1.8990068218419E-2,-3.2529748770505E-2,-2.1841717175414E-2,-5.2838357969930E-5,
    -4.7184321073267E-4,-3.0001780793026E-4,4.7661393906987E-5,-4.4141845330846E-6,
    -7.2694996297594E-16,-3.1679644845054E-5,-2.8270797985312E-6,-8.5205128120103E-10,
    -2.2425281908000E-6,-6.5171222895601E-7,-1.4341729937924E-13,-4.0516996860117E-7,
    -1.2734301741641E-9,-1.7424871230634E-10,-6.8762131295531E-19,1.4478307828521E-20,
    2.6335781662795E-23,-1.1947622640071E-23,1.8228094581404E-24,-9.3537087292458E-26};
    protected final static double p_ = 16.53; // MPa
    protected final static double t_ = 1386;  // K

    // Internal function γ(π,τ): reduced Gibbs Free Energy = g/RT
    static double gamma(double pi, double tau) {
        double gamma = 0;

        for (int k=0; k<34; k++) {
            gamma += n[k] * Math.pow((7.1-pi), i[k]) * Math.pow((tau-1.222), j[k]);
        }

        return gamma;
    }

    // Internal function γ_π(π,τ) (partial differential wrt π)
    static double gamma_pi(double pi, double tau) {
        double gamma_pi = 0;

        for (int k=0; k<34; k++) {
            gamma_pi += -n[k] * i[k] * Math.pow((7.1-pi), (i[k]-1)) * Math.pow((tau-1.222), j[k]);
        }

        return gamma_pi;
    }

    // Internal function γ_τ(π,τ) (partial differential wrt τ)
    static double gamma_tau(double pi, double tau) {
        double gamma_tau = 0;

        for (int k=0; k<34; k++) {
            gamma_tau += n[k] * j[k] * Math.pow((7.1-pi), i[k]) * Math.pow((tau-1.222), (j[k]-1));
        }

        return gamma_tau;
    }

    // Internal function γ_πτ(π,τ) (partial differential wrt π and τ)
    static double gamma_pi_tau(double pi, double tau) {
        double gamma_pi_tau = 0;

        for (int k=0; k<34; k++) {
            gamma_pi_tau += -n[k] * i[k] * j[k] * Math.pow((7.1-pi), (i[k]-1)) * Math.pow((tau-1.222), (j[k]-1));
        }

        return gamma_pi_tau;
    }

    // Internal function γ_ππ(π,τ) (2nd-order partial differential wrt π)
    static double gamma_pi2(double pi, double tau) {
        double gamma_pi2 = 0;

        for (int k=0; k<34; k++) {
            gamma_pi2 += n[k] * i[k] * (i[k]-1) * Math.pow((7.1-pi), (i[k]-2)) * Math.pow((tau-1.222), j[k]);
        }

        return gamma_pi2;
    }

    // Internal function γ_ττ(π,τ) (2nd-order partial differential wrt τ)
    static double gamma_tau2(double pi, double tau) {
        double gamma_tau2 = 0;

        for (int k=0; k<34; k++) {
            gamma_tau2 += n[k] * j[k] * (j[k]-1) * Math.pow((7.1-pi), i[k]) * Math.pow((tau-1.222), (j[k]-2));
        }

        return gamma_tau2;
    }

    // Gibbs Free Energy g(p,t)
    public static double g(double p, double t) {
        double pi = p/p_;
        double tau = t_/t;

        return gamma(pi,tau) * SteamCalc.R * t;
    }

    // Specific volume v(p,t)
    public static double v(double p, double t) {
        double pi = p/p_;
        double tau = t_/t;

        return pi * gamma_pi(pi,tau) * SteamCalc.R * t / (1000*p); // must adjust by factor of 10^3. Dunno why.
    }

    // Specific internal energy u(p,t)
    public static double u(double p, double t) {
        double pi = p/p_;
        double tau = t_/t;

        return (tau * gamma_tau(pi,tau) - pi * gamma_pi(pi,tau)) * SteamCalc.R * t;
    }

    // Specific entropy s(p,t)
    public static double s(double p, double t) {
        double pi = p/p_;
        double tau = t_/t;

        return (tau * gamma_tau(pi,tau) - gamma(pi,tau)) * SteamCalc.R;
    }

    // Specific enthalpy h(p,t)
    public static double h(double p, double t) {
        double pi = p/p_;
        double tau = t_/t;

        return tau * gamma_tau(pi,tau) * SteamCalc.R * t;
    }

    // Specific isobaric heat capacity cp(p,t)
    public static double cp(double p, double t) {
        double pi = p/p_;
        double tau = t_/t;

        return -tau*tau * gamma_tau2(pi,tau) * SteamCalc.R;
    }

    // Specific isochoric heat capacity cv(p,t)
    public static double cv(double p, double t) {
        double pi = p/p_;
        double tau = t_/t;

        return (-tau*tau * gamma_tau2(pi,tau)
                + Math.pow((gamma_pi(pi,tau) - tau * gamma_pi_tau(pi,tau)), 2)
                /gamma_pi2(pi,tau))* SteamCalc.R;
    }

    // Speed of sound w(p,t)
    public static double w(double p, double t) {
        double pi = p/p_;
        double tau = t_/t;

        double a = gamma_pi(pi,tau);
        double b = a - tau * gamma_pi_tau(pi, tau);
        double c = a*a/(b*b/(tau*tau*gamma_tau2(pi,tau)) - gamma_pi2(pi,tau));
//        return gamma_pi(pi,tau) * Math.sqrt(SteamCalc.R*t
//                / ((Math.pow((gamma_pi(pi,tau) - tau * gamma_pi_tau(pi,tau)), 2)
//                /tau/tau/gamma_tau2(pi,tau)) - gamma_pi2(pi,tau)));
        return Math.sqrt(10 * SteamCalc.R * t * c); // for some reason, must adjust by factor of 10. Dunno why.
    }

    // Backward equation T(p,h)
    public static double t_ph(double p, double h) {
        int i_[] = {0,0,0,0,0,0,1,1,1,1,1,1,1,2,2,3,3,4,5,6};
        int j_[] = {0,1,2,6,22,32,0,1,2,3,4,10,32,10,32,10,32,32,32,32};
        double n_[] = {-2.3872489924521E2,4.0421188637945E2,1.1349746881718E2,
        -5.8457616048039,-1.5285482413140E-4,-1.0866707695377E-6,-1.3391744872602E1,
        4.3211039183559E1,-5.4010067170506E1,3.0535892203916E1,-6.5964749423638,
        9.3965400878363E-3,1.1573647505340E-7,-2.5858641282073E-5,-4.0644363084799E-9,
        6.6456186191635E-8,8.0670734103027E-11,-9.3477771213947E-13,
        5.8265442020601E-15,-1.5020185953503E-17};
        double p__ = 1;    // MPa
        double t__ = 1;    // K
        double h_ = 2500; // kJ/kg
        double pi = p/p__;
        double eta = h/h_;
        double tau = 0;

        for (int k=0; k<20; k++) {
            tau += n_[k] * Math.pow(pi, i_[k]) * Math.pow(eta+1, j_[k]);
        }
        return tau*t__;
    }

    // Backward equation T(p,s)
    public static double t_ps(double p, double s) {
        int i_[] = {0,0,0,0,0,0,1,1,1,1,1,1,2,2,2,2,2,3,3,4};
        int j_[] = {0,1,2,3,11,31,0,1,2,3,12,31,0,1,2,9,31,10,32,32};
        double n_[] = {1.7478268058307E2,3.4806930892873E1,6.5292584978455,
        0.33039981775489,-1.9281382923196E-7,-2.4909197244573E-23,-0.26107636489332,
        0.22592965981586,-6.4256463395226E-2,7.8876289270526E-3,3.5672110607366E-10,
        1.7332496994895E-24,5.6608900654837E-4,-3.2635483139717E-4,4.4778286690632E-5,
        -5.1322156908507E-10,-4.2522657042207E-26,2.6400441360689E-13,
        7.8124600459723E-29,-3.0732199903668E-31};
        double p__ = 1; // MPa
        double t__ = 1; // K
        double s_ = 1;  // kJ/kgK
        double pi = p/p__;
        double sigma = s/s_;
        double tau = 0;

        for (int k=0; k<20; k++) {
            tau += n_[k] * Math.pow(pi, i_[k]) * Math.pow(sigma+2, j_[k]);
        }
        return tau*t__;
    }

    public static double p_hs(double h, double s) {
        int i_[] = {0,0,0,0,0,0,0,0,1,1,1,1,2,2,2,3,4,4,5};
        int j_[] = {0,1,2,4,5,6,8,14,0,1,4,6,0,1,10,4,1,4,0};
        double n_[] = {-0.691997014660582,-18.3612548787560,-9.28332409297335,
        65.9639569909906,-16.2060388912024,450.620017338667,854.680678224170,
        6075.23214001162,32.6487682621856,-26.9408844582931,-319.947848334300,
        -928.354307043320,30.3634537455249,-65.0540422444146,-4309.91316516130,
        -747.512324096068,730.000345529245,1142.84032569021,-436.407041874559};
        double p__ = 100; // MPa
        double h__ = 3400; // kJ/kg
        double s__ = 7.6; // kJ/(kg K)
        double eta = h/h__;
        double sigma = s/s__;
        double pi = 0;

        for (int k=0; k<19; k++) {
            pi += n_[k] * Math.pow(eta+0.05, i_[k]) * Math.pow(sigma+0.05, j_[k]);
        }

        return pi * p__;
    }

    public static double t_hs(double h, double s) {
        return t_ph(p_hs(h,s), h);
    }
}
