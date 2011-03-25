/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package com.spayc;

/**
 *
 * @author Kehinde
 */
class Region2 {
    // Shared constants for region 2
    // Ideal gas component
    final static int jo[] = {0,1,-5,-4,-3,-2,-1,2,3};
    final static double no[] = {-9.6927686500217,1.0086655968018E1,-5.6087911283020E-3,
    7.1452738081455E-2,-0.40710498223928,1.4240819171444,-4.3839511319450,
    -0.28408632460772,2.1268463753307E-2};
    final static double no1 = -9.6937268393049; // no1 and no2 to be used in
    final static double no2 = 1.0087275970006E1; // metastable region (saturated steam)
    final static double p_ = 1; // p* (MPa)
    final static double t_ = 540; // T* (K)
    // Residual gas component
    final static int i[] = {1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,5,6,6,6,7,7,7,8,
    8,9,10,10,10,16,16,18,20,20,20,21,22,23,24,24,24};
    final static int j[] = {0,1,2,3,6,1,2,4,7,36,0,1,3,6,35,1,2,3,7,3,16,35,0,11,
    25,8,36,13,4,10,14,29,50,57,20,35,48,21,53,39,26,40,58};
    final static double n[] = {-1.7731742473213E-3,-1.7834862292358E-2,-4.5996013696365E-2,
    -5.7581259083432E-2,-5.0325278727930E-2,-3.3032641670203E-5,-1.8948987516315E-4,
    -3.9392777243355E-3,-4.3797295650573E-2,-2.6674547914087E-5,2.0481737692309E-8,
    4.3870667284435E-7,-3.2277677238570E-5,-1.5033924542148E-3,-4.0668253562649E-2,
    -7.8847309559367E-10,1.2790717852285E-8,4.8225372718507E-7,2.2922076337661E-6,
    -1.6714766451061E-11,-2.1171472321355E-3,-2.3895741934104E1,-5.9059564324270E-18,
    -1.2621808899101E-6,-3.8946842435739E-2,1.1256211360459E-11,-8.2311340897998,
    1.9809712802088E-8,1.0406965210174E-19,-1.0234747095929E-13,-1.0018179379511E-9,
    -8.0882908646985E-11,0.10693031879409,-0.33662250574171,8.9185845355421E-25,
    3.0629316876232E-13,-4.2002467698208E-6,-5.9056029685639E-26,3.7826947613457E-6,
    -1.2768608934681E-15,7.3087610595061E-29,5.5414715350778E-17,-9.4369707241210E-7};
    // Residual gas component for metastable region (saturated steam)
    final static int i_m[] = {1,1,1,1,2,2,2,3,3,4,4,5,5};
    final static int j_m[] = {0,2,5,11,1,7,16,4,16,7,10,9,10};
    final static double n_m[] = {-7.3362260186506E-3,-8.8223831943146E-2,
    -7.2334555213245E-2,-4.0813178534455E-3,2.0097803380207E-3,-5.3045921898642E-2,
    -7.6190409086970E-3,-6.3498037657313E-3,-8.6043093028588E-2,7.5321581522770E-3,
    -7.9238375446139E-3,-2.2888160778447E-4,-2.6456501482810E-3};
    // Internal function γo(π,τ): reduced Gibbs Free Energy = g/RT (ideal gas part)
    static double gammao(double pi, double tau, boolean meta) {
        double gamma = Math.log(pi);
        gamma += (meta? no1: no[1-1]) * Math.pow(tau, jo[1-1]);
        gamma += (meta? no2: no[2-1]) * Math.pow(tau, jo[2-1]);

        for (int k=2; k<9; k++) {
            gamma += no[k] * Math.pow(tau, jo[k]);
        }
        return gamma;
    }

    // Internal function γr(π,τ): reduced Gibbs Free Energy = g/RT (residual gas part)
    static double gammar(double pi, double tau) {
        double gamma = 0;

        for (int k=0; k<43; k++) {
            gamma += n[k] * Math.pow(pi, i[k]) * Math.pow(tau-0.5, j[k]);
        }
        return gamma;
    }
    // metastable version of above (for saturated vapor)
    static double gammarm(double pi, double tau) {
        double gamma = 0;

        for (int k=0; k<13; k++) {
            gamma += n_m[k] * Math.pow(pi, i_m[k]) * Math.pow(tau-0.5, j_m[k]);
        }
        return gamma;
    }

    static double gammao_pi(double pi, double tau, boolean meta) {
        return 1/pi;
    }

    static double gammao_pi2(double pi, double tau, boolean meta) {
        return -1/(pi*pi);
    }

    static double gammao_tau(double pi, double tau, boolean meta) {
        double gamma = (meta? no1: no[1-1]) * jo[1-1] * Math.pow(tau, jo[1-1] - 1);
        gamma += (meta? no2: no[2-1]) * jo[2-1] * Math.pow(tau, jo[2-1] - 1);

        for (int k=2; k<9; k++) {
            gamma += no[k] * jo[k] * Math.pow(tau, jo[k] - 1);
        }
        return gamma;
    }

    static double gammao_tau2(double pi, double tau, boolean meta) {
        double gamma = (meta? no1: no[1-1]) * jo[1-1] * (jo[1-1] - 1) * Math.pow(tau, jo[1-1] - 2);
        gamma += (meta? no2: no[2-1]) * jo[2-1] * (jo[2-1] - 1) * Math.pow(tau, jo[2-1] - 2);

        for (int k=2; k<9; k++) {
            gamma += no[k] * jo[k] * (jo[k] - 1) * Math.pow(tau, jo[k] - 2);
        }
        return gamma;
    }

    static double gammao_pi_tau(double pi, double tau, boolean meta) {
        return 0;
    }

    static double gammar_pi(double pi, double tau) {
        double gamma = 0;

        for (int k=0; k<43; k++) {
            gamma += n[k] * i[k] * Math.pow(pi, i[k] - 1) * Math.pow(tau-0.5, j[k]);
        }
        return gamma;
    }

    static double gammar_pi2(double pi, double tau) {
        double gamma = 0;

        for (int k=0; k<43; k++) {
            gamma += n[k] * i[k] * (i[k] - 1) * Math.pow(pi, i[k] - 2) * Math.pow(tau-0.5, j[k]);
        }
        return gamma;
    }

    static double gammar_tau(double pi, double tau) {
        double gamma = 0;

        for (int k=0; k<43; k++) {
            gamma += n[k] * Math.pow(pi, i[k]) * j[k] * Math.pow(tau-0.5, j[k] - 1);
        }
        return gamma;
    }

    static double gammar_tau2(double pi, double tau) {
        double gamma = 0;

        for (int k=0; k<43; k++) {
            gamma += n[k] * Math.pow(pi, i[k]) * j[k] * (j[k] - 1) * Math.pow(tau-0.5, j[k] - 2);
        }
        return gamma;
    }

    static double gammar_pi_tau(double pi, double tau) {
        double gamma = 0;

        for (int k=0; k<43; k++) {
            gamma += n[k] * i[k] * Math.pow(pi, i[k] - 1) * j[k] * Math.pow(tau-0.5, j[k] - 1);
        }
        return gamma;
    }

    static double gammarm_pi(double pi, double tau) {
        double gamma = 0;

        for (int k=0; k<13; k++) {
            gamma += n_m[k] * i_m[k] * Math.pow(pi, i_m[k] - 1) * Math.pow(tau-0.5, j_m[k]);
        }
        return gamma;
    }

    static double gammarm_pi2(double pi, double tau) {
        double gamma = 0;

        for (int k=0; k<13; k++) {
            gamma += n_m[k] * i_m[k] * (i_m[k] - 1) * Math.pow(pi, i_m[k] - 2) * Math.pow(tau-0.5, j_m[k]);
        }
        return gamma;
    }

    static double gammarm_tau(double pi, double tau) {
        double gamma = 0;

        for (int k=0; k<13; k++) {
            gamma += n_m[k] * Math.pow(pi, i_m[k]) * j_m[k] * Math.pow(tau-0.5, j_m[k] - 1);
        }
        return gamma;
    }

    static double gammarm_tau2(double pi, double tau) {
        double gamma = 0;

        for (int k=0; k<13; k++) {
            gamma += n_m[k] * Math.pow(pi, i_m[k]) * j_m[k] * (j_m[k] - 1) * Math.pow(tau-0.5, j_m[k] - 2);
        }
        return gamma;
    }

    static double gammarm_pi_tau(double pi, double tau) {
        double gamma = 0;

        for (int k=0; k<13; k++) {
            gamma += n_m[k] * i_m[k] * Math.pow(pi, i_m[k] - 1) * j_m[k] * Math.pow(tau-0.5, j_m[k] - 1);
        }
        return gamma;
    }

    // Gibbs Free Energy g(p, t)
    public static double g(double p, double t) {
        double pi = p/p_;
        double tau = t_/t;
        return SteamCalc.R * t * (gammao(pi, tau, false) + gammar(pi, tau));
    }

    // Gibbs Free Energy for saturated vapor
    public static double gsat(double psat, double tsat) {
        double pi = psat/p_;
        double tau = t_/tsat;
        return SteamCalc.R * tsat * (gammao(pi, tau, true) + gammarm(pi, tau));
    }

    // Specific volume v(p,t)
    public static double v(double p, double t) {
        double pi = p/p_;
        double tau = t_/t;
        return (SteamCalc.R * t * 0.001/p) * pi * (gammao_pi(pi, tau, false) + gammar_pi(pi, tau));
    }

    // Specific volume v(p,t) for saturated steam
    public static double vsat(double psat, double tsat) {
        double pi = psat/p_;
        double tau = t_/tsat;
        return (SteamCalc.R * tsat * 0.001/psat) * pi * (gammao_pi(pi, tau, true) + gammarm_pi(pi, tau));
    }

    // Specific internal energy u(p,t)
    public static double u(double p, double t) {
        double pi = p/p_;
        double tau = t_/t;
        return (SteamCalc.R * t) * (tau * (gammao_tau(pi, tau, false) + gammar_tau(pi, tau))
                - pi * (gammao_pi(pi, tau, false) + gammar_pi(pi, tau)));
    }

    // Specific internal energy u(p,t) for saturated steam
    public static double usat(double psat, double tsat) {
        double pi = psat/p_;
        double tau = t_/tsat;
        return (SteamCalc.R * tsat) * (tau * (gammao_tau(pi, tau, true) + gammarm_tau(pi, tau))
                - pi * (gammao_pi(pi, tau, true) + gammarm_pi(pi, tau)));
    }

    // Specific entropy s(p,t)
    public static double s(double p, double t) {
        double pi = p/p_;
        double tau = t_/t;
        return (SteamCalc.R) * (tau * (gammao_tau(pi, tau, false) + gammar_tau(pi, tau))
                - (gammao(pi, tau, false) + gammar(pi, tau)));
    }

    // Specific entropy s(p,t) for saturated steam
    public static double ssat(double psat, double tsat) {
        double pi = psat/p_;
        double tau = t_/tsat;
        return (SteamCalc.R) * (tau * (gammao_tau(pi, tau, true) + gammarm_tau(pi, tau))
                - (gammao(pi, tau, true) + gammarm(pi, tau)));
    }

    // Specific enthalpy h(p,t)
    public static double h(double p, double t) {
        double pi = p/p_;
        double tau = t_/t;
        return (SteamCalc.R * t) * tau * (gammao_tau(pi, tau, false) + gammar_tau(pi, tau));
    }

    // Specific enthalpy h(p,t) for saturated steam
    public static double hsat(double psat, double tsat) {
        double pi = psat/p_;
        double tau = t_/tsat;
        return (SteamCalc.R * tsat) * tau * (gammao_tau(pi, tau, true) + gammarm_tau(pi, tau));
    }

    // Specific isobaric heat capacity cp(p,t)
    public static double cp(double p, double t) {
        double pi = p/p_;
        double tau = t_/t;
        return (SteamCalc.R) * -tau*tau * (gammao_tau2(pi, tau, false) + gammar_tau2(pi, tau));
    }

    // Specific isobaric heat capacity cp(p,t) for saturated steam
    public static double cpsat(double psat, double tsat) {
        double pi = psat/p_;
        double tau = t_/tsat;
        return (SteamCalc.R) * -tau*tau * (gammao_tau2(pi, tau, true) + gammarm_tau2(pi, tau));
    }

    // Specific isochoric heat capacity cv(p,t)
    public static double cv(double p, double t) {
        double pi = p/p_;
        double tau = t_/t;
        double a = 1 + pi * gammar_pi(pi, tau) - tau*pi *gammar_pi_tau(pi, tau);
        return (SteamCalc.R) * (-tau*tau * (gammao_tau2(pi, tau, false) + gammar_tau2(pi, tau))
                - a*a / (1 - pi*pi * gammar_pi2(pi, tau)));
    }

    // Specific isochoric heat capacity cv(p,t) for saturated steam
    public static double cvsat(double psat, double tsat) {
        double pi = psat/p_;
        double tau = t_/tsat;
        double a = 1 + pi * gammarm_pi(pi, tau) - tau*pi *gammarm_pi_tau(pi, tau);
        return (SteamCalc.R) * (-tau*tau * (gammao_tau2(pi, tau, true) + gammarm_tau2(pi, tau))
                - a*a / (1 - pi*pi * gammarm_pi2(pi, tau)));
    }

    // Speed of sound w(p,t)
    public static double w(double p, double t) {
        double pi = p/p_;
        double tau = t_/t;
        double a = pi * gammar_pi(pi, tau);
        double b = 1 + a - tau*pi * gammar_pi_tau(pi, tau);
        double c = (1 + a)*(1 + a) / ((1 - pi*pi * gammar_pi2(pi, tau))
                + b*b / (tau*tau * (gammao_tau2(pi, tau, false) + gammar_tau2(pi, tau))));
        return Math.sqrt(1000 * c * SteamCalc.R * t); // still don't know why, but must multiply with 1000
    }

    // Speed of sound w(p,t) for saturated steam
    public static double wsat(double psat, double tsat) {
        double pi = psat/p_;
        double tau = t_/tsat;
        double a = pi * gammarm_pi(pi, tau);
        double b = 1 + a - tau*pi * gammarm_pi_tau(pi, tau);
        double c = (1 + a)*(1 + a) / ((1 - pi*pi * gammarm_pi2(pi, tau))
                + b*b / (tau*tau * (gammao_tau2(pi, tau, true) + gammarm_tau2(pi, tau))));
        return Math.sqrt(1000*c * SteamCalc.R * tsat);
    }

    // B2bc equation for determining subregions 2b,2c
    // for backward equation T(p,h)
    static double b2bc(double h) {
        double n1 = 905.84278514723;
        double n2 = -0.67955786399241;
        double n3 = 1.2809002730136E-4;
        double n4 = 2.6526571908428E3;
        double n5 = 4.5257578905948;
        double p__ = 1;
        double h__ = 1;

        double eta = h/h__;
        double pi = n1 + n2*eta + n3*eta*eta;
        return pi * p__;
    }

    // Backward equation T(p, h)
    public static double t_ph(double p, double h) {
        if (p <= 4) return t_ph2a(p, h);
        double pp = b2bc(h);
        if (p <= pp) return t_ph2b(p, h);
        return t_ph2c(p, h);
    }

    // T(p,h) for region 2a
    static double t_ph2a(double p, double h) {
        int i_r[] = {0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,4,4,4,5,5,5,6,6,7};
        int j_r[] = {0,1,2,3,7,20,0,1,2,3,7,9,11,18,44,0,2,7,36,38,40,42,44,24,44,12,32,44,32,36,42,34,44,28};
        double n_r[] = {1.0898952318288E3,8.4951654495535E2,-1.0781748091826E2,33.153654801263,
        -7.4232016790248,11.765048724356,1.8445749355790,-4.1792700549624,6.2478196935812,
        -17.344563108114,-200.58176862096,271.96065473796,-455.11318285818,3.0919688604755E3,
        2.5226640357872E5,-6.1707422868339E-3,-0.31078046629583,11.670873077107,
        1.2812798404046E8,-9.8554909623276E8,2.8224546973002E9,-3.5948971410703E9,
        1.7227349913197E9,-1.3551334240775E4,1.2848734664650E7,1.3865724283226,
        2.3598832556514E5,-1.3105236545054E7,7.3999835474766E3,-5.5196697030060E5,
        3.7154085996233E6,1.9127729239660E4,-4.1535164835634E5,-62.459855192507};
        double p__ = 1; // MPa
        double h__ = 2000; // kJ/kg
        double t__ = 1; // K
        double theta = 0;
        double pi = p/p__;
        double eta = h/h__;

        for (int k=0; k<34; k++) {
            theta += n_r[k] * Math.pow(pi, i_r[k]) * Math.pow(eta-2.1, j_r[k]);
        }

        return theta * t__;
    }

    // T(p,h) for region 2b
    static double t_ph2b(double p, double h) {
        int i_r[] = {0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,4,4,
        5,5,5,6,7,7,9,9};
        int j_r[] = {0,1,2,12,18,24,28,40,0,2,6,12,18,24,28,40,2,8,18,40,1,2,12,
        24,2,12,18,24,28,40,18,24,40,28,2,28,1,40};
        double n_r[] = {1.4895041079516E3,743.07798314034,-97.708318797837,
        2.4742464705674,-0.63281320016026,1.1385952129658,-0.47811863648625,
        8.5208123431544E-3,0.93747147377932,3.3593118604916,3.3809355601454,
        0.16844539671904,0.73875745236695,-0.47128737436186,0.15020273139707,
        -2.1764114219750E-3,-2.1810755324761E-2,-0.10829784403677,-4.6333324635812E-2,
        7.1280351959551E-5,1.1032831789999E-4,1.8955248387902E-4,3.0891541160537E-3,
        1.3555504554949E-3,2.8640237477456E-7,-1.0779857357512E-5,-7.6462712454814E-5,
        1.4052392818316E-5,-3.1083814331434E-5,-1.0302738212103E-6,2.8217281635040E-7,
        1.2704902271945E-6,7.3803353468292E-8,-1.1030139238909E-8,-8.1456365207833E-14,
        -2.5180545682962E-11,-1.7565233969407E-18,8.6934156344163E-15};
        double p__ = 1; // MPa
        double h__ = 2000; // kJ/kg
        double t__ = 1; // K
        double theta = 0;
        double pi = p/p__;
        double eta = h/h__;

        for (int k=0; k<38; k++) {
            theta += n_r[k] * Math.pow(pi-2, i_r[k]) * Math.pow(eta-2.6, j_r[k]);
        }

        return theta * t__;
    }

    // T(p,h) for region 2c
    static double t_ph2c(double p, double h) {
        int i_r[] = {-7,-7,-6,-6,-5,-5,-2,-2,-1,-1,0,0,1,1,2,6,6,6,6,6,6,6,6};
        int j_r[] = {0,4,0,2,0,2,0,1,0,2,0,1,4,8,4,0,1,4,10,12,16,20,22};
        double n_r[] = {-3.2368398555242E12,7.3263350902181E12,3.5825089945447E11,
        -5.8340131851590E11,-1.0783068217470E10,2.0825544563171E10,6.1074783564516E5,
        8.5977722535580E5,-2.5745723604170E4,3.1081088422714E4,1.2082315865936E3,
        482.19755109255,3.7966001272486,-10.842984880077,-4.5364172676660E-2,
        1.4559115658698E-13,1.1261597407230E-12,-1.7804982240686E-11,1.2324579690832E-7,
        -1.1606921130984E-6,2.7846367088554E-5,-5.9270038474176E-4,1.2918582991878E-3};
        double p__ = 1; // MPa
        double h__ = 2000; // kJ/kg
        double t__ = 1; // K
        double theta = 0;
        double pi = p/p__;
        double eta = h/h__;

        for (int k=0; k<23; k++) {
            theta += n_r[k] * Math.pow(pi+25, i_r[k]) * Math.pow(eta-1.8, j_r[k]);
        }

        return theta * t__;
    }

    // Backward equation T(p, s)
    public static double t_ps(double p, double s) {
        if (p <= 4) return t_ps2a(p, s);
        if (s >= 5.85) return t_ps2b(p, s);
        return t_ps2c(p, s);
    }

    // T(p,s) for region 2a
    static double t_ps2a(double p, double s) {
        double i_r[] = {-1.5,-1.5,-1.5,-1.5,-1.5,-1.5,-1.25,-1.25,-1.25,-1,-1,-1,
        -1,-1,-1,-0.75,-0.75,-0.5,-0.5,-0.5,-0.5,-0.25,-0.25,-0.25,-0.25,0.25,
        0.25,0.25,0.25,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.75,0.75,0.75,0.75,1,1,1.25,
        1.25,1.5,1.5};
        int j_r[] = {-24,-23,-19,-13,-11,-10,-19,-15,-6,-26,-21,-17,-16,-9,-8,
        -15,-14,-26,-13,-9,-7,-27,-25,-11,-6,1,4,8,11,0,1,5,6,10,14,16,0,4,9,17,
        7,18,3,15,5,18};
        double n_r[] = {-3.9235983861984E5,5.1526573827270E5,4.0482443161048E4,
        -321.93790923902,96.961424218694,-22.867846371773,-4.4942914124357E5,
        -5.0118336020166E3,0.35684463560015,4.4235335848190E4,-1.3673388811708E4,
        4.2163260207864E5,2.2516925837475E4,474.42144865646,-149.31130797647,
        -1.9781126320452E5,-2.3554399470760E4,-1.9070616302076E4,5.5375669883164E4,
        3.8293691437363E3,-603.91860580567,1.9363102620331E3,4.2660643698610E3,
        -5.9780638872718E3,-704.01463926862,338.36784107553,20.862786635187,
        3.3834172656196E-2,-4.3124428414893E-5,166.53791356412,-139.86292055898,
        -0.78849547999872,7.2132411753872E-2,-5.9754839398283E-3,-1.2141358953904E-5,
        2.3227096733871E-7,-10.538463566194,2.0718925496502,-7.2193155260427E-2,
        2.0749887081120E-7,-1.8340657911379E-2,2.9036272348696E-7,0.21037527893619,
        2.5681239729999E-4,-1.2799002933781E-2,-8.2198102652018E-6};
        double p__ = 1; // MPa
        double s__ = 2; // kJ/kg K
        double t__ = 1; // K
        double pi = p/p__;
        double sigma = s/s__;
        double theta = 0;

        for (int k=0; k<46; k++) {
            theta += n_r[k] * Math.pow(pi, i_r[k]) * Math.pow(sigma-2, j_r[k]);
        }

        return theta * t__;
    }

    // T(p,s) for region 2b
    static double t_ps2b(double p, double s) {
        int i_r[] = {-6,-6,-5,-5,-4,-4,-4,-3,-3,-3,-3,-2,-2,-2,-2,-1,-1,-1,-1,-1,
        0,0,0,0,0,0,0,1,1,1,1,1,1,2,2,2,3,3,3,4,4,5,5,5};
        int j_r[] = {0,11,0,11,0,1,11,0,1,11,12,0,1,6,10,0,1,5,8,9,0,1,2,4,5,6,9,
        0,1,2,3,7,8,0,1,5,0,1,3,0,1,0,1,2};
        double n_r[] = {3.1687665083497E5,20.864175881858,-3.9859399803599E5,
        -21.816058518877,2.2369785194242E5,-2.7841703445817E3,9.9207436071480,
        -7.5197512299157E4,2.9708605951158E3,-3.4406878548526,0.38815564249115,
        1.7511295085750E4,-1.4237112854449E3,1.0943803364167,0.89971619308495,
        -3.3759740098958E3,471.62885818355,-1.9188241993679,0.41078580492196,
        -0.33465378172097,1.3870034777505E3,-406.63326195838,41.727347159610,
        2.1932549434532,-1.0320050009077,0.35882943516703,5.2511453726066E-3,
        12.838916450705,-2.8642437219381,0.56912683664855,-9.9962954584931E-2,
        -3.2632037778459E-3,2.3320922576723E-4,-0.15334809857450,2.9072288239902E-2,
        3.7534702741167E-4,1.7296691702411E-3,-3.8556050844504E-4,-3.5017712292608E-5,
        -1.4566393631492E-5,5.6420857267269E-6,4.1286150074605E-8,-2.0684671118824E-8,
        1.6409393674725E-9};
        double p__ = 1; // MPa
        double s__ = 0.7853; // kJ/kg K
        double t__ = 1; // K
        double pi = p/p__;
        double sigma = s/s__;
        double theta = 0;

        for (int k=0; k<44; k++) {
            theta += n_r[k] * Math.pow(pi, i_r[k]) * Math.pow(10-sigma, j_r[k]);
        }

        return theta * t__;
    }

    // T(p,s) for region 2c
    static double t_ps2c(double p, double s) {
        int i_r[] = {-2,-2,-1,0,0,0,0,1,1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,7,7,7,
        7,7};
        int j_r[] = {0,1,0,0,1,2,3,0,1,3,4,0,1,2,0,1,5,0,1,4,0,1,2,0,1,0,1,3,4,5};
        double n_r[] = {909.68501005365,2.4045667088420E3,-591.62326387130,
        541.45404128074,-270.98308411192,979.76525097926,-469.66772959435,
        14.399274604723,-19.104204230429,5.3299167111971,-21.252975375934,
        -0.31147334413760,0.60334840894623,-4.2764839702509E-2,5.8185597255259E-3,
        -1.4597008284753E-2,5.6631175631027E-3,-7.6155864584577E-5,2.2440342919332E-4,
        -1.2561095013413E-5,6.3323132660934E-7,-2.0541989675375E-6,3.6405370390082E-8,
        -2.9759897789215E-9,1.0136618529763E-8,5.9925719692351E-12,-2.0677870105164E-11,
        -2.0874278181886E-11,1.0162166825089E-10,-1.6429828281347E-10};
        double p__ = 1; // MPa
        double s__ = 2.9251; // kJ/kg K
        double t__ = 1; // K
        double pi = p/p__;
        double sigma = s/s__;
        double theta = 0;

        for (int k=0; k<30; k++) {
            theta += n_r[k] * Math.pow(pi, i_r[k]) * Math.pow(2-sigma, j_r[k]);
        }

        return theta * t__;
    }

    // Backward equation p(h, s)
    public static double p_hs(double h, double s) {
        if (h <= h2ab(s)) return p_hs2a(h, s);
        if (s >= 5.85) return p_hs2b(h, s);
        return p_hs2c(h, s);
    }

    // Backward equation T(h, s)
    public static double t_hs(double h, double s) {
        if (h <= h2ab(s)) return t_ph2a(p_hs2a(h, s), h);
        if (s >= 5.85) return t_ph2b(p_hs2b(h, s), h);
        return t_ph2c(p_hs2c(h, s), h);
    }

    // Boundary equation between Region 2a and 2b
    static double h2ab(double s) {
        double n1 = -3498.98083432139;
        double n2 = 2575.60716905876;
        double n3 = -421.073558227969;
        double n4 = 27.6349063799944;
        double s__ = 1; // kJ/(kg K)
        double h__ = 1; // kJ/kg
        double sigma = s/s__;
        double eta = n1 + (n2 + (n3 + n4 * sigma) * sigma) * sigma;
        return eta * h__;
    }

    // Backward equation p(h, s) for region 2a
    static double p_hs2a(double h, double s) {
        int i_[] = {0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,3,3,3,3,3,4,5,5,6,7};
        int j_[] = {1,3,6,16,20,22,0,1,2,3,5,6,10,16,20,22,3,16,20,0,2,3,6,16,16,3,16,3,1};
        double n_[] = {-1.82575361923032E-2,-0.125229548799536,0.592290437320145,
        6.04769706185122,238.624965444474,-298.639090222922,5.12250813040750E-2,
        -0.437266515606486,0.413336902999504,-5.16468254574773,-5.57014838445711,
        12.8555037824478,11.4144108953290,-119.504225652714,-2847.77985961560,
        4317.57846408006,1.12894040802650,1974.09186206319,1516.12444706087,
        1.41324451421235E-2,0.585501282219601,-2.97258075863012,5.94567314847319,
        -6236.56565798905,9659.86235133332,6.81500934948134,-6332.07286824489,
        -5.58919224465760,4.00645798472063E-2};
        double p__ = 4; // MPa
        double h__ = 4200; // kJ/kg
        double s__ = 12; // kJ/(kg K)
        double eta = h/h__;
        double sigma = s/s__;
        double pi = 0;

        for (int k=0; k<29; k++) {
            pi += n_[k] * Math.pow(eta-0.5, i_[k]) * Math.pow(sigma-1.2, j_[k]);
        }

        return pi*pi*pi*pi * p__;
    }

    // Backward equation p(h, s) for region 2b
    static double p_hs2b(double h, double s) {
        int i_[] = {0,0,0,0,0,1,1,1,1,1,1,2,2,2,3,3,3,3,4,4,5,5,6,6,6,7,7,8,8,8,
        8,12,14};
        int j_[] = {0,1,2,4,8,0,1,2,3,5,12,1,6,18,0,1,7,12,1,16,1,12,1,8,18,1,16,
        1,3,14,18,10,16};
        double n_[] = {8.01496989929495E-2,-0.543862807146111,0.337455597421283,
        8.90555451157450,313.840736431485,0.797367065977789,-1.21616973556240,
        8.72803386937477,-16.9769781757602,-186.552827328416,9.51159274344237E4,
        -18.9168510120494,-4334.07037194840,5.43212633012715E8,0.144793408386013,
        128.024559637516,-6.72309534071268E4,3.36972380095287E7,-586.634196762720,
        -2.21403224769889E10,1716.06668708389,-5.70817595806302E8,-3121.09693178482,
        -2.07841384633010E6,3.05605946157786E12,3221.57004314333,3.26810259797295E11,
        -1441.04158934487,410.694867802691,1.09077066873024E11,-2.47964654258893E13,
        1.88801906865134E9,-1.23651009018773E14};
        double p__ = 100; // MPa
        double h__ = 4100; // kJ/kg
        double s__ = 7.9; // kJ/(kg K)
        double eta = h/h__;
        double sigma = s/s__;
        double pi = 0;

        for (int k=0; k<33; k++) {
            pi += n_[k] * Math.pow(eta-0.6, i_[k]) * Math.pow(sigma-1.01, j_[k]);
        }

        return pi*pi*pi*pi * p__;
    }

    // Backward equation p(h, s) for region 2c
    static double p_hs2c(double h, double s) {
        int i_[] = {0,0,0,0,0,0,1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,5,5,5,5,6,6,10,12,16};
        int j_[] = {0,1,2,3,4,8,0,2,5,8,14,2,3,7,10,18,0,5,8,16,18,18,1,4,6,14,8,18,
        7,7,10};
        double n_[] = {0.112225607199012,-3.39005953606712,-32.0503911730094,
        -197.597305104900,-407.693861553446,1.32943775222331E4,1.70846839774007,
        37.3694198142245,3581.44365815434,4.23014446424664E5,-7.51071025760063E8,
        52.3446127607898,-228.351290812417,-9.60652417056937E5,-8.07059292526074E7,
        1.62698017225669E12,0.772465073604171,4.63929973837746E4,-1.37317885134128E7,
        1.70470392630512E12,-2.51104628187308E13,3.17748830835520E13,53.8685623675312,
        -5.53089094625169E4,-1.02861522421405E6,2.04249418756234E12,2.73918446626977E8,
        -2.63963146312685E15,-1.07890854108088E9,-2.96492620980124E10,-1.11754907323424E15};
        double p__ = 100;
        double h__ = 3500;
        double s__ = 5.9;
        double eta = h/h__;
        double sigma = s/s__;
        double pi = 0;

        for (int k=0; k<31; k++) {
            pi += n_[k] * Math.pow(eta-0.7, i_[k]) * Math.pow(sigma-1.1, j_[k]);
        }

        return pi*pi*pi*pi * p__;
    }
}
