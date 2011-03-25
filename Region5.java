/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package com.spayc;

/**
 *
 * @author Kehinde
 */
class Region5 {
    // Shared constants for Region5: High Temperature region
    final static double p_ = 1; // MPa
    final static double t_ = 1000; // K
    // Ideal gas component
    final static int jo[] = {0,1,-3,-2,-1,2};
    final static double no[] = {-13.179983674201,6.8540841634434,
    -2.4805148933466E-2,0.36901534980333,-3.1161318213925,-0.32961626538917};
    // Residual component
    final static int ir[] = {1,1,1,2,2,3};
    final static int jr[] = {1,2,3,3,9,7};
    final static double nr[] = {1.5736404855259E-3,9.0153761673944E-4,
    -5.0270077677648E-3,2.2440037409485E-6,-4.1163275453471E-6,3.7919454822955E-8};

    /** Gibbs Free Energy G(p,T)
     *
     * @param p
     * @param t
     * @return
     */
    public static double g(double p, double t) {
        double pi = p/p_;
        double tau = t_/t;

        return (gammao(pi, tau) + gammar(pi, tau)) * SteamCalc.R * t;
    }

    /** Specific volume v(p, T)
     *
     * @param p
     * @param t
     * @return
     */
    public static double v(double p, double t) {
        double pi = p/p_;
        double tau = t_/t;

        return pi * (gammao_pi(pi, tau) + gammar_pi(pi, tau)) * SteamCalc.R * t/(p*1000);
    }

    /** Specific internal energy u(p,T)
     *
     * @param p
     * @param t
     * @return
     */
    public static double u(double p, double t) {
        double pi = p/p_;
        double tau = t_/t;

        return SteamCalc.R * t * (tau * (gammao_tau(pi,tau) + gammar_tau(pi,tau))
                - pi * (gammao_pi(pi,tau) + gammar_pi(pi,tau)));
    }

    /** Specific entropy s(p,T)
     *
     * @param p
     * @param t
     * @return
     */
    public static double s(double p, double t) {
        double pi = p/p_;
        double tau = t_/t;

        return SteamCalc.R * (tau * (gammao_tau(pi, tau) + gammar_tau(pi, tau))
                - (gammao(pi, tau) + gammar(pi, tau)));
    }

    /** Specific enthalpy h(p, T)
     *
     * @param p
     * @param t
     * @return
     */
    public static double h(double p, double t) {
        double pi = p/p_;
        double tau = t_/t;

        return SteamCalc.R * t * tau * (gammao_tau(pi, tau) + gammar_tau(pi, tau));
    }

    /** Specific isobaric heat capacity Cp(p,T)
     * 
     * @param p
     * @param t
     * @return
     */
    public static double cp(double p, double t) {
        double pi = p/p_;
        double tau = t_/t;

        return SteamCalc.R * -tau*tau * (gammao_tau2(pi, tau) + gammar_tau2(pi, tau));
    }

    /** Specific isochoric heat capacity Cv(p,T)
     * 
     * @param p
     * @param t
     * @return
     */
    public static double cv(double p, double t) {
        double pi = p/p_;
        double tau = t_/t;
        double a = 1 + pi * gammar_pi(pi, tau) - tau * pi * gammar_pi_tau(pi, tau);

        return cp(p, t) - SteamCalc.R * a*a/(1 - pi*pi * gammar_pi2(pi, tau));
    }

    /** Speed of sound w(p, T)
     *
     * @param p
     * @param t
     * @return
     */
    public static double w(double p, double t) {
        double pi = p/p_;
        double tau = t_/t;
        double a = pi * gammar_pi(pi, tau);
        double b = 1 + a - tau * pi * gammar_pi_tau(pi, tau);

        return Math.sqrt(1000* SteamCalc.R * t * (1 + 2*a + a*a)/((1 - pi*pi*gammar_pi2(pi, tau))
                + b*b/(tau*tau*(gammao_tau2(pi, tau) + gammar_tau2(pi, tau)))));
    }

    static double gammao(double pi, double tau) {
        double gammao = Math.log(pi);

        for (int k=0; k<6; k++) {
            gammao += no[k] * Math.pow(tau, jo[k]);
        }

        return gammao;
    }

    static double gammao_pi(double pi, double tau) {
        return 1/pi;
    }

    static double gammao_pi2(double pi, double tau) {
        return -1/(pi*pi);
    }

    static double gammao_tau(double pi, double tau) {
        double gammao_tau = 0;

        for (int k=0; k<6; k++) {
            gammao_tau += no[k] * jo[k] * Math.pow(tau, jo[k]-1);
        }

        return gammao_tau;
    }

    static double gammao_tau2(double pi, double tau) {
        double gammao_tau2 = 0;

        for (int k=0; k<6; k++) {
            gammao_tau2 += no[k] * jo[k] * (jo[k]-1) * Math.pow(tau, jo[k]-2);
        }

        return gammao_tau2;
    }

    static double gammao_pi_tau(double pi, double tau) { return 0; }
    
    static double gammar(double pi, double tau) {
        double gammar = 0;

        for (int k=0; k<6; k++) {
            gammar += nr[k] * Math.pow(pi, ir[k]) * Math.pow(tau, jr[k]);
        }

        return gammar;
    }

    static double gammar_pi(double pi, double tau) {
        double gammar = 0;

        for (int k=0; k<6; k++) {
            gammar += nr[k] * ir[k] * Math.pow(pi, ir[k]-1) * Math.pow(tau, jr[k]);
        }

        return gammar;
    }

    static double gammar_pi2(double pi, double tau) {
        double gammar = 0;

        for (int k=0; k<6; k++) {
            gammar += nr[k] * ir[k] * (ir[k]-1) * Math.pow(pi, ir[k]-2) * Math.pow(tau, jr[k]);
        }

        return gammar;
    }

    static double gammar_tau(double pi, double tau) {
        double gammar = 0;

        for (int k=0; k<6; k++) {
            gammar += nr[k] * jr[k] * Math.pow(pi, ir[k]) * Math.pow(tau, jr[k]-1);
        }

        return gammar;
    }

    static double gammar_tau2(double pi, double tau) {
        double gammar = 0;

        for (int k=0; k<6; k++) {
            gammar += nr[k] * jr[k] * (jr[k]-1) * Math.pow(pi, ir[k]) * Math.pow(tau, jr[k]-2);
        }

        return gammar;
    }

    static double gammar_pi_tau(double pi, double tau) {
        double gammar = 0;

        for (int k=0; k<6; k++) {
            gammar += nr[k] * ir[k] * jr[k] * Math.pow(pi, ir[k]-1) * Math.pow(tau, jr[k]-1);
        }

        return gammar;
    }

}
