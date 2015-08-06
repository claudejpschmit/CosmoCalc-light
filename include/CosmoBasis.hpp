#pragma once

#include <string>
#include <map>
#include <cmath>
#include <boost/math/constants/constants.hpp>

using namespace std;

class CosmoBasis {
    
    public:
        
        // ------------ Functions -------------- //
        
        CosmoBasis();     
        CosmoBasis(map<string,double> params);             
        ~CosmoBasis();      
       
        void generate_params(map<string,double> params);
        map<string, double> give_current_params();
        map<string, double> give_fiducial_params();
        double sph_bessel_camb(unsigned int l, double x);
    protected:
  
        void check_params();

        double E(double z);
        double Z(double z);

        // ------------ Variables -------------- //

        /// \brief params contains the working parameters.
        map<string, double> current_params;
        /**
         * \brief fiducial_params contains the fiducial parameters, 
         *        which are not meant to be changed.
         */
        map<string, double> fiducial_params;
        
        
        /// \brief The CMB temperature today [K].
        double T_CMB;
        /// \brief The Hubble constant in standard units today [km s^-1 Mpc^-1].
        double H_0;
        /// \brief The Hubble parameter today.
        double h;
        /// \brief The relative baryon density today.
        double O_b;
        /// \brief The relative Cold Dark Matter density today.
        double O_cdm;
        /// \brief The relative photon density today.
        double O_gamma;
        /// \brief The relative relativistic neutrino density today.
        double O_nu_rel;
        /// \brief The non-relative relativistic neutrino density today.
        double O_nu;
        /// \brief The relative radiation density today.
        double O_R;
        /// \brief The relative curvature density today.
        double O_k;
        /// \brief The relative matter density today.
        double O_M;
        /// \brief The total relative density today.
        double O_tot;
        /// \brief The relative Vacuum density today.
        double O_V;
        /// \brief The hubble distance today [Mpc].
        double D_H;
        /// \brief The hubble time.
        double t_H;
        
        /// \brief Cosmological Bias factor..
        double b_bias;
        /// \brief Scales at horizon crossing [Mpc^-1].
        double k_eq;
        /// \brief T_* = hc/k_b Lambda_21cm [K].
        double T_star;

        // ------------ Constants -------------- //
        
        const double pi = boost::math::constants::pi<double>();
        /// \brief The speed of light in [m/s].
        const double c = 299792458.0;
        /// \brief The Boltzmann Constant in [m^2 kg s^-2 K^-1].
        const double k_b = 1.3806488*pow(10,-23);
        /// \brief The Baryon mass in [kg].
        const double m_b = 1.674927351*pow(10,-27);
        /// \brief The Electron mass in [kg].
        const double m_e = 9.10938291*pow(10,-31);
        /// \brief The Charge of an electron in [Coulomb].
        const double e = 1.60217657*pow(10,-19);
        /// \brief Plancks Constant in [m^2 kg s^-1].
        const double h_planck = 6.62606957*pow(10,-34);
        /// \brief The Gravitational Constant in [m^3 kg^-1 s^2].
        const double G = 6.67384*pow(10,-11);
        
        /// \brief Spontaneous decay rate of the spin flip transition [s^-1].
        const double A_10 = 2.85*pow(10,-15);
        /// \brief Variance of the matter fluctuations today 
        //         smoothed on a scale of 8 h^-1 Mpc.
        const double sigma_8 = 0.8;

        /// \brief Width of the Reionization regime.
        const double delta_z_rei = 4;
        /// \brief Center of Reionization.
        const double z_rei = 10;
        /// \brief Redshift for Recombination..
        const double z_CMB = 1100;

        /// \brief Beta factor.
        const double beta = 0.7;
        
};
