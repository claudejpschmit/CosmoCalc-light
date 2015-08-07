#pragma once

#include "CosmoBasis.hpp"
#include <string>
#include <map>
#include <cmath>
#include <vector>
#include "interpolation.h"
#include "CAMB_interface.hpp"
#include "Global21cmInterface.hpp"
#include "Helper.hpp"
#include "LevinIntegrator.hpp"

using namespace std;
using namespace alglib;

class CosmoCalc : public CosmoBasis {
    
    public:
        
        CosmoCalc(map<string, double> params, int *Pk_index, int *Tb_index, int *q_index);
        ~CosmoCalc();
        void update_G21(map<string, double> params, int *Tb_index);
        void update_Pk_interpolator_direct(map<string, double> params, int *Pk_index);

        void update_q(map<string, double> params, int *q_index);
        double q_interp(double z, int q_index);
        double r_interp(double z);
        double Hf_interp(double z);

        double Cl_new(int l, double k1, double k2, double k_low,\
                double k_high, int n_levin, int Pk_index, int Tb_index, int q_index);
        double Cl_noise(int l, double k1, double k2);
 
        double Pk_interp(double k, double z, int Pk_index);
        double Tb_interp(double z, int Tb_index);

        double Cl(int l, double k1, double k2, double k_low, double k_high,\
                int Pk_index, int Tb_index, int q_index);
        
    protected:

        
        // ------------ Variables -------------- //
        int zsteps_Ml, Pk_steps;
        int lmin_bess;
        double zmin_Ml, zmax_Ml, stepsize_Ml, prefactor_Ml, k_stepsize;
             
        spline1dinterpolant q_p_interp, H_f_interp;

        vector<Pk_interpolator> Pks, Pks_full;
        vector<Tb_interpolator> Tbs, Tbs_full;
        vector<q_interpolator> qs;

        CAMB_CALLER *CAMB;
        Global21cmInterface *G21;
        //Levin *LEVIN;
       
};
