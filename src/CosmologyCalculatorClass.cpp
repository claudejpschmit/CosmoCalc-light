#include "CosmologyCalculatorClass.hpp"
#include <string>
#include <iostream>
#include "Integrator.hpp"
#include <fstream>

CosmoCalc::CosmoCalc(map<string, double> params, int *Pk_index, int *Tb_index, int *q_index)
    :
        CosmoBasis(params)
{
    cout << "... Beginning to build CosmoCalc ..." << endl;
    this->prefactor_Ml = 2*this->b_bias*this->c/pi;
    this->zmin_Ml = this->fiducial_params["zmin"];
    this->zmax_Ml = this->fiducial_params["zmax"];
    this->zsteps_Ml = this->fiducial_params["zsteps"];
    this->stepsize_Ml = (this->zmax_Ml - this->zmin_Ml)/(double)this->zsteps_Ml;
    this->Pk_steps = this->fiducial_params["Pk_steps"];
    this->k_stepsize = this->fiducial_params["k_stepsize"];

    //generate object that is the CAMB interface.
    CAMB = new CAMB_CALLER;

    cout << "... precalculating Ml dependencies ..." << endl;
    this->update_q(fiducial_params, q_index);
    this->prefactor_Ml = 2*this->b_bias * this->c / this->pi;
    cout << "... Dependencies calculated ..." << endl;
    cout << "... Initializing Pk interpolator ..." << endl;
    this->update_Pk_interpolator_direct(this->fiducial_params, Pk_index);
    cout << "... Pks calculated ..." << endl;
    cout << "... generating 21cm interface ..." << endl;
    G21 = new Global21cmInterface();
    this->update_G21(fiducial_params, Tb_index);
    cout << "... 21cm interface built ..." << endl;
    cout << "... CosmoCalc built ..." << endl;
}

double CosmoCalc::Cl(int l, double k1, double k2, double k_low, double k_high, int Pk_index, int Tb_index, int q_index)
{
    return this->Cl_new(l,k1,k2,k_low,k_high,8, Pk_index, Tb_index, q_index);
    //return (k1+k2) * k1;
}

double CosmoCalc::Cl_noise(int l, double k1, double k2)
{
    //TODO: integrand needs to be corrected.
    auto integrand = [&](double z)
    {
        double r;
        r = r_interp(z);
        double jl = sph_bessel_camb(l,k1*r);
        double hub = Hf_interp(z)*1000.0;
        return r*r*jl/hub; 
    };

    if (k1==k2) {
        // in mK
        double Tsys = 700000;
        double fcover = 1.0;
        double lmax = 5000;
        double tau = 365.25*24*60*60;
        double prefactor = 2.0 *pi*c*c * Tsys*Tsys/(fcover*fcover * fiducial_params["df"] * lmax * lmax * tau);
        double integral = integrate_simps(integrand, this->zmin_Ml, this->zmax_Ml,\
                this->zsteps_Ml);
        return prefactor * integral * integral;
    } else {
        return 0.0;
    }
}

CosmoCalc::~CosmoCalc()
{
    delete CAMB;
    delete G21;
}

double CosmoCalc::D_C(double z)
{
    return this->D_H * this->Z(z);
}

double CosmoCalc::H(double z)
{
    return this->H_0 * this->E(z);
}

void CosmoCalc::update_q(map<string,double> params, int *q_index)
{
    bool do_calc = true;
    for (unsigned int i = 0; i < qs.size(); ++i) {
        
        if (params["ombh2"] == qs[i].ombh2 && params["omnuh2"] == qs[i].omnuh2 &&\
                params["omch2"] == qs[i].omch2 && params["omk"] == qs[i].omk &&\
                params["hubble"] == qs[i].hubble && params["T_CMB"] == qs[i].t_cmb) {

            do_calc = false;
            *q_index = i;
            break;
        }
    }

    if (do_calc) {
        q_interpolator interp;
        interp.ombh2 = params["ombh2"];
        interp.omnuh2 = params["omnuh2"];
        interp.omch2 = params["omch2"];
        interp.omk = params["omk"];
        interp.hubble = params["hubble"];
        interp.t_cmb = params["T_CMB"];
        // TODO: Do this in a way that works with parallelism....
        //UPDATE D_C to use the above parameters.

        double T_CMB2 = params["T_CMB"];
        double H_02 = params["hubble"];
        double h2 = H_02 / 100.0;
        double O_b2 = params["ombh2"] / pow(h2,2);
        double O_cdm2 = params["omch2"] / pow(h2,2);
        double O_nu2 = params["omnuh2"] / pow(h2,2);
        double O_gamma2 = pow(pi,2) * pow(T_CMB2/11605.0,4) / (15.0*8.098*pow(10,-11)*pow(h2,2));
        double O_nu_rel2 = O_gamma2 * 3.0 * 7.0/8.0 * pow(4.0/11.0, 4.0/3.0);
        double O_R2 = O_gamma2 + O_nu_rel2;
        double O_k2 = params["omk"];
        double O_M2 = O_b2 + O_cdm2 + O_nu2;
        double O_tot2 = 1.0 - O_k2;
        double O_V2 = O_tot2 - O_M2 - O_R2;
        double D_H2 = c / (1000.0 * H_02);

        real_1d_array xs, ys, hs;
        xs.setlength(this->zsteps_Ml+1);
        ys.setlength(this->zsteps_Ml+1);
        hs.setlength(this->zsteps_Ml+1);
        double z;
        for (int n = 0; n <= this->zsteps_Ml; ++n) {
            z = this->zmin_Ml + n * this->stepsize_Ml;
            xs[n] = z;
    
            auto integrand = [&](double x){
                return 1/sqrt(O_V2 + O_R2 * pow(1+z,4) + O_M2 * pow(1+z,3) + O_k2 * pow(1+z,2));
            };
            double Z = integrate(integrand, 0.0, z, 1000, simpson());

            ys[n] = D_H2 * Z;
            hs[n] = H_02 * sqrt(O_V2 + O_R2 * pow(1+z,4) + O_M2 * pow(1+z,3) + O_k2 * pow(1+z,2));
        }
        spline1dinterpolant interpolator, interpolator_Hf;
        spline1dbuildlinear(xs,ys,interpolator);
        spline1dbuildlinear(xs,hs,interpolator_Hf);
 
        interp.h = h2;
        interp.interpolator = interpolator;
        interp.interpolator_Hf = interpolator_Hf;
        
        qs.push_back(interp);
        *q_index = qs.size() - 1;
    }
}

double CosmoCalc::Hf_interp(double z)
{
    return spline1dcalc(qs[0].interpolator_Hf,z);
}

double CosmoCalc::q_interp(double z, int q_index)
{
    return spline1dcalc(qs[q_index].interpolator,z);
}

double CosmoCalc::r_interp(double z)
{
    return spline1dcalc(qs[0].interpolator,z);
}

void CosmoCalc::update_G21(map<string,double> params, int *Tb_index)
{
    bool do_calc = true;
    for (unsigned int i = 0; i < Tbs.size(); ++i) {
        if (params["ombh2"] == Tbs[i].ombh2 && params["omnuh2"] == Tbs[i].omnuh2 &&\
                params["omch2"] == Tbs[i].omch2 && params["omk"] == Tbs[i].omk &&\
                params["hubble"] == Tbs[i].hubble && params["sigma8"] == Tbs[i].s8 &&\
                params["T_CMB"] == Tbs[i].T_CMB && params["n_s"] == Tbs[i].n_s &&\
                params["fstar"] == Tbs[i].fstar && params["fesc"] == Tbs[i].fesc &&\
                params["nion"] == Tbs[i].nion && params["fx"] == Tbs[i].fx &&\
                params["flya"] == Tbs[i].flya) {

            do_calc = false;
            *Tb_index = i;
            break;
        }
    }
    if (do_calc) {
        Tb_interpolator interp;
        interp.ombh2 = params["ombh2"];
        interp.omnuh2 = params["omnuh2"];
        interp.omch2 = params["omch2"];
        interp.omk = params["omk"];
        interp.hubble = params["hubble"];
        interp.s8 = params["sigma8"];
        interp.T_CMB = params["T_CMB"];
        interp.n_s = params["n_s"];
        interp.fstar = params["fstar"];
        interp.fesc = params["fesc"];
        interp.nion = params["nion"];
        interp.fx = params["fx"];
        interp.flya = params["flya"];

        cout << "G21 is being updated" << endl;
        G21->updateGlobal21cm(params);
        vector<double> vz, vTb;
        G21->getTb(&vz, &vTb);

        real_1d_array g21_z, g21_Tb;
        g21_z.setlength(vz.size());
        g21_Tb.setlength(vTb.size());

        for (unsigned int i = 0; i < vz.size(); i++){
            g21_z[i] = vz[i];
        }
        for (unsigned int i = 0; i < vTb.size(); i++){
            g21_Tb[i] = vTb[i];
        }

        spline1dinterpolant interpolator;
        spline1dbuildcubic(g21_z, g21_Tb, interpolator);
        interp.interpolator = interpolator;

        Tbs.push_back(interp);
        *Tb_index = Tbs.size() - 1;
        
        cout << "G21 update done" << endl;
    }
}

double CosmoCalc::Tb_interp(double z, int Tb_index)
{
    return spline1dcalc(Tbs[Tb_index].interpolator,z) * 1000.0;
}

void CosmoCalc::update_Pk_interpolator_direct(map<string, double> params, int *Pk_index)
{
    bool do_calc = true;
    for (unsigned int i = 0; i < Pks.size(); ++i) {
        if (params["ombh2"] == Pks[i].ombh2 && params["omnuh2"] == Pks[i].omnuh2 &&\
                params["omch2"] == Pks[i].omch2 && params["omk"] == Pks[i].omk &&\
                params["hubble"] == Pks[i].hubble) {

            do_calc = false;
            *Pk_index = i;
            break;
        }
    }


    if (do_calc) {
        Pk_interpolator interp;
        interp.ombh2 = params["ombh2"];
        interp.omnuh2 = params["omnuh2"];
        interp.omch2 = params["omch2"];
        interp.omk = params["omk"];
        interp.hubble = params["hubble"];

        CAMB->call(params);    
        vector<double> vk = CAMB->get_k_values();
        vector<vector<double>> Pz = CAMB->get_Pz_values();

        double z_stepsize = (params["zmax"] - params["zmin"])/(params["Pk_steps"] - 1);
        vector<double> vz, vP;
        for (unsigned int i = 0; i < Pz.size(); ++i) {
            vz.push_back(params["zmin"] + i * z_stepsize);
            vP.insert(vP.end(), Pz[i].begin(), Pz[i].end());
        }

        real_1d_array matterpowerspectrum_k, matterpowerspectrum_z, matterpowerspectrum_P;
        matterpowerspectrum_k.setlength(vk.size());
        matterpowerspectrum_z.setlength(vz.size());
        matterpowerspectrum_P.setlength(vP.size());
        for (unsigned int i = 0; i < vk.size(); i++){
            matterpowerspectrum_k[i] = vk[i];
        }
        for (unsigned int i = 0; i < vP.size(); i++){
            matterpowerspectrum_P[i] = vP[i];
        }
        for (unsigned int i = 0; i < vz.size(); i++){
            matterpowerspectrum_z[i] = vz[i];
        }

        spline2dinterpolant interpolator;
        spline2dbuildbilinearv(matterpowerspectrum_k, vk.size(),matterpowerspectrum_z, vz.size(),\
                matterpowerspectrum_P, 1, interpolator);
        interp.interpolator = interpolator;

        Pks.push_back(interp);
        *Pk_index = Pks.size() - 1;
    }
}
double CosmoCalc::Pk_interp(double k, double z, int Pk_index)
{
    return spline2dcalc(Pks[Pk_index].interpolator, k, z);
}

double CosmoCalc::Cl_new(int l, double k1, double k2, double k_low,\
        double k_high, int n_levin, int Pk_index, int Tb_index, int q_index)
{
    double a;
    double low;
    double hhh = pow(qs[q_index].h,3);

    auto integrand1 = [&](double z)
    {
        double r,q;
        r = r_interp(z);
        q = q_interp(z, q_index);

        if (l < 1000){
            low = (double)l/(1.2*q);
            a = (double)(l+1000)/(1.5*q);
        } else {
            low = (double)l/(q);
            a = (double)(l+1000)/(1.5*q);
        }
        double lower_kappa_bound;
        if (low > k_low)
            lower_kappa_bound = low;
        else
            lower_kappa_bound = k_low;

        int steps = (int)((a - lower_kappa_bound)/this->k_stepsize);
        if (steps % 2 == 1)
            ++steps;

        auto integrand2 = [&](double zp)
        {
            double rp,qp;
            rp = r_interp(zp);
            qp = q_interp(zp, q_index);

            auto integrand3 = [&](double kappa)
            {
                double sP = sqrt(this->Pk_interp(kappa*qs[q_index].h,z, Pk_index)/hhh);
                double sPp = sqrt(this->Pk_interp(kappa*qs[q_index].h,zp, Pk_index)/hhh);
                return kappa*kappa * sP * sPp * this->sph_bessel_camb(l,kappa*q) *\
                    this->sph_bessel_camb(l,kappa*qp);
            };
            double integral = integrate_simps(integrand3, lower_kappa_bound, a, steps);
            double integral3;
            
            LEVIN = new Levin(a, k_high);

            auto foo = [&](double kappa)
            {
                double sP1 = sqrt(this->Pk_interp(kappa*qs[q_index].h,z, Pk_index)/hhh);
                double sPp1 = sqrt(this->Pk_interp(kappa*qs[q_index].h,zp, Pk_index)/hhh);
                return kappa*kappa * sP1 * sPp1;     
            };

            
            if (z == zp)
                integral3 = LEVIN->integrate_2sphj_1r(foo,q,l,n_levin);
            else
                integral3 = LEVIN->integrate_2sphj_2r(foo,q,qp,l,n_levin);

            delete LEVIN;
            integral3 += integral;
            return rp*rp / (Hf_interp(zp)*1000.0) * this->Tb_interp(zp, Tb_index) *\
                this->sph_bessel_camb(l,k2*rp) * integral3;
        };
        double integral2 = integrate_simps(integrand2, this->zmin_Ml, this->zmax_Ml, this->zsteps_Ml);
        return r*r / (Hf_interp(z)*1000.0) * this->Tb_interp(z, Tb_index) *\
            this->sph_bessel_camb(l,k1*r) * integral2;
    };
    double integral1 = integrate_simps(integrand1,this->zmin_Ml, this->zmax_Ml, this->zsteps_Ml);
    return pow(this->prefactor_Ml,2) * integral1;
}
