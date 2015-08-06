#pragma once

#include "CosmologyCalculatorClass.hpp"
#include <armadillo>
#include <fstream>
#include <string>

using namespace arma;

class Fisher {
    public:
        Fisher(map<string, double> params, string Fl_filename);
        ~Fisher();

        void update_Model(map<string, double> new_params, int *Pk_index, int *Tb_index, int *q_index);
        mat compute_Cl(int l, int Pk_index, int Tb_index, int q_index, vector<double> krange);
        vector<vector<double>> Cl_derivative_matrix(int l, string param_key, int *Pk_index,\
                int *Tb_index, int *q_index, map<string,double> params, vector<double> krange);
        double compute_Fl(int l, string param_key1, string param_key2,\
                int *Pk_index, int *Tb_index, int *q_index, map<string,double> params);
        double F(string param_key1, string param_key2);

    private:
        vector<double> give_kmodes(int l, double kmax, int steps);
        
        ofstream Fl_file;
        CosmoCalc *CALC;
        map<string, double> current_params, fiducial_params, var_params;
        double kmin, kmax;
};
