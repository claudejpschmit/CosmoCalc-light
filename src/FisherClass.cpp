#include "FisherClass.hpp"
#include <fstream>
#include <time.h>

Fisher::Fisher(map<string, double> params, string Fl_filename)
{
    cout << "... Beginning to build FisherClass ..." << endl;
    int Pk_index, Tb_index, q_index;
    CALC = new CosmoCalc(params, &Pk_index, &Tb_index, &q_index);
    this->current_params = CALC->give_current_params();
    this->fiducial_params = CALC->give_fiducial_params();
    kmin = this->fiducial_params["kmin"];
    kmax = this->fiducial_params["kmax"];
    string model_params_keys[] = {"ombh2", "omch2", "omnuh2", "omk", "hubble"};
    for (int i = 0; i < 5; ++i) {
        string key = model_params_keys[i];
        if (current_params[key] == 0.0)
            var_params.insert(pair<string,double>(key,0.0001));
        else
            var_params.insert(pair<string,double>(key,current_params[key]/100));
    }
    
    Fl_file.open(Fl_filename);

    cout << "... Fisher built ..." << endl;
}

Fisher::~Fisher()
{
    Fl_file.close();
    delete CALC;
}

void Fisher::update_Model(map<string, double> new_params, int *Pk_index, int *Tb_index, int *q_index)
{
    this->CALC->update_q(new_params, q_index);
    this->CALC->update_Pk_interpolator_direct(new_params, Pk_index);
    this->CALC->update_G21(new_params, Tb_index);
}

mat Fisher::compute_Cl(int l, int Pk_index, int Tb_index, int q_index, vector<double> krange)
{
    mat Cl = randu<mat>(krange.size(),krange.size());
    for (unsigned int i = 0; i < krange.size(); ++i) {
        double k1 = krange[i];
        for (unsigned int j = i; j < krange.size(); ++j) {
            double k2 = krange[j];
            double res = this->CALC->Cl(l, k1, k2, this->kmin, this->kmax, Pk_index, Tb_index, q_index);
            Cl(i,j) = res;
            Cl(j,i) = res;
        }
    }

    return Cl;
}

vector<vector<double>> Fisher::Cl_derivative_matrix(int l, string param_key, int *Pk_index, int *Tb_index, int *q_index, vector<double> krange)
{
    map<string,double> working_params = fiducial_params;
    double h = this->var_params[param_key];
    //double x = working_params[param_key];
    double x = fiducial_params[param_key];

    vector<vector<double>> res, f1matrix, f2matrix, f3matrix, f4matrix;
    vector<double> row;
    working_params[param_key] = x + 2 * h;
    this->update_Model(working_params, Pk_index, Tb_index, q_index);
    for (unsigned int i = 0; i < krange.size(); ++i) {
        double k1 = krange[i];
        row.clear();
        for (unsigned int j = 0; j < krange.size(); ++j) {
            double k2 = krange[j];
            double cl = this->CALC->Cl(l, k1, k2, this->kmin, this->kmax, *Pk_index, *Tb_index, *q_index);
            row.push_back(cl);
        }
        f1matrix.push_back(row);
    }
    working_params[param_key] = x + h;

    this->update_Model(working_params, Pk_index, Tb_index, q_index);
    for (unsigned int i = 0; i < krange.size(); ++i) {
        double k1 = krange[i];
        row.clear();
        for (unsigned int j = 0; j < krange.size(); ++j) {
            double k2 = krange[j];
            row.push_back(this->CALC->Cl(l, k1, k2, this->kmin, this->kmax, *Pk_index, *Tb_index, *q_index));
        }
        f2matrix.push_back(row);
    }

    working_params[param_key] = x - h;
    this->update_Model(working_params, Pk_index, Tb_index, q_index);
    for (unsigned int i = 0; i < krange.size(); ++i) {
        double k1 = krange[i];
        row.clear();
        for (unsigned int j = 0; j < krange.size(); ++j) {
            double k2 = krange[j];
            row.push_back(this->CALC->Cl(l, k1, k2, this->kmin, this->kmax, *Pk_index, *Tb_index, *q_index));
        }
        f3matrix.push_back(row);
    }

    working_params[param_key] = x - 2 * h;
    this->update_Model(working_params, Pk_index, Tb_index, q_index);
    for (unsigned int i = 0; i < krange.size(); ++i) {
        double k1 = krange[i];
        row.clear();
        for (unsigned int j = 0; j < krange.size(); ++j) {
            double k2 = krange[j];
            row.push_back(this->CALC->Cl(l, k1, k2, this->kmin, this->kmax, *Pk_index, *Tb_index, *q_index));
        }
        f4matrix.push_back(row);
    }
    
    working_params[param_key] = x;
    this->update_Model(working_params, Pk_index, Tb_index, q_index);

    double num;
    for (unsigned int i = 0; i < krange.size(); ++i) {
        row.clear();
        for (unsigned int j = 0; j < krange.size(); ++j) {
            num = -f1matrix[i][j] + 8*f2matrix[i][j] - 8*f3matrix[i][j] +\
                  f4matrix[i][j];
            num = num / (12.0 * h);    
            row.push_back(num);
        }
        res.push_back(row);
    }

    return res;
}

double Fisher::compute_Fl(int l, string param_key1, string param_key2, int *Pk_index, int *Tb_index, int *q_index)
{
    vector<vector<double>> Cl_alpha, Cl_beta;
    //This determines the size of the Cl matrices.
    int ksteps_Cl = 4;
    vector<double> krange = give_kmodes(l, this->fiducial_params["kmax"], ksteps_Cl); 
    //for (int i = 0; i < krange.size(); i++)
    //    cout << krange[i] << " ";
    //cout << endl;

    mat Cl = randu<mat>(krange.size(),krange.size());
    mat Cl_inv = Cl;

    cout << "... derivative matrix calulation started" << endl;
    Cl_alpha = this->Cl_derivative_matrix(l, param_key1, Pk_index, Tb_index, q_index, krange);
    if (param_key1 == param_key2)
        Cl_beta = Cl_alpha;
    else
        Cl_beta = this->Cl_derivative_matrix(l, param_key2, Pk_index, Tb_index, q_index, krange);

    cout << "-> The derivative matrices are done for l = " << l << endl;
    cout << "... The Cl and Cl_inv matrices will be calculated for l = " << l << endl;

    Cl = compute_Cl(l, *Pk_index, *Tb_index, *q_index, krange);
    //cout << Cl << endl;
    Cl_inv = Cl.i();
    
    cout << "-> Cl & Cl_inv are done for l = " << l << endl;
    mat Cl_a, Cl_b;
    Cl_a = randu<mat> (krange.size(), krange.size());
    Cl_b = Cl_a;
    for (unsigned int i = 0; i < krange.size(); ++i) {
        for (unsigned int j = 0; j < krange.size(); ++j) {
            Cl_a(i,j) = Cl_alpha[i][j];
            Cl_b(i,j) = Cl_beta[i][j];

        }
    } 
   
    mat product = Cl_a * Cl_inv;
    product = product * Cl_b;
    product = product * Cl_inv;

    return 0.5 * trace(product);
}

double Fisher::F(string param_key1, string param_key2)
{
    int Pk_index, Tb_index, q_index;
    int lmax = 15;
    double sum = 0;
    // IMPORTANT! l has to start at 1 since Nl_bar has j_(l-1) in it!

    //first calculation outside so that all the dependencies are calculated.
    int l0 = 10;
    cout << "Computation of Fl starts for l = " << l0 << endl;
    double fl = this->compute_Fl(l0, param_key1, param_key2, &Pk_index, &Tb_index, &q_index);
    cout << "fl with l = " << l0 << " is: " << fl << endl;
    Fl_file << l0 << " " << fl << endl;
    sum += (2*l0 + 1) * fl;
    
    // The following line parallelizes the code
    #pragma omp parallel private(Pk_index, Tb_index, q_index) 
    {
        #pragma omp for reduction (+:sum)
        for (int l = l0+1; l <= lmax; ++l) {
            int Pk_index, Tb_index, q_index;
            cout << "Computation of Fl starts for l = " << l << endl;
            double fl = this->compute_Fl(l, param_key1, param_key2, &Pk_index, &Tb_index, &q_index);
            cout << "fl with l = " << l << " is: " << fl << endl;
            Fl_file << l << " " << fl << endl;
            sum += (2*l + 1) * fl;
        }
    }
    return sum;
}

vector<double> Fisher::give_kmodes(int l, double k_max, int steps)
{
    double k_min; 
    if (l < 100) {
        k_min = (double)l/20000.0;
        if (k_min < 0.0001)
            k_min = 0.0001;
    } else if (l < 500) {
        k_min = (double)l/15000.0;
    } else {
        k_min = (double)l/10000.0;
    }
    double stepsize = (k_max - k_min)/(double)steps;
    vector<double> range;
    for (int i = 0; i <= steps; ++i)
    {
        range.push_back(k_min + i * stepsize); 
    }
    return range;
}
