#include "CAMB_interface.hpp"
#include <fstream>
#include <iostream>
#include <vector>
#include <sstream>
#include <time.h>

CAMB_CALLER::CAMB_CALLER()
    :
        run_first_time(true)
{
    ifstream params_ini_file;
    string line;
    params_ini_file.open("CAMB/params.ini");

    while (params_ini_file.good()) {
        getline(params_ini_file, line);
        file_content.push_back(line);
    }

    parameter_names[0] = "ombh2";
    parameter_names[1] = "omch2";
    parameter_names[2] = "omnuh2";
    parameter_names[3] = "omk";
    parameter_names[4] = "hubble";
    parameter_names[5] = "transfer_num_redshifts";
    parameter_names[6] = "transfer_redshift(1)";
    parameter_names[7] = "transfer_matterpower(1)";


}

CAMB_CALLER::~CAMB_CALLER()
{
}
//TODO: !!!!!!! build in something that checks whether a certain parameter
//set has already been computed and then skips the CAMB call and returns the 
//relevant Pks, or even better, depending on how many times CAMB actually needs 
//to run, I could store a list of Pk interpolators in memory and then return 
//a position in that list where the right interpolator can be found... This 
//would need CosmoCalc to store the Pk interpolators with reference to the
//corresponding parameter values that were used to calculate it.
void CAMB_CALLER::call(map<string, double> params)
{
    cout << "CAMB is called with an ombh2 value of " << params["ombh2"] << endl;
    update_params_ini(params);

    //call camb with new_params.ini
    system("./CAMB/camb CAMB/new_params.ini");

    //recovering Power spectrum.
    read_matterpower_files(params["Pk_steps"]);    
}

void CAMB_CALLER::call_full(map<string, double> params)
{
    cout << "CAMB is called with an ombh2 value of " << params["ombh2"] << endl;
    update_params_ini_full(params);

    //call camb with new_params.ini
    system("./CAMB/camb CAMB/new_params.ini");

    //recovering Power spectrum.
    read_matterpower_files(params["zmax_interp"]+1);    
}


void CAMB_CALLER::update_params_ini_full(map<string, double> params)
{
    int n_redshifts = params["zmax_interp"] + 1;
    double zmin = 0;
    double zmax = params["zmax_interp"];
    double stepsize_z = (zmax - zmin)/(double)(n_redshifts - 1);

    int found_n_params = 0;
    for (int i = 0; i < file_content.size(); ++i) {
        if (file_content[i].find("output_root =") != string::npos) {
            file_content[i] = "output_root = CAMB/test";
        }
        for (int j = 0; j < 8; j++) {
            //pos[j] = file_content[i].find(parameter_names[j]);
            if (file_content[i].find(parameter_names[j]) != string::npos) {
                found_n_params += 1;
                string pn = parameter_names[j];
                stringstream val;
                if (pn == "transfer_num_redshifts")
                    val << n_redshifts;
                else if (pn == "transfer_redshift(1)")
                    val << zmax;
                else if (pn =="transfer_matterpower(1)") 
                    val << "matterpower_1.dat";
                else
                    val << params[pn];

                string new_parameter = parameter_names[j] + " = " + val.str();
                file_content[i] = new_parameter;

                break;
            }
        }
        if (found_n_params == 8){
            break;
        }
    }
    if (run_first_time) {
        for (int i = 1; i < n_redshifts; ++i){
            stringstream line1, line2;
            line1 << "transfer_redshift(" << i+1 << ") = " << (zmax - i * stepsize_z);
            line2 << "transfer_matterpower(" << i+1 << ") = matterpower_" << i+1 << ".dat";
            file_content.push_back(line1.str());
            file_content.push_back(line2.str());
        }
        run_first_time = false;
    }
    create_output_file();
}


void CAMB_CALLER::update_params_ini(map<string, double> params)
{
    int n_redshifts = params["Pk_steps"];
    double zmin = params["zmin"];
    double zmax = params["zmax"];
    double stepsize_z = (zmax - zmin)/(double)(n_redshifts - 1);

    int found_n_params = 0;
    for (int i = 0; i < file_content.size(); ++i) {
        if (file_content[i].find("output_root =") != string::npos) {
            file_content[i] = "output_root = CAMB/test";
        }
        for (int j = 0; j < 8; j++) {
            //pos[j] = file_content[i].find(parameter_names[j]);
            if (file_content[i].find(parameter_names[j]) != string::npos) {
                found_n_params += 1;
                string pn = parameter_names[j];
                stringstream val;
                if (pn == "transfer_num_redshifts")
                    val << n_redshifts;
                else if (pn == "transfer_redshift(1)")
                    val << zmax;
                else if (pn =="transfer_matterpower(1)") 
                    val << "matterpower_1.dat";
                else
                    val << params[pn];

                string new_parameter = parameter_names[j] + " = " + val.str();
                file_content[i] = new_parameter;

                break;
            }
        }
        if (found_n_params == 8){
            break;
        }
    }
    if (run_first_time) {
        for (int i = 1; i < n_redshifts; ++i){
            stringstream line1, line2;
            line1 << "transfer_redshift(" << i+1 << ") = " << (zmax - i * stepsize_z);
            line2 << "transfer_matterpower(" << i+1 << ") = matterpower_" << i+1 << ".dat";
            file_content.push_back(line1.str());
            file_content.push_back(line2.str());
        }
        run_first_time = false;
    }
    create_output_file();
}

void CAMB_CALLER::create_output_file()
{
    ofstream output;
    output.open("CAMB/new_params.ini");
    for (int i = 0; i < file_content.size(); i++) {
        output << file_content[i] << endl;
    }
}

void CAMB_CALLER::read_matterpower_files(int nfiles)
{
    ifstream file;
    int file_index = nfiles;
    Pz_values.clear();
    k_values.clear();
    while (file_index >= 1){
        string filename;
        filename = "CAMB/test_matterpower_" + to_string(file_index) + ".dat";
        file.open(filename);
        double k, P;
        vector<double> P_values;
        if (file_index == 1) {
            while(file >> k >> P) {
                k_values.push_back(k);
                P_values.push_back(P);
            }
        } else {
            while(file >> k >> P) {
                P_values.push_back(P);
            }
        }

        Pz_values.push_back(P_values);
        file.close();
        --file_index;
    }
}

vector<vector<double>> CAMB_CALLER::get_Pz_values()
{
    return Pz_values;
}
vector<double> CAMB_CALLER::get_k_values()
{
    return k_values;
}
