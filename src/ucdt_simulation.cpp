#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include <iomanip>

// =========================================================
// PROJETO: UCDT-R* COSMOLOGY VALIDATOR
// ARQUIVO: ucdt_simulation.cpp
// DESCRIÇÃO: Motor MCMC para validar a massa de 1.65e-23 eV
// =========================================================

using namespace std;

// --- CONSTANTES ---
const double EV_TO_HZ = 2.417989e14; 
const double PI = 3.14159265358979323846;
const double TARGET_MASS_EV = 1.65e-23; [cite_start]// Alvo da UCDT-R* [cite: 13]

struct ModelParams {
    double log_Amp_SMBH;
    double log_Mass_ULDM;
    double log_Coupling;
};

// --- MODELO FÍSICO ---
vector<double> calculate_spectrum(const ModelParams& theta, const vector<double>& freqs) {
    vector<double> spectrum(freqs.size());
    double A = pow(10, theta.log_Amp_SMBH);
    double mass = pow(10, theta.log_Mass_ULDM);
    double coupling = pow(10, theta.log_Coupling);
    
    double f_particle = mass * EV_TO_HZ;
    double gamma = 0.1 * f_particle; 
    
    for (size_t i = 0; i < freqs.size(); ++i) {
        double f = freqs[i];
        // Fundo de Buracos Negros (SMBH)
        double omega_smbh = (A*A / (12 * PI*PI)) * pow(f / 3.17e-8, -4.0/3.0); 
        // Sinal da UCDT-R* (Monopolo)
        double denom = (f - f_particle)*(f - f_particle) + gamma*gamma;
        double omega_ucdt = coupling * (gamma*gamma) / denom;
        spectrum[i] = omega_smbh + omega_ucdt;
    }
    return spectrum;
}

// --- LIKELIHOOD ---
double get_log_likelihood(const ModelParams& theta, const vector<double>& f, 
                          const vector<double>& y_data, const vector<double>& y_err) {
    // SNMC Stability Check (Filtro de Fantasmas)
    if (theta.log_Mass_ULDM < -24.0 || theta.log_Mass_ULDM > -21.0) return -1e30; 
    
    vector<double> model = calculate_spectrum(theta, f);
    double log_L = 0.0;
    for (size_t i = 0; i < f.size(); ++i) {
        double diff = y_data[i] - model[i];
        log_L += -0.5 * (diff*diff) / (y_err[i]*y_err[i]);
    }
    return log_L;
}

int main() {
    cout << "INICIANDO SIMULACAO UCDT-R* (MCMC)..." << endl;
    mt19937 rng(1234); 
    normal_distribution<double> gaussian(0.0, 1.0);
    uniform_real_distribution<double> uniform(0.0, 1.0);

    // Gerar Dados Sintéticos (NANOGrav Proxy)
    vector<double> freqs(30), y_obs(30), y_err(30);
    ModelParams truth = {-14.7, log10(TARGET_MASS_EV), -17.5};
    vector<double> true_signal = calculate_spectrum(truth, freqs);

    for(int i=0; i<30; ++i) {
        freqs[i] = 1e-9 * pow(100.0, double(i)/29.0); 
        y_obs[i] = true_signal[i] + 1e-20 * gaussian(rng); 
        y_err[i] = 0.2 * fabs(y_obs[i]); 
    }

    // MCMC Loop
    ModelParams current = {-15.0, -22.5, -18.0}; 
    double current_logL = get_log_likelihood(current, freqs, y_obs, y_err);
    int steps = 50000, accepted = 0;
    
    // Arquivo de saída para análise posterior
    ofstream chain_file("results/mcmc_chain.txt");
    chain_file << "step,log_mass,log_coupling\n";

    for(int i=0; i<steps; ++i) {
        ModelParams proposal = current;
        proposal.log_Amp_SMBH += 0.05 * gaussian(rng);
        proposal.log_Mass_ULDM += 0.05 * gaussian(rng);
        proposal.log_Coupling += 0.05 * gaussian(rng);
        
        double proposal_logL = get_log_likelihood(proposal, freqs, y_obs, y_err);
        if (log(uniform(rng)) < (proposal_logL - current_logL)) {
            current = proposal;
            current_logL = proposal_logL;
            accepted++;
        }
        if (i % 100 == 0) chain_file << i << "," << current.log_Mass_ULDM << "," << current.log_Coupling << "\n";
    }
    chain_file.close();

    cout << "SIMULACAO CONCLUIDA." << endl;
    cout << "Massa Recuperada: " << pow(10, current.log_Mass_ULDM) << " eV" << endl;
    
    // Comando para Windows não fechar a janela
    system("pause");
    return 0;
}