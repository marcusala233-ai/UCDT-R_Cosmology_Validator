#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <stdexcept>
#include <string>

// ==========================================
// 1. ESTRUTURAS DE DADOS
// ==========================================

struct ModelConstants {
    double Lambda_UV; 
    double f_a;       
    double G5;        
    double phi_0;     
};

struct UniverseState {
    double t;      
    double a;      
    double phi;    
    double d_phi;  
};

// ==========================================
// 2. FÍSICA: TENSOR E POTENCIAL
// ==========================================

class EnergyMomentumTensor5D {
private:
    ModelConstants constants;

public:
    EnergyMomentumTensor5D(ModelConstants cons) : constants(cons) {}

    // Potencial pNGB: Gera a Energia Escura
    double potential_pNGB(double phi) const {
        return std::pow(constants.Lambda_UV, 4) * (1.0 - std::cos(phi / constants.f_a));
    }

    // Derivada do Potencial: Gera a Força de Restauração
    double d_potential_pNGB(double phi) const {
        return (std::pow(constants.Lambda_UV, 4) / constants.f_a) * std::sin(phi / constants.f_a);
    }
};

// ==========================================
// 3. SOLVER DE FRIEDMANN (EVOLUÇÃO)
// ==========================================

class FriedmannSolver {
private:
    ModelConstants constants;
    EnergyMomentumTensor5D physics;
    const double PI = 3.14159265359;

public:
    FriedmannSolver(ModelConstants cons) : constants(cons), physics(cons) {}

    // Calcula as derivadas (velocidades de mudança)
    std::vector<double> get_derivatives(const UniverseState& state) {
        // Densidade de Matéria "Velha" (apenas para background)
        double rho_m = 0.0 * std::pow(state.a, -3.0); 
        
        // Densidade do Campo UCDT-R
        double V_phi = physics.potential_pNGB(state.phi);
        double rho_phi = 0.5 * state.d_phi * state.d_phi + V_phi;
        double rho_total = rho_m + rho_phi;
        
        // Hubble Parameter
        double H_squared = (8.0 * PI * constants.G5 / 3.0) * rho_total;
        double H = std::sqrt(std::abs(H_squared));

        // Klein-Gordon Equation
        double dV = physics.d_potential_pNGB(state.phi);
        double friction = 3.0 * H * state.d_phi;
        double phi_ddot = -friction - dV;

        return { state.a * H, state.d_phi, phi_ddot };
    }

    // Integrador Numérico RK4
    UniverseState step_rk4(UniverseState current, double dt) {
        UniverseState next = current;
        auto k1 = get_derivatives(current);
        
        UniverseState s2 = current;
        s2.t += 0.5*dt; s2.a += 0.5*dt*k1[0]; s2.phi += 0.5*dt*k1[1]; s2.d_phi += 0.5*dt*k1[2];
        auto k2 = get_derivatives(s2);

        UniverseState s3 = current;
        s3.t += 0.5*dt; s3.a += 0.5*dt*k2[0]; s3.phi += 0.5*dt*k2[1]; s3.d_phi += 0.5*dt*k2[2];
        auto k3 = get_derivatives(s3);

        UniverseState s4 = current;
        s4.t += dt; s4.a += dt*k3[0]; s4.phi += dt*k3[1]; s4.d_phi += dt*k3[2];
        auto k4 = get_derivatives(s4);

        next.t += dt;
        next.a += (dt/6.0) * (k1[0] + 2*k2[0] + 2*k3[0] + k4[0]);
        next.phi += (dt/6.0) * (k1[1] + 2*k2[1] + 2*k3[1] + k4[1]);
        next.d_phi += (dt/6.0) * (k1[2] + 2*k2[2] + 2*k3[2] + k4[2]);
        return next;
    }
    
    // Calcula Equation of State (w)
    double get_w_parameter(const UniverseState& s) {
        double V = physics.potential_pNGB(s.phi);
        double Kinetic = 0.5 * s.d_phi * s.d_phi;
        
        double Pressure = Kinetic - V;
        double Density = Kinetic + V;

        if (Density <= 1e-9) return -1.0; // Evita divisão por zero no vácuo
        return Pressure / Density;
    }
};

// ==========================================
// 4. MAIN: SIMULAÇÃO TURBO
// ==========================================

int main() {
    std::cout << "=== UCDT-R* UNIFICATION DEMO (TURBO) ===" << std::endl;
    std::cout << "Objetivo: Ver a queda do campo e o nascimento da Materia Escura\n" << std::endl;

    ModelConstants params;
    params.G5 = 0.1;         
    params.Lambda_UV = 1.5;  
    params.f_a = 0.5; // Declive acentuado para forçar a queda
    
    FriedmannSolver solver(params);

    UniverseState u;
    u.t = 0.0;
    u.a = 0.1;   
    u.phi = 2.0; // Começando alto no potencial
    u.d_phi = 0.0;

    // CONFIGURAÇÃO TURBO
    double dt = 0.05;  // Passos de tempo maiores
    int steps = 300;   // Mais passos para ver a historia completa

    std::cout << std::left << std::setw(8) << "Tempo" 
              << std::setw(12) << "Fator(a)" 
              << std::setw(12) << "Campo(phi)" 
              << std::setw(12) << "w_EOS" 
              << std::setw(15) << "  Status" << std::endl;
    std::cout << "------------------------------------------------------------" << std::endl;

    for(int i=0; i <= steps; i++) {
        double w = solver.get_w_parameter(u);
        
        std::string status = "[TRANSICAO]";
        
        // Classificação Ajustada
        if (w < -0.7) status = "[DARK ENERGY]"; 
        else if (std::abs(w) < 0.4) status = "[DARK MATTER]"; // Janela de oscilacao
        else status = "[QUEDA LIVRE]"; // Fase cinetica (w=1)

        // Imprimir apenas a cada 10 passos OU se mudou para Dark Matter
        if (i % 10 == 0 || status == "[DARK MATTER]") {
            std::cout << std::left << std::setw(8) << std::fixed << std::setprecision(2) << u.t 
                      << std::setw(12) << std::setprecision(4) << u.a 
                      << std::setw(12) << std::setprecision(4) << u.phi 
                      << std::setw(12) << std::setprecision(3) << w 
                      << status << std::endl;
        }

        u = solver.step_rk4(u, dt);
        
        // Fail-safe
        if (u.a > 1000.0) break;
    }
    
    return 0;
}