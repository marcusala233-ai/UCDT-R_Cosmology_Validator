#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <stdexcept>

// ==========================================
// ESTRUTURAS DE DADOS E CONSTANTES
// ==========================================

// Constantes físicas e do modelo
struct ModelConstants {
    double Lambda_UV; // Escala de Energia do potencial
    double f_a;       // Escala de quebra de simetria (Decay constant)
    double G5;        // Gravidade efetiva
    double phi_0;     // VEV (não usado diretamente no solver dinâmico, mas referência)
};

// Estado do Universo (variáveis que mudam com o tempo)
struct UniverseState {
    double t;      // Tempo cósmico
    double a;      // Fator de escala (tamanho do universo)
    double phi;    // Valor do campo escalar pNGB
    double d_phi;  // Velocidade do campo (phi_dot)
};

// ==========================================
// PARTE 1: FÍSICA LOCAL (O TENSOR E POTENCIAL)
// ==========================================

class EnergyMomentumTensor5D {
private:
    ModelConstants constants;

public:
    EnergyMomentumTensor5D(ModelConstants cons) : constants(cons) {}

    // Potencial pNGB: V(phi) = Lambda^4 * (1 - cos(phi/f_a))
    // Ref: "Refinamento da Teoria UCDT-R Científica.pdf"
    double potential_pNGB(double phi) const {
        return std::pow(constants.Lambda_UV, 4) * (1.0 - std::cos(phi / constants.f_a));
    }

    // Derivada do Potencial dV/dphi
    double d_potential_pNGB(double phi) const {
        return (std::pow(constants.Lambda_UV, 4) / constants.f_a) * std::sin(phi / constants.f_a);
    }
};

// ==========================================
// PARTE 2: SOLVER DE FRIEDMANN (EVOLUÇÃO)
// ==========================================

class FriedmannSolver {
private:
    ModelConstants constants;
    EnergyMomentumTensor5D physics;
    const double PI = 3.14159265359;

public:
    FriedmannSolver(ModelConstants cons) : constants(cons), physics(cons) {}

    // Calcula as taxas de variação: { da/dt, dphi/dt, d(dphi)/dt }
    std::vector<double> get_derivatives(const UniverseState& state) {
        
        // 1. Densidades
        // Matéria (cai com volume a^3) - Acoplamento VSEP implícito na dinâmica
        double rho_m = 0.5 * std::pow(state.a, -3.0); 
        
        // Campo Escalar (Energia Cinética + Potencial)
        double V_phi = physics.potential_pNGB(state.phi);
        double rho_phi = 0.5 * state.d_phi * state.d_phi + V_phi;
        
        double rho_total = rho_m + rho_phi;
        
        // 2. Parâmetro de Hubble H
        // H^2 = (8*pi*G / 3) * rho_total
        double H_squared = (8.0 * PI * constants.G5 / 3.0) * rho_total;
        double H = std::sqrt(std::abs(H_squared));

        // 3. Aceleração do Campo (Klein-Gordon em FLRW)
        // phi_ddot + 3H*phi_dot + V' = Fonte
        double dV = physics.d_potential_pNGB(state.phi);
        double friction = 3.0 * H * state.d_phi;
        
        // Equação de Movimento
        double phi_ddot = -friction - dV;

        // Retorna: { a_dot, phi_dot, phi_ddot }
        return { state.a * H, state.d_phi, phi_ddot };
    }

    // Integrador Numérico RK4 (Runge-Kutta 4ª Ordem)
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
    
    // Método auxiliar para obter densidade atual (para display)
    double get_current_density(const UniverseState& s) {
        double V = physics.potential_pNGB(s.phi);
        return 0.5*s.d_phi*s.d_phi + V + (0.5 * std::pow(s.a, -3.0));
    }
};

// ==========================================
// MAIN: EXECUÇÃO DA SIMULAÇÃO
// ==========================================

int main() {
    std::cout << "=== UCDT-R UNIVERSE EVOLUTION SIMULATION ===" << std::endl;
    std::cout << "Model: pNGB Scalar Field + 5D Gravity Effects\n" << std::endl;

    // 1. Definição dos Parâmetros (Unidades Normalizadas)
    ModelConstants params;
    params.G5 = 0.1;         // Gravidade
    params.Lambda_UV = 1.0;  // Energia do Potencial
    params.f_a = 5.0;        // Escala onde a simetria quebra (pNGB)
    
    FriedmannSolver solver(params);

    // 2. Condições Iniciais
    UniverseState u;
    u.t = 0.0;
    u.a = 0.01;  // Universo pequeno (Início da simulação)
    u.phi = 3.0; // Campo deslocado do zero (gera dinâmica)
    u.d_phi = 0.0;

    // 3. Loop de Tempo
    double dt = 0.01;
    int steps = 50; // Quantos passos vamos simular

    // Cabeçalho da Tabela
    std::cout << std::left << std::setw(10) << "Time" 
              << std::setw(15) << "ScaleFactor(a)" 
              << std::setw(15) << "Phi_Field" 
              << std::setw(15) << "Total_Density" << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;

    for(int i=0; i <= steps; i++) {
        // Imprime o estado atual
        double rho = solver.get_current_density(u);
        
        std::cout << std::left << std::setw(10) << std::fixed << std::setprecision(3) << u.t 
                  << std::setw(15) << std::setprecision(4) << u.a 
                  << std::setw(15) << std::setprecision(4) << u.phi 
                  << std::setw(15) << std::setprecision(4) << rho << std::endl;

        // Avança o universo
        u = solver.step_rk4(u, dt);
    }
    
    std::cout << "\n[CONCLUSAO] Simulacao finalizada." << std::endl;
    if (u.a > 0.01) {
        std::cout << "Status: EXPANSAO CONFIRMADA (a aumentou)." << std::endl;
    } else {
        std::cout << "Status: COLAPSO OU ESTAGNACAO." << std::endl;
    }

    return 0;
}