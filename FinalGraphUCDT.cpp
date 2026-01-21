#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <stdexcept>
#include <string>

// ==========================================
// 1. ESTRUTURAS FÍSICAS (Idêntico ao V5)
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
// 2. FÍSICA: UCDT-R* (Idêntico ao V5)
// ==========================================

class EnergyMomentumTensor5D {
private:
    ModelConstants constants;

public:
    EnergyMomentumTensor5D(ModelConstants cons) : constants(cons) {}

    double potential_pNGB(double phi) const {
        return std::pow(constants.Lambda_UV, 4) * (1.0 - std::cos(phi / constants.f_a));
    }

    double d_potential_pNGB(double phi) const {
        return (std::pow(constants.Lambda_UV, 4) / constants.f_a) * std::sin(phi / constants.f_a);
    }
};

class FriedmannSolver {
private:
    ModelConstants constants;
    EnergyMomentumTensor5D physics;
    const double PI = 3.14159265359;

public:
    FriedmannSolver(ModelConstants cons) : constants(cons), physics(cons) {}

    std::vector<double> get_derivatives(const UniverseState& state) {
        double rho_m = 0.0; 
        double V_phi = physics.potential_pNGB(state.phi);
        double rho_phi = 0.5 * state.d_phi * state.d_phi + V_phi;
        double rho_total = rho_m + rho_phi;
        
        double H_squared = (8.0 * PI * constants.G5 / 3.0) * rho_total;
        double H = std::sqrt(std::abs(H_squared));

        double dV = physics.d_potential_pNGB(state.phi);
        double friction = 3.0 * H * state.d_phi;
        double phi_ddot = -friction - dV;

        return { state.a * H, state.d_phi, phi_ddot };
    }

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
    
    double get_w_parameter(const UniverseState& s) {
        double V = physics.potential_pNGB(s.phi);
        double Kinetic = 0.5 * s.d_phi * s.d_phi;
        double Pressure = Kinetic - V;
        double Density = Kinetic + V;
        if (Density <= 1e-9) return -1.0;
        return Pressure / Density;
    }
};

// ==========================================
// 3. MAIN: VISUALIZADOR GRÁFICO (ASCII ART)
// ==========================================

void draw_bar(double w) {
    // Mapeia w de [-1.1 a +1.1] para uma barra de 50 caracteres
    // w = -1 (Energia Escura) -> Esquerda
    // w = 0  (Matéria Escura) -> Meio
    // w = +1 (Cinética)       -> Direita
    
    int width = 40;
    int position = (int)((w + 1.1) / 2.2 * width); 
    if (position < 0) position = 0;
    if (position >= width) position = width - 1;

    std::cout << " [";
    for (int i = 0; i < width; i++) {
        if (i == position) std::cout << "O"; // A "partícula" se movendo
        else if (i == width/2) std::cout << "|"; // A marca de Matéria Escura (0)
        else std::cout << "-";
    }
    std::cout << "] ";
}

int main() {
    std::cout << "\n=== UCDT-R* VISUALIZADOR EM TEMPO REAL ===" << std::endl;
    std::cout << "Legenda do Grafico:" << std::endl;
    std::cout << "Esquerda [-1] = Energia Escura (Inflacao)" << std::endl;
    std::cout << "Meio     [ 0] = Materia Escura (Alvo)" << std::endl;
    std::cout << "Direita  [+1] = Energia Cinetica Pura\n" << std::endl;

    ModelConstants params;
    params.G5 = 0.1; params.Lambda_UV = 1.5; params.f_a = 0.5; // Turbo Settings
    FriedmannSolver solver(params);
    
    UniverseState u;
    u.t = 0.0; u.a = 0.1; u.phi = 2.0; u.d_phi = 0.0;

    double dt = 0.05;
    int steps = 300;

    for(int i=0; i <= steps; i++) {
        double w = solver.get_w_parameter(u);
        
        // Imprime a cada 5 passos para criar animação fluida
        if (i % 5 == 0) {
            std::cout << "t=" << std::fixed << std::setw(5) << std::setprecision(2) << u.t << " ";
            draw_bar(w); 
            
            if (w < -0.7) std::cout << "DE (Inflacao)";
            else if (std::abs(w) < 0.4) std::cout << "** DM (Materia) **";
            else std::cout << "Oscilacao";
            
            std::cout << std::endl;
        }

        u = solver.step_rk4(u, dt);
        if (u.a > 1000.0) break;
    }
    
    return 0;
}