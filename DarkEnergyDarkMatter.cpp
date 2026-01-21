#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <stdexcept>

// ==========================================
// ESTRUTURAS E CONSTANTES
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
// FÍSICA: TENSOR E POTENCIAL
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

// ==========================================
// SOLVER DE FRIEDMANN (EVOLUÇÃO)
// ==========================================

class FriedmannSolver {
private:
    ModelConstants constants;
    EnergyMomentumTensor5D physics;
    const double PI = 3.14159265359;

public:
    FriedmannSolver(ModelConstants cons) : constants(cons), physics(cons) {}

    std::vector<double> get_derivatives(const UniverseState& state) {
        double rho_m = 0.5 * std::pow(state.a, -3.0); 
        
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
    
    // Calcula Equation of State (w) = Pressão / Densidade
    double get_w_parameter(const UniverseState& s) {
        double V = physics.potential_pNGB(s.phi);
        double Kinetic = 0.5 * s.d_phi * s.d_phi;
        
        double Pressure = Kinetic - V;
        double Density = Kinetic + V;

        if (Density == 0) return 0.0;
        return Pressure / Density;
    }
};

// ==========================================
// MAIN
// ==========================================

int main() {
    std::cout << "=== UCDT-R* UNIFICATION TEST ===" << std::endl;
    std::cout << "Monitorando Parametro w (w=-1: DE, w=0: DM)\n" << std::endl;

    ModelConstants params;
    params.G5 = 0.1;         
    params.Lambda_UV = 1.0;  
    params.f_a = 5.0;        
    
    FriedmannSolver solver(params);

    UniverseState u;
    u.t = 0.0;
    u.a = 0.01;  
    u.phi = 3.0; // Inicia longe do mínimo (potencial domina)
    u.d_phi = 0.0;

    double dt = 0.01;
    int steps = 50; 

    std::cout << std::left << std::setw(8) << "Tempo" 
              << std::setw(12) << "Fator(a)" 
              << std::setw(12) << "Campo(phi)" 
              << std::setw(12) << "w_EOS" 
              << std::setw(15) << "  Status" << std::endl;
    std::cout << "------------------------------------------------------------" << std::endl;

    for(int i=0; i <= steps; i++) {
        double w = solver.get_w_parameter(u);
        
        std::string status = "Misto";
        if (w < -0.9) status = "[DARK ENERGY]";
        else if (w > -0.1 && w < 0.1) status = "[DARK MATTER]";
        else if (w > 0.8) status = "[KINETIC]";

        std::cout << std::left << std::setw(8) << std::fixed << std::setprecision(2) << u.t 
                  << std::setw(12) << std::setprecision(4) << u.a 
                  << std::setw(12) << std::setprecision(4) << u.phi 
                  << std::setw(12) << std::setprecision(3) << w 
                  << status << std::endl;

        u = solver.step_rk4(u, dt);
    }
    
    return 0;
}