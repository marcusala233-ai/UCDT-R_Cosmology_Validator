#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <string>

// ==========================================
// 1. ESTRUTURAS FÍSICAS
// ==========================================

struct ModelConstants {
    double Lambda_UV; 
    double f_a;       
    double G5;        
};

struct UniverseState {
    double t;      
    double a;      
    double phi;    
    double d_phi;
};

// ==========================================
// 2. TENSOR DE ENERGIA-MOMENTO
// ==========================================

class EnergyMomentumTensor5D {
private:
    ModelConstants constants;

public:
    EnergyMomentumTensor5D(ModelConstants cons) : constants(cons) {}

    double potential(double phi) const {
        return std::pow(constants.Lambda_UV, 4) * (1.0 - std::cos(phi / constants.f_a));
    }

    double d_potential(double phi) const {
        return (std::pow(constants.Lambda_UV, 4) / constants.f_a) * std::sin(phi / constants.f_a);
    }
};

// ==========================================
// 3. ANALISADOR DE CICLOS (NOVO!)
// ==========================================

class CycleAnalyzer {
private:
    EnergyMomentumTensor5D physics;
    
    // Acumuladores
    double sum_P;
    double sum_rho;
    double sum_cs2_inst; // Velocidade som instantanea
    int count;

public:
    CycleAnalyzer(ModelConstants cons) : physics(cons), sum_P(0), sum_rho(0), sum_cs2_inst(0), count(0) {}

    // Adiciona um ponto de dados ao ciclo atual
    void accumulate(const UniverseState& state) {
        double V = physics.potential(state.phi);
        double K = 0.5 * state.d_phi * state.d_phi;
        
        double P = K - V;
        double rho = K + V;
        
        sum_P += P;
        sum_rho += rho;
        
        // Cs2 instantaneo (apenas para estatística, não determina física macro)
        if (rho > 1e-10) sum_cs2_inst += (P/rho); 
        
        count++;
    }

    // Fecha o ciclo e retorna médias
    struct CycleResult {
        double w_avg;    // Equation of State Média
        double rho_avg;  // Densidade Média
        std::string regime;
        bool is_structure_forming;
    };

    CycleResult compute_and_reset() {
        if (count == 0) return {0,0,"ERRO", false};

        double avg_P = sum_P / count;
        double avg_rho = sum_rho / count;
        double w_avg = avg_P / avg_rho;

        // Reset
        sum_P = 0; sum_rho = 0; sum_cs2_inst = 0; count = 0;

        // Classificação do Regime Baseada na Média
        std::string regime;
        bool forming = false;

        if (w_avg < -0.3) {
            regime = "[ENERGIA ESCURA] (Repulsivo)";
        } else if (std::abs(w_avg) < 0.15) {
            // Este é o alvo! w ~ 0
            regime = "** MATERIA ESCURA ** (Aglomeravel)";
            forming = true;
        } else {
            regime = "[RADIACAO/CINETICA] (Quente)";
        }

        return {w_avg, avg_rho, regime, forming};
    }
    
    int get_count() const { return count; }
};

// ==========================================
// 4. SOLVER
// ==========================================

class FriedmannSolver {
private:
    ModelConstants constants;
    EnergyMomentumTensor5D physics;
    const double PI = 3.14159265359;

public:
    FriedmannSolver(ModelConstants cons) : constants(cons), physics(cons) {}

    std::vector<double> get_derivatives(const UniverseState& s) {
        double V = physics.potential(s.phi);
        double rho = 0.5 * s.d_phi * s.d_phi + V;
        double H = std::sqrt((8.0 * PI * constants.G5 / 3.0) * rho);
        double dV = physics.d_potential(s.phi);
        double friction = 3.0 * H * s.d_phi;
        
        return { s.a * H, s.d_phi, -friction - dV };
    }

    UniverseState step_rk4(UniverseState curr, double dt) {
        UniverseState next = curr;
        auto k1 = get_derivatives(curr);
        
        UniverseState s2 = curr; s2.t += 0.5*dt; s2.a += 0.5*dt*k1[0]; s2.phi += 0.5*dt*k1[1]; s2.d_phi += 0.5*dt*k1[2];
        auto k2 = get_derivatives(s2);
        
        UniverseState s3 = curr; s3.t += 0.5*dt; s3.a += 0.5*dt*k2[0]; s3.phi += 0.5*dt*k2[1]; s3.d_phi += 0.5*dt*k2[2];
        auto k3 = get_derivatives(s3);
        
        UniverseState s4 = curr; s4.t += dt; s4.a += dt*k3[0]; s4.phi += dt*k3[1]; s4.d_phi += dt*k3[2];
        auto k4 = get_derivatives(s4);

        next.t += dt;
        next.a += (dt/6.0)*(k1[0]+2*k2[0]+2*k3[0]+k4[0]);
        next.phi += (dt/6.0)*(k1[1]+2*k2[1]+2*k3[1]+k4[1]);
        next.d_phi += (dt/6.0)*(k1[2]+2*k2[2]+2*k3[2]+k4[2]);
        return next;
    }
};

// ==========================================
// 5. MAIN
// ==========================================

int main() {
    std::cout << "=== UCDT-R* CYCLE AVERAGED STRUCTURE TEST ===" << std::endl;
    std::cout << "Metodo: Media movel de w = <P> / <rho> para eliminar ruido de oscilacao.\n" << std::endl;

    ModelConstants params;
    params.G5 = 0.1; params.Lambda_UV = 1.5; params.f_a = 0.5; // Turbo
    
    FriedmannSolver solver(params);
    CycleAnalyzer analyzer(params);

    UniverseState u;
    u.t = 0.0; u.a = 0.1; u.phi = 2.0; u.d_phi = 0.0;

    double dt = 0.02; // Passo pequeno para capturar a oscilação
    int total_steps = 1000;
    int cycle_window = 40; // Quantos passos formam um ciclo de média

    std::cout << std::left << std::setw(10) << "Tempo(fim)" 
              << std::setw(15) << "<w> (EOS)" 
              << "Regime Efetivo (Macro)" << std::endl;
    std::cout << "-------------------------------------------------------------" << std::endl;

    for(int i=0; i <= total_steps; i++) {
        // 1. Acumular dados do passo atual
        analyzer.accumulate(u);

        // 2. Se fechou o ciclo, calcular média e imprimir
        if (analyzer.get_count() >= cycle_window) {
            auto res = analyzer.compute_and_reset();
            
            // Só imprime se mudou algo ou a cada N ciclos para não poluir
            std::cout << std::left << std::setw(10) << std::fixed << std::setprecision(2) << u.t 
                      << std::setw(15) << std::setprecision(4) << res.w_avg 
                      << res.regime << std::endl;
        }

        // 3. Avançar física
        u = solver.step_rk4(u, dt);
        if (u.a > 200.0) break;
    }

    return 0;
}