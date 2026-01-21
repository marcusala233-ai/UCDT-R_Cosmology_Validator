#ifndef ENERGY_MOMENTUM_TENSOR_5D_HPP
#define ENERGY_MOMENTUM_TENSOR_5D_HPP

#include <vector>
#include <cmath>
#include <stdexcept>

// Estrutura para segurar as constantes cosmológicas e do modelo UCDT-R
struct ModelConstants {
    double Lambda_UV; // Escala de Energia Ultravioleta (eV)
    double f_a;       // Constante de decaimento da simetria pNGB (eV)
    double G5;        // Constante gravitacional 5D
    double phi_0;     // Valor esperado no vácuo (VEV)
};

// Estrutura para o estado do campo em um ponto do grid
struct FieldState {
    double phi;           // Valor do campo escalar
    std::vector<double> d_phi; // Derivadas parciais (0..4)
    // Nota: d_phi[4] é a derivada na direção temporal extra (tau)
};

class EnergyMomentumTensor5D {
private:
    ModelConstants constants;

public:
    EnergyMomentumTensor5D(ModelConstants cons) : constants(cons) {}

    // ========================================================================
    // Potencial pNGB (Pseudo Nambu-Goldstone Boson)
    // Implementa a correção de naturalidade descrita na refatoração UCDT-R*
    // ========================================================================
    double potential_pNGB(double phi) const {
        // V(phi) = Lambda^4 * (1 - cos(phi / f_a))
        return std::pow(constants.Lambda_UV, 4) * (1.0 - std::cos(phi / constants.f_a));
    }

    // ========================================================================
    // Cálculo das Componentes do Tensor T_AB
    // ========================================================================
    std::vector<std::vector<double>> compute(
        const std::vector<std::vector<double>>& metric,
        const std::vector<std::vector<double>>& inverse_metric,
        const FieldState& state,
        const std::vector<std::vector<double>>& fluid_T
    ) {
        std::vector<std::vector<double>> T_total(5, std::vector<double>(5, 0.0));

        // 1. Calcular o termo cinético contraído: (nabla phi)^2
        double kinetic_term = 0.0;
        for (int c = 0; c < 5; ++c) {
            for (int d = 0; d < 5; ++d) {
                kinetic_term += inverse_metric[c][d] * state.d_phi[c] * state.d_phi[d];
            }
        }

        // 2. Calcular o valor do Potencial V(phi)
        double V_val = potential_pNGB(state.phi);

        // 3. Montar o Tensor
        for (int a = 0; a < 5; ++a) {
            for (int b = 0; b < 5; ++b) {
                // Termo A: d_A phi * d_B phi
                double term_A = state.d_phi[a] * state.d_phi[b];

                // Termo B: Lagrangiana * Métrica
                double lagrangian = 0.5 * kinetic_term + V_val;
                double term_B = metric[a][b] * lagrangian;

                // T_AB_scalar = Term_A - Term_B
                double T_scalar = term_A - term_B;

                // T_AB_Total = T_scalar + T_matter
                T_total[a][b] = T_scalar + fluid_T[a][b];
            }
        }
        
        return T_total;
    }
};

#endif // ENERGY_MOMENTUM_TENSOR_5D_HPP