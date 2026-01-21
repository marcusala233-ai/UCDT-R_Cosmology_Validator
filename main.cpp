#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include "EnergyMomentumTensor5D.hpp"

// Função auxiliar para imprimir matrizes de forma legível
void print_tensor(const std::vector<std::vector<double>>& tensor, const std::string& name) {
    std::cout << "--- Tensor: " << name << " ---" << std::endl;
    for (const auto& row : tensor) {
        std::cout << "[ ";
        for (double val : row) {
            std::cout << std::setw(12) << std::scientific << val << " ";
        }
        std::cout << "]" << std::endl;
    }
    std::cout << "--------------------------" << std::endl;
}

int main() {
    std::cout << "=== UCDT-R COSMOLOGY VALIDATOR ===" << std::endl;
    std::cout << "Validacao da Condicao SNMC (Compactificacao Estabilizada)\n" << std::endl;

    // 1. Configuração dos Parâmetros (baseados no refinamento pNGB)
    ModelConstants constants;
    constants.Lambda_UV = 1.0e27; 
    constants.f_a = 1.0e28;       
    constants.G5 = 1.0;           
    constants.phi_0 = 1.0e28;     

    EnergyMomentumTensor5D tensorCalc(constants);

    // 2. Métrica 5D (Assinatura 3,2: -, +, +, +, -)
    std::vector<std::vector<double>> metric(5, std::vector<double>(5, 0.0));
    std::vector<std::vector<double>> inv_metric(5, std::vector<double>(5, 0.0));

    // Diagonal: t, x, y, z, tau
    metric[0][0] = -1.0; 
    metric[1][1] = 1.0;  
    metric[2][2] = 1.0;  
    metric[3][3] = 1.0;  
    metric[4][4] = -1.0; // Tempo extra negativo

    for(int i=0; i<5; i++) inv_metric[i][i] = 1.0 / metric[i][i];

    // 3. Estado do Campo (Teste)
    FieldState state;
    state.phi = 0.5 * constants.phi_0; 
    state.d_phi = {1.0e5, 0.0, 0.0, 0.0, 1.0e-2}; // Pequena variação no tempo extra

    // 4. Fluido de Matéria (Vácuo para este teste)
    std::vector<std::vector<double>> fluid_T(5, std::vector<double>(5, 0.0));

    // 5. Cálculo e Verificação
    try {
        auto T_total = tensorCalc.compute(metric, inv_metric, state, fluid_T);
        print_tensor(T_total, "T_AB (Total Energy-Momentum)");

        double T_44 = T_total[4][4];
        std::cout << "\n>>> VERIFICACAO SNMC (T_44) <<<" << std::endl;
        std::cout << "Valor em T_44 (Densidade na dimensao extra): " << T_44 << std::endl;
        
        if (T_44 < 0) {
            std::cout << "[ALERTA] T_44 Negativo. Risco de instabilidade fantasma." << std::endl;
        } else {
            std::cout << "[SUCESSO] T_44 Estavel. Condicao SNMC respeitada." << std::endl;
        }

    } catch (const std::exception& e) {
        std::cerr << "Erro: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}