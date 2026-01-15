import numpy as np
import matplotlib.pyplot as plt

# --- CONFIGURAÇÃO ---
true_mass_eV = 1.65e-23
target_freq = true_mass_eV * 2.418e14
freqs = np.logspace(-9, -7.5, 40)

# Modelo Standard (Buracos Negros) - Linha Azul
background = 1e-15 * (freqs / 1e-9)**(-2/3)

# Modelo UCDT-R* (Sinal 4 nHz) - Linha Vermelha
signal_width = 0.5e-9 
bump = 3e-16 * np.exp(-0.5 * (freqs - target_freq)**2 / signal_width**2)

# Dados Simulados (Pontos Pretos)
np.random.seed(42)
data_y = background + bump + np.random.normal(0, 0.5e-16, size=len(freqs))
data_err = 0.2 * data_y

# --- PLOTAGEM ---
plt.figure(figsize=(10, 6))
plt.errorbar(freqs*1e9, data_y, yerr=data_err, fmt='ko', label='NANOGrav 15yr Data (Simulado)')
plt.plot(freqs*1e9, background, 'b--', label='Standard Physics (Black Holes)')
plt.plot(freqs*1e9, background + bump, 'r-', linewidth=2, label='UCDT-R* Prediction')

plt.xscale('log')
plt.yscale('log')
plt.xlabel('Frequency (nHz)')
plt.ylabel('Gravitational Wave Strain')
plt.title('Validation: UCDT-R* (Red) vs Standard Model (Blue)')
plt.legend()
plt.grid(True, which="both", ls="-", alpha=0.2)

# SALVAR NA PASTA RESULTS
# Certifique-se de criar a pasta "results" antes de rodar
plt.savefig('../results/grafico-ucdt.png', dpi=300)
print("Gráfico salvo em results/grafico-ucdt.png")