

# UCDT-R* Cosmology Validator
### Unified Geometric Unification (Refactored) - 5D pNGB Model

[![DOI](https://zenodo.org/badge/1134625974.svg)](https://doi.org/10.5281/zenodo.18290985)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)

![UCDT-R Theory Visual Compendium](ucdt_visual_compendium.png)

> **Status do Projeto:** ✅ Unificação de Fundo Validada (Background Dynamics) | 🚧 Análise de Estrutura em Progresso (Jeans Instability)

## 🌌 O Compêndio Teórico UCDT-R*

A **UCDT-R* (Refactored)** é uma extensão da Relatividade Geral definida em um hiperespaço de assinatura $(3,2)$, projetada para resolver o problema da Energia Escura e Matéria Escura através de um único mecanismo geométrico-escalar, livre das instabilidades (fantasmas) e do fine-tuning extremo do modelo original.

### 1. O Mecanismo de Unificação (Potencial pNGB)
Diferente do modelo polinomial instável ($\lambda \phi^4$), a UCDT-R* utiliza um potencial protegido por simetria de deslocamento (shift symmetry):

$$V(\phi) = \Lambda_{UV}^4 \left[ 1 - \cos\left(\frac{\phi}{f_a}\right) \right]$$

* **Fase Inflacionária (Topo):** Quando $\phi \approx 0$ (ou $2\pi n$), o campo rola lentamente (Slow-Roll), gerando pressão negativa $P \approx -\rho$, mimetizando **Energia Escura**.
* **Fase Material (Fundo):** Quando o campo cai no poço de potencial, ele oscila com frequência $m_{eff} \gg H$. A média temporal da pressão se anula $\langle P \rangle \approx 0$, comportando-se como **Matéria Escura Fria (CDM)**.

### 2. Geometria 5D Estabilizada (SNMC)
O modelo opera em uma variedade $\mathcal{M}^{(3,2)}$. Para evitar modos fantasmas (ghosts) típicos de dois tempos, impomos a **Condição de Compactificação Não-Mínima Estabilizada (SNMC)** no tensor de energia-momento 5D:

$$T_{44}^{(5)} \ge 0 \quad \text{(Condição de Estabilidade no Código)}$$

Isto garante que a dimensão extra $\tau$ não viole a unitariedade no setor efetivo 4D.

---

## 💻 Resultados da Simulação (Prova de Conceito)

O validador numérico (`FriedmannSolver` via Runge-Kutta 4) confirmou a transição de fase utilizando os parâmetros "Turbo" ($\Lambda_{UV}=1.5, f_a=0.5$).

### Evidência de Unificação (Log de Saída)
Abaixo, o registro da Equação de Estado ($w = P/\rho$) durante a evolução cósmica simulada:

```text
t= 0.00  [-O------------------|-------------------] DE (Inflacao/Energia Escura)
t= 0.50  [--------O-----------|-------------------] Queda do Potencial
t= 0.75  [--------------------|O------------------] ** DM (Materia Escura) **
t= 1.00  [--------------------|----------------O--] Oscilacao Cinetica
t= 1.50  [--------------------O-------------------] ** DM (Materia Escura) **
...
(Média temporal em t > 1.0 converge para w = 0)

📊 Estrutura do Repositório
• /src: Código fonte das equações cosmológicas e likelihoods.
• /data: Datasets observacionais (Pantheon+, SH0ES, etc).
• /plots: Gráficos gerados (Corner plots e Hubble Diagrams).
• CITATION.cff: Arquivo de metadados para citação acadêmica.

🚀 Como Reproduzir os Resultados
Requer compilador C++ e Python (opcional para plotagem HD).

1. Compilar e Rodar a Simulação
g++ main.cpp -o ucdt_validator
./ucdt_validator

2. Gerar Gráficos de Alta Resolução
python plot_ucdt.py
sso gerará o arquivo ucdt_final_proof.png com a evolução detalhada de $a(t)$, $\phi(t)$ e $w(t)$.

🔮 Próximos Passos:
Estrutura em Grande EscalaA próxima fase do desenvolvimento focará na Instabilidade de Jeans.Objetivo: Provar que a "sound speed" efetiva ($c_s^2$) cai para zero na fase de oscilação.
Módulo: StructureValidator.cpp (Em desenvolvimento).

🤝 Citação
Se você utilizar este software ou a teoria UCDT-R em sua pesquisa, por favor, cite utilizando o DOI arquivado no Zenodo (botão acima) ou o arquivo CITATION.cff presente neste repositório.

📝 Licença
Distribuído sob a licença MIT. Veja o arquivo LICENSE para mais informações.
Autor: Marcus Ala Pedreira Roriz
Pesquisador em Física Teórica e Cosmologia Computacional
