

# UCDT-R Cosmology Validator
### Validação Cosmológica da Teoria da Unificação de Campos Quadridimensionais (Restrita)

[![DOI](https://zenodo.org/badge/1134625974.svg)](https://doi.org/10.5281/zenodo.18290985)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)

## 🔭 Sobre o Projeto

Este repositório contém o código-fonte e os dados utilizados para a validação estatística e cosmológica do modelo **UCDT-R (Unification of Quadridimensional Fields Theory - Restricted)**.

O objetivo principal deste projeto é testar a viabilidade de um modelo de universo baseado em uma variedade espaço-temporal de 5 dimensões, confrontando suas previsões teóricas com dados observacionais recentes da cosmologia (Supernovas Tipo Ia, BAO e CMB).

## 📐 Fundamentação Teórica

A UCDT-R propõe uma extensão da Relatividade Geral onde os campos fundamentais são unificados através de uma geometria pentadimensional. Diferente dos modelos padrão $\Lambda$CDM, a UCDT-R sugere que a aceleração cósmica pode emergir naturalmente da geometria extra, sem a necessidade exclusiva de uma Constante Cosmológica exótica.

Este validador foca em:
1.  **Métrica 5D:** Análise da evolução do fator de escala $a(t)$ sob as equações de campo da UCDT.
2.  **Parâmetros Cosmológicos:** Restrição dos parâmetros livres do modelo ($\Omega_m$, $\Omega_k$, e os parâmetros extras da UCDT).
3.  **Teste de Ajuste:** Comparação da Luminosidade-Distância ($d_L$) prevista pelo modelo contra o catálogo Pantheon+.

## ⚙️ Metodologia e Algoritmos

O núcleo deste validador utiliza **Inferência Bayesiana** via simulações de Monte Carlo via Cadeias de Markov (MCMC).

* **Linguagem:** Python 3.10+
* **Amostragem MCMC:** `emcee` (The MCMC Hammer)
* **Visualização de Posteriores:** `corner.py`
* **Cálculo Numérico:** `NumPy`, `SciPy`, `Astropy`

## 🚀 Instalação e Uso

Para reproduzir os resultados ou testar o modelo com novos dados:

```bash
# Clone este repositório
git clone [https://github.com/marcusala233-ai/UCDT-R_Cosmology_Validator.git](https://github.com/marcusala233-ai/UCDT-R_Cosmology_Validator.git)

# Entre na pasta
cd UCDT-R_Cosmology_Validator

# Instale as dependências
pip install -r requirements.txt

📊 Estrutura do Repositório
• /src: Código fonte das equações cosmológicas e likelihoods.
• /data: Datasets observacionais (Pantheon+, SH0ES, etc).
• /plots: Gráficos gerados (Corner plots e Hubble Diagrams).
• CITATION.cff: Arquivo de metadados para citação acadêmica.

🤝 Citação
Se você utilizar este software ou a teoria UCDT-R em sua pesquisa, por favor, cite utilizando o DOI arquivado no Zenodo (botão acima) ou o arquivo CITATION.cff presente neste repositório.

📝 Licença
Distribuído sob a licença MIT. Veja o arquivo LICENSE para mais informações.
Autor: Marcus Ala Pedreira Roriz
Pesquisador em Física Teórica e Cosmologia Computacional
