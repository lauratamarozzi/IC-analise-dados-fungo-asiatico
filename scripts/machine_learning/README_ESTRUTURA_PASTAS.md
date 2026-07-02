# 📁 ESTRUTURA DE PASTAS RECOMENDADA

## Organização do Projeto

```
ic-analise-dados-ferrugem-asiatica/
│
├── scripts/
│   ├── data_processing/                    ← Scripts de processamento (já existem)
│   │   ├── subtracao_grafos_v3.py
│   │   └── teste_fisher_fdr.py
│   │
│   └── machine_learning/                   ← NOVA SUBPASTA (todo ML aqui!)
│       ├── comparacao_manual_vs_deep.py    ← Script principal
│       ├── requirements_ml_comparison.txt  ← Dependências
│       ├── README_COMPARACAO_ML.md         ← Documentação
│       ├── COMO_RODAR.txt
│       └── SUMARIO_ENTREGA.txt
│
├── dados_limpos/                           ← Dados brutos
│   ├── balanceado/
│   │   └── K8108_TEs_balanceado.fasta
│   └── K8108_genes_nao_TEs.fasta
│
└── resultados/                             ← Saídas dos scripts
    ├── excel/                              (subtracao_grafos_v3.py)
    │   └── features_kmers_v3_bidirecional.xlsx
    │
    ├── csv/                                (teste_fisher_fdr.py)
    │   ├── kmers_completo_k15.csv
    │   ├── teste_fisher_resultados.csv
    │   └── teste_fisher_candidatos_fortes.csv
    │
    ├── graficos/                           (teste_fisher_fdr.py)
    │   ├── distribuicao_pvalores.png
    │   ├── vulcao_fisher.png
    │   └── top20_kmers_candidato_forte.png
    │
    └── machine_learning/                   ← NOVA SUBPASTA (comparacao_manual_vs_deep.py)
        ├── ml_manual/
        │   ├── modelo_random_forest.pkl
        │   ├── features_importancia.png
        │   ├── matriz_confusao_manual.png
        │   └── roc_curve_manual.png
        │
        ├── ml_deep/
        │   ├── modelo_cnn_1d.h5
        │   ├── historico_treinamento.png
        │   ├── matriz_confusao_deep.png
        │   └── roc_curve_deep.png
        │
        └── comparacao/                     ← LEIA ISTO PRIMEIRO!
            ├── relatorio_comparativo.txt
            ├── comparacao_metricas.csv
            ├── comparacao_visual.png
            ├── matriz_confusao_manual.png
            ├── matriz_confusao_deep.png
            ├── roc_curve_manual.png
            └── roc_curve_deep.png
```

---

## ✅ PASSO-A-PASSO PARA ORGANIZAR

### **Passo 1: Criar a subpasta `machine_learning/` em `scripts/`**

```bash
# Navegue até a pasta scripts
cd scripts

# Crie a subpasta
mkdir machine_learning

# Ou no Windows (PowerShell):
New-Item -ItemType Directory -Name machine_learning
```

### **Passo 2: Copiar os 5 arquivos para `scripts/machine_learning/`**

Copie esses arquivos (que você baixou) para **`scripts/machine_learning/`**:

- ✅ `comparacao_manual_vs_deep.py`
- ✅ `requirements_ml_comparison.txt`
- ✅ `README_COMPARACAO_ML.md`
- ✅ `COMO_RODAR.txt`
- ✅ `SUMARIO_ENTREGA.txt`

Após: sua pasta deve ficar assim:

```
scripts/machine_learning/
├── comparacao_manual_vs_deep.py
├── requirements_ml_comparison.txt
├── README_COMPARACAO_ML.md
├── COMO_RODAR.txt
└── SUMARIO_ENTREGA.txt
```

### **Passo 3: Criar a subpasta `data_processing/` em `scripts/`** (opcional, para organizar melhor)

Se quiser organizar melhor, pode mover os scripts antigos:

```bash
cd scripts

# Crie a subpasta
mkdir data_processing

# Mova os scripts antigos (ou copie-os)
# Windows PowerShell:
Move-Item subtracao_grafos_v3.py data_processing/
Move-Item teste_fisher_fdr.py data_processing/

# Ou Linux/Mac:
mv subtracao_grafos_v3.py data_processing/
mv teste_fisher_fdr.py data_processing/
```

---

## 🚀 COMO RODAR

### **Opção A: De dentro da subpasta `machine_learning/`**

```bash
# Navegue até a subpasta
cd scripts/machine_learning

# Instale dependências (primeira vez)
pip install -r requirements_ml_comparison.txt

# Rode o script
python comparacao_manual_vs_deep.py
```

### **Opção B: De qualquer lugar (se quiser)**

```bash
# De onde estiver
cd scripts/machine_learning
python comparacao_manual_vs_deep.py
```

---

## 📊 ONDE VER OS RESULTADOS

Após rodar, abra:

```
resultados/machine_learning/comparacao/relatorio_comparativo.txt
```

Ou visualize os gráficos em:

```
resultados/machine_learning/
├── ml_manual/features_importancia.png
├── ml_deep/historico_treinamento.png
└── comparacao/comparacao_visual.png
```

---

## ⚡ TROUBLESHOOTING

### **Erro: "File not found: teste_fisher_candidatos_fortes.csv"**

**Causa:** O script não conseguiu encontrar os arquivos de entrada.

**Solução:**
1. Verifique que você está rodando de `scripts/machine_learning/`
2. Verifique que `resultados/csv/teste_fisher_candidatos_fortes.csv` existe
3. Se não existir, rode antes:
   ```bash
   cd scripts/data_processing
   python teste_fisher_fdr.py
   ```

### **Erro: "ModuleNotFoundError"**

**Solução:**
```bash
cd scripts/machine_learning
pip install -r requirements_ml_comparison.txt
```

### **Script muito lento**

Normal em CPU (~15 minutos). Se quiser acelerar, instale GPU:
```bash
pip install tensorflow[and-cuda]
```

---

## 📝 RESUMO RÁPIDO

```bash
# 1. Organize a estrutura de pastas (conforme descrito acima)

# 2. Copie os 5 arquivos para scripts/machine_learning/

# 3. Navegue até lá
cd scripts/machine_learning

# 4. Instale dependências
pip install -r requirements_ml_comparison.txt

# 5. Rode o script
python comparacao_manual_vs_deep.py

# 6. Leia o resultado
# Abra: resultados/machine_learning/comparacao/relatorio_comparativo.txt
```

---

## ✨ VANTAGENS DESTA ESTRUTURA

✅ **Melhor organização:**
  - Data processing scripts em uma subpasta
  - Machine Learning scripts em outra subpasta
  - Fácil encontrar o que precisa

✅ **Fácil reutilizar:**
  - Copie a subpasta `machine_learning/` se quiser usar em outro projeto
  - Tudo que precisa está junto

✅ **Escalável:**
  - Se quiser adicionar mais modelos de ML depois, fica fácil
  - Basta adicionar outro script na subpasta `machine_learning/`

✅ **Profissional:**
  - Reflete boas práticas de organização em projetos reais
  - Útil para quando virar código de produção

---

**Pronto? Crie a subpasta `machine_learning/`, copie os arquivos e rode!** 🚀