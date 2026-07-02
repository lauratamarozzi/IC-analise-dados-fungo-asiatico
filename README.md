# XAI-Genomics: Explicabilidade em Elementos de Transposição

Projeto de Iniciação Científica (IC) focado na aplicação de **Inteligência Artificial Explicável (XAI)** e de **Engenharia de Dados** para a caracterização matemática de Elementos de Transposição (TEs) no genoma do fungo *Phakopsora pachyrhizi* (ferrugem asiática da soja). 

Desenvolvido junto ao **CyberLab** (Grupo de Pesquisa em Segurança Cibernética e Inteligência Artificial) da **Universidade Estadual de Londrina (UEL)**.

---

## O Desafio Biológico e Computacional

A ferrugem asiática da soja figura entre os patógenos agrícolas de maior impacto. A compreensão de seu genoma — especificamente dos Elementos de Transposição, que consistem em sequências de DNA capazes de se mover e induzir mutações — é fundamental para o desenvolvimento de defesas genéticas. 

A identificação clássica desses elementos baseia-se comumente em softwares de alinhamento com bancos de dados externos ou em algoritmos de modelagem preditiva opaca ("Caixa-Preta" ou *Black-Box*). Tais métodos indicam a localização das anomalias, porém falham em elucidar a regra estrutural subjacente. O presente projeto visa construir um pipeline computacional capaz de realizar a engenharia reversa desses padrões biológicos.

---

## Evolução Metodológica e Redirecionamento

O desenvolvimento científico em dados pressupõe processos iterativos. O desenho arquitetural deste estudo foi submetido a uma rigorosa validação de viabilidade:

### Abordagem Inicial: A Hipótese GRAMEP
* **Proposta Inicial:** Utilizar o software GRAMEP para o treinamento de um modelo de *Machine Learning* capaz de prever TEs diretamente na sequência de DNA.
* **Obstáculo Técnico:** Identificou-se uma limitação crítica referente à verdade de campo (*Ground Truth*). Os dados genômicos do isolado K8108 disponíveis já consistiam em anotações geradas por algoritmos preditivos prévios (pipeline REPET/TEdenovo).
* **Conclusão:** O treinamento de um modelo opaco a partir de dados sintetizados por outro sistema de mesma natureza resultaria na mera replicação de vieses algorítmicos, sem assegurar a extração de um aprendizado biológico autêntico.

### Redirecionamento: Inteligência Artificial Explicável (XAI)
A metodologia foi reestruturada visando a auditoria algorítmica e a explicabilidade estrutural. Em detrimento da simples predição de coordenadas dos TEs, o sistema emprega a **Teoria dos Grafos** conjugada a algoritmos interpretáveis (*White-Box*, como Árvores de Decisão e *Random Forest*). O intuito é demonstrar matematicamente quais assinaturas geométricas (K-mers) determinam a classificação de uma sequência como um Elemento de Transposição.

---

## Pipeline de Engenharia de Dados

O processamento analítico divide-se nas seguintes etapas fundamentais:

### 1. Extração de Características (Janela Deslizante)
* **Processamento de Complexidade $O(N)$:** Leitura computacionalmente eficiente de arquivos no formato FASTA através do uso de iteradores, prevenindo o esgotamento de memória (*Out of Memory* - OOM).
* **Fragmentação Sequencial:** A cadeia de DNA é segmentada em blocos (K-mers) adotando as resoluções de **K=4, 6, 8 e 15**.
* **Passo de Leitura (Step=1):** Substituição da leitura tradicional baseada em códons (step=3) pela leitura sobreposta contínua (step=1). Este ajuste metodológico assegura a captura de todas as arestas do grafo funcional, fundamentado no fato de que os TEs apresentam altas taxas de mutação e frequentemente se inserem em regiões não codificantes.

### 2. Análise Diferencial Global (Subtração de Grafos)
* Mapeamento estatístico das frequências de transição de nucleotídeos referentes à população de TEs em contraste com a população de Controle (genes regulares).
* Subtração de estruturas de dicionários (Tabela Hash) para a quantificação do nível de exclusividade de cada formação geométrica.
* Geração de gráficos de dispersão (*Scatter Plots*) a fim de fornecer provas visuais a respeito da separabilidade matemática entre as classes.

### 3. Matriz Tabular e Tratamento de Esparsidade
* Conversão de dados biológicos não estruturados (textuais) para conjuntos de dados tabulares estruturados (CSV).
* Elaboração de matrizes esparsas individuais. Nesse modelo, cada linha representa um Elemento de Transposição distinto e cada coluna denota uma característica estrutural quantificada, adotando o valor numérico **0** nos casos de ausência natural do K-mer.
* **Seleção e Normalização de Características:** Preparação do conjunto de dados para os modelos do Scikit-Learn mediante a aplicação do método `MinMaxScaler`, mitigando possíveis vieses oriundos de K-mers com frequências absolutas desproporcionais.

### 4. Benchmarking Computacional
* Aferição ininterrupta da utilização de recursos de hardware empregando as bibliotecas nativas `tracemalloc` e `time`. Esta etapa garante o monitoramento rigoroso das complexidades de tempo e espaço, validando a viabilidade e a escalabilidade da solução analítica em arquiteturas de hardware convencionais.

---

## Stack Tecnológico

* **Linguagem:** Python 3.13+
* **Engenharia de Dados e Análise:** Pandas, NumPy, Collections (DefaultDict / Counter)
* **Bioinformática:** Biopython (`SeqIO`)
* **Machine Learning:** Scikit-Learn (Decision Trees, Random Forest, Normalização)
* **Visualização e Monitoramento:** Matplotlib, Rich, Tracemalloc
