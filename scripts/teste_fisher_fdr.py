"""
teste_fisher_fdr.py
===================

Testa TODOS os k-mers para significância estatística via Teste Exato de Fisher,
depois aplica correção FDR para não cair em armadilha de comparações múltiplas.

O Teste Exato de Fisher monta uma tabela de contingência 2×2:

                    | k-mer presente | k-mer ausente
    ─────────────────────────────────────────────
    Sequências TE   |       a        |      b
    Sequências Ctrl |       c        |      d

E calcula: qual a probabilidade de essa distribuição ocorrer por ACASO?

Saída: CSV com p-valores brutos e ajustados (FDR). K-mers com p-ajustado < 0.05
são estatisticamente significativos — candidatos para a etapa de ML.

DEPENDÊNCIAS
─────────────
pip install biopython pandas scipy openpyxl rich
(NÃO precisa mais de "statsmodels" — a correção FDR Benjamini-Hochberg foi
reimplementada localmente neste arquivo, veja benjamini_hochberg().
Isso corrige o ModuleNotFoundError: No module named 'statsmodels'.)

USO
───
python teste_fisher_fdr.py

ENTRADA
───────
- resultados/excel/features_kmers_v3_bidirecional.xlsx, gerado pelo
  subtracao_grafos_v3.py (aba "completo_k15", com TODOS os k-mers)

SAÍDA (organizada em subpastas dentro de resultados/, mesma estrutura
do subtracao_grafos_v3.py)
─────
resultados/csv/teste_fisher_resultados.csv → p-valores brutos e ajustados
resultados/graficos/
   - distribuicao_pvalores.png (histograma bruto vs. ajustado)
   - vulcao_fisher.png (delta_TE vs. -log10(p_ajustado))
   - top20_kmers_significativos.png
"""

import os
import gc
import logging
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact
from Bio import SeqIO
from collections import defaultdict
from rich.logging import RichHandler
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, MofNCompleteColumn

logging.basicConfig(level="INFO", format="%(message)s", handlers=[RichHandler(show_time=False)])
log = logging.getLogger("rich")


# ─────────────────────────────────────────────────────────────────────────────
# CORREÇÃO FDR (BENJAMINI-HOCHBERG) — implementação local
# ─────────────────────────────────────────────────────────────────────────────
# O script original dependia de `statsmodels.stats.multitest.multipletests`,
# mas esse pacote não vem instalado por padrão (era o ModuleNotFoundError que
# você recebeu) e é uma dependência pesada só para esta única função. Como o
# algoritmo Benjamini-Hochberg é simples, foi reimplementado aqui — assim o
# script roda só com biopython/pandas/scipy/matplotlib/rich, que você já tem.

def benjamini_hochberg(p_values, alpha=0.05):
    """
    Reimplementação do método FDR Benjamini-Hochberg (equivalente a
    `multipletests(p_values, alpha=alpha, method="fdr_bh")` do statsmodels).

    Para m testes com p-valores ordenados p(1) ≤ p(2) ≤ ... ≤ p(m):
        p_ajustado(i) = min_{j >= i} [ p(j) * m / j ]
    (o "min" à direita garante monotonicidade, exatamente como no BH clássico)

    Parâmetros
    ----------
    p_values : array-like de p-valores brutos (pode conter NaN)
    alpha    : limiar de significância

    Retorna
    -------
    reject       : array de bool, True onde p_ajustado <= alpha
    p_adjusted   : array de p-valores ajustados, na MESMA ordem da entrada
    """
    p_values = np.asarray(p_values, dtype=float)
    n = len(p_values)

    # Trata NaN como não-significativo (p ajustado = 1.0)
    nan_mask = np.isnan(p_values)
    p_values_safe = np.where(nan_mask, 1.0, p_values)

    ordem = np.argsort(p_values_safe)
    ranks = np.empty(n, dtype=int)
    ranks[ordem] = np.arange(1, n + 1)  # rank 1-indexado na ordem crescente

    p_ordenado = p_values_safe[ordem]
    m = n
    # p(j) * m / j, calculado na ordem crescente
    ajustado_ordenado = p_ordenado * m / np.arange(1, m + 1)

    # Impõe monotonicidade: da direita (maior p) para a esquerda,
    # cada valor não pode ser maior que o próximo à direita.
    ajustado_ordenado = np.minimum.accumulate(ajustado_ordenado[::-1])[::-1]
    ajustado_ordenado = np.clip(ajustado_ordenado, 0, 1)

    # Devolve para a ordem original de entrada
    p_adjusted = ajustado_ordenado[ranks - 1]
    p_adjusted = np.where(nan_mask, np.nan, p_adjusted)

    reject = np.where(nan_mask, False, p_adjusted <= alpha)

    return reject, p_adjusted


# ─────────────────────────────────────────────────────────────────────────────
# EXTRAÇÃO E TESTE ESTATÍSTICO
# ─────────────────────────────────────────────────────────────────────────────

def contar_sequencias_com_kmer(arquivo_fasta, k_mers_set):
    """
    Para cada sequência do FASTA, calcula quais k-mers estão presentes
    (presença = 1, ausência = 0).

    Retorna um dicionário: {kmer: set(índices de sequências onde ocorre)}
    """
    kmer_to_seqs = defaultdict(set)
    seq_index = 0

    with open(arquivo_fasta, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seq = str(record.seq).upper()
            if "N" in seq:
                seq_index += 1
                continue

            # Extrai todos os k-mers dessa sequência
            k = len(next(iter(k_mers_set)))  # pega o tamanho de um k-mer qualquer
            for i in range(0, len(seq) - k + 1):
                kmer = seq[i: i + k]
                if kmer in k_mers_set:
                    kmer_to_seqs[kmer].add(seq_index)

            seq_index += 1

    return dict(kmer_to_seqs), seq_index


def aplicar_teste_fisher(df_kmers, arquivo_tes, arquivo_ctrl):
    """
    Para cada k-mer da lista, realiza o Teste Exato de Fisher:

    Tabela 2×2:
    ┌─────────────────────┬──────────────┬──────────────┐
    │                     │ k-mer presnt │ k-mer ausent │
    ├─────────────────────┼──────────────┼──────────────┤
    │ Sequências TE       │      a       │      b       │
    │ Sequências Ctrl     │      c       │      d       │
    └─────────────────────┴──────────────┴──────────────┘

    p-valor = Pr(tabela observada | hipótese nula de independência)
    """
    k = len(df_kmers.iloc[0]["kmer"])

    log.info(f"🔍  Contando presença/ausência de k-mers nos FASTs...")

    with Progress(SpinnerColumn(),
                  TextColumn("[progress.description]{task.description}"),
                  BarColumn(), MofNCompleteColumn()) as prg:
        t1 = prg.add_task("Processando TEs...", total=1)
        kmer_to_seqs_te, n_seqs_te = contar_sequencias_com_kmer(
            arquivo_tes, set(df_kmers["kmer"])
        )
        prg.update(t1, advance=1)

        t2 = prg.add_task("Processando Controle...", total=1)
        kmer_to_seqs_ctrl, n_seqs_ctrl = contar_sequencias_com_kmer(
            arquivo_ctrl, set(df_kmers["kmer"])
        )
        prg.update(t2, advance=1)

    log.info(f"   TE:   {n_seqs_te:,} sequências")
    log.info(f"   Ctrl: {n_seqs_ctrl:,} sequências")

    # Monta tabelas 2×2 e calcula p-valores
    p_values = []

    for _, row in df_kmers.iterrows():
        kmer = row["kmer"]

        # Conta sequências com o k-mer
        seqs_com_kmer_te   = len(kmer_to_seqs_te.get(kmer, set()))
        seqs_com_kmer_ctrl = len(kmer_to_seqs_ctrl.get(kmer, set()))

        # Tabela 2×2
        a = seqs_com_kmer_te
        b = n_seqs_te - a
        c = seqs_com_kmer_ctrl
        d = n_seqs_ctrl - c

        # Teste Exato de Fisher (two-tailed)
        odds_ratio, p_valor = fisher_exact([[a, b], [c, d]], alternative="two-sided")

        p_values.append({
            "kmer": kmer,
            "seqs_TE_com_kmer": a,
            "seqs_TE_sem_kmer": b,
            "seqs_Ctrl_com_kmer": c,
            "seqs_Ctrl_sem_kmer": d,
            "odds_ratio": odds_ratio,
            "p_value_bruto": p_valor,
        })

    df_resultado = pd.DataFrame(p_values)

    # Aplica FDR (Benjamini-Hochberg) — implementação local, sem statsmodels
    reject, p_adjust = benjamini_hochberg(df_resultado["p_value_bruto"], alpha=0.05)

    df_resultado["p_value_ajustado"] = p_adjust
    df_resultado["significativo"] = reject

    return df_resultado


# ─────────────────────────────────────────────────────────────────────────────
# VISUALIZAÇÃO
# ─────────────────────────────────────────────────────────────────────────────

def plotar_distribuicao_pvalores(df, pasta_destino):
    """Histograma dos p-valores brutos e ajustados."""
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # p-valores brutos
    axes[0].hist(df["p_value_bruto"], bins=50, color="steelblue", edgecolor="black", alpha=0.7)
    axes[0].axvline(0.05, color="red", linestyle="--", linewidth=2, label="α = 0.05")
    axes[0].set_xlabel("p-valor bruto", fontsize=11)
    axes[0].set_ylabel("Frequência", fontsize=11)
    axes[0].set_title("Distribuição de p-valores brutos\n(antes de correção FDR)", fontsize=12)
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    # p-valores ajustados
    axes[1].hist(df["p_value_ajustado"], bins=50, color="coral", edgecolor="black", alpha=0.7)
    axes[1].axvline(0.05, color="red", linestyle="--", linewidth=2, label="α = 0.05")
    axes[1].set_xlabel("p-valor ajustado (FDR)", fontsize=11)
    axes[1].set_ylabel("Frequência", fontsize=11)
    axes[1].set_title("Distribuição de p-valores ajustados\n(após correção FDR)", fontsize=12)
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    caminho = os.path.join(pasta_destino, "distribuicao_pvalores.png")
    plt.savefig(caminho, dpi=300, bbox_inches="tight")
    plt.close()
    log.info(f"   📊 {os.path.basename(caminho)}")


def plotar_vulcao(df, pasta_destino):
    """Vulcão: delta_TE vs. -log10(p_ajustado). Marca k-mers significativos."""
    df_plot = df.copy()
    df_plot["-log10(p_ajustado)"] = -np.log10(df_plot["p_value_ajustado"] + 1e-300)

    fig, ax = plt.subplots(figsize=(10, 7))

    # Não-significativos (cinza)
    ns = df_plot[~df_plot["significativo"]]
    ax.scatter(ns["delta_TE"] * 1000, ns["-log10(p_ajustado)"],
               s=30, color="#cccccc", alpha=0.5, label="Não significativo (p > 0.05)")

    # Significativos (colorido)
    sig = df_plot[df_plot["significativo"]]
    ax.scatter(sig["delta_TE"] * 1000, sig["-log10(p_ajustado)"],
               s=80, c=sig["delta_TE"] * 1000, cmap="RdYlBu_r",
               alpha=0.8, edgecolors="w", linewidth=0.7,
               label="Significativo (p ≤ 0.05)")

    # Linhas de referência
    ax.axhline(-np.log10(0.05), color="red", linestyle="--", alpha=0.5, linewidth=1.5)
    ax.axvline(0, color="black", linestyle="-", alpha=0.3, linewidth=0.8)

    # Anotações dos top 10
    for _, row in df_plot.nlargest(10, "delta_TE").iterrows():
        ax.annotate(row["kmer"],
                    (row["delta_TE"] * 1000, row["-log10(p_ajustado)"]),
                    xytext=(5, 5), textcoords="offset points",
                    fontsize=8, fontweight="bold", color="#333")

    ax.set_xlabel("delta_TE (frequência relativa × 1000)", fontsize=11)
    ax.set_ylabel("-log10(p-valor ajustado)", fontsize=11)
    ax.set_title("Vulcão: Magnitude vs. Significância\n(K-mers discriminativos TE vs. Controle)", fontsize=12)
    ax.legend(loc="upper left")
    ax.grid(True, alpha=0.3)
    plt.tight_layout()

    caminho = os.path.join(pasta_destino, "vulcao_fisher.png")
    plt.savefig(caminho, dpi=300, bbox_inches="tight")
    plt.close()
    log.info(f"   📊 {os.path.basename(caminho)}")


def plotar_top_significativos(df, pasta_destino, top_n=20):
    """Barras horizontais dos top N k-mers mais significativos."""
    sig = df[df["significativo"]].nlargest(top_n, "delta_TE").sort_values("delta_TE")

    if sig.empty:
        log.warning("   ⚠️  Nenhum k-mer significativo encontrado.")
        return

    fig, ax = plt.subplots(figsize=(10, max(5, len(sig) * 0.4)))

    bars = ax.barh(sig["kmer"], sig["delta_TE"] * 1000, color="steelblue", edgecolor="white")
    ax.bar_label(bars, fmt="%.6f", padding=3, fontsize=8)
    ax.set_xlabel("delta_TE (frequência relativa × 1000)", fontsize=10)
    ax.set_title(f"Top {top_n} K-mers estatisticamente significativos", fontsize=12, pad=10)
    ax.grid(axis="x", linestyle="--", alpha=0.4)
    plt.tight_layout()

    caminho = os.path.join(pasta_destino, f"top{top_n}_kmers_significativos.png")
    plt.savefig(caminho, dpi=300, bbox_inches="tight")
    plt.close()
    log.info(f"   📊 {os.path.basename(caminho)}")


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

def main():
    caminho_script = os.path.dirname(os.path.abspath(__file__))
    pasta_dados = os.path.join(caminho_script, "..", "dados_limpos")
    pasta_result = os.path.join(caminho_script, "..", "resultados")

    # Mesma organização em subpastas usada pelo subtracao_grafos_v3.py
    pasta_excel    = os.path.join(pasta_result, "excel")
    pasta_csv      = os.path.join(pasta_result, "csv")
    pasta_graficos = os.path.join(pasta_result, "graficos")
    for pasta in (pasta_result, pasta_excel, pasta_csv, pasta_graficos):
        os.makedirs(pasta, exist_ok=True)

    arquivo_tes = os.path.join(pasta_dados, "balanceado", "K8108_TEs_balanceado.fasta")
    arquivo_ctrl = os.path.join(pasta_dados, "K8108_genes_nao_TEs.fasta")
    excel_input = os.path.join(pasta_excel, "features_kmers_v3_bidirecional.xlsx")

    # Valida inputs
    for arq in [arquivo_tes, arquivo_ctrl, excel_input]:
        if not os.path.exists(arq):
            log.error(f"❌  Arquivo não encontrado: {arq}")
            if arq == excel_input:
                log.error("   → Rode subtracao_grafos_v3.py antes deste script.")
            return

    # Lê a aba COMPLETA do Excel (todos os k-mers, não só o top 15) para
    # testar estatisticamente todo o universo de candidatos.
    log.info("📖 Lendo k-mers do Excel (aba completa)...")
    df_kmers = pd.read_excel(excel_input, sheet_name="completo_k15")
    log.info(f"   {len(df_kmers):,} k-mers para testar")

    # Aplica Teste de Fisher
    log.info("📊  Aplicando Teste Exato de Fisher + FDR...")
    df_resultado = aplicar_teste_fisher(df_kmers, arquivo_tes, arquivo_ctrl)

    # Merge com dados originais (delta_TE, etc)
    df_resultado = df_resultado.merge(
        df_kmers[["kmer", "delta_TE", "freq_TE_rel", "freq_Ctrl_rel"]],
        on="kmer", how="left"
    )

    # Ordena por p_value_ajustado
    df_resultado = df_resultado.sort_values("p_value_ajustado")

    # Resumo no terminal
    n_sig = df_resultado["significativo"].sum()
    log.info(f"\n{'─'*60}")
    log.info(f"  RESUMO ESTATÍSTICO")
    log.info(f"{'─'*60}")
    log.info(f"  Total de k-mers testados: {len(df_resultado):,}")
    log.info(f"  Significativos (p ≤ 0.05): {n_sig}")
    if n_sig > 0:
        log.info(f"  Taxa de descobertas: {100*n_sig/len(df_resultado):.1f}%")

    # Top 10 significativos
    if n_sig > 0:
        log.info(f"\n{'─'*60}")
        log.info(f"  TOP 10 K-MERS SIGNIFICATIVOS")
        log.info(f"{'─'*60}")
        for _, row in df_resultado[df_resultado["significativo"]].head(10).iterrows():
            log.info(f"  {row['kmer']}  "
                     f"p={row['p_value_ajustado']:.2e}  "
                     f"ΔTE={row['delta_TE']*1000:.6f}‰")

    # Salva CSV (todos os k-mers testados, com p-valor bruto e ajustado)
    output_csv = os.path.join(pasta_csv, "teste_fisher_resultados.csv")
    df_resultado.to_csv(output_csv, index=False)
    log.info(f"\n✅  CSV salvo: {output_csv}")

    # Gráficos
    log.info("📈  Gerando gráficos...")
    plotar_distribuicao_pvalores(df_resultado, pasta_graficos)
    plotar_vulcao(df_resultado, pasta_graficos)
    plotar_top_significativos(df_resultado, pasta_graficos, top_n=20)

    log.info(f"\n🎉  Concluído! Próximo passo: ML com os k-mers significativos")


if __name__ == "__main__":
    main()