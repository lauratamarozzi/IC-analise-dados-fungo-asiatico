"""
subtracao_grafos_v3.py
======================

O QUE MUDOU EM RELAÇÃO À V2
-----------------------------
A v2 tinha um `if diff_rel > 0` que DESCARTAVA silenciosamente todos os
k-mers onde o Controle é mais frequente que o TE. Isso significava que
a planilha gerada só mostrava metade da informação: assinaturas de TE.

Esta v3 corrige isso de duas formas complementares:

  1. SUBTRAÇÃO COMPLETA BIDIRECIONAL
     - Calcula   ΔTE  = freq_TE  − freq_Ctrl  (positivo → TE > Ctrl)
     - Calcula   ΔCtrl = freq_Ctrl − freq_TE   (positivo → Ctrl > TE)
     - Mantém TODOS os k-mers com suporte mínimo, sem filtro de sinal.

  2. EXCEL COM 3 ABAS
     ┌─────────────────────────────────────────────────────┐
     │  Aba "top_TE"   → top N k-mers com ΔTE  mais alto  │
     │                   = assinaturas dos Transposons      │
     │                                                      │
     │  Aba "top_Ctrl" → top N k-mers com ΔCtrl mais alto  │
     │                   = assinaturas dos genes controle   │
     │                                                      │
     │  Aba "completo" → lista COMPLETA ordenada por ΔTE   │
     │                   (valores negativos incluídos)       │
     └─────────────────────────────────────────────────────┘

     Com TOP_N=15 você obtém exatamente os "15 primeiros de TE" e os
     "15 primeiros de Controle" que eram o objetivo original — mas a aba
     "completo" (e o CSV espelho dela) sempre traz TODOS os k-mers que
     passaram o filtro de suporte mínimo, não só os 15 primeiros.

  3. ORGANIZAÇÃO EM SUBPASTAS
     A partir desta versão, "resultados/" deixa de receber tudo solto e
     passa a ter:
       resultados/excel/     → o .xlsx com as 3 abas por k
       resultados/csv/       → kmers_completo_k{k}.csv (espelho da aba
                                "completo", para abrir sem precisar do Excel)
       resultados/graficos/  → todos os .png (scatter e barras)

DEPENDÊNCIAS
-------------
pip install biopython pandas openpyxl matplotlib rich

USO
---
python subtracao_grafos_v3.py

PARÂMETROS AJUSTÁVEIS (bloco main())
--------------------------------------
  K_TAMANHOS  → lista de valores de k a testar (default: [15])
  FREQ_MINIMA → suporte mínimo de ocorrências em AO MENOS UMA das classes
  TOP_N       → quantos k-mers exibir no resumo e salvar nas abas de topo
"""

import os
import gc
import logging
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
from collections import defaultdict
from rich.logging import RichHandler
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, MofNCompleteColumn

logging.basicConfig(level="INFO", format="%(message)s", handlers=[RichHandler(show_time=False)])
log = logging.getLogger("rich")


# ─────────────────────────────────────────────────────────────────────────────
# EXTRAÇÃO DE K-MERS
# ─────────────────────────────────────────────────────────────────────────────

def extrair_kmers(arquivo_fasta: str, k: int, step: int = 1):
    """
    Percorre um FASTA e conta todos os k-mers por janela deslizante.

    Retorna:
        contagem   (dict)  : {kmer: n_ocorrências}
        total_kmers (int)  : soma total de k-mers extraídos (necessário
                             para calcular frequência relativa depois).

    Sequências com N são puladas para não poluir os k-mers com bases
    indefinidas.
    """
    contagem = defaultdict(int)
    total_kmers = 0

    if not os.path.exists(arquivo_fasta):
        log.error(f"❌  Arquivo não encontrado: {arquivo_fasta}")
        return contagem, total_kmers

    with open(arquivo_fasta, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seq = str(record.seq).upper()
            if "N" in seq:
                continue
            for i in range(0, len(seq) - k + 1, step):
                kmer = seq[i: i + k]
                contagem[kmer] += 1
                total_kmers += 1

    return contagem, total_kmers


# ─────────────────────────────────────────────────────────────────────────────
# SUBTRAÇÃO BIDIRECIONAL
# ─────────────────────────────────────────────────────────────────────────────

def calcular_subtracao_bidirecional(
    mapa_te, total_te,
    mapa_ctrl, total_ctrl,
    freq_minima: int = 5,
):
    """
    Compara TE vs Controle por FREQUÊNCIA RELATIVA nas DUAS direções:

        delta_TE   = freq_TE_rel  − freq_Ctrl_rel   → positivo = TE > Ctrl
        delta_Ctrl = freq_Ctrl_rel − freq_TE_rel    → positivo = Ctrl > TE

    Não há mais filtro de sinal — todos os k-mers com suporte mínimo em
    ao menos uma das classes são mantidos no resultado.

    O filtro FREQ_MINIMA exige que o k-mer apareça pelo menos `freq_minima`
    vezes em AO MENOS UMA das classes.  K-mers com apenas 1-2 ocorrências
    em ambas as classes são ruído amostral, não assinatura estrutural.

    Parâmetros
    ----------
    mapa_te / mapa_ctrl : dict {kmer: contagem_absoluta}
    total_te / total_ctrl : int  (total de k-mers da classe)
    freq_minima : int  (suporte mínimo para manter o k-mer)

    Retorna
    -------
    DataFrame com colunas:
        kmer, freq_TE_abs, freq_Ctrl_abs,
        freq_TE_rel, freq_Ctrl_rel,
        delta_TE, delta_Ctrl
    """
    todos_kmers = set(mapa_te.keys()) | set(mapa_ctrl.keys())
    registros = []

    for kmer in todos_kmers:
        freq_te_abs   = mapa_te.get(kmer, 0)
        freq_ctrl_abs = mapa_ctrl.get(kmer, 0)

        # Descarta k-mers raros nas duas classes (ruído amostral)
        if freq_te_abs < freq_minima and freq_ctrl_abs < freq_minima:
            continue

        freq_te_rel   = freq_te_abs   / total_te   if total_te   else 0.0
        freq_ctrl_rel = freq_ctrl_abs / total_ctrl if total_ctrl else 0.0

        delta_te   = freq_te_rel   - freq_ctrl_rel   # TE − Ctrl
        delta_ctrl = freq_ctrl_rel - freq_te_rel     # Ctrl − TE  (= −delta_te)

        registros.append({
            "kmer"          : kmer,
            "freq_TE_abs"   : freq_te_abs,
            "freq_Ctrl_abs" : freq_ctrl_abs,
            "freq_TE_rel"   : freq_te_rel,
            "freq_Ctrl_rel" : freq_ctrl_rel,
            "delta_TE"      : delta_te,    # positivo → enriquecido em TE
            "delta_Ctrl"    : delta_ctrl,  # positivo → enriquecido em Controle
        })

    return pd.DataFrame(registros)


# ─────────────────────────────────────────────────────────────────────────────
# VISUALIZAÇÃO
# ─────────────────────────────────────────────────────────────────────────────

def plotar_scatter_bidirecional(df: pd.DataFrame, k: int, pasta_destino: str, top_n: int = 50):
    """
    Scatter plot mostrando TODOS os k-mers selecionados (com sinal),
    destacando os top_n por delta_TE (em verde/amarelo) e os top_n por
    delta_Ctrl (em azul), para visualizar as duas direções de enriquecimento.
    """
    if df.empty:
        log.warning(f"⚠️  DataFrame vazio para K={k}; scatter não gerado.")
        return

    # Converte para ‰ para os eixos ficarem legíveis
    x = df["freq_Ctrl_rel"] * 1000
    y = df["freq_TE_rel"]   * 1000

    fig, ax = plt.subplots(figsize=(11, 9))

    # Plota todos os k-mers em cinza claro como pano de fundo
    ax.scatter(x, y, s=15, color="#cccccc", alpha=0.4, label="Todos os k-mers", zorder=1)

    # Destaca top_n por delta_TE (enriquecidos em TE)
    top_te = df.nlargest(top_n, "delta_TE")
    sz_te  = (top_te["delta_TE"] / top_te["delta_TE"].max()) * 350 + 40
    sc_te  = ax.scatter(
        top_te["freq_Ctrl_rel"] * 1000,
        top_te["freq_TE_rel"]   * 1000,
        s=sz_te, c=top_te["delta_TE"], cmap="YlOrRd",
        alpha=0.85, edgecolors="w", linewidth=0.7,
        label=f"Top {top_n} enriquecidos em TE", zorder=3,
    )

    # Destaca top_n por delta_Ctrl (enriquecidos no Controle)
    top_ctrl = df.nlargest(top_n, "delta_Ctrl")
    sz_ctrl  = (top_ctrl["delta_Ctrl"] / top_ctrl["delta_Ctrl"].max()) * 350 + 40
    sc_ctrl  = ax.scatter(
        top_ctrl["freq_Ctrl_rel"] * 1000,
        top_ctrl["freq_TE_rel"]   * 1000,
        s=sz_ctrl, c=top_ctrl["delta_Ctrl"], cmap="Blues",
        alpha=0.85, edgecolors="w", linewidth=0.7,
        label=f"Top {top_n} enriquecidos em Controle", zorder=3,
    )

    # Linha neutra TE = Ctrl
    max_val = max(x.max(), y.max()) * 1.05
    ax.plot([0, max_val], [0, max_val], color="red", linestyle="--",
            alpha=0.5, linewidth=1.2, label="Linha neutra (TE = Ctrl)")

    # Anotações para os 10 melhores de cada grupo
    for _, row in top_te.head(10).iterrows():
        ax.annotate(row["kmer"],
                    (row["freq_Ctrl_rel"] * 1000, row["freq_TE_rel"] * 1000),
                    xytext=(6, 6), textcoords="offset points",
                    fontsize=7.5, color="#8B0000", fontweight="bold")

    for _, row in top_ctrl.head(10).iterrows():
        ax.annotate(row["kmer"],
                    (row["freq_Ctrl_rel"] * 1000, row["freq_TE_rel"] * 1000),
                    xytext=(6, -12), textcoords="offset points",
                    fontsize=7.5, color="#00008B", fontweight="bold")

    plt.colorbar(sc_te,   ax=ax, label="delta_TE (enriquecimento em TE)",   pad=0.01)
    plt.colorbar(sc_ctrl, ax=ax, label="delta_Ctrl (enriquecimento em Controle)", pad=0.07)

    ax.set_title(f"K-mers K={k}: enriquecimento bidirecional TE vs. Controle\n"
                 f"(frequências relativas em ‰)", fontsize=13, pad=12)
    ax.set_xlabel("Frequência relativa no Controle (‰)", fontsize=11)
    ax.set_ylabel("Frequência relativa nos TEs (‰)", fontsize=11)
    ax.legend(loc="upper left", fontsize=9)
    ax.grid(True, linestyle="--", alpha=0.35)
    plt.tight_layout()

    caminho = os.path.join(pasta_destino, f"scatter_bidirecional_k{k}.png")
    plt.savefig(caminho, dpi=300, bbox_inches="tight")
    plt.close()
    log.info(f"   📊 Gráfico salvo: {os.path.basename(caminho)}")


def plotar_barras_top(df: pd.DataFrame, coluna_delta: str, label: str,
                      k: int, pasta_destino: str, top_n: int = 15, cor: str = "steelblue"):
    """
    Gráfico de barras horizontais dos top_n k-mers por coluna_delta.
    Facilita a leitura visual dos candidatos mais discriminativos.
    """
    if df.empty:
        return

    subset = df.nlargest(top_n, coluna_delta)[["kmer", coluna_delta]].reset_index(drop=True)
    subset = subset.sort_values(coluna_delta)  # menor embaixo, maior no topo

    fig, ax = plt.subplots(figsize=(9, max(4, top_n * 0.45)))
    bars = ax.barh(subset["kmer"], subset[coluna_delta] * 1000, color=cor, edgecolor="white")
    ax.bar_label(bars, fmt="%.4f", padding=3, fontsize=8)
    ax.set_xlabel(f"{coluna_delta} × 1000 (diferença de frequência relativa em ‰)", fontsize=10)
    ax.set_title(f"Top {top_n} k-mers — {label} (K={k})", fontsize=12, pad=10)
    ax.grid(axis="x", linestyle="--", alpha=0.4)
    plt.tight_layout()

    nome_arquivo = f"barras_top{top_n}_{coluna_delta.lower()}_k{k}.png"
    caminho = os.path.join(pasta_destino, nome_arquivo)
    plt.savefig(caminho, dpi=300, bbox_inches="tight")
    plt.close()
    log.info(f"   📊 Gráfico salvo: {os.path.basename(caminho)}")


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

def main():
    # ── Parâmetros ajustáveis ──────────────────────────────────────────────
    K_TAMANHOS  = [15]   # Adicione outros valores se quiser comparar
    PASSO       = 1      # Step da janela deslizante
    FREQ_MINIMA = 5      # Suporte mínimo (em ao menos 1 classe)
    TOP_N       = 15     # Quantos k-mers exibir no resumo e nas abas de topo
    # ──────────────────────────────────────────────────────────────────────

    caminho_script = os.path.dirname(os.path.abspath(__file__))
    pasta_dados    = os.path.join(caminho_script, "..", "dados_limpos")
    pasta_result   = os.path.join(caminho_script, "..", "resultados")

    # ── Subpastas de organização dentro de "resultados" ────────────────────
    # excel/    → planilhas .xlsx (top_TE, top_Ctrl, completo)
    # csv/      → a lista COMPLETA de k-mers também em CSV, sem precisar
    #             abrir o Excel para achar a aba certa
    # graficos/ → todos os .png gerados
    pasta_excel    = os.path.join(pasta_result, "excel")
    pasta_csv      = os.path.join(pasta_result, "csv")
    pasta_graficos = os.path.join(pasta_result, "graficos")
    for pasta in (pasta_result, pasta_excel, pasta_csv, pasta_graficos):
        os.makedirs(pasta, exist_ok=True)

    # Usa o FASTA de TE já balanceado (gerado por balancear_volume_te.py)
    arquivo_tes  = os.path.join(pasta_dados, "balanceado", "K8108_TEs_balanceado.fasta")
    arquivo_ctrl = os.path.join(pasta_dados, "K8108_genes_nao_TEs.fasta")

    for arq, nome in [(arquivo_tes, "K8108_TEs_balanceado.fasta"),
                      (arquivo_ctrl, "K8108_genes_nao_TEs.fasta")]:
        if not os.path.exists(arq):
            log.error(f"❌  {nome} não encontrado em {os.path.dirname(arq)}")
            if "balanceado" in arq:
                log.error("   → Rode balancear_volume_te.py antes deste script.")
            return

    excel_path = os.path.join(pasta_excel, "features_kmers_v3_bidirecional.xlsx")

    with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
        for k in K_TAMANHOS:
            log.info("=" * 60)
            log.info(f"🚀  K={k} | Step={PASSO} | Freq. mínima={FREQ_MINIMA} | Top-N={TOP_N}")

            # 1. Extração de k-mers ─────────────────────────────────────────
            with Progress(SpinnerColumn(),
                          TextColumn("[progress.description]{task.description}"),
                          BarColumn(), MofNCompleteColumn()) as prg:
                t1 = prg.add_task("Extraindo k-mers de TEs (balanceado)...", total=1)
                mapa_te, total_te = extrair_kmers(arquivo_tes, k, PASSO)
                prg.update(t1, advance=1)

                t2 = prg.add_task("Extraindo k-mers de Controle...       ", total=1)
                mapa_ctrl, total_ctrl = extrair_kmers(arquivo_ctrl, k, PASSO)
                prg.update(t2, advance=1)

            log.info(f"   TE:      {total_te:>12,} k-mers  |  {len(mapa_te):>8,} únicos")
            log.info(f"   Ctrl:    {total_ctrl:>12,} k-mers  |  {len(mapa_ctrl):>8,} únicos")
            if total_ctrl:
                log.info(f"   Razão TE/Ctrl: {total_te / total_ctrl:.3f}x  "
                         f"(~1.0 = balanceado ✓)")

            # 2. Subtração bidirecional ─────────────────────────────────────
            log.info(f"📐  Calculando subtração bidirecional...")
            df = calcular_subtracao_bidirecional(
                mapa_te, total_te, mapa_ctrl, total_ctrl, FREQ_MINIMA
            )

            if df.empty:
                log.warning(f"⚠️  Nenhum k-mer passou o filtro para K={k}.")
                continue

            # Ordena pelo delta_TE para a aba principal
            df_sorted = df.sort_values("delta_TE", ascending=False).reset_index(drop=True)

            # 3. Resumo no terminal ─────────────────────────────────────────
            log.info(f"\n{'─'*55}")
            log.info(f"  TOP {TOP_N} ENRIQUECIDOS EM TE   (delta_TE mais positivo)")
            log.info(f"{'─'*55}")
            top_te = df_sorted.head(TOP_N)
            for _, row in top_te.iterrows():
                log.info(f"  {row['kmer']}  "
                         f"TE={row['freq_TE_rel']*1000:.4f}‰  "
                         f"Ctrl={row['freq_Ctrl_rel']*1000:.4f}‰  "
                         f"ΔTE={row['delta_TE']*1000:.4f}‰")

            log.info(f"\n{'─'*55}")
            log.info(f"  TOP {TOP_N} ENRIQUECIDOS EM CONTROLE  (delta_Ctrl mais positivo)")
            log.info(f"{'─'*55}")
            top_ctrl = df.nlargest(TOP_N, "delta_Ctrl")
            for _, row in top_ctrl.sort_values("delta_Ctrl", ascending=False).iterrows():
                log.info(f"  {row['kmer']}  "
                         f"Ctrl={row['freq_Ctrl_rel']*1000:.4f}‰  "
                         f"TE={row['freq_TE_rel']*1000:.4f}‰  "
                         f"ΔCtrl={row['delta_Ctrl']*1000:.4f}‰")

            # 4. Salva no Excel ─────────────────────────────────────────────
            # Aba 1: top N enriquecidos em TE
            top_te.to_excel(writer, sheet_name=f"top{TOP_N}_TE_k{k}", index=False)

            # Aba 2: top N enriquecidos em Controle
            top_ctrl.sort_values("delta_Ctrl", ascending=False).to_excel(
                writer, sheet_name=f"top{TOP_N}_Ctrl_k{k}", index=False
            )

            # Aba 3: lista completa ordenada por delta_TE (inclui valores negativos)
            # ESTA aba já contém TODOS os k-mers que passaram FREQ_MINIMA,
            # não só o TOP_N — é a lista completa que você pediu.
            df_sorted.to_excel(writer, sheet_name=f"completo_k{k}", index=False)

            # 4b. CSV dedicado com a lista COMPLETA (fora do Excel) ─────────
            # Facilita abrir/filtrar sem precisar navegar entre abas.
            csv_completo = os.path.join(pasta_csv, f"kmers_completo_k{k}.csv")
            df_sorted.to_csv(csv_completo, index=False)

            log.info(f"\n   ✅  {len(df_sorted):,} k-mers salvos (TODOS, não só o top {TOP_N}) | "
                     f"ΔTE varia de {df_sorted['delta_TE'].min()*1000:.4f}‰ "
                     f"a {df_sorted['delta_TE'].max()*1000:.4f}‰")
            log.info(f"   📄  Lista completa também em CSV: {os.path.basename(csv_completo)}")

            # 5. Gráficos ───────────────────────────────────────────────────
            log.info("📈  Gerando gráficos...")
            plotar_scatter_bidirecional(df_sorted, k, pasta_graficos, top_n=TOP_N)
            plotar_barras_top(df_sorted, "delta_TE",   "Enriquecido em TE",      k, pasta_graficos, TOP_N, cor="#C0392B")
            plotar_barras_top(df_sorted, "delta_Ctrl", "Enriquecido em Controle",k, pasta_graficos, TOP_N, cor="#2980B9")

            # 6. Libera memória ─────────────────────────────────────────────
            del mapa_te, mapa_ctrl, df, df_sorted
            gc.collect()

    log.info(f"\n🎉  Concluído!")
    log.info(f"   📊 Excel:    {excel_path}")
    log.info(f"   📄 CSV completo: {pasta_csv}")
    log.info(f"   📈 Gráficos: {pasta_graficos}")


if __name__ == "__main__":
    main()