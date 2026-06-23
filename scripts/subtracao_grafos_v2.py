"""
subtracao_grafos_v2.py

O QUE MUDOU EM RELAÇÃO À VERSÃO ORIGINAL
-----------------------------------------
1. NORMALIZAÇÃO: a subtração agora usa FREQUÊNCIA RELATIVA
   (contagem_do_kmer / total_de_kmers_da_classe) em vez de contagem bruta.
   Isso é o que corrige matematicamente o viés de volume: mesmo que o pool
   de TEs tenha mais nucleotídeos totais que o Controle, a pergunta deixa
   de ser "quantas vezes esse k-mer aparece no total?" (que favorece quem
   tem mais dados) e passa a ser "que FRAÇÃO dos k-mers dessa classe é
   este?" — uma medida invariante ao tamanho do pool.
   (Mesma lógica do CPM/TMM usado no artigo da Francismar et al., 2023,
   seção "TE analysis", para normalizar contagens de TE entre bibliotecas
   de tamanhos diferentes.)

2. FILTRO DE SUPORTE MÍNIMO: com k grande (ex: k=15), a maioria dos k-mers
   observados aparece só 1-2 vezes — isso é ruído amostral, não assinatura
   estrutural recorrente. O parâmetro FREQ_MINIMA descarta esses casos.

3. SUPORTE A K=15 e EXPORTAÇÃO PARA EXCEL (uma aba por valor de K).

4. Por padrão, lê o FASTA de TE já balanceado por balancear_volume_te.py
   (K8108_TEs_balanceado.fasta), não o FASTA bruto de 93%.

DEPENDÊNCIA EXTRA
------------------
pip install openpyxl

USO
---
python subtracao_grafos_v2.py
"""

import os
import logging
import gc
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
from collections import defaultdict
from rich.logging import RichHandler
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, MofNCompleteColumn

logging.basicConfig(level="INFO", format="%(message)s", handlers=[RichHandler(show_time=False)])
log = logging.getLogger("rich")


def extrair_kmers_grafo(arquivo_fasta, k, step=1):
    """
    Mesma lógica da v1, mas agora também retorna o TOTAL de k-mers
    extraídos da classe (necessário para normalizar por frequência
    relativa depois).
    """
    contagem = defaultdict(int)
    total_kmers = 0
    try:
        with open(arquivo_fasta, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                seq = str(record.seq).upper()
                if 'N' in seq:
                    continue
                for i in range(0, len(seq) - k + 1, step):
                    kmer = seq[i: i + k]
                    contagem[kmer] += 1
                    total_kmers += 1
    except FileNotFoundError:
        log.error(f"Arquivo não encontrado: {arquivo_fasta}")
    return contagem, total_kmers


def calcular_diferenca_grafos(mapa_tes, total_te, mapa_controle, total_ctrl, freq_minima=5):
    """
    Compara TE vs Controle por FREQUÊNCIA RELATIVA (proporção dentro de
    cada classe), não contagem bruta.
    """
    features = []
    todos_kmers = set(mapa_tes.keys()) | set(mapa_controle.keys())

    for kmer in todos_kmers:
        freq_te_abs = mapa_tes.get(kmer, 0)
        freq_ctrl_abs = mapa_controle.get(kmer, 0)

        # Descarta k-mers com pouquíssimas observações nas duas classes:
        # com suporte tão baixo, "exclusividade" é coincidência estatística,
        # não padrão estrutural recorrente.
        if freq_te_abs < freq_minima and freq_ctrl_abs < freq_minima:
            continue

        freq_te_rel = freq_te_abs / total_te if total_te else 0
        freq_ctrl_rel = freq_ctrl_abs / total_ctrl if total_ctrl else 0
        diff_rel = freq_te_rel - freq_ctrl_rel

        if diff_rel > 0:
            features.append({
                "kmer": kmer,
                "freq_TE_abs": freq_te_abs,
                "freq_Ctrl_abs": freq_ctrl_abs,
                "freq_TE_rel": freq_te_rel,
                "freq_Ctrl_rel": freq_ctrl_rel,
                "diff_relativa": diff_rel,
            })

    return pd.DataFrame(features)


def plotar_top_features(df, k_atual, pasta_destino, top_n=50):
    if df.empty:
        log.warning(f"⚠️  Nenhuma feature exclusiva encontrada para K={k_atual} "
                     f"(considere reduzir FREQ_MINIMA).")
        return

    df_top = df.head(top_n)
    plt.figure(figsize=(10, 8))

    tamanhos = (df_top["diff_relativa"] / df_top["diff_relativa"].max()) * 600 + 50

    # Frequências em ‰ (por mil k-mers) só para deixar os eixos legíveis
    x = df_top["freq_Ctrl_rel"] * 1000
    y = df_top["freq_TE_rel"] * 1000

    scatter = plt.scatter(
        x=x, y=y, s=tamanhos, c=df_top["diff_relativa"],
        cmap="viridis", alpha=0.7, edgecolors="w", linewidth=1
    )

    plt.title(f"Dispersão de K-mers (K={k_atual}) — frequência relativa normalizada",
              fontsize=13, pad=15)
    plt.xlabel("Frequência relativa no Controle (‰ dos k-mers da classe)", fontsize=12)
    plt.ylabel("Frequência relativa nos TEs (‰ dos k-mers da classe)", fontsize=12)

    max_val = max(x.max(), y.max())
    plt.plot([0, max_val], [0, max_val], color='red', linestyle='--',
              alpha=0.5, label='Linha Neutra')

    for _, row in df.head(10).iterrows():
        plt.annotate(
            row['kmer'],
            (row['freq_Ctrl_rel'] * 1000, row['freq_TE_rel'] * 1000),
            xytext=(8, 8), textcoords='offset points',
            fontsize=8, fontweight='bold', color="#333333"
        )

    cbar = plt.colorbar(scatter)
    cbar.set_label('Diferença relativa (exclusividade)', rotation=270, labelpad=15)
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.4)
    plt.tight_layout()

    caminho_grafico = os.path.join(pasta_destino, f"scatter_features_k{k_atual}_normalizado.png")
    plt.savefig(caminho_grafico, dpi=300, bbox_inches='tight')
    plt.close()


def main():
    # Adicione outros valores (ex: [8, 12, 15]) se quiser comparar o
    # efeito de k no mesmo gráfico/planilha.
    K_TAMANHOS = [15]
    PASSO = 1
    FREQ_MINIMA = 5  # suporte mínimo por k-mer; ajuste conforme o tamanho do dataset

    caminho_script = os.path.dirname(os.path.abspath(__file__))
    pasta_dados = os.path.join(caminho_script, "..", "dados_limpos")
    resultados = os.path.join(caminho_script, "..", "resultados")
    os.makedirs(resultados, exist_ok=True)

    # Usa o FASTA já balanceado por balancear_volume_te.py
    arquivo_tes = os.path.join(pasta_dados, "balanceado", "K8108_TEs_balanceado.fasta")
    arquivo_ctrl = os.path.join(pasta_dados, "K8108_genes_nao_TEs.fasta")

    if not os.path.exists(arquivo_tes):
        log.error(f"❌ {arquivo_tes} não encontrado.")
        log.error("   Rode balancear_volume_te.py antes deste script.")
        return

    excel_path = os.path.join(resultados, "features_kmers_normalizado.xlsx")
    with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
        for k in K_TAMANHOS:
            log.info("=" * 55)
            log.info(f"🚀 Iniciando análise para K={k} (Step={PASSO}, freq. mínima={FREQ_MINIMA})")

            with Progress(SpinnerColumn(), TextColumn("[progress.description]{task.description}"),
                          BarColumn(), MofNCompleteColumn()) as prg:
                t1 = prg.add_task("Processando Transposons (balanceado)...", total=1)
                grafo_tes, total_te = extrair_kmers_grafo(arquivo_tes, k, PASSO)
                prg.update(t1, advance=1)

                t2 = prg.add_task("Processando Controle...", total=1)
                grafo_ctrl, total_ctrl = extrair_kmers_grafo(arquivo_ctrl, k, PASSO)
                prg.update(t2, advance=1)

            log.info(f"   Total de k-mers extraídos -> TE: {total_te:,} | Controle: {total_ctrl:,} "
                      f"(razão: {total_te / total_ctrl:.2f}x)" if total_ctrl else "")

            log.info(f"📊 Subtraindo grafos (frequência relativa) para K={k}...")
            df_features = calcular_diferenca_grafos(grafo_tes, total_te, grafo_ctrl, total_ctrl, FREQ_MINIMA)
            df_features = df_features.sort_values(by="diff_relativa", ascending=False)

            sheet_name = f"k{k}"
            df_features.to_excel(writer, sheet_name=sheet_name, index=False)

            log.info("📈 Gerando gráfico de dispersão...")
            plotar_top_features(df_features, k, resultados)

            n_unicos = len(set(grafo_tes) | set(grafo_ctrl))
            log.info(f"✅ K={k} concluído: {len(df_features)} k-mers exclusivos "
                      f"(suporte >= {FREQ_MINIMA}) de {n_unicos:,} k-mers únicos observados.")

            del grafo_tes, grafo_ctrl, df_features
            gc.collect()

    log.info(f"🎉 Excel salvo em: {excel_path}")


if __name__ == "__main__":
    main()
