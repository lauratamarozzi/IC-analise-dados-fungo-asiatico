import os
import logging
import gc
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
from collections import defaultdict
from rich.logging import RichHandler
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, MofNCompleteColumn

# Configuração de Logging
logging.basicConfig(level="INFO", format="%(message)s", handlers=[RichHandler(show_time=False)])
log = logging.getLogger("rich")

def extrair_kmers_grafo(arquivo_fasta, k, step=1):
    """
    Transforma a sequência em um 'grafo' de frequências de transição.
    Passo 1 captura todas as sobreposições possíveis (arestas).
    """
    contagem = defaultdict(int)
    try:
        with open(arquivo_fasta, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                seq = str(record.seq).upper()
                if 'N' in seq: continue
                
                for i in range(0, len(seq) - k + 1, step):
                    kmer = seq[i : i + k]
                    contagem[kmer] += 1
    except FileNotFoundError:
        log.error(f"Arquivo não encontrado: {arquivo_fasta}")
    return contagem

def calcular_diferenca_grafos(mapa_tes, mapa_controle):
    features_exclusivas = []
    todos_kmers = set(mapa_tes.keys()) | set(mapa_controle.keys())
    
    for kmer in todos_kmers:
        freq_te = mapa_tes.get(kmer, 0)
        freq_ctrl = mapa_controle.get(kmer, 0)
        
        diff = freq_te - freq_ctrl
        
        if diff > 0:
            features_exclusivas.append({
                "kmer": kmer,
                "freq_TE": freq_te,
                "freq_Ctrl": freq_ctrl,
                "diff_exclusiva": diff
            })
            
    return pd.DataFrame(features_exclusivas)

def plotar_top_features(df, k_atual, pasta_destino, top_n=50):
    df_top = df.head(top_n)
    plt.figure(figsize=(10, 8))
    
    tamanhos = (df_top["diff_exclusiva"] / df_top["diff_exclusiva"].max()) * 600 + 50
    
    scatter = plt.scatter(
        x=df_top["freq_Ctrl"], 
        y=df_top["freq_TE"], 
        s=tamanhos, 
        c=df_top["diff_exclusiva"], 
        cmap="viridis",
        alpha=0.7,
        edgecolors="w",
        linewidth=1
    )
    
    plt.title(f"Dispersão de K-mers (K={k_atual}): Transposons vs. Controle", fontsize=14, pad=15)
    plt.xlabel("Frequência no Controle (Não-TEs)", fontsize=12)
    plt.ylabel("Frequência nos Transposons (TEs)", fontsize=12)
    
    max_val = max(df_top["freq_Ctrl"].max(), df_top["freq_TE"].max())
    plt.plot([0, max_val], [0, max_val], color='red', linestyle='--', alpha=0.5, label='Linha Neutra')
    
    for i, row in df.head(10).iterrows():
        plt.annotate(
            row['kmer'],
            (row['freq_Ctrl'], row['freq_TE']),
            xytext=(8, 8),
            textcoords='offset points',
            fontsize=9,
            fontweight='bold',
            color="#333333"
        )
        
    cbar = plt.colorbar(scatter)
    cbar.set_label('Magnitude de Exclusividade', rotation=270, labelpad=15)
    
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.4)
    plt.tight_layout()
    
    caminho_grafico = os.path.join(pasta_destino, f"scatter_features_k{k_atual}.png")
    plt.savefig(caminho_grafico, dpi=300, bbox_inches='tight')
    plt.close()

def main():
    # LISTA DE Ks PARA TESTAR (Loop) E NOVO PASSO
    K_TAMANHOS = [15] 
    PASSO = 1
    
    caminho_script = os.path.dirname(os.path.abspath(__file__))
    pasta_dados = os.path.join(caminho_script, "..", "dados_limpos")
    resultados = os.path.join(caminho_script, "..", "resultados")
    
    if not os.path.exists(resultados):
        os.makedirs(resultados)
    
    arquivo_tes = os.path.join(pasta_dados, "K8108_TEs_limpos.fasta")
    arquivo_ctrl = os.path.join(pasta_dados, "K8108_genes_nao_TEs.fasta")
    
    # Executa o pipeline para cada tamanho de K na lista
    for k in K_TAMANHOS:
        log.info(f"==================================================")
        log.info(f"🚀 Iniciando análise para K={k} (Step={PASSO})")
        
        with Progress(SpinnerColumn(), TextColumn("[progress.description]{task.description}"), BarColumn(), MofNCompleteColumn()) as prg:
            t1 = prg.add_task("Processando Transposons...", total=1)
            grafo_tes = extrair_kmers_grafo(arquivo_tes, k, PASSO)
            prg.update(t1, advance=1)
            
            t2 = prg.add_task("Processando Controle...", total=1)
            grafo_ctrl = extrair_kmers_grafo(arquivo_ctrl, k, PASSO)
            prg.update(t2, advance=1)

        log.info(f"📊 Subtraindo grafos para K={k}...")
        df_features = calcular_diferenca_grafos(grafo_tes, grafo_ctrl)
        df_features = df_features.sort_values(by="diff_exclusiva", ascending=False)
        
        output_csv = os.path.join(resultados, f"features_importancia_k{k}.csv")
        df_features.to_csv(output_csv, index=False)
        
        log.info(f"📈 Gerando gráfico de dispersão...")
        plotar_top_features(df_features, k, resultados)
        
        log.info(f"✅ Análise do K={k} concluída. Foram analisadas {len(df_features)} transições.")
        
        # Limpeza forçada de memória RAM para não travar o PC do laboratório
        del grafo_tes
        del grafo_ctrl
        del df_features
        gc.collect()

    log.info("🎉 Todos os Ks foram processados com sucesso!")

if __name__ == "__main__":
    main()