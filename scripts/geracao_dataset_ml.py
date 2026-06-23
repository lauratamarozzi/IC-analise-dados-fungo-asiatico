import os
import pandas as pd
from Bio import SeqIO
from collections import Counter
import logging
from rich.logging import RichHandler
from rich.progress import track

# Configuração de Logging
logging.basicConfig(level="INFO", format="%(message)s", handlers=[RichHandler(show_time=False)])
log = logging.getLogger("rich")

def processar_sequencias_individuais(arquivo_fasta, k, step, rotulo_classe):
    """
    Lê o FASTA e conta os K-mers separadamente para CADA sequência.
    Retorna uma lista de dicionários (uma linha por sequência).
    """
    dados_matriz = []
    
    try:
        registros = list(SeqIO.parse(arquivo_fasta, "fasta"))
        log.info(f"Extraindo K-mers (K={k}) de {len(registros)} sequências em: {os.path.basename(arquivo_fasta)}")
        
        for record in track(registros, description=f"Processando {rotulo_classe}..."):
            seq = str(record.seq).upper()
            if 'N' in seq: continue  # Ignora sequências com nucleotídeos não identificados
            
            # Conta os k-mers apenas desta sequência usando o passo 1
            kmers_desta_seq = [seq[i : i + k] for i in range(0, len(seq) - k + 1, step)]
            contagem = Counter(kmers_desta_seq)
            
            # Cria a "linha" da tabela com o ID, o Rótulo e as contagens
            linha = {
                "ID_Sequencia": record.id,
                "Classe": rotulo_classe  # "TE" ou "Controle"
            }
            linha.update(contagem) # Adiciona os k-mers como colunas
            
            dados_matriz.append(linha)
            
    except FileNotFoundError:
        log.error(f"Arquivo não encontrado: {arquivo_fasta}")
        
    return dados_matriz

def main():
    K_TAMANHOS = [4, 6] 
    PASSO = 1      
    
    caminho_script = os.path.dirname(os.path.abspath(__file__))
    pasta_dados = os.path.join(caminho_script, "..", "dados_limpos")
    pasta_resultados = os.path.join(caminho_script, "..", "resultados")
    
    if not os.path.exists(pasta_resultados):
        os.makedirs(pasta_resultados)
    
    arquivo_tes = os.path.join(pasta_dados, "K8108_TEs_limpos.fasta")
    arquivo_ctrl = os.path.join(pasta_dados, "K8108_genes_nao_TEs.fasta")
    
    for k in K_TAMANHOS:
        log.info(f"\n==================================================")
        log.info(f"🤖 CONSTRUINDO DATASET DE MACHINE LEARNING | K={k}")
        
        # 1. Processa os TEs (Classe 1)
        dados_tes = processar_sequencias_individuais(arquivo_tes, k, PASSO, rotulo_classe="TE")
        
        # 2. Processa o Controle (Classe 0)
        dados_ctrl = processar_sequencias_individuais(arquivo_ctrl, k, PASSO, rotulo_classe="Controle")
        
        # 3. Junta tudo em um único DataFrame do Pandas
        log.info("Unindo dados e preenchendo K-mers ausentes com 0...")
        df_completo = pd.DataFrame(dados_tes + dados_ctrl)
        
        # Preenche os NaN (k-mers que não apareceram em determinada sequência) com zero
        df_completo = df_completo.fillna(0)
        
        # Move a coluna 'Classe' para o final para ficar mais organizado
        colunas = [c for c in df_completo.columns if c not in ['ID_Sequencia', 'Classe']]
        df_completo = df_completo[['ID_Sequencia'] + colunas + ['Classe']]
        
        # 4. Salva a Matriz Pronta para o Machine Learning
        output_csv = os.path.join(pasta_resultados, f"dataset_ml_k{k}.csv")
        df_completo.to_csv(output_csv, index=False)
        
        linhas, colunas_totais = df_completo.shape
        log.info(f"✅ Dataset K={k} salvo! Matriz final: {linhas} sequências x {colunas_totais - 2} K-mers (features).")

if __name__ == "__main__":
    main()