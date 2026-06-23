import os
import logging
import gc
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from collections import defaultdict
import pandas as pd
from rich.logging import RichHandler
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, MofNCompleteColumn

# Configuração de Logging
logging.basicConfig(level="INFO", format="%(message)s", handlers=[RichHandler(show_time=False)])
log = logging.getLogger("rich")

def dividir_sequencia_em_tres(record):
    """
    Divide uma sequência em 3 partes iguais.
    
    Exemplo:
        Sequência de 9.000bp:
        Cabeça: posições 0    → 3.000bp
        Meio:   posições 3000 → 6.000bp
        Fim:    posições 6000 → 9.000bp
    
    Para sequências com tamanho não divisível por 3, o 'meio' absorve
    os nucleotídeos extras (comportamento do Python com divisão inteira).
    """
    seq_str = str(record.seq).upper()
    tamanho = len(seq_str)
    
    # Calcular pontos de corte
    corte1 = tamanho // 3
    corte2 = 2 * (tamanho // 3)
    
    cabeca = seq_str[:corte1]
    meio   = seq_str[corte1:corte2]
    fim    = seq_str[corte2:]
    
    # Criar novos records com IDs descritivos
    record_cabeca = SeqRecord(
        Seq(cabeca),
        id=f"{record.id}_cabeca",
        description=f"regiao=cabeca tamanho={len(cabeca)}bp original={tamanho}bp"
    )
    record_meio = SeqRecord(
        Seq(meio),
        id=f"{record.id}_meio",
        description=f"regiao=meio tamanho={len(meio)}bp original={tamanho}bp"
    )
    record_fim = SeqRecord(
        Seq(fim),
        id=f"{record.id}_fim",
        description=f"regiao=fim tamanho={len(fim)}bp original={tamanho}bp"
    )
    
    return record_cabeca, record_meio, record_fim

def processar_e_dividir_tes(arquivo_tes, pasta_saida, tamanho_minimo=300):
    """
    Lê o arquivo de TEs e divide cada sequência em 3 regiões.
    
    Filtra sequências muito pequenas (< tamanho_minimo bp) pois
    dividir em 3 geraria fragmentos biologicamente inúteis.
    """
    
    log.info(f"📂 Lendo arquivo: {arquivo_tes}")
    
    lista_cabecas = []
    lista_meios   = []
    lista_fins     = []
    
    # Estatísticas
    stats = {
        'total_lidas': 0,
        'descartadas_tamanho': 0,
        'processadas': 0,
        'tamanhos_originais': [],
        'tamanhos_cabeca': [],
        'tamanhos_meio':   [],
        'tamanhos_fim':    []
    }
    
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        MofNCompleteColumn()
    ) as prg:
        
        # Primeiro, contar total para a barra de progresso
        total = sum(1 for _ in SeqIO.parse(arquivo_tes, "fasta"))
        task = prg.add_task("Dividindo sequências em 3 regiões...", total=total)
        
        for record in SeqIO.parse(arquivo_tes, "fasta"):
            stats['total_lidas'] += 1
            tamanho = len(record.seq)
            
            # Filtrar sequências muito pequenas
            if tamanho < tamanho_minimo:
                stats['descartadas_tamanho'] += 1
                prg.update(task, advance=1)
                continue
            
            # Dividir em 3 regiões
            cabeca, meio, fim = dividir_sequencia_em_tres(record)
            
            lista_cabecas.append(cabeca)
            lista_meios.append(meio)
            lista_fins.append(fim)
            
            # Coletar estatísticas
            stats['tamanhos_originais'].append(tamanho)
            stats['tamanhos_cabeca'].append(len(cabeca.seq))
            stats['tamanhos_meio'].append(len(meio.seq))
            stats['tamanhos_fim'].append(len(fim.seq))
            stats['processadas'] += 1
            
            prg.update(task, advance=1)
    
    # Salvar os 3 arquivos FASTA
    log.info("💾 Salvando arquivos FASTA das 3 regiões...")
    
    arquivo_cabeca = os.path.join(pasta_saida, "K8108_TEs_cabeca.fasta")
    arquivo_meio   = os.path.join(pasta_saida, "K8108_TEs_meio.fasta")
    arquivo_fim    = os.path.join(pasta_saida, "K8108_TEs_fim.fasta")
    
    SeqIO.write(lista_cabecas, arquivo_cabeca, "fasta")
    SeqIO.write(lista_meios,   arquivo_meio,   "fasta")
    SeqIO.write(lista_fins,    arquivo_fim,    "fasta")
    
    return stats, arquivo_cabeca, arquivo_meio, arquivo_fim

def gerar_relatorio_estatistico(stats, pasta_saida):
    """Gera relatório com estatísticas das regiões geradas"""
    
    def media(lista):
        return sum(lista) / len(lista) if lista else 0
    
    log.info("\n" + "="*55)
    log.info("📊 RELATÓRIO DE DIVISÃO DAS SEQUÊNCIAS")
    log.info("="*55)
    log.info(f"Total de sequências lidas:        {stats['total_lidas']}")
    log.info(f"Descartadas (muito pequenas):     {stats['descartadas_tamanho']}")
    log.info(f"Sequências processadas:           {stats['processadas']}")
    log.info("")
    log.info(f"{'Região':<12} {'Qtd':>6} {'Média (bp)':>12} {'Min (bp)':>10} {'Max (bp)':>10}")
    log.info("-"*55)
    log.info(f"{'Original':<12} {stats['processadas']:>6} {media(stats['tamanhos_originais']):>12.1f} {min(stats['tamanhos_originais']):>10} {max(stats['tamanhos_originais']):>10}")
    log.info(f"{'Cabeça':<12} {len(stats['tamanhos_cabeca']):>6} {media(stats['tamanhos_cabeca']):>12.1f} {min(stats['tamanhos_cabeca']):>10} {max(stats['tamanhos_cabeca']):>10}")
    log.info(f"{'Meio':<12} {len(stats['tamanhos_meio']):>6} {media(stats['tamanhos_meio']):>12.1f} {min(stats['tamanhos_meio']):>10} {max(stats['tamanhos_meio']):>10}")
    log.info(f"{'Fim':<12} {len(stats['tamanhos_fim']):>6} {media(stats['tamanhos_fim']):>12.1f} {min(stats['tamanhos_fim']):>10} {max(stats['tamanhos_fim']):>10}")
    log.info("="*55)
    
    # Salvar relatório em arquivo
    relatorio_path = os.path.join(pasta_saida, "relatorio_divisao_regioes.txt")
    with open(relatorio_path, 'w', encoding='utf-8') as f:
        f.write("RELATÓRIO DE DIVISÃO DE TEs EM REGIÕES\n")
        f.write("Projeto: GRAMEP - P. pachyrhizi XAI\n")
        f.write("Autora: Laura Tamarozzi\n")
        f.write("="*55 + "\n\n")
        f.write(f"Total lidas:         {stats['total_lidas']}\n")
        f.write(f"Descartadas:         {stats['descartadas_tamanho']}\n")
        f.write(f"Processadas:         {stats['processadas']}\n\n")
        f.write(f"Média original:      {media(stats['tamanhos_originais']):.1f} bp\n")
        f.write(f"Média cabeça:        {media(stats['tamanhos_cabeca']):.1f} bp\n")
        f.write(f"Média meio:          {media(stats['tamanhos_meio']):.1f} bp\n")
        f.write(f"Média fim:           {media(stats['tamanhos_fim']):.1f} bp\n")
    
    log.info(f"📋 Relatório salvo em: {relatorio_path}")

def main():
    
    log.info("🧬 DIVISÃO DE TEs EM REGIÕES - P. PACHYRHIZI")
    log.info("🔬 Metodologia: Cabeça / Meio / Fim (1/3 cada)")
    log.info("👩‍🔬 Laura Tamarozzi")
    log.info("="*55)
    
    # Caminhos
    caminho_script = os.path.dirname(os.path.abspath(__file__))
    pasta_dados    = os.path.join(caminho_script, "..", "dados_limpos")
    pasta_saida    = os.path.join(caminho_script, "..", "dados_limpos", "regioes_tes")
    
    # Criar pasta de saída
    if not os.path.exists(pasta_saida):
        os.makedirs(pasta_saida)
        log.info(f"📁 Pasta criada: {pasta_saida}")
    
    # Arquivo de TEs limpos (resultado da limpeza anterior)
    arquivo_tes = os.path.join(pasta_dados, "K8108_TEs_limpos.fasta")
    
    # Verificar existência
    if not os.path.exists(arquivo_tes):
        log.error(f"❌ Arquivo não encontrado: {arquivo_tes}")
        log.error("Verifique se a limpeza dos dados foi executada primeiro.")
        return
    
    # Processar e dividir
    stats, arq_cabeca, arq_meio, arq_fim = processar_e_dividir_tes(
        arquivo_tes,
        pasta_saida,
        tamanho_minimo=300  # TEs menores que 300bp serão descartados
    )
    
    # Gerar relatório
    gerar_relatorio_estatistico(stats, pasta_saida)
    
    # Confirmação final
    log.info("\n✅ DIVISÃO CONCLUÍDA!")
    log.info("📁 ARQUIVOS GERADOS:")
    log.info(f"   🔴 Cabeça (5'): {arq_cabeca}")
    log.info(f"   🟡 Meio (central): {arq_meio}")
    log.info(f"   🔵 Fim (3'): {arq_fim}")
    log.info("\n🚀 PRÓXIMO PASSO:")
    log.info("   Rodar o subtracao_grafos.py apontando para cada região")
    log.info("   e comparar individualmente com os não-TEs.")

if __name__ == "__main__":
    main()
