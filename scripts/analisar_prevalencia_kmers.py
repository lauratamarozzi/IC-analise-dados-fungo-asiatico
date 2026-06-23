"""
analisar_prevalencia_kmers.py

PERGUNTA QUE ESSE SCRIPT RESPONDE
----------------------------------
A contagem bruta de um k-mer (ex: 658 ocorrências de "TTTTATGTGCCCCAA")
NÃO diz se isso vem de:

  (a) UMA única cópia de TE com uma repetição em tandem gigante desse
      motivo (ex: um microssatélite expandido em 1 elemento só) — nesse
      caso é uma curiosidade local, não uma "regra" da família; ou

  (b) MUITAS cópias diferentes de TE que compartilham esse motivo uma
      vez cada — esse sim é candidato a assinatura estrutural recorrente
      e biologicamente interessante.

Os dois cenários produzem o MESMO total de ocorrências na sua planilha,
mas têm significados completamente diferentes. Este script calcula, para
uma lista de k-mers de interesse:

  - n_sequencias_distintas : em quantas sequências de TE o k-mer aparece
                              pelo menos 1 vez (prevalência na família)
  - ocorrencia_maxima_em_1_sequencia : a maior contagem do k-mer DENTRO
                              de uma única sequência (alto = expansão
                              local em tandem)
  - concentracao_pct       : que fração do total de ocorrências vem dessa
                              única sequência mais "carregada"

USO
---
Edite KMERS_DE_INTERESSE abaixo (cole os k-mers que quer investigar — por
exemplo, os 15 mais frequentes da sua planilha) e rode o script.
"""

import os
import logging
from collections import Counter
from Bio import SeqIO
import pandas as pd
from rich.logging import RichHandler

logging.basicConfig(level="INFO", format="%(message)s", handlers=[RichHandler(show_time=False)])
log = logging.getLogger("rich")

# Cole aqui os k-mers que quer investigar
KMERS_DE_INTERESSE = [
    "AAAAAAAAAAAAAAA",
    "TTTTTTTTTTTTTTT",
    "ATATATATATATATA",
    "TATATATATATATAT",
    "AAAAAAATAAAAAAA",
    "AAAAAAAATAAAAAA",
    "TTTTATGTGCCCCAA",
    "TTTATGTGCCCCAAA",
    "AAAAAATAAAAAAAA",
    "TATGTGCCCCAAAGA",
    "TTATGTGCCCCAAAG",
    "CTTTTATGTGCCCCA",
    "AAAAATAAAAAAAAA",
    "CCATCAGAAAAAAAA",
    "TAAAAAAAAAAAAAA",
]


def diversidade_simbolica(kmer):
    """
    Quantos nucleotídeos distintos compõem o k-mer.
    1-2 = repetição simples (homopolímero/microssatélite, ex: poly-A, AT).
    3-4 = motivo mais específico, candidato mais forte a "regra estrutural".
    """
    return len(set(kmer))


def contar_ocorrencias_sobrepostas(seq, kmer):
    """
    Conta ocorrências PERMITINDO sobreposição (mesma convenção de
    step=1 usada em extrair_kmers_grafo). str.count() do Python NÃO
    permite sobreposição e subestimaria k-mers repetitivos como
    homopolímeros — por isso esse contador manual é necessário aqui.
    """
    count = 0
    start = 0
    while True:
        idx = seq.find(kmer, start)
        if idx == -1:
            break
        count += 1
        start = idx + 1
    return count


def analisar(arquivo_fasta, kmers):
    distribuicao = {k: Counter() for k in kmers}  # kmer -> {id_seq: contagem}

    for record in SeqIO.parse(arquivo_fasta, "fasta"):
        seq = str(record.seq).upper()
        for kmer in kmers:
            n = contar_ocorrencias_sobrepostas(seq, kmer)
            if n > 0:
                distribuicao[kmer][record.id] = n

    linhas = []
    for kmer in kmers:
        distrib = distribuicao[kmer]
        ocorrencias_totais = sum(distrib.values())
        n_seqs = len(distrib)
        max_em_uma_seq = max(distrib.values()) if distrib else 0
        concentracao = round(100 * max_em_uma_seq / ocorrencias_totais, 1) if ocorrencias_totais else 0

        linhas.append({
            "kmer": kmer,
            "diversidade_simbolica": diversidade_simbolica(kmer),
            "ocorrencias_totais": ocorrencias_totais,
            "n_sequencias_distintas": n_seqs,
            "ocorrencia_maxima_em_1_sequencia": max_em_uma_seq,
            "concentracao_pct": concentracao,
        })

    return pd.DataFrame(linhas)


def classificar_kmer(row):
    """
    Classifica automaticamente cada k-mer com base nos dados calculados.
    Retorna uma string de diagnóstico para o relatório.
    """
    if row["diversidade_simbolica"] <= 2:
        tipo = "Homopolímero/Microssatélite (ex: poly-A, AT-repeat)"
    else:
        tipo = "Motivo composto (sequência específica)"

    if row["concentracao_pct"] >= 70:
        distribuicao = "EXPANSÃO LOCAL — concentrado em poucas cópias"
        interpretacao = "Curiosidade pontual. Não generalizar como regra da família."
    elif row["n_sequencias_distintas"] >= 50:
        distribuicao = "ASSINATURA RECORRENTE — espalhado por muitas cópias"
        interpretacao = "Candidato forte a assinatura estrutural. Vale BLASTar."
    else:
        distribuicao = "INTERMEDIÁRIO — distribuição moderada"
        interpretacao = "Investigar em mais detalhe antes de concluir."

    return tipo, distribuicao, interpretacao


def salvar_relatorio_txt(df, pasta_relatorios, arquivo_fasta_usado):
    """Gera um relatório legível em .txt na pasta relatorios."""
    os.makedirs(pasta_relatorios, exist_ok=True)
    caminho = os.path.join(pasta_relatorios, "relatorio_prevalencia_kmers.txt")

    separador = "=" * 65

    with open(caminho, "w", encoding="utf-8") as f:
        f.write(separador + "\n")
        f.write("RELATÓRIO DE PREVALÊNCIA DE K-MERS DE INTERESSE\n")
        f.write("Projeto: Caracterização de TEs - P. pachyrhizi (K8108)\n")
        f.write("Autora:  Laura Tamarozzi\n")
        f.write(separador + "\n\n")
        f.write(f"Arquivo analisado : {arquivo_fasta_usado}\n")
        f.write(f"K-mers investigados: {len(df)}\n\n")

        f.write(separador + "\n")
        f.write("LEGENDA DAS COLUNAS\n")
        f.write(separador + "\n")
        f.write("  diversidade_simbolica          : quantos nucleotídeos distintos\n")
        f.write("                                   compõem o k-mer (1-4).\n")
        f.write("                                   1-2 = homopolímero/microssatélite.\n")
        f.write("                                   3-4 = motivo específico.\n")
        f.write("  ocorrencias_totais             : total de vezes que o k-mer\n")
        f.write("                                   aparece em todo o FASTA de TEs.\n")
        f.write("  n_sequencias_distintas         : em quantas cópias de TE diferentes\n")
        f.write("                                   o k-mer aparece ao menos 1 vez.\n")
        f.write("                                   (alto = assinatura recorrente)\n")
        f.write("  ocorrencia_maxima_em_1_sequencia: a maior contagem do k-mer dentro\n")
        f.write("                                   de uma ÚNICA sequência.\n")
        f.write("                                   (alto = expansão local em tandem)\n")
        f.write("  concentracao_pct               : % das ocorrências totais que vêm\n")
        f.write("                                   dessa única sequência mais carregada.\n")
        f.write("                                   (>70% = provavelmente ruído local)\n\n")

        f.write(separador + "\n")
        f.write("ANÁLISE INDIVIDUAL\n")
        f.write(separador + "\n\n")

        for _, row in df.iterrows():
            tipo, distribuicao, interpretacao = classificar_kmer(row)
            f.write(f"K-mer              : {row['kmer']}\n")
            f.write(f"Tipo               : {tipo}\n")
            f.write(f"Diversidade        : {row['diversidade_simbolica']} nucleotídeo(s) distinto(s)\n")
            f.write(f"Ocorrências totais : {row['ocorrencias_totais']}\n")
            f.write(f"Sequências distintas: {row['n_sequencias_distintas']}\n")
            f.write(f"Máx. em 1 sequência: {row['ocorrencia_maxima_em_1_sequencia']}\n")
            f.write(f"Concentração       : {row['concentracao_pct']}%\n")
            f.write(f"Distribuição       : {distribuicao}\n")
            f.write(f"Interpretação      : {interpretacao}\n")
            f.write("-" * 65 + "\n\n")

        f.write(separador + "\n")
        f.write("RESUMO\n")
        f.write(separador + "\n\n")

        n_assinatura = sum(
            1 for _, r in df.iterrows()
            if r["concentracao_pct"] < 70 and r["n_sequencias_distintas"] >= 50
        )
        n_local = sum(1 for _, r in df.iterrows() if r["concentracao_pct"] >= 70)
        n_intermediario = len(df) - n_assinatura - n_local

        f.write(f"  Candidatos a assinatura estrutural recorrente : {n_assinatura}\n")
        f.write(f"  Prováveis expansões locais (ruído/tandem)     : {n_local}\n")
        f.write(f"  Casos intermediários (investigar mais)        : {n_intermediario}\n\n")
        f.write("  Os candidatos a assinatura recorrente são os mais indicados\n")
        f.write("  para consulta em BLAST/Repbase e inclusão como features\n")
        f.write("  prioritárias no treinamento da Árvore de Decisão.\n\n")
        f.write(separador + "\n")

    return caminho


def main():
    caminho_script = os.path.dirname(os.path.abspath(__file__))
    pasta_dados = os.path.join(caminho_script, "..", "dados_limpos")

    arquivo_tes = os.path.join(pasta_dados, "balanceado", "K8108_TEs_balanceado.fasta")
    if not os.path.exists(arquivo_tes):
        log.warning("FASTA balanceado não encontrado, usando o FASTA original de TEs.")
        arquivo_tes = os.path.join(pasta_dados, "K8108_TEs_limpos.fasta")

    log.info(f"📂 Analisando prevalência em: {arquivo_tes}")
    df = analisar(arquivo_tes, KMERS_DE_INTERESSE)
    df = df.sort_values("n_sequencias_distintas", ascending=False)

    # Salvar CSV na pasta resultados
    resultados = os.path.join(caminho_script, "..", "resultados")
    os.makedirs(resultados, exist_ok=True)
    saida_csv = os.path.join(resultados, "prevalencia_kmers_interesse.csv")
    df.to_csv(saida_csv, index=False)

    # Salvar relatório .txt na pasta relatorios
    pasta_relatorios = os.path.join(caminho_script, "..", "relatorios")
    saida_txt = salvar_relatorio_txt(df, pasta_relatorios, arquivo_tes)

    log.info(f"\n{df.to_string(index=False)}")
    log.info(f"\n✅ CSV salvo em      : {saida_csv}")
    log.info(f"✅ Relatório salvo em: {saida_txt}")


if __name__ == "__main__":
    main()
