"""
balancear_volume_te.py

OBJETIVO
--------
Corrigir o viés identificado pelos orientadores: o pool de TEs tem muito
mais nucleotídeos totais que o pool de Controle (já que TEs ~= 93% do
genoma), então qualquer comparação de CONTAGEM ABSOLUTA de k-mers favorece
artificialmente o TE, mesmo sem nenhuma assinatura estrutural real.

ESTRATÉGIA
----------
1. Mede o volume total (pb) do grupo Controle (K8108_genes_nao_TEs.fasta).
2. Define esse volume como o ALVO de pb a serem usados do lado do TE.
3. Divide esse alvo igualmente entre as 3 regiões estruturais do TE
   (cabeça / meio / fim, geradas pelo dividir_tes_regioes.py), preservando
   representação posicional equilibrada.
4. Em cada região, embaralha as sequências (seed fixa = reprodutibilidade)
   e seleciona um SUBCONJUNTO DE FRAGMENTOS INTEIROS (nunca corta uma
   sequência no meio) até atingir o alvo de pb daquela região.
5. Salva um único FASTA balanceado (K8108_TEs_balanceado.fasta) +
   um relatório CSV com os percentuais reais usados de cada região.

PRÉ-REQUISITO
-------------
Rodar antes: dividir_tes_regioes.py (gera regioes_tes/K8108_TEs_{cabeca,meio,fim}.fasta)

USO
---
python balancear_volume_te.py
"""

import os
import random
import logging
from Bio import SeqIO
import pandas as pd
from rich.logging import RichHandler

logging.basicConfig(level="INFO", format="%(message)s", handlers=[RichHandler(show_time=False)])
log = logging.getLogger("rich")

# Seed fixa: documentar este valor no TCC para garantir reprodutibilidade
# do experimento (qualquer pessoa que rode o script com SEED=42 obtém
# exatamente o mesmo FASTA balanceado).
SEED = 42
random.seed(SEED)


def contar_bp_total(arquivo_fasta):
    """Soma o total de pares de base e número de sequências de um FASTA."""
    total_bp = 0
    n_seqs = 0
    for record in SeqIO.parse(arquivo_fasta, "fasta"):
        total_bp += len(record.seq)
        n_seqs += 1
    return total_bp, n_seqs


def subamostrar_regiao(arquivo_regiao, alvo_bp, nome_regiao):
    """
    Embaralha os fragmentos da região e acumula sequências inteiras
    (sem cortar nenhuma no meio) até atingir o volume-alvo em pb.

    Selecionar fragmentos inteiros (em vez de truncar) preserva a
    integridade biológica de cada subsequência cabeça/meio/fim.
    """
    registros = list(SeqIO.parse(arquivo_regiao, "fasta"))
    random.shuffle(registros)

    selecionados = []
    acumulado = 0
    for rec in registros:
        if acumulado >= alvo_bp:
            break
        selecionados.append(rec)
        acumulado += len(rec.seq)

    pct_usado = 100 * len(selecionados) / len(registros) if registros else 0
    log.info(
        f"  {nome_regiao:<8} {len(selecionados):>5}/{len(registros):<5} sequências "
        f"({pct_usado:5.1f}%)  ->  {acumulado:>10,} pb  (alvo: {alvo_bp:,} pb)"
    )
    return selecionados, acumulado, pct_usado


def main():
    caminho_script = os.path.dirname(os.path.abspath(__file__))
    pasta_dados = os.path.join(caminho_script, "..", "dados_limpos")
    pasta_regioes = os.path.join(pasta_dados, "regioes_tes")
    pasta_saida = os.path.join(pasta_dados, "balanceado")
    os.makedirs(pasta_saida, exist_ok=True)

    arquivo_ctrl = os.path.join(pasta_dados, "K8108_genes_nao_TEs.fasta")
    regioes = {
        "cabeca": os.path.join(pasta_regioes, "K8108_TEs_cabeca.fasta"),
        "meio":   os.path.join(pasta_regioes, "K8108_TEs_meio.fasta"),
        "fim":    os.path.join(pasta_regioes, "K8108_TEs_fim.fasta"),
    }

    for nome, caminho in regioes.items():
        if not os.path.exists(caminho):
            log.error(f"❌ Arquivo da região '{nome}' não encontrado: {caminho}")
            log.error("   Rode dividir_tes_regioes.py antes deste script.")
            return

    log.info("📏 Medindo volume total do grupo Controle (referência do balanceamento)...")
    bp_ctrl, n_ctrl = contar_bp_total(arquivo_ctrl)
    log.info(f"   Controle: {n_ctrl} sequências, {bp_ctrl:,} pb totais\n")

    # Alvo: ~1 pb de TE para cada 1 pb de Controle (razão 1:1), dividido
    # igualmente entre as 3 regiões estruturais. Se os orientadores
    # pedirem outra razão (ex: 2:1), basta multiplicar bp_ctrl abaixo.
    RAZAO_TE_CONTROLE = 1.0
    alvo_total = bp_ctrl * RAZAO_TE_CONTROLE
    alvo_por_regiao = int(alvo_total // 3)

    log.info(f"🎯 Alvo de subamostragem: {alvo_por_regiao:,} pb por região "
              f"(total ≈ {int(alvo_total):,} pb, razão TE/Controle = {RAZAO_TE_CONTROLE})\n")

    relatorio = []
    todos_selecionados = []
    for nome, caminho in regioes.items():
        selecionados, acumulado, pct = subamostrar_regiao(caminho, alvo_por_regiao, nome)
        todos_selecionados.extend(selecionados)
        relatorio.append({
            "regiao": nome,
            "pb_alvo": alvo_por_regiao,
            "pb_obtido": acumulado,
            "n_sequencias_usadas": len(selecionados),
            "percentual_do_total_disponivel": round(pct, 2),
        })

    arquivo_te_balanceado = os.path.join(pasta_saida, "K8108_TEs_balanceado.fasta")
    SeqIO.write(todos_selecionados, arquivo_te_balanceado, "fasta")

    bp_te_final = sum(len(r.seq) for r in todos_selecionados)
    log.info(f"\n✅ FASTA balanceado salvo em: {arquivo_te_balanceado}")
    log.info(f"   TE balanceado: {len(todos_selecionados)} sequências, {bp_te_final:,} pb")
    log.info(f"   Controle:      {n_ctrl} sequências, {bp_ctrl:,} pb")
    log.info(f"   Razão final TE/Controle: {bp_te_final / bp_ctrl:.3f}  (alvo era {RAZAO_TE_CONTROLE})")

    df_rel = pd.DataFrame(relatorio)
    caminho_relatorio = os.path.join(pasta_saida, "relatorio_balanceamento.csv")
    df_rel.to_csv(caminho_relatorio, index=False)
    log.info(f"   Relatório salvo em: {caminho_relatorio}")
    log.info(f"   Seed usada: {SEED}  (citar no TCC para reprodutibilidade)")
    log.info("\n🚀 PRÓXIMO PASSO: rodar subtracao_grafos_v2.py apontando para "
              "K8108_TEs_balanceado.fasta")


if __name__ == "__main__":
    main()
