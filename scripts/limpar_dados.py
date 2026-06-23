#!/usr/bin/env python3
"""
Script de limpeza para os dados do projeto GRAMEP
Baseado na análise real dos headers do Phakopsora pachyrhizi
Autora: Laura Tamarozzi
Orientadores: Daniel Kaster, Fabrício Lopes, Matheus Pimenta
"""

from Bio import SeqIO
import os
from datetime import datetime
from collections import Counter

def criar_pastas():
    """Cria as pastas necessárias se não existirem"""
    pastas = ['dados_limpos', 'analises', 'relatorios']
    for pasta in pastas:
        if not os.path.exists(pasta):
            os.makedirs(pasta)
            print(f"📁 Pasta criada: {pasta}")

def limpar_tes(arquivo_entrada, arquivo_saida_limpo, arquivo_saida_removidos):
    """
    Limpeza baseada no padrão real dos dados P. pachyrhizi
    
    CRITÉRIOS ESPECÍFICOS PARA SEUS DADOS:
    - Manter: RLX, DTX, RIX, DHX, DMX, DXX, RSX (famílias conhecidas)
    - Remover: noCat (não categorizados)
    - Manter preferencialmente: comp (completos)
    - Avaliar: incomp (incompletos) - manter só se >= 500bp
    """
    
    print(f"\n🧹 LIMPEZA DOS TEs - P. PACHYRHIZI")
    print("=" * 55)
    
    # Famílias aceitas baseadas na sua análise
    familias_aceitas = ['rlx', 'dtx', 'rix', 'dhx', 'dmx', 'dxx', 'rsx', 'rxx']
    
    # Categorias a remover
    categorias_rejeitar = ['nocat']  # Não categorizados
    
    sequencias_limpas = []
    sequencias_removidas = []
    
    estatisticas = {
        'total': 0,
        'mantidas': 0,
        'removidas': 0,
        'familias_mantidas': Counter(),
        'motivos_remocao': Counter()
    }
    
    print("🔍 Processando sequências...")
    
    for record in SeqIO.parse(arquivo_entrada, "fasta"):
        estatisticas['total'] += 1
        header = record.description.lower()
        tamanho = len(record.seq)
        
        # Extrair família (primeira parte do header)
        familia = header.split('-')[0] if '-' in header else header.split('_')[0]
        
        # Verificar se é uma família aceita
        familia_aceita = any(fam in familia for fam in familias_aceitas)
        
        # Verificar se é categoria rejeitada
        categoria_rejeitada = any(cat in header for cat in categorias_rejeitar)
        
        # Verificar completude para elementos incompletos
        is_completo = 'comp' in header
        is_incompleto = 'incomp' in header
        
        # Tamanho mínimo para elementos incompletos
        tamanho_adequado = True
        if is_incompleto:
            tamanho_adequado = tamanho >= 500  # Mínimo 500bp para incompletos
        else:
            tamanho_adequado = tamanho >= 100   # Mínimo 100bp para outros
        
        # Decisão final
        manter = (familia_aceita and 
                 not categoria_rejeitada and 
                 tamanho_adequado)
        
        if manter:
            sequencias_limpas.append(record)
            estatisticas['mantidas'] += 1
            estatisticas['familias_mantidas'][familia] += 1
        else:
            sequencias_removidas.append(record)
            estatisticas['removidas'] += 1
            
            # Documentar motivo da remoção
            if not familia_aceita and not categoria_rejeitada:
                estatisticas['motivos_remocao']['familia_nao_reconhecida'] += 1
            elif categoria_rejeitada:
                estatisticas['motivos_remocao']['nao_categorizado'] += 1
            elif not tamanho_adequado:
                estatisticas['motivos_remocao']['muito_pequeno'] += 1
    
    # Salvar arquivos
    SeqIO.write(sequencias_limpas, arquivo_saida_limpo, "fasta")
    SeqIO.write(sequencias_removidas, arquivo_saida_removidos, "fasta")
    
    # Relatório detalhado
    print(f"\n📊 RELATÓRIO DE LIMPEZA:")
    print(f"   Total analisado: {estatisticas['total']}")
    print(f"   ✅ Mantidas: {estatisticas['mantidas']} ({estatisticas['mantidas']/estatisticas['total']*100:.1f}%)")
    print(f"   ❌ Removidas: {estatisticas['removidas']} ({estatisticas['removidas']/estatisticas['total']*100:.1f}%)")
    
    print(f"\n🏷️ FAMÍLIAS MANTIDAS:")
    for familia, count in estatisticas['familias_mantidas'].most_common():
        percentual = count / estatisticas['mantidas'] * 100 if estatisticas['mantidas'] > 0 else 0
        print(f"   {familia.upper()}: {count} ({percentual:.1f}%)")
    
    print(f"\n🎯 MOTIVOS DE REMOÇÃO:")
    for motivo, count in estatisticas['motivos_remocao'].items():
        print(f"   {motivo}: {count}")
    
    print(f"\n💾 ARQUIVOS SALVOS:")
    print(f"   ✅ TEs limpos: {arquivo_saida_limpo}")
    print(f"   🗑️ TEs removidos: {arquivo_saida_removidos}")
    
    return estatisticas

def preparar_dados_negativos(arquivo_cds, arquivo_saida, max_sequencias=None):
    """Prepara genes (não-TEs) para usar como dados negativos"""
    
    print(f"\n🧬 PREPARANDO DADOS NEGATIVOS (GENES)")
    print("=" * 40)
    
    sequencias_genes = []
    tamanhos = []
    
    for record in SeqIO.parse(arquivo_cds, "fasta"):
        # Filtrar por tamanho razoável
        if 300 <= len(record.seq) <= 10000:
            sequencias_genes.append(record)
            tamanhos.append(len(record.seq))
    
    # Limitar quantidade se especificado
    if max_sequencias and len(sequencias_genes) > max_sequencias:
        import random
        sequencias_genes = random.sample(sequencias_genes, max_sequencias)
        print(f"🎲 Selecionados aleatoriamente: {max_sequencias}")
    
    # Salvar
    SeqIO.write(sequencias_genes, arquivo_saida, "fasta")
    
    print(f"📊 ESTATÍSTICAS DOS GENES:")
    print(f"   Total selecionado: {len(sequencias_genes)}")
    if tamanhos:
        print(f"   Tamanho médio: {sum(tamanhos)/len(tamanhos):.1f} bp")
        print(f"   Faixa de tamanho: {min(tamanhos)} - {max(tamanhos)} bp")
    print(f"   Arquivo salvo: {arquivo_saida}")

def validar_limpeza(arquivo_limpo):
    """Valida o arquivo limpo usando os padrões específicos do P. pachyrhizi"""
    
    print(f"\n🔬 VALIDAÇÃO DO ARQUIVO LIMPO")
    print("=" * 45)
    
    if not os.path.exists(arquivo_limpo):
        print(f"❌ Arquivo não encontrado: {arquivo_limpo}")
        return
    
    familias_encontradas = []
    tamanhos = []
    completos = 0
    incompletos = 0
    
    for record in SeqIO.parse(arquivo_limpo, "fasta"):
        header = record.description.lower()
        tamanhos.append(len(record.seq))
        
        # Extrair família
        familia = header.split('-')[0] if '-' in header else header.split('_')[0]
        familias_encontradas.append(familia.upper())
        
        # Contar completos vs incompletos
        if 'comp' in header and 'incomp' not in header:
            completos += 1
        elif 'incomp' in header:
            incompletos += 1
    
    if len(tamanhos) == 0:
        print("❌ Nenhuma sequência no arquivo limpo!")
        return
    
    # Estatísticas gerais
    contagem_familias = Counter(familias_encontradas)
    
    print(f"📈 ESTATÍSTICAS FINAIS:")
    print(f"   Total de TEs limpos: {len(tamanhos)}")
    print(f"   Tamanho médio: {sum(tamanhos)/len(tamanhos):.1f} bp")
    print(f"   Faixa de tamanho: {min(tamanhos)} - {max(tamanhos)} bp")
    print(f"   Elementos completos: {completos}")
    print(f"   Elementos incompletos: {incompletos}")
    
    print(f"\n🏷️ DISTRIBUIÇÃO DAS FAMÍLIAS:")
    for familia, count in contagem_familias.most_common():
        percentual = count / len(familias_encontradas) * 100
        print(f"   {familia}: {count} ({percentual:.1f}%)")
    
    # Correlação com dados esperados (baseado no artigo)
    print(f"\n🔬 INTERPRETAÇÃO BIOLÓGICA:")
    rlx_percent = contagem_familias.get('RLX', 0) / len(familias_encontradas) * 100
    dtx_percent = contagem_familias.get('DTX', 0) / len(familias_encontradas) * 100
    
    print(f"   RLX (LTR retrotransposons): {rlx_percent:.1f}%")
    print(f"   DTX (DNA transposons): {dtx_percent:.1f}%")
    print(f"   Proporção completos/incompletos: {completos}:{incompletos}")
    
    if rlx_percent > dtx_percent:
        print(f"   ✅ Padrão consistente: LTR > DNA transposons (esperado pelo artigo)")
    else:
        print(f"   ⚠️ Padrão inesperado: DNA transposons > LTR")

def balancear_datasets(arquivo_tes_limpo, arquivo_genes, proporcao=1.0):
    """Balanceia os datasets para treinamento"""
    
    print(f"\n⚖️ BALANCEANDO DATASETS PARA TREINAMENTO")
    print("=" * 45)
    
    # Contar TEs limpos
    tes_count = sum(1 for _ in SeqIO.parse(arquivo_tes_limpo, "fasta"))
    genes_count = sum(1 for _ in SeqIO.parse(arquivo_genes, "fasta"))
    
    print(f"📊 CONTAGENS ATUAIS:")
    print(f"   TEs limpos: {tes_count}")
    print(f"   Genes disponíveis: {genes_count}")
    
    # Calcular quantidade ideal
    genes_target = int(tes_count * proporcao)
    
    if genes_target <= genes_count:
        print(f"✅ Selecionaremos {genes_target} genes para balancear")
        
        # Selecionar genes aleatoriamente
        import random
        genes_list = list(SeqIO.parse(arquivo_genes, "fasta"))
        genes_selecionados = random.sample(genes_list, genes_target)
        
        # Salvar dataset balanceado
        arquivo_genes_balanceado = "dados_limpos/K8108_genes_balanceado.fasta"
        SeqIO.write(genes_selecionados, arquivo_genes_balanceado, "fasta")
        
        print(f"💾 Dataset balanceado salvo: {arquivo_genes_balanceado}")
        print(f"📈 Proporção final TEs:Genes = {tes_count}:{genes_target} (1:{proporcao})")
    else:
        print(f"⚠️ Poucos genes disponíveis. Usando todos os {genes_count} genes.")

def gerar_relatorio():
    """Gera relatório específico para os dados do P. pachyrhizi"""
    
    arquivo_relatorio = "relatorios/relatorio_limpeza.txt"
    
    with open(arquivo_relatorio, 'w', encoding='utf-8') as f:
        f.write("RELATÓRIO DE LIMPEZA - P. PACHYRHIZI\n")
        f.write("=" * 55 + "\n")
        f.write(f"Data: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Autora: Laura Tamarozzi\n\n")
        
        f.write("PADRÃO DE NOMENCLATURA IDENTIFICADO:\n")
        f.write("- RLX: LTR retrotransposons (família principal)\n")
        f.write("- DTX: DNA transposons (segunda família)\n")
        f.write("- RIX, DHX, DMX, DXX, RSX: Outras famílias\n")
        f.write("- noCat: Não categorizados (REMOVIDOS)\n\n")
        
        f.write("CRITÉRIOS DE LIMPEZA APLICADOS:\n")
        f.write("- Manter apenas famílias reconhecidas\n")
        f.write("- Remover elementos 'noCat'\n")
        f.write("- Tamanho mínimo: 100bp (500bp para incompletos)\n\n")
        
        f.write("ARQUIVOS GERADOS:\n")
        f.write("- dados_limpos/K8108_TEs_limpos.fasta\n")
        f.write("- dados_limpos/K8108_genes_nao_TEs.fasta\n")
        f.write("- dados_limpos/K8108_genes_balanceado.fasta\n")
        f.write("- dados_limpos/K8108_TEs_removidos.fasta\n\n")
    
    print(f"📋 Relatório salvo: {arquivo_relatorio}")

def main():
    """Função principal com limpeza"""
    
    print("🧹 LIMPEZA - PROJETO GRAMEP P. PACHYRHIZI")
    print("🔬 Baseada na análise real dos dados")
    print("👩‍🔬 Laura Tamarozzi")
    print("=" * 65)
    print(f"⏰ Início: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print()
    
    # Criar pastas
    criar_pastas()
    
    # Arquivos
    arquivo_tes = "dados_originais/K8108_TEs_anotados.fa"
    arquivo_genes = "dados_originais/K8108_nao_TEs_genes.fna"
    
    # ETAPA 1: Limpeza dos TEs
    print("1️⃣ LIMPEZA DOS TEs")
    estatisticas = limpar_tes(
        arquivo_tes,
        "dados_limpos/K8108_TEs_limpos.fasta",
        "dados_limpos/K8108_TEs_removidos.fasta"
    )
    
    # ETAPA 2: Preparar dados negativos
    print("\n2️⃣ PREPARANDO DADOS NEGATIVOS")
    preparar_dados_negativos(arquivo_genes, "dados_limpos/K8108_genes_nao_TEs.fasta")
    
    # ETAPA 3: Validar resultados
    print("\n3️⃣ VALIDANDO RESULTADOS")
    validar_limpeza("dados_limpos/K8108_TEs_limpos.fasta")
    
    # ETAPA 4: Balancear datasets
    print("\n4️⃣ BALANCEANDO DATASETS")
    balancear_datasets(
        "dados_limpos/K8108_TEs_limpos.fasta",
        "dados_limpos/K8108_genes_nao_TEs.fasta"
    )
    
    # ETAPA 5: Gerar relatório
    print("\n5️⃣ GERANDO RELATÓRIO")
    gerar_relatorio()
    
    # Finalização
    print("\n" + "=" * 65)
    print("✅ LIMPEZA CONCLUÍDA!")
    print("\n📁 ARQUIVOS PRONTOS PARA O GRAMEP:")
    print("   🎯 TEs limpos: dados_limpos/K8108_TEs_limpos.fasta")
    print("   🎯 Genes balanceados: dados_limpos/K8108_genes_balanceado.fasta")
    print("\n🚀 PRÓXIMO PASSO: Treinar o GRAMEP com esses arquivos!")
    print(f"⏰ Fim: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

if __name__ == "__main__":
    main()