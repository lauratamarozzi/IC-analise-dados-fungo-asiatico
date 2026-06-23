#!/usr/bin/env python3
"""
Script para analisar os headers dos TEs e entender o padrão de nomenclatura
Isso vai nos ajudar a ajustar os critérios de limpeza
"""

from Bio import SeqIO
from collections import Counter
import re

def analisar_headers_detalhado(arquivo_fasta):
    """Analisa os headers em detalhes para entender a nomenclatura"""
    
    print("🔍 ANALISANDO HEADERS DOS TEs EM DETALHES")
    print("=" * 50)
    
    headers_completos = []
    palavras_encontradas = []
    prefixos = []
    
    # Ler todas as sequências
    for i, record in enumerate(SeqIO.parse(arquivo_fasta, "fasta")):
        header_completo = record.description
        headers_completos.append(header_completo)
        
        # Extrair palavras do header
        palavras = re.findall(r'\b[A-Za-z]+\b', header_completo.lower())
        palavras_encontradas.extend(palavras)
        
        # Extrair prefixo (primeira parte antes do hífen)
        if '-' in header_completo:
            prefixo = header_completo.split('-')[0]
            prefixos.append(prefixo)
        
        # Mostrar primeiros 10 headers completos
        if i < 10:
            print(f"{i+1:2d}: {header_completo}")
    
    print(f"\nTotal de headers analisados: {len(headers_completos)}")
    
    # Analisar palavras mais comuns
    contador_palavras = Counter(palavras_encontradas)
    print(f"\n📊 PALAVRAS MAIS COMUNS NOS HEADERS:")
    for palavra, freq in contador_palavras.most_common(20):
        print(f"   {palavra}: {freq}")
    
    # Analisar prefixos
    contador_prefixos = Counter(prefixos)
    print(f"\n🏷️ PREFIXOS MAIS COMUNS:")
    for prefixo, freq in contador_prefixos.most_common(10):
        print(f"   {prefixo}: {freq}")
    
    # Procurar padrões de classificação
    print(f"\n🔎 PROCURANDO PADRÕES DE CLASSIFICAÇÃO:")
    
    # Palavras relacionadas a TEs que podem estar nos headers
    te_keywords = [
        'gypsy', 'copia', 'ltr', 'tir', 'line', 'sine', 'helitron', 
        'cacta', 'mariner', 'hat', 'mudr', 'retrotransposon', 'transposon',
        'repeat', 'element', 'dhx', 'comp', 'reversed', 'map'
    ]
    
    for keyword in te_keywords:
        count = sum(1 for header in headers_completos if keyword.lower() in header.lower())
        if count > 0:
            print(f"   '{keyword}': {count} headers")
    
    # Verificar se há algum padrão de família nos headers
    print(f"\n🧬 ANALISANDO POSSÍVEIS CLASSIFICAÇÕES:")
    for header in headers_completos[:20]:  # Primeiros 20
        # Tentar extrair partes que podem indicar classificação
        parts = re.split(r'[-_]', header)
        print(f"   Partes: {parts[:5]}...")  # Primeiras 5 partes

def main():
    arquivo = "dados_originais/K8108_TEs_anotados.fa"
    analisar_headers_detalhado(arquivo)

if __name__ == "__main__":
    main()