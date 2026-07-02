"""
comparacao_manual_vs_deep.py
============================

Compara duas estratégias de ML para detectar TEs:

  1. MANUAL: Random Forest com features pré-extraídas
     - Input: k-mers, frequência, diversidade simbólica, delta_TE
     - Vantagem: interpretável, rápido, precisa pouco dado
     - Processamento: paralelo (um modelo por k-mer)

  2. DEEP: CNN 1D com autofeaturization
     - Input: sequência DNA bruta (A/C/G/T codificados)
     - Vantagem: descobre padrões automáticos, potencialmente mais preciso
     - Processamento: a rede aprende features sozinha

DEPENDÊNCIAS
─────────────
pip install biopython pandas numpy matplotlib scikit-learn tensorflow rich openpyxl

USO
───
python comparacao_manual_vs_deep.py

SAÍDA
─────
resultados/
  ├─ ml_manual/
  │   ├─ modelo_random_forest.pkl
  │   ├─ metricas_manual.csv
  │   └─ features_importancia.png
  │
  ├─ ml_deep/
  │   ├─ modelo_cnn_1d.h5
  │   ├─ historico_treinamento.png
  │   └─ features_aprendidas.png
  │
  └─ comparacao/
      ├─ comparacao_metricas.csv
      ├─ comparacao_visual.png
      ├─ matriz_confusao_manual.png
      ├─ matriz_confusao_deep.png
      └─ relatorio_comparativo.txt
"""

import os
import gc
import logging
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from collections import defaultdict

# Machine Learning
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (
    accuracy_score, precision_score, recall_score, f1_score,
    roc_auc_score, confusion_matrix, roc_curve, auc
)
import pickle

# Deep Learning
import tensorflow as tf
from tensorflow.keras import Sequential, layers
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.callbacks import EarlyStopping, ReduceLROnPlateau

# Bioinformatics
from Bio import SeqIO

# Logging
from rich.logging import RichHandler
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, MofNCompleteColumn
from rich.table import Table
from rich.console import Console

logging.basicConfig(level="INFO", format="%(message)s", handlers=[RichHandler(show_time=False)])
log = logging.getLogger("rich")
console = Console()


# ─────────────────────────────────────────────────────────────────────────────
# PREPARAÇÃO DE DADOS
# ─────────────────────────────────────────────────────────────────────────────

def carregar_sequencias_fasta(arquivo_fasta, label):
    """
    Carrega sequências do FASTA com label (TE ou Controle).
    Filtra sequências com N e muito curtas.
    """
    sequencias = []
    labels = []
    
    with open(arquivo_fasta, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seq = str(record.seq).upper()
            
            # Filtra
            if "N" in seq or len(seq) < 100:
                continue
            
            sequencias.append(seq)
            labels.append(label)
    
    return sequencias, labels


def codificar_dna_numerico(sequencia, max_len=2000):
    """
    Codifica sequência DNA em números (A=1, C=2, G=3, T=4, pad=0).
    Padroniza tamanho com truncamento ou padding.
    """
    mapping = {'A': 1, 'C': 2, 'G': 3, 'T': 4}
    
    # Converte para números
    numerica = np.array([mapping.get(base, 0) for base in sequencia])
    
    # Padroniza tamanho
    if len(numerica) > max_len:
        numerica = numerica[:max_len]  # Trunca se maior
    else:
        numerica = np.pad(numerica, (0, max_len - len(numerica)), 'constant')  # Padding
    
    return numerica


def extrair_features_manuais_sequencia(sequencia, kmer_fortes_set):
    """
    Para uma sequência, calcula:
    - quantos k-mers "candidato_forte" ela contém
    - frequência de aparição de cada k-mer
    
    Retorna um vetor de features para passar ao Random Forest.
    """
    features = {}
    k = len(next(iter(kmer_fortes_set)))  # kmer_fortes_set é um set, não dict
    
    for i in range(len(sequencia) - k + 1):
        kmer = sequencia[i:i+k]
        if kmer in kmer_fortes_set:
            if kmer not in features:
                features[kmer] = 0
            features[kmer] += 1
    
    # Normaliza por tamanho da sequência
    total_kmers = len(sequencia) - k + 1
    for kmer in features:
        features[kmer] = features[kmer] / total_kmers if total_kmers > 0 else 0
    
    return features


# ─────────────────────────────────────────────────────────────────────────────
# MODELO 1: MANUAL FEATURES + RANDOM FOREST
# ─────────────────────────────────────────────────────────────────────────────

def treinar_modelo_manual(X_train, X_test, y_train, y_test, pasta_destino):
    """
    Treina Random Forest com features manuais pré-extraídas.
    """
    log.info("🤖 Treinando Random Forest (features manuais)...")
    
    # Treina modelo
    modelo = RandomForestClassifier(
        n_estimators=100,
        max_depth=15,
        random_state=42,
        n_jobs=-1,
        class_weight='balanced'  # Importante para dados desbalanceados
    )
    modelo.fit(X_train, y_train)
    
    # Predições
    y_pred = modelo.predict(X_test)
    y_proba = modelo.predict_proba(X_test)[:, 1]
    
    # Métricas
    metricas = {
        'accuracy': accuracy_score(y_test, y_pred),
        'precision': precision_score(y_test, y_pred),
        'recall': recall_score(y_test, y_pred),
        'f1': f1_score(y_test, y_pred),
        'roc_auc': roc_auc_score(y_test, y_proba),
        'n_features': X_train.shape[1],
        'tempo_treino': 'instant'
    }
    
    log.info(f"   ✅ Acurácia: {metricas['accuracy']:.4f}")
    log.info(f"   ✅ Precisão: {metricas['precision']:.4f}")
    log.info(f"   ✅ Recall: {metricas['recall']:.4f}")
    log.info(f"   ✅ F1-Score: {metricas['f1']:.4f}")
    log.info(f"   ✅ ROC-AUC: {metricas['roc_auc']:.4f}")
    
    # Salva modelo
    with open(os.path.join(pasta_destino, "modelo_random_forest.pkl"), "wb") as f:
        pickle.dump(modelo, f)
    
    # Feature importance
    feature_importance = pd.DataFrame({
        'feature_index': range(len(modelo.feature_importances_)),
        'importance': modelo.feature_importances_
    }).sort_values('importance', ascending=False)
    
    # Visualização
    fig, ax = plt.subplots(figsize=(10, 6))
    top_n = min(20, len(feature_importance))
    top_features = feature_importance.head(top_n)
    ax.barh(top_features['feature_index'].astype(str), top_features['importance'], color='steelblue')
    ax.set_xlabel('Importância', fontsize=11)
    ax.set_title('Top 20 Features — Random Forest', fontsize=12)
    plt.tight_layout()
    plt.savefig(os.path.join(pasta_destino, "features_importancia.png"), dpi=300, bbox_inches='tight')
    plt.close()
    
    # Matriz de confusão
    cm = confusion_matrix(y_test, y_pred)
    fig, ax = plt.subplots(figsize=(8, 6))
    im = ax.imshow(cm, interpolation='nearest', cmap='Blues')
    ax.figure.colorbar(im, ax=ax)
    ax.set(xticks=np.arange(cm.shape[1]), yticks=np.arange(cm.shape[0]),
           xticklabels=['Controle', 'TE'], yticklabels=['Controle', 'TE'],
           ylabel='Real', xlabel='Predito')
    for i in range(cm.shape[0]):
        for j in range(cm.shape[1]):
            ax.text(j, i, format(cm[i, j], 'd'), ha="center", va="center",
                   color="white" if cm[i, j] > cm.max() / 2 else "black")
    ax.set_title('Matriz de Confusão — Random Forest')
    plt.tight_layout()
    plt.savefig(os.path.join(pasta_destino, "matriz_confusao_manual.png"), dpi=300, bbox_inches='tight')
    plt.close()
    
    # Curva ROC
    fpr, tpr, _ = roc_curve(y_test, y_proba)
    roc_auc = auc(fpr, tpr)
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(fpr, tpr, color='darkorange', lw=2, label=f'ROC curve (AUC = {roc_auc:.3f})')
    ax.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--', label='Acaso (AUC = 0.5)')
    ax.set_xlim([0.0, 1.0])
    ax.set_ylim([0.0, 1.05])
    ax.set_xlabel('Taxa de Falsos Positivos')
    ax.set_ylabel('Taxa de Verdadeiros Positivos')
    ax.set_title('Curva ROC — Random Forest')
    ax.legend(loc="lower right")
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(pasta_destino, "roc_curve_manual.png"), dpi=300, bbox_inches='tight')
    plt.close()
    
    return modelo, metricas, y_pred, y_proba


# ─────────────────────────────────────────────────────────────────────────────
# MODELO 2: AUTOFEATURIZATION COM CNN 1D
# ─────────────────────────────────────────────────────────────────────────────

def construir_cnn_1d(max_len=2000, num_classes=2):
    """
    Constrói arquitetura CNN 1D para sequências DNA.
    """
    modelo = Sequential([
        # Input: sequência numérica codificada (max_len,)
        # Embedding não é necessário pois já numerizamos; começamos com Conv1D direto
        
        # Conv block 1: aprende motivos pequenos (k-mers de 5-10 bp)
        layers.Conv1D(64, kernel_size=5, activation='relu', padding='same', input_shape=(max_len, 1)),
        layers.BatchNormalization(),
        layers.MaxPooling1D(pool_size=2),
        layers.Dropout(0.2),
        
        # Conv block 2: aprende motivos médios (k-mers de 10-20 bp)
        layers.Conv1D(128, kernel_size=10, activation='relu', padding='same'),
        layers.BatchNormalization(),
        layers.MaxPooling1D(pool_size=2),
        layers.Dropout(0.2),
        
        # Conv block 3: aprende motivos maiores (padrões de 20+ bp)
        layers.Conv1D(256, kernel_size=15, activation='relu', padding='same'),
        layers.BatchNormalization(),
        layers.MaxPooling1D(pool_size=4),
        layers.Dropout(0.3),
        
        # Reduz dimensionalidade
        layers.GlobalMaxPooling1D(),
        
        # Dense layers para classificação
        layers.Dense(256, activation='relu'),
        layers.Dropout(0.4),
        layers.Dense(128, activation='relu'),
        layers.Dropout(0.3),
        
        # Output
        layers.Dense(num_classes, activation='softmax')
    ])
    
    return modelo


def treinar_modelo_deep(X_train, X_test, y_train, y_test, pasta_destino):
    """
    Treina CNN 1D com autofeaturization.
    """
    log.info("🧠 Treinando CNN 1D (autofeaturization)...")
    
    # Reshape para (amostras, tamanho_sequência, 1 canal)
    X_train_cnn = X_train.reshape(X_train.shape[0], X_train.shape[1], 1)
    X_test_cnn = X_test.reshape(X_test.shape[0], X_test.shape[1], 1)
    
    # Constrói modelo
    modelo = construir_cnn_1d(max_len=X_train.shape[1])
    modelo.compile(
        optimizer=Adam(learning_rate=0.001),
        loss='sparse_categorical_crossentropy',
        metrics=['accuracy', tf.keras.metrics.Precision(), tf.keras.metrics.Recall()]
    )
    
    # Callbacks
    early_stop = EarlyStopping(monitor='val_loss', patience=5, restore_best_weights=True)
    reduce_lr = ReduceLROnPlateau(monitor='val_loss', factor=0.5, patience=3, min_lr=1e-6)
    
    # Treina
    log.info("   Treinando (isso pode levar alguns minutos)...")
    historico = modelo.fit(
        X_train_cnn, y_train,
        validation_split=0.2,
        epochs=30,
        batch_size=32,
        callbacks=[early_stop, reduce_lr],
        verbose=0
    )
    
    # Predições
    y_pred_proba = modelo.predict(X_test_cnn, verbose=0)
    y_pred = np.argmax(y_pred_proba, axis=1)
    y_proba = y_pred_proba[:, 1]
    
    # Métricas
    metricas = {
        'accuracy': accuracy_score(y_test, y_pred),
        'precision': precision_score(y_test, y_pred),
        'recall': recall_score(y_test, y_pred),
        'f1': f1_score(y_test, y_pred),
        'roc_auc': roc_auc_score(y_test, y_proba),
        'n_features': 'automático (CNN aprendeu)',
        'tempo_treino': f'{len(historico.history["loss"])} epochs'
    }
    
    log.info(f"   ✅ Acurácia: {metricas['accuracy']:.4f}")
    log.info(f"   ✅ Precisão: {metricas['precision']:.4f}")
    log.info(f"   ✅ Recall: {metricas['recall']:.4f}")
    log.info(f"   ✅ F1-Score: {metricas['f1']:.4f}")
    log.info(f"   ✅ ROC-AUC: {metricas['roc_auc']:.4f}")
    
    # Salva modelo
    modelo.save(os.path.join(pasta_destino, "modelo_cnn_1d.h5"))
    
    # Histórico de treinamento
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    
    axes[0, 0].plot(historico.history['loss'], label='Train Loss')
    axes[0, 0].plot(historico.history['val_loss'], label='Val Loss')
    axes[0, 0].set_xlabel('Epoch')
    axes[0, 0].set_ylabel('Loss')
    axes[0, 0].set_title('Loss ao longo do treinamento')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)
    
    axes[0, 1].plot(historico.history['accuracy'], label='Train Accuracy')
    axes[0, 1].plot(historico.history['val_accuracy'], label='Val Accuracy')
    axes[0, 1].set_xlabel('Epoch')
    axes[0, 1].set_ylabel('Accuracy')
    axes[0, 1].set_title('Acurácia ao longo do treinamento')
    axes[0, 1].legend()
    axes[0, 1].grid(True, alpha=0.3)
    
    axes[1, 0].plot(historico.history['precision'], label='Train Precision')
    axes[1, 0].plot(historico.history['val_precision'], label='Val Precision')
    axes[1, 0].set_xlabel('Epoch')
    axes[1, 0].set_ylabel('Precision')
    axes[1, 0].set_title('Precisão ao longo do treinamento')
    axes[1, 0].legend()
    axes[1, 0].grid(True, alpha=0.3)
    
    axes[1, 1].plot(historico.history['recall'], label='Train Recall')
    axes[1, 1].plot(historico.history['val_recall'], label='Val Recall')
    axes[1, 1].set_xlabel('Epoch')
    axes[1, 1].set_ylabel('Recall')
    axes[1, 1].set_title('Recall ao longo do treinamento')
    axes[1, 1].legend()
    axes[1, 1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(pasta_destino, "historico_treinamento.png"), dpi=300, bbox_inches='tight')
    plt.close()
    
    # Matriz de confusão
    cm = confusion_matrix(y_test, y_pred)
    fig, ax = plt.subplots(figsize=(8, 6))
    im = ax.imshow(cm, interpolation='nearest', cmap='Blues')
    ax.figure.colorbar(im, ax=ax)
    ax.set(xticks=np.arange(cm.shape[1]), yticks=np.arange(cm.shape[0]),
           xticklabels=['Controle', 'TE'], yticklabels=['Controle', 'TE'],
           ylabel='Real', xlabel='Predito')
    for i in range(cm.shape[0]):
        for j in range(cm.shape[1]):
            ax.text(j, i, format(cm[i, j], 'd'), ha="center", va="center",
                   color="white" if cm[i, j] > cm.max() / 2 else "black")
    ax.set_title('Matriz de Confusão — CNN 1D')
    plt.tight_layout()
    plt.savefig(os.path.join(pasta_destino, "matriz_confusao_deep.png"), dpi=300, bbox_inches='tight')
    plt.close()
    
    # Curva ROC
    fpr, tpr, _ = roc_curve(y_test, y_proba)
    roc_auc = auc(fpr, tpr)
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(fpr, tpr, color='darkorange', lw=2, label=f'ROC curve (AUC = {roc_auc:.3f})')
    ax.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--', label='Acaso (AUC = 0.5)')
    ax.set_xlim([0.0, 1.0])
    ax.set_ylim([0.0, 1.05])
    ax.set_xlabel('Taxa de Falsos Positivos')
    ax.set_ylabel('Taxa de Verdadeiros Positivos')
    ax.set_title('Curva ROC — CNN 1D')
    ax.legend(loc="lower right")
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(pasta_destino, "roc_curve_deep.png"), dpi=300, bbox_inches='tight')
    plt.close()
    
    return modelo, metricas, y_pred, y_proba


# ─────────────────────────────────────────────────────────────────────────────
# COMPARAÇÃO E RELATÓRIO
# ─────────────────────────────────────────────────────────────────────────────

def gerar_relatorio_comparativo(metricas_manual, metricas_deep, pasta_destino):
    """
    Cria visualização e relatório comparativo dos dois modelos.
    """
    # Tabela comparativa
    df_comparacao = pd.DataFrame({
        'Métrica': ['Acurácia', 'Precisão', 'Recall', 'F1-Score', 'ROC-AUC'],
        'Random Forest (Manual)': [
            metricas_manual['accuracy'],
            metricas_manual['precision'],
            metricas_manual['recall'],
            metricas_manual['f1'],
            metricas_manual['roc_auc']
        ],
        'CNN 1D (Deep)': [
            metricas_deep['accuracy'],
            metricas_deep['precision'],
            metricas_deep['recall'],
            metricas_deep['f1'],
            metricas_deep['roc_auc']
        ]
    })
    
    # Calcula diferenças
    df_comparacao['Diferença (Deep - Manual)'] = (
        df_comparacao['CNN 1D (Deep)'] - df_comparacao['Random Forest (Manual)']
    )
    df_comparacao['Vencedor'] = df_comparacao['Diferença (Deep - Manual)'].apply(
        lambda x: '🟦 Deep' if x > 0.01 else ('🟩 Manual' if x < -0.01 else '🟪 Empate')
    )
    
    # Salva CSV
    df_comparacao.to_csv(os.path.join(pasta_destino, "comparacao_metricas.csv"), index=False)
    
    # Visualização em barras
    fig, ax = plt.subplots(figsize=(11, 7))
    
    x = np.arange(len(df_comparacao))
    width = 0.35
    
    bars1 = ax.bar(x - width/2, df_comparacao['Random Forest (Manual)'], width,
                   label='Random Forest (Manual)', color='steelblue', alpha=0.8)
    bars2 = ax.bar(x + width/2, df_comparacao['CNN 1D (Deep)'], width,
                   label='CNN 1D (Deep)', color='coral', alpha=0.8)
    
    ax.set_ylabel('Score', fontsize=11)
    ax.set_title('Comparação: Feature Engineering Manual vs. Autofeaturization', fontsize=12, pad=15)
    ax.set_xticks(x)
    ax.set_xticklabels(df_comparacao['Métrica'])
    ax.legend(fontsize=10)
    ax.set_ylim([0, 1.1])
    ax.grid(axis='y', linestyle='--', alpha=0.3)
    
    # Adiciona valores nas barras
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{height:.3f}', ha='center', va='bottom', fontsize=8)
    
    plt.tight_layout()
    plt.savefig(os.path.join(pasta_destino, "comparacao_visual.png"), dpi=300, bbox_inches='tight')
    plt.close()
    
    # Relatório em texto
    relatorio = f"""
╔══════════════════════════════════════════════════════════════════════════════╗
║                 RELATÓRIO COMPARATIVO DE MODELOS DE ML                       ║
║                    Manual Feature Engineering vs. Deep Learning              ║
╚══════════════════════════════════════════════════════════════════════════════╝

┌─────────────────────────────────────────────────────────────────────────────┐
│ MODELO 1: RANDOM FOREST (FEATURE ENGINEERING MANUAL)                        │
└─────────────────────────────────────────────────────────────────────────────┘

  • Algoritmo:          Random Forest (100 árvores)
  • Features:           {metricas_manual['n_features']} k-mers candidatos_fortes
  • Tempo de treino:    {metricas_manual['tempo_treino']}
  
  Desempenho:
  ├─ Acurácia:         {metricas_manual['accuracy']:.4f} ({metricas_manual['accuracy']*100:.2f}%)
  ├─ Precisão:         {metricas_manual['precision']:.4f}
  ├─ Recall:           {metricas_manual['recall']:.4f}
  ├─ F1-Score:         {metricas_manual['f1']:.4f}
  └─ ROC-AUC:          {metricas_manual['roc_auc']:.4f}

  Interpretação:
  ✓ Cada feature (k-mer) tem um peso (importância) calculado
  ✓ Fácil explicar por quê o modelo tomou uma decisão
  ✓ Requer feature engineering prévio (nós já fizemos com Fisher + FDR)
  ✓ Rápido de treinar e predizer

┌─────────────────────────────────────────────────────────────────────────────┐
│ MODELO 2: CNN 1D (AUTOFEATURIZATION - DEEP LEARNING)                        │
└─────────────────────────────────────────────────────────────────────────────┘

  • Algoritmo:          Convolutional Neural Network 1D
  • Features:           Aprendidas automaticamente pela rede
  • Arquitetura:        3 camadas conv + GlobalMaxPooling + 2 camadas dense
  • Tempo de treino:    {metricas_deep['tempo_treino']}
  
  Desempenho:
  ├─ Acurácia:         {metricas_deep['accuracy']:.4f} ({metricas_deep['accuracy']*100:.2f}%)
  ├─ Precisão:         {metricas_deep['precision']:.4f}
  ├─ Recall:           {metricas_deep['recall']:.4f}
  ├─ F1-Score:         {metricas_deep['f1']:.4f}
  └─ ROC-AUC:          {metricas_deep['roc_auc']:.4f}

  Interpretação:
  ✓ Aprende automaticamente quais padrões DNA são importantes
  ✓ Convolução 1D procura motivos (k-mers) de diferentes tamanhos
  ✓ Mais flexível: não precisa de feature engineering prévio
  ✗ Mais "caixa preta" — difícil entender decisões específicas
  ✗ Mais lento de treinar, precisa de GPU para sequências gigantescas

┌─────────────────────────────────────────────────────────────────────────────┐
│ COMPARAÇÃO DIRETA                                                            │
└─────────────────────────────────────────────────────────────────────────────┘

"""
    
    # Adiciona tabela ao relatório
    relatorio += df_comparacao.to_string(index=False)
    
    # Análise
    melhor_acuracia = 'Deep' if metricas_deep['accuracy'] > metricas_manual['accuracy'] else 'Manual'
    melhor_f1 = 'Deep' if metricas_deep['f1'] > metricas_manual['f1'] else 'Manual'
    
    diferenca_acuracia = abs(metricas_deep['accuracy'] - metricas_manual['accuracy'])
    
    relatorio += f"""

┌─────────────────────────────────────────────────────────────────────────────┐
│ ANÁLISE E CONCLUSÕES                                                         │
└─────────────────────────────────────────────────────────────────────────────┘

  Melhor Acurácia:     {melhor_acuracia} ({diferenca_acuracia*100:.2f}% de diferença)
  Melhor F1-Score:     {melhor_f1}
  
"""
    
    if diferenca_acuracia < 0.05:
        relatorio += f"""
  📊 RESULTADO: Modelos têm desempenho SIMILAR
  
  Ambas estratégias funcionam bem neste dataset. Recomendação:
  
  → Use Random Forest (Manual) se:
    • Precisa explicar decisões (TCC, publicação)
    • Quer treinar rápido
    • Tem pouco dado ou poder computacional limitado
    • Quer interpretabilidade máxima
    
  → Use CNN 1D (Deep) se:
    • Vai trabalhar em múltiplos organismos (escalabilidade)
    • Quer máxima acurácia
    • Tem GPU disponível
    • Pode aceitar "caixa preta"
"""
    elif metricas_deep['accuracy'] > metricas_manual['accuracy']:
        relatorio += f"""
  📊 RESULTADO: Deep Learning (CNN 1D) é MELHOR
  
  A autofeaturization descobriu padrões que manual engineering não capturou.
  Recomendação: Use CNN 1D para produção, mas documente Random Forest
  como baseline nos papers/TCC.
"""
    else:
        relatorio += f"""
  📊 RESULTADO: Feature Engineering Manual (Random Forest) é MELHOR
  
  Neste caso, a abordagem manual foi mais eficiente. Muito comum!
  Recomendação: Use Random Forest para seu TCC/publicação.
"""
    
    relatorio += f"""

┌─────────────────────────────────────────────────────────────────────────────┐
│ PRÓXIMOS PASSOS                                                              │
└─────────────────────────────────────────────────────────────────────────────┘

  1. Validação em dados novos (sequências de outro isolado/espécie)
  2. Análise de features importantes (qual k-mer/ padrão mais discriminativo?)
  3. Ablation study (remover features uma a uma, ver queda de performance)
  4. Documentar no TCC: "Comparamos duas estratégias de ML..."

═══════════════════════════════════════════════════════════════════════════════
"""
    
    with open(os.path.join(pasta_destino, "relatorio_comparativo.txt"), "w") as f:
        f.write(relatorio)
    
    log.info(relatorio)


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

def main():
    caminho_script = os.path.dirname(os.path.abspath(__file__))
    # Sobe 2 níveis: de scripts/machine_learning/ para a raiz do projeto
    pasta_raiz = os.path.dirname(os.path.dirname(caminho_script))
    pasta_dados    = os.path.join(pasta_raiz, "dados_limpos")
    pasta_result   = os.path.join(pasta_raiz, "resultados")
    
    # Subpastas organizadas em resultados/machine_learning/
    pasta_ml       = os.path.join(pasta_result, "machine_learning")
    pasta_ml_manual = os.path.join(pasta_ml, "ml_manual")
    pasta_ml_deep   = os.path.join(pasta_ml, "ml_deep")
    pasta_comparacao = os.path.join(pasta_ml, "comparacao")
    
    for pasta in [pasta_result, pasta_ml, pasta_ml_manual, pasta_ml_deep, pasta_comparacao]:
        os.makedirs(pasta, exist_ok=True)
    
    # Arquivos de entrada (saídas dos scripts anteriores)
    arquivo_tes = os.path.join(pasta_dados, "balanceado", "K8108_TEs_balanceado.fasta")
    arquivo_ctrl = os.path.join(pasta_dados, "K8108_genes_nao_TEs.fasta")
    excel_input = os.path.join(pasta_result, "excel", "features_kmers_v3_bidirecional.xlsx")
    csv_candidatos = os.path.join(pasta_result, "csv", "teste_fisher_candidatos_fortes.csv")
    
    # Valida inputs
    for arq in [arquivo_tes, arquivo_ctrl, csv_candidatos]:
        if not os.path.exists(arq):
            log.error(f"❌  Arquivo não encontrado: {arq}")
            return
    
    log.info("=" * 80)
    log.info("🚀 COMPARAÇÃO: FEATURE ENGINEERING MANUAL vs. AUTOFEATURIZATION")
    log.info("=" * 80)
    
    # ── 1. Carrega sequências ──────────────────────────────────────────────
    log.info("\n📖 Carregando sequências...")
    seqs_te, labels_te = carregar_sequencias_fasta(arquivo_tes, 1)
    seqs_ctrl, labels_ctrl = carregar_sequencias_fasta(arquivo_ctrl, 0)
    
    log.info(f"   TE:       {len(seqs_te):,} sequências")
    log.info(f"   Controle: {len(seqs_ctrl):,} sequências")
    
    # Combina
    todas_seqs = seqs_te + seqs_ctrl
    todos_labels = labels_te + labels_ctrl
    
    # ── 2. Carrega k-mers candidatos (para modelo manual) ────────────────
    log.info("\n📊 Carregando k-mers candidatos fortes...")
    df_candidatos = pd.read_csv(csv_candidatos)
    kmer_fortes = set(df_candidatos['kmer'].values)
    log.info(f"   {len(kmer_fortes):,} k-mers candidatos fortes")
    
    # ── 3. Prepara dados para MODELO MANUAL ──────────────────────────────
    log.info("\n🔧 Preparando features manuais...")
    X_manual = []
    with Progress(SpinnerColumn(),
                  TextColumn("[progress.description]{task.description}"),
                  BarColumn(), MofNCompleteColumn()) as prg:
        task = prg.add_task("Extraindo k-mers por sequência...", total=len(todas_seqs))
        for seq in todas_seqs:
            features_dict = extrair_features_manuais_sequencia(seq, kmer_fortes)
            feature_vector = np.array([features_dict.get(kmer, 0) for kmer in kmer_fortes])
            X_manual.append(feature_vector)
            prg.update(task, advance=1)
    
    X_manual = np.array(X_manual)
    log.info(f"   Shape: {X_manual.shape} (amostras × features)")
    
    # Split treino/teste
    X_train_manual, X_test_manual, y_train, y_test = train_test_split(
        X_manual, todos_labels, test_size=0.2, random_state=42, stratify=todos_labels
    )
    
    # Normaliza (importante para RF com features de escalas diferentes)
    scaler = StandardScaler()
    X_train_manual = scaler.fit_transform(X_train_manual)
    X_test_manual = scaler.transform(X_test_manual)
    
    # ── 4. Prepara dados para MODELO DEEP ──────────────────────────────
    log.info("\n🔧 Preparando sequências para CNN (codificação numérica)...")
    X_deep = []
    with Progress(SpinnerColumn(),
                  TextColumn("[progress.description]{task.description}"),
                  BarColumn(), MofNCompleteColumn()) as prg:
        task = prg.add_task("Codificando DNA...", total=len(todas_seqs))
        for seq in todas_seqs:
            seq_numerica = codificar_dna_numerico(seq, max_len=2000)
            X_deep.append(seq_numerica)
            prg.update(task, advance=1)
    
    X_deep = np.array(X_deep)
    log.info(f"   Shape: {X_deep.shape} (amostras × tamanho_sequência)")
    
    # Split (mesmo split que manual para comparação justa)
    indices = np.arange(len(todas_seqs))
    idx_train, idx_test = train_test_split(indices, test_size=0.2, random_state=42,
                                           stratify=todos_labels)
    X_train_deep = X_deep[idx_train]
    X_test_deep = X_deep[idx_test]
    y_train = np.array(todos_labels)[idx_train]
    y_test = np.array(todos_labels)[idx_test]
    
    # ── 5. TREINA MODELO MANUAL ──────────────────────────────────────────
    log.info("\n" + "="*80)
    modelo_manual, metricas_manual, y_pred_manual, y_proba_manual = treinar_modelo_manual(
        X_train_manual, X_test_manual, y_train, y_test, pasta_ml_manual
    )
    
    # ── 6. TREINA MODELO DEEP ────────────────────────────────────────────
    log.info("\n" + "="*80)
    modelo_deep, metricas_deep, y_pred_deep, y_proba_deep = treinar_modelo_deep(
        X_train_deep, X_test_deep, y_train, y_test, pasta_ml_deep
    )
    
    # ── 7. COMPARAÇÃO ────────────────────────────────────────────────────
    log.info("\n" + "="*80)
    log.info("📊 Gerando relatório comparativo...")
    gerar_relatorio_comparativo(metricas_manual, metricas_deep, pasta_comparacao)
    
    log.info(f"\n✅ Modelos treinados e salvos!")
    log.info(f"   Manual:      {pasta_ml_manual}")
    log.info(f"   Deep:        {pasta_ml_deep}")
    log.info(f"   Comparação:  {pasta_comparacao}")
    log.info(f"\n📂 Estrutura: resultados/machine_learning/")
    log.info(f"   ├── ml_manual/")
    log.info(f"   ├── ml_deep/")
    log.info(f"   └── comparacao/  ← LEIA: relatorio_comparativo.txt")


if __name__ == "__main__":
    main()