#!/bin/bash

# Script para fazer o blast das proteínas referência para proteôma de Leishmania amazonensis

# Lembrando que é necessário fazer a indexação do proteoma (fasta) de referência antes de rodar esse script
# makeblastdb \
#   -in (onde ta o fasta do proteoma de referência)
#   -dbtype prot
#   -out (pasta de output)/(nome do arquivo) Ex: db_lla/Leishmania_amazonensis
#   -title "(título do proteoma)" Ex: "Leishmania amazonensis proteome"

echo "Script para realizar BLAST das proteínas referência encontradas."
echo ""

# Diretórios utilizados
INPUT_DIR="fastas_ref"
OUTPUT_DIR="blast_results"
DB_DIR="db_lla/Leishmania_amazonensis"

# Variáveis utilizadas
EVALUE="1e-5"
THREADS=4
MAX_HITS=5

# Criando o diretório de output, caso não exista
echo "Criando diretório de output."
mkdir -p $OUTPUT_DIR
echo ""

# Loop para fazer o blast com todos os fastas de referência
for arquivo in $INPUT_DIR/*.fasta; do
    
    nome=$(basename "$arquivo" .fasta)
    echo "Rodando BLAST para: $nome"

    mamba run -n blast_env blastp \
        -query "$arquivo" \
        -db "$DB_DIR" \
        -out "$OUTPUT_DIR/${nome}_blast_result.tsv" \
        -outfmt "6 qseqid sseqid pident qcovs length evalue bitscore stitle" \
        -evalue $EVALUE \
        -num_threads $THREADS \
        -max_target_seqs $MAX_HITS

    echo "BLAST para $nome terminado. Salvo em: $OUTPUT_DIR/${nome}_blast_result.tsv"
    echo ""
done

echo "BLASTp concluido para todos os arquivos!"
echo ""
echo "Os arquivos salvos contém os seguintes itens, em sequência:"
echo "ID_query; ID_hit; %ident; %cobertura_query; comprimento_align; e-value; qualidade_align; info_hit"

# Removendo arquivos vazios
for arquivo in $OUTPUT_DIR/*.tsv; do
    if [ ! -s $arquivo ]; then
        rm $arquivo
    fi
done

