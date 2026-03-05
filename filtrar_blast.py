import pandas as pd
import os, glob

# Variáveis para diretórios
results_dir = "blast_results/"
output_dir = "blast_filtrado/"
os.makedirs(output_dir, exist_ok=True)

# Variáveis filtro
ident_min = 30 # % identidade mínima
cover_min = 70 # % cobertura de alinhamento mínima
e_val_max = 1e-5 # e-value máximo

# Colunas dos tsv de resultados do blast
colunas = ["qseqid","sseqid","pident","qcovs","length","evalue","bitscore","stitle"]

# Função para limpar os IDs
def limpar_id(id_sujo):
    # Inicialmente assumindo que não tem id unipro
    uniprot_id = None
    
    # Separar o id, caso esteja com muita info
    if "|" in id_sujo:
        partes = id_sujo.split("|")

        # Seleciona somente o ID por si só
        id_limpo = partes[1]

        # identificadores de ID uniprot
        if partes[0] in ("sp", "tr"):
            uniprot_id = partes[1]
    else:
        id_limpo = id_sujo # Caso o id ja esteja "limpo"
    
    return id_limpo, uniprot_id

# armazenando o nome dos arquivos de resultado blast
arquivos_tsv = sorted(glob.glob(f"{results_dir}*.tsv"))

# Loop para operar sobre os arquivos
for tsv in arquivos_tsv:
    # Extraindo apenas o nome do composto
    nome_base = os.path.basename(tsv).replace("_proteinas_blast_result.tsv", "")
    print(f'\nProcessando resultados do BLAST para {nome_base}.')

    df = pd.read_csv(tsv, sep="\t", names=colunas) # Ler o tsv

    print(f"    Hits brutos: {len(df)}")

    # Aplicando o filtro
    df_filtrado = df[
        (df["pident"] >= ident_min) &
        (df["qcovs"] >= cover_min) &
        (df["evalue"] <= e_val_max)
    ].copy()

    print(f"    Hits pós filtro: {len(df_filtrado)}")

    if df_filtrado.empty:
        print(" Nenhum hit passou no filtro.")
        continue

    # Selecionando apenas os melhores hits individuais de acordo com o bitscore 
    df_melhores = (
        df_filtrado.sort_values("bitscore", ascending=False)
        .drop_duplicates(subset="qseqid", keep="first")
        .reset_index(drop=True)
    )

    print(f"    Melhores hits únicos: {len(df_melhores)}")

    # Limpando os indices, e adicionando tais colunas
    df_melhores.insert(3, "id_limpo_hit", df_melhores["sseqid"].apply(lambda x: limpar_id(x)[0]))
    df_melhores.insert(4, "id_uniprot_hit", df_melhores["sseqid"].apply(lambda x: limpar_id(x)[1]))

    # Renomeando colunas
    df_melhores = df_melhores.rename(columns={
        "qseqid": "id_query",
        "sseqid": "id_bruto_hit",
        "pident": "identidade",
        "qcovs": "cobertura_query",
        "length": "tamanho_align",
        "stitle": "info_hit" 
    })

    # Reordenando colunas
    df_melhores = df_melhores[[
        "id_query",
        "id_bruto_hit",
        "id_limpo_hit",
        "id_uniprot_hit",
        "identidade",
        "cobertura_query",
        "tamanho_align",
        "evalue",
        "bitscore",
        "info_hit"
    ]]

    # Escrevendo arquivos de saída
    saida = f"{output_dir}{nome_base}_blast_filtrado.csv"
    df_melhores.to_csv(saida, index=False)
    print(f"    Resultados do BLAST filtrados salvos em: {saida}")
