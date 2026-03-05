from Bio import SeqIO
import pandas as pd
import os, glob

# Variáveis diretório e caminho
input_dir = "blast_filtrado/"
path_proteoma = "proteoma_leish/Leishmania_amazonensis_AnnotatedProteins.fasta"
output_dir = "fastas_hits/"

os.makedirs(output_dir, exist_ok=True)

# Carregando o proteoma
print("Carregando o proteoma:")
proteoma = SeqIO.to_dict(SeqIO.parse(path_proteoma, "fasta"))
print(f"-> {len(proteoma)} sequêncais carregadas.")

# Função para extrair as sequências, a partir do ID (neste caso ID do Trypdb, presente no proteoma carregado)
def busca_id_proteoma(id, proteoma_usado):
    
    # Busca diretamente o id no proteoma
    if id in proteoma_usado:
        return proteoma_usado[id]

    # Se não encontrar o id, tenta encontrar parte dele (id do trypdb pode aparecer como sujo ou parte do orig)
    for chave in proteoma_usado:
        if id in chave or chave in id:
            return proteoma_usado[chave]

    return None # Caso não encontre de jeito nenhum


# Extraindo os arquivos de input
arquivos_csv = sorted(glob.glob(f"{input_dir}*.csv"))

# Loop para executar as ações para cada arquivo input
for arquivo in arquivos_csv:
    
    # Extraindo o nome do arquivo
    nome_base = os.path.basename(arquivo).replace("_blast_filtrado.csv", "")
    print(f"\nProcessando {nome_base}.\n")

    # Abrir o arquivo como df
    df = pd.read_csv(arquivo)

    # Listas importantes
    sequencias_extraidas = []
    nao_encontradas = []
    resultados_alphafold = []

    # Loop para cada linha (proteína) em cada arquivo
    for _, row in df.iterrows():

        id_limpo = row["id_limpo_hit"] # Extraindo o ID

        sequencia = busca_id_proteoma(id_limpo, proteoma) # Buscando o FASTA relacionado ao ID

        if sequencia is None:
            print(f"    Não encontrado no proteoma: {id_limpo}")
            nao_encontradas.append(id_limpo)
            continue

        # Adicionando a sequência encontrada à lista
        sequencias_extraidas.append(sequencia)
        str_sequencia = str(sequencia.seq)

        print(f"    {id_limpo} encontrada! ({len(str_sequencia)} aa)")

    # Salvar FASTA com as sequencias extraidas
    fasta_saida = f"{output_dir}{nome_base}_sequencias_hits.fasta"
    SeqIO.write(sequencias_extraidas, fasta_saida, "fasta")
    print(f"\n    FASTA salvo em: {fasta_saida} ({len(sequencias_extraidas)} sequências).")

