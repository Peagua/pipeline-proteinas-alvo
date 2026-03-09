import pandas as pd
import glob,os,time
from Bio import SeqIO
from rcsbapi.search import SeqSimilarityQuery
from rcsbapi.data import DataQuery

input_dir = "blast_filtrado/"
path_proteoma = "proteoma_leish/Leishmania_amazonensis_AnnotatedProteins.fasta"
output_dir = "fastas_melhores_hits/"

os.makedirs(output_dir, exist_ok=True)

# Carregando o proteoma
proteoma = SeqIO.to_dict(SeqIO.parse(path_proteoma, "fasta"))

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

# Quantas proteínas triar
top_melhores = int(input("Quantas proteínas quer retornar como melhores hits? "))

# Puxando os arquivos
arquivos = sorted(glob.glob(f"{input_dir}*.csv"))

for arquivo in arquivos:

    # Extraindo nome e abrindo o arquivo como df
    nome = os.path.basename(arquivo).replace("_blast_filtrado.csv", "")
    df = pd.read_csv(arquivo)

    # Estabelecendo o score para triagem
    # Valores dos pesos podem estar sujeito a mudanças, de acordo com a prioridade
    df["score_triagem"]=(
        df["identidade"]*0.5 + df["cobertura_query"]*0.3 + (df["bitscore"]/df["bitscore"].max())*100*0.2
    ) 

    # Criando df ordenado pelo score da triagem
    df_ordenado = df.sort_values("score_triagem", ascending=False)

    print(f"\n{nome} - Top {top_melhores}")
    print(df_ordenado[[
        "id_query", "id_limpo_hit", "identidade",
        "cobertura_query", "bitscore", "score_triagem"
    ]].head(top_melhores).to_string(index=False)) # Printa apenas os Top melhores 

    # Agora, para colocar as sequências das top melhores em um único fasta por composto
    sequencias = []

    # Iterando cada ID das top melhores
    for id_hit in df_ordenado["id_limpo_hit"].head(top_melhores):
        # Usando a funçãozinha
        sequence = busca_id_proteoma(id_hit, proteoma)
        sequencias.append(sequence)

    # Escrevendo o fasta
    fasta_saida = f"{output_dir}{nome}_melhores_hits.fasta"
    SeqIO.write(sequencias, fasta_saida, "fasta")
    print(f"\n    {top_melhores} Melhores FASTAS salvos em {fasta_saida}")

    # Procurando por proteínas catalogadas no PDB que são similares às melhores sequências
    # Isso é principalmente para ter uma estrutura para se basear para os próximos passos
    # E também pra dar nome às proteínas que selecionei
    ## Como trabalho com LLa muito provavelmente não vou encontrar proteínas dessa espécie catalogadas no PDB
    ### Mas similares já serve pra dar um norte
    sequencias_busca = []

    # Pega APENAS a sequência dentro do arquivo fasta, e transforma em STR
    # Será usado pelo API do PDB
    for seq_record in SeqIO.parse(fasta_saida, "fasta"):
        sequencias_busca.append(str(seq_record.seq))
    
    print(f"\nIDs PDB para proteínas similares às Top {top_melhores} para {nome}:")

    # Iterando para cada sequência
    for i in range(len(sequencias_busca)):

        # Estabelecendo o query para a busca no API
        # O identity cutoff é opcional e subjetivo de acordo com o objetivo do estudo
        query = SeqSimilarityQuery(
            value=sequencias_busca[i],
            identity_cutoff=0.3
        )

        # Como resultado, quero o polymer_entity (ID PDB)
        results = list(query(return_type="polymer_entity"))
        
        # Agora, se for encontrado, realiza busca de metadados (para extrair o nome da proteína)
        if len(results)>0:
            entrada = results[0].replace("_1", "") # Limpao ID

            # Outro API de busca do PDB, agora para extrair metadados
            data_query = DataQuery(
                input_type="entries",
                input_ids=[entrada],
                return_data_list=["struct.title"] # Retorana elementos do título
            )

            # Executando o API de busca
            metadados = data_query.exec()
            
            # Agora, printando no terminal o ID PDB e nome da proteína
            entry = metadados["data"]["entries"][0]

            id_pdb = entry['rcsb_id']
            nome_proteina = entry['struct']['title']
            url_pdb = f"https://www.rcsb.org/structure/{id_pdb}"

            print(f"\n    {i+1}ª - {id_pdb} - {nome_proteina}\n        URL: {url_pdb}")
        
        time.sleep(0.4) # Respeitando o tempo do API


