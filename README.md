# Pipeline para encontrar proteínas em Leishmania

Projeto para encontrar proteínas potencialmente interativas com certas moléculas (é pro meu TCC).
Muitas estruturas nesse projeto seguem padrões para utilizar na minha pesquisa.
Não necessáriamente é reprodutível para todos tipos de dados.

## Workflow
### Encontrando proteínas alvo de estudo sobre moleculas identidade ou similares
Inicialmente, deve-se procurar o smiles da molécula desejada no [PubChem](https://pubchem.ncbi.nlm.nih.gov), e buscar nas abas de "Identity" e "Similarity" pela opção de "Proteins" em "Linked Data". Essa opção mostra proteínas alvo dos estudos envolvendo a molécula identidade (se existir) ou similiares ao smiles forncecido.

### Extraindo as sequências das proteínas alvo referência
Em seguida, deve-se baixar o .csv das proteínas, contendo as informações de "protid", "Taxonomy" e "Protein_Accession". Após baixar os arquivos, deve-se preparar eles para o script [script_fasta_referencia.py](script_fasta_referencia.py).
- Dentro de "proteinas_triagem", criar um diretório com o nome identificador do composto (Ex.: RMS)
- Mover os arquivos baixados para a pasta criada
- Os arquivos devem seguir um identificador, se a origem foi de moléculas identidade ou similar (Ex.: RMS_I.csv)

A partir dai, é só executar o script, alterando a lista de compostos com os nomes dos identificadores de cada csv de cada molécula testada.

### BLAST com proteoma de Leishmania
Após obter os fastas das referêncais, parte-se para o BLAST em proteoma de Leishmania amazonensis, com o objetivo de extrair os ortólogos em Leishmania.

O repositório já possui o fasta do proteoma, e o database indexado. Entretanto, se não houver, deve-se realizar essa etapa (Baixar o fasta com o proteoma anotado, e indexação, seguindo instruções dentro do script do blast).

Em seguida, basta rodar o script [blast_script.sh](blast_script.sh), que então todo o processo será realizado, com os outputs em "blast_results", no formato .tsv.

**IMPORTANTE:** Eu criei esse script utilizando um ambiente virtual com a ferramenta do blast instalado, utilizando o mamba e miniconda3. Então é recomendado fazer o mesmo (ou alterar o script né). 

### Filtrando os resultados do BLAST
Com os resultados do blast em mãos, é necessário filtrar os melhores hits.

Os filtros que utilizei foram: Identidade >= 30%; Cobertura de alinhamento >= 70%; E-value <= 1.10^-5.

O script que realiza esta ação é o [filtrar_blast.py](filtrar_blast.py).

Após aplicar os filtros, os melhores hits ÚNICOS são extraidos, e aramzenados em arquivos .csv no diretório blast_filtrado.

### Obtendo FASTAS dos Hits
Tendo as IDs dos hits, nos arquivos .csv, o próximo passo é rodar o [script_fastas_hits.py](script_fastas_hits.py).

Este script realiza a busca dos IDs de cada hit dentro do proteoma, no formato FASTA, de Leishmania amazonensis, contido na pasta [proteoma_leish](proteoma_leish/). O script retorna um arquivo .fasta para cada arquivo .csv filtrado do BLAST.

Em seguida, deve-se seguir para a seleção das proteínas para a próxima etapa, de modelagem e docking.