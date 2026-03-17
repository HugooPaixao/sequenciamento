# Compressões Huffman e RLE

A empresa de biotecnologia Poxim Tech possui um sistema de diagnóstico para doenças genéticas, comparando a sequência de DNA com genes conhecidos.

A sequência de DNA é composta somente pelos símbolos A, C, G e T. Uma doença genética pode possuir vários genes associados, com sequências de tamanho entre 100 e 1000 símbolos, denotadas exclusivamente por letras maiúsculas e números com tamanho entre 4 e 8 caracteres.

Para tratar os efeitos da mutação nos genes, é feita a busca por combinações que possuam o tamanho mínimo de subcadeia, com pelo menos 90% de compatibilidade para manifestação da doença. No diagnóstico será calculada a probabilidade de manifestação da doença de acordo com a quantidade de genes detectados na sequência de DNA.

## Exemplo

Diagnóstico da doença CTRLF4 com genes `AATTTGGCCC` e `GGGGGGGGGG`.

* DNA: `AAAAAAAAAATTTTTTTTTTTGGGGGGGGG`
* Tamanho da subcadeia: 3
* AATTTGGCCC: 5 combinações = 50%
* GGGGGGGGGG: 9 combinações = 90%
* Chance de 50% de ocorrência da doença CTRLF4

## Formato de arquivo de entrada

```text
[#Tamanho da subcadeia]
[B0 ... BN-1]
[#Número de doenças]
[Código0] [#Genes0] [G00] ... [G0i-1]
...
[CódigoM-1] [#GenesM-1] [GM-10] ... [GM-1j-1]
```

## Exemplo de entrada

```text
3
AAAATTTCGTTAAATTTGAACATAGGGATA
4
ABCDE 3 AAA AAT AAAG
XY1WZ2AB 1 TTTTTTGGGG
H1N1 4 ACTG AACCGGTT AATAAT AAAAAAAAGA
HUEBR 1 CATAGGGATT
```

## Formato de arquivo de saída

É feita a ordenação estável em ordem decrescente dos resultados, utilizando como critério a probabilidade de ocorrência da doença e fazendo o arredondamento dos percentuais para fins de comparação.

## Exemplo de saída

```text
XY1WZ2AB ->100%
HUEBR ->100%
ABCDE ->67%
H1N1 ->25%
```

---

> Todas as entradas e saídas seguem rigorosamente os formatos das descrições.
> O algoritmo foi implementado para fins acadêmicos.
