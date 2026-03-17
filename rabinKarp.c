#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct d {
    char cod[9];
    char** genes;
    int qtdGenes;
    int porcentagem;
    int indice;
} Doenca;

typedef struct No {
    char* chave;
    unsigned long valHash;
    struct No* prox;
} No;

No** tabela;
int TAMANHO_HASH;

unsigned long valCaractere(char c) {
    return (unsigned long)c;
}

unsigned long potenciaMod(unsigned long base, int exp, unsigned long mod) {
    unsigned long resultado = 1;
    unsigned long b = base % mod;
    while (exp > 0) {
        if (exp & 1)
            resultado = (resultado * b) % mod;
        b = (b * b) % mod;
        exp >>= 1;
    }
    return resultado;
}

unsigned long indiceHash(unsigned long h) {
    return h % TAMANHO_HASH;
}

void inserir(const char* s, int K, unsigned long valHash) {
    unsigned long indice = indiceHash(valHash);

    No* no = (No*)malloc(sizeof(No));
    no->chave = strndup(s, K);
    no->valHash = valHash;
    no->prox = tabela[indice];
    tabela[indice] = no;
}

int existeHash(unsigned long valHash, const char* s, int K) {
    unsigned long indice = indiceHash(valHash);
    No* no = tabela[indice];
    while (no) {
        if (no->valHash == valHash && strncmp(no->chave, s, K) == 0) {
            return 1;
        }
        no = no->prox;
    }
    return 0;
}

void inicializar(int tamanhoDNA, int K) {
    int numSubs = tamanhoDNA - K + 1;
    if (numSubs <= 0)
        TAMANHO_HASH = 1;
    else
        TAMANHO_HASH = numSubs * 2 + 1;

    tabela = (No**)calloc(TAMANHO_HASH, sizeof(No*));
}

void liberarHash() {
    for (int i = 0; i < TAMANHO_HASH; i++) {
        No* no = tabela[i];
        while (no) {
            No* tmp = no;
            no = no->prox;
            free(tmp->chave);
            free(tmp);
        }
    }

    free(tabela);
}

void rabinKarp(const char* DNA, int tamanhoDNA, int K, unsigned long base, unsigned long modulo) {
    if (tamanhoDNA < K || modulo == 0)
        return;

    unsigned long h = potenciaMod(base, K - 1, modulo);
    unsigned long t = 0;

    for (int i = 0; i < K; i++)
        t = (base * t + valCaractere(DNA[i])) % modulo;

    inserir(DNA, K, t);

    for (int j = 0; j < tamanhoDNA - K; j++) {
        t = (base * ((t + modulo - (valCaractere(DNA[j]) * h) % modulo) % modulo) + valCaractere(DNA[j + K])) % modulo;
        inserir(DNA+j+1, K, t);
    }
}

int geneAtivo(const char* gene, int K, unsigned long base, unsigned long modulo) {
    int tamanhoGene = (int)strlen(gene);
    int numSub = tamanhoGene - K + 1;

    if (numSub <= 0 || modulo == 0)
        return 0;

    unsigned long h = potenciaMod(base, K - 1, modulo);
    unsigned long t = 0;
    for (int i = 0; i < K; i++)
        t = (base * t + valCaractere(gene[i])) % modulo;

    int encontradas = 0;
    if (existeHash(t, gene, K))
        encontradas++;

    for (int s = 0; s < numSub - 1; s++) {
        t = (base * ((t + modulo - (valCaractere(gene[s]) * h) % modulo) % modulo) + valCaractere(gene[s + K])) % modulo;
        if (existeHash(t, gene + s + 1, K))
            encontradas++;
    }

    float fracao = (float)encontradas / (float)numSub;
    return fracao >= 0.90; // >= 90% de semelhnca
}

void countingSort(Doenca* v, int n) {
    int maxPorcentagem = 0;

    for (int i = 0; i < n; i++) {
        if (v[i].porcentagem > maxPorcentagem)
            maxPorcentagem = v[i].porcentagem;
    }

    int* count = (int*)calloc(maxPorcentagem + 1, sizeof(int));
    Doenca* saida = (Doenca*)malloc(n * sizeof(Doenca));
    int* inicio = (int*)malloc((maxPorcentagem + 1) * sizeof(int));

    for (int i = 0; i < n; i++)
        count[v[i].porcentagem]++;

    int pos = 0;
    for (int p = maxPorcentagem; p >= 0; p--) {
        inicio[p] = pos;
        pos += count[p];
    }

    for (int i = 0; i < n; i++) {
        int p = v[i].porcentagem;
        saida[inicio[p]++] = v[i];
    }

    for (int i = 0; i < n; i++)
        v[i] = saida[i];

    free(count);
    free(saida);
    free(inicio);
}

int main(int argc, char* argv[]) {
    printf("#ARGS = %i\n", argc);
    printf("PROGRAMA = %s\n", argv[0]);
    printf("ARG1 = %s, ARG2 = %s\n", argv[1], argv[2]);

    unsigned long base = 257; //2^8 +1
    unsigned long modulo = 1000000009; // primo suficinetemente grande

    FILE* input = fopen(argv[1], "r");
    FILE* output = fopen(argv[2], "w");

    int K, M;
    fscanf(input, "%d", &K);

    char buffer[1000001];
    fscanf(input, "%1000000s", buffer);

    char* DNA = strdup(buffer);
    int tamanhoDNA = (int)strlen(DNA);

    fscanf(input, "%d", &M);
    Doenca* doencas = (Doenca*)malloc(M * sizeof(Doenca));

    for (int i = 0; i < M; i++) {
        fscanf(input, "%s %d", doencas[i].cod, &doencas[i].qtdGenes);
        doencas[i].genes = (char**)malloc(doencas[i].qtdGenes * sizeof(char*));
        doencas[i].indice = i;
        for (int j = 0; j < doencas[i].qtdGenes; j++) {
            char bufferGene[1001];
            fscanf(input, "%1000s", bufferGene);
            doencas[i].genes[j] = strdup(bufferGene);
        }
    }

    inicializar(tamanhoDNA, K);
    rabinKarp(DNA, tamanhoDNA, K, base, modulo);

    for (int i = 0; i < M; i++) {
        int ativos = 0;
        for (int g = 0; g < doencas[i].qtdGenes; g++) {
            if (geneAtivo(doencas[i].genes[g], K, base, modulo))
                ativos++;
        }
        if (doencas[i].qtdGenes > 0)
            doencas[i].porcentagem = (int)round((100.0 * ativos / doencas[i].qtdGenes));
        else
            doencas[i].porcentagem = 0;
    }

    countingSort(doencas, M);

    for (int i = 0; i < M; i++)
        fprintf(output, "%s->%d%%\n", doencas[i].cod, doencas[i].porcentagem);

    for (int i = 0; i < M; i++) {
        for (int j = 0; j < doencas[i].qtdGenes; j++)
            free(doencas[i].genes[j]);

        free(doencas[i].genes);
    }

    free(doencas);
    free(DNA);
    liberarHash();

    fclose(input);
    fclose(output);

    return 0;
}
