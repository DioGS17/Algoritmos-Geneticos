import random

class Balsa():
    def __init__(self):
        # Atributos de construção da balsa
        self.camada = 3             # Nº de camadas (eixo z)
        self.coluna = 3             # Nº de colunas (eixo x)
        self.linha = 4              # Nº de linhas  (eixo y)
        self.posicoes = []          # Matriz tridimensional binária. 0 = sem container; 1 = com container
        
        # Número de containers
        self.num_containers = 0

        # Atributos para calcular a estabilidade
        self.CM0 = [2, 2.5, 0]      # Coordenadas [x,y,z] do centro de massa(CM) da balsa incialmente vazia
        self.CM = [2, 2.5, 0]       # Coordenadas [x,y,z] do CM deslocado

        # Inicializa a matriz posicoes
        for i in range(self.camada):
            camada_coord = []
            for j in range(self.linha):
                x_coord = []
                for k in range(self.coluna):
                    x_coord.append(0)
                camada_coord.append(x_coord)
            self.posicoes.append(camada_coord)
        
    def print_posicoes(self):
        for i in range(self.camada):
            print(f"Camada {i}:")
            for j in range(self.linha):
                print(f"{j}: {self.posicoes[i][j]}")
            print()
    
    def distancia_cm(self):
        x_cm = 0
        y_cm = 0
        z_cm = 0
        m = 0
        for i in range(self.camada):
            for j in range(self.linha):
                for k in range(self.coluna):
                    if self.posicoes[i][j][k] == 1:
                        x_cm += k+1
                        y_cm += j+1
                        z_cm += i+1
                        m += 1

        x_cm = x_cm/m
        y_cm = y_cm/m
        z_cm = z_cm/m

        self.CM = [x_cm, y_cm, z_cm]
        d = ((x_cm - self.CM0[0])**2 + (y_cm - self.CM0[1])**2 + (z_cm - self.CM0[2])**2)**(1/2)

        return d
    
    def colocar_container(self, camada:int, linha:int, coluna:int):
        # Impede containers flutuantes ou sobrepostos
        if (camada > 0 and self.posicoes[camada-1][linha][coluna] == 0) or (self.posicoes[camada][linha][coluna] == 1):
            return 0
        
        self.posicoes[camada][linha][coluna] = 1
        self.num_containers += 1
        return 1


class Individuo():
    def __init__(self, geracao:int):
        self.geracao = geracao
        self.nota_fitness = 0
        self.cromossomo = []
        self.num_containers = 0
        
    def inicializa_cromossomo(self):
        for _ in range(30):
            gene = [
                # Coordenadas do container na balsa
                random.randint(0, 2),   # Deck (0, 1, 2)    (eixo z)
                random.randint(0, 3),   # Row (0, 1, 2, 3)  (eixo y)
                random.randint(0, 2)    # Column (0, 1, 2)  (eixo x)
            ]
            self.cromossomo.append(gene)
    
    def verifica_repeticoes(self):
        repeticoes = 0
        for i in range(30):
            for j in range(i+1, 30):
                if self.cromossomo[i] == self.cromossomo[j]:
                    repeticoes += 1
        return repeticoes

    def fitness(self, balsa):
        for i in range(30):
            n1 = balsa.colocar_container(self.cromossomo[i][0], self.cromossomo[i][1], self.cromossomo[i][2])
            if n1 == 1:
                self.nota_fitness += 1/(balsa.distancia_cm()+1)
        

class AlgoritmoGenetico():
    def __init__(self, tam_populacao:int, taxa_mutacao:float):
        self.tam_populacao = tam_populacao
        self.populacao = []
        self.taxa_mutacao = taxa_mutacao
        self.elites = []
        self.melhor_solucao = None

    def inicializa_populacao(self):
        for i in range(self.tam_populacao):
            ind = Individuo(1)
            ind.inicializa_cromossomo()
            self.populacao.append(ind)
        
        for ind in self.populacao:
            balsa = Balsa()
            ind.fitness(balsa)

    def ordena_populacao(self):
        self.populacao = sorted(self.populacao, key=lambda ind: ind.nota_fitness, reverse=True)
        self.melhor_solucao = self.populacao[0]
        
    def soma_avaliacoes(self):
        soma = sum(ind.nota_fitness for ind in self.populacao)
        return soma
    
    def seleciona_pai(self, soma_avaliacao):
        pai = -1
        valor_sorteado = random.random() * soma_avaliacao
        soma = 0
        i = 0
        while i < len(self.populacao) and soma < valor_sorteado:
            soma += self.populacao[i].nota_fitness
            pai += 1
            i += 1
        return pai
    
    def crossover(self, pai1:Individuo, pai2:Individuo, geracao:int):
        filho1 = Individuo(geracao)
        filho2 = Individuo(geracao)

        corte = round(random.random() * 30)

        filho1.cromossomo = pai1.cromossomo[0:corte] + pai2.cromossomo[corte:]
        filho2.cromossomo = pai2.cromossomo[0:corte] + pai1.cromossomo[corte:]
        return filho1, filho2
    
    def mutacao(self, filho:Individuo):
        if random.random() < self.taxa_mutacao:
            i1 = random.randint(0, 29)
            i2 = random.randint(0, 29)
            filho.cromossomo[i1], filho.cromossomo[i2] = filho.cromossomo[i2], filho.cromossomo[i1]
    
    def evolui(self, geracao:int):
        self.elite = self.populacao[:2]
        nova_populacao = []
        soma_avaliacao = self.soma_avaliacoes()
        
        while len(nova_populacao) < self.tam_populacao:
            pai1 = self.seleciona_pai(soma_avaliacao)
            pai2 = self.seleciona_pai(soma_avaliacao)
            
            filho1, filho2 = self.crossover(self.populacao[pai1], self.populacao[pai2], geracao)
            
            self.mutacao(filho1)
            self.mutacao(filho2)
            nova_populacao.append(filho1)
            nova_populacao.append(filho2)
        
        self.populacao = nova_populacao
        for individuo in self.populacao:
            balsa = Balsa()
            individuo.fitness(balsa)
            individuo.num_containers = balsa.num_containers
        
        self.ordena_populacao()
        self.populacao[-2:] = self.elite
        self.ordena_populacao()

        print(f"\nGeração {geracao}: {self.melhor_solucao.nota_fitness}; Nº de containers: {self.melhor_solucao.num_containers}\n")


tam_populacao = 1000
n_geracoes = 1000
ag = AlgoritmoGenetico(tam_populacao, 0.05)
ag.inicializa_populacao()
ag.ordena_populacao()
print(f"\nGeração 1: {ag.melhor_solucao.nota_fitness}\n")

for i in range(1, n_geracoes):
    ag.evolui(i+1)

print(f"\nGeração: {ag.melhor_solucao.geracao}")
print(f"Melhor solução:")
i = 1
for gene in ag.melhor_solucao.cromossomo:
    print(f"Container {i}: Camada: {gene[0]}; Linha: {gene[1]}; Coluna: {gene[2]}")
    i += 1

print()

balsa = Balsa()
for gene in ag.melhor_solucao.cromossomo:
    balsa.colocar_container(gene[0], gene[1], gene[2])

balsa.print_posicoes()
d = balsa.distancia_cm()

print(f"\nCM: {balsa.CM}; Deslocamento: {d}")