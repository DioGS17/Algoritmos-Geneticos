import random
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

class Individuo():
    def __init__(self):
        self.nota_fitness = 0
        self.cromossomo = []            # 7 barras retas e 13 areas
        self.areas = []                 # Lista de áreas (mm^2)
        self.comp = []                  # Lista de comprimentos (m)
        self.lista_pontos = []
        self.massa_especifica = 7870    # Kg/m^3
        self.elasticidade = 200         # GPa
        self.delta = 0                  # Deformação no ponto C (mm)
        self.massa = 0.0                # Massa (Kg)
        self.fo = 0.0                   # Função objetivo

        for _ in range(7):  # adiciona os tamanhos das barras e calcula e adiciona diagonais
            gene = random.randint(1, 3)
            self.cromossomo.append(gene)

        for _ in range(13): # adiciona as areas ao cromossomo 
            num = random.randint(3, 5)
            self.cromossomo.append(num*100)


    def define_pontos(self):
        ponto_a = [0,0]
        ponto_b = [self.cromossomo[0], 0]
        ponto_c = [ponto_b[0]+self.cromossomo[1], 0]
        ponto_d = [ponto_c[0]+self.cromossomo[2], 0]
        ponto_e = [ponto_d[0]+self.cromossomo[3], 0]
        ponto_f= [ponto_b[0], self.cromossomo[4]]
        ponto_h= [ponto_c[0], self.cromossomo[5]]
        ponto_g= [ponto_d[0], self.cromossomo[6]]
        self.lista_pontos.append(ponto_a)
        self.lista_pontos.append(ponto_b)
        self.lista_pontos.append(ponto_c)
        self.lista_pontos.append(ponto_d)
        self.lista_pontos.append(ponto_e)
        self.lista_pontos.append(ponto_f)
        self.lista_pontos.append(ponto_g)
        self.lista_pontos.append(ponto_h)
    

    def calcula_massa_barras(self):
        massa_total = 0
        for i in range(13):
            massa_total += self.comp[i] * (self.areas[i]*10**(-6)) * self.massa_especifica
        return massa_total


    def analisa_base(self):
        # A base da treliça deve ter comprimento minimo de 8m
        soma_base = sum(self.cromossomo[:4])
        if soma_base < 8:
            return 0.0001
        return 1


    def funcao_delta(self):
        self.comp = self.cromossomo[:7]
        self.areas = self.cromossomo[7:]

        # Forças
        f = [0, 70, 10, 40, 10]
        # f = [10, 80, 20, 20, 10]
        carga = 1                       # Carga unitária
        fN = [0 for i in range(13)]     # Forças internas
        fn = [0 for i in range(13)]     # Forças virtuais

        # Ângulos
        alpha = np.arctan(self.comp[4]/self.comp[0])
        alpha2 = np.arctan(self.comp[6]/self.comp[3])
        beta = np.arctan(self.comp[0]/self.comp[4])
        beta2 = np.arctan(self.comp[3]/self.comp[6])
        gama = np.arctan(self.comp[1]/self.comp[4])
        gama2 = np.arctan(self.comp[2]/self.comp[6])
        theta = None            
        theta2 = None
        omega = np.arctan(self.comp[4]/self.comp[1])
        mi = np.arctan(self.comp[6]/self.comp[2])
        

        # Cálculo das reações
        # Reação em E
        rEy = (self.comp[0]*f[1] + sum(self.comp[:2])*f[2] + sum(self.comp[:3])*f[3] + sum(self.comp[:4])*f[4])/sum(self.comp[:4])

        # Reação em A
        rAy = sum(f) - rEy

        # Cálculo das forças internas
        # Nó A
        fN[4] = (f[0] - rAy)/np.sin(alpha)
        fN[0] = -fN[4]*np.cos(alpha)

        # Nó B
        fN[1] = fN[0]
        fN[5] = 0

        # Nó E
        fN[10] = (f[4] - rEy)/np.sin(alpha2)
        fN[3] = -fN[10]*np.cos(alpha2)
        # print(f"fN11 = {fN[10]}, fN4 = {fN[3]}")

        # Nó D
        fN[2] = fN[3]
        fN[9] = 0

        # Nó F
        if self.comp[5] > self.comp[4]:
            theta = np.arctan((self.comp[5] - self.comp[4])/self.comp[1])
            fN[11] = (fN[4]*(np.cos(beta)*np.tan(gama)+np.sin(beta)) + f[1]*np.tan(gama))/(np.sin(theta)*np.tan(gama) + np.cos(theta))
            fN[6] = (fN[4]*np.sin(beta) - fN[11]*np.cos(theta))/np.sin(gama)
        elif self.comp[5] == self.comp[4]:
            theta = np.arctan(self.comp[5]/self.comp[1])
            fN[6] = (-fN[4]*np.cos(beta) - fN[5] - f[1])/np.cos(gama)
            fN[11] = -fN[6]*np.cos(theta) + fN[4]*np.sin(beta)
        elif self.comp[5] < self.comp[4]:
            theta = np.arctan((self.comp[4] - self.comp[5])/self.comp[1])
            fN[11] = (fN[4]*(np.sin(beta)+np.cos(beta)*np.tan(gama)) + f[1]*np.tan(gama))/(np.cos(theta) - np.sin(theta)*np.tan(gama))
            fN[6] = (fN[4]*np.sin(beta) - fN[11]*np.cos(theta))/np.sin(gama)

        # Nó G
        if self.comp[5] > self.comp[6]:
            theta2 = np.arctan((self.comp[5] - self.comp[6])/self.comp[2])
            fN[12] = (fN[10]*(np.sin(beta2)+np.cos(beta2)*np.tan(gama2))+f[3]*np.tan(gama2))/(np.cos(theta2)+np.sin(theta2)*np.tan(gama2))
            fN[8] = (-fN[12]*np.cos(theta2)+fN[10]*np.sin(beta2))/np.sin(gama2)
        elif self.comp[5] == self.comp[6]:
            theta2 = np.arctan(self.comp[5]/self.comp[2])
            fN[8] = (-fN[9] - fN[10]*np.cos(beta2) - f[3])/np.cos(gama2)
            fN[12] = -fN[8]*np.sin(gama2) + fN[10]*np.sin(beta2)
        elif self.comp[5] < self.comp[6]:
            theta2 = np.arctan((self.comp[6] - self.comp[5])/self.comp[2])
            fN[12] = (fN[10]*(np.sin(beta2)+np.cos(beta2)*np.tan(gama2))+f[3]*np.tan(gama2))/(np.cos(theta2)-np.sin(theta2)*np.tan(gama2))
            fN[8] = (-fN[12]*np.cos(theta2)+fN[10]*np.sin(beta2))/np.sin(gama2)

        # Nó B
        fN[7] = -fN[6]*np.sin(omega) - fN[8]*np.sin(mi)

        # Cálculo das reações virtuais
        # Reação em E
        rnEy = (sum(self.comp[:2]) * carga)/sum(self.comp[:4])

        # Reação em A
        rnAy = carga - rnEy

        # Cálculo da forças virtuais
        # Nó A
        fn[4] = -rnAy/np.sin(alpha)
        fn[0] = -fn[4]*np.cos(alpha)

        # Nó B
        fn[1] = fn[0]
        fn[5] = 0

        # Nó E
        fn[10] = -rnEy/np.sin(alpha2)
        fn[3] = -fn[10]*np.cos(alpha2)

        # Nó D
        fn[2] = fn[3]
        fn[9] = 0

        # Nó F
        if self.comp[5] > self.comp[4]:
            fn[11] = (fn[4]*(np.cos(beta)*np.tan(gama)+np.sin(beta)))/(np.sin(theta)*np.tan(gama)+np.cos(theta))
            fn[6] = (fn[4]*np.sin(beta) - fn[11]*np.cos(theta))/np.sin(gama)
        elif self.comp[5] == self.comp[4]:
            fn[6] = (-fn[4]*np.cos(beta))/np.cos(gama)
            fn[11] = -fn[6]*np.cos(theta) + fn[4]*np.sin(beta)
        elif self.comp[5] < self.comp[4]:
            fn[11] = (fn[4]*(np.sin(beta)+np.cos(beta)*np.tan(gama)))/(np.cos(theta) - np.sin(theta)*np.tan(gama))
            fn[6] = (fn[4]*np.sin(beta) - fn[11]*np.cos(theta))/np.sin(gama)

        # Nó G
        if self.comp[5] > self.comp[6]:
            fn[12] = (fn[10]*(np.sin(beta2) + np.cos(beta2)*np.tan(gama2)))/(np.cos(theta2)+np.sin(theta2)*np.tan(gama2))
            fn[8] = (-fn[12]*np.cos(theta2) + fn[10]*np.sin(beta2))/np.sin(gama2)
        elif self.comp[5] == self.comp[6]:
            fn[8] = (-fn[10]*np.cos(beta2))/np.cos(gama2)
            fn[12] = -fn[8]*np.sin(gama2) + fn[10]*np.sin(beta2)
        elif self.comp[5] < self.comp[6]:
            fn[12] = (fn[10]*(np.sin(beta2)+np.cos(beta2)*np.tan(gama2)))/(np.cos(theta2) - np.sin(theta2)*np.tan(gama2))
            fn[8] = (-fn[12]*np.cos(theta2) + fn[10]*np.sin(beta2))/np.sin(gama2)

        # Nó B
        fn[7] = -fn[6]*np.sin(omega) - fn[8]*np.sin(mi) + carga

        # Inserindo as diagonais
        l5 = (self.comp[4]**2 + self.comp[0]**2)**(0.5)
        l7 = (self.comp[4]**2 + self.comp[1]**2)**(0.5)
        l9 = (self.comp[6]**2 + self.comp[2]**2)**(0.5)
        l11 = (self.comp[6]**2 + self.comp[3]**2)**(0.5)
        l12 = ((self.comp[4] - self.comp[5])**2 + self.comp[1]**2)**(0.5)
        l13 = ((self.comp[6] - self.comp[5])**2 + self.comp[2]**2)**(0.5)

        self.comp.insert(4, l5)
        self.comp.insert(6, l7)
        self.comp.insert(8, l9)
        self.comp.append(l11)
        self.comp.append(l12)
        self.comp.append(l13)

        # Deformação
        delta = [0 for i in range(13)]
        for i in range(len(delta)):
            delta[i] = (fN[i]*fn[i]*self.comp[i])/(self.areas[i]*self.elasticidade)

        # print("\nResultados da função delta")
        # for i in range(len(delta)):
        #     print(f"Barra {i+1}: fN = {fN[i]:.3f} kN; fn = {fn[i]:.3f} kN; L = {self.comp[i]:.2f} m; a = {self.areas[i]} mm^2; d = {delta[i]:.5f} m")

        # print(f"{sum(delta*10**3):.3f} mm")

        return sum(delta)*10**3


    def fitness(self):
        self.delta = self.funcao_delta()
        self.massa = self.calcula_massa_barras()

        self.fo = 1/(0.1*self.massa + self.delta)
        self.nota_fitness = self.fo * self.analisa_base()
    
    def crossover(self, outro_individuo):
        corte = round(random.random() * 20)
        cromossomo1 = []
        cromossomo2 = []
        cromossomo1 = outro_individuo.cromossomo[0:corte] + self.cromossomo[corte:]
        cromossomo2 = self.cromossomo[0:corte] + outro_individuo.cromossomo[corte:]       

        filho1 = Individuo()
        filho2 = Individuo()

        filho1.cromossomo = cromossomo1
        filho2.cromossomo = cromossomo2
        return filho1, filho2
        
    def mutacao(self, taxa_mutacao):
        if random.random() < taxa_mutacao:
            n = random.randint(0, 6)
            self.cromossomo[n] = random.randint(1, 3)

        if random.random() < taxa_mutacao:
            n = random.randint(7, 19)
            self.cromossomo[n] = random.randint(3, 5)*100
        return self

    

class AlgoritmoGenetico():
    def __init__(self, tamanho_populacao):
        self.tamanho_populacao = tamanho_populacao
        self.populacao = []
        self.melhor_solucao = None
        self.elite = []
        self.lista_solucoes = []
        self.medias = []
        self.std_deviations = []

        self.indice_melhor_solucao = 0
        self.o_top = None

    def inicializa_populacao(self):
        for i in range(self.tamanho_populacao):
            self.populacao.append(Individuo())
        self.o_top = self.populacao[0]
    
    def ordena_populacao(self):
        self.populacao = sorted(self.populacao, key=lambda ind: ind.nota_fitness, reverse=True)
        self.melhor_solucao = self.populacao[0]
        if self.melhor_solucao.nota_fitness>self.o_top.nota_fitness:
            self.o_top = self.melhor_solucao
            self.indice_melhor_solucao = GENERATION_I
    

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
    
    def get_std_deviation(self, media):
        var = 0
        for i in self.populacao:
            var += ((i.nota_fitness - media) * (i.nota_fitness - media))/(self.tamanho_populacao)
        return math.sqrt(var)

    def evolui(self, taxa_mutacao):
        self.elite = self.populacao[:2]
        nova_populacao = []
        soma_avaliacao = self.soma_avaliacoes()
        self.medias.append(soma_avaliacao/self.tamanho_populacao)
        self.std_deviations.append(self.get_std_deviation(soma_avaliacao/self.tamanho_populacao))
        
        while len(nova_populacao) < self.tamanho_populacao:
            pai1 = self.seleciona_pai(soma_avaliacao)
            pai2 = self.seleciona_pai(soma_avaliacao)
            
            filho1, filho2 = self.populacao[pai1].crossover(self.populacao[pai2])
            
            nova_populacao.append(filho1.mutacao(taxa_mutacao))
            nova_populacao.append(filho2.mutacao(taxa_mutacao))
        
        self.populacao = nova_populacao
        for individuo in self.populacao:
            individuo.fitness()
        
        self.ordena_populacao()
        self.populacao[-2:] = self.elite
        self.ordena_populacao()
        print(f"Geração {GENERATION_I+1}: {self.melhor_solucao.nota_fitness}")

        self.lista_solucoes.append(self.melhor_solucao.nota_fitness)

            
    def mostra_melhor(self):
        print(f"\nGeração: {self.indice_melhor_solucao}")
        print(f"\nDelta: {ag.o_top.delta:.2f} mm")
        print(f"\nMassa: {ag.o_top.massa:.1f} kg")
        print(f"\nFunção objetivo: {ag.o_top.fo}")
        print(f"\nCromossomo: {ag.melhor_solucao.cromossomo}")
    
    def define_cor(self, index):
        if self.o_top.areas[index] == 300:
            return '#00008B'
        elif self.o_top.areas[index] == 400:
            return '#00E813'
        elif self.o_top.areas[index] == 500:
            return '#EC0015'

    def desenhar_trelica(self):
        self.o_top.define_pontos()
        lista_pontos = self.o_top.lista_pontos

        plt.figure(figsize=(8,8))
        plt.title('Treliça Geração: {}'.format(self.indice_melhor_solucao))
        
        plt.axis([-0.5, 8, -0.5, 8])
        
        blue_line = mlines.Line2D([], [], color='#00008B', label='A1 = 3x10^-4 m²')
        red_line = mlines.Line2D([], [], color='#EC0015', label='A3 = 5x10^-4 m²')
        green_line = mlines.Line2D([], [], color='#00E813', label='A2 = 4x10^-4 m²')
        plt.legend(handles=[blue_line, green_line, red_line, ])
        plt.text(0.5, 6.5, f'Delta: {self.o_top.delta:.2f} mm')
        plt.text(0.5, 7, f'Massa Total: {self.o_top.massa:.1f} Kg')
        plt.text(0.5, 7.5, f'Função Objetivo: {self.o_top.fo:.4f}/Kg*mm')
        plt.text(1, 5.5, "Guilherme Montenegro, Diogo Gomes, Caio Bertoldo,\n Leonardo Abinader, Rodolfo Simoes")
        
        plt.grid(color='w')
        plt.plot((lista_pontos[1][0],lista_pontos[0][0]), (lista_pontos[1][1],lista_pontos[0][1]), color=self.define_cor(0), linewidth=2)
        plt.plot((lista_pontos[2][0],lista_pontos[1][0]), (lista_pontos[2][1],lista_pontos[1][1]), color=self.define_cor(1), linewidth=2)
        plt.plot((lista_pontos[3][0],lista_pontos[2][0]), (lista_pontos[3][1],lista_pontos[2][1]), color=self.define_cor(2), linewidth=2)
        plt.plot((lista_pontos[4][0],lista_pontos[3][0]), (lista_pontos[4][1],lista_pontos[3][1]), color=self.define_cor(3), linewidth=2)

        
        plt.plot((lista_pontos[5][0],lista_pontos[0][0]), (lista_pontos[5][1],lista_pontos[0][1]), color=self.define_cor(4), linewidth=2)
        plt.plot((lista_pontos[5][0],lista_pontos[1][0]), (lista_pontos[5][1],lista_pontos[1][1]), color=self.define_cor(5), linewidth=2)
        plt.plot((lista_pontos[5][0],lista_pontos[2][0]), (lista_pontos[5][1],lista_pontos[2][1]), color=self.define_cor(6), linewidth=2)
        plt.plot((lista_pontos[2][0],lista_pontos[7][0]), (lista_pontos[2][1],lista_pontos[7][1]), color=self.define_cor(7), linewidth=2)
        plt.plot((lista_pontos[2][0],lista_pontos[6][0]), (lista_pontos[2][1],lista_pontos[6][1]), color=self.define_cor(8), linewidth=2)
        plt.plot((lista_pontos[3][0],lista_pontos[6][0]), (lista_pontos[3][1],lista_pontos[6][1]), color=self.define_cor(9), linewidth=2)
        plt.plot((lista_pontos[6][0],lista_pontos[4][0]), (lista_pontos[6][1],lista_pontos[4][1]), color=self.define_cor(10), linewidth=2)
        plt.plot((lista_pontos[5][0],lista_pontos[7][0]), (lista_pontos[5][1],lista_pontos[7][1]), color=self.define_cor(11), linewidth=2)
        plt.plot((lista_pontos[7][0],lista_pontos[6][0]), (lista_pontos[7][1],lista_pontos[6][1]), color=self.define_cor(12), linewidth=2)

        plt.show()
        plt.close()
    
    def desenha_graficos(self):
        plt.plot(self.lista_solucoes)
        plt.title('Evolução do fitness')
        plt.xlabel('Geração')
        plt.ylabel('Fitness')
        plt.show()


        x = []
        for i in range (1, GENERATION):
            x.append(i)

        # Criação do gráfico de erro
        plt.figure(figsize=(10, 6))
        plt.errorbar(x, self.medias, yerr=self.std_deviations, fmt='o', ecolor='black', capsize=5, linestyle='-', color='red', alpha=0.5)
        plt.plot(self.lista_solucoes)

        # Adição de rótulos e título
        plt.xlabel('Geração')
        plt.ylabel('Valores')
        plt.title('Fitness x (Média Fitness ± Desvio Padrão)')
        plt.grid(True)

        # Exibição do gráfico
        plt.show()



# Parâmetros do AG
POPULATION_SIZE = 750
MUTATE_RATE = 0.05
GENERATION = 250
GENERATION_I = 0

# # Geração 1
ag = AlgoritmoGenetico(POPULATION_SIZE)
ag.inicializa_populacao()
ag.ordena_populacao()
GENERATION_I += 1
print(f"\nGeração {GENERATION_I}: {ag.melhor_solucao.nota_fitness}")

print("\n=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n")

# Geração Z
while GENERATION_I < GENERATION:
    ag.evolui(MUTATE_RATE)
    GENERATION_I += 1
    print("\n=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n")

print("-=-= Melhor solução =-=-")
ag.mostra_melhor()
ag.desenhar_trelica()
ag.desenha_graficos()