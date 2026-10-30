import networkx as nx
import math
import random
from random import choice
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from datetime import datetime

def SEIR_deletion(m,
                  n,
                  i0,
                  time_total,
                  step,
                  lamb,
                  beta,
                  gamma,
                  sigma,
                  miu,
                  pi_sin,
                  pi_iin,
                  pi_rin,
                  pi_sout,
                  pi_iout,
                  pi_rout):

    time_series = []  # 存储时间点
    S_series = []    # 存储S数量
    E_series = []    # 存储E数量
    I_series = []    # 存储I数量
    R_series = []    # 存储R数量

    global count_in
    global count_out
    count_in = 0
    count_out = 0

    kinds = ["S", "E", "I", "R"]  # 修改为SEIR四种状态
    status = [0]

    continous_time = {kinds[0]: {}, kinds[1]: {}, kinds[2]: {}, kinds[3]: {}}  # 增加E状态的统计

    graph = nx.generators.random_graphs.watts_strogatz_graph(n, m, .7)

    node_kind = {}
    node_out = {}
    time_seq = []

    I_degree = {}
    total_degree = {}
    k_i = {}


    def record_time():
        total_node = len(graph.nodes)
        # 记录各状态实际数量
        S_count = sum(1 for i in graph.nodes if node_kind[i] == kinds[0])
        E_count = sum(1 for i in graph.nodes if node_kind[i] == kinds[1])
        I_count = sum(1 for i in graph.nodes if node_kind[i] == kinds[2])
        R_count = sum(1 for i in graph.nodes if node_kind[i] == kinds[3])

        # 存储时间序列数据
        S_series.append(S_count)
        E_series.append(E_count)
        I_series.append(I_count)
        R_series.append(R_count)
        time_series.append(len(time_series))  # 记录时间点

        # 原有的记录
        time_seq.append((get_total_rate(kinds[0]) * total_node,
                         get_total_rate(kinds[1]) * total_node,
                         get_total_rate(kinds[2]) * total_node,
                         get_total_rate(kinds[3]) * total_node))


    def exp_distri(lamb):
        return math.log(1 - random.random()) * (-1 / lamb)


    def get_rate(target, kind):
        neighbors = graph.neighbors(target)
        total_n = len(neighbors)
        kind_n = 0
        for i in neighbors:
            if node_kind[i] == kind:
                kind_n += 1
        return 1.0 * kind_n / total_n


    def get_total_rate(kind):
        total_n = len(graph.nodes)
        kind_n = 0
        for i in graph.nodes:
            if node_kind[i] == kind:
                kind_n += 1
        return 1.0 * kind_n / total_n


    def get_time_interval(node):
        i_total = get_total_rate(kinds[2])  # I状态的比例
        s_total = get_total_rate(kinds[0])
        e_total = get_total_rate(kinds[1])  # E状态的比例
        r_total = get_total_rate(kinds[3])

        if node_kind[node] == kinds[0]:  # S状态
            k = len(graph.edges) * 2.0 / len(graph.nodes())
            rate_s_serv = beta * k * i_total * len(graph.nodes) + pi_sout * miu
            if rate_s_serv > 0:
                temp = exp_distri(rate_s_serv)
            else:
                temp = 1000
            if temp < node_out[node] or node_out[node] == 0:
                node_out[node] = temp

        elif node_kind[node] == kinds[1]:  # E状态
            rate_e_serv = sigma + pi_iout * miu
            temp = exp_distri(rate_e_serv)
            if temp < node_out[node] or node_out[node] == 0:
                node_out[node] = temp

        elif node_kind[node] == kinds[2]:  # I状态
            rate_i_serv = gamma + pi_iout * miu
            temp = exp_distri(rate_i_serv)
            if temp < node_out[node] or node_out[node] == 0:
                node_out[node] = temp

        elif node_kind[node] == kinds[3]:  # R状态
            rate_r_serv = pi_rout * miu  # R状态只有移除的可能
            temp = exp_distri(rate_r_serv)
            if temp < node_out[node] or node_out[node] == 0:
                node_out[node] = temp


    def new_node():
        node_new = 0
        while graph.has_node(node_new):
            node_new += 1
        k = 0
        while k < m:
            dead = choice(list(graph.nodes()))
            k += 1
            if dead != node_new:
                graph.add_edge(dead, node_new)
            if graph.degree(node_new) < k + 1:
                k = k - 1

        temp_r = random.random()
        if temp_r <= pi_sin:
            node_kind[node_new] = kinds[0]  # S
        elif temp_r <= pi_sin + pi_iin:
            node_kind[node_new] = kinds[2]  # I (保持原有的分布逻辑)
        else:
            node_kind[node_new] = kinds[3]  # R

        node_out[node_new] = 10000
        get_time_interval(node_new)


    def remove_node(min_node):
        global count_out
        neis = list(graph.neighbors(min_node))
        i = int(random.randint(0, len(neis) - 1))
        selected_one = neis[i]
        for j in neis:
            if j == selected_one:
                continue
            if j in graph.neighbors(selected_one):
                continue
            else:
                graph.add_edge(j, selected_one)
        graph.remove_node(min_node)
        node_out.pop(min_node)
        node_kind.pop(min_node)
        count_out += 1


    def initialize():
        for i in graph.nodes:
            node_kind[i] = kinds[0]
        for i in range(i0):
            i_first = random.choice(list(graph.nodes()))
            node_kind[i_first] = kinds[2]  # 初始感染者为I状态

        for i in graph.nodes:
            if node_kind[i] == kinds[0]:
                k = graph.degree(i)
                i_total = get_total_rate(kinds[2])
                node_out[i] = exp_distri(beta * k * i_total * len(graph.nodes))
            elif node_kind[i] == kinds[2]:
                node_out[i] = exp_distri(gamma)


    def simulation(t, graph):
        time_left = t
        count_in = 0
        count_out = 0

        last_time = time_total / 4
        next_node = exp_distri(lamb)

        while time_left > 0:
            min_node = 10000
            min_time = 1000
            for i in node_out.keys():
                if node_out[i] < min_time:
                    min_node = i
                    min_time = node_out[i]
            is_new_node = False
            if next_node < min_time:
                is_new_node = True

            for i in node_out.keys():
                if is_new_node:
                    node_out[i] -= next_node
                elif i == min_node:
                    node_out[i] = 0
                else:
                    node_out[i] -= min_time

            if not is_new_node:
                time2 = time_left - min_time
                next_node -= min_time
            else:
                time2 = time_left - next_node
                next_node = exp_distri(lamb)

            decrete_time = int(time_left)
            while decrete_time > time2:
                record_time()
                decrete_time -= step

            time_left = time2

            i_total = get_total_rate(kinds[2])
            if not is_new_node:
                if node_kind[min_node] == kinds[0]:  # S状态
                    k = len(graph.edges) * 2.0 / len(graph.nodes())
                    if random.random() <= beta * k * i_total * len(graph.nodes) / (
                            beta * k * i_total * len(graph.nodes) + pi_sout * miu):
                        node_kind[min_node] = kinds[1]  # S->E
                        get_time_interval(min_node)
                    else:
                        remove_node(min_node)

                elif node_kind[min_node] == kinds[1]:  # E状态
                    if random.random() <= sigma / (sigma + pi_iout * miu):
                        node_kind[min_node] = kinds[2]  # E->I
                        get_time_interval(min_node)
                    else:
                        remove_node(min_node)

                elif node_kind[min_node] == kinds[2]:  # I状态
                    if random.random() <= gamma / (gamma + pi_iout * miu):
                        node_kind[min_node] = kinds[3]  # I->R
                        get_time_interval(min_node)
                    else:
                        remove_node(min_node)

                elif node_kind[min_node] == kinds[3]:  # R状态
                    remove_node(min_node)  # R状态只能被移除

            if is_new_node:
                new_node()
                count_in += 1

            if time_left < last_time:
                for temp in range(len(kinds)):
                    Is = sum([1 for i in node_kind.keys() if node_kind[i] == kinds[temp]])
                    if Is in continous_time[kinds[temp]].keys():
                        continous_time[kinds[temp]][Is] += last_time - time_left
                    else:
                        continous_time[kinds[temp]][Is] = last_time - time_left
                last_time = time_left

        # print("in", count_in)
        # print("out", count_out)

    print("SEIR with deletion BEGINS")
    initialize()
    simulation(time_total, graph)
    return S_series, E_series, I_series, R_series, time_series
#plot_time_evolution()

m = 4
n = 1000
i0 = 1

step = 1
lamb = 1.5
#alpha = .4
beta = .004
gamma = 0.9
sigma = 0.7
miu = 0.025

pi_sin = 0.85
pi_iin = 0.10
pi_rin = 0.05

pi_sout = 0.90
pi_iout = 0.05
pi_rout = 0.05

df = pd.read_csv(r'D:\论文\project 1\code_stochastic\data2\regions\London\epidemic\cases.csv')
# df['date'] = pd.to_datetime(df['date'], format='%Y/%m/%d')
#real_infections = df.loc[(df['date']>=datetime(2020,3,11)) & (df['date']<=datetime(2020,4,30))]['value_daily']
#seven_day_avg = real_infections.rolling(window=7, min_periods=1).mean().values
real_infections = df.loc[(df['year_week']>='2020-44') & (df['year_week']<='2021-06')]['cases'].values[5:]/8796500

time_total = len(real_infections)-1

S_series, E_series, I_series, R_series, time_series = SEIR_deletion(m,
              n,
              i0,
              time_total,
              step,
              lamb,
              beta,
              gamma,
              sigma,
              miu,
              pi_sin,
              pi_iin,
              pi_rin,
              pi_sout,
              pi_iout,
              pi_rout)

print(len(time_series))
print(len(real_infections))


plt.figure(figsize=(10, 6))
# plt.plot(time_series, S_series, 'b-', label='Susceptible', linewidth=2)
# plt.plot(time_series, E_series, 'y-', label='Exposed', linewidth=2)
plt.plot(time_series, np.array(I_series)/1000, 'r-', label='Infected', linewidth=2)
plt.scatter(time_series, real_infections)
# plt.plot(time_series, R_series, 'g-', label='Recovered', linewidth=2)

plt.xlabel('Time Steps', fontsize=12)
plt.ylabel('Number of Individuals', fontsize=12)
plt.title('SEIR Model Time Evolution', fontsize=14)
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend(fontsize=10)
plt.tight_layout()
plt.show()