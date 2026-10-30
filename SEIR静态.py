import networkx as nx
import math
import random
from random import choice
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np

# 保持原有参数
m = 4
n = 1000
i0 = 1
time_total = 100

step = 1
lamb = 1.5
alpha = .4
beta = .05
gamma = .7
sigma = 0.2  # 新增参数：E->I的转化率
miu = 0.025

pi_sin = 0.85
pi_iin = 0.10
pi_rin = 0.05

pi_sout = 0.90
pi_iout = 0.05
pi_rout = 0.05

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

            elif node_kind[min_node] == kinds[1]:  # E状态
                if random.random() <= sigma:
                    node_kind[min_node] = kinds[2]  # E->I
                    get_time_interval(min_node)

            elif node_kind[min_node] == kinds[2]:  # I状态
                if random.random() <= gamma / (gamma + pi_iout * miu):
                    node_kind[min_node] = kinds[3]  # I->R
                    get_time_interval(min_node)

        if time_left < last_time:
            for temp in range(len(kinds)):
                Is = sum([1 for i in node_kind.keys() if node_kind[i] == kinds[temp]])
                if Is in continous_time[kinds[temp]].keys():
                    continous_time[kinds[temp]][Is] += last_time - time_left
                else:
                    continous_time[kinds[temp]][Is] = last_time - time_left
            last_time = time_left


def plot_time_evolution():
    plt.figure(figsize=(10, 6))
    plt.plot(time_series, S_series, 'b-', label='Susceptible', linewidth=2)
    plt.plot(time_series, E_series, 'y-', label='Exposed', linewidth=2)
    plt.plot(time_series, I_series, 'r-', label='Infected', linewidth=2)
    plt.plot(time_series, R_series, 'g-', label='Recovered', linewidth=2)

    plt.xlabel('Time Steps', fontsize=12)
    plt.ylabel('Number of Individuals', fontsize=12)
    plt.title('SEIR Model Time Evolution', fontsize=14)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend(fontsize=10)
    plt.tight_layout()
    plt.show()

# 运行模拟
print("BEGINS")
initialize()
simulation(time_total, graph)
plot_time_evolution()
# 绘图和统计分析
sir_pi = [{} for i in range(20)]

for k in range(len(kinds)):
    sir_pi[k][0] = [0 for i in range(1 + max(continous_time[kinds[k]].keys()))]
    for i in continous_time[kinds[k]].keys():
        sir_pi[k][0][i] = continous_time[kinds[k]][i]

for k in range(len(kinds)):
    for i in sir_pi[k].keys():
        total_is = sum(sir_pi[k][i])
        sir_pi[k][i] = [j * 1.0 / total_is for j in sir_pi[k][i]]

print(continous_time)
S = []
E = []
I = []
R = []
PS = []
PE = []
PI = []
PR = []

for i in range(len(sir_pi[0][0])):
    if sir_pi[0][0][i] != 0:
        S.append(i)
        PS.append(sir_pi[0][0][i])

for i in range(len(sir_pi[1][0])):
    if sir_pi[1][0][i] != 0:
        E.append(i)
        PE.append(sir_pi[1][0][i])

for i in range(len(sir_pi[2][0])):
    if sir_pi[2][0][i] != 0:
        I.append(i)
        PI.append(sir_pi[2][0][i])

for i in range(len(sir_pi[3][0])):
    if sir_pi[3][0][i] != 0:
        R.append(i)
        PR.append(sir_pi[3][0][i])

# 计算期望值
ES = sum([S[i] * PS[i] for i in range(0, len(S))])
EE = sum([E[i] * PE[i] for i in range(0, len(E))])
EI = sum([I[i] * PI[i] for i in range(0, len(I))])
ER = sum([R[i] * PR[i] for i in range(0, len(R))])

print("ES:", ES)
print("EE:", EE)
print("EI:", EI)
print("ER:", ER)

# 绘制最终结果图
plt.figure(figsize=(10, 6))
plt.scatter(S, PS, marker="o", color="blue", label="S")
plt.scatter(E, PE, marker="s", color="yellow", label="E")
plt.scatter(I, PI, marker="+", color="red", label="I")
plt.scatter(R, PR, marker="^", color="green", label="R")
plt.xlabel("N")
plt.ylabel("P(N)")
plt.legend(frameon=False)
plt.show()