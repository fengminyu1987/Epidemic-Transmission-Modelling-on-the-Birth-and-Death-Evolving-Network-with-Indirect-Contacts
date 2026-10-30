import networkx as nx
import math
import random
from random import choice

from matplotlib.backends.backend_pdf import PdfPages

m=4
n=1000
i0=1
time_total=4

step=1
lamb=1.5
alpha=.4
beta=.05
gamma=.7
miu=0.025

pi_sin=0.85
pi_iin=0.10
pi_rin=0.05

pi_sout=0.90
pi_iout=0.05
pi_rout=0.05

global count_in
global count_out
count_in=0
count_out=0

kinds=["S","I","R"]
status=[0]

continous_time = {kinds[0]: {}, kinds[1]: {}, kinds[2]: {}}

#use WS network
graph=nx.generators.random_graphs.watts_strogatz_graph(n,m,.7)

node_kind={}

node_out= {}

time_seq=[]


I_degree={}
total_degree={}
k_i={}

def record_time():
	total_node=len(graph.nodes)
	time_seq.append((get_total_rate(kinds[0])*total_node,get_total_rate(kinds[1])*total_node,get_total_rate(kinds[2])*total_node))

# exponential distribution
def exp_distri(lamb):
	return math.log(1 - random.random()) * (-1 / lamb)

# the rate of neighbor of special kind
def get_rate(target,kind):
	neighbors=graph.neighbors(target)
	total_n=len(neighbors)
	kind_n=0
	for i in neighbors:
		if node_kind[i]==kind:
			kind_n+=1
	return 1.0*kind_n/total_n

# the rate of special kind
def get_total_rate(kind):
	total_n=len(graph.nodes)
	kind_n=0
	for i in graph.nodes:
		if node_kind[i]==kind:
			kind_n+=1
	return 1.0*kind_n/total_n

############要改#####################
# the stay time of node
def get_time_interval(node):
	i_total = get_total_rate(kinds[1])
	s_total = get_total_rate(kinds[0])
	r_total = get_total_rate(kinds[2])
	if node_kind[node]==kinds[0]:
		k = len(graph.edges)*2.0/len(graph.nodes())
		rate_s_serv=beta * k * i_total * len(graph.nodes)+pi_sout*miu################注意问一下！！！len(graph.nodes)怎么回事##
		if rate_s_serv>0:
			temp = exp_distri(rate_s_serv)
		else:temp =1000
		if temp<node_out[node] or node_out[node]==0:
			node_out[node]=temp

	if node_kind[node]==kinds[1]:
		rate_i_serv=(gamma+pi_iout*miu)
		temp = exp_distri(rate_i_serv)
		if temp<node_out[node] or node_out[node]==0:
			node_out[node]=temp
	if node_kind[node]==kinds[2]:
		rate_r_serv=(alpha+pi_rout*miu)
		temp = exp_distri(rate_r_serv)
		if temp<node_out[node] or node_out[node]==0:
			node_out[node]=temp

# create a new node and add it to the graph
def new_node():
	node_new=0
	while graph.has_node(node_new):
		node_new+=1
	k=0
	while k<m:
		dead = choice(list(graph.nodes()))
		k+=1
		if dead!=node_new:
			graph.add_edge(dead,node_new)
		# 如果连接重复,重roll
		if graph.degree(node_new) < k + 1:
			k = k - 1
	temp_r=random.random()
	if temp_r<=pi_sin:
		node_kind[node_new]=kinds[0]
	elif temp_r<=pi_sin+pi_iin:
		node_kind[node_new] = kinds[1]
	else:
		node_kind[node_new] = kinds[2]
	node_out[node_new]=10000
	get_time_interval(node_new)

def remove_node(min_node):
	global count_out
	neis=list(graph.neighbors(min_node))
	i=int(random.randint(0,len(neis)-1))
	selected_one=neis[i]#要继承的那个节点
	for j in neis:
		if j==selected_one:
			continue
		if j in graph.neighbors(selected_one):
			continue
		else:
			graph.add_edge(j,selected_one)
	graph.remove_node(min_node)
	node_out.pop(min_node)
	node_kind.pop(min_node)
	count_out += 1



# initialize the graph
def initialize():
	for i in graph.nodes:
		node_kind[i]=kinds[0]
	#initial i
	for i in range(i0):
		i_first = random.choice(list(graph.nodes()))
		node_kind[i_first]=kinds[1]
	s_total=get_total_rate(kinds[0])
	i_total = get_total_rate(kinds[1])
	total=len(graph.nodes)
	#initial time
	for i in graph.nodes:
		if node_kind[i]==kinds[0]:
			k=graph.degree(i)
			node_out[i]=exp_distri(beta*k*i_total*len(graph.nodes))
		elif node_kind[i]==kinds[1]:
			node_out[i] = exp_distri(gamma)##############i_total*gamma? 好像不用 这里是对每个节点的out


# when the simulation reach the stable situation

# simulation
def simulation(t,graph):
	time_left=t
	#decrete_time = int(time_left)
	count_in = 0
	count_out = 0

	last_time=time_total/4
	# CREATE NEW NODE
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
			#record_time()
			decrete_time -= step
			#print("\r%d " % decrete_time, end=time.asctime())
		time_left = time2

		###########要改##########################
		i_total = get_total_rate(kinds[1])
		if not is_new_node:
			if node_kind[min_node] == kinds[0]:
				k = len(graph.edges) * 2.0 / len(graph.nodes())
				if random.random() <= beta * k * i_total * len(graph.nodes) / (
						beta * k * i_total * len(graph.nodes) + pi_sout * miu):
					node_kind[min_node] = kinds[1]
					get_time_interval(min_node)
				else:
					remove_node(min_node)
			elif node_kind[min_node] == kinds[1]:
				if random.random() <= gamma / (gamma + pi_iout * miu):
					node_kind[min_node] = kinds[2]
					get_time_interval(min_node)
				else:
					remove_node(min_node)
			elif node_kind[min_node] == kinds[2]:
				if random.random() <= alpha / (alpha + pi_rout * miu):
					node_kind[min_node] = kinds[0]
					get_time_interval(min_node)
				else:
					remove_node(min_node)

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
	# when_stable()

	print("in", count_in)
	print("out", count_out)

print("BEGINS")

initialize()

print("node_degree ",graph.degree)
print("node_out ",node_out)
print("node_kind ",node_kind)


sir_pi=[{} for i in range(20)]

status[0]=0
graph=nx.generators.random_graphs.watts_strogatz_graph(n,m,.3)
initialize()
simulation(time_total,graph)
for k in range(len(kinds)):
	sir_pi[k][0] = [0 for i in range(1 + max(continous_time[kinds[k]].keys()))]
	for i in continous_time[kinds[k]].keys():
		sir_pi[k][0][i]=continous_time[kinds[k]][i]

#continous_time={kinds[0]:{},kinds[1]:{},kinds[2]:{}}


print(len(graph.nodes))
print("node_degree ",graph.degree)
print("node_out ",node_out)
print("node_kind ",node_kind)


import matplotlib.pyplot as plt
import numpy as np

for k in range(len(kinds)):
	for i in sir_pi[k].keys():
		total_is=sum(sir_pi[k][i])
		sir_pi[k][i]=[j * 1.0 / total_is for j in sir_pi[k][i]]


font = {'family': 'serif',
        'weight': 'normal',
        'size': 18,
        }


print(continous_time)
S=[]
I=[]
R=[]
PS=[]
PI=[]
PR=[]
for i in range(len(sir_pi[0][0])):
	if sir_pi[0][0][i]!=0:
		S.append(i)
		PS.append(sir_pi[0][0][i])

for i in range(len(sir_pi[1][0])):
	if sir_pi[1][0][i]!=0:
		I.append(i)
		PI.append(sir_pi[1][0][i])

for i in range(len(sir_pi[2][0])):
	if sir_pi[2][0][i]!=0:
		R.append(i)
		PR.append(sir_pi[2][0][i])

ES_list=[S[i]*PS[i] for i in range(0,len(S))]
ES=sum(ES_list)

EI_list=[I[i]*PI[i] for i in range(0,len(I))]
EI=sum(EI_list)

ER_list=[R[i]*PR[i] for i in range(0,len(R))]
ER=sum(ER_list)

print("ES:", ES)
print("EI:", EI)
print("ER:", ER)


fig = plt.figure()
plt.subplots_adjust(left=0.2, bottom=0.2, right=0.9)
plt.scatter(S, PS, marker="o", color="blue",label="S")
plt.scatter(I, PI, marker="+",color="red",label="I")
plt.scatter(R, PR, marker="^",color="green",label="R")
#plt.scatter(range(0,len(sir_pi[1][0])), sir_pi[1][0], marker="+",color="red")
#plt.scatter(range(0,len(sir_pi[2][0])), sir_pi[2][0], marker="^",color="green")
plt.xlabel("N")
plt.ylabel("P(N)")
plt.legend(frameon=False)
plt.show()


'''
plt.xlabel("I")
plt.ylabel("P(I)")
plt.legend(frameon=False)
plt.show()




plt.xlabel("R")
plt.ylabel("P(R)")
plt.legend(frameon=False)
plt.show()
'''

#
# fig = plt.figure()
# plt.subplots_adjust(left=0.2, bottom=0.2, right=0.9)
# plt.scatter(range(len(sir_pi[0])), sir_pi[0], color="green", label=r"$\gamma$=%.2f" % gammas[0])
# plt.scatter(range(len(sir_pi[3])), sir_pi[3], color="blue", marker="^", label=r"$\gamma$=%.2f" % gammas[1])
# plt.scatter(range(len(sir_pi[4])), sir_pi[4], color="red", marker="+", label=r"$\gamma$=%.2f" % gammas[2])
# plt.xlabel("R",fontdict=font)
# plt.ylabel("P(R)",fontdict=font)
# plt.legend(frameon=False)
# with PdfPages('8.pdf') as pdf:
#      # As many times as you like, create a figure fig and save it:
#      pdf.savefig(fig)
#
#
# fig = plt.figure()
# plt.subplots_adjust(left=0.2, bottom=0.2, right=0.9)
# plt.scatter(range(len(sir_pi[5])), sir_pi[5], color="green", label=r"$\lambda$=%.2f" % lambdas[0])
# plt.scatter(range(len(sir_pi[6])), sir_pi[6], color="blue", marker="^", label=r"$\lambda$=%.2f" % lambdas[1])
# plt.scatter(range(len(sir_pi[7])), sir_pi[7], color="red", marker="+", label=r"$\lambda$=%.2f" % lambdas[2])
# plt.xlabel("R",fontdict=font)
# plt.ylabel("P(R)",fontdict=font)
# plt.legend(frameon=False)
# with PdfPages('9.pdf') as pdf:
#      # As many times as you like, create a figure fig and save it:
#      pdf.savefig(fig)