import networkx as nx
import math
import random
from random import choice
import time
import numpy as np
import matplotlib.pyplot as plt

m=4
n=1000
i0=1
time_left=400
step=1
lamb=1.5
alpha=.4
beta=.0025
gamma=.7
miu=0.025
#u=1

pi_sin=0.85
pi_iin=0.10
pi_rin=0.05

pi_sout=0.90
pi_iout=0.05
pi_rout=0.05

#f=math.exp(-np.random.poisson(u))

#retention=0
kinds=["S","I","R"]

global count_in
global count_out
count_in=0
count_out=0
continous_time={kinds[0]:{},kinds[1]:{},kinds[2]:{}}
#use WS network
graph=nx.generators.random_graphs.watts_strogatz_graph(n,m,.3)

time_seq=[]
node_kind={}

node_out= {}

I_degree={}
total_degree={}
k_i={}

#record time sequnce
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

print("BEGINS")

initialize()

print("node_degree ",graph.degree)
print("node_out ",node_out)
print("node_kind ",node_kind)

# CREATE NEW NODE
next_node=exp_distri(lamb)

print(time.asctime())
print("BEGIN")


decrete_time=int(time_left)

last_time=500

while time_left>0:
	# TIME PASS
	min_node=10000
	min_time=1000
	for i in node_out.keys():
		if node_out[i]<min_time:
			min_node=i
			min_time=node_out[i]
	is_new_node=False
	if next_node<min_time:
		is_new_node=True

	for i in node_out.keys():
		if is_new_node:
			node_out[i] -=next_node
		elif i==min_node :
			node_out[i]=0
		else:
			node_out[i]-=min_time

	if not is_new_node:
		time2=time_left-min_time
		next_node-=min_time
	else:
		time2=time_left-next_node
		next_node=exp_distri(lamb)

	if time_left<decrete_time:
		decrete_time=int(time_left)
	while decrete_time>time2:
		record_time()
		decrete_time-=step
		print("\r%d "%decrete_time,end=time.asctime())
	time_left=time2

    ###########要改##########################
	i_total=get_total_rate(kinds[1])
	if not is_new_node:
		if node_kind[min_node]==kinds[0]:
			k = len(graph.edges) * 2.0 / len(graph.nodes())
			if random.random()<=beta * k * i_total * len(graph.nodes)/(beta * k * i_total * len(graph.nodes)+pi_sout*miu):
				node_kind[min_node]=kinds[1]
				get_time_interval(min_node)
			else:
				remove_node(min_node)
		elif node_kind[min_node]==kinds[1]:
			if random.random()<=gamma/(gamma+pi_iout*miu):
				node_kind[min_node]=kinds[2]
				get_time_interval(min_node)
			else:
				remove_node(min_node)
		elif node_kind[min_node]==kinds[2]:
			if random.random()<=alpha/(alpha+pi_rout*miu):
				node_kind[min_node]=kinds[0]
				get_time_interval(min_node)
			else:
				remove_node(min_node)

	if is_new_node:
		new_node()
		count_in+=1

	if time_left < last_time:
		for temp in range(len(kinds)):
			Is = sum([1 for i in node_kind.keys() if node_kind[i] == kinds[temp]])
			if Is in continous_time[kinds[temp]].keys():
				continous_time[kinds[temp]][Is] += last_time - time_left
			else:
				continous_time[kinds[temp]][Is] = last_time - time_left
		last_time = time_left
	#when_stable()

sir_pi=[{} for i in range(len(kinds))]

for k in range(len(kinds)):
	sir_pi[k] = [0 for i in range(1 + max(continous_time[kinds[k]].keys()))]
	for i in continous_time[kinds[k]].keys():
		sir_pi[k][i]=continous_time[kinds[k]][i]

for k in range(len(kinds)):
	total_is=sum(sir_pi[k])
	sir_pi[k]=[j * 1.0 / total_is for j in sir_pi[k]]


# for k in range(len(kinds)):
# 	fig = plt.figure()
# 	plt.subplots_adjust(left=0.2, bottom=0.2, right=0.9)
# 	plt.scatter(range(len(sir_pi[k])), sir_pi[k], color="green")
# 	plt.xlabel("%s"%kinds[k],fontdict=font)
# 	plt.ylabel("P(%s)"%kinds[k],fontdict=font)
# 	plt.legend(frameon=False)
# 	with PdfPages('%d.pdf'%(j*3+k)) as pdf:
# 		     # As many times as you like, create a figure fig and save it:
# 		pdf.savefig(fig)


print(time_seq)
print(len(time_seq))
print(len(graph.nodes))
print("node_degree ",graph.degree)
print("node_out ",node_out)
print("node_kind ",node_kind)
print("in",count_in)
print("out",count_out)

import matplotlib.pyplot as plt

fig = plt.figure(figsize=(10, 6))
plt.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.12)

plt.plot(range(0,30*len(time_seq),30)[:5000], [i[0] for i in time_seq], "b-",linewidth=2.5,markeredgewidth=.01, label="S")
plt.plot(range(0,30*len(time_seq),30)[:5000], [i[1] for i in time_seq], "r-",linewidth=4.5,markeredgewidth=.5, label="I")
plt.plot(range(0,30*len(time_seq),30)[:5000], [i[2] for i in time_seq], 'g-',linewidth=3,markeredgewidth=.2, label="R")
plt.xlabel("t",fontsize=20)
plt.ylabel("N",fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.legend(fontsize=16,frameon=False,loc="upper right")

name='heritage_Nt,b'+str(beta)+',g'+str(gamma)+',a'+str(alpha)+'_'+str(lamb)+'_'+str(miu)+'t'+str(time_left)

plt.savefig("%s.pdf"%name)

import json
file=open("%s.json"%name,mode="w")
file.write(json.dumps(time_seq))
file.close()

plt.show()