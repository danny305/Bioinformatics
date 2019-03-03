
import pydot
from pprint import pprint
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
#%matplotlib inline





#### Markov Model ####
"""Describes a stochastic process where the assumed probability of future state(s) depends only on the current 
process state and not on any of the preceding states."""

#create state space and initial state probabilities

states = ['sleeping', 'eating', 'pooping']
pi = [0.35, 0.35, 0.3]
state_space = pd.Series(pi, index=states, name='states')
print(state_space)
print(state_space.sum())


# create transition matrix: transition probability of changing states given a state
# matrix is size (M x M) where M is number of states.

q_df = pd.DataFrame(columns=states, index=states)
q_df.loc[states[0]] = [0.4, 0.2, 0.4]
q_df.loc[states[1]] = [0.45, 0.45, 0.1]
q_df.loc[states[2]] = [0.45, 0.25, 0.3]


print(q_df)

print(q_df.sum(axis=1))
#
#creates a dictionary with a tuple as a key and the value being the transition between the 2 states in the tuple
def _get_markov_edges(Q):
    edges = {(idx,col):Q.loc[idx, col] for col in Q.columns for idx in Q.index}
    return edges
#
# edges_wts = _get_markov_edges(q_df)
#
# pprint(edges_wts)
#
#
# #create graph object
# G = nx.MultiDiGraph()
#
# #nodes correspond to states
# G.add_nodes_from(states)
# print(f'Nodes:\n{G.nodes()}\n')
#
#
# #edges represent the transition probability
#
# for (frm,to),prob in edges_wts.items():
#     G.add_edge(frm, to, weight=prob, label=prob)
#
# pprint(G.edges(data=True))
#
# pos = nx.drawing.nx_pydot.graphviz_layout(G, prog='dot')
# nx.draw_networkx(G, pos)
#
#
# nx.drawing.nx_pydot.write_dot(G,'pet_dog_markov.dot')
#
# (graph,) = pydot.graph_from_dot_file('pet_dog_markov.dot')
# graph.write_png('./class_hw/Bioinformatics/pet_dog_markov.png')






#### Hidden Markov Model ####

# create state space and initial state probabilites
# we assume they are equiprobable.

hidden_states = ['healthy', 'sick']
pi = [0.5, 0.5]
state_space = pd.Series(pi, index=hidden_states, name='states')

#create hidden transition matrix
#a or alpha = transition probability matrix of changing between states given a state
# matrix is size (M x M) where M is number of states

a_df = pd.DataFrame(columns=hidden_states, index=hidden_states)
a_df.loc[hidden_states[0]] = [0.7, 0.3]
a_df.loc[hidden_states[1]] = [0.4, 0.6]

a = a_df.values


print(a_df)

# create matrix of observable (emission) probabilities
# b or beta = observation probabilities given a state
# matrix size (M x O) M: # of hidden states O: # of different observations.

obs_states = ['sleeping', 'eating', 'pooping']

b_df = pd.DataFrame(columns=obs_states, index=hidden_states)
b_df.loc[hidden_states[0]] = [0.2, 0.6, 0.2]
b_df.loc[hidden_states[1]] = [0.4, 0.1, 0.5]

b = b_df.values

print(b_df)


# Create graph edges

hide_edges_wts = _get_markov_edges(a_df)
pprint(hide_edges_wts)

emit_edges_wts = _get_markov_edges(b_df)
pprint(emit_edges_wts)


# Create graph object
G = nx.MultiDiGraph()


# add nodes correspond to the hidden states
G.add_nodes_from(hidden_states)


# Add edges for transition probabilities between hidden probabilities
for (frm,to),prob in hide_edges_wts.items():
    G.add_edge(frm, to, weight=prob, label=prob)


# Add edges for emission probabilities for every hidden state
for (h_state,em_state), prob in emit_edges_wts.items():
    G.add_edge(h_state, em_state, weight=prob, label=prob)


pprint(G.edges(data=True))

# pos = nx.drawing.nx_pydot.graphviz_layout(G, prog='neato')
# nx.draw_networkx(G, pos)
#
#
# nx.drawing.nx_pydot.write_dot(G,'pet_dog_hm.dot')
#
# (graph,) = pydot.graph_from_dot_file('pet_dog_hm.dot')
# graph.write_png('pet_dog_hm.png')



# Now we will discern the health of your dog over time given a sequence of observation behaviors

# observations are encoded numerically

obs_map = {'sleeping': 0, 'eating': 1, 'pooping': 2}
obs_code = np.array([1, 1, 2, 1, 0, 1, 2, 1, 0, 2, 2, 0, 1, 0, 1])

inv_obs_map = {val: key for key,val in obs_map.items()}
obs_seq = [inv_obs_map[obs] for obs in obs_code]

print(pd.DataFrame(np.column_stack([obs_code, obs_seq]), columns=['Obs_Code', 'Obs_Seq']))


# define Viterbi algorithm for shortest path
# code adapted from Stephen Marsland's, Machine Learning An Algorithmic Perspective, Vol. 2
# https://github.com/alexsosn/MarslandMLAlgo/blob/master/Ch16/HMM.py


def viterbi(pi, a, b, obs):
    """pi = initial probabilites,
    a = hidden states transition matrix,
    b = emission state probability matrix given a hidden state,
    obs =  sequence of observations"""

    print('PI: ', pi)
    print(f'A: {a}')
    print('B\n',b) #(2,3) matrix
    nStates =np.shape(b)[0]  # 2 states
    print('nStates\n', nStates)
    T = np.shape(obs)[0]     # 15 sequence of observation length (15,)
    print('obs\n', obs)
    print('T\n', T)

    # initialize blank path
    path = np.zeros(T,dtype=int)  # 1-D array with 15 zeros
    print(f'path: {path}')

    # delta --> highest probability of any path that reaches state i
    delta = np.zeros((nStates, T))   # (2,15) 2-D array of all zeros
    print(f'delta: {delta}')

    # phi --> argmax by time step for each state
    phi = np.zeros((nStates, T))

    # initialize delta and phi
    print('pi: ', pi)
    print('b column 0: ', b[0,obs[0]])
    delta[:,0] = pi * b[:, obs[0]]
    print(f'delta initialized: {delta[:,0]}')
    phi[:, 0] = 0

    # forward algorithm
    print('\nStart Walk Forward\n')
    for time in range(1,T):
        for state in range(nStates):
            # max returns the actual largest value
            print(f'state: {state}, time {time}')
            print(np.max(delta[:, time-1] * a[:, state]) * b[state, obs[time]])
            """ 
            b multiplication is not included in the max calculation () bc its just a scalar value. The max is decided by 
            when you multiply the prev state vector by the transition vector a.
            """
            delta[state,time] = np.max(delta[:, time-1] * a[:, state]) * b[state, obs[time]]
            #
            phi[state, time] = np.argmax(delta[:, time-1] * a[:, state]) #* b[state, obs[time]]
            print(f's={state} and t={time}: phi[{state}, {time}] = {phi[state,time]}')



    # find optimal path
    print('-'*50)
    print(delta)
    print(phi)


    print('Start Backtrace\n')
    print(delta[:, T-1])
    """
    We initialize the backtrace with the hidden state with the biggest delta at the end (the state with the 
    highest probability at the last position) then given this state we move backwards through the 
    probability using the 
    """
    path[T-1] = np.argmax(delta[:,T-1])

    print(f'T: {T}')
    for t in range(T-2, -1, -1):
        print(f't+1: {t+1}')
        print('path[t+1]', path[t+1])
        print(f'phi[{path[t+1]},{t+1}]',phi[path[t+1],[t+1]])
        path[t] = phi[path[t+1], [t+1]]
        print(f'set path[{t}] = to {phi[path[t+1], [t+1]]}')
        print(f'path[{t}] = {path[t]}')

    return path, delta, phi

path, delta, phi = viterbi(pi,a, b, obs_code)

print('\nsingle best state path: \n', path)
print('\ndelta:\n', delta)
print('\nphi:\n', phi)