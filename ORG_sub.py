import numpy as np
import math
import random
import pandas as pd
from GraphRicciCurvature.OllivierRicci import OllivierRicci
import networkx as nx
import networkit as nk
from typing import Union
from typing import List, Tuple
from littleballoffur.sampler import Sampler
from littleballoffur.backend import NetworKitBackEnd, NetworkXBackEnd


class ORGSampler(Sampler):
    r"""An implementation of graph expoloration subsampling by ORG-sub1 algorithm.
    The random walker calculates the edges' OR curvature adptaively, then goes to
    to the next edge with the most different OR curvature. The sampled graph is always connected.  
    For details about the algorithm, please see the Algorithm 1 in manuscript.
    Args:
        number_of_nodes (int): Number of nodes. Default is 100.
        seed (int): Random seed. Default is 42.
    """

    def __init__(self, number_of_nodes: int = 100,seed: int = 42,alpha=0.5):
        self.number_of_nodes = number_of_nodes
        self.seed = seed
        self._set_seed()
        self.alpha=alpha
        
    def _create_initial_node_set(self, graph, start_edge):
        """
        Choosing an initial node.
        """
        self._edges = list(graph.edges())
        self._edges_np = np.array(self._edges)
        self._orc=OllivierRicci(graph, alpha=alpha, verbose="TRACE")
        self._targets = list() 
        self._current_edge = random.choice(
            range(self.backend.get_number_of_edges(graph))
        )

        self._sampled_edges=set([self._current_edge])
        self._current_curvature = list(self._orc.compute_ricci_curvature_edges([(self._edges[self._current_edge][0],self._edges[self._current_edge][1])]).values())[0]
        self._sampled_nodes = set([self._edges[self._current_edge][0]])
        self._sampled_nodes.add(self._edges[self._current_edge][1])

        self._targets1 = np.union1d(np.where(self._edges_np[:,0] == self._edges[self._current_edge][0]),np.where(self._edges_np[:,1] == self._edges[self._current_edge][0])) # [i for i in range(len(edges))]
        self._targets2 = np.union1d(np.where(self._edges_np[:,0] == self._edges[self._current_edge][1]),np.where(self._edges_np[:,1] == self._edges[self._current_edge][1]))
        self._targets = np.union1d(self._targets1,self._targets2)
        self._targets = list(set(self._targets.tolist()).difference(self._sampled_edges))
        self._pre_curvature = self._current_curvature

            
    def _do_a_step(self, graph):
        """
        Doing a single random walk step.
        """

        if len(self._targets)==0:
            new_edge=random.choice(list(set([i for i in range(len(self._edges))]).difference(self._sampled_edges)))
            self._current_edge = new_edge
            self._current_curvature = list(self._orc.compute_ricci_curvature_edges([(edges[self._current_edge][0],edges[self._current_edge][1])]).values())[0]
        else:
            self._target_edges=[self._edges[i] for i in self._targets]
            self._target_curvatures=list(self._orc.compute_ricci_curvature_edges(self._target_edges).values())
            self._edgecurvatures=np.column_stack((self._targets, self._target_curvatures))

            self._curvatures = abs(self._edgecurvatures[:,1]-self._pre_curvature)
            newidx=int(np.where(self._curvatures==np.amax(self._curvatures))[0][0])
            new_edge = int(self._edgecurvatures[newidx,0])
            self._current_edge = new_edge
            self._current_curvature = self._edgecurvatures[newidx,1]
        self._sampled_edges.add(self._current_edge)
        
        self._sampled_nodes.add(self._edges[self._current_edge][0])
        self._sampled_nodes.add(self._edges[self._current_edge][1])

        self._targets1 = np.union1d(np.where(self._edges_np[:,0] == self._edges[self._current_edge][0]),np.where(self._edges_np[:,1] == self._edges[self._current_edge][0])) # [i for i in range(len(edges))]
        self._targets2 = np.union1d(np.where(self._edges_np[:,0] == self._edges[self._current_edge][1]),np.where(self._edges_np[:,1] == self._edges[self._current_edge][1]))
        self._targets = np.union1d(self._targets1,self._targets2)
        self._targets = list(set(self._targets.tolist()).difference(self._sampled_edges))
        
        self._pre_curvature = self._current_curvature

    def sample(
        self, graph, start_edge= None) :
        """
        Sampling nodes by the proposed ORG-sub1 algorithm.
        Arg types:
            * **graph** *(NetworkX or NetworKit graph)* - The graph to be sampled from.
            * **start_node** *(int, optional)* - The start node.
        Return types:
            * **new_graph** *(NetworkX or NetworKit graph)* - The graph of sampled edges.
        """
        self._deploy_backend(graph)
        self._check_number_of_nodes(graph)
        self._create_initial_node_set(graph, start_edge)

        while len(self._sampled_nodes) < self.number_of_nodes:
            print("A new step!")
            self._do_a_step(graph)
        new_graph = self.backend.get_subgraph(graph, self._sampled_nodes)
        return new_graph

############################## An example ################################

seed=20
ncom1=20
ncom2=20
ncom3=20
ncom4=20
sizes = [ncom1, ncom2, ncom3, ncom4]
pin=0.8
pout=0.1
prop=0.1
alpha=0.5
probs = [[pin, pout, pout, pout], [pout, pin, pout, pout], [pout,pout,pin,pout],[pout,pout,pout,pin]]

######################## Generate a SBM model ##################

G = nx.stochastic_block_model(sizes, probs, seed=seed)
number_of_nodes=int(prop*G.number_of_nodes())

sampler = ORGSampler(number_of_nodes = number_of_nodes,seed=seed,alpha=alpha)
sub_G = sampler.sample(G)
