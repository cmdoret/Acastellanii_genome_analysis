import sys
import cooler
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import markov_clustering as mc
import networkx as nx

clr = cooler.Cooler(sys.argv[1])

chroms = clr.chroms()[:]

inter = np.zeros((chroms.shape[0], chroms.shape[0]))

for i1, c1 in enumerate(chroms['name']):
    for i2, c2 in enumerate(chroms['name']):
        if c1 != c2:
            
            inter[i1, i2] = clr.matrix(balance=False, sparse=True).fetch(c1, c2).mean()

plt.imshow(inter)
plt.show()

markov_inter = mc.run_mcl(inter)
clusters = mc.get_clusters(markov_inter)
mc.draw_graph(markov_inter, clusters, labels={num: name for num, name in enumerate(chroms['name'])})
plt.show()
