"""
Calculating network metrics

Useful sources:
  * https://www.python-course.eu/graphs_python.php

@author: Scott Campit
"""

import networkx as nx
import xml.etree.cElementTree as ET

# First construct a dictionary, where each unique biomass metabolite is the key, and the values are a list of reactions associated with the reaction. The keys will be the nodes, and the edges will be the reactions

path = r'/mnt/c/Users/scampit/Desktop/eMetOncoFit/models/''
model = 'RECON1.xml'
gssm = {}


gem = nx.Graph()
