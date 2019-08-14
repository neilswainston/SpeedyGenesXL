'''
AutoGenes (c) University of Liverpool 2019

AutoGenes is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author: neilswainston
'''
# pylint: disable=invalid-name
from synbiochem.utils.graph_utils import add_edge, add_vertex, Graph


class GraphWriter():
    '''GraphWriter base class.'''

    def __init__(self, output_name):
        self.__output_name = output_name
        self.__graph = Graph()
        self.__initialised = False

    def get_graph(self):
        '''Gets graph.'''
        if not self.__initialised:
            self._initialise()
            self.__initialised = True

        return self.__graph

    def get_output_name(self):
        '''Gets output name.'''
        return self.__output_name

    def _initialise(self):
        '''Initialise graph.'''

    def _add_vertex(self, name, kwds):
        '''Add vertex to graph.'''
        return add_vertex(self.__graph, name, kwds)

    def _add_edge(self, vertex_from, vertex_to, kwds):
        '''Add edge to graph.'''
        return add_edge(self.__graph, vertex_from, vertex_to, kwds)
