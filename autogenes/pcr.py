'''
AutoGenes (c) University of Liverpool 2019

AutoGenes is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author: neilswainston
'''
# pylint: disable=invalid-name
# pylint: disable=too-few-public-methods
# pylint: disable=too-many-arguments
from autogenes.graph_writer import GraphWriter


class PcrWriter(GraphWriter):
    '''Class for generating PCR worklist graphs.'''

    def __init__(self, comps_vol, wt_primer_vol, mut_primer_vol, total_vol,
                 output_name):
        self._comps_vol = comps_vol
        self.__wt_primer_vol = wt_primer_vol
        self.__mut_primer_vol = mut_primer_vol
        self.__total_vol = total_vol
        GraphWriter.__init__(self, output_name)

    def _add_pcr(self, pcr_id, pcr_comps_ids, primer_ids):
        '''Add PCR reaction to worklist graph.'''
        mm_vol = self.__total_vol
        pcr = self._add_vertex(pcr_id, {'is_reagent': False})

        for pcr_comps_id in pcr_comps_ids:
            pcr_comps = self._add_vertex(pcr_comps_id,
                                         {'is_reagent': False})

            self._add_edge(pcr_comps, pcr, {'Volume': self._comps_vol})
            mm_vol -= self._comps_vol

        # Add outer oligos:
        for primer_id in primer_ids:
            primer = self._add_vertex(primer_id[0],
                                      {'is_reagent': False})

            primer_vol = self.__mut_primer_vol if primer_id[1] \
                else self.__wt_primer_vol

            self._add_edge(primer, pcr, {'Volume': primer_vol})

            mm_vol -= primer_vol

        mm = self._add_vertex('mm', {'is_reagent': True})
        self._add_edge(mm, pcr, {'Volume': mm_vol})
