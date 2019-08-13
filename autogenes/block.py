'''
AutoGenes (c) University of Liverpool 2019

AutoGenes is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author: neilswainston
'''
# pylint: disable=too-few-public-methods
# pylint: disable=too-many-arguments
from collections import Counter

from autogenes import get_dil_oligo_id, get_block_id, get_pos_muts
from autogenes.graph_writer import GraphWriter
from autogenes.pcr import PcrWriter


class InnerBlockPoolWriter(GraphWriter):
    '''Class for generating pooled inner block worklist graphs.'''

    def __init__(self, designs, wt_oligo_vol, mut_oligo_vol, output_name):
        self.__designs = designs
        self.__wt_oligo_vol = wt_oligo_vol
        self.__mut_oligo_vol = mut_oligo_vol
        GraphWriter.__init__(self, output_name)

    def _initialise(self):
        block_ids = []

        for design in self.__designs:
            for block_idx, block in enumerate(design):
                block_id = get_block_id(block_idx, block)

                if block_id not in block_ids:
                    inner_pool_id = block_id + '_ib'

                    inner_pool = self._add_vertex(inner_pool_id,
                                                  {'is_reagent': False})

                    # Pool *inner* oligos:
                    for oligo_id in block[1:-1]:
                        dil_id, is_mut = get_dil_oligo_id(oligo_id)
                        oligo = self._add_vertex(dil_id,
                                                 {'is_reagent': False})

                        self._add_edge(oligo, inner_pool,
                                       {'Volume': self.__mut_oligo_vol
                                        if is_mut else self.__wt_oligo_vol})

                    block_ids.append(block_id)


class BlockPcrWriter(PcrWriter):
    '''Class for generating block PCR worklist graphs.'''

    def __init__(self, designs, comps_vol, wt_primer_vol, mut_primer_vol,
                 total_vol, output_name):
        self.__designs = designs
        PcrWriter.__init__(self, comps_vol, wt_primer_vol, mut_primer_vol,
                           total_vol, output_name)

    def _initialise(self):
        block_ids = []

        for design in self.__designs:
            for block_idx, block in enumerate(design):
                block_id = get_block_id(block_idx, block)

                if block_id not in block_ids:
                    pcr_comps_ids = [block_id + '_ib']
                    pcr_id = block_id + '_b'
                    primer_ids = [get_dil_oligo_id(block[idx])
                                  for idx in [0, -1]]

                    self._add_pcr(pcr_id, pcr_comps_ids, primer_ids)

                block_ids.append(block_id)


class BlockPoolWriter(GraphWriter):
    '''Class for generating pooled block worklist graphs.'''

    def __init__(self, designs, min_vol, max_vol, is_mut, output_name):
        self.__designs = designs
        self.__min_vol = min_vol
        self.__max_vol = max_vol
        self.__is_mut = is_mut
        GraphWriter.__init__(self, output_name)

    def _initialise(self):
        pool_steps = {}

        for design in self.__designs:
            for block_idx, block in enumerate(design):
                block_id = get_block_id(block_idx, block)
                pcr_id = block_id + '_b'
                pool_id = '_'.join([str(val) if val > 0 else 'wt'
                                    for val in get_pos_muts(block_id)]) + '_p'

                pool_steps[pcr_id] = pool_id

        pool_counter = Counter(list(pool_steps.values()))
        max_pool_vol = pool_counter.most_common(1)[0][1] * self.__min_vol

        for pcr_id, pool_id in pool_steps.items():
            if (self.__is_mut and 'wt' not in pool_id) \
                    or (not self.__is_mut and 'wt' in pool_id):
                pcr = self._add_vertex(pcr_id, {'is_reagent': False})
                pool = self._add_vertex(pool_id, {'is_reagent': False})
                vol = min(self.__max_vol,
                          (1 / pool_counter[pool_id] * max_pool_vol))
                self._add_edge(pcr, pool, {'Volume': vol})
