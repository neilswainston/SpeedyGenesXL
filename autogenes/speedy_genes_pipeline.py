'''
AutoGenes (c) University of Liverpool 2019

AutoGenes is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author: neilswainston
'''
from collections import defaultdict
import itertools
import os
import sys
from time import gmtime, strftime

from synbiochem import utils

from autogenes import pipeline, worklist
from autogenes.block import InnerBlockPoolWriter, \
    BlockPcrWriter, BlockPoolWriter
from autogenes.dilution import WtOligoDilutionWriter
from autogenes.gene import CombiGenePcrWriter
from autogenes.pool import MutOligoPoolWriter


def run(plate_dir, max_mutated, n_blocks, out_dir_parent, exp_name):
    '''run method.'''
    assert len(exp_name) < 6

    dte = strftime("%y%m%d", gmtime())

    input_plates = pipeline.get_input_plates(plate_dir)
    oligos, mutant_oligos, primers = _read_plates(input_plates)
    designs = _combine(oligos, mutant_oligos, max_mutated, n_blocks)

    writers = [
        WtOligoDilutionWriter(oligos + primers, designs, 20, 20, 200,
                              exp_name + '-wt-dil'),
        MutOligoPoolWriter(mutant_oligos, 10, exp_name + '-mut-pl'),
        InnerBlockPoolWriter(designs, 2.5, 5, exp_name + '-templ'),
        BlockPcrWriter(designs, 1.2, 1.5, 3, 25, exp_name + '-pcr1'),
        BlockPoolWriter(designs, 2, 25, False, exp_name + '-wt-bk'),
        BlockPoolWriter(designs, 2, 25, True, exp_name + '-mut-bk'),
        CombiGenePcrWriter(designs, 4, 1.5, 1.5, 3, 25,
                           [['5-primer_dil', False], ['28_dil', False]],
                           exp_name + '-pcr2')

    ]

    out_dir_name = os.path.join(out_dir_parent, dte + exp_name)

    pipeline.run(writers, input_plates,
                 parent_out_dir_name=out_dir_name)

    worklist.format_worklist(out_dir_name)


def _read_plates(input_plates):
    '''Read plates.'''
    oligos = utils.sort(
        [obj['id'] for obj in input_plates['wt'].get_all().values()
         if obj['id'].isdigit()])

    primers = [obj['id'] for obj in input_plates['wt'].get_all().values()
               if not obj['id'].isdigit()]

    mutant_oligos = defaultdict(list)

    for obj in input_plates['mut'].get_all().values():
        mutant_oligos[obj['parent']].append(obj['id'])

    return oligos, mutant_oligos, primers


def _combine(oligos, mutant_oligos, max_mutated, n_blocks):
    '''Design combinatorial assembly.'''

    # Assertion sanity checks:
    assert len(oligos) % 2 == 0
    assert len(oligos) / n_blocks >= 2
    assert mutant_oligos if max_mutated > 0 else True

    designs = []

    # Get combinations:
    for n_mutated in range(max_mutated + 1):
        designs.extend(_get_combis(oligos, mutant_oligos, n_mutated, n_blocks))

    return designs


def _get_combis(oligos, mutant_oligos, n_mutated, n_blocks):
    '''Get combinations.'''
    designs = []

    for combi in itertools.combinations(list(mutant_oligos), n_mutated):
        design = list(oligos)

        for wt_id in combi:
            design[design.index(wt_id)] = wt_id + 'm'

        block_lengths = [0] * n_blocks

        for idx in itertools.cycle(range(0, n_blocks)):
            block_lengths[idx] = block_lengths[idx] + 2

            if sum(block_lengths) == len(design):
                break

        idx = 0
        blocks = []

        for val in block_lengths:
            blocks.append(design[idx: idx + val])
            idx = idx + val

        designs.append(blocks)

    return designs


def main(args):
    '''main method.'''
    import cProfile

    prf = cProfile.Profile()
    prf.enable()

    run(args[0], int(args[1]), int(args[2]), args[3], args[4])

    prf.disable()

    prf.print_stats(sort='cumtime')


if __name__ == '__main__':
    main(sys.argv[1:])
