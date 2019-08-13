'''
AutoGenes (c) University of Liverpool 2019

AutoGenes is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author: neilswainston
'''


def get_dil_oligo_id(oligo_id):
    '''Get diluted oligo id.'''
    is_mut = oligo_id[-1] == 'm'
    return oligo_id if is_mut else oligo_id + '_dil', is_mut


def get_block_id(block_idx, block):
    '''Get block id.'''
    mutations = [oligo[:-1] for oligo in block if oligo[-1] == 'm']

    return str(block_idx + 1) + '_' + \
        ('&'.join(mutations) if mutations else 'wt')


def get_design_id(design):
    '''Get design id.'''
    return '-'.join([get_block_id(block_idx, block)
                     for block_idx, block in enumerate(design)])


def get_pos_muts(block_id):
    '''Parse block id to get position and number of mutations.'''
    tokens = block_id.split('_')
    return int(tokens[0]), \
        0 if tokens[1] == 'wt' else (tokens[1].count('&') + 1)


def get_primers(designs):
    '''Get primers.'''
    primers = [[block[0], block[-1]]
               for design in designs
               for block in design]

    return list({primer for pair in primers for primer in pair})
