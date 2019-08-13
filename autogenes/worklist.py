'''
AutoGenes (c) University of Liverpool 2019

AutoGenes is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author: neilswainston
'''
# pylint: disable=invalid-name
# pylint: disable=not-an-iterable
# pylint: disable=too-few-public-methods
# pylint: disable=too-many-instance-attributes
# pylint: disable=too-many-locals
# pylint: disable=ungrouped-imports
# pylint: disable=unsubscriptable-object
# pylint: disable=wrong-import-order
from collections import defaultdict
from operator import itemgetter
import os
import re

from scipy.spatial.distance import cityblock
from synbiochem.utils.graph_utils import get_roots

from autogenes import plate, smart_sort_opt
import pandas as pd

_VALUES_RENAME = {('src_plate', 'dest_plate'):
                  {('reagents'): 'MastermixTrough'}}

_COLUMNS_RENAME = {'src_name': 'ComponentName',
                   'src_plate': 'SourcePlateBarcode',
                   'src_well': 'SourcePlateWell',
                   'dest_plate': 'DestinationPlateBarcode',
                   'dest_well': 'DestinationPlateWell'}

_COLUMNS_ORDER = ['Volume',
                  'SourcePlateBarcode',
                  'SourcePlateWell',
                  'DestinationPlateBarcode',
                  'DestinationPlateWell',
                  'ComponentName']


class WorklistGenerator():
    '''Class to generate worklists.'''

    def __init__(self, graph):
        self.__graph = graph
        self.__worklist = None
        self.__input_plates = {}
        self.__plate_names = {'reagents': 'reagents',
                              'output': 'output'}
        self.__added_comps = {}

    def get_worklist(self, input_plates=None, plate_names=None):
        '''Gets worklist and required plates.'''
        required_plates = {}

        if not self.__worklist:
            self.__create_worklist(input_plates, plate_names)

        worklists = []

        for dest_plate, worklist in self.__worklist.groupby('dest_plate'):
            worklist.name = dest_plate
            worklists.append(worklist)
            required_plates[dest_plate] = self.__input_plates[dest_plate]
            required_plates.update({nme: self.__input_plates[nme]
                                    for nme in worklist['src_plate'].unique()})

        return worklists, required_plates

    def __create_worklist(self, input_plates, plate_names):
        '''Creates worklist and plates.'''
        data = []

        if input_plates:
            self.__input_plates.update(input_plates)

        if plate_names:
            self.__plate_names.update(plate_names)

        for root in get_roots(self.__graph):
            self.__traverse(root, 0, data)

        self.__worklist = pd.DataFrame(data)

        self.__write_input_plates()
        self.__add_locations()

    def __write_input_plates(self):
        '''Writes input_plates from worklist.'''
        # Write input plate:
        if 'src_well_fixed' not in self.__worklist:
            self.__worklist['src_well_fixed'] = None
        else:
            self.__worklist['src_well_fixed'] = \
                self.__worklist['src_well_fixed'].where(
                    (pd.notnull(self.__worklist['src_well_fixed'])), None)

        if 'dest_well_fixed' not in self.__worklist:
            self.__worklist['dest_well_fixed'] = None
        else:
            self.__worklist['dest_well_fixed'] = \
                self.__worklist['dest_well_fixed'].where(
                    (pd.notnull(self.__worklist['dest_well_fixed'])), None)

        inpt = \
            self.__worklist.loc[self.__worklist['src_is_input']
                                ][['src_name', 'src_well_fixed']].values

        for val in inpt:
            self.__add_component(val[0],
                                 'input',
                                 False,
                                 val[1])

        # Write reagents plate:
        reags = \
            self.__worklist.loc[self.__worklist['src_is_reagent']
                                ][['src_name', 'src_well_fixed']].values

        for val in sorted(reags, key=itemgetter(0)):
            self.__add_component(val[0],
                                 self.__plate_names['reagents'],
                                 True,
                                 val[1])

        # Write intermediates:
        intrm = self.__worklist[~(self.__worklist['src_is_input']) &
                                ~(self.__worklist['src_is_reagent'])]

        for _, row in intrm.sort_values('level', ascending=False).iterrows():
            self.__add_component(row['src_name'],
                                 row['level'],
                                 False,
                                 row['src_well_fixed'])

        # Write products:
        for _, row in self.__worklist.iterrows():
            if row['level'] == 0:
                self.__add_component(row['dest_name'],
                                     self.__plate_names['output'],
                                     False,
                                     row['dest_well_fixed'])

    def __add_locations(self):
        '''Add locations to worklist.'''
        locations = self.__worklist.apply(lambda row: self.__get_location(
            row['src_name'], row['dest_name']), axis=1)

        loc_df = locations.apply(pd.Series)
        loc_df.index = self.__worklist.index
        loc_df.columns = ['src_plate',
                          'src_well',
                          'src_idx',
                          'src_row',
                          'src_col',
                          'src_plate_size',
                          'src_pipette_idx',
                          'dest_plate',
                          'dest_well',
                          'dest_idx',
                          'dest_row',
                          'dest_col',
                          'dest_plate_size',
                          'dest_pipette_idx']

        self.__worklist = optimise(pd.concat([self.__worklist, loc_df],
                                             axis=1))

    def __get_location(self, src_name, dest_name):
        '''Get location.'''
        srcs = self.__added_comps[src_name]
        dests = self.__added_comps[dest_name]

        shortest_dist = float('inf')
        opt_pair = None

        for src_plt, src_wells in srcs.items():
            for dest_plt, dest_wells in dests.items():
                for src_well in src_wells:
                    for dest_well in dest_wells:
                        src_ind = plate.get_indices(src_well)
                        dest_ind = plate.get_indices(dest_well)
                        dist = cityblock(src_ind, dest_ind)

                        if dist < shortest_dist:
                            src_idx = self.__input_plates[src_plt].get_idx(
                                *src_ind)
                            dest_idx = self.__input_plates[dest_plt].get_idx(
                                *dest_ind)

                            shortest_dist = dist
                            opt_pair = [src_plt, src_well, src_idx,
                                        *plate.get_indices(src_well),
                                        *self.__get_pipette_idx(src_plt,
                                                                src_idx),
                                        dest_plt, dest_well, dest_idx,
                                        *plate.get_indices(dest_well),
                                        *self.__get_pipette_idx(dest_plt,
                                                                dest_idx)]
        return opt_pair

    def __get_pipette_idx(self, src_plt, src_idx):
        '''Get pipetting index (supporting 96 and 384 well plates.'''
        plt = self.__input_plates[src_plt]
        return plt.size(), src_idx % 2 if plt.size() == 384 else 0

    def __traverse(self, dest, level, data):
        '''Traverse tree.'''
        for src in dest.predecessors():
            opr = src[1].copy()

            for key, val in src[0].attributes().items():
                opr['src_' + key] = val

            for key, val in dest.attributes().items():
                opr['dest_' + key] = val

            opr['level'] = level
            opr['src_is_input'] = \
                not src[0].indegree() and not src[0].attributes()['is_reagent']

            data.append(opr)
            self.__traverse(src[0], level + 1, data)

    def __add_component(self, component, plate_id, is_reagent, well_name):
        '''Add component.'''
        if component not in self.__added_comps:
            try:
                (wells, plt) = plate.add_component({'id': component},
                                                   plate_id,
                                                   is_reagent,
                                                   self.__input_plates,
                                                   well_name)

                self.__added_comps[component] = {plt.get_name(): wells}
            except KeyError:
                # Occurs when plate is full:
                self.__add_component(component,
                                     self.__get_next_plate_id(
                                         plate_id),
                                     is_reagent,
                                     well_name)

    def __get_next_plate_id(self, full_plate_id):
        '''Get next plate id.'''
        grps = re.match(r'(.*)~(\d+)', full_plate_id)

        if grps:
            new_plate_id = grps[1] + '~' + str(int(grps[2]) + 1)
        else:
            new_plate_id = full_plate_id + '~2'

        if new_plate_id not in self.__input_plates:
            # Add new plate:
            full_plate = self.__input_plates[full_plate_id]
            rows, cols = full_plate.shape()

            new_plate = plate.Plate(new_plate_id,
                                    rows=rows, cols=cols,
                                    col_ord=full_plate.get_col_order(),
                                    properties=full_plate.get_properties())

            self.__input_plates[new_plate_id] = new_plate

        return new_plate_id


def optimise(df, optimiser=smart_sort_opt):
    '''Optimise.'''
    optimised_dfs = []
    cols = ['level',
            'src_is_reagent',
            # 'src_plate',
            'dest_plate']

    for _, group_df in df.groupby(cols):
        if group_df['src_is_reagent'].all():
            group_df = group_df.copy()
            group_df.loc[:, 'src_name'] = pd.Categorical(group_df['src_name'],
                                                         ['water',
                                                          'buffer', 'ladder',
                                                          'mm', 'mm_dig',
                                                          'mm_lcr',
                                                          'ampligase'])
            for _, subgroup_df in group_df.groupby('src_name'):
                if not subgroup_df.empty:
                    optimised_dfs.append(optimiser.optimise(subgroup_df))
        else:
            optimised_dfs.append(optimiser.optimise(group_df))

    optimised_df = pd.concat(optimised_dfs)

    return optimised_df.sort_values(cols,
                                    ascending=[False, False, True])


def to_csv(wrklst, out_dir_name='.'):
    '''Export worklist as csv file.'''
    path = os.path.abspath(os.path.join(out_dir_name,
                                        str(wrklst.name) + '_worklist.csv'))
    wrklst.to_csv(path, encoding='utf-8', index=False)


def format_worklist(dir_name):
    '''Rename columns to SYNBIOCHEM-specific headers.'''
    dfs = []
    dir_dfs = defaultdict(list)

    for(dirpath, _, filenames) in os.walk(dir_name):
        for filename in filenames:
            if filename.endswith('worklist.csv'):
                filepath = os.path.join(dirpath, filename)
                df = pd.read_csv(filepath)
                _rename_values(df)
                _rename_cols(df)
                df = _reorder_cols(df)
                # df.to_csv(filepath, encoding='utf-8', index=False)
                dfs.append(df)
                dir_dfs[dirpath].append(df)
                os.remove(filepath)

    for dirpath, dfs in dir_dfs.items():
        filepath = os.path.join(dirpath, 'worklist.csv')
        worklist_df = pd.concat(dfs)
        worklist_df = worklist_df[_COLUMNS_ORDER + ['dest_name']]
        worklist_df.to_csv(filepath, index=False)

    return dfs


def _rename_values(df):
    '''Rename values.'''
    for columns, replacement in _VALUES_RENAME.items():
        for to_replace, value in replacement.items():
            df[list(columns)] = df[list(columns)].replace(to_replace, value)


def _rename_cols(df):
    '''Rename cols.'''
    df.rename(columns=_COLUMNS_RENAME, inplace=True)


def _reorder_cols(df):
    '''Reorder cols.'''
    columns = _COLUMNS_ORDER + \
        sorted([col for col in df.columns if col not in _COLUMNS_ORDER])

    return df.reindex(columns, axis=1)
