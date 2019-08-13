'''
AutoGenes (c) University of Liverpool 2019

AutoGenes is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author: neilswainston
'''
# pylint: disable=invalid-name
# pylint: disable=too-many-arguments
import itertools
import math
import os

import pandas as pd


class Plate():
    '''Class to represent a well plate.'''

    def __init__(self, name, rows=8, cols=12, col_ord=False, properties=None,
                 plate=None):

        if not properties:
            properties = ['id']

        assert 'id' in properties

        perms = list(itertools.product(properties, list(range(1, cols + 1))))
        columns = pd.MultiIndex.from_tuples(perms)

        if plate is not None:
            self.__plate = plate
        else:
            self.__plate = pd.DataFrame(index=[chr(r + ord('A'))
                                               for r in range(0, rows)],
                                        columns=columns)
        self.__plate.name = name
        self.__col_ord = col_ord
        self.__next = 0

    def get_name(self):
        '''Get name.'''
        return self.__plate.name

    def get_col_order(self):
        '''Get column order.'''
        return self.__col_ord

    def get_properties(self):
        '''Get properties.'''
        return list(self._Plate__plate.columns.levels[0])

    def shape(self):
        '''Get plate shape.'''
        return self.__plate['id'].shape

    def size(self):
        '''Get plate size.'''
        rows, cols = self.shape()
        return rows * cols

    def set(self, obj, row, col):
        '''Set object at a given row, col.'''
        self.__next = max(self.__next, self.get_idx(row, col) + 1)

        for key, val in obj.items():
            self.__plate.loc[:, (key, col + 1)][row] = val

    def get(self, row, col):
        '''Get object at a given row, col.'''
        keys = self.get_properties()

        return {key: self._Plate__plate.loc[:, (key, col + 1)][row]
                for key in keys
                if _is_value(self._Plate__plate.loc[:, (key, col + 1)][row])}

    def get_all(self):
        '''Get all objects.'''
        rows, cols = self.shape()
        return {get_well_name(row, col): self.get(row, col)
                for row in range(rows) for col in range(cols)
                if self.get(row, col)}

    def get_by_well(self, well_name):
        '''Get by well, e.g. by C12.'''
        row, col = get_indices(well_name)
        return self.get(row, col)

    def add(self, obj, well_name=None):
        '''Adds an object to the next well.'''
        if well_name:
            row, col = get_indices(well_name)

            for key, val in obj.items():
                self.__plate.loc[:, (key, col + 1)][row] = val
            return well_name

        # else:
        return self.__set(obj, self.__next)

    def add_line(self, obj):
        '''Adds a line of objects (row or col) in next empty line.'''
        if self.__col_ord:
            line_len = len(self.__plate.columns)
        else:
            line_len = len(self.__plate.index)

        start = ((self.__next + line_len - 1) // line_len) * line_len

        for idx in range(start, start + line_len):
            row, col = self.get_row_col(idx)
            self.set(obj, row, col)

    def find(self, src_terms):
        '''Finds an object.'''
        wells = [[str(row) + str(col)
                  for col in self.__plate[key].columns
                  for row in self.__plate[key].index[
            self.__plate[key][col] == val]]
            for key, val in src_terms.items()]

        if len(wells) > 1:
            wells = list(set(wells[0]).intersection(*wells))
        else:
            wells = wells[0]

        return wells

    def get_row_col(self, idx):
        '''Map idx to well.'''
        return get_row_col(idx, self.__plate.shape, self.__col_ord)

    def get_idx(self, row, col):
        '''Map row, col to idx.'''
        return get_idx(row, col, self.__plate.shape, self.__col_ord)

    def to_csv(self, out_dir_name='.'):
        '''Export plate to csv.'''
        if not os.path.exists(out_dir_name):
            os.makedirs(out_dir_name)

        filepath = os.path.abspath(os.path.join(out_dir_name,
                                                str(self.__plate.name) +
                                                '.csv'))
        self.__plate.to_csv(filepath, encoding='utf-8')

    def __set(self, obj, idx):
        '''Sets an object in the given well.'''
        row, col = self.get_row_col(idx)
        self.set(obj, row, col)
        return get_well_name(row, col)

    def __repr__(self):
        return self.__plate.__repr__()


def get_row_col(idx, shape=(8, 12), col_ord=False):
    '''Map idx to well.'''
    rows, cols = shape

    if col_ord:
        return int(idx / cols), int(idx % cols)

    return int(idx % rows), int(idx / rows)


def get_idx(row, col, shape=(8, 12), col_ord=False):
    '''Map idx to well, column ordered.'''
    rows, cols = shape

    if col_ord:
        return row * cols + col

    return col * rows + row


def get_indices(well_name):
    '''Get indices from well name.'''
    return ord(well_name[0]) - ord('A'), int(well_name[1:]) - 1


def get_well_name(row, col):
    '''Get well name from indices.'''
    return str(chr(row + ord('A'))) + str(col + 1)


def find(plates, obj):
    '''Find object in plates.'''
    found = {}

    for plt in plates.values():
        wells = plt.find(obj)

        if wells:
            found[plt.get_name()] = wells

    return found


def add_component(component, plate_id, is_reagent, plates, well_name):
    '''Add a component to a plate.'''
    for plate in plates.values():
        wells = plate.find(component)

        if wells:
            return wells, plate

    if plate_id not in plates:
        plate = Plate(plate_id)
        plates[plate_id] = plate
    else:
        plate = plates[plate_id]

    if is_reagent:
        if well_name:
            return [plate.add(component, well_name)], plate
        # else:
        plate.add_line(component)
        return add_component(component, plate_id, is_reagent, plates,
                             well_name)

    return [plate.add(component, well_name)], plate


def from_table(df, name):
    '''Generate Plate from tabular data.'''
    # df = pd.read_csv(filename, dtype={'id': object, 'parent': object})
    df['id'] = df['id'].astype(str)

    if 'parent' in df.columns.values:
        df['parent'] = df['parent'].astype(str)

    # 96 or 384?
    if len(df) > 96:
        rows, cols = 16, 24
    else:
        rows, cols = 8, 12

    props = list(df.columns[df.columns != 'well'])

    plt = Plate(name.split('.')[0], rows=rows, cols=cols, properties=props)

    for _, row in df.iterrows():
        dct = row.to_dict()
        well = dct.pop('well')
        plt.add(dct, well)

    return plt


def from_plate(df, name):
    '''Generate Plate from tabular data.'''
    df = df.astype(object)
    df.set_index('Unnamed: 0', inplace=True)
    df.columns = pd.MultiIndex.from_product([['id'], df.columns.astype(int)])
    return Plate(name.split('.')[0], plate=df)


def _is_value(val):
    '''Return boolean depending on whether value is None or NaN.'''
    return bool(val and not (isinstance(val, float) and math.isnan(val)))
