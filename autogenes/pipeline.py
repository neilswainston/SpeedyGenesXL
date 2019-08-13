'''
AutoGenes (c) University of Liverpool 2019

AutoGenes is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author: neilswainston
'''
# pylint: disable=invalid-name
# pylint: disable=wrong-import-order
import os
import shutil

from autogenes import plate, worklist
import pandas as pd


def get_input_plates(dir_name):
    '''Get input plates.'''
    input_plates = {}

    for(dirpath, _, filenames) in os.walk(dir_name):
        for filename in filenames:
            if filename[-4:] == '.csv':
                df = pd.read_csv(os.path.join(dirpath, filename))
                _, name = os.path.split(filename)

                if 'well' in df.columns.values:
                    plt = plate.from_table(df, name)

                else:
                    plt = plate.from_plate(df, name)

                input_plates[plt.get_name()] = plt

    return input_plates


def run(wrtrs, input_plates=None, plate_names=None,
        parent_out_dir_name='.'):
    '''Run pipeline.'''
    if not plate_names:
        plate_names = {}

    parent_out_dir = os.path.abspath(parent_out_dir_name)

    if os.path.exists(parent_out_dir):
        shutil.rmtree(parent_out_dir)

    for idx, writers in enumerate(wrtrs):
        if isinstance(writers, list):
            for wrt_idx, writer in enumerate(writers):
                input_plates.update(_run_writer(writer,
                                                str(idx + 1) + '_' +
                                                str(wrt_idx + 1),
                                                input_plates,
                                                plate_names,
                                                parent_out_dir))
        else:
            input_plates.update(_run_writer(writers,
                                            str(idx + 1),
                                            input_plates,
                                            plate_names,
                                            parent_out_dir))


def _run_writer(writer, name, input_plates, plate_names,
                parent_out_dir):
    '''Run a writer.'''
    out_dir = os.path.join(parent_out_dir, name)
    os.makedirs(out_dir)

    worklist_gen = worklist.WorklistGenerator(writer.get_graph())
    plate_names['output'] = writer.get_output_name()
    wrklsts, plates = worklist_gen.get_worklist(input_plates, plate_names)

    for plt in plates.values():
        plt.to_csv(os.path.join(out_dir, 'plates'))

    for wrklst in wrklsts:
        worklist.to_csv(wrklst, out_dir)

    summary_df = _summarise(wrklsts)
    summary_df.to_csv(os.path.join(out_dir, 'input_summary.csv'), index=False)

    return plates


def _summarise(worklists):
    '''Summarise worklists.'''
    vols = {}
    dest_plates = set()

    for wrklst in worklists:
        # Add destination plates:
        dest_plates.update(wrklst['dest_plate'].values)

        # Add volumes:
        for keys, group_df in wrklst.groupby(['src_plate',
                                              'src_well',
                                              'src_name']):

            if keys not in vols:
                vols[keys] = 0

            vols[keys] += group_df['Volume'].sum()

    out = [list(key) + [value, None] for key, value in vols.items()]

    for dest_plate in dest_plates:
        out.append([None, None, None, None, dest_plate])

    return pd.DataFrame(out,
                        columns=['src_plate', 'src_well', 'src_name',
                                 'total_volume', 'dest_plate'])
