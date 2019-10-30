# mz_explorer - .mzXML and .mzML m/z data extractor
# Written by Harlan R. Williams
#
# References:
# https://www.ibm.com/developerworks/xml/library/x-hiperfparse/

import argparse
from base64 import b64decode
import csv
from os import path
import numpy as np
import re
import sys
import struct
import zlib

from lxml import etree
# import pandas


# Opening x.mzML
# Scanning...
# Total run time: 24.99 minutes
# Options:
# e - Export selection to a csv file
# h - Display options
# r - Select an m/z range to limit scans to
# t - Select a time range to pull scans from
#
# t
# Enter the range in minutes (none for whole file) > 5 - 6
# Selected: 34 scans
#
# r
# Enter the range to limit m/z values to > 1000 - 2000
# Done
#
# e
# Enter a filename (none for default) >
# Wrote x.csv


parser = argparse.ArgumentParser(description='Extract LC/MS spectral data from .mzML and .mzXML files')
parser.add_argument('filepath', type=str, help='Filepath of an .mzML or .mzXML file')


def main():
    args = parser.parse_args()
    filepath = args.filepath

    if not path.isfile(filepath):
        print('Error: File does not exist')
        sys.exit(1)

    print('Opening {}'.format(filepath))
    filetype = filepath.split('.')[-1]
    name = filepath.split('.')[0].split('/')[-1]

    time_start, time_stop = 0, 0

    # create a list of scan IDs and tell the user its length
    if filetype == 'mzXML':
        lns = '{http://sashimi.sourceforge.net/schema_revision/mzXML_3.2}'
        number_scans = 0
        context = etree.iterparse(filepath, events=('end',), tag=(lns + 'msRun'))
        for event, elem in context:
            if elem.tag == lns + 'msRun':
                number_scans = int(elem.attrib.get('scanCount'))
                time_start = float(elem.attrib.get('startTime').split('T')[-1].split('S')[0])
                time_stop = float(elem.attrib.get('endTime').split('T')[-1].split('S')[0])
                total_time = time_stop - time_start
                print('Found {} scans'.format(number_scans))
                print('Total run time: {}'.format(total_time), 'seconds')
                break
            else:
            	break
            while elem.getprevious() != None:
                del elem.getparent()[0]
    else:
        print('Error: Filetype is not of type .mzXML nor .mzML')
        sys.exit(1)

    repeat = True
    count = 0
    while repeat is True:
        # prompt the user for an index or range of indices to extract spectrum data from
        indices = None
        while indices is None:
            indices = get_range_from_user("Enter a time range to pull data from, " \
                                          f"or 'q' to exit ({time_start} - {time_stop}) > ",
                                          low=time_start, high=time_stop)
            if indices is None:  # error parsing input
                print('Error: Out of range')
            elif indices is False:  # user entered 'q' to quit
                sys.exit(1)
        lim_mz_range = False
        mz_range = []
        if input('Limit m/z to a range? (y/n)> ') == 'y':
            lim_mz_range = True
            mz_range = get_range_from_user('Enter the range > ')

        write_every_row = False
        if input('Write m/z for each row? (y/n)> ') == 'y':
            write_every_row = True
        root = []
        tolerance = 0.001  # floats +/- tolerance are considered equal
        has_set_root = False
        first_row_written = False
        with open('{}_{}.csv'.format(name, count), 'w') as csvfile:
            writer = csv.writer(csvfile)
            if filetype == 'mzXML':
                lns = '{http://sashimi.sourceforge.net/schema_revision/mzXML_3.2}'
                context = etree.iterparse(filepath, events=('end',), tag=(lns + 'scan'))
                for event, elem in context:
                    retention_time = float(elem.attrib.get('retentionTime').split('T')[-1].split('S')[0])
                    if retention_time >= indices[0] and retention_time <= indices[1]:
                        peaks = elem.getchildren()[0]
                        mz_array, int_array = grab_data_mzxml(peaks)
                        if lim_mz_range and mz_range is not None:
                            mz_indices = find_mz_range_indices(mz_array, mz_range, tolerance)
                            mz_array = mz_array[int(mz_indices[0]):int(mz_indices[1])]
                            int_array = int_array[int(mz_indices[0]):int(mz_indices[1])]
                        expected_delta_mz = find_expected_delta_mz(mz_array)
                        print(expected_delta_mz)
                        adj_mz_array, adj_int_array = [], []

                        first_mz = mz_array[0]
                        lower_mz_offset = round((first_mz - mz_range[0]) / expected_delta_mz)
                        for i in range(lower_mz_offset):
                            adj_mz_array.append(mz_range[0] + expected_delta_mz * i)
                            adj_int_array.append(0)
                        # for i in range(1, len(mz_array)):
                        #    delta_mz = mz_array[i] - mz_array[i - 1]
                        #    if abs(delta_mz - expected_delta_mz) > tolerance:
                        #        mz_offset = round(delta_mz / expected_delta_mz)
                        #        for j in range(mz_offset):
                        #            adj_mz_array.append(mz_array[i] + expected_delta_mz * j)
                        #            adj_int_array.append(0)
                        #    adj_mz_array.append(mz_array[i])
                        #    adj_int_array.append(int_array[i])
                        for i in range(1, len(mz_array) - 1):
                            adj_mz_array.append(mz_array[i])
                            adj_int_array.append(int_array[i])
                            if int_array[i] == 0:
                                delta_mz_next = mz_array[i + 1] - mz_array[i]
                                delta_mz_prev = mz_array[i] - mz_array[i - 1]
                                if abs(delta_mz_next - delta_mz_prev) > tolerance:
                                    mz_offset = round(delta_mz_next / delta_mz_prev)
                                    for j in range(1, mz_offset):
                                        adj_mz_array.append(mz_array[i] + delta_mz_prev * j)
                                        adj_int_array.append(0)

                        if not first_row_written or write_every_row:
                            first_row_written = True
                            writer.writerow(adj_mz_array)
                        writer.writerow(adj_int_array)
                    elem.clear()
                    if retention_time > indices[1]:
                        break
                    while elem.getprevious() is not None:
                        del elem.getparent()[0]
        repeat = input('Obtain data from another range? (y/n)> ') == 'y'


def get_range_from_user(prompt, low=-1, high=-1):
    indices = input(prompt)
    if indices == 'q':
        return False

    indices = [float(x) for x in indices.split('-')]

    print(indices)

    if len(indices) == 2:
        indices.sort()
    elif len(indices) == 1:
        indices += indices
    else:
        return None

    indices = tuple(indices)
    # print(indices)
    # if low > 0 and high > 0:  # if range checking
    #     if min(indices) <= low or max(indices) > high:
    #         return None
    return indices


def align_mz_int_arrays_to_root(root, mz_array, int_array, tolerance):
    pivot_tolerance = 50
    middle_pivot = len(root) // 2
    possible_pivots = range(middle_pivot - pivot_tolerance, middle_pivot + pivot_tolerance)

    found_pivot = False
    for pivot_ind in possible_pivots:  # Try a range of possible pivots incase a pivot is missing
        shift_by = [i for i, number in enumerate(mz_array) if abs(number - root[pivot_ind]) < tolerance]
        if len(shift_by) == 1:
            shift_by = shift_by[0] - pivot_ind
            found_pivot = True
            break
    if not found_pivot:
    	raise ValueError('No suitable pivot found.')

    print(shift_by)

    if shift_by > 0:
        mz_array = mz_array[shift_by:] + (shift_by + 1) * [float('nan')]
        int_array = int_array[shift_by:] + (shift_by + 1) * [float('nan')]
    elif shift_by < 0:
        mz_array = abs(shift_by) * [float('nan')] + mz_array[abs(shift_by):]
        int_array = abs(shift_by) * [float('nan')] + int_array[abs(shift_by):]
    return mz_array, int_array


def grab_data_mzxml(peaks):
    float_size_bytes = 0
    compressed = False

    if peaks.attrib.get('compressionType') == 'zlib':
        compressed = True

    if peaks.attrib.get('precision') == '32':
        fmt = 'f'
        float_size_bytes = 4
    elif peaks.attrib.get('precision') == '64':
        fmt = 'd'
        float_size_bytes = 8
    else:
        return None  # couldn't determine float size

    data = b64decode(bytes(peaks.text, encoding='ascii'))
    if compressed:
        data = zlib.decompress(data)
    data = list(struct.unpack('!{}{}'.format((str(len(data) // float_size_bytes)), fmt), data))
    mz_array = data[0::2]
    int_array = data[1::2]  # data in the .mzXML format is stored as alternating m/z ratio and intensity data

    # need to add check for byte order
    return mz_array, int_array

# find the array indices which correspond to a certain range of m/z values
def find_mz_range_indices(mz_array, mz_range, tolerance=0.001):
    assert mz_range[0] < mz_range[1]
    mz_indices = []
    found_lower = False
    found_upper = False
    for ind, val in zip(range(len(mz_array)), mz_array):
        if val > mz_range[0] and not found_lower:
            mz_indices.append(ind)
            found_lower = True
        if val > mz_range[1] and not found_upper:
            mz_indices.append(ind)
            found_upper = True
    return mz_indices

# takes the first 100 values of an m/z array and finds an "expected" delta m/z
from math import log10

def find_expected_delta_mz(mz_array, max_diff=1.25):
    assert max_diff < 2

    if len(mz_array) < 100:
        value_range = len(mz_array)
    else:
        value_range = 100
    delta_array = []
    for i in range(1, value_range):
        delta_array.append(mz_array[i] - mz_array[i - 1])

    min_delta = min(delta_array)
    acc = 0
    vals = 0

    for i in delta_array:
        if i < min_delta * max_diff:
            acc += i
            vals += 1

    avg = acc / vals
    return avg

# def grab_data_mzml(binary_data_list):
#     mz_array, int_array = [], []
#     float_size_bytes = 0
#     compressed = False
#     binary_data_elems = binary_data_list.getchildren()

#     for elem in binary_data_elems:
#         if elem[0].attrib.get('name') == '32-bit float':
#             float_size_bytes = 4
#         elif elem[0].attrib.get('name') == '64-bit float':
#             float_size_bytes = 8
#         else:
#             return None  # couldn't determine float size

#         if elem[1].attrib.get('name') == 'zlib compression':
#             compressed = True

#         data = b64decode(bytes(elem[3].text, encoding='ascii'))
#         if compressed:
#             data = zlib.decompress(data)
#         floats = list(struct.unpack('f' * (len(data) // float_size_bytes), data))
#         if elem[2].attrib.get('name') == 'm/z array':
#             mz_array = floats
#         elif elem[2].attrib.get('name') == 'intensity array':
#             int_array = floats

#     return mz_array, int_array


if __name__ == '__main__':
    main()
