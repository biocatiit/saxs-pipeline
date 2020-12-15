# coding: utf-8
#
#    Project: BioCAT SAXS pipeline
#             https://github.com/biocatiit/saxs-pipeline
#
#
#    Principal author:       Jesse Hopkins
#
#    This is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This software is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this software.  If not, see <http://www.gnu.org/licenses/>.

import numpy as np

def text_round(value, round_to):
    value = float(value)

    low_bound = 1./(10.**(round_to))

    if round_to > 1:
        high_bound = 1000./(10.**(round_to-1))
    else:
        high_bound = 1000

    if value < low_bound or value > high_bound:
            value = str(np.format_float_scientific(value, round_to, trim='0',
            exp_digits=1))
    else:
        value = str(round(value, round_to))

    return value

def rotate_list(my_list):
    rot_list = []

    dim1 = len(my_list)
    dim2 = len(my_list[0])

    for j in range(dim2):
        data = []

        for k in range(dim1):
            data.append(my_list[k][j])

        rot_list.append(data)

    return rot_list
