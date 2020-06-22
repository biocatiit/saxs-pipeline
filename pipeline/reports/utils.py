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
