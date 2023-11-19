#!/bin/env python3
"""
Opus plot reader. Expects one spectrum in the file.
Ondrej Chvala <ochvala@utexas.edu>
"""
import re


def integrate_opus(plt_file_name: str) -> float:
    """
    Integrates the first spectrum in an OPUS .plt file
    """
    integral: float = 0
    match_positive_number = re.compile(' +[0-9]+.?[0-9]*(?:[Ee][+-]? *[0-9]+)?')
    debug: int = 0

    with open(plt_file_name, 'r') as f:
        i_line: int = 0
        prev_x: float = -1.0
        prev_y: float = -1.0
        # is_bin_read: bool = False   # flag if the full bin was read in
        for line in f.read().splitlines():
            i_line += 1
            if i_line <= 6:
                continue
            m = re.findall(match_positive_number, line)
            if len(m) != 2:  # expecting only two numbers
                break
            x: float = float(m[0])
            y: float = float(m[1])
            if y != prev_y:  # read the first line of bin
                if debug > 1:
                    print(i_line, line)
                prev_x = x
                prev_y = y
            else:
                if debug > 1:
                    print(i_line, line, x, prev_x, y, prev_y, integral)
                if prev_x >= x:
                    raise ValueError(f"Bins mis-formatted at line {i_line}")
                integral += y * (x - prev_x)

    return integral
