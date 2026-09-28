#!/bin/env python3
"""
Opus plot reader. Expects one spectrum in the file.
Ondrej Chvala <ochvala@utexas.edu>
"""


def integrate_opus(plt_file_name: str) -> float:
    """
    Integrates the first spectrum in an OPUS .plt file.
    After 6 header lines, OPUS writes each histogram bin as two points, (E_low, y) and (E_high, y).
    The integral is sum(y * (E_high - E_low)), e.g. particles/s for units=intensity [1/(s MeV)] vs E [MeV].
    Adjacent bins with equal y are separate bins. Reading stops at the first line that is not a pair of numbers,
    such as the time label of the next spectrum. An unterminated last point adds nothing.
    """
    integral: float = 0

    with open(plt_file_name, 'r') as f:
        i_line: int = 0
        bin_low: tuple[float, float] | None = None  # (x, y) of the open bin, None between bins
        for line in f.read().splitlines():
            i_line += 1
            if i_line <= 6:
                continue
            tokens: list[str] = line.split()
            if len(tokens) != 2:  # expecting only two numbers
                break
            try:
                x: float = float(tokens[0])
                y: float = float(tokens[1])
            except ValueError:
                break
            if bin_low is None:  # the first point of a bin
                bin_low = (x, y)
                continue
            x_low, y_low = bin_low
            if y != y_low or x < x_low:
                raise ValueError(f"Bins mis-formatted at line {i_line}")
            integral += y * (x - x_low)
            bin_low = None

    return integral
