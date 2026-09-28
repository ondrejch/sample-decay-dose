"""Regenerate ``data/isotopic_data.json`` from the NIST "Atomic Weights and Isotopic Compositions" listing.

Usage::

    python -m sample_decay_dose.download_NIST_nuclide_data            # parse the shipped data/aw.html
    python -m sample_decay_dose.download_NIST_nuclide_data --source URL_OR_FILE --output FILE

Input and output default to the package's ``data`` directory, resolved from this file, so the result is the
same from any working directory.
"""
import argparse
import json
import re
from pathlib import Path

DATA_DIR = Path(__file__).resolve().parent / 'data'
DEFAULT_SOURCE = DATA_DIR / 'aw.html'
DEFAULT_OUTPUT = DATA_DIR / 'isotopic_data.json'
NIST_URL = "https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&ascii=ascii2&isotype=all"

# Deuterium and tritium are listed as separate symbols by NIST and belong to hydrogen. The NIST listing still
# uses the temporary IUPAC names Uup (Z=115) and Uus (Z=117); their official names are Mc and Ts.
SYMBOL_MAP = {'d': 'h', 't': 'h', 'uup': 'mc', 'uus': 'ts'}


def _read_source(source) -> str:
    """Return the text of a local file, or of an http(s) URL."""
    source = str(source)
    if source.startswith(('http://', 'https://')):
        import requests  # Imported here so that parsing a local file does not need it.
        response = requests.get(source, timeout=60)
        response.raise_for_status()
        return response.text
    if source.startswith('file://'):
        from urllib.parse import urlparse
        from urllib.request import url2pathname
        source = url2pathname(urlparse(source).path)
    return Path(source).read_text(encoding='utf-8')


def parse_nist_text(text: str) -> dict:
    """Parse the NIST 'Key = Value' records into {symbol: {mass_number: {'mass': ..., 'abundance': ...}}}."""
    isotopic_data = {}
    current_record = {}

    # Helper to process and save a completed record
    def process_record(record):
        if not record or "Atomic Symbol" not in record or "Mass Number" not in record:
            return

        symbol = record["Atomic Symbol"].lower()
        symbol = SYMBOL_MAP.get(symbol, symbol)

        if symbol not in isotopic_data:
            isotopic_data[symbol] = {}

        mass_number = int(record["Mass Number"])

        # Default abundance to 0 if not present
        abundance_str = record.get("Isotopic Composition", "0")
        mass_str = record.get("Relative Atomic Mass", "0")

        # Remove uncertainties from mass, e.g., "1.0078250322(1)" -> "1.0078250322"
        mass = float(re.sub(r"\(.*\)", "", mass_str))

        # Abundance can be a range, e.g., "[0.999816,0.999974]". We take the average.
        if '[' in abundance_str:
            try:
                low, high = map(float, abundance_str.strip('[]').split(','))
                abundance = (low + high) / 2.0
            except ValueError:
                abundance = 0.0
        elif not abundance_str.strip():  # Handle empty abundance
            abundance = 0.0
        else:
            # Also remove uncertainties from abundance
            abundance = float(re.sub(r"\(.*\)", "", abundance_str))

        isotopic_data[symbol][mass_number] = {"mass": mass, "abundance": abundance}

    # Regex to capture "Key = Value" lines
    record_regex = re.compile(r"^\s*([^=]+?)\s*=\s*(.*)")

    for line in text.splitlines():
        # A new record starts with "Atomic Number"
        if line.strip().startswith("Atomic Number"):
            # Process the previous record before starting a new one
            process_record(current_record)
            current_record = {}

        match = record_regex.match(line)
        if match:
            key, value = match.groups()
            current_record[key.strip()] = value.strip()

    # Process the last record in the file
    process_record(current_record)
    return isotopic_data


def download_and_parse_nist_data(source=DEFAULT_SOURCE, output_path=DEFAULT_OUTPUT) -> dict:
    """
    Reads the NIST isotopic data (by default the shipped copy in data/aw.html; NIST_URL gives the live listing),
    parses it, and saves it as JSON. The default output is the package's data/isotopic_data.json, which
    sample_decay_dose.data loads at import.
    """
    print(f"Reading data from: {source}")
    text = _read_source(source)

    print("Parsing data...")
    isotopic_data = parse_nist_text(text)

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    print(f"Saving parsed data to {output_path}...")
    with open(output_path, 'w', encoding='utf-8') as f:
        json.dump(isotopic_data, f, indent=4)

    print("Done.")
    return isotopic_data


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument('--source', default=str(DEFAULT_SOURCE),
                        help=f"Local file or http(s) URL of the NIST listing (default: {DEFAULT_SOURCE}; "
                             f"live listing: {NIST_URL}).")
    parser.add_argument('--output', default=str(DEFAULT_OUTPUT), help=f"Output JSON (default: {DEFAULT_OUTPUT}).")
    args = parser.parse_args(argv)
    download_and_parse_nist_data(args.source, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
