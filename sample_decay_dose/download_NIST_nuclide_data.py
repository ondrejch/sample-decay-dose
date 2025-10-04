import requests
from requests_file import FileAdapter

import re
import json
from pathlib import Path


def download_and_parse_nist_data():
    """
    Downloads isotopic data from NIST (or uses a local file),
    parses it, and saves it as a JSON file.
    """
    # The original NIST URL seems to be unstable.
    # The code is adapted to use a local file provided by the user.
    # url = "https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&ascii=ascii2&isotype=all"

    # Using the local file as specified.
    # Ensure 'aw_shortened.html' is in the correct path relative to the script.
    # The path provided by the user was absolute, making it relative for portability.
    local_file_path = Path(__file__).parent.parent / 'sample_decay_dose/data/aw.html'
    url = local_file_path.as_uri()

    print(f"Reading data from local file: {url}")
    s = requests.Session()
    s.mount('file://', FileAdapter())

    try:
        response = s.get(url)
        response.raise_for_status()
        print("File read successfully.")
    except requests.exceptions.RequestException as e:
        print(f"Error reading local file: {e}")
        return

    lines = response.text.splitlines()

    isotopic_data = {}
    current_record = {}

    # Helper to process and save a completed record
    def process_record(record):
        if not record or "Atomic Symbol" not in record or "Mass Number" not in record:
            return

        symbol_map = {'d': 'h', 't': 'h'}
        symbol = record["Atomic Symbol"].lower()
        symbol = symbol_map.get(symbol, symbol)  # Map D and T to H

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

    print("Parsing data...")
    # Regex to capture "Key = Value" lines
    record_regex = re.compile(r"^\s*([^=]+?)\s*=\s*(.*)")

    for line in lines:
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

    # Create the data directory if it doesn't exist
    data_dir = Path("data")
    data_dir.mkdir(exist_ok=True)

    output_path = data_dir / "isotopic_data.json"

    print(f"Saving parsed data to {output_path}...")
    with open(output_path, 'w', encoding='utf-8') as f:
        json.dump(isotopic_data, f, indent=4)

    print("Done.")


if __name__ == "__main__":
    download_and_parse_nist_data()
