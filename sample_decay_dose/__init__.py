import json
import os


def _load_isotopic_data():
    """
    Loads isotopic data from the JSON file and makes it available
    at the package level.
    """
    # The path is constructed relative to this __init__.py file.
    # It navigates up to the repository root, then into the 'data' directory.
    data_path = os.path.join(os.path.dirname(__file__), 'data/isotopic_data.json')
    try:
        with open(data_path, 'r', encoding='utf-8') as f:
            data = json.load(f)
            # JSON keys are strings by default; this converts the mass numbers back to integers.
            return {elem: {int(mass_num): props for mass_num, props in isotopes.items()} for elem, isotopes in
                data.items()}
    except (FileNotFoundError, json.JSONDecodeError) as e:
        # Provide a clear error message if the data file is missing or corrupt.
        raise RuntimeError(f"Could not load isotopic data from {data_path}. "
                           f"Please ensure 'download_NIST_nuclide_data.py' has been run successfully. Error: {e}")


# This line executes when the 'sample_decay_dose' package is imported,
# loading the data once and making it available globally within the package.
ISOTOPIC_DATA = _load_isotopic_data()
