#!/bin/env python3
"""
Th-EIRENE with rapid Pa removal analysis script
Ondrej Chvala <ochvala@utexas.edu>
"""
import os
import json5
import numpy as np
import scipy
from sample_decay_dose import SampleDose

f71_file_name: str = 'ThEIRENE.f71'
decay_f71_file_name: str = 'decay-u233.f71'
decay_inp_file_name: str = 'decay-u233.inp'
fuel_volume: float = 1.06806e+07  # [cm^3]
MTiHM: dict = {' 5.00': 10.4159200656709, '19.75': 10.2647529843048}
fuel_type: dict = {' 5.00': 'LEU+Th', '19.75': 'HALEU+Th'}
flux_per_MW: dict = {' 5.00': 1.3070e+11, '19.75': 1.2441e+11}
seconds_per_day: float = 60.0 * 60.0 * 24.0

SampleDose.ATOM_DENS_MINIMUM = 0.0 # 1e-20
# fuel_salt_at_dens_BOC: dict = SampleDose.get_burned_nuclide_atom_dens(f71_file_name, 1)


def decay_Pa_origen_deck(lib_pos: int = None) -> str:
    """ Protactinium decay Origen deck """
    assert lib_pos is not None
    from sample_decay_dose.SampleDose import NOW
    origen_output = f'''
=shell
cp -r ${{INPDIR}}/ThEIRENE.f71 .
end

=origen
' {NOW} 
options{{
    digits=6
}}
bounds {{
    neutron="scale.rev13.xn200g47v7.1"
    gamma="scale.rev13.xn200g47v7.1"
    beta=[100L 1.0e7 1.0e-3]
}}
case {{
    gamma=yes
    neutron=yes
    beta=yes
    lib {{ % decay only library
        file="end7dec"
    }}
    mat {{
      load {{ 
            file="ThEIRENE.f71" pos={lib_pos}
        }}

    }}
    time {{
        units=YEARS
        start=0
        t=[500L 0.001 5]

    }}
    save {{
        file="{decay_f71_file_name}"
    }}
}}
end
'''
    return origen_output


runs: dict = {
    # '20-20MWth/05.00pct_1y': {'power': 20.0, 'enr': ' 5.00'},
    '20-20MWth/05.00pct_1y_gas_removal': {'power': 20.0, 'enr': ' 5.00'},
    # '20-20MWth/19.75pct_1y': {'power': 20.0, 'enr': '19.75'},
    '20-20MWth/19.75pct_1y_gas_removal': {'power': 20.0, 'enr': '19.75'},
    '21-01MWth/05.00pct_1y': {'power': 1.0, 'enr': ' 5.00'},
    '21-01MWth/05.00pct_1y_gas_removal': {'power': 1.0, 'enr': ' 5.00'},
    '21-01MWth/19.75pct_1y': {'power': 1.0, 'enr': '19.75'},
    '21-01MWth/19.75pct_1y_gas_removal': {'power': 1.0, 'enr': '19.75'},
    '22-100MWth/05.00pct_1y': {'power': 100.0, 'enr': ' 5.00'},
    # '22-100MWth/19.75pct_1y': {'power': 100.0, 'enr': '19.75'},
}

os.chdir('/home/o/ThEIRENE/03-Triton-moveiso/')
cwd: str = os.getcwd()
for my_path, my_run in runs.items():
    my_f71_file_name = os.path.join(my_path, f71_file_name)
    thermal_power: float = my_run['power']
    enrichment: str = my_run['enr']
    print(my_f71_file_name)
    df = SampleDose.get_f71_nuclide_case(my_f71_file_name, 'gram', [13])
    dt = df.transpose()
    for col in dt.columns:  # remove ' around nuclide names
        dt.rename(columns={col: col.replace("'", "")}, inplace=True)
    dt.set_index('time', inplace=True)
    dt = dt.multiply(MTiHM[enrichment])
    u233_series = dt['U233'].iloc[(np.where(dt.index > 1e7)[0][0]):]
    slope, intercept, r, p, se = scipy.stats.linregress(u233_series.index, u233_series)
    u233_1y_grams: float = u233_series.iloc[-1]
    runs[my_path]['u233 g per day'] = slope * seconds_per_day
    runs[my_path]['u233 g in 1y'] = u233_1y_grams
    runs[my_path]['u233 g per MWday'] = slope * seconds_per_day / thermal_power
    runs[my_path]['u233 g in 1MWy'] = u233_1y_grams / thermal_power
    print(f'U233: {runs[my_path]["u233 g per day"]:.3e} g/day, '
          f'total year 1: {runs[my_path]["u233 g in 1y"]:.3e} g')
    print(f'U233: {runs[my_path]["u233 g per MWday"]:.3e} g/day, '
          f'total year 1: {runs[my_path]["u233 g in 1MWy"]:.3e} g/MW')

    # ORIGEN decay of Pa-233 box
    last_Pa_pos: int = SampleDose.get_last_position_for_case(my_f71_file_name, 13)
    # print(last_Pa_pos)
    with open(os.path.join(my_path, decay_inp_file_name), 'w') as f:    # write ORIGEN decay deck
        f.write(decay_Pa_origen_deck(last_Pa_pos))
    os.chdir(my_path)
    SampleDose.run_scale(decay_inp_file_name)       # run ORIGEN decay in the correct directory
    df = SampleDose.get_f71_nuclide_case(decay_f71_file_name, 'becq', [1])
    # print(df, os.getcwd(), os.path.exists(decay_f71_file_name), decay_f71_file_name)
    os.chdir(cwd)

    dt = df.transpose()
    for col in dt.columns:  # remove ' around nuclide names
        dt.rename(columns={col: col.replace("'", "")}, inplace=True)
    dt.set_index('time', inplace=True)
    dt = dt.multiply(MTiHM[enrichment])
    ac225_series = dt['Ac225'].iloc[20:]
    slope, intercept, r, p, se = scipy.stats.linregress(ac225_series.index, ac225_series)
    ac225_5y_Bqs: float = ac225_series.iloc[-1]
    runs[my_path]['ac225 Bq per day'] = slope * seconds_per_day
    runs[my_path]['ac225 Bq in 1y'] = ac225_5y_Bqs
    runs[my_path]['ac225 Bq per MWday'] = slope * seconds_per_day / thermal_power
    runs[my_path]['ac225 Bq in 5y/MW'] = ac225_5y_Bqs / thermal_power
    print(f'Ac225: {runs[my_path]["ac225 Bq per day"]:.3e} Bq/day, '
          f'total year 1: {runs[my_path]["ac225 Bq in 1y"]:.3e} Bq')
    print(f'Ac225: {runs[my_path]["ac225 Bq per MWday"]:.3e} Bq/day, '
          f'total year 1: {runs[my_path]["ac225 Bq in 5y/MW"]:.3e} Bq/MW')

with open('u233-results.json', 'w') as f:
    json5.dump(runs, f, indent=4)
