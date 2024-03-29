#!/bin/env python3
"""
Thorium slug irradiation in SCALE/Origen
Ondrej Chvala <ochvala@utexas.edu>
"""

scale_file: str = 'irr--workday-2y_dec-5y.inp'
volume: float = 12.87  # OD=0.5in, H=4in
th_metal_at_dens: float = 0.0303653  # atom density of Th metal
f71_name: str = 'th-wd.f71'
flux: float = 2e12
weeks: int = 105

header: str = f'''=origen
options{{
    digits=6
}}
bounds {{
    neutron="scale.rev13.xn200g47v7.1"
    gamma="scale.rev13.xn200g47v7.1"
    beta=[100L 1.0e7 1.0e-3]
}}
case(week1) {{
    lib {{
        file="/opt/scale6.3_data/arplibs/irt4m8t_e19.f33" 
        pos=1 
    }}
    mat{{
       iso=[th-232={th_metal_at_dens}]
        units=ATOMS-PER-BARN-CM
        volume={volume}
    }}
    time {{
        units=DAYS
        start=0
        t=[0.166667 0.333334 1  1.166667 1.333334 2  2.166667 2.333334 3  3.166667 3.333334 4  4.166667 4.333334 5 7]
    }}
    flux=[2R {flux} 0 2R {flux} 0 2R {flux} 0 2R {flux} 0 2R {flux} 0 0]
    save {{
        file="{f71_name}"
    }}
}}
'''


def case(week_number: int) -> str:
    out: str = f'''case(week{week_number}) {{
    lib {{
        file="/opt/scale6.3_data/arplibs/irt4m8t_e19.f33" 
        pos=1 
    }}
    time {{
        units=DAYS
        start=0
        t=[0.166667 0.333334 1  1.166667 1.333334 2  2.166667 2.333334 3  3.166667 3.333334 4  4.166667 4.333334 5 7]
    }}
    flux=[2R {flux} 0 2R {flux} 0 2R {flux} 0 2R {flux} 0 2R {flux} 0 0]
    save {{
        file="{f71_name}"
    }}
}}
'''
    return out


footer: str = '''case(decay) {
    gamma=yes
    neutron=yes
    beta=yes
    lib {
        file="end7dec"
    }
    time {
        units=DAYS
        t=[70I 730.5 2556.68]
    }
    save {
        file="th-wd.f71"
    }
}
end
'''

deck_str: str = header
for i in range(weeks-2):
    nw: int = i + 2
    deck_str += case(nw)
deck_str += footer

with open(scale_file, 'w') as f:
    f.write(deck_str)
