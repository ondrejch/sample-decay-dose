import os
import re
import numpy as np
import json5
import matplotlib.pyplot as plt
# from leadcell import sample_mass
sample_mass: float = 0.1                # 0.1 g sample

cwd: str = os.getcwd()
particles = {'2': 'Gamma, contact dose', '6': 'Gamma, 30cm handling dose'}
data = {'2': 'slategrey', '6': 'crimson'}

LABEL = 'fs'
labels = {'fs': ['decay time', 'days', f'Burned fuel salt sample {sample_mass:.1f} g'], }

dose = {}  # doses [mrem/h]
errd = {}  # stdev of doses
r = {}

for d in data.keys():
    with open(os.path.join(cwd, 'responses.json')) as fin:
        r[d] = json5.load(fin)
        dose[d] = np.array([v[d]['value'] for _, v in r[d].items()], float) * 1e3  # mrem/h !!!
        errd[d] = np.array([v[d]['stdev'] for _, v in r[d].items()], float) * 1e3  # mrem/h !!!

xlist = list(r[list(data.keys())[0]].keys())
x = np.array(xlist, float)  # x coordinate - sample masses
closest_to_1month: str = min(xlist, key=lambda t: abs(float(t) - 30.0))
idx_1month: int = list(xlist).index(closest_to_1month)

# Plots!
plt.close('all')
plt.xscale('linear')
plt.yscale('linear')
plt.grid()
plt.title(labels[LABEL][2])
plt.xlabel(f'Sample {labels[LABEL][0]} [{labels[LABEL][1]}]')
plt.ylabel('Dose [mrem/h]')

for d in data.keys():
    my_title = f'{particles[d]}, at {float(closest_to_1month):.1f} days = {dose[d][idx_1month]:.3f} mrem/h'
    print(my_title)
    plt.errorbar(x, dose[d], errd[d], ls='none', color=f'{data[d]}', capsize=0.8)
    plt.scatter(x, dose[d], color=f'{data[d]}', s=5, label=my_title)

plt.legend()
plt.tight_layout()
label_file_name = labels[LABEL][0].replace(' ', '_')
plt.savefig(f'dose_burned-salt_{label_file_name}.png', dpi=1000)

plt.xscale('log')
plt.yscale('log')
plt.savefig(f'dose_burned-salt_{label_file_name}-loglog.png', dpi=1000)
