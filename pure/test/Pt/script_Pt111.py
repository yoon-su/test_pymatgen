# Import statements 
from pymatgen import Structure, Lattice, MPRester, Molecule
from pymatgen.analysis.adsorption import * 
from pymatgen.core.surface import generate_all_slabs 
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer 
from matplotlib import pyplot as plt 
#% matplotlib inline 
# Note that you must provide your own API Key, which can 
# be accessed via the Dashboard at materialsproject.org 
mpr = MPRester()
fcc_pt = Structure.from_spacegroup("Fm-3m", Lattice.cubic(3.92), ["Pt", "Pt"], [[0, 0, 0], [0.5, 0.5, 0.5]])
slabs = generate_all_slabs(fcc_pt, max_index=1, min_slab_size=8.0, min_vacuum_size=10.0)
pt_111 = [slab for slab in slabs if slab.miller_index==(1,1,1)][0]
asf_pt_111 = AdsorbateSiteFinder(pt_111) 
ads_sites = asf_pt_111.find_adsorption_sites() 
#print(ads_sites) 
assert len(ads_sites) == 4

fig = plt.figure() 
ax = fig.add_subplot(111) 
plot_slab(pt_111, ax, adsorption_sites=True)


fig = plt.figure() 
ax = fig.add_subplot(111) 
adsorbate = Molecule("H", [[0, 0, 0]]) 
ads_structs = asf_pt_111.generate_adsorption_structures(adsorbate, repeat=[1, 1, 1])
plot_slab(ads_structs[1], ax, adsorption_sites=False, decay=0.09)

fig = plt.figure() 
ax = fig.add_subplot(111) 
plot_slab(ads_structs[2], ax, adsorption_sites=False, decay=0.09)
fig = plt.figure() 
ax = fig.add_subplot(111) 
plot_slab(ads_structs[3], ax, adsorption_sites=False, decay=0.09)
plt.show()
