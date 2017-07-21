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
fcc_ni = Structure.from_spacegroup("Fm-3m", Lattice.cubic(3.5), ["Ni", "Ni"], [[0, 0, 0], [0.5, 0.5, 0.5]])
slabs = generate_all_slabs(fcc_ni, max_index=1, min_slab_size=8.0, min_vacuum_size=10.0)
ni_100 = [slab for slab in slabs if slab.miller_index==(1,0,0)][0]
asf_ni_100 = AdsorbateSiteFinder(ni_100) 
ads_sites = asf_ni_100.find_adsorption_sites() 
#print(ads_sites) 
assert len(ads_sites) == 4
fig = plt.figure() 
ax = fig.add_subplot(111) 
plot_slab(ni_100, ax, adsorption_sites=True)

fig = plt.figure() 
ax = fig.add_subplot(111) 
adsorbate = Molecule("H", [[0, 0, 0]]) 
ads_structs = asf_ni_100.generate_adsorption_structures(adsorbate, repeat=[1, 1, 1])
plot_slab(ads_structs[0], ax, adsorption_sites=False, decay=0.09)
plt.show()
