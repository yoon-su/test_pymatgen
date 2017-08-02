#!/usr/bin/env python
import time, sys, os
import numpy as np

version = 20170802
nt = time.localtime()
now_time = "%s_%s_%s_%s:%s:%s" % (nt[0],nt[1],nt[2],nt[3],nt[4],nt[5])
usage    = ' Usage: %s [element1],[element2], lattice constant,  h,k,l \n ' % sys.argv[0]
foottext = '\n Thank you\n## Yoon Su Shim (KAIST, Graduate School of EEWS) <yoonsushim@kaist.ac.kr>'
print("## Creating images of absorbate on slab (fcc) using 'pymatgen'")
print("## Version : %s \n" % version)
print("## Updated: combine with ASE for molecule database and A2P2")

print(now_time)

# Check input file
if len(sys.argv) != 7:
	if len(sys.argv) == 2:
		filename = str(sys.argv[1])
		print('Reading from structure file ...')
	else:
	        print(usage)
	        print(foottext)
	        sys.exit(1)

if len(sys.argv) == 7:
	element1  = sys.argv[1]
	element2  = sys.argv[2]
	a0       = float(sys.argv[3])
	MI       = (int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6]))

# Import statements 
from pymatgen import Structure, Lattice, MPRester, Molecule
from pymatgen.analysis.adsorption import * 
from pymatgen.core.surface import generate_all_slabs 
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer 
from matplotlib import pyplot as plt 
from pymatgen.io.vasp import Poscar
from pymatgen.core import IStructure

#% matplotlib inline 
# Note that you must provide your own API Key, which can 
# be accessed via the Dashboard at materialsproject.org 
#mpr = MPRester()


#lattice = Lattice.cubic(a0)
if len(sys.argv) ==7:
	#a0_x = 3.3
	#a0_y = 4.5
	#lattice = Lattice([[a0_x,0,0],[0,a0_y,0],[0,0, a0_x]]) # To describe uniaxial strain
	lattice = Lattice.cubic(a0)
	#structure = Structure(lattice, [element1,element1,element1,element2], [[0, 0, 0], [0.5, 0.5, 0.0],[0, 0.5,0.5],[0.5,0,0.5]])
	structure = Structure.from_spacegroup("Fm-3m", lattice, [element1, element2], [[0, 0, 0], [0.5, 0.5, 0.5]])
	slabs = generate_all_slabs(structure, max_index=1, min_slab_size=8.0, min_vacuum_size=10.0)
	metal_ml = [slab for slab in slabs if slab.miller_index==MI][0]

if len(sys.argv) ==2: 
	metal_ml = IStructure.from_file(filename)
	print(metal_ml)#

asf_metal_ml = AdsorbateSiteFinder(metal_ml) 
ads_sites = asf_metal_ml.find_adsorption_sites() 
assert len(ads_sites) == 4

fig = plt.figure() 
ax = fig.add_subplot(111) 
plot_slab(metal_ml, ax, adsorption_sites=True)
plt.savefig('slab_%s%s_site.png' % (element1,element2))

from ase.build import molecule as asemole
h  = Molecule("H", [[0,0,0]])
o  = Molecule("O", [[0,0,0]])
o2 = Molecule("OO", [[0, 0, 0],[0,0,1.2]]) 
oh = Molecule("OH", [[0, 0, 0], [-0.793, 0.384, 0.422]])
ooh = Molecule("OOH", [[0, 0, 0], [-1.067, -0.403, 0.796], [-0.696, -0.272, 1.706]])

h2o2_ase      = asemole('H2O2') ## from "ASE" databases
h2o2_position = h2o2_ase.get_positions()
h2o2_species  = h2o2_ase.get_chemical_symbols()
h2o2_sp       = h2o2_species[0]+h2o2_species[1]+h2o2_species[2]+h2o2_species[3]
h2o2          = Molecule(h2o2_sp,h2o2_position)

h2o_ase      = asemole('H2O') ## from "ASE" databases
h2o_position = h2o_ase.get_positions()
h2o_species  = h2o_ase.get_chemical_symbols()
h2o_sp       = h2o_species[0]+h2o_species[1]+h2o_species[2]
h2o          = Molecule(h2o_sp,h2o_position)
adsorbates =[h,o,o2,oh,ooh,h2o2,h2o] 

#adsorbate = Molecule("OO", [[0, 0, 0],[0,0,1.2]]) 
for adsorbate in adsorbates:
        print(adsorbate)
        ads_structs = asf_metal_ml.generate_adsorption_structures(adsorbate, repeat=[2, 2, 1])

        i=0
        while (i < len(ads_structs)):
        	fig = plt.figure() 
        	ax = fig.add_subplot(111) 
        	plot_slab(ads_structs[i], ax, adsorption_sites=False, decay=0.09)
                NAME = [SPECIES.name for SPECIES  in adsorbate.species]
                nnn=''
                for nn in NAME:
                 nnn += nn
                print(nnn)
        	plt.savefig('slab_%s%s_absorb%s_site%i.png' % (element1,element2,nnn,i))
        	Poscar(ads_structs[i]).write_file('POSCAR_%s%s_absorb%s_site%i.vasp' % (element1,element2,nnn,i))
        	i+=1
        #plt.show()

os.system('python a2p2.py')
os.system('mkdir temp')
os.system('mv POSCAR* *.png temp/')
os.system('mv temp %s%s_%i%i%i' % (element1,element2,MI[0],MI[1],MI[2]))

#import glob
#dirs = glob.glob('dir_*')
#for dir in dirs:
#        os.chdir(dir)
#        os.system('python ../a2p2.py')
#        os.chdir('..')
