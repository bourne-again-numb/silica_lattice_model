# silica_lattice_model
Silica Lattice Model

v2.8
first neighbor hard-sphere TAA-TAA repulsions
second neighbor TAA-ISIO(SI-) attractions

This model represents TMA model.


TAA => 	hardsphere first neighbor repulsions with every molecule,
	hardsphere second neighbor repulsions TAA-TAA
	hardsphere third neighbor repulsions TAA-TAA	

SI => 	second neighbor hardsphere repulsion between SI-SI
	third neighbor harsphere repulsion between SI-SI

a change in the ring counting algorithm explicitly moving from Si-O-Si, to

check for the 3 and 4 membered rings. 30 % faster than v2.6

input of nTEOS,nTAAOH,nH2O instead of [SI],[SN],[TAA] concentrations

generation of three q distributions,one for SN, one for SI and one for both


