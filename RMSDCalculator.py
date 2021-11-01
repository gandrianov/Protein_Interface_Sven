import PDBContainer

def get_center_of_mass(atoms):
	avg_x, avg_y, avg_z = 0.0 , 0.0 , 0.0

	for atom_name, atom_pdb in atoms.items():
		avg_x += atom_pdb.x
		avg_y += atom_pdb.y
		avg_z += atom_pdb.z

	avg_x/= len(atoms)
	avg_y/= len(atoms)
	avg_z/= len(atoms)

	return PDBContainer.PDBAtom(avg_x, avg_y, avg_z)

def get_RMSD(lig_coords, protein_coords):
	
	best_lig_atom = None
	best_prot_atom = None
	RMSD = 99999

	for lig_name, lig_coord in lig_coords.items():
		for prot_name, prot_coord in protein_coords.items():
			dx = (lig_coord.x - prot_coord.x) ** 2
			dy = (lig_coord.y - prot_coord.y) ** 2
			dz = (lig_coord.z - prot_coord.z) ** 2

			d = (dx + dy + dz) ** 0.5

			# print(d, prot_name)

			if d < RMSD:
				best_lig_atom = lig_name
				best_prot_atom = prot_name
				RMSD = d

	return RMSD, best_prot_atom, best_lig_atom



if __name__ == '__main__':

    container =  PDBContainer.PDBContainer("/Users/svenmiller/PDBrenum/mmCIF/6pai.cif")

    ligand_coords = container.get_ligand_coords()
    protein_coords = container.get_protein_coords()

    for lig_chain, value in ligand_coords.items():
    	for ligand_name, atom_coords in value.items():
    		# center_mass = {"Center":get_center_of_mass(atom_coords)}
    		for protein_chain, prot_coords in protein_coords.items():

    			print("====" * 20)

    			# rmsd, atom = get_RMSD(center_mass, prot_coords)
    			rmsd, patom, latom = get_RMSD(atom_coords, prot_coords)
    			print(lig_chain, ligand_name, protein_chain, rmsd, len(atom_coords), patom, latom)
    			

    			print("====" * 20)

    		# 	print(lig_chain, ligand_name, protein_chain, rmsd, len(atom_coords), atom)




	# ligand_center_mass = get_center_of_mass(ligand_coords)