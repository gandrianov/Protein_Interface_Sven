import PDBContainer
import RMSDCalculator
import multiprocessing

def calculate_RMSD(filename):

	try:
		container = PDBContainer.PDBContainer(filename)
		
		ligand_coords = container.get_ligand_coords()
		ligand_center_mass = RMSDCalculator.get_center_of_mass(ligand_coords)

		protein_coords = container.get_protein_coords()

		return RMSDCalculator.get_RMSD(ligand_center_mass, protein_coords)
	except Exception as e:
		print(e)
		return None


if __name__ == "__main__":

	# test for one file
	filename = ""
	rmsd = calculate_RMSD(filename)
	print(rmsd)

	# test for multiple files in parallel
	filenames = glob.glob("")[:100] # take 100 for testing

	with multiprocessing.Pool(4) as pool:
		rmsds = pool.map(calculate_RMSD, filenames)

		with open("output.txt", "w") as fwr:
			for fname, score in zip(filenames, rmsds):
				fwr.write(fname, score, "\n")

			fwr.close()