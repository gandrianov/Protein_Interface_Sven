import sys
import PDBContainer
import RMSDCalculator
import multiprocessing
import gzip, glob
import gc
import tqdm


def read_cif(filename):

    with open(filename, "rb") as f:
        if "gz" in filename:
            f = str(gzip.decompress(f.read()))
        else:
            f = str(f.read())

    return f

def KeywordFilter(string):
    return 'polypeptide' in string and 'Homo sapiens' in string and "polydeoxyribonucleotide" not in string
    
def calculate_RMSD(filename, num_atom=15, rmsd_threshold=10.0):

    try:
        output = {}

        string = read_cif(filename)

        if KeywordFilter(string) is False:
            return {}

        container = PDBContainer.PDBContainer(pdb_string=string)
        
        protein_coords = container.get_protein_coords()
        ligand_coords = container.get_ligand_coords()

        for chain, value in ligand_coords.items():
            for lig_name, atoms in value.items():
                # if lig_name in ["HOH", "ATP" ,"Br","Zn"]:
                #     continue
                # if len(atoms) > num_atom:
                ligand_center_mass = RMSDCalculator.get_center_of_mass(atoms)
                ligand_center_mass = {"Center":ligand_center_mass}
                for protein_chain, prot_coord in protein_coords.items():
                    rmsd = RMSDCalculator.get_RMSD(ligand_center_mass, prot_coord)
                    if rmsd[0] < rmsd_threshold:
                        rmsd = list(rmsd)
                        rmsd[0] = str(rmsd[0])
                        output[f"{lig_name} {chain} {protein_chain}"] = rmsd

        del container, protein_coords, ligand_coords
        collected = gc.collect()
        return output

    except Exception as e:
        print(e, filename)
        return {}


if __name__ == "__main__":

    # filenames = glob.glob("/Users/svenmiller/PDBrenum/mmCIF/*.cif") # take 100 for testing

    filenames = ["/Users/svenmiller/PDBrenum/mmCIF/6tz6.cif.gz"]

    with multiprocessing.Pool(8) as pool:
        # rmsds = pool.map(calculate_RMSD, filenames)
        rmsds = list(tqdm.tqdm(pool.imap(calculate_RMSD, filenames), total=len(filenames)))

        fwr = open(sys.argv[1], "w")
        fwr.write("Filename\tLig:LigChain:ProtChain\tDistance\tProtAtomName\tLigAtomName\n")

        for i, rmsd in enumerate(rmsds):
            # if len(rmsd) >= 2:
            for key, distance in rmsd.items():
                if float(distance[0]) < 10.0:            
                    fwr.write("\t".join([filenames[i], key] + distance) + "\n")

        fwr.close()


    

    # # test for one file
    # filename = "/Users/svenmiller/PDBrenum/mmCIF/6pai.cif"
    # file = read_cif(filename)
    # rmsd = calculate_RMSD(file)

    # fwr = open(sys.argv[1], "w")
    # fwr.write("Filename\tLig:LigChain:ProtChain\tDistance\tProtAtomName\tLigAtomName\n")

    # if len(rmsd) >= 2:
    #     for key, distance in rmsd.items():
    #         if float(distance[0]) < 10.0:            
    #             fwr.write("\t".join([filename, key] + distance) + "\n")

    # fwr.close()

    # print(rmsd)

    # test for multiple files in parallel
    

    # with multiprocessing.Pool(4) as pool:
    #   rmsds = pool.map(calculate_RMSD, filenames)

    #   with open("output.txt", "w") as fwr:
    #       for fname, score in zip(filenames, rmsds):
    #           fwr.write(fname, score, "\n")

    #       fwr.close()