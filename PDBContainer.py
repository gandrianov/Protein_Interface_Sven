class PDBAtom:
    def __init__(self, x,y,z):
        self.x=x
        self.y=y
        self.z=z


class PDBContainer:
    def __init__(self, filename):
        self.rows = self.read_file(filename)
        self.protein_rows = self.parse_protein_rows(self.rows) # <chain_id>:<list of coordinates>
        self.protein_coordinates = {key:self.parse_atom(value) for key, value in self.protein_rows.items()}
        self.ligand_rows = self.parse_ligand_rows(self.rows)
        # self.ligand_coordinates = {for }

        self.ligand_coordinates = {}

        for key, value in self.ligand_rows.items():
            for key2, value2 in value.items():
                data = self.parse_atom(value2)
                if key not in self.ligand_coordinates:
                    self.ligand_coordinates[key] = {}
                if key2 not in self.ligand_coordinates[key]:
                    self.ligand_coordinates[key][key2] = data

        # print(self.ligand_coordinates)

    def read_file(self, filename):
        with open(filename, 'r') as f:
            rows= f.read().split('\n')
        
        return rows

    def parse_protein_rows(self, rows):

        protein_atom = {}
        for row in rows:
            if row.startswith('ATOM'):
                chain_name = [e for e in row.split(" ") if e != ""]
                chain_name = chain_name[6] + chain_name[-3]
                if chain_name not in protein_atom:
                    protein_atom[chain_name] = [row]
                else:
                    protein_atom[chain_name].append(row)
                
        return protein_atom

    def parse_ligand_rows(self, rows):
        hets = {}
        for row in rows:
            if row.startswith('HETATM'):
                chain_name = [e for e in row.split(" ") if e != ""]
                chain_name = chain_name[6] + chain_name[-3]
                ligand_name = [e for e in row.split(" ") if e != ""][5]
                # print(ligand_name)

                if chain_name not in hets:
                    hets[chain_name] = {}
                if ligand_name not in hets[chain_name]:
                    hets[chain_name][ligand_name] = []
                
                hets[chain_name][ligand_name].append(row)

        return hets

    def parse_atom(self, rows):
        coordinates = {}
        for row in rows: 
            x,y,z = [e for e in row.split(" ") if e != ""][10:13]
            x,y,z = float(x), float(y), float(z)
            atom_name = [e for e in row.split(" ") if e != ""][1:4:2]
            atom_name = "_".join(atom_name)
            # print(atom_name,x,y,z)
            atom=PDBAtom(x,y,z)
            coordinates[atom_name]=atom

        return coordinates

    def get_ligand_coords(self):
        return self.ligand_coordinates

    def get_protein_coords(self):
        return self.protein_coordinates


if __name__ == '__main__':

    # atom=PDBAtom(0.1, 0.2, 0.3)

    container =  PDBContainer("/Users/svenmiller/PDBrenum/mmCIF/6pai.cif")
    # print(container.protein_rows.keys())