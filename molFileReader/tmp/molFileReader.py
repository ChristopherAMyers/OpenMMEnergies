import numpy as np
from numpy import pi
import os, sys

class _Components:
    def __init__(self):
        self.atoms = []
        self.coords = []
        self.comment = ""
        self.n_atoms = 0

class _GroComps():
    def __init__(self):
        self.coords = []
        self.comment = ""
        self.n_atoms = 0
        self.res_nums = []
        self.res_names = []
        self.vels = []
        self.names = []
        self.box = []
        self.atom_num = []
        self.time = None

class _QC_Components:
    def __init__(self):
        self.atoms = np.empty(0)
        self.coords = np.empty((0, 3))
        self.charge = 0
        self.multiplicity = 0
        self.n_atoms = 0

class QC(_Components):
    def __init__(self):
        self.total_charge = 0
        self.total_multiplicity = 1
        self.fragments = list([_QC_Components()])
        self.fragments.clear()
        self.n_frags = 0
        self.n_atoms = 0
        self.other_lines = []

    def import_qc(self, fileLoc):
        if os.path.isfile(fileLoc):
            with open(fileLoc, 'r') as file:
                read_mol = False
                line = file.readline()
                found_mol = False
                other_lines = []
                while line:
                    sp = line.split()
                    if len(sp) == 0:
                        other_lines.append(line)
                        line = file.readline()
                        continue

                    #   found end of molecule section
                    if sp[0].lower() == '$end' and found_mol:
                        read_mol = False
                        line = file.readline()
                        continue

                    #   separator signals a new fragment
                    if sp[0] in ['-'*2, '-'*3, '-'*4]:

                        self.fragments.append(_QC_Components())
                        sp = file.readline().split()

                        self.fragments[-1].charge = int(sp[0])
                        self.fragments[-1].multiplicity = int(sp[1])
                        line = file.readline()
                        sp = line.split()
                        read_mol = True
                    
                    #   record coordinates and element and add to current fragment
                    if read_mol:
                        coords = np.array([[float(x) for x in sp[1:4]]])
                        self.fragments[-1].coords = np.append(self.fragments[-1].coords, coords, axis=0)
                        self.fragments[-1].atoms = np.append(self.fragments[-1].atoms, sp[0])
                        self.n_atoms += 1
                        line = file.readline()

                    #   start of molecule section 
                    elif sp[0].lower() == '$molecule':
                        found_mol = True
                        line = file.readline()
                        sp = line.split()
                        if len(sp) != 2:
                            print(" Wrong format for charge and multiplicity line")
                            exit()
                        self.total_charge = int(sp[0])
                        self.total_multiplicity = int(sp[1])
                        read_mol = True

                        line = file.readline()
                        sp = line.split()
                        #   if next line is not a separator, we are not using multiple fragments
                        if sp[0] not in ['-'*2, '-'*3, '-'*4]:
                            read_mol = True
                            self.fragments.append(_QC_Components())
                            self.fragments[-1].charge = self.total_charge
                            self.fragments[-1].multiplicity = self.total_multiplicity
                        


                    else:
                        other_lines.append(line)
                        line = file.readline()

            self.n_frags = len(self.fragments)
            self.other_lines = other_lines
                

        else:
            raise FileNotFoundError(' QC file not found. Terminating program')

    def input_lines(self):
        lines = []
        lines.append('$molecule\n')
        lines.append('{:d} {:d}\n'.format(self.total_charge, self.total_multiplicity))
        for n in range(self.n_frags):
            frag = self.fragments[n]
            if self.n_frags > 1:
                lines.append('---\n')
                lines.append('{:d} {:d}'.format(frag.charge, frag.multiplicity))
            for i, coord in enumerate(frag.coords):
                lines.append('{:2s}  {:15.8f}  {:15.8f}  {:15.8f}\n'.format(frag.atoms[i], *(coord)))
        lines.append('$end\n')

        return lines
            
    

class XYZ(_Components):
    def __init__(self):
        super(XYZ, self).__init__()
        self.frames = list([_Components()])
        self.frames.clear()

    def set_frame_num(self, idx):
        self.coords = self.frames[idx].coords
        self.atoms= self.frames[idx].atoms
        self.n_atoms = self.frames[idx].n_atoms
        self.comment = self.frames[idx].comment

    def import_xyz(self, fileLoc):
        if os.path.isfile(fileLoc):
            with open(fileLoc, 'r') as file:
                count = 2
                numAtoms = 0
                for line in file.readlines():
                    
                    #   if reached the number of atoms, add to frames
                    if count == numAtoms + 2:
                        if len(self.frames) != 0:
                            if self.frames[-1].n_atoms != len(self.frames[-1].atoms):
                                print("WARNING: Found " + str(len(self.frames[-1].atoms)) + " atoms in XYZ file:")
                                print("         " + fileLoc)
                                print("         but should be " + str(self.frames[-1].n_atoms))
                                self.frames[-1].n_atoms = len(self.frames[-1].atoms)
                        count = 0

                    #   exit import if a blank line is found
                    if line not in ['\n', '\r\n']:
                        if count == 0:
                            self.frames.append(_Components())
                            line = line.split()
                            self.frames[-1].n_atoms = int(line[0])
                            numAtoms = int(line[0])
                        elif count == 1:
                            self.frames[-1].comment = line.strip('\n')
                        else:
                            line = line.split()
                            self.frames[-1].atoms.append(line[0])
                            self.frames[-1].coords.append([float(line[1]), float(line[2]), float(line[3])])
                        count += 1
                    else:
                        break
            

            self.set_frame_num(0)
            
            self.atoms = np.array(self.atoms)
            self.coords = np.array(self.coords)
        else:
            print("ERROR: File ", fileLoc, " not found!")
            print("Program will now terminate")
            exit()


    def add_atom(self, coord, atom):
        #if not isinstance(self.coords, list):
        #    print("ERROR: Molecule coordinates have been convered form a list")
        #    exit()
        #if not isinstance(coord, list):
        #    print("ERROR: new coordinate must be of type list")
        #    exit()
        self.coords.append(list(coord))
        self.atoms.append(atom)
        self.n_atoms += 1

    def add_frame(self, coords = np.array([]), comment = None):
        self.frames.append(_Components())
        if len(coords) == 0:
            self.frames[-1].coords = self.coords
        else:
            self.frames[-1].coords = np.copy(coords)
        self.frames[-1].atoms = self.atoms
        self.frames[-1].n_atoms = self.n_atoms
        if comment == None:
            self.frames[-1].comment = self.comment
        else:
            self.frames[-1].comment = comment
    
    def recenter(self, vec):
        '''
        re-center coordinates so that the vector 'vec' is at the origin
        '''
        self.coords = center(self.coords, vec)

    def rotate_plane(self, vec1, vec2):
        '''
        rotate coordinates so that 'vec1' and 'vec2' are in the xy-plane
        '''
        self.coords = orient_plane(self.coords, vec1, vec2)

    def rotate_plane_xz(self, vec1, vec2):
        '''
        rotate coordinates so that 'vec1' and 'vec2' are in the xz-plane
        '''
        self.coords = orient_plane_xz(self.coords, vec1, vec2)

    def rotate_plane_yz(self, vec1, vec2):
        '''
        rotate coordinates so that 'vec1' and 'vec2' are in the xz-plane
        '''
        self.coords = orient_plane_yz(self.coords, vec1, vec2)

    def write(self, fileName):
        '''
        write coordinates to a new XYZ file
        '''
        for n in range(len(self.frames)):
            if n == 0:
                write_xyz(self.frames[n].atoms, self.frames[n].coords, fileName, comment=self.frames[n].comment)
            else:
                write_xyz(self.frames[n].atoms, self.frames[n].coords, fileName, filemode='a', comment=self.frames[n].comment)
        #write_xyz(self.atoms, self.coords, fileName)
        

class GRO(_Components):
    def __init__(self):
        super(GRO, self).__init__()
        self.frames = list([_GroComps()])
        self.frames.clear()

    def set_frame_num(self, idx):
        self.coords = self.frames[idx].coords
        self.n_atoms = self.frames[idx].n_atoms
        self.comment = self.frames[idx].comment
        self.res_names = self.frames[idx].res_names
        self.res_nums = self.frames[idx].res_nums
        self.names = self.frames[idx].names
        self.box = self.frames[idx].box
        self.time = self.frames[idx].time
        self.atom_num = self.frames[idx].atom_num

    def _parse_add_line(self, line):
        n = -1
        self.frames[n].res_nums.append(int(line[:5]))
        line = line[5:]

        self.frames[n].res_names.append(line[:5].split()[0])
        line = line[5:]

        self.frames[n].names.append(line[:5].split()[0])
        line = line[5:]

        self.frames[n].atom_num.append(int(line[:5]))
        line = line[5:]

        x = float(line[:8])*10
        line = line[8:]
        y = float(line[:8])*10
        line = line[8:]
        z = float(line[:8])*10
        line = line[8:]
        self.frames[n].coords.append([x, y, z])

        if len(line) > 0 and line not in ['\n', '\r\n']:
            vx = float(line[:8])
            line = line[8:]
            vy = float(line[:8])
            line = line[8:]
            vz = float(line[:8])
            line = line[8:]
            self.frames[n].vels.append([vx, vy, vz])
        else:
            self.frames[n].vels.append([0.0, 0.0, 0.0])

        

    def import_gro(self, fileLoc):
        if os.path.isfile(fileLoc):
            with open(fileLoc, 'r') as file:
                count = 3
                numAtoms = 0
                for line in file.readlines():
                    
                    #   if reached the number of atoms, add to frames
                    if count == numAtoms + 3:
                        if len(self.frames) != 0:
                            if self.frames[-1].n_atoms != len(self.frames[-1].coords):
                                print("WARNING: Found " + str(len(self.frames[-1].coords)) + " atoms in XYZ file:")
                                print("         " + fileLoc)
                                print("         but should be " + str(self.frames[-1].n_atoms))
                                self.frames[-1].n_atoms = len(self.coords[-1].atoms)
                        count = 0

                    #   exit import if a blank line is found
                    if line not in ['\n', '\r\n']:
                        if count == numAtoms + 2:
                            line = line.split()
                            line = [float(x) for x in line]
                            self.frames[-1].box = line[:3]
                        elif count == 0:
                            self.frames.append(_GroComps())
                            line = line.split('t=')
                            self.frames[-1].comment = line[0]
                            if len(line) == 2:
                                self.frames[-1].time = float(line[1])
                        elif count == 1:
                            line = line.split()
                            self.frames[-1].n_atoms = int(line[0])
                            numAtoms = int(line[0])
                        
                        else:
                            self._parse_add_line(line)
                        count += 1
                    else:
                        break
            

            self.set_frame_num(0)
            
            self.atoms = np.array(self.atoms)
            self.coords = np.array(self.coords)
        else:
            print("ERROR: File ", fileLoc, " not found!")
            print("Program will now terminate")
            exit()

def center(coord_list, centCoord):
    '''
    re-centers atomic coordinates 'coord_list' by 'centCoord'.
    'coord_list' can be a list or a numpy array, but will
    return a numpy array
    '''
    newCoords = coord_list.copy()
    for n in np.arange(0, len(newCoords)):
        newCoords[n] -= centCoord
    return np.array(newCoords)

def rotate_x(coords, phi):
    '''
    rotation about x-axis by angle 'phi'
    '''
    newCoords = coords.copy()
    dim = len(coords)

    cosP = np.cos(phi)
    sinP = np.sin(phi)
    rotX = np.array([
        [1, 0, 0],
        [0, cosP, -sinP],
        [0, sinP, cosP]
    ])
    for n in np.arange(0, dim):
        newCoords[n] = np.matmul(rotX, coords[n])
    return newCoords

def rotate_y(coords, phi):
    '''
    rotation about y-axis by angle 'phi'
    '''
    newCoords = coords.copy()
    dim = len(coords)

    #normXY = np.sqrt(vec[0]**2 + vec[1]**2)
    #cosP = vec[0]/normXY
    #sinP = vec[1]/normXY
    cosP = np.cos(phi)
    sinP = np.sin(phi)
    rotYm = np.array([
        [cosP, 0, sinP],
        [0, 1, 0],
        [-sinP, 0, cosP]
    ])
    for n in np.arange(0, dim):
        newCoords[n] = np.matmul(rotYm, coords[n])
    return newCoords

def rotate_z(coords, phi):
    '''
    rotation about z-axis by angle 'phi'
    '''
    newCoords = coords.copy()
    dim = len(coords)

    #cosT = vec[3]/np.linalg.norm(vec)
    #sinT = np.sqrt(1 - cosT**2
    cosT = np.cos(phi)
    sinT = np.sin(phi)
    rotZm = np.array([
        [cosT, -sinT, 0],
        [sinT, cosT, 0],
        [0, 0, 1]
    ])
    for n in np.arange(0, dim):
        newCoords[n] = np.matmul(rotZm, coords[n])

    return newCoords

def get_phi_spherical(vec):
    x, y, z = np.array(vec)
    if abs(x) <= 1E-8 and y > 0:
        phi = pi/2
    elif abs(x) <= 1E-8 and y < 0:
        phi = 3*pi/2
    else:
        phi = np.arctan2(y, x)

    return phi

def orient_plane(coords, vec1, vec2):
    '''
    reorients 'coords' so that 'vec1' is along x-axis
    and 'vec2' is in xy-plane
    '''
    allCoords = np.append(coords, np.array([vec1]), axis=0)
    allCoords = np.append(allCoords, np.array([vec2]), axis = 0)

    phi = np.arctan2(allCoords[-2][1], allCoords[-2][0])
    allCoords = rotate_z(allCoords, -phi)

    theta = np.pi*0.5 - np.arccos(allCoords[-2][2] / np.linalg.norm(allCoords[-2]))
    allCoords = rotate_y(allCoords, theta)

    theta = np.pi*0.5 - np.arccos(allCoords[-1][2] / np.sqrt(allCoords[-1][1]**2 + allCoords[-1][2]**2))
    tmp = rotate_x(allCoords, -theta)
    if abs(tmp[-1][2]) > 1E-8:
        allCoords = rotate_x(allCoords, theta)
    else:
        allCoords = tmp

    return allCoords[0:-2]

def orient_plane_xz(coords, vec1, vec2):
    '''
    reorients 'coords' so that 'vec1' is along z-axis
    and 'vec2' is in xz-plane
    '''

    allCoords = np.append(coords, np.array([vec1]), axis=0)
    allCoords = np.append(allCoords, np.array([vec2]), axis = 0)

    phi = get_phi_spherical(allCoords[-2])
    allCoords = rotate_z(allCoords, -phi + pi*0.5)

    theta = np.arccos(allCoords[-2][2] / np.linalg.norm(allCoords[-2]))
    allCoords = rotate_x(allCoords, -theta + pi)

    phi = get_phi_spherical(allCoords[-1])
    allCoords = rotate_z(allCoords, -phi)

    return allCoords[0:-2]

def orient_plane_yz(coords, vec1, vec2):
    '''
    reorients 'coords' so that 'vec1' is along z-axis
    and 'vec2' is in yz-plane
    '''

    allCoords = np.append(coords, np.array([vec1]), axis=0)
    allCoords = np.append(allCoords, np.array([vec2]), axis = 0)

    phi = get_phi_spherical(allCoords[-2])
    allCoords = rotate_z(allCoords, -phi + pi*0.5)

    theta = np.arccos(allCoords[-2][2] / np.linalg.norm(allCoords[-2]))
    allCoords = rotate_x(allCoords, -theta + pi)

    phi = get_phi_spherical(allCoords[-1])
    allCoords = rotate_z(allCoords, -phi + pi/2)

    return allCoords[0:-2]

def import_qmol(fileLoc):
    '''
    import Q-Chem molecule
    '''
    mols = []
    coords = []
    atoms = []
    parse = False
    count = 0
    fragCount = 0
    fragChg = 0
    totalChg = 0
    with open(fileLoc, 'r') as file:
        for line in file.readlines():
            sp = line.split()
            if len(sp) > 0:
                if not(parse) and sp[0].lower() == "$molecule":
                    parse = True
                    count = 1
                    fragCount = 1
                elif parse and sp[0].lower() == "$end":
                    parse = False
                    mols.append([fragChg, atoms, coords])
                elif parse:
                    if count == 1:
                        totalChg = int(sp[0])
                        fragChg = totalChg
                        count = 2
                        fragCount = 2
                    elif count > 1 and sp[0] == "--":
                        if len(coords) != 0:
                            mols.append([fragChg, atoms, coords])
                        coords = []
                        atoms = []
                        fragCount = 1
                    elif fragCount == 1:
                        fragChg = int(sp[0])
                        fragCount += 1
                    elif len(sp) == 4:
                        atoms.append(sp[0])
                        coords.append([float(x) for x in sp[1:4]])

    print("Number of fragments: ", len(mols))
    for n in range(0, len(mols)):
        print("Fragment 1")
        print("\tNumber of atoms: {:3d}".format(len(mols[n][1])))
        print("\tTotal charge:    {:3d}".format(mols[n][0]))
    return mols

def write_xyz(atoms, coords, xyzFile, filemode = 'w', comment="Generated xyz file"):
    '''
    write 'atoms' and 'coords' to a
    new XYZ file 'xyzFile'
    '''
    n_atoms = len(atoms)
    if n_atoms == 0:
        return
    with open(xyzFile, filemode) as file:
        file.write(str(int(n_atoms)) + "\n")
        if comment == "":
            file.write("Generated xyz file\n")
        else:
            file.write(comment + " \n")
        for n in range(0, n_atoms):
            print(coords)
            file.write("{:3s}  {:13.8f}  {:13.8f}  {:13.8f}\n"
                .format(atoms[n], coords[n][0], coords[n][1], coords[n][2]))

def write_qcin(atoms, coords, charges, qmolFile, comment = ""):
    '''
    write a new Q-Chem file
    '''
    atoms1, atoms2 = atoms
    coords1, coords2 = coords
    chg1, chg2 = charges
    with open(qmolFile, 'w') as file:

        if comment != "":
            file.write("$comment\n")
            file.write(comment)
            file.write("\n$end\n\n")


        file.write("$molecule\n")
        file.write(str(chg1 + chg2) + " 1\n")

        file.write("-- Frag 1\n")
        file.write(str(chg1) + " 1\n")
        for n in range(0, len(atoms1)):
            file.write("{:3s}  {:13.8f}  {:13.8f}  {:13.8f}\n"
                .format(atoms1[n], coords1[n][0], coords1[n][1], coords1[n][2]))

        file.write("-- Frag 2\n")
        file.write(str(chg2) + " 1\n")
        for n in range(0, len(atoms2)):
            file.write("{:3s}  {:13.8f}  {:13.8f}  {:13.8f}\n"
                .format(atoms2[n], coords2[n][0], coords2[n][1], coords2[n][2]))

        file.write("$end\n")
        file.write("\n\n$rem\n")
        file.write("  METHOD   PBE\n")
        file.write("  DFT_D    D3\n")
        file.write("  JOBTYPE   BSSE\n")
        file.write("  BASIS   pcseg-2\n")
        file.write("$end\n")


if __name__ == "__main__":
    import_qmol(sys.argv[1])
