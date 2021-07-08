import numpy as np
import sys
from math import acos, cos, sin

def parse_qc_lines(qc_lines):
    '''
        Get Q-Chem input sections and their inlcuded options

        Parameters
        ----------
        qc_lines: str
            file location of input file

        Returns
        -------
        dict
            parsed input options. Each key is the section tile and it's
            contents are a list of lists of the options provided
    '''
    qm_lines = {}

    reading_sec = None
    for line in qc_lines:
        line = line.strip()
        if "$" in line:
            if "$end" in line:
                reading_sec = None
            else:
                reading_sec = str.lower(line[1:])
                qm_lines[reading_sec] = []
        
        elif reading_sec is not None:
            
            qm_lines[reading_sec].append(line.replace('=', '').split())

    return qm_lines

def get_torsion_angle(A, B, C, D):
    """
        Calculate the torsion angle given 4 atom
        numpy vectors A, B, C, D

        Parameters
        ----------
        A, B, C, D: np.array
            4 vectors that define the torsion angle

        Returns
        -------
        float
            torsion angle of the 4 vector
        
        Math taken from
        https://charmm-gui.org/?doc=lecture&module=scientific&lesson=9
    """

    C_ba = B - A
    C_bc = B - C
    C_cb = C - B
    C_cd = C - D

    cross_1 = np.cross(C_ba, C_bc)
    cross_2 = np.cross(C_cb, C_cd)
    cross_1 /= np.linalg.norm(cross_1)
    cross_2 /= np.linalg.norm(cross_2)
    angle = acos(np.dot(cross_1, cross_2))

    C_ijk = np.cross(C_ba, C_cd)
    if np.dot(C_ijk, C_bc) < 0:
        return angle
    else:
        return -angle

def rodrigues_rotation(vec, rot_k, angle):
    """
        Rotates a vector 'vec' about the vector 'rot_k' by 'angle'

        Parameters
        ----------
        vec: np.array
            vector in R3 to rotate
        rot_k: np.array
            unit vector in R3 to rotate about
        angle: float
            angle in radians to rotate

        Returns
        -------
        np.array
            rotated vector
    """
    cos_t = cos(angle)
    sin_t = sin(angle)
    dot = np.dot(rot_k, vec)
    cross = np.cross(rot_k, vec)
    return vec*cos_t + cross*sin_t + rot_k*dot*(1 - cos_t)

if __name__ == '__main__':\
    #   import Q-Chem input file
    qm_lines = parse_qc_lines(open(sys.argv[1], 'r').readlines())

    #   search for $scan section and assume torsion scan is requested
    scan_lines = qm_lines['scan']
    scan_info = scan_lines[0]
    idx_numbers = [int(x) - 1 for x in scan_info[1:5]]
    start_angle = float(scan_info[5])
    end_angle =   float(scan_info[6])
    interval =    float(scan_info[7])
    scan_idx_list = [int(x[0]) - 1 for x in scan_lines[1:]]
    print(" Start angle: {:10.2f}".format(start_angle))
    print(" End angle:   {:10.2f}".format(end_angle))
    print(" Interval:    {:10.2f}".format(interval))
    print(" No. of atoms to rotate: {:3d}".format(len(scan_idx_list)))

    #   get initial torsion angle
    coords = np.array([np.array(x[1:4], dtype=float) for x in qm_lines['molecule'][1:] ])
    atoms = np.array([np.array(x[0]) for x in qm_lines['molecule'][1:] ])
    idx_i, idx_j, idx_k, idx_l = idx_numbers
    init_angle = get_torsion_angle(coords[idx_i], coords[idx_j], coords[idx_k], coords[idx_l])*180/np.pi
    
    print("\n Initial torsion angle: {:8.2f} degrees".format(init_angle))
    
    #   atoms will be rotated about this vector direction
    rot_vec = coords[idx_k] - coords[idx_j]
    rot_vec /= np.linalg.norm(rot_vec)

    #   xyz file of all rotations
    coord_file = open('all_coords.xyz', 'w')

    #   loop over angles and create input file
    with open('run.in', 'w') as file:
        angles = np.arange(start_angle, end_angle + interval, interval)
        print(" There are a total of {:d} angles".format(len(angles)))
        for n, angle in enumerate(angles):
            radians = angle*np.pi/180.0

            #   recenter by atom j
            new_coords = np.copy(coords) - coords[idx_j]
            #   rotate the specified atoms 
            for idx in scan_idx_list:
                new_angle = (angle - init_angle)*np.pi/180.0
                new_coords[idx] = rodrigues_rotation(new_coords[idx], rot_vec, new_angle)
            #   put back to original center
            new_coords += coords[idx_j]

            #   write header for coord file
            coord_file.write('{:d}\n'.format(len(coords)))
            coord_file.write('Angle = {:10.6f} \n'.format(angle))

            #   add current angle as a comment
            file.write('\n\n$comment\n')
            file.write('    Angle = {:10.6f} \n'.format(angle))
            file.write('$end \n\n')

            #   write rem options
            file.write('$rem\n')
            if n > 0:
                #   all other angles uses the previous solved orbitals as an initial guess
                file.write("    {:20s}{:20s}\n".format('scf_guess', 'read'))
            for line in qm_lines['rem']:
                file.write("    {:20s}{:20s}\n".format(line[0], line[1]))
            file.write('$end\n\n')

            #   write molecule section
            file.write('$molecule\n')
            charge, mult = qm_lines['molecule'][0]
            file.write('{:3s}{:2s}\n'.format(charge, mult))
            for n, coord in enumerate(new_coords):
                file.write("{:2s}  {:15.8f}  {:15.8f}  {:15.8f}\n".format(atoms[n], *tuple(coord)))
                coord_file.write("{:2s}  {:15.8f}  {:15.8f}  {:15.8f}\n".format(atoms[n], *tuple(coord)))
            file.write('$end\n\n')

            #   write other Q-Chem sections
            for key in qm_lines.keys():
                if key in ['rem', 'molecule', 'scan']: continue
                file.write('${:s}\n'.format(key))
                for line in qm_lines[key]:
                    file.write("{:20s}{:20s}\n".format(line[0], line[1]))
                file.write('$end\n\n')

            #   flag for consecutive jobs
            if n != len(angles) - 1:
                file.write('\n@@@\n\n')

    coord_file.close()
    print(" Q-Chem input file written to run.in")
    print(" All coordinates printed to all_coords.xyz")
    print("")
    print(" Done. Good Luck!")
