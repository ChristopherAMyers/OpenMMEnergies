Create a Q-Chem job that produces a potential energy scan for torsion angle rotations.\

#### Input file formatFormat
This script takes in a single Q-Chem template file with an altered version of the $scan section.\
The first line of $scan follows the Q-Chem default for torsion scans (atom1, atom2, atom3, atom4, angle_start, angle_end, interval).\
The following lines are single valued and include the indicies of the atoms which to rotate about.\
\
See the *example_template.in* file as a model for this format. An example of the generated outputs are *example_run.in* and *example_coords.xyz*
