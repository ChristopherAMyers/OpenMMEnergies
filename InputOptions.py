from distutils.util import strtobool


class InputOptions(object):
    '''
        Input file with options that control the functionality
        of energy.py. The format follows that of Q-Chem, with
        the $rem section contianing the main progrma options.

        Comments use a ! character, sections start with a $ followed 
        by a secion title, and all sections conclude with with $end. 
        For example, the main options section would look like this:
        
        $rem
            optimize True       ! optimize coordinates in pdb file. Only prints energy for first frame.
            print_eda           ! don't print energy decomposition, only print total energy
        $end
    '''
    def __init__(self, input_file=None):
        self.density_chg = False
        self.optimize = False
        self.print_eda = True
        self.print_forces = False
        self.opt_mode = 'openmm'
        self.opt_freeze_main = False
        self.opt_freeze_drude = False
        self.nonbonded_eda = False
        self.nonbonded_res_only = False
        self.inter_res_only = False
        self.fragments = {}
        self.compute_forces = False

        if input_file is not None:
            self.read(input_file)


    def read(self, input_file):
        """ Read input file options

            Parameters
            ----------
            input_file: str
                location of the input file
        """
        input_lines = {'rem': [], 'frags': []}
        reading_sec = None
        for line in open(input_file, 'r'):     
            line = line.strip()
            if "$" in line:
                if "$end" in line:
                    reading_sec = None
                else:
                    reading_sec = str.lower(line[1:])
                    input_lines[reading_sec] = []
            
            elif reading_sec is not None:
                input_lines[reading_sec].append(line.replace('=', '').split())


        for line in input_lines['rem']:
            option = line[0].lower()
            value = line[1].lower()
            if option == 'nonbonded_eda':               self.nonbonded_eda = strtobool(value)
            if option == 'nonbonded_res_only':          self.nonbonded_res_only = strtobool(value)
            if option == 'inter_res_only':              self.inter_res_only = strtobool(value)
            if option == 'density_chg':                 self.density_chg = strtobool(value)
            if option == 'optimize':                    self.optimize = strtobool(value)
            if option == 'print_eda':                   self.print_eda = strtobool(value)
            if option == 'frag_opt':                    self.frag_opt = strtobool(value)
            if option == 'opt_freeze_main':             self.opt_freeze_main = strtobool(value)
            if option == 'opt_freeze_drude':            self.opt_freeze_drude = strtobool(value)
            if option == 'opt_mode':
                if value == 'openmm':                   self.opt_mode = 'openmm'
                elif value == 'bfgs':                   self.opt_mode = 'bfgs'
                else:
                    raise ValueError('Invalid value for rem option "opt_mode"')
            if option == 'compute_forces':              self.compute_forces = int(value)

        fragments = {}
        for line in input_lines['frags']:
            name = line[0].lower()
            range_syntex = line[1].lower()

            if name not in fragments:
                fragments[name] = set()


            if "-" in range_syntex:
                sp = range_syntex.strip().split('-')
                start = int(sp[0])
                if len(sp) >= 2:
                    end = int(sp[1])
                else:
                    exit("Error: fragment range must have a START and END id.")
            else:
                start = end = int(range_syntex)

            for i in range(start, end + 1):
                fragments[name].add(i)

        if len(fragments) == 1:
            exit("Error: At least two fragments must be defined")
        elif len(fragments) >= 2:
            for name in fragments:
                fragments[name] = tuple(sorted(fragments[name]))
            self.fragments = fragments



            

        