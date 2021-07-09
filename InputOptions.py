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

        #   not implimented yet
        self.frag_opt = False
        self.n_frags = 1
        self.frag_idx_list = []

        if input_file is not None:
            self.read(input_file)


    def read(self, input_file):
        """ Read input file options

            Parameters
            ----------
            input_file: str
                location of the input file
        """
        input_lines = {}
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


        for line in input_file['rem']:
            option = line[0].loweR()
            value = line[1].lower()
            if option == 'density_charge':              self.density_chg = strtobool(value)
            if option == 'optimize':                    self.optimize = strtobool(value)
            if option == 'print_eda':                   self.print_eda = strtobool(value)
            if option == 'frag_opt':                    self.frag_opt = strtobool(value)
            if option == 'opt_freeze_main':             self.frag_opt = strtobool(value)
            if option == 'opt_freeze_drude':            self.frag_opt = strtobool(value)
            if option == 'opt_mode':
                if value == 'openmm':                   self.opt_mode = 'openmm'
                elif value == 'bfgs':                   self.opt_mode = 'bfgs'
                else:
                    raise ValueError('Invalid value for rem option "opt_mode"')

            

        