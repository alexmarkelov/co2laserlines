'''
wavelength module includes 2 classes:
Line for storing information about emission line of 12C16O2
molecule.
Branch contains number of the Lines related to one of the four
emission branches of CO2 molecule.

Wavelength of CO2 emission line can be obtained from Vitterman
book or calculated.
'''

import molconst12c16o2

VITTERMANN_FILE_NAME = 'CO2laserlinesVitterman.dat'
LIGHT_VELOCITY = 299792458
branch_list = ['P', 'R']
center_lines_list = [10, 9]
mol_const_dict = {
                  '10': tuple(molconst12c16o2.CO210MKM),
                  '9': tuple(molconst12c16o2.CO29MKM)
                 }

jnum_const_dict = {
                    'P': 'const_branch_P',
                    'R': 'const_branch_R'
                  }
unit_selection = {
                  'cm-1': 'lambda x: x',
                  'mkm': 'convert_rcm_in_mkm',
                  'hz': 'convert_rcm_in_hz'
                  }


source_selection = {
                    'calc': 'calculate_wavelength',
                    'vitt': 'get_wavelength_vitt'
                   }


vitt_const = {
              '10': 0,
              '9': 2,
              'P': 0,
              'R': 1
             }


def unit_checker(unit):
    ''''Verify the wavelength unit.'''
    if unit in unit_selection.keys():
        return True
    return False


def line_parameter_checker(center: int, branch: str,
                           jnum: int) -> bool:
    '''Check the correctness of line parameter: center = (9, 10),
    branch = ('P', 'R'), jnum = range(2, 58, 2)
    '''
    if center in center_lines_list and branch in branch_list and \
       jnum in range(2, 60, 2):
        return True
    return False


def const_branch_P(j_num: int) -> int:
    return -j_num


def const_branch_R(j_num: int) -> int:
    return j_num + 1


def get_line_name(center: int, branch: str, jnum: int) -> str:
    '''Get line name from line parameter(center, branch, jnum)'''
    line_name = str(center) + branch + str(jnum)
    return line_name


def get_line_center_branch(line) -> str:
    '''Get str consist of line center and branch'''
    center_branch = str(line.center) + line.branch
    return center_branch


def convert_rcm_in_mkm(wavelength_rcm: float) -> float:
    '''Convert cm-1 in mkm'''
    wavelength_in_mkm = 10000.0/float(wavelength_rcm)
    return wavelength_in_mkm


def convert_rcm_in_hz(wavelength_rcm: float) -> float:
    '''Convert wavelength[cm-1] in frequency[hz]'''
    frequency_hz = LIGHT_VELOCITY*wavelength_rcm*100
    return frequency_hz


def get_vitt_num(center: int, branch: str, jnum: int) -> int:
    '''Get position in list from file with wavelengths'''
    number = vitt_const[str(center)] + vitt_const[branch] + (jnum/2 - 1)*4
    return int(number)


def get_wavelength_vitt(center: int, branch: str, jnum: int, *args,
                        file_name: str = VITTERMANN_FILE_NAME) -> float:
    '''Get wavelength from file'''
    with open(file_name, 'r') as lines_Vitterman:
        list_wavelength = lines_Vitterman.read().split()
        num = get_vitt_num(center, branch, jnum)
        wavelength = float(list_wavelength[num].replace(',', '.'))
    return wavelength


def get_branch_const(center, branch, jnum) -> tuple:
    '''Get constants for calculation wavelength'''
    constants = [eval(jnum_const_dict[branch])(jnum),
                     *mol_const_dict[str(center)]]
    branch_const = tuple(constants)
    return branch_const


def wavelength_calculator(line_const: tuple) -> float:
    '''Calculator of wavelength using tuple of constants'''
    wavelength = 0
    line_const_list = list(line_const)
    jnum = line_const_list.pop(0)
    for number, constanta in enumerate(line_const_list):
        wavelength += constanta*jnum**number
    return wavelength


def calculate_wavelength(center: int=10, branch: str='P',
                         jnum: int=20) -> float:
    '''Calculate wavelength of CO2 laser emission line by using formula
    from paper "Infrared energy levels and intensities of carbon dioxide-II."
    L.S Rothman, L. D. P. Young, J. Quant. Spect. Radiat. Transfer. 1981
    '''
    if line_parameter_checker(center, branch, jnum):
        branch_const = get_branch_const(center, branch, jnum)
        wavelength = wavelength_calculator(branch_const)
        return wavelength
    return None


class Line:
    '''
    Line is a class for description vibrational-rotational transition
    in CO2 molecula. Line includes number parameters of the line,
    center - central wavelength of transition(10 or 9), branch - branch
    of the transition('P' or 'R'), jnum - rotational number of lower
    energy level, unit - unit of wavelength, source - calculated or vitterman.
    '''
    def __init__(self, center: int=10, branch: str='P', jnum: int=20, *args,
                 unit: str='cm-1', source: str='calc'):
        self.__center = center
        self.__branch = branch
        self.__jnum = jnum
        self.__unit = unit
        self.__name = str(center) + branch + str(jnum)
        self.__wavelength = eval(source_selection[source])(self.__center,
                                                self.__branch, self.__jnum)
        self.__frequency = convert_rcm_in_hz(self.__wavelength)
        self.unit_check_and_calc()

    @property
    def center(self):
        return self.__center

    @property
    def branch(self):
        return self.__branch

    @property
    def jnum(self):
        return self.__jnum

    @property
    def unit(self):
        return self.__unit

    @unit.setter
    def unit(self, unit: str='mkm'):
        if unit_checker(unit):
            self.__unit = unit
            self.__wavelength = calculate_wavelength([self.__name],
                                            sel_unit=self.__unit)[self.__name]
            return self.__unit
        return None

    @property
    def wavelength(self):
        return self.__wavelength

    @property
    def name(self):
        return self.__name

    @property
    def frequency(self):
        return self.__frequency

    def unit_check_and_calc(self):
        if self.__unit != 'cm-1':
            if unit_checker(self.__unit):
                self.__wavelength = eval(unit_selection[self.__unit])(
                    self.__wavelength)
        pass


class Branch:
    '''Branch is list of Lines'''
    def __init__(self, center: int=10, branch: str='P', jnum_min: int=2,
                 jnum_max: int=58, *args, unit: str='mkm', source: str='calc'):
        self.__center = center
        self.__branch = branch
        self.__jnum_min = jnum_min
        self.__jnum_max = jnum_max
        self.__unit = unit
        self.__source = source
        self.__lines = list()
        self.__get_lines()

    def __get_lines(self):
        for jnum in range(self.__jnum_min, self.__jnum_max + 2, 2):
            self.__lines.append(Line(self.__center, self.__branch, jnum,
                                     unit=self.__unit, source=self.__source))
        pass

    def __iter__(self):
        return iter(self.__lines)

    def __getitem__(self, item):
        return self.__lines[item]


def compare_wavelength_source():
    '''Compare wavelength from Vitterman and from calculation'''
    calc = [Branch(center, branch, 2, 58, unit='cm-1', source='calc') for
            center in center_lines_list for branch in branch_list]
    vitt = [Branch(center, branch, 2, 58, unit='cm-1', source='vitt') for
            center in center_lines_list for branch in branch_list]
    print('Compare wavelength from two sources', '\n', 'Line', '-'*10,
          'df[hz]')
    for br_num, branch in enumerate(calc):
        for line_num, line in enumerate(branch):
            print(line.name, '-> df = ',
                  line.frequency - vitt[br_num][line_num].frequency)
        print('\n', '-'*50, '\n')
    pass


if __name__ == '__main__':
    compare_wavelength_source()
