import molconst12c16o2

VITTERMANN_FILE_NAME = 'CO2laserlinesVitterman.dat'
LIGHT_VELOCITY = 299792458
branch_list = ['P', 'R']
center_lines_list = ['10', '9']
mol_const_dict = {
                  '10': tuple(molconst12c16o2.CO210MKM),
                  '9': tuple(molconst12c16o2.CO29MKM)
                 }

jnum_const_dict = {
                    'P': 'const_branch_P',
                    'R': 'const_branch_R'
                  }
unit_selection = {
                  'cm-1': '(lambda x: x)',
                  'mkm': 'convert_wavelength_cm_in_mkm',
                  'hz': 'convert_wavelength_cm_in_hz'
                  }


def const_branch_P(j_num: int) -> int:
    return -j_num


def const_branch_R(j_num: int) -> int:
    return j_num + 1


def get_line_jnum(line_name: str) -> int:
    jnum = line_name.split('_')[-1]
    return int(jnum)


def get_line_branch(line_name: str) -> str:
    branch = line_name.split('_')[1]
    return branch


def get_line_center(line_name: str) -> str:
    center = line_name.split('_')[0]
    return center


def get_line_name(center: int, branch: str, jnum: int) -> str:
    line_name = str(center) + '_' + branch + '_' + str(jnum)
    return line_name


def convert_wavelength_cm_in_mkm(wavelength_in_cm: float) -> float:
    wavelength_in_mkm = 10000.0/float(wavelength_in_cm)
    return wavelength_in_mkm


def convert_wavelength_cm_in_hz(wavelength_in_cm: float) -> float:
    wavelength_in_hz = LIGHT_VELOCITY*wavelength_in_cm*100
    return wavelength_in_hz


def get_CO2_emission_lines_one_branch(
        center: int = 10, branch: str = 'P',
        jnum_min: int = 2, jnum_max: int = 58) -> list:
    assert jnum_min % 2 == 0 and (
        jnum_max % 2 == 0 and jnum_min >= 2 and jnum_max <= 58)
    CO2_emission_lines = [get_line_name(center, branch, jnum) for
                          jnum in range(jnum_min, jnum_max + 2, 2)]
    return CO2_emission_lines


def get_CO2_emission_lines(jnum_min: int = 2, jnum_max: int = 58) -> list:
    CO2emission_lines = list()
    for center in center_lines_list:
        for branch in branch_list:
            CO2emission_lines += get_CO2_emission_lines_one_branch(center,
                                                    branch, jnum_min, jnum_max)
    return CO2emission_lines


def get_lines_Vittermann(file_name: str = VITTERMANN_FILE_NAME) -> list:
    wavelength_center_lines = list()
    with open(file_name, 'r') as lines_Vittermann:
        for line in lines_Vittermann:
            wavelength_center_lines.append(list(
                line.replace(',', '.').split()))
    converted_lines = list()
    for lines in zip(*wavelength_center_lines):
        converted_lines += lines
    return converted_lines


def get_lines_Vitt_dict(wavelength_convert_func = lambda x: float(x)) -> dict:
    emission_lines = get_CO2_emission_lines()
    CO2laser_line_dict = dict()
    wavelength_center_lines = get_lines_Vittermann()
    for number, line in enumerate(emission_lines):
        CO2laser_line_dict[line] = wavelength_convert_func(
        wavelength_center_lines[number])
    return CO2laser_line_dict


def get_branch_const_dict(emission_lines: list) -> dict:
    branch_const_dict = dict()
    for line in emission_lines:
        center = get_line_center(line)
        branch = get_line_branch(line)
        jnum = get_line_jnum(line)
        constants = [eval(jnum_const_dict[branch])(jnum),
                     *mol_const_dict[center]]
        branch_const_dict[line] = tuple(constants)
    return branch_const_dict


def wavelength_calculator(line_const: tuple) -> float:
    wavelength = 0
    line_const_list = list(line_const)
    jnum = line_const[0]
    line_const_list.pop(0)
    for number, constanta in enumerate(line_const_list):
        wavelength = wavelength + constanta*jnum**number
    return wavelength


def calculate_wavelength(emission_lines: list, *args, sel_unit = 'cm-1') -> dict:
    branch_const_dict = get_branch_const_dict(emission_lines)
    lines_wavelength_calc_dict = dict()
    for line in emission_lines:
        lines_wavelength_calc_dict[line] = eval(unit_selection[sel_unit])(
            wavelength_calculator(branch_const_dict[line]))
    return lines_wavelength_calc_dict


def compare_all_lines_wavelength(jnum_min: int = 2, jnum_max: int = 58):
    CO2laser_lines_Vitt_dict = get_lines_Vitt_dict()
    emission_lines = get_CO2_emission_lines(jnum_min, jnum_max)
    CO2laser_lines_calc = calculate_wavelength(emission_lines)
    line = 'line'
    wavelength = 'WLVit'
    wcalc = 'WLcalc'
    wdiff = 'difference'
    print(f'{line:<8} -> {wavelength:<10} -> {wcalc:<10} -> {wdiff:<10}')
    for line, wavelength in CO2laser_lines_calc.items():
        wVitt = CO2laser_lines_Vitt_dict[line]
        wcalc = wavelength
        wdiff = wVitt - wcalc
        print(f'{line:<8} -> {wVitt:<10} -> {wcalc:<10} -> {wdiff:<10}')
    pass


class CO2lines:
    def __init__(self, *args, center = 10,  branch = 'P',
                 jnum_min = 2, jnum_max = 58, unit = 'cm-1', **kwards):
        self.unit = unit
        self.emission_lines = get_CO2_emission_lines_one_branch(
            center, branch, jnum_min, jnum_max)
        self.CO2lines_dict = calculate_wavelength(self.emission_lines,
                                                  sel_unit = self.unit)

    def __getitem__(self, line):
        return self.CO2lines_dict[line]

    def __iter__(self):
        return iter(self.CO2lines_dict)

    def printlines(self):
        for line, wavelength in self.CO2lines_dict.items():
            line = line.replace(r'_', '')
            print(f' {line}  ->  {wavelength:<12} {self.unit}')
        pass


if __name__ == '__main__':
    p10 = CO2lines(center = 10, branch = 'P', jnum_min = 12, jnum_max = 12,
                   unit = 'cm-1')
    p10.printlines()
    print(p10['10_P_12'])
    for item in p10:
        print(item, '    ', p10[item])
    compare_all_lines_wavelength()
