import os
import pathlib
import shutil
import subprocess

operators = [['ge ', '>='],
             ['le ', '<='],
             ['lt ', '<'],
             ['gt ', '>'],
             ['ne ', '!='],
             ['eq ', '='],
             ['contains '],
             ['datestartswith ']]


def split_filter_part(filter_part):
    for operator_type in operators:
        for operator in operator_type:
            if operator in filter_part:
                name_part, value_part = filter_part.split(operator, 1)
                name = name_part[name_part.find('{') + 1: name_part.rfind('}')]

                value_part = value_part.strip()
                v0 = value_part[0]
                if v0 == value_part[-1] and v0 in ("'", '"', '`'):
                    value = value_part[1: -1].replace('\\' + v0, v0)
                else:
                    try:
                        value = float(value_part)
                    except ValueError:
                        value = value_part

                # word operators need spaces after them in the filter string,
                # but we don't want these later
                return name, operator_type[0].strip(), value

    return [None] * 3


def gen_fchk(chk):
    if not os.path.exists(chk.replace('.chk','.fchk')):
        subprocess.run(["formchk", f'{chk}'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
def gen_cube(fchk,mo):
    filename = f'{fchk}{mo}.cube'
    if os.path.exists(filename):
        return True
    input_ = f'''5
    4
    {mo}
    1 
    2'''
    input_ = bytes(input_, 'utf-8')
    subprocess.run(["Multiwfn", f'{fchk}'], input=input_, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    shutil.move('MOvalue.cub',filename)
    return True

def gen_hole_electron_cube(fchk,value):
    keys = ['Sm index','Sr index','Ghost-hunter index']
    if not value:
        return ''
    filename = f'{fchk}{value}EH.cube'
    if os.path.exists(filename):
        with open(filename + '.txt', "r") as text_file:
            return text_file.read()
    input_ = f'''18
    1
    {fchk.replace('.fchk','.log')}
    {value}
    1 
    2
    1
    13
    '''
    # input_ = bytes(input_, 'utf-8')
    r = subprocess.run(["Multiwfn", f'{fchk}'], input=input_, capture_output=True, text= True)
    out ='\n'
    for i in r.stdout.split('\n'):
        for key in keys:
            if key in i:
                out += i
                out +='\n'
    with open(filename + '.txt', "w") as text_file:
        text_file.write(out)
    shutil.move('transdens.cub',filename)
    return out