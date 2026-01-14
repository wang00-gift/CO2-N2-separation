import MDAnalysis as mda

def select_area(file_tpr, file_xtc, atom_name, left, right):
    """
    选择区域
    :param file_tpr:
    :param file_xtc:
    :param left:
    :param right:
    :return:
    """
    u = mda.Universe(file_tpr, file_xtc)
    atom = []
    for se in u.trajectory[-1:]:
        atom = u.select_atoms(f'name {atom_name} and prop z < {right} and prop z > {left}')
    return len(atom)


def main(num_file, atom_name, left, right):
    for i in range(1, num_file+1):
        number = select_area(str(i)+'.tpr', str(i)+'.xtc', atom_name, left, right)
        with open(f"0il-{atom_name}-area-num.txt", "a+") as w:
            w.write(str(number) + '\n')


if __name__ == '__main__':
    num_file = 200
    atom_name = 'C'
    left = 123
    right = 153
    main(num_file, atom_name, left, right)
