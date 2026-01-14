import os
import time
import numpy as np


def Mdrun(file_mdp, file_gro, file_top, file_out):
    os.system('gmx grompp -f ' + file_mdp + ' -c ' + file_gro+ ' -r ' + file_gro +  ' -p ' + file_top + ' -o ' + file_out + '.tpr')
    os.system('gmx mdrun -v -deffnm ' + file_out +' -ntomp 18' +' -ntmpi 1')     


def DelGMXGro(file_gro, atom1_name, atom2_name, local_zr, local_zl):
    '''
    Read single frame .gro file,
    :param file_gro:
    :return: box, number of atom1 and atom2, moles. moles: list consisted of Molecule class.
    '''
    with open(file_gro, 'r') as f:
        lines = [line.rstrip() for line in f.readlines()]
    totalmoles = []
    openedfile = open(file_gro,'r')
    natoms = len(openedfile.readlines()) - 3
    atoms = lines[2: 3 + natoms]
    n = 0
    atom1_lnum = 0
    atom2_lnum = 0
    atom1_rnum = 0
    atom2_rnum = 0
    atom1_mnum = 0
    atom2_mnum = 0
    while n < natoms:
        typename = atoms[n][5:10].strip()  
        coordinate = list(map(float, atoms[n][20:44].split()))  
        space = list(map(float, atoms[n][44:69].split()))
        symbol = atoms[n][10:15].strip()
        coordinates = []
        spaces = []
        mole = []
        if symbol != atom1_name and symbol != atom2_name or coordinate[2] <= local_zr:
            if symbol == atom1_name and   local_zl < coordinate[2] <= local_zr:
                atom1_mnum += 1
            elif symbol == atom2_name and   local_zl < coordinate[2] <= local_zr:
                atom2_mnum += 1
            coordinates.append(coordinate)
            atomlist = []
            x = coordinate[0]
            y = coordinate[1]
            z = coordinate[2]
            sx = space[0]
            sy = space[1]
            sz = space[2]
            atomlist.append(typename)
            atomlist.append((atoms[n][10:15]).strip())
            atomlist.append(x)
            atomlist.append(y)
            atomlist.append(z)
            atomlist.append(sx)
            atomlist.append(sy)
            atomlist.append(sz)
            mole.append(atomlist)
            n += 1
            totalmoles.append(mole)
            if symbol == atom2_name and coordinate[2] <= local_zl:
                atom2_lnum += 1
            elif symbol == atom1_name and coordinate[2] <= local_zl:
                atom1_lnum += 1
        elif symbol == atom1_name and coordinate[2] > local_zr:
            atom1_rnum += 1
            n += 5
        elif symbol == atom2_name and coordinate[2] > local_zr:
            atom2_rnum += 1
            n += 3
        else:
            n += 3
    boxv = np.array(list(map(float, lines[-1].split())))
    natom = len(totalmoles)
    f.close()
    return natom, boxv, totalmoles, atom1_lnum, atom2_lnum, atom1_rnum, atom2_rnum, atom1_mnum, atom2_mnum


def PackmolGro(file_inp, atom1_inum, atom2_inum):
    pack_infos = []
    with open(file_inp, 'r') as f:
        lines2 = [line.rsplit() for line in f.readlines()]
        lines2[7][1] = str(atom1_inum)
        lines2[11][1] = str(atom2_inum)
        if atom1_inum > 0 and atom2_inum > 0:
            for lines in lines2:
                pack_infos.append(lines)
        elif atom1_inum > 0 and atom2_inum <= 0:
            for lines in lines2[:10]:
                pack_infos.append(lines)
        elif atom1_inum <= 0 and atom2_inum > 0:
            for lines in lines2[0:6]:
                pack_infos.append(lines)
            for lines in lines2[10:]:
                pack_infos.append(lines)
        elif atom1_inum <= 0 and atom2_inum <= 0:
            for lines in lines2[0:6]:
                pack_infos.append(lines)
    return pack_infos


def Top(file_top, atom1_tnum, atom2_tnum):
    with open(file_top, 'r') as f:
        info_top = [line.rsplit() for line in f.readlines()]
        info_top[-2][-1] = str(atom1_tnum)
        info_top[-1][-1] = str(atom2_tnum)
    return info_top


def Newgro(file_gro1, file_gro2):
    with open(file_gro1, 'r') as f1:
        with open(file_gro2, 'r') as f2:
            infos = []
            new_infos = []
            n = 2
            lines1 = [lines.rstrip() for lines in f1.readlines()]
            lines2 = [lines.rstrip() for lines in f2.readlines()]
            natom1 = len(lines1)
            natom2 = len(lines2) - 1
            while n < natom2:
                lines1[n] = lines2[n]
                n += 1
            lines1[-1] = lines2[-1]
            atoms = lines1[2:natom1 - 1]
            infos.append(lines1[0])
            infos.append(lines1[1])
            print(infos)
            j = 0
            while j < natom1 - 3:
                typename = atoms[j][5:8].strip()
                if typename == 'COF' or typename == 'BR' or typename == 'NSC' or typename == 'NSA':                 
                # 膜的残基名
                    infos.append(atoms[j])
                    j += 1
                elif typename == 'CO2':          
                    infos.append(atoms[j])
                    infos.append(atoms[j + 1])
                    infos.append(atoms[j + 2])
                    infos.append(atoms[j + 3])
                    infos.append(atoms[j + 4])
                    j += 5
                else:
                    j += 1
            j = 0
            while j < natom1 - 3:
                typename = atoms[j][5:8].strip()
                if typename == 'N2':             
                    infos.append(atoms[j])
                    infos.append(atoms[j + 1])
                    infos.append(atoms[j + 2])
                    j += 3
                else:
                    j += 1
            k = 2
            new_infos.append(lines1[0])
            new_infos.append(lines1[1])
            while k < natom1 - 1:
                aa = k - 1
                bb = ('%5d' % aa)
                new_info = infos[k][:15] + str(bb) + infos[k][20:]
                new_infos.append(new_info)
                k += 1
            new_infos.append(lines2[-1])
    return new_infos


def main( file_gro, file_inp, file_mdp, file_top,
          local_zr, local_zl, atom1_name, atom2_name, gasnum):
    ''' ------------------------------------------------------------------------------------- '''
    natom, boxv, totalmoles, atom1_lnum, atom2_lnum, atom1_rnum, atom2_rnum, atom1_mnum, atom2_mnum = \
        DelGMXGro(file_gro, atom1_name, atom2_name, local_zr, local_zl)
    gro_file = open("del.gro", 'w+')
    now_time = time.strftime("%a %b %d %H:%M:%S %Y", time.localtime())
    gro_file.write("%s\n" % ('Generated on ' + str(now_time)))
    gro_file.write("%5d\n" % natom)
    n = 1
    for info in totalmoles:
        for info1 in info:
            resid = n
            resname = info1[0]
            atomtype = info1[1]
            x = info1[2]
            y = info1[3]
            z = info1[4]
            sx = info1[5]
            sy = info1[6]
            sz = info1[7]
            gro_file.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n" %
                               (resid, resname, atomtype, resid, x, y, z, sx, sy, sz))
        n += 1
    gro_file.write("%10.5f%10.5f%10.5f\n" % (boxv[0], boxv[1], boxv[2]))
    gro_file.close()
    os.system('gmx editconf -f del.gro -o del.pdb')

    # print(atom1_rnum, atom2_rnum)

    atom1_inum = gasnum - atom1_lnum                
    atom2_inum = gasnum - atom2_lnum
    lines2 = PackmolGro(file_inp, atom1_inum, atom2_inum)
    space = ' '
    with open('md.inp', 'a') as w:
        for line2 in lines2:
            info = space.join(line2)
            w.write(info + '\n')                    
    w.close()

    os.system('packmol < md.inp')
    os.system('gmx editconf -f md.pdb -o md.gro')
    new_infos = Newgro('md.gro', 'del.gro')
    atom1_tnum = 0
    atom2_tnum = 0
    with open('mdrun.gro', 'w+') as w1:
        for new_info in new_infos:
            if new_info[12:15].strip() == atom1_name:
                atom1_tnum += 1
            elif new_info[12:15].strip() == atom2_name:
                atom2_tnum += 1
            else:
                pass
            w1.write(new_info + '\n')

    info_top = Top(file_top, atom1_tnum, atom2_tnum)
    top_file = open('md.top', 'a')
    for infos_tops in info_top:
        top_file.write(space.join(infos_tops) + '\n')
    top_file.close()                              

    Mdrun(file_mdp, 'mdrun.gro', 'md.top', '1')
    with open('gas_rum.txt', 'w+') as w2:

        a = str(atom1_rnum)
        b = str(atom2_rnum)
        w2.write(a + '     ' + b + '\n')


if __name__ == '__main__':
    file_gro = 'start.gro'                                    
    file_inp = 'start.inp'
    file_mdp = 'md.mdp'
    file_top = 'topol.top'
    local_zr = 16.7                                    
    local_zl = 10.65                                      
    atom1_name = 'C'                       
    atom2_name = 'N11'
    gasnum = 100                                         
    main(file_gro, file_inp, file_mdp, file_top,
         local_zr, local_zl, atom1_name, atom2_name, gasnum)
