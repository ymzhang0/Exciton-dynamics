import h5py
import re
import numpy as np
import os
import time
import math
import spglib
import multiprocessing
from math import pi
import matplotlib.pyplot as plt
"""
Setting and Signing
"""
target_file = "5-epw55-1.out"
#target_file = 'sample'
yambo_set_up = 'r_setup'
path = './' ## the current path for python
elph_file = target_file  ## file name of elphon matrix, to be read

#exciton_file = 'o-2D_WR_WC.exc_qpt1_weights_at_'  ##  file name of exciton, to be read
exciton_file = 'o-2D_WR_WC.exc_qpt1_weights_at_'
#exciton_file = 'o-3D_BSE.exc_qpt1_weights_at_'
#exciton_energy = 'o-GWBSE.exc_qpt1_E_sorted'
exciton_energy = 'o-2D_WR_WC.exc_qpt1_E_sorted'
#exciton_file = 'o-3D_BSE.exc_qpt1_weights_at_'

#f = h5py.File('out_elph.h5', 'w')  ## data for EPW electron-phonon matrix, to be wirtten
iita = 0.0000000000000000001j
save = 0#!!!!Do you need to rewrite h5 file!

#Yes: 1    No: 0

"""
Initial and Test (time, smearing)..
"""

def time_cost(n):
    start = time.perf_counter()
    for i in range(n):
#        get_value_from_elph(8,8,1,1,1)
        get_value_from_elph_qkmnl(8,8,1,1,1)
#    out_epw2h5_qkmnl(elph_file)
    end = time.perf_counter()
    print('time cost: ', end - start)

def smearing_test(smi, smf, interval, T, S, q=2):
    res = open('dataT', 'w')
    res.close()
    Y = []
    x= []
    for enu, i in enumerate(np.linspace(smi, smf, interval)):
        print('\n\n#######################\nThe number of smearing test:', enu + 1,'/' , interval,'\n#######################')
        Y = Y + [sum_k(smearing=i, T=T, S=S ,q=q)]
        x = x + [i]
        res = open('dataT', 'a')
        res.write(('%f %f\n' % (x[-1], Y[-1])))
        res.close()
#    Y = str2int(Y)
#    for i in range(len(Y)):
#        Y[i] = Y[i] * (2 * pi) * 0.001 * 0.001 * e_val / h_bar / (10 ** (15))
#    zz = [x]
#    zz.append(Y)
#    result = np.array(zz)
#    xT = result.transpose()
#    for i in xT:
#        print(i)
#        res.write('%f  %f \n' % (i[0], i[1]))
#    res.close()
    plt.title('Smearing Converence TEST')
    plt.plot(x, Y)
    plt.show()

def initial_option():
    print(
        '\nBBBBBBB      HH       HH\nBB     B     HH       HH\nBB    BB     HH       HH\nBBBBBB       HH H H H HH\nBB    BB     HH       HH\nBB     B     HH       HH\nBBBBBBB      HH       HH  o\n')
    global f
    if os.path.isfile('out_elph.h5'):
        print("out_elph.h5 exists! ")
        if save == 1:
            os.system('rm out_elph.h5')
            print("old out_elph.h5 has been deleted")
            f = h5py.File('out_elph.h5', 'w')
            out_epw2h5(elph_file)
            qkmnl_map(elph_file)
            print('matching kpoints between Yambo and epw:')
            Kmaps()
            print('matching successfully')
        else:
            f = h5py.File('out_elph.h5', 'r')
            qkmnl_map(elph_file)
            print('matching kpoints between Yambo and epw:')
            Kmaps()
            print('matching successfully')
            print('We will use old h5 file')
    else:
        print("NO h5file exists! We will create a new h5 for you")
        f = h5py.File('out_elph.h5', 'w')
        out_epw2h5(elph_file)
        qkmnl_map(elph_file)
        Kmaps()
    t=0.1
    time.sleep(t)

def initial():
    print(
        '\nBBBBBBB      HH       HH\nBB     B     HH       HH\nBB    BB     HH       HH\nBBBBBB       HH H H H HH\nBB    BB     HH       HH\nBB     B     HH       HH\nBBBBBBB      HH       HH  o\n')
    global f
    f = h5py.File('out_elph.h5', 'w')
    try:
        if os.path.isfile('out_elph.h5'):
            print("out_elph.h5 exists! ")
            os.system('rm out_elph.h5')
            print("old out_elph.h5 has been deleted")
        else:
            print("Initialization has been finished\nPython will run in 1s")
    except:
        print("Initialization has been finished\nPython will run in 1s")
    t=0.1
    time.sleep(t)
    out_epw2h5(elph_file)
    qkmnl_map(elph_file)

"""
Reading and Writing
"""

def sign(x):
    if x >= 0:
        return -1
    else:
        return 1

def location(coordinates):
    idx = (mesh[1] * coordinates[1] + mesh[1] * (1 + sign(coordinates[1]))/2)*\
          mesh[0] + mesh[0] * coordinates[0] + mesh[0] * (1 + sign(coordinates[0]))/2

    return np.around(idx).astype(int)

def T_rev_switch(coordinates,switch = 0):
    if switch == 0:
        if ((0<coordinates) & (coordinates<0.5)).all():
            return coordinates
        else:
            translation = np.where(coordinates<0,1,0)
            return coordinates + translation
    if switch == 1:
        if ((0<np.array(coordinates)) & (np.array(coordinates)<0.5)).all():
            return coordinates
        else:
            translation = np.where(coordinates > 0.5, -1, 0)
            return coordinates + translation

def ibzkpts(path):
    ibzkpts_rlu = []
    ibzkpts_cc = []
    i = 0
    num_ibz = 0
    for line in open(path + yambo_set_up,'r'):
        if re.search(r'^\s+K-points', line):
            num_ibz = int(re.findall(r'\d+', line)[0])
            break
    for line in open(path+yambo_set_up,'r'):
        if i < num_ibz :
            if re.search(r'^\sQ\s\[\d+\].*iku',line):
                temp = np.array(re.findall(r'\d+\.?\d+',line)[-4:])
                ibzkpts_rlu.append(temp.astype(float))
                i = i + 1
        else:
            if i < 2*num_ibz:
                if re.search(r'^\sQ\s\[\d+\].*cc',line):
                    temp = np.array(re.findall(r'\d+\.?\d+',line)[-4:])
                    ibzkpts_cc.append(temp.astype(float))
                    i = i + 1
            else:
                break
    return ibzkpts_rlu

#def rbzkpts(path):
#    rbzkpts_cc = []
#    i = 0
#    for line in open(path + elph_file,'r'):
#        if re.search(r'number of k points=',line):
#            num_rbz = int(re.findall(r'\d+', line)[0])
#            break
#    for line in open(path + elph_file,'r'):
#        if i < num_rbz :
#            if re.search(r'^\s+k\(\s+\d+\)', line):
#                    temp = np.array(re.findall(r'\d+\.?\d+', line)[-4:])
#                    rbzkpts_cc.append(temp.astype(float))
#                    i = i + 1
#        else:
#            break
#     kmesh = [60,60,1]
#     num_k = kmesh[0]*kmesh[1]*kmesh[2]
#
#     x = np.linspace(0.0,1.0,kmesh[0]+1)[0:-1]
#     y = np.linspace(0.0,1.0,kmesh[1]+1)[0:-1]
#     z = np.linspace(0.0,1.0,kmesh[2]+1)[0:-1]
#     weight = np.array([1.0 /num_k]*num_k).reshape(kmesh)
#
#     xv,yv,zv = np.meshgrid(x,y,z)
#     rbzkpts_cc = np.stack((xv.T,yv.T,zv.T,weight.T), axis=3).reshape(num_k,4)
#
#     return rbzkpts_cc

def Kmaps(iq=2):
    # global kmaps
    # kmaps = {}
    # ibzkpts_cc = ibzkpts(path)
    # rbzkpts_cc = []
    # for i in rbzkpts_cc:
    #     diff = ibzkpts_cc[ibz] - i
    #     diff_pos = list(map(abs, list(diff[:-1])))
    #     norm.append(np.linalg.norm(diff_pos))
    # return norm.index(min(norm))
    # global kmaps
    # global k_num
    # k_num = '?'
    # for line in open(target_file, 'r'):
    #     if re.search("Using uniform k-mesh:", line):
    #         k_num = int(line.split()[3])*int(line.split()[4])*int(line.split()[5])
    #         break
    #
    # kmaps = {}
    # ibzkpts_rlu = np.array(ibzkpts(path))[:,:-1]
    # rbzkpts_rlu = []
    #
    # for ik in range(k_num+1)[1:]:
    #     rbzkpts_rlu.append(f['q%.0fk%.0f' % (iq, ik)]['q_k_coordinate'][1])
    #
    # rbzkpts_rlu = np.array(rbzkpts_rlu)
    # for ibzkpt in range(len(ibzkpts_rlu)):
    #     norm = []
    #     for rbzkpt in range(len(rbzkpts_rlu)):
    #         diff = ibzkpts_rlu[ibzkpt] - rbzkpts_rlu[rbzkpt]
    #         norm.append(np.linalg.norm(diff))
    #     kmaps[ibzkpt] = norm.index(min(norm)) + 1
    global unfolded_mapping_to_yambo
    global epw_index
    lattice = [[1.0000,0.0000,0.0000],[-0.5000,0.8660,0.0000],[0.0000,0.0000,3.5000]]
    positions = [[0.00000,0.00000,1.75000],[-0.00000,0.57735,1.34996],[0.50000,0.28868,2.15004]]
    numbers = [1, 2, 2]
    cell = (lattice, np.dot(positions, np.linalg.inv(lattice)), numbers)
    epw_grid = []
    for ik in range(1,k_num+1):
        epw_grid.append(f['q%.0fk%.0f' % (iq, ik)]['q_k_coordinate'][1])
    epw_grid = np.array(epw_grid)
    unfolded_index = []
    unfolded_mapping_to_yambo = []
    ibz_grid = np.array(ibzkpts(path))[:,:-1]
    yambo_index = []
    epw_index = []
    mapping, grid = spglib.get_ir_reciprocal_mesh(mesh, cell, is_shift=[0, 0, 0], symprec=1e-3)
    for i in ibz_grid:
        yambo_index.append(location(i))

    irreducible_index = mapping[np.array(yambo_index)]

    for i in irreducible_index:
        unfolded_index = np.append(unfolded_index, np.argwhere(mapping == i)).astype(int)
        unfolded_mapping_to_yambo = np.append(unfolded_mapping_to_yambo, len(np.argwhere(mapping == i))).astype(int)
    unfolded_mapping_to_yambo = np.repeat(np.arange(len(yambo_index)), unfolded_mapping_to_yambo)+1
    unfolded_coordinates = grid[unfolded_index] / mesh

    for i in unfolded_coordinates:
        norm = []
        for j in epw_grid:
            diff = T_rev_switch(i, 0) - j
            norm.append(np.linalg.norm(diff))
        epw_index.append(norm.index(min(norm)))
    epw_index = np.array(epw_index) + 1
    np.savetxt("unfolded_mapping_to_yambo",unfolded_mapping_to_yambo)
    np.savetxt("epw_index",np.array(epw_index))
def out_epw2h5_qkmnl(target_file):
    contain = ['ini']
    q_num='?'
    k_num='?'
    j = 0
    for line in open(target_file, 'r'):
        if re.search("Using uniform q-mesh", line):
            q_num = int(line.split()[3])*int(line.split()[4])*int(line.split()[5])
        if re.search("Using uniform k-mesh:", line):
            k_num = int(line.split()[3])*int(line.split()[4])*int(line.split()[5])
        if re.search("iq =", line):
            iq = int(line.split()[2])
            q = str2int(line.split()[4:7])
        if re.search("ik =", line):
            ik = int(line.split()[2])
            k = str2int(line.split()[4:7])
            if contain[0] == 'ini':
                print('Loop Start')
            else:
                if q_num != '?':
                    if k_num != '?':
                        print("q: %s/%s   k: %.1f" % (iq, q_num, ik*100/k_num), "%" )
                    else:
                        print("q: %s/%s   k: %.1f" % (iq, q_num, ik * 100 ))
                else:
                    if k_num != '?':
                        print("q: %s/?   k: %.1f" % (iq, ik/k_num), "%" )
                    else:
                        print('Warning: We do not know thw q and k number for now')
        if re.match('.*\.\d\d\d\d.*\.\d\d\d\d    [0-9]\.[0-9]+E[+-][0-9]+$',line):
            temp = line.split()
            grp = f.create_group('q%.0fk%.0fm%sn%simode%s' % (iq, ik, temp[0], temp[1], temp[2]))
#           grp.create_dataset('m', data=arrT[0])
#           grp.create_dataset('n', data=arrT[1])
#           grp.create_dataset('imode', data=arrT[2])
            grp.create_dataset('enk', data=float(temp[3]))
            grp.create_dataset('enk+q', data=float(temp[4]))
            grp.create_dataset('omega[meV]', data=float(temp[5]))
            grp.create_dataset('|g|', data=float(temp[6]))
            grp.create_dataset('q_k_coordinate', data=[q, k])
        if j == 0 or j == 1 or j == 2:
            if re.search('------------------------------------------------------------------------------', line):
                j = j + 1
                if j == 1:
                    contain = []
                if j == 2:
                    del contain[0]
                    arr = np.array(contain)
                    arrT = arr.transpose()
#                for i in range(len(arrT[0])):
#                    grp = f.create_group('q%.0fk%.0fm%.0fn%.0fimode%.0f' % (iq, ik, arrT[0][i], arrT[1][i], arrT[2][i]))
#                    grp.create_dataset('m', data=arrT[0])
#                    grp.create_dataset('n', data=arrT[1])
#                    grp.create_dataset('imode', data=arrT[2])
#                    grp.create_dataset('enk', data=arrT[3][i])
#                    grp.create_dataset('enk+q', data=arrT[4][i])
#                    grp.create_dataset('omega[meV]', data=arrT[5][i])
#                    grp.create_dataset('|g|', data=arrT[6][i])
#                    grp.create_dataset('q_k_coordinate', data=[q, k])
#                    grp.create_dataset('min_max_band', data=[min(list(arrT[0])),max(arrT[1])])
            pass
        if j ==1:
            contain.append(str2int(line.split()))
    band_max = max(arrT[0])
    band_min = min(arrT[0])
    mode_num = max(arrT[2])
    f.create_dataset('min_max_nmode_band', data=[band_min, band_max, mode_num])
    print('\nTransfer has been down, out.h5 file has been successfully generated')
# Input: the output file of electron-matrix; Output: out_elph.h5
# !!!! This could not test the q = 0

def out_epw2h5(target_file):
    contain = ['ini']
    global mesh
    global k_num
    q_num = '?'
    k_num = '?'
    j = 0
    for line in open(target_file, 'r'):
        if re.search("Using uniform q-mesh", line):
            q_num = int(line.split()[3])*int(line.split()[4])*int(line.split()[5])
        if re.search("Using uniform k-mesh:", line):
            k_num = int(line.split()[3])*int(line.split()[4])*int(line.split()[5])
            mesh = np.array(line.split()[3:6]).astype(int)
        if re.search("iq =", line):
            iq = int(line.split()[2])
            q = str2int(line.split()[4:7])
        if re.search("ik =", line):
            ik = int(line.split()[2])
            k = str2int(line.split()[4:7])
            if contain[0] == 'ini':
                print('Loop Start')
            else:
                if q_num != '?':
                    if k_num != '?':
                        print("q: %s/%s   k: %.1f" % (iq, q_num, ik * 100 / k_num), "%")
                    else:
                        print("q: %s/%s   k: %.1f" % (iq, q_num, ik ))
                else:
                    if k_num != '?':
                        print("q: %s/?   k: %.1f" % (iq, ik *100 / k_num), "%")
                    else:
                        print('Warning: We do not know thw q and k number for now')
        if re.search('------------------------------------------------------------------------------', line):
            j = j + 1
            if j%2 == 1:
                contain = []
            else:
                del contain[0]
                arr = np.array(contain)
                arrT = arr.transpose()
                grp = f.create_group('q%.0fk%.0f' % (iq, ik ))
                grp.create_dataset('m', data=arrT[0])
                grp.create_dataset('n', data=arrT[1])
                grp.create_dataset('imode', data=arrT[2])
                grp.create_dataset('enk', data=arrT[3])
                grp.create_dataset('enk+q', data=arrT[4])
                grp.create_dataset('omega[meV]', data=arrT[5])
                grp.create_dataset('|g|', data=arrT[6])
                grp.create_dataset('q_k_coordinate', data=[q,k])
            pass
#        if j ==2:
#            global qkmnl_2_num_dic
#           qkmnl_2_num_dic = {}
#           for i in range(len(arrT[0])):
#                qkmnl_2_num_dic['m%.0fn%.0fimode%.0f' % (arrT[0][i], arrT[1][i], arrT[2][i])] = i
        contain.append(str2int(line.split()))
    print('\nTransfer has been down, out.h5 file has been successfully generated')

def qkmnl_map(target_file):
    contain = ['ini']
    q_num = '?'
    k_num  = '?'
    j = 0
    for line in open(target_file, 'r'):
        if re.search("Using uniform q-mesh", line):
            q_num = int(line.split()[3])*int(line.split()[4])*int(line.split()[5])
        if re.search("Using uniform k-mesh:", line):
            k_num = int(line.split()[3])*int(line.split()[4])*int(line.split()[5])
        if re.search("iq =", line):
            iq = int(line.split()[2])
            q = str2int(line.split()[4:7])
        if re.search("ik =", line):
            ik = int(line.split()[2])
            k = str2int(line.split()[4:7])
        if re.search('------------------------------------------------------------------------------', line):
            j = j + 1
            if j%2 == 1:
                contain = []
            else:
                del contain[0]
                arr = np.array(contain)
                arrT = arr.transpose()
            pass
        if j ==2:
            global qkmnl_2_num_dic
            qkmnl_2_num_dic = {}
            for i in range(len(arrT[0])):
                qkmnl_2_num_dic['m%.0fn%.0fimode%.0f' % (arrT[0][i], arrT[1][i], arrT[2][i])] = i
            break
        contain.append(str2int(line.split()))
    print('\nTransfer has been down, out.h5 file has been successfully generated')

def val_band_num(path):
    for line in open(path):
        if re.search(r'\sElectrons\s+',line):
            return int(float(re.findall(r'\d+\.?\d+', line)[0])/2)
            break

def fermi_level(path):
    for line in open(path):
        if re.search('Fermi', line):
            return float(line.split()[5])
            break

def get_exicton_energy(index,path):
    energy = {}
    for line in open(path):
        if not re.search('^#', line):
            energy[int(float(line.split()[-1]))] = [float(line.split()[0]), float(line.split()[1])]
    return energy[index][0], energy[index][1]

def get_value_from_elph_qkmnl(m, n, imode, k, q=2):
    kq_targe = 'q%.0fk%.0fm%.0fn%.0fimode%.0f' % (q, k,m, n, imode)
    [g_h5, omega_h5, enk_h5, enk_q_h5, K_value, Q_value] = [f[kq_targe + '/|g|'][()], f[kq_targe + '/omega[meV]'][()],\
    f[kq_targe + '/enk'][()], f[kq_targe + '/enk+q'][()], f[kq_targe + '/q_k_coordinate'][1], f[kq_targe + '/q_k_coordinate'][0]]
#        else:
#            print('the index of elph is out of range')
    return enk_h5, enk_q_h5, omega_h5, g_h5, K_value, Q_value
# Input: m, n, k, q; Output: |g|
# !!!! This could test the q = 0

def get_value_from_elph(m, n, imode, k, q=2):
    kq_targe = 'q%.0fk%.0f' % (q, k)
    g_h5 = 'Check your process, Eorro: get_value_from_elph'
    ind = qkmnl_2_num_dic['m%sn%simode%s' % (m, n, imode)]
#    for i in range(len(f[kq_targe+'/m'][:])):
#        m_h5 = f[kq_targe + '/m'][i]
#        n_h5 = f[kq_targe + '/n'][i]
#        imode_h5 = f[kq_targe + '/imode'][i]
#        if [m_h5, n_h5, imode_h5] == [m, n, imode]:
    [g_h5, omega_h5, enk_h5, enk_q_h5, K_value, Q_value] = [f[kq_targe + '/|g|'][ind], f[kq_targe + '/omega[meV]'][ind],\
    f[kq_targe + '/enk'][ind], f[kq_targe + '/enk+q'][ind], f[kq_targe + '/q_k_coordinate'][1], f[kq_targe + '/q_k_coordinate'][0]]
#        else:
#            print('the index of elph is out of range')
    return enk_h5, enk_q_h5, omega_h5, g_h5, K_value, Q_value
# Input: m, n, k, q; Output: |g|
# !!!! This could test the q = 0

def get_skip(path):
    z = 0
    for line in open(path + elph_file, 'r'):
        if re.search('excluded bands', line):
            return int(line.split()[-1].split(')')[0])
            z = 1
            break
    if z == 0:
#        print('there is no skip band information')
        print('Warning: there is no skip band value in out file,we will use 0 for calcuaiton')
        return 0

def get_value_from_exction(m, n, S, k, skip):
    A_per_ex = np.loadtxt(path+exciton_file+'%d'%S)
    A_per_ex = np.delete(A_per_ex, [3, 5], axis=1)
    for trans in A_per_ex:
        tran = np.array(trans[0:4])
        tran = tran.astype(int)
        if (tran.astype(int) == [m + skip, n + skip, unfolded_mapping_to_yambo[k], unfolded_mapping_to_yambo[k]]).all():
            return trans[-2]
            break
        else:
            pass
#            print("Warning: unmatched electron-band number!")
    return 0

# Input: m, n, q; Output: AcvS
# !!!! This could test the q = 0

"""
Tools and Maths 
"""
def str2int(list):
    new_list = []
    for i in list:
        try:
            new_list = new_list + [float(i)]
        except:
            pass
    return new_list
# Input: ['num1', 'num2', ...] ; Output: [num1, num2, num3, ...]

def BE(omega, T, mode_num):
    u = 0
    K_B = 8.617343e-05  ### Dimension: eV/K
    return (mode_num / (np.exp((omega * 0.001 - u) / (K_B * T)) - 1))

def FD(omega, u_f, T, mode_num):
    K_B = 8.617343e-05  ### Dimension: eV/K
    return (mode_num / (np.exp((omega * 0.001 - u_f) / (K_B * T)) + 1))

def L_Hospital(m, n, imode, k, T, q_num=3):
    g_set = []
    q_set = []
    omega_set = []
    for ii in range(1, q_num + 1):
        [enk, enk_q, omega, g, K_va, Q_va] = get_value_from_elph(m, n, imode, k, ii)
        q_set.append(vec2len(Q_va))
        g_set.append(g)
        omega_set.append(omega)
    result = np.polyfit(g_set, q_set, 1)[0] / np.polyfit(omega_set, q_set, 1)[0]
    return result

def vec2len(vec=[0,0,0]):
    return math.sqrt(vec[0]**2+vec[1]**2+vec[2]**2)

#def Gauss(x, sigma):
#    return np.exp(-1 * (x ** 2) / (2 * (sigma ** 2))) / (math.sqrt(2 * np.pi) * sigma)

def Gauss(x, sigma):
    if abs(x) < sigma:
        return 1
    else:
        return 0

"""
Sum and Integration for Exciton-Phonon Interaction
"""

def sum_g_cv(smearing, T, m, n, imode, k, skip ,q=2):
    sum_c = 0
    sum_v = 0
    val_bands = val_band_num(path+'r_setup') - skip
    mode_num = int(max(f['q%.0fk%.0f/imode' % (q, epw_index[k])][:]))
    max_band = int(max(f['q%.0fk%.0f/m' % (q, epw_index[k])][:]))
    min_band = int(min(f['q%.0fk%.0f/m' % (q, epw_index[k])][:]))
#    min_band = int(f['min_max_nmode_band'][0])
#    max_band = int(f['min_max_nmode_band'][1])
#    mode_num = int(f['min_max_nmode_band'][2])
#    if imode > 3:
#       print('    SUM for imode: %s' %imode)
    for i in range(val_bands+1, max_band+1):
        [enk_h5, enk_q_h5, omega, g, K_value, Q_value] = get_value_from_elph(n, i, imode, epw_index[k], q)
        sum_c = sum_c + (g**2)*(BE(omega, T, mode_num) * Gauss(enk_h5 - enk_q_h5 + omega*0.001, smearing) + \
                                (BE(omega, T, mode_num)+1) * Gauss(enk_h5 - enk_q_h5 - omega*0.001, smearing))
        #if n==5 and i==5 and imode==6:        
#        print('k k+q omea:', enk_h5, enk_q_h5, omega*0.001 ,'x:',enk_h5 - enk_q_h5 + omega*0.001,'Gauss:', Gauss(enk_h5 - enk_q_h5 + omega*0.001, smearing))
#        print('      sum_c:', sum_c)
    for i in range(min_band, val_bands + 1):
        [enk_h5, enk_q_h5, omega, g, K_value,Q_value] = get_value_from_elph(m, i, imode, epw_index[k], q)
        sum_v = sum_v + (g ** 2) * ((BE(omega, T, mode_num) + 1) * Gauss(enk_h5 - enk_q_h5 + omega*0.001, smearing) + \
                                    BE(omega,T, mode_num) * Gauss(enk_h5 - enk_q_h5 - omega*0.001, smearing))
        #if m==5 and i==5 and imode==6:
#        print('    CV_SUM:  Gauss:', Gauss(enk_h5 - enk_q_h5 + omega*0.001, smearing), '   g:',g,'   BE:',BE(omega,T,mode_num), '   sum_v:', sum_v)
#        print('      sum_v:', sum_v)
#           print('        SUM progress:  %s/%.0f' % (val_bands-i, val_bands - min_band))
#    else:
#       print('    SUM for imode: %s' % imode)
#        for i in range(val_bands+1, max_band+1):
#            [enk_h5, enk_q_h5, omega, g, K_value, Q_value] = get_value_from_elph(m, i, imode, search_kpt(k)+1, q)
#            sum_c = sum_c + 0
#           print('        SUM progress:  %s/%.0f' %(i-val_bands, max_band-val_bands))
#        for i in range(min_band, val_bands + 1):
#            [enk_h5, enk_q_h5, omega, g, K_value, Q_value] = get_value_from_elph(n, i, imode, search_kpt(k)+1, q)
#            sum_v = sum_v + 0
#           print('        SUM progress:  %s/%.0f' % (val_bands - i, val_bands - min_band))
    return sum_c + sum_v

def sum_g_l(smearing, T, m, n, k, skip, imode=0 ,q=2):
    if imode ==0:
        sum_l = 0
        mode_num = int(max(f['q%.0fk%.0f/imode' % (q, epw_index[k])][:]))
    #    mode_num = int(f['min_max_nmode_band'][2])
        for i in range(1, mode_num+1):
            sum_l = sum_l + sum_g_cv(smearing, T, m, n, i, k, skip, q=q)
    #        print('sum_l: ', sum_l)
        return sum_l
    else:
        return sum_g_cv(smearing, T, m, n, imode, k, skip, q=q)


def sum_vc(smearing, T, S, k, skip, imode=0, q=2):
    A_vc = 0
    val_bands = val_band_num(path+'r_setup') - skip
    max_band = int(max(f['q%.0fk%.0f/m' % (q, epw_index[k])][:]))
    min_band = int(min(f['q%.0fk%.0f/m' % (q, epw_index[k])][:]))
#    min_band = int(f['min_max_nmode_band'][0])
#    max_band = int(f['min_max_nmode_band'][1])
#    print('min_band, max_band, val_band:', min_band, max_band, val_bands)
    for i in range(min_band, val_bands + 1):
        for j in range(val_bands+1, max_band+1):
            x = sum_g_l(smearing, T, i, j, k, skip, imode=imode, q=q) * get_value_from_exction(i, j, S, k, skip)
            A_vc = A_vc + x
#            if i ==6 and j ==6:
#    print('  A_vc: ', A_vc, 'sum_g_l:',x)
#            print('the SUM for C: %s/%s' % (i - val_bands +2, max_band - val_bands))
#        print('the SUM for V: %s/%s' %(val_bands - i+1, val_bands - min_band+1))
    return A_vc

def sum_k(smearing,T,S,imode=0,q=2):
    h_bar = 6.62607015 * 10 ** (-34) / (2 * pi)
    e_val = 1.6 * 10 ** (-19)
    skip = get_skip(path)
    print('the skip band is:', skip)
    scattering_rate = 0
    ibzkpts_rlu = ibzkpts(path)
    co = (2 * pi) * 0.001 * 0.001 * e_val / h_bar / (10 ** (15))
    for i in range(len(unfolded_mapping_to_yambo)):
        # weight = int((ibzkpts_rlu[i][-1]) * k_num)
        val = sum_vc(smearing,T,S,i, skip, imode=imode, q=q) * co
        # scattering_rate = scattering_rate + val * weight
        scattering_rate = scattering_rate + val
        print('The sum for K progress: %s/%s' % (i,len(unfolded_mapping_to_yambo)),'   scattering rate:', scattering_rate, '[fs-1]   value for this k: %.f:' % val,'[fs-1]')
#        print('Imaginary Self-Eergy:', ImSE)
#        print('\nBBBBBBB      HH       HH\nBB     B     HH       HH\nBB    BB     HH       HH\nBBBBBB       HH H H H HH\nBB    BB     HH       HH\nBB     B     HH       HH\nBBBBBBB      HH       HH  o\n')

    return scattering_rate

def cal_exph(smearing,T, S_f, n_eig, S_i=1, imode=0, q=2):
    res = open('Result_of_exciton_%s_%s' % (S_i, S_f), 'w')
    res.close()
    Y = []
    x = []
    if imode == -1:
        for i in range(S_i, S_f, round((S_f-S_i)/n_eig)):
            x = []
            temp = ''
            print('\nBBBBBBB      HH       HH\nBB     B     HH       HH\nBB    BB     HH       HH\nBBBBBB       HH H H H HH\nBB    BB     HH       HH\nBB     B     HH       HH\nBBBBBBB      HH       HH  o\n')
            print('\n\n#######################\nThe number of exciton progress:',
                  int((i - S_i) / (round((S_f - S_i) / n_eig)) + 1), '/', n_eig,
                  '\n#######################')
            mode_num = int(max(f['q1k1/imode'][:]))
            try:
                x = x + [get_exicton_energy(i, exciton_energy)[0]]
                print('we find energy information')
            except:
                x = x + [' ']
                print('Warning: we do not fine exciton energy information')
            for lam in range(1, mode_num+1):
                print('\nthe number of exciton', int((i - S_i) / (round((S_f - S_i) / n_eig)) + 1), '/', n_eig,'   imode: %s' % lam)
                x = x + [sum_k(smearing=smearing, T=T, S=i, imode=lam, q=q)]
            Y.append(x)
            res = open('Result_of_exciton_%s_%s' % (S_i, S_f), 'a')
            for el in Y[-1]:
                temp = temp + " " +str(el)
            res.write(temp+'\n')
            res.close()
#            Y.append(x)
#        for ele in Y:
#            temp = ''
#            for el in ele:
#                temp = temp + str(el)
#            res.write(temp+'\n')
#        res.close()
    else:
        res = open('Result_of_exciton_%s_%s' % (S_i, S_f), 'w')
        for i in range(S_i, S_f, round((S_f-S_i)/n_eig)):
            print('\nBBBBBBB      HH       HH\nBB     B     HH       HH\nBB    BB     HH       HH\nBBBBBB       HH H H H HH\nBB    BB     HH       HH\nBB     B     HH       HH\nBBBBBBB      HH       HH  o\n')
            print('\n\n#######################\nThe number of exciton progress:',
                  int((i - S_i) / (round((S_f - S_i) / n_eig)) + 1), '/', n_eig,
                  '\n#######################')
            Y = Y + [sum_k(smearing=smearing, T=T, S=i, imode=imode, q=q)]
            try:
                [ex_energy, ex_weight] = get_exicton_energy(i, exciton_energy)
                x = x + [ex_energy]
                print('We find exciton energy')
            except:
                x = x + [i]
                print('there is no exciton energy found')
        Y =str2int(Y)
        zz = [x]
        zz.append(Y)
        result = np.array(zz)
        xT = result.transpose()
        for i in xT:
            print(i)
            res.write('%f  %f \n' % (i[0], i[1]))
        res.close()
# if Imode=0: it calculates all phonon modes, and this is default set
# if Imode=[1-9]+: it calculates the #imode
# if imode=-1: it calculates all  imode respectively
# start with S_i, end but not include S_f, n_eig means the number of exciton you want to calculate
"""
Sum and Integration for high-order of electron-Phonon Interaction
"""
def q_search(k_kp, q_list):
    norm = []
    for i in q_list:
        dif = np.array(k_kp) - i
        norm.append(vec2len(dif))
    return norm.index(min(norm))+1
# input the k_kp and q_list set
# output nearby q index (real)

def q_k_list():
    for line in open(target_file, 'r'):
        if re.search("Using uniform q-mesh", line):
            q_num = int(line.split()[3])*int(line.split()[4])
        if re.search("Using uniform k-mesh:", line):
            k_num = int(line.split()[3])*int(line.split()[4])
    q_num = 1
    q_list = []
    k_list = []
    for i in range(1,q_num+1):
        q_list.append(f['q%.0fk1/q_k_coordinate' % i][0])
    for j in range(1,k_num+1):
        k_list.append(f['q1k%.0f/q_k_coordinate' % j][1])
    return q_list, k_list

def sum_pi_ss(T, fermi, m, n, k, imode_p, k_p,omega,q_p, mode_num):
    pi_ss = 0
    e_u_k = get_value_from_elph(m, m, imode_p, k + 1, q=1)[0]
    e_up_kp = get_value_from_elph(n, n, imode_p, k_p + 1, q=1)[0]
    omega_qp_imodep = get_value_from_elph(m, n, imode_p, k + 1, q=q_p)[2]
    for s in [-1, 1]:
        for s_p in [-1, 1]:
            pi_ss = pi_ss + ((FD(e_u_k,fermi,T,1)-FD(e_up_kp-s*s_p*omega_qp_imodep,fermi,T,1))/(e_u_k-(e_up_kp-s*s_p*omega_qp_imodep)))*((s*(BE(s*omega_qp_imodep,T, mode_num)+FD(s_p*e_up_kp,fermi,T,1)))/(omega*(omega+iita+s_p*(e_u_k-e_up_kp+s*omega_qp_imodep))))

def sum_pi_2(T, imode, omega, segment=1, inde=1):
    rbzkpts_cc = rbzkpts(path)
    fermi = fermi_level('./'+elph_file)
    pi_2 = 0
    max_band = int(max(f['q%.0fk%.0f/m' % (1, 1)][:]))
    min_band = int(min(f['q%.0fk%.0f/m' % (1, 1)][:]))
    mode_num = int(max(f['q%.0fk%.0f/imode' % (1, 1)][:]))
#    min_band = int(f['min_max_nmode_band'][0])
#    max_band = int(f['min_max_nmode_band'][1])
#    mode_num = int(f['min_max_nmode_band'][2])
    [q_list, k_list] = q_k_list()
    period = len(rbzkpts_cc)//segment
    if inde != segment:
        k_i = (inde-1)*period
        k_f = inde*period
    else:
        k_i = (inde-1)*period
        k_f = len(rbzkpts_cc)
    for i_k in range(k_i, k_f):
        print('the sum of k:',i_k - k_i,'/',period,'segment:', inde)
        for i_u in range(min_band, max_band + 1):
#            print('the sum of u:', i_u-min_band, '/', max_band+1-min_band)
            for i_u_p in range(min_band, max_band + 1):
                print('    the sum of u_p:', i_u_p-min_band, '/', max_band+1-min_band)
                for i_imode_p in range(1, mode_num + 1):
                    print('        the sum of imode_p:', i_imode_p, '/', mode_num)
                    for i_k_p in range(len(rbzkpts_cc)):
                        g_uu_imode_k_0 = get_value_from_elph(i_u, i_u, imode, i_k+1, q=1)[3]
                        g_upup_imode_kp_0 = get_value_from_elph(i_u_p, i_u_p, imode, i_k_p+1, q=1)[3]
                        diff = k_list[i_k] - k_list[i_k_p]
                        q_p = q_search(diff, q_list)
                        if q_p == 1:
                            continue
                        else:
                            g_uup_imodep_k_qp = get_value_from_elph(i_u, i_u_p, i_imode_p, i_k+1, q=q_p)[3]
                            co = (g_uu_imode_k_0 ** 2) * (1 - (g_upup_imode_kp_0/g_uu_imode_k_0)) * (g_uup_imodep_k_qp ** 2)
                            pi_2 = pi_2 + (sum_pi_ss(T, fermi, i_u, i_u_p, i_k, i_imode_p, i_k_p, omega, q_p, mode_num) * co)
    print('the pi-2 has been successfully calculated with T: %s, imode: %s, omega: %s' % (T, imode, omega))
    print('pi_2:', pi_2)
    return pi_2

"""
Multi_process
"""

if __name__ == '__main__':
    pass
    initial_option()
#    smearing_test(0.001,0.1, 100, 300, 3)
    cal_exph(smearing=0.001, T=300, S_f=16, n_eig=10, S_i=6, imode=-1)
#    sum_k(smearing=0.001, T=300, S=3, imode=0)
#    sum_pi_2(T=10, imode=4, omega=0.001)





