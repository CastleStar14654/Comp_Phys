import matplotlib.pyplot as plt
import struct

def main():
    ys = [1, -1]
    xs = [0, 0]
    x = 0
    size = struct.calcsize('d')
    # for seed in (14654, 35641, 6244212, 19780503):
    #     with open(f'output/ising_2.1_{seed}.out', 'br') as file:
    #         while bt := file.read(size):
    #             y = struct.unpack('d', bt)[0]
    #             ys.append(y)
    #             xs.append(x)
    #             x += 32**2
    #     plt.plot(xs, ys)
    #     plt.tight_layout(1.08,None,None,(0.04,0,1,1))
    #     plt.ylabel(r'$\bar{s}$')
    #     plt.xlabel(r'iterate count $n$')
    #     # plt.savefig(f'figures/ising_2.1_{seed}.eps'); plt.clf()
    #     plt.show()
    #     ys.clear(); ys.extend((-1, 1))
    #     xs.clear(); xs.extend((0, 0))
    #     x = 0

    # =================================================================

    d = {}
    key_lst_label = {
        'energy':'$U$',
        'C_H':'$C_H$',
        'M':'$M$',
        # 'M_abs':r'|$M$|',
        'chi':r'$\chi$',
        'chi_abs':r'|$\chi$|'
    }
    key_lst = ['#', 'beta', 'energy', 'C_H', 'M', 'M_abs', 'chi', 'chi_abs']
    # with open('2_Ising_model_2.2.txt', encoding='UTF-16 LE') as file:
    #     file.readline(); file.readline(); file.readline()
    #     d[32] = ([{key: [] for key in key_lst[2:]} for i in range(4)], [])
    #     lst = d[32][0]
    #     betas = d[32][1]
    #     while line := file.readline():
    #         line = line.split()
    #         number = int(line[0])
    #         if number == 0:
    #             betas.append(float(line[1]))
    #         for i in range(2, len(key_lst)):
    #             lst[number][key_lst[i]].append(float(line[i]))

    # for key in key_lst_label:
    #     for network in d[32][0]:
    #         plt.plot(d[32][1], network[key], 'C0.')
    #     plt.ylabel(key_lst_label[key])
    #     plt.xlabel(r'$\beta J$')
    #     plt.tight_layout(1.08,None,None,(0.02,0,1,1))
    #     plt.savefig(f'figures/ising_2.2_{key}.eps'); plt.clf()


    for L in (16, 24, 32, 40, 48, 56, 64):
        with open(f'2_Ising_model_2.3_{L}.txt', encoding='UTF-16 LE') as file:
            file.readline(); file.readline(); file.readline(); file.readline()
            d[L] = ([{key: [] for key in key_lst[2:]} for i in range(4)], [])
            lst = d[L][0]
            betas = d[L][1]
            while line := file.readline():
                line = line.split()
                number = int(line[0])
                if number == 0:
                    betas.append(float(line[1]))
                for i in range(2, len(key_lst)):
                    lst[number][key_lst[i]].append(float(line[i]))

    # for key in key_lst_label:
    #     for L in (16, 24, 32, 40, 48, 56, 64):
    #         plt.plot(ave_d[L][1], ave_d[L][0][key], label=L)
    #     plt.legend()
    #     plt.ylabel(key_lst_label[key])
    #     plt.xlabel(r'$\beta J$')
    #     plt.tight_layout(1.08,None,None,(0.02,0,1,1))
    #     plt.show()

    for key in ('chi', 'chi_abs'):
        for L in (16, 24, 32, 40, 48, 56, 64):
            for network in d[L][0]:
                plt.plot(d[L][1], network[key], 'C0.')
            plt.ylabel(key_lst_label[key])
            plt.xlabel(r'$\beta J$')
            plt.tight_layout(1.08,None,None,(0.02,0,1,1))
            # plt.show()
            plt.savefig(f'figures/ising_2.3_{L}_{key}.eps'); plt.clf()




    # for key in key_lst_label:
    #     count = 0
    #     for L in (16, 24, 32, 40, 48, 56, 64):
    #         for network in d[L][0]:
    #             a = plt.plot(d[L][1], network[key], f'C{count}.')
    #         # a[0].set_label(L)
    #         count += 1
    #         plt.title(L)
    #         plt.ylabel(key_lst_label[key])
    #         plt.xlabel(r'$\beta J$')
    #         plt.tight_layout(1.08,None,None,(0.02,0,1,1))
    #         plt.show()

    # for key in key_lst_label:
    #     count = 0
    #     for L in (16, 24, 32, 40, 48, 56, 64):
    #         for network in d[L][0]:
    #             a = plt.plot(d[L][1], network[key], f'C{count}.')
    #         a[0].set_label(L)
    #         count += 1
    #     plt.legend()
    #     plt.ylabel(key_lst_label[key])
    #     plt.xlabel(r'$\beta J$')
    #     plt.tight_layout(1.08,None,None,(0.02,0,1,1))
    #     plt.show()

    # for key in key_lst_label:
    #     count = 0
    #     for L in (32, 40, 48):
    #         for network in d[L][0]:
    #             a = plt.plot(d[L][1], network[key], f'C{count}.')
    #         a[0].set_label(L)
    #         count += 1
    #     plt.legend()
    #     plt.ylabel(key_lst_label[key])
    #     plt.xlabel(r'$\beta J$')
    #     plt.tight_layout(1.08,None,None,(0.02,0,1,1))
    #     plt.show()

    # for key in key_lst_label:
    #     count = 0
    #     for L in (48, 56, 64):
    #         for network in d[L][0]:
    #             a = plt.plot(d[L][1], network[key], f'C{count}.')
    #         a[0].set_label(L)
    #         count += 1
    #     plt.legend()
    #     plt.ylabel(key_lst_label[key])
    #     plt.xlabel(r'$\beta J$')
    #     plt.tight_layout(1.08,None,None,(0.02,0,1,1))
    #     plt.show()

    # =================================================================

    # for L in (16, 24, 40, 48, 56, 64):
    #     for j in range(4):
    #         with open(f'output/ising_2.3_{L}_{j}.out', 'br') as file:
    #             while bt := file.read(size):
    #                 y = struct.unpack('d', bt)[0]
    #                 ys.append(y)
    #                 xs.append(x)
    #                 x += L**2
    #         plt.plot(xs, ys)
    #         plt.tight_layout(1.08,None,None,(0.04,0,1,0.96))
    #         plt.title(f'{L} ' + str(j))
    #         # plt.savefig(f'figures/{name}_weighted_{i}.eps'); plt.clf()
    #         plt.show()
    #         ys.clear(); ys.extend((-1, 1))
    #         xs.clear(); xs.extend((0, 0))
    #         x = 0

    # =================================================
    '''
    hs = []
    ss = []
    for N in (512, 1024, 2048):
        for M in (512, 1024, 2048):
            with open(f'output/ising_2.5_{M}_{N}.out', 'br') as file:
                while bth := file.read(size):
                    bts = file.read(size)
                    h = struct.unpack('d', bth)[0]
                    s = struct.unpack('d', bts)[0]
                    hs.append(h)
                    ss.append(s)
            plt.plot(hs, ss, '.', label=f"{M=}")
            hs.clear()
            ss.clear()
        plt.tight_layout(1.08,None,None,(0.03,0.03,1,1))
        # plt.title(f'{M=} {N=}')
        plt.legend()
        plt.ylabel(r'$\bar{s}$')
        plt.xlabel(r'$H$')
        plt.savefig(f'figures/ising_2.5_N_{N}.eps'); plt.clf()
        # plt.show()

    for M, N in ((512, 2048), (1024, 1024), (2048, 512)):
        with open(f'output/ising_2.5_{M}_{N}.out', 'br') as file:
            while bth := file.read(size):
                bts = file.read(size)
                h = struct.unpack('d', bth)[0]
                s = struct.unpack('d', bts)[0]
                hs.append(h)
                ss.append(s)
        plt.plot(hs, ss, '.', label=f"{M=}, {N=}")
        hs.clear()
        ss.clear()
    plt.tight_layout(1.08,None,None,(0.03,0.03,1,1))
    # plt.title(f'{M=} {N=}')
    plt.legend()
    plt.ylabel(r'$\bar{s}$')
    plt.xlabel(r'$H$')
    # plt.savefig(f'figures/ising_2.5_MN.eps'); plt.clf()
    plt.show()

    for M, N in ((512, 1024), (1024, 512)):
        with open(f'output/ising_2.5_{M}_{N}.out', 'br') as file:
            while bth := file.read(size):
                bts = file.read(size)
                h = struct.unpack('d', bth)[0]
                s = struct.unpack('d', bts)[0]
                hs.append(h)
                ss.append(s)
        plt.plot(hs, ss, '.', label=f"{M=}, {N=}")
        hs.clear()
        ss.clear()
    plt.tight_layout(1.08,None,None,(0.03,0.03,1,1))
    # plt.title(f'{M=} {N=}')
    plt.legend()
    plt.ylabel(r'$\bar{s}$')
    plt.xlabel(r'$H$')
    # plt.savefig(f'figures/ising_2.5_MN_2small.eps'); plt.clf()
    plt.show()

    for M, N in ((1024, 2048), (2048, 1024)):
        with open(f'output/ising_2.5_{M}_{N}.out', 'br') as file:
            while bth := file.read(size):
                bts = file.read(size)
                h = struct.unpack('d', bth)[0]
                s = struct.unpack('d', bts)[0]
                hs.append(h)
                ss.append(s)
        plt.plot(hs, ss, '.', label=f"{M=}, {N=}")
        hs.clear()
        ss.clear()
    plt.tight_layout(1.08,None,None,(0.03,0.03,1,1))
    # plt.title(f'{M=} {N=}')
    plt.legend()
    plt.ylabel(r'$\bar{s}$')
    plt.xlabel(r'$H$')
    # plt.savefig(f'figures/ising_2.5_MN_2big.eps'); plt.clf()
    plt.show()
    '''


if __name__ == '__main__':
    main()
