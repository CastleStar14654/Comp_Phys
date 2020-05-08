import matplotlib.pyplot as plt
# from matplotlib import rc


def main():
    # problem_1 #
    xs = []
    ys = []
    for name in ('chebyshev', 'hermite'):
        for i in (2, 16, 128, 1024):
            xs.clear()
            ys.clear()
            with open(f'figures/{name}_weighted_{i}') as file:
                for line in file:
                    x, y = map(float, line.split('\t'))
                    xs.append(x)
                    ys.append(y)
            plt.plot(xs, ys)
            plt.tight_layout()
            plt.savefig(f'figures/{name}_weighted_{i}.eps'); plt.clf()
    # end:problem_1 #

    xs.clear()
    ys.clear()
    for name in ('legendre',):
        for i in (16, 64, 256, 1024):
            xs.clear()
            ys.clear()
            with open(f'figures/{name}_zeroes_{i}') as file:
                x = 0
                for line in file:
                    y = float(line.rstrip())
                    xs.append(x)
                    ys.append(y)
                    x += 1
            plt.plot(xs, ys, '.')
            plt.tight_layout()
            plt.savefig(f'figures/{name}_zeroes_{i}.eps'); plt.clf()

    for name in ('legendre',):
        for n, i in zip((16, 64, 256, 1024), (6, 28, 110, 400)):
            xs.clear()
            ys.clear()
            with open(f'figures/{name}_poly_{n}_{i}') as file:
                for line in file:
                    x, y = map(float, line.split('\t'))
                    xs.append(x)
                    ys.append(y)
            plt.plot(xs, ys)
            plt.tight_layout()
            plt.savefig(f'figures/{name}_poly_{n}_{i}.eps'); plt.clf()

if __name__ == '__main__':
    main()
