import random
import matplotlib.pyplot as plt
from matplotlib import rc


def main():
    # gx #
    def g(x):
        return x*(1-x)
    # end:gx #

    # problem_1 #
    r = 0.5
    for x_0 in (random.random() for i in range(10)):
        plt.plot(list(logistic_generator(g, r, x_0, 11)), '.--')
    plt.xlabel('$n$'); plt.ylabel('$x^*$')
    plt.tight_layout()
    plt.savefig(f'p1_{r}.eps'); plt.clf()
    r = 1.5
    for x_0 in (random.random() for i in range(10)):
        plt.plot(list(logistic_generator(g, r, x_0, 15)), '.--')
    plt.xlabel('$n$'); plt.ylabel('$x^*$')
    plt.tight_layout()
    plt.savefig(f'p1_{r}.eps'); plt.clf()
    # end:problem_1 #

    # problem_2 #
    N = 50
    rs = [i/N for i in range(3*N)]
    xs = []
    for r in rs:
        x_0 = random.random()
        it = logistic_generator(g, r, x_0,1000)
        for i in range(999):
            next(it)
        xs.append(next(it))
    plt.plot(rs, xs, '.--')
    plt.xlabel('$r$'); plt.ylabel('$x^*$')
    plt.tight_layout()
    plt.savefig(f'p2_{r}.eps'); plt.clf()
    # end:problem_2 #

    # problem_3 #
    r = 3.1
    for x_0 in (random.random() for i in range(3)):
        plt.plot(list(logistic_generator(g, r, x_0,50)), '.--')
    plt.xlabel('$n$'); plt.ylabel('$x^*$')
    plt.tight_layout()
    plt.savefig(f'p3_{r}.eps'); plt.clf()
    # end:problem_3 #

    # problem_4 #
    N = 100
    rs += [i/N for i in range(3*N, int(3.44*N))]
    for r in rs:
        x_0 = random.random()
        total = 1000
        it = logistic_generator(g, r, x_0,total)
        for i in range(total-2):
            next(it)
        obj = list(set(it))
        plt.plot([r]*len(obj),obj, '.', c='C0')
    plt.xlabel('$r$'); plt.ylabel('$x^*$')
    plt.tight_layout()
    plt.savefig(f'p4_{r}.eps'); plt.clf()
    # end:problem_4 #

    # problem_5 #
    r = 3.5
    plt.plot(list(logistic_generator(g, r, random.random(), 250)), '.')
    plt.xlabel(f'$n$, $r={r}$'); plt.ylabel('$x^*$')
    plt.tight_layout()
    plt.savefig(f'p5_{r}.eps'); plt.clf()
    r = 3.55
    plt.plot(list(logistic_generator(g, r, random.random(), 250)), '.')
    plt.xlabel(f'$n$, $r={r}$'); plt.ylabel('$x^*$')
    plt.tight_layout()
    plt.savefig(f'p5_{r}.eps'); plt.clf()
    # mid:problem_5 #
    N = 100
    rs = [i/N for i in range(4*N)]
    speeds = []
    for r in rs:
        speeds.append(converge_speed(g, r))
    plt.plot(rs, speeds, '.--')
    plt.xlabel('$r$'); plt.ylabel('Converging Speed')
    plt.tight_layout()
    plt.savefig('p5_speed_line.eps'); plt.clf()
    # end:problem_5 #

    # problem_6 #
    N = 300
    total = 10000
    ry = 4
    count = 0
    for px, py in zip(
        [0, 1, 3, 3.449, 3.5439, 3.5643],
        [0, 0, 2/3, 0.8499, 0.88402, 0.89076]
    ):
        rs = [px + i*(ry-px)/N for i in range(N)]
        for r in rs:
            x_0 = random.random()
            it = logistic_generator(g, r, x_0, total)
            for i in range(total-128):
                next(it)
            obj = list(set(it))
            plt.plot([r]*len(obj), obj, '.', c='C0')
        plt.xlabel('$r$' + (f', from $T={2**(count-1)}$ to {2**count}' if count else ''))
        plt.ylabel('$x^*$')
        margin = (max(obj) - py)/50
        plt.ylim(top=max(obj)+margin, bottom=py-margin)
        plt.tight_layout()
        plt.savefig(f'p6_{ry}_{count}.eps'); plt.clf()
        count += 1
    ry = 3.56994592
    count = 0
    for px, py in zip(
        [0, 1, 3, 3.449, 3.5439, 3.5643, 3.56874, 3.569684,
            3.569888, 3.5699334, 3.56994271],
        [0, 0, 2/3, 0.8499, 0.88402, 0.89076, 0.892132, 0.8924125,
            0.8924706, 0.8924831, 0.89248563]
    ):
        rs = [px + i*(ry-px)/N for i in range(N)]
        for r in rs:
            x_0 = random.random()
            it = logistic_generator(g, r, x_0, total)
            for i in range(total-2*1024):
                next(it)
            obj = list(set(it))
            plt.plot([r]*len(obj), obj, '.', c='C0')
        plt.xlabel('$r$' + (f', from $T={2**(count-1)}$ to {2**count}' if count else ''))
        plt.ylabel('$x^*$')
        margin = (max(obj) - py)/50
        plt.ylim(top=max(obj)+margin, bottom=py-margin)
        plt.tight_layout()
        plt.savefig(f'p6_{ry}_{count}.eps'); plt.clf()
        count += 1
    # end:problem_6 #

    # problem_7
    rs = [0, 1, 3, 3.449, 3.5439, 3.5643, 3.56874, 3.569684,
          3.569888, 3.5699334, 3.56994271]
    delta = [n - p for n, p in zip(rs[2:], rs[1:-1])]
    Fs = [p/n for n, p in zip(delta[1:], delta[:-1])]
    print(sum(Fs)/len(Fs))
    print(delta)
    print(Fs)


# log_gen #
def logistic_generator(g, r, x_0, n=50):
    '''g: function, r, x_0, n -> generator of Logistic f(x) = r*g(x)
    the iteration will repeat `n' times
    '''
    for i in range(n):
        yield x_0
        x_0 = r*g(x_0)
# end:log_gen #

# converge_speed #
def converge_speed(g, r, n=30000, skip=5):
    '''g: function, r -> converge speed of Logistic f(x) = r*g(x)
    for series with cycle other than 1, this speed is defined as that
    of each cycle.
    n: calculate the series till n
    '''
    # find cycle T
    it = logistic_generator(g, r, random.random(), n)
    lst = list(it)
    T = len(set(lst[-512:]))
    # find x*
    x_star = lst[skip+((n-skip)//T-1)*T]
    # get sub series
    res = []
    old_error = lst[skip] - x_star
    for i in range(1, 50):
        error = lst[skip+i*T] - x_star
        if error:
            res.append(abs(error/old_error))
            old_error = error
        else:
            break
    res = sum(res)/len(res) if res else 0
    return res if res <= 1 else 1
# end:converge_speed #

if __name__ == '__main__':
    main()
