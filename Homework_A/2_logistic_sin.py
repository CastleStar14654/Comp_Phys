import random
import matplotlib.pyplot as plt
from matplotlib import rc
from math import sin, pi


def main():
    # gx #
    def g(x):
        return sin(x*pi)
    # end:gx #

    # problem_1 #
    r = 0.1
    for x_0 in (random.random() for i in range(10)):
        plt.plot(list(logistic_generator(g, r, x_0, 11)), '.--')
    plt.xlabel('$n$'); plt.ylabel('$x^*$')
    plt.tight_layout()
    plt.savefig(f'p9.1_{r}.eps'); plt.clf()
    # plt.show()
    r = 0.5
    for x_0 in (random.random() for i in range(10)):
        plt.plot(list(logistic_generator(g, r, x_0, 15)), '.--')
    plt.xlabel('$n$'); plt.ylabel('$x^*$')
    plt.tight_layout()
    plt.savefig(f'p9.1_{r}.eps'); plt.clf()
    # plt.show()
    # end:problem_1 #

    # # problem_2 #
    N = 150
    rs = [i/N for i in range(int(0.72*N))]
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
    plt.savefig(f'p9.2_{r}.eps'); plt.clf()
    # plt.show()
    # end:problem_2 #

    # problem_3 #
    r = 0.82
    for x_0 in (random.random() for i in range(3)):
        plt.plot(list(logistic_generator(g, r, x_0,50)), '.--')
    plt.xlabel('$n$'); plt.ylabel('$x^*$')
    plt.tight_layout()
    plt.savefig(f'p9.3_{r}.eps'); plt.clf()
    # plt.show()
    # end:problem_3 #

    # problem_4 #
    N = 150
    rs += [i/N for i in range(int(0.72*N), int(0.831*N))]
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
    plt.savefig(f'p9.4_{r}.eps'); plt.clf()
    # plt.show()
    # end:problem_4 #

    # problem_5 #
    r = 0.84
    plt.plot(list(logistic_generator(g, r, random.random(), 250)), '.')
    plt.xlabel(f'$n$, $r={r}$'); plt.ylabel('$x^*$')
    plt.tight_layout()
    plt.savefig(f'p9.5_{r}.eps'); plt.clf()
    # plt.show()
    r = 0.86
    plt.plot(list(logistic_generator(g, r, random.random(), 250)), '.')
    plt.xlabel(f'$n$, $r={r}$'); plt.ylabel('$x^*$')
    plt.tight_layout()
    plt.savefig(f'p9.5_{r}.eps'); plt.clf()
    # plt.show()
    # mid:problem_5 #
    N = 400
    rs = [i/N for i in range(1*N)]
    speeds = []
    for r in rs:
        speeds.append(converge_speed(g, r))
    plt.plot(rs, speeds, '.--')
    plt.xlabel('$r$'); plt.ylabel('Converging Speed')
    plt.tight_layout()
    plt.savefig('p9.5_speed_line.eps'); plt.clf()
    # plt.show()
    # end:problem_5 #

    # problem_6 #
    N = 300
    total = 10000
    ry = 1
    count = 0
    for px, py in zip(
        [0, 0.317, 0.719, 0.8328, 0.85648],
        [0, 0, 0.645, 0.8205, 0.85650]
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
        plt.savefig(f'p9.6_{ry}_{count}.eps'); plt.clf()
        # plt.show()
        count += 1
    ry = 0.86557934
    count = 0
    for px, py in zip(
        [0, 0.317, 0.719, 0.8328, 0.85848, 0.864062,
         0.865252, 0.8655092, 0.8655637],
        [0, 0, 0.645, 0.8205, 0.85650, 0.863745,
         0.865201, 0.865501, 0.8655625]
    ):
        rs = [px + i*(ry-px)/N for i in range(N)]
        for r in rs:
            x_0 = random.random()
            it = logistic_generator(g, r, x_0, total)
            for i in range(total-2*1024):
                next(it)
            obj = list(set(it))
            plt.plot([r]*len(obj), obj, '.', c='C0')
        plt.xlabel(
            '$r$' + (f', from $T={2**(count-1)}$ to {2**count}' if count else ''))
        plt.ylabel('$x^*$')
        margin = (max(obj) - py)/50
        plt.ylim(top=max(obj)+margin, bottom=py-margin)
        plt.tight_layout()
        plt.savefig(f'p9.6_{ry}_{count}.eps'); plt.clf()
        # plt.show()
        count += 1
    end: problem_6

    # problem_7
    rs = [0, 0.317, 0.719, 0.8328, 0.85848,
          0.864062, 0.865252, 0.8655092, 0.8655637]
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
