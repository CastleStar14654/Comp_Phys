#include <random>
#include <iostream>
#include <vector>
#include <array>
#include <bitset>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <iterator>
#include <functional>
#include <bit>
#include <thread>

using namespace std;

/*beg:Ising_dec*/
template <size_t L, size_t HIS_SIZE = 21013>
class Ising_2D
{
private:
    array<bitset<L>, L> data;
    ranlux48 ran;
    double _beta;
    double _H;
    array<double, HIS_SIZE> spin_his;        // history of spin()
    array<double, HIS_SIZE> hamiltonian_his; // history of hamiltonian()
    size_t history_data_count{0ULL};         // count of data sent to history
    size_t skip_at_start;
    size_t start_skip_counter{0};
    const size_t ave_history_skip{47};

    double ave(function<double(double, double)> func_spin_hamiltonian) const;
    double delta_hamiltonian(size_t x, size_t y) const; // change of hamiltonian if flip (x, y)

public:
    Ising_2D() {}
    Ising_2D(ranlux48::result_type s) : ran{s} {}
    void set_seed(ranlux48::result_type s) { ran.seed(s); }

    // spin at x, y
    bool at(size_t x, size_t y) const { return data[x][y]; }
    double beta() const { return _beta; }
    double mag_field() const { return _H; }
    /* set `_beta` and `_H`, then iterate `times` */
    void iterate(double beta, double H, size_t times, size_t start_skip,
                 bool random_start = true, size_t output_skip = 0, ostream *p_os = &cout);
    // beta * hamiltonian
    double hamiltonian() const;
    // space ensemble average beta * hamiltonian
    double energy() const;
    // \partial energy() / \partial T
    double c_h() const;
    // space average spin
    double spin() const;
    // space ensemble average spin
    double mag_sus() const;
    double mag_sus_abs() const;
    // \partial mag_sus() / \partial H
    double chi() const;
    double chi_abs() const;
};

template <size_t L, size_t HIS_SIZE>
ostream &operator<<(ostream &os, const Ising_2D<L, HIS_SIZE> &ising);
/*end:Ising_dec*/

// ================================================

int main(int argc, const char **argv)
{
    default_random_engine ran{};
    cout << setprecision(10);
    /*p2.1*/
    if (argc > 1 && argv[1][0] == '-' && argv[1][1] == '1')
    {
        cout << "=======================\n"
             << "PROBLEM 2.1" << endl;
        // note: default seed is 19780503
        for (ranlux48::result_type seed : {35641, 19780503, 14654, 6244212})
        {
            Ising_2D<32> ising_32{seed};
            double beta{.5};
            double H{0.};
            ofstream ofs{"output/ising_2.1_" + to_string(seed) + ".out", ios_base::binary};
            ising_32.iterate(beta, H, 32 * 32 * 2048 * 2, 32 * 32 * 512, true, 32 * 32, &ofs);
        }
    }
    /*p2.2*/
    if (argc > 1 && argv[1][0] == '-' && argv[1][1] == '2')
    {
        cout << "=======================\n"
             << "PROBLEM 2.2\n"
             << "#\tbeta\tenergy\t\tC_H\t\tM\t\tM_abs\t\tchi\t\tchi_abs" << endl;
        // note: default seed is 19780503
        double H{0.};
        double beta{0.00};
        vector<Ising_2D<32> *> ps_ising_32{};
        vector<ofstream> ofss;
        for (size_t i = 0; i < 4; i++)
        {
            ps_ising_32.push_back(new Ising_2D<32>());
            ps_ising_32.back()->iterate(beta, H, 0, 0, true);
            ofss.emplace_back("output/ising_2.2_" + to_string(i) + ".out", ios_base::binary);
        }
        for (; beta < 1.01; beta += 0.01)
        {
            vector<thread> threads{};
            for (size_t i = 0; i < 4; i++)
            {
                ps_ising_32[i]->set_seed(ran());
                threads.emplace_back(
                    ps_ising_32[i]->iterate, ps_ising_32[i],
                    beta, H, 32 * 32 * 2048, 32 * 32 * 1024, false, 32 * 32, &ofss[i]);
            }
            for (auto &th : threads)
            {
                th.join();
            }
            for (size_t i = 0; i < 4; i++)
            {
                cout << i << '\t'
                     << beta << '\t'
                     << ps_ising_32[i]->energy() << '\t'
                     << ps_ising_32[i]->c_h() << '\t'
                     << ps_ising_32[i]->mag_sus() << '\t'
                     << ps_ising_32[i]->mag_sus_abs() << '\t'
                     << ps_ising_32[i]->chi() << '\t'
                     << ps_ising_32[i]->chi_abs() << endl;
            }
        }
        for (auto &&p : ps_ising_32)
        {
            delete p;
        }
    }
    /*p2.3*/
    if (argc > 1 && argv[1][0] == '-' && argv[1][1] == '3')
    {
        cout << "=======================\n"
             << "PROBLEM 2.3\n";
        if (argc > 2 && argv[2][0] == '1' && argv[2][1] == '6')
        {
            cout << "L=16 ------------------\n"
                 << "#\tbeta\tenergy\t\tC_H\t\tM\t\tM_abs\t\tchi\t\tchi_abs" << endl;
            double H{0.};
            double beta{0.30};
            vector<Ising_2D<16, 3407> *> ps_ising_16{};
            vector<ofstream> ofss;
            for (size_t i = 0; i < 4; i++)
            {
                ps_ising_16.push_back(new Ising_2D<16, 3407>());
                ps_ising_16.back()->iterate(beta, H, 0, 0, true);
                ofss.emplace_back("output/ising_2.3_16_" + to_string(i) + ".out", ios_base::binary);
            }
            for (; beta < 0.501; beta += 0.008)
            {
                vector<thread> threads{};
                for (size_t i = 0; i < 4; i++)
                {
                    ps_ising_16[i]->set_seed(ran());
                    threads.emplace_back(
                        ps_ising_16[i]->iterate, ps_ising_16[i],
                        beta + 0.002 * i, H, 16 * 16 * 2048, 16 * 16 * 1024, false, 16 * 16, &ofss[i]);
                }
                for (auto &th : threads)
                {
                    th.join();
                }
                for (size_t i = 0; i < 4; i++)
                {
                    cout << i << '\t'
                         << beta + 0.002 * i << '\t'
                         << ps_ising_16[i]->energy() << '\t'
                         << ps_ising_16[i]->c_h() << '\t'
                         << ps_ising_16[i]->mag_sus() << '\t'
                         << ps_ising_16[i]->mag_sus_abs() << '\t'
                         << ps_ising_16[i]->chi() << '\t'
                         << ps_ising_16[i]->chi_abs() << endl;
                }
            }
            for (auto &&p : ps_ising_16)
            {
                delete p;
            }
        }
        if (argc > 2 && argv[2][0] == '2' && argv[2][1] == '4')
        {
            cout << "L=24 ------------------\n"
                 << "#\tbeta\tenergy\t\tC_H\t\tM\t\tM_abs\t\tchi\t\tchi_abs" << endl;
            double H{0.};
            double beta{0.30};
            vector<Ising_2D<24, 4591> *> ps_ising_24{};
            vector<ofstream> ofss;
            for (size_t i = 0; i < 4; i++)
            {
                ps_ising_24.push_back(new Ising_2D<24, 4591>());
                ps_ising_24.back()->iterate(beta, H, 0, 0, true);
                ofss.emplace_back("output/ising_2.3_24_" + to_string(i) + ".out", ios_base::binary);
            }
            for (; beta < 0.501; beta += 0.008)
            {
                vector<thread> threads{};
                for (size_t i = 0; i < 4; i++)
                {
                    ps_ising_24[i]->set_seed(ran());
                    threads.emplace_back(
                        ps_ising_24[i]->iterate, ps_ising_24[i],
                        beta + 0.002 * i, H, 24 * 24 * 2048, 24 * 24 * 1024, false, 24 * 24, &ofss[i]);
                }
                for (auto &th : threads)
                {
                    th.join();
                }
                for (size_t i = 0; i < 4; i++)
                {
                    cout << i << '\t'
                         << beta + 0.002 * i<< '\t'
                         << ps_ising_24[i]->energy() << '\t'
                         << ps_ising_24[i]->c_h() << '\t'
                         << ps_ising_24[i]->mag_sus() << '\t'
                         << ps_ising_24[i]->mag_sus_abs() << '\t'
                         << ps_ising_24[i]->chi() << '\t'
                         << ps_ising_24[i]->chi_abs() << endl;
                }
            }
            for (auto &&p : ps_ising_24)
            {
                delete p;
            }
        }
        if (argc > 2 && argv[2][0] == '3' && argv[2][1] == '2')
        {
            cout << "L=32 ------------------\n"
                 << "#\tbeta\tenergy\t\tC_H\t\tM\t\tM_abs\t\tchi\t\tchi_abs" << endl;
            double H{0.};
            double beta{0.35};
            vector<Ising_2D<32> *> ps_ising_32{};
            vector<ofstream> ofss;
            for (size_t i = 0; i < 4; i++)
            {
                ps_ising_32.push_back(new Ising_2D<32>());
                ps_ising_32.back()->iterate(beta, H, 0, 0, true);
                ofss.emplace_back("output/ising_2.3_32_" + to_string(i) + ".out", ios_base::binary);
            }
            for (; beta < 0.501; beta += 0.008)
            {
                vector<thread> threads{};
                for (size_t i = 0; i < 4; i++)
                {
                    ps_ising_32[i]->set_seed(ran());
                    threads.emplace_back(
                        ps_ising_32[i]->iterate, ps_ising_32[i],
                        beta + 0.002 * i, H, 32 * 32 * 3072, 32 * 32 * 1024, false, 32 * 32, &ofss[i]);
                }
                for (auto &th : threads)
                {
                    th.join();
                }
                for (size_t i = 0; i < 4; i++)
                {
                    cout << i << '\t'
                         << beta + 0.002 * i<< '\t'
                         << ps_ising_32[i]->energy() << '\t'
                         << ps_ising_32[i]->c_h() << '\t'
                         << ps_ising_32[i]->mag_sus() << '\t'
                         << ps_ising_32[i]->mag_sus_abs() << '\t'
                         << ps_ising_32[i]->chi() << '\t'
                         << ps_ising_32[i]->chi_abs() << endl;
                }
            }
            for (auto &&p : ps_ising_32)
            {
                delete p;
            }
        }
        if (argc > 2 && argv[2][0] == '4' && argv[2][1] == '0')
        {
            cout << "L=40 ------------------\n"
                 << "#\tbeta\tenergy\t\tC_H\t\tM\t\tM_abs\t\tchi\t\tchi_abs" << endl;
            double H{0.};
            double beta{0.35};
            vector<Ising_2D<40, 26573> *> ps_ising_40{};
            vector<ofstream> ofss;
            for (size_t i = 0; i < 4; i++)
            {
                ps_ising_40.push_back(new Ising_2D<40, 26573>());
                ps_ising_40.back()->iterate(beta, H, 0, 0, true);
                ofss.emplace_back("output/ising_2.3_40_" + to_string(i) + ".out", ios_base::binary);
            }
            for (; beta < 0.481; beta += 0.006)
            {
                vector<thread> threads{};
                for (size_t i = 0; i < 4; i++)
                {
                    ps_ising_40[i]->set_seed(ran());
                    threads.emplace_back(
                        ps_ising_40[i]->iterate, ps_ising_40[i],
                        beta + 0.0015 * i, H, 40 * 40 * 4096, 40 * 40 * 2048, false, 40 * 40, &ofss[i]);
                }
                for (auto &th : threads)
                {
                    th.join();
                }
                for (size_t i = 0; i < 4; i++)
                {
                    cout << i << '\t'
                         << beta + 0.0015 * i << '\t'
                         << ps_ising_40[i]->energy() << '\t'
                         << ps_ising_40[i]->c_h() << '\t'
                         << ps_ising_40[i]->mag_sus() << '\t'
                         << ps_ising_40[i]->mag_sus_abs() << '\t'
                         << ps_ising_40[i]->chi() << '\t'
                         << ps_ising_40[i]->chi_abs() << endl;
                }
            }
            for (auto &&p : ps_ising_40)
            {
                delete p;
            }
        }
        if (argc > 2 && argv[2][0] == '4' && argv[2][1] == '8')
        {
            cout << "L=48 ------------------\n"
                 << "#\tbeta\tenergy\t\tC_H\t\tM\t\tM_abs\t\tchi\t\tchi_abs" << endl;
            double H{0.};
            double beta{0.35};
            vector<Ising_2D<48, 48649> *> ps_ising_48{};
            vector<ofstream> ofss;
            for (size_t i = 0; i < 4; i++)
            {
                ps_ising_48.push_back(new Ising_2D<48, 48649>());
                ps_ising_48.back()->iterate(beta, H, 0, 0, true);
                ofss.emplace_back("output/ising_2.3_48_" + to_string(i) + ".out", ios_base::binary);
            }
            for (; beta < 0.481; beta += 0.006)
            {
                vector<thread> threads{};
                for (size_t i = 0; i < 4; i++)
                {
                    ps_ising_48[i]->set_seed(ran());
                    threads.emplace_back(
                        ps_ising_48[i]->iterate, ps_ising_48[i],
                        beta + 0.0015 * i, H, 48 * 48 * 4096, 48 * 48 * 2048, false, 48 * 48, &ofss[i]);
                }
                for (auto &th : threads)
                {
                    th.join();
                }
                for (size_t i = 0; i < 4; i++)
                {
                    cout << i << '\t'
                         << beta + 0.0015 * i<< '\t'
                         << ps_ising_48[i]->energy() << '\t'
                         << ps_ising_48[i]->c_h() << '\t'
                         << ps_ising_48[i]->mag_sus() << '\t'
                         << ps_ising_48[i]->mag_sus_abs() << '\t'
                         << ps_ising_48[i]->chi() << '\t'
                         << ps_ising_48[i]->chi_abs() << endl;
                }
            }
            for (auto &&p : ps_ising_48)
            {
                delete p;
            }
        }
        if (argc > 2 && argv[2][0] == '5' && argv[2][1] == '6')
        {
            cout << "L=56 ------------------\n"
                 << "#\tbeta\tenergy\t\tC_H\t\tM\t\tM_abs\t\tchi\t\tchi_abs" << endl;
            double H{0.};
            double beta{0.38};
            vector<Ising_2D<56, 49993> *> ps_ising_56{};
            vector<ofstream> ofss;
            for (size_t i = 0; i < 4; i++)
            {
                ps_ising_56.push_back(new Ising_2D<56, 49993>());
                ps_ising_56.back()->iterate(beta, H, 0, 0, true);
                ofss.emplace_back("output/ising_2.3_56_" + to_string(i) + ".out", ios_base::binary);
            }
            for (; beta < 0.481; beta += 0.004)
            {
                vector<thread> threads{};
                for (size_t i = 0; i < 4; i++)
                {
                    ps_ising_56[i]->set_seed(ran());
                    threads.emplace_back(
                        ps_ising_56[i]->iterate, ps_ising_56[i],
                        beta + 0.001 * i, H, 56 * 56 * 4096, 56 * 56 * 2048, false, 56 * 56, &ofss[i]);
                }
                for (auto &th : threads)
                {
                    th.join();
                }
                for (size_t i = 0; i < 4; i++)
                {
                    cout << i << '\t'
                         << beta + 0.001 * i<< '\t'
                         << ps_ising_56[i]->energy() << '\t'
                         << ps_ising_56[i]->c_h() << '\t'
                         << ps_ising_56[i]->mag_sus() << '\t'
                         << ps_ising_56[i]->mag_sus_abs() << '\t'
                         << ps_ising_56[i]->chi() << '\t'
                         << ps_ising_56[i]->chi_abs() << endl;
                }
            }
            for (auto &&p : ps_ising_56)
            {
                delete p;
            }
        }
        if (argc > 2 && argv[2][0] == '6' && argv[2][1] == '4')
        {
            cout << "L=64 ------------------\n"
                 << "#\tbeta\tenergy\t\tC_H\t\tM\t\tM_abs\t\tchi\t\tchi_abs" << endl;
            double H{0.};
            double beta{0.38};
            vector<Ising_2D<64, 49993> *> ps_ising_64{};
            vector<ofstream> ofss;
            for (size_t i = 0; i < 4; i++)
            {
                ps_ising_64.push_back(new Ising_2D<64, 49993>());
                ps_ising_64.back()->iterate(beta, H, 0, 0, true);
                ofss.emplace_back("output/ising_2.3_64_" + to_string(i) + ".out", ios_base::binary);
            }
            for (; beta < 0.461; beta += 0.004)
            {
                vector<thread> threads{};
                for (size_t i = 0; i < 4; i++)
                {
                    ps_ising_64[i]->set_seed(ran());
                    threads.emplace_back(
                        ps_ising_64[i]->iterate, ps_ising_64[i],
                        beta + 0.001 * i, H, 64 * 64 * 4096, 64 * 64 * 2048, false, 64 * 64, &ofss[i]);
                }
                for (auto &th : threads)
                {
                    th.join();
                }
                for (size_t i = 0; i < 4; i++)
                {
                    cout << i << '\t'
                         << beta + 0.001 * i<< '\t'
                         << ps_ising_64[i]->energy() << '\t'
                         << ps_ising_64[i]->c_h() << '\t'
                         << ps_ising_64[i]->mag_sus() << '\t'
                         << ps_ising_64[i]->mag_sus_abs() << '\t'
                         << ps_ising_64[i]->chi() << '\t'
                         << ps_ising_64[i]->chi_abs() << endl;
                }
            }
            for (auto &&p : ps_ising_64)
            {
                delete p;
            }
        }
    }
    /*p2.4*/
    if (argc > 1 && argv[1][0] == '-' && argv[1][1] == '4')
    {
        cout << "=======================\n"
             << "PROBLEM 2.4" << endl;
        double H{0.};
        Ising_2D<32> ising_32{};
        for (double beta : {0.2, 0.42, 0.6})
        {
            cout << "beta=" << beta << " ----------------" << endl;
            ising_32.iterate(beta, H, 32 * 32 * 2048, 32 * 32 * 512, true, 0, nullptr);
            cout << ising_32 << endl;
        }
    }
    /*p2.5*/
    if (argc > 1 && argv[1][0] == '-' && argv[1][1] == '5')
    {
        cout << "=======================\n"
             << "PROBLEM 2.5" << endl;
        double beta{0.5};
        double H0{1.};
        Ising_2D<32, 10> ising_32{};
        for (size_t N : {512, 1024, 2048})
        {
            for (size_t M : {512, 1024, 2048})
            {
                cout << "M=" << M << ", N=" << N << endl;
                ofstream ofs{"output/ising_2.5_" + to_string(M) + "_" + to_string(N) + ".out",
                             ios_base::binary};
                ising_32.set_seed(ran());
                ising_32.iterate(beta, H0, 32 * 32 * 512, 32 * 32 * 512, true, 0, nullptr);
                for (size_t k = 0; k < 3 * N; k++)
                {
                    double H{H0 * cos(k * 6.283185307179586 / N)};
                    ising_32.iterate(beta, H, M, M, false);
                    double spin{ising_32.spin()};
                    ofs.write(static_cast<char *>((void *)(&H)), sizeof(double));
                    ofs.write(static_cast<char *>((void *)(&spin)), sizeof(double));
                }
            }
        }
    }
    /*pend*/
}

// ================================================

/*beg:average*/
template <size_t L, size_t HIS_SIZE>
double Ising_2D<L, HIS_SIZE>::ave(function<double(double, double)> func_spin_hamiltonian) const
{
    if (history_data_count >= HIS_SIZE)
    {
        double sum{};
        double compensation{};
        for (size_t i = 0; i < HIS_SIZE; i++)
        {
            double to_add{func_spin_hamiltonian(spin_his[i], hamiltonian_his[i])};
            to_add -= compensation;
            double temp{sum + to_add};
            compensation = temp - sum - to_add;
            sum = temp;
        }
        return sum / HIS_SIZE;
    }
    else
    {
        return numeric_limits<double>::quiet_NaN();
    }
}
/*end:average*/

/*beg:delta_h*/
template <size_t L, size_t HIS_SIZE>
double Ising_2D<L, HIS_SIZE>::delta_hamiltonian(size_t x, size_t y) const
{
    double res{-2.};
    size_t pos_count{0}; // number of positive around it
    pos_count += data[x][(y + 1) % L];
    pos_count += data[x][(L - 1 + y) % L];
    pos_count += data[(x + 1) % L][y];
    pos_count += data[(L - 1 + x) % L][y];
    res *= pos_count * 2. - 4.;
    res -= 2. * mag_field();
    if (data[x][y])
        res *= -1;
    return res;
}
/*end:delta_h*/

/*beg:iterate*/
template <size_t L, size_t HIS_SIZE>
void Ising_2D<L, HIS_SIZE>::iterate(double beta, double H, size_t times, size_t start_skip,
                                    bool random_start, size_t output_skip, ostream *p_os)
{
    _beta = beta;
    _H = H;
    if (random_start)
        for (size_t i = 0; i < L; i++)
        {
            data[i].reset();
            for (size_t j = 0; j < L; j++)
                if (uniform_int_distribution<>(0, 1)(ran))
                {
                    data[i].set(j);
                }
        }

    history_data_count = 0ULL;
    start_skip_counter = 0ULL;
    skip_at_start = start_skip;

    for (size_t n = 0; n < times; n++)
    {
        if (start_skip_counter < skip_at_start)
        {
            start_skip_counter++;
        }
        else if (uniform_int_distribution<size_t>(1, ave_history_skip)(ran) == 1)
        {
            history_data_count++;
            size_t index{history_data_count < HIS_SIZE
                             ? history_data_count
                             : uniform_int_distribution<size_t>(0, history_data_count)(ran)};
            if (index < HIS_SIZE)
            {
                hamiltonian_his[index] = hamiltonian();
                spin_his[index] = spin();
            }
        }

        if (output_skip && !(n % output_skip))
        {
            double current_spin{spin()};
            p_os->write(static_cast<char *>((void *)(&current_spin)), sizeof(double));
        }

        size_t i{uniform_int_distribution<size_t>(0, L - 1)(ran)};
        size_t j{uniform_int_distribution<size_t>(0, L - 1)(ran)};
        if (uniform_real_distribution<>(0., 1. + exp(_beta * delta_hamiltonian(i, j)))(ran) < 1.)
        // go to new position
        {
            data[i].flip(j);
        }
    }
}
/*end:iterate*/

/*beg:hamiltonian*/
template <size_t L, size_t HIS_SIZE>
double Ising_2D<L, HIS_SIZE>::hamiltonian() const
{
    double res{-mag_field() * spin() * static_cast<double>(L * L)};
    size_t diff_count{};
    for (size_t i = 0; i < L - 1; i++)
    {
        diff_count += (data[i] ^ data[i + 1]).count();
    }
    diff_count += (data[L - 1] ^ data[0]).count();
    for (size_t i = 0; i < L; i++)
    {
        bitset<L> temp{data[i] << 1};
        temp[0] = data[i][L - 1];
        diff_count += (data[i] ^ temp).count();
    }
    res -= 2. * (static_cast<double>(L * L) - diff_count);
    return res;
}
/*end:hamiltonian*/

/*beg:energy*/
template <size_t L, size_t HIS_SIZE>
double Ising_2D<L, HIS_SIZE>::energy() const
{
    static auto func{
        [](double, double h) -> double {
            return h;
        }};
    return ave(func) / static_cast<double>(L * L);
}
/*end:energy*/

/*beg:c_h*/
template <size_t L, size_t HIS_SIZE>
double Ising_2D<L, HIS_SIZE>::c_h() const
{
    double ene{energy()};
    static auto func{
        [](double, double h) -> double {
            return h * h;
        }};
    constexpr double size{static_cast<double>(L * L)};
    return beta() * beta() * (ave(func) / size - size * ene * ene);
}
/*end:c_h*/

/*beg:spin*/
template <size_t L, size_t HIS_SIZE>
double Ising_2D<L, HIS_SIZE>::spin() const
{
    size_t count{0};
    for (auto &&bs : data)
    {
        count += bs.count();
    }
    return (2. / (L * L)) * count - 1.;
}
/*end:spin*/

/*beg:mag_sus*/
template <size_t L, size_t HIS_SIZE>
double Ising_2D<L, HIS_SIZE>::mag_sus() const
{
    static auto func{
        [](double s, double) -> double {
            return s;
        }};
    return ave(func);
}
/*end:mag_sus*/

template <size_t L, size_t HIS_SIZE>
double Ising_2D<L, HIS_SIZE>::mag_sus_abs() const
{
    static auto func{
        [](double s, double) -> double {
            return abs(s);
        }};
    return ave(func);
}

/*beg:chi*/
template <size_t L, size_t HIS_SIZE>
double Ising_2D<L, HIS_SIZE>::chi() const
{
    double mag{mag_sus()};
    static auto func{
        [](double s, double) -> double {
            return s * s;
        }};
    constexpr double size{static_cast<double>(L * L)};
    return beta() * size * (ave(func) - mag * mag);
}
/*end:chi*/

template <size_t L, size_t HIS_SIZE>
double Ising_2D<L, HIS_SIZE>::chi_abs() const
{
    double mag{mag_sus_abs()};
    static auto func{
        [](double s, double) -> double {
            return s * s;
        }};
    constexpr double size{static_cast<double>(L * L)};
    return beta() * size * (ave(func) - mag * mag);
}

template <size_t L, size_t HIS_SIZE>
ostream &operator<<(ostream &os, const Ising_2D<L, HIS_SIZE> &ising)
{
    for (size_t i = 0; i < L - 1; i++)
    {
        for (size_t j = 0; j < L; j++)
        {
            os << (ising.at(i, j) ? '#' : ' ') << ' ';
        }
        os << '\n';
    }
    if (L >= 1)
        for (size_t j = 0; j < L; j++)
        {
            os << (ising.at(L - 1, j) ? '#' : ' ') << ' ';
        }
    return os;
}
