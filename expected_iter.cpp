#include <bits/stdc++.h>
using namespace std;

constexpr size_t N = 1 << 6, F = 1, M = 2;

namespace nat
{
    vector<uint32_t> div[N];
    uint32_t phi[N];

    void init()
    {
        auto r = vector<uint32_t>(N);

        for (size_t i = 2; i < N; ++i)
            if (!r[i])
            {
                r[i] = i;
                for (size_t j = i * i; j < N; j += i)
                    r[j] = i;
            }

        for (size_t i = 1; i < N; ++i)
            for (size_t j = i; j < N; j += i)
                div[j].push_back(i);

        phi[1] = 1;
        for (size_t i = 2; i < N; ++i)
        {
            size_t q = 1, n = i;
            while (!(n % r[i]))
            {
                q *= r[i];
                n /= r[i];
            }

            if (n == 1)
                phi[i] = (r[i] - 1) * (q / r[i]);
            else
                phi[i] = phi[q] * phi[i / q];
        }
    }

    double gcd_prob(uint32_t k, uint32_t d) { return (double)phi[k / d] / k; }
}

double expected_m1(uint32_t k, size_t factor = 0, uint32_t d_sum = 0)
{
    if (factor == F)
        return 1.0 / sqrt(d_sum);

    double e = 0.0;
    for (auto d : nat::div[k])
        e += expected_m1(k, factor + 1, d_sum + d) * nat::gcd_prob(k, d);
    return e;
}

uint32_t _k[M];
double inv_lg_k[M];

double expected_f1(size_t machine = 0, double d_div_lg_sum = 0.0)
{
    if (machine == M)
        return 1.0 / sqrt(d_div_lg_sum);

    double e = 0.0;
    for (auto d : nat::div[_k[machine]])
        e += expected_f1(machine + 1, d_div_lg_sum +
                                          (double)d * inv_lg_k[machine]) *
             nat::gcd_prob(_k[machine], d);
    return e;
}

int main()
{
    nat::init();
    printf("%-8s%-8sexpected iterations / sqrt((p - 1) ln 2)\n", "k_0", "k_1");

    for (uint32_t k1 = 1; k1 < N; ++k1)
    {
        _k[0] = k1;
        inv_lg_k[0] = 1.0 / log2(2.0 * k1);
        inv_lg_k[0] *= inv_lg_k[0];
        for (uint32_t k2 = k1 + 1; k2 < N; ++k2)
        {
            _k[1] = k2;
            inv_lg_k[1] = 1.0 / log2(2.0 * k2);
            inv_lg_k[1] *= inv_lg_k[1];
            printf("%-8" PRIu32 "%-8" PRIu32 "%lf\n", k1, k2, expected_f1());
        }
    }
}