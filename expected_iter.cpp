#include <bits/stdc++.h>
using namespace std;

constexpr size_t N = 1 << 4, F = 2, M = 6;

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

uint32_t k[M];
double inv_lg_k[M];

double expected(
    size_t machine = 0, size_t factor = 0, uint32_t d_sum = 0,
    double d_div_lg_sum = 0.0)
{
    if (machine == M)
        return 1.0 / sqrt(d_div_lg_sum);
    if (factor == F)
        return expected(
            machine + 1, 0, 0, d_div_lg_sum + d_sum * inv_lg_k[machine]);

    double e = 0.0;
    for (auto d : nat::div[k[machine]])
        e += expected(machine, factor + 1, d_sum + d, d_div_lg_sum) *
             nat::gcd_prob(k[machine], d);
    return e;
}

void iterate_k_cartesian_prod(size_t machine = 0)
{
    if (machine == M)
    {
        for (size_t i = 0; i < M; ++i)
            printf("%-6" PRIu32, k[i]);
        printf("%lf\n", expected());
        return;
    }

    size_t start = !machine ? 1 : k[machine - 1];
    for (size_t i = start; i < N; ++i)
    {
        k[machine] = i;
        inv_lg_k[machine] = 1.0 / log2(2.0 * i);
        inv_lg_k[machine] *= inv_lg_k[machine];
        iterate_k_cartesian_prod(machine + 1);
    }
}

int main()
{
    nat::init();
    for (size_t i = 0; i < M; ++i)
        printf("k%-5zu", i);
    printf("expected iterations / sqrt((p - 1) ln 2)\n");
    iterate_k_cartesian_prod();
}