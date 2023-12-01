#include <bits/stdc++.h>
using namespace std;

constexpr size_t N = 1 << 12, F = 2;

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
                for (size_t j = i; j < N; j += i)
                    r[j] = i;
            }

        for (size_t i = 1; i < N; ++i)
            for (size_t j = i; j < N; j += i)
                div[j].push_back(i);

        phi[1] = 1;
        for (size_t i = 2; i < N; ++i)
        {
            for (auto d : div[i])
                printf("%" PRIu32 " ", d);
            printf("\n");
            size_t q = 1, n = i;
            while (!(n % r[i]))
            {
                q *= r[i];
                n /= r[i];
            }

            if (q == i)
                phi[i] = (r[i] - 1) * (q / r[i]);
            else
                phi[i] = phi[q] * phi[i / q];

            // printf("phi(%zu) = %" PRIu32 "\n", i, phi[i]);
        }
    }

    double gcd_prob(uint32_t k, uint32_t d)
    {
        assert(!(k % d));
        return (double)phi[k / d] / (double)k;
    }
}

double expected(uint32_t k, size_t i = 0, uint32_t d_sum = 0)
{
    if (i == F)
        return 1.0 / sqrt((double)d_sum);

    double e = 0.0;
    for (auto d : nat::div[k])
        e += expected(k, i + 1, d_sum + d) * nat::gcd_prob(k, d);
    return e;
}

int main()
{
    nat::init();
    printf("%-8cexpected iterations / sqrt((p - 1) ln 2)\n", 'k');
    for (uint32_t k = 1; k < N; ++k)
        printf("%-8" PRIu32 "%lf\n", k, expected(k) * log2((double)k * 2.0));
}