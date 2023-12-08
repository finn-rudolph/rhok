#include <bits/stdc++.h>
using namespace std;

constexpr size_t N = 1000, M = 4, F = 2, P = 50;

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

double p_mu_lambda[N][N];

// the assigned k and the gcd
uint32_t k[M], h[M];
double lg_k[M];

double mu_lambda(size_t machine = 0, size_t factor = 0, double curr_min = DBL_MAX)
{
    if (machine == M)
        return curr_min;
    if (factor == F)
        return mu_lambda(machine + 1, 0, curr_min);

    double r = 0.0, prob_so_far = 0.0;
    for (size_t mu_plus_lambda = 1; mu_plus_lambda <= P / (h[machine] - 1); ++mu_plus_lambda)
    {
        if (mu_plus_lambda * lg_k[machine] >= curr_min)
        {
            assert(prob_so_far < 1.0);
            return r + (1.0 - prob_so_far) * mu_lambda(machine, factor + 1, curr_min);
        }

        prob_so_far += p_mu_lambda[P / (h[machine] - 1)][mu_plus_lambda] * mu_plus_lambda;
        r += p_mu_lambda[P / (h[machine] - 1)][mu_plus_lambda] * mu_plus_lambda *
             mu_lambda(machine, factor + 1, mu_plus_lambda * lg_k[machine]);
    }
    return r;
}

double expected(size_t machine = 0, size_t factor = 0)
{
    if (machine == M)
        return mu_lambda();
    if (factor == F)
        return expected(machine + 1);

    double e = 0.0;
    for (auto d : nat::div[k[machine]])
    {
        h[machine] = 2 * d;
        e += expected(machine, factor + 1) * nat::gcd_prob(k[machine], d);
    }
    return e;
}

void iterate_k_cartesian_prod(size_t machine = 0)
{
    if (machine == M)
    {
        for (size_t i = 0; i < M; ++i)
            printf("%-6" PRIu32, k[i]);
        printf("%lf\n", expected());
        fflush(stdout);
        return;
    }

    size_t start = !machine ? 1 : k[machine - 1];
    for (size_t i = start; i <= 12; ++i)
    {
        k[machine] = i;
        lg_k[machine] = log2(2.0 * i);
        iterate_k_cartesian_prod(machine + 1);
    }
}

int main()
{
    nat::init();

    for (size_t q = 1; q < N; ++q)
    {
        p_mu_lambda[q][1] = 1.0 / q;
        for (size_t i = 2; i < N; ++i)
        {
            p_mu_lambda[q][i] = p_mu_lambda[q][i - 1] * ((double)(q - i + 1) / q);
        }
    }

    cerr << "commencing" << endl;
    iterate_k_cartesian_prod();
}