#include <cmath>
#include <vector>

#include "corax/corax.h"
#include "corax/model/gamma.h"
#include "gtest/gtest.h"

using std::vector;

double gamma_cdf(double alpha, double x) {
    return corax_compute_incomplete_gamma(alpha * x, alpha, std::lgamma(alpha));
}

double compute_wasserstein1(double alpha, const vector<double> &rate_weights, const vector<double> &rates, double step_size = 5e-5, double max_quantile = 1 - 5e-4) {
    assert(step_size > 0);
    assert(max_quantile > 0 && max_quantile <= 1);
    assert(rate_weights.size() == rates.size());
    // Assert that rates are sorted in ascending order
    for (auto it = rates.cbegin() + 1; it != rates.cend(); ++it) {
        assert(*(it - 1) <= *it);
    }

    size_t steps = 0;
    double acc = 0;
    auto it_rate = rates.cbegin();
    auto it_weight = rate_weights.cbegin();
    double discrete_rate = *it_rate;
    double discrete_weight = *it_weight;
    double rate = 0;
    double cdf_val = gamma_cdf(alpha, rate);

    while (cdf_val < max_quantile) {
        if (cdf_val > discrete_weight) {
            ++it_rate;
            ++it_weight;
            discrete_rate = *it_rate;
            discrete_weight += *it_weight;
        }
        acc += abs(rate - discrete_rate);

        cdf_val = gamma_cdf(alpha, rate);
        rate += step_size;
        ++steps;
    }

    return acc / static_cast<double>(steps);
}

template<typename T>
std::string vec2string(const vector<T> &v, const std::string &sep = "/") {
    std::string result;

    for (auto i = 0UL; i < v.size(); ++i) {
        result += std::to_string(v.at(i));

        if (i + 1 < v.size()) {
            result += sep;
        }
    }

    return result;
}

struct GammaTest : public testing::TestWithParam<std::tuple<double, unsigned int>> {};

TEST_P(GammaTest, lloyd_max_discretization) {
    const auto alpha = std::get<0>(GetParam());
    const auto ncategories = std::get<1>(GetParam());

    vector<double> rate_weights(ncategories, 0), rates(ncategories, 0);

    int ret = corax_compute_gamma_cats(alpha, ncategories, rates.data(), rate_weights.data(), CORAX_GAMMA_RATES_WEIGHTS_LLOYD_MAX);
    ASSERT_EQ(ret, CORAX_SUCCESS);
    const double wasserstein_lloyd_max = compute_wasserstein1(alpha, rate_weights, rates);
    RecordProperty("Lloyd-Max-Wasserstein1", wasserstein_lloyd_max);
    RecordProperty("Lloyd-Max-Rates", vec2string(rates));
    RecordProperty("Lloyd-Max-Weights", vec2string(rate_weights));

    ret = corax_compute_gamma_cats(alpha, ncategories, rates.data(), rate_weights.data(), CORAX_GAMMA_RATES_MEAN);
    ASSERT_EQ(ret, CORAX_SUCCESS);
    const double wasserstein_mean = compute_wasserstein1(alpha, rate_weights, rates);
    RecordProperty("Mean-Wasserstein1", wasserstein_mean);
    RecordProperty("Mean-Rates", vec2string(rates));
    RecordProperty("Mean-Weights", vec2string(rate_weights));

    // For less than 4 categories and on alpha values around 2-ish, the wasserstein distance is worse for lloyd max than for G4m.
    // TODO: determine why
    if (ncategories >= 4) {
        // Lloyd-Max should always fit the gamma distribution better than mean categories
        EXPECT_LT(wasserstein_lloyd_max, wasserstein_mean);

        if (alpha <= 2) {
            // For small values of α, we expect substantial improvement
            EXPECT_LT(wasserstein_lloyd_max / wasserstein_mean, 0.75);
        }
    }
}

INSTANTIATE_TEST_SUITE_P(GammaTests, GammaTest,
        testing::Combine(testing::Values(CORAX_OPT_MIN_ALPHA, 0.1, 0.3, 0.7, 1, 1.5, 3, 10, 20, 50, 150),
                         testing::Values(2, 3, 4, 8, 10)));
