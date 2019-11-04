#ifndef UTIL_H
#define UTIL_H

#include <tuple>
#include <unordered_map>

namespace ptmscoring {

    class LogMath {
        public:
            float log_sum(float, float);
            float log_bin_coef(size_t, size_t);
    };

    class BinomialDist {
        LogMath lmath;
        float log_prob_success, log_prob_fail;
        std::unordered_map<uint32_t, float> cache;
        public:
            BinomialDist();
            BinomialDist(float);
            float log_pmf(size_t, size_t);
            float log_pvalue(size_t, size_t);
            float log10_pvalue(size_t, size_t);
            void initialize_cache(size_t, float);
    };
}

#endif
