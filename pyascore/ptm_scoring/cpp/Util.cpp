#include <iostream>
#include <cmath>
#include "Util.h"

namespace ptmscoring {

    uint32_t bit_concat(uint16_t a, uint16_t b) {
        return (uint32_t) a << 16 | b;
    }

    /////////////////////////////
    // LogMath Implementations //
    /////////////////////////////

    float LogMath::log_sum (float a, float b) {
        if (std::isinf(a)) {
            return b;
        } else if (std::isinf(b)){
            return a;
        } else {
            float max_val = std::max(a, b);
            float add_this = std::log(std::exp(a - max_val) + std::exp(b - max_val));
            return max_val + add_this;
        }
    }

    float LogMath::log_bin_coef (size_t k, size_t n) {
        float coef = 0;
        k = std::min(n-k, k);

        for (size_t mult = n - k + 1; mult <= n; mult++) {
            coef += std::log(mult);
        }

        for (size_t mult = 2; mult <= k; mult++) {
            coef -= std::log(mult);
        }

        return coef;
    }

    //////////////////////////////////
    // BinomialDist Implementations //
    //////////////////////////////////

    BinomialDist::BinomialDist() {
        log_prob_success = std::log(.5);
        log_prob_fail = std::log(.5);
    }

    BinomialDist::BinomialDist (float prob) {
        log_prob_success = std::log(prob);
        log_prob_fail = std::log(1. - prob);
    }

    float BinomialDist::log_pmf (size_t successes, size_t trials) {
        return lmath.log_bin_coef(successes, trials) + successes * log_prob_success + (trials - successes) * log_prob_fail;
    }

    float BinomialDist::log_pvalue (size_t successes, size_t trials) {
        if (successes > trials) {
            std::cout << "ERROR -- Successes more than trials:" << successes << " > " << trials << std::endl;
            throw 10;
        } else if (successes == 0) { return 0.; }

        size_t starting_k = successes;
        while(cache.count(bit_concat(starting_k, trials)) == 0 and starting_k < trials + 1) {
            starting_k++;
        }

        if (starting_k == trials + 1) {cache[bit_concat(trials + 1, trials)] = -INFINITY;}

        for (starting_k--; starting_k >= successes; starting_k--) {
            cache[bit_concat(starting_k, trials)] = lmath.log_sum(cache.at(bit_concat(starting_k + 1, trials)), 
                                                                      log_pmf(starting_k, trials));
        }
        return cache.at(bit_concat(successes, trials));
    }

    float BinomialDist::log10_pvalue(size_t successes, size_t trials) {
        return std::log10(std::exp(1)) * log_pvalue(successes, trials);
    }

}
