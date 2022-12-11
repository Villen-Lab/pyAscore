#include <iostream>
#include <algorithm>
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
        while((cache.count(bit_concat(starting_k, trials)) == 0) && (starting_k < trials + 1)) {
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

    /////////////////////////////////
    // PowerSetSum Implementations //
    /////////////////////////////////

    PowerSetSum::PowerSetSum() {
       this->max_depth = 0;
       this->pos = 0;
       this->sums.push_back(0.);
    }

    PowerSetSum::PowerSetSum(const std::vector<float>& target, size_t max_depth) {
        reset(target, max_depth);
    }

    void PowerSetSum::initializeSums(const std::vector<float>& target, size_t start, size_t depth) {
       float base_val = sums.back();
       for (; start < target.size(); start++) {
           sums.push_back(base_val + target.at(start));
           if (start < target.size() - 1 && depth < max_depth - 1) {       
               initializeSums(target, start + 1, depth + 1);
           }
       }
    }

    void PowerSetSum::deduplicateSums() {
        std::sort(sums.begin(), sums.end());

        auto front_iter = sums.begin();
        auto back_iter = sums.begin();
        size_t final_size = 1;
        for (front_iter++; front_iter != sums.end(); front_iter++) {
           if (*front_iter != *back_iter) {
               final_size++;
               back_iter++;
               *back_iter = *front_iter;
           } 
        }
        sums.resize(final_size);
    }
           
    void PowerSetSum::reset() {
        this->pos = 0;
    }

    void PowerSetSum::reset(const std::vector<float>& target, size_t max_depth) {
       if (target.size() < max_depth) {
           this->max_depth = target.size();
       } else {
           this->max_depth = max_depth;
       }

       sums.clear();
       this->pos = 0;
       this->sums.push_back(0.);
       initializeSums(target, 0, 0);
       deduplicateSums();    
    }

    size_t PowerSetSum::getPos() {
        return pos;
    }

    bool PowerSetSum::hasNext() {
        return pos < sums.size() - 1;
    }

    void PowerSetSum::next() {
        if (!hasNext()) {
             throw 40;
        }
        pos++;
    }

    float PowerSetSum::getSum() {
        return sums.at(pos);
    }
}
