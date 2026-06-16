#ifndef HAP_LIKELIHOOD_H
#define HAP_LIKELIHOOD_H

#include <string>
#include <cmath>
#include <algorithm>
#include <vector>
#include "util.h"

namespace HapLL {
// --- Core Model Parameters ---
    constexpr double MU_INS = 1.0e-9;
    constexpr double MU_DEL = 2.5e-9;
    constexpr double BETA   = 1.15;
    constexpr double GAMMA  = 0.1;
    constexpr double ALPHA  = 0.9;
    constexpr double P_BG   = 0.25;

    // --- Expanded Container for Indel Repeat Analysis ---
    struct RepeatInfo {
        std::string unit = "";
        int count_in_indel = 0;    // Number of units inside the deleted/inserted sequence
        int count_in_flanking = 0; // Number of units found in left + right flanking regions
        int total_template_count = 0; // Total unit copies present in the original template allele
        int total_length = 1;      // L parameter used for context scaling
    };

    /**
     * @brief Computes Levenshtein distance normalized similarity (0.0 to 1.0)
     */
    inline double calculate_similarity(const std::string& s1, const std::string& s2) {
        if (s1.empty() && s2.empty()) return 1.0;
        if (s1.empty() || s2.empty()) return 0.0;
        
        int m = s1.length(), n = s2.length();
        std::vector<std::vector<int>> dp(m + 1, std::vector<int>(n + 1, 0));
        
        for (int i = 0; i <= m; ++i) dp[i][0] = i;
        for (int j = 0; j <= n; ++j) dp[0][j] = j;
        
        for (int i = 1; i <= m; ++i) {
            for (int j = 1; j <= n; ++j) {
                int cost = (s1[i - 1] == s2[j - 1]) ? 0 : 1;
                dp[i][j] = std::min({dp[i - 1][j] + 1, dp[i][j - 1] + 1, dp[i - 1][j - 1] + cost});
            }
        }
        return std::max(0.0, 1.0 - (static_cast<double>(dp[m][n]) / std::max(m, n)));
    }

    /**
     * @brief Finds the smallest sub-repeat unit within the indel sequence itself.
     * Example: "ATAATA" -> returns "ATA"
     */
    inline std::string find_indel_base_unit(const std::string& indel_seq) {
        int n = indel_seq.length();
        for (int m = 1; m <= n; ++m) {
            if (n % m == 0) {
                std::string sub = indel_seq.substr(0, m);
                bool is_composite = true;
                for (int i = m; i < n; i += m) {
                    if (indel_seq.substr(i, m) != sub) {
                        is_composite = false;
                        break;
                    }
                }
                if (is_composite) return sub;
            }
        }
        return indel_seq;
    }

    /**
     * @brief Counts continuous backward repeats of the unit in the left flank
     */
    inline int count_left_repeats(const std::string& lt_seq, const std::string& unit) {
        int m = unit.length();
        int idx = lt_seq.length();
        int count = 0;
        while (idx >= m) {
            if (lt_seq.substr(idx - m, m) == unit) {
                count++;
                idx -= m;
            } else {
                break;
            }
        }
        return count;
    }

    /**
     * @brief Counts continuous forward repeats of the unit in the right flank
     */
    inline int count_right_repeats(const std::string& rt_seq, const std::string& unit) {
        int m = unit.length();
        int n = rt_seq.length();
        int idx = 0;
        int count = 0;
        while (idx + m <= n) {
            if (rt_seq.substr(idx, m) == unit) {
                count++;
                idx += m;
            } else {
                break;
            }
        }
        return count;
    }

    inline double calc_log_likelihood(bool is_insertion, int k, int L, 
                                      const std::string& ins_seq, const std::string& ref_prev) {
        double log_mu = std::log(is_insertion ? MU_INS : MU_DEL);
        double log_phi = BETA * std::max(0, L - 1);
        double log_psi = std::log(1.0 - GAMMA) + (k - 1) * std::log(GAMMA);
        
        double log_omega = 0.0;
        if (is_insertion) {
            double sim = calculate_similarity(ins_seq, ref_prev);
            double prob_dup = ALPHA * sim;
            double prob_rand = (1.0 - ALPHA) * std::pow(P_BG, k);
            log_omega = std::log(prob_dup + prob_rand);
        }
        return log_mu + log_phi + log_psi + log_omega;
    }

    /**
     * @brief Main evaluation pipeline factoring in internal repeat structural breakdown
     */
    inline double evaluate_variant(const Variant& v, RepeatInfo& out_repeat) {
        if (v.ref.length() == v.alt.length()) {
            return -1000.0; // Bypass SNVs
        }

        bool is_insertion = v.alt.length() > v.ref.length();
        int k = std::abs(static_cast<int>(v.ref.length()) - static_cast<int>(v.alt.length()));

        // 1. Isolate the exact altered core sequence segment
        std::string indel_seq = is_insertion ? v.alt.substr(v.ref.length()) 
                                             : v.ref.substr(v.alt.length());

        // 2. Resolve internal pattern composition (e.g., "ATAATA" -> "ATA")
        out_repeat.unit = find_indel_base_unit(indel_seq);
        out_repeat.count_in_indel = indel_seq.length() / out_repeat.unit.length();

        // 3. Quantify how often this unit blankets the structural neighborhood
        int left_flank_count = count_left_repeats(v.sample_lt_seq, out_repeat.unit);
        int right_flank_count = count_right_repeats(v.sample_rt_seq, out_repeat.unit);
        out_repeat.count_in_flanking = left_flank_count + right_flank_count;

        // 4. Calculate total copies on the unmodified template strand (L derivation)
        if (is_insertion) {
            // For insertions, the unmutated template only contained the flanking units
            out_repeat.total_template_count = out_repeat.count_in_flanking;
        } else {
            // For deletions, the unmutated template contained flanks AND the deleted chunk
            out_repeat.total_template_count = out_repeat.count_in_flanking + out_repeat.count_in_indel;
        }
        
        out_repeat.total_length = std::max(1, out_repeat.total_template_count * static_cast<int>(out_repeat.unit.length()));

        // 5. Setup sequence matching profile for insertion models
        std::string ins_seq = "";
        std::string ref_prev = "";
        if (is_insertion) {
            ins_seq = indel_seq;
            if (v.sample_lt_seq.length() >= static_cast<size_t>(k)) {
                ref_prev = v.sample_lt_seq.substr(v.sample_lt_seq.length() - k);
            } else {
                ref_prev = v.sample_lt_seq;
            }
        }

        return calc_log_likelihood(is_insertion, k, out_repeat.total_length, ins_seq, ref_prev);
    }
}

#endif 
