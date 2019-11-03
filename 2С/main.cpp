#include <iostream>
#include <vector>
#include <iostream>
#include <cmath>

const int alphabet = 256;

class SuffixArrayFounder {
public:
    explicit SuffixArrayFounder(const std::string &s) : size(s.size()), suffix_arrays(s.size()),
                                                        equivalent_clases(s.size(), std::vector<long long>(s.size())) {
        std::vector<long long> counting(alphabet, 0);
        for (long long i = 0; i < s.size(); ++i) {
            ++counting[s[i]];
        }
        for (long long i = 1; i < alphabet; ++i) {
            counting[i] += counting[i - 1];
        }
        for (long long i = s.size() - 1; i >= 0; --i) {
            suffix_arrays[--counting[s[i]]] = i;
        }
        long long amount_of_classes = 1;
        equivalent_clases[0][suffix_arrays[0]] = 0;
        for (long long i = 1; i < s.size(); ++i) {
            if (s[suffix_arrays[i]] != s[suffix_arrays[i - 1]]) { ++amount_of_classes; }
            equivalent_clases[0][suffix_arrays[i]] = amount_of_classes - 1;
        }
        for (long long j = 0; (1 << j) < s.size(); ++j) {
            std::vector<long long> permitation_in_the_second_order(s.size());
            for (long long i = 0; i < s.size(); ++i) {
                permitation_in_the_second_order[i] = suffix_arrays[i] - (1 << j);
                if (permitation_in_the_second_order[i] < 0) {
                    permitation_in_the_second_order[i] += s.size();
                }
            }
            counting = std::vector<long long>(amount_of_classes, 0);
            for (long long i = 0; i < s.size(); ++i) {
                counting[equivalent_clases[j][permitation_in_the_second_order[i]]]++;
            }
            for (long long i = 1; i < amount_of_classes; ++i) {
                counting[i] += counting[i - 1];
            }
            for (long long i = s.size() - 1; i >= 0; --i) {
                suffix_arrays[--counting[equivalent_clases[j][permitation_in_the_second_order[i]]]] = permitation_in_the_second_order[i];
            }
            equivalent_clases[j + 1][0] = 0;
            amount_of_classes = 1;
            for (long long i = 1; i < s.size(); ++i) {
                long long mid1 = (suffix_arrays[i] + (1 << j)) % s.size(), mid2 =
                        (suffix_arrays[i - 1] + (1 << j)) % s.size();
                if ((equivalent_clases[j][suffix_arrays[i]] != equivalent_clases[j][suffix_arrays[i - 1]])
                    || (equivalent_clases[j][mid1] != equivalent_clases[j][mid2])) {
                    amount_of_classes++;
                }
                equivalent_clases[j + 1][suffix_arrays[i]] = amount_of_classes - 1;
            }
        }
    }

    std::vector<long long>& show_suffix_array() { return suffix_arrays; }

    std::vector<std::vector<long long>> & show_equivalent_classes() { return equivalent_clases; }

    long long show_size() const { return size; }

    long long LCP(long long i, long long j) {
        long long ans = 0;
        for (long long k = floor(log2(size)); k >= 0; --k) {
            if (equivalent_clases[k][i] == equivalent_clases[k][j]) {
                ans += (1 << k);
                i += (1 << k);
                j += (1 << k);
            }
        }
        return ans;
    }

    std::vector<long long> LCPs() {
        std::vector<long long> lcps(size - 1);
        lcps[0] = 0;
        for (long long i = 1; i < size - 1; ++i) { lcps[i] = LCP(suffix_arrays[i], suffix_arrays[i + 1]); }
        return lcps;
    }

private:
    long long size;
    std::vector<long long> suffix_arrays;
    std::vector<std::vector<long long>> equivalent_clases;
};

bool Are_not_in_one_string(long long index1, long long index2, long long size_of_first) {
    return ((index1 > size_of_first && index2 < size_of_first) || (index2 > size_of_first && index1 < size_of_first));
}

std::string common_k_stats(const std::string &s, const std::string &t, long long k) {
    std::string str = s + '$' + t + '#';
    std::string ans = "-1";
    SuffixArrayFounder S = SuffixArrayFounder(str);
    std::vector<long long> suffix_array = S.show_suffix_array();
    std::vector<long long> Lcps = S.LCPs();
    bool k_stats_exists = false;
    long long minimum_common_substr = 0;
    for (long long i = 2; i < suffix_array.size() - 1; ++i) {
        if (Are_not_in_one_string(suffix_array[i], suffix_array[i + 1], s.size()) && Lcps[i] != 0) {
            long long cur_substr = std::max(Lcps[i] - minimum_common_substr, (long long) 0);
            if (k > cur_substr) {
                k -= cur_substr;
                minimum_common_substr = Lcps[i];
            } else {
                ans = str.substr(suffix_array[i], k + minimum_common_substr);
                break;
            }
        }
        minimum_common_substr = std::min(minimum_common_substr, Lcps[i]);
    }
    return ans;
}


int main() {
    std::string s, t;
    std::cin >> s >> t;
    long long k;
    std::cin >> k;
    std::cout << common_k_stats(s, t, k);
}