#include <iostream>
#include <vector>
#include <iostream>
#include <cmath>

const int alphabet = 256;

class SuffixArrayFounder {
public:
    explicit SuffixArrayFounder(const std::string &s) : size(s.size()), suffix_arrays(s.size()),
                                                        equivalent_clases(s.size()) {
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
        equivalent_clases[suffix_arrays[0]] = 0;
        for (long long i = 1; i < s.size(); ++i) {
            if (s[suffix_arrays[i]] != s[suffix_arrays[i - 1]]) { ++amount_of_classes; }
            equivalent_clases[suffix_arrays[i]] = amount_of_classes - 1;
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
                counting[equivalent_clases[permitation_in_the_second_order[i]]]++;
            }
            for (long long i = 1; i < amount_of_classes; ++i) {
                counting[i] += counting[i - 1];
            }
            for (long long i = s.size() - 1; i >= 0; --i) {
                suffix_arrays[--counting[equivalent_clases[permitation_in_the_second_order[i]]]] = permitation_in_the_second_order[i];
            }
            std:: vector <long long> new_equivalent_clases(s.size());
            new_equivalent_clases[0] = 0;
            amount_of_classes = 1;
            for (long long i = 1; i < s.size(); ++i) {
                long long mid1 = (suffix_arrays[i] + (1 << j)) % s.size(), mid2 =
                        (suffix_arrays[i - 1] + (1 << j)) % s.size();
                if ((equivalent_clases[suffix_arrays[i]] != equivalent_clases[suffix_arrays[i - 1]])
                    || (equivalent_clases[mid1] != equivalent_clases[mid2])) {
                    amount_of_classes++;
                }
                new_equivalent_clases[suffix_arrays[i]] = amount_of_classes - 1;
            }
            equivalent_clases = new_equivalent_clases;
        }
    }

    std::vector<long long>& show_suffix_array()  { return suffix_arrays; }

    std::vector<long long> &show_equivalent_classes(){ return equivalent_clases; }

    long long show_size() const { return size; }

    std::vector<long long> Kasai_algorithm(std::string &s) {
        std::vector<long long > pos(size);//inverted suf
        std::vector<long long > lcp(size);
        for (long long i = 0; i < size; i++) {
            pos[suffix_arrays[i]] = i;
        }
        long long k = 0;
        for (long long i = 0; i < size; ++i) {
            if (k > 0) {
                k--;
            }
            if (pos[i] == size - 1) {
                lcp[size - 1] = -1;
                k = 0;
                continue;
            } else {
                int j = suffix_arrays[pos[i] + 1];
                while (std::max(i + k, j + k) < size && s[i+k]==s[j+k] ){
                    k++;
                }
                lcp[pos[i]]=k;
            }
        }
        return lcp;
    }

private:
    long long size;
    std::vector<long long> suffix_arrays;
    std::vector <long long> equivalent_clases;
};

bool Are_not_in_one_string(long long index1, long long index2, long long size_of_first) {
    return ((index1 > size_of_first && index2 < size_of_first) || (index2 > size_of_first && index1 < size_of_first));
}

std::string common_k_stats(const std::string &s, const std::string &t, long long k) {
    std::string str = s + '$' + t + '#';
    std::string ans = "-1";
    SuffixArrayFounder S = SuffixArrayFounder(str);
    std::vector<long long> suffix_array = S.show_suffix_array();
    std::vector<long long> Lcps = S.Kasai_algorithm(str);
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