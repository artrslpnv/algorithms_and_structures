#include <iostream>
#include <vector>
#include <iostream>
#include <cmath>

const long alphabet = 256;

class SuffixArrayFounder {
public:
    explicit SuffixArrayFounder(const std::string &s) : size(s.size()), suffix_arrays(s.size()),
                                                        equivalent_clases(s.size(), std::vector<long>(s.size())) {
        std::vector<long> counting(alphabet, 0);
        for (long i = 0; i < s.size(); ++i) {
            ++counting[s[i]];
        }
        for (long i = 1; i < alphabet; ++i) {
            counting[i] += counting[i - 1];
        }
        for (long i = s.size() - 1; i >= 0; --i) {
            suffix_arrays[--counting[s[i]]] = i;
        }
        long amount_of_classes = 1;
        equivalent_clases[0][suffix_arrays[0]] = 0;
        for (long i = 1; i < s.size(); ++i) {
            if (s[suffix_arrays[i]] != s[suffix_arrays[i - 1]]) { ++amount_of_classes; }
            equivalent_clases[0][suffix_arrays[i]] = amount_of_classes - 1;
        }
        for (long j = 0; (1 << j) < s.size(); ++j) {
            std::vector<long> permitation_in_the_second_order(s.size());
            for (long i = 0; i < s.size(); ++i) {
                permitation_in_the_second_order[i] = suffix_arrays[i] - (1 << j);
                if (permitation_in_the_second_order[i] < 0) {
                    permitation_in_the_second_order[i] += s.size();
                }
            }
            counting = std::vector<long>(amount_of_classes, 0);
            for (long i = 0; i < s.size(); ++i) {
                counting[equivalent_clases[j][permitation_in_the_second_order[i]]]++;
            }
            for (long i = 1; i < amount_of_classes; ++i) {
                counting[i] += counting[i - 1];
            }
            for (long i = s.size() - 1; i >= 0; --i) {
                suffix_arrays[--counting[equivalent_clases[j][permitation_in_the_second_order[i]]]] = permitation_in_the_second_order[i];
            }
            equivalent_clases[j + 1][0] = 0;
            amount_of_classes = 1;
            for (long i = 1; i < s.size(); ++i) {
                long mid1 = (suffix_arrays[i] + (1 << j)) % s.size(), mid2 =
                        (suffix_arrays[i - 1] + (1 << j)) % s.size();
                if ((equivalent_clases[j][suffix_arrays[i]] != equivalent_clases[j][suffix_arrays[i - 1]])
                    || (equivalent_clases[j][mid1] != equivalent_clases[j][mid2])) {
                    amount_of_classes++;
                }
                equivalent_clases[j + 1][suffix_arrays[i]] = amount_of_classes - 1;
            }
        }
    }

    std::vector<long> show_suffix_array() const { return suffix_arrays; }

    std::vector<std::vector<long>> show_equivalent_classes() const { return equivalent_clases; }

    long show_size() const { return size; }

    long LCP(long i, long j) {
        long ans = 0;
        for (long k = floor(log2(size)); k >= 0; --k) {
            if (equivalent_clases[k][i] == equivalent_clases[k][j]) {
                ans += (1 << k);
                i += (1 << k);
                j += (1 << k);
            }
        }
        return ans;
    }

private:
    long size;
    std::vector<long> suffix_arrays;
    std::vector<std::vector<long>> equivalent_clases;
};

long long CountDifferentSubstr(std::string &s) {
    SuffixArrayFounder S = SuffixArrayFounder(s + "$");
    long long ans = 0;
    for (long i = 1; i < S.show_size() - 1; ++i) { // c 1 потому что игнорируем 0 элемент
        ans += (S.show_size() - 1 - S.show_suffix_array()[i] -
                S.LCP(S.show_suffix_array()[i], S.show_suffix_array()[i + 1]));
    }
    ans += S.show_size() - 1 - S.show_suffix_array()[s.size()];
    return ans;
}


int main() {
    std::string s;
    std::cin >> s;
    std::cout << CountDifferentSubstr(s);
}