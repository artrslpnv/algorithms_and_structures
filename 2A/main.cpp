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

long long CountDifferentSubstr(std::string &s) {
    SuffixArrayFounder S = SuffixArrayFounder(s);
    std::vector <long long > lcps=S.Kasai_algorithm(s);
    long long ans = 0;
    for (long long i = 1; i < S.show_size(); ++i) { // c 1 потому что игнорируем 0 элемент
        ans += (S.show_size() - 1 - S.show_suffix_array()[i]);
    }
    for (long long i=0;i<s.size()-1;++i){
        ans-=lcps[i];
    }
    return ans;
}


int main() {
    std::string s;
    std::cin >> s;
    s=s+'$';
    std::cout << CountDifferentSubstr(s);
}