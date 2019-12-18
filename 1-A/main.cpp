#include <iostream>
#include <vector>
#include <string>

void iteration(int &previos_pf, std::vector<int> &prefix, std::string &s, std::string &temp, int i, int &j) {
    j = previos_pf;
    while (j > 0 and s[i] != temp[j]) {
        j = prefix[j - 1];
    }
    if (temp[j] == s[i]) {
        ++j;
    }
}

void find_entries(std::string &s, std::string &temp, std::vector<int> &prefix, std::vector<int> &entries) {
    int size = temp.length();
    int previos_pf = 0;
    for (int i = 1; i < size; ++i) {
        previos_pf = prefix[i - 1];
        int j = 0;
        iteration(previos_pf, prefix, temp, temp, i, j);
        prefix[i] = j;
    }
    previos_pf = 0;
    for (int i = 0; i < s.size(); ++i) {
        int j = 0;
        iteration(previos_pf, prefix, s, temp, i, j);
        if (temp[i] != s[j]) { j = 0; }
        previos_pf = j;
        if (j == size) {
            entries.push_back(1 + i - size);
        }
    }
}

int main() {
    std::string temp, str;
    std::cin >> temp >> str;
    std::vector<int> prefix_for_temp(temp.length());
    prefix_for_temp[0] = 0;
    std::vector<int> enrty;
    find_entries(str, temp, prefix_for_temp, enrty);
    for (auto i : enrty) { std::cout << i << " "; }
}
