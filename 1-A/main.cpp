#include <iostream>
#include <vector>
#include <string>

void prefix_function(std::string &s, std::string &temp, std::vector<int> &prefix, std::vector<int> &entries) {
    int size = temp.length();
    for (int i = 1; i < size; ++i) {
        int index = prefix[i - 1];
        while (index > 0 && temp[i] != temp[index]) {
            index = prefix[index - 1];
        }
        if (temp[i] == temp[index]) { ++index; }
        prefix[i] = index;
    }
    int previos_pf = 0;

    for (int i = 0; i < s.size(); ++i) {
        int j = previos_pf;
        while (j > 0 and s[i] != temp[j]) {
            j = prefix[j - 1];
        }
        if (temp[j] != s[i]) {
            j = 0;
        } else { ++j; }
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
    prefix_function(str, temp, prefix_for_temp, enrty);
    for (auto i : enrty) { std::cout << i << " "; }
}