
#include <iostream>
#include <vector>
#include <string>

std::string lexographic_minimum(std::vector<int> &prefix) {
    if (!(prefix.size())) {
        return "";
    }
    std::string str = "a";
    for (int i = 1; i < prefix.size(); ++i) {
        if (prefix[i] == 0) {
            std::vector<bool> alphabet(26);
            alphabet.assign(26, false);
            int previos_prefix_func_value = prefix[i - 1];
            while (previos_prefix_func_value != 0) {
                alphabet[str[previos_prefix_func_value] - 'a'] = true;
                previos_prefix_func_value = prefix[previos_prefix_func_value - 1];
            }
            char character_to_put = 'a' + 1;
            while (true) {
                if (!alphabet[character_to_put - 'a']) {
                    break;
                }
                ++character_to_put;
            }
            str = str + character_to_put;
        } else {
            int symbol_index = prefix[i] - 1;
            str += str[symbol_index];
        }
    }
    return str;
}

int main() {
    std::vector<int> prefix_function;
    int x = 0;
    while (std::cin >> x) {
        prefix_function.push_back(x);
    }
    for (auto u : lexographic_minimum(prefix_function)) { std::cout << u; }
}
