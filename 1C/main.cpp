#include <iostream>
#include <vector>


const int chars_amount = 26;

class AhoCorasickAutomat {
public:
    explicit AhoCorasickAutomat(std::string &pattern_template);

    std::vector<int> find_entries(std::string &text);
private:
    struct Node;

    int get_vertice_connected_by_suffix_link(int index);

    int find_link(int index, char character);

    void adding_of_pattern( std::pair<int, int> &positions_of_smaller_patterns, int pattern_index);

    void position_of_subputterns_finding( std::string &pattern_template);

    std::vector<Node> Bor;
    std::vector<std::pair<int, int>> subpattern_positions_;
    std::string _pattern_template;
};

struct AhoCorasickAutomat::Node {
    explicit Node(int parent, char char_was_added_with);
    explicit Node();
    char char_was_added_with;
    int parent;
    int suf_link;
    bool is_terminal;
    std:: vector <int> edges;
    std :: vector <int> automat_step;
    std::vector<int> connected_patterns;
};

AhoCorasickAutomat::Node::Node()
{edges.assign(chars_amount,-1);
    char_was_added_with = 0;
    automat_step.assign(chars_amount,-1);
    parent = -1;
    suf_link = -1;

}

AhoCorasickAutomat::Node::Node(int _parent, char _char_was_added_with)
{edges.assign(chars_amount,-1);
    char_was_added_with = _char_was_added_with;
    automat_step.assign(chars_amount,-1);
    parent = _parent;
    suf_link = -1;

}

void AhoCorasickAutomat::adding_of_pattern( std::pair<int, int> &positions_of_smaller_patterns, int pattern_index) {
    int number_of_current_Node = 0;
    for (int i = positions_of_smaller_patterns.first; i <= positions_of_smaller_patterns.second; i++) {
        char character = _pattern_template[i] - 'a';
        if (Bor[number_of_current_Node].edges[character] == -1) {
            Bor.push_back(Node(number_of_current_Node, character));
            Bor[number_of_current_Node].edges[character] = Bor.size() - 1;
        }
        number_of_current_Node = Bor[number_of_current_Node].edges[character];
    }
    Bor[number_of_current_Node].is_terminal = true;
    Bor[number_of_current_Node].connected_patterns.push_back(pattern_index);
}



AhoCorasickAutomat::AhoCorasickAutomat( std::string &pattern_template)
{Bor.assign(1, Node());
    Bor[0].suf_link = 0;
    _pattern_template=pattern_template;
    position_of_subputterns_finding(pattern_template);
    for (int i = 0; i < subpattern_positions_.size(); i++) {
        adding_of_pattern(subpattern_positions_[i], i);
    }
}


void AhoCorasickAutomat::position_of_subputterns_finding( std::string &pattern_template) {
    std::pair<int, int> current_subpattern_pos;
    if (pattern_template[0]!='?') {
        current_subpattern_pos.first = 0;
    }
    if (_pattern_template[1] == '?') {if( _pattern_template[0]!='?'){
            current_subpattern_pos.second = 0;
            subpattern_positions_.push_back(current_subpattern_pos);
        }}
    for (int i = 1; i < _pattern_template.length() - 1; i++) {
        if (pattern_template[i - 1] == '?' && pattern_template[i]!='?') {
            current_subpattern_pos.first = i;
        }
        if (pattern_template[i + 1] == '?' && pattern_template[i]!='?') {
            current_subpattern_pos.second = i;
            subpattern_positions_.push_back(current_subpattern_pos);
        }
    }
    if (_pattern_template[pattern_template.length() - 2] == '?' ) {
        if(pattern_template[pattern_template.length() - 1]!='?'){
            current_subpattern_pos.first = pattern_template.length() - 1;
        }
    }
    if (_pattern_template[_pattern_template.length() - 1]!='?') {
        current_subpattern_pos.second = pattern_template.length() - 1;
        subpattern_positions_.push_back(current_subpattern_pos);
    }
}
int AhoCorasickAutomat::get_vertice_connected_by_suffix_link(int index) {
    if (Bor[index].suf_link == -1) {
        if (Bor[index].parent != 0) {
            Bor[index].suf_link = find_link(get_vertice_connected_by_suffix_link(Bor[index].parent), Bor[index].char_was_added_with);
        }
        else  { Bor[index].suf_link = 0; }
    }
    return Bor[index].suf_link;
}

int AhoCorasickAutomat::find_link(int index, char character) {
    if (Bor[index].automat_step[character] == -1) {
        if (Bor[index].edges[character] == -1 && index!=0){Bor[index].automat_step[character] = find_link(get_vertice_connected_by_suffix_link(index), character);}
        else if (Bor[index].edges[character] != -1) {
            Bor[index].automat_step[character] = Bor[index].edges[character];
        }
        else {
            Bor[index].automat_step[character] = 0;
        }
    }
    return Bor[index].automat_step[character];
}

std::vector<int> AhoCorasickAutomat::find_entries( std::string &text) {
    int vertice = 0;
    std::vector<int> entries(text.length());
    std::vector<int> answer;
    for (int i = 0; i < text.length(); i++) {
        vertice = find_link(vertice, text[i] - 'a');
        int u = vertice;
        if (Bor[u].is_terminal) {
            for (int j = 0; j < Bor[u].connected_patterns.size(); j++) {
                int begin_ind = i + subpattern_positions_[Bor[u].connected_patterns[j]].first;
                begin_ind -= subpattern_positions_[Bor[u].connected_patterns[j]].second;
                if ((begin_ind >= subpattern_positions_[Bor[u].connected_patterns[j]].first)
                    && (begin_ind - subpattern_positions_[Bor[u].connected_patterns[j]].first +
                        _pattern_template.length() < text.length() + 1)) {
                    entries[begin_ind - subpattern_positions_[Bor[u].connected_patterns[j]].first]++; }}}
        u = get_vertice_connected_by_suffix_link(u);
        while (u != 0) {
            if (Bor[u].is_terminal) {
                for (int j = 0; j < Bor[u].connected_patterns.size(); j++) {
                    int begin_ind = i + subpattern_positions_[Bor[u].connected_patterns[j]].first;
                    begin_ind-=subpattern_positions_[Bor[u].connected_patterns[j]].second;
                    if ((begin_ind >= subpattern_positions_[Bor[u].connected_patterns[j]].first)
                        && (begin_ind - subpattern_positions_[Bor[u].connected_patterns[j]].first +
                            _pattern_template.length()< text.length()+ 1)) {
                        entries[begin_ind - subpattern_positions_[Bor[u].connected_patterns[j]].first]++;
                    }}}u = get_vertice_connected_by_suffix_link(u);
        }
    }
    for (int i = 0; i < entries.size(); ++i) {
        if (entries[i] == subpattern_positions_.size()) {
            answer.push_back(i);
        }
    }
    return answer;
}

int main() {
    std::string pattern, text;
    std::cin >> pattern >> text;
    AhoCorasickAutomat AhoCorasickAutomat(pattern);
    std::vector<int> entries = AhoCorasickAutomat.find_entries(text);
    for (int i = 0;i<entries.size() ;++i) {
        std::cout << entries[i] << " ";
    }
    return 0;
}
