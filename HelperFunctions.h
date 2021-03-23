

#ifndef SPEEDPPR_HELPERFUNCTIONS_H
#define SPEEDPPR_HELPERFUNCTIONS_H

#define MSG(...)  { std::cout << #__VA_ARGS__ << ":" << (__VA_ARGS__) << std::endl; }


#include <string>
#include <vector>
#include <cstring>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <chrono>
#include <fstream>
#include "BasicDefinition.h"

extern double getCurrentTime();


inline uint32_t parse_integer(const std::string &_str, size_t &_end) {
    uint32_t rtn = 0;
    for (_end = 0; !isdigit(_str[_end]);) { ++_end; }
    for (; isdigit(_str[_end]); ++_end) {
        rtn *= 10;
        rtn += _str[_end] - '0';
    }
    return rtn;
}

inline uint32_t parse_integer(const std::string &_str, const size_t &_start, size_t &_end) {
    uint32_t rtn = 0;
    for (_end = _start; !isdigit(_str[_end]);) { ++_end; }
    for (; isdigit(_str[_end]); ++_end) {
        rtn *= 10;
        rtn += _str[_end] - '0';
    }
    return rtn;
}

inline double parse_double(const std::string &_str, const size_t &_start, size_t &_end) {
    size_t rtn = 0, scale = 1;
    for (_end = _start; isspace(_str[_end]); ++_end);
    if (!isdigit(_str[_end])) { printf("Error in parsing double.\n"); }
    for (; isdigit(_str[_end]); ++_end) {
        rtn *= 10;
        rtn += _str[_end] - '0';
    }
    for (; isspace(_str[_end]);  ++_end);
    if (_str[_end++] != '.') {
        printf("Error in parsing double. Expecting Dot. \n");
        exit(1);
    }
    for (; isdigit(_str[_end]); ++_end, scale *= 10) {
        rtn *= 10;
        rtn += _str[_end] - '0';
    }
    return rtn / (double) scale;
}


inline void
show(const std::vector<PageRankScoreType> &pi, const unsigned int &_numOfVertices, unsigned int &_num_to_show) {
    if (_num_to_show > pi.size()) {
        _num_to_show = pi.size();
    }
    std::priority_queue<PageRankScoreType,
            std::vector<PageRankScoreType>,
            std::greater<> > top_k_queue;
    for (VertexIdType id = 0; id < _numOfVertices; ++id) {
        if (top_k_queue.size() < _num_to_show) {
            top_k_queue.push(pi[id]);
        } else if (top_k_queue.top() < pi[id]) {
            top_k_queue.pop();
            top_k_queue.push(pi[id]);
        }
    }
    printf("Top Page Rank Scores.\n");
    std::vector<PageRankScoreType> top_k_list;
    while (top_k_queue.empty() == false) {
        top_k_list.emplace_back(top_k_queue.top());
        top_k_queue.pop();
    }
    while (top_k_list.empty() == false) {
        printf("%.13f\n", top_k_list.back());
        top_k_list.pop_back();
    }
}


template<class FLOAT_TYPE>
inline void
save_answer(const std::vector<FLOAT_TYPE> &_means, const VertexIdType &_numOfVertices, const std::string &_file_name) {
    std::ofstream file(_file_name);
    if (file.is_open()) {
        std::vector<IdScorePair<FLOAT_TYPE>> pairs(_numOfVertices, 0);
        for (VertexIdType id = 0; id < _numOfVertices; ++id) {
            pairs[id] = {id, _means[id]};
        }
        std::sort(pairs.begin(), pairs.end(), IdScorePairComparatorGreater<FLOAT_TYPE>());
        file.precision(std::numeric_limits<FLOAT_TYPE>::max_digits10);
        for (const auto &pair : pairs) {
            if (pair.score > 0) {
                file << pair.id << "\t" << pair.score << "\n";
            }
        }
        file.close();
    } else {
        printf("ERROR IN " __FILE__ " LINE %u\n", __LINE__);
        printf("FILE NOT EXISTS.\n%s\n", _file_name.c_str());
        exit(1);
    }
}


struct Param {
    std::string graph_file;
    std::string query_file;
    std::string answer_folder;
    std::string index_file;
    std::string meta_file;
    std::string graph_binary_file;
    std::string algorithm = "PowItr";
    std::string output_folder;
    std::string estimation_folder;
    double epsilon = 0.5;
    double alpha = 0.2;
    double l1_error = 0.000'000'01;
    unsigned int num_top_k = 1;
    unsigned int query_size = 0;
    bool with_idx = false;
    bool is_top_k = false;
    bool is_undirected_graph = false;
    bool specified_l1_error = false;
    bool output_estimations = false;
};

extern Param param;

extern Param parseArgs(int nargs, char **args);

template<class T>
inline void show_vector(const std::string &_header, const std::vector<T> &_vec) {
    if (_vec.empty()) {
        std::cout << "Empty Vector." << std::endl;
    } else {
        std::cout << std::endl << _header;
        bool identical = true;
        const T &elem = _vec.front();
        std::for_each(_vec.begin(), _vec.end(), [&](const T &e) { identical &= (e == elem); });
        if (identical) {
            std::cout << "\tSize of the Vector: " << _vec.size() << "\t Value of Each Element: " << elem;
        } else {
            std::cout << std::endl;
            std::copy(begin(_vec), end(_vec), std::ostream_iterator<T>(std::cout, "\t"));
        }
        std::cout << std::endl;
    }
}

#endif //SPEEDPPR_HELPERFUNCTIONS_H
