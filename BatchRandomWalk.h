

#ifndef SPEEDPPR_BATCHRANDOMWALK_H
#define SPEEDPPR_BATCHRANDOMWALK_H


#include <vector>
#include <numeric>
#include <stack>
#include <random>
#include <queue>
#include "BasicDefinition.h"
#include "Graph.h"
#include "MyQueue.h"
#include "MyRandom.h"

template<class FLOAT_TYPE>
class Alias {
private:
    const unsigned int size;
    const std::vector<VertexIdType> &first;
    std::vector<VertexIdType> second;
    std::vector<FLOAT_TYPE> probability;
public:
    Alias(const std::vector<VertexIdType> &_active_ids, const std::vector<FLOAT_TYPE> &_active_residuals) :
            size(_active_ids.size()),
            first(_active_ids),
            second(_active_ids.size(), 0),
            probability(_active_residuals) {
        const FLOAT_TYPE sum = std::accumulate(_active_residuals.begin(), _active_residuals.end(), 0.0);
        std::stack<VertexIdType, std::vector<VertexIdType >> small;
        std::stack<VertexIdType, std::vector<VertexIdType >> big;
        const FLOAT_TYPE size_over_sum = size / sum;
        for (VertexIdType id = 0; id < size; ++id) {
            probability[id] *= size_over_sum;
            if (probability[id] > 1) {
                big.push(id);
            } else {
                small.push(id);
            }
        }
        while (!small.empty() && !big.empty()) {
            const VertexIdType small_id = small.top();
            small.pop();
            const VertexIdType big_id = big.top();
            second[small_id] = first[big_id];
            probability[big_id] -= (1 - probability[small_id]);
            if (probability[big_id] < 1) {
                small.push(big_id);
                big.pop();
            }
        }
    }


    inline VertexIdType generate_random_id() const {
        const unsigned int bucket_id = MinimalStandardGenerator::uniform_int(size);
//        const unsigned int bucket_id = MyRandom::rand_int(0, size - 1);
        // First or second
        return MinimalStandardGenerator::bias_coin_is_head(probability[bucket_id]) ? first[bucket_id]
                                                                                   : second[bucket_id];
    }

};


class WalkCache {
private:
    const Graph &graph;
    std::vector<VertexIdType> walks;
    std::vector<VertexIdType> start_indices;
public:
    static const uint32_t num_zero_walks;

public:

    explicit WalkCache(const Graph &_graph) :
            graph(_graph),
            walks(_graph.getNumOfVertices() * num_zero_walks
                  + _graph.getNumOfEdges() + _graph.get_num_dead_end(), 0),
            start_indices(_graph.getNumOfVertices(), 0) {
    }

    void generate() {
        const VertexIdType num_vertices = graph.getNumOfVertices();
        for (VertexIdType sid = 0, index = 0; sid < num_vertices; ++sid) {
            if (sid % 500000 == 0) { std::cout << sid << " vertices processed.\n"; }
            start_indices[sid] = index;
            for (uint32_t j = 0; j < num_zero_walks; ++j) {
                VertexIdType current_id = sid;
                while (SFMT64::bias_coin_is_tail(graph.get_alpha())) {
                    const VertexIdType &idx_start = graph.get_neighbor_list_start_pos(current_id);
                    const VertexIdType &idx_end = graph.get_neighbor_list_start_pos(current_id + 1);
                    const VertexIdType degree = idx_end - idx_start;
                    const VertexIdType shift = SFMT64::uniform_int(degree);
                    const VertexIdType nid = graph.getOutNeighbor(idx_start + shift);
                    current_id = nid;
                }
                walks[index++] = current_id;
            }
            const VertexIdType &sid_idx_start = graph.get_neighbor_list_start_pos(sid);
            const VertexIdType &sid_idx_end = graph.get_neighbor_list_start_pos(sid + 1);
            const VertexIdType sid_degree = sid_idx_end - sid_idx_start;
            for (uint32_t j = 0; j < sid_degree; ++j) {
                const VertexIdType sid_shift = SFMT64::uniform_int(sid_degree);
                VertexIdType current_id = graph.getOutNeighbor(sid_idx_start + sid_shift);
                while (SFMT64::bias_coin_is_tail(graph.get_alpha())) {
                    // continues to walk
                    const VertexIdType &idx_start = graph.get_neighbor_list_start_pos(current_id);
                    const VertexIdType &idx_end = graph.get_neighbor_list_start_pos(current_id + 1);
                    const VertexIdType degree = idx_end - idx_start;
                    const VertexIdType shift = SFMT64::uniform_int(degree);
                    const VertexIdType nid = graph.getOutNeighbor(idx_start + shift);
                    current_id = nid;
                }
                walks[index++] = current_id;
            }
        }
    }

    void save(const std::string &_filename) const {
        const auto start = getCurrentTime();
        if (std::FILE *f = std::fopen(_filename.c_str(), "wb")) {
            std::fwrite(walks.data(), sizeof walks[0], walks.size(), f);
            std::fclose(f);
        } else {
            printf("WalkCache::save; File Not Exists.\n");
        }
        const auto end = getCurrentTime();
        printf("Time Used For Saving : %.2f\n", end - start);
    }

    void load(const std::string &_filename) {
        const auto start = getCurrentTime();
        walks.clear();
        //TODO check here
        walks.resize(graph.getNumOfVertices() * num_zero_walks
                     + graph.getNumOfEdges() + graph.get_num_dead_end(), 0);
        assert(walks.size() ==
               graph.getNumOfVertices() * num_zero_walks + graph.get_neighbor_list_start_pos(graph.getNumOfVertices()));
        if (std::FILE *f = std::fopen(_filename.c_str(), "rb")) {
            MSG(walks.size())
            size_t rtn = std::fread(walks.data(), sizeof walks[0], walks.size(), f);
            printf("Returned Value of fread: %zu\n", rtn);
            std::fclose(f);
            start_indices.clear();
            start_indices.resize(graph.getNumOfVertices(), 0);
            for (VertexIdType prev_id = 0, id = 1; id < graph.getNumOfVertices(); ++prev_id, ++id) {
                start_indices[id] =
                        start_indices[prev_id]
                        + std::max(1u, graph.original_out_degree(prev_id))
                        + num_zero_walks;
            }
        } else {
            printf("WalkCache::load; File Not Exists.\n");
            exit(1);
        }
        const auto end = getCurrentTime();
        printf("Time Used For Loading Cache : %.2f\n", end - start);
    }

    void show() {
        show_vector("The cache walk vector", walks);
    }

    inline const VertexIdType &get_zero_hop_start_index(const VertexIdType &_vid) const {
        assert(_vid < graph.getNumOfVertices());
        return start_indices[_vid];
    }

    inline VertexIdType get_one_hop_start_index(const VertexIdType &_vid) const {
#ifdef DEBUG_MODE
//        if (_vid >= graph.getNumOfVertices()) {
//            printf("Error in " __FILE__ " line %d\n", __LINE__);
//            exit(1);
//        }
#endif
        assert(_vid < graph.getNumOfVertices());
        return start_indices[_vid] + num_zero_walks;
    }

    inline const VertexIdType &get_walk(const VertexIdType &_index) const {
        assert(_index < walks.size());
        return walks[_index];
    }
};

#endif //SPEEDPPR_BATCHRANDOMWALK_H
