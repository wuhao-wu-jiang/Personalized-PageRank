

#ifndef SPEEDPPR_GRAPH_H
#define SPEEDPPR_GRAPH_H

#ifndef USE_REVERSE_POS
#define USE_REVERSE_POS
#endif

#include <string>
#include <vector>
#include <algorithm>
#include <cassert>
#include <fstream>
#include "BasicDefinition.h"
#include "HelperFunctions.h"


class Graph {

private:
    VertexIdType numOfVertices = 0;
    EdgeSizeType numOfEdges = 0;
    VertexIdType num_deadend_vertices = 0;
    VertexIdType sid = 0;
    VertexIdType dummy_id = 0;
    PageRankScoreType alpha = 0.2;
    EdgeSizeType max_size_edge_list = 0;
    std::vector<VertexIdType> out_degrees;
    std::vector<VertexIdType> in_degrees;
    std::vector<VertexIdType> start_pos_in_out_neighbor_lists;
    std::vector<VertexIdType> start_pos_in_appearance_pos_lists;
    std::vector<VertexIdType> out_neighbors_lists;
    std::vector<VertexIdType> appearance_pos_lists;
    std::vector<VertexIdType> deadend_vertices;

public:

    inline size_t get_num_dead_end() const {
        return deadend_vertices.size();
    }

    inline void set_dummy_out_degree_zero() {
        out_degrees[dummy_id] = 0;
        start_pos_in_out_neighbor_lists[dummy_id + 1] = start_pos_in_out_neighbor_lists[dummy_id];
    }

    inline void set_dummy_neighbor(const VertexIdType &_id) {
        out_degrees[dummy_id] = 1;
        start_pos_in_out_neighbor_lists[dummy_id + 1] = start_pos_in_out_neighbor_lists[dummy_id] + 1;
        out_neighbors_lists[start_pos_in_out_neighbor_lists[dummy_id]] = _id;
    }

    inline void reset_set_dummy_neighbor() {
        out_degrees[dummy_id] = 0;
        out_neighbors_lists[start_pos_in_out_neighbor_lists[dummy_id]] = dummy_id;
        set_dummy_out_degree_zero();
    }

    inline const VertexIdType &get_dummy_id() const {
        return dummy_id;
    }

    inline const VertexIdType &get_sid() const {
        return sid;
    }

    inline const PageRankScoreType &get_alpha() const {
        return alpha;
    }

    inline void set_alpha(const PageRankScoreType _alpha = 0.2) {
        alpha = _alpha;
    }

    inline void fill_dead_end_neighbor_with_id(const VertexIdType &_id) {
        for (VertexIdType index = 0; index < num_deadend_vertices; ++index) {
            const VertexIdType &id = deadend_vertices[index];
            const VertexIdType &start = start_pos_in_out_neighbor_lists[id];
            out_neighbors_lists[start] = _id;
        }
    }


    inline void change_in_neighbors_adj(const VertexIdType &_sid, const VertexIdType &_target) {
        const VertexIdType &idx_start = start_pos_in_appearance_pos_lists[_sid];
        const VertexIdType &idx_end = start_pos_in_appearance_pos_lists[_sid + 1];
        for (VertexIdType index = idx_start; index < idx_end; ++index) {
            out_neighbors_lists[appearance_pos_lists[index]] = _target;
        }
    }

    inline void restore_neighbors_adj(const VertexIdType &_sid) {
        const VertexIdType &idx_start = start_pos_in_appearance_pos_lists[_sid];
        const VertexIdType &idx_end = start_pos_in_appearance_pos_lists[_sid + 1];
        for (VertexIdType index = idx_start; index < idx_end; ++index) {
            out_neighbors_lists[appearance_pos_lists[index]] = _sid;
        }
    }

    inline void set_source_and_alpha(const VertexIdType _sid, const PageRankScoreType _alpha) {
        sid = _sid;
        alpha = _alpha;
//        fill_dead_end_neighbor_with_id(_sid);
    }

    Graph() = default;

    ~Graph() = default;

    inline const VertexIdType &getNumOfVertices() const {
        return numOfVertices;
    }

    /**
     * @param _vid
     * @return return the original out degree
     */
    inline const VertexIdType &original_out_degree(const VertexIdType &_vid) const {
        assert(_vid < numOfVertices);
        return out_degrees[_vid];
    }

    inline const VertexIdType &get_neighbor_list_start_pos(const VertexIdType &_vid) const {
        assert(_vid < numOfVertices + 2);
        return start_pos_in_out_neighbor_lists[_vid];
    }


    inline const VertexIdType &getOutNeighbor(const VertexIdType &_index) const {
//        if (_index >= start_pos_in_out_neighbor_lists[dummy_id + 1]) {
//            MSG("Time to check " __FILE__)
//            MSG(__LINE__)
//        }
        assert(_index < start_pos_in_out_neighbor_lists[dummy_id + 1]);
        return out_neighbors_lists[_index];
    }

    inline const EdgeSizeType &getNumOfEdges() const {
        return numOfEdges;
    }


    void read_binary(const std::string &_attribute_file,
                     const std::string &_graph_file) {
        {
            std::string line;
            std::ifstream attribute_file(_attribute_file);
            if (attribute_file.is_open()) {
                std::getline(attribute_file, line);
                size_t start1 = line.find_first_of('=');
                numOfVertices = std::stoul(line.substr(start1 + 1));
                std::getline(attribute_file, line);
                size_t start2 = line.find_first_of('=');
                numOfEdges = std::stoul(line.substr(start2 + 1));
                dummy_id = numOfVertices;
                printf("The Number of Vertices: %d\n", numOfVertices);
                printf("The Number of Edges: %d\n", numOfEdges);
                attribute_file.close();
            } else {
                printf(__FILE__ "; LINE %d; File Not Exists.\n", __LINE__);
                std::cout << _attribute_file << std::endl;
                exit(1);
            }
        }
        const auto start = getCurrentTime();
        // create temporary graph
        std::vector<Edge> edges(numOfEdges);
        if (std::FILE *f = std::fopen(_graph_file.c_str(), "rb")) {
            size_t rtn = std::fread(edges.data(), sizeof edges[0], edges.size(), f);
            printf("Returned Value of fread: %zu\n", rtn);
            std::fclose(f);
        } else {
            printf("Graph::read; File Not Exists.\n");
            std::cout << _graph_file << std::endl;
            exit(1);
        }
        const auto end = getCurrentTime();
        printf("Time Used For Loading BINARY : %.2f\n", end - start);

        // read the edges
        // the ids must be in the range from [0 .... the number of vertices - 1];
        numOfEdges = 0;
        out_degrees.clear();
        out_degrees.resize(numOfVertices + 2, 0);
        in_degrees.clear();
        in_degrees.resize(numOfVertices + 2, 0);
        for (auto &edge : edges) {
            const VertexIdType &from_id = edge.from_id;
            const VertexIdType &to_id = edge.to_id;
            // remove self loop
            if (from_id != to_id) {
                //the edge read is a directed one
                ++out_degrees[from_id];
                ++in_degrees[to_id];
                ++numOfEdges;
            }
        }
        /* final count */
//        printf("%d-th Directed Edge Processed.\n", numOfEdges);

        // sort the adj list
//        for (auto &neighbors : matrix) {
//            std::sort(neighbors.begin(), neighbors.end());
//        }

        // process the dead_end
        uint32_t degree_max = 0;
        deadend_vertices.clear();
        for (VertexIdType i = 0; i < numOfVertices; ++i) {
            if (out_degrees[i] == 0) {
                deadend_vertices.emplace_back(i);
            }
            degree_max = std::max(degree_max, out_degrees[i]);
        }
        num_deadend_vertices = deadend_vertices.size();
        printf("The number of dead end vertices:%d\n", num_deadend_vertices);

        // process pos_list list
        start_pos_in_appearance_pos_lists.clear();
        start_pos_in_appearance_pos_lists.resize(numOfVertices + 2, 0);
        for (VertexIdType i = 0, j = 1; j < numOfVertices; ++i, ++j) {
            start_pos_in_appearance_pos_lists[j] = start_pos_in_appearance_pos_lists[i] + in_degrees[i];
        }
        start_pos_in_appearance_pos_lists[numOfVertices] = numOfEdges;

        // process out list
        start_pos_in_out_neighbor_lists.clear();
        start_pos_in_out_neighbor_lists.resize(numOfVertices + 2, 0);
        for (VertexIdType current_id = 0, next_id = 1; next_id < numOfVertices + 1; ++current_id, ++next_id) {
            start_pos_in_out_neighbor_lists[next_id] =
                    start_pos_in_out_neighbor_lists[current_id] + std::max(out_degrees[current_id], 1u);
        }
        // process dummy vertex
        assert(start_pos_in_out_neighbor_lists[numOfVertices] == numOfEdges + deadend_vertices.size());
        out_degrees[dummy_id] = 0;
        start_pos_in_out_neighbor_lists[numOfVertices + 1] = start_pos_in_out_neighbor_lists[numOfVertices];
        ////////////////////////////////////////////////////////////

        // compute the positions
        std::vector<VertexIdType> out_positions_to_fill(start_pos_in_out_neighbor_lists.begin(),
                                                        start_pos_in_out_neighbor_lists.end());
        // fill the edge list
        out_neighbors_lists.clear();
        out_neighbors_lists.resize(numOfEdges + num_deadend_vertices + degree_max, 0);
        uint32_t edges_processed = 0;
        uint32_t msg_gap = std::max(1u, numOfEdges / 10);
        std::vector<std::pair<VertexIdType, VertexIdType>> position_pair;
        position_pair.reserve(numOfEdges);
        for (auto &edge : edges) {
            const VertexIdType &from_id = edge.from_id;
            const VertexIdType &to_id = edge.to_id;
            // remove self loop
            if (from_id != to_id) {
                VertexIdType &out_position = out_positions_to_fill[from_id];
                assert(out_position < out_positions_to_fill[from_id + 1]);
                out_neighbors_lists[out_position] = to_id;
                position_pair.emplace_back(to_id, out_position);
                ++out_position;
                if (++edges_processed % msg_gap == 0) {
                    printf("%u edges processed.\n", edges_processed);
                }
            }
        }
        edges.clear();
        MSG(edges_processed);
        printf("%s\n", std::string(30, '-').c_str());
#ifdef USE_REVERSE_POS
        std::vector<VertexIdType> in_positions_to_fill(start_pos_in_appearance_pos_lists.begin(),
                                                       start_pos_in_appearance_pos_lists.end());
        in_positions_to_fill[numOfVertices] = numOfEdges;
        const double time_sort_start = getCurrentTime();
        std::sort(position_pair.begin(), position_pair.end(), std::less<>());
        const double time_sort_end = getCurrentTime();
//        MSG(time_sort_end - time_sort_start);
        printf("%s\n", std::string(30, '-').c_str());
        appearance_pos_lists.clear();
        appearance_pos_lists.resize(numOfEdges + num_deadend_vertices + degree_max, 0);
        uint32_t in_pos_pair = 0;
        for (const auto &pair : position_pair) {
            const VertexIdType &to_id = pair.first;
            const VertexIdType &pos = pair.second;
            VertexIdType &in_position = in_positions_to_fill[to_id];
            assert(in_position < in_positions_to_fill[to_id + 1]);
            appearance_pos_lists[in_position] = pos;
            ++in_position;
            if (++in_pos_pair % msg_gap == 0) {
                MSG(in_pos_pair);
            }
        }
#endif

        // fill the dummy ids
        for (const VertexIdType &id : deadend_vertices) {
            out_neighbors_lists[out_positions_to_fill[id]++] = dummy_id;
        }

#ifdef DEBUG_MODE
//        show();
        // sanity check
//        std::memcpy(in_positions_to_fill.data(), start_pos_in_appearance_pos_lists.data(),
//                    sizeof(VertexIdType) * numOfVertices);
//        for (uint32_t index = 0; index < numOfEdges + num_deadend_vertices; ++index) {
//            const VertexIdType id = out_neighbors_lists[index];
//            if (id != numOfVertices) {
//                VertexIdType &in_position = in_positions_to_fill[id];
//                if (appearance_pos_lists[in_position] != index) {
//                    printf("Error in " __FILE__  ", line %d, in position error.\n", __LINE__);
//                    MSG(index);
//                    MSG(id);
//                    MSG(in_position);
//                    exit(1);
//                }
//                ++in_position;
//            }
//        }
#endif
        const double time_end = getCurrentTime();
        printf("Graph Build Finished. TIME: %.4f\n", time_end - start);
        printf("%s\n", std::string(110, '-').c_str());

    }

    void show() const {
        // we need to show the dummy
        const VertexIdType num_to_show = std::min(numOfVertices + 1, 50u);
        // show the first elements
        show_vector("The Out Degrees of The Vertices:",
                    std::vector<VertexIdType>(out_degrees.data(), out_degrees.data() + num_to_show));
        show_vector("The Start Positions of The Vertices in Out Neighbor Lists:",
                    std::vector<VertexIdType>(start_pos_in_out_neighbor_lists.data(),
                                              start_pos_in_out_neighbor_lists.data() + num_to_show));
        show_vector("The In Degrees of The Vertices:",
                    std::vector<VertexIdType>(in_degrees.data(), in_degrees.data() + num_to_show));
        show_vector("The Start Positions of The Vertices in Appearance List:",
                    std::vector<VertexIdType>(start_pos_in_appearance_pos_lists.data(),
                                              start_pos_in_appearance_pos_lists.data() + num_to_show));
        // assume that the number of vertices >= the number of edges; otherwise, there is a potential bug here.
        show_vector("Out Neighbor Lists:",
                    std::vector<VertexIdType>(out_neighbors_lists.data(),
                                              out_neighbors_lists.data() +
                                              std::min(numOfEdges + num_deadend_vertices, 50u)));
        show_vector("The Appearance Positions of Vertices in the Out Neighbor Lists:",
                    std::vector<VertexIdType>(appearance_pos_lists.data(),
                                              appearance_pos_lists.data() + std::min(numOfEdges, 50u)));
//        show_vector("The adj list of the middel vertex", matrix[numOfVertices / 2]);
        printf("The position the id appears in outNeighbor List:\n");
        for (VertexIdType id = 0; id < numOfVertices; ++id) {
            const VertexIdType &idx_start = start_pos_in_appearance_pos_lists[id];
            const VertexIdType &idx_end = start_pos_in_appearance_pos_lists[id + 1];
            printf("Id:%u;\tPositions: ", id);
            for (VertexIdType index = idx_start; index < idx_end; ++index) {
                printf("%u, ", appearance_pos_lists[index]);
            }
            printf("\n");
        }
        show_vector("Dead End Vertices List:",
                    std::vector<VertexIdType>(deadend_vertices.data(),
                                              deadend_vertices.data() +
                                              std::min(num_deadend_vertices, 50u)));
        printf("\n%s\n", std::string(120, '=').c_str());
    }
};


#endif //SPEEDPPR_GRAPH_H
