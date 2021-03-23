

#ifndef SPEEDPPR_CLEANGRAPH_H
#define SPEEDPPR_CLEANGRAPH_H

#ifndef DEEP_CLEAN
#define DEEP_CLEAN
#endif

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <limits>
#include <unordered_set>
#include <algorithm>
#include "BasicDefinition.h"
#include "HelperFunctions.h"

class CleanGraph {
    VertexIdType numOfVertices = 0;
    EdgeSizeType numOfEdges = 0;
public:

    void duplicate_edges(const std::string &_input_file,
                         const std::string &_output_file) {
        ////////////////////////////////////////////////////////////
        std::ifstream inf(_input_file.c_str());
        if (!inf.is_open()) {
            printf("CleanGraph::duplicate_edges; File not exists.\n");
            printf("%s\n", _input_file.c_str());
            exit(1);
        }
        ////////////////////////////////////////////////////////////
        // status indicator
        printf("\nReading Input Graph\n");

        std::string line;
        /**
         * skip the headers, we assume the headers are the comments that
         * begins with '#'
         */
        while (std::getline(inf, line) && line[0] == '#') {}
        if (line.empty() || !isdigit(line[0])) {
            printf("Error in CleanGraph::duplicate_edges.Raw File Format Error.\n");
            printf("%s\n", line.c_str());
            exit(1);
        }
        ////////////////////////////////////////////////////////////////////////////////
        // create temporary graph
        std::vector<Edge> edges;
        numOfEdges = 0;
        /**
         * read the raw file
         */
        size_t num_lines = 0;
        // process the first line=
        {
            VertexIdType fromId, toID;
            ++num_lines;
            size_t end = 0;
            fromId = std::stoul(line, &end);
            toID = std::stoul(line.substr(end));
            edges.emplace_back(fromId, toID);
            // duplicate the edge
            edges.emplace_back(toID, fromId);
        }
        // read the edges
        for (VertexIdType fromId, toID; inf >> fromId >> toID;) {
            edges.emplace_back(fromId, toID);
            // duplicate the edge
            edges.emplace_back(toID, fromId);
            if (++num_lines % 5000000 == 0) { printf("%zu Valid Lines Read.\n", num_lines); }
        }
        ////////////////////////////////////////////////////////////////////////////////

        // close the file
        inf.close();
        /* final count */
        printf("%zu Lines Read.\n", num_lines);
        printf("Finish Reading.\n");
        numOfEdges = edges.size();
        printf("There are %d Edges After Duplications.\n", numOfEdges);

        // write the file

        if (FILE *file = fopen(_output_file.c_str(), "w")) {
            for (const auto &edge : edges) {
                fprintf(file, "%u\t%u\n", edge.from_id, edge.to_id);
            }
            printf("%d Edges Written. Duplication Finished.\n", numOfEdges);
            printf("\n%s\n", std::string(80, '=').c_str());
            fclose(file);
        } else {
            printf("CleanGraph::duplicate_edges; Output File Can Not Be Openned.\n");
            printf("%s\n", _output_file.c_str());
        }
    }

    /**
     * Each vertex is then given a label between 0 to (#vertices - 1).
     * Currently We don't deal with parallel edges.
     * @param _input_file
     * @param _is_directed
     * @param _output_folder
     */
    void clean_graph(const std::string &_input_file,
                     const std::string &_output_folder) {
        ////////////////////////////////////////////////////////////
        std::ifstream inf(_input_file.c_str());
        if (!inf.is_open()) {
            printf("CleanGraph::clean_graph; File not exists.\n");
            printf("%s\n", _input_file.c_str());
            exit(1);
        }
        ////////////////////////////////////////////////////////////
        // status indicator
        printf("\nReading Input Graph\n");

        std::string line;
        /**
         * skip the headers, we assume the headers are the comments that
         * begins with '#'
         */
        while (std::getline(inf, line) && line[0] == '#') {}
        if (line.empty() || !isdigit(line[0])) {
            printf("Error in CleanGraph::clean_graph. Raw File Format Error.\n");
            printf("%s\n", line.c_str());
            exit(1);
        }
        ////////////////////////////////////////////////////////////////////////////////
        // create temporary graph
        std::vector<Edge> edges;
        numOfEdges = 0;
        /**
         * read the raw file
         */
        size_t num_lines = 0;
        // process the first line
        {
            VertexIdType fromId, toID;
            ++num_lines;
            size_t end = 0;
            fromId = std::stoul(line, &end);
            toID = std::stoul(line.substr(end));
            // remove self-loops
#ifdef DEEP_CLEAN
            if (fromId != toID) { edges.emplace_back(fromId, toID); }
#else
            edges.emplace_back(fromId, toID);
#endif
        }
        // read the edges
        for (VertexIdType fromId, toID; inf >> fromId >> toID;) {
#ifdef DEEP_CLEAN
            if (fromId != toID) { edges.emplace_back(fromId, toID); }
#else
            edges.emplace_back(fromId, toID);
#endif
            if (++num_lines % 5000000 == 0) { printf("%zu Valid Lines Read.\n", num_lines); }
        }
        ////////////////////////////////////////////////////////////////////////////////

        // close the file
        inf.close();
        /* final count */
        printf("%zu Lines Read.\n", num_lines);
        numOfEdges = edges.size();
        printf("%d-th Non-Self Loop Edges.\n", numOfEdges);
        printf("Finish Reading.\n");
        printf("\n%s\n", std::string(30, '-').c_str());

        // find the maximum id
        size_t id_max = 0;
        size_t id_min = std::numeric_limits<uint32_t>::max();
        for (const auto &pair : edges) {
            id_max = std::max(id_max, (size_t) std::max(pair.from_id, pair.to_id));
            id_min = std::min(id_min, (size_t) std::min(pair.from_id, pair.to_id));
        }
        printf("Maximum ID: %zu\n", id_max);
        printf("Minimum ID: %zu\n", id_min);
        if (id_max >= std::numeric_limits<uint32_t>::max()) {
            printf("Warning. Change VertexIdType First.\n");
            exit(1);
        }
        const VertexIdType one_plus_id_max = id_max + 1;
        std::vector<VertexIdType> out_degree(one_plus_id_max, 0);
        std::vector<VertexIdType> in_degree(one_plus_id_max, 0);
        // compute the degrees.
        for (const auto &edge : edges) {
            ++out_degree[edge.from_id];
            ++in_degree[edge.to_id];
        }
        // count the number of dead-end vertices
        uint32_t original_dead_end_num = 0;
        uint32_t num_isolated_points = 0;
        uint32_t max_degree = 0;
        for (VertexIdType id = 0; id < one_plus_id_max; ++id) {
            if (out_degree[id] == 0) {
                ++original_dead_end_num;
                if (in_degree[id] == 0) {
                    ++num_isolated_points;
                }
            }
            // compute maximum out degree
            max_degree = std::max(out_degree[id], max_degree);
        }
        printf("The number of dead end vertices: %u\n", original_dead_end_num);
        printf("The number of Isolated Points: %u\n", num_isolated_points);
        printf("The maximum out degree is: %u\n", max_degree);

        // re-label the vertices.
#ifdef DEEP_CLEAN
        std::vector<std::vector<VertexIdType>> out_adj_list(one_plus_id_max, std::vector<VertexIdType>());
        for (const auto &edge : edges) {
            out_adj_list[edge.from_id].push_back(edge.to_id);
        }
        numOfVertices = 0;
        std::vector<bool> visited(one_plus_id_max, false);
        std::vector<VertexIdType> id_map(one_plus_id_max, 0);
        for (VertexIdType index = 0; index < one_plus_id_max; ++index) {
            // skip singleton and visited vertices
            if (!visited[index] && out_degree[index] > 0) {
                std::vector<VertexIdType> current_level({index});
                visited[index] = true;
                std::vector<VertexIdType> next_level;
                while (!current_level.empty()) {
                    for (const VertexIdType &current_id : current_level) {
                        id_map[current_id] = numOfVertices;
                        ++numOfVertices;
                        for (const VertexIdType &nid : out_adj_list[current_id]) {
                            if (!visited[nid]) {
                                next_level.push_back(nid);
                                visited[nid] = true;
                            }
                        }
                    }
                    current_level.swap(next_level);
                    next_level.clear();
                }
            }
        }
        printf("Relabeling Finish.\n");
        printf("The number of valid vertices: %u\n", numOfVertices);
        // re-write the edges ids
        for (auto &edge : edges) {
            edge.from_id = id_map[edge.from_id];
            edge.to_id = id_map[edge.to_id];
        }
#else
        // we assume the vertice ids are in the arrange of 0 ... numOfVertices - 1
        numOfVertices = one_plus_id_max;
#endif

        // sort the edges
        std::sort(edges.begin(), edges.end());

        // Write the attribute file
        numOfEdges = edges.size();
        std::string attribute_file = _output_folder + '/' + "attribute.txt";
        if (std::FILE *file = std::fopen(attribute_file.c_str(), "w")) {
            std::fprintf(file, "n=%d\nm=%d\n", numOfVertices, numOfEdges);
            std::fclose(file);
        } else {
            printf("Graph::clean_graph; File Not Exists.\n");
            printf("%s\n", attribute_file.c_str());
            exit(1);
        }

        // write the graph in txt
        std::string graph_txt_file = _output_folder + '/' + "graph.txt";
        if (std::FILE *file = std::fopen(graph_txt_file.c_str(), "w")) {
            fprintf(file, "%u\n", numOfVertices);
            for (const auto &edge : edges) {
                fprintf(file, "%u\t%u\n", edge.from_id, edge.to_id);
            }
            printf("Writing Txt Finished.\n");
            std::fclose(file);
        } else {
            printf("Graph::clean_graph; File Not Exists.\n");
            printf("%s\n", graph_txt_file.c_str());
            exit(1);
        }

        // write the graph in binary
        std::string graph_bin_file = _output_folder + '/' + "graph.bin";
        if (std::FILE *file = std::fopen(graph_bin_file.c_str(), "wb")) {
            std::fwrite(edges.data(), sizeof edges[0], edges.size(), file);
            printf("Writing Binary Finished.\n");
            std::fclose(file);
        } else {
            printf("Graph::clean_graph; File Not Exists.\n");
            printf("%s\n", graph_bin_file.c_str());
            exit(1);
        }
        printf("%s\n", std::string(110, '=').c_str());
    }

};


#endif //SPEEDPPR_CLEANGRAPH_H
