

#ifndef SPEEDPPR_QUERYGENERATOR_H
#define SPEEDPPR_QUERYGENERATOR_H


#include <vector>
#include "BasicDefinition.h"
#include "HelperFunctions.h"
#include "MyRandom.h"
#include "Graph.h"

struct QueryGenerator {
    static void generate(const Graph &_graph, const std::string &_file_name) {
        //generate random query node
        std::vector<VertexIdType> source_vertices;
        while (source_vertices.size() < param.query_size) {
            const VertexIdType sid = SFMT64::uniform_int(_graph.getNumOfVertices());
            if (_graph.original_out_degree(sid) != 0) {
                source_vertices.push_back(sid);
            }
        }
        if (FILE *file = fopen(_file_name.c_str(), "w")) {
            for (const auto &sid : source_vertices) {
                fprintf(file, "%u\n", sid);
            }
            fclose(file);
        } else {
            printf("ERROR IN " __FILE__ " LINE %u. FILE CAN'T BE OPENED \n", __LINE__);
            printf("%s\n", _file_name.c_str());
            exit(1);
        }
    }
};

#endif //SPEEDPPR_QUERYGENERATOR_H
