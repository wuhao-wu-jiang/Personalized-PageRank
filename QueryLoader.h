

#ifndef SPEEDPPR_QUERYLOADER_H
#define SPEEDPPR_QUERYLOADER_H


#include <vector>
#include <string>
#include <fstream>
#include "BasicDefinition.h"
#include "HelperFunctions.h"

class QueryLoader {
    std::vector<VertexIdType> vertices;
public:
    QueryLoader(const std::string &_query_file_name, const uint32_t &_query_size) {
        std::ifstream file(_query_file_name.c_str());
        if (file.good() == false) {
            printf("Query File Not Exists.\n");
            exit(1);
        }
        for (VertexIdType sid; (file >> sid) && vertices.size() < _query_size;) {
            vertices.emplace_back(sid);
        }
        if (vertices.empty()) {
            printf("Error in QueryLoader::QueryLoader. Empty Query File\n");
        } else if (vertices.size() > _query_size) {
            printf("Error in QueryLoader::QueryLoader.\n");
        } else {
            printf("The number of query points loaded: %zu\n", vertices.size());
        }
        file.close();
    }

    const std::vector<VertexIdType> &source_vertices() const {
        return vertices;
    }

};


#endif //SPEEDPPR_QUERYLOADER_H
