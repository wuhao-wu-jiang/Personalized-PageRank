
#ifndef SPEEDPPR_BASICDEFINITION_H
#define SPEEDPPR_BASICDEFINITION_H

#include <queue>

#ifndef DEBUG_MODE
#define DEBUG_MODE
#endif

#ifndef FAST_READ
#define FAST_READ
#endif

typedef unsigned int VertexIdType;
typedef unsigned int EdgeSizeType;
typedef double PageRankScoreType;


template<class FLOAT_TYPE>
struct IdScorePair {
    VertexIdType id = 0;
    FLOAT_TYPE score = 0;


    IdScorePair(const VertexIdType &_id = 0, const FLOAT_TYPE &_score = 0) :
            id(_id), score(_score) {}
};

template<class FLOAT_TYPE>
struct IdScorePairComparatorGreater {
    // Compare 2 IdScorePair objects using name
    bool operator()(const IdScorePair<FLOAT_TYPE> &pair1, const IdScorePair<FLOAT_TYPE> &pair2) {
        return pair1.score > pair2.score || pair1.score == pair2.score && pair1.id < pair2.id;
    }
};

template<class FLOAT_TYPE>
struct IdScorePairComparatorLess {
    // Compare 2 IdScorePair objects using name
    bool operator()(const IdScorePair<FLOAT_TYPE> &pair1, const IdScorePair<FLOAT_TYPE> &pair2) {
        return pair1.score < pair2.score || pair1.score == pair2.score && pair1.id < pair2.id;
    }
};

typedef std::priority_queue<IdScorePair<float>, std::vector<IdScorePair<float> >, IdScorePairComparatorGreater<float> > IdScorePairMaxQueue_float;

struct Edge {
    VertexIdType from_id;
    VertexIdType to_id;

    Edge() : from_id(0), to_id(0) {}

    Edge(const VertexIdType &_from, const VertexIdType &_to) :
            from_id(_from), to_id(_to) {}

    bool operator<(const Edge &_edge) const {
        return from_id < _edge.from_id || (from_id == _edge.from_id && to_id < _edge.to_id);
    }
};

#endif //SPEEDPPR_BASICDEFINITION_H
