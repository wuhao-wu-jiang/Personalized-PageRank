

#ifndef SPEEDPPR_MYQUEUE_H
#define SPEEDPPR_MYQUEUE_H

#include <vector>
#include <cmath>
#include <iostream>
#include <cassert>
#include "BasicDefinition.h"

class MyQueue {
private:
    const VertexIdType mask;
    std::vector<VertexIdType> queue;
    VertexIdType num = 0;
    VertexIdType idx_front = 0;
    VertexIdType idx_last_plus_one = 0;
private:
    static inline VertexIdType compute_queue_size(const VertexIdType &_numOfVertices) {
        return (1u) << (uint32_t) ceil(log2(_numOfVertices + 2u));
    }

public:
    explicit MyQueue(const VertexIdType &_numOfVertices) :
            mask(compute_queue_size(_numOfVertices) - 1),
            queue(mask + 2u, 0) {
#ifdef DEBUG_MODE
//        std::cout << "Queue Size: " << queue.size() << " Mask: " << mask << std::endl;
#endif
    }

    inline void clear() {
        idx_front = 0;
        idx_last_plus_one = 0;
        num = 0;
    }


    inline const VertexIdType &size() const {
        return num;
    }

    inline const VertexIdType &front() const {
        return queue[idx_front];
    }

    inline void pop() {
        --num;
        ++idx_front;
        idx_front &= mask;
    }

    inline void push(const VertexIdType &_elem) {
        ++num;
        queue[idx_last_plus_one] = _elem;
        ++idx_last_plus_one;
        idx_last_plus_one &= mask;
    }

    inline bool empty() const {
#ifdef DEBUG_MODE
        assert(num != 0 || idx_last_plus_one == idx_front);
        assert(num == 0 || idx_last_plus_one != idx_front);
#endif
        return idx_last_plus_one == idx_front;
    }
};

struct FwdPushStructure {
    // reserve one slot for the dummy vertex
    MyQueue active_vertices;
    // reserve one slot for the dummy vertex
    std::vector<bool> is_active;

    explicit FwdPushStructure(const VertexIdType &numOfVertices) :
            active_vertices(numOfVertices + 1),
            is_active(numOfVertices + 1, false) {}
};

#endif //SPEEDPPR_MYQUEUE_H
