

#ifndef SPEEDPPR_ITERATEDMETHODS_H
#define SPEEDPPR_ITERATEDMETHODS_H

#ifndef SHOW_TIME_POWER
#define SHOW_TIME_POWER
#endif

//#ifndef SHOW_DEBUG_POWER
//#define SHOW_DEBUG_POWER
//#endif

#ifndef SHOW_SERIES_POWER
#define SHOW_SERIES_POWER
#endif

#include <vector>
#include <cassert>
#include <cmath>
#include <numeric>
#include <queue>
#include "BasicDefinition.h"
#include "Graph.h"
#include <thread>
#include "MyQueue.h"


class IteratedMethods {

private:
    const double alpha;
    const double one_minus_alpha;
    const double l1_error;
    const VertexIdType numOfVertices;
    const Graph &graph;
//    std::vector<ExecutionLog> logs;

private:
    unsigned int computeNumOfIterations(const PageRankScoreType &_l1_error) const {
        assert(_l1_error <= 1);
        assert(alpha <= 1);
        return static_cast<uint32_t>(log(2.0 / _l1_error) / log(1.0 / (1.0 - alpha)));
    }

public:

    IteratedMethods(const Graph &_graph,
                    const PageRankScoreType _l1_error) :
            alpha(_graph.get_alpha()),
            one_minus_alpha(1 - _graph.get_alpha()),
            l1_error(_l1_error),
            numOfVertices(_graph.getNumOfVertices()),
            graph(_graph) {
    }

private:
    void
    ground_truth_iteration(const VertexIdType &sid, std::vector<PageRankScoreType> &pi,
                           std::vector<PageRankScoreType> &residuals) const {
///////////////////////////////////////////////////////////////////////////////////////////////////
        const double time_start = getCurrentTime();
        std::cout.precision(std::numeric_limits<double>::max_digits10);
        assert(pi.size() >= numOfVertices + 1);
        assert(residuals.size() >= numOfVertices + 1);
        // we must assume that pi and residual has size >= numOfVertices + 1
        std::fill(pi.begin(), pi.end(), 0);
        std::fill(residuals.begin(), residuals.end(), 0);
        residuals[sid] = 1;
        size_t number_of_pushes = 0;
        double real_r_sum = one_minus_alpha;
        const int num_epoch = 8;
        for (int epoch = 0, iter = 0; epoch < num_epoch && real_r_sum > l1_error; ++epoch) {
            double l1_error_this_epoch = pow(l1_error, (1.0 + epoch) / num_epoch);
            const double threshold = l1_error_this_epoch / graph.getNumOfEdges() / 2;
            uint32_t num_iter = computeNumOfIterations(l1_error_this_epoch);
            for (unsigned int k = 0; k < num_iter && real_r_sum > l1_error_this_epoch; ++iter, ++k) {
                for (VertexIdType id = 0, next_id = 1; id < numOfVertices; ++id, ++next_id) {
                    const VertexIdType &idx_start = graph.get_neighbor_list_start_pos(id);
                    const VertexIdType &idx_end = graph.get_neighbor_list_start_pos(next_id);
                    const VertexIdType degree = idx_end - idx_start;
                    if (residuals[id] > threshold * degree) {
                        const double alpha_residual = alpha * residuals[id];
                        pi[id] += alpha_residual;
                        const double increment = (residuals[id] - alpha_residual) / degree;
                        residuals[id] = 0;
                        number_of_pushes += degree;
                        for (unsigned int j = idx_start; j < idx_end; ++j) {
                            const VertexIdType &nid = graph.getOutNeighbor(j);
                            residuals[nid] += increment;
                        }
                    }
                }
                const double r_sum = std::accumulate(residuals.begin(), residuals.end(), 0.0);
                real_r_sum = one_minus_alpha * (r_sum - residuals[numOfVertices]) / (1 - residuals[numOfVertices]);
                {
                    const double time_used = getCurrentTime() - time_start;
                    printf("Sid: %u \tTime Used: %.1f\tr_sum: %.19f\n", sid, time_used, real_r_sum);
                }
            }
        }
        const double dummy_scale_factor = 1.0 / (1 - residuals[numOfVertices]);
        {
            double r_sum = std::accumulate(residuals.begin(), residuals.end(), 0.0);
            real_r_sum = one_minus_alpha * (r_sum - residuals[numOfVertices]) * dummy_scale_factor;
            const double time_used = getCurrentTime() - time_start;
            printf("Sid: %u \tTime Used: %.1f\tr_sum: %.19f\n", sid, time_used, real_r_sum);
        }
        for (VertexIdType id = 0; id < numOfVertices; ++id) {
            pi[id] *= dummy_scale_factor;
            residuals[id] *= dummy_scale_factor;
            pi[id] += residuals[id] * alpha;
        }
    }

    void run_ground_truth(const std::vector<VertexIdType> &_sids) const {
        std::vector<PageRankScoreType> pi(numOfVertices + 1, 0);
        std::vector<PageRankScoreType> residuals(numOfVertices + 1, 0);
        for (const auto &sid : _sids) {
            const double time_start = getCurrentTime();
            ground_truth_iteration(sid, pi, residuals);
            const double time_end = getCurrentTime();
            printf("-----------------\nTime Used for Sid:%d\t%.2f\n", sid, time_end - time_start);
            save_answer(pi, numOfVertices, param.answer_folder + "/" + std::to_string(sid) + ".txt");
        }
    }

public:

    void
    forward(const VertexIdType &sid, std::vector<PageRankScoreType> &pi,
            std::vector<PageRankScoreType> &residuals) const {
        long long number_of_pushes = 0;
        long long previous_number_of_pushes = 0;
        if (pi.size() < numOfVertices + 1) {
            pi.resize(numOfVertices + 1, 0);
            residuals.resize(numOfVertices + 1, 0);
        }
        std::fill(pi.begin(), pi.end(), 0);
        std::fill(residuals.begin(), residuals.end(), 0);
        residuals[sid] = 1;
        double r_sum = 1.0;
        MyQueue active_vertices(numOfVertices);
        std::vector<bool> is_active(numOfVertices, false);
        active_vertices.push(sid);
        is_active[sid] = true;
        // this prevents dummy vertex from entering the queue
        is_active[numOfVertices] = true;
//        const double threshold_to_reject = l1_error;
        const double threshold_to_reject = l1_error / graph.getNumOfEdges();
        const double time_start = getCurrentTime();
        while (!active_vertices.empty() && r_sum > l1_error) {
#ifdef DEBUG_MODE
            if (number_of_pushes - previous_number_of_pushes >= graph.getNumOfEdges()) {
                previous_number_of_pushes = number_of_pushes;
                const size_t num_iter = number_of_pushes / graph.getNumOfEdges();
                const double time_current = getCurrentTime();
                printf("#Iter:%s%lu\tr_sum:%.12f\tTime Used:%.4f\t#Pushes:%llu\n",
                       (num_iter < 10 ? "0" : ""), num_iter, r_sum, time_current - time_start, number_of_pushes);
            }
#endif
            const VertexIdType id = active_vertices.front();
            active_vertices.pop();
            is_active[id] = false;
            const VertexIdType &idx_start = graph.get_neighbor_list_start_pos(id);
            const VertexIdType &idx_end = graph.get_neighbor_list_start_pos(id + 1);
            const VertexIdType degree = idx_end - idx_start;
            if (residuals[id] > degree * threshold_to_reject) {
                const double alpha_residual = alpha * residuals[id];
                pi[id] += alpha_residual;
                r_sum -= alpha_residual;
                const double increment = (residuals[id] - alpha_residual) / degree;
                residuals[id] = 0;
                number_of_pushes += degree;
                for (unsigned int j = idx_start; j < idx_end; ++j) {
                    const VertexIdType &nid = graph.getOutNeighbor(j);
                    residuals[nid] += increment;
                    if (!is_active[nid]) {
                        active_vertices.push(nid);
                        is_active[nid] = true;
                    }
                }
            }
        }
        const double dummy_scale_factor = 1.0 / (1 - residuals[graph.get_dummy_id()]);
#ifdef DEBUG_MODE
//        MSG(dummy_scale_factor)
#endif
        for (VertexIdType id = 0; id < numOfVertices; ++id) {
            pi[id] *= dummy_scale_factor;
            residuals[id] *= dummy_scale_factor;
        }
#ifdef DEBUG_MODE
        {
            const double time_current = getCurrentTime();
            const size_t num_iter = number_of_pushes / graph.getNumOfEdges();
            printf("#Iter:%s%lu\tr_sum:%.12f\tTime Used:%.4f\t#Pushes:%llu\n",
                   (num_iter < 10 ? "0" : ""), num_iter, std::accumulate(residuals.begin(), residuals.end(), 0.0),
                   time_current - time_start, number_of_pushes);
        }
#endif
    }


    void
    naive_power_iteration(const VertexIdType &sid, std::vector<PageRankScoreType> &pi,
                          std::vector<PageRankScoreType> &residuals) const {
        if (pi.size() < numOfVertices + 1) {
            pi.resize(numOfVertices + 1, 0);
            residuals.resize(numOfVertices + 1, 0);
        }
        std::fill(pi.begin(), pi.end(), 0);
        std::fill(residuals.begin(), residuals.end(), 0);
        residuals[sid] = 1;
        double r_sum = 1.0;
        size_t number_of_pushes = 0;
        std::vector<PageRankScoreType> new_residuals(numOfVertices + 1, 0);
        printf("#Iter:00\tr_sum:%.12f\tTime Used:%.4f\t#Pushes:0\n", r_sum, 0.0);
        const double time_start = getCurrentTime();
///////////////////////////////////////////////////////////////////////////////////////////////////
        for (uint32_t num_iter = 0; r_sum > l1_error; ++num_iter) {
            for (VertexIdType id = 0, next_id = 1; id < numOfVertices; ++id, ++next_id) {
                const VertexIdType &idx_start = graph.get_neighbor_list_start_pos(id);
                const VertexIdType &idx_end = graph.get_neighbor_list_start_pos(next_id);
                const VertexIdType degree = idx_end - idx_start;
                number_of_pushes += degree;
                const double alpha_residual = alpha * residuals[id];
                pi[id] += alpha_residual;
                r_sum -= alpha_residual;
                const double increment = (residuals[id] - alpha_residual) / degree;
                residuals[id] = 0;
                for (unsigned int j = idx_start; j < idx_end; ++j) {
                    const VertexIdType &nid = graph.getOutNeighbor(j);
                    new_residuals[nid] += increment;
                }

            }
            residuals.swap(new_residuals);
            const double time_current = getCurrentTime();
            printf("#Iter:%s%u\tr_sum:%.12f\tTime Used:%.4f\t#Pushes:%zu\n",
                   (num_iter < 10 ? "0" : ""), num_iter, r_sum, time_current - time_start, number_of_pushes);
        }
#ifdef DEBUG_MODE
        if (residuals[graph.get_dummy_id()] != 0) {
            printf("ERROR IN " __FILE__ " LINE %u. FORGET TO DEAL WITH DEAD END.\n", __LINE__);
            exit(1);
        }
        do {
            const double time_current = getCurrentTime();
            const size_t num_iter = number_of_pushes / graph.getNumOfEdges();
            printf("#Iter:%s%lu\tr_sum:%.12f\tTime Used:%.4f\t#Pushes:%lu\n",
                   (num_iter < 10 ? "0" : ""), num_iter, std::accumulate(residuals.begin(), residuals.end(), 0.0),
                   time_current - time_start, number_of_pushes);
        } while (false);
#endif
    }

    void
    forward_iteration(const VertexIdType &sid, std::vector<PageRankScoreType> &pi,
                      std::vector<PageRankScoreType> &residuals,
                      FwdPushStructure &fwdPushStructure) {
        // don't use this function for multi-thread computing
        const double time_start = getCurrentTime();
        if (pi.size() < numOfVertices + 1) {
            pi.resize(numOfVertices + 1, 0);
            residuals.resize(numOfVertices + 1, 0);
        }
        std::fill(pi.begin(), pi.end(), 0);
        std::fill(residuals.begin(), residuals.end(), 0);
///////////////////////////////////////////////////////////////////////////////////////////////////
        residuals[sid] = 1;
        double r_sum = 1.0;
        double threshold_to_reject = 1.0 / graph.getNumOfEdges();
        size_t number_of_pushes = 0;
        size_t previous_number_of_pushes = 0;
///////////////////////////////////////////////////////////////////////////////////////////////////
//        // reserve one slot for the dummy vertex
//        MyQueue active_vertices(numOfVertices + 1);
//        // reserve one slot for the dummy vertex
//        std::vector<bool> is_active(numOfVertices + 1, false);
        auto &active_vertices = fwdPushStructure.active_vertices;
        auto &is_active = fwdPushStructure.is_active;
        active_vertices.clear();
        std::fill(is_active.begin(), is_active.end(), false);
        //////////////////////////////////////////////////////////
        active_vertices.push(sid);
        is_active[sid] = true;
        // this prevents dummy vertex from entering the queue
        is_active[numOfVertices] = true;
        // threshold to power iteration
        uint32_t switch_size = (numOfVertices >> 2u);
        while (!active_vertices.empty() && active_vertices.size() <= switch_size) {
            const VertexIdType id = active_vertices.front();
            active_vertices.pop();
            is_active[id] = false;
            const VertexIdType &idx_start = graph.get_neighbor_list_start_pos(id);
            const VertexIdType &idx_end = graph.get_neighbor_list_start_pos(id + 1);
            const VertexIdType degree = idx_end - idx_start;
#ifdef DEBUG_MODE
//            assert(degree != 0);
#endif
            // we need to skip the dummy vertex, which has degree = 0 by default
            if (residuals[id] > degree * threshold_to_reject) {
                const double alpha_residual = alpha * residuals[id];
                pi[id] += alpha_residual;
                r_sum -= alpha_residual;
                const double increment = (residuals[id] - alpha_residual) / degree;
                residuals[id] = 0;
                number_of_pushes += degree;
                for (unsigned int j = idx_start; j < idx_end; ++j) {
                    const VertexIdType &nid = graph.getOutNeighbor(j);
                    residuals[nid] += increment;
                    if (!is_active[nid]) {
                        active_vertices.push(nid);
                        is_active[nid] = true;
                    }
                }
            }
#ifdef SHOW_SERIES_POWER
            if (number_of_pushes - previous_number_of_pushes >= graph.getNumOfEdges()) {
                previous_number_of_pushes = number_of_pushes;
                const size_t num_iter = number_of_pushes / graph.getNumOfEdges();
                const double time_used = getCurrentTime() - time_start;
//                logs.emplace_back(r_sum, time_used, number_of_pushes);
                printf("#Iter:%s%lu\tr_sum:%.12f\tTime Used:%.4f\t#Pushes:%zu\n",
                       (num_iter < 10 ? "0" : ""), num_iter, r_sum, time_used, number_of_pushes);
            }
#endif
        }
#ifdef SHOW_SERIES_POWER
        //        {
        //            const double time_current = getCurrentTime();
        //            printf("r_sum:%.12f\tTime Used:%.4f\t#Pushes:%zu\n", r_sum, time_current - time_start, number_of_pushes);
        //        }
#endif
///////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef SHOW_SERIES_POWER
        const double time_epoch_start = getCurrentTime();
#endif
        const int num_epoch = 8;
        double real_r_sum = one_minus_alpha * (r_sum - residuals[numOfVertices]) / (1 - residuals[numOfVertices]);
        for (int epoch = 0, iter = 0; epoch < num_epoch && real_r_sum > l1_error; ++epoch) {
            double l1_error_this_epoch = pow(l1_error, (1.0 + epoch) / num_epoch);
            const double threshold = l1_error_this_epoch / graph.getNumOfEdges() / 2;
            uint32_t num_iter = computeNumOfIterations(l1_error_this_epoch);
            for (unsigned int k = 0; k < num_iter && real_r_sum > l1_error_this_epoch; ++iter, ++k) {
//                // id < numOfVertices skip the dummy
                for (VertexIdType id = 0, next_id = 1; id < numOfVertices; ++id, ++next_id) {
                    const VertexIdType &idx_start = graph.get_neighbor_list_start_pos(id);
                    const VertexIdType &idx_end = graph.get_neighbor_list_start_pos(next_id);
                    const VertexIdType degree = idx_end - idx_start;
                    if (residuals[id] > threshold * degree) {
                        const double alpha_residual = alpha * residuals[id];
                        pi[id] += alpha_residual;
                        r_sum -= alpha_residual;
                        const double increment = (residuals[id] - alpha_residual) / degree;
                        residuals[id] = 0;
                        number_of_pushes += degree;
                        for (unsigned int j = idx_start; j < idx_end; ++j) {
                            const VertexIdType &nid = graph.getOutNeighbor(j);
                            residuals[nid] += increment;
                        }
                    }
#ifdef SHOW_SERIES_POWER
                    if (number_of_pushes - previous_number_of_pushes >= graph.getNumOfEdges()) {
                        previous_number_of_pushes = number_of_pushes;
                        const size_t num_round = number_of_pushes / graph.getNumOfEdges();
                        real_r_sum =
                                one_minus_alpha * (r_sum - residuals[numOfVertices]) / (1 - residuals[numOfVertices]);
                        const double time_used = getCurrentTime() - time_start;
                        printf("#Iter:%s%lu\tr_sum:%.12f\tTime Used:%.4f\t#Pushes:%zu\n",
                               (num_round < 10 ? "0" : ""), num_round, real_r_sum, time_used,
                               number_of_pushes);
                    }
#endif
                }
                real_r_sum = one_minus_alpha * (r_sum - residuals[numOfVertices]) / (1 - residuals[numOfVertices]);
#ifdef SHOW_SERIES_POWER
                //                const double time_current = getCurrentTime();
                //                printf("r_sum:%.12f\tTime Used:%.4f\t#Pushes:%zu\n", real_r_sum, time_current - time_start,
                //                       number_of_pushes);
#endif
            }
#ifndef SHOW_SERIES_POWER
            //            const double time_end = getCurrentTime();
            //            printf("Sid: %d\tThe Sum of Residuals:%.16f\tTime Used:%.4f\t#Iter:%u\t#Pushes:%zu\n", sid, real_r_sum,
            //                   time_end - time_epoch_start, iter, number_of_pushes);
#endif
        }
        const double dummy_scale_factor = 1.0 / (1 - residuals[numOfVertices]);
#ifdef DEBUG_MODE
        {
            real_r_sum =
                    one_minus_alpha * (r_sum - residuals[numOfVertices]) * dummy_scale_factor;
            const double time_used = getCurrentTime() - time_start;
            printf("#Iter:**\tr_sum:%.12f\tTime Used:%.4f\t#Pushes:%zu\n", real_r_sum, time_used, number_of_pushes);
        }
#endif
        for (VertexIdType id = 0; id < numOfVertices; ++id) {
            pi[id] *= dummy_scale_factor;
            residuals[id] *= dummy_scale_factor;
            pi[id] += residuals[id] * alpha;
        }
    }

    void
    run3(const VertexIdType &sid, std::vector<PageRankScoreType> &pi, std::vector<PageRankScoreType> &residuals) const {
        if (pi.size() < numOfVertices + 1) {
            pi.resize(numOfVertices + 1, 0);
            residuals.resize(numOfVertices + 1, 0);
        }
        std::fill(pi.begin(), pi.end(), 0);
        std::fill(residuals.begin(), residuals.end(), 0);
        // reserve one slot for the dummy vertex
        MyQueue active_vertices(numOfVertices + 1);
        // reserve one slot for the dummy vertex
        std::vector<bool> is_active(numOfVertices + 1, false);
        active_vertices.push(sid);
        is_active[sid] = true;
        // this prevents dummy vertex from entering the queue
        is_active[numOfVertices] = true;
        double num_walks = ceil(numOfVertices * (2 + (2.0 / 3.0) * 0.5) * log(numOfVertices) / (0.5 * 0.5));
        double r_sum = num_walks;
        residuals[sid] = num_walks;
        ////////////////////////////////////////////////////////
        const auto avg_deg = static_cast<PageRankScoreType>(graph.getNumOfEdges() / (double) graph.getNumOfVertices());
        const VertexIdType queue_threshold = (numOfVertices / avg_deg * 4);
        ////////////////////////////////////////////////////////
        size_t number_of_pushes = 0;
        uint32_t num_active = 0;
        const PageRankScoreType increment_factor = std::exp(1.0);
//        const PageRankScoreType increment_factor = 9;
        ////////////////////////////////////////////////////////
        const double time_forward_start = getCurrentTime();
        ////////////////////////////////////////////////////////
        const uint32_t initial_size = std::max(num_walks / (1000 * log(numOfVertices)), 1.0);
        const uint32_t step_size = std::max(powf(initial_size, 1.0 / 3.0), 2.0f);
        for (uint32_t scale_factor = initial_size;
             scale_factor >= 1 && active_vertices.size() < queue_threshold;) {
//            const PageRankScoreType scale_factor_over_one_minus_alpha = scale_factor / one_minus_alpha;
            while (!active_vertices.empty() && active_vertices.size() < queue_threshold) {
                const VertexIdType id = active_vertices.front();
                active_vertices.pop();
                is_active[id] = false;
                const PageRankScoreType residual = residuals[id];
                const VertexIdType &idx_start = graph.get_neighbor_list_start_pos(id);
                const VertexIdType &idx_end = graph.get_neighbor_list_start_pos(id + 1);
                const PageRankScoreType degree_f = idx_end - idx_start;
                const PageRankScoreType one_minus_alpha_residual = one_minus_alpha * residual;
                if (one_minus_alpha_residual >= degree_f * scale_factor) {
                    const PageRankScoreType alpha_residual = residual - one_minus_alpha_residual;
                    pi[id] += alpha_residual;
                    r_sum -= alpha_residual;
                    residuals[id] = 0;
                    const PageRankScoreType increment = one_minus_alpha_residual / degree_f;
#ifdef SHOW_SERIES_POWER
                    number_of_pushes += degree_f;
#endif
                    for (uint32_t j = idx_start; j < idx_end; ++j) {
                        const VertexIdType &nid = graph.getOutNeighbor(j);
                        residuals[nid] += increment;
                        if (!is_active[nid]) {
                            active_vertices.push(nid);
                            is_active[nid] = true;
                        }
                    }
                }
            }
            scale_factor /= step_size;
            if (active_vertices.empty()) {
                for (VertexIdType id = 0; id < numOfVertices; ++id) {
                    if (one_minus_alpha * residuals[id] >= scale_factor) {
                        active_vertices.push(id);
                        is_active[id] = true;
                    }
                }
            }
        }
#ifdef SHOW_SERIES_POWER
        {
            const double time_used_current = getCurrentTime() - time_forward_start;
            const PageRankScoreType one_over_num_walk = 1.0 / num_walks;
            const PageRankScoreType real_r_sum =
                    (r_sum - residuals[numOfVertices]) * one_over_num_walk /
                    (1 - residuals[numOfVertices] * one_over_num_walk);
            printf("r_sum:%.12f\tTime Used:%.4f\t#Pushes:%zu\n", real_r_sum, time_used_current, number_of_pushes);
        }
#endif
#ifdef SHOW_TIME_POWER
        //        const double time_forward_end = getCurrentTime();
        //        MSG(time_forward_end - time_forward_start);
#endif
        ////////////////////////////////////////////////////////
        num_active = active_vertices.size();
        for (bool more_pushes = true; more_pushes;) {
#ifdef SHOW_TIME_POWER
            const double time_round_start = getCurrentTime();
#endif
            for (; num_active > queue_threshold;) {
                num_active = 0;
                for (VertexIdType id = 0, next_id = 1, degree, idx_start,
                             idx_end = graph.get_neighbor_list_start_pos(id);
                     id < numOfVertices; ++id, ++next_id) {
                    idx_start = idx_end;
                    idx_end = graph.get_neighbor_list_start_pos(next_id);
                    degree = idx_end - idx_start;
                    const PageRankScoreType &residual = residuals[id];
                    if (residual >= degree) {
                        const PageRankScoreType alpha_residual = alpha * residual;
                        pi[id] += alpha_residual;
                        r_sum -= alpha_residual;
                        const PageRankScoreType increment = (residual - alpha_residual) / degree;
                        residuals[id] = 0;
                        num_active += degree;
                        for (uint32_t j = idx_start; j < idx_end; ++j) {
                            const VertexIdType &nid = graph.getOutNeighbor(j);
                            residuals[nid] += increment;
                        }
                    }
                }
#ifdef SHOW_SERIES_POWER
                number_of_pushes += num_active;
#endif
#ifdef SHOW_DEBUG_POWER
                printf("Active Size In Linear Push:%d\n", num_active);
#endif
            }

#ifdef SHOW_TIME_POWER
            const double time_power_end = getCurrentTime();
            const double _time_for_power_ = time_power_end - time_round_start;
            MSG(_time_for_power_);
#endif
            active_vertices.clear();
            std::fill(is_active.begin(), is_active.end(), false);
            // this prevents dummy vertex from entering the queue
            is_active[numOfVertices] = true;
            for (VertexIdType id = 0; id < numOfVertices; ++id) {
                if (residuals[id] >= 1) {
                    active_vertices.push(id);
                    is_active[id] = true;
                }
            }
            while (!active_vertices.empty()) {
                const VertexIdType id = active_vertices.front();
                active_vertices.pop();
                is_active[id] = false;
                const VertexIdType &idx_start = graph.get_neighbor_list_start_pos(id);
                const VertexIdType &idx_end = graph.get_neighbor_list_start_pos(id + 1);
                const auto degree_f = static_cast<PageRankScoreType>(idx_end - idx_start);
                const PageRankScoreType &residual = residuals[id];
                if (residual >= degree_f) {
                    const PageRankScoreType alpha_residual = alpha * residual;
                    pi[id] += alpha_residual;
                    r_sum -= alpha_residual;
                    const PageRankScoreType increment = (residual - alpha_residual) / degree_f;
                    residuals[id] = 0;
#ifdef SHOW_SERIES_POWER
                    number_of_pushes += degree_f;
#endif
                    for (uint32_t j = idx_start; j < idx_end; ++j) {
                        const VertexIdType &nid = graph.getOutNeighbor(j);
                        residuals[nid] += increment;
                        if (!is_active[nid]) {
                            active_vertices.push(nid);
                            is_active[nid] = true;
                        }
                    }
                }
            }
#ifdef SHOW_TIME_POWER
            const double time_forward_2_end = getCurrentTime();
            const double time_for_forward = time_forward_2_end - time_power_end;
            MSG(time_for_forward);
#endif
//            const PageRankScoreType r_sum = std::accumulate(residuals.begin(), residuals.end(), 0.0);
            const double time_used_current = getCurrentTime() - time_forward_start;
            const PageRankScoreType one_over_num_walk = 1.0 / num_walks;
            const PageRankScoreType real_r_sum =
                    one_minus_alpha * (r_sum - residuals[numOfVertices]) * one_over_num_walk /
                    (1 - residuals[numOfVertices] * one_over_num_walk);
#ifdef SHOW_SERIES_POWER
            printf("r_sum:%.12f\tTime Used:%.4f\t#Pushes:%zu\n", real_r_sum, time_used_current, number_of_pushes);

#endif
#ifdef SHOW_DEBUG_POWER
            std::cout << "The number of operations: " << number_of_pushes << std::endl;
#endif
            num_active = 0;
            if (real_r_sum > l1_error) {
                more_pushes = true;
                num_walks *= increment_factor;
                for (VertexIdType id = 0; id < numOfVertices; ++id) {
                    if (pi[id] || residuals[id]) {
                        pi[id] *= increment_factor;
                        residuals[id] *= increment_factor;
                        if (residuals[id] >= 1) {
                            ++num_active;
                        }
                    }
                }
                pi[numOfVertices] *= increment_factor;
                residuals[numOfVertices] *= increment_factor;
                r_sum *= increment_factor;
#ifdef SHOW_TIME_POWER
                //                const double time_update_end = getCurrentTime();
                //                const double time_for_update_ = time_update_end - time_forward_2_end;
                //                MSG(time_for_update_);
#endif
            } else {
                more_pushes = false;
            }
#ifdef SHOW_TIME_POWER
            const double time_current = getCurrentTime();
            const double time_this_round_ = time_current - time_round_start;
            const double time_used_so_far = time_current - time_forward_start;
            MSG(time_this_round_)
            MSG(time_used_so_far)
            printf("---\n");
#endif
        }
        const PageRankScoreType one_over_num_walks = (1.0 / num_walks);
        MSG(pi[numOfVertices])
        MSG(residuals[graph.get_dummy_id()] * one_over_num_walks)
        const double dummy_scale_factor = 1.0 / (1 - residuals[graph.get_dummy_id()] * one_over_num_walks);
        printf("Sid: %d\tDummy Scale Factor: %.12f\n", sid, dummy_scale_factor);
        for (VertexIdType id = 0; id < numOfVertices; ++id) {
            pi[id] *= (dummy_scale_factor * one_over_num_walks);
            residuals[id] *= (dummy_scale_factor * one_over_num_walks);
        }
        for (VertexIdType id = 0; id < numOfVertices; ++id) {
            pi[id] += residuals[id] * alpha;
        }
    }

private:
    static void run_ground_truth_iteration(IteratedMethods *pIteratedMethods, const std::vector<VertexIdType> _sids) {
        pIteratedMethods->run_ground_truth(_sids);
    }

public:

    void multi_thread_ground_truth_iteration(const std::vector<VertexIdType> &_sids,
                                             uint32_t _num_threads = std::thread::hardware_concurrency()) {
        uint32_t num_of_source = _sids.size();
        if (num_of_source == 0) {
            printf("Error in multi_thread_ground_truth_iteration. Zero Source Ids Provided.\n");
            exit(1);
        }
        _num_threads = std::min(_num_threads, num_of_source);
        printf("The number of threads: %d\n", _num_threads);
        std::vector<std::thread> threads;

        uint32_t batch_size = std::ceil(num_of_source / (float) _num_threads);
        for (uint32_t thread_id = 0, start = 0, end;
             thread_id < _num_threads && start < num_of_source; ++thread_id) {
            end = std::min(start + batch_size, num_of_source);
            std::vector<VertexIdType> thread_sids(_sids.begin() + start, _sids.begin() + end);
            threads.emplace_back(IteratedMethods::run_ground_truth_iteration, this, thread_sids);
            start += batch_size;
        }
        for (auto &thread : threads) {
            thread.join();
        }
    }
};


#endif //SPEEDPPR_ITERATEDMETHODS_H
