/**
 * @file global.h, part of the project TransOMS under Apache License 2.0
 * @author jdlph (jdlph@hotmail.com) and xzhou99 (xzhou74@asu.edu)
 * @brief Alias, constants, and forward declarations of some classes.
 *
 * @copyright Copyright (c) 2023 Peiheng Li, Ph.D. and Xuesong (Simon) Zhou, Ph.D.
 */

#ifndef GUARD_GLOBAL_H
#define GUARD_GLOBAL_H

#define MLC_DEQUE

#include <random>
#include <string>
#include <tuple>

namespace transoms
{
using size_type = unsigned int;
// origin zone no, destination zone no, demand period no, agent type no
using ColumnVecKey = std::tuple<size_type, size_type, unsigned short, unsigned short>;

// some constants
constexpr unsigned short CHUNK = 256;
constexpr unsigned short COORD_X = 91;
constexpr unsigned short COORD_Y = 181;
constexpr unsigned short BACKWAVE_SPEED = 12;
constexpr unsigned short JAM_DENSITY = 200;
constexpr unsigned short MINUTES_IN_HOUR = 60;
constexpr unsigned short SECONDS_IN_MINUTE = 60;
constexpr unsigned short SECONDS_IN_HOUR = 3600;

constexpr double EPSILON = 0.00001;
constexpr double MIN_COL_VOL = 0.1;

const std::string ALL_MODES {"all"};
const std::string AT_DEFAULT_NAME {"passenger"};
const std::string AT_LEGACY_NAME {"auto"};

// forward declarations of supply classes to be used in demand.h
class Column;
class SpecialEvent;

// forward declarations of demand classes to be used in supply.h
class Agent;
class AgentType;
class DemandPeriod;

#ifndef MLC_DEQUE
// a struct for heap-Dijkstra's algorithm
struct HeapNode {
    HeapNode(size_type node_no_, double cost_) : node_no {node_no_}, cost {cost_}
    {
    }

    friend bool operator<(const HeapNode& h1, const HeapNode& h2)
    {
        return h1.cost > h2.cost;
    }

    size_type node_no;
    double cost;
};
#endif

template<typename T>
T uniform(T lb, T ub)
{
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<> dis(lb, ub);

    return dis(gen);
}

} // namespace transoms

#endif