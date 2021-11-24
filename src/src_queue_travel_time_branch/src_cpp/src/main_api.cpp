/* Portions Copyright 2010-2021 Xuesong Zhou and Peiheng Li
 *
 * If you help write or modify the code, please also list your names here.
 * The reason of having Copyright info here is to ensure all the modified version, as a whole, under the GPL
 * and further prevent a violation of the GPL.
 *
 * More about "How to use GNU licenses for your own software"
 * http://www.gnu.org/licenses/gpl-howto.html
 */

// Peiheng, 02/03/21, remove them later after adopting better casting
#pragma warning(disable : 4305 4267 4018)
// stop warning: "conversion from 'int' to 'float', possible loss of data"
#pragma warning(disable: 4244)



#ifdef _WIN32
#include "pch.h"
#endif

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cstring>
#include <cstdio>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <functional>
#include <stack>
#include <list>
#include <vector>
#include <map>
#include <omp.h>
#include "config.h"
#include "utils.h"

using std::max;
using std::min;
using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::map;
using std::ifstream;
using std::ofstream;
using std::istringstream;

// some basic parameters setting

//Pls make sure the _MAX_K_PATH > Agentlite.cpp's g_number_of_column_generation_iterations+g_reassignment_number_of_K_paths and the _MAX_ZONE remain the same with .cpp's defination
constexpr auto _MAX_LABEL_COST = 1.0e+15;

constexpr auto _MAX_AGNETTYPES = 10; //because of the od demand store format,the MAX_demandtype must >=g_DEMANDTYPES.size()+1;
constexpr auto _MAX_TIMEPERIODS = 5; // time period set to 4: mid night, morning peak, mid-day and afternoon peak;

constexpr auto _MAX_MEMORY_BLOCKS = 100;

constexpr auto _MAX_LINK_SIZE_IN_A_PATH = 5000;		// lu
constexpr auto _MAX_LINK_SIZE_FOR_A_NODE = 200;

constexpr auto _MAX_TIMESLOT_PerPeriod = 100; // max 96 15-min slots per day
constexpr auto _MAX_TIMEINTERVAL_PerDay = 300; // max 96*3 5-min slots per day
constexpr auto _MAX_DAY_PerYear = 360; // max 96*3 5-min slots per day

constexpr auto _default_saturation_flow_rate = 1530;

constexpr auto MIN_PER_TIMESLOT = 15;

extern class CTMCLink;
extern std::map<int, int> g_dayDataMap;
/* make sure we change the following two parameters together*/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
constexpr auto number_of_seconds_per_interval = 1;
constexpr auto number_of_interval_per_min = 60;

/* number_of_seconds_per_interval should satisify the ratio of 60/number_of_seconds_per_interval is an integer*/

// Linear congruential generator
constexpr auto LCG_a = 17364;
constexpr auto LCG_c = 0;
constexpr auto LCG_M = 65521;  // it should be 2^32, but we use a small 16-bit number to save memory

constexpr auto STRING_LENGTH_PER_LINE = 50000;

char str_buffer[STRING_LENGTH_PER_LINE];

void g_createNode2ZoneMappingInQEM();
// FILE* g_pFileOutputLog = nullptr;

template <typename T>
T** AllocateDynamicArray(int nRows, int nCols)
{
    T** dynamicArray;

    dynamicArray = new (std::nothrow) T*[nRows];

    if (!dynamicArray)
    {
        dtalog.output() << "Error: insufficient memory.";
        cout << "Error: insufficient memory.";
        g_ProgramStop();
    }

    for (int i = 0; i < nRows; ++i)
    {
        dynamicArray[i] = new (std::nothrow) T[nCols];

        if (!dynamicArray[i])
        {
            dtalog.output() << "Error: insufficient memory.";
            cout <<  "Error: insufficient memory.";
            g_ProgramStop();
        }
    }

    return dynamicArray;
}

template <typename T>
void DeallocateDynamicArray(T** dArray, int nRows, int nCols)
{
    if (!dArray)
        return;

    for (int x = 0; x < nRows; ++x)
        delete[] dArray[x];

    delete[] dArray;
}

template <typename T>
T*** Allocate3DDynamicArray(int nX, int nY, int nZ)
{
    T*** dynamicArray = new (std::nothrow) T**[nX];

    if (!dynamicArray)
    {
        dtalog.output() << "Error: insufficient memory.";
        cout << "Error: insufficient memory.";
        g_ProgramStop();
    }

    for (int x = 0; x < nX; ++x)
    {
        if (x % 1000 == 0)
        {
            dtalog.output() << "allocating 3D memory for " << x << endl;
            cout << "allocating 3D memory for " << x << endl;
        }

        dynamicArray[x] = new (std::nothrow) T*[nY];

        if (!dynamicArray[x])
        {
            dtalog.output() << "Error: insufficient memory.";
            cout << "Error: insufficient memory.";
            g_ProgramStop();
        }

        for (int y = 0; y < nY; ++y)
        {
            dynamicArray[x][y] = new (std::nothrow) T[nZ];
            if (!dynamicArray[x][y])
            {
                dtalog.output() << "Error: insufficient memory.";
                cout << "Error: insufficient memory.";
                g_ProgramStop();
            }
        }
    }

    for (int x = 0; x < nX; ++x)
        for (int y = 0; y < nY; ++y)
            for (int z = 0; z < nZ; ++z)
                dynamicArray[x][y][z] = 0;

    return dynamicArray;
}

template <typename T>
void Deallocate3DDynamicArray(T*** dArray, int nX, int nY)
{
    if (!dArray)
        return;

    for (int x = 0; x < nX; ++x)
    {
        for (int y = 0; y < nY; ++y)
            delete[] dArray[x][y];

        delete[] dArray[x];
    }

    delete[] dArray;
}

template <typename T>
T**** Allocate4DDynamicArray(int nM, int nX, int nY, int nZ)
{
    T**** dynamicArray = new (std::nothrow) T***[nX];

    if (!dynamicArray)
    {
        dtalog.output() << "Error: insufficient memory.";
        cout << "Error: insufficient memory.";
        g_ProgramStop();
    }

    for (int m = 0; m < nM; ++m)
    {
        if (m % 1000 == 0)
        { 
            dtalog.output() << "allocating 4D memory for " << m << " zones" << endl;
            cout << "allocating 4D memory for " << m << " zones" << endl;
        }

        dynamicArray[m] = new (std::nothrow) T**[nX];

        if (!dynamicArray[m])
        {
            dtalog.output() << "Error: insufficient memory.";
            cout << "Error: insufficient memory.";
            g_ProgramStop();
        }

        for (int x = 0; x < nX; ++x)
        {
            dynamicArray[m][x] = new (std::nothrow) T*[nY];

            if (!dynamicArray[m][x])
            {
                dtalog.output() << "Error: insufficient memory.";
                cout << "Error: insufficient memory.";
                g_ProgramStop();
            }

            for (int y = 0; y < nY; ++y)
            {
                dynamicArray[m][x][y] = new (std::nothrow) T[nZ];
                if (!dynamicArray[m][x][y])
                {
                    dtalog.output() << "Error: insufficient memory.";
                    cout << "Error: insufficient memory.";
                    g_ProgramStop();
                }
            }
        }
    }

    return dynamicArray;
}

template <typename T>
void Deallocate4DDynamicArray(T**** dArray, int nM, int nX, int nY)
{
    if (!dArray)
        return;

    for (int m = 0; m < nM; ++m)
    {
        for (int x = 0; x < nX; ++x)
        {
            for (int y = 0; y < nY; ++y)
                delete[] dArray[m][x][y];

            delete[] dArray[m][x];
        }

        delete[] dArray[m];
    }

    delete[] dArray;
}


__int64 g_GetCellID(double x, double y, double grid_resolution)
{
    __int64 xi;
    xi = floor(x / grid_resolution);

    __int64 yi;
    yi = floor(y / grid_resolution);

    __int64 x_code, y_code, code;
    x_code = fabs(xi) * grid_resolution * 1000000000000;
    y_code = fabs(yi)* grid_resolution * 100000;
    code = x_code + y_code;
    return code;
};



class CDemand_Period {
public:
    CDemand_Period() : demand_period{ 0 }, starting_time_slot_no{ 0 }, ending_time_slot_no{ 0 }, demand_period_id{ 1 }, time_period{"0700_0800"}
    {
        
    }

    int get_time_horizon_in_min()
    {
        return (ending_time_slot_no - starting_time_slot_no) * 15;
    }

    string demand_period;
    int starting_time_slot_no;
    int ending_time_slot_no;
    string time_period;
    int demand_period_id;
};

class CAgent_type {
public:
    CAgent_type() : agent_type_no{ 1 }, value_of_time{ 1 }, agent_type{ "auto" }, PCE{ 1.0 }, avg_speed{ 60 }, trip_time_budget_in_min{ 25 }, trip_ratio{ 1.0 }
    {

    }
    float avg_speed;
    float trip_ratio;
    int agent_type_no;
    // dollar per hour
    float value_of_time;
    // link type, product consumption equivalent used, for travel time calculation
    float PCE;
    float trip_time_budget_in_min;
    string agent_type;
};

class CLinkType
{
public:
    CLinkType() : link_type{ 1 }, number_of_links{ 0 }, traffic_flow_code{ 0 }
    {
    }


    int link_type;
    int number_of_links;
    int traffic_flow_code;

    string link_type_name;
    string type_code;
};

class CColumnPath {
public:
    CColumnPath() : path_node_vector{ nullptr }, path_link_vector{ nullptr }, path_seq_no{ 0 },
        path_switch_volume{ 0 }, path_volume{ 0 }, path_travel_time{ 0 }, path_distance{ 0 }, path_toll{ 0 },
        path_gradient_cost{ 0 }, path_gradient_cost_difference{ 0 }, path_gradient_cost_relative_difference{ 0 }
    {
    }

    ~CColumnPath()
    {
        if (m_node_size >= 1)
        {
            delete[] path_node_vector;
            delete[] path_link_vector;
        }
    }

    // Peiheng, 02/02/21, consider using move later
    void AllocateVector(int node_size, const int* node_vector, int link_size, const int* link_vector, bool backwardflag = false)
    {
        m_node_size = node_size;
        m_link_size = link_size;
        // dynamic array
        path_node_vector = new int[node_size];
        path_link_vector = new int[link_size];

        if(backwardflag)
        {
            // copy backward
            for (int i = 0; i < m_node_size; ++i)
                path_node_vector[i] = node_vector[m_node_size - 1 - i];

            for (int i = 0; i < m_link_size; ++i)
                path_link_vector[i] = link_vector[m_link_size - 1 - i];
        }
        else
        {
            // copy forward
            for (int i = 0; i < m_node_size; ++i)
                path_node_vector[i] = node_vector[i];

            for (int i = 0; i < m_link_size; ++i)
                path_link_vector[i] = link_vector[i];
        }
    }

    // Peiheng, 02/02/21, consider using move later
    void AllocateVector(const std::vector<int>& node_vector, const std::vector<int>& link_vector, bool backwardflag = false)
    {
        m_node_size = node_vector.size();
        m_link_size = link_vector.size();
        // dynamic array
        path_node_vector = new int[m_node_size];
        path_link_vector = new int[m_link_size];

        if(backwardflag)
        {
            // copy backward
            for (int i = 0; i < m_node_size; ++i)
                path_node_vector[i] = node_vector[m_node_size - 1 - i];

            for (int i = 0; i < m_link_size; ++i)
                path_link_vector[i] = link_vector[m_link_size - 1 - i];
        }
        else
        {
            // copy forward
            for (int i = 0; i < m_node_size; ++i)
                path_node_vector[i] = node_vector[i];

            for (int i = 0; i < m_link_size; ++i)
                path_link_vector[i] = link_vector[i];
        }
    }

    int* path_node_vector;
    int* path_link_vector;

    int path_seq_no;
    // path volume
    double path_volume;
	double path_switch_volume;
	double path_travel_time;
	double path_distance;
	double path_toll;
    // first order graident cost.
	double path_gradient_cost;
    // first order graident cost - least gradient cost
	double path_gradient_cost_difference;
    // first order graident cost - least gradient cost
	double path_gradient_cost_relative_difference;

    int m_node_size;
    int m_link_size;
    std::vector<int> agent_simu_id_vector;
};

class CAgentPath {
public:
    CAgentPath() : path_id{ 0 }, node_sum{ -1 }, travel_time{ 0 }, distance{ 0 }, volume{ 0 }
    {
    }

    int path_id;
    int node_sum;
    float travel_time;
    float distance;
    float volume;

    int o_node_no;
    int d_node_no;

    std::vector <int> path_link_sequence;
};

class CColumnVector {

public:
    // this is colletion of unique paths
    CColumnVector() : cost{ 0 }, time{ 0 }, distance{ 0 }, od_volume{ 0 }, bfixed_route{ false }
    {
    }

    float cost;
    float time;
    float distance;
    // od volume
    double od_volume;
    bool bfixed_route;
    // first key is the sum of node id;. e.g. node 1, 3, 2, sum of those node ids is 6, 1, 4, 2 then node sum is 7.
    // Peiheng, 02/02/21, potential memory leak, fix it
    std::map <int, CColumnPath> path_node_sequence_map;
};

class CAgent_Column {
public:
    CAgent_Column() : cost{ 0 }
    {
    }

    float cost;
    float volume;
    float travel_time;
    float distance;

    int agent_id;
    int o_zone_id;
    int d_zone_id;
    int o_node_id;
    int d_node_id;

    string agent_type;
    string demand_period;

    vector<int> path_node_vector;
    vector<int> path_link_vector;
    vector<float> path_time_vector;
};

// event structure in this "event-based" traffic simulation

class DTAVehListPerTimeInterval
{
public:
    std::vector<int> m_AgentIDVector;
};

std::map<int, DTAVehListPerTimeInterval> g_AgentTDListMap;

class CAgent_Simu
{
public:
    CAgent_Simu() : agent_vector_seq_no{ -1 }, path_toll{ 0 }, departure_time_in_min{0}, m_bGenereated{ false }, m_bCompleteTrip{ false },
        m_Veh_LinkArrivalTime_in_simu_interval{ nullptr }, m_Veh_LinkDepartureTime_in_simu_interval{ nullptr }
    {
    }

    ~CAgent_Simu()
    {
        DeallocateMemory();
    }

    void AllocateMemory()
    {
        if (!m_Veh_LinkArrivalTime_in_simu_interval)
        {
            m_current_link_seq_no = 0;
            m_Veh_LinkArrivalTime_in_simu_interval = new int[path_link_seq_no_vector.size()];
            m_Veh_LinkDepartureTime_in_simu_interval = new int[path_link_seq_no_vector.size()];

            for (int i = 0; i < path_link_seq_no_vector.size(); ++i)
            {
                m_Veh_LinkArrivalTime_in_simu_interval[i] = -1;
                m_Veh_LinkDepartureTime_in_simu_interval[i] = -1;
            }

            m_path_link_seq_no_vector_size = path_link_seq_no_vector.size();
            departure_time_in_simu_interval = int(departure_time_in_min * 60.0 / number_of_seconds_per_interval + 0.5);  // round off
        }
    }

    void DeallocateMemory()
    {
        // Peiheng, 02/02/21, it is OK to delete nullptr
        if (m_Veh_LinkArrivalTime_in_simu_interval)
            delete[] m_Veh_LinkArrivalTime_in_simu_interval;
        if (m_Veh_LinkDepartureTime_in_simu_interval)
            delete[] m_Veh_LinkDepartureTime_in_simu_interval;

        m_Veh_LinkArrivalTime_in_simu_interval = nullptr;
        m_Veh_LinkDepartureTime_in_simu_interval = nullptr;
    }

    float GetRandomRatio()
    {
        // Peiheng, 02/02/21, m_RandomSeed is uninitialized
        //m_RandomSeed is automatically updated.
        m_RandomSeed = (LCG_a * m_RandomSeed + LCG_c) % LCG_M;

        return float(m_RandomSeed) / LCG_M;
    }

    int agent_vector_seq_no;
    float path_toll;
    float departure_time_in_min;
    bool m_bGenereated;
    bool m_bCompleteTrip;
    int* m_Veh_LinkArrivalTime_in_simu_interval;
    int* m_Veh_LinkDepartureTime_in_simu_interval;

    int agent_service_type;
    int demand_type;
    int agent_id;

    int m_current_link_seq_no;
    int m_path_link_seq_no_vector_size;

    int departure_time_in_simu_interval;
    float arrival_time_in_min;

    unsigned int m_RandomSeed;

    // external input
    std::vector<int> path_link_seq_no_vector;
    // internal output
    std::vector<float> time_seq_no_vector;
    std::vector<int> path_timestamp_vector;
    std::vector<int> path_node_id_vector;
};

vector<CAgent_Simu*> g_agent_simu_vector;

class Assignment {
public:
    // default is UE
    Assignment() : assignment_mode{ 0 }, g_number_of_memory_blocks{ 8 }, g_number_of_threads{ 1 }, g_link_type_file_loaded{ true }, g_agent_type_file_loaded{ false },
        total_demand_volume{ 0.0 }, g_origin_demand_array{ nullptr }, g_column_pool{ nullptr }, g_number_of_in_memory_simulation_intervals{ 500 },
        g_number_of_column_generation_iterations{ 20 }, g_number_of_demand_periods{ 24 }, g_number_of_links{ 0 }, g_number_of_timing_arcs{ 0 },
        g_number_of_nodes{ 0 }, g_number_of_zones{ 0 }, g_number_of_agent_types{ 0 }, g_reassignment_tau0{ 999 }, debug_detail_flag{ 1 }, g_subarea_mode{ 0 }, b_unit_mile_mode{false},
        m_GridResolution {0.01}, accessibility_output {1}

    {
    }

    ~Assignment()
    {
        if (g_column_pool)
            Deallocate4DDynamicArray(g_column_pool, g_number_of_zones, g_number_of_zones, g_number_of_agent_types);

        if (g_origin_demand_array)
            Deallocate3DDynamicArray(g_origin_demand_array, g_number_of_zones, g_number_of_agent_types);

        DeallocateLinkMemory4Simulation();
    }

    void InitializeDemandMatrix(int number_of_zones, int number_of_agent_types, int number_of_time_periods)
    {
        total_demand_volume = 0.0;
        g_number_of_zones = number_of_zones;
        g_number_of_agent_types = number_of_agent_types;

        g_column_pool = Allocate4DDynamicArray<CColumnVector>(number_of_zones, number_of_zones, max(1, number_of_agent_types), number_of_time_periods);
        g_origin_demand_array = Allocate3DDynamicArray<float>(number_of_zones, max(1, number_of_agent_types), number_of_time_periods);

        for (int i = 0; i < number_of_zones; ++i)
        {
            for (int at = 0;at < number_of_agent_types; ++at)
            {
                for (int tau = 0;tau < g_number_of_demand_periods; ++tau)
                    g_origin_demand_array[i][at][tau] = 0.0;
            }
        }

        for (int i = 0; i < number_of_agent_types; ++i)
        {
            for (int tau = 0;tau < g_number_of_demand_periods;++tau)
                total_demand[i][tau] = 0.0;
        }

        g_DemandGlobalMultiplier = 1.0f;
    }

    int get_in_memory_time(int t)
    {
        return t % g_number_of_in_memory_simulation_intervals;
    }

    void STTrafficSimulation();
    //OD demand estimation estimation
    void Demand_ODME(int OD_updating_iterations);
    void Mapping_TMC_Identification();

    std::map<string,int> m_TMClink_map;
    std::map<string, int> m_TMC_corridor_map;

    bool Map_TMC_Reading();

    void AllocateLinkMemory4Simulation();
    void DeallocateLinkMemory4Simulation();

    double m_GridResolution;
    bool b_unit_mile_mode;

    int g_subarea_mode;  // 1: subarea cuts: 2: focusing approach to generate aggregated zones outside subarea with reduced od demand
    std::vector<double> m_subarea_vec_x, m_subarea_vec_y;
    std::map<int,int> m_subarea_node_flag_map;  // inside the subarea
    std::map<int, int> m_subarea_node_id_map;   //to be include in the subarea network
    std::map<int, int> m_subarea_boundary_node_map;   //to be include in the subarea network

    int assignment_mode;
    int g_number_of_memory_blocks;
    int g_number_of_threads;
    int accessibility_output;

    bool g_link_type_file_loaded;
    bool g_agent_type_file_loaded;

    float total_demand_volume;
    float*** g_origin_demand_array;
    CColumnVector**** g_column_pool;

    // the data horizon in the memory
    int g_number_of_in_memory_simulation_intervals;
    int g_number_of_column_generation_iterations;
    int g_number_of_demand_periods;

    int g_number_of_links;
    int g_number_of_timing_arcs;
    int g_number_of_nodes;
    int g_number_of_zones;
    int g_number_of_agent_types;

    // Peiheng, 02/06/21, useless members
    int g_reassignment_tau0;
    int debug_detail_flag;

    // hash table, map external node number to internal node sequence no.
    std::map<int, int> g_node_id_to_seq_no_map;
    // from integer to integer map zone_id to zone_seq_no
    std::map<int, int> g_zoneid_to_zone_seq_no_mapping;
    std::map<string, int> g_link_id_map;

    std::vector<CDemand_Period> g_DemandPeriodVector;
    int g_LoadingStartTimeInMin;
    int g_LoadingEndTimeInMin;



    std::vector<CAgent_type> g_AgentTypeVector;
    void AddAgentType(string agent_type_str, float VOT, float PCE, float AvgSpeed, float TripTimeBudgetInMin, float TripRatio)
    {
        CAgent_type agent_type;
        agent_type.agent_type_no = g_AgentTypeVector.size();
        agent_type.avg_speed = AvgSpeed;
        agent_type.agent_type = agent_type_str;
        agent_type.value_of_time = VOT;
        agent_type.PCE = PCE;
        agent_type.trip_time_budget_in_min = TripTimeBudgetInMin;
        agent_type.trip_ratio = TripRatio;
        agent_type_2_seqno_mapping[agent_type.agent_type] = g_AgentTypeVector.size();
        g_AgentTypeVector.push_back(agent_type);
        g_number_of_agent_types = g_AgentTypeVector.size();

    }

    std::map<int, int> zone_id_to_centriod_node_no_mapping;  // this is an one-to-one mapping
    std::map<int, int> zone_id_2_node_no_mapping;  // this is used to mark if this zone_id has been identified or not
    std::map<int, _int64> zone_id_2_cell_id_mapping;  // this is used to mark if this zone_id has been identified or not
    std::map<_int64, int> cell_id_mapping;  // this is used to mark if this cell_id has been identified or not

    std::map<int, CLinkType> g_LinkTypeMap;


    std::map<string, int> demand_period_to_seqno_mapping;
    std::map<string, int> agent_type_2_seqno_mapping;

    float total_demand[_MAX_AGNETTYPES][_MAX_TIMEPERIODS];
    float g_DemandGlobalMultiplier;

    // used in ST Simulation
    float** m_LinkOutFlowCapacity;

    // in min
    float** m_LinkTDWaitingTime;
    // number of simulation time intervals
    int** m_LinkTDTravelTime;

    float** m_LinkCumulativeArrival;
    float** m_LinkCumulativeDeparture;

    int g_start_simu_interval_no;
    int g_number_of_simulation_intervals;
    // is shorter than g_number_of_simulation_intervals
    int g_number_of_loading_intervals;
    // the data horizon in the memory in min
    int g_number_of_simulation_horizon_in_min;
};

Assignment assignment;

class CVDF_Period
{
public:
    CVDF_Period() : m{ 0.5 }, VOC{ 0 }, gamma{ 3.47f }, mu{ 1000 }, PHF{ 3 },
		alpha{ 0.15f }, beta{ 4 }, rho{ 1 }, preload{ 0 }, penalty{ 0 }, marginal_base{ 1 },
        starting_time_slot_no{ 0 }, ending_time_slot_no{ 0 },
        cycle_length{ 29 }, red_time{ 0 }, t0{ 0 }, t3{ 0 }
    {
		for (int at = 0; at < _MAX_AGNETTYPES; at++)
		{
			toll[at] = 0;
            pce[at] = 0;
		}

		for (int t = 0; t < _MAX_TIMESLOT_PerPeriod; ++t)
        {
            Queue[t] = 0;
            waiting_time[t] = 0;
            arrival_rate[t] = 0;
            discharge_rate[t] = 0;
            travel_time[t] = 0;
        }
    }



    ~CVDF_Period()
    {
    }

    float get_waiting_time(int relative_time_slot_no)
    {
        if (relative_time_slot_no >=0 && relative_time_slot_no < _MAX_TIMESLOT_PerPeriod)
            return waiting_time[relative_time_slot_no];
        else
            return 0;
    }

    float PerformSignalVDF(float hourly_per_lane_volume, float red, float cycle_length)
    {
        float lambda = hourly_per_lane_volume;
        float mu = _default_saturation_flow_rate; //default saturation flow ratesa
        float s_bar = 1.0 / 60.0 * red * red / (2*cycle_length); // 60.0 is used to convert sec to min
        float uniform_delay = s_bar / max(1 - lambda / mu, 0.1f);

        return uniform_delay;
    }

	double PerformBPR(double volume)
    {
        // take nonnegative values
        volume = max(0.0, volume);

        // Peiheng, 02/02/21, useless block
        if (volume > 1.0)
        {
            int debug = 1;
        }

        VOC = volume / max(0.00001, capacity);
        avg_travel_time = FFTT + FFTT * alpha * pow((volume+preload) / max(0.00001, capacity), beta);
		total_travel_time = (volume + preload)*avg_travel_time;
        marginal_base = FFTT * alpha * beta*pow((volume + preload) / max(0.00001, capacity), beta - 1);

        return avg_travel_time;
        // volume --> avg_traveltime
    }

    // input period based volume
	double PerformBPR_X(double volume)
    {
        bValidQueueData = false;
        congestion_period_P = 0;

        float FFTT_in_hour = FFTT / 60.0;
        // convert avg_travel_time from unit of min to hour
        float avg_travel_time_in_hour = avg_travel_time / 60.0;

        // Step 1: Initialization
        int L = ending_time_slot_no - starting_time_slot_no;  // in 15 min slot
        if (L >= _MAX_TIMESLOT_PerPeriod - 1)
            return 0;

        for (int t = starting_time_slot_no; t <= ending_time_slot_no; ++t)
        {
            Queue[t] = 0;
            waiting_time[t] = 0;
            arrival_rate[t] = 0;
            discharge_rate[t]= mu/2.0;
            travel_time[t] = FFTT_in_hour;
        }

        // avg_travel time should be per min
        // avg_travel_time = (FFTT_in_hour + avg_waiting_time)*60.0;

        //int L = ending_time_slot_no - starting_time_slot_no;  // in 15 min slot
        float mid_time_slot_no = starting_time_slot_no + L / 2.0;  // t1;  // we can discuss thi
        // Case 1: fully uncongested region
        if (volume <= L * mu / 2)
        {
            // still keep 0 waiting time for all time period
            congestion_period_P = 0;
        }
        else
        {
            // partially congested region
            //if (volume > L * mu / 2 ) // Case 2
            float P = 0;
            // unit: hour  // volume / PHF  is D, D/mu = congestion period
            congestion_period_P = volume / PHF / mu;
            P = congestion_period_P * 4; //unit: 15 time slot

            t0 = max(0.0, mid_time_slot_no - P / 2.0);
            t3 = min((double)_MAX_TIMESLOT_PerPeriod-1, mid_time_slot_no + P / 2.0);

            // we need to recalculate gamma coefficient based on D/C ratio. based on assumption of beta = 4;
            // average waiting time based on eq. (32)
            // https://www.researchgate.net/publication/348488643_Introduction_to_connection_between_point_queue_model_and_BPR
            // w= gamma*power(D/mu,4)/(80*mu)
            // gamma = w*80*mu/power(D/mu,4)
            // avg_travel_time has been calculated based on standard BPR function, thus, the condition is that, PerformBPR() should be called before PerformBPR_X

            float uncongested_travel_time_in_hour = FFTT_in_hour;
            float w = avg_travel_time_in_hour - uncongested_travel_time_in_hour; //unit: hour
            // mu is read from the external link.csv file
            float D_over_mu = congestion_period_P;

            gamma = w * 80 * mu / pow(D_over_mu, 4);

            // derive t0 and t3 based on congestion duration p
            int t2 = m * (t3 - t0) + t0;
            //tt_relative is relative time
            for (int tt_relative = 0; tt_relative <= L; ++tt_relative)
            {
                //absolute time index
                int time_abs = starting_time_slot_no + tt_relative;
                if (time_abs < t0)
                {
                    //first uncongested phase with mu/2 as the approximate flow rates
                    waiting_time[time_abs] = 0;  // per hour
                    arrival_rate[time_abs] = mu / 2;
                    discharge_rate[time_abs] = mu / 2.0;
                    travel_time[time_abs] = FFTT_in_hour;  // per hour
                }

                if (time_abs >= t0 && time_abs <= t3 && t3 > t0)
                {
                    float t_ph = time_abs / (60.0/ MIN_PER_TIMESLOT);
                    float est_t0h = t0 / (60.0 / MIN_PER_TIMESLOT);
                    float t2_ph = t2 / (60.0 / MIN_PER_TIMESLOT);
                    float est_t3h = t3 / (60.0 / MIN_PER_TIMESLOT);

                    //second congested phase based on the gamma calculated from the dynamic eq. (32)
                    Queue[time_abs] = 1 / (4.0 ) * gamma * (t_ph - est_t0h) * (t_ph - est_t0h) * (t_ph - est_t3h) * (t_ph - est_t3h);
                    // unit is hour
                    waiting_time[time_abs] = 1 / (4.0*mu) *gamma *(t_ph - est_t0h)*(t_ph - est_t0h) * (t_ph - est_t3h)*(t_ph - est_t3h);
                    arrival_rate[time_abs] = gamma * (t_ph - est_t0h)*(t_ph - t2_ph)*(t_ph - est_t3h) + mu;
                    discharge_rate[time_abs] = mu;
                    travel_time[time_abs] = FFTT_in_hour + waiting_time[time_abs]; // per hour
                    bValidQueueData = true;
                }

                if (time_abs > t3)
                {
                    //third uncongested phase with mu/2 as the approximate flow rates
                    waiting_time[time_abs] = 0;
                    arrival_rate[time_abs] = mu / 2;
                    discharge_rate[time_abs] = mu / 2.0;
                    travel_time[time_abs] = FFTT_in_hour;
                }
                // avg_waiting_time = gamma / (120 * mu)*pow(P, 4.0) *60.0;// avg_waiting_time  should be per min
                // dtalog() << avg_waiting_time << endl;
                // avg_travel_time = (FFTT_in_hour + avg_waiting_time)*60.0; // avg_travel time should be per min
            }
        }

        return avg_travel_time;
    }

	double m;
    // we should also pass uncongested_travel_time as length/(speed_at_capacity)
	double VOC;
    //updated BPR-X parameters
	double gamma;
	double mu;
    //peak hour factor
	double PHF;
    //standard BPR parameter
	double alpha;
	double beta;
    double preload;


	double toll[_MAX_AGNETTYPES];
    double pce[_MAX_AGNETTYPES];

	double penalty; 
	string allowed_uses;


	double rho;
	double marginal_base;
    // in 15 min slot
    int starting_time_slot_no;
    int ending_time_slot_no;

    float cycle_length;
    float red_time;
    int t0, t3;

    bool bValidQueueData;
    string period;

	double capacity;
	double FFTT;

	double congestion_period_P;
    // inpput
	double volume;

    //output
    double avg_delay;
	double avg_travel_time = 0;
	double avg_waiting_time = 0;
	double total_travel_time = 0;

    // t starting from starting_time_slot_no if we map back to the 24 hour horizon
    float Queue[_MAX_TIMESLOT_PerPeriod];
    float waiting_time[_MAX_TIMESLOT_PerPeriod];
    float arrival_rate[_MAX_TIMESLOT_PerPeriod];

    float discharge_rate[_MAX_TIMESLOT_PerPeriod];
    float travel_time[_MAX_TIMESLOT_PerPeriod];
};


class CLink
{
public:
    // construction
	CLink() :main_node_id{ -1 }, obs_count{ -1 }, upper_bound_flag{ 0 }, est_count_dev{ 0 }, free_speed{0},
        BWTT_in_simulation_interval{ 100 }, zone_seq_no_for_outgoing_connector{ -1 }, number_of_lanes{ 1 }, lane_capacity{ 1999 },
        length{ 1 }, free_flow_travel_time_in_min{ 1 }, link_spatial_capacity{ 100 },
        service_arc_flag{ false }, Scenario_evaluation_flag{ false }, traffic_flow_code{ 0 }, spatial_capacity_in_vehicles{ 999999 }, link_type{ 2 }, FT{ 1 }, AT{ 1 }, s3_m{ 4 }, tmc_road_order {0}
    {
        for (int tau = 0; tau < _MAX_TIMEPERIODS; ++tau)
        {
            flow_volume_per_period[tau] = 0;
            queue_length_perslot[tau] = 0;
            travel_time_per_period[tau] = 0;
            TDBaseTT[tau] = 0;
            TDBaseCap[tau] = 0;
            TDBaseFlow[tau] = 0;
            TDBaseQueue[tau] = 0;
            //cost_perhour[tau] = 0;
            Scenario_STA_VOC_Ratio[tau] = 1;

            for(int at = 0; at < _MAX_AGNETTYPES; ++at)
                volume_per_period_per_at[tau][at] = 0;
        }
    }

    ~CLink()
    {
    }

    // Peiheng, 02/05/21, useless block
    void free_memory()
    {
    }

    void CalculateTD_VDFunction();

    float get_VOC_ratio(int tau)
    {
        return (flow_volume_per_period[tau] + TDBaseFlow[tau]) / max(0.00001, TDBaseCap[tau]);
    }

    float get_speed(int tau)
    {
        return length / max(travel_time_per_period[tau], 0.0001) * 60;  // per hour
    }

    void calculate_marginal_cost_for_agent_type(int tau, int agent_type_no, float PCE_agent_type)
    {
        // volume * dervative
        // BPR_term: volume * FFTT * alpha * (beta) * power(v/c, beta-1),

        travel_marginal_cost_per_period[tau][agent_type_no] = VDF_period[tau].marginal_base * PCE_agent_type;
    }

    float get_generalized_first_order_gradient_cost_of_second_order_loss_for_agent_type(int tau, int agent_type_no)
    {
        // *60 as 60 min per hour
        float generalized_cost = travel_time_per_period[tau] + VDF_period[tau].penalty + VDF_period[tau].toll[agent_type_no] / assignment.g_AgentTypeVector[agent_type_no].value_of_time * 60;

        // system optimal mode or exterior panalty mode
        if (assignment.assignment_mode == 4)
            generalized_cost += travel_marginal_cost_per_period[tau][agent_type_no];

        return generalized_cost;
    }

    int main_node_id;

    float obs_count;
    int upper_bound_flag;
    float est_count_dev;

    int BWTT_in_simulation_interval;
    int zone_seq_no_for_outgoing_connector;  // 

    int number_of_lanes;
    double lane_capacity;
	double length;
	double free_flow_travel_time_in_min;
	double free_speed;

	double cost;
	double link_spatial_capacity;

    bool service_arc_flag;
    int traffic_flow_code;
    int spatial_capacity_in_vehicles;

    // 1. based on BPR.

    int link_seq_no;
    string link_id;
    string name;
    string geometry;
    bool AllowAgentType(string agent_type, int tau)
    {
        if (VDF_period[tau].allowed_uses.size() == 0 || VDF_period[tau].allowed_uses == "all")  // if the allowed_uses is empty then all types are allowed.
            return true;
        else
        {
            if (VDF_period[tau].allowed_uses.find(agent_type) != string::npos)  // otherwise, only an agent type is listed in this "allowed_uses", then this agent type is allowed to travel on this link
                return true;
            else
            {
                return false;
            }


        }
    }

    int from_node_seq_no;
    int to_node_seq_no;
    __int64 from_node_cell_id;
    __int64 to_node_cell_id;

    int link_type;
    string link_type_name;
    string link_type_code;
    string movement_str;
    string TMC_code;
    string tmc_corridor_name;
    int tmc_corridor_id;
  
    int tmc_road_order;

    int tmc_road_sequence;
    string tmc_road, tmc_direction, tmc_intersection;
    float tmc_reference_speed;
    float tmc_mean_speed;
    float VDF_STA_speed[5];
    float VDF_STA_VOC[5];
    float VDF_STA_volume[5];

    float Scenario_STA_VOC_Ratio[5];
    bool Scenario_evaluation_flag;

    float tmc_volume;
    GDPoint TMC_from, TMC_to;
    float TMC_highest_speed;
    int FT;
    int AT;
    float PCE;
    float fftt;
    float vc; // critical speed;
    float kc; // critical density;
    float s3_m; // m factor in s3 model

    void UpdateKC(float free_speed_value)
    {  
        int updateKc_method = 0; // 0: HCM
//        float speed_ratio = free_speed / max(1, speed);
        kc = lane_capacity * pow(2, 2 / s3_m) / max(1, free_speed_value);

        if (kc > 50)
            kc = 50;

        if (updateKc_method==0 && free_speed_value >= 55 && lane_capacity >=1500) // free flow speed > 55, treat as freeway facility type
        {
            kc = 45;  // 45 vehicles per mile per lane based on HCM
            vc = lane_capacity / kc;
            s3_m = 2 * log(2) / log(free_speed_value / vc);
        }else  // non freeway facility, capacity is partially determined by g/c ratio
        {
            vc = free_speed_value * 0.7;
            s3_m = 2 * log(2) / log(free_speed_value / vc);
        }
        TMC_highest_speed = free_speed_value;


       

    }

    double get_volume_from_speed(float speed, float free_speed_value)
    {
        //test data free_speed = 55.0f; 
        //speed = 52;
        //kc = 23.14167648;

        if (speed < 0)
            return -1;

        double speed_ratio = free_speed_value / max(1,speed);
        if (speed_ratio <= 1.00001)
            speed_ratio = 1.00001;

     /*   float volume = 0;*/
        double ratio_difference = pow(speed_ratio, s3_m / 2) - 1;

        double ratio_difference_final = max(ratio_difference,0.00000001);

        double volume = speed * kc* pow(ratio_difference_final, 1/s3_m);

        return volume;

    }

    

    CVDF_Period VDF_period[_MAX_TIMEPERIODS];

	double TDBaseTT[_MAX_TIMEPERIODS];
	double TDBaseCap[_MAX_TIMEPERIODS];
	double TDBaseFlow[_MAX_TIMEPERIODS];
	double TDBaseQueue[_MAX_TIMEPERIODS];

    int type;

    //static
    //float flow_volume;
    //float travel_time;

    double flow_volume_per_period[_MAX_TIMEPERIODS];
	double  volume_per_period_per_at[_MAX_TIMEPERIODS][_MAX_AGNETTYPES];

	double  queue_length_perslot[_MAX_TIMEPERIODS];  // # of vehicles in the vertical point queue
	double travel_time_per_period[_MAX_TIMEPERIODS];
	double  travel_marginal_cost_per_period[_MAX_TIMEPERIODS][_MAX_AGNETTYPES];

    int number_of_periods;

    //std::vector <SLinkMOE> m_LinkMOEAry;
    //beginning of simulation data

    //toll related link
    //int m_TollSize;
    //Toll *pTollVector;  // not using SLT here to avoid issues with OpenMP

    // for discrete event simulation
    // link-in queue of each link
    std::list<int> EntranceQueue;
    // link-out queue of each link
    std::list<int> ExitQueue;
};

class GridLinkSet
{
public:
    __int64 cell_id;
    std::vector<int> m_LinkNoVector;
};

class CTMCLink
{
public:

    CTMCLink()
    {
        bWithSensorSpeedData = false;
        for (int t = 0; t < _MAX_TIMEINTERVAL_PerDay; ++t)
        {
            speed_sum[t] = -1;
            avg_speed[t] = -1;
            est_speed[t] = -1;
            speed_lowest[t] = 99999;
            count[t] = 0;
        }


        //for (int d = 0; d < _MAX_DAY_PerYear; d++)
        //{
        //    for (int t = 0; t < _MAX_TIMEINTERVAL_PerDay; t++)
        //    {
        //        speed_day[d][t] = -1;
        //    }
        //}

    }

    ~CTMCLink()
    {
    }


    float GetHighestSpeed()
    {
        float highest_speed = 0;

        for (int t_in_min = 6 * 60; t_in_min < 20 * 60; t_in_min += 5)
        {

            float avg_speed = RecordAvgSpeed(t_in_min);

            if (avg_speed > highest_speed)
                highest_speed = avg_speed;
        }

        return highest_speed;
    }


    float ScanCongestionDuration(int peak_no, int starting_time_in_hour, int ending_time_in_hour, float &FD_vc, float& obs_t0_in_hour, float& obs_t3_in_hour,
        CLink* pLink, float& assignment_volume, float& assignment_VMT, float& assignment_VHT, float& assignment_VDT, float& assignment_VCDT, float& congestion_demand,  float& congestion_dc_ratio, float& mean_congestion_flow_rate_mu, float& mean_congestion_speed,
        float& highest_speed, float& mean_speed, float& t2_speed, float& t2_queue, float& gamma)
    {

        obs_t0_in_hour = -1;
        obs_t3_in_hour = -1;
        est_QHF[peak_no] = 1;
        est_QDF_n[peak_no] = 1;
        est_QDF_s[peak_no] = 1;

        assignment_VHT = 0;
        assignment_VMT = 0;
        assignment_VDT = 0;

        int obs_t0_in_interval = starting_time_in_hour * 12;
        int obs_t3_in_interval = ending_time_in_hour * 12;

        float total_speed_value = 0;
        float total_speed_count = 0;

        highest_speed = 0;


        for (int t_in_min = 6 * 60; t_in_min < 20 * 60; t_in_min += 5)
        {

            float avg_speed = RecordAvgSpeed(t_in_min);

            if (avg_speed > highest_speed)
                highest_speed = avg_speed;
        }


        for (int t_in_min = starting_time_in_hour * 60; t_in_min < ending_time_in_hour * 60; t_in_min += 5)
        {

            float avg_speed = RecordAvgSpeed(t_in_min);

            total_speed_value += avg_speed;
            total_speed_count++;
        }


        // use the measured speed to update 
        pLink->UpdateKC(highest_speed);
        FD_vc = pLink->vc;

        mean_speed = total_speed_value / max(1, total_speed_count);

        int t_mid = (starting_time_in_hour + ending_time_in_hour) * 12 / 2;
        est_QHF[peak_no] = max(0.1, ending_time_in_hour - starting_time_in_hour); // number of hours in peak period
        assignment_volume = 0;
        
        int t_lowest_speed = -1;
        double lowest_speed = 99999;


        if(avg_speed[t_mid] > 0)
        {
            float VHT_sum = 0;
            float VMT_sum = 0;
            float VDT_sum = 0;
            float VCDT_sum = 0; // congestion_delay
            for (int t_in_min = starting_time_in_hour * 60; t_in_min < ending_time_in_hour * 60; t_in_min += 5)
            {

                float avg_speed = RecordAvgSpeed(t_in_min);
                float volume = pLink->get_volume_from_speed(avg_speed, highest_speed);
                float VHT = volume / 12 * pLink->length / avg_speed;
                float VMT = volume / 12 * pLink->length;
                float VDT = volume / 12 * max(0,pLink->length / avg_speed- pLink->length / (max(1,highest_speed))) ;  //base is free-flow travel time 
                float VCDT = volume / 12 * max(0, pLink->length / avg_speed - pLink->length / (max(1, pLink->vc)));  // base is vc; speed at capacity
                VHT_sum += VHT * pLink->number_of_lanes;
                VMT_sum += VMT * pLink->number_of_lanes;
                VDT_sum += VDT * pLink->number_of_lanes;
                VCDT_sum += VCDT * pLink->number_of_lanes;
                assignment_volume += volume / 12; // 12 5-min interval per hour


                if (avg_speed < lowest_speed)
                {
                    t_lowest_speed = t_in_min / 5;
                    lowest_speed = avg_speed;
                }
                 


            }

            assignment_VHT = VHT_sum;
            assignment_VMT = VMT_sum;
            assignment_VDT = VDT_sum;
            assignment_VCDT = VCDT_sum;
        }

        if (avg_speed[t_mid] < 0)
            return 0;

        if (avg_speed[t_mid] >= FD_vc )
        {
            // condition 1: the congestion in the middle 
            //not met
            // condition 2: 
            if (lowest_speed <= FD_vc)
            {
                t_mid = t_lowest_speed; // update t_mid using lowest speed time index.
            }
            else
            {
                return 0;
            }
        }


        // below is congested period analysis

        bool b_t0_found_flag = false;
        bool b_t3_found_flag = false;

        if (avg_speed[t_mid] < FD_vc)  // exit if the speed is higher than vc
        {
            obs_t3_in_interval = t_mid;
            obs_t3_in_hour = t_mid * 1.0 / 12.0;
            b_t3_found_flag = true;

            obs_t0_in_interval = t_mid;
            b_t0_found_flag = true;
            obs_t0_in_hour = t_mid * 1.0 / 12.0;

        }

        for (int t = t_mid + 1; t < 24 * 12; t += 1)  // move forward from mid t
        {
            if (avg_speed[t] < 1)
            {
                int i_no_data = 1;

//                dtalog.output() << " Error: TMC_link " << pLink->TMC_code.c_str() << " has no data at timt t " << t << g_time_coding(t * 5).c_str() << endl;;
                break;
            }

            obs_t3_in_interval = t;
            obs_t3_in_hour = t * 1.0 / 12.0;
            b_t3_found_flag = true;

            if (avg_speed[t] > FD_vc)  // exit if the speed is higher than vc
            {
                break;
            }
        }


        for (int t = t_mid - 1; t >= 0 * 12; t -= 1) // move backward from t_mid 
        {
            if (avg_speed[t] < 1)
            {
                int i_no_data = 1;
                dtalog.output() << " Error: TMC_link " << pLink->TMC_code.c_str() << " has no data at timt t =" << t*5 << " min" << endl;;
                break;
            }
            if (avg_speed[t] > FD_vc)
            {

                break;
            }

            obs_t0_in_interval = t;
            b_t0_found_flag = true;
            obs_t0_in_hour = t * 1.0 / 12.0;

        }


        if (b_t0_found_flag == false || b_t3_found_flag == false)  // two boundaries are found
            return 0;


       congestion_demand = 0;
       total_speed_count = 0;  // initial values
       total_speed_value = 0; 

        for (int t = obs_t0_in_interval; t <= obs_t3_in_interval; t += 1) // move between congestion duration per interval
        {
            if(avg_speed[t] >=1)
            {
            total_speed_value += avg_speed[t];
            float volume = pLink->get_volume_from_speed(avg_speed[t], highest_speed);
            congestion_demand += volume/12; // 12 5-min interval per hour
            total_speed_count ++;
            }

        }




        // test
        float test_volume = pLink->get_volume_from_speed(avg_speed[t_mid], highest_speed);

        mean_congestion_flow_rate_mu = congestion_demand / max(1, total_speed_count)*12;
        float    obs_P_in_hour = (obs_t3_in_hour - obs_t0_in_hour);  // congestion duration P

        congestion_demand = mean_congestion_flow_rate_mu * obs_P_in_hour;  // recalculate D in case there are missing data. 
        
        congestion_dc_ratio = congestion_demand / max(1, pLink->lane_capacity);  // unit: demand: # of vehicles, lane_capacity # of vehicles per hours: dc ratio has a unit of hour, but it is different from P

        mean_congestion_speed = total_speed_value / max(1, total_speed_count);

        t2_speed = avg_speed[t_mid];  // if we use a pure second order model, we should consider t2= 2/3(t3-t0)+ t0
        float t2_waiting_time = 1 / t2_speed - 1 / FD_vc;  // assume 1 mile road, unit: 1 hour
        t2_queue = t2_waiting_time* mean_congestion_flow_rate_mu;
        gamma = t2_queue * 12/ pow((obs_t3_in_hour - obs_t0_in_hour) /2.0, 4);

        
    // calibration
        if (obs_P_in_hour > 5) // setup upper bound
            obs_P_in_hour = 5;

        if(obs_P_in_hour > 1)  // only for congested hour
        {
        est_QHF[peak_no] = max(1,assignment_volume / max(pLink->lane_capacity, congestion_demand));
        est_QDF_n[peak_no] = max(1,log(obs_P_in_hour) / log(congestion_dc_ratio));
        float Vmin_ratio = FD_vc / max(1,t2_speed);
        est_QDF_s[peak_no] = max(0.0001, log(Vmin_ratio) / log(obs_P_in_hour));
        est_gamma[peak_no] = gamma;
        }


        return obs_P_in_hour;
    }

    string tmc_code;

    int matching_link_no;
    bool bWithSensorSpeedData;

    // construction
    void AddSpeedData(int day_no, int time_in_min, float speed_value)
    {
        bWithSensorSpeedData = true;
        int t = time_in_min / 5;

        if (t >= 0 && t < _MAX_TIMEINTERVAL_PerDay)
        {
            //            speed_day[day_no][t] = speed_value;

            if (speed_sum[t] < 0 || count[t] ==0)
                speed_sum[t] = 0; 

            speed_sum[t] += speed_value;
            count[t] += 1;

            if (speed_value < speed_lowest[t])
                speed_lowest[t] = speed_value;

        }

    }

    float RecordAvgSpeed(int time_in_min)
    {
        int t = time_in_min / 5;

        if (count[t] == 0)
            return -1;

        if (t >= 0 && t < _MAX_TIMEINTERVAL_PerDay)
        {

            avg_speed[t] =  speed_sum[t] / max(1, count[t]);
            return avg_speed[t];
        }

    }

    float GetAvgSpeed(int time_in_min)
    {
        int t = time_in_min / 5;
        if (t >= 0 && t < _MAX_TIMEINTERVAL_PerDay)
        {

            return speed_sum[t] / max(1, count[t]);

        }

    }

    float GetAvgSpeed_15min(int time_in_min)
    {
        int t = time_in_min / 5;

        float total_speed_value = 0;
        int total_speed_count = 0;

        for (int tt = 0; tt < 3; tt++)
        {

            if (t + tt >= 0 && t + tt < _MAX_TIMEINTERVAL_PerDay)
            {
                total_speed_value += speed_sum[t + tt] / max(1, count[t + tt]);
                total_speed_count++;
            }
        }

        return total_speed_value / max(1, total_speed_count);

    }

    float GetAvgHourlySpeed(int time_in_min)
    {
        int t = time_in_min / 5;

        float total_speed_value = 0;
        int total_speed_count = 0;

        for(int tt=0; tt<12;tt++)
        {
            
            if (t+tt >= 0 && t+tt < _MAX_TIMEINTERVAL_PerDay)
            {
                total_speed_value+= speed_sum[t+tt] / max(1, count[t+tt]);
                total_speed_count++;
            }
        }

        return total_speed_value / max(1, total_speed_count);
    }

    float GetDaySpeed(int day_no, int time_in_min)
    {
        //int t = time_in_min / 5;
        //if (t >= 0 && t < _MAX_TIMEINTERVAL_PerDay && day_no >= 0 && day_no < _MAX_DAY_PerYear)
        //{

        //    return speed_day[day_no][t];

        //}

        return -1;
    }

    float GetLowestSpeed(int time_in_min)
    {
        int t = time_in_min / 5;
        if (t >= 0 && t < _MAX_TIMEINTERVAL_PerDay)
        {

            return speed_lowest[t];

        }

    }


    float reference_speed;

    float avg_speed[_MAX_TIMEINTERVAL_PerDay];
    float speed_sum[_MAX_TIMEINTERVAL_PerDay];
    //float speed_day[_MAX_DAY_PerYear][_MAX_TIMEINTERVAL_PerDay];

    float est_speed[_MAX_TIMEINTERVAL_PerDay];
    float est_per_hour_per_lane[_MAX_TIMEINTERVAL_PerDay];

    float pred_speed[_MAX_TIMEINTERVAL_PerDay];
    float pred_per_hour_per_lane[_MAX_TIMEINTERVAL_PerDay];

    // AM 0: noon 1, PM: 2

    float est_QHF[_MAX_TIMEPERIODS];
    float est_QDF_n[_MAX_TIMEPERIODS];  // P= (D/C)^n
    float est_QDF_s[_MAX_TIMEPERIODS];  // vc/vt2= (P)^s, P Vmin slope

    // first layer of varaibles 
    float est_D[_MAX_TIMEPERIODS];
    float est_DCRatio[_MAX_TIMEPERIODS];

    // 2nd layer of varaibles 
    float est_P[_MAX_TIMEPERIODS];
    float est_mu[_MAX_TIMEPERIODS];
    float est_vcd[_MAX_TIMEPERIODS];
    float est_vt2[_MAX_TIMEPERIODS];

    float est_gamma[_MAX_TIMEPERIODS];
    float est_t0[_MAX_TIMEPERIODS];
    float est_t3[_MAX_TIMEPERIODS];

    float GetEstSpeed(int time_in_min)
    {
        int t = time_in_min / 5;
        if (t >= 0 && t < _MAX_TIMEINTERVAL_PerDay)
        {
            if(est_speed[t]>=1)
                return est_speed[t];
        }
        return 0;
    }
    float GetEstHourlySpeed(int time_in_min)
    {
        int t = time_in_min / 5;
        float total_speed_value = 0;
        int total_speed_count = 0;

        for (int tt = 0; tt < 12; tt++)
        {

            if (t + tt >= 0 && t + tt < _MAX_TIMEINTERVAL_PerDay)
            {
                if (est_speed[t+tt] >= 1)
                {
                    total_speed_value += est_speed[t + tt];
                    total_speed_count++;
                }
            }
        }

        return total_speed_value / max(1, total_speed_count);
    }
    float GetPredSpeed(int time_in_min)
    {
        int t = time_in_min / 5;
        if (t >= 0 && t < _MAX_TIMEINTERVAL_PerDay)
        {
            if (pred_speed[t] >= 1)
                return pred_speed[t];
        }
        return 0;
    }
    float GetPredHourlySpeed(int time_in_min)
    {
        int t = time_in_min / 5;
        float total_speed_value = 0;
        int total_speed_count = 0;

        for (int tt = 0; tt < 12; tt++)
        {

            if (t + tt >= 0 && t + tt < _MAX_TIMEINTERVAL_PerDay)
            {
                if (pred_speed[t + tt] >= 1)
                {
                    total_speed_value += pred_speed[t + tt];
                    total_speed_count++;
                }
            }
        }

        return total_speed_value / max(1, total_speed_count);
    }

    float PerformEstimation(bool bPredictionFlag, CLink* pLink, int peak_no, float assign_period_start_time_in_hour, float assign_period_end_time_in_hour, float t2, float V, float laneCapacity,  float vf, float vcd, float vct, float& DTASpeed, float& DTAP, float& TDAvgSpeedDiff)
    {
        int p = peak_no;

        TDAvgSpeedDiff = -1;

        est_D[p] = V / est_QHF[p];
        double coef_a = pow(pLink->kc, pLink->s3_m);
        double coef_b = pow(pLink->kc, pLink->s3_m) * pow(vf, pLink->s3_m / 2.0);
        double coef_c = pow(est_D[p], pLink->s3_m);

        DTASpeed = pow((coef_b + pow(coef_b * coef_b - 4.0 * coef_a * coef_c, 0.5)) / (2.0 * coef_a), 2.0 / pLink->s3_m);    //under uncongested condition

        est_DCRatio[p] = est_D[p] / laneCapacity;
        est_P[p] = pow(est_DCRatio[p], est_QDF_n[p]);
        est_mu[p] = pow(est_DCRatio[p], est_QDF_n[p] * (-1)) * est_D[p];
        est_vt2[p] = vcd / pow(est_P[p], est_QDF_s[p]); //est_vt2=vc/P^s   

        DTAP = est_P[p];
        // setup up default values
        for (int t_in_min = assign_period_start_time_in_hour * 60; t_in_min < assign_period_end_time_in_hour * 60; t_in_min += 5)
        {
            int t_interval = t_in_min / 5;

            if(bPredictionFlag)
                pred_speed[t_interval] = DTASpeed;
            else
                est_speed[t_interval] = DTASpeed;
        }

        if (DTAP < 1) // for P < 1, we will skip the creation of curve below vc. 
            return DTASpeed;

        float est_q_t2 = 1.0 / est_vt2[p] * est_mu[p];
        est_gamma[p] = est_q_t2 * 4 / pow(est_P[p] / 2, 4);  // because q_tw = w*mu =1/4 * gamma (P/2)^4, => 1/vt2 * mu = 1/4 * gamma  * (P/2)^4
        
        est_t0[p] = t2 - 0.5 * est_P[p];
        est_t3[p] = t2 + 0.5 * est_P[p];;
        
        //if (est_P[p] > 2 && est_DCRatio[p] > 2)
        //{
        //    int idebug = 1;
        //}


        for (int t_in_min = assign_period_start_time_in_hour * 60; t_in_min < assign_period_end_time_in_hour * 60; t_in_min += 5)
        {
            float t = t_in_min / 60.0;  // t in hour
            float td_queue;

            if (est_t0[p] <= t && t <= est_t3[p])
                td_queue = 0.25 * est_gamma[p] * pow((t - est_t0[p]), 2) * pow((t - est_t3[p]), 2);
            else
                td_queue = 0;

            float td_speed;
            if (est_P[p] < 0.25)
            {
                if (t < est_t0[p])
                {
                    td_speed = vf - ((vf - vct) / (est_t0[p] - assign_period_start_time_in_hour)) * (t - assign_period_start_time_in_hour);
                }
                else if (t > est_t3[p])
                {
                    td_speed = vf - ((vf - vct) / (assign_period_end_time_in_hour - est_t3[p])) * (assign_period_end_time_in_hour - t);
                }
                else
                    td_speed = vcd;
            }
            else // est_P > 0.25
            {
                if (t < est_t0[p])
                    td_speed = vf - ((vf - vcd) / (est_t0[p] - assign_period_start_time_in_hour)) * (t - assign_period_start_time_in_hour);
                else if (t > est_t3[p])
                    td_speed = vf - ((vf - vcd) / (assign_period_end_time_in_hour - est_t3[p]))* (assign_period_end_time_in_hour - t);
                else

                    td_speed = 1 / ( (td_queue / est_mu[p]) + (1 / vcd) );

//                cout << "td_queue t" << t << " =  " << td_queue << ", speed =" << td_speed << endl;
            }

            int t_interval = t_in_min / 5;

            if (bPredictionFlag)
                pred_speed[t_interval] = td_speed;
            else
                est_speed[t_interval] = td_speed;

        }

        if (bPredictionFlag == false)
        {
            // tally the mean DTA speed
            float est_speed_total = 0;
            int est_speed_count = 0;
            float total_speed_abs_diff = 0;
            for (int t_in_min = assign_period_start_time_in_hour * 60; t_in_min < assign_period_end_time_in_hour * 60; t_in_min += 5)
            {
                int t_interval = t_in_min / 5;
                est_speed_total += est_speed[t_interval];
                est_speed_count++;

                if (avg_speed[t_interval] >= 1)  // with data
                {
                    total_speed_abs_diff += fabs(est_speed[t_interval] - avg_speed[t_interval]);
                }
            }

            TDAvgSpeedDiff = total_speed_abs_diff / max(1, est_speed_count);
            DTASpeed = est_speed_total / max(1, est_speed_count);
        return est_vt2[p];
        }
        else
        {
            // tally the mean DTA speed
            float pred_speed_total = 0;
            int pred_speed_count = 0;
            float total_speed_abs_diff = 0;
            for (int t_in_min = assign_period_start_time_in_hour * 60; t_in_min < assign_period_end_time_in_hour * 60; t_in_min += 5)
            {
                int t_interval = t_in_min / 5;
                pred_speed_total += pred_speed[t_interval];
                pred_speed_count++;

                if (est_speed[t_interval] >= 1)  // with data
                {
                    total_speed_abs_diff += fabs(pred_speed[t_interval] - est_speed[t_interval]);
                }
            }

            TDAvgSpeedDiff = total_speed_abs_diff / max(1, pred_speed_count);
            DTASpeed = pred_speed_total / max(1, pred_speed_count);
            return est_vt2[p];

        }
    
    }

    float speed_lowest[_MAX_TIMEINTERVAL_PerDay];

    float volume_per_hour_per_lane[_MAX_TIMEINTERVAL_PerDay];
    float count[_MAX_TIMEINTERVAL_PerDay];

};
class CServiceArc
{
public:
    CServiceArc() : cycle_length{ 0 }, red_time{ 0 }, capacity{ 0 }, travel_time_delta{ 0 }
    {
    }

    int cycle_length;
    float red_time;
    float capacity;
    int travel_time_delta;

    int link_seq_no;
    int starting_time_no;
    int ending_time_no;
    int time_interval_no;
};

class CNode
{
public:
    CNode() : zone_id{ -1 }, zone_org_id{ -1 }, prohibited_movement_size{ 0 }, node_seq_no{ -1 }, production{ 0 }, attraction{ 0 }, is_activity_node{ 0 }, ctrl_type{ 0 }
    {
    }

    //int accessible_node_count;

    int zone_id;
    __int64 cell_id;
    // original zone id for non-centriod nodes
    int  zone_org_id;
    int prohibited_movement_size;
    // sequence number
    int node_seq_no;

    //external node number
    int node_id;
    int ctrl_type;

    double x;
    double y;

    int is_activity_node;
    double production;
    double attraction;

    std::vector<int> m_outgoing_link_seq_no_vector;
    std::vector<int> m_incoming_link_seq_no_vector;

    std::vector<int> m_to_node_seq_no_vector;
    std::map<int, int> m_to_node_2_link_seq_no_map;

    std::map<string, int> m_prohibited_movement_string_map;
};


std::vector<CNode> g_node_vector;
std::vector<CLink> g_link_vector;
std::vector<CTMCLink> g_TMC_vector;
std::vector<CServiceArc> g_service_arc_vector;

class CODMatrix
{
public:
    std::map <int, float> value_map;
    std::map <int, float> distance_map;
    std::map <int, double> disutility_map;

};

class COZone
{
public:
    COZone() : obs_production{ 0 }, obs_attraction{ 0 },
        est_production{ 0 }, est_attraction{ 0 },
        est_production_dev{ 0 }, est_attraction_dev{ 0 }
    {
    }

    _int64 cell_id;
    double cell_x;
    double cell_y;

    float obs_production;
    float obs_attraction;

    float gravity_production[_MAX_AGNETTYPES];
    float gravity_attraction[_MAX_AGNETTYPES];

    float gravity_est_production[_MAX_AGNETTYPES];
    float gravity_est_attraction[_MAX_AGNETTYPES];

    float est_production;
    float est_attraction;

    float est_production_dev;
    float est_attraction_dev;

    // 0, 1,
    int zone_seq_no;
    // external zone id // this is origin zone
    int zone_id;
    int node_seq_no;

    float obs_production_upper_bound_flag;
    float obs_attraction_upper_bound_flag;

    CODMatrix m_ODMatrix[_MAX_AGNETTYPES][_MAX_TIMEPERIODS];
    CODMatrix m_ODAccessibilityMatrix[_MAX_AGNETTYPES][_MAX_TIMEPERIODS];

    std::vector<int> m_activity_node_vector;


};

std::vector<COZone> g_zone_vector;
class CAGBMAgent
{
public:
    int agent_id;
    int income;
    int gender;
    int vehicle;
    int purpose;
    int flexibility;
    float preferred_arrival_time;
    float travel_time_in_min;
    float free_flow_travel_time;
    int from_zone_seq_no;
    int to_zone_seq_no;
    int type;
    int time_period;
    int k_path;
    float volume;
    float arrival_time_in_min;
};

std::vector<CAGBMAgent> g_agbmagent_vector;

struct CNodeForwardStar{
    CNodeForwardStar() : OutgoingLinkNoArray{ nullptr }, OutgoingNodeNoArray{ nullptr }, OutgoingLinkSize{ 0 }
    {
    }

    // Peiheng, 03/22/21, release memory
    ~CNodeForwardStar()
    {
        delete[] OutgoingLinkNoArray;
        delete[] OutgoingNodeNoArray;
    }

    int* OutgoingLinkNoArray;
    int* OutgoingNodeNoArray;
    int OutgoingLinkSize;
};

class NetworkForSP  // mainly for shortest path calculation
{
public:
    NetworkForSP() : temp_path_node_vector_size{ _MAX_LINK_SIZE_IN_A_PATH }, m_value_of_time{ 10 }, bBuildNetwork{ false }, m_memory_block_no{ 0 }
    {
    }

    ~NetworkForSP()
    {
        if (m_SENodeList)  //1
            delete[] m_SENodeList;

        if (m_node_status_array)  //2
            delete[] m_node_status_array;

        if (m_label_time_array)  //3
            delete[] m_label_time_array;

        if (m_label_distance_array)  //4
            delete[] m_label_distance_array;

        if (m_node_predecessor)  //5
            delete[] m_node_predecessor;

        if (m_link_predecessor)  //6
            delete[] m_link_predecessor;

        if (m_node_label_cost)  //7
            delete[] m_node_label_cost;

        if (m_link_flow_volume_array)  //8
            delete[] m_link_flow_volume_array;

        if (m_link_genalized_cost_array) //9
            delete[] m_link_genalized_cost_array;

        if (m_link_outgoing_connector_zone_seq_no_array) //10
            delete[] m_link_outgoing_connector_zone_seq_no_array;

        if (m_link_FFTT_array) //9
            delete[] m_link_FFTT_array;

        // Peiheng, 03/22/21, memory release on OutgoingLinkNoArray and OutgoingNodeNoArray
        // is taken care by the destructor of CNodeForwardStar
        if (NodeForwardStarArray)
            delete[] NodeForwardStarArray;
    }

    int temp_path_node_vector_size;
    float m_value_of_time;
    bool bBuildNetwork;
    int m_memory_block_no;

    //node seq vector for each ODK
    int temp_path_node_vector[_MAX_LINK_SIZE_IN_A_PATH];
    //node seq vector for each ODK
    int temp_path_link_vector[_MAX_LINK_SIZE_IN_A_PATH];

    bool m_bSingleSP_Flag;

    // assigned nodes for computing
    std::vector<int> m_origin_node_vector;
    std::vector<int> m_origin_zone_seq_no_vector;

    int m_tau; // assigned nodes for computing
    int m_agent_type_no; // assigned nodes for computing

    CNodeForwardStar* NodeForwardStarArray;

    int m_threadNo;  // internal thread number

    int m_ListFront; // used in coding SEL
    int m_ListTail;  // used in coding SEL
    int* m_SENodeList; // used in coding SEL

    // label cost for shortest path calcuating
    double* m_node_label_cost;
    // time-based cost
	double* m_label_time_array;
    // distance-based cost
	double* m_label_distance_array;

    // predecessor for nodes
    int* m_node_predecessor;
    // update status
    int* m_node_status_array;
    // predecessor for this node points to the previous link that updates its label cost (as part of optimality condition) (for easy referencing)
    int* m_link_predecessor;

	double* m_link_flow_volume_array;

	double* m_link_genalized_cost_array;
    double* m_link_FFTT_array;

    int* m_link_outgoing_connector_zone_seq_no_array;

    // major function 1:  allocate memory and initialize the data
    void AllocateMemory(int number_of_nodes, int number_of_links)
    {
        NodeForwardStarArray = new CNodeForwardStar[number_of_nodes];

        m_SENodeList = new int[number_of_nodes];  //1

        m_LinkBasedSEList = new int[number_of_links];  //1;  // dimension: number of links

        m_node_status_array = new int[number_of_nodes];  //2
        m_label_time_array = new double[number_of_nodes];  //3
        m_label_distance_array = new double[number_of_nodes];  //4
        m_node_predecessor = new int[number_of_nodes];  //5
        m_link_predecessor = new int[number_of_nodes];  //6
        m_node_label_cost = new double[number_of_nodes];  //7

        m_link_flow_volume_array = new double[number_of_links];  //8

        m_link_genalized_cost_array = new double[number_of_links];  //9
        m_link_outgoing_connector_zone_seq_no_array = new int[number_of_links]; //10
    }

    void UpdateGeneralizedLinkCost(int agent_type_no)
    {
        for (int i = 0; i < g_link_vector.size(); ++i)
        {
            CLink* pLink = &(g_link_vector[i]);

            if(pLink!=NULL)
            {
            m_link_genalized_cost_array[i] = pLink->travel_time_per_period[m_tau] + pLink->VDF_period[m_tau].penalty + pLink->VDF_period[m_tau].toll[agent_type_no] / m_value_of_time * 60;  // *60 as 60 min per hour
            m_link_FFTT_array[i] = pLink->free_flow_travel_time_in_min;
            }
          //route_choice_cost 's unit is min
        }
    }

    void BuildNetwork(Assignment* p_assignment)
    {
        if (bBuildNetwork)
            return;

        int m_outgoing_link_seq_no_vector[_MAX_LINK_SIZE_FOR_A_NODE];
        int m_to_node_seq_no_vector[_MAX_LINK_SIZE_FOR_A_NODE];

        for (int i = 0; i < g_link_vector.size(); ++i)
        {
            CLink* pLink = &(g_link_vector[i]);
            m_link_outgoing_connector_zone_seq_no_array[i] = pLink->zone_seq_no_for_outgoing_connector;
        }

        for (int i = 0; i < p_assignment->g_number_of_nodes; ++i) //Initialization for all non-origin nodes
        {
            int outgoing_link_size = 0;

            for (int j = 0; j < g_node_vector[i].m_outgoing_link_seq_no_vector.size(); ++j)
            {

                int link_seq_no = g_node_vector[i].m_outgoing_link_seq_no_vector[j];
                // only predefined allowed agent type can be considered
                if (g_link_vector[link_seq_no].AllowAgentType(p_assignment->g_AgentTypeVector[m_agent_type_no].agent_type, m_tau))
                    if(g_link_vector[link_seq_no].AllowAgentType (p_assignment->g_AgentTypeVector[m_agent_type_no].agent_type,m_tau))

                {
                    m_outgoing_link_seq_no_vector[outgoing_link_size] = link_seq_no;
                    m_to_node_seq_no_vector[outgoing_link_size] = g_node_vector[i].m_to_node_seq_no_vector[j];

                    outgoing_link_size++;

                    if (outgoing_link_size >= _MAX_LINK_SIZE_FOR_A_NODE)
                    {
                        dtalog.output() << " Error: node id: " << g_node_vector[i].node_id << " zone id:" 
                            << g_node_vector[i].zone_org_id << " has outgoing_link_size >= _MAX_LINK_SIZE_FOR_A_NODE" << endl;
                        cout << " Error: outgoing_link_size >= _MAX_LINK_SIZE_FOR_A_NODE" << endl;
                        outgoing_link_size = _MAX_LINK_SIZE_FOR_A_NODE - 1;
                        break; // continue with the last within bound size 

//                        g_ProgramStop();
                    }
                }
            }

            int node_seq_no = g_node_vector[i].node_seq_no;
            NodeForwardStarArray[node_seq_no].OutgoingLinkSize = outgoing_link_size;

            if(outgoing_link_size >= 1)
            {
                NodeForwardStarArray[node_seq_no].OutgoingLinkNoArray = new int[outgoing_link_size];
                NodeForwardStarArray[node_seq_no].OutgoingNodeNoArray = new int[outgoing_link_size];
            }

            for (int j = 0; j < outgoing_link_size; ++j)
            {
                NodeForwardStarArray[node_seq_no].OutgoingLinkNoArray[j] = m_outgoing_link_seq_no_vector[j];
                NodeForwardStarArray[node_seq_no].OutgoingNodeNoArray[j] = m_to_node_seq_no_vector[j];
            }
        }

        // after dynamic arrays are created for forward star
        if (dtalog.debug_level() == 2)
        {
            dtalog.output() << "add outgoing link data into dynamic array" << endl;
            cout << "add outgoing link data into dynamic array" << endl;

            for (int i = 0; i < g_node_vector.size(); ++i)
            {
                if (g_node_vector[i].zone_org_id > 0) // for each physical node
                { // we need to make sure we only create two way connectors between nodes and zones
                    dtalog.output() << "node id= " << g_node_vector[i].node_id << " with zone id " << g_node_vector[i].zone_org_id << "and "
                                    << NodeForwardStarArray[i].OutgoingLinkSize << " outgoing links." << endl;
                    cout << "node id= " << g_node_vector[i].node_id << " with zone id " << g_node_vector[i].zone_org_id << "and "
                        << NodeForwardStarArray[i].OutgoingLinkSize << " outgoing links." << endl;

                    for (int j = 0; j < NodeForwardStarArray[i].OutgoingLinkSize; j++)
                    {
                        int link_seq_no = NodeForwardStarArray[i].OutgoingLinkNoArray[j];
                        dtalog.output() << "  outgoing node = " << g_node_vector[g_link_vector[link_seq_no].to_node_seq_no].node_id << endl;
                        cout << "  outgoing node = " << g_node_vector[g_link_vector[link_seq_no].to_node_seq_no].node_id << endl;
                    }
                }
                else
                {
                    if (dtalog.debug_level() == 3)
                    {
                        dtalog.output() << "node id= " << g_node_vector[i].node_id << " with "
                                        << NodeForwardStarArray[i].OutgoingLinkSize << " outgoing links." << endl;

                        cout << "node id= " << g_node_vector[i].node_id << " with "
                            << NodeForwardStarArray[i].OutgoingLinkSize << " outgoing links." << endl;
                        for (int j = 0; j < NodeForwardStarArray[i].OutgoingLinkSize; ++j)
                        {
                            int link_seq_no = NodeForwardStarArray[i].OutgoingLinkNoArray[j];
                            dtalog.output() << "  outgoing node = " << g_node_vector[g_link_vector[link_seq_no].to_node_seq_no].node_id << endl;
                            cout << "  outgoing node = " << g_node_vector[g_link_vector[link_seq_no].to_node_seq_no].node_id << endl;
                        }
                    }
                }
            }
        }

        m_value_of_time = p_assignment->g_AgentTypeVector[m_agent_type_no].value_of_time;
        bBuildNetwork = true;
    }

    // SEList: scan eligible List implementation: the reason for not using STL-like template is to avoid overhead associated pointer allocation/deallocation
    inline void SEList_clear()
    {
        m_ListFront = -1;
        m_ListTail = -1;
    }

    inline void SEList_push_front(int node)
    {
        if (m_ListFront == -1)  // start from empty
        {
            m_SENodeList[node] = -1;
            m_ListFront = node;
            m_ListTail = node;
        }
        else
        {
            m_SENodeList[node] = m_ListFront;
            m_ListFront = node;
        }
    }

    inline void SEList_push_back(int node)
    {
        if (m_ListFront == -1)  // start from empty
        {
            m_ListFront = node;
            m_ListTail = node;
            m_SENodeList[node] = -1;
        }
        else
        {
            m_SENodeList[m_ListTail] = node;
            m_SENodeList[node] = -1;
            m_ListTail = node;
        }
    }

    inline bool SEList_empty()
    {
        return(m_ListFront == -1);
    }

    //	inline int SEList_front()
    //	{
    //		return m_ListFront;
    //	}

    //	inline void SEList_pop_front()
    //	{
    //		int tempFront = m_ListFront;
    //		m_ListFront = m_SENodeList[m_ListFront];
    //		m_SENodeList[tempFront] = -1;
    //	}

    //major function: update the cost for each node at each SP tree, using a stack from the origin structure

    double backtrace_shortest_path_tree(Assignment& assignment, int iteration_number, int o_node_index);

    //major function 2: // time-dependent label correcting algorithm with double queue implementation
    float optimal_label_correcting(int processor_id, Assignment* p_assignment, int iteration_k, int o_node_index, int d_node_no = -1, bool pure_travel_time_cost = false)
    {
        // g_debugging_flag = 1;
        int SE_loop_count = 0;

        if (iteration_k == 0)
            BuildNetwork(p_assignment);  // based on agent type and link type

		int agent_type = m_agent_type_no; // assigned nodes for computing
		UpdateGeneralizedLinkCost(agent_type);

        int origin_node = m_origin_node_vector[o_node_index]; // assigned nodes for computing
        int origin_zone_no = m_origin_zone_seq_no_vector[o_node_index]; // assigned nodes for computing


        if (p_assignment->g_number_of_nodes >= 1000 && origin_zone_no %97 == 0)
        {
            dtalog.output() << "label correcting for zone " << origin_zone_no <<  " in processor " << processor_id <<  endl;
        }

        if (dtalog.debug_level() >= 2)
        {
            dtalog.output() << "SP iteration k =  " << iteration_k << ": origin node: " << g_node_vector[origin_node].node_id << endl;
        }

        int number_of_nodes = p_assignment->g_number_of_nodes;
        //Initialization for all non-origin nodes
        for (int i = 0; i < number_of_nodes; ++i)
        {
            // not scanned
            m_node_status_array[i] = 0;
            m_node_label_cost[i] = _MAX_LABEL_COST;
            // pointer to previous NODE INDEX from the current label at current node and time
            m_link_predecessor[i] = -1;
            // pointer to previous NODE INDEX from the current label at current node and time
            m_node_predecessor[i] = -1;
            // comment out to speed up comuting
            ////m_label_time_array[i] = 0;
            ////m_label_distance_array[i] = 0;
        }

        // int internal_debug_flag = 0;
        if (NodeForwardStarArray[origin_node].OutgoingLinkSize == 0)
            return 0;

        //Initialization for origin node at the preferred departure time, at departure time, cost = 0, otherwise, the delay at origin node
        m_label_time_array[origin_node] = 0;
        m_node_label_cost[origin_node] = 0.0;
        //Mark:	m_label_distance_array[origin_node] = 0.0;

        // Peiheng, 02/05/21, duplicate initialization, remove them later
        // pointer to previous NODE INDEX from the current label at current node and time
        m_link_predecessor[origin_node] = -1;
        // pointer to previous NODE INDEX from the current label at current node and time
        m_node_predecessor[origin_node] = -1;

        SEList_clear();
        SEList_push_back(origin_node);

        int from_node, to_node;
        int link_sqe_no;
		double new_time = 0;
		double new_distance = 0;
		double new_to_node_cost = 0;
        int tempFront;
        int log_path = 0;

        while (!(m_ListFront == -1))   //SEList_empty()
        {
            // from_node = SEList_front();
            // SEList_pop_front();  // remove current node FromID from the SE list

            from_node = m_ListFront;//pop a node FromID for scanning
            tempFront = m_ListFront;
            m_ListFront = m_SENodeList[m_ListFront];
            m_SENodeList[tempFront] = -1;

            m_node_status_array[from_node] = 2;

            if (log_path >= 2)
            {
                dtalog.output() << "SP:scan SE node: " << g_node_vector[from_node].node_id << " with "
                                << NodeForwardStarArray[from_node].OutgoingLinkSize  << " outgoing link(s). "<< endl;
            }
            //scan all outbound nodes of the current node

            int pred_link_seq_no = m_link_predecessor[from_node];

            // for each link (i,j) belong A(i)
            for (int i = 0; i < NodeForwardStarArray[from_node].OutgoingLinkSize; ++i)
            {
                to_node = NodeForwardStarArray[from_node].OutgoingNodeNoArray[i];
                link_sqe_no = NodeForwardStarArray[from_node].OutgoingLinkNoArray[i];

                if (log_path >= 2)
                    dtalog.output() << "SP:  checking outgoing node " << g_node_vector[to_node].node_id << endl;

                // if(map (pred_link_seq_no, link_sqe_no) is prohibitted )
                //     then continue; //skip this is not an exact solution algorithm for movement

                if (g_node_vector[from_node].prohibited_movement_size >= 1)
                {
                    if (pred_link_seq_no >= 0)
                    {
                        string	movement_string;
                        string ib_link_id = g_link_vector[pred_link_seq_no].link_id;
                        string ob_link_id = g_link_vector[link_sqe_no].link_id;
                        movement_string = ib_link_id + "->" + ob_link_id;

                        if (g_node_vector[from_node].m_prohibited_movement_string_map.find(movement_string) != g_node_vector[from_node].m_prohibited_movement_string_map.end())
                        {
                            dtalog.output() << "prohibited movement " << movement_string << " will not be used " << endl;
                            cout << "prohibited movement " << movement_string << " will not be used " << endl;
                            continue;
                        }
                    }
                }

                //remark: the more complicated implementation can be found in paper Shortest Path Algorithms In Transportation Models: Classical and Innovative Aspects
                //	A note on least time path computation considering delays and prohibitions for intersection movements

                if (m_link_outgoing_connector_zone_seq_no_array[link_sqe_no] >= 0)
                {
                    if(m_link_outgoing_connector_zone_seq_no_array[link_sqe_no] != origin_zone_no)
                    {
                        // filter out for an outgoing connector with a centriod zone id different from the origin zone seq no
                        continue;
                    }
                }

                //very important: only origin zone can access the outbound connectors,
                //the other zones do not have access to the outbound connectors

                // Mark				new_time = m_label_time_array[from_node] + pLink->travel_time_per_period[tau];
                // Mark				new_distance = m_label_distance_array[from_node] + pLink->length;

                if (m_link_genalized_cost_array[link_sqe_no] > 60 && m_link_genalized_cost_array[link_sqe_no]/max(0.01,m_link_FFTT_array[link_sqe_no]) > 10)
                {   // avoid using over congested links
                    continue;
                }
                new_to_node_cost = m_node_label_cost[from_node] + m_link_genalized_cost_array[link_sqe_no];

                if (log_path >=2)
                {
                    dtalog.output() << "SP:  checking from node " << g_node_vector[from_node].node_id
                                    << "  to node" << g_node_vector[to_node].node_id << " cost = " << new_to_node_cost << endl;
                }

                if (new_to_node_cost < m_node_label_cost[to_node]) // we only compare cost at the downstream node ToID at the new arrival time t
                {
                    if (log_path >=2)
                    {
                        dtalog.output() << "SP:  updating node: " << g_node_vector[to_node].node_id << " current cost:" << m_node_label_cost[to_node]
                                        << " new cost " << new_to_node_cost << endl;
                    }

                    // update cost label and node/time predecessor
                    // m_label_time_array[to_node] = new_time;
                    // m_label_distance_array[to_node] = new_distance;
                    m_node_label_cost[to_node] = new_to_node_cost;
                    // pointer to previous physical NODE INDEX from the current label at current node and time
                    m_node_predecessor[to_node] = from_node;
                    // pointer to previous physical NODE INDEX from the current label at current node and time
                    m_link_predecessor[to_node] = link_sqe_no;

                    if (log_path >=2)
                    {
                        dtalog.output() << "SP: add node " << g_node_vector[to_node].node_id << " new cost:" << new_to_node_cost
                                        << " into SE List " << g_node_vector[to_node].node_id << endl;
                    }

                    // deque updating rule for m_node_status_array
                    if (m_node_status_array[to_node] == 0)
                    {
                        ///// SEList_push_back(to_node);
                        ///// begin of inline block
                        if (m_ListFront == -1)  // start from empty
                        {
                            m_ListFront = to_node;
                            m_ListTail = to_node;
                            m_SENodeList[to_node] = -1;
                        }
                        else
                        {
                            m_SENodeList[m_ListTail] = to_node;
                            m_SENodeList[to_node] = -1;
                            m_ListTail = to_node;
                        }
                        ///// end of inline block

                        m_node_status_array[to_node] = 1;
                    }

                    if (m_node_status_array[to_node] == 2)
                    {
                        /////SEList_push_front(to_node);
                        ///// begin of inline block
                        if (m_ListFront == -1)  // start from empty
                        {
                            m_SENodeList[to_node] = -1;
                            m_ListFront = to_node;
                            m_ListTail = to_node;
                        }
                        else
                        {
                            m_SENodeList[to_node] = m_ListFront;
                            m_ListFront = to_node;
                        }
                        ///// end of inline block

                        m_node_status_array[to_node] = 1;
                    }
                }
            }
        }

        if (log_path>=1)
        {
            dtalog.output() << "SPtree at iteration k = " << iteration_k <<  " origin node: "
                            << g_node_vector[origin_node].node_id  << endl;

            dtalog.output() << "origin node, origin zone no,origin zone id,SP node,label cost,node_pred_id,x_coord,y_coord" << endl;

            //Initialization for all non-origin nodes
            for (int i = 0; i < p_assignment->g_number_of_nodes; ++i)
            {
                int node_pred_id = -1;
                int node_pred_no = m_node_predecessor[i];

                if (node_pred_no >= 0)
                    node_pred_id = g_node_vector[node_pred_no].node_id;

                if(m_node_label_cost[i] < 9999)
                {
                    dtalog.output()  << origin_node << "," << origin_zone_no << "," << g_node_vector[i].zone_org_id << "," << g_node_vector[i].node_id << "," << m_node_label_cost[i] << ","
                                     << node_pred_id << "," << std::setprecision(9) << g_node_vector[i].x << "," << std::setprecision(9) << g_node_vector[i].y << endl;
                }
            }
        }

        if (d_node_no >= 1)
            return m_node_label_cost[d_node_no];
        else
            return 0;  // one to all shortest pat
    }

    ////////// link based SE List

    int m_LinkBasedSEListFront;
    int m_LinkBasedSEListTail;
    int* m_LinkBasedSEList;  // dimension: number of links

    // SEList: Scan List implementation: the reason for not using STL-like template is to avoid overhead associated pointer allocation/deallocation
    void LinkBasedSEList_clear()
    {
        m_LinkBasedSEListFront = -1;
        m_LinkBasedSEListTail = -1;
    }

    void LinkBasedSEList_push_front(int link)
    {
        if (m_LinkBasedSEListFront == -1)  // start from empty
        {
            m_LinkBasedSEList[link] = -1;
            m_LinkBasedSEListFront = link;
            m_LinkBasedSEListTail = link;
        }
        else
        {
            m_LinkBasedSEList[link] = m_LinkBasedSEListFront;
            m_LinkBasedSEListFront = link;
        }
    }

    void LinkBasedSEList_push_back(int link)
    {
        if (m_LinkBasedSEListFront == -1)  // start from empty
        {
            m_LinkBasedSEListFront = link;
            m_LinkBasedSEListTail = link;
            m_LinkBasedSEList[link] = -1;
        }
        else
        {
            m_LinkBasedSEList[m_LinkBasedSEListTail] = link;
            m_LinkBasedSEList[link] = -1;
            m_LinkBasedSEListTail = link;
        }
    }

    bool LinkBasedSEList_empty()
    {
        return(m_LinkBasedSEListFront == -1);
    }

    int LinkBasedSEList_front()
    {
        return m_LinkBasedSEListFront;
    }

    void LinkBasedSEList_pop_front()
    {
        int tempFront = m_LinkBasedSEListFront;
        m_LinkBasedSEListFront = m_LinkBasedSEList[m_LinkBasedSEListFront];
        m_LinkBasedSEList[tempFront] = -1;
    }
};


std::vector<NetworkForSP*> g_NetworkForSP_vector;
NetworkForSP g_RoutingNetwork;


void g_ReadDemandFileBasedOnDemandFileList(Assignment& assignment)
{
    //	fprintf(g_pFileOutputLog, "number of zones =,%lu\n", g_zone_vector.size());

    assignment.InitializeDemandMatrix(g_zone_vector.size(), assignment.g_AgentTypeVector.size(), assignment.g_DemandPeriodVector.size());

    float total_demand_in_demand_file = 0;

    CCSVParser parser;
    dtalog.output() << endl;
    dtalog.output() << "Step 1.8: Reading file section [demand_file_list] in setting.csv..." << endl;
    
    cout << "Step 1.8: Reading file section [demand_file_list] in setting.csv..." << endl;

    parser.IsFirstLineHeader = false;

    if (parser.OpenCSVFile("settings.csv", false))
    {
        while (parser.ReadRecord_Section())
        {
            if(parser.SectionName == "[demand_file_list]")
            {
                int file_sequence_no = 1;

                string format_type = "null";

                int demand_format_flag = 0;

                if (!parser.GetValueByFieldName("file_sequence_no", file_sequence_no))
                    break;

                // skip negative sequence no
                if (file_sequence_no <= -1)
                    continue;

                string file_name, demand_period, agent_type;
                parser.GetValueByFieldName("file_name", file_name);
                parser.GetValueByFieldName("demand_period", demand_period);
                parser.GetValueByFieldName("format_type", format_type);

                if (format_type.find("null") != string::npos)  // skip negative sequence no
                {
                    dtalog.output() << "Please provide format_type in section [demand_file_list.]" << endl;
                    cout << "Please provide format_type in section [demand_file_list.]" << endl;
                    g_ProgramStop();
                }

                parser.GetValueByFieldName("agent_type", agent_type);

                int agent_type_no = 0;
                int demand_period_no = 0;

                if (assignment.demand_period_to_seqno_mapping.find(demand_period) != assignment.demand_period_to_seqno_mapping.end())
                    demand_period_no = assignment.demand_period_to_seqno_mapping[demand_period];
                else
                {
                    dtalog.output() << "Error: demand period in section [demand_file_list]" << demand_period << "cannot be found." << endl;
                    cout << "Error: demand period in section [demand_file_list]" << demand_period << "cannot be found." << endl;
                    g_ProgramStop();
                }

                bool b_multi_agent_list = false;

                if (agent_type == "multi_agent_list")
                    b_multi_agent_list = true;
                else
                {
                    if (assignment.agent_type_2_seqno_mapping.find(agent_type) != assignment.agent_type_2_seqno_mapping.end())
                        agent_type_no = assignment.agent_type_2_seqno_mapping[agent_type];
                    else
                    {
                        dtalog.output() << "Error: agent_type in agent_type " << agent_type << "cannot be found." << endl;
                        cout << "Error: agent_type in agent_type " << agent_type << "cannot be found." << endl;
                        g_ProgramStop();
                    }
                }

                if (demand_period_no > _MAX_TIMEPERIODS)
                {
                    dtalog.output() << "demand_period_no should be less than settings in demand_period section. Please change the parameter settings in the source code." << endl;
                    cout << "demand_period_no should be less than settings in demand_period section. Please change the parameter settings in the source code." << endl;
                    g_ProgramStop();
                }

                if (format_type.find("column") != string::npos)  // or muliti-column
                {
                    bool bFileReady = false;
                    int error_count = 0;

                    // read the file formaly after the test.
                    FILE* st;
                    fopen_ss(&st, file_name.c_str(), "r");
                    if (st)
                    {
                        bFileReady = true;
                        int line_no = 0;

                        while (true)
                        {
                            int origin_zone = (int) (g_read_float(st));
                            int destination_zone =(int) g_read_float(st);
                            float demand_value = g_read_float(st);

                            if (origin_zone <= -1)
                            {
                                if (line_no == 1 && !feof(st))  // read only one line, but has not reached the end of the line
                                {
                                    dtalog.output() << endl << "Error: Only one line has been read from file. Are there multiple columns of demand type in file " << file_name << " per line?" << endl;
                                    cout << endl << "Error: Only one line has been read from file. Are there multiple columns of demand type in file " << file_name << " per line?" << endl;
                                    g_ProgramStop();
                                }
                                break;
                            }

                            if (assignment.g_zoneid_to_zone_seq_no_mapping.find(origin_zone) == assignment.g_zoneid_to_zone_seq_no_mapping.end())
                            {
                                if(error_count < 10)
                                {
                                    dtalog.output() << endl << "Warning: origin zone " << origin_zone << "  has not been defined in node.csv" << endl;
                                    cout << endl << "Warning: origin zone " << origin_zone << "  has not been defined in node.csv" << endl;
                                }

                                error_count++;
                                 // origin zone has not been defined, skipped.
                                continue;
                            }

                            if (assignment.g_zoneid_to_zone_seq_no_mapping.find(destination_zone) == assignment.g_zoneid_to_zone_seq_no_mapping.end())
                            {
                                if (error_count < 10)
                                {
                                    dtalog.output() << endl << "Warning: destination zone " << destination_zone << "  has not been defined in node.csv" << endl;
                                    cout << endl << "Warning: destination zone " << destination_zone << "  has not been defined in node.csv" << endl;
                                }

                                error_count++;
                                // destination zone has not been defined, skipped.
                                continue;
                            }

                            int from_zone_seq_no = 0;
                            int to_zone_seq_no = 0;
                            from_zone_seq_no = assignment.g_zoneid_to_zone_seq_no_mapping[origin_zone];
                            to_zone_seq_no = assignment.g_zoneid_to_zone_seq_no_mapping[destination_zone];

                            // encounter return
                            if (demand_value < -99)
                                break;

                            assignment.total_demand[agent_type_no][demand_period_no] += demand_value;
                            assignment.g_column_pool[from_zone_seq_no][to_zone_seq_no][agent_type_no][demand_period_no].od_volume += demand_value;
                            assignment.total_demand_volume += demand_value;
                            assignment.g_origin_demand_array[from_zone_seq_no][agent_type_no][demand_period_no] += demand_value;

                            // we generate vehicles here for each OD data line
                            if (line_no <= 5)  // read only one line, but has not reached the end of the line
                                dtalog.output() << "o_zone_id:" << origin_zone << ", d_zone_id: " << destination_zone << ", value = " << demand_value << endl;

                            line_no++;
                        }  // scan lines

                        fclose(st);

                        dtalog.output() << "total_demand_volume is " << assignment.total_demand_volume << endl << endl;
                    }
                    else
                    {
                        // open file
                        dtalog.output() << "Error: File " << file_name << " cannot be opened.\n It might be currently used and locked by EXCEL." << endl;
                        g_ProgramStop();
                    }
                }else if (format_type.find("geo") != string::npos)  // or muliti-column
                {
                    bool bFileReady = false;
                    int error_count = 0;

                    // read the file formaly after the test.
                    CCSVParser parser;
                    if (parser.OpenCSVFile(file_name, false))
                    {
                        bFileReady = true;
                        int line_no = 0;

                        while (parser.ReadRecord())
                        {
                            int origin_zone = -1;
                            int destination_zone = -1;
                            float demand_value = 0;
                            parser.GetValueByFieldName("o_zone_id", origin_zone);
                            parser.GetValueByFieldName("d_zone_id", destination_zone);
                            parser.GetValueByFieldName("volume", demand_value);


                            if (origin_zone <= -1)
                            {
                                if (line_no == 1)  // read only one line, but has not reached the end of the line
                                {
                                    dtalog.output() << endl << "Error: Only one line has been read from file. Are there multiple columns of demand type in file " << file_name << " per line?" << endl;
                                    cout << endl << "Error: Only one line has been read from file. Are there multiple columns of demand type in file " << file_name << " per line?" << endl;
                                    g_ProgramStop();
                                }
                                break;
                            }

                            if (assignment.g_zoneid_to_zone_seq_no_mapping.find(origin_zone) == assignment.g_zoneid_to_zone_seq_no_mapping.end())
                            {
                                if (error_count < 10)
                                {
                                    dtalog.output() << endl << "Warning: origin zone " << origin_zone << "  has not been defined in node.csv" << endl;
                                    cout << endl << "Warning: origin zone " << origin_zone << "  has not been defined in node.csv" << endl;
                                }

                                error_count++;
                                // origin zone has not been defined, skipped.
                                continue;
                            }

                            if (assignment.g_zoneid_to_zone_seq_no_mapping.find(destination_zone) == assignment.g_zoneid_to_zone_seq_no_mapping.end())
                            {
                                if (error_count < 10)
                                {
                                    dtalog.output() << endl << "Warning: destination zone " << destination_zone << "  has not been defined in node.csv" << endl;
                                    cout << endl << "Warning: destination zone " << destination_zone << "  has not been defined in node.csv" << endl;
                                }

                                error_count++;
                                // destination zone has not been defined, skipped.
                                continue;
                            }

                            int from_zone_seq_no = 0;
                            int to_zone_seq_no = 0;
                            from_zone_seq_no = assignment.g_zoneid_to_zone_seq_no_mapping[origin_zone];
                            to_zone_seq_no = assignment.g_zoneid_to_zone_seq_no_mapping[destination_zone];

                            // encounter return
                            if (demand_value < -99)
                                break;

                            assignment.total_demand[agent_type_no][demand_period_no] += demand_value;
                            assignment.g_column_pool[from_zone_seq_no][to_zone_seq_no][agent_type_no][demand_period_no].od_volume += demand_value;
                            assignment.total_demand_volume += demand_value;
                            assignment.g_origin_demand_array[from_zone_seq_no][agent_type_no][demand_period_no] += demand_value;

                            // we generate vehicles here for each OD data line
                            if (line_no <= 5)  // read only one line, but has not reached the end of the line
                                dtalog.output() << "o_zone_id:" << origin_zone << ", d_zone_id: " << destination_zone << ", value = " << demand_value << endl;

                            line_no++;
                        }  // scan lines
                        parser.CloseCSVFile();


                        dtalog.output() << "total_demand_volume is " << assignment.total_demand_volume << endl << endl;
                    }
                    else
                    {
                        // open file
                        dtalog.output() << "Error: File " << file_name << " cannot be opened.\n It might be currently used and locked by EXCEL." << endl;
                        g_ProgramStop();
                    }
                }


                else if (format_type.compare("agent_csv") == 0 || format_type.compare("routing_policy") == 0)
                {
                    CCSVParser parser;
                    if (parser.OpenCSVFile(file_name, false))
                    {
                        int total_demand_in_demand_file = 0;
                        // read agent file line by line,

                        int agent_id, o_zone_id, d_zone_id;
                        string agent_type, demand_period;

                        std::vector <int> node_sequence;

                        while (parser.ReadRecord())
                        {
                            total_demand_in_demand_file++;
                            if (total_demand_in_demand_file % 1000 == 0)
                                dtalog.output() << "demand_volume is " << total_demand_in_demand_file << endl;

                            parser.GetValueByFieldName("agent_id", agent_id);
                            parser.GetValueByFieldName("o_zone_id", o_zone_id);
                            parser.GetValueByFieldName("d_zone_id", d_zone_id);

                            CAgentPath agent_path_element;

                            int o_node_id;
                            int d_node_id;
                            float routing_ratio = 0;

                            parser.GetValueByFieldName("path_id", agent_path_element.path_id);
                            parser.GetValueByFieldName("o_node_id", o_node_id);
                            parser.GetValueByFieldName("d_node_id", d_node_id);

                            agent_path_element.o_node_no = assignment.g_node_id_to_seq_no_map[o_node_id];
                            agent_path_element.d_node_no = assignment.g_node_id_to_seq_no_map[d_node_id];

                            int from_zone_seq_no = 0;
                            int to_zone_seq_no = 0;
                            from_zone_seq_no = assignment.g_zoneid_to_zone_seq_no_mapping[o_zone_id];
                            to_zone_seq_no = assignment.g_zoneid_to_zone_seq_no_mapping[d_zone_id];

                            if (format_type.compare("agent_csv") == 0)
                            {
                                parser.GetValueByFieldName("volume", agent_path_element.volume);

                                assignment.total_demand[agent_type_no][demand_period_no] += agent_path_element.volume;
                                assignment.g_column_pool[from_zone_seq_no][to_zone_seq_no][agent_type_no][demand_period_no].od_volume += agent_path_element.volume;
                                assignment.total_demand_volume += agent_path_element.volume;
                                assignment.g_origin_demand_array[from_zone_seq_no][agent_type_no][demand_period_no] += agent_path_element.volume;
                            }

                            if(format_type.compare("routing_policy") == 0)
                            {
                                parser.GetValueByFieldName("ratio", routing_ratio);

                                float ODDemandVolume = assignment.g_column_pool[from_zone_seq_no][to_zone_seq_no][agent_type_no][demand_period_no].od_volume;

                                if (ODDemandVolume <= 0.001)
                                {
                                    dtalog.output() << "ODDemandVolume <= 0.001 for OD pair" << o_zone_id << "->" << d_zone_id  << "in routing policy file " << file_name.c_str() << ". Please check" <<  endl;
                                    g_ProgramStop();
                                }

                                agent_path_element.volume = routing_ratio * ODDemandVolume;
                                //assignment.g_origin_demand_array[from_zone_seq_no][agent_type_no][demand_period_no] should be loaded first
                            }

                            //apply for both agent csv and routing policy
                            assignment.g_column_pool[from_zone_seq_no][to_zone_seq_no][agent_type_no][demand_period_no].bfixed_route = true;

                            bool bValid = true;

                            string path_node_sequence;
                            parser.GetValueByFieldName("node_sequence", path_node_sequence);

                            std::vector<int> node_id_sequence;

                            g_ParserIntSequence(path_node_sequence, node_id_sequence);

                            std::vector<int> node_no_sequence;
                            std::vector<int> link_no_sequence;

                            int node_sum = 0;
                            for (int i = 0; i < node_id_sequence.size(); ++i)
                            {
                                if (assignment.g_node_id_to_seq_no_map.find(node_id_sequence[i]) == assignment.g_node_id_to_seq_no_map.end())
                                {
                                    bValid = false;
                                    //has not been defined
                                    continue;
                                    // warning
                                }

                                int internal_node_seq_no = assignment.g_node_id_to_seq_no_map[node_id_sequence[i]];  // map external node number to internal node seq no.
                                node_no_sequence.push_back(internal_node_seq_no);

                                node_sum += internal_node_seq_no;
                                if (i >= 1)
                                {
                                    // check if a link exists
                                    int link_seq_no = -1;
                                    // map external node number to internal node seq no.
                                    int prev_node_seq_no = assignment.g_node_id_to_seq_no_map[node_id_sequence[i - 1]];
                                    int current_node_no = node_no_sequence[i];

                                    if (g_node_vector[prev_node_seq_no].m_to_node_2_link_seq_no_map.find(current_node_no) != g_node_vector[prev_node_seq_no].m_to_node_2_link_seq_no_map.end())
                                    {
                                        link_seq_no = g_node_vector[prev_node_seq_no].m_to_node_2_link_seq_no_map[node_no_sequence[i]];
                                        link_no_sequence.push_back(link_seq_no);
                                    }
                                    else
                                        bValid = false;
                                }
                            }

                            if (bValid)
                            {
                                agent_path_element.node_sum = node_sum; // pointer to the node sum based path node sequence;
                                agent_path_element.path_link_sequence = link_no_sequence;

                                CColumnVector* pColumnVector = &(assignment.g_column_pool[from_zone_seq_no][to_zone_seq_no][agent_type_no][demand_period_no]);

                                 // we cannot find a path with the same node sum, so we need to add this path into the map,
                                if (pColumnVector->path_node_sequence_map.find(node_sum) == pColumnVector->path_node_sequence_map.end())
                                {
                                    // add this unique path
                                    int path_count = pColumnVector->path_node_sequence_map.size();
                                    pColumnVector->path_node_sequence_map[node_sum].path_seq_no = path_count;
                                    pColumnVector->path_node_sequence_map[node_sum].path_volume = 0;
                                    //assignment.g_column_pool[m_origin_zone_seq_no][destination_zone_seq_no][agent_type][tau].time = m_label_time_array[i];
                                    //assignment.g_column_pool[m_origin_zone_seq_no][destination_zone_seq_no][agent_type][tau].path_node_sequence_map[node_sum].path_distance = m_label_distance_array[i];
                                    pColumnVector->path_node_sequence_map[node_sum].path_toll = 0;

                                    pColumnVector->path_node_sequence_map[node_sum].AllocateVector(node_no_sequence, link_no_sequence, false);
                                }

                                pColumnVector->path_node_sequence_map[node_sum].path_volume += agent_path_element.volume;
                            }
                        }
                    }
                    else
                    {
                        //open file
                        dtalog.output() << "Error: File " << file_name << " cannot be opened.\n It might be currently used and locked by EXCEL." << endl;
                        g_ProgramStop();
                    }
                }
                else
                {
                    dtalog.output() << "Error: format_type = " << format_type << " is not supported. Currently STALite supports format such as column and agent_csv." << endl;
                    g_ProgramStop();
                }
            }
        }
    }
    //else
    {  // default

        bool bFileReady = false;
        int error_count = 0;


        int demand_period_no = 0;

        string file_name;
        file_name = "demand_geo.csv";
        // read the file formaly after the test.
        CCSVParser parser;
        if (parser.OpenCSVFile(file_name, false))
        {
            bFileReady = true;
            int line_no = 0;

            while (parser.ReadRecord())
            {
                int origin_zone = -1;
                int destination_zone = -1;
                float demand_value = 0;
                string agent_type;
                parser.GetValueByFieldName("agent_type", agent_type);
                parser.GetValueByFieldName("o_zone_id", origin_zone);
                parser.GetValueByFieldName("d_zone_id", destination_zone);
                parser.GetValueByFieldName("volume", demand_value);

                int agent_type_no = assignment.agent_type_2_seqno_mapping[agent_type];
                
                if (origin_zone <= -1)
                {
                    if (line_no == 1)  // read only one line, but has not reached the end of the line
                    {
                        dtalog.output() << endl << "Error: Only one line has been read from file. Are there multiple columns of demand type in file " << file_name << " per line?" << endl;
                        cout << endl << "Error: Only one line has been read from file. Are there multiple columns of demand type in file " << file_name << " per line?" << endl;
                        g_ProgramStop();
                    }
                    break;
                }

                if (assignment.g_zoneid_to_zone_seq_no_mapping.find(origin_zone) == assignment.g_zoneid_to_zone_seq_no_mapping.end())
                {
                    if (error_count < 10)
                    {
                        dtalog.output() << endl << "Warning: origin zone " << origin_zone << "  has not been defined in node.csv" << endl;
                        cout << endl << "Warning: origin zone " << origin_zone << "  has not been defined in node.csv" << endl;
                    }

                    error_count++;
                    // origin zone has not been defined, skipped.
                    continue;
                }

                if (assignment.g_zoneid_to_zone_seq_no_mapping.find(destination_zone) == assignment.g_zoneid_to_zone_seq_no_mapping.end())
                {
                    if (error_count < 10)
                    {
                        dtalog.output() << endl << "Warning: destination zone " << destination_zone << "  has not been defined in node.csv" << endl;
                        cout << endl << "Warning: destination zone " << destination_zone << "  has not been defined in node.csv" << endl;
                    }

                    error_count++;
                    // destination zone has not been defined, skipped.
                    continue;
                }

                int from_zone_seq_no = 0;
                int to_zone_seq_no = 0;
                from_zone_seq_no = assignment.g_zoneid_to_zone_seq_no_mapping[origin_zone];
                to_zone_seq_no = assignment.g_zoneid_to_zone_seq_no_mapping[destination_zone];

                // encounter return
                if (demand_value < -99)
                    break;

                assignment.total_demand[agent_type_no][demand_period_no] += demand_value;
                assignment.g_column_pool[from_zone_seq_no][to_zone_seq_no][agent_type_no][demand_period_no].od_volume += demand_value;
                assignment.total_demand_volume += demand_value;
                assignment.g_origin_demand_array[from_zone_seq_no][agent_type_no][demand_period_no] += demand_value;

                // we generate vehicles here for each OD data line
                if (line_no <= 5)  // read only one line, but has not reached the end of the line
                    dtalog.output() << "o_zone_id:" << origin_zone << ", d_zone_id: " << destination_zone << ", value = " << demand_value << endl;

                line_no++;
            }  // scan lines
            parser.CloseCSVFile();


            dtalog.output() << "total_demand_volume is " << assignment.total_demand_volume << endl << endl;
        }
    }
}

void g_ReadOutputFileConfiguration(Assignment& assignment)
{
    dtalog.output() << "Step 1.9: Reading file section [output_file_configuration] in setting.csv..." << endl;

    cout << "Step 1.8: Reading file section [output_file_configuration] in setting.csv..." << endl;

    CCSVParser parser;
    parser.IsFirstLineHeader = false;
    if (parser.OpenCSVFile("settings.csv", false))
    {
        while (parser.ReadRecord_Section())
        {
            if (parser.SectionName == "[output_file_configuration]")
            {
                parser.GetValueByFieldName("accessibility_output", assignment.accessibility_output, false, false);
            }
        }

    parser.CloseCSVFile();
    }
}

void g_AddNewVirtualConnectorLink(int internal_from_node_seq_no, int internal_to_node_seq_no, int zone_seq_no = -1)
{
    // create a link object
    CLink link;

    link.from_node_seq_no = internal_from_node_seq_no;
    link.to_node_seq_no = internal_to_node_seq_no;
    link.link_seq_no = assignment.g_number_of_links;
    link.to_node_seq_no = internal_to_node_seq_no;
    //virtual connector
    link.link_type = -1;
    link.link_id = "virtual" + std::to_string(link.link_seq_no);

    //only for outgoing connectors
    link.zone_seq_no_for_outgoing_connector = zone_seq_no;

    //BPR
    link.traffic_flow_code = 0;

    link.spatial_capacity_in_vehicles = 99999;
    link.number_of_lanes = 10;
    link.lane_capacity = 999999;
    link.link_spatial_capacity = 99999;
    link.length = 0.00001;

    for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
    {
        //setup default values
        link.VDF_period[tau].capacity = 99999;
        // 60.0 for 60 min per hour
        link.VDF_period[tau].FFTT = 0;
        link.VDF_period[tau].alpha = 0;
        link.VDF_period[tau].beta = 0;

        link.TDBaseTT[tau] = 0;
        link.TDBaseCap[tau] = 99999;
        link.travel_time_per_period[tau] = 0;
    }

    // add this link to the corresponding node as part of outgoing node/link
    g_node_vector[internal_from_node_seq_no].m_outgoing_link_seq_no_vector.push_back(link.link_seq_no);
    // add this link to the corresponding node as part of outgoing node/link
    g_node_vector[internal_to_node_seq_no].m_incoming_link_seq_no_vector.push_back(link.link_seq_no);
     // add this link to the corresponding node as part of outgoing node/link
    g_node_vector[internal_from_node_seq_no].m_to_node_seq_no_vector.push_back(link.to_node_seq_no);
    // add this link to the corresponding node as part of outgoing node/link
    g_node_vector[internal_from_node_seq_no].m_to_node_2_link_seq_no_map[link.to_node_seq_no] = link.link_seq_no;

    g_link_vector.push_back(link);

    assignment.g_number_of_links++;
}

//source: https://stackoverflow.com/questions/217578/how-can-i-determine-whether-a-2d-point-is-within-a-polygon/2922778#2922778 

int g_PtInPoly(std::vector<double> vertx, std::vector<double> verty, double testx, double testy)
{
    int nvert = vertx.size();
    int i, j, c = 0;
    for (i = 0, j = nvert - 1; i < nvert; j = i++) {
        if (((verty[i] > testy) != (verty[j] > testy)) &&
            (testx < (vertx[j] - vertx[i]) * (testy - verty[i]) / (verty[j] - verty[i]) + vertx[i]))
            c = !c;
    }
    return c;
}

void g_ReadSubarea(Assignment& assignment)
{
    bool bNodeNonExistError = false;
    int lineno = 0;

    CCSVParser parser;

    if (parser.OpenCSVFile("subarea.csv", false))
    {
        int i = 0;
        while (parser.ReadRecord())
        {
            int subarea_id;
            assignment.g_subarea_mode = 1;   // focusing approach

            if (parser.GetValueByFieldName("subarea_id", subarea_id) == false)
            {

            }

            string name;


            std::vector<CCoordinate> CoordinateVector;
            string geo_string;

            if (parser.GetValueByFieldName("geometry", geo_string))
            {
                // overwrite when the field "geometry" exists
                CGeometry geometry(geo_string);
                CoordinateVector = geometry.GetCoordinateList();

                int si;
                for (si = 0; si < CoordinateVector.size(); si++)
                {

                    assignment.m_subarea_vec_x.push_back(CoordinateVector[si].X);
                    assignment.m_subarea_vec_y.push_back(CoordinateVector[si].Y);

                }

            }
            break; // only once record
            lineno++;
        }
        parser.CloseCSVFile();
    }
    else
    {
        return;
    }
        dtalog.output() << "subarea contains # of shape points = " << assignment.m_subarea_vec_x.size() << endl;
//        m_ZoneDataLoadingStatus.Format("%d zone info records are loaded from file %s.", lineno, lpszFileName);


        CCSVParser parser_node;

        if (parser_node.OpenCSVFile("node.csv", true))
        {
            while (parser_node.ReadRecord())  // if this line contains [] mark, then we will also read field headers.
            {
                int node_id;
                if (!parser_node.GetValueByFieldName("node_id", node_id))
                    continue;
                double x, y;
                parser_node.GetValueByFieldName("x_coord", x, true, false);
                parser_node.GetValueByFieldName("y_coord", y, true, false);

                int b_in_subarea = g_PtInPoly(assignment.m_subarea_vec_x, assignment.m_subarea_vec_y, x, y);
                if (b_in_subarea == 1)
                    assignment.m_subarea_node_flag_map[node_id] = 1;
            }

            parser_node.CloseCSVFile();
            dtalog.output() << "# of within_subarea node = " << assignment.m_subarea_node_flag_map.size() << endl;
       }


        CCSVParser parser_link;

        if (parser_link.OpenCSVFile("link.csv", true))
        {
            while (parser_link.ReadRecord())  // if this line contains [] mark, then we will also read field headers.
            {
                int from_node_id;
                if (!parser_link.GetValueByFieldName("from_node_id", from_node_id))
                    continue;

                int to_node_id;
                if (!parser_link.GetValueByFieldName("to_node_id", to_node_id))
                    continue;
                
                int from_node_flag = 0;
                int to_node_flag = 0;

                if(assignment.g_subarea_mode==1)
                {
                    if (assignment.m_subarea_node_flag_map.find(from_node_id) != assignment.m_subarea_node_flag_map.end())
                        from_node_flag = 1;

                    if (assignment.m_subarea_node_flag_map.find(to_node_id) != assignment.m_subarea_node_flag_map.end())
                        to_node_flag = 1;

                    if (from_node_flag && to_node_flag)
                    {
                        assignment.m_subarea_node_id_map[from_node_id] = 1;
                        assignment.m_subarea_node_id_map[to_node_id] = 1;
                    }
                    if (from_node_flag==1 && to_node_flag==0)
                    {
                        assignment.m_subarea_node_id_map[from_node_id] = 1;
                        assignment.m_subarea_node_id_map[to_node_id] = 1;
                        assignment.m_subarea_boundary_node_map[to_node_id]= 1;
                        dtalog.output() << " boundary out node = " << to_node_id << endl;
                    }
                    if (from_node_flag == 0 && to_node_flag == 1)
                    {
                        assignment.m_subarea_node_id_map[from_node_id] = 1;
                        assignment.m_subarea_node_id_map[to_node_id] = 1;
                        assignment.m_subarea_boundary_node_map[from_node_id] = 1;
                        dtalog.output() << " boundary in node = " << from_node_id << endl;
                    }
                }
                //this link will be included in the subarea

             }

            parser_link.CloseCSVFile();
            }

        if (assignment.g_subarea_mode == 1)
        { 
            dtalog.output() << "# of cross_subarea node = " << assignment.m_subarea_boundary_node_map.size() << endl;
        }


}
void g_OutputTMCFiles(bool bReadingDataReady)
{
    FILE* g_pFileTMCLink = fopen("TMC_link.csv", "w");

    if (g_pFileTMCLink != NULL)
    {
        fprintf(g_pFileTMCLink, "link_id,tmc,tmc_corridor_name,tmc_corridor_id,tmc_road_order,tmc_road_sequence,tmc_road,tmc_direction,tmc_intersection,tmc_highest_speed,link_no,from_node_id,to_node_id,link_type,");
        fprintf(g_pFileTMCLink, "link_type_code, FT, AT, nlanes,length, free_speed, tmc_reference_speed, tmc_mean_speed, tmc_volume, capacity, kc, ");
        fprintf(g_pFileTMCLink, "AM_vc,AM_QHF,AM_QDF_n,AM_QDF_s,AM_t0, AM_t3, AM_P, AM_Assign_V,AM_Assign_VMT,AM_Assign_VHT, AM_Assign_VDT, AM_Assign_VCDT,AM_D, AM_DC_ratio, AM_mu, AM_vu, AM_vf_reference, AM_v_mean, AM_t2_v, AM_t2_queue, AM_gamma, AM_DTASpeed1, AM_DTAP1, AM_DTATDSpdDiff,");
        fprintf(g_pFileTMCLink, "MD_vc,MD_QHF,MD_QDF_n,MD_QDF_s,MD_t0, MD_t3, MD_P, MD_Assign_V,MD_Assign_VMT,MD_Assign_VHT, MD_Assign_VDT, MD_Assign_VCDT,MD_D, MD_DC_ratio, MD_mu, MD_vu, MD_vf_reference, MD_v_mean, MD_t2_v, MD_t2_queue, MD_gamma, MD_DTASpeed1, MD_DTAP1, MD_DTATDSpdDiff,");
        fprintf(g_pFileTMCLink, "PM_vc,PM_QHF,PM_QDF_n,PM_QDF_s,PM_t0, PM_t3, PM_P, PM_Assign_V,PM_Assign_VMT,PM_Assign_VHT, PM_Assign_VDT, PM_Assign_VCDT,PM_D, PM_DC_ratio, PM_mu, PM_vu, PM_vf_reference, PM_v_mean, PM_t2_v, PM_t2_queue, PM_gamma,");
        fprintf(g_pFileTMCLink, "geometry, tmc_geometry, STA_speed1,STA_VOC1,STA_volume1,STA_speed2,STA_VOC2,STA_volume2,STA_speed3,STA_VOC3,STA_volume3,STA_speed4,STA_VOC4,STA_volume4,");

        //if(bReadingDataReady)
        {
            for (int t = 6 * 60; t < 20 * 60; t += 60)
            {
                int hour = t / 60;
                int minute = t - hour * 60;
                fprintf(g_pFileTMCLink, "vh%02d,", hour, minute);
            }
            for (int t = 6 * 60; t < 20 * 60; t += 60)
            {
                int hour = t / 60;
                int minute = t - hour * 60;

                fprintf(g_pFileTMCLink, "evh%02d,", hour);
            }

            for (int t = 6 * 60; t < 20 * 60; t += 60)
            {
                int hour = t / 60;
                int minute = t - hour * 60;

                fprintf(g_pFileTMCLink, "evh%02ddiff,", hour);
            }

            fprintf(g_pFileTMCLink, "evhMAE,evhMAPE,evhRMSE,");

            for (int t = 6 * 60; t < 20 * 60; t += 60)
            {
                int hour = t / 60;
                int minute = t - hour * 60;
                fprintf(g_pFileTMCLink, "vhr%02d,", hour, minute);
            }
            for (int t = 6 * 60; t < 20 * 60; t += 60)
            {
                int hour = t / 60;
                int minute = t - hour * 60;

                fprintf(g_pFileTMCLink, "evhr%02d,", hour);
            }
            fprintf(g_pFileTMCLink, "STA_AM_AE,STA_MD_AE,STA_PM_AE,STA_MAE,STA_MAPE,STA_RMSE,");


            for (int t = 6 * 60; t < 20 * 60; t += 5)
        {
            int hour = t / 60;
            int minute = t - hour * 60;

            fprintf(g_pFileTMCLink, "v%02d:%02d,", hour, minute);
        }

            for (int t = 0 * 60; t < 24 * 60; t += 15)
            {
                int hour = t / 60;
                int minute = t - hour * 60;

                fprintf(g_pFileTMCLink, "v%02d:%02d,", hour, minute);
            }
            
            for (int t = 6 * 60; t < 20 * 60; t += 15)
        {
            int hour = t / 60;
            int minute = t - hour * 60;

            fprintf(g_pFileTMCLink, "vr%02d:%02d,", hour, minute);
        }


        for (int t = 6 * 60; t < 20 * 60; t += 15)
        {
            int hour = t / 60;
            int minute = t - hour * 60;

            fprintf(g_pFileTMCLink, "q%02d:%02d,", hour, minute);
        }
        for (int t = 6 * 60; t < 20 * 60; t += 15)
        {
            int hour = t / 60;
            int minute = t - hour * 60;

            fprintf(g_pFileTMCLink, "q%02d:%02d,", hour, minute);
        }
        /// <summary>
        /// / etimation 
        /// </summary>
        for (int t = 6 * 60; t < 20 * 60; t += 15)
        {
            int hour = t / 60;
            int minute = t - hour * 60;

            fprintf(g_pFileTMCLink, "ev%02d:%02d,", hour, minute);
        }
        }


        fprintf(g_pFileTMCLink, "\n");

        std::map<_int64, int> TMC_long_id_mapping;  // this is used to mark if this cell_id has been identified or not

        //// sort data records
        for (int i = 0; i < g_link_vector.size(); i++)
        {
            if (g_link_vector[i].link_id == "5666")
            {
             int idebug = 1;
            }

            if (g_link_vector[i].TMC_code.size() > 0)
            {
               

                _int64 TMC_long_key = (g_link_vector[i].tmc_corridor_id * 10000 + g_link_vector[i].tmc_road_sequence) * 10 + g_link_vector[i].link_seq_no;
                TMC_long_id_mapping[TMC_long_key] = g_link_vector[i].link_seq_no;
            }
        }


        std::map<_int64, int>::iterator it;

        for (it = TMC_long_id_mapping.begin(); it != TMC_long_id_mapping.end(); ++it)
        {
            int i = it->second;
            float highest_speed = g_link_vector[i].free_speed;
        
            if (assignment.m_TMClink_map.find(g_link_vector[i].TMC_code) != assignment.m_TMClink_map.end())
            {
                int tmc_index = assignment.m_TMClink_map[g_link_vector[i].TMC_code];

                highest_speed = g_TMC_vector[tmc_index].GetHighestSpeed();
                //    if (g_TMC_vector[tmc_index].bWithSensorSpeedData == false)  // with reading speed data only
                //        continue;
            }


                float free_speed = g_link_vector[i].free_speed;
                g_link_vector[i].UpdateKC(g_link_vector[i].free_speed);


                fprintf(g_pFileTMCLink, "%s,%s,%s,%d,%d,%d,%s,%s,%s,%.2f,%d,%d,%d,%d,%s,%d,%d,%d,%f,%f,%f,%f,%f,%f,%f,",
                    g_link_vector[i].link_id.c_str(),
                    g_link_vector[i].TMC_code.c_str(),
                    g_link_vector[i].tmc_corridor_name.c_str(),
                    g_link_vector[i].tmc_corridor_id,
                    g_link_vector[i].tmc_road_order,
                    g_link_vector[i].tmc_road_sequence,
                    g_link_vector[i].tmc_road.c_str(),
                    g_link_vector[i].tmc_direction.c_str(),
                    g_link_vector[i].tmc_intersection.c_str(),
                    highest_speed,
                    g_link_vector[i].link_seq_no,
                    g_node_vector[g_link_vector[i].from_node_seq_no].node_id,
                    g_node_vector[g_link_vector[i].to_node_seq_no].node_id,
                    g_link_vector[i].link_type,
                    g_link_vector[i].link_type_code.c_str(),
                    g_link_vector[i].FT,
                    g_link_vector[i].AT,
                    g_link_vector[i].number_of_lanes,
                    g_link_vector[i].length,
                    g_link_vector[i].free_speed,
                    g_link_vector[i].tmc_reference_speed,
                    g_link_vector[i].tmc_mean_speed,
                    g_link_vector[i].tmc_volume,
                    g_link_vector[i].lane_capacity,
                    g_link_vector[i].kc

                );


                if (assignment.m_TMClink_map.find(g_link_vector[i].TMC_code) == assignment.m_TMClink_map.end())
                {
                    fprintf(g_pFileTMCLink, "\n");
                    continue;

                }
                float FD_vc = 0;
                float obs_t0_in_hour = -1;
                float obs_t3_in_hour = -1;
                float mean_congestion_speed = 0;
                float congestion_dc_ratio = 0;
                float obs_P_in_hour = 0;
                float mean_congestion_mu = 0;
                highest_speed = 0;
                float mean_speed = 0;
                float congestion_D = 0;
                float assignment_D = 0;
                float assignment_VHT = 0;
                float assignment_VMT = 0;
                float assignment_VDT = 0;
                float assignment_VCDT = 0;

                float t2_speed = 0;
                float t2_queue = 0;
                float gamma = 0;
                float assign_period_start_time_in_hour = 6;
                float assign_period_end_time_in_hour = 11;
                int starting_time_in_hour = 6;
                int ending_time_in_hour = 10;
                float t2 = 7.5;
                float DTASpeed = -1;
                float DTAP = -1;
                float TDAvgSpeedDiff=-1;

                float STA_AM_MAE = -1;
                float STA_MD_MAE = -1;
                float STA_PM_MAE = -1;
                float STA_AM_APE = -1;
                float STA_MD_APE = -1;
                float STA_PM_APE = -1;
                // P D analysis


                if (assignment.m_TMClink_map.find(g_link_vector[i].TMC_code) != assignment.m_TMClink_map.end())
                {
                    int tmc_index = assignment.m_TMClink_map[g_link_vector[i].TMC_code];

                    if (g_link_vector[i].TMC_code == "115-04402")
                    {
                        int idebug = 1;
                    }
                    int peak_no = 0;  //AM

                    assign_period_start_time_in_hour = 6;
                    assign_period_end_time_in_hour = 10;
                    starting_time_in_hour = 6;
                    ending_time_in_hour = 10;


                    CLink* pLink = &(g_link_vector[i]);

                    if(pLink->link_id == "213")
                        int debug_flag = 1;
                    
                    obs_P_in_hour = g_TMC_vector[tmc_index].ScanCongestionDuration(peak_no, starting_time_in_hour,
                        ending_time_in_hour, FD_vc, obs_t0_in_hour, obs_t3_in_hour, pLink, assignment_D, assignment_VMT, assignment_VHT, assignment_VDT, assignment_VCDT,congestion_D,
                        congestion_dc_ratio, mean_congestion_mu, mean_congestion_speed, highest_speed, mean_speed, t2_speed, t2_queue, gamma);

                    float V = assignment_D;
                    float laneCapacity = pLink->lane_capacity;
                    float vf = highest_speed;
                    float vcd = pLink->vc;
                    float vct = pLink->vc;
                    t2 = 7.75;  // calibrated based on 101 data

                    g_TMC_vector[tmc_index].PerformEstimation(false,pLink,peak_no,assign_period_start_time_in_hour, assign_period_end_time_in_hour, t2, V, laneCapacity, vf, vcd, vct, DTASpeed, DTAP, TDAvgSpeedDiff);

                    STA_AM_MAE = fabs(g_link_vector[i].VDF_STA_speed[1]-mean_speed);
                    STA_AM_APE = STA_AM_MAE/max(1,mean_speed);

                    float QHF = 3;
                    if(congestion_D>1 && obs_P_in_hour >= 1)
                        QHF = V / congestion_D;

                    fprintf(g_pFileTMCLink, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,",
                        FD_vc, 
                        QHF, g_TMC_vector[tmc_index].est_QDF_n [peak_no], g_TMC_vector[tmc_index].est_QDF_s[peak_no],
                        obs_t0_in_hour, obs_t3_in_hour, obs_P_in_hour, assignment_D, assignment_VMT, assignment_VHT, assignment_VDT, assignment_VCDT,
                        congestion_D, congestion_dc_ratio, mean_congestion_mu, mean_congestion_speed,
                        highest_speed, mean_speed,
                        t2_speed, t2_queue, gamma, DTASpeed, DTAP, TDAvgSpeedDiff);

                    ///////////////////////////////////////  noon
                    peak_no = 1;

                    assign_period_start_time_in_hour = 10;
                    assign_period_end_time_in_hour = 14;
                    starting_time_in_hour = 12;
                    ending_time_in_hour = 13;
                    t2 = 12.5;  //12 noon

                    obs_P_in_hour = g_TMC_vector[tmc_index].ScanCongestionDuration(peak_no, starting_time_in_hour,
                        ending_time_in_hour, FD_vc, obs_t0_in_hour, obs_t3_in_hour, pLink, assignment_D, assignment_VMT, assignment_VHT, assignment_VDT, assignment_VCDT,congestion_D,
                        congestion_dc_ratio, mean_congestion_mu, mean_congestion_speed, highest_speed, mean_speed, t2_speed, t2_queue, gamma);

                    V = assignment_D;
                    vcd = pLink->vc;
                    vct = pLink->vc;

                    g_TMC_vector[tmc_index].PerformEstimation(false,pLink, peak_no, assign_period_start_time_in_hour, assign_period_end_time_in_hour, t2, V, laneCapacity, vf, vcd, vct,DTASpeed, DTAP, TDAvgSpeedDiff);
                    STA_MD_MAE = fabs(g_link_vector[i].VDF_STA_speed[1] - mean_speed);
                    STA_MD_APE = STA_MD_MAE / max(1, mean_speed);

                    QHF = 3;
                    if (congestion_D > 1 && obs_P_in_hour >= 1)
                        QHF = V / congestion_D;

                    fprintf(g_pFileTMCLink, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,",
                        FD_vc,
                        QHF, g_TMC_vector[tmc_index].est_QDF_n[peak_no], g_TMC_vector[tmc_index].est_QDF_s[peak_no],
                        obs_t0_in_hour, obs_t3_in_hour, obs_P_in_hour, assignment_D, assignment_VMT, assignment_VHT, assignment_VDT, assignment_VCDT,
                        congestion_D, congestion_dc_ratio, mean_congestion_mu, mean_congestion_speed,
                        highest_speed, mean_speed,
                        t2_speed, t2_queue, gamma, DTASpeed, DTAP, TDAvgSpeedDiff);

                    ///////////////////////////////////////  PM
                    peak_no = 2;

                    assign_period_start_time_in_hour = 14;
                    assign_period_end_time_in_hour = 20;
                    starting_time_in_hour = 15;
                    ending_time_in_hour = 18;
                    t2 = 17;  //5pm

                    obs_P_in_hour = g_TMC_vector[tmc_index].ScanCongestionDuration(peak_no, starting_time_in_hour,
                        ending_time_in_hour, FD_vc, obs_t0_in_hour, obs_t3_in_hour, pLink, assignment_D, assignment_VMT, assignment_VHT, assignment_VDT, assignment_VCDT,congestion_D,
                        congestion_dc_ratio, mean_congestion_mu, mean_congestion_speed, highest_speed, mean_speed, t2_speed, t2_queue, gamma);

                    V = assignment_D;
                    vcd = pLink->vc;
                    vct = pLink->vc;

                    g_TMC_vector[tmc_index].PerformEstimation(false,pLink, peak_no, assign_period_start_time_in_hour, assign_period_end_time_in_hour, t2, V, laneCapacity, vf, vcd, vct, DTASpeed, DTAP, TDAvgSpeedDiff);
                    STA_PM_MAE = fabs(g_link_vector[i].VDF_STA_speed[3] - mean_speed);
                    STA_PM_APE = STA_PM_MAE / max(1, mean_speed);
                    
                    QHF = 3;
                    if (congestion_D > 1 && obs_P_in_hour>=1)
                        QHF = max(1, V / congestion_D);

                    fprintf(g_pFileTMCLink, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,",
                        FD_vc, 
                        QHF, g_TMC_vector[tmc_index].est_QDF_n[peak_no], g_TMC_vector[tmc_index].est_QDF_s[peak_no],
                        obs_t0_in_hour, obs_t3_in_hour, obs_P_in_hour, assignment_D, assignment_VMT, assignment_VHT, assignment_VDT, assignment_VCDT,
                        congestion_D, congestion_dc_ratio, mean_congestion_mu, mean_congestion_speed,
                        highest_speed, mean_speed,
                        t2_speed, t2_queue, gamma, DTASpeed, DTAP, TDAvgSpeedDiff);


                }
                else
                {
                // no data AM

                fprintf(g_pFileTMCLink, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,",
                    FD_vc, -1, -1, -1, obs_t0_in_hour, obs_t3_in_hour, obs_P_in_hour, assignment_D, assignment_VMT, assignment_VHT, assignment_VDT, assignment_VCDT,
                    congestion_D,congestion_dc_ratio, mean_congestion_mu,mean_congestion_speed,
                    highest_speed, mean_speed,
                    t2_speed, t2_queue, gamma, DTASpeed, DTAP, TDAvgSpeedDiff);
                // no data MD
                fprintf(g_pFileTMCLink, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,",
                    FD_vc, -1, -1, -1, obs_t0_in_hour, obs_t3_in_hour, obs_P_in_hour, assignment_D, assignment_VMT, assignment_VHT, assignment_VDT, assignment_VCDT,
                    congestion_D, congestion_dc_ratio, mean_congestion_mu, mean_congestion_speed,
                    highest_speed, mean_speed,
                    t2_speed, t2_queue, gamma, DTASpeed, DTAP, TDAvgSpeedDiff);

                // no data PM
                fprintf(g_pFileTMCLink, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,",
                    FD_vc, -1, -1, -1, obs_t0_in_hour, obs_t3_in_hour, obs_P_in_hour, assignment_D, assignment_VMT, assignment_VHT, assignment_VDT, assignment_VCDT,
                    congestion_D, congestion_dc_ratio, mean_congestion_mu, mean_congestion_speed,
                    highest_speed, mean_speed,
                    t2_speed, t2_queue, gamma, DTASpeed, DTAP, TDAvgSpeedDiff);
                }
                // vc, t0, t3, P, D, mu,

                fprintf(g_pFileTMCLink, "\"%s\",", g_link_vector[i].geometry.c_str());

                fprintf(g_pFileTMCLink, "\"LINESTRING (");

                fprintf(g_pFileTMCLink, "%f %f,", g_link_vector[i].TMC_from.x, g_link_vector[i].TMC_from.y);
                fprintf(g_pFileTMCLink, "%f %f,", g_link_vector[i].TMC_to.x, g_link_vector[i].TMC_to.y);
                fprintf(g_pFileTMCLink, ")\"");
                fprintf(g_pFileTMCLink, ",");

               /* if (bReadingDataReady)*/
                {
                if (assignment.m_TMClink_map.find(g_link_vector[i].TMC_code) != assignment.m_TMClink_map.end())
                {
                    int tmc_index = assignment.m_TMClink_map[g_link_vector[i].TMC_code];

                    for (int tau = 1; tau <= 4; tau++)
                    {
                        fprintf(g_pFileTMCLink, "%f,%f,%f,",
                            g_link_vector[i].VDF_STA_speed[tau], g_link_vector[i].VDF_STA_VOC[tau], g_link_vector[i].VDF_STA_volume[tau]);
                    }

                    double ObsSpeed[25], EstSpeed[25], EstSpeedDiff[25];
                    double MAE_total = 0;
                    double MAPE_total = 0;
                    double RMSE_total = 0;
                    int count_total = 0;

                    for (int t = 6 * 60; t < 20 * 60; t += 60)
                    {
                        int hour = t / 60;
                        ObsSpeed[hour] = g_TMC_vector[tmc_index].GetAvgHourlySpeed(t);
                        fprintf(g_pFileTMCLink, "%.1f,", ObsSpeed[hour]);
                    }
                    for (int t = 6 * 60; t < 20 * 60; t += 60)
                    {
                        int hour = t / 60;
                        EstSpeed[hour] = g_TMC_vector[tmc_index].GetEstHourlySpeed(t);
                        if (EstSpeed[hour] > 1)  // valid
                        {
                            EstSpeedDiff[hour] = fabs(EstSpeed[hour] - ObsSpeed[hour]);
                            MAE_total += fabs(EstSpeedDiff[hour]);
                            MAPE_total += fabs(EstSpeedDiff[hour]) / max(1, ObsSpeed[hour]);
                            RMSE_total += EstSpeedDiff[hour] * EstSpeedDiff[hour];
                            count_total += 1;
                        }
                        else
                            EstSpeedDiff[hour] = -1;
                    }
                    for (int t = 6 * 60; t < 20 * 60; t += 60)
                    {
                        int hour = t / 60;
                        EstSpeed[hour] = g_TMC_vector[tmc_index].GetEstHourlySpeed(t);

                        fprintf(g_pFileTMCLink, "%.1f,", EstSpeed[hour]);
                    }
                    for (int t = 6 * 60; t < 20 * 60; t += 60)
                    {
                        int hour = t / 60;
                        fprintf(g_pFileTMCLink, "%.1f,", EstSpeedDiff[hour]);
                    }                   
                    
                    double MSE_total = RMSE_total / max(1, count_total);

                    fprintf(g_pFileTMCLink, "%.2f,%.2f,%.2f,", MAE_total/max(1, count_total), MAPE_total / max(1, count_total) * 100, pow(MSE_total,0.5));

                    /// <summary>
                    ///
                    /// </summary>
                    /// <param name="bReadingDataReady"></param>
                    for (int t = 6 * 60; t < 20 * 60; t += 60)
                    {
                        int hour = t / 60;
                        ObsSpeed[hour] = g_TMC_vector[tmc_index].GetAvgHourlySpeed(t);

                        float speed_ratio = ObsSpeed[hour] / max(1, g_link_vector[i].TMC_highest_speed);

                        if (speed_ratio >= 1)
                            speed_ratio = 1;
                        fprintf(g_pFileTMCLink, "%.3f,", speed_ratio);
                    }
                    for (int t = 6 * 60; t < 20 * 60; t += 60)
                    {
                        int hour = t / 60;
                        EstSpeed[hour] = g_TMC_vector[tmc_index].GetEstHourlySpeed(t);
                        float speed_ratio = EstSpeed[hour] / max(1, g_link_vector[i].TMC_highest_speed);

                        if (speed_ratio >= 1)
                            speed_ratio = 1;
                        fprintf(g_pFileTMCLink, "%.3f,", speed_ratio);
                    }

                    double STA_RMSE = pow((STA_AM_MAE * STA_AM_MAE + STA_MD_MAE * STA_MD_MAE + STA_PM_MAE * STA_PM_MAE) / 3, 0.5);
                    fprintf(g_pFileTMCLink, "%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,", STA_AM_MAE, STA_MD_MAE, STA_PM_MAE, (STA_AM_MAE + STA_MD_MAE + STA_PM_MAE) / 3, (STA_AM_APE + STA_MD_APE + STA_PM_APE) / 3*100, STA_RMSE);

                    for (int t = 6 * 60; t < 20 * 60; t += 5)
                    {

                        fprintf(g_pFileTMCLink, "%.1f,", g_TMC_vector[tmc_index].GetAvgSpeed(t));
                    }

                    for (int t = 0 * 60; t < 24 * 60; t += 15)
                    {

                        fprintf(g_pFileTMCLink, "%.1f,", g_TMC_vector[tmc_index].GetAvgSpeed_15min(t));
                    }

                    for (int t = 6 * 60; t < 20 * 60; t += 15)
                    {
                        float speed_ratio = g_TMC_vector[tmc_index].GetAvgSpeed_15min(t) / max(1, g_link_vector[i].TMC_highest_speed);

                        if (speed_ratio >= 1)
                            speed_ratio = 1;    

                        fprintf(g_pFileTMCLink, "%.2f,", speed_ratio);
                    }

                    for (int t = 6 * 60; t < 20 * 60; t += 15)
                    {
                        float speed = g_TMC_vector[tmc_index].GetAvgSpeed(t);
                        double volume = g_link_vector[i].get_volume_from_speed(speed, g_link_vector[i].TMC_highest_speed);

                        fprintf(g_pFileTMCLink, "%.2f,", volume);
                    }

                    for (int t = 6 * 60; t < 20 * 60; t += 15)
                    {
                        float speed = g_TMC_vector[tmc_index].GetAvgSpeed_15min(t);
                        double volume = g_link_vector[i].get_volume_from_speed(speed, g_link_vector[i].TMC_highest_speed);

                        fprintf(g_pFileTMCLink, "%.2f,", volume);
                    }

                    // estimation 
                    for (int t = 6 * 60; t < 20 * 60; t += 15)
                    {

                        fprintf(g_pFileTMCLink, "%.1f,", g_TMC_vector[tmc_index].GetEstSpeed(t));
                    }

                //    for (int t = 6 * 60; t < 20 * 60; t += 5)
                //    {
                //        float est_speed = g_TMC_vector[tmc_index].GetEstSpeed(t);
                //        float speed_ratio = 1;

                //        if(est_speed >1)  // valid result
                //        speed_ratio = g_TMC_vector[tmc_index].GetEstSpeed(t) / max(1, g_link_vector[i].TMC_highest_speed);

                //        if (speed_ratio >= 1)
                //            speed_ratio = 1;


                //        fprintf(g_pFileTMCLink, "%.2f,", speed_ratio);
                //    }


                }
                }
                fprintf(g_pFileTMCLink, "\n");

            


        }

        fclose(g_pFileTMCLink);
    }
    else
    {
        dtalog.output() << "Error: File TMC_link.csv cannot be opened.\n It might be currently used and locked by EXCEL." << endl;
        g_ProgramStop();

    }


}
void g_OutputTMCScenarioFiles()
{

    CCSVParser parser_scenario;

    if (parser_scenario.OpenCSVFile("scenario.csv", true))
    {
        while (parser_scenario.ReadRecord())  // if this line contains [] mark, then we will also read field headers.
        {
            //string tmc;

            //for (int i = 0; i < g_link_vector.size(); i++)
            //{
            //    if (g_link_vector[i].link_type >= 0 && g_link_vector[i].TMC_code.size() > 0)
            //    {
            //        if (g_link_vector[i].TMC_code.compare(tmc)!=0)
            //        {
            //            g_link_vector[i].Scenario_evaluation_flag = true;

            //            parser_scenario.GetValueByFieldName("QDF_voc_ratio1", g_link_vector[i].Scenario_STA_VOC_Ratio[1], false);
            //            parser_scenario.GetValueByFieldName("QDF_voc_ratio2", g_link_vector[i].Scenario_STA_VOC_Ratio[2], false);
            //            parser_scenario.GetValueByFieldName("QDF_voc_ratio3", g_link_vector[i].Scenario_STA_VOC_Ratio[3], false);
            //            parser_scenario.GetValueByFieldName("QDF_voc_ratio4", g_link_vector[i].Scenario_STA_VOC_Ratio[4], false);



            //        }
            //    }
            //}

            int from_node_id;
            if (!parser_scenario.GetValueByFieldName("from_node_id", from_node_id))
                continue;

            int to_node_id;
            if (!parser_scenario.GetValueByFieldName("to_node_id", to_node_id))
                continue;


            if (assignment.g_subarea_mode >= 1)  // subarea handling 
            {
                if (assignment.g_node_id_to_seq_no_map.find(from_node_id) == assignment.g_node_id_to_seq_no_map.end())
                {
                    continue; //has not been defined
                }

                if (assignment.g_node_id_to_seq_no_map.find(to_node_id) == assignment.g_node_id_to_seq_no_map.end())
                {
                    continue; //has not been defined
                }

                if (assignment.m_subarea_node_id_map.find(from_node_id) == assignment.m_subarea_node_id_map.end())
                {
                    continue; //has not been defined
                }

                if (assignment.m_subarea_node_id_map.find(to_node_id) == assignment.m_subarea_node_id_map.end())
                {
                    continue; //has not been defined
                }
            }

            if (assignment.g_node_id_to_seq_no_map.find(from_node_id) == assignment.g_node_id_to_seq_no_map.end())
            {
                dtalog.output() << "Error: from_node_id " << from_node_id << " in file scenario.csv is not defined in node.csv." << endl;
                continue; //has not been defined
            }

            if (assignment.g_node_id_to_seq_no_map.find(to_node_id) == assignment.g_node_id_to_seq_no_map.end())
            {
                dtalog.output() << "Error: to_node_id " << to_node_id << " in file scenario.csv is not defined in node.csv." << endl;
                continue; //has not been defined
            }


            int internal_from_node_seq_no = assignment.g_node_id_to_seq_no_map[from_node_id];  // map external node number to internal node seq no.
            int internal_to_node_seq_no = assignment.g_node_id_to_seq_no_map[to_node_id];

            if (g_node_vector[internal_from_node_seq_no].m_to_node_2_link_seq_no_map.find(internal_to_node_seq_no) != g_node_vector[internal_from_node_seq_no].m_to_node_2_link_seq_no_map.end())
            {
                int link_seq_no = g_node_vector[internal_from_node_seq_no].m_to_node_2_link_seq_no_map[internal_to_node_seq_no];
                int i = link_seq_no;
                            g_link_vector[i].Scenario_evaluation_flag = true;

                            parser_scenario.GetValueByFieldName("QDF_voc_ratio1", g_link_vector[i].Scenario_STA_VOC_Ratio[1], false);
                            parser_scenario.GetValueByFieldName("QDF_voc_ratio2", g_link_vector[i].Scenario_STA_VOC_Ratio[2], false);
                            parser_scenario.GetValueByFieldName("QDF_voc_ratio3", g_link_vector[i].Scenario_STA_VOC_Ratio[3], false);
                            parser_scenario.GetValueByFieldName("QDF_voc_ratio4", g_link_vector[i].Scenario_STA_VOC_Ratio[4], false);


               }
            }
        
    }

    FILE* g_pFileTMCLink = fopen("TMC_link_scenario.csv", "w");

    if (g_pFileTMCLink != NULL)
    {
        fprintf(g_pFileTMCLink, "link_id,tmc,tmc_corridor_name,tmc_corridor_id,tmc_road_order,tmc_road_sequence,tmc_road,tmc_direction,tmc_intersection,link_no,from_node_id,to_node_id,link_type,");
        fprintf(g_pFileTMCLink, "link_type_code, FT, AT, nlanes, free_speed, tmc_reference_speed, tmc_mean_speed, tmc_volume, capacity, kc, ");
        fprintf(g_pFileTMCLink, "AM_vc, AM_t0, AM_t3, AM_P, AM_Assign_V, AM_D,AM_DC_ratio_base,AM_DC_ratio_future, AM_mu, AM_vu, AM_vf_reference, AM_v_mean, AM_t2_v, AM_t2_queue, AM_gamma, AM_DTASpeed1, AM_DTAP1, AM_DTATDSpdDiff,");
        fprintf(g_pFileTMCLink, "PM_vc, PM_t0,PM_t3,PM_P,PM_Assign_V,PM_D,PM_DC_ratio_base,PM_DC_ratio_future,PM_mu,PM_vu,PM_vf_reference,PM_v_mean, PM_t2_v, PM_t2_queue, PM_gamma, PM_DTASpeed1, PM_DTAP1, PM_DTATDSpdDiff,");
        fprintf(g_pFileTMCLink, ",geometry,tmc_geometry,reading_road_name,QDM_voc_ratio1,QDM_voc_ratio2,QDM_voc_ratio3,");

        fprintf(g_pFileTMCLink, "|,");

            for (int t = 6 * 60; t < 20 * 60; t += 60)
            {
                int hour = t / 60;
                int minute = t - hour * 60;

                fprintf(g_pFileTMCLink, "vh%02d,", hour, minute);
            }
            fprintf(g_pFileTMCLink, "|,");
            for (int t = 6 * 60; t < 20 * 60; t += 60)
            {
                int hour = t / 60;

                fprintf(g_pFileTMCLink, "pvh%02d,", hour);
            }
            fprintf(g_pFileTMCLink, "|,");

            for (int t = 6 * 60; t < 20 * 60; t += 60)
            {
                int hour = t / 60;

                fprintf(g_pFileTMCLink, "pvh%02ddiff,", hour);
            }
            fprintf(g_pFileTMCLink, "|,");

            fprintf(g_pFileTMCLink, "pvhAE,pvhAPE,");

            fprintf(g_pFileTMCLink, "|,");
            // travel time 
            for (int t = 6 * 60; t < 20 * 60; t += 60)
            {
                int hour = t / 60;

                fprintf(g_pFileTMCLink, "tth%02d,", hour);
            }
            fprintf(g_pFileTMCLink, "|,");
            for (int t = 6 * 60; t < 20 * 60; t += 60)
            {
                int hour = t / 60;

                fprintf(g_pFileTMCLink, "ptth%02d,", hour);
            }
            fprintf(g_pFileTMCLink, "|,");

            // travel time 
            for (int t = 6 * 60; t < 20 * 60; t += 60)
            {
                int hour = t / 60;

                fprintf(g_pFileTMCLink, "vtth%02d,", hour);
            }
            fprintf(g_pFileTMCLink, "|,");
            for (int t = 6 * 60; t < 20 * 60; t += 60)
            {
                int hour = t / 60;

                fprintf(g_pFileTMCLink, "ptth%02d,", hour);
            }
            fprintf(g_pFileTMCLink, "|,");

            // 5 min prediction 

            for (int t = 6 * 60; t < 20 * 60; t += 5)
            {
                int hour = t / 60;
                int minute = t - hour * 60;

                fprintf(g_pFileTMCLink, "tt%02d:%02d,", hour, minute);
            }

        fprintf(g_pFileTMCLink, "\n");

        std::map<_int64, int> TMC_long_id_mapping;  // this is used to mark if this cell_id has been identified or not

        //// sort data records
        for (int i = 0; i < g_link_vector.size(); i++)
        {
            if (g_link_vector[i].link_type >= 0 && g_link_vector[i].TMC_code.size() > 0)
            {
                _int64 TMC_long_key = (g_link_vector[i].tmc_corridor_id * 10000 + g_link_vector[i].tmc_road_sequence) * 1000000 + g_link_vector[i].link_seq_no;
                TMC_long_id_mapping[TMC_long_key] = g_link_vector[i].link_seq_no;
            }
        }


        std::map<_int64, int>::iterator it;

        for (it = TMC_long_id_mapping.begin(); it != TMC_long_id_mapping.end(); ++it)
        {
            int i = it->second;

            if (g_link_vector[i].link_type >= 0 && g_link_vector[i].TMC_code.size() > 0)
            {

                if (assignment.m_TMClink_map.find(g_link_vector[i].TMC_code) != assignment.m_TMClink_map.end())
                {
                    int tmc_index = assignment.m_TMClink_map[g_link_vector[i].TMC_code];
                    if (g_TMC_vector[tmc_index].bWithSensorSpeedData == false)  // with reading speed data only
                        continue;
                }



                float free_speed = g_link_vector[i].free_speed;
                g_link_vector[i].UpdateKC(g_link_vector[i].free_speed);

                fprintf(g_pFileTMCLink, "%s,%s,%s,%d,%d,%d,%s,%s,%s,%d,%d,%d,%d,%s,%d,%d,%d,%f,%f,%f,%f,%f,%f,",
                    g_link_vector[i].link_id.c_str(),
                    g_link_vector[i].TMC_code.c_str(),
                    g_link_vector[i].tmc_corridor_name.c_str(),
                    g_link_vector[i].tmc_corridor_id,
                    g_link_vector[i].tmc_road_order,
                    g_link_vector[i].tmc_road_sequence,
                    g_link_vector[i].tmc_road.c_str(),
                    g_link_vector[i].tmc_direction.c_str(),
                    g_link_vector[i].tmc_intersection.c_str(),
                    g_link_vector[i].link_seq_no,
                    g_node_vector[g_link_vector[i].from_node_seq_no].node_id,
                    g_node_vector[g_link_vector[i].to_node_seq_no].node_id,
                    g_link_vector[i].link_type,
                    g_link_vector[i].link_type_code.c_str(),
                    g_link_vector[i].FT,
                    g_link_vector[i].AT,
                    g_link_vector[i].number_of_lanes,
                    g_link_vector[i].free_speed,
                    g_link_vector[i].tmc_reference_speed,
                    g_link_vector[i].tmc_mean_speed,
                    g_link_vector[i].tmc_volume,
                    g_link_vector[i].lane_capacity,
                    g_link_vector[i].kc

                    //g_link_vector[i].VDF_period[0].FFTT,
                    //g_link_vector[i].VDF_period[0].capacity,
                    //g_link_vector[i].VDF_period[0].alpha,
                    //g_link_vector[i].VDF_period[0].beta,
                );

                float FD_vc = 0;
                float obs_t0_in_hour = -1;
                float obs_t3_in_hour = -1;
                float mean_congestion_speed = 0;
                float congestion_dc_ratio = 0;
                float obs_P_in_hour = 0;
                float mean_congestion_mu = 0;
                float highest_speed = 0;
                float mean_speed = 0;
                float congestion_D = 0;
                float assignment_D = 0;
                float assignment_VMT = 0;
                float assignment_VHT = 0;
                float assignment_VDT = 0;
                float assignment_VCDT = 0;
                float t2_speed = 0;
                float t2_queue = 0;
                float gamma = 0;
                float assign_period_start_time_in_hour = 6;
                float assign_period_end_time_in_hour = 11;
                int starting_time_in_hour = 6;
                int ending_time_in_hour = 9;
                float t2 = 7.5;
                float DTASpeed = -1;
                float DTAP = -1;
                float TDAvgSpeedDiff = -1;

                // P D analysis
                if (assignment.m_TMClink_map.find(g_link_vector[i].TMC_code) != assignment.m_TMClink_map.end())
                {
                    int tmc_index = assignment.m_TMClink_map[g_link_vector[i].TMC_code];

                    int peak_no = 0;  //AM

                    assign_period_start_time_in_hour = 6;
                    assign_period_end_time_in_hour = 10;
                    starting_time_in_hour = 6;
                    ending_time_in_hour = 9;

                    congestion_dc_ratio = 0;

                    CLink* pLink = &(g_link_vector[i]);

                    int from_node_id = g_node_vector[g_link_vector[i].from_node_seq_no].node_id;
                    int to_node_id = g_node_vector[g_link_vector[i].to_node_seq_no].node_id;

                    if (from_node_id == 19229 && to_node_id == 19233)
                    {
                        int idebug = 1;
                    }

                    obs_P_in_hour = g_TMC_vector[tmc_index].ScanCongestionDuration(peak_no, starting_time_in_hour,
                        ending_time_in_hour, FD_vc, obs_t0_in_hour, obs_t3_in_hour, pLink, assignment_D, assignment_VMT, assignment_VHT, assignment_VDT, assignment_VCDT,congestion_D,
                        congestion_dc_ratio, mean_congestion_mu, mean_congestion_speed, highest_speed, mean_speed, t2_speed, t2_queue, gamma);


                    float V = assignment_D ;  // apply demand flow change: AM
                    float laneCapacity = pLink->lane_capacity;
                    float vf = highest_speed;
                    float vcd = pLink->vc;
                    float vct = pLink->vc;
                    t2 = 7.5;

                    if (pLink->Scenario_evaluation_flag && g_link_vector[i].TMC_code.compare("115 + 04402") != 0)
                    {
                        int i_debug = 1;
                    }


                    if (pLink->Scenario_STA_VOC_Ratio[1] > 1.0001)
                    {
                        int i_debug = 1;
                    }

                    float AM_DC_ratio_base, AM_DC_ratio_future;
                    g_TMC_vector[tmc_index].PerformEstimation(false, pLink, peak_no, assign_period_start_time_in_hour, assign_period_end_time_in_hour, t2, V,
                        laneCapacity, vf, vcd, vct, DTASpeed, DTAP, TDAvgSpeedDiff);
                    
                    AM_DC_ratio_base = g_TMC_vector[tmc_index].est_DCRatio[peak_no];

                    V = assignment_D * pLink->Scenario_STA_VOC_Ratio[1];  // apply demand flow change: AM

                    g_TMC_vector[tmc_index].PerformEstimation(true,pLink, peak_no, assign_period_start_time_in_hour, assign_period_end_time_in_hour, t2, V,
                        laneCapacity, vf, vcd, vct, DTASpeed, DTAP, TDAvgSpeedDiff);

                    AM_DC_ratio_future = g_TMC_vector[tmc_index].est_DCRatio[peak_no];

                    fprintf(g_pFileTMCLink, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,",
                        FD_vc, obs_t0_in_hour, obs_t3_in_hour, obs_P_in_hour, assignment_D,
                        congestion_D, AM_DC_ratio_base, AM_DC_ratio_future, g_TMC_vector[tmc_index].est_mu[peak_no], mean_congestion_speed,
                        highest_speed, mean_speed,
                        g_TMC_vector[tmc_index].est_vt2[peak_no], t2_queue, gamma, DTASpeed, DTAP, TDAvgSpeedDiff);

                    ///////////////////////////////////////  noon
                    peak_no = 1;

                    assign_period_start_time_in_hour = 10;
                    assign_period_end_time_in_hour = 14;
                    starting_time_in_hour = 12;
                    ending_time_in_hour = 13;
                    t2 = 12;  //12 noon


                    obs_P_in_hour = g_TMC_vector[tmc_index].ScanCongestionDuration(peak_no, starting_time_in_hour,
                        ending_time_in_hour, FD_vc, obs_t0_in_hour, obs_t3_in_hour, pLink, assignment_D, assignment_VMT, assignment_VHT, assignment_VDT, assignment_VCDT,congestion_D,
                        congestion_dc_ratio, mean_congestion_mu, mean_congestion_speed, highest_speed, mean_speed, t2_speed, t2_queue, gamma);

                    V = assignment_D;
                    vcd = pLink->vc;
                    vct = pLink->vc;

                    g_TMC_vector[tmc_index].PerformEstimation(false,pLink, peak_no, assign_period_start_time_in_hour, assign_period_end_time_in_hour, t2, V, laneCapacity, vf, vcd, vct, DTASpeed, DTAP, TDAvgSpeedDiff);
                    V = assignment_D * pLink->Scenario_STA_VOC_Ratio[2];  // apply demand flow change noon
                    g_TMC_vector[tmc_index].PerformEstimation(true, pLink, peak_no, assign_period_start_time_in_hour, assign_period_end_time_in_hour, t2, V, laneCapacity, vf, vcd, vct, DTASpeed, DTAP, TDAvgSpeedDiff);

                    ///////////////////////////////////////  PM
                    peak_no = 2;

                    assign_period_start_time_in_hour = 14;
                    assign_period_end_time_in_hour = 20;
                    starting_time_in_hour = 15;
                    ending_time_in_hour = 19;
                    t2 = 17;  //5pm

                    obs_P_in_hour = g_TMC_vector[tmc_index].ScanCongestionDuration(peak_no, starting_time_in_hour,
                        ending_time_in_hour, FD_vc, obs_t0_in_hour, obs_t3_in_hour, pLink, assignment_D, assignment_VMT, assignment_VHT, assignment_VDT, assignment_VCDT,congestion_D,
                        congestion_dc_ratio, mean_congestion_mu, mean_congestion_speed, highest_speed, mean_speed, t2_speed, t2_queue, gamma);

                    V = assignment_D;
                    vcd = pLink->vc;
                    vct = pLink->vc;

                    float PM_DC_ratio_base, PM_DC_ratio_future;

                    g_TMC_vector[tmc_index].PerformEstimation(false,pLink, peak_no, assign_period_start_time_in_hour, assign_period_end_time_in_hour, t2, V, laneCapacity, vf, vcd, vct, DTASpeed, DTAP, TDAvgSpeedDiff);

                    PM_DC_ratio_base = g_TMC_vector[tmc_index].est_DCRatio[peak_no];
                    V = assignment_D * pLink->Scenario_STA_VOC_Ratio[3];  // apply demand flow change noon
                    g_TMC_vector[tmc_index].PerformEstimation(true, pLink, peak_no, assign_period_start_time_in_hour, assign_period_end_time_in_hour, t2, V, laneCapacity, vf, vcd, vct, DTASpeed, DTAP, TDAvgSpeedDiff);
                    PM_DC_ratio_future = g_TMC_vector[tmc_index].est_DCRatio[peak_no];
                    fprintf(g_pFileTMCLink, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,",
                        FD_vc, obs_t0_in_hour, obs_t3_in_hour, obs_P_in_hour, assignment_D,
                        congestion_D, PM_DC_ratio_base, PM_DC_ratio_future, g_TMC_vector[tmc_index].est_mu[peak_no], mean_congestion_speed,
                        highest_speed, mean_speed,
                        g_TMC_vector[tmc_index].est_vt2[peak_no], t2_queue, gamma, DTASpeed, DTAP, TDAvgSpeedDiff);
                }
                else
                {
                    // no data

                    fprintf(g_pFileTMCLink, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,",
                        FD_vc, obs_t0_in_hour, obs_t3_in_hour, obs_P_in_hour, assignment_D,
                        congestion_D, congestion_dc_ratio, mean_congestion_mu, mean_congestion_speed,
                        highest_speed, mean_speed,
                        t2_speed, t2_queue, gamma, DTASpeed, DTAP, TDAvgSpeedDiff);
                }
                // vc, t0, t3, P, D, mu,

                fprintf(g_pFileTMCLink, "\"%s\",", g_link_vector[i].geometry.c_str());

                fprintf(g_pFileTMCLink, "\"LINESTRING (");

                fprintf(g_pFileTMCLink, "%f %f,", g_link_vector[i].TMC_from.x, g_link_vector[i].TMC_from.y);
                fprintf(g_pFileTMCLink, "%f %f,", g_link_vector[i].TMC_to.x, g_link_vector[i].TMC_to.y);
                fprintf(g_pFileTMCLink, ")\"");
                fprintf(g_pFileTMCLink, ",,");
                fprintf(g_pFileTMCLink, "%f,%f,%f,", g_link_vector[i].Scenario_STA_VOC_Ratio[1], g_link_vector[i].Scenario_STA_VOC_Ratio[2], g_link_vector[i].Scenario_STA_VOC_Ratio[3]);
                    

                    if (assignment.m_TMClink_map.find(g_link_vector[i].TMC_code) != assignment.m_TMClink_map.end())
                    {
                        int tmc_index = assignment.m_TMClink_map[g_link_vector[i].TMC_code];

                        double PredSpeed[25], EstSpeed[25], PredSpeedDiff[25];
                        double PredTravelTime[25], EstTravelTime[25], PredTravelTimeDiff[25];
                        double MAE_total = 0;
                        double TravelTime_total = 0;
                        double MAPE_total = 0;
                        int count_total = 0;

                        for (int t = 6 * 60; t < 20 * 60; t += 60)
                        {
                            int hour = t / 60;
                            EstSpeed[hour] = g_TMC_vector[tmc_index].GetEstHourlySpeed(t);
                            fprintf(g_pFileTMCLink, "%.1f,", EstSpeed[hour]);
                        }
                        fprintf(g_pFileTMCLink, ",");
                        for (int t = 6 * 60; t < 20 * 60; t += 60)
                        {
                            int hour = t / 60;
                            PredSpeed[hour] = g_TMC_vector[tmc_index].GetPredHourlySpeed(t);
                            PredSpeedDiff[hour] = PredSpeed[hour] - EstSpeed[hour];
                            MAE_total += PredSpeedDiff[hour];
                            MAPE_total += fabs(PredSpeedDiff[hour]) / max(1, PredSpeed[hour]);

                            count_total += 1;

                            fprintf(g_pFileTMCLink, "%.1f,", PredSpeed[hour]);
                        }
                        fprintf(g_pFileTMCLink, ",");
                        for (int t = 6 * 60; t < 20 * 60; t += 60)
                        {
                            int hour = t / 60;
                            fprintf(g_pFileTMCLink, "%.1f,", PredSpeedDiff[hour]);
                        }
                        fprintf(g_pFileTMCLink, ",");

                        fprintf(g_pFileTMCLink, "%.2f,%.2f,", MAE_total / max(1, count_total), MAPE_total / max(1, count_total) * 100);
                        fprintf(g_pFileTMCLink, ",");

                        // travel time, to be computed for corridor level travel time
                        for (int t = 6 * 60; t < 20 * 60; t += 60)
                        {
                            int hour = t / 60;
                            EstSpeed[hour] = g_TMC_vector[tmc_index].GetEstHourlySpeed(t);
                            EstTravelTime[hour] = g_link_vector[i].length / max(1, EstSpeed[hour]) * 60;

                            fprintf(g_pFileTMCLink, "%.2f,", EstTravelTime[hour]);
                        }
                        fprintf(g_pFileTMCLink, ",");
                        for (int t = 6 * 60; t < 20 * 60; t += 60)
                        {
                            int hour = t / 60;
                            PredSpeed[hour] = g_TMC_vector[tmc_index].GetPredHourlySpeed(t);
                            PredTravelTime[hour] = g_link_vector[i].length / max(1, PredSpeed[hour]) * 60;
                            count_total += 1;

                            fprintf(g_pFileTMCLink, "%.2f,", PredTravelTime[hour]);
                        }

                        fprintf(g_pFileTMCLink, ",");

                        // 
                        for (int t = 6 * 60; t < 20 * 60; t += 5)
                        {
                            fprintf(g_pFileTMCLink, "%.1f,", g_link_vector[i].length / max(1, g_TMC_vector[tmc_index].pred_speed[t]) * 60);
                        }
                }
                fprintf(g_pFileTMCLink, "\n");

            }


        }

        fclose(g_pFileTMCLink);
    }
    else
    {
        dtalog.output() << "Error: File TMC_link_scenario.csv cannot be opened.\n It might be currently used and locked by EXCEL." << endl;
        g_ProgramStop();

    }
    
}
void g_OutputModelFiles()
{
    FILE* g_pFileModelNode = fopen("model_node.csv", "w");

    if (g_pFileModelNode != NULL)
    {
        fprintf(g_pFileModelNode, "node_id,node_no,activity_node_flag,zone_id,ctrl_type,production,attraction,x_coord,y_coord\n");
        for (int i = 0; i < g_node_vector.size(); i++)
        {
            if(g_node_vector[i].node_id >=0)
            {
                fprintf(g_pFileModelNode, "%d,%d,%d,%d,%d,%f,%f,%f,%f\n",
                g_node_vector[i].node_id,
                g_node_vector[i].node_seq_no,
                g_node_vector[i].is_activity_node,
                g_node_vector[i].zone_org_id,
                g_node_vector[i].ctrl_type,
                g_node_vector[i].production,
                g_node_vector[i].attraction,
                g_node_vector[i].x,
                g_node_vector[i].y
            );

            }

        }

        fclose(g_pFileModelNode);
    }
    else
    {
    dtalog.output() << "Error: File model_node.csv cannot be opened.\n It might be currently used and locked by EXCEL." << endl;
    g_ProgramStop();


    }

    FILE* g_pFileModelLink = fopen("model_link.csv", "w");

    if (g_pFileModelLink != NULL)
    {
        fprintf(g_pFileModelLink, "link_id,link_no,from_node_id,to_node_id,link_type,link_type_name,FT,AT,nlanes,free_speed,capacity,kc,geometry\n");

        //VDF_fftt1,VDF_cap1,VDF_alpha1,VDF_beta1
        for (int i = 0; i < g_link_vector.size(); i++)
        {
            if(g_link_vector[i].link_type >=0)
            {
                fprintf(g_pFileModelLink, "%s,%d,%d,%d,%d,%s,%d,%d,%d,%f,%f,%f,",
                    g_link_vector[i].link_id.c_str(),
                    g_link_vector[i].link_seq_no,
                    g_node_vector[g_link_vector[i].from_node_seq_no].node_id,
                    g_node_vector[g_link_vector[i].to_node_seq_no].node_id,
                    g_link_vector[i].link_type,
                    g_link_vector[i].link_type_name.c_str(),
                    g_link_vector[i].FT,
                    g_link_vector[i].AT,
                    g_link_vector[i].number_of_lanes,
                    g_link_vector[i].free_speed,
                    g_link_vector[i].lane_capacity,
                    g_link_vector[i].kc
                    //g_link_vector[i].VDF_period[0].FFTT,
                    //g_link_vector[i].VDF_period[0].capacity,
                    //g_link_vector[i].VDF_period[0].alpha,
                    //g_link_vector[i].VDF_period[0].beta,
                );

                fprintf(g_pFileModelLink, "\"%s\",\n",  g_link_vector[i].geometry.c_str());
            }
            
              
        }

        fclose(g_pFileModelLink);
    }
    else
    {
        dtalog.output() << "Error: File model_link.csv cannot be opened.\n It might be currently used and locked by EXCEL." << endl;
        g_ProgramStop();

    }



}

double g_CalculateP2PDistanceInLonglatFromLatitudeLongitude(double p1x, double p1y, double p2x, double p2y)
{
    double distance = sqrt((p2y - p1y) * (p2y - p1y) + (p2x - p1x) * (p2x - p1x));
    return distance;
}


double g_CalculateP2PDistanceInMeterFromLatitudeLongitude(double p1x, double p1y, double p2x, double p2y)
{
    double PI = 3.1415926;
    double Equatorial_Radius = 3963.19059*1609; // unit: mile-> meter
    double toradians = 3.1415926 / 180.0;
    double todeg = 180.0 / PI;

    double p2lat = p2x * toradians;
    double p2lng = p2y * toradians;

    double p1lat = p1x * toradians;
    double p1lng = p1y * toradians;

    double distance = acos(sin(p1lat) * sin(p2lat) + cos(p1lat) * cos(p2lat) * cos(p2lng - p1lng)) * Equatorial_Radius;  // unit: mile
    return distance;
}

void g_TripGeneration(Assignment& assignment)
{
    // accessibility 
    for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); ++tau)
    {
        for (int at = 0; at < assignment.g_AgentTypeVector.size(); ++at)
        {

           for (int orig = 0; orig < g_zone_vector.size(); ++orig)  // o
            {
                    for (int dest = 0; dest < g_zone_vector.size(); ++dest)  // d
                    {
                        if(orig!=dest)
                        {
                        float distance_in_meter = g_CalculateP2PDistanceInMeterFromLatitudeLongitude(g_zone_vector[orig].cell_x, g_zone_vector[orig].cell_y, g_zone_vector[dest].cell_x, g_zone_vector[dest].cell_y);
                        g_zone_vector[orig].m_ODAccessibilityMatrix[at][tau].distance_map[dest] = distance_in_meter/1000;
                        float travel_time_in_min = distance_in_meter/1000 / assignment.g_AgentTypeVector[at].avg_speed * 60;

                         g_zone_vector[orig].m_ODAccessibilityMatrix[at][tau].value_map[dest] = travel_time_in_min;
                        }
                        
                    }
                
            }
        }
    }

    int out_of_bound_log_count = 0;
    int trip_accessibility_log_count = 0;
    int trip_distribution_log_count = 0;

    for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); ++tau)
    {
        for (int at = 0; at < assignment.g_AgentTypeVector.size(); ++at)
        {
            for (int orig = 0; orig < g_zone_vector.size(); ++orig)  // o
                                                                     // gravity model;
            {
                if (g_zone_vector[orig].gravity_production[at] > 0.00001)
                {

                    float total_attraction_utility = 0;
                    int count = 0;

                    for (int d = 0; d < g_zone_vector.size(); ++d)
                    {
                        if (orig != d)
                        {

                            g_zone_vector[d].gravity_est_attraction[at] = 0;

                            float cut_off = assignment.g_AgentTypeVector[at].trip_time_budget_in_min ;
                            if (g_zone_vector[d].gravity_attraction[at] > 0)
                            {
                                //double disutility = g_zone_vector[orig].m_ODAccessibilityMatrix[at][tau].value_map[d] * beta;
                                 double exp_disutility = g_zone_vector[d].gravity_attraction[at];
                                g_zone_vector[orig].m_ODAccessibilityMatrix[at][tau].disutility_map[d] = exp_disutility;

                                if (g_zone_vector[orig].m_ODAccessibilityMatrix[at][tau].value_map[d] < cut_off)
                                {
                                    if(trip_accessibility_log_count <=100)
                                    {
                                    dtalog.output() << "agent type: " << assignment.g_AgentTypeVector[at].agent_type.c_str() << ", o: " << orig << ",d:" << d <<
                                        ", gc distance = " << g_zone_vector[orig].m_ODAccessibilityMatrix[at][tau].distance_map[d] <<
                                        ", travel time =" << g_zone_vector[orig].m_ODAccessibilityMatrix[at][tau].value_map[d] <<
                                        ",value = " << exp_disutility << endl;
                                    }
                                    total_attraction_utility += exp_disutility;
                                    trip_accessibility_log_count++;
                                    count++;
                                }
                                else
                                {
                                    if(out_of_bound_log_count < 10)
                                    {
                                    dtalog.output() << "out of bound: agent type: " << assignment.g_AgentTypeVector[at].agent_type.c_str() << ",o: " << orig << ",d:" << d <<
                                        ", gc distance = " << g_zone_vector[orig].m_ODAccessibilityMatrix[at][tau].distance_map[d] <<
                                        ", travel time =" << g_zone_vector[orig].m_ODAccessibilityMatrix[at][tau].value_map[d] <<
                                        ",value = " << exp_disutility << endl;
                                    }

                                    out_of_bound_log_count++;

                                }
                            }
                        }
                    }

                    dtalog.output() << "o: " << orig << ", total_attraction_utility =" << total_attraction_utility << endl;

                    if(count>0)
                    {
                        for (int dest = 0; dest < g_zone_vector.size(); ++dest)  // d
                        {
                            if (orig != dest )
                            {
                                if (g_zone_vector[dest].gravity_attraction[at] > 0)
                                {
                                    float cut_off = assignment.g_AgentTypeVector[at].trip_time_budget_in_min;
                                    if (g_zone_vector[orig].m_ODAccessibilityMatrix[at][tau].value_map[dest] < cut_off)
                                    {

                                   double ratio = g_zone_vector[orig].m_ODAccessibilityMatrix[at][tau].disutility_map[dest] / total_attraction_utility;
                                

                                   g_zone_vector[orig].m_ODMatrix[at][tau].value_map[dest] = g_zone_vector[orig].gravity_production[at]*ratio;

                                   if (trip_distribution_log_count < 100)
                                   {
                                       dtalog.output() << "agent type: " << assignment.g_AgentTypeVector[at].agent_type.c_str() << ", o: " << orig << ",d:" << dest << ", ratio =" << ratio <<
                                           ",trip = " << g_zone_vector[orig].m_ODMatrix[at][tau].value_map[dest] << endl;
                                   }
                                   trip_distribution_log_count++;

                                g_zone_vector[dest].gravity_est_attraction[at] += g_zone_vector[orig].m_ODMatrix[at][tau].value_map[dest];
                                    }
                                }

                            }
                        }
                    }
                }

            }
        }
    }
}
double g_RandomGenerateActivityNodes(Assignment& assignment)
{

    int activity_node_count = 0;
    for (int i = 0; i < g_node_vector.size(); i++)
    {

        if (g_node_vector[i].is_activity_node >= 1)
        {
            activity_node_count++;
        }
    }


    if (activity_node_count <= 1)
    {
        activity_node_count = 0;
        int sampling_rate = 10;

        for (int i = 0; i < g_node_vector.size(); i++)
        {

            if (i % sampling_rate==0)
            {
                g_node_vector[i].is_activity_node = 10;//random generation
                activity_node_count++;
            }
        }

        if (activity_node_count <= 1)
        {
            activity_node_count = 0;
            sampling_rate = 2;

            for (int i = 0; i < g_node_vector.size(); i++)
            {

                if (i % sampling_rate==0)
                {
                    g_node_vector[i].is_activity_node = 10;//random generation
                    activity_node_count++;
                }
            }
            // still no activity nodes, define all nodes as activity nodes
            if (activity_node_count <= 1)
            {
                activity_node_count = 0;

                for (int i = 0; i < g_node_vector.size(); i++)
                {

                    g_node_vector[i].is_activity_node = 10;//random generation
                    activity_node_count++;
                }
            }
        }


    }


// calculate avg near by distance; 
    double total_near_by_distance = 0;
    activity_node_count = 0; 
    for (int i = 0; i < g_node_vector.size(); i++)
    {
        double min_near_by_distance = 100;
        if (g_node_vector[i].is_activity_node)
        {
            activity_node_count ++;
            for (int j = 0; j < g_node_vector.size(); j++)
            {
                if (i != j && g_node_vector[j].is_activity_node)
                {
                    double near_by_distance = g_CalculateP2PDistanceInLonglatFromLatitudeLongitude(g_node_vector[i].x, g_node_vector[i].y, g_node_vector[j].x, g_node_vector[j].y);

                    if (near_by_distance < min_near_by_distance)
                        min_near_by_distance = near_by_distance;

                }

            }

            total_near_by_distance += min_near_by_distance;
            activity_node_count++;
        }
    }

    double nearby_distance = total_near_by_distance/max(1, activity_node_count);
    return nearby_distance;

}

void g_GridZoneGeneration(Assignment& assignment)
{
    dtalog.output() << "Step 1.4.1: QEM mode for creating node 2 zone mapping" << endl;

    double activity_nearbydistance = g_RandomGenerateActivityNodes(assignment);
    // initialization of grid rectangle boundary
    double left = 100000000;
    double right = -100000000;
    double top = -1000000000;
    double  bottom = 1000000000;

    for (int i = 0; i < g_node_vector.size(); i++)
    {
        // exapnd the grid boundary according to the nodes
        left = min(left, g_node_vector[i].x);
        right = max(right, g_node_vector[i].x);
        top = max(top, g_node_vector[i].y);
        bottom = min(bottom, g_node_vector[i].y);

    }

    int grid_size = 8;

    if (g_node_vector.size() > 3000)
        grid_size = 10;
    if (g_node_vector.size() > 10000)
        grid_size = 20;
    if (g_node_vector.size() > 40000)
        grid_size = 30;

    double temp_resolution = (((right - left) / grid_size + (top - bottom) / grid_size)) / 2.0;

    if (activity_nearbydistance*4 < temp_resolution)
    {
        temp_resolution = activity_nearbydistance * 4;
       
    }


    vector<double> ResolutionVector;

    ResolutionVector.push_back(0.00005);
    ResolutionVector.push_back(0.0001);
    ResolutionVector.push_back(0.0002);
    ResolutionVector.push_back(0.00025);
    ResolutionVector.push_back(0.0005);
    ResolutionVector.push_back(0.00075);
    ResolutionVector.push_back(0.001);
    ResolutionVector.push_back(0.002);
    ResolutionVector.push_back(0.0025);
    ResolutionVector.push_back(0.005);
    ResolutionVector.push_back(0.0075);
    ResolutionVector.push_back(0.01);
    ResolutionVector.push_back(0.02);
    ResolutionVector.push_back(0.025);
    ResolutionVector.push_back(0.05);
    ResolutionVector.push_back(0.075);
    ResolutionVector.push_back(0.1);
    ResolutionVector.push_back(0.2);
    ResolutionVector.push_back(0.25);
    ResolutionVector.push_back(0.5);
    ResolutionVector.push_back(0.75);
    ResolutionVector.push_back(1);
    ResolutionVector.push_back(2);
    ResolutionVector.push_back(2.5);
    ResolutionVector.push_back(5);
    ResolutionVector.push_back(7.5); 
    ResolutionVector.push_back(10);
    ResolutionVector.push_back(20);
    ResolutionVector.push_back(25);
    ResolutionVector.push_back(50);
    ResolutionVector.push_back(75);

    double ClosestResolution = 1;

    if (temp_resolution < ResolutionVector[0])
        temp_resolution = ResolutionVector[0];

    for (unsigned int i = 0; i < ResolutionVector.size() - 1; i++)
    {
        if ((temp_resolution > ResolutionVector[i] + 0.000001) && temp_resolution < ResolutionVector[i + 1])
        {
            temp_resolution = ResolutionVector[i + 1]; // round up
            break;

        }
    }

    assignment.m_GridResolution = temp_resolution;

    assignment.zone_id_2_node_no_mapping.clear();
    dtalog.output() << "Step 1.4.2: Grid Resolution " << assignment.m_GridResolution << endl;

    int activity_node_count = 0;
    for (int i = 0; i < g_node_vector.size(); i++)
    {

        if (g_node_vector[i].is_activity_node >= 1)
        {

            if (g_node_vector[i].node_id == 966)
            {
                int itest = 1;
            }
            __int64 cell_id = g_GetCellID(g_node_vector[i].x, g_node_vector[i].y, assignment.m_GridResolution);
            int zone_id;

            if (assignment.cell_id_mapping.find(cell_id) == assignment.cell_id_mapping.end())  // create a cell
            {
                //create zone
                assignment.cell_id_mapping[cell_id] = g_node_vector[i].node_id;


                dtalog.output() << "Step 1.2: creating cell " << cell_id << " using node id " << g_node_vector[i].node_id << endl;

                zone_id = assignment.cell_id_mapping[cell_id]; // which is the node id when a cell is created. 
                if (assignment.zone_id_2_node_no_mapping.find(zone_id) == assignment.zone_id_2_node_no_mapping.end()) // create a zone 
                {
                    dtalog.output() << "Step 1.2: creating zone " << zone_id << " using node id " << g_node_vector[i].node_id << endl;
                    //create zone
                    assignment.zone_id_2_node_no_mapping[zone_id] = i;
                    assignment.zone_id_2_cell_id_mapping[zone_id] = cell_id;
                    g_node_vector[i].zone_org_id = zone_id;
                }
            }
            else
            {
                zone_id = assignment.cell_id_mapping[cell_id]; // which is the node id when a cell is created. 
                // for physcial nodes because only centriod can have valid zone_id.
                g_node_vector[i].zone_org_id = zone_id;

            }

            activity_node_count++;


        }
    }
    
    dtalog.output() << "Step 1.4.3: creating " << assignment.zone_id_2_node_no_mapping.size() << " zones." << " # of activity nodes =" << activity_node_count << endl;

}

void g_ReadInputData(Assignment& assignment)
{

    g_ReadSubarea(assignment);
    assignment.g_LoadingStartTimeInMin = 99999;
    assignment.g_LoadingEndTimeInMin = 0;

    //step 0:read demand period file
    CCSVParser parser_demand_period;
	dtalog.output() << "_____________" << endl;
	dtalog.output() << "Step 1: Reading input data" << endl;
	dtalog.output() << "_____________" << endl;

    dtalog.output() << "Step 1.1: Reading section [demand_period] in setting.csv..." << endl;

    parser_demand_period.IsFirstLineHeader = false;
    if (parser_demand_period.OpenCSVFile("settings.csv", false))
    {
        while (parser_demand_period.ReadRecord_Section())
        {
            if (parser_demand_period.SectionName == "[demand_period]")
            {
                CDemand_Period demand_period;

                if (!parser_demand_period.GetValueByFieldName("demand_period_id", demand_period.demand_period_id))
                    break;

                if (!parser_demand_period.GetValueByFieldName("demand_period", demand_period.demand_period))
                {
                    dtalog.output() << "Error: Field demand_period in file demand_period cannot be read." << endl;
                    g_ProgramStop();
                }

                vector<float> global_minute_vector;

                if (!parser_demand_period.GetValueByFieldName("time_period", demand_period.time_period))
                {
                    dtalog.output() << "Error: Field time_period in file demand_period cannot be read." << endl;
                    g_ProgramStop();
                }

                //input_string includes the start and end time of a time period with hhmm format
                global_minute_vector = g_time_parser(demand_period.time_period); //global_minute_vector incldue the starting and ending time
                if (global_minute_vector.size() == 2)
                {
                    demand_period.starting_time_slot_no = global_minute_vector[0] / MIN_PER_TIMESLOT;
                    demand_period.ending_time_slot_no = global_minute_vector[1] / MIN_PER_TIMESLOT;

                    if (global_minute_vector[0] < assignment.g_LoadingStartTimeInMin)
                        assignment.g_LoadingStartTimeInMin = global_minute_vector[0];

                    if (global_minute_vector[1] > assignment.g_LoadingEndTimeInMin)
                        assignment.g_LoadingEndTimeInMin = global_minute_vector[1];

                    //g_fout << global_minute_vector[0] << endl;
                    //g_fout << global_minute_vector[1] << endl;
                }

                assignment.demand_period_to_seqno_mapping[demand_period.demand_period] = assignment.g_DemandPeriodVector.size();
                assignment.g_DemandPeriodVector.push_back(demand_period);
            }
        }

        parser_demand_period.CloseCSVFile();

        if(assignment.g_DemandPeriodVector.size() == 0)
        {
            dtalog.output() << "Error:  Section demand_period has no information." << endl;
            g_ProgramStop();
        }
    }
    else
    {
        dtalog.output() << "Error: File settings.csv cannot be opened.\n Continue to use default values." << endl;

        CDemand_Period demand_period;
        demand_period.demand_period_id = 1;
        demand_period.demand_period = "am";
        assignment.demand_period_to_seqno_mapping[demand_period.demand_period] = 0;
        demand_period.starting_time_slot_no = 7*60/ MIN_PER_TIMESLOT;
        demand_period.ending_time_slot_no = 8*60 / MIN_PER_TIMESLOT;

        assignment.g_LoadingStartTimeInMin = 7 * 60;
        assignment.g_LoadingEndTimeInMin = 8*60;

        //g_fout << global_minute_vector[0] << endl;
        //g_fout << global_minute_vector[1] << endl;

        assignment.demand_period_to_seqno_mapping[demand_period.demand_period] = assignment.g_DemandPeriodVector.size();
        assignment.g_DemandPeriodVector.push_back(demand_period);

    }

 
    dtalog.output() << "number of demand periods = " << assignment.g_DemandPeriodVector.size() << endl;

    assignment.g_number_of_demand_periods = assignment.g_DemandPeriodVector.size();
    //step 1:read demand type file

    dtalog.output() << "Step 1.2: Reading section [link_type] in setting.csv..." << endl;

    CCSVParser parser_link_type;
    parser_link_type.IsFirstLineHeader = false;
    if (parser_link_type.OpenCSVFile("settings.csv", false))
    {
        // create a special link type as virtual connector
        CLinkType element_vc;
        // -1 is for virutal connector
        element_vc.link_type = -1;
        element_vc.type_code = "c";
        element_vc.traffic_flow_code = 0;
        assignment.g_LinkTypeMap[element_vc.link_type] = element_vc;
        //end of create special link type for virtual connectors

        int line_no = 0;

        while (parser_link_type.ReadRecord_Section())
        {
            if (parser_link_type.SectionName == "[link_type]")
            {
                CLinkType element;

                if (!parser_link_type.GetValueByFieldName("link_type", element.link_type))
                {
                    if (line_no == 0)
                    {
                        dtalog.output() << "Error: Field link_type cannot be found in file link_type.csv." << endl;
                        g_ProgramStop();
                    }
                    else
                    {
                        // read empty line
                        break;
                    }
                }

                if (assignment.g_LinkTypeMap.find(element.link_type) != assignment.g_LinkTypeMap.end())
                {
                    dtalog.output() << "Error: Field link_type " << element.link_type << " has been defined more than once in file link_type.csv." << endl;
                    g_ProgramStop();
                }

                string traffic_flow_code_str;
                parser_link_type.GetValueByFieldName("type_code", element.type_code, true);
                parser_link_type.GetValueByFieldName("traffic_flow_code", traffic_flow_code_str);

                // by default bpr
                element.traffic_flow_code = 0;

                if (traffic_flow_code_str == "point_queue")
                    element.traffic_flow_code = 1;

                if (traffic_flow_code_str == "spatial_queue")
                    element.traffic_flow_code = 2;

                if (traffic_flow_code_str == "kw")
                    element.traffic_flow_code = 3;

                dtalog.output() << "important: traffic_flow_code on link type " << element.link_type  << " is " << element.traffic_flow_code  << endl;


                assignment.g_LinkTypeMap[element.link_type] = element;
                line_no++;
            }
        }

        parser_link_type.CloseCSVFile();
    }
    else
    {
        CLinkType element_vc;
        // -1 is for virutal connector
        element_vc.link_type = 1;
        element_vc.type_code = "f";
        element_vc.traffic_flow_code = 0;
        assignment.g_LinkTypeMap[element_vc.link_type] = element_vc;

        // -1 is for virutal connector
        element_vc.link_type = 2;
        element_vc.type_code = "f";
        element_vc.traffic_flow_code = 0;
        assignment.g_LinkTypeMap[element_vc.link_type] = element_vc;

        element_vc.link_type = 3;
        element_vc.type_code = "f";
        element_vc.traffic_flow_code = 0;
        assignment.g_LinkTypeMap[element_vc.link_type] = element_vc;

        element_vc.link_type = 4;
        element_vc.type_code = "f";
        element_vc.traffic_flow_code = 0;
        assignment.g_LinkTypeMap[element_vc.link_type] = element_vc;


        // -1 is for virutal connector
        element_vc.link_type = -1;
        element_vc.type_code = "c";
        element_vc.traffic_flow_code = 0;
        assignment.g_LinkTypeMap[element_vc.link_type] = element_vc;
        //end of create special link type for virtual connectors

    }

    dtalog.output() << "number of link types = " << assignment.g_LinkTypeMap.size() << endl;

    CCSVParser parser_agent_type;
    dtalog.output() << "Step 1.3: Reading section [agent_type] in setting.csv..." << endl;

    parser_agent_type.IsFirstLineHeader = false;
    if (parser_agent_type.OpenCSVFile("settings.csv", false))
    {
        assignment.g_AgentTypeVector.clear();
        while (parser_agent_type.ReadRecord_Section())
        {
            if(parser_agent_type.SectionName == "[agent_type]")
            {
                CAgent_type agent_type;
                agent_type.agent_type_no = assignment.g_AgentTypeVector.size();

                if (!parser_agent_type.GetValueByFieldName("agent_type", agent_type.agent_type))
                    break;

                parser_agent_type.GetValueByFieldName("VOT", agent_type.value_of_time, false, false);

                // scan through the map with different node sum for different paths
                parser_agent_type.GetValueByFieldName("PCE", agent_type.PCE, false, false);
                assignment.agent_type_2_seqno_mapping[agent_type.agent_type] = assignment.g_AgentTypeVector.size();

                assignment.g_AgentTypeVector.push_back(agent_type);
                assignment.g_number_of_agent_types = assignment.g_AgentTypeVector.size();
            }
        }
        parser_agent_type.CloseCSVFile();

        if (assignment.g_AgentTypeVector.size() == 0 )
            dtalog.output() << "Error: Section agent_type does not contain information." << endl;
    }
    else
    {
        assignment.AddAgentType("auto", 10, 1, 40, 25,0.6);
        assignment.AddAgentType("walk", 10, 0.5, 2, 10,0.1);
        assignment.AddAgentType("bike", 10, 0.5, 5,15,0.1);
        assignment.AddAgentType("truck", 10, 0.5, 40, 100,0.1);
        assignment.AddAgentType("cav", 10, 0.5, 40, 30,0.1);
    }

    if (assignment.g_AgentTypeVector.size() >= _MAX_AGNETTYPES)
    {
        dtalog.output() << "Error: agent_type = " << assignment.g_AgentTypeVector.size() << " in section agent_type is too large. " << "_MAX_AGNETTYPES = " << _MAX_AGNETTYPES << "Please contact program developers!";
        g_ProgramStop();
    }

    dtalog.output() << "number of agent typess = " << assignment.g_AgentTypeVector.size() << endl;

    assignment.g_number_of_nodes = 0;
    assignment.g_number_of_links = 0;  // initialize  the counter to 0

    int internal_node_seq_no = 0;
    // step 3: read node file


    std::map<int, int> zone_id_production;
    std::map<int, int> zone_id_attraction;

    CCSVParser parser;


    dtalog.output() << "Step 1.4: Reading node data in node.csv..."<< endl;

    if (parser.OpenCSVFile("node.csv", true))
    {
        while (parser.ReadRecord())  // if this line contains [] mark, then we will also read field headers.
        {
            int node_id;
            if (!parser.GetValueByFieldName("node_id", node_id))
                continue;

            if (node_id == 202904)
            {
                int iBebug = 1;
            }

            if (assignment.g_subarea_mode == 1)  // windowing approach == 1, skip outside nodes
            {
                if (assignment.m_subarea_node_id_map.find(node_id) == assignment.m_subarea_node_id_map.end())
                    continue; // skip this non-subarea nodes
            }

            if (assignment.g_node_id_to_seq_no_map.find(node_id) != assignment.g_node_id_to_seq_no_map.end())
            {
                //has been defined
                continue;
            }

            assignment.g_node_id_to_seq_no_map[node_id] = internal_node_seq_no;

            // create a node object
            CNode node;
            node.node_id = node_id;
            node.node_seq_no = internal_node_seq_no;

            int zone_id = -1;
            parser.GetValueByFieldName("zone_id", zone_id);
            if (zone_id >= 1)
            {
                node.is_activity_node = 1;  // from zone
            }

            parser.GetValueByFieldName("x_coord", node.x,true, false);
            parser.GetValueByFieldName("y_coord", node.y,true, false);
            parser.GetValueByFieldName("ctrl_type", node.ctrl_type,false);
            int POI_id = 0;
            parser.GetValueByFieldName("POI_id", POI_id, false, false);

            if (POI_id >= 1)
            {
                node.is_activity_node = 3;  // from POI
            }

            string activity_type;
            parser.GetValueByFieldName("activity_type", activity_type,false);
            {
                //if(activity_type=="residential")
                //{
                //node.is_activity_node = 4;
                //}
            }


            if (assignment.g_subarea_mode >= 1)
            {
                int is_boundary = 0;
                parser.GetValueByFieldName("is_boundary", is_boundary, false, false);

                if (is_boundary >= 1)
                {
                    node.is_activity_node = 2;  // from boundary

                    if (zone_id <= 1)
                    {
                        zone_id = node_id + 100000;
                    }
                }

                if (assignment.m_subarea_boundary_node_map.find(node_id) != assignment.m_subarea_boundary_node_map.end())
                {
                    node.is_activity_node = 2; // dynamically define is_boundary flag as 1 based on subarea
                    zone_id = node_id+100000;
                }

            }

            // this is an activity node // we do not allow zone id of zero
            if(zone_id>=1)
            {
                // for physcial nodes because only centriod can have valid zone_id.
                node.zone_org_id = zone_id;
                if (assignment.zone_id_2_node_no_mapping.find(zone_id) == assignment.zone_id_2_node_no_mapping.end())
                {
                    //create zone
                    assignment.zone_id_2_node_no_mapping[zone_id] = internal_node_seq_no;
                }

                // for od calibration, I think we don't need to implement for now
                if (assignment.assignment_mode == 3 || assignment.assignment_mode == 10)
                {
                    parser.GetValueByFieldName("production", node.production,false);
                    parser.GetValueByFieldName("attraction", node.attraction, false);

                }
            }

            /*node.x = x;
            node.y = y;*/
            internal_node_seq_no++;

            // push it to the global node vector
            g_node_vector.push_back(node);
            assignment.g_number_of_nodes++;

            if (assignment.g_number_of_nodes % 5000 == 0)
                dtalog.output() << "reading " << assignment.g_number_of_nodes << " nodes.. " << endl;
        }

        dtalog.output() << "number of nodes = " << assignment.g_number_of_nodes << endl;

    	// fprintf(g_pFileOutputLog, "number of nodes =,%d\n", assignment.g_number_of_nodes);
        parser.CloseCSVFile();
    }

    /// creat cell struture

    //if (assignment.assignment_mode == 10)  // quick estimation mode
    //{        // initialize zone vector
    //    g_GridZoneGeneration(assignment);

    //}
        // special handling if # of zones is less than the number of memory blocks
        if (assignment.zone_id_2_node_no_mapping.size() < assignment.g_number_of_memory_blocks)
            assignment.g_number_of_memory_blocks = max(1, assignment.zone_id_2_node_no_mapping.size() / 2); // reduce the memeory size requirements


        std::map<int, int> waring_message_link_type_map;
        // initialize zone vector
        dtalog.output() << "Step 1.5: Initializing O-D zone vector..." << endl;

        std::map<int, int>::iterator it;

        for (it = assignment.zone_id_2_node_no_mapping.begin(); it != assignment.zone_id_2_node_no_mapping.end(); ++it)
        {
            COZone ozone;

            if (it->first == 966)
            {
                int itest = 1;
            }
            // for each zone, we have to also create centriod
            ozone.zone_id = it->first;  // zone_id
            ozone.cell_id = assignment.zone_id_2_cell_id_mapping[it->first];
            ozone.zone_seq_no = g_zone_vector.size();
            ozone.cell_x = g_node_vector[it->second].x;
            ozone.cell_y = g_node_vector[it->second].y;

            dtalog.output() << "create zone id = " << ozone.zone_id << " with representive node id " << it->second << ",x = " << g_node_vector[it->second].x << ",y=" <<
                ozone.cell_y << endl;


            assignment.g_zoneid_to_zone_seq_no_mapping[ozone.zone_id] = ozone.zone_seq_no;  // create the zone id to zone seq no mapping

            // create a centriod
            CNode node;
            // very large number as a special id
            node.node_id = -1 * (ozone.zone_id) - 1000000;
            node.node_seq_no = g_node_vector.size();
            assignment.g_node_id_to_seq_no_map[node.node_id] = node.node_seq_no;
            node.zone_id = ozone.zone_id;
            // push it to the global node vector
            g_node_vector.push_back(node);
            assignment.g_number_of_nodes++;

            ozone.node_seq_no = node.node_seq_no;
            // this should be the only one place that defines this mapping
            assignment.zone_id_to_centriod_node_no_mapping[ozone.zone_id] = node.node_seq_no;
            // add element into vector
            g_zone_vector.push_back(ozone);
        }


        dtalog.output() << "number of zones = " << g_zone_vector.size() << endl;
        // step 4: read link file

        CCSVParser parser_link;

        double lane_capacity_vector[20];
        double free_speed_vector[20];
        double nlanes_vector[20];


        for (int link_type_i = 0; link_type_i < 20; link_type_i++)
        {
            lane_capacity_vector[link_type_i] = 1800;
            free_speed_vector[link_type_i] = 70;
            nlanes_vector[link_type_i] = 2;
        }

        lane_capacity_vector[1] = 1700;  //motorway	
        lane_capacity_vector[2] = 1700;  //truck
        lane_capacity_vector[3] = 900;  //primary
        lane_capacity_vector[4] = 700;  //secondary
        lane_capacity_vector[5] = 500;  //tertiary
        lane_capacity_vector[6] = 550;  //residential
        lane_capacity_vector[11] = 500;  //unclassified

        free_speed_vector[1] = 70;  //motorway	
        free_speed_vector[2] = 70;  //truck
        free_speed_vector[3] = 40;  //primary
        free_speed_vector[4] = 35;  //secondary
        free_speed_vector[5] = 35;  //tertiary
        free_speed_vector[6] = 35;  //residential
        free_speed_vector[11] = 35;  //unclassified

        nlanes_vector[1] = 3;  //motorway	
        nlanes_vector[2] = 2;  //truck
        nlanes_vector[3] = 3;  //primary
        nlanes_vector[4] = 2;  //secondary
        nlanes_vector[5] = 2;  //tertiary
        nlanes_vector[6] = 1;  //residential
        nlanes_vector[11] = 1;  //unclassified

    
        dtalog.output() << "Step 1.6: Reading link data in link.csv... " << endl;
        if (parser_link.OpenCSVFile("link.csv", true))
        {
            while (parser_link.ReadRecord())  // if this line contains [] mark, then we will also read field headers.
            {
                int from_node_id;
                if (!parser_link.GetValueByFieldName("from_node_id", from_node_id))
                    continue;

                int to_node_id;
                if (!parser_link.GetValueByFieldName("to_node_id", to_node_id))
                    continue;

                string linkID;
                parser_link.GetValueByFieldName("link_id", linkID);
                // add the to node id into the outbound (adjacent) node list

                if (assignment.g_subarea_mode >= 1)  // subarea handling 
                {
                    if (assignment.g_node_id_to_seq_no_map.find(from_node_id) == assignment.g_node_id_to_seq_no_map.end())
                    {
                        continue; //has not been defined
                    }

                    if (assignment.g_node_id_to_seq_no_map.find(to_node_id) == assignment.g_node_id_to_seq_no_map.end())
                    {
                        continue; //has not been defined
                    }

                    if (assignment.m_subarea_node_id_map.find(from_node_id) == assignment.m_subarea_node_id_map.end())
                    {
                        continue; //has not been defined
                    }

                    if (assignment.m_subarea_node_id_map.find(to_node_id) == assignment.m_subarea_node_id_map.end())
                    {
                        continue; //has not been defined
                    }
                }

                if (assignment.g_node_id_to_seq_no_map.find(from_node_id) == assignment.g_node_id_to_seq_no_map.end())
                {
                    dtalog.output() << "Error: from_node_id " << from_node_id << " in file link.csv is not defined in node.csv." << endl;
                    continue; //has not been defined
                }

                if (assignment.g_node_id_to_seq_no_map.find(to_node_id) == assignment.g_node_id_to_seq_no_map.end())
                {
                    dtalog.output() << "Error: to_node_id " << to_node_id << " in file link.csv is not defined in node.csv." << endl;
                    continue; //has not been defined
                }


                if (linkID.size() > 0 && assignment.g_link_id_map.find(linkID) != assignment.g_link_id_map.end())
                    dtalog.output() << "Error: link_id " << linkID.c_str() << " has been defined more than once. Please check link.csv." << endl;

                int internal_from_node_seq_no = assignment.g_node_id_to_seq_no_map[from_node_id];  // map external node number to internal node seq no.
                int internal_to_node_seq_no = assignment.g_node_id_to_seq_no_map[to_node_id];

                // create a link object
                CLink link;

                link.from_node_seq_no = internal_from_node_seq_no;
                link.to_node_seq_no = internal_to_node_seq_no;
                link.link_seq_no = assignment.g_number_of_links;
                link.to_node_seq_no = internal_to_node_seq_no;
                link.link_id = linkID;

                assignment.g_link_id_map[link.link_id] = 1;

                string movement_str;
                parser_link.GetValueByFieldName("movement_str", movement_str, false);
                parser_link.GetValueByFieldName("geometry", link.geometry, false);

                if(link.geometry.size() == 0)
                {
                    link.geometry = "LINESTRING (" + std::to_string(+g_node_vector[internal_from_node_seq_no].x) + " " + std::to_string(g_node_vector[internal_from_node_seq_no].y) + "," +
                        std::to_string(g_node_vector[internal_to_node_seq_no].x) + " " + std::to_string(g_node_vector[internal_to_node_seq_no].y) + "),";

                }

                parser_link.GetValueByFieldName("name", link.name, false);

                // and valid
                if (movement_str.size() > 0)
                {
                    int main_node_id = -1;


                    link.movement_str = movement_str;
                    link.main_node_id = main_node_id;
                }

                // Peiheng, 05/13/21, if setting.csv does not have corresponding link type or the whole section is missing, set it as 2 (i.e., Major arterial)
                int link_type = 2;
                parser_link.GetValueByFieldName("link_type", link_type, false);
                parser_link.GetValueByFieldName("link_type_name", link.link_type_name,false);
                parser_link.GetValueByFieldName("link_type_code", link.link_type_code, false);

                string TMC_code;

                parser_link.GetValueByFieldName("tmc", link.TMC_code, false);

                if(link.TMC_code.size()>0)
                {
                parser_link.GetValueByFieldName("tmc_corridor_name", link.tmc_corridor_name, false);
                link.tmc_corridor_id = 1;
                link.tmc_road_sequence = 1;
                parser_link.GetValueByFieldName("tmc_corridor_id", link.tmc_corridor_id, false);
                parser_link.GetValueByFieldName("tmc_road_sequence", link.tmc_road_sequence, false);
                }

                
                
                for (int tau = 0; tau <= 4; tau++)
                {
                    link.VDF_STA_speed[tau] = -1;
                    link.VDF_STA_VOC[tau] = -1;
                    link.VDF_STA_volume[tau] = -1;
                }

                parser_link.GetValueByFieldName("VDF_STA_speed1", link.VDF_STA_speed[1], false);
                parser_link.GetValueByFieldName("VDF_STA_speed2", link.VDF_STA_speed[2], false);
                parser_link.GetValueByFieldName("VDF_STA_speed3", link.VDF_STA_speed[3], false);
                parser_link.GetValueByFieldName("VDF_STA_speed4", link.VDF_STA_speed[4], false);

                parser_link.GetValueByFieldName("VDF_STA_VOC1", link.VDF_STA_VOC[1], false);
                parser_link.GetValueByFieldName("VDF_STA_VOC2", link.VDF_STA_VOC[2], false);
                parser_link.GetValueByFieldName("VDF_STA_VOC3", link.VDF_STA_VOC[3], false);
                parser_link.GetValueByFieldName("VDF_STA_VOC4", link.VDF_STA_VOC[4], false);

                parser_link.GetValueByFieldName("VDF_STA_volume1", link.VDF_STA_volume[1], false);
                parser_link.GetValueByFieldName("VDF_STA_volume2", link.VDF_STA_volume[2], false);
                parser_link.GetValueByFieldName("VDF_STA_volume3", link.VDF_STA_volume[3], false);
                parser_link.GetValueByFieldName("VDF_STA_volume4", link.VDF_STA_volume[4], false);

                parser_link.GetValueByFieldName("FT", link.FT, false);
                parser_link.GetValueByFieldName("AT", link.AT, false);

                if (assignment.g_LinkTypeMap.find(link_type) == assignment.g_LinkTypeMap.end() && waring_message_link_type_map.find(link_type) == waring_message_link_type_map.end())
                {
                    dtalog.output() << "link type " << link_type << " in link.csv is not defined for link " << from_node_id << "->" << to_node_id << " in link_type.csv" << endl;
                    waring_message_link_type_map[link_type] = 1;
                    // link.link_type has been taken care by its default constructor
                    //g_ProgramStop();
                }
                else
                {
                    // link type should be defined in settings.csv
                    link.link_type = link_type;
                }

                if (assignment.g_LinkTypeMap[link.link_type].type_code == "c")  // suggestion: we can move "c" as "connector" in allowed_uses
                {
                    if (g_node_vector[internal_from_node_seq_no].zone_org_id >= 0)
                    {
                        int zone_org_id = g_node_vector[internal_from_node_seq_no].zone_org_id;
                        if (assignment.g_zoneid_to_zone_seq_no_mapping.find(zone_org_id) != assignment.g_zoneid_to_zone_seq_no_mapping.end())
                            link.zone_seq_no_for_outgoing_connector = assignment.g_zoneid_to_zone_seq_no_mapping[zone_org_id];
                    }
                }

                double length = 1.0; // km or mile
                double free_speed = 1.0;
                double k_jam = 200;
                double bwtt_speed = 12;  //miles per hour
                int number_of_lanes = 1;

                double lane_capacity = 1800;
                float loading_ratio = 0.2;


                if (link_type < 20)
                {
                    free_speed = free_speed_vector[link_type];
                    number_of_lanes = nlanes_vector[link_type];
                    lane_capacity = lane_capacity_vector[link_type];
                }

                parser_link.GetValueByFieldName("length", length);
                parser_link.GetValueByFieldName("free_speed", free_speed);

                link.free_speed = free_speed;
                free_speed = max(0.1, free_speed);


                parser_link.GetValueByFieldName("lanes", number_of_lanes);
                parser_link.GetValueByFieldName("capacity", lane_capacity);


                //loading ratio
                loading_ratio = 0.01;
                if (lane_capacity > 2500)
                    loading_ratio = 0.1; // connectors
                else if (lane_capacity > 1500)
                    loading_ratio = 0.4;
                else if (lane_capacity > 1000)
                    loading_ratio = 0.3;
                else if (lane_capacity > 800)
                    loading_ratio = 0.2;


                // use half of capacity as reference production and attraction
                if (g_node_vector[internal_from_node_seq_no].production < 0.01)  // no value
                {
                    float valid_lane_capacity = min(2000, lane_capacity);
                    float global_loading_factor = 0.5;  // for ODME
                    g_node_vector[internal_from_node_seq_no].production = valid_lane_capacity * number_of_lanes * loading_ratio* global_loading_factor;
                }
                if (g_node_vector[internal_to_node_seq_no].attraction < 0.01)  // no value
                {
                    float valid_lane_capacity = min(2000, lane_capacity);

                    g_node_vector[internal_to_node_seq_no].attraction = valid_lane_capacity * number_of_lanes * loading_ratio;

                }

                link.free_flow_travel_time_in_min = length / 1609 / free_speed * 60;  // assumption: length is in meter

                link.traffic_flow_code = assignment.g_LinkTypeMap[link.link_type].traffic_flow_code;

                //spatial queue and kinematic wave
                if (link.traffic_flow_code >= 2)
                    link.spatial_capacity_in_vehicles = max(1.0, length * number_of_lanes * k_jam);

                // kinematic wave
                if (link.traffic_flow_code == 3)
                    link.BWTT_in_simulation_interval = length / bwtt_speed * 3600 / number_of_seconds_per_interval;

                // Peiheng, 02/03/21, useless block
                if (linkID == "10")
                    int i_debug = 1;

                char VDF_field_name[20];

                for (int at = 0; at < assignment.g_AgentTypeVector.size(); at++)
                {
                    double pce_at = 1; // default
                    sprintf(VDF_field_name, "VDF_pce%s", assignment.g_AgentTypeVector[at].agent_type.c_str());

                    parser_link.GetValueByFieldName(VDF_field_name, pce_at, false, true);

                    for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
                    {
                        link.VDF_period[tau].pce[at] = pce_at;
                    }

                }
                    link.UpdateKC(link.free_speed);

                for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
                {
                    //setup default values
                    link.VDF_period[tau].capacity = lane_capacity * number_of_lanes;
                    link.VDF_period[tau].FFTT = length / 1609 / free_speed * 60.0;  // 60.0 for 60 min per hour // assumption: length is in meter
                    link.VDF_period[tau].alpha = 0.15;
                    link.VDF_period[tau].beta = 4;
                    link.VDF_period[tau].preload = 0;

                    for (int at = 0; at < assignment.g_AgentTypeVector.size(); at++)
                    {
                        link.VDF_period[tau].toll[at] = 0;
                    }


                    link.VDF_period[tau].starting_time_slot_no = assignment.g_DemandPeriodVector[tau].starting_time_slot_no;
                    link.VDF_period[tau].ending_time_slot_no = assignment.g_DemandPeriodVector[tau].ending_time_slot_no;

                    int demand_period_id = assignment.g_DemandPeriodVector[tau].demand_period_id;
                    sprintf(VDF_field_name, "VDF_fftt%d", demand_period_id);
                    parser_link.GetValueByFieldName(VDF_field_name, link.VDF_period[tau].FFTT, false, false);  // FFTT should be per min

                    sprintf(VDF_field_name, "VDF_cap%d", demand_period_id);
                    parser_link.GetValueByFieldName(VDF_field_name, link.VDF_period[tau].capacity, false, false);  // capacity should be per period per link (include all lanes)

                    sprintf(VDF_field_name, "VDF_alpha%d", demand_period_id);
                    parser_link.GetValueByFieldName(VDF_field_name, link.VDF_period[tau].alpha, false, false);

                    sprintf(VDF_field_name, "VDF_beta%d", demand_period_id);
                    parser_link.GetValueByFieldName(VDF_field_name, link.VDF_period[tau].beta, false, false);

                    sprintf(VDF_field_name, "VDF_allowed_uses%d", demand_period_id);
                    parser_link.GetValueByFieldName(VDF_field_name, link.VDF_period[tau].allowed_uses, false);

                    sprintf(VDF_field_name, "VDF_preload%d", demand_period_id);
                    parser_link.GetValueByFieldName(VDF_field_name, link.VDF_period[tau].preload, false, false);

                    for (int at = 0; at < assignment.g_AgentTypeVector.size(); at++)
                    {
                        sprintf(VDF_field_name, "VDF_toll%s%d", assignment.g_AgentTypeVector[at].agent_type.c_str(), demand_period_id);
                        parser_link.GetValueByFieldName(VDF_field_name, link.VDF_period[tau].toll[at], false, false);

                        if (link.VDF_period[tau].toll[at] > 0.001)
                        {
                            dtalog.output() << "link " << from_node_id << "->" << to_node_id << " has a toll of " << link.VDF_period[tau].toll[at] << " for agent type "
                                << assignment.g_AgentTypeVector[at].agent_type.c_str() << " at demand period " << demand_period_id << endl;
                        }
                    }

                    sprintf(VDF_field_name, "VDF_penalty%d", demand_period_id);
                    parser_link.GetValueByFieldName(VDF_field_name, link.VDF_period[tau].penalty, false, false);

                    sprintf(VDF_field_name, "VDF_PHF%d", demand_period_id);
                    parser_link.GetValueByFieldName(VDF_field_name, link.VDF_period[tau].PHF, false, false);

                    sprintf(VDF_field_name, "VDF_mu%d", demand_period_id);
                    parser_link.GetValueByFieldName(VDF_field_name, link.VDF_period[tau].mu, false, false);  // mu should be per hour per link, so that we can calculate congestion duration and D/mu in BPR-X

                    //sprintf(VDF_field_name, "VDF_gamma%d", demand_period_id);  // remove gamma
                    //parser_link.GetValueByFieldName(VDF_field_name, link.VDF_period[tau].gamma);
                }

                // for each period

                float default_cap = 1000;
                float default_BaseTT = 1;

                // setup default value
                for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
                {
                    link.TDBaseTT[tau] = default_BaseTT;
                    link.TDBaseCap[tau] = default_cap;
                }

                //link.m_OutflowNumLanes = number_of_lanes;//visum lane_cap is actually link_cap

                link.number_of_lanes = number_of_lanes;
                link.lane_capacity = lane_capacity;
                link.link_spatial_capacity = k_jam * number_of_lanes * length;

                link.length = length;
                for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
                    link.travel_time_per_period[tau] = length / free_speed * 60;

                // min // calculate link cost based length and speed limit // later we should also read link_capacity, calculate delay

                //int sequential_copying = 0;
                //
                //parser_link.GetValueByFieldName("sequential_copying", sequential_copying);

                g_node_vector[internal_from_node_seq_no].m_outgoing_link_seq_no_vector.push_back(link.link_seq_no);  // add this link to the corresponding node as part of outgoing node/link
                g_node_vector[internal_to_node_seq_no].m_incoming_link_seq_no_vector.push_back(link.link_seq_no);  // add this link to the corresponding node as part of outgoing node/link

                g_node_vector[internal_from_node_seq_no].m_to_node_seq_no_vector.push_back(link.to_node_seq_no);  // add this link to the corresponding node as part of outgoing node/link
                g_node_vector[internal_from_node_seq_no].m_to_node_2_link_seq_no_map[link.to_node_seq_no] = link.link_seq_no;  // add this link to the corresponding node as part of outgoing node/link

                g_link_vector.push_back(link);

                assignment.g_number_of_links++;

                if (assignment.g_number_of_links % 10000 == 0)
                    dtalog.output() << "reading " << assignment.g_number_of_links << " links.. " << endl;
            }

            parser_link.CloseCSVFile();
        }
        // we now know the number of links
        dtalog.output() << "number of links = " << assignment.g_number_of_links << endl;



        // after we read the physical links
        // we create virtual connectors
        for (int i = 0; i < g_node_vector.size(); ++i)
        {
            if (g_node_vector[i].zone_org_id >= 0) // for each physical node
            { // we need to make sure we only create two way connectors between nodes and zones

                int internal_from_node_seq_no, internal_to_node_seq_no, zone_seq_no;

                internal_from_node_seq_no = g_node_vector[i].node_seq_no;
                internal_to_node_seq_no = assignment.zone_id_to_centriod_node_no_mapping[g_node_vector[i].zone_org_id];
                zone_seq_no = assignment.g_zoneid_to_zone_seq_no_mapping[g_node_vector[i].zone_org_id];

                // incomming virtual connector
                g_AddNewVirtualConnectorLink(internal_from_node_seq_no, internal_to_node_seq_no, -1);
                // outgoing virtual connector
                g_AddNewVirtualConnectorLink(internal_to_node_seq_no, internal_from_node_seq_no, zone_seq_no);
            }
        }

        dtalog.output() << "number of links + connectors =" << assignment.g_number_of_links << endl;

        if (dtalog.debug_level() == 2)
        {
            for (int i = 0; i < g_node_vector.size(); ++i)
            {
                if (g_node_vector[i].zone_org_id > 0) // for each physical node
                {
                    // we need to make sure we only create two way connectors between nodes and zones
                    dtalog.output() << "node id= " << g_node_vector[i].node_id << " with zone id " << g_node_vector[i].zone_org_id << "and "
                        << g_node_vector[i].m_outgoing_link_seq_no_vector.size() << " outgoing links." << endl;
                    for (int j = 0; j < g_node_vector[i].m_outgoing_link_seq_no_vector.size(); ++j)
                    {
                        int link_seq_no = g_node_vector[i].m_outgoing_link_seq_no_vector[j];
                        dtalog.output() << "  outgoing node = " << g_node_vector[g_link_vector[link_seq_no].to_node_seq_no].node_id << endl;
                    }
                }
                else
                {
                    if (dtalog.debug_level() == 3)
                    {
                        dtalog.output() << "node id= " << g_node_vector[i].node_id << " with " << g_node_vector[i].m_outgoing_link_seq_no_vector.size() << " outgoing links." << endl;
                        for (int j = 0; j < g_node_vector[i].m_outgoing_link_seq_no_vector.size(); ++j)
                        {
                            int link_seq_no = g_node_vector[i].m_outgoing_link_seq_no_vector[j];
                            dtalog.output() << "  outgoing node = " << g_node_vector[g_link_vector[link_seq_no].to_node_seq_no].node_id << endl;
                        }
                    }
                }
            }
        }

        assignment.Mapping_TMC_Identification();
        bool ReadingDataReady = assignment.Map_TMC_Reading();
        g_OutputTMCFiles(ReadingDataReady);

        g_OutputTMCScenarioFiles();

        CCSVParser parser_movement;
        int prohibited_count = 0;

        if (parser_movement.OpenCSVFile("movement.csv", false))  // not required
        {
            while (parser_movement.ReadRecord())
            {
                string ib_link_id;
                int node_id = 0;
                string ob_link_id;
                int prohibited_flag = 0;

                if (!parser_movement.GetValueByFieldName("node_id", node_id))
                    break;

                if (assignment.g_node_id_to_seq_no_map.find(node_id) == assignment.g_node_id_to_seq_no_map.end())
                {
                    dtalog.output() << "Error: node_id " << node_id << " in file movement.csv is not defined in node.csv." << endl;
                    //has not been defined
                    continue;
                }

                parser_movement.GetValueByFieldName("ib_link_id", ib_link_id);
                parser_movement.GetValueByFieldName("ob_link_id", ob_link_id);

                if (assignment.g_link_id_map.find(ib_link_id) != assignment.g_link_id_map.end())
                    dtalog.output() << "Error: ib_link_id " << ib_link_id.c_str() << " has not been defined in movement.csv. Please check link.csv." << endl;

                if (assignment.g_link_id_map.find(ob_link_id) != assignment.g_link_id_map.end())
                    dtalog.output() << "Error: ob_link_id " << ob_link_id.c_str() << " has not been defined in movement.csv. Please check link.csv." << endl;

                float penalty = 0;
                parser_movement.GetValueByFieldName("penalty", penalty);

                if (penalty >= 99)
                {
                    string	movement_string;
                    movement_string = ib_link_id + "->" + ob_link_id;

                    int node_no = assignment.g_node_id_to_seq_no_map[node_id];
                    g_node_vector[node_no].prohibited_movement_size++;
                    g_node_vector[node_no].m_prohibited_movement_string_map[movement_string] = 1;

                    prohibited_count++;
                }
            }

            dtalog.output() << "Step XX: Reading movement.csv data with " << prohibited_count << " prohibited records." << endl;
            parser_movement.CloseCSVFile();
        }

        // tally total zone based production and attraction for QEM and ODME mode

        if (assignment.assignment_mode == 10 || assignment.assignment_mode == 3)
        {
            std::map<int, int>::iterator it;
            for (int at = 0; at < assignment.g_AgentTypeVector.size(); ++at)
            { 
                for (int d = 0; d < g_zone_vector.size(); ++d)
                {
                    g_zone_vector[d].gravity_production[at] = 0;
                    g_zone_vector[d].gravity_attraction[at] = 0;

                }

            }

                for (int i = 0; i < g_node_vector.size(); i++)
            {
                if (g_node_vector[i].zone_org_id > 0)
                {
                    int zone_seq_no = assignment.g_zoneid_to_zone_seq_no_mapping[g_node_vector[i].zone_org_id];
                    g_zone_vector[zone_seq_no].m_activity_node_vector.push_back(g_node_vector[i].node_id);

                    for (int at = 0; at < assignment.g_AgentTypeVector.size(); ++at)
                    {

                        g_zone_vector[zone_seq_no].gravity_attraction[at] += g_node_vector[i].attraction * assignment.g_AgentTypeVector[at].trip_ratio;
                        g_zone_vector[zone_seq_no].gravity_production[at] += g_node_vector[i].production * assignment.g_AgentTypeVector[at].trip_ratio;
                    }
                }
             
            }



            FILE* g_pFileZone = nullptr;
            g_pFileZone = fopen("model_zone.csv", "w");

            if (g_pFileZone == NULL)
            {
                cout << "File model_zone.csv cannot be opened." << endl;
                g_ProgramStop();
            }
            else
            {


                fprintf(g_pFileZone, "zone_id,node_vector,cell_id,");

                for (int at = 0; at < assignment.g_AgentTypeVector.size(); ++at)
                {
                    fprintf(g_pFileZone,"%s_production, %s_attraction,", assignment.g_AgentTypeVector[at].agent_type.c_str(), assignment.g_AgentTypeVector[at].agent_type.c_str());
                }

               fprintf(g_pFileZone, ",geometry\n");

                std::map<int, int>::iterator it;

                for (int d = 0; d < g_zone_vector.size(); ++d)
                {

                    int x_i = floor(g_zone_vector[d].cell_x / assignment.m_GridResolution);
                    int y_i = floor(g_zone_vector[d].cell_y / assignment.m_GridResolution);

                    double x_coord_left = x_i * assignment.m_GridResolution;
                    double y_coord_bottom = y_i * assignment.m_GridResolution;
                    double x_coord_right = x_coord_left + assignment.m_GridResolution;
                    double y_coord_top = y_coord_bottom + assignment.m_GridResolution;

                    fprintf(g_pFileZone, "%d,", g_zone_vector[d].zone_id);
                    for (int an = 0; an < g_zone_vector[d].m_activity_node_vector.size(); ++an)
                    {
                        fprintf(g_pFileZone, "%d;", g_zone_vector[d].m_activity_node_vector[an]);
                    }


                    fprintf(g_pFileZone, ",%jd",  g_zone_vector[d].cell_id);

                    for (int at = 0; at < assignment.g_AgentTypeVector.size(); ++at)
                    {
                    fprintf(g_pFileZone, "%f,%f,", g_zone_vector[d].gravity_production[at], g_zone_vector[d].gravity_attraction[at]);
                    }


                    fprintf(g_pFileZone, "\"LINESTRING (");

                    fprintf(g_pFileZone, "%f %f,", x_coord_left, y_coord_top);
                    fprintf(g_pFileZone, "%f %f,", x_coord_right, y_coord_top);
                    fprintf(g_pFileZone, "%f %f,", x_coord_right, y_coord_bottom);
                    fprintf(g_pFileZone, "%f %f,", x_coord_left, y_coord_bottom);
                    fprintf(g_pFileZone, "%f %f,", x_coord_left, y_coord_top);
                    fprintf(g_pFileZone, ")\"");
                    fprintf(g_pFileZone, "\n");
                }
                fclose(g_pFileZone);
            }
        }

        if (assignment.assignment_mode == 10)
        {
            g_TripGeneration(assignment);

            dtalog.output() << "writing demand_geo.csv.." << endl;

            FILE* g_pFileODMatrix = nullptr;
            fopen_ss(&g_pFileODMatrix, "demand_geo.csv", "w");

            if (!g_pFileODMatrix)
            {
                dtalog.output() << "File demand_geo.csv cannot be opened." << endl;
                g_ProgramStop();
            }
            else
            {

                fprintf(g_pFileODMatrix, "demand_period,time_period,agent_type,o_zone_id,d_zone_id,volume,geometry\n");
                int demand_writing_log_count = 0;
                // reset the estimated production and attraction
                for (int orig = 0; orig < g_zone_vector.size(); ++orig)  // o
                {
                    for (int dest = 0; dest < g_zone_vector.size(); ++dest)  // d
                    {
                        for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); ++tau)
                        {
                            for (int at = 0; at < assignment.g_AgentTypeVector.size(); ++at)
                            {
                                if (g_zone_vector[orig].gravity_production[at] >= 0)
                                {
                                    if (g_zone_vector[dest].gravity_attraction[at] > 0)
                                    {
                                        float value = 0;
                                        if (g_zone_vector[orig].m_ODMatrix[at][tau].value_map.find(dest) != g_zone_vector[orig].m_ODMatrix[at][tau].value_map.end())
                                        {
                                            value = g_zone_vector[orig].m_ODMatrix[at][tau].value_map[dest];
                                        }

                                        if (value > 0.000001)
                                        {
                                            fprintf(g_pFileODMatrix, "%s,%s,%s,%d,%d,%.4f,", assignment.g_DemandPeriodVector[tau].demand_period.c_str(), assignment.g_DemandPeriodVector[tau].time_period.c_str(),

                                                assignment.g_AgentTypeVector[at].agent_type.c_str(), g_zone_vector[orig].zone_id, g_zone_vector[dest].zone_id, value);

                                            if (demand_writing_log_count < 100)
                                            {
                                                dtalog.output() << "orig= " << g_zone_vector[orig].zone_id << " dest= " << g_zone_vector[dest].zone_id << ":" << value << endl;
                                            }
                                            demand_writing_log_count++;

                                            fprintf(g_pFileODMatrix, "\"LINESTRING (");

                                            fprintf(g_pFileODMatrix, "%f %f,", g_zone_vector[orig].cell_x, g_zone_vector[orig].cell_y);
                                            fprintf(g_pFileODMatrix, "%f %f,", g_zone_vector[dest].cell_x, g_zone_vector[dest].cell_y);
                                            fprintf(g_pFileODMatrix, ")\"");
                                            fprintf(g_pFileODMatrix, "\n");
                                        }

                                    }
                                }

                            }
                        }
                    }
                }

                fclose(g_pFileODMatrix);
            }


            ////////////////////////////
            dtalog.output() << "writing demand_mtx.csv.." << endl;


            fopen_ss(&g_pFileODMatrix, "demand_mtx.csv", "w");

            if (!g_pFileODMatrix)
            {
                dtalog.output() << "File demand_mtx.csv cannot be opened." << endl;
                g_ProgramStop();
            }
            else
            {
                for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); ++tau)
                {
                    for (int at = 0; at < assignment.g_AgentTypeVector.size(); ++at)
                    {

                        fprintf(g_pFileODMatrix, "demand_period,agent_type,OD,");
                        for (int d = 0; d < g_zone_vector.size(); ++d)
                        {
                            fprintf(g_pFileODMatrix, "%d,", g_zone_vector[d].zone_id);
                        }


                        fprintf(g_pFileODMatrix, "subtotal_est,subtotal_target,ratio_diff\n");

                        // reset the estimated production and attraction
                        for (int orig = 0; orig < g_zone_vector.size(); ++orig)  // o
                        {
                            fprintf(g_pFileODMatrix, "%s,%s,", assignment.g_DemandPeriodVector[tau].demand_period.c_str(),                       
                                assignment.g_AgentTypeVector[at].agent_type.c_str());
                            float total_production = 0;

                            fprintf(g_pFileODMatrix, "%d,", g_zone_vector[orig].zone_id);

                            for (int dest = 0; dest < g_zone_vector.size(); ++dest)  // d
                            {
                                float value = 0;
                                if (g_zone_vector[orig].gravity_production[at] >= 0 && g_zone_vector[dest].gravity_attraction[at] > 0)
                                {
                                    value = g_zone_vector[orig].m_ODMatrix[at][tau].value_map[dest];
                                }
                                total_production += value;
                                fprintf(g_pFileODMatrix, "%f,", value);

                            }
                            float percentage_difference = (total_production - g_zone_vector[orig].gravity_production[at]) / max(0.001, g_zone_vector[orig].gravity_production[at]);
                            fprintf(g_pFileODMatrix, "%f,%f,%f\n", total_production, g_zone_vector[orig].gravity_production[at], percentage_difference);

                        }

                        fprintf(g_pFileODMatrix, "est,");

                        for (int dest = 0; dest < g_zone_vector.size(); ++dest)  // d
                        {
                            fprintf(g_pFileODMatrix, "%f,", g_zone_vector[dest].gravity_est_attraction[at]);
                        }
                        fprintf(g_pFileODMatrix, "\n");
                        fprintf(g_pFileODMatrix, "target,");

                        for (int dest = 0; dest < g_zone_vector.size(); ++dest)  // d
                        {
                            fprintf(g_pFileODMatrix, "%f,", g_zone_vector[dest].gravity_attraction[at]);
                        }

                        fprintf(g_pFileODMatrix, "\n");
                        fprintf(g_pFileODMatrix, "ratio_diff,");

                        for (int dest = 0; dest < g_zone_vector.size(); ++dest)  // d
                        {
                            float percentage_difference = (g_zone_vector[dest].gravity_est_attraction[at] - g_zone_vector[dest].gravity_attraction[at]) / max(0.001, g_zone_vector[dest].gravity_attraction[at]);

                            fprintf(g_pFileODMatrix, "%f,", percentage_difference);
                        }
                        fprintf(g_pFileODMatrix, "\n");
                    }
                }
            

            fclose(g_pFileODMatrix);
        }


        ////////////////////set reference half of link capacity as link volume

        //for (int l = 0; l < g_link_vector.size(); ++l)
        //{
        //    if(g_link_vector[l].lane_capacity >=800 && g_link_vector[l].lane_capacity < 2000) // remark: link_type = -1 is virtual connector
        //    {
        //        g_link_vector[l].obs_count = g_link_vector[l].lane_capacity * g_link_vector[l].number_of_lanes/2;
        //    }
        //}

        for (int z = 0; z < g_zone_vector.size(); ++z)  // d
        {
            g_zone_vector[z].obs_production = -1;
            g_zone_vector[z].obs_attraction = -1;  // invalidate the data, we will focus on link count data first 
        }
    


        assignment.assignment_mode = 3; // // end of QEM model, reset to ODME;
       } 

       g_OutputModelFiles();
       
       

}

void g_reload_service_arc_data(Assignment& assignment)
{
    dtalog.output() << "Step 1.7: Reading service arc in timing.csv..." << endl;

    CCSVParser parser_service_arc;
    if (parser_service_arc.OpenCSVFile("timing.csv", false))
    {
        while (parser_service_arc.ReadRecord())  // if this line contains [] mark, then we will also read field headers.
        {
            int from_node_id = 0;
            if (!parser_service_arc.GetValueByFieldName("from_node_id", from_node_id))
            {
                dtalog.output() << "Error: from_node_id in file timing.csv is not defined." << endl;
                continue;
            }

            int to_node_id = 0;
            if (!parser_service_arc.GetValueByFieldName("to_node_id", to_node_id))
            {
                continue;
            }

            if (assignment.g_node_id_to_seq_no_map.find(from_node_id) == assignment.g_node_id_to_seq_no_map.end())
            {
                dtalog.output() << "Error: from_node_id " << from_node_id << " in file timing.csv is not defined in node.csv." << endl;
                //has not been defined
                continue;
            }
            if (assignment.g_node_id_to_seq_no_map.find(to_node_id) == assignment.g_node_id_to_seq_no_map.end())
            {
                dtalog.output() << "Error: to_node_id " << to_node_id << " in file timing.csv is not defined in node.csv." << endl;
                    //has not been defined
                continue;
            }

            // map external node number to internal node seq no.
            int internal_from_node_seq_no = assignment.g_node_id_to_seq_no_map[from_node_id];
            int internal_to_node_seq_no = assignment.g_node_id_to_seq_no_map[to_node_id];

            // create a link object
            CServiceArc service_arc;

            if (g_node_vector[internal_from_node_seq_no].m_to_node_2_link_seq_no_map.find(internal_to_node_seq_no) != g_node_vector[internal_from_node_seq_no].m_to_node_2_link_seq_no_map.end())
            {
                service_arc.link_seq_no = g_node_vector[internal_from_node_seq_no].m_to_node_2_link_seq_no_map[internal_to_node_seq_no];
                g_link_vector[service_arc.link_seq_no].service_arc_flag = true;
            }
            else
            {
                dtalog.output() << "Error: Link " << from_node_id << "->" << to_node_id << " in file timing.csv is not defined in link.csv." << endl;
                continue;
            }

            string time_period;
            if (!parser_service_arc.GetValueByFieldName("time_window", time_period))
            {
                dtalog.output() << "Error: Field time_window in file timing.csv cannot be read." << endl;
                g_ProgramStop();
                break;
            }

            vector<float> global_minute_vector;

            //input_string includes the start and end time of a time period with hhmm format
            global_minute_vector = g_time_parser(time_period); //global_minute_vector incldue the starting and ending time
            if (global_minute_vector.size() == 2)
            {
                if (global_minute_vector[0] < assignment.g_LoadingStartTimeInMin)
                    global_minute_vector[0] = assignment.g_LoadingStartTimeInMin;

                if (global_minute_vector[0] > assignment.g_LoadingEndTimeInMin)
                    global_minute_vector[0] = assignment.g_LoadingEndTimeInMin;

                if (global_minute_vector[1] < assignment.g_LoadingStartTimeInMin)
                    global_minute_vector[1] = assignment.g_LoadingStartTimeInMin;

                if (global_minute_vector[1] > assignment.g_LoadingEndTimeInMin)
                    global_minute_vector[1] = assignment.g_LoadingEndTimeInMin;

                if (global_minute_vector[1] < global_minute_vector[0])
                    global_minute_vector[1] = global_minute_vector[0];

                //this could contain sec information.
                service_arc.starting_time_no = (global_minute_vector[0] - assignment.g_LoadingStartTimeInMin) * 60 / number_of_seconds_per_interval;
                service_arc.ending_time_no = (global_minute_vector[1] - assignment.g_LoadingStartTimeInMin) * 60 / number_of_seconds_per_interval;
            }
            else
                continue;

            float time_interval = 0;
            parser_service_arc.GetValueByFieldName("time_interval", time_interval, false, false);
            service_arc.time_interval_no = max(1.0, time_interval * 60.0 / number_of_seconds_per_interval);

            int travel_time_delta_in_min = 0;
            parser_service_arc.GetValueByFieldName("travel_time_delta", travel_time_delta_in_min, false, false);
            service_arc.travel_time_delta =
                max((int) g_link_vector[service_arc.link_seq_no].free_flow_travel_time_in_min * 60 / number_of_seconds_per_interval,
                travel_time_delta_in_min * 60 / number_of_seconds_per_interval);

            service_arc.travel_time_delta = max(1, service_arc.travel_time_delta);

            // capacity in the space time arcs
            float capacity = 1;
            parser_service_arc.GetValueByFieldName("capacity", capacity);
            service_arc.capacity = max(0.0f, capacity);

            // capacity in the space time arcs
            parser_service_arc.GetValueByFieldName("cycle_length", service_arc.cycle_length);

            // capacity in the space time arcs
            parser_service_arc.GetValueByFieldName("red_time", service_arc.red_time);

            for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
            {
                    // to do: we need to consider multiple periods in the future, Xuesong Zhou, August 20, 2020.
                g_link_vector[service_arc.link_seq_no].VDF_period[tau].red_time = service_arc.red_time;
                g_link_vector[service_arc.link_seq_no].VDF_period[tau].cycle_length = service_arc.cycle_length;
            }

            g_service_arc_vector.push_back(service_arc);
            assignment.g_number_of_timing_arcs++;

            if (assignment.g_number_of_timing_arcs % 10000 == 0)
                dtalog.output() << "reading " << assignment.g_number_of_timing_arcs << " timing_arcs.. " << endl;
        }

        parser_service_arc.CloseCSVFile();
    }

    dtalog.output() << endl;
    dtalog.output() << "Step 1.8: Reading file section [demand_file_list] in setting.csv..." << endl;

    CCSVParser parser;
    parser.IsFirstLineHeader = false;

    if (parser.OpenCSVFile("settings.csv", false))
    {
        while (parser.ReadRecord_Section())
        {
            if (parser.SectionName == "[capacity_scenario]")
            {
                int from_node_id = 0;
                if (!parser.GetValueByFieldName("from_node_id", from_node_id))
                {
                    dtalog.output() << "Error: from_node_id in file timing.csv is not defined." << endl;
                    continue;
                }

                int to_node_id = 0;
                if (!parser.GetValueByFieldName("to_node_id", to_node_id))
                    continue;

                if (assignment.g_node_id_to_seq_no_map.find(from_node_id) == assignment.g_node_id_to_seq_no_map.end())
                {
                    dtalog.output() << "Error: from_node_id " << from_node_id << " in file timing.csv is not defined in node.csv." << endl;
                    //has not been defined
                    continue;
                }
                if (assignment.g_node_id_to_seq_no_map.find(to_node_id) == assignment.g_node_id_to_seq_no_map.end())
                {
                    dtalog.output() << "Error: to_node_id " << to_node_id << " in file timing.csv is not defined in node.csv." << endl;
                    //has not been defined
                    continue;
                }

                // create a link object
                CServiceArc service_arc;
                // map external node number to internal node seq no.
                int internal_from_node_seq_no = assignment.g_node_id_to_seq_no_map[from_node_id];
                int internal_to_node_seq_no = assignment.g_node_id_to_seq_no_map[to_node_id];

                if (g_node_vector[internal_from_node_seq_no].m_to_node_2_link_seq_no_map.find(internal_to_node_seq_no) != g_node_vector[internal_from_node_seq_no].m_to_node_2_link_seq_no_map.end())
                {
                    service_arc.link_seq_no = g_node_vector[internal_from_node_seq_no].m_to_node_2_link_seq_no_map[internal_to_node_seq_no];
                    g_link_vector[service_arc.link_seq_no].service_arc_flag = true;
                }
                else
                {
                    dtalog.output() << "Error: Link " << from_node_id << "->" << to_node_id << " in file timing.csv is not defined in link.csv." << endl;
                    continue;
                }

                string time_period;
                if (!parser.GetValueByFieldName("time_window", time_period))
                {
                    dtalog.output() << "Error: Field time_window in file timing.csv cannot be read." << endl;
                    g_ProgramStop();
                    break;
                }

                vector<float> global_minute_vector;

                //input_string includes the start and end time of a time period with hhmm format
                global_minute_vector = g_time_parser(time_period); //global_minute_vector incldue the starting and ending time
                if (global_minute_vector.size() == 2)
                {
                    if (global_minute_vector[0] < assignment.g_LoadingStartTimeInMin)
                        global_minute_vector[0] = assignment.g_LoadingStartTimeInMin;

                    if (global_minute_vector[0] > assignment.g_LoadingEndTimeInMin)
                        global_minute_vector[0] = assignment.g_LoadingEndTimeInMin;

                    if (global_minute_vector[1] < assignment.g_LoadingStartTimeInMin)
                        global_minute_vector[1] = assignment.g_LoadingStartTimeInMin;

                    if (global_minute_vector[1] > assignment.g_LoadingEndTimeInMin)
                        global_minute_vector[1] = assignment.g_LoadingEndTimeInMin;

                    if (global_minute_vector[1] < global_minute_vector[0])
                        global_minute_vector[1] = global_minute_vector[0];

                    //this could contain sec information.
                    service_arc.starting_time_no = (global_minute_vector[0] - assignment.g_LoadingStartTimeInMin) * 60 / number_of_seconds_per_interval;
                    service_arc.ending_time_no = (global_minute_vector[1] - assignment.g_LoadingStartTimeInMin) * 60 / number_of_seconds_per_interval;
                }
                else
                    continue;

                float time_interval = 0;
                parser.GetValueByFieldName("time_interval", time_interval, false, false);
                service_arc.time_interval_no = max(1.0, time_interval * 60.0 / number_of_seconds_per_interval);

                int travel_time_delta_in_min = 0;
                parser.GetValueByFieldName("travel_time_delta", travel_time_delta_in_min, false, false);
                service_arc.travel_time_delta =
                    max((int)g_link_vector[service_arc.link_seq_no].free_flow_travel_time_in_min * 60 / number_of_seconds_per_interval,
                        travel_time_delta_in_min * 60 / number_of_seconds_per_interval);

                service_arc.travel_time_delta = max(1, service_arc.travel_time_delta);

                // capacity in the space time arcs
                float capacity = 1;
                parser.GetValueByFieldName("capacity", capacity);
                service_arc.capacity = max(0.0f, capacity);

                g_service_arc_vector.push_back(service_arc);
                assignment.g_number_of_timing_arcs++;

                dtalog.output() << "reading " << assignment.g_number_of_timing_arcs << " capacity reduction scenario.. " << endl;
            }
        }

        parser.CloseCSVFile();
    }
    // we now know the number of links
    dtalog.output() << "number of timing records = " << assignment.g_number_of_timing_arcs << endl << endl;
}

void g_reset_link_volume_in_master_program_without_columns(int number_of_links, int iteration_index, bool b_self_reducing_path_volume)
{
    int number_of_demand_periods = assignment.g_number_of_demand_periods;

    if(iteration_index == 0)
    {
        for (int i = 0; i < number_of_links; ++i)
        {
            for (int tau = 0; tau < number_of_demand_periods; ++tau)
            {
                // used in travel time calculation
                g_link_vector[i].flow_volume_per_period[tau] = 0;
            }
        }
    }
    else
    {
        for (int i = 0; i < number_of_links; ++i)
        {
            for (int tau = 0; tau < number_of_demand_periods; ++tau)
            {
                if (b_self_reducing_path_volume)
                {
                    // after link volumn "tally", self-deducting the path volume by 1/(k+1) (i.e. keep k/(k+1) ratio of previous flow)
                    // so that the following shortes path will be receiving 1/(k+1) flow
                    g_link_vector[i].flow_volume_per_period[tau] = g_link_vector[i].flow_volume_per_period[tau] * (float(iteration_index) / float(iteration_index + 1));
                }
            }
        }
    }
}

void g_reset_and_update_link_volume_based_on_columns(int number_of_links, int iteration_index, bool b_self_reducing_path_volume)
{
    for (int i = 0; i < number_of_links; ++i)
    {
        for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
        {
            // used in travel time calculation
            g_link_vector[i].flow_volume_per_period[tau] = 0;
            // reserved for BPR-X
            g_link_vector[i].queue_length_perslot[tau] = 0;

            for (int at = 0; at < assignment.g_AgentTypeVector.size(); ++at)
                g_link_vector[i].volume_per_period_per_at[tau][at] = 0;
        }
    }

    if(iteration_index >= 0)
    {
        for (int at = 0; at < assignment.g_AgentTypeVector.size(); ++at)  //m
        {
//#pragma omp parallel for

            std::map<int, CColumnPath>::iterator it;
            int zone_size = g_zone_vector.size();
            int tau_size = assignment.g_DemandPeriodVector.size();

            float link_volume_contributed_by_path_volume;

            int link_seq_no;
            float PCE_ratio;
            int nl;

            std::map<int, CColumnPath>::iterator it_begin;
            std::map<int, CColumnPath>::iterator it_end;

            int column_vector_size;
            CColumnVector* p_column;

            for (int orig = 0; orig < zone_size; ++orig)  // o
            {
                for (int dest = 0; dest < zone_size; ++dest) //d
                {
                    for (int tau = 0; tau < tau_size; ++tau)  //tau
                    {
                        p_column = &(assignment.g_column_pool[orig][dest][at][tau]);
                        if (p_column->od_volume > 0)
                        {
                            column_vector_size = p_column->path_node_sequence_map.size();

                            it_begin = p_column->path_node_sequence_map.begin();
                            it_end = p_column->path_node_sequence_map.end();
                            for (it = it_begin ; it != it_end; ++it)
                            {
                                link_volume_contributed_by_path_volume = it->second.path_volume;  // assign all OD flow to this first path

                                // add path volume to link volume
                                for (nl = 0; nl < it->second.m_link_size; ++nl)  // arc a
                                {
                                    link_seq_no = it->second.path_link_vector[nl];

                                    // MSA updating for the existing column pools
                                    // if iteration_index = 0; then update no flow discount is used (for the column pool case)
                                    PCE_ratio = g_link_vector[link_seq_no].VDF_period[tau].pce[at];  // updated on 08/16/2021 for link dependent and agent type dependent pce factor mainly for trucks 
                                    #pragma omp critical
                                    {
                                        g_link_vector[link_seq_no].flow_volume_per_period[tau] += link_volume_contributed_by_path_volume * PCE_ratio;
                                        g_link_vector[link_seq_no].volume_per_period_per_at[tau][at] += link_volume_contributed_by_path_volume;  // pure volume, not consider PCE
                                    }
                                }

                                // this  self-deducting action does not agents with fixed routing policies.
                                if(!p_column->bfixed_route && b_self_reducing_path_volume)
                                {
                                    //after link volumn "tally", self-deducting the path volume by 1/(k+1) (i.e. keep k/(k+1) ratio of previous flow) so that the following shortes path will be receiving 1/(k+1) flow
                                    it->second.path_volume = it->second.path_volume * (float(iteration_index) / float(iteration_index + 1));
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

double update_link_travel_time_and_cost()
{
    if (assignment.assignment_mode == 2)
    {
        //compute the time-dependent delay from simulation
        //for (int l = 0; l < g_link_vector.size(); l++)
        //{
        //	float volume = assignment.m_LinkCumulativeDeparture[l][assignment.g_number_of_simulation_intervals - 1];  // link flow rates
        //	float waiting_time_count = 0;

        //for (int tt = 0; tt < assignment.g_number_of_simulation_intervals; tt++)
        //{
        //	waiting_time_count += assignment.m_LinkTDWaitingTime[l][tt/number_of_interval_per_min];   // tally total waiting cou
        //}

        //for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); tau++)
        //{
        //	float travel_time = g_link_vector[l].free_flow_travel_time_in_min  + waiting_time_count* number_of_seconds_per_interval / max(1, volume) / 60;
        //	g_link_vector[l].travel_time_per_period[tau] = travel_time;

        //}
    }

#pragma omp parallel for
    for (int i = 0; i < g_link_vector.size(); ++i)
    {
        // step 1: travel time based on VDF
        g_link_vector[i].CalculateTD_VDFunction();

        for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); ++tau)
        {
            for (int at = 0; at < assignment.g_AgentTypeVector.size(); ++at)
            {
                float PCE_agent_type = assignment.g_AgentTypeVector[at].PCE;

                // step 2: marginal cost for SO
                g_link_vector[i].calculate_marginal_cost_for_agent_type(tau, at, PCE_agent_type);

                //if (g_debug_level  >= 3 && assignment.assignment_mode >= 2 && assignment.g_pFileDebugLog != NULL)
                //	fprintf(assignment.g_pFileDebugLog, "Update link cost: link %d->%d: tau = %d, at = %d, travel_marginal =  %.3f\n",

                //		g_node_vector[g_link_vector[l].from_node_seq_no].node_id,
                //		g_node_vector[g_link_vector[l].to_node_seq_no].node_id,
                //		tau, at,
                //		g_link_vector[l].travel_marginal_cost_per_period[tau][at]);
            }
        }
    }

	double total_network_travel_time = 0;
	for (int i = 0; i < g_link_vector.size(); ++i)
	{
		for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); ++tau)
		{
			total_network_travel_time += g_link_vector[i].VDF_period[tau].total_travel_time;
		}

	}
	return total_network_travel_time;
}

// changes here are also for odmes, don't need to implement the changes in this function for now
double g_reset_and_update_link_volume_based_on_ODME_columns(int number_of_links, int iteration_no)
{
    float total_gap = 0;
    float sub_total_gap_link_count = 0;
    float sub_total_gap_P_count = 0;
    float sub_total_gap_A_count = 0;

    // reset the link volume
    for (int i = 0; i < number_of_links; ++i)
    {
        for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
        {
            // used in travel time calculation
            g_link_vector[i].flow_volume_per_period[tau] = 0;
        }
    }

    // reset the estimated production and attraction

            for (int orig = 0; orig < g_zone_vector.size(); ++orig)  // o
        {
            g_zone_vector[orig].est_attraction = 0;
            g_zone_vector[orig].est_production = 0;
        }

    for (int at = 0; at < assignment.g_AgentTypeVector.size(); ++at)  //m
    {
        //#pragma omp parallel for

        std::map<int, CColumnPath>::iterator it;
        int zone_size = g_zone_vector.size();
        int tau_size = assignment.g_DemandPeriodVector.size();

        float link_volume_contributed_by_path_volume;

        int link_seq_no;
        float PCE_ratio;
        int nl;

        std::map<int, CColumnPath>::iterator it_begin;
        std::map<int, CColumnPath>::iterator it_end;

        int column_vector_size;
        CColumnVector* p_column;

        for (int orig = 0; orig < zone_size; ++orig)  // o
        {
            for (int dest = 0; dest < zone_size; ++dest) //d
            {
                for (int tau = 0; tau < tau_size; ++tau)  //tau
                {
                    p_column = &(assignment.g_column_pool[orig][dest][at][tau]);
                    if (p_column->od_volume > 0)
                    {
                        // continuous: type 0
                        column_vector_size = p_column->path_node_sequence_map.size();

                        it_begin = p_column->path_node_sequence_map.begin();
                        it_end = p_column->path_node_sequence_map.end();
                        for (it = it_begin; it != it_end; ++it)  // path k
                        {
                            link_volume_contributed_by_path_volume = it->second.path_volume;  // assign all OD flow to this first path

                            g_zone_vector[orig].est_production += it->second.path_volume;
                            g_zone_vector[dest].est_attraction += it->second.path_volume;

                            // add path volume to link volume
                            for (nl = 0; nl < it->second.m_link_size; ++nl)  // arc a
                            {
                                link_seq_no = it->second.path_link_vector[nl];

                                // MSA updating for the existing column pools
                                // if iteration_index = 0; then update no flow discount is used (for the column pool case)
                                PCE_ratio = 1;
                                //#pragma omp critical
                                {
                                    g_link_vector[link_seq_no].flow_volume_per_period[tau] += link_volume_contributed_by_path_volume * PCE_ratio;
                                    g_link_vector[link_seq_no].volume_per_period_per_at[tau][at] += link_volume_contributed_by_path_volume;  // pure volume, not consider PCE
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    total_gap = 0;
    sub_total_gap_link_count = 0;
    sub_total_gap_P_count = 0;
    sub_total_gap_A_count = 0;
    int total_link_count = 0; 
    int total_prod_count = 0;
    int total_attr_count = 0;

    // calcualte deviation for each measurement type
    for (int i = 0; i < number_of_links; ++i)
    {
        g_link_vector[i].CalculateTD_VDFunction();

        if (g_link_vector[i].obs_count >= 1)  // with data
        {
            int tau = 0;
            g_link_vector[i].est_count_dev = g_link_vector[i].flow_volume_per_period[tau] - g_link_vector[i].obs_count;
                //dtalog.output() << "link " << g_node_vector [g_link_vector[i].from_node_seq_no].node_id
                //                << "->" << g_node_vector[g_link_vector[i].to_node_seq_no].node_id
                //                << " , obs: " << g_link_vector[i].obs_count << " ,est: " << g_link_vector[i].flow_volume_per_period[tau]
                //                << " , dev:" << g_link_vector[i].est_count_dev << endl;
            total_gap += abs(g_link_vector[i].est_count_dev);
            total_link_count += 1;
            sub_total_gap_link_count += abs(g_link_vector[i].est_count_dev)/ max(1,g_link_vector[i].obs_count);
        }
    }

    for (int orig = 0; orig < g_zone_vector.size(); ++orig)  // o
    {
        //if (g_zone_vector[orig].obs_attraction >= 1)  // with observation
        //{
        //    g_zone_vector[orig].est_attraction_dev = g_zone_vector[orig].est_attraction[at] - g_zone_vector[orig].obs_attraction;
        //        //dtalog.output() << "zone " << g_zone_vector[orig].zone_id << "A: obs:" << g_zone_vector[orig].obs_attraction
        //        //                << ",est:," << g_zone_vector[orig].est_attraction << ",dev:," << g_zone_vector[orig].est_attraction_dev << endl;

        //    total_gap += abs(g_zone_vector[orig].est_attraction_dev);
        //    sub_total_gap_A_count += g_zone_vector[orig].est_attraction_dev / g_zone_vector[orig].obs_attraction;
        //    total_attr_count++;
        //}

        //if (g_zone_vector[orig].obs_production >= 1)  // with observation
        //{
        //    g_zone_vector[orig].est_production_dev = g_zone_vector[orig].est_production - g_zone_vector[orig].obs_production;

        //        //dtalog.output() << "zone " << g_zone_vector[orig].zone_id << "P: obs:" << g_zone_vector[orig].obs_production
        //        //                << ",est:," << g_zone_vector[orig].est_production << ",dev:," << g_zone_vector[orig].est_production_dev << endl;

        //    total_gap += abs(g_zone_vector[orig].est_production_dev);
        //    sub_total_gap_P_count += g_zone_vector[orig].est_production_dev / g_zone_vector[orig].obs_production;
        //    total_prod_count++;
        //}
    }

    dtalog.output() << "ODME #" << iteration_no/*<< " total abs gap= " << total_gap*/
                    << " ,%_gap_link: " << fabs(sub_total_gap_link_count)/max(1, total_link_count) *100
                    << " ,%_gap_P: " << sub_total_gap_P_count/ max(1,total_prod_count) *100
                    << " ,%_gap_A: " << sub_total_gap_A_count/max(1, total_attr_count) * 100 << endl;
    double gap = fabs(sub_total_gap_link_count / max(1, total_link_count));
    return gap; 
}

void g_update_gradient_cost_and_assigned_flow_in_column_pool(Assignment& assignment, int inner_iteration_number)
{
    float total_gap = 0;
    float total_relative_gap = 0;
    float total_gap_count = 0;

    // we can have a recursive formulat to reupdate the current link volume by a factor of k/(k+1),
    // and use the newly generated path flow to add the additional 1/(k+1)
    g_reset_and_update_link_volume_based_on_columns(g_link_vector.size(), inner_iteration_number, false);

    // step 4: based on newly calculated path volumn, update volume based travel time, and update volume based resource balance, update gradie
    update_link_travel_time_and_cost();
    // step 0

    //step 1: calculate shortest path at inner iteration of column flow updating
#pragma omp parallel for
    for (int orig = 0; orig < g_zone_vector.size(); ++orig)  // o
    {
        CColumnVector* p_column;
        std::map<int, CColumnPath>::iterator it, it_begin, it_end;
        int column_vector_size;

        float least_gradient_cost = 999999;
        int least_gradient_cost_path_seq_no = -1;
        int least_gradient_cost_path_node_sum_index = -1;
        int path_seq_count = 0;

        double path_toll = 0;
		double path_gradient_cost = 0;
		double path_distance = 0;
		double path_travel_time = 0;
        int link_seq_no;

		double link_travel_time;
		double total_switched_out_path_volume = 0;

		double step_size = 0;
		double previous_path_volume = 0;

        for (int dest = 0; dest < g_zone_vector.size(); ++dest) //d
        {
            for (int at = 0; at < assignment.g_AgentTypeVector.size(); ++at)  //m
            {
                for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); ++tau)  //tau
                {
                    p_column = &(assignment.g_column_pool[orig][dest][at][tau]);
                    if (p_column->od_volume > 0)
                    {
                        column_vector_size = p_column->path_node_sequence_map.size();

                        // scan through the map with different node sum for different paths
                        /// step 1: update gradient cost for each column path
                        //if (o = 7 && d == 15)
                        //{

                        //	if (assignment.g_pFileDebugLog != NULL)
                        //		fprintf(assignment.g_pFileDebugLog, "CU: iteration %d: total_gap=, %f,total_relative_gap,%f,\n", inner_iteration_number, total_gap, total_gap / max(0.00001, total_gap_count));
                        //}
                        least_gradient_cost = 999999;
                        least_gradient_cost_path_seq_no = -1;
                        least_gradient_cost_path_node_sum_index = -1;
                        path_seq_count = 0;

                        it_begin = p_column->path_node_sequence_map.begin();
                        it_end = p_column->path_node_sequence_map.end();
                        for (it = it_begin; it != it_end; ++it)
                        {
                            path_toll = 0;
                            path_gradient_cost = 0;
                            path_distance = 0;
                            path_travel_time = 0;

                            for (int nl = 0; nl < it->second.m_link_size; ++nl)  // arc a
                            {
                                link_seq_no = it->second.path_link_vector[nl];
                                path_toll += g_link_vector[link_seq_no].VDF_period[tau].toll[at];
                                path_distance += g_link_vector[link_seq_no].length;
                                link_travel_time = g_link_vector[link_seq_no].travel_time_per_period[tau];
                                path_travel_time += link_travel_time;

                                path_gradient_cost += g_link_vector[link_seq_no].get_generalized_first_order_gradient_cost_of_second_order_loss_for_agent_type(tau, at);
                            }

                            it->second.path_toll = path_toll;
                            it->second.path_travel_time = path_travel_time;
                            it->second.path_gradient_cost = path_gradient_cost;

                            if (column_vector_size == 1)  // only one path
                            {
                                total_gap_count += (it->second.path_gradient_cost * it->second.path_volume);
                                break;
                            }

                            if (path_gradient_cost < least_gradient_cost)
                            {
                                least_gradient_cost = path_gradient_cost;
                                least_gradient_cost_path_seq_no = it->second.path_seq_no;
                                least_gradient_cost_path_node_sum_index = it->first;
                            }
                        }

                        if (column_vector_size >= 2)
                        {
                            // step 2: calculate gradient cost difference for each column path
                            total_switched_out_path_volume = 0;
                            for (it = it_begin; it != it_end; ++it)
                            {
                                if (it->second.path_seq_no != least_gradient_cost_path_seq_no)  //for non-least cost path
                                {
                                    it->second.path_gradient_cost_difference = it->second.path_gradient_cost - least_gradient_cost;
                                    it->second.path_gradient_cost_relative_difference = it->second.path_gradient_cost_difference / max(0.0001f, least_gradient_cost);

                                    total_gap += (it->second.path_gradient_cost_difference * it->second.path_volume);
                                    total_gap_count += (it->second.path_gradient_cost * it->second.path_volume);

                                    step_size = 1.0 / (inner_iteration_number + 2) * p_column->od_volume;

                                    previous_path_volume = it->second.path_volume;

                                    //recall that it->second.path_gradient_cost_difference >=0
                                    // step 3.1: shift flow from nonshortest path to shortest path
                                    it->second.path_volume = max(0.0, it->second.path_volume - step_size * it->second.path_gradient_cost_relative_difference);

                                    //we use min(step_size to ensure a path is not switching more than 1/n proportion of flow
                                    it->second.path_switch_volume = (previous_path_volume - it->second.path_volume);
                                    total_switched_out_path_volume += (previous_path_volume - it->second.path_volume);
                                }
                            }

                            //step 3.2 consider least cost path, receive all volume shifted from non-shortest path
                            if (least_gradient_cost_path_seq_no != -1)
                            {
                                p_column->path_node_sequence_map[least_gradient_cost_path_node_sum_index].path_volume += total_switched_out_path_volume;
                                total_gap_count += (p_column->path_node_sequence_map[least_gradient_cost_path_node_sum_index].path_gradient_cost *
                                    p_column->path_node_sequence_map[least_gradient_cost_path_node_sum_index].path_volume);
                            }
                        }
                    }
                }
            }
        }
    }
}

void g_column_pool_optimization(Assignment& assignment, int column_updating_iterations)
{
    // column_updating_iterations is internal numbers of column updating
    for (int n = 0; n < column_updating_iterations; ++n)
    {
        dtalog.output() << "Current iteration number: " << n << endl;
        g_update_gradient_cost_and_assigned_flow_in_column_pool(assignment, n);

        if(dtalog.debug_level() >=3)
        {
            for (int i = 0; i < g_link_vector.size(); ++i)
            {
                dtalog.output() << "link: " << g_node_vector[g_link_vector[i].from_node_seq_no].node_id << "-->"
                                << g_node_vector[g_link_vector[i].to_node_seq_no].node_id << ", "
                                << "flow count:" << g_link_vector[i].flow_volume_per_period[0] << endl;
            }
        }
    }
}



void g_output_accessibility()
{

    if (assignment.accessibility_output == 0)
        return; 

    dtalog.output() << "writing accessibility.csv.." << endl;

    float path_time_vector[_MAX_LINK_SIZE_IN_A_PATH];
    FILE* g_pFileODMOE = nullptr;
    fopen_ss(&g_pFileODMOE, "accessibility.csv", "w");

    if (!g_pFileODMOE)
    {
        dtalog.output() << "File accessibility.csv cannot be opened." << endl;
        g_ProgramStop();
    }

    fprintf(g_pFileODMOE, "agent_id,o_zone_id,d_zone_id,path_id,agent_type,demand_period,volume,toll,travel_time,great_circle_distance_in_meter,network_distance,geometry\n");

    int count = 1;

    clock_t start_t, end_t;
    start_t = clock();
    clock_t iteration_t;

    int buffer_len;

    int agent_type_size = assignment.g_AgentTypeVector.size();
    int zone_size = g_zone_vector.size();
    int demand_period_size = assignment.g_DemandPeriodVector.size();

    CColumnVector* p_column;

    float path_toll = 0;
    float path_distance = 0;
    float path_travel_time = 0;
    float time_stamp = 0;

    std::map<int, CColumnPath>::iterator it, it_begin, it_end;

    dtalog.output() << "writing data for " << zone_size << "  zones " << endl;

    for (int orig = 0; orig < zone_size; ++orig)
    {
        if (g_zone_vector[orig].zone_id % 100 == 0)
            dtalog.output() << "o zone id =  " << g_zone_vector[orig].zone_id << endl;

        for (int at = 0; at < agent_type_size; ++at)
        {
            for (int dest = 0; dest < zone_size; ++dest)
            {
                for (int tau = 0; tau < demand_period_size; ++tau)
                {
                    p_column = &(assignment.g_column_pool[orig][dest][at][tau]);
                    if (p_column->od_volume > 0)
                    {
                        time_stamp = (assignment.g_DemandPeriodVector[tau].starting_time_slot_no + assignment.g_DemandPeriodVector[tau].ending_time_slot_no) / 2.0 * MIN_PER_TIMESLOT;

                        // scan through the map with different node sum for different continuous paths
                        it_begin = p_column->path_node_sequence_map.begin();
                        it_end = p_column->path_node_sequence_map.end();

                        for (it = it_begin; it != it_end; ++it)
                        {
                            if (count % 100000 == 0)
                            {
                                end_t = clock();
                                iteration_t = end_t - start_t;
                                dtalog.output() << "writing " << count / 1000 << "K agents with CPU time " << iteration_t / 1000.0 << " s" << endl;
                            }

                            path_toll = 0;
                            path_distance = 0;
                            path_travel_time = 0;
                            path_time_vector[0] = time_stamp;

                            for (int nl = 0; nl < it->second.m_link_size; ++nl)  // arc a
                            {
                                int link_seq_no = it->second.path_link_vector[nl];
                                path_toll += g_link_vector[link_seq_no].VDF_period[tau].toll[at];
                                path_distance += g_link_vector[link_seq_no].length;
                                float link_travel_time = g_link_vector[link_seq_no].travel_time_per_period[tau];
                                path_travel_time += link_travel_time;
                                time_stamp += link_travel_time;
                                path_time_vector[nl + 1] = time_stamp;
                            }

                            int virtual_link_delta = 1;
                            // fixed routes have physical nodes always, without virtual connectors
                            if (p_column->bfixed_route)
                                virtual_link_delta = 0;

                            // assignment_mode = 1, path flow mode
                            if (assignment.assignment_mode == 1 || assignment.assignment_mode == 3)
                            {
                                float great_circule_distance = g_zone_vector[orig].m_ODAccessibilityMatrix[at][tau].distance_map[dest]*1000;
                                buffer_len = 0;
                                buffer_len = sprintf(str_buffer, "%d,%d,%d,%d,%s,%s,%.2f,%.1f,%.4f,%.4f,%.4f,",
                                    count,
                                    g_zone_vector[orig].zone_id,
                                    g_zone_vector[dest].zone_id,
                                    it->second.path_seq_no,
                                    assignment.g_AgentTypeVector[at].agent_type.c_str(),
                                    assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
                                    it->second.path_volume,
                                    path_toll,
                                    path_travel_time,
                                    great_circule_distance,
                                    path_distance);


                                buffer_len += sprintf(str_buffer + buffer_len, "\"LINESTRING (");

                                for (int ni = 0 + virtual_link_delta; ni < it->second.m_node_size - virtual_link_delta; ++ni)
                                {
                                    buffer_len += sprintf(str_buffer + buffer_len, "%f %f", g_node_vector[it->second.path_node_vector[ni]].x,
                                        g_node_vector[it->second.path_node_vector[ni]].y);

                                    if (ni != it->second.m_node_size - virtual_link_delta - 1)
                                        buffer_len += sprintf(str_buffer + buffer_len, ", ");
                                }

                                buffer_len += sprintf(str_buffer + buffer_len, ")\"\n");
                                fprintf(g_pFileODMOE, "%s", str_buffer);
                                count++;
                            }

                        }
                    }
                }
            }
        }
    }
    fclose(g_pFileODMOE);
    /////////////////////////////////////////////////////////////////////////////////////////////////////
}
void g_output_simulation_result(Assignment& assignment)
{
    dtalog.output() << "writing link_performance.csv.." << endl;

    int b_debug_detail_flag = 0;
    FILE* g_pFileLinkMOE = nullptr;

    fopen_ss(&g_pFileLinkMOE,"link_performance.csv", "w");
    if (!g_pFileLinkMOE)
    {
        dtalog.output() << "File link_performance.csv cannot be opened." << endl;
        g_ProgramStop();
    }
    else
    {
        if (assignment.assignment_mode <= 1 || assignment.assignment_mode == 3)  //ODME
        {
            // Option 2: BPR-X function
            fprintf(g_pFileLinkMOE, "link_id,link_type,link_type_name,name,from_node_id,to_node_id,nlanes,length,lane_hour_capacity,time_period,VDF_capacity,volume,VOC,travel_time,speed,queue,density,geometry,");
             //ODME
            fprintf(g_pFileLinkMOE, "tmc_corridor_name,tmc_road_sequence,");

            if (assignment.assignment_mode == 3)
                fprintf(g_pFileLinkMOE, "obs_count,dev");

            fprintf(g_pFileLinkMOE, "notes\n");


                    for (int i = 0; i < g_link_vector.size(); ++i)
                    {
                        // virtual connectors
                        if (g_link_vector[i].link_type == -1)
                            continue;

                        for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
                        {
                            float speed = g_link_vector[i].free_speed;

                            if (g_link_vector[i].VDF_period[tau].avg_travel_time > 0.001f)
                                speed = g_link_vector[i].free_speed / (g_link_vector[i].VDF_period[tau].avg_travel_time / g_link_vector[i].VDF_period[tau].FFTT);
                            fprintf(g_pFileLinkMOE, "%s,%d,%s,%s,%d,%d,%d,%.3f,%.3f,%s,%.3f,%.3f,%.3f,%.3f,%.3f,0,0,\"%s\",",
                                g_link_vector[i].link_id.c_str(),
                                g_link_vector[i].link_type,
                                g_link_vector[i].link_type_name.c_str(),
                                g_link_vector[i].name.c_str(),
                                g_node_vector[g_link_vector[i].from_node_seq_no].node_id,
                                g_node_vector[g_link_vector[i].to_node_seq_no].node_id,
                                g_link_vector[i].number_of_lanes,
                                g_link_vector[i].length,
                                g_link_vector[i].lane_capacity,
                                assignment.g_DemandPeriodVector[tau].time_period.c_str(),
                                g_link_vector[i].VDF_period[tau].capacity,
                                g_link_vector[i].flow_volume_per_period[tau],
                                g_link_vector[i].VDF_period[tau].VOC,
                                g_link_vector[i].VDF_period[tau].avg_travel_time,
                                speed,  /* 60.0 is used to convert min to hour */
                                g_link_vector[i].geometry.c_str());

                            fprintf(g_pFileLinkMOE, "%s,%d,", g_link_vector[i].tmc_corridor_name.c_str(), g_link_vector[i].tmc_road_sequence);


                            if (assignment.assignment_mode == 3)  //ODME
                            {
                                if (g_link_vector[i].obs_count >= 1) //ODME
                                    fprintf(g_pFileLinkMOE, "%.1f,%.1f,", g_link_vector[i].obs_count, g_link_vector[i].est_count_dev);
                                else
                                    fprintf(g_pFileLinkMOE, ",,,");
                            }
                            fprintf(g_pFileLinkMOE, "period-based\n");


                            // print out for BPR-X
                            bool b_print_out_for_BPR_X = false;
                            if (b_print_out_for_BPR_X)
                            {
                                // skip the printout for the nonqueued link or invalid queue data
                                if (g_link_vector[i].VDF_period[tau].t0 == g_link_vector[i].VDF_period[tau].t3 || !g_link_vector[i].VDF_period[tau].bValidQueueData)
                                    continue;

                                int start_time_slot_no = max(g_link_vector[i].VDF_period[tau].starting_time_slot_no, g_link_vector[i].VDF_period[tau].t0);
                                int end_time_slot_no = min(g_link_vector[i].VDF_period[tau].ending_time_slot_no, g_link_vector[i].VDF_period[tau].t3);
                                //tt here is absolute time index
                                for (int tt = start_time_slot_no; tt <= end_time_slot_no; ++tt)
                                {
                                    // 15 min per interval
                                    int time = tt * MIN_PER_TIMESLOT;

                                    float speed = g_link_vector[i].length / (max(0.001f, g_link_vector[i].VDF_period[tau].travel_time[tt]));
                                    float V_mu_over_V_f_ratio = 0.5; // to be calibrated.
                                    float physical_queue = g_link_vector[i].VDF_period[tau].Queue[tt] / (1 - V_mu_over_V_f_ratio);  // per lane
                                    float density = g_link_vector[i].VDF_period[tau].discharge_rate[tt] / max(0.001f, speed);
                                    if (density > 150)  // 150 as kjam.
                                        density = 150;

                                    fprintf(g_pFileLinkMOE, "%s,%d,%d,%s_%s,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,\"%s\",",
                                        g_link_vector[i].link_id.c_str(),
                                        g_node_vector[g_link_vector[i].from_node_seq_no].node_id,
                                        g_node_vector[g_link_vector[i].to_node_seq_no].node_id,
                                        g_time_coding(time).c_str(), g_time_coding(time + MIN_PER_TIMESLOT).c_str(),
                                        g_link_vector[i].VDF_period[tau].discharge_rate[tt],
                                        g_link_vector[i].VDF_period[tau].travel_time[tt] * 60,  /*convert per hour to min*/
                                        speed,
                                        g_link_vector[i].VDF_period[tau].VOC,
                                        physical_queue,
                                        density,
                                        g_link_vector[i].geometry.c_str());

                                    fprintf(g_pFileLinkMOE, "slot-based\n");
                                }
                            }
                        }
                }  // for each link
                    fclose(g_pFileLinkMOE);
        }
        if (assignment.assignment_mode == 2 || assignment.assignment_mode == 3)  // space time based simulation // ODME
        {
            FILE* g_pFileTDLinkMOE = nullptr;
            fopen_ss(&g_pFileTDLinkMOE, "link_TD_performance.csv", "w");

            // Option 2: BPR-X function
            fprintf(g_pFileTDLinkMOE, "link_id,from_node_id,to_node_id,time_period,volume,CA,CD,density,queue,travel_time,waiting_time_in_sec,speed,");
            fprintf(g_pFileTDLinkMOE, "notes\n");

            //Initialization for all nodes
            for (int i = 0; i < g_link_vector.size(); ++i)
            {
                // virtual connectors
                if (g_link_vector[i].link_type == -1)
                    continue;

                // first loop for time t
                for (int t = 0; t < assignment.g_number_of_simulation_intervals; ++t)
                {
                    int sampling_time_interval = 60;
                    if (t % (sampling_time_interval / number_of_seconds_per_interval) == 0)
                    {
                        int time_in_min = t / (60 / number_of_seconds_per_interval);  //relative time

                        float volume = 0;
                        float queue = 0;
                        float waiting_time_in_sec = 0;
                        int arrival_rate = 0;
                        float avg_waiting_time_in_sec = 0;

                        float travel_time = (float)(assignment.m_LinkTDTravelTime[i][t / number_of_interval_per_min]) * number_of_seconds_per_interval / sampling_time_interval + avg_waiting_time_in_sec / 60.0;;
                        float speed = g_link_vector[i].length / (g_link_vector[i].free_flow_travel_time_in_min / 60.0);
                        float virtual_arrival = 0;

                        if (time_in_min >= 1)
                        {
                            volume = assignment.m_LinkCumulativeDeparture[i][t] - assignment.m_LinkCumulativeDeparture[i][t - sampling_time_interval / number_of_seconds_per_interval];

                            if (t - assignment.m_LinkTDTravelTime[i][t / number_of_interval_per_min] >= 0)  // if we have waiting time
                                virtual_arrival = assignment.m_LinkCumulativeArrival[i][t - assignment.m_LinkTDTravelTime[i][t / number_of_interval_per_min]];

                            queue = virtual_arrival - assignment.m_LinkCumulativeDeparture[i][t];
                            //							waiting_time_in_min = queue / (max(1, volume));

                            float waiting_time_count = 0;

                            waiting_time_in_sec = assignment.m_LinkTDWaitingTime[i][t / number_of_interval_per_min] * number_of_seconds_per_interval;

                            arrival_rate = assignment.m_LinkCumulativeArrival[i][t + sampling_time_interval / number_of_seconds_per_interval] - assignment.m_LinkCumulativeArrival[i][t];
                            avg_waiting_time_in_sec = waiting_time_in_sec / max(1, arrival_rate);

                            travel_time = (float)(assignment.m_LinkTDTravelTime[i][t / number_of_interval_per_min]) * number_of_seconds_per_interval / sampling_time_interval + avg_waiting_time_in_sec / 60.0;
                            speed = g_link_vector[i].length / (max(0.00001f, travel_time) / sampling_time_interval);
                        }

                        float density = (assignment.m_LinkCumulativeArrival[i][t] - assignment.m_LinkCumulativeDeparture[i][t]) / (g_link_vector[i].length * g_link_vector[i].number_of_lanes);

                        fprintf(g_pFileTDLinkMOE, "%s,%d,%d,%s_%s,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,",
                                g_link_vector[i].link_id.c_str(),
                                g_node_vector[g_link_vector[i].from_node_seq_no].node_id,
                                g_node_vector[g_link_vector[i].to_node_seq_no].node_id,
                                g_time_coding(assignment.g_LoadingStartTimeInMin + time_in_min).c_str(),
                                g_time_coding(assignment.g_LoadingStartTimeInMin + time_in_min + 1).c_str(),
                                volume,
                                assignment.m_LinkCumulativeArrival[i][t],
                                assignment.m_LinkCumulativeDeparture[i][t],
                                density,
                                queue,
                                travel_time,
                                avg_waiting_time_in_sec,
                                speed);

                        fprintf(g_pFileTDLinkMOE, "simulation-based\n");
                    }
                }  // for each time t
            }  // for each link l
            fclose(g_pFileTDLinkMOE);
        
        }//assignment mode 2 as simulation

        //for (int l = 0; l < g_link_vector.size(); l++) //Initialization for all nodes
        //{
        //		if (g_link_vector[l].link_type == -1)  // virtual connectors
        //	continue;
        //	for (int tau = 0; tau < assignment.g_number_of_demand_periods; tau++)
        //	{
        //
        //			int starting_time = g_link_vector[l].VDF_period[tau].starting_time_slot_no;
        //			int ending_time = g_link_vector[l].VDF_period[tau].ending_time_slot_no;

        //			for (int t = 0; t <= ending_time - starting_time; t++)
        //			{
        //				fprintf(g_pFileTDLinkMOE, "%s,%s,%s,%d,%.3f,%.3f,%.3f,%.3f,%s\n",

        //					g_link_vector[l].link_id.c_str(),
        //					g_node_vector[g_link_vector[l].from_node_seq_no].node_id.c_str(),
        //					g_node_vector[g_link_vector[l].to_node_seq_no].node_id.c_str(),
        //					t,
        //					g_link_vector[l].VDF_period[tau].discharge_rate[t],
        //					g_link_vector[l].VDF_period[tau].travel_time[t],
        //					g_link_vector[l].length / g_link_vector[l].VDF_period[tau].travel_time[t] * 60.0,
        //					g_link_vector[l].VDF_period[tau].congestion_period_P,
        //					"timeslot-dependent");
        //			}

        //		}

        //}


    }
    g_output_accessibility();

    float path_time_vector[_MAX_LINK_SIZE_IN_A_PATH];
    int ODMOE = 1;
    FILE* g_pFileODMOE = nullptr;
    int count = 0;
    if (ODMOE == 1)
    {        //////////////////////////////////////////////////
        dtalog.output() << "writing od_performance.csv.." << endl;

        clock_t start_t, end_t;
        start_t = clock();
        clock_t iteration_t;

        g_pFileODMOE = nullptr;
        fopen_ss(&g_pFileODMOE, "od_performance.csv", "w");

        if (!g_pFileODMOE)
        {
            dtalog.output() << "File od_performance.csv cannot be opened." << endl;
            g_ProgramStop();
        }

        fprintf(g_pFileODMOE, "agent_id,o_zone_id,d_zone_id,agent_type,demand_period,volume,toll,travel_time,distance,geometry\n");


        int buffer_len;

        int agent_type_size = assignment.g_AgentTypeVector.size();
        int zone_size = g_zone_vector.size();
        int demand_period_size = assignment.g_DemandPeriodVector.size();

        CColumnVector* p_column;

        float path_toll = 0;
        float path_distance = 0;
        float path_travel_time = 0;
        float time_stamp = 0;

        std::map<int, CColumnPath>::iterator it, it_begin, it_end;

        dtalog.output() << "writing data for " << zone_size << "  zones " << endl;

        for (int orig = 0; orig < zone_size; ++orig)
        {
            if (g_zone_vector[orig].zone_id % 100 == 0)
                dtalog.output() << "o zone id =  " << g_zone_vector[orig].zone_id << endl;

            for (int at = 0; at < agent_type_size; ++at)
            {
                for (int dest = 0; dest < zone_size; ++dest)
                {
                    for (int tau = 0; tau < demand_period_size; ++tau)
                    {
                        p_column = &(assignment.g_column_pool[orig][dest][at][tau]);
                        if (p_column->od_volume > 0)
                        {
                            time_stamp = (assignment.g_DemandPeriodVector[tau].starting_time_slot_no + assignment.g_DemandPeriodVector[tau].ending_time_slot_no) / 2.0 * MIN_PER_TIMESLOT;

                            // scan through the map with different node sum for different continuous paths
                            it_begin = p_column->path_node_sequence_map.begin();
                            it_end = p_column->path_node_sequence_map.end();

                            float od_volume = 0;
                            float od_toll = 0;
                            float od_distance = 0;
                            float od_travel_time = 0;

                            for (it = it_begin; it != it_end; ++it)
                            {
                                if (count % 100000 == 0)
                                {
                                    end_t = clock();
                                    iteration_t = end_t - start_t;
                                    dtalog.output() << "writing " << count / 1000 << "K agents with CPU time " << iteration_t / 1000.0 << " s" << endl;
                                }

                                path_toll = 0;
                                path_distance = 0;
                                path_travel_time = 0;
                                path_time_vector[0] = time_stamp;

                                for (int nl = 0; nl < it->second.m_link_size; ++nl)  // arc a
                                {
                                    int link_seq_no = it->second.path_link_vector[nl];
                                    path_toll += g_link_vector[link_seq_no].VDF_period[tau].toll[at];
                                    path_distance += g_link_vector[link_seq_no].length;
                                    float link_travel_time = g_link_vector[link_seq_no].travel_time_per_period[tau];
                                    path_travel_time += link_travel_time;
                                    time_stamp += link_travel_time;
                                    path_time_vector[nl + 1] = time_stamp;
                                }

                                od_volume += it->second.path_volume;
                                od_toll += path_toll;
                                od_distance += path_distance;
                                od_travel_time += path_travel_time;

                            }
                            // assignment_mode = 1, path flow mode
                            if (assignment.assignment_mode == 1 || assignment.assignment_mode == 3)
                            {
                                fprintf(g_pFileODMOE, "%d,%d,%d,%s,%s,%.2f,%.1f,%.4f,%.4f,",
                                    count,
                                    g_zone_vector[orig].zone_id,
                                    g_zone_vector[dest].zone_id,
                                    assignment.g_AgentTypeVector[at].agent_type.c_str(),
                                    assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
                                    od_volume,
                                    od_toll / max(0.0001, od_volume),
                                    od_distance / max(0.0001, od_volume),
                                    od_distance / max(0.0001, od_volume));

                                fprintf(g_pFileODMOE, "\"LINESTRING (");

                                fprintf(g_pFileODMOE, "%f %f,", g_zone_vector[orig].cell_x, g_zone_vector[orig].cell_y);
                                fprintf(g_pFileODMOE, "%f %f,", g_zone_vector[dest].cell_x, g_zone_vector[dest].cell_y);
                                fprintf(g_pFileODMOE, ")\"");
                                fprintf(g_pFileODMOE, "\n");

                                count++;
                            }


                        }
                    }
                }
            }
        }
        fclose(g_pFileODMOE);
    }

    if (assignment.assignment_mode == 0)
    {
        fopen_ss(&g_pFileODMOE, "path.csv", "w");
        fclose(g_pFileODMOE);
    }
    else if(assignment.assignment_mode >= 1)
    {
        dtalog.output() << "writing path.csv.." << endl;


        FILE* g_pFileODMOE = nullptr;
        fopen_ss(&g_pFileODMOE,"path.csv", "w");

        if (!g_pFileODMOE)
        {
            dtalog.output() << "File path.csv cannot be opened." << endl;
            g_ProgramStop();
        }

        fprintf(g_pFileODMOE, "agent_id,o_zone_id,d_zone_id,path_id,agent_type,demand_period,volume,toll,travel_time,distance,node_sequence,link_sequence,time_sequence,time_decimal_sequence,link_travel_time_sequence,geometry\n");

        int count = 1;

        clock_t start_t, end_t;
        start_t = clock();
        clock_t iteration_t;

        int buffer_len;

        int agent_type_size = assignment.g_AgentTypeVector.size();
        int zone_size = g_zone_vector.size();
        int demand_period_size = assignment.g_DemandPeriodVector.size();

        CColumnVector* p_column;

        float path_toll = 0;
        float path_distance = 0;
        float path_travel_time = 0;
        float time_stamp = 0;

        std::map<int, CColumnPath>::iterator it, it_begin, it_end;

        dtalog.output() << "writing data for " << zone_size << "  zones " << endl;

        for (int orig = 0; orig < zone_size; ++orig)
        {
            if(g_zone_vector[orig].zone_id % 100 == 0)
                dtalog.output() << "o zone id =  " << g_zone_vector[orig].zone_id << endl;

            for (int at = 0; at < agent_type_size; ++at)
            {
                for (int dest = 0; dest < zone_size; ++dest)
                {
                    for (int tau = 0; tau < demand_period_size; ++tau)
                    {
                        p_column = &(assignment.g_column_pool[orig][dest][at][tau]);
                        if (p_column->od_volume > 0)
                        {
                            time_stamp = (assignment.g_DemandPeriodVector[tau].starting_time_slot_no + assignment.g_DemandPeriodVector[tau].ending_time_slot_no) / 2.0 * MIN_PER_TIMESLOT;

                            // scan through the map with different node sum for different continuous paths
                            it_begin = p_column->path_node_sequence_map.begin();
                            it_end = p_column->path_node_sequence_map.end();

                            for (it = it_begin;it != it_end; ++it)
                            {
                                if (count%100000 ==0)
                                {
                                    end_t = clock();
                                    iteration_t = end_t - start_t;
                                    dtalog.output() << "writing " << count/1000 << "K agents with CPU time " << iteration_t / 1000.0 << " s" << endl;
                                }

                                path_toll = 0;
                                path_distance = 0;
                                path_travel_time = 0;
                                path_time_vector[0] = time_stamp;

                                for (int nl = 0; nl < it->second.m_link_size; ++nl)  // arc a
                                {
                                    int link_seq_no = it->second.path_link_vector[nl];
                                    path_toll += g_link_vector[link_seq_no].VDF_period[tau].toll[at];
                                    path_distance += g_link_vector[link_seq_no].length;
                                    float link_travel_time = g_link_vector[link_seq_no].travel_time_per_period[tau];
                                    path_travel_time += link_travel_time;
                                    time_stamp += link_travel_time;
                                    path_time_vector[nl+1] = time_stamp;
                                }

                                int virtual_link_delta = 1;
                                    // fixed routes have physical nodes always, without virtual connectors
                                if (p_column->bfixed_route)
                                    virtual_link_delta = 0;

                                // assignment_mode = 1, path flow mode
                                if(assignment.assignment_mode == 1 || assignment.assignment_mode == 3)
                                {
                                    buffer_len = 0;
                                    buffer_len = sprintf(str_buffer, "%d,%d,%d,%d,%s,%s,%.2f,%.1f,%.4f,%.4f,",
                                                         count,
                                                         g_zone_vector[orig].zone_id,
                                                         g_zone_vector[dest].zone_id,
                                                         it->second.path_seq_no,
                                                         assignment.g_AgentTypeVector[at].agent_type.c_str(),
                                                         assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
                                                         it->second.path_volume,
                                                         path_toll,
                                                         path_travel_time,
                                                         path_distance);

                                    /* Format and print various data */
                                    for (int ni = 0+ virtual_link_delta; ni <it->second.m_node_size- virtual_link_delta; ++ni)
                                        buffer_len += sprintf(str_buffer + buffer_len, "%d;", g_node_vector[it->second.path_node_vector[ni]].node_id);
                                    buffer_len += sprintf(str_buffer+ buffer_len, ",");

                                    for (int nl = 0 + virtual_link_delta; nl < it->second.m_link_size - virtual_link_delta; ++nl)
                                    {
                                        int link_seq_no = it->second.path_link_vector[nl];
                                        buffer_len += sprintf(str_buffer + buffer_len, "%s;", g_link_vector[link_seq_no].link_id.c_str());
                                    }
                                    buffer_len += sprintf(str_buffer + buffer_len, ",");

                                    for (int nt = 0 + virtual_link_delta; nt < it->second.m_link_size+1 - virtual_link_delta; ++nt)
                                        buffer_len += sprintf(str_buffer + buffer_len, "%s;", g_time_coding(path_time_vector[nt]).c_str());
                                    buffer_len += sprintf(str_buffer + buffer_len, ",");

                                    for (int nt = 0 + virtual_link_delta; nt < it->second.m_link_size+1 - virtual_link_delta; ++nt)
                                        buffer_len += sprintf(str_buffer + buffer_len, "%.2f;", path_time_vector[nt]);
                                    buffer_len += sprintf(str_buffer + buffer_len, ",");

                                    for (int nt = 0 + virtual_link_delta; nt < it->second.m_link_size - virtual_link_delta; ++nt)
                                        buffer_len += sprintf(str_buffer + buffer_len, "%.2f;", path_time_vector[nt+1]- path_time_vector[nt]);
                                    buffer_len += sprintf(str_buffer + buffer_len, ",");

                                    if (buffer_len >= STRING_LENGTH_PER_LINE - 1)
                                    {
                                        dtalog.output() << "Error: buffer_len >= STRING_LENGTH_PER_LINE." << endl;
                                        g_ProgramStop();
                                    }

                                    buffer_len += sprintf(str_buffer + buffer_len,  "\"LINESTRING (");

                                    for (int ni = 0 + virtual_link_delta; ni < it->second.m_node_size - virtual_link_delta; ++ni)
                                    {
                                        buffer_len += sprintf(str_buffer + buffer_len, "%f %f", g_node_vector[it->second.path_node_vector[ni]].x,
                                                              g_node_vector[it->second.path_node_vector[ni]].y);

                                        if (ni != it->second.m_node_size - virtual_link_delta - 1)
                                            buffer_len += sprintf(str_buffer + buffer_len, ", ");
                                    }

                                    buffer_len += sprintf(str_buffer + buffer_len, ")\"\n");
                                    fprintf(g_pFileODMOE, "%s", str_buffer);
                                    count++;
                                }
                              
                            }
                        }
                    }
                }
            }
        }
        fclose(g_pFileODMOE);


/////////////////////////////////////////////////////////////////////////////////////////////////////

        if (assignment.assignment_mode >= 2)  // 2 and 3
        {
        dtalog.output() << "writing trajectory.csv.." << endl;

        float path_time_vector[_MAX_LINK_SIZE_IN_A_PATH];
        FILE* g_pFileTrajectory = nullptr;
        fopen_ss(&g_pFileTrajectory, "trajectory.csv", "w");

        if (!g_pFileTrajectory)
        {
            dtalog.output() << "File trajectory.csv cannot be opened." << endl;
            g_ProgramStop();
        }

        fprintf(g_pFileTrajectory, "agent_id,o_zone_id,d_zone_id,path_id,agent_type,demand_period,volume,toll,travel_time,distance,node_sequence,link_sequence,time_sequence,time_decimal_sequence,link_travel_time_sequence,geometry\n");

        int count = 1;

        clock_t start_t, end_t;
        start_t = clock();
        clock_t iteration_t;

        int buffer_len;

        int agent_type_size = assignment.g_AgentTypeVector.size();
        int zone_size = g_zone_vector.size();
        int demand_period_size = assignment.g_DemandPeriodVector.size();

        CColumnVector* p_column;

        float path_toll = 0;
        float path_distance = 0;
        float path_travel_time = 0;
        float time_stamp = 0;

        std::map<int, CColumnPath>::iterator it, it_begin, it_end;

        dtalog.output() << "writing data for " << zone_size << "  zones " << endl;

        for (int orig = 0; orig < zone_size; ++orig)
        {
            if (g_zone_vector[orig].zone_id % 100 == 0)
                dtalog.output() << "o zone id =  " << g_zone_vector[orig].zone_id << endl;

            for (int at = 0; at < agent_type_size; ++at)
            {
                for (int dest = 0; dest < zone_size; ++dest)
                {
                    for (int tau = 0; tau < demand_period_size; ++tau)
                    {
                        p_column = &(assignment.g_column_pool[orig][dest][at][tau]);
                        if (p_column->od_volume > 0)
                        {
                            time_stamp = (assignment.g_DemandPeriodVector[tau].starting_time_slot_no + assignment.g_DemandPeriodVector[tau].ending_time_slot_no) / 2.0 * MIN_PER_TIMESLOT;

                            // scan through the map with different node sum for different continuous paths
                            it_begin = p_column->path_node_sequence_map.begin();
                            it_end = p_column->path_node_sequence_map.end();

                            for (it = it_begin; it != it_end; ++it)
                            {
                                if (count % 100000 == 0)
                                {
                                    end_t = clock();
                                    iteration_t = end_t - start_t;
                                    dtalog.output() << "writing " << count / 1000 << "K agents with CPU time " << iteration_t / 1000.0 << " s" << endl;
                                }

                                path_toll = 0;
                                path_distance = 0;
                                path_travel_time = 0;
                                path_time_vector[0] = time_stamp;

                                for (int nl = 0; nl < it->second.m_link_size; ++nl)  // arc a
                                {
                                    int link_seq_no = it->second.path_link_vector[nl];
                                    path_toll += g_link_vector[link_seq_no].VDF_period[tau].toll[at];
                                    path_distance += g_link_vector[link_seq_no].length;
                                    float link_travel_time = g_link_vector[link_seq_no].travel_time_per_period[tau];
                                    path_travel_time += link_travel_time;
                                    time_stamp += link_travel_time;
                                    path_time_vector[nl + 1] = time_stamp;
                                }

                                int virtual_link_delta = 1;
                                // fixed routes have physical nodes always, without virtual connectors
                                if (p_column->bfixed_route)
                                    virtual_link_delta = 0;

                              {
                                    // assignment_mode = 2, simulated agent flow mode

                                    for (int vi = 0; vi < it->second.agent_simu_id_vector.size(); ++vi)
                                    {
                                        buffer_len = 0;
                                        // some bugs for output link performances before
                                        buffer_len = sprintf(str_buffer, "%d,%d,%d,%d,%s,%s,1,%.1f,%.4f,%.4f,",
                                            count,
                                            g_zone_vector[orig].zone_id,
                                            g_zone_vector[dest].zone_id,
                                            it->second.path_seq_no,
                                            assignment.g_AgentTypeVector[at].agent_type.c_str(),
                                            assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
                                            path_toll,
                                            path_travel_time,
                                            path_distance);

                                        /* Format and print various data */

                                        for (int ni = 0 + virtual_link_delta; ni < it->second.m_node_size - virtual_link_delta; ++ni)
                                            buffer_len += sprintf(str_buffer + buffer_len, "%d;", g_node_vector[it->second.path_node_vector[ni]].node_id);
                                        buffer_len += sprintf(str_buffer + buffer_len, ",");

                                        for (int nl = 0 + virtual_link_delta; nl < it->second.m_link_size - virtual_link_delta; ++nl)
                                        {
                                            int link_seq_no = it->second.path_link_vector[nl];
                                            buffer_len += sprintf(str_buffer + buffer_len, "%s;", g_link_vector[link_seq_no].link_id.c_str());
                                        }
                                        buffer_len += sprintf(str_buffer + buffer_len, ",");

                                        int agent_simu_id = it->second.agent_simu_id_vector[vi];
                                        CAgent_Simu* pAgentSimu = g_agent_simu_vector[agent_simu_id];
                                        for (int nt = 0 + virtual_link_delta; nt < it->second.m_link_size + 1 - virtual_link_delta; ++nt)
                                        {
                                            float time_in_min = 0;

                                            if (nt < it->second.m_link_size - virtual_link_delta)
                                                time_in_min = assignment.g_LoadingStartTimeInMin + pAgentSimu->m_Veh_LinkArrivalTime_in_simu_interval[nt] * number_of_seconds_per_interval / 60.0;
                                            else
                                                time_in_min = assignment.g_LoadingStartTimeInMin + pAgentSimu->m_Veh_LinkDepartureTime_in_simu_interval[nt - 1] * number_of_seconds_per_interval / 60.0;  // last link in the path

                                            path_time_vector[nt] = time_in_min;
                                        }

                                        for (int nt = 0 + virtual_link_delta; nt < it->second.m_link_size + 1 - virtual_link_delta; ++nt)
                                            buffer_len += sprintf(str_buffer + buffer_len, "%s;", g_time_coding(path_time_vector[nt]).c_str());
                                        buffer_len += sprintf(str_buffer + buffer_len, ",");

                                        for (int nt = 0 + virtual_link_delta; nt < it->second.m_link_size + 1 - virtual_link_delta; ++nt)
                                            buffer_len += sprintf(str_buffer + buffer_len, "%.2f;", path_time_vector[nt]);
                                        buffer_len += sprintf(str_buffer + buffer_len, ",");

                                        for (int nt = 0 + virtual_link_delta; nt < it->second.m_link_size - virtual_link_delta; ++nt)
                                            buffer_len += sprintf(str_buffer + buffer_len, "%.2f;", path_time_vector[nt + 1] - path_time_vector[nt]);
                                        buffer_len += sprintf(str_buffer + buffer_len, "\n");

                                        if (buffer_len >= STRING_LENGTH_PER_LINE - 1)
                                        {
                                            dtalog.output() << "Error: buffer_len >= STRING_LENGTH_PER_LINE." << endl;
                                            g_ProgramStop();
                                        }

                                        fprintf(g_pFileTrajectory, "%s", str_buffer);
                                        count++;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        fclose(g_pFileTrajectory);

       }
    }
}

void g_output_TD_queue_result(Assignment& assignment)
{
    dtalog.output() << "writing link_performance.csv.." << endl;

    int b_debug_detail_flag = 0;
    FILE* g_pFileLinkMOE = nullptr;

    fopen_ss(&g_pFileLinkMOE, "link_performance.csv", "w");
    if (!g_pFileLinkMOE)
    {
        dtalog.output() << "File link_performance.csv cannot be opened." << endl;
        g_ProgramStop();
    }
    else
    {
        if (assignment.assignment_mode <= 1 || assignment.assignment_mode == 3)  //ODME
        {
            // Option 2: BPR-X function
            fprintf(g_pFileLinkMOE, "link_id,link_type,link_type_name,name,from_node_id,to_node_id,nlanes,length,lane_hour_capacity,time_period,VDF_capacity,volume,VOC,travel_time,speed,queue,density,geometry,");
            //ODME
            if (assignment.assignment_mode == 3)
                fprintf(g_pFileLinkMOE, "obs_count,dev");

            fprintf(g_pFileLinkMOE, "notes\n");


            for (int i = 0; i < g_link_vector.size(); ++i)
            {
                // virtual connectors
                if (g_link_vector[i].link_type == -1)
                    continue;

                for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
                {
                    float speed = g_link_vector[i].free_speed;

                    if (g_link_vector[i].VDF_period[tau].avg_travel_time > 0.001f)
                        speed = g_link_vector[i].free_speed / (g_link_vector[i].VDF_period[tau].avg_travel_time / g_link_vector[i].VDF_period[tau].FFTT);
                    fprintf(g_pFileLinkMOE, "%s,%d,%s,%s,%d,%d,%d,%.3f,%.3f,%s,%.3f,%.3f,%.3f,%.3f,%.3f,0,0,\"%s\",",
                        g_link_vector[i].link_id.c_str(),
                        g_link_vector[i].link_type,
                        g_link_vector[i].link_type_name.c_str(),
                        g_link_vector[i].name.c_str(),
                        g_node_vector[g_link_vector[i].from_node_seq_no].node_id,
                        g_node_vector[g_link_vector[i].to_node_seq_no].node_id,
                        g_link_vector[i].number_of_lanes,
                        g_link_vector[i].length,
                        g_link_vector[i].lane_capacity,
                        assignment.g_DemandPeriodVector[tau].time_period.c_str(),
                        g_link_vector[i].VDF_period[tau].capacity,
                        g_link_vector[i].flow_volume_per_period[tau],
                        g_link_vector[i].VDF_period[tau].VOC,
                        g_link_vector[i].VDF_period[tau].avg_travel_time,
                        speed,  /* 60.0 is used to convert min to hour */
                        g_link_vector[i].geometry.c_str());

                    if (assignment.assignment_mode == 3)  //ODME
                    {
                        if (g_link_vector[i].obs_count >= 1) //ODME
                            fprintf(g_pFileLinkMOE, "%.1f,%.1f,", g_link_vector[i].obs_count, g_link_vector[i].est_count_dev);
                        else
                            fprintf(g_pFileLinkMOE, ",,,");
                    }
                    fprintf(g_pFileLinkMOE, "period-based\n");


                    // print out for BPR-X
                    bool b_print_out_for_BPR_X = false;
                    if (b_print_out_for_BPR_X)
                    {
                        // skip the printout for the nonqueued link or invalid queue data
                        if (g_link_vector[i].VDF_period[tau].t0 == g_link_vector[i].VDF_period[tau].t3 || !g_link_vector[i].VDF_period[tau].bValidQueueData)
                            continue;

                        int start_time_slot_no = max(g_link_vector[i].VDF_period[tau].starting_time_slot_no, g_link_vector[i].VDF_period[tau].t0);
                        int end_time_slot_no = min(g_link_vector[i].VDF_period[tau].ending_time_slot_no, g_link_vector[i].VDF_period[tau].t3);
                        //tt here is absolute time index
                        for (int tt = start_time_slot_no; tt <= end_time_slot_no; ++tt)
                        {
                            // 15 min per interval
                            int time = tt * MIN_PER_TIMESLOT;

                            float speed = g_link_vector[i].length / (max(0.001f, g_link_vector[i].VDF_period[tau].travel_time[tt]));
                            float V_mu_over_V_f_ratio = 0.5; // to be calibrated.
                            float physical_queue = g_link_vector[i].VDF_period[tau].Queue[tt] / (1 - V_mu_over_V_f_ratio);  // per lane
                            float density = g_link_vector[i].VDF_period[tau].discharge_rate[tt] / max(0.001f, speed);
                            if (density > 150)  // 150 as kjam.
                                density = 150;

                            fprintf(g_pFileLinkMOE, "%s,%d,%d,%s_%s,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,\"%s\",",
                                g_link_vector[i].link_id.c_str(),
                                g_node_vector[g_link_vector[i].from_node_seq_no].node_id,
                                g_node_vector[g_link_vector[i].to_node_seq_no].node_id,
                                g_time_coding(time).c_str(), g_time_coding(time + MIN_PER_TIMESLOT).c_str(),
                                g_link_vector[i].VDF_period[tau].discharge_rate[tt],
                                g_link_vector[i].VDF_period[tau].travel_time[tt] * 60,  /*convert per hour to min*/
                                speed,
                                g_link_vector[i].VDF_period[tau].VOC,
                                physical_queue,
                                density,
                                g_link_vector[i].geometry.c_str());

                            fprintf(g_pFileLinkMOE, "slot-based\n");
                        }
                    }
                }
            }  // for each link
            fclose(g_pFileLinkMOE);
        }
        if (assignment.assignment_mode == 2 || assignment.assignment_mode == 3)  // space time based simulation // ODME
        {
            FILE* g_pFileTDLinkMOE = nullptr;
            fopen_ss(&g_pFileTDLinkMOE, "link_TD_performance.csv", "w");

            // Option 2: BPR-X function
            fprintf(g_pFileTDLinkMOE, "link_id,from_node_id,to_node_id,time_period,volume,CA,CD,density,queue,travel_time,waiting_time_in_sec,speed,");
            fprintf(g_pFileTDLinkMOE, "notes\n");

            //Initialization for all nodes
            for (int i = 0; i < g_link_vector.size(); ++i)
            {
                // virtual connectors
                if (g_link_vector[i].link_type == -1)
                    continue;

                // first loop for time t
                for (int t = 0; t < assignment.g_number_of_simulation_intervals; ++t)
                {
                    int sampling_time_interval = 60;
                    if (t % (sampling_time_interval / number_of_seconds_per_interval) == 0)
                    {
                        int time_in_min = t / (60 / number_of_seconds_per_interval);  //relative time

                        float volume = 0;
                        float queue = 0;
                        float waiting_time_in_sec = 0;
                        int arrival_rate = 0;
                        float avg_waiting_time_in_sec = 0;

                        float travel_time = (float)(assignment.m_LinkTDTravelTime[i][t / number_of_interval_per_min]) * number_of_seconds_per_interval / sampling_time_interval + avg_waiting_time_in_sec / 60.0;;
                        float speed = g_link_vector[i].length / (g_link_vector[i].free_flow_travel_time_in_min / 60.0);
                        float virtual_arrival = 0;

                        if (time_in_min >= 1)
                        {
                            volume = assignment.m_LinkCumulativeDeparture[i][t] - assignment.m_LinkCumulativeDeparture[i][t - sampling_time_interval / number_of_seconds_per_interval];

                            if (t - assignment.m_LinkTDTravelTime[i][t / number_of_interval_per_min] >= 0)  // if we have waiting time
                                virtual_arrival = assignment.m_LinkCumulativeArrival[i][t - assignment.m_LinkTDTravelTime[i][t / number_of_interval_per_min]];

                            queue = virtual_arrival - assignment.m_LinkCumulativeDeparture[i][t];
                            //							waiting_time_in_min = queue / (max(1, volume));

                            float waiting_time_count = 0;

                            waiting_time_in_sec = assignment.m_LinkTDWaitingTime[i][t / number_of_interval_per_min] * number_of_seconds_per_interval;

                            arrival_rate = assignment.m_LinkCumulativeArrival[i][t + sampling_time_interval / number_of_seconds_per_interval] - assignment.m_LinkCumulativeArrival[i][t];
                            avg_waiting_time_in_sec = waiting_time_in_sec / max(1, arrival_rate);

                            travel_time = (float)(assignment.m_LinkTDTravelTime[i][t / number_of_interval_per_min]) * number_of_seconds_per_interval / sampling_time_interval + avg_waiting_time_in_sec / 60.0;
                            speed = g_link_vector[i].length / (max(0.00001f, travel_time) / sampling_time_interval);
                        }

                        float density = (assignment.m_LinkCumulativeArrival[i][t] - assignment.m_LinkCumulativeDeparture[i][t]) / (g_link_vector[i].length * g_link_vector[i].number_of_lanes);

                        fprintf(g_pFileTDLinkMOE, "%s,%d,%d,%s_%s,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,",
                            g_link_vector[i].link_id.c_str(),
                            g_node_vector[g_link_vector[i].from_node_seq_no].node_id,
                            g_node_vector[g_link_vector[i].to_node_seq_no].node_id,
                            g_time_coding(assignment.g_LoadingStartTimeInMin + time_in_min).c_str(),
                            g_time_coding(assignment.g_LoadingStartTimeInMin + time_in_min + 1).c_str(),
                            volume,
                            assignment.m_LinkCumulativeArrival[i][t],
                            assignment.m_LinkCumulativeDeparture[i][t],
                            density,
                            queue,
                            travel_time,
                            avg_waiting_time_in_sec,
                            speed);

                        fprintf(g_pFileTDLinkMOE, "simulation-based\n");
                    }
                }  // for each time t
            }  // for each link l
            fclose(g_pFileTDLinkMOE);

        }//assignment mode 2 as simulation

        //for (int l = 0; l < g_link_vector.size(); l++) //Initialization for all nodes
        //{
        //		if (g_link_vector[l].link_type == -1)  // virtual connectors
        //	continue;
        //	for (int tau = 0; tau < assignment.g_number_of_demand_periods; tau++)
        //	{
        //
        //			int starting_time = g_link_vector[l].VDF_period[tau].starting_time_slot_no;
        //			int ending_time = g_link_vector[l].VDF_period[tau].ending_time_slot_no;

        //			for (int t = 0; t <= ending_time - starting_time; t++)
        //			{
        //				fprintf(g_pFileTDLinkMOE, "%s,%s,%s,%d,%.3f,%.3f,%.3f,%.3f,%s\n",

        //					g_link_vector[l].link_id.c_str(),
        //					g_node_vector[g_link_vector[l].from_node_seq_no].node_id.c_str(),
        //					g_node_vector[g_link_vector[l].to_node_seq_no].node_id.c_str(),
        //					t,
        //					g_link_vector[l].VDF_period[tau].discharge_rate[t],
        //					g_link_vector[l].VDF_period[tau].travel_time[t],
        //					g_link_vector[l].length / g_link_vector[l].VDF_period[tau].travel_time[t] * 60.0,
        //					g_link_vector[l].VDF_period[tau].congestion_period_P,
        //					"timeslot-dependent");
        //			}

        //		}

        //}


    }
    g_output_accessibility();

    if (assignment.assignment_mode == 0)
    {
        FILE* g_pFileODMOE = nullptr;
        fopen_ss(&g_pFileODMOE, "path.csv", "w");
        fclose(g_pFileODMOE);
    }
    else if (assignment.assignment_mode >= 1)
    {
        dtalog.output() << "writing path.csv.." << endl;

        float path_time_vector[_MAX_LINK_SIZE_IN_A_PATH];
        FILE* g_pFileODMOE = nullptr;
        fopen_ss(&g_pFileODMOE, "path.csv", "w");

        if (!g_pFileODMOE)
        {
            dtalog.output() << "File path.csv cannot be opened." << endl;
            g_ProgramStop();
        }

        fprintf(g_pFileODMOE, "agent_id,o_zone_id,d_zone_id,path_id,agent_type,demand_period,volume,toll,travel_time,distance,node_sequence,link_sequence,time_sequence,time_decimal_sequence,link_travel_time_sequence,geometry\n");

        int count = 1;

        clock_t start_t, end_t;
        start_t = clock();
        clock_t iteration_t;

        int buffer_len;

        int agent_type_size = assignment.g_AgentTypeVector.size();
        int zone_size = g_zone_vector.size();
        int demand_period_size = assignment.g_DemandPeriodVector.size();

        CColumnVector* p_column;

        float path_toll = 0;
        float path_distance = 0;
        float path_travel_time = 0;
        float time_stamp = 0;

        std::map<int, CColumnPath>::iterator it, it_begin, it_end;

        dtalog.output() << "writing data for " << zone_size << "  zones " << endl;

        for (int orig = 0; orig < zone_size; ++orig)
        {
            if (g_zone_vector[orig].zone_id % 100 == 0)
                dtalog.output() << "o zone id =  " << g_zone_vector[orig].zone_id << endl;

            for (int at = 0; at < agent_type_size; ++at)
            {
                for (int dest = 0; dest < zone_size; ++dest)
                {
                    for (int tau = 0; tau < demand_period_size; ++tau)
                    {
                        p_column = &(assignment.g_column_pool[orig][dest][at][tau]);
                        if (p_column->od_volume > 0)
                        {
                            time_stamp = (assignment.g_DemandPeriodVector[tau].starting_time_slot_no + assignment.g_DemandPeriodVector[tau].ending_time_slot_no) / 2.0 * MIN_PER_TIMESLOT;

                            // scan through the map with different node sum for different continuous paths
                            it_begin = p_column->path_node_sequence_map.begin();
                            it_end = p_column->path_node_sequence_map.end();

                            for (it = it_begin; it != it_end; ++it)
                            {
                                if (count % 100000 == 0)
                                {
                                    end_t = clock();
                                    iteration_t = end_t - start_t;
                                    dtalog.output() << "writing " << count / 1000 << "K agents with CPU time " << iteration_t / 1000.0 << " s" << endl;
                                }

                                path_toll = 0;
                                path_distance = 0;
                                path_travel_time = 0;
                                path_time_vector[0] = time_stamp;

                                for (int nl = 0; nl < it->second.m_link_size; ++nl)  // arc a
                                {
                                    int link_seq_no = it->second.path_link_vector[nl];
                                    path_toll += g_link_vector[link_seq_no].VDF_period[tau].toll[at];
                                    path_distance += g_link_vector[link_seq_no].length;
                                    float link_travel_time = g_link_vector[link_seq_no].travel_time_per_period[tau];
                                    path_travel_time += link_travel_time;
                                    time_stamp += link_travel_time;
                                    path_time_vector[nl + 1] = time_stamp;
                                }

                                int virtual_link_delta = 1;
                                // fixed routes have physical nodes always, without virtual connectors
                                if (p_column->bfixed_route)
                                    virtual_link_delta = 0;

                                // assignment_mode = 1, path flow mode
                                if (assignment.assignment_mode == 1 || assignment.assignment_mode == 3)
                                {
                                    buffer_len = 0;
                                    buffer_len = sprintf(str_buffer, "%d,%d,%d,%d,%s,%s,%.2f,%.1f,%.4f,%.4f,",
                                        count,
                                        g_zone_vector[orig].zone_id,
                                        g_zone_vector[dest].zone_id,
                                        it->second.path_seq_no,
                                        assignment.g_AgentTypeVector[at].agent_type.c_str(),
                                        assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
                                        it->second.path_volume,
                                        path_toll,
                                        path_travel_time,
                                        path_distance);

                                    /* Format and print various data */
                                    for (int ni = 0 + virtual_link_delta; ni < it->second.m_node_size - virtual_link_delta; ++ni)
                                        buffer_len += sprintf(str_buffer + buffer_len, "%d;", g_node_vector[it->second.path_node_vector[ni]].node_id);
                                    buffer_len += sprintf(str_buffer + buffer_len, ",");

                                    for (int nl = 0 + virtual_link_delta; nl < it->second.m_link_size - virtual_link_delta; ++nl)
                                    {
                                        int link_seq_no = it->second.path_link_vector[nl];
                                        buffer_len += sprintf(str_buffer + buffer_len, "%s;", g_link_vector[link_seq_no].link_id.c_str());
                                    }
                                    buffer_len += sprintf(str_buffer + buffer_len, ",");

                                    for (int nt = 0 + virtual_link_delta; nt < it->second.m_link_size + 1 - virtual_link_delta; ++nt)
                                        buffer_len += sprintf(str_buffer + buffer_len, "%s;", g_time_coding(path_time_vector[nt]).c_str());
                                    buffer_len += sprintf(str_buffer + buffer_len, ",");

                                    for (int nt = 0 + virtual_link_delta; nt < it->second.m_link_size + 1 - virtual_link_delta; ++nt)
                                        buffer_len += sprintf(str_buffer + buffer_len, "%.2f;", path_time_vector[nt]);
                                    buffer_len += sprintf(str_buffer + buffer_len, ",");

                                    for (int nt = 0 + virtual_link_delta; nt < it->second.m_link_size - virtual_link_delta; ++nt)
                                        buffer_len += sprintf(str_buffer + buffer_len, "%.2f;", path_time_vector[nt + 1] - path_time_vector[nt]);
                                    buffer_len += sprintf(str_buffer + buffer_len, ",");

                                    if (buffer_len >= STRING_LENGTH_PER_LINE - 1)
                                    {
                                        dtalog.output() << "Error: buffer_len >= STRING_LENGTH_PER_LINE." << endl;
                                        g_ProgramStop();
                                    }

                                    buffer_len += sprintf(str_buffer + buffer_len, "\"LINESTRING (");

                                    for (int ni = 0 + virtual_link_delta; ni < it->second.m_node_size - virtual_link_delta; ++ni)
                                    {
                                        buffer_len += sprintf(str_buffer + buffer_len, "%f %f", g_node_vector[it->second.path_node_vector[ni]].x,
                                            g_node_vector[it->second.path_node_vector[ni]].y);

                                        if (ni != it->second.m_node_size - virtual_link_delta - 1)
                                            buffer_len += sprintf(str_buffer + buffer_len, ", ");
                                    }

                                    buffer_len += sprintf(str_buffer + buffer_len, ")\"\n");
                                    fprintf(g_pFileODMOE, "%s", str_buffer);
                                    count++;
                                }

                            }
                        }
                    }
                }
            }
        }
        fclose(g_pFileODMOE);

        int ODMOE = 1;
        if (ODMOE == 1)
        {        //////////////////////////////////////////////////
            dtalog.output() << "writing od_performance.csv.." << endl;


            g_pFileODMOE = nullptr;
            fopen_ss(&g_pFileODMOE, "od_performance.csv", "w");

            if (!g_pFileODMOE)
            {
                dtalog.output() << "File od_performance.csv cannot be opened." << endl;
                g_ProgramStop();
            }

            fprintf(g_pFileODMOE, "agent_id,o_zone_id,d_zone_id,agent_type,demand_period,volume,toll,travel_time,distance,geometry\n");

            count = 1;


            int buffer_len;

            int agent_type_size = assignment.g_AgentTypeVector.size();
            int zone_size = g_zone_vector.size();
            int demand_period_size = assignment.g_DemandPeriodVector.size();

            CColumnVector* p_column;

            float path_toll = 0;
            float path_distance = 0;
            float path_travel_time = 0;
            float time_stamp = 0;

            std::map<int, CColumnPath>::iterator it, it_begin, it_end;

            dtalog.output() << "writing data for " << zone_size << "  zones " << endl;

            for (int orig = 0; orig < zone_size; ++orig)
            {
                if (g_zone_vector[orig].zone_id % 100 == 0)
                    dtalog.output() << "o zone id =  " << g_zone_vector[orig].zone_id << endl;

                for (int at = 0; at < agent_type_size; ++at)
                {
                    for (int dest = 0; dest < zone_size; ++dest)
                    {
                        for (int tau = 0; tau < demand_period_size; ++tau)
                        {
                            p_column = &(assignment.g_column_pool[orig][dest][at][tau]);
                            if (p_column->od_volume > 0)
                            {
                                time_stamp = (assignment.g_DemandPeriodVector[tau].starting_time_slot_no + assignment.g_DemandPeriodVector[tau].ending_time_slot_no) / 2.0 * MIN_PER_TIMESLOT;

                                // scan through the map with different node sum for different continuous paths
                                it_begin = p_column->path_node_sequence_map.begin();
                                it_end = p_column->path_node_sequence_map.end();

                                float od_volume = 0;
                                float od_toll = 0;
                                float od_distance = 0;
                                float od_travel_time = 0;

                                for (it = it_begin; it != it_end; ++it)
                                {
                                    if (count % 100000 == 0)
                                    {
                                        end_t = clock();
                                        iteration_t = end_t - start_t;
                                        dtalog.output() << "writing " << count / 1000 << "K agents with CPU time " << iteration_t / 1000.0 << " s" << endl;
                                    }

                                    path_toll = 0;
                                    path_distance = 0;
                                    path_travel_time = 0;
                                    path_time_vector[0] = time_stamp;

                                    for (int nl = 0; nl < it->second.m_link_size; ++nl)  // arc a
                                    {
                                        int link_seq_no = it->second.path_link_vector[nl];
                                        path_toll += g_link_vector[link_seq_no].VDF_period[tau].toll[at];
                                        path_distance += g_link_vector[link_seq_no].length;
                                        float link_travel_time = g_link_vector[link_seq_no].travel_time_per_period[tau];
                                        path_travel_time += link_travel_time;
                                        time_stamp += link_travel_time;
                                        path_time_vector[nl + 1] = time_stamp;
                                    }

                                    od_volume += it->second.path_volume;
                                    od_toll += path_toll;
                                    od_distance += path_distance;
                                    od_travel_time += path_travel_time;

                                }
                                // assignment_mode = 1, path flow mode
                                if (assignment.assignment_mode == 1 || assignment.assignment_mode == 3)
                                {
                                    fprintf(g_pFileODMOE, "%d,%d,%d,%s,%s,%.2f,%.1f,%.4f,%.4f,",
                                        count,
                                        g_zone_vector[orig].zone_id,
                                        g_zone_vector[dest].zone_id,
                                        assignment.g_AgentTypeVector[at].agent_type.c_str(),
                                        assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
                                        od_volume,
                                        od_toll / max(0.0001, od_volume),
                                        od_distance / max(0.0001, od_volume),
                                        od_distance / max(0.0001, od_volume));

                                    fprintf(g_pFileODMOE, "\"LINESTRING (");

                                    fprintf(g_pFileODMOE, "%f %f,", g_zone_vector[orig].cell_x, g_zone_vector[orig].cell_y);
                                    fprintf(g_pFileODMOE, "%f %f,", g_zone_vector[dest].cell_x, g_zone_vector[dest].cell_y);
                                    fprintf(g_pFileODMOE, ")\"");
                                    fprintf(g_pFileODMOE, "\n");

                                    count++;
                                }


                            }
                        }
                    }
                }
            }
            fclose(g_pFileODMOE);
        }

        /////////////////////////////////////////////////////////////////////////////////////////////////////

        if (assignment.assignment_mode >= 2)  // 2 and 3
        {
            dtalog.output() << "writing trajectory.csv.." << endl;

            float path_time_vector[_MAX_LINK_SIZE_IN_A_PATH];
            FILE* g_pFileTrajectory = nullptr;
            fopen_ss(&g_pFileTrajectory, "trajectory.csv", "w");

            if (!g_pFileTrajectory)
            {
                dtalog.output() << "File trajectory.csv cannot be opened." << endl;
                g_ProgramStop();
            }

            fprintf(g_pFileTrajectory, "agent_id,o_zone_id,d_zone_id,path_id,agent_type,demand_period,volume,toll,travel_time,distance,node_sequence,link_sequence,time_sequence,time_decimal_sequence,link_travel_time_sequence,geometry\n");

            int count = 1;

            clock_t start_t, end_t;
            start_t = clock();
            clock_t iteration_t;

            int buffer_len;

            int agent_type_size = assignment.g_AgentTypeVector.size();
            int zone_size = g_zone_vector.size();
            int demand_period_size = assignment.g_DemandPeriodVector.size();

            CColumnVector* p_column;

            float path_toll = 0;
            float path_distance = 0;
            float path_travel_time = 0;
            float time_stamp = 0;

            std::map<int, CColumnPath>::iterator it, it_begin, it_end;

            dtalog.output() << "writing data for " << zone_size << "  zones " << endl;

            for (int orig = 0; orig < zone_size; ++orig)
            {
                if (g_zone_vector[orig].zone_id % 100 == 0)
                    dtalog.output() << "o zone id =  " << g_zone_vector[orig].zone_id << endl;

                for (int at = 0; at < agent_type_size; ++at)
                {
                    for (int dest = 0; dest < zone_size; ++dest)
                    {
                        for (int tau = 0; tau < demand_period_size; ++tau)
                        {
                            p_column = &(assignment.g_column_pool[orig][dest][at][tau]);
                            if (p_column->od_volume > 0)
                            {
                                time_stamp = (assignment.g_DemandPeriodVector[tau].starting_time_slot_no + assignment.g_DemandPeriodVector[tau].ending_time_slot_no) / 2.0 * MIN_PER_TIMESLOT;

                                // scan through the map with different node sum for different continuous paths
                                it_begin = p_column->path_node_sequence_map.begin();
                                it_end = p_column->path_node_sequence_map.end();

                                for (it = it_begin; it != it_end; ++it)
                                {
                                    if (count % 100000 == 0)
                                    {
                                        end_t = clock();
                                        iteration_t = end_t - start_t;
                                        dtalog.output() << "writing " << count / 1000 << "K agents with CPU time " << iteration_t / 1000.0 << " s" << endl;
                                    }

                                    path_toll = 0;
                                    path_distance = 0;
                                    path_travel_time = 0;
                                    path_time_vector[0] = time_stamp;

                                    for (int nl = 0; nl < it->second.m_link_size; ++nl)  // arc a
                                    {
                                        int link_seq_no = it->second.path_link_vector[nl];
                                        path_toll += g_link_vector[link_seq_no].VDF_period[tau].toll[at];
                                        path_distance += g_link_vector[link_seq_no].length;
                                        float link_travel_time = g_link_vector[link_seq_no].travel_time_per_period[tau];
                                        path_travel_time += link_travel_time;
                                        time_stamp += link_travel_time;
                                        path_time_vector[nl + 1] = time_stamp;
                                    }

                                    int virtual_link_delta = 1;
                                    // fixed routes have physical nodes always, without virtual connectors
                                    if (p_column->bfixed_route)
                                        virtual_link_delta = 0;

                                    {
                                        // assignment_mode = 2, simulated agent flow mode

                                        for (int vi = 0; vi < it->second.agent_simu_id_vector.size(); ++vi)
                                        {
                                            buffer_len = 0;
                                            // some bugs for output link performances before
                                            buffer_len = sprintf(str_buffer, "%d,%d,%d,%d,%s,%s,1,%.1f,%.4f,%.4f,",
                                                count,
                                                g_zone_vector[orig].zone_id,
                                                g_zone_vector[dest].zone_id,
                                                it->second.path_seq_no,
                                                assignment.g_AgentTypeVector[at].agent_type.c_str(),
                                                assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
                                                path_toll,
                                                path_travel_time,
                                                path_distance);

                                            /* Format and print various data */

                                            for (int ni = 0 + virtual_link_delta; ni < it->second.m_node_size - virtual_link_delta; ++ni)
                                                buffer_len += sprintf(str_buffer + buffer_len, "%d;", g_node_vector[it->second.path_node_vector[ni]].node_id);
                                            buffer_len += sprintf(str_buffer + buffer_len, ",");

                                            for (int nl = 0 + virtual_link_delta; nl < it->second.m_link_size - virtual_link_delta; ++nl)
                                            {
                                                int link_seq_no = it->second.path_link_vector[nl];
                                                buffer_len += sprintf(str_buffer + buffer_len, "%s;", g_link_vector[link_seq_no].link_id.c_str());
                                            }
                                            buffer_len += sprintf(str_buffer + buffer_len, ",");

                                            int agent_simu_id = it->second.agent_simu_id_vector[vi];
                                            CAgent_Simu* pAgentSimu = g_agent_simu_vector[agent_simu_id];
                                            for (int nt = 0 + virtual_link_delta; nt < it->second.m_link_size + 1 - virtual_link_delta; ++nt)
                                            {
                                                float time_in_min = 0;

                                                if (nt < it->second.m_link_size - virtual_link_delta)
                                                    time_in_min = assignment.g_LoadingStartTimeInMin + pAgentSimu->m_Veh_LinkArrivalTime_in_simu_interval[nt] * number_of_seconds_per_interval / 60.0;
                                                else
                                                    time_in_min = assignment.g_LoadingStartTimeInMin + pAgentSimu->m_Veh_LinkDepartureTime_in_simu_interval[nt - 1] * number_of_seconds_per_interval / 60.0;  // last link in the path

                                                path_time_vector[nt] = time_in_min;
                                            }

                                            for (int nt = 0 + virtual_link_delta; nt < it->second.m_link_size + 1 - virtual_link_delta; ++nt)
                                                buffer_len += sprintf(str_buffer + buffer_len, "%s;", g_time_coding(path_time_vector[nt]).c_str());
                                            buffer_len += sprintf(str_buffer + buffer_len, ",");

                                            for (int nt = 0 + virtual_link_delta; nt < it->second.m_link_size + 1 - virtual_link_delta; ++nt)
                                                buffer_len += sprintf(str_buffer + buffer_len, "%.2f;", path_time_vector[nt]);
                                            buffer_len += sprintf(str_buffer + buffer_len, ",");

                                            for (int nt = 0 + virtual_link_delta; nt < it->second.m_link_size - virtual_link_delta; ++nt)
                                                buffer_len += sprintf(str_buffer + buffer_len, "%.2f;", path_time_vector[nt + 1] - path_time_vector[nt]);
                                            buffer_len += sprintf(str_buffer + buffer_len, "\n");

                                            if (buffer_len >= STRING_LENGTH_PER_LINE - 1)
                                            {
                                                dtalog.output() << "Error: buffer_len >= STRING_LENGTH_PER_LINE." << endl;
                                                g_ProgramStop();
                                            }

                                            fprintf(g_pFileTrajectory, "%s", str_buffer);
                                            count++;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            fclose(g_pFileTrajectory);

        }
    }
}
void g_output_simulation_result_for_signal_api(Assignment& assignment)
{
    bool movement_str_flag = false;
    //Initialization for all nodes
    for (int i = 0; i < g_link_vector.size(); ++i)
    {
        if (g_link_vector[i].movement_str.length() >= 1)
            movement_str_flag = true;
    }

    // no movement_str data
    if (!movement_str_flag)
        return;

    dtalog.output() << "writing link_performance_sig.csv.." << endl;

    int b_debug_detail_flag = 0;
    FILE* g_pFileLinkMOE = nullptr;

    fopen_ss(&g_pFileLinkMOE, "link_performance_sig.csv", "w");
    if (!g_pFileLinkMOE)
    {
        dtalog.output() << "File link_performance_sig.csv cannot be opened." << endl;
        g_ProgramStop();
    }

    if (assignment.assignment_mode <= 1)
    {
        // Option 2: BPR-X function
        fprintf(g_pFileLinkMOE, "link_id,from_node_id,to_node_id,demand_period,time_period,movement_str,main_node_id,volume,travel_time,speed,VOC,");

        fprintf(g_pFileLinkMOE, "notes\n");

        //Initialization for all nodes
        for (int i = 0; i < g_link_vector.size(); ++i)
        {
            for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
            {
                if (g_link_vector[i].movement_str.length() >= 1)
                {
                    float speed = g_link_vector[i].length / (max(0.001, g_link_vector[i].VDF_period[tau].avg_travel_time) / 60.0);
                    fprintf(g_pFileLinkMOE, "%s,%d,%d,%s,%s,%s,%d,%.3f,%.3f,%.3f,%.3f,",
                        g_link_vector[i].link_id.c_str(),
                        g_node_vector[g_link_vector[i].from_node_seq_no].node_id,
                        g_node_vector[g_link_vector[i].to_node_seq_no].node_id,
                        assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
                        assignment.g_DemandPeriodVector[tau].time_period.c_str(),
                        g_link_vector[i].movement_str.c_str(),
                        g_link_vector[i].main_node_id,
                        g_link_vector[i].flow_volume_per_period[tau],
                        g_link_vector[i].VDF_period[tau].avg_travel_time,
                        speed,  /* 60.0 is used to convert min to hour */
                        g_link_vector[i].VDF_period[tau].VOC);

                    fprintf(g_pFileLinkMOE, "period-based\n");
                }
            }
        }
    }

    if (assignment.assignment_mode == 2 || assignment.assignment_mode == 3)  // space time based simulation
    {
        // Option 2: BPR-X function
        fprintf(g_pFileLinkMOE, "link_id,from_node_id,to_node_id,time_period,volume,CA,CD,queue,travel_time,waiting_time_in_sec,speed,");
        fprintf(g_pFileLinkMOE, "notes\n");

        for (int i = 0; i < g_link_vector.size(); ++i) //Initialization for all nodes
        {
            // Peiheng, 02/02/21, useless block
            if (g_link_vector[i].link_id == "146")
                int idebug = 1;

            for (int t = 0; t < assignment.g_number_of_simulation_intervals; ++t)  // first loop for time t
            {
                if (t % (60 / number_of_seconds_per_interval) == 0)
                {
                    int time_in_min = t / (60 / number_of_seconds_per_interval);  //relative time

                    float volume = 0;
                    float queue = 0;
                    float waiting_time_in_sec = 0;
                    int arrival_rate = 0;
                    float avg_waiting_time_in_sec = 0;

                    float travel_time = 0;
                    float speed = g_link_vector[i].length / (g_link_vector[i].free_flow_travel_time_in_min / 60.0);
                    float virtual_arrival = 0;

                    if (time_in_min >= 1)
                    {
                        volume = assignment.m_LinkCumulativeDeparture[i][t] - assignment.m_LinkCumulativeDeparture[i][t - 60 / number_of_seconds_per_interval];

                        if (t - assignment.m_LinkTDTravelTime[i][t/ number_of_interval_per_min] >= 0)
                            virtual_arrival = assignment.m_LinkCumulativeArrival[i][t - assignment.m_LinkTDTravelTime[i][t/ number_of_interval_per_min]];

                        queue = virtual_arrival - assignment.m_LinkCumulativeDeparture[i][t];
                        // waiting_time_in_min = queue / (max(1, volume));

                        float waiting_time_count = 0;
                        for (int tt = t; tt < t + 60 / number_of_seconds_per_interval; ++tt)
                            waiting_time_count += assignment.m_LinkTDWaitingTime[i][tt / number_of_interval_per_min];

                        if (waiting_time_count >= 1)
                        {
                            waiting_time_in_sec = waiting_time_count * number_of_seconds_per_interval;

                            arrival_rate = assignment.m_LinkCumulativeArrival[i][t + 60 / number_of_seconds_per_interval] - assignment.m_LinkCumulativeArrival[i][t];
                            avg_waiting_time_in_sec = waiting_time_in_sec / max(1, arrival_rate);
                        }
                        else
                            avg_waiting_time_in_sec = 0;

                        travel_time = (float)(assignment.m_LinkTDTravelTime[i][t/ number_of_interval_per_min]) * number_of_seconds_per_interval / 60.0 + avg_waiting_time_in_sec / 60.0;
                        speed = g_link_vector[i].length / (max(0.00001f, travel_time) / 60.0);
                    }

                    fprintf(g_pFileLinkMOE, "%s,%d,%d,%s_%s,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,",
                            g_link_vector[i].link_id.c_str(),
                            g_node_vector[g_link_vector[i].from_node_seq_no].node_id,
                            g_node_vector[g_link_vector[i].to_node_seq_no].node_id,
                            g_time_coding(assignment.g_LoadingStartTimeInMin + time_in_min).c_str(),
                            g_time_coding(assignment.g_LoadingStartTimeInMin + time_in_min + 1).c_str(),
                            volume,
                            assignment.m_LinkCumulativeArrival[i][t],
                            assignment.m_LinkCumulativeDeparture[i][t],
                            queue,
                            travel_time,
                            avg_waiting_time_in_sec,
                            speed);

                    fprintf(g_pFileLinkMOE, "simulation-based\n");
                }
            }  // for each time t
        }  // for each link l
    }//assignment mode 2 as simulation

    //for (int l = 0; l < g_link_vector.size(); l++) //Initialization for all nodes
    //{
    //	for (int tau = 0; tau < assignment.g_number_of_demand_periods; tau++)
    //	{
    //
    //			int starting_time = g_link_vector[l].VDF_period[tau].starting_time_slot_no;
    //			int ending_time = g_link_vector[l].VDF_period[tau].ending_time_slot_no;

    //			for (int t = 0; t <= ending_time - starting_time; t++)
    //			{
    //				fprintf(g_pFileLinkMOE, "%s,%s,%s,%d,%.3f,%.3f,%.3f,%.3f,%s\n",

    //					g_link_vector[l].link_id.c_str(),
    //					g_node_vector[g_link_vector[l].from_node_seq_no].node_id.c_str(),
    //					g_node_vector[g_link_vector[l].to_node_seq_no].node_id.c_str(),
    //					t,
    //					g_link_vector[l].VDF_period[tau].discharge_rate[t],
    //					g_link_vector[l].VDF_period[tau].travel_time[t],
    //					g_link_vector[l].length / g_link_vector[l].VDF_period[tau].travel_time[t] * 60.0,
    //					g_link_vector[l].VDF_period[tau].congestion_period_P,
    //					"timeslot-dependent");
    //			}

    //		}

    //}

    fclose(g_pFileLinkMOE);
}

//***
// major function 1:  allocate memory and initialize the data
// void AllocateMemory(int number_of_nodes)
//
//major function 2: // time-dependent label correcting algorithm with double queue implementation
//int optimal_label_correcting(int origin_node, int destination_node, int departure_time, int g_debugging_flag, FILE* g_pFileDebugLog, Assignment& assignment, int time_period_no = 0, int agent_type = 1, float VOT = 10)

//	//major function: update the cost for each node at each SP tree, using a stack from the origin structure
//int tree_cost_updating(int origin_node, int departure_time, int g_debugging_flag, FILE* g_pFileDebugLog, Assignment& assignment, int time_period_no = 0, int agent_type = 1)

//***

// The one and only application object

int g_number_of_CPU_threads()
{
    int number_of_threads = omp_get_max_threads();
    int max_number_of_threads = 4000;

    if (number_of_threads > max_number_of_threads)
        number_of_threads = max_number_of_threads;

    return number_of_threads;
}

void g_assign_computing_tasks_to_memory_blocks(Assignment& assignment)
{
    //fprintf(g_pFileDebugLog, "-------g_assign_computing_tasks_to_memory_blocks-------\n");
    // step 2: assign node to thread
    dtalog.output() << "Step 2: Assigning computing tasks to memory blocks..." << endl;

    NetworkForSP* PointerMatrx[_MAX_AGNETTYPES][_MAX_TIMEPERIODS][_MAX_MEMORY_BLOCKS];

    for (int at = 0; at < assignment.g_AgentTypeVector.size(); ++at)
    {
        for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); ++tau)
        {
            //assign all nodes to the corresponding thread
            for (int z = 0; z < g_zone_vector.size(); ++z)
            {
                //fprintf(g_pFileDebugLog, "%f\n",g_origin_demand_array[zone_seq_no][at][tau]);
                if (z < assignment.g_number_of_memory_blocks)
                {
                    NetworkForSP* p_NetworkForSP = new NetworkForSP();

                    p_NetworkForSP->m_origin_node_vector.push_back(g_zone_vector[z].node_seq_no);
                    p_NetworkForSP->m_origin_zone_seq_no_vector.push_back(z);

                    p_NetworkForSP->m_agent_type_no = at;
                    p_NetworkForSP->m_tau = tau;
                    p_NetworkForSP->AllocateMemory(assignment.g_number_of_nodes, assignment.g_number_of_links);

                    PointerMatrx[at][tau][z] = p_NetworkForSP;

                    g_NetworkForSP_vector.push_back(p_NetworkForSP);
                }
                else  // zone seq no is greater than g_number_of_memory_blocks
                {
                    //get the corresponding memory block seq no
                    // take residual of memory block size to map a zone no to a memory block no.
                    int memory_block_no = z % assignment.g_number_of_memory_blocks;
                    NetworkForSP* p_NetworkForSP = PointerMatrx[at][tau][memory_block_no];
                    p_NetworkForSP->m_origin_node_vector.push_back(g_zone_vector[z].node_seq_no);
                    p_NetworkForSP->m_origin_zone_seq_no_vector.push_back(z);
                }
            }
        }
    }

    dtalog.output() << "There are " << g_NetworkForSP_vector.size() << " networks in memory." << endl;
}

void g_deallocate_computing_tasks_from_memory_blocks()
{
    //fprintf(g_pFileDebugLog, "-------g_assign_computing_tasks_to_memory_blocks-------\n");
    // step 2: assign node to thread
    for (int n = 0; n < g_NetworkForSP_vector.size(); ++n)
    {
        NetworkForSP* p_NetworkForSP = g_NetworkForSP_vector[n];
        delete p_NetworkForSP;
    }
}

//void g_reset_link_volume_for_all_processors()
//{
//#pragma omp parallel for
//    for (int ProcessID = 0; ProcessID < g_NetworkForSP_vector.size(); ++ProcessID)
//    {
//        NetworkForSP* pNetwork = g_NetworkForSP_vector[ProcessID];
//        //Initialization for all non-origin nodes
//        int number_of_links = assignment.g_number_of_links;
//        for (int i = 0; i < number_of_links; ++i)
//            pNetwork->m_link_flow_volume_array[i] = 0;
//    }
//}


void g_reset_link_volume_for_all_processors()
{
	int number_of_memory_blocks = min((int)g_NetworkForSP_vector.size(), assignment.g_number_of_memory_blocks);

#pragma omp parallel for
	for (int ProcessID = 0; ProcessID < number_of_memory_blocks; ++ProcessID)
	{
		for (int blk = 0; blk < assignment.g_AgentTypeVector.size()*assignment.g_DemandPeriodVector.size(); ++blk)
		{
			NetworkForSP* pNetwork = g_NetworkForSP_vector[blk*assignment.g_number_of_memory_blocks + ProcessID];
			//Initialization for all non-origin nodes
			int number_of_links = assignment.g_number_of_links;
			for (int i = 0; i < number_of_links; ++i)
				pNetwork->m_link_flow_volume_array[i] = 0;
		}

	}
}


void g_fetch_link_volume_for_all_processors()
{
    for (int ProcessID = 0; ProcessID < g_NetworkForSP_vector.size(); ++ProcessID)
    {
        NetworkForSP* pNetwork = g_NetworkForSP_vector[ProcessID];

        for (int i = 0; i < g_link_vector.size(); ++i)
            g_link_vector[i].flow_volume_per_period[pNetwork->m_tau] += pNetwork->m_link_flow_volume_array[i];
    }
    // step 1: travel time based on VDF
}

//major function: update the cost for each node at each SP tree, using a stack from the origin structure
double NetworkForSP::backtrace_shortest_path_tree(Assignment& assignment, int iteration_number_outterloop, int o_node_index)
{
	double total_origin_least_cost = 0;
    int origin_node = m_origin_node_vector[o_node_index]; // assigned no
    int m_origin_zone_seq_no = m_origin_zone_seq_no_vector[o_node_index]; // assigned no

    //if (assignment.g_number_of_nodes >= 100000 && m_origin_zone_seq_no % 100 == 0)
    //{
    //	g_fout << "backtracing for zone " << m_origin_zone_seq_no << endl;
    //}

    int departure_time = m_tau;
    int agent_type = m_agent_type_no;

    if (g_node_vector[origin_node].m_outgoing_link_seq_no_vector.size() == 0)
        return 0;

    // given,  m_node_label_cost[i]; is the gradient cost , for this tree's, from its origin to the destination node i'.

    //	fprintf(g_pFileDebugLog, "------------START: origin:  %d  ; Departure time: %d  ; demand type:  %d  --------------\n", origin_node + 1, departure_time, agent_type);
    float k_path_prob = 1;
    k_path_prob = float(1) / float(iteration_number_outterloop + 1);  //XZ: use default value as MSA, this is equivalent to 1/(n+1)
    // MSA to distribute the continuous flow
    // to do, this is for each nth tree.

    //change of path flow is a function of cost gap (the updated node cost for the path generated at the previous iteration -m_node_label_cost[i] at the current iteration)
    // current path flow - change of path flow,
    // new propability for flow staying at this path
    // for current shortest path, collect all the switched path from the other shortest paths for this ODK pair.
    // check demand flow balance constraints

    int num = 0;
    int number_of_nodes = assignment.g_number_of_nodes;
    int number_of_links = assignment.g_number_of_links;
    int l_node_size = 0;  // initialize local node size index
    int l_link_size = 0;
    int node_sum = 0;

    float path_travel_time = 0;
    float path_distance = 0;

    int current_node_seq_no = -1;  // destination node
    int current_link_seq_no = -1;
    int destination_zone_seq_no;
    double ODvolume, volume;
    CColumnVector* pColumnVector;

    for (int i = 0; i < number_of_nodes; ++i)
    {
        if (g_node_vector[i].zone_id >= 1)
        {
            if (i == origin_node) // no within zone demand
                continue;
            //fprintf(g_pFileDebugLog, "--------origin  %d ; destination node: %d ; (zone: %d) -------\n", origin_node + 1, i+1, g_node_vector[i].zone_id);
            //fprintf(g_pFileDebugLog, "--------iteration number outterloop  %d ;  -------\n", iteration_number_outterloop);
            destination_zone_seq_no = assignment.g_zoneid_to_zone_seq_no_mapping[g_node_vector[i].zone_id];

            pColumnVector = &(assignment.g_column_pool[m_origin_zone_seq_no][destination_zone_seq_no][agent_type][m_tau]);

            if (pColumnVector->bfixed_route) // with routing policy, no need to run MSA for adding new columns
                continue;

            ODvolume = pColumnVector->od_volume;
            volume = ODvolume * k_path_prob;
            // this is contributed path volume from OD flow (O, D, k, per time period

            if (ODvolume > 0.000001)
            {
                l_node_size = 0;  // initialize local node size index
                l_link_size = 0;
                node_sum = 0;

                path_travel_time = 0;
                path_distance = 0;

                current_node_seq_no = i;  // destination node
                current_link_seq_no = -1;

				total_origin_least_cost += ODvolume * m_node_label_cost[current_node_seq_no];
                // backtrace the sp tree from the destination to the root (at origin)
                while (current_node_seq_no >= 0 && current_node_seq_no < number_of_nodes)
                {
                    temp_path_node_vector[l_node_size++] = current_node_seq_no;

                    if (l_node_size >= temp_path_node_vector_size)
                    {
                        dtalog.output() << "Error: l_node_size >= temp_path_node_vector_size" << endl;
                        g_ProgramStop();
                    }

                    // this is valid node
                    if (current_node_seq_no >= 0 && current_node_seq_no < number_of_nodes)
                    {
                        node_sum += current_node_seq_no;
                        current_link_seq_no = m_link_predecessor[current_node_seq_no];

                        // fetch m_link_predecessor to mark the link volume
                        if (current_link_seq_no >= 0 && current_link_seq_no < number_of_links)
                        {
                            temp_path_link_vector[l_link_size++] = current_link_seq_no;

                            // pure link based computing mode
                            if (assignment.assignment_mode == 0)
                            {
                                // this is critical for parallel computing as we can write the volume to data
                                m_link_flow_volume_array[current_link_seq_no]+= volume;
                            }

                            //path_travel_time += g_link_vector[current_link_seq_no].travel_time_per_period[tau];
                            //path_distance += g_link_vector[current_link_seq_no].length;
                        }
                    }
                    current_node_seq_no = m_node_predecessor[current_node_seq_no];  // update node seq no
                }
                //fprintf(g_pFileDebugLog, "\n");

                // we obtain the cost, time, distance from the last tree-k
                if(assignment.assignment_mode >=1) // column based mode
                {
                    // we cannot find a path with the same node sum, so we need to add this path into the map,
                    if (pColumnVector->path_node_sequence_map.find(node_sum) == assignment.g_column_pool[m_origin_zone_seq_no][destination_zone_seq_no][agent_type][m_tau].path_node_sequence_map.end())
                    {
                        // add this unique path
                        int path_count = pColumnVector->path_node_sequence_map.size();
                        pColumnVector->path_node_sequence_map[node_sum].path_seq_no = path_count;
                        pColumnVector->path_node_sequence_map[node_sum].path_volume = 0;
                        //assignment.g_column_pool[m_origin_zone_seq_no][destination_zone_seq_no][agent_type][tau].time = m_label_time_array[i];
                        //assignment.g_column_pool[m_origin_zone_seq_no][destination_zone_seq_no][agent_type][tau].path_node_sequence_map[node_sum].path_distance = m_label_distance_array[i];
                        pColumnVector->path_node_sequence_map[node_sum].path_toll = m_node_label_cost[i];

                        pColumnVector->path_node_sequence_map[node_sum].AllocateVector(
                            l_node_size,
                            temp_path_node_vector,
                            l_link_size,
                            temp_path_link_vector,true);
                    }

                    pColumnVector->path_node_sequence_map[node_sum].path_volume += volume;
                }
            }
        }
    }
	return total_origin_least_cost;
}

void  CLink::CalculateTD_VDFunction()
{
    for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
    {
        float starting_time_slot_no = assignment.g_DemandPeriodVector[tau].starting_time_slot_no;
        float ending_time_slot_no = assignment.g_DemandPeriodVector[tau].ending_time_slot_no;

        // signalized with red_time > 1
        if (this->movement_str.length() > 1 && VDF_period[tau].red_time > 1)
        {
            // arterial streets with the data from sigal API
            float hourly_per_lane_volume = flow_volume_per_period[tau] / (max(1.0f, (ending_time_slot_no - starting_time_slot_no)) * 15 / 60 / number_of_lanes);
            float red_time = VDF_period[tau].red_time;
            float cycle_length = VDF_period[tau].cycle_length;
            //we dynamically update cycle length, and green time/red time, so we have dynamically allocated capacity and average delay
            travel_time_per_period[tau] = VDF_period[tau].PerformSignalVDF(hourly_per_lane_volume, red_time, cycle_length);
            travel_time_per_period[tau] += VDF_period[tau].PerformBPR(flow_volume_per_period[tau]);
            // update capacity using the effective discharge rates, will be passed in to the following BPR function
            VDF_period[tau].capacity = (1 - red_time / cycle_length) * _default_saturation_flow_rate * number_of_lanes;
        }

        // either non-signalized or signalized with red_time < 1 and cycle_length < 30
        if (this->movement_str.length() == 0 || (this->movement_str.length() > 1 && VDF_period[tau].red_time < 1 && VDF_period[tau].cycle_length < 30))
        {
            travel_time_per_period[tau] = VDF_period[tau].PerformBPR(flow_volume_per_period[tau]);
            VDF_period[tau].PerformBPR_X(flow_volume_per_period[tau]);  // only for freeway segments
        }
    }
}

double network_assignment(int assignment_mode, int iteration_number, int column_updating_iterations, int number_of_memory_blocks)
{
    int signal_updating_iterations = 0;

    // k iterations for column generation
    assignment.g_number_of_column_generation_iterations = iteration_number;
    // 0: link UE: 1: path UE, 2: Path SO, 3: path resource constraints
    assignment.assignment_mode = assignment_mode;

    assignment.g_number_of_memory_blocks = min(max(1, number_of_memory_blocks), _MAX_MEMORY_BLOCKS);


    if (assignment.assignment_mode == 0)
        column_updating_iterations = 0;

    // step 1: read input data of network / demand tables / Toll
    g_ReadInputData(assignment);
    g_reload_service_arc_data(assignment);
    g_ReadDemandFileBasedOnDemandFileList(assignment);

    //step 2: allocate memory and assign computing tasks
    g_assign_computing_tasks_to_memory_blocks(assignment); // static cost based label correcting

    // definte timestamps
    clock_t start_t, end_t, total_t;
    start_t = clock();
    clock_t iteration_t, cumulative_lc, cumulative_cp, cumulative_lu;

    //step 3: column generation stage: find shortest path and update path cost of tree using TD travel time
    dtalog.output() << endl;
    dtalog.output() << "Step 3: Column Generation for Traffic Assignment..." << endl;
    dtalog.output() << "Total Column Generation iteration: " << assignment.g_number_of_column_generation_iterations << endl;
    for (int iteration_number = 0; iteration_number < assignment.g_number_of_column_generation_iterations; iteration_number++)
    {
        dtalog.output() << endl;
        dtalog.output() << "Current iteration number:" << iteration_number << endl;
        end_t = clock();
        iteration_t = end_t - start_t;
        dtalog.output() << "Current CPU time: " << iteration_t / 1000.0 << " s" << endl;

        // commment out for DLL version
        // if (signal_updating_iterations >=1 && iteration_number >= signal_updating_iterations)
        // {
        //     g_fout << "use SignalAPI to recalibrate signal timing at iteration " << iteration_number << endl;
        //     SignalAPI(iteration_number, assignment_mode, 0);
        //     g_reload_service_arc_data(assignment);
        // }

        // step 3.1 update travel time and resource consumption
        clock_t start_t_lu = clock();

        double total_system_travel_time = 0;
        double total_least_system_travel_time = 0;
        // initialization at beginning of shortest path
        total_system_travel_time = update_link_travel_time_and_cost();

        if (assignment.assignment_mode == 0)
        {
            //fw
            g_reset_link_volume_in_master_program_without_columns(g_link_vector.size(), iteration_number, true);
            g_reset_link_volume_for_all_processors();
        }
        else
        {
            // we can have a recursive formulat to reupdate the current link volume by a factor of k/(k+1),
            //  and use the newly generated path flow to add the additional 1/(k+1)
            g_reset_and_update_link_volume_based_on_columns(g_link_vector.size(), iteration_number, true);
        }

        if (dtalog.debug_level() >= 3)
        {
            dtalog.output() << "Results:" << endl;
            for (int i = 0; i < g_link_vector.size(); ++i) {
                dtalog.output() << "link: " << g_node_vector[g_link_vector[i].from_node_seq_no].node_id << "-->"
                    << g_node_vector[g_link_vector[i].to_node_seq_no].node_id << ", "
                    << "flow count:" << g_link_vector[i].flow_volume_per_period[0] << endl;
            }
        }

        end_t = clock();
        iteration_t = end_t - start_t_lu;
        // g_fout << "Link update with CPU time " << iteration_t / 1000.0 << " s; " << (end_t - start_t) / 1000.0 << " s" << endl;

        //****************************************//
        //step 3.2 computng block for continuous variables;

        clock_t start_t_lc = clock();
        clock_t start_t_cp = clock();

        cumulative_lc = 0;
        cumulative_cp = 0;
        cumulative_lu = 0;

        int number_of_memory_blocks = min((int)g_NetworkForSP_vector.size(), assignment.g_number_of_memory_blocks);

#pragma omp parallel for  // step 3: C++ open mp automatically create n threads., each thread has its own computing thread on a cpu core
        //for (int ProcessID = 0; ProcessID < g_NetworkForSP_vector.size(); ++ProcessID)
        //{
        //    int agent_type_no = g_NetworkForSP_vector[ProcessID]->m_agent_type_no;

        //    for (int o_node_index = 0; o_node_index < g_NetworkForSP_vector[ProcessID]->m_origin_node_vector.size(); ++o_node_index)
        //    {
        //        start_t_lc = clock();
        //        g_NetworkForSP_vector[ProcessID]->optimal_label_correcting(ProcessID, &assignment, iteration_number, o_node_index);
        //        end_t = clock();
        //        cumulative_lc += end_t - start_t_lc;

        //        start_t_cp = clock();
        //        g_NetworkForSP_vector[ProcessID]->backtrace_shortest_path_tree(assignment, iteration_number, o_node_index);
        //        end_t = clock();
        //        cumulative_cp += end_t - start_t_cp;
        //    }
        //    // perform one to all shortest path tree calculation
        //}

        for (int ProcessID = 0; ProcessID < number_of_memory_blocks; ++ProcessID)
        {
            for (int blk = 0; blk < assignment.g_AgentTypeVector.size() * assignment.g_DemandPeriodVector.size(); ++blk)
            {
                NetworkForSP* pNetwork = g_NetworkForSP_vector[blk * assignment.g_number_of_memory_blocks + ProcessID];

                for (int o_node_index = 0; o_node_index < pNetwork->m_origin_node_vector.size(); ++o_node_index)
                {
                    start_t_lc = clock();
                    pNetwork->optimal_label_correcting(ProcessID, &assignment, iteration_number, o_node_index);


                    end_t = clock();
                    cumulative_lc += end_t - start_t_lc;

                    start_t_cp = clock();
                    double total_origin_least_travel_time = pNetwork->backtrace_shortest_path_tree(assignment, iteration_number, o_node_index);


#pragma omp critical
                    {
                        total_least_system_travel_time += total_origin_least_travel_time;
                    }
                    end_t = clock();
                    cumulative_cp += end_t - start_t_cp;
                }
            }

        }

        // link based computing mode, we have to collect link volume from all processors.
        if (assignment.assignment_mode == 0)
            g_fetch_link_volume_for_all_processors();

        // g_fout << "LC with CPU time " << cumulative_lc / 1000.0 << " s; " << endl;
        // g_fout << "column generation with CPU time " << cumulative_cp / 1000.0 << " s; " << endl;

        //****************************************//

        // last iteraion before performing signal timing updating
        if (signal_updating_iterations >= 1 && iteration_number >= signal_updating_iterations - 1)
            g_output_simulation_result_for_signal_api(assignment);

        double relative_gap = (total_system_travel_time - total_least_system_travel_time) / max(0.00001, total_least_system_travel_time);
        dtalog.output() << "iteration: " << iteration_number << ",systemTT: " << total_system_travel_time << ", least system TT:" <<
            total_least_system_travel_time << ",gap=" << relative_gap << endl;
    }
    dtalog.output() << endl;

    // step 1.8: column updating stage: for given column pool, update volume assigned for each column
    dtalog.output() << "Step 4: Column Pool Updating" << endl;
    dtalog.output() << "Total Column Pool Updating iteration: " << column_updating_iterations << endl;
    start_t = clock();
    g_column_pool_optimization(assignment, column_updating_iterations);
    dtalog.output() << endl;

    // post route assignment aggregation
    if (assignment.assignment_mode != 0)
    {
        // we can have a recursive formulat to reupdate the current link volume by a factor of k/(k+1),
        // and use the newly generated path flow to add the additional 1/(k+1)
        g_reset_and_update_link_volume_based_on_columns(g_link_vector.size(), iteration_number, false);
    }
    else
        g_reset_link_volume_in_master_program_without_columns(g_link_vector.size(), iteration_number, false);

    // initialization at the first iteration of shortest path
    update_link_travel_time_and_cost();

    if (assignment.assignment_mode == 3)
    {
        dtalog.output() << "Step 6: O-D estimation for traffic assignment.." << endl;
        assignment.Demand_ODME(column_updating_iterations);
        dtalog.output() << endl;
    }

    if (assignment.assignment_mode == 2
        || assignment.assignment_mode == 3
        || assignment.assignment_mode == 10)
    {
        dtalog.output() << "Step 5: Simulation for traffic assignment.." << endl;
        assignment.STTrafficSimulation();
        dtalog.output() << endl;
    }



    end_t = clock();
    total_t = (end_t - start_t);
    dtalog.output() << "Done!" << endl;

    dtalog.output() << "CPU Running Time for column pool updating: " << total_t / 1000.0 << " s" << endl;

    start_t = clock();

    //step 5: output simulation results of the new demand
    g_output_simulation_result(assignment);
    g_output_TD_queue_result(assignment);


    end_t = clock();
    total_t = (end_t - start_t);
    dtalog.output() << "Output for assignment with " << assignment.g_number_of_column_generation_iterations << " iterations. Traffic assignment completes!" << endl;
    dtalog.output() << "CPU Running Time for outputting simulation results: " << total_t / 1000.0 << " s" << endl;

    dtalog.output() << "free memory.." << endl;
    g_node_vector.clear();

    for (int i = 0; i < g_link_vector.size(); ++i)
        g_link_vector[i].free_memory();
    g_link_vector.clear();

    dtalog.output() << "done." << endl;

    return 1;
}

void Assignment::AllocateLinkMemory4Simulation()
{
    g_number_of_simulation_intervals = (g_LoadingEndTimeInMin - g_LoadingStartTimeInMin + 60) * 60 /number_of_seconds_per_interval;
    g_number_of_loading_intervals = (g_LoadingEndTimeInMin - g_LoadingStartTimeInMin) * 60 / number_of_seconds_per_interval;

    g_number_of_simulation_horizon_in_min = (int)(g_number_of_simulation_intervals / number_of_interval_per_min +1);
    // add + 120 as a buffer
    g_number_of_in_memory_simulation_intervals = g_number_of_simulation_intervals;

    dtalog.output() << "allocate 2D dynamic memory LinkOutFlowCapacity..." << endl;

    m_LinkOutFlowCapacity = AllocateDynamicArray <float>(g_number_of_links, g_number_of_simulation_intervals);  //1
    // discharge rate per simulation time interval
    dtalog.output() << "allocate 2D dynamic memory m_LinkCumulativeArrival..." << endl;
    m_LinkCumulativeArrival = AllocateDynamicArray <float>(g_number_of_links, g_number_of_simulation_intervals);  //2

    dtalog.output() << "allocate 2D dynamic memory m_LinkCumulativeDeparture..." << endl;
    m_LinkCumulativeDeparture = AllocateDynamicArray <float>(g_number_of_links, g_number_of_simulation_intervals);  //3

    dtalog.output() << "allocate 2D dynamic memory m_LinkTDTravelTime..." << endl;
    m_LinkTDTravelTime = AllocateDynamicArray <int>(g_number_of_links, g_number_of_simulation_horizon_in_min); //4

    dtalog.output() << "allocate 2D dynamic memory m_LinkTDWaitingTime..." << endl;
    m_LinkTDWaitingTime = AllocateDynamicArray <float>(g_number_of_links, g_number_of_simulation_horizon_in_min); //5

    dtalog.output() << "initializing time dependent capacity data..." << endl;

    unsigned int RandomSeed = 101;
    float residual;
    float random_ratio = 0;

#pragma omp parallel for
    for (int i = 0; i < g_number_of_links; ++i)
    {
        float cap_count = 0;
        float discharge_rate = g_link_vector[i].lane_capacity * g_link_vector[i].number_of_lanes / 3600.0 * number_of_seconds_per_interval;
        float discharge_rate_after_loading = 10 * g_link_vector[i].lane_capacity * g_link_vector[i].number_of_lanes / 3600.0 * number_of_seconds_per_interval;

        for (int t = 0; t < g_number_of_simulation_intervals; ++t)
        {
            m_LinkTDTravelTime[i][t/ number_of_interval_per_min] = max(1, (int)(g_link_vector[i].free_flow_travel_time_in_min * 60 / number_of_seconds_per_interval));

            if (t >= g_number_of_loading_intervals)
                m_LinkOutFlowCapacity[i][t] = discharge_rate_after_loading;
                /* 10 times of capacity to discharge all flow */
            else
                m_LinkOutFlowCapacity[i][t] = discharge_rate;

            m_LinkTDWaitingTime[i][t / number_of_interval_per_min] = 0;

            residual = m_LinkOutFlowCapacity[i][t] - (int)(m_LinkOutFlowCapacity[i][t]);
            //RandomSeed is automatically updated.
            RandomSeed = (LCG_a * RandomSeed + LCG_c) % LCG_M;
            random_ratio = float(RandomSeed) / LCG_M;

            if (random_ratio < residual)
                m_LinkOutFlowCapacity[i][t] = (int)(m_LinkOutFlowCapacity[i][t]) +1;
            else
                m_LinkOutFlowCapacity[i][t] = (int)(m_LinkOutFlowCapacity[i][t]);

            cap_count += m_LinkOutFlowCapacity[i][t];

            // convert per hour capacity to per second and then to per simulation interval
            m_LinkCumulativeArrival[i][t] = 0;
            m_LinkCumulativeDeparture[i][t] = 0;
        }
    }

    // for each link with  for link type code is 's'
    //reset the time-dependent capacity to zero
    //go through the records defined in service_arc file
    //only enable m_LinkOutFlowCapacity[l][t] for the timestamp in the time window per time interval and reset the link TD travel time.

    for (unsigned li = 0; li < g_link_vector.size(); ++li)
    {
        if(g_link_vector[li].service_arc_flag && assignment.g_LinkTypeMap[g_link_vector[li].link_type].type_code != "f")
        {
            // reset for signalized links (not freeway links as type code != 'f' for the case of freeway workzones)
            // only for the loading period
            for (int t = 0; t < g_number_of_loading_intervals; ++t)
                m_LinkOutFlowCapacity[li][t] = 0;
        }
    }

    for (int si = 0; si < g_service_arc_vector.size(); ++si)
    {
        CServiceArc* pServiceArc = &(g_service_arc_vector[si]);
        int l = pServiceArc->link_seq_no;

        int number_of_cycles = (g_LoadingEndTimeInMin - g_LoadingStartTimeInMin) * 60 / max(1, pServiceArc->cycle_length );  // unit: seconds;

        for(int cycle_no = 0; cycle_no < number_of_cycles; ++cycle_no)
        {
            int count = 0;
            // relative time horizon
            for (int t = cycle_no * pServiceArc->cycle_length + pServiceArc->starting_time_no; t <= cycle_no * pServiceArc->cycle_length + pServiceArc->ending_time_no; t += pServiceArc->time_interval_no)
                count++;

            // relative time horizon
            for (int t = cycle_no * pServiceArc->cycle_length + pServiceArc->starting_time_no; t <= cycle_no * pServiceArc->cycle_length + pServiceArc->ending_time_no; t+= pServiceArc->time_interval_no)
            {
                // active capacity for this space time arc
                m_LinkOutFlowCapacity[l][t] = pServiceArc->capacity/max(1, count);

                residual = m_LinkOutFlowCapacity[l][t] - (int)(m_LinkOutFlowCapacity[l][t]);
                //RandomSeed is automatically updated.
                RandomSeed = (LCG_a * RandomSeed + LCG_c) % LCG_M;
                random_ratio = float(RandomSeed) / LCG_M;

                if (random_ratio < residual)
                    m_LinkOutFlowCapacity[l][t] = (int)(m_LinkOutFlowCapacity[l][t]) + 1;
                else
                    m_LinkOutFlowCapacity[l][t] = (int)(m_LinkOutFlowCapacity[l][t]);

                // enable time-dependent travel time
                m_LinkTDTravelTime[l][t/ number_of_interval_per_min] = pServiceArc->travel_time_delta;
                m_LinkTDWaitingTime[l][t / number_of_interval_per_min] = 0;
            }
        }
    }

    dtalog.output() << "End of initializing time dependent capacity data." << endl;
}

void Assignment::DeallocateLinkMemory4Simulation()
{
	// g_fout << "deallocate 2D dynamic memory m_LinkOutFlowCapacity..." << endl;
    if(m_LinkOutFlowCapacity)
        DeallocateDynamicArray(m_LinkOutFlowCapacity , g_number_of_links, g_number_of_simulation_intervals);  //1
	// g_fout << "deallocate 2D dynamic memory m_LinkCumulativeArrival..." << endl;
    if(m_LinkCumulativeArrival)
        DeallocateDynamicArray(m_LinkCumulativeArrival, g_number_of_links, g_number_of_simulation_intervals);  //2
	// g_fout << "deallocate 2D dynamic memory m_LinkCumulativeDeparture..." << endl;
    if (m_LinkCumulativeDeparture)
        DeallocateDynamicArray(m_LinkCumulativeDeparture, g_number_of_links, g_number_of_simulation_intervals); //3
	// g_fout << "deallocate 2D dynamic memory m_LinkTDTravelTime..." << endl;
    if (m_LinkTDTravelTime)
        DeallocateDynamicArray(m_LinkTDTravelTime, g_number_of_links, g_number_of_simulation_horizon_in_min); //3
	// g_fout << "deallocate 2D dynamic memory m_LinkTDWaitingTime..." << endl;
    if (m_LinkTDWaitingTime)
        DeallocateDynamicArray(m_LinkTDWaitingTime, g_number_of_links, g_number_of_simulation_horizon_in_min); //4
}

void Assignment::STTrafficSimulation()
{
    //given p_agent->path_link_seq_no_vector path link sequence no for each agent
    int TotalCumulative_Arrival_Count = 0;
    int TotalCumulative_Departure_Count = 0;

    clock_t start_t;
    start_t = clock();

    AllocateLinkMemory4Simulation();

    int agent_type_size = g_AgentTypeVector.size();
    int zone_size = g_zone_vector.size();
    int demand_period_size = g_DemandPeriodVector.size();

    CColumnVector* p_column;
    float path_toll = 0;
    float path_distance = 0;
    float path_travel_time = 0;
    float time_stamp = 0;

    std::map<int, CColumnPath>::iterator it, it_begin, it_end;

    for (int orig = 0; orig < zone_size; ++orig)
    {
        if (orig % 100 == 0)
            dtalog.output() << "generating " << g_agent_simu_vector.size()/1000 << " K agents for "  << orig << "  zones " << endl;

        for (int at = 0; at < agent_type_size; ++at)
        {
            for (int dest = 0; dest < zone_size; ++dest)
            {
                for (int tau = 0; tau < demand_period_size; ++tau)
                {
                    p_column = &(assignment.g_column_pool[orig][dest][at][tau]);
                    if (p_column->od_volume > 0)
                    {
                        // scan through the map with different node sum for different continuous paths
                        it_begin = p_column->path_node_sequence_map.begin();
                        it_end = p_column->path_node_sequence_map.end();

                        for (it = it_begin; it != it_end; ++it)
                        {
                            path_toll = 0;
                            path_distance = 0;
                            path_travel_time = 0;

                            int VehicleSize = (it->second.path_volume +0.5);

                            for(int v = 0; v < VehicleSize; ++v)
                            {
                                if (it->second.path_volume < 1)
                                    time_stamp = it->second.path_volume * (assignment.g_DemandPeriodVector[tau].ending_time_slot_no - assignment.g_DemandPeriodVector[tau].starting_time_slot_no) *MIN_PER_TIMESLOT ;
                                else
                                    time_stamp = v *1.0 / it->second.path_volume * (assignment.g_DemandPeriodVector[tau].ending_time_slot_no - assignment.g_DemandPeriodVector[tau].starting_time_slot_no) * MIN_PER_TIMESLOT;

                                if (it->second.m_link_size == 0)   // only load agents with physical path
                                    continue;

                                CAgent_Simu* pAgent = new CAgent_Simu();
                                pAgent->agent_id = g_agent_simu_vector.size();
                                pAgent->departure_time_in_min = time_stamp;

                                it->second.agent_simu_id_vector.push_back(pAgent->agent_id);

                                int simulation_time_intervalNo = (int)(pAgent->departure_time_in_min * 60/ number_of_seconds_per_interval);
                                g_AgentTDListMap[simulation_time_intervalNo].m_AgentIDVector.push_back(pAgent->agent_id);

                                for (int nl = 0; nl < it->second.m_link_size; ++nl)  // arc a
                                {
                                    int link_seq_no = it->second.path_link_vector[nl];
                                    pAgent->path_link_seq_no_vector.push_back(link_seq_no);
                                }
                                pAgent->AllocateMemory();

                                int FirstLink = pAgent->path_link_seq_no_vector[0];

                                pAgent->m_Veh_LinkArrivalTime_in_simu_interval[0] = simulation_time_intervalNo;
                                pAgent->m_Veh_LinkDepartureTime_in_simu_interval[0] = pAgent->m_Veh_LinkArrivalTime_in_simu_interval[0] + m_LinkTDTravelTime[FirstLink][simulation_time_intervalNo/ number_of_interval_per_min];

                                g_agent_simu_vector.push_back(pAgent);
                            }
                        }
                    }
                }
            }
        }
    }

    dtalog.output() << "number of simulation zones:" << zone_size << endl;

    int current_active_agent_id = 0;
    // the number of threads is redifined.
    int number_of_threads = omp_get_max_threads();
    // first loop for time t
    for (int t = 0; t < g_number_of_simulation_intervals; ++t)
    {
        int link_size = g_link_vector.size();
        for (int li = 0; li < link_size; ++li)
        {
            CLink* pLink = &(g_link_vector[li]);
            if (t >= 1)
            {
                m_LinkCumulativeDeparture[li][t] = m_LinkCumulativeDeparture[li][t - 1];
                m_LinkCumulativeArrival[li][t] = m_LinkCumulativeArrival[li][t - 1];
            }
        }

        int number_of_simu_interval_per_min = 60 / number_of_seconds_per_interval;
        if (t % (number_of_simu_interval_per_min*10) == 0)
            dtalog.output() << "simu time= " << t / number_of_simu_interval_per_min << " min, CA = " << TotalCumulative_Arrival_Count << " CD=" << TotalCumulative_Departure_Count << endl;

        if (g_AgentTDListMap.find(t) != g_AgentTDListMap.end())
        {
            for (int Agent_v = 0; Agent_v < g_AgentTDListMap[t].m_AgentIDVector.size(); ++Agent_v)
            {
                int agent_id = g_AgentTDListMap[t].m_AgentIDVector[Agent_v];

                CAgent_Simu* p_agent = g_agent_simu_vector[agent_id];
                p_agent->m_bGenereated = true;
                int FirstLink = p_agent->path_link_seq_no_vector[0];
                m_LinkCumulativeArrival[FirstLink][t] += 1;

                g_link_vector[FirstLink].EntranceQueue.push_back(p_agent->agent_id);
                TotalCumulative_Arrival_Count++;
            }
        }

#pragma omp parallel for    // parallel computing for each link
        for (int li = 0; li < link_size; ++li)
        {
            CLink* pLink = &(g_link_vector[li]);

            // if there are Agents in the entrance queue
            while (pLink->EntranceQueue.size() > 0)
            {
                int agent_id = pLink->EntranceQueue.front();
                pLink->EntranceQueue.pop_front();
                pLink->ExitQueue.push_back(agent_id);
                CAgent_Simu* p_agent = g_agent_simu_vector[agent_id];
                int arrival_time = p_agent->m_Veh_LinkArrivalTime_in_simu_interval[p_agent->m_current_link_seq_no];
                p_agent->m_Veh_LinkDepartureTime_in_simu_interval[p_agent->m_current_link_seq_no] = arrival_time + m_LinkTDTravelTime[li][arrival_time/ number_of_interval_per_min];
            }
        }

        int node_size = g_node_vector.size();

#pragma omp parallel for  // parallel computing for each node
        for (int node = 0; node < node_size; ++node)  // for each node
        {
            // for each incoming link
            for (int i = 0; i < g_node_vector[node].m_incoming_link_seq_no_vector.size(); ++i)
            {
                int incoming_link_index = (i + t)% (g_node_vector[node].m_incoming_link_seq_no_vector.size());
                int link = g_node_vector[node].m_incoming_link_seq_no_vector[incoming_link_index];  // we will start with different first link from the incoming link list,
                // equal change, regardless of # of lanes or main line vs. ramp, but one can use service arc, to control the effective capacity rates, e.g. through a metered ramp, to
                // allow mainline to use the remaining flow
                CLink* pLink = &(g_link_vector[link]);

                // check if the current link has sufficient capacity
                // most critical and time-consuming task, check link outflow capacity
                while (m_LinkOutFlowCapacity[link][t] >= 1 && pLink->ExitQueue.size() >=1)
                {
                    int agent_id = pLink->ExitQueue.front();
                    CAgent_Simu* p_agent = g_agent_simu_vector[agent_id];

                    if (p_agent->m_Veh_LinkDepartureTime_in_simu_interval[p_agent->m_current_link_seq_no] > t)
                    {
                        // the future departure time on this link is later than the current time
                        break;
                    }

                    if (p_agent->m_current_link_seq_no == p_agent->path_link_seq_no_vector.size() - 1)
                    {
                        // end of path
                        pLink->ExitQueue.pop_front();
                        p_agent->m_bCompleteTrip = true;
                        m_LinkCumulativeDeparture[link][t] += 1;

                        #pragma omp critical
                        {
                            TotalCumulative_Departure_Count += 1;
                        }

                    }
                    else
                    {
                        // not complete the trip. move to the next link's entrance queue
                        int next_link_seq_no = p_agent->path_link_seq_no_vector[p_agent->m_current_link_seq_no + 1];
                        CLink* pNextLink = &(g_link_vector[next_link_seq_no]);

                        // spatial queue
                        if (pNextLink->traffic_flow_code == 2)
                        {
                            int current_vehicle_count = m_LinkCumulativeArrival[next_link_seq_no][t - 1] - m_LinkCumulativeDeparture[next_link_seq_no][t - 1];
                            if(current_vehicle_count > pNextLink->spatial_capacity_in_vehicles)
                            {
                                // spatical queue in the next link is blocked, break the while loop from here, as a first in first out queue.
								// g_fout << "spatical queue in the next link is blocked on link seq  " <<  g_node_vector[pNextLink->from_node_seq_no].node_id  << " -> " << g_node_vector[pNextLink->to_node_seq_no].node_id  <<endl;
                                break;
                            }
                        }

                        // kinematic wave
                        if (pNextLink->traffic_flow_code == 3)
                        {
                            int lagged_time_stamp = max(0, t - 1 - pNextLink->BWTT_in_simulation_interval);
                            int current_vehicle_count = m_LinkCumulativeArrival[next_link_seq_no][t - 1] - m_LinkCumulativeDeparture[next_link_seq_no][lagged_time_stamp];
                            if (current_vehicle_count > pNextLink->spatial_capacity_in_vehicles)
                            {
                                // spatical queue in the next link is blocked, break the while loop from here, as a first in first out queue.
								// g_fout << "spatical queue in the next link is blocked on link seq  " <<  g_node_vector[pNextLink->from_node_seq_no].node_id  << " -> " << g_node_vector[pNextLink->to_node_seq_no].node_id  <<endl;
                                break;
                            }
                        }

                        pLink->ExitQueue.pop_front();
                        pNextLink->EntranceQueue.push_back(agent_id);
                        p_agent->m_Veh_LinkDepartureTime_in_simu_interval[p_agent->m_current_link_seq_no] = t;
                        p_agent->m_Veh_LinkArrivalTime_in_simu_interval[p_agent->m_current_link_seq_no + 1] = t;

                        float travel_time = p_agent->m_Veh_LinkDepartureTime_in_simu_interval[p_agent->m_current_link_seq_no] - p_agent->m_Veh_LinkArrivalTime_in_simu_interval[p_agent->m_current_link_seq_no];
                        //for each waited vehicle
                        float waiting_time = travel_time - m_LinkTDTravelTime[link][p_agent->m_Veh_LinkArrivalTime_in_simu_interval[p_agent->m_current_link_seq_no]/ number_of_interval_per_min];

                        m_LinkTDWaitingTime[link][p_agent->m_Veh_LinkArrivalTime_in_simu_interval[p_agent->m_current_link_seq_no] / number_of_interval_per_min] += waiting_time;

                        m_LinkCumulativeDeparture[link][t] += 1;
                        m_LinkCumulativeArrival[next_link_seq_no][t] += 1;
                    }

                    //move
                    p_agent->m_current_link_seq_no += 1;
                    m_LinkOutFlowCapacity[link][t] -= 1;
                }
            }
        } // conditions
    }  // departure time events
}

bool Assignment::Map_TMC_Reading()
{
    // step 1: read measurement.csv
    CCSVParser parser_measurement;


    if (parser_measurement.OpenCSVFile("Reading.csv", true))
    {
        int count = 0;
        while (parser_measurement.ReadRecord())  // if this line contains [] mark, then we will also read field headers.
        {
            string tmc;

            parser_measurement.GetValueByFieldName("tmc_code", tmc);
            string measurement_tstamp;

            float tmc_reference_speed = 0;
            float speed = 0;
            bool bMatched_flag = false;
            parser_measurement.GetValueByFieldName("measurement_tstamp", measurement_tstamp, false);

            float global_time;
            int day_of_week_flag = 0; 
            int day_of_year = 0;

            if (measurement_tstamp.size() < 18)
                continue; // skip empty lines 

            global_time = g_measurement_tstamp_parser(measurement_tstamp, day_of_week_flag, day_of_year);
            
            if (day_of_week_flag == 0 && day_of_week_flag == 6)
                continue; 

            parser_measurement.GetValueByFieldName("speed", speed, false);
           
            parser_measurement.GetValueByFieldName("reference_speed", tmc_reference_speed, false);

            string ROADNAME;
            parser_measurement.GetValueByFieldName("ROADNAME", ROADNAME, false);

            if (assignment.m_TMClink_map.find(tmc) == assignment.m_TMClink_map.end())
            {
                CTMCLink tmc_link;
                tmc_link.tmc_code = tmc;
                assignment.m_TMClink_map[tmc] = g_TMC_vector.size();
                g_TMC_vector.push_back(tmc_link);

            }
            
            int index = assignment.m_TMClink_map[tmc];

            g_TMC_vector[index].AddSpeedData(day_of_year,global_time,speed);



            if (count % 100000==0)
            {
                dtalog.output() << "reading " << count/100000 << "00k TMC data items" << endl;

            }
            count++;
        }
        parser_measurement.CloseCSVFile();

            dtalog.output() << "reading data for " << g_TMC_vector.size() << " TMC links." << endl;
        return true;
    }
    else
    {

    }

    return false;
  
}

// updates for OD re-generations
void Assignment::Demand_ODME(int OD_updating_iterations)
{
    // step 1: read measurement.csv
    CCSVParser parser_measurement;

    if (parser_measurement.OpenCSVFile("measurement.csv", true))
    {
        while (parser_measurement.ReadRecord())  // if this line contains [] mark, then we will also read field headers.
        {
            string measurement_type = "link";
            parser_measurement.GetValueByFieldName("measurement_type", measurement_type);

            int upper_bound_flag = 0;
            parser_measurement.GetValueByFieldName("upper_bound_flag", upper_bound_flag);

            if (measurement_type == "link")
            {
                int from_node_id;
                if (!parser_measurement.GetValueByFieldName("from_node_id", from_node_id))
                    continue;

                int to_node_id;
                if (!parser_measurement.GetValueByFieldName("to_node_id", to_node_id))
                    continue;

                // add the to node id into the outbound (adjacent) node list
                if (g_node_id_to_seq_no_map.find(from_node_id) == assignment.g_node_id_to_seq_no_map.end())
                {
                    dtalog.output() << "Error: from_node_id " << from_node_id << " in file measurement.csv is not defined in node.csv." << endl;
                    //has not been defined
                    continue;
                }
                if (g_node_id_to_seq_no_map.find(to_node_id) == assignment.g_node_id_to_seq_no_map.end())
                {
                    dtalog.output() << "Error: to_node_id " << to_node_id << " in file measurement.csv is not defined in node.csv." << endl;
                    //has not been defined
                    continue;
                }

                float count = -1;
                parser_measurement.GetValueByFieldName("count", count);

                // map external node number to internal node seq no.
                int internal_from_node_seq_no = assignment.g_node_id_to_seq_no_map[from_node_id];
                int internal_to_node_seq_no = assignment.g_node_id_to_seq_no_map[to_node_id];

                if (g_node_vector[internal_from_node_seq_no].m_to_node_2_link_seq_no_map.find(internal_to_node_seq_no) != g_node_vector[internal_from_node_seq_no].m_to_node_2_link_seq_no_map.end())
                {
                    int link_seq_no = g_node_vector[internal_from_node_seq_no].m_to_node_2_link_seq_no_map[internal_to_node_seq_no];

                    g_link_vector[link_seq_no].obs_count = count;
                    g_link_vector[link_seq_no].upper_bound_flag = upper_bound_flag;
                }
                else
                {
                    dtalog.output() << "Error: Link " << from_node_id << "->" << to_node_id << " in file timing.csv is not defined in link.csv." << endl;
                    continue;
                }
            }

            if (measurement_type == "production")
            {
                int o_zone_id;
                if (!parser_measurement.GetValueByFieldName("o_zone_id", o_zone_id))
                    continue;

                if (g_zoneid_to_zone_seq_no_mapping.find(o_zone_id) != g_zoneid_to_zone_seq_no_mapping.end())
                {
                    float obs_production = -1;
                    if (parser_measurement.GetValueByFieldName("count", obs_production))
                    {
                        g_zone_vector[g_zoneid_to_zone_seq_no_mapping[o_zone_id]].obs_production = obs_production;
                        g_zone_vector[g_zoneid_to_zone_seq_no_mapping[o_zone_id]].obs_production_upper_bound_flag = upper_bound_flag;
                    }
                }
            }

            if (measurement_type == "attraction")
            {
                int o_zone_id;
                if (!parser_measurement.GetValueByFieldName("d_zone_id", o_zone_id))
                    continue;

                if (g_zoneid_to_zone_seq_no_mapping.find(o_zone_id) != g_zoneid_to_zone_seq_no_mapping.end())
                {
                    float obs_attraction = -1;
                    if (parser_measurement.GetValueByFieldName("count", obs_attraction))
                    {
                        g_zone_vector[g_zoneid_to_zone_seq_no_mapping[o_zone_id]].obs_attraction = obs_attraction;
                        g_zone_vector[g_zoneid_to_zone_seq_no_mapping[o_zone_id]].obs_attraction_upper_bound_flag = upper_bound_flag;
                    }
                }
            }
        }

        parser_measurement.CloseCSVFile();
    }




   // step 1: input the measurements of
    // Pi
    // Dj
    // link l

    // step 2: loop for adjusting OD demand
    double prev_gap = 9999999;
    for (int s = 0; s < OD_updating_iterations; ++s)
    {
        float total_gap = 0;
        float total_relative_gap = 0;
        float total_gap_count = 0;
        //step 2.1
        // we can have a recursive formulat to reupdate the current link volume by a factor of k/(k+1),
        // and use the newly generated path flow to add the additional 1/(k+1)
        double gap = g_reset_and_update_link_volume_based_on_ODME_columns(g_link_vector.size(),s);

        double gap_improvement = gap - prev_gap;

        if (s>=1 && gap_improvement > 0.001)  // convergency criterion
            break;

        prev_gap = gap;

        //step 2.2: based on newly calculated path volumn, update volume based travel time, and update volume based measurement error/deviation

        //step 3: calculate shortest path at inner iteration of column flow updating
//#pragma omp parallel for
        for (int orig = 0; orig < g_zone_vector.size(); ++orig)  // o
        {
            CColumnVector* p_column;
            std::map<int, CColumnPath>::iterator it, it_begin, it_end;
            int column_vector_size;
            int path_seq_count = 0;

            float path_toll = 0;
            float path_gradient_cost = 0;
            float path_distance = 0;
            float path_travel_time = 0;

            int link_seq_no;

            float total_switched_out_path_volume = 0;

            float step_size = 0;
            float previous_path_volume = 0;

            for (int dest = 0; dest < g_zone_vector.size(); ++dest) //d
            {
                for (int at = 0; at < assignment.g_AgentTypeVector.size(); ++at)  //at
                {
                    string agent_type = assignment.g_AgentTypeVector[at].agent_type;

                    // only focus on auto mode based ODME now
                    if (agent_type == "bike" || agent_type == "cav" || agent_type == "ev" || agent_type == "truck")
                        continue; 

                    for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); ++tau)  //tau
                    {
                        p_column = &(assignment.g_column_pool[orig][dest][at][tau]);
                        if (p_column->od_volume > 0 && !p_column->bfixed_route)
                        {
                            column_vector_size = p_column->path_node_sequence_map.size();
                            path_seq_count = 0;

                            it_begin = p_column->path_node_sequence_map.begin();
                            it_end = p_column->path_node_sequence_map.end();
                            int i = 0;
                            for (it = it_begin; it != it_end; ++it, ++i) // for each k
                            {
                                path_gradient_cost = 0;
                                path_distance = 0;
                                path_travel_time = 0;

                                // step 3.1 origin production flow gradient

                                // est_production_dev = est_production - obs_production;
                                // requirement: when equality flag is 1,
                                if (g_zone_vector[orig].obs_production > 0)
                                {
                                    if(g_zone_vector[orig].obs_production_upper_bound_flag ==0)
                                        path_gradient_cost += g_zone_vector[orig].est_production_dev;

                                    if (g_zone_vector[orig].obs_production_upper_bound_flag == 1 && g_zone_vector[orig].est_production_dev>0)
                                        /*only if est_production is greater than obs value , otherwise, do not apply*/
                                        path_gradient_cost += g_zone_vector[orig].est_production_dev;
                                }

                                // step 3.2 destination attraction flow gradient

                                if (g_zone_vector[dest].obs_attraction > 0)
                                {
                                    if (g_zone_vector[orig].obs_attraction_upper_bound_flag == 0)
                                        path_gradient_cost += g_zone_vector[dest].est_attraction_dev;

                                    if (g_zone_vector[orig].obs_attraction_upper_bound_flag == 1 && g_zone_vector[dest].est_attraction_dev > 0)
                                        path_gradient_cost += g_zone_vector[dest].est_attraction_dev;
                                }

                                float est_count_dev = 0;
                                for (int nl = 0; nl < it->second.m_link_size; ++nl)  // arc a
                                {
                                    // step 3.3 link flow gradient
                                    link_seq_no = it->second.path_link_vector[nl];
                                    if(g_link_vector[link_seq_no].obs_count >= 1)
                                    {
                                        if (g_link_vector[link_seq_no].upper_bound_flag == 0)
                                        {
                                            path_gradient_cost += g_link_vector[link_seq_no].est_count_dev;
                                            est_count_dev += g_link_vector[link_seq_no].est_count_dev;
                                        }

                                        if (g_link_vector[link_seq_no].upper_bound_flag == 1 && g_link_vector[link_seq_no].est_count_dev > 0)
                                        {
                                            path_gradient_cost += g_link_vector[link_seq_no].est_count_dev;
                                            est_count_dev += g_link_vector[link_seq_no].est_count_dev;
                                        }
                                    }
                                }

                                it->second.path_gradient_cost = path_gradient_cost;

                                step_size = 0.01; // very small size 
                                float prev_path_volume = it->second.path_volume;
                                float change = step_size * it->second.path_gradient_cost;

                                float change_lower_bound = it->second.path_volume * 0.05 * (-1);
                                float change_upper_bound = it->second.path_volume * 0.05 ;

                                // reset
                                if (change < change_lower_bound)
                                    change = change_lower_bound;

                                // reset
                                if (change > change_upper_bound)
                                    change = change_upper_bound;

                                it->second.path_volume = max(1.0, it->second.path_volume - change);

                                if (dtalog.log_odme() == 1)
                                {
                                    dtalog.output() << "OD " << orig << "-> " << dest << " path id:" << i << ", prev_vol"
                                                    << prev_path_volume << ", gradient_cost = " << it->second.path_gradient_cost
                                                    << prev_path_volume << ", gradient_cost = " << it->second.path_gradient_cost
                                                    << prev_path_volume << ", gradient_cost = " << it->second.path_gradient_cost
                                                    << prev_path_volume << ", gradient_cost = " << it->second.path_gradient_cost
                                                    << prev_path_volume << ", gradient_cost = " << it->second.path_gradient_cost
                                                    << " link," << g_link_vector[link_seq_no].est_count_dev
                                                    << " P," << g_zone_vector[orig].est_production_dev
                                                    << " A," << g_zone_vector[orig].est_attraction_dev
                                                    << "proposed change = " << step_size * it->second.path_gradient_cost
                                                    << "proposed change = " << step_size * it->second.path_gradient_cost
                                                    << "proposed change = " << step_size * it->second.path_gradient_cost
                                                    << "proposed change = " << step_size * it->second.path_gradient_cost
                                                    << "proposed change = " << step_size * it->second.path_gradient_cost
                                                    << "actual change = " << change
                                                    <<"new vol = " << it->second.path_volume <<endl;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    //if (assignment.g_pFileDebugLog != NULL)
    //	fprintf(assignment.g_pFileDebugLog, "CU: iteration %d: total_gap=, %f,total_relative_gap=, %f,\n", s, total_gap, total_gap / max(0.0001, total_gap_count));
}




void Assignment::Mapping_TMC_Identification()
{
    // step 1: read measurement.csv
    CCSVParser parser_measurement;
    string geo_string;
    CGeometry geometry(geo_string);


    //double left = 100000000;
    //double right = -100000000;
    //double top = -1000000000;
    //double  bottom = 1000000000;

    //for (int i = 0; i < g_node_vector.size(); i++)
    //{
    //    // exapnd the grid boundary according to the nodes
    //    left = min(left, g_node_vector[i].x);
    //    right = max(right, g_node_vector[i].x);
    //    top = max(top, g_node_vector[i].y);
    //    bottom = min(bottom, g_node_vector[i].y);

    //}

    //int grid_size = 30;

    //double temp_resolution = (((right - left) / grid_size + (top - bottom) / grid_size)) / 2.0;

    //assignment.m_GridResolution = temp_resolution;
    std::map <__int64, GridLinkSet > Grid_link_id_map;

    for (int i = 0; i < g_link_vector.size(); ++i)
    {
        CLink* pLink = &(g_link_vector[i]);

        pLink->from_node_cell_id = g_GetCellID(g_node_vector[g_link_vector[i].from_node_seq_no].x, g_node_vector[g_link_vector[i].from_node_seq_no].y, assignment.m_GridResolution);
        pLink->to_node_cell_id = g_GetCellID(g_node_vector[g_link_vector[i].to_node_seq_no].x, g_node_vector[g_link_vector[i].to_node_seq_no].y, assignment.m_GridResolution);

        Grid_link_id_map[pLink->from_node_cell_id].m_LinkNoVector.push_back(i);

        if (pLink->from_node_cell_id != pLink->to_node_cell_id)  // different cell ids
        {
            Grid_link_id_map[pLink->to_node_cell_id].m_LinkNoVector.push_back(i);
        }
    }


    std::map<string, GDPoint > matching_TMC_pt0;
    std::map<string, GDPoint > matching_TMC_pt1;
    std::map<string, string > matching_TMC_road;
    std::map<string, int > matching_TMC_flag;
    std::map<string, int > TMC_corridor_id_map;

    dtalog.output() << "Mapping_TMC_Identification" << endl;
    cout << "Mapping_TMC_Identification";

    int count = 0;
    if (parser_measurement.OpenCSVFile("TMC_Identification.csv", true))
    {
        while (parser_measurement.ReadRecord())  // if this line contains [] mark, then we will also read field headers.
        {
            string tmc;

            parser_measurement.GetValueByFieldName("tmc", tmc);



            if (count % 100 == 0 || tmc.size() == 0)
            {

                dtalog.output() << "Mapping " << count << " TMC" << tmc << endl;
                cout << "Mapping " << count << " TMC" << tmc << endl;
                count++;
            }


            if (tmc == "112+06308")
            {
                int idebug = 1;
            }
            string tmc_road, tmc_direction, tmc_intersection, tmc_corridor_name, tmc_state;

            int corridor_id = -1;
            int road_order = 0;

            int road_sequence = 0;
            float tmc_reference_speed = 0;
            float tmc_mean_speed = 0;
            float tmc_volume = 0;

            parser_measurement.GetValueByFieldName("tmc_corridor_name", tmc_corridor_name, false);


            parser_measurement.GetValueByFieldName("road", tmc_road, false);
            parser_measurement.GetValueByFieldName("direction", tmc_direction, false);
            parser_measurement.GetValueByFieldName("state", tmc_state, false);
            parser_measurement.GetValueByFieldName("road_order", road_order, false);

            

            if (tmc_corridor_name.size() == 0)
                tmc_corridor_name = tmc_road + "_" + tmc_direction + "_" + tmc_state;

            if (corridor_id == -1)
            {
                corridor_id = TMC_corridor_id_map[tmc_corridor_name];
            }

            string link_type_code;
            parser_measurement.GetValueByFieldName("link_type_code", link_type_code, false);
            double start_latitude;
            double start_longitude;
            double end_latitude;
            double end_longitude;


            if (!parser_measurement.GetValueByFieldName("start_latitude", start_latitude, true, false))
                continue;

            if (!parser_measurement.GetValueByFieldName("start_longitude", start_longitude, true, false))
                continue;

            if (!parser_measurement.GetValueByFieldName("end_latitude", end_latitude, true, false))
                continue;

            if (!parser_measurement.GetValueByFieldName("end_longitude", end_longitude, true, false))
                continue;

            GDPoint pt0, pt1, FromPt, ToPt;
            pt0.x = start_longitude;
            pt0.y = start_latitude;

            pt1.x = end_longitude;
            pt1.y = end_latitude;

            __int64 from_pt_cell_id = g_GetCellID(pt0.x, pt0.y, assignment.m_GridResolution);
            __int64 to_pt_cell_id = g_GetCellID(pt1.x, pt1.y, assignment.m_GridResolution);


            bool bMatched_flag = false;
            if (Grid_link_id_map.find(from_pt_cell_id) == Grid_link_id_map.end())
                continue;

            std::vector<int> m_LinkNoVector;

            for (int k = 0; k < Grid_link_id_map[from_pt_cell_id].m_LinkNoVector.size(); k++)
            {
                m_LinkNoVector.push_back(Grid_link_id_map[from_pt_cell_id].m_LinkNoVector[k]);  // merge other link no ids
            }

            if (from_pt_cell_id != to_pt_cell_id)  // different cell ids
            {
                if (Grid_link_id_map.find(to_pt_cell_id) != Grid_link_id_map.end())
                {
                    for (int k = 0; k < Grid_link_id_map[to_pt_cell_id].m_LinkNoVector.size(); k++)
                    {
                        m_LinkNoVector.push_back(Grid_link_id_map[to_pt_cell_id].m_LinkNoVector[k]);  // merge other link no ids
                    }
                }

            }

            for (int k = 0; k < m_LinkNoVector.size(); ++k)
            {
                int i = m_LinkNoVector[k];
                CLink* pLink = &(g_link_vector[i]);

                if (pLink->TMC_code.size() > 0)  // TMC code has been predefined, skip mapping process
                {
                    matching_TMC_flag[pLink->TMC_code] = 1;
                    continue;
                }

                if (tmc_reference_speed >= 20 && fabs(tmc_reference_speed - pLink->free_speed) >= 10 && pLink->lane_capacity < 2500)  // prevent wrong type maching
                    continue;

                if (from_pt_cell_id != pLink->from_node_cell_id && to_pt_cell_id != pLink->to_node_cell_id)  // test if both cell ids of GPS points and link are different
                    continue;

                if (link_type_code.size() > 0)  // key condition for link type code based matching 
                {
                    if (pLink->link_type_code.size() > 0 && link_type_code.compare(pLink->link_type_code) != 0)
                        continue;


                }
                else
                {
                    // use geometry based matching by default 
                }

                // condition 1: if TMC_link_type code and the data content of link_type code exist in link.csv, we need to have a perfect string matching 
                // condition 2: if the data content link_type code does not exist in link.csv, we will use geometry based matching 


                int from_node_id = g_node_vector[g_link_vector[i].from_node_seq_no].node_id;
                int to_node_id = g_node_vector[g_link_vector[i].to_node_seq_no].node_id;

                if (pLink->link_type >= 0 && pLink->lane_capacity < 2500)  // prevent connectors
                {

                    FromPt.x = g_node_vector[pLink->from_node_seq_no].x;
                    FromPt.y = g_node_vector[pLink->from_node_seq_no].y;

                    ToPt.x = g_node_vector[pLink->to_node_seq_no].x;
                    ToPt.y = g_node_vector[pLink->to_node_seq_no].y;

                    if (geometry.g_GetTwoPoint2LineIntersectionFlag(&pt0, &pt1, &FromPt, &ToPt) == true)
                    {
                        geometry.g_GetTwoPoint2LineIntersectionFlag(&pt0, &pt1, &FromPt, &ToPt);

                        pLink->TMC_code = tmc;
                        pLink->tmc_road = tmc_road;
                        pLink->tmc_road_order = road_order;
                        pLink->tmc_corridor_name = tmc_corridor_name;
                        pLink->tmc_corridor_id = corridor_id;
                        pLink->tmc_road_sequence = road_sequence;
                        pLink->tmc_direction = tmc_direction;
                        pLink->tmc_intersection = tmc_intersection;
                        pLink->tmc_reference_speed = tmc_reference_speed;
                        pLink->tmc_mean_speed = tmc_mean_speed;
                        pLink->tmc_volume = tmc_volume;


                        pLink->TMC_from = pt0;
                        pLink->TMC_to = pt1;
                        bMatched_flag = true;
                    }
                }
            }  //end for all links

            // stage 2 matching: If a TMC has not found a matching link, I will find the closest link (regardless of TMC2Link distance) and mark the TMC for this link. 
            if (bMatched_flag == false)
            {

                double L2LDistance_minimum = 9999999;
                int link_no_minimum = -1;

                for (int k = 0; k < m_LinkNoVector.size(); ++k)
                {
                    int i = m_LinkNoVector[k];
                    CLink* pLink = &(g_link_vector[i]);

                    if (pLink->TMC_code.size() > 0)  // TMC code has been predefined, skip mapping process
                    {
                        matching_TMC_flag[pLink->TMC_code] = 1;
                        continue;
                    }

                    if (tmc_reference_speed >= 20 && fabs(tmc_reference_speed - pLink->free_speed) >= 10 && pLink->lane_capacity < 2500)  // prevent wrong type maching
                        continue;

                    if (from_pt_cell_id != pLink->from_node_cell_id && to_pt_cell_id != pLink->to_node_cell_id)  // test if both cell ids of GPS points and link are different
                        continue;

                    if (link_type_code.size() > 0)  // key condition for link type code based matching 
                    {
                        if (pLink->link_type_code.size() > 0 && link_type_code.compare(pLink->link_type_code) != 0)
                            continue;
                    }
                    else
                    {
                        // use geometry based matching by default 
                    }

                    // condition 1: if TMC_link_type code and the data content of link_type code exist in link.csv, we need to have a perfect string matching 
                    // condition 2: if the data content link_type code does not exist in link.csv, we will use geometry based matching 


                    int from_node_id = g_node_vector[g_link_vector[i].from_node_seq_no].node_id;
                    int to_node_id = g_node_vector[g_link_vector[i].to_node_seq_no].node_id;

                    if (pLink->link_type >= 0 && pLink->lane_capacity < 2500)  // prevent connectors
                    {

                        FromPt.x = g_node_vector[pLink->from_node_seq_no].x;
                        FromPt.y = g_node_vector[pLink->from_node_seq_no].y;

                        ToPt.x = g_node_vector[pLink->to_node_seq_no].x;
                        ToPt.y = g_node_vector[pLink->to_node_seq_no].y;



                        double L2LDistance = 999999;
                        geometry.g_GetTwoPointSegment2LineIntersectionFlag(&pt0, &pt1, &FromPt, &ToPt, L2LDistance);
                        {
                            if (L2LDistance < L2LDistance_minimum && L2LDistance < 1)
                            {
                                L2LDistance_minimum = L2LDistance;
                                link_no_minimum = i;
                            }


                        }
                    }
                }  //end for all links

                if (link_no_minimum >= 0 && link_no_minimum < g_link_vector.size())
                {

                    CLink* pLink = &(g_link_vector[link_no_minimum]);
                    pLink->TMC_code = tmc;
                    pLink->tmc_road = tmc_road;
                    pLink->tmc_road_order = road_order;
                    pLink->tmc_corridor_name = tmc_corridor_name;
                    pLink->tmc_corridor_id = corridor_id;
                    pLink->tmc_road_sequence = road_sequence;
                    pLink->tmc_direction = tmc_direction;
                    pLink->tmc_intersection = tmc_intersection;
                    pLink->tmc_reference_speed = tmc_reference_speed;
                    pLink->tmc_mean_speed = tmc_mean_speed;
                    pLink->tmc_volume = tmc_volume;

                    //if (pLink->link_id == "819")
                    //{
                    //    int idebug = 1;
                    //    double L2LDistance = 999999;
                    //    geometry.g_GetTwoPointSegment2LineIntersectionFlag(&pt0, &pt1, &FromPt, &ToPt, L2LDistance);
                    //}


                    pLink->TMC_from = pt0;
                    pLink->TMC_to = pt1;

                    bMatched_flag = true;

                }
            }


            if (bMatched_flag == false)
            {  // not matched:

                //dtalog.output() << "Warning:  not matched TMC code: ." << tmc << endl;
                //cout << "Warning:  not matched TMC code:";
                matching_TMC_flag[tmc] = 0;
            }
            else
            {
                matching_TMC_flag[tmc] = 1;

            }
            matching_TMC_pt0[tmc] = pt0;
            matching_TMC_pt1[tmc] = pt1;
            matching_TMC_road[tmc] = tmc_road;

        }
    }


    parser_measurement.CloseCSVFile();


    FILE* g_pFileTMCLink = fopen("TMC_mapping.csv", "w");

    if (g_pFileTMCLink != NULL)
    {
        fprintf(g_pFileTMCLink, "tmc,matching_flag,road,v_reference,geometry\n");

        std::map<string, GDPoint>::iterator it;

        for (it = matching_TMC_pt0.begin(); it != matching_TMC_pt0.end(); ++it)
        {
            fprintf(g_pFileTMCLink, "\"%s\",%d,", it->first.c_str(), matching_TMC_flag[it->first]);
            fprintf(g_pFileTMCLink, "\"%s\",", matching_TMC_road[it->first].c_str());
            fprintf(g_pFileTMCLink, "\"LINESTRING (");

            fprintf(g_pFileTMCLink, "%f %f,", matching_TMC_pt0[it->first].x, matching_TMC_pt0[it->first].y);
            fprintf(g_pFileTMCLink, "%f %f,", matching_TMC_pt1[it->first].x, matching_TMC_pt1[it->first].y);
            fprintf(g_pFileTMCLink, ")\"");
            fprintf(g_pFileTMCLink, "\n");
        }
        fclose(g_pFileTMCLink);
    }
    else
    {
        dtalog.output() << "Error: File TMC_mapping.csv cannot be opened.\n It might be currently used and locked by EXCEL." << endl;
        g_ProgramStop();

    }
}
