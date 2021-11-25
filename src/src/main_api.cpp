/* Portions Copyright 2019 Xuesong Zhou and Peiheng Li
 
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
constexpr auto _INFO_ZONE_ID = 100000;

constexpr auto _MAX_AGNETTYPES = 10; //because of the od demand store format,the MAX_demandtype must >=g_DEMANDTYPES.size()+1;
constexpr auto _MAX_TIMEPERIODS = 5; // time period set to 4: mid night, morning peak, mid-day and afternoon peak;
constexpr auto _MAX_MEMORY_BLOCKS = 100;

constexpr auto _MAX_LINK_SIZE_IN_A_PATH = 5000;		// lu
constexpr auto _MAX_LINK_SIZE_FOR_A_NODE = 200;

constexpr auto _MAX_TIMESLOT_PerPeriod = 100; // max 96 15-min slots per day
constexpr auto _default_saturation_flow_rate = 1530;

constexpr auto MIN_PER_TIMESLOT = 15;

/* make sure we change the following two parameters together*/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
constexpr auto number_of_seconds_per_interval = 0.25;  // consistent with the cell length of 7 meters
constexpr auto number_of_simu_interval_reaction_time = 4;  // reaction time as 1 second, 4 simu intervals, CAV: 0.5 seconds

constexpr auto number_of_simu_intervals_in_min = 240; // 60/0.25 number_of_seconds_per_interval

/* number_of_seconds_per_interval should satisify the ratio of 60/number_of_seconds_per_interval is an integer*/

// Linear congruential generator
constexpr auto LCG_a = 17364;
constexpr auto LCG_c = 0;
constexpr auto LCG_M = 65521;  // it should be 2^32, but we use a small 16-bit number to save memory

// FILE* g_pFileOutputLog = nullptr;

template <typename T>
T* Allocate1DDynamicArray(int nRows)
{
    T* dynamicVector;

    dynamicVector = new (std::nothrow) T[nRows]();

    if (dynamicVector == NULL)
    {
        exit(1);

    }
    return dynamicVector;
}

template <typename T>
void Deallocate1DDynamicArray(T* dVector, int nRows)
{
    if (!dVector)
        return;
    delete[] dVector;
}


template <typename T>
T** Allocate2DDynamicArray(int nRows, int nCols)
{
    T** dynamicArray;

    dynamicArray = new (std::nothrow) T*[nRows];

    if (!dynamicArray)
    {
        dtalog.output() << "Error: insufficient memory.";
        g_ProgramStop();
    }

    for (int i = 0; i < nRows; ++i)
    {
        dynamicArray[i] = new (std::nothrow) T[nCols];

        if (!dynamicArray[i])
        {
            dtalog.output() << "Error: insufficient memory.";
            g_ProgramStop();
        }
    }

    return dynamicArray;
}

template <typename T>
void Deallocate2DDynamicArray(T** dArray, int nRows, int nCols)
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
        g_ProgramStop();
    }

    for (int x = 0; x < nX; ++x)
    {
        if (x % 1000 == 0)
        {
            dtalog.output() << "allocating 3D memory for " << x << endl;
        }

        dynamicArray[x] = new (std::nothrow) T*[nY];

        if (!dynamicArray[x])
        {
            dtalog.output() << "Error: insufficient memory.";
            g_ProgramStop();
        }

        for (int y = 0; y < nY; ++y)
        {
            dynamicArray[x][y] = new (std::nothrow) T[nZ];
            if (!dynamicArray[x][y])
            {
                dtalog.output() << "Error: insufficient memory.";
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
        g_ProgramStop();
    }

    for (int m = 0; m < nM; ++m)
    {
        if (m % 1000 == 0)
            dtalog.output() << "allocating 4D memory for " << m << " zones" << endl;

        dynamicArray[m] = new (std::nothrow) T**[nX];

        if (!dynamicArray[m])
        {
            dtalog.output() << "Error: insufficient memory.";
            g_ProgramStop();
        }

        for (int x = 0; x < nX; ++x)
        {
            dynamicArray[m][x] = new (std::nothrow) T*[nY];

            if (!dynamicArray[m][x])
            {
                dtalog.output() << "Error: insufficient memory.";
                g_ProgramStop();
            }

            for (int y = 0; y < nY; ++y)
            {
                dynamicArray[m][x][y] = new (std::nothrow) T[nZ];
                if (!dynamicArray[m][x][y])
                {
                    dtalog.output() << "Error: insufficient memory.";
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
    y_code = fabs(yi) * grid_resolution * 100000;
    code = x_code + y_code;
    return code;
};

string g_GetCellCode(double x, double y, double grid_resolution, double left, double top)
{
    std::string s("ABCDEFGHIJKLMNOPQRSTUVWXYZ");
    std::string str_letter;
    std::string code;

    __int64 xi;
    xi = floor(x / grid_resolution) - floor(left / grid_resolution);

    __int64 yi;
    yi = ceil(top / grid_resolution) - floor(y / grid_resolution);

    int digit = (int)(xi / 26);
    if (digit >= 1)
        str_letter = s.at(digit % s.size());

    int reminder = xi - digit * 26;
    str_letter += s.at(reminder % s.size());

    std::string num_str = std::to_string(yi);

    code = str_letter + num_str;
        
    return code;

}

class CDemand_Period {
public:
    CDemand_Period() : demand_period{ 0 }, starting_time_slot_no{ 0 }, ending_time_slot_no{ 0 }, m_RandomSeed {101}
    {
    }

    int get_time_horizon_in_min()
    {
        return (ending_time_slot_no - starting_time_slot_no) * 15;
    }

    unsigned int m_RandomSeed;

    float GetRandomRatio()
    {
        //m_RandomSeed is automatically updated.
        m_RandomSeed = (LCG_a * m_RandomSeed + LCG_c) % LCG_M;

        return float(m_RandomSeed) / LCG_M;
    }

    void compute_cumulative_profile(int starting_slot_no, int ending_slot_no)
    {
        float total_ratio = 0;
        for (int s = starting_slot_no; s < ending_slot_no; s++)
        {
            total_ratio += departure_time_ratio[s];
        }
        if (total_ratio < 0.000001)
            total_ratio = 0.000001;

        cumulative_departure_time_ratio[starting_slot_no] = 0;
        float cumulative_ratio = 0;
        for (int s = starting_slot_no; s < ending_slot_no; s++)
        {
             cumulative_ratio += departure_time_ratio[s]/ total_ratio;
             cumulative_departure_time_ratio[s] = cumulative_ratio;
             dtalog.output() << "cumulative profile ratio at slot  " << s << " = "  << cumulative_departure_time_ratio[s] << endl;
        }
        dtalog.output() << "final cumulative profile ratio" << cumulative_departure_time_ratio[ending_slot_no - 1] << endl;

    }
    int get_time_slot_no()
    {
        float r = GetRandomRatio();
         for (int s = starting_time_slot_no; s < ending_time_slot_no; s++)
        {
             if (r < cumulative_departure_time_ratio[s])
                 return s;
        }
         return starting_time_slot_no;  // first time slot as the default value
    }

    float departure_time_ratio[_MAX_TIMESLOT_PerPeriod];
    float cumulative_departure_time_ratio[_MAX_TIMESLOT_PerPeriod];

    string demand_period;
    int starting_time_slot_no;
    int ending_time_slot_no;
    string time_period;
    int demand_period_id;
};

class CAgent_type {
public:
    CAgent_type() : agent_type_no{ 1 }, value_of_time{ 1 }, time_headway_in_sec {1}
    {
    }

    int agent_type_no;
    // dollar per hour
    float value_of_time;
    // link type, product consumption equivalent used, for travel time calculation
    float PCE;
    float time_headway_in_sec;
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
        path_gradient_cost{ 0 }, path_gradient_cost_difference{ 0 }, path_gradient_cost_relative_difference{ 0 }, subarea_output_flag{1}, measurement_flag {0}
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
    int subarea_output_flag;
    int measurement_flag;
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
    std::map <int, bool> diverted_vehicle_map;


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
   std::vector <int> path_link_sequence;
};

class CColumnVector {

public:
    // this is colletion of unique paths
    CColumnVector() : cost{ 0 }, time{ 0 }, distance{ 0 }, od_volume{ 0 }, bfixed_route{ false }, m_passing_sensor_flag{ -1 }, information_type{ 0 }
    {
    }

    float cost;
    float time;
    float distance;
    // od volume
    double od_volume;
    bool bfixed_route;
    int information_type;
    

    int m_passing_sensor_flag;
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
    CAgent_Simu() : agent_vector_seq_no{ -1 }, path_toll{ 0 }, departure_time_in_min{ 0 }, m_bGenereated{ false }, m_bCompleteTrip{ false },
        path_travel_time_in_min{ 0 }, path_distance{ 0 }, diversion_flag{ 0 }, time_headway{ number_of_simu_interval_reaction_time }, PCE_unit_size{ 1 }
    {
    }

    ~CAgent_Simu()
    {
       
    }

    void AllocateMemory()
    {

            m_current_link_seq_no = 0;

            m_Veh_LinkArrivalTime_in_simu_interval.reserve(path_link_seq_no_vector.size());
            m_Veh_LinkDepartureTime_in_simu_interval.reserve(path_link_seq_no_vector.size());

            for (int i = 0; i < path_link_seq_no_vector.size(); ++i)
            {
                m_Veh_LinkArrivalTime_in_simu_interval.push_back(-1);
                m_Veh_LinkDepartureTime_in_simu_interval.push_back(-1);
            }

            m_path_link_seq_no_vector_size = path_link_seq_no_vector.size();
            departure_time_in_simu_interval = int(departure_time_in_min * 60.0 / number_of_seconds_per_interval + 0.5);  // round off
        
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
    int diversion_flag;
    bool m_bGenereated;
    bool m_bCompleteTrip;

    int agent_service_type;
    int demand_type;
    int agent_id;

    // for column pool index
    int at;
    int tau;
    int dest;

    int m_current_link_seq_no;
    int m_path_link_seq_no_vector_size;

    int departure_time_in_simu_interval;
    float arrival_time_in_min;
    float path_travel_time_in_min;
    float path_distance;

    unsigned int m_RandomSeed;

    // external input
    std::vector<int> path_link_seq_no_vector;
    std::vector<int>  m_Veh_LinkArrivalTime_in_simu_interval;
    std::vector<int>  m_Veh_LinkDepartureTime_in_simu_interval;
    int time_headway;  // in terms of simulation interval 
    int PCE_unit_size;  // the number of units: 
};

vector<CAgent_Simu*> g_agent_simu_vector;

class Assignment {
public:
    // default is UE
    Assignment() : assignment_mode{ 0 }, g_number_of_memory_blocks{ 8 }, g_number_of_threads{ 1 }, g_link_type_file_loaded{ true }, g_agent_type_file_loaded{ false },
        total_demand_volume{ 0.0 },  g_column_pool{ nullptr }, g_number_of_in_memory_simulation_intervals{ 500 },
        g_number_of_column_generation_iterations{ 20 }, g_number_of_demand_periods{ 24 },g_number_of_links{ 0 }, g_number_of_timing_arcs{ 0 },
        g_number_of_nodes{ 0 }, g_number_of_zones{ 0 }, g_number_of_agent_types{ 0 },  debug_detail_flag{ 1 }, path_output{ 0 }, trajectory_output{ 1 }, major_path_volume_threshold{ 6 }, trajectory_sampling_rate{ 1.0 }, td_link_performance_sampling_interval_in_min{ 60 }, td_link_performance_sampling_interval_hd_in_min{ 15 }, trajectory_diversion_only{ 0 }, m_GridResolution{ 0.01 }
    {
    }

    ~Assignment()
    {
        if (g_column_pool)
            Deallocate4DDynamicArray(g_column_pool, g_number_of_zones, g_number_of_zones, g_number_of_agent_types);


        DeallocateLinkMemory4Simulation();
    }

    void InitializeDemandMatrix(int number_of_zones, int number_of_agent_types, int number_of_time_periods)
    {
        total_demand_volume = 0.0;
        g_number_of_zones = number_of_zones;
        g_number_of_agent_types = number_of_agent_types;

        g_column_pool = Allocate4DDynamicArray<CColumnVector>(number_of_zones, number_of_zones, max(1, number_of_agent_types), number_of_time_periods);

        for (int i = 0; i < number_of_zones; ++i)
        {
            g_origin_demand_array[i] = 0.0;
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
    void STMesoTrafficSimulation();

    //OD demand estimation estimation
    void Demand_ODME(int OD_updating_iterations);
    void AllocateLinkMemory4Simulation();
    void UpdateRTPath(CAgent_Simu* pAgent);
    bool RTSP_RealTimeShortestPathFinding(int time_slot_no, int simu_interval_t);
    void DeallocateLinkMemory4Simulation();

    double m_GridResolution;
    int assignment_mode;
    int g_number_of_memory_blocks;
    int g_number_of_threads;
    int path_output;
    int trajectory_output;
    float trajectory_sampling_rate;
    int trajectory_diversion_only;
    int td_link_performance_sampling_interval_in_min;
    float td_link_performance_sampling_interval_hd_in_min;

    float major_path_volume_threshold;

    bool g_link_type_file_loaded;
    bool g_agent_type_file_loaded;

    float total_demand_volume;
    std::map<int, float> g_origin_demand_array;
    CColumnVector**** g_column_pool;

    // the data horizon in the memory
    int g_number_of_in_memory_simulation_intervals;
    int g_number_of_column_generation_iterations;
    int g_number_of_demand_periods;


    std::map<_int64, int> cell_id_mapping;  // this is used to mark if this cell_id has been identified or not;

    int g_number_of_links;
    int g_number_of_timing_arcs;
    int g_number_of_nodes;
    int g_number_of_zones;
    int g_number_of_agent_types;

    std::map<int, int> node_seq_no_2_info_zone_seq_no_mapping;  // this is used to mark if this zone_id has been identified or not
    std::map<int, int> zone_seq_no_2_info_mapping;  // this is used to mark if this zone_id has been identified or not

    int debug_detail_flag;

    // hash table, map external node number to internal node sequence no.
    std::map<int, int> g_node_id_to_seq_no_map;
    std::map<string, int> g_mvmt_key_to_link_no_map;
    // from integer to integer map zone_id to zone_seq_no
    std::map<int, int> g_zoneid_to_zone_seq_no_mapping;
    std::map<string, int> g_link_id_map;

    std::map<int, double> zone_id_X_mapping;
    std::map<int, double> zone_id_Y_mapping;

    std::vector<CDemand_Period> g_DemandPeriodVector;
    int g_LoadingStartTimeInMin;
    int g_LoadingEndTimeInMin;

    std::vector<CAgent_type> g_AgentTypeVector;
    std::map<int, CLinkType> g_LinkTypeMap;

    std::map<string, int> demand_period_to_seqno_mapping;
    std::map<string, int> agent_type_2_seqno_mapping;

    float total_demand[_MAX_AGNETTYPES][_MAX_TIMEPERIODS];
    float g_DemandGlobalMultiplier;

    // used in ST Simulation
    float** m_LinkOutFlowCapacity;  // per second interval for simplicity
    int** m_LinkOutFlowState;  // per second interval for simplicity


    // in min
    float** m_LinkTDWaitingTime;
    std::vector<float> m_LinkTotalWaitingTimeVector;;
    // number of simulation time intervals

    float** m_LinkCumulativeArrivalVector;
    float** m_LinkCumulativeDepartureVector;

    float* m_LinkCACount;  // CA, assign this value to m_LinkCumulativeArrivalVector at a given time in min
    float* m_LinkCDCount;  // CD

    int g_start_simu_interval_no;
    int g_number_of_simulation_intervals;
       // is shorter than g_number_of_simulation_intervals
    int g_number_of_loading_intervals_in_sec;
    // the data horizon in the memory in min
    int g_number_of_intervals_in_min;

    int g_number_of_intervals_in_sec;

};

Assignment assignment;

class CVDF_Period
{
public:
    CVDF_Period() : m{ 0.5 }, VOC{ 0 }, gamma{ 3.47f }, mu{ 1000 }, PHF{ 3 },
		alpha{ 0.15f }, beta{ 4 }, rho{ 1 }, preload{ 0 }, penalty{ 0 }, marginal_base{ 1 },
        starting_time_slot_no{ 0 }, ending_time_slot_no{ 0 },
        cycle_length{ -1 }, red_time{ 0 }, effective_green_time{ 0 }, t0{ 0 }, t3{ 0 }, start_green_time{ -1 }, end_green_time{ -1 }
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

        if (cycle_length > 1 )
        {
            if(volume > 10)
                avg_travel_time += cycle_length/60.0 / 2.0; // add additional random delay due to signalized intersection by 
            else
            {
                int debug = 1;
            }
        }
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
                    float t0_ph = t0 / (60.0 / MIN_PER_TIMESLOT);
                    float t2_ph = t2 / (60.0 / MIN_PER_TIMESLOT);
                    float t3_ph = t3 / (60.0 / MIN_PER_TIMESLOT);

                    //second congested phase based on the gamma calculated from the dynamic eq. (32)
                    Queue[time_abs] = 1 / (4.0 ) * gamma * (t_ph - t0_ph) * (t_ph - t0_ph) * (t_ph - t3_ph) * (t_ph - t3_ph);
                    // unit is hour
                    waiting_time[time_abs] = 1 / (4.0*mu) *gamma *(t_ph - t0_ph)*(t_ph - t0_ph) * (t_ph - t3_ph)*(t_ph - t3_ph);
                    arrival_rate[time_abs] = gamma * (t_ph - t0_ph)*(t_ph - t2_ph)*(t_ph - t3_ph) + mu;
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
    float effective_green_time;
    int start_green_time;
    int end_green_time;

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
        length{ 1 }, free_flow_travel_time_in_min{ 1 }, link_spatial_capacity{ 100 }, capacity_reduction_flag {0} ,
        timing_arc_flag{ false }, traffic_flow_code{ 0 }, spatial_capacity_in_vehicles{ 999999 }, link_type{ 2 }, subarea_id{ -1 }, RT_flow_volume{ 0 },
        cell_type{ -1 }, saturation_flow_rate { 1800}, TD_link_reduction_start_time_slot_no {99999}
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

            for(int at = 0; at < _MAX_AGNETTYPES; ++at)
                volume_per_period_per_at[tau][at] = 0;
        }

        for (int tau = 0; tau < _MAX_TIMESLOT_PerPeriod; ++tau)
        {
            RT_travel_time_vector[tau] = -1;
            RT_speed_vector[tau] = -1;

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
    void CalculateTD_RTVDFunction();


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
    int zone_seq_no_for_outgoing_connector;

    int number_of_lanes;
    double lane_capacity;
    double saturation_flow_rate;

    float TD_link_capacity[_MAX_TIMESLOT_PerPeriod];
    int TD_link_reduction_start_time_slot_no;
    std::map <int, bool> TD_link_closure_map;

	double length;
	double free_flow_travel_time_in_min;
	double free_speed;

	double cost;
	double link_spatial_capacity;

    bool timing_arc_flag;
    int traffic_flow_code;
    int spatial_capacity_in_vehicles;
    int time_to_be_released;

    // 1. based on BPR.

    int link_seq_no;
    int capacity_reduction_flag;
    string link_id;
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
    int link_type;

    int cell_type;
    string mvmt_txt_id;
    string path_code_str;
    string tmc_corridor_name;
    string link_type_name;
    string link_type_code;

    float PCE;
    float fftt;

    CVDF_Period VDF_period[_MAX_TIMEPERIODS];

	double TDBaseTT[_MAX_TIMEPERIODS];
	double TDBaseCap[_MAX_TIMEPERIODS];
	double TDBaseFlow[_MAX_TIMEPERIODS];
	double TDBaseQueue[_MAX_TIMEPERIODS];

    int type;

    //static
    //float flow_volume;
    //float travel_time;

    int subarea_id;
    double flow_volume_per_period[_MAX_TIMEPERIODS];
    double RT_flow_volume;
    double background_flow_volume_per_period[_MAX_TIMEPERIODS];

	double  volume_per_period_per_at[_MAX_TIMEPERIODS][_MAX_AGNETTYPES];

	double  queue_length_perslot[_MAX_TIMEPERIODS];  // # of vehicles in the vertical point queue
	double travel_time_per_period[_MAX_TIMEPERIODS];
    double RT_travel_time;

    float RT_travel_time_vector[_MAX_TIMESLOT_PerPeriod];  // for each 15 min time slot
    float RT_speed_vector[_MAX_TIMESLOT_PerPeriod];  // for each 15 min time slot


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
    int current_driving_AgentID;
    int win_count;
    int lose_count;

};


class CSignalTiming
{
public:
    CSignalTiming() : cycle_length{ 60 }, red_time{ 0 }, saturation_flow_rate{ 1800 }, VDF_capacity{ 1800 }, start_green_time{ 0 }, end_green_time{ 9999 }
    {
    }

    int cycle_length;
    float red_time;  // for VDF
    float saturation_flow_rate;
    float VDF_capacity;  // for VDF

    int link_seq_no;
    int start_green_time;
    int end_green_time;
};

class CNode
{
public:
    CNode() : zone_id{ -1 }, zone_org_id{ -1 }, prohibited_movement_size{ 0 }, node_seq_no{ -1 }, subarea_id {-1}, is_activity_node{ 0 }, is_information_zone{ 0 }
    {
    }

    //int accessible_node_count;

    int zone_id;
    __int64 cell_id;
    string cell_str;
    // original zone id for non-centriod nodes
    int zone_org_id;
    int subarea_id;
    int prohibited_movement_size;
    // sequence number
    int node_seq_no;

    //external node number
    int node_id;

    int is_activity_node;
    int is_information_zone;

    double x;
    double y;

    std::vector<int> m_outgoing_link_seq_no_vector;
    std::vector<int> m_incoming_link_seq_no_vector;

    std::vector<int> m_to_node_seq_no_vector;
    std::map<int, int> m_to_node_2_link_seq_no_map;

    std::map<string, int> m_prohibited_movement_string_map;
};

class CInfoCell {
public: 
    __int64 cell_id;
    string cell_str;
    std::vector<GDPoint> m_ShapePoints;

    void CreateCell(double x, double y, double grid_resolution)
    {

        __int64 xi;
        xi = floor(x / grid_resolution);

        __int64 yi;
        yi = floor(y / grid_resolution);

        double left, right, top, bottom;
        
        left = xi * grid_resolution;
        right = (xi + 1) * grid_resolution;

        top = (yi +1)* grid_resolution;
        bottom = (yi) * grid_resolution;

        GDPoint	pt0, pt1, pt2, pt3, pt4;

        pt0.x = left; 	pt0.y = top;
        pt1.x = right; 	pt1.y = top;
        pt2.x = right; 	pt2.y = bottom;
        pt3.x = left; 	pt3.y = bottom;
        pt4.x = left; 	pt4.y = top;

        
        m_ShapePoints.push_back(pt0);
        m_ShapePoints.push_back(pt1);
        m_ShapePoints.push_back(pt2);
        m_ShapePoints.push_back(pt3);
        m_ShapePoints.push_back(pt4);

    }

};

class CTMC_Corridor_Info {
public:
    CTMC_Corridor_Info() 
    {
    }

    void reset()
    {
        total_VMT = 0;
        total_VHT = 0;
        total_VDT = 0;
        lowest_speed = 9999;
        highest_speed = -1;
        link_count = 0;
    }

    double get_avg_speed()
    {
        return total_VMT / max(0.001,total_VHT);  //miles per hour
    }
    double total_VMT;
    double total_VHT;
    double total_VDT;

    double avg_speed;
    double lowest_speed;
    double highest_speed;
    double total_congestion_duration;
    int link_count;

    std::map<int, int> road_sequence_map;

};

std::vector<CNode> g_node_vector;
std::vector<CLink> g_link_vector;
std::vector<CSignalTiming> g_signal_timing_arc_vector;
std::map<string, CTMC_Corridor_Info> g_tmc_corridor_vector;
std::map<string, CInfoCell> g_info_cell_map;


class COZone
{
public:
    COZone() : obs_production{ -1 }, obs_attraction{ -1 },
        est_production{ -1 }, est_attraction{ -1 },
        est_production_dev{ 0 }, est_attraction_dev{ -1 }, total_demand{ 0 }, b_real_time_information {false}
    {
    }

    bool b_real_time_information;

    float obs_production;
    float obs_attraction;

    float est_production;
    float est_attraction;

    float est_production_dev;
    float est_attraction_dev;

    // 0, 1,
    int zone_seq_no;
    // external zone id // this is origin zone
    int zone_id;
    int node_seq_no;
    float total_demand;
    float obs_production_upper_bound_flag;
    float obs_attraction_upper_bound_flag;

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
    NetworkForSP() : temp_path_node_vector_size{ _MAX_LINK_SIZE_IN_A_PATH }, m_value_of_time{ 10 }, bBuildNetwork{ false }, m_memory_block_no{ 0 }, m_agent_type_no{ 0 }, m_tau{ 0 }, b_real_time_information{ false }
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

    bool b_real_time_information; 
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
            m_link_genalized_cost_array[i] = pLink->travel_time_per_period[m_tau] + pLink->VDF_period[m_tau].penalty + pLink->VDF_period[m_tau].toll[agent_type_no] / m_value_of_time * 60 + pLink->RT_travel_time;  // *60 as 60 min per hour
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
                if(g_link_vector[link_seq_no].AllowAgentType (p_assignment->g_AgentTypeVector[m_agent_type_no].agent_type,m_tau))
                {
                    m_outgoing_link_seq_no_vector[outgoing_link_size] = link_seq_no;
                    m_to_node_seq_no_vector[outgoing_link_size] = g_node_vector[i].m_to_node_seq_no_vector[j];

                    outgoing_link_size++;

                    if (outgoing_link_size >= _MAX_LINK_SIZE_FOR_A_NODE)
                    {
                        dtalog.output() << " Error: outgoing_link_size >= _MAX_LINK_SIZE_FOR_A_NODE" << endl;
                        g_ProgramStop();
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

            for (int i = 0; i < g_node_vector.size(); ++i)
            {
                if (g_node_vector[i].zone_org_id > 0) // for each physical node
                { // we need to make sure we only create two way connectors between nodes and zones
                    dtalog.output() << "node id= " << g_node_vector[i].node_id << " with zone id " << g_node_vector[i].zone_org_id << "and "
                                    << NodeForwardStarArray[i].OutgoingLinkSize << " outgoing links." << endl;

                    for (int j = 0; j < NodeForwardStarArray[i].OutgoingLinkSize; j++)
                    {
                        int link_seq_no = NodeForwardStarArray[i].OutgoingLinkNoArray[j];
                        dtalog.output() << "  outgoing node = " << g_node_vector[g_link_vector[link_seq_no].to_node_seq_no].node_id << endl;
                    }
                }
                else
                {
                    if (dtalog.debug_level() == 3)
                    {
                        dtalog.output() << "node id= " << g_node_vector[i].node_id << " with "
                                        << NodeForwardStarArray[i].OutgoingLinkSize << " outgoing links." << endl;

                        for (int j = 0; j < NodeForwardStarArray[i].OutgoingLinkSize; ++j)
                        {
                            int link_seq_no = NodeForwardStarArray[i].OutgoingLinkNoArray[j];
                            dtalog.output() << "  outgoing node = " << g_node_vector[g_link_vector[link_seq_no].to_node_seq_no].node_id << endl;
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
        int origin_zone = m_origin_zone_seq_no_vector[o_node_index]; // assigned nodes for computing


        if (p_assignment->g_number_of_nodes >= 1000 && origin_zone%97 == 0)
            dtalog.output() << "label correcting for zone " << origin_zone <<  " in processor " << processor_id <<  endl;

        if (dtalog.debug_level() >= 2)
            dtalog.output() << "SP iteration k =  " << iteration_k << ": origin node: " << g_node_vector[origin_node].node_id << endl;

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
        while (!(m_ListFront == -1))   //SEList_empty()
        {
            // from_node = SEList_front();
            // SEList_pop_front();  // remove current node FromID from the SE list

            from_node = m_ListFront;//pop a node FromID for scanning
            tempFront = m_ListFront;
            m_ListFront = m_SENodeList[m_ListFront];
            m_SENodeList[tempFront] = -1;

            m_node_status_array[from_node] = 2;

            if (dtalog.log_path() >= 2)
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

                if (dtalog.log_path() >= 2)
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
                            continue;
                        }
                    }
                }

                //remark: the more complicated implementation can be found in paper Shortest Path Algorithms In Transportation Models: Classical and Innovative Aspects
                //	A note on least time path computation considering delays and prohibitions for intersection movements

                if (m_link_outgoing_connector_zone_seq_no_array[link_sqe_no] >= 0)
                {
                    if(m_link_outgoing_connector_zone_seq_no_array[link_sqe_no] != origin_zone)
                    {
                        // filter out for an outgoing connector with a centriod zone id different from the origin zone seq no
                        continue;
                    }
                }

                //very important: only origin zone can access the outbound connectors,
                //the other zones do not have access to the outbound connectors

                // Mark				new_time = m_label_time_array[from_node] + pLink->travel_time_per_period[tau];
                // Mark				new_distance = m_label_distance_array[from_node] + pLink->length;
                float additional_cost = 0;

                if (g_link_vector[link_sqe_no].RT_travel_time > 1)  // used in real time routing only
                {
                    additional_cost = g_link_vector[link_sqe_no].RT_travel_time;

                    //if (g_link_vector[link_sqe_no].RT_travel_time > 999)
                    //    continue; //skip this link due to closure
                }

                   
                new_to_node_cost = m_node_label_cost[from_node] + m_link_genalized_cost_array[link_sqe_no] + additional_cost;

                if (dtalog.log_path())
                {
                    dtalog.output() << "SP:  checking from node " << g_node_vector[from_node].node_id
                                    << "  to node" << g_node_vector[to_node].node_id << " cost = " << new_to_node_cost << endl;
                }

                if (new_to_node_cost < m_node_label_cost[to_node]) // we only compare cost at the downstream node ToID at the new arrival time t
                {
                    if (dtalog.log_path())
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

                    if (dtalog.log_path())
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

        if (dtalog.log_path())
        {
            dtalog.output() << "SPtree at iteration k = " << iteration_k <<  " origin node: "
                            << g_node_vector[origin_node].node_id  << endl;

            //Initialization for all non-origin nodes
            for (int i = 0; i < p_assignment->g_number_of_nodes; ++i)
            {
                int node_pred_id = -1;
                int node_pred_no = m_node_predecessor[i];

                if (node_pred_no >= 0)
                    node_pred_id = g_node_vector[node_pred_no].node_id;

                if(m_node_label_cost[i] < 9999)
                {
                    dtalog.output() << "SP node: " << g_node_vector[i].node_id << " label cost " << m_node_label_cost[i] << "time "
                                    << m_label_time_array[i] << "node_pred_id " << node_pred_id << endl;
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
std::vector<NetworkForSP*> g_NetworkForRTSP_vector;
NetworkForSP g_RoutingNetwork;


void g_ReadDemandFileBasedOnDemandFileList(Assignment& assignment)
{
    //	fprintf(g_pFileOutputLog, "number of zones =,%lu\n", g_zone_vector.size());

    assignment.InitializeDemandMatrix(g_zone_vector.size(), assignment.g_AgentTypeVector.size(), assignment.g_DemandPeriodVector.size());

    float total_demand_in_demand_file = 0;

    CCSVParser parser;
    dtalog.output() << endl;
    dtalog.output() << "Step 1.8: Reading file section [demand_file_list] in setting.csv..." << endl;
    parser.IsFirstLineHeader = false;

    if (parser.OpenCSVFile("settings.csv", false))
    {
        while (parser.ReadRecord_Section())
        {
            if (parser.SectionName == "[demand_file_list]")
            {
                int file_sequence_no = 1;

                string format_type = "null";

                int demand_format_flag = 0;

                if (!parser.GetValueByFieldName("file_sequence_no", file_sequence_no))
                    break;

                // skip negative sequence no
                if (file_sequence_no <= -1)
                    continue;

                double loading_scale_factor = 1.0;
                string file_name, demand_period, agent_type;
                parser.GetValueByFieldName("file_name", file_name);
                parser.GetValueByFieldName("demand_period", demand_period);
                parser.GetValueByFieldName("format_type", format_type);
                parser.GetValueByFieldName("loading_scale_factor", loading_scale_factor, false);



                if (format_type.find("null") != string::npos)  // skip negative sequence no
                {
                    dtalog.output() << "Please provide format_type in section [demand_file_list.]" << endl;
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
                        g_ProgramStop();
                    }
                }

                if (demand_period_no > _MAX_TIMEPERIODS)
                {
                    dtalog.output() << "demand_period_no should be less than settings in demand_period section. Please change the parameter settings in the source code." << endl;
                    g_ProgramStop();
                }

                if (format_type.find("column") != string::npos)  // or muliti-column
                {
                    bool bFileReady = false;
                    int error_count = 0;
                    int critical_OD_count = 0;
                    double critical_OD_volume = 0;

                    // read the file formaly after the test.
                    FILE* st;
                    fopen_ss(&st, file_name.c_str(), "r");
                    if (st)
                    {
                        bFileReady = true;
                        int line_no = 0;

                        dtalog.output() << endl << "VIODP  o,d,volume,geometry" << endl;

                        while (true)
                        {
                            int origin_zone = (int)(g_read_float(st));
                            int destination_zone = (int)g_read_float(st);
                            float demand_value = g_read_float(st);

                            if (origin_zone <= -1)
                            {
                                if (line_no == 1 && !feof(st))  // read only one line, but has not reached the end of the line
                                {
                                    dtalog.output() << endl << "Error: Only one line has been read from file. Are there multiple columns of demand type in file " << file_name << " per line?" << endl;
                                    g_ProgramStop();
                                }
                                break;
                            }

                            if (assignment.g_zoneid_to_zone_seq_no_mapping.find(origin_zone) == assignment.g_zoneid_to_zone_seq_no_mapping.end())
                            {
                                if (error_count < 10)
                                    dtalog.output() << endl << "Warning: origin zone " << origin_zone << "  has not been defined in node.csv" << endl;

                                error_count++;
                                // origin zone has not been defined, skipped.
                                continue;
                            }

                            if (assignment.g_zoneid_to_zone_seq_no_mapping.find(destination_zone) == assignment.g_zoneid_to_zone_seq_no_mapping.end())
                            {
                                if (error_count < 10)
                                    dtalog.output() << endl << "Warning: destination zone " << destination_zone << "  has not been defined in node.csv" << endl;

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

                            demand_value *= loading_scale_factor;
                            if (demand_value >= 5)
                            {
                                critical_OD_volume += demand_value ;
                                critical_OD_count += 1;
                                //dtalog.output() << origin_zone << "," << destination_zone << "," << demand_value << "," << "\"LINESTRING( " <<
                                //    assignment.zone_id_X_mapping[origin_zone] << " " << assignment.zone_id_Y_mapping[origin_zone] << "," <<
                                //    assignment.zone_id_X_mapping[destination_zone] << " " << assignment.zone_id_Y_mapping[destination_zone] << ")\" " << endl;

                            }

                            assignment.total_demand[agent_type_no][demand_period_no] += demand_value;
                            assignment.g_column_pool[from_zone_seq_no][to_zone_seq_no][agent_type_no][demand_period_no].od_volume += demand_value;
                            assignment.total_demand_volume += demand_value;
                            assignment.g_origin_demand_array[from_zone_seq_no] += demand_value;

                            // we generate vehicles here for each OD data line
                            if (line_no <= 5)  // read only one line, but has not reached the end of the line
                                dtalog.output() << "o_zone_id:" << origin_zone << ", d_zone_id: " << destination_zone << ", value = " << demand_value << endl;

                            line_no++;
                        }  // scan lines

                        fclose(st);

                        dtalog.output() << "total demand volume is " << assignment.total_demand_volume << endl;
                        dtalog.output() << "crtical demand volume has " << critical_OD_count << " OD pairs in size," << critical_OD_volume << ", " << ", account for " << critical_OD_volume / max(0.1, assignment.total_demand_volume) * 100 << "%%" << endl;

                        dtalog.output() << "crtical OD zones volume has " << critical_OD_count << " OD pairs in size," << critical_OD_volume << ", " << ", account for " << critical_OD_volume / max(0.1, assignment.total_demand_volume) * 100 << "%%" << endl;


                        std::map<int, float>::iterator it;
                        int count_zone_demand = 0;
                        for (it = assignment.g_origin_demand_array.begin(); it != assignment.g_origin_demand_array.end(); ++it)
                        {
                            //if (it->second > 5)
                            //{
                            //    dtalog.output() << "o_zone " << it->first << ", d_zone=," << it->second << endl;
                            //    count_zone_demand++;
                            //}
                        }
                        dtalog.output() << "There are  " << count_zone_demand << " zones with positive demand" << endl;

                    }
                    else
                    {
                        // open file
                        dtalog.output() << "Error: File " << file_name << " cannot be opened.\n It might be currently used and locked by EXCEL." << endl;
                        g_ProgramStop();
                    }
                }
                else if (format_type.compare("path") == 0)
                {

                    int path_counts = 0;
                    float sum_of_path_volume = 0;
                    CCSVParser parser;
                    if (parser.OpenCSVFile(file_name, false))
                    {
                        int total_path_in_demand_file = 0;
                        // read agent file line by line,

                        int agent_id, o_zone_id, d_zone_id;
                        string agent_type, demand_period;

                        std::vector <int> node_sequence;

                        while (parser.ReadRecord())
                        {
                            total_path_in_demand_file++;
                            if (total_path_in_demand_file % 1000 == 0)
                                dtalog.output() << "total_path_in_demand_file is " << total_path_in_demand_file << endl;

                            parser.GetValueByFieldName("agent_id", agent_id);
                            parser.GetValueByFieldName("o_zone_id", o_zone_id);
                            parser.GetValueByFieldName("d_zone_id", d_zone_id); 

                            CAgentPath agent_path_element;

                            agent_path_element.path_id = 0;
                            parser.GetValueByFieldName("path_id", agent_path_element.path_id, false);


                            int from_zone_seq_no = 0;
                            int to_zone_seq_no = 0;
                            from_zone_seq_no = assignment.g_zoneid_to_zone_seq_no_mapping[o_zone_id];
                            to_zone_seq_no = assignment.g_zoneid_to_zone_seq_no_mapping[d_zone_id];

                            if (format_type.compare("path") == 0)
                            {
                                double volume = 0;
                                parser.GetValueByFieldName("volume", volume);
                                volume *= loading_scale_factor;
                                agent_path_element.volume = volume;
                                path_counts++;
                                sum_of_path_volume += agent_path_element.volume;

                                assignment.total_demand[agent_type_no][demand_period_no] += agent_path_element.volume;
                                assignment.g_column_pool[from_zone_seq_no][to_zone_seq_no][agent_type_no][demand_period_no].od_volume += agent_path_element.volume;
                                assignment.total_demand_volume += agent_path_element.volume;
                                assignment.g_origin_demand_array[from_zone_seq_no] += agent_path_element.volume;
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
                        dtalog.output() << "total_demand_volume loaded from path file is " << sum_of_path_volume << " with " << path_counts << "paths." << endl;

                    }
                    else
                    {
                        //open file
                        dtalog.output() << "Error: File " << file_name << " cannot be opened.\n It might be currently used and locked by EXCEL." << endl;
                        g_ProgramStop();
                    }
                }
                else if (format_type.compare("matrix") == 0)
                {
                    bool bFileReady = false;
                    int error_count = 0;
                    int critical_OD_count = 0;
                    double critical_OD_volume = 0;

                    vector<int> LineIntegerVector;

                    CCSVParser parser;
                    parser.IsFirstLineHeader = false;
                    if (parser.OpenCSVFile(file_name, true))
                    {
                        int control_type_code;
                        int i = 0;
                        if (parser.ReadRecord())
                        {
                            parser.ConvertLineStringValueToIntegers();
                            LineIntegerVector = parser.LineIntegerVector;
                        }
                        parser.CloseCSVFile();
                    }

                    int number_of_zones = LineIntegerVector.size();


                    bFileReady = false;
                    int i;

                    FILE* st;
                    fopen_s(&st, file_name.c_str(), "r");
                    if (st != NULL)
                    {
                        // read the first line
                        g_read_a_line(st);

                        cout << "number of zones to be read = " << number_of_zones << endl;

                        //test if a zone has been defined. 
                        for (int destination_zone_index = 0; destination_zone_index < number_of_zones; destination_zone_index++)
                        {
                            int zone = LineIntegerVector[destination_zone_index];
                            if (assignment.g_zoneid_to_zone_seq_no_mapping.find(zone) == assignment.g_zoneid_to_zone_seq_no_mapping.end())
                            {
                                if (error_count < 10)
                                    dtalog.output() << endl << "Warning: destination zone " << zone << "  has not been defined in node.csv" << endl;

                                error_count++;
                                // destination zone has not been defined, skipped.
                                continue;
                            }

                        }


                        int line_no = 0;
                        for (int origin_zone_index = 0; origin_zone_index < number_of_zones; origin_zone_index++)
                        {
                            int origin_zone = (int)(g_read_float(st)); // read the origin zone number

                            if (assignment.g_zoneid_to_zone_seq_no_mapping.find(origin_zone) == assignment.g_zoneid_to_zone_seq_no_mapping.end())
                            {
                                if (error_count < 10)
                                    dtalog.output() << endl << "Warning: destination zone " << origin_zone << "  has not been defined in node.csv" << endl;

                                error_count++;
                                // destination zone has not been defined, skipped.
                                continue;
                            }

                            cout << "Reading file no." << file_sequence_no << " " << file_name << " at zone " << origin_zone << " ... " << endl;

                            for (int destination_zone_index = 0; destination_zone_index < number_of_zones; destination_zone_index++)
                            {
                                int destination_zone = LineIntegerVector[destination_zone_index];

                                float demand_value = g_read_float(st);

                                int from_zone_seq_no = 0;
                                int to_zone_seq_no = 0;
                                from_zone_seq_no = assignment.g_zoneid_to_zone_seq_no_mapping[origin_zone];
                                to_zone_seq_no = assignment.g_zoneid_to_zone_seq_no_mapping[destination_zone];

                                // encounter return
                                if (demand_value < -99)
                                    break;

                                demand_value *= loading_scale_factor;
                                if (demand_value >= 1)
                                {
                                    critical_OD_volume += demand_value ;
                                    critical_OD_count += 1;
                                    //dtalog.output() << origin_zone << "," << destination_zone << "," << demand_value << "," << "\"LINESTRING( " <<
                                    //    assignment.zone_id_X_mapping[origin_zone] << " " << assignment.zone_id_Y_mapping[origin_zone] << "," <<
                                    //    assignment.zone_id_X_mapping[destination_zone] << " " << assignment.zone_id_Y_mapping[destination_zone] << ")\" " << endl;

                                }

                                assignment.total_demand[agent_type_no][demand_period_no] += demand_value;
                                assignment.g_column_pool[from_zone_seq_no][to_zone_seq_no][agent_type_no][demand_period_no].od_volume += demand_value;
                                assignment.total_demand_volume += demand_value;
                                assignment.g_origin_demand_array[from_zone_seq_no] += demand_value;

                                // we generate vehicles here for each OD data line
                                if (line_no <= 5)  // read only one line, but has not reached the end of the line
                                    dtalog.output() << "o_zone_id:" << origin_zone << ", d_zone_id: " << destination_zone << ", value = " << demand_value << endl;

                                line_no++;
                            }  // scan lines

                        }

                            fclose(st);

                            dtalog.output() << "total demand volume is " << assignment.total_demand_volume << endl;
                            dtalog.output() << "crtical demand volume has " << critical_OD_count << " OD pairs in size," << critical_OD_volume << ", " << ", account for " << critical_OD_volume / max(0.1, assignment.total_demand_volume) * 100 << "%%" << endl;

                            dtalog.output() << "crtical OD zones volume has " << critical_OD_count << " OD pairs in size," << critical_OD_volume << ", " << ", account for " << critical_OD_volume / max(0.1, assignment.total_demand_volume) * 100 << "%%" << endl;


                            std::map<int, float>::iterator it;
                            int count_zone_demand = 0;
                            for (it = assignment.g_origin_demand_array.begin(); it != assignment.g_origin_demand_array.end(); ++it)
                            {
                                if (it->second > 0.001)
                                {
                                    dtalog.output() << "o_zone " << it->first << ", demand=," << it->second << endl;
                                    count_zone_demand++;
                                }
                            }
                            dtalog.output() << "There are  " << count_zone_demand << " zones with positive demand" << endl;

                        
                    } //end reading file
                    else
                    {
                        // open file
                        dtalog.output() << "Error: File " << file_name << " cannot be opened.\n It might be currently used and locked by EXCEL." << endl;
                        g_ProgramStop();
                    }

                }
                else{
                    dtalog.output() << "Error: format_type = " << format_type << " is not supported. Currently STALite supports format such as column, matrix and path." << endl;
                    g_ProgramStop();
                }
            }
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
                parser.GetValueByFieldName("path_output", assignment.path_output, false, false);
                parser.GetValueByFieldName("major_path_volume_threshold", assignment.major_path_volume_threshold, false, false);
                parser.GetValueByFieldName("trajectory_output", assignment.trajectory_output, false, false);
                parser.GetValueByFieldName("trajectory_sampling_rate", assignment.trajectory_sampling_rate, false, false);
                parser.GetValueByFieldName("trajectory_diversion_only", assignment.trajectory_diversion_only, false, false);
                parser.GetValueByFieldName("td_link_performance_sampling_interval_in_min", assignment.td_link_performance_sampling_interval_in_min, false, false);
                parser.GetValueByFieldName("td_link_performance_sampling_interval_hd_in_min", assignment.td_link_performance_sampling_interval_hd_in_min, false, false);

                dtalog.output() << "td_link_performance_sampling_interval_in_min= " << assignment.td_link_performance_sampling_interval_in_min << " min" << endl;
                dtalog.output() << "td_link_performance_sampling_interval_hd_in_min= " << assignment.td_link_performance_sampling_interval_hd_in_min << " min" << endl;

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


double g_CalculateP2PDistanceInLonglatFromLatitudeLongitude(double p1x, double p1y, double p2x, double p2y)
{
    double distance = sqrt((p2y - p1y) * (p2y - p1y) + (p2x - p1x) * (p2x - p1x));
    return distance;
}


double g_CalculateP2PDistanceInMeterFromLatitudeLongitude(double p1x, double p1y, double p2x, double p2y)
{
    double PI = 3.1415926;
    double Equatorial_Radius = 3963.19059 * 1609; // unit: mile-> meter
    double toradians = 3.1415926 / 180.0;
    double todeg = 180.0 / PI;

    double p2lat = p2x * toradians;
    double p2lng = p2y * toradians;

    double p1lat = p1x * toradians;
    double p1lng = p1y * toradians;

    double distance = acos(sin(p1lat) * sin(p2lat) + cos(p1lat) * cos(p2lat) * cos(p2lng - p1lng)) * Equatorial_Radius;  // unit: mile
    return distance;
}

double g_CheckActivityNodes(Assignment& assignment)
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

            if (i % sampling_rate == 0)
            {
                g_node_vector[i].is_activity_node = 10;//random generation
                activity_node_count++;
            }
        }

        //if (activity_node_count <= 1)
        //{
        //    activity_node_count = 0;
        //    sampling_rate = 2;

        //    for (int i = 0; i < g_node_vector.size(); i++)
        //    {

        //        if (i % sampling_rate == 0)
        //        {
        //            g_node_vector[i].is_activity_node = 10;//random generation
        //            activity_node_count++;
        //        }
        //    }
        //     still no activity nodes, define all nodes as activity nodes
        //    if (activity_node_count <= 1)
        //    {
        //        activity_node_count = 0;

        //        for (int i = 0; i < g_node_vector.size(); i++)
        //        {

        //            g_node_vector[i].is_activity_node = 10;//random generation
        //            activity_node_count++;
        //        }
        //    }
        //}


    }


    // calculate avg near by distance; 
    double total_near_by_distance = 0;
    activity_node_count = 0;
    for (int i = 0; i < g_node_vector.size(); i++)
    {
        double min_near_by_distance = 100;
        if (g_node_vector[i].is_activity_node)
        {
            activity_node_count++;
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

    double nearby_distance = total_near_by_distance / max(1, activity_node_count);
    return nearby_distance;

}


void g_OutputModelFiles(int mode)
{
    if (mode == 1)
    {
        FILE* g_pFileModelNode = fopen("model_node.csv", "w");

        if (g_pFileModelNode != NULL)
        {
            fprintf(g_pFileModelNode, "node_id,node_no,activity_node_flag,zone_id,cell_id,cell_code,info_zone_flag,x_coord,y_coord\n");
            for (int i = 0; i < g_node_vector.size(); i++)
            {
                if (g_node_vector[i].node_id >= 0)
                {
                    fprintf(g_pFileModelNode, "%d,%d,%d,%d,%ld,%s,%d,%f,%f\n",
                        g_node_vector[i].node_id,
                        g_node_vector[i].node_seq_no,
                        g_node_vector[i].is_activity_node,
                        g_node_vector[i].zone_org_id,
                        g_node_vector[i].cell_id,
                        g_node_vector[i].cell_str.c_str(),
                        g_node_vector[i].is_information_zone,
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

    }

    if (mode == 2)
    {
        FILE* g_pFileModelLink = fopen("model_link.csv", "w");

        if (g_pFileModelLink != NULL)
        {
            fprintf(g_pFileModelLink, "link_id,link_no,from_node_id,to_node_id,link_type,link_type_name,nlanes,free_speed,capacity,capacity_reduction,geometry\n");

            //VDF_fftt1,VDF_cap1,VDF_alpha1,VDF_beta1
            for (int i = 0; i < g_link_vector.size(); i++)
            {
                if (g_link_vector[i].link_type >= 0)
                {
                    fprintf(g_pFileModelLink, "%s,%d,%d,%s,%d,%d,%d,%f,%f,%d",
                        g_link_vector[i].link_id.c_str(),
                        g_link_vector[i].link_seq_no,
                        g_node_vector[g_link_vector[i].from_node_seq_no].node_id,
                        g_node_vector[g_link_vector[i].to_node_seq_no].node_id,
                        g_link_vector[i].link_type,
                        g_link_vector[i].link_type_name.c_str(),
                        g_link_vector[i].number_of_lanes,
                        g_link_vector[i].free_speed,
                        g_link_vector[i].lane_capacity
                        //g_link_vector[i].VDF_period[0].FFTT,
                        //g_link_vector[i].VDF_period[0].capacity,
                        //g_link_vector[i].VDF_period[0].alpha,
                        //g_link_vector[i].VDF_period[0].beta,
                    );

                    fprintf(g_pFileModelLink, "\"%s\",\n", g_link_vector[i].geometry.c_str());
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

    if (mode == 3)  // cell
    {
        FILE* g_pFileZone = nullptr;
        g_pFileZone = fopen("model_cell.csv", "w");

        if (g_pFileZone == NULL)
        {
            cout << "File model_cell.csv cannot be opened." << endl;
            g_ProgramStop();
        }
        else
        {


            fprintf(g_pFileZone, "cell_code,geometry\n");

            std::map<string, CInfoCell>::iterator it;

            for (it = g_info_cell_map.begin(); it != g_info_cell_map.end(); ++it)
            {
 
                fprintf(g_pFileZone, "%s,", it->first.c_str());
                fprintf(g_pFileZone, "\"LINESTRING (");

                for (int s = 0; s < it->second.m_ShapePoints.size(); s++)
                {
                    fprintf(g_pFileZone, "%f %f,", it->second.m_ShapePoints[s].x, it->second.m_ShapePoints[s].y);
                }

                fprintf(g_pFileZone, ")\"");
                fprintf(g_pFileZone, "\n");
            }
            fclose(g_pFileZone);
        }
    }
}

void g_InfoGridGeneration(Assignment& assignment)
{
    dtalog.output() << "Step QEM mode for creating node 2 zone mapping" << endl;

    double activity_nearbydistance = g_CheckActivityNodes(assignment);
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

    //if (activity_nearbydistance * 4 < temp_resolution)
    //{
    //    temp_resolution = activity_nearbydistance * 4;

    //}


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

    for (int i = 0; i < g_node_vector.size(); i++)
    {
        g_node_vector[i].cell_id = g_GetCellID(g_node_vector[i].x, g_node_vector[i].y, assignment.m_GridResolution);

        g_node_vector[i].cell_str = g_GetCellCode(g_node_vector[i].x, g_node_vector[i].y, assignment.m_GridResolution, left, top);

        if (g_info_cell_map.find(g_node_vector[i].cell_str) == g_info_cell_map.end())
        {
            CInfoCell cell;

            cell.cell_str = g_node_vector[i].cell_str;
            cell.CreateCell(g_node_vector[i].x, g_node_vector[i].y, assignment.m_GridResolution);
            g_info_cell_map[g_node_vector[i].cell_str] = cell;
        }
    }


 /*   assignment.zone_id_2_node_no_mapping.clear();*/
    dtalog.output() << "Step 1.4.2: Grid Resolution " << assignment.m_GridResolution << endl;

    //int activity_node_count = 0;
    //for (int i = 0; i < g_node_vector.size(); i++)
    //{

    //    if (g_node_vector[i].is_activity_node >= 1)
    //    {

    //        //if (g_node_vector[i].node_id == 966)
    //        //{
    //        //    int itest = 1;
    //        //}
    //        __int64 cell_id = g_GetCellID(g_node_vector[i].x, g_node_vector[i].y, assignment.m_GridResolution);
    //        int zone_id;

    //        //if (assignment.cell_id_mapping.find(cell_id) == assignment.cell_id_mapping.end())  // create a cell
    //        //{
    //        //    create zone
    //        //    assignment.cell_id_mapping[cell_id] = g_node_vector[i].node_id;


    //        //    dtalog.output() << "Step 1.2: creating cell " << cell_id << " using node id " << g_node_vector[i].node_id << endl;

    //        //    zone_id = assignment.cell_id_mapping[cell_id]; // which is the node id when a cell is created. 
    //        //    if (assignment.zone_id_2_node_no_mapping.find(zone_id) == assignment.zone_id_2_node_no_mapping.end()) // create a zone 
    //        //    {
    //        //        dtalog.output() << "Step 1.2: creating zone " << zone_id << " using node id " << g_node_vector[i].node_id << endl;
    //        //        create zone
    //        //        assignment.zone_id_2_node_no_mapping[zone_id] = i;
    //        //        assignment.zone_id_2_cell_id_mapping[zone_id] = cell_id;
    //        //        g_node_vector[i].zone_org_id = zone_id;
    //        //    }
    //        //}
    //        //else
    //        //{
    //        //    zone_id = assignment.cell_id_mapping[cell_id]; // which is the node id when a cell is created. 
    //        //     for physcial nodes because only centriod can have valid zone_id.
    //        //    g_node_vector[i].zone_org_id = zone_id;

    //        //}

    //        activity_node_count++;


    //    }
    //}

    //dtalog.output() << "Step 1.4.3: creating " << assignment.zone_id_2_node_no_mapping.size() << " zones." << " # of activity nodes =" << activity_node_count << endl;
}
void g_InfoZoneGeneration(Assignment& assignment)
{
    dtalog.output() << "Step QEM mode for creating node 2 zone mapping" << endl;

    double activity_nearbydistance = g_CheckActivityNodes(assignment);
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

    //if (activity_nearbydistance * 4 < temp_resolution)
    //{
    //    temp_resolution = activity_nearbydistance * 4;

    //}


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

    for (int i = 0; i < g_node_vector.size(); i++)
    {
        g_node_vector[i].cell_id = g_GetCellID(g_node_vector[i].x, g_node_vector[i].y, assignment.m_GridResolution);
    }


    /*   assignment.zone_id_2_node_no_mapping.clear();*/
    dtalog.output() << "Step 1.4.2: Grid Resolution " << assignment.m_GridResolution << endl;

    //int activity_node_count = 0;
    //for (int i = 0; i < g_node_vector.size(); i++)
    //{

    //    if (g_node_vector[i].is_activity_node >= 1)
    //    {

    //        //if (g_node_vector[i].node_id == 966)
    //        //{
    //        //    int itest = 1;
    //        //}
    //        __int64 cell_id = g_GetCellID(g_node_vector[i].x, g_node_vector[i].y, assignment.m_GridResolution);
    //        int zone_id;

    //        //if (assignment.cell_id_mapping.find(cell_id) == assignment.cell_id_mapping.end())  // create a cell
    //        //{
    //        //    create zone
    //        //    assignment.cell_id_mapping[cell_id] = g_node_vector[i].node_id;


    //        //    dtalog.output() << "Step 1.2: creating cell " << cell_id << " using node id " << g_node_vector[i].node_id << endl;

    //        //    zone_id = assignment.cell_id_mapping[cell_id]; // which is the node id when a cell is created. 
    //        //    if (assignment.zone_id_2_node_no_mapping.find(zone_id) == assignment.zone_id_2_node_no_mapping.end()) // create a zone 
    //        //    {
    //        //        dtalog.output() << "Step 1.2: creating zone " << zone_id << " using node id " << g_node_vector[i].node_id << endl;
    //        //        create zone
    //        //        assignment.zone_id_2_node_no_mapping[zone_id] = i;
    //        //        assignment.zone_id_2_cell_id_mapping[zone_id] = cell_id;
    //        //        g_node_vector[i].zone_org_id = zone_id;
    //        //    }
    //        //}
    //        //else
    //        //{
    //        //    zone_id = assignment.cell_id_mapping[cell_id]; // which is the node id when a cell is created. 
    //        //     for physcial nodes because only centriod can have valid zone_id.
    //        //    g_node_vector[i].zone_org_id = zone_id;

    //        //}

    //        activity_node_count++;


    //    }
    //}

    //dtalog.output() << "Step 1.4.3: creating " << assignment.zone_id_2_node_no_mapping.size() << " zones." << " # of activity nodes =" << activity_node_count << endl;
}



void g_ReadInputData(Assignment& assignment)
{
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

                    if (assignment.g_LoadingEndTimeInMin < assignment.g_LoadingStartTimeInMin)
                    {
                        assignment.g_LoadingEndTimeInMin = assignment.g_LoadingStartTimeInMin + 1; // in case user errror
                    }
                    //g_fout << global_minute_vector[0] << endl;
                    //g_fout << global_minute_vector[1] << endl;


                    char time_interval_field_name[20];

                    for (int s = max(0,demand_period.starting_time_slot_no-1); s <= demand_period.ending_time_slot_no; s++)
                    {
                        demand_period.cumulative_departure_time_ratio[s] = 0;
                    }
                     

                    for (int s = demand_period.starting_time_slot_no; s <= demand_period.ending_time_slot_no; s++)
                    {
                        sprintf(time_interval_field_name, "time_interval%d", s);
                        demand_period.departure_time_ratio[s] = 1.0 / (demand_period.ending_time_slot_no - demand_period.starting_time_slot_no + 1);
                        parser_demand_period.GetValueByFieldName(time_interval_field_name, demand_period.departure_time_ratio[s],false);
                    }

                    demand_period.compute_cumulative_profile(demand_period.starting_time_slot_no, demand_period.ending_time_slot_no);

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
        dtalog.output() << "Error: File settings.csv cannot be opened.\n It might be currently used and locked by EXCEL." << endl;
        g_ProgramStop();
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
                parser_agent_type.GetValueByFieldName("headway", agent_type.time_headway_in_sec, false, false);
                assignment.agent_type_2_seqno_mapping[agent_type.agent_type] = assignment.g_AgentTypeVector.size();

                assignment.g_AgentTypeVector.push_back(agent_type);
                assignment.g_number_of_agent_types = assignment.g_AgentTypeVector.size();
            }
        }
        parser_agent_type.CloseCSVFile();

        if (assignment.g_AgentTypeVector.size() == 0 )
            dtalog.output() << "Error: Section agent_type does not contain information." << endl;
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

    std::map<int, int> zone_id_to_centriod_node_id_mapping;  // this is an one-to-one mapping
    std::map<int, int> zone_id_mapping;  // this is used to mark if this zone_id has been identified or not
    std::map<int, int> info_zone_id_2_node_id_mapping;  // this is used to mark if this zone_id has been identified or not



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

            int subarea_id = -1;
            parser.GetValueByFieldName("subarea_id", subarea_id,false);
            node.subarea_id = subarea_id;
            // this is an activity node // we do not allow zone id of zero
            if(zone_id>=1)
            {
                // for physcial nodes because only centriod can have valid zone_id.
                node.zone_org_id = zone_id;
                if (zone_id_mapping.find(zone_id) == zone_id_mapping.end())
                {
                    //create zone
                    zone_id_mapping[zone_id] = node_id;

                    assignment.zone_id_X_mapping[zone_id] = node.x;
                    assignment.zone_id_Y_mapping[zone_id] = node.y;
                }

                // for od calibration, I think we don't need to implement for now
                if (assignment.assignment_mode == 5)
                {
                    float production = 0;
                    float attraction = 0;
                    parser.GetValueByFieldName("production", production);
                    parser.GetValueByFieldName("attraction", attraction);

                    zone_id_production[zone_id] = production;
                    zone_id_attraction[zone_id] = attraction;
                }
            }
            int info_zone_id = -1;
            parser.GetValueByFieldName("info_zone_id", info_zone_id,false);

            if(info_zone_id>=1)
            {
            info_zone_id = _INFO_ZONE_ID + node_id;
            info_zone_id_2_node_id_mapping[info_zone_id] = node_id;
            node.is_information_zone = 1;
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
        dtalog.output() << "number of zones = " << zone_id_mapping.size() << endl;
        dtalog.output() << "number of info zones = " << info_zone_id_2_node_id_mapping.size() << endl;

    	// fprintf(g_pFileOutputLog, "number of nodes =,%d\n", assignment.g_number_of_nodes);
        parser.CloseCSVFile();
    }

    g_InfoGridGeneration(assignment);
    //g_InfoZoneMapping(assignment);
    g_OutputModelFiles(1);  // node
    g_OutputModelFiles(3);  // info cell
    // initialize zone vector
    dtalog.output() << "Step 1.5: Initializing O-D zone vector..." << endl;

    std::map<int, int>::iterator it;

    for (it = zone_id_mapping.begin(); it != zone_id_mapping.end(); ++it)
    {
        COZone ozone;

        // for each zone, we have to also create centriod
        ozone.zone_id = it->first;  // zone_id
        ozone.zone_seq_no = g_zone_vector.size();
        ozone.obs_production = zone_id_production[it->first];
        ozone.obs_attraction = zone_id_attraction[it->first];

        assignment.g_zoneid_to_zone_seq_no_mapping[ozone.zone_id] = ozone.zone_seq_no;  // create the zone id to zone seq no mapping

        // create a centriod
        CNode node;
        // very large number as a special id
        node.node_id = -1* ozone.zone_id;
        node.node_seq_no = g_node_vector.size();
        assignment.g_node_id_to_seq_no_map[node.node_id] = node.node_seq_no;
        node.zone_id = ozone.zone_id;
        // push it to the global node vector
        g_node_vector.push_back(node);
        assignment.g_number_of_nodes++;

        ozone.node_seq_no = node.node_seq_no;
        // this should be the only one place that defines this mapping
        zone_id_to_centriod_node_id_mapping[ozone.zone_id] = node.node_id;
        // add element into vector
        g_zone_vector.push_back(ozone);
    }

    // information zone 
    for (it = info_zone_id_2_node_id_mapping.begin(); it != info_zone_id_2_node_id_mapping.end(); ++it)
    {
        COZone ozone;

        // for each zone, we have to also create centriod
        ozone.zone_id = it->first;  // zone_id
        ozone.zone_seq_no = g_zone_vector.size();

        assignment.g_zoneid_to_zone_seq_no_mapping[ozone.zone_id] = ozone.zone_seq_no;  // create the zone id to zone seq no mapping
        ozone.node_seq_no = assignment.g_node_id_to_seq_no_map[it->second];
        // this should be the only one place that defines this mapping
        zone_id_to_centriod_node_id_mapping[ozone.zone_id] = it->second;
        // add element into vector
        ozone.b_real_time_information = true;

        //establish mapping relationship, for future use in RT simulation
        assignment.node_seq_no_2_info_zone_seq_no_mapping[ozone.node_seq_no] = ozone.zone_seq_no;
        assignment.zone_seq_no_2_info_mapping[ozone.zone_seq_no] = true;

        g_zone_vector.push_back(ozone);
    }

    dtalog.output() << "adding "<< info_zone_id_2_node_id_mapping.size() << "info zones " << endl;

   // gravity model.
    if (assignment.assignment_mode == 5)
    {
        dtalog.output() << "writing demand.csv.." << endl;

        FILE* g_pFileODMatrix = nullptr;
        fopen_ss(&g_pFileODMatrix, "demand.csv", "w");

        if (!g_pFileODMatrix)
        {
            dtalog.output() << "File demand.csv cannot be opened." << endl;
            g_ProgramStop();
        }
        else
        {
            fprintf(g_pFileODMatrix, "o_zone_id,d_zone_id,volume\n");

            float total_attraction = 0;

            for (int d = 0; d < g_zone_vector.size(); ++d)
            {
                if(g_zone_vector[d].obs_attraction>0)
                    total_attraction += g_zone_vector[d].obs_attraction;
            }

            // reset the estimated production and attraction
            for (int orig = 0; orig < g_zone_vector.size(); ++orig)  // o
            {
                if(g_zone_vector[orig].obs_production >=0)
                {
                    for (int dest = 0; dest < g_zone_vector.size(); ++dest)  // d
                    {
                        if (g_zone_vector[dest].obs_attraction > 0)
                        {
                            float value = g_zone_vector[orig].obs_production * g_zone_vector[dest].obs_attraction / max(0.0001f, total_attraction);
                            fprintf(g_pFileODMatrix, "%d,%d,%.4f,\n", g_zone_vector[orig].zone_id, g_zone_vector[dest].zone_id, value);
                            dtalog.output() << "orig= " << g_zone_vector[orig].zone_id << " dest= " << g_zone_vector[dest].zone_id << ":" <<  value << endl;
                        }
                    }
                }
            }

            fclose(g_pFileODMatrix);
        }
    }

    dtalog.output() << "number of zones = " << g_zone_vector.size() << endl;
    // step 4: read link file

    CCSVParser parser_link;

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

            if (assignment.g_link_id_map.find(linkID) != assignment.g_link_id_map.end())
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
            parser_link.GetValueByFieldName("mvmt_txt_id", movement_str, false);
            int cell_type = -1;
            if(parser_link.GetValueByFieldName("cell_type", cell_type, false)==true)
                link.cell_type = cell_type;




            parser_link.GetValueByFieldName("geometry", link.geometry,false);
            parser_link.GetValueByFieldName("path_code", link.path_code_str, false);
            parser_link.GetValueByFieldName("tmc_corridor_name", link.tmc_corridor_name, false);
            parser_link.GetValueByFieldName("link_type_name", link.link_type_name, false);
            parser_link.GetValueByFieldName("link_type_code", link.link_type_code, false);

            // and valid
            if (movement_str.size() > 0)
            {
                int main_node_id = -1;


                link.mvmt_txt_id = movement_str;
                link.main_node_id = main_node_id;
            }

            // Peiheng, 05/13/21, if setting.csv does not have corresponding link type or the whole section is missing, set it as 2 (i.e., Major arterial)
            int link_type = 2;
            parser_link.GetValueByFieldName("link_type", link_type, false);

            if (assignment.g_LinkTypeMap.find(link_type) == assignment.g_LinkTypeMap.end())
            {
                dtalog.output() << "link type " << link_type << " in link.csv is not defined for link " << from_node_id << "->"<< to_node_id << " in link_type.csv" << endl;
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

			double lane_capacity = 1800;
            parser_link.GetValueByFieldName("length", length);
            if (length < 0.007)
            {
                length = 0.007;  // minimum length
            }
            parser_link.GetValueByFieldName("free_speed", free_speed);

            if (free_speed <= 0.1)
                free_speed = 60;

            free_speed = max(0.1, free_speed);

			link.free_speed = free_speed;
            


            int number_of_lanes = 1;
            parser_link.GetValueByFieldName("lanes", number_of_lanes);
            parser_link.GetValueByFieldName("capacity", lane_capacity);
            
            link.free_flow_travel_time_in_min = length / free_speed * 60;
            link.traffic_flow_code = assignment.g_LinkTypeMap[link.link_type].traffic_flow_code;

            //spatial queue and kinematic wave
            link.spatial_capacity_in_vehicles = max(1.0,length * number_of_lanes * k_jam);

            // kinematic wave
            if (link.traffic_flow_code == 3)
                link.BWTT_in_simulation_interval = length / bwtt_speed *3600/ number_of_seconds_per_interval;

            // Peiheng, 02/03/21, useless block
            if (linkID == "10")
                int i_debug = 1;

            char VDF_field_name[20];

            for (int at = 0; at < assignment.g_AgentTypeVector.size(); at++)
            {
                double pce_at = 1; // default
                sprintf(VDF_field_name, "VDF_pce%s", assignment.g_AgentTypeVector[at].agent_type.c_str());

                parser_link.GetValueByFieldName(VDF_field_name, pce_at, false, true);

                if (pce_at > 1.001)  // log
                {
                    //dtalog.output() << "link " << from_node_id << "->" << to_node_id << " has a pce of " << pce_at << " for agent type "
                    //    << assignment.g_AgentTypeVector[at].agent_type.c_str() << endl;
                }


                for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
                {
                    link.VDF_period[tau].pce[at] = pce_at;
                }

            }

            // for traffic simulation
            for (int s = 0; s < _MAX_TIMESLOT_PerPeriod; s++)
            {
                link.TD_link_capacity[s] = lane_capacity * number_of_lanes;;
            }

            for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
            {
                //setup default values
                link.VDF_period[tau].capacity = lane_capacity * number_of_lanes;
                link.VDF_period[tau].FFTT = length / free_speed * 60.0;  // 60.0 for 60 min per hour
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

					if(link.VDF_period[tau].toll[at] >0.001)
					{ 
					dtalog.output() << "link " << from_node_id << "->" << to_node_id << " has a toll of " << link.VDF_period[tau].toll[at] << " for agent type "
						<< assignment.g_AgentTypeVector[at].agent_type.c_str() << " at demand period " << demand_period_id <<  endl;
					}
				}

				sprintf(VDF_field_name, "VDF_penalty%d", demand_period_id);
				parser_link.GetValueByFieldName(VDF_field_name, link.VDF_period[tau].penalty, false, false);

                if (link.cell_type >= 1) // micro lane-changing arc
                {
                    // additinonal min: 24 seconds 0.4 min
                    link.VDF_period[tau].penalty += 0.4;
                }

				sprintf(VDF_field_name, "VDF_PHF%d", demand_period_id);
                parser_link.GetValueByFieldName(VDF_field_name, link.VDF_period[tau].PHF, false, false);

                sprintf(VDF_field_name, "VDF_mu%d", demand_period_id);
                parser_link.GetValueByFieldName(VDF_field_name, link.VDF_period[tau].mu, false, false);  // mu should be per hour per link, so that we can calculate congestion duration and D/mu in BPR-X

                parser_link.GetValueByFieldName("cycle_length", link.VDF_period[tau].cycle_length,false,false);

                if(link.VDF_period[tau].cycle_length>=1)
                {
                    link.timing_arc_flag = true;

                parser_link.GetValueByFieldName("start_green_time", link.VDF_period[tau].start_green_time);
                parser_link.GetValueByFieldName("end_green_time", link.VDF_period[tau].end_green_time);
                parser_link.GetValueByFieldName("red_time", link.VDF_period[tau].red_time,false);
                parser_link.GetValueByFieldName("green_time", link.VDF_period[tau].effective_green_time,false);
                }

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
            link.link_spatial_capacity = k_jam * number_of_lanes*length;

            link.length = max(0.00001,length);
            for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
                link.travel_time_per_period[tau] = length / free_speed * 60;

            // min // calculate link cost based length and speed limit // later we should also read link_capacity, calculate delay

            //int sequential_copying = 0;
            //
            //parser_link.GetValueByFieldName("sequential_copying", sequential_copying);

            g_node_vector[internal_from_node_seq_no].m_outgoing_link_seq_no_vector.push_back(link.link_seq_no);  // add this link to the corresponding node as part of outgoing node/link
            g_node_vector[internal_to_node_seq_no].m_incoming_link_seq_no_vector.push_back(link.link_seq_no);  // add this link to the corresponding node as part of outgoing node/link

            g_node_vector[internal_from_node_seq_no].m_to_node_seq_no_vector .push_back(link.to_node_seq_no);  // add this link to the corresponding node as part of outgoing node/link
            g_node_vector[internal_from_node_seq_no].m_to_node_2_link_seq_no_map[link.to_node_seq_no] = link.link_seq_no;  // add this link to the corresponding node as part of outgoing node/link

            g_link_vector.push_back(link);

            string mvmt_key;
            parser_link.GetValueByFieldName("mvmt_key", mvmt_key, false);
            if(mvmt_key.size()>4) // main_node_id _ movement code
            {
            assignment.g_mvmt_key_to_link_no_map[mvmt_key] = assignment.g_number_of_links;
            }

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
            int node_id = zone_id_to_centriod_node_id_mapping[g_node_vector[i].zone_org_id];
            internal_to_node_seq_no = assignment.g_node_id_to_seq_no_map[node_id];
            zone_seq_no = assignment.g_zoneid_to_zone_seq_no_mapping[g_node_vector[i].zone_org_id];

            // incomming virtual connector
            g_AddNewVirtualConnectorLink(internal_from_node_seq_no, internal_to_node_seq_no, -1);
            // outgoing virtual connector
            g_AddNewVirtualConnectorLink(internal_to_node_seq_no, internal_from_node_seq_no, zone_seq_no);
        }
    }

    dtalog.output() << "number of links =" << assignment.g_number_of_links << endl;

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

    CCSVParser parser_movement;
    int prohibited_count = 0;

    if (parser_movement.OpenCSVFile("movement.csv", false))  // not required
    {
        while (parser_movement.ReadRecord())
        {
            string ib_link_id;
            int node_id = 0;
            string ob_link_id;
            int prohibited_flag  = 0;

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
}


//
//void g_reload_timing_arc_data(Assignment& assignment)
//{
//    dtalog.output() << "Step 1.7: Reading service arc in timing.csv..." << endl;
//
//    CCSVParser parser_timing_arc;
//    if (parser_timing_arc.OpenCSVFile("timing.csv", false))
//    {
//        while (parser_timing_arc.ReadRecord())  // if this line contains [] mark, then we will also read field headers.
//        {
//            string mvmt_key;
//            if (!parser_timing_arc.GetValueByFieldName("mvmt_key", mvmt_key))
//            {
//                dtalog.output() << "Error: mvmt_key in file timing.csv is not defined." << endl;
//                continue;
//            }
//            // create a link object
//            CSignalTiming timing_arc;
//
//            if (assignment.g_mvmt_key_to_link_no_map.find(mvmt_key) == assignment.g_mvmt_key_to_link_no_map.end())
//            {
//                dtalog.output() << "Error: mvmt_key " << mvmt_key << " in file timing.csv is not defined in link.csv." << endl;
//                //has not been defined
//                continue;
//            }
//            else
//            {
//                timing_arc.link_seq_no = assignment.g_mvmt_key_to_link_no_map[mvmt_key];
//                g_link_vector[timing_arc.link_seq_no].timing_arc_flag = true;
//            }
//
//            string time_period;
//            if (!parser_timing_arc.GetValueByFieldName("time_window", time_period))
//            {
//                dtalog.output() << "Error: Field time_window in file timing.csv cannot be read." << endl;
//                g_ProgramStop();
//                break;
//            }
//
//            vector<float> global_minute_vector;
//
//            //input_string includes the start and end time of a time period with hhmm format
//            global_minute_vector = g_time_parser(time_period); //global_minute_vector incldue the starting and ending time
//            if (global_minute_vector.size() == 2)
//            {
//                if (global_minute_vector[0] < assignment.g_LoadingStartTimeInMin)
//                    global_minute_vector[0] = assignment.g_LoadingStartTimeInMin;
//
//                if (global_minute_vector[0] > assignment.g_LoadingEndTimeInMin)
//                    global_minute_vector[0] = assignment.g_LoadingEndTimeInMin;
//
//                if (global_minute_vector[1] < assignment.g_LoadingStartTimeInMin)
//                    global_minute_vector[1] = assignment.g_LoadingStartTimeInMin;
//
//                if (global_minute_vector[1] > assignment.g_LoadingEndTimeInMin)
//                    global_minute_vector[1] = assignment.g_LoadingEndTimeInMin;
//
//                if (global_minute_vector[1] < global_minute_vector[0])
//                    global_minute_vector[1] = global_minute_vector[0];
//
//            }
//            else
//                continue;
//
//            float time_interval = 0;
//
//
//            // capacity in the space time arcs
//            float capacity = 1;
//            parser_timing_arc.GetValueByFieldName("capacity", capacity);
//            timing_arc.VDF_capacity = max(0.0f, capacity);
//
//            // capacity in the space time arcs
//            parser_timing_arc.GetValueByFieldName("cycle_length", timing_arc.cycle_length);
//
//            // capacity in the space time arcs
//            parser_timing_arc.GetValueByFieldName("red_time", timing_arc.red_time);
//            parser_timing_arc.GetValueByFieldName("start_green_time", timing_arc.start_green_time);
//            parser_timing_arc.GetValueByFieldName("end_green_time", timing_arc.end_green_time);
//
//            for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
//            {
//                    // to do: we need to consider multiple periods in the future, Xuesong Zhou, August 20, 2020.
//                g_link_vector[timing_arc.link_seq_no].VDF_period[tau].red_time = timing_arc.red_time;
//                g_link_vector[timing_arc.link_seq_no].VDF_period[tau].cycle_length = timing_arc.cycle_length;
//            }
//
//
//            g_signal_timing_arc_vector.push_back(timing_arc);
//            assignment.g_number_of_timing_arcs++;
//
//            if (assignment.g_number_of_timing_arcs % 10000 == 0)
//                dtalog.output() << "reading " << assignment.g_number_of_timing_arcs << " timing_arcs.. " << endl;
//        }
//
//        parser_timing_arc.CloseCSVFile();
//    }
//
//    dtalog.output() << endl;
//    dtalog.output() << "Step 1.8: Reading file section [demand_file_list] in setting.csv..." << endl;
//    // we now know the number of links
//    dtalog.output() << "number of timing records = " << assignment.g_number_of_timing_arcs << endl << endl;
//}

void g_load_scenario_data(Assignment& assignment)
{
    dtalog.output() << "Step 2.0: Reading scenario data..." << endl;

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
                CSignalTiming timing_arc;
                // map external node number to internal node seq no.
                int internal_from_node_seq_no = assignment.g_node_id_to_seq_no_map[from_node_id];
                int internal_to_node_seq_no = assignment.g_node_id_to_seq_no_map[to_node_id];

                if (g_node_vector[internal_from_node_seq_no].m_to_node_2_link_seq_no_map.find(internal_to_node_seq_no) != g_node_vector[internal_from_node_seq_no].m_to_node_2_link_seq_no_map.end())
                {
                    timing_arc.link_seq_no = g_node_vector[internal_from_node_seq_no].m_to_node_2_link_seq_no_map[internal_to_node_seq_no];
                    g_link_vector[timing_arc.link_seq_no].timing_arc_flag = true;
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

               }
                else
                    continue;

               // capacity in the space time arcs
                float capacity = 1;
                parser.GetValueByFieldName("capacity", capacity);
                timing_arc.VDF_capacity = max(0.0f, capacity);

                g_link_vector[timing_arc.link_seq_no].capacity_reduction_flag = 1;

                int info_zone_id = 1;
                parser.GetValueByFieldName("info_zone", info_zone_id,false);


                for (int s = global_minute_vector[0] / 15; s <= global_minute_vector[1] / 15; s++)
                {
                    g_link_vector[timing_arc.link_seq_no].TD_link_capacity[s] = capacity;
                    if (capacity < 1)
                    {
                        g_link_vector[timing_arc.link_seq_no].TD_link_closure_map[s] = true;
                    }
                    
                    
                    if (s < g_link_vector[timing_arc.link_seq_no].TD_link_reduction_start_time_slot_no)
                        g_link_vector[timing_arc.link_seq_no].TD_link_reduction_start_time_slot_no = s;

                }


                g_signal_timing_arc_vector.push_back(timing_arc);
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
            CColumnVector* p_column_pool;

            for (int orig = 0; orig < zone_size; ++orig)  // o
            {
                for (int dest = 0; dest < zone_size; ++dest) //d
                {
                    for (int tau = 0; tau < tau_size; ++tau)  //tau
                    {
                        p_column_pool = &(assignment.g_column_pool[orig][dest][at][tau]);
                        if (p_column_pool->od_volume > 0)
                        {

                            column_vector_size = p_column_pool->path_node_sequence_map.size();

                            it_begin = p_column_pool->path_node_sequence_map.begin();
                            it_end = p_column_pool->path_node_sequence_map.end();
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
                                    //#pragma omp critical
                                    {
                                        g_link_vector[link_seq_no].flow_volume_per_period[tau] += link_volume_contributed_by_path_volume * PCE_ratio;
                                        g_link_vector[link_seq_no].volume_per_period_per_at[tau][at] += link_volume_contributed_by_path_volume;  // pure volume, not consider PCE
                                    }
                                }

                                // this  self-deducting action does not agents with fixed routing policies.
                                if(!p_column_pool->bfixed_route && b_self_reducing_path_volume)
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
        //	float volume = assignment.m_LinkCumulativeDepartureVector[l][assignment.g_number_of_simulation_intervals - 1];  // link flow rates
        //	float waiting_time_count = 0;

        //for (int tt = 0; tt < assignment.g_number_of_simulation_intervals; tt++)
        //{
        //	waiting_time_count += assignment.m_LinkTDWaitingTime[l][tt/number_of_simu_intervals_in_min];   // tally total waiting cou
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
double g_reset_and_update_link_volume_based_on_ODME_columns(int number_of_links, int iteration_no, double& system_gap)
{
    float total_gap = 0;
    float sub_total_gap_link_count = 0;
    float sub_total_system_gap_count = 0;
    system_gap = 0;
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
        CColumnVector* p_column_pool;

        for (int orig = 0; orig < zone_size; ++orig)  // o
        {
            for (int dest = 0; dest < zone_size; ++dest) //d
            {
                for (int tau = 0; tau < tau_size; ++tau)  //tau
                {
                    p_column_pool = &(assignment.g_column_pool[orig][dest][at][tau]);
                    if (p_column_pool->od_volume > 0)
                    {
                        // continuous: type 0
                        column_vector_size = p_column_pool->path_node_sequence_map.size();

                        it_begin = p_column_pool->path_node_sequence_map.begin();
                        it_end = p_column_pool->path_node_sequence_map.end();
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

    int total_link_count = 0;

    // calcualte deviation for each measurement type
    for (int i = 0; i < number_of_links; ++i)
    {
        g_link_vector[i].CalculateTD_VDFunction();

        if (g_link_vector[i].obs_count >= 1)  // with data
        {
            int tau = 0;

                g_link_vector[i].est_count_dev = g_link_vector[i].flow_volume_per_period[tau] + g_link_vector[i].VDF_period[tau].preload - g_link_vector[i].obs_count;

            if (dtalog.debug_level() == 2)
            {
                dtalog.output() << "link " << g_node_vector [g_link_vector[i].from_node_seq_no].node_id
                                << "->" << g_node_vector[g_link_vector[i].to_node_seq_no].node_id
                                << "obs:, " << g_link_vector[i].obs_count << "est:, " << g_link_vector[i].flow_volume_per_period[tau]
                                << "dev:," << g_link_vector[i].est_count_dev << endl;
            }
            if(g_link_vector[i].upper_bound_flag==0)
            { 
            total_gap += abs(g_link_vector[i].est_count_dev);
            sub_total_gap_link_count += fabs(g_link_vector[i].est_count_dev / g_link_vector[i].obs_count);
            sub_total_system_gap_count += g_link_vector[i].est_count_dev / g_link_vector[i].obs_count;
            }
            else
            {  // upper bound constraints 
                if(g_link_vector[i].est_count_dev>0)
                { 
                total_gap += abs(g_link_vector[i].est_count_dev);
                sub_total_gap_link_count += fabs(g_link_vector[i].est_count_dev / g_link_vector[i].obs_count);
                sub_total_system_gap_count += g_link_vector[i].est_count_dev / g_link_vector[i].obs_count;
                }
            }
            total_link_count += 1;
        }
    }

    for (int orig = 0; orig < g_zone_vector.size(); ++orig)  // o
    {
        if (g_zone_vector[orig].obs_attraction >= 1)  // with observation
        {
            g_zone_vector[orig].est_attraction_dev = g_zone_vector[orig].est_attraction - g_zone_vector[orig].obs_attraction;

            if (dtalog.debug_level() == 2)
            {
                dtalog.output() << "zone " << g_zone_vector[orig].zone_id << "A: obs:" << g_zone_vector[orig].obs_attraction
                                << ",est:," << g_zone_vector[orig].est_attraction << ",dev:," << g_zone_vector[orig].est_attraction_dev << endl;
            }

            total_gap += abs(g_zone_vector[orig].est_attraction_dev);
            sub_total_gap_A_count += g_zone_vector[orig].est_attraction_dev / g_zone_vector[orig].obs_attraction;
        }

        if (g_zone_vector[orig].obs_production >= 1)  // with observation
        {
            g_zone_vector[orig].est_production_dev = g_zone_vector[orig].est_production - g_zone_vector[orig].obs_production;

            if (dtalog.debug_level() == 2)
            {
                dtalog.output() << "zone " << g_zone_vector[orig].zone_id << "P: obs:" << g_zone_vector[orig].obs_production
                                << ",est:," << g_zone_vector[orig].est_production << ",dev:," << g_zone_vector[orig].est_production_dev << endl;
            }

            total_gap += abs(g_zone_vector[orig].est_production_dev);
            sub_total_gap_P_count += g_zone_vector[orig].est_production_dev / g_zone_vector[orig].obs_production;
        }
    }

    dtalog.output() << "ODME #" << iteration_no/*<< " total abs gap= " << total_gap*/
        << " ,%link_MAPE: " << (sub_total_gap_link_count) / max(1, total_link_count) * 100 <<
        " ,%system_MPE: " << (sub_total_system_gap_count) / max(1, total_link_count) * 100 <<  endl;
    double gap = sub_total_gap_link_count / max(1, total_link_count);
    system_gap = sub_total_system_gap_count / max(1, total_link_count);

    return gap;
}

void g_update_gradient_cost_and_assigned_flow_in_column_pool(Assignment& assignment, int inner_iteration_number)
{
    double total_system_cost_gap = 0;
    float total_relative_gap = 0;
    double total_system_travel_cost = 0;

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
        CColumnVector* p_column_pool;
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
                    p_column_pool = &(assignment.g_column_pool[orig][dest][at][tau]);
                    if (p_column_pool->od_volume > 0)
                    {
                        column_vector_size = p_column_pool->path_node_sequence_map.size();

                        // scan through the map with different node sum for different paths
                        /// step 1: update gradient cost for each column path

                        least_gradient_cost = 999999;
                        least_gradient_cost_path_seq_no = -1;
                        least_gradient_cost_path_node_sum_index = -1;
                        path_seq_count = 0;

                        it_begin = p_column_pool->path_node_sequence_map.begin();
                        it_end = p_column_pool->path_node_sequence_map.end();
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
                                total_system_travel_cost += (it->second.path_gradient_cost * it->second.path_volume);
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

                                    total_system_cost_gap += (it->second.path_gradient_cost_difference * it->second.path_volume);
                                    total_system_travel_cost += (it->second.path_gradient_cost * it->second.path_volume);

                                    step_size = 1.0 / (inner_iteration_number + 2) * p_column_pool->od_volume;

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
                                p_column_pool->path_node_sequence_map[least_gradient_cost_path_node_sum_index].path_volume += total_switched_out_path_volume;
                                total_system_travel_cost += (p_column_pool->path_node_sequence_map[least_gradient_cost_path_node_sum_index].path_gradient_cost *
                                    p_column_pool->path_node_sequence_map[least_gradient_cost_path_node_sum_index].path_volume);
                            }
                        }
                    }
                }
            }
        }
    }

    dtalog.output() << "column updating: iteration= " << inner_iteration_number << ", total_gap=" << total_system_cost_gap
        << ",total_relative_gap=" << total_system_cost_gap / max(0.00001, total_system_travel_cost) << endl;

}

void g_column_pool_optimization(Assignment& assignment, int column_updating_iterations)
{
    // column_updating_iterations is internal numbers of column updating
    for (int n = 0; n < column_updating_iterations; ++n)
    {
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

void g_output_assignment_result(Assignment& assignment)
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

            // Option 2: BPR-X function
            fprintf(g_pFileLinkMOE, "link_id,from_node_id,to_node_id,link_type_name,link_type_code,time_period,volume,travel_time,speed,speed_ratio,VOC,capacity,queue,density,geometry,");

            //ODME
            if (assignment.assignment_mode == 3)
                fprintf(g_pFileLinkMOE, "obs_count,upper_type,dev,");

            fprintf(g_pFileLinkMOE, "notes\n");

            //Initialization for all nodes
            for (int i = 0; i < g_link_vector.size(); ++i)
            {
                // virtual connectors
                if (g_link_vector[i].link_type == -1)
                    continue;

                for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
                {
                    float speed = g_link_vector[i].free_speed;  // default speed 

                    if (g_link_vector[i].VDF_period[tau].avg_travel_time > 0.001f)
                        speed = g_link_vector[i].length / (g_link_vector[i].VDF_period[tau].avg_travel_time / 60.0);

                    float speed_ratio = speed/max(1,g_link_vector[i].free_speed);  // default speed 

                    fprintf(g_pFileLinkMOE, "%s,%d,%d,%s,%s,%s,%.3f,%.3f,%.3f,%.3f,%.1f,%.3f,0,0,\"%s\",",
                        g_link_vector[i].link_id.c_str(),
                        g_node_vector[g_link_vector[i].from_node_seq_no].node_id,
                        g_node_vector[g_link_vector[i].to_node_seq_no].node_id,
                        g_link_vector[i].link_type_name.c_str(),
                        g_link_vector[i].link_type_code.c_str(),
                        assignment.g_DemandPeriodVector[tau].time_period.c_str(),
                        g_link_vector[i].flow_volume_per_period[tau] + g_link_vector[i].VDF_period[tau].preload,
                        g_link_vector[i].VDF_period[tau].avg_travel_time,
                        speed,  /* 60.0 is used to convert min to hour */
                        speed_ratio,
                        g_link_vector[i].VDF_period[tau].VOC,
                        g_link_vector[i].VDF_period[tau].capacity,
                        g_link_vector[i].geometry.c_str());

                    if (assignment.assignment_mode == 3)  //ODME
                    {
                        if (g_link_vector[i].obs_count >= 1) //ODME
                            fprintf(g_pFileLinkMOE, "%.1f,%d,%.1f,", g_link_vector[i].obs_count, g_link_vector[i].upper_bound_flag, g_link_vector[i].est_count_dev);
                        else
                            fprintf(g_pFileLinkMOE, ",,");
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
        
        }
        fclose(g_pFileLinkMOE);
    }

    if (assignment.assignment_mode == 0 || assignment.path_output ==0)  //LUE
    {
        FILE* g_pFilePathMOE = nullptr;
        fopen_ss(&g_pFilePathMOE, "path.csv", "w");
        fclose(g_pFilePathMOE);

    }
    else if(assignment.assignment_mode >= 1)  //UE mode, or ODME, DTA
    {
        dtalog.output() << "writing path.csv.." << endl;

        float path_time_vector[_MAX_LINK_SIZE_IN_A_PATH];
        FILE* g_pFilePathMOE = nullptr;
        fopen_ss(&g_pFilePathMOE,"path.csv", "w");

        if (!g_pFilePathMOE)
        {
            dtalog.output() << "File path.csv cannot be opened." << endl;
            g_ProgramStop();
        }

       fprintf(g_pFilePathMOE, "agent_id,o_zone_id,d_zone_id,path_id,information_type,agent_type,demand_period,volume,assign_volume,subarea_flag,sensor_flag,toll,travel_time,VDF_travel_time,distance,path_code,node_sequence,link_sequence,time_sequence,geometry\n");

        int count = 1;

        clock_t start_t, end_t;
        start_t = clock();
        clock_t iteration_t;

        
        int agent_type_size = assignment.g_AgentTypeVector.size();
        int zone_size = g_zone_vector.size();
        int demand_period_size = assignment.g_DemandPeriodVector.size();

        CColumnVector* p_column_pool;

        float path_toll = 0;
        float path_distance = 0;
        float path_travel_time = 0;
        float path_delay = 0;
        float path_FF_travel_time = 0;
        float time_stamp = 0;

        std::map<int, CColumnPath>::iterator it, it_begin, it_end;

        if (assignment.major_path_volume_threshold > 0.00001)  // performing screening of path flow pattern
        {

            //initialization 
            bool b_subarea_mode = false;

            int number_of_links = g_link_vector.size();
            for (int i = 0; i < number_of_links; ++i)
            {
                for (int tau = 0; tau < demand_period_size; ++tau)
                {
                    // used in travel time calculation
                    g_link_vector[i].background_flow_volume_per_period[tau] = 0;
                }

                if (g_node_vector[g_link_vector[i].from_node_seq_no].subarea_id >= 1 && g_node_vector[g_link_vector[i].to_node_seq_no].node_id >= 1)
                {
             
                g_link_vector[i].subarea_id = g_node_vector[g_link_vector[i].from_node_seq_no].subarea_id;
                b_subarea_mode = true;
                }
                else
                    g_link_vector[i].subarea_id = 0;

            }


            /// <summary>  screening the path flow pattern
            for (int orig = 0; orig < zone_size; ++orig)
            {

                for (int at = 0; at < agent_type_size; ++at)
                {
                    for (int dest = 0; dest < zone_size; ++dest)
                    {
                        for (int tau = 0; tau < demand_period_size; ++tau)
                        {
                            p_column_pool = &(assignment.g_column_pool[orig][dest][at][tau]);
                            if (p_column_pool->od_volume > 0 )
                            {
                                // scan through the map with different node sum for different continuous paths
                                it_begin = p_column_pool->path_node_sequence_map.begin();
                                it_end = p_column_pool->path_node_sequence_map.end();

                                for (it = it_begin; it != it_end; ++it)
                                {
                                    int subarea_output_flag = 0;
                                    if (b_subarea_mode == true)
                                    {
                                        int insubarea_flag = 0;

                                        for (int nl = 0; nl < it->second.m_link_size; ++nl)  // arc a
                                        {
                                            int link_seq_no = it->second.path_link_vector[nl];

                                            if (g_link_vector[link_seq_no].subarea_id >= 1)
                                            {
                                                insubarea_flag = 1;
                                                break;
                                            }

                                        }
                                        // 
                                        if (insubarea_flag && it->second.path_volume > assignment.major_path_volume_threshold)
                                        {
                                            subarea_output_flag = 1;
                                        }

                                    }
                                    else
                                    {
                                        if (it->second.path_volume > assignment.major_path_volume_threshold)
                                            subarea_output_flag = 1;

                                    }
                                    if (subarea_output_flag==0)
                                    {
                                        it->second.subarea_output_flag = 0;  // disable the output of this column into path.csv

                                        for (int nl = 0; nl < it->second.m_link_size; ++nl)  // arc a
                                        {
                                            int link_seq_no = it->second.path_link_vector[nl];
                                            g_link_vector[link_seq_no].background_flow_volume_per_period[tau] += it->second.path_volume;
                                        }
                                    }

                                }
                            }
                        }
                        /// </summary>
                        /// <param name="assignment"></param>
                    }

                }
            }

            /// output background_link_volume.csv
            dtalog.output() << "writing link_performance.csv.." << endl;

            int b_debug_detail_flag = 0;
            FILE* g_pFileLinkMOE = nullptr;

            fopen_ss(&g_pFileLinkMOE, "link_background_volume.csv", "w");
            if (!g_pFileLinkMOE)
            {
                dtalog.output() << "File link_background_volume.csv cannot be opened." << endl;
                g_ProgramStop();
            }
            else
            {
                    fprintf(g_pFileLinkMOE, "link_id,from_node_id,to_node_id,from_cell_code,time_period,volume,background_volume,major_path_volume,ratio_of_major_path_flow,geometry,");

                    fprintf(g_pFileLinkMOE, "notes\n");

                    //Initialization for all nodes
                    for (int i = 0; i < g_link_vector.size(); ++i)
                    {
                        // virtual connectors
                        if (g_link_vector[i].link_type == -1)
                            continue;

                        for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
                        {
                            double volume = g_link_vector[i].flow_volume_per_period[tau] + g_link_vector[i].VDF_period[tau].preload;
                            double major_path_link_volume = g_link_vector[i].flow_volume_per_period[tau] + g_link_vector[i].VDF_period[tau].preload - g_link_vector[i].background_flow_volume_per_period[tau];
                            double ratio = major_path_link_volume / max(volume,0.000001);

                            if (volume < 0.0000001)
                                ratio = -1;
                            fprintf(g_pFileLinkMOE, "%s,%d,%d,%s,%s,%.3f,%.3f,%.3f,%.3f,\"%s\",",
                                g_link_vector[i].link_id.c_str(),
                                g_node_vector[g_link_vector[i].from_node_seq_no].node_id,
                                g_node_vector[g_link_vector[i].to_node_seq_no].node_id,
                                g_node_vector[g_link_vector[i].from_node_seq_no].cell_str.c_str(),
                                assignment.g_DemandPeriodVector[tau].time_period.c_str(),
                                g_link_vector[i].flow_volume_per_period[tau] + g_link_vector[i].VDF_period[tau].preload,
                                g_link_vector[i].background_flow_volume_per_period[tau],
                                major_path_link_volume,
                                ratio,
                                g_link_vector[i].geometry.c_str());
                            fprintf(g_pFileLinkMOE, "\n");

                        }

                    }

                    fclose(g_pFileLinkMOE);
                }
        

        } // end of path flow pattern screening 
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
                        p_column_pool = &(assignment.g_column_pool[orig][dest][at][tau]);
                        if (p_column_pool->od_volume > 0 || assignment.zone_seq_no_2_info_mapping.find(orig) != assignment.zone_seq_no_2_info_mapping.end())
                        {
                            
                            int information_type = 0;

                            if (assignment.zone_seq_no_2_info_mapping.find(orig) != assignment.zone_seq_no_2_info_mapping.end())
                            {
                                information_type = 1;
                                continue; // too many output 
                            }

                            time_stamp = (assignment.g_DemandPeriodVector[tau].starting_time_slot_no + assignment.g_DemandPeriodVector[tau].ending_time_slot_no) / 2.0 * MIN_PER_TIMESLOT;

                            // scan through the map with different node sum for different continuous paths
                            it_begin = p_column_pool->path_node_sequence_map.begin();
                            it_end = p_column_pool->path_node_sequence_map.end();

                            for (it = it_begin;it != it_end; ++it)
                            {
                                if (it->second.subarea_output_flag == 0)
                                    continue; 

                                if (count%100000 ==0)
                                {
                                    end_t = clock();
                                    iteration_t = end_t - start_t;
                                    dtalog.output() << "writing " << count/1000 << "K agents with CPU time " << iteration_t / 1000.0 << " s" << endl;
                                }

                                path_toll = 0;
                                path_distance = 0;
                                path_travel_time = 0;
                                path_FF_travel_time = 0;

                                path_time_vector[0] = time_stamp;
                                string path_code_str;

                                for (int nl = 0; nl < it->second.m_link_size; ++nl)  // arc a
                                {

                                    int link_seq_no = it->second.path_link_vector[nl];
                                    if (g_link_vector[link_seq_no].link_type >=0)
                                    {
                                    path_toll += g_link_vector[link_seq_no].VDF_period[tau].toll[at];
                                    path_distance += g_link_vector[link_seq_no].length;
                                    float link_travel_time = g_link_vector[link_seq_no].travel_time_per_period[tau];
                                    path_travel_time += link_travel_time;
                                    path_FF_travel_time += g_link_vector[link_seq_no].VDF_period[tau].FFTT;
                                    time_stamp += link_travel_time;
                                    path_time_vector[nl + 1] = time_stamp;

                                    path_code_str += g_link_vector[link_seq_no].path_code_str;
                                    }
                                    else
                                    {
                                        int virtual_link = 0;
                                    }
                                }

                                double total_agent_path_travel_time = 0;

                                for (int vi = 0; vi < it->second.agent_simu_id_vector.size(); ++vi)
                                {
                                    int agent_simu_id = it->second.agent_simu_id_vector[vi];
                                    CAgent_Simu* pAgentSimu = g_agent_simu_vector[agent_simu_id];
                                    total_agent_path_travel_time += pAgentSimu->path_travel_time_in_min;
                                }

                                double final_path_travel_time = path_travel_time;  // by default

                                if (it->second.agent_simu_id_vector.size() > 1)  // with simulated agents
                                {
                                    final_path_travel_time = total_agent_path_travel_time / it->second.agent_simu_id_vector.size();
                                }



                                int virtual_first_link_delta = 1;
                                int virtual_last_link_delta = 1;
                                // fixed routes have physical nodes always, without virtual connectors
                                if (p_column_pool->bfixed_route)
                                {
                                
                                    virtual_first_link_delta = 0;
                                    virtual_last_link_delta = 1;
                                }
                                if (information_type == 1 && it->second.path_volume < 0.1)
                                    continue;

                                if(information_type ==1)  // information diversion start from physical nodes
                                {
                                    virtual_first_link_delta = 0;
                                    virtual_last_link_delta = 1;

                                    it->second.path_volume = it->second.diverted_vehicle_map.size();
                                }


                                // assignment_mode = 1, path flow mode
                                if(assignment.assignment_mode >= 1 )
                                {

                                    fprintf(g_pFilePathMOE, "%d,%d,%d,%d,%d,%s,%s,%.2f,%d,%d,%d,%.1f,%.1f,%.4f,%.4f,%s,",
                                                         count,
                                                         g_zone_vector[orig].zone_id,
                                                         g_zone_vector[dest].zone_id,
                                                         it->second.path_seq_no,
                                                         information_type,
                                                         assignment.g_AgentTypeVector[at].agent_type.c_str(),
                                                         assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
                                                         it->second.path_volume,
                                                         it->second.agent_simu_id_vector.size(),
                                                         it->second.subarea_output_flag,
                                                         it->second.measurement_flag,
                                                         path_toll,
                                                         final_path_travel_time,
                                                         path_travel_time,
                                                         //path_FF_travel_time,
                                                         //final_path_travel_time- path_FF_travel_time,
                                                         path_distance, path_code_str.c_str());

                                    /* Format and print various data */
                                    for (int ni = 0+ virtual_first_link_delta; ni <it->second.m_node_size- virtual_last_link_delta; ++ni)
                                        fprintf(g_pFilePathMOE, "%d;", g_node_vector[it->second.path_node_vector[ni]].node_id);

                                    fprintf(g_pFilePathMOE, ",");
                                    int link_seq_no;
                                    for (int nl = 0 + virtual_first_link_delta; nl < it->second.m_link_size - virtual_last_link_delta; ++nl)
                                    {
                                        link_seq_no = it->second.path_link_vector[nl];
                                        fprintf(g_pFilePathMOE, "%s;", g_link_vector[link_seq_no].link_id.c_str());
                                    }
                                    fprintf(g_pFilePathMOE, ",");

                                    for (int nt = 0 + virtual_first_link_delta; nt < it->second.m_link_size+1 - virtual_last_link_delta; ++nt)
                                        fprintf(g_pFilePathMOE, "%s;", g_time_coding(path_time_vector[nt]).c_str());

                                    fprintf(g_pFilePathMOE, ",");

                                    //for (int nt = 0 + virtual_first_link_delta; nt < it->second.m_link_size+1 - virtual_last_link_delta; ++nt)
                                    //    fprintf(g_pFilePathMOE, "%.2f;", path_time_vector[nt]);

                                    //fprintf(g_pFilePathMOE, ",");

                                    //for (int nt = 0 + virtual_first_link_delta; nt < it->second.m_link_size - virtual_last_link_delta; ++nt)
                                    //    fprintf(g_pFilePathMOE, "%.2f;", path_time_vector[nt+1]- path_time_vector[nt]);

                                    //fprintf(g_pFilePathMOE, ",");

                                    fprintf(g_pFilePathMOE,  "\"LINESTRING (");

                                    for (int ni = 0 + virtual_first_link_delta; ni < it->second.m_node_size - virtual_last_link_delta; ++ni)
                                    {
                                        fprintf(g_pFilePathMOE, "%f %f", g_node_vector[it->second.path_node_vector[ni]].x,
                                                              g_node_vector[it->second.path_node_vector[ni]].y);

                                        if (ni != it->second.m_node_size - virtual_last_link_delta - 1)
                                            fprintf(g_pFilePathMOE, ", ");
                                    }

                                    fprintf(g_pFilePathMOE, ")\"\n");
                                    count++;
                                }
                               
                            }
                        }
                    }
                }
            }
        }
        fclose(g_pFilePathMOE);
    }
}

void g_output_TD_link_performance(Assignment& assignment, int output_mode = 1)
{
    dtalog.output() << "writing TD_link_performance.csv.." << endl;

    int b_debug_detail_flag = 0;
    FILE* g_pFileLinkMOE = nullptr;

    string file_name = "TD_link_performance.csv";

    if (output_mode == 0)
    {
        file_name = "TD_link_performance_hd.csv";
    }

    if (output_mode == 2)
    {
        file_name = "TD_link_performance_horizon.csv";
    }

    fopen_ss(&g_pFileLinkMOE, file_name.c_str(), "w");

    if (!g_pFileLinkMOE)
    {
        dtalog.output() << "File " << file_name.c_str() << "cannot be opened." << endl;
        g_ProgramStop();
    }
    else
    {

        // Option 2: BPR-X function
        fprintf(g_pFileLinkMOE, "link_id,tmc_corridor_name,link_type_name,from_node_id,to_node_id,from_cell_code,lanes,length,free_speed,FFTT,time_period,volume,CA,CD,density,queue_length,discharge_rate,TD_free_flow_travel_time,waiting_time_in_sec,RT_speed,speed,geometry,");
        fprintf(g_pFileLinkMOE, "notes\n");


        //Initialization for all nodes
        for (int i = 0; i < g_link_vector.size(); ++i)
        {
            // virtual connectors
            if (g_link_vector[i].link_type == -1)
                continue;

            // first loop for time t
            for (int t = 1; t < assignment.g_number_of_intervals_in_min; ++t)
            {
                if (assignment.td_link_performance_sampling_interval_in_min < 1)
                    assignment.td_link_performance_sampling_interval_in_min = 1;

                int sampling_time_interval = assignment.td_link_performance_sampling_interval_in_min;

                if (output_mode == 0)
                    sampling_time_interval = assignment.td_link_performance_sampling_interval_hd_in_min*60; // min by min

                if (output_mode == 2)
                    sampling_time_interval = assignment.g_number_of_loading_intervals_in_sec/60; // simulation horizon

                if (t % (sampling_time_interval ) == 0)
                {
                    int time_in_min = t ;  //relative time

                    float volume = 0;
                    float queue = 0;
                    float waiting_time_in_sec = 0;
                    int arrival_rate = 0;
                    float avg_waiting_time_in_sec = 0;

                    float travel_time = (float)(g_link_vector[i].free_flow_travel_time_in_min + avg_waiting_time_in_sec / 60.0);
                    float speed = g_link_vector[i].length / (g_link_vector[i].free_flow_travel_time_in_min / 60.0);
                    float virtual_arrival = 0;

                    float discharge_rate = g_link_vector[i].lane_capacity * g_link_vector[i].number_of_lanes * assignment.td_link_performance_sampling_interval_in_min / 60.0;

                    if (time_in_min >= 1)
                    {
                        volume = assignment.m_LinkCumulativeDepartureVector[i][t] - assignment.m_LinkCumulativeDepartureVector[i][t - sampling_time_interval ];


                        queue = assignment.m_LinkCumulativeArrivalVector[i][t] - assignment.m_LinkCumulativeDepartureVector[i][t];
                        //							waiting_time_in_min = queue / (max(1, volume));

                        float waiting_time_count = 0;

                        waiting_time_in_sec = assignment.m_LinkTDWaitingTime[i][t] * number_of_seconds_per_interval;

                        if (output_mode==2)
                        {
                            waiting_time_in_sec = assignment.m_LinkTotalWaitingTimeVector[i];
                            avg_waiting_time_in_sec = waiting_time_in_sec / max(1, arrival_rate);
                        }

                        arrival_rate = assignment.m_LinkCumulativeArrivalVector[i][t ] - assignment.m_LinkCumulativeArrivalVector[i][t- sampling_time_interval];
                        avg_waiting_time_in_sec = waiting_time_in_sec / max(1, arrival_rate);

                        travel_time = (float)(g_link_vector[i].free_flow_travel_time_in_min + avg_waiting_time_in_sec / 60.0);
                        speed = g_link_vector[i].length / (max(0.00001f, travel_time / 60.0));
                    }

                    if (speed >= 1000)
                    {
                        int i_debug = 1;
                    }
                    float density = (assignment.m_LinkCumulativeArrivalVector[i][t] - assignment.m_LinkCumulativeDepartureVector[i][t]) / (g_link_vector[i].length * g_link_vector[i].number_of_lanes);

                    if (output_mode == 0)
                    {
                        if(assignment.m_LinkCumulativeArrivalVector[i][t] < 1000)
                            continue; //skip
                    }

                    fprintf(g_pFileLinkMOE, "%s,%s,%s,%d,%d,%s,%d,%.3f,%.3f,%.3f,%s_%s,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,\"%s\",",
                        g_link_vector[i].link_id.c_str(),
                        g_link_vector[i].tmc_corridor_name.c_str(),
                        g_link_vector[i].link_type_name.c_str(),

                        g_node_vector[g_link_vector[i].from_node_seq_no].node_id,
                        g_node_vector[g_link_vector[i].to_node_seq_no].node_id,
                        g_node_vector[g_link_vector[i].from_node_seq_no].cell_str.c_str(),
                        g_link_vector[i].number_of_lanes,
                        g_link_vector[i].length,
                        g_link_vector[i].free_speed,
                        g_link_vector[i].free_flow_travel_time_in_min,

                        g_time_coding(assignment.g_LoadingStartTimeInMin + time_in_min - sampling_time_interval).c_str(),
                        g_time_coding(assignment.g_LoadingStartTimeInMin + time_in_min).c_str(),
                        volume,
                        assignment.m_LinkCumulativeArrivalVector[i][t],
                        assignment.m_LinkCumulativeDepartureVector[i][t],
                        density,
                        queue,
                        discharge_rate,
                        travel_time,
                        avg_waiting_time_in_sec,
                        speed,
                        g_link_vector[i].geometry.c_str());

                    fprintf(g_pFileLinkMOE, "simulation-based\n");
                }
            }  // for each time t
        }  // for each link l
        fclose(g_pFileLinkMOE);
    }//assignment mode 2 as simulation

}

void g_output_TD_link_state(Assignment& assignment, int output_mode = 1)
{
    dtalog.output() << "writing TD_link_state.csv.." << endl;

    int b_debug_detail_flag = 0;
    FILE* g_pFileLinkMOE = nullptr;

    string file_name = "TD_link_state.csv";

    fopen_ss(&g_pFileLinkMOE, file_name.c_str(), "w");

    if (!g_pFileLinkMOE)
    {
        dtalog.output() << "File " << file_name.c_str() << "cannot be opened." << endl;
        g_ProgramStop();
    }
    else
    {

        // Option 2: BPR-X function
        fprintf(g_pFileLinkMOE, "link_id,from_node_id,to_node_id,time_period,duration_in_sec,state,state_name\n");


        //Initialization for all nodes
        for (unsigned li = 0; li < g_link_vector.size(); ++li)
        {
            if (g_link_vector[li].timing_arc_flag)
            {
                // reset for signalized links (not freeway links as type code != 'f' for the case of freeway workzones)
                // only for the loading period

                int t = 0;
                int last_t = t;
                int current_state = assignment.m_LinkOutFlowState[li][t];
                while (t < assignment.g_number_of_loading_intervals_in_sec - 1)
                {
                    int next_state = assignment.m_LinkOutFlowState[li][t + 1];

                    if (next_state == current_state && t < assignment.g_number_of_loading_intervals_in_sec-2)
                    {
                        // do nothing 
                    }   
                    else
                    {  // change of state
                        string state_str;

                        if (current_state == 1)
                            state_str = "g";

                        if (current_state == 0)
                            state_str = "r";

                        fprintf(g_pFileLinkMOE, "%s,%d,%d,%s_%s,%d,%d,%s\n",
                            g_link_vector[li].link_id.c_str(),
                            g_node_vector[g_link_vector[li].from_node_seq_no].node_id,
                            g_node_vector[g_link_vector[li].to_node_seq_no].node_id,
                            g_time_coding(assignment.g_LoadingStartTimeInMin + last_t / 60.0).c_str(),
                            g_time_coding(assignment.g_LoadingStartTimeInMin + (t+1) / 60.0).c_str(),
                            t+1 - last_t,
                            current_state,
                            state_str.c_str());

                        last_t = t+1;
                        current_state = assignment.m_LinkOutFlowState[li][t+1];

                    }
                    t++;
                }


            }

        }
        fclose(g_pFileLinkMOE);

    }
}

void g_output_TD_link_capacity(Assignment& assignment, int output_mode = 1)
{
    dtalog.output() << "writing TD_link_capacity.csv.." << endl;

    int b_debug_detail_flag = 0;
    FILE* g_pFileLinkMOE = nullptr;

    string file_name = "TD_link_capacity.csv";

    fopen_ss(&g_pFileLinkMOE, file_name.c_str(), "w");

    if (!g_pFileLinkMOE)
    {
        dtalog.output() << "File " << file_name.c_str() << "cannot be opened." << endl;
        g_ProgramStop();
    }
    else
    {

        // Option 2: BPR-X function
        fprintf(g_pFileLinkMOE, "link_id,from_node_id,to_node_id,time_stamp,capacity,state,state_name\n");


        //Initialization for all nodes
        for (unsigned li = 0; li < g_link_vector.size(); ++li)
        {
            if (g_link_vector[li].flow_volume_per_period[0] > 1)
            {
                // reset for signalized links (not freeway links as type code != 'f' for the case of freeway workzones)
                // only for the loading period

                int t = 0;
                int last_t = t;
                while (t < assignment.g_number_of_loading_intervals_in_sec - 1)
                {
                    int current_state = -1;
                    
                    if(g_link_vector[li].timing_arc_flag == 1)
                        current_state =  assignment.m_LinkOutFlowState[li][t];

                    float current_cap = assignment.m_LinkOutFlowCapacity[li][t];

                        string state_str;

                        if (current_state == 1)
                            state_str = "g";

                        if (current_state == 0)
                            state_str = "r";

                        fprintf(g_pFileLinkMOE, "%s,%d,%d,T%s,%.3f,%d,%s\n",
                            g_link_vector[li].link_id.c_str(),
                            g_node_vector[g_link_vector[li].from_node_seq_no].node_id,
                            g_node_vector[g_link_vector[li].to_node_seq_no].node_id,
                            g_time_coding(assignment.g_LoadingStartTimeInMin + t / 60.0).c_str(),
                            current_cap,
                            current_state,
                            state_str.c_str());

                    t++;
                }


            }

        }
        fclose(g_pFileLinkMOE);

    }
}
void g_output_simulation_agents(Assignment& assignment)
{
    if (assignment.assignment_mode == 0 || assignment.trajectory_output == 0)  //LUE
    {
        FILE* g_pFilePathMOE = nullptr;
        fopen_ss(&g_pFilePathMOE, "agent.csv", "w");
        fclose(g_pFilePathMOE);
    }
    else if (assignment.assignment_mode >= 1)  //UE mode, or ODME, DTA
    {
        dtalog.output() << "writing agent.csv.." << endl;

        float path_time_vector[_MAX_LINK_SIZE_IN_A_PATH];
        FILE* g_pFileAgent = nullptr;
        fopen_ss(&g_pFileAgent, "agent.csv", "w");

        if (!g_pFileAgent)
        {
            dtalog.output() << "File agent.csv cannot be opened." << endl;
            g_ProgramStop();
        }

        fprintf(g_pFileAgent, "agent_id,o_zone_id,d_zone_id,OD_index,path_id,#_of_links,diversion_flag,agent_type,demand_period,volume,toll,travel_time,distance,speed,departure_time_in_min,arrival_time_in_min,departure_time_slot_no,\n");

        int count = 1;

        clock_t start_t, end_t;
        start_t = clock();
        clock_t iteration_t;

        int agent_type_size = assignment.g_AgentTypeVector.size();
        int zone_size = g_zone_vector.size();
        int demand_period_size = assignment.g_DemandPeriodVector.size();

        CColumnVector* p_column_pool;

        float path_toll = 0;
        float path_distance = 0;
        float path_travel_time = 0;
        float time_stamp = 0;

        if (assignment.trajectory_sampling_rate < 0.01)
            assignment.trajectory_sampling_rate = 0.01;
        int sampling_step = 100 / int(100 * assignment.trajectory_sampling_rate + 0.5);

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
                        p_column_pool = &(assignment.g_column_pool[orig][dest][at][tau]);
                        if (p_column_pool->od_volume > 0)
                        {

                            // scan through the map with different node sum for different continuous paths
                            it_begin = p_column_pool->path_node_sequence_map.begin();
                            it_end = p_column_pool->path_node_sequence_map.end();

                            for (it = it_begin; it != it_end; ++it)
                            {
                                if (count % 100000 == 0)
                                {
                                    end_t = clock();
                                    iteration_t = end_t - start_t;
                                    dtalog.output() << "writing " << count / 1000 << "K agents with CPU time " << iteration_t / 1000.0 << " s" << endl;
                                }

                                if (count % sampling_step != 0)
                                    continue;


                                path_toll = 0;
                                path_distance = 0;
                                path_travel_time = 0;
                                path_time_vector[0] = time_stamp;



                                // assignment_mode = 1, path flow mode
                                {
                                    // assignment_mode = 2, simulated agent flow mode //DTA simulation 

                                    for (int vi = 0; vi < it->second.agent_simu_id_vector.size(); ++vi)
                                    {
                                        int agent_simu_id = it->second.agent_simu_id_vector[vi];
                                        CAgent_Simu* pAgentSimu = g_agent_simu_vector[agent_simu_id];
                                        time_stamp = assignment.g_LoadingStartTimeInMin + pAgentSimu->departure_time_in_min;

                                        float departure_time_in_min = time_stamp ;
                                        float arrival_time_in_min = assignment.g_LoadingStartTimeInMin + pAgentSimu->arrival_time_in_min;
                                        int departure_time_in_slot_no = time_stamp / MIN_PER_TIMESLOT;
                                        float speed = pAgentSimu->path_distance / max(0.001,pAgentSimu->path_travel_time_in_min) * 60;

                                        fprintf(g_pFileAgent, "%d,%d,%d,%d->%d,%d,%d,%d,%s,%s,1,%.1f,%.4f,%.4f,%.4f,%.4f,%.4f,%d",
                                            pAgentSimu->agent_id,
                                            g_zone_vector[orig].zone_id,
                                            g_zone_vector[dest].zone_id,
                                            g_zone_vector[orig].zone_id,
                                            g_zone_vector[dest].zone_id,
                                            it->second.path_seq_no,
                                            it->second.m_link_size,
                                            pAgentSimu->diversion_flag,
                                            assignment.g_AgentTypeVector[at].agent_type.c_str(),
                                            assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
                                            path_toll,
                                            pAgentSimu->path_travel_time_in_min,
                                            pAgentSimu->path_distance, speed, 
                                            departure_time_in_min,
                                            arrival_time_in_min,
                                            departure_time_in_slot_no);

                                        count++;

                                        fprintf(g_pFileAgent, "\n");
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        fclose(g_pFileAgent);
    }
}

void g_output_simulation_result(Assignment& assignment) 
{
    g_output_TD_link_performance(assignment, 1);
    g_output_TD_link_performance(assignment, 2);
    if (assignment.assignment_mode == 2)  //DTA mode
    {
    g_output_TD_link_state(assignment, 1);
     }
//    g_output_TD_link_capacity(assignment, 1);

    g_output_simulation_agents(assignment);

    if (assignment.assignment_mode == 0 || assignment.trajectory_output==0)  //LUE
    {
        FILE* g_pFilePathMOE = nullptr;
        fopen_ss(&g_pFilePathMOE, "trajectory.csv", "w");
        fclose(g_pFilePathMOE);
    }
    else if (assignment.assignment_mode >= 1)  //UE mode, or ODME, DTA
    {
        dtalog.output() << "writing trajectory.csv.." << endl;

        float path_time_vector[_MAX_LINK_SIZE_IN_A_PATH];
        FILE* g_pFilePathMOE = nullptr;
        fopen_ss(&g_pFilePathMOE, "trajectory.csv", "w");

        if (!g_pFilePathMOE)
        {
            dtalog.output() << "File trajectory.csv cannot be opened." << endl;
            g_ProgramStop();
        }

        fprintf(g_pFilePathMOE, "agent_id,o_zone_id,d_zone_id,path_id,diversion_flag,agent_type,PCE_unit,demand_period,volume,toll,travel_time,distance,node_sequence,link_sequence,time_sequence,waiting_time_in_simu_interval,geometry\n");

        int count = 1;

        clock_t start_t, end_t;
        start_t = clock();
        clock_t iteration_t;

        int agent_type_size = assignment.g_AgentTypeVector.size();
        int zone_size = g_zone_vector.size();
        int demand_period_size = assignment.g_DemandPeriodVector.size();

        CColumnVector* p_column_pool;

        float path_toll = 0;
        float path_distance = 0;
        float path_travel_time = 0;
        float time_stamp = 0;

        if (assignment.trajectory_sampling_rate < 0.01)
            assignment.trajectory_sampling_rate = 0.01;
        int sampling_step = 100/int(100*assignment.trajectory_sampling_rate+0.5);

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
                        p_column_pool = &(assignment.g_column_pool[orig][dest][at][tau]);
                        if (p_column_pool->od_volume > 0)
                        {

                            // scan through the map with different node sum for different continuous paths
                            it_begin = p_column_pool->path_node_sequence_map.begin();
                            it_end = p_column_pool->path_node_sequence_map.end();

                            for (it = it_begin; it != it_end; ++it)
                            {
                                if (count % 100000 == 0)
                                {
                                    end_t = clock();
                                    iteration_t = end_t - start_t;
                                    dtalog.output() << "writing " << count / 1000 << "K agents with CPU time " << iteration_t / 1000.0 << " s" << endl;
                                }

                                if (count % sampling_step != 0)
                                    continue;


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

                                    time_stamp += link_travel_time;
                                    path_time_vector[nl + 1] = time_stamp;
                                }

                                int virtual_link_delta = 1;
                                int virtual_begin_link_delta = 1;
                                int virtual_end_link_delta = 1;
                                // fixed routes have physical nodes always, without virtual connectors
                                if (p_column_pool->bfixed_route)
                                {
                                    virtual_begin_link_delta = 0;
                                    virtual_end_link_delta = 1;

                                }

                                // assignment_mode = 1, path flow mode
                                {
                                    // assignment_mode = 2, simulated agent flow mode //DTA simulation 

                                    for (int vi = 0; vi < it->second.agent_simu_id_vector.size(); ++vi)
                                    {


                                        int agent_simu_id = it->second.agent_simu_id_vector[vi];
                                        CAgent_Simu* pAgentSimu = g_agent_simu_vector[agent_simu_id];




                                        if (pAgentSimu->agent_id == 81)
                                        {
                                            int idebug = 1;
                                        }
                                        if (assignment.trajectory_diversion_only == 1 && pAgentSimu->diversion_flag == 0)  // diversion flag only, then we skip the non-diversion path
                                            continue;

                                        time_stamp = assignment.g_LoadingStartTimeInMin + pAgentSimu->departure_time_in_min;
                                        for (int nt = 0 + virtual_link_delta; nt < pAgentSimu->path_link_seq_no_vector.size() + 1 - virtual_link_delta; ++nt)
                                        {
                                            float time_in_min = 0;

                                            if (nt < pAgentSimu->path_link_seq_no_vector.size() - virtual_link_delta)
                                                time_in_min = assignment.g_LoadingStartTimeInMin + pAgentSimu->m_Veh_LinkArrivalTime_in_simu_interval[nt] * number_of_seconds_per_interval / 60.0;
                                            else
                                                time_in_min = assignment.g_LoadingStartTimeInMin + pAgentSimu->m_Veh_LinkDepartureTime_in_simu_interval[nt - 1] * number_of_seconds_per_interval / 60.0;  // last link in the path

                                            path_time_vector[nt] = time_in_min;
                                        }

                                        float vehicle_travel_time = pAgentSimu->path_travel_time_in_min;

                                        
                                        // some bugs for output link performances before
                                        fprintf(g_pFilePathMOE, "%d,%d,%d,%d,%d,%s,%d,%s,1,%.1f,%.4f,%.4f,",
                                            pAgentSimu->agent_id,
                                            g_zone_vector[orig].zone_id,
                                            g_zone_vector[dest].zone_id,
                                            it->second.path_seq_no,
                                            pAgentSimu->diversion_flag,
                                            assignment.g_AgentTypeVector[at].agent_type.c_str(),
                                            pAgentSimu->PCE_unit_size,
                                            assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
                                            path_toll,
                                            vehicle_travel_time,
                                            path_distance);

                                        /* Format and print various data */

                                        for (int ni = 0 + virtual_link_delta; ni < pAgentSimu->path_link_seq_no_vector.size(); ++ni)
                                        {
                                            int node_id = g_node_vector[g_link_vector[pAgentSimu->path_link_seq_no_vector[ni]].from_node_seq_no].node_id;
                                            fprintf(g_pFilePathMOE, "%d;", node_id);

                                        }
                              
                                        fprintf(g_pFilePathMOE, ",");

                                        //for (int nl = 0 + virtual_link_delta; nl < pAgentSimu->path_link_seq_no_vector.size() - virtual_link_delta; ++nl)
                                        //{
                                        //    int link_seq_no = it->second.path_link_vector[nl];
                                        //    fprintf(g_pFilePathMOE, "%s;", g_link_vector[link_seq_no].link_id.c_str());
                                        //}
                                        fprintf(g_pFilePathMOE, ",");




                                        for (int nt = 0 + virtual_link_delta; nt < pAgentSimu->path_link_seq_no_vector.size() + 1 - virtual_link_delta; ++nt)
                                            fprintf(g_pFilePathMOE, "%s;", g_time_coding(path_time_vector[nt]).c_str());
                                        
                                        fprintf(g_pFilePathMOE, ",");

                                        //// waiting time in simu interval
                                        //int waiting_time_in_simu_interaval = 0;
                                        //for (int nt = 0 + virtual_link_delta; nt < pAgentSimu->path_link_seq_no_vector.size() - virtual_link_delta; ++nt)
                                        //{
                                        //    int link_seq_no = it->second.path_link_vector[nt];

                                        //    waiting_time_in_simu_interaval = (path_time_vector[nt + 1] - path_time_vector[nt] - g_link_vector[link_seq_no].free_flow_travel_time_in_min) * number_of_simu_intervals_in_min;
                                        //    fprintf(g_pFilePathMOE, "%d;", waiting_time_in_simu_interaval);

                                        //        
                                        //}

                                        fprintf(g_pFilePathMOE, ",");


                                        fprintf(g_pFilePathMOE, "\"LINESTRING (");


                                        for (int ni = 0 + virtual_link_delta; ni < pAgentSimu->path_link_seq_no_vector.size(); ++ni)
                                        {
                                        
                                           int node_no = g_link_vector[pAgentSimu->path_link_seq_no_vector[ni]].from_node_seq_no;
                                            fprintf(g_pFilePathMOE, "%f %f", g_node_vector[node_no].x,
                                                g_node_vector[node_no].y);

                                            if (ni != pAgentSimu->path_link_seq_no_vector.size() - 1)
                                                fprintf(g_pFilePathMOE, ", ");
                                        }

                                       
                                        

                                        fprintf(g_pFilePathMOE, ")\"\n");

                                        count++;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        fclose(g_pFilePathMOE);
    }

    int b_trace_file = false;

    // output trace file
    if (assignment.assignment_mode == 0 || assignment.trajectory_output == 0)  //LUE
    {
        FILE* g_pFilePathMOE = nullptr;
        fopen_ss(&g_pFilePathMOE, "trace.csv", "w");
        fclose(g_pFilePathMOE);
    }
    else if (assignment.assignment_mode >= 1)  //UE mode, or ODME, DTA
    {
        dtalog.output() << "writing trace.csv.." << endl;

        float path_time_vector[_MAX_LINK_SIZE_IN_A_PATH];
        FILE* g_pFilePathMOE = nullptr;
        fopen_ss(&g_pFilePathMOE, "trace.csv", "w");

        if (!g_pFilePathMOE)
        {
            dtalog.output() << "File trace.csv cannot be opened." << endl;
            g_ProgramStop();
        }

        fprintf(g_pFilePathMOE, "agent_id,seq_no,node_id,timestamp,timestamp_in_min,trip_time_in_min,travel_time_in_sec,waiting_time_in_simu_interval,x_coord,y_coord\n");

        int count = 1;

        clock_t start_t, end_t;
        start_t = clock();
        clock_t iteration_t;

        int agent_type_size = assignment.g_AgentTypeVector.size();
        int zone_size = g_zone_vector.size();
        int demand_period_size = assignment.g_DemandPeriodVector.size();

        CColumnVector* p_column_pool;

        float path_toll = 0;
        float path_distance = 0;
        float path_travel_time = 0;
        float time_stamp = 0;

        if (assignment.trajectory_sampling_rate < 0.01)
            assignment.trajectory_sampling_rate = 0.01;
        int sampling_step = 1 ;

        std::map<int, CColumnPath>::iterator it, it_begin, it_end;

        dtalog.output() << "writing data for " << zone_size << "  zones " << endl;

        for (int orig = 0; orig < zone_size; ++orig)
        {
            for (int at = 0; at < agent_type_size; ++at)
            {
                for (int dest = 0; dest < zone_size; ++dest)
                {
                    for (int tau = 0; tau < demand_period_size; ++tau)
                    {
                        p_column_pool = &(assignment.g_column_pool[orig][dest][at][tau]);
                        if (p_column_pool->od_volume > 0)
                        {

                            // scan through the map with different node sum for different continuous paths
                            it_begin = p_column_pool->path_node_sequence_map.begin();
                            it_end = p_column_pool->path_node_sequence_map.end();

                            for (it = it_begin; it != it_end; ++it)
                            {
                                if (count % 100000 == 0)
                                {
                                    end_t = clock();
                                    iteration_t = end_t - start_t;
                                    dtalog.output() << "writing " << count / 1000 << "K agents with CPU time " << iteration_t / 1000.0 << " s" << endl;
                                }

                                if (count % sampling_step != 0)
                                    continue;

                                if (count >= 1000)
                                    break;

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

                                // Xuesong: 11/20/2021, need to check again 
                                // fixed routes have physical nodes always, without virtual connectors
                                //if (p_column_pool->bfixed_route)
                                //    virtual_link_delta = 0;

                                // assignment_mode = 1, path flow mode
                                {
                                    // assignment_mode = 2, simulated agent flow mode //DTA simulation 

                                    for (int vi = 0; vi < it->second.agent_simu_id_vector.size(); ++vi)
                                    {

                                        int agent_simu_id = it->second.agent_simu_id_vector[vi];
                                        CAgent_Simu* pAgentSimu = g_agent_simu_vector[agent_simu_id];


                                        if (pAgentSimu->agent_id == 81)
                                        {
                                            int idebug = 1;
                                        }
                                        //if (assignment.trajectory_diversion_only == 1 && pAgentSimu->diversion_flag == 0)  // diversion flag only, then we skip the non-diversion path
                                        //     continue;

                                        for (int nt = 0 + virtual_link_delta; nt < pAgentSimu->path_link_seq_no_vector.size() + 1 - virtual_link_delta; ++nt)
                                        {
                                            float time_in_min = 0;

                                            if (nt < pAgentSimu->path_link_seq_no_vector.size() - virtual_link_delta)
                                                time_in_min = assignment.g_LoadingStartTimeInMin + pAgentSimu->m_Veh_LinkArrivalTime_in_simu_interval[nt] * number_of_seconds_per_interval / 60.0;
                                            else
                                                time_in_min = assignment.g_LoadingStartTimeInMin + pAgentSimu->m_Veh_LinkDepartureTime_in_simu_interval[nt - 1] * number_of_seconds_per_interval / 60.0;  // last link in the path

                                            path_time_vector[nt] = time_in_min;
                                        }


                                        for (int nt = 0 + virtual_link_delta; nt < pAgentSimu->path_link_seq_no_vector.size() - virtual_link_delta; ++nt)
                                        {
                                            float trip_time_in_min = path_time_vector[nt] - path_time_vector[0 + virtual_link_delta];

                                            int node_id = g_node_vector[g_link_vector[pAgentSimu->path_link_seq_no_vector[nt]].from_node_seq_no].node_id;

                                            if (nt >= pAgentSimu->path_link_seq_no_vector.size() - virtual_link_delta)
                                            {
                                                node_id = g_node_vector[g_link_vector[pAgentSimu->path_link_seq_no_vector[nt]].to_node_seq_no].node_id;

                                            }
                                            int link_seq_no = it->second.path_link_vector[nt];
                                            float travel_time_in_sec = (path_time_vector[nt + 1] - path_time_vector[nt])*60;
                                            int waiting_time_in_simu_interaval = (path_time_vector[nt + 1] - path_time_vector[nt] - g_link_vector[link_seq_no].free_flow_travel_time_in_min) * number_of_simu_intervals_in_min;

                                            int node_no = g_link_vector[pAgentSimu->path_link_seq_no_vector[nt]].from_node_seq_no;

                                            fprintf(g_pFilePathMOE, "%d,%d,%d,T%s,%.5f,%f,%.4f,%d,%f,%f\n",
                                                pAgentSimu->agent_id, nt, node_id, 
                                                g_time_coding(path_time_vector[nt]).c_str(),
                                                path_time_vector[nt],
                                                trip_time_in_min,
                                                travel_time_in_sec,
                                                waiting_time_in_simu_interaval,
                                                g_node_vector[node_no].x, g_node_vector[node_no].y);

                                        }
                                        count++;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        fclose(g_pFilePathMOE);
    }

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

    NetworkForSP* PointerMatrx[_MAX_AGNETTYPES][_MAX_TIMEPERIODS][_MAX_MEMORY_BLOCKS] = {NULL};
    NetworkForSP* RTPointerMatrx[_MAX_AGNETTYPES][_MAX_TIMEPERIODS][_MAX_MEMORY_BLOCKS] = { NULL };

    int computing_zone_count = 0;

    for (int at = 0; at < assignment.g_AgentTypeVector.size(); ++at)
    {
        for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); ++tau)
        {
            //assign all nodes to the corresponding thread
            for (int z = 0; z < g_zone_vector.size(); ++z)
            {
                if (g_zone_vector[z].b_real_time_information == true)  // do not consider zones creatd for real time information  here
                    continue;

                 if (z < assignment.g_number_of_memory_blocks)
                {
                    NetworkForSP* p_NetworkForSP = new NetworkForSP();

                        p_NetworkForSP->m_origin_node_vector.push_back(g_zone_vector[z].node_seq_no);
                        p_NetworkForSP->m_origin_zone_seq_no_vector.push_back(z);

                        p_NetworkForSP->m_agent_type_no = at;
                        p_NetworkForSP->m_tau = tau;
                        computing_zone_count++;
                    
                    p_NetworkForSP->AllocateMemory(assignment.g_number_of_nodes, assignment.g_number_of_links);

                    PointerMatrx[at][tau][z] = p_NetworkForSP;

                    g_NetworkForSP_vector.push_back(p_NetworkForSP);
                }
                else  // zone seq no is greater than g_number_of_memory_blocks
                {
                     if (assignment.g_origin_demand_array[z] > 0.001)
                     {
                         //get the corresponding memory block seq no
                        // take residual of memory block size to map a zone no to a memory block no.
                         int memory_block_no = z % assignment.g_number_of_memory_blocks;
                         NetworkForSP* p_NetworkForSP = PointerMatrx[at][tau][memory_block_no];
                         p_NetworkForSP->m_origin_node_vector.push_back(g_zone_vector[z].node_seq_no);
                         p_NetworkForSP->m_origin_zone_seq_no_vector.push_back(z);
                         computing_zone_count++;
                     }
                }
            }
        }
        
    }

    int info_computing_zone_no = 0;
    for (int at = 0; at < assignment.g_AgentTypeVector.size(); ++at)
    {
        for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); ++tau)
        {
            //assign all nodes to the corresponding thread
            for (int z = 0; z < g_zone_vector.size(); ++z)
            {
                if (g_zone_vector[z].b_real_time_information == false)
                    continue; 
                // consider zones created for real time information 

                    int memory_block_no = info_computing_zone_no % assignment.g_number_of_memory_blocks;
                    NetworkForSP* p_NetworkForSP = RTPointerMatrx[at][tau][memory_block_no];
                    if (p_NetworkForSP == NULL)
                    {  // exception handling as the information zone id does not start from 0
                        int ibebug = 1;
                        p_NetworkForSP = new NetworkForSP();

                        p_NetworkForSP->m_origin_node_vector.push_back(g_zone_vector[z].node_seq_no);
                        p_NetworkForSP->m_origin_zone_seq_no_vector.push_back(z);

                        p_NetworkForSP->m_agent_type_no = at;
                        p_NetworkForSP->m_tau = tau;
                        computing_zone_count++;

                        p_NetworkForSP->AllocateMemory(assignment.g_number_of_nodes, assignment.g_number_of_links);

                        RTPointerMatrx[at][tau][memory_block_no] = p_NetworkForSP;

                        g_NetworkForRTSP_vector.push_back(p_NetworkForSP);
                        info_computing_zone_no++;
                    }
                    else
                    {
                    p_NetworkForSP->m_origin_node_vector.push_back(g_zone_vector[z].node_seq_no);
                    p_NetworkForSP->m_origin_zone_seq_no_vector.push_back(z);
                    info_computing_zone_no++;
                    }
            }
        }

    }
    dtalog.output() << "There are " << g_NetworkForSP_vector.size() << " SP networks in memory." << endl;
    dtalog.output() << "There are " << g_NetworkForRTSP_vector.size() << " RTSP networks in memory." << endl;
    dtalog.output() << "There are " << computing_zone_count << " zones to be computed in CPU." << endl;

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

            if (ODvolume > 0.000001 || assignment.zone_seq_no_2_info_mapping.find(m_origin_zone_seq_no)!= assignment.zone_seq_no_2_info_mapping.end())
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
    RT_travel_time = 0; // reset RT_travel time for each end of simulation iteration 
    for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
    {
        float starting_time_slot_no = assignment.g_DemandPeriodVector[tau].starting_time_slot_no;
        float ending_time_slot_no = assignment.g_DemandPeriodVector[tau].ending_time_slot_no;

        travel_time_per_period[tau] = VDF_period[tau].PerformBPR(flow_volume_per_period[tau]);

        // signalized with red_time > 1
        if (this->timing_arc_flag && this->mvmt_txt_id.length() > 1 && VDF_period[tau].red_time > 1)
        {
            // arterial streets with the data from sigal API
            float hourly_per_lane_volume = flow_volume_per_period[tau] / (max(1.0f, (ending_time_slot_no - starting_time_slot_no)) * 15 / 60 / number_of_lanes);
            float red_time = VDF_period[tau].red_time;
            float cycle_length = VDF_period[tau].cycle_length;
            //we dynamically update cycle length, and green time/red time, so we have dynamically allocated capacity and average delay
            travel_time_per_period[tau] = VDF_period[tau].PerformSignalVDF(hourly_per_lane_volume, red_time, cycle_length);

            // update capacity using the effective discharge rates, will be passed in to the following BPR function
            VDF_period[tau].capacity = (1 - red_time / cycle_length) * _default_saturation_flow_rate * number_of_lanes;
        }

        // either non-signalized or signalized with red_time < 1 and cycle_length < 30
        ////if (this->mvmt_txt_id.length() == 0 || (this->mvmt_txt_id.length() > 1 && VDF_period[tau].red_time < 1 && VDF_period[tau].cycle_length < 30))
        ////{
        ////    travel_time_per_period[tau] = VDF_period[tau].PerformBPR(flow_volume_per_period[tau]);
        ////    //VDF_period[tau].PerformBPR_X(flow_volume_per_period[tau]);  // only for freeway segments
        ////}
    }
}

void  CLink::CalculateTD_RTVDFunction()
{
    
//
//    for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
//    {
//  //       signalized with red_time > 1
//        if (this->movement_str.length() > 1)
//        {
////             arterial streets with the data from sigal API
//            float hourly_per_lane_volume = RT_flow_volume / (max(1.0f, (ending_time_slot_no - starting_time_slot_no)) * 15 / 60 / number_of_lanes);
//            float red_time = VDF_period[tau].red_time;
//            float cycle_length = VDF_period[tau].cycle_length;
//    //        we dynamically update cycle length, and green time/red time, so we have dynamically allocated capacity and average delay
//            travel_time_per_period[tau] = VDF_period[tau].PerformSignalVDF(hourly_per_lane_volume, red_time, cycle_length);
//            travel_time_per_period[tau] += VDF_period[tau].PerformBPR(RT_flow_volume);
//      //       update capacity using the effective discharge rates, will be passed in to the following BPR function
//            VDF_period[tau].capacity = (1 - red_time / cycle_length) * _default_saturation_flow_rate * number_of_lanes;
//        }
//
////         either non-signalized or signalized with red_time < 1 and cycle_length < 30
//        if (this->movement_str.length() == 0 || (this->movement_str.length() > 1 && VDF_period[tau].red_time < 1 && VDF_period[tau].cycle_length < 30))
//        {
//            travel_time_per_period[tau] = VDF_period[tau].PerformBPR(RT_flow_volume);
//        }
//    }
}
double network_assignment(int assignment_mode, int iteration_number, int column_updating_iterations, int ODME_iterations, int number_of_memory_blocks)
{
    int signal_updating_iterations = 0;

    // k iterations for column generation
    assignment.g_number_of_column_generation_iterations = iteration_number;
    // 0: link UE: 1: path UE, 2: Path SO, 3: path resource constraints
    assignment.assignment_mode = assignment_mode;

	assignment.g_number_of_memory_blocks = min( max(1, number_of_memory_blocks), _MAX_MEMORY_BLOCKS);
	

    if (assignment.assignment_mode == 0)
        column_updating_iterations = 0;

    // step 1: read input data of network / demand tables / Toll
    g_ReadInputData(assignment);
    //g_reload_timing_arc_data(assignment); // no need to load timing data from timing.csv
    g_load_scenario_data(assignment);
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
			for (int blk = 0; blk < assignment.g_AgentTypeVector.size()*assignment.g_DemandPeriodVector.size(); ++blk)
			{
                int network_copy_no = blk* assignment.g_number_of_memory_blocks + ProcessID;
                if (network_copy_no >= g_NetworkForSP_vector.size())
                    continue;

				NetworkForSP* pNetwork = g_NetworkForSP_vector[network_copy_no];

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

    if (assignment.assignment_mode == 2)
    {
        dtalog.output() << "Step 5: Simulation for traffic assignment.." << endl;
        assignment.STTrafficSimulation();
        dtalog.output() << endl;
    }

    if (assignment.assignment_mode == 3)
    {
        dtalog.output() << "Step 6: O-D estimation for traffic assignment.." << endl;
        assignment.Demand_ODME(ODME_iterations);
        dtalog.output() << endl;
    }

    end_t = clock();
    total_t = (end_t - start_t);
    dtalog.output() << "Done!" << endl;

    dtalog.output() << "CPU Running Time for column pool updating: " << total_t / 1000.0 << " s" << endl;

    start_t = clock();

    //step 5: output simulation results of the new demand
    g_ReadOutputFileConfiguration(assignment);
    g_output_assignment_result(assignment);
    g_output_simulation_result(assignment);

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
bool Assignment::RTSP_RealTimeShortestPathFinding(int time_slot_no, int simu_time_interval_t)
{
clock_t start_t_lc = clock();
clock_t start_t_cp = clock();

bool bRealTimeInformationActicvated = false;
for (int i = 0; i < g_link_vector.size(); ++i)
{
    CLink* pLink = &(g_link_vector[i]);
    int slot_no = time_slot_no;
    float RT_capacity = pLink->TD_link_capacity[slot_no];

    if (time_slot_no >= pLink->TD_link_reduction_start_time_slot_no)  // later than the earliest start time slot
    {
        bRealTimeInformationActicvated = true;
    }
}
if (bRealTimeInformationActicvated == false)  // as long as there are incident activated
return false;

for (int i = 0; i < g_link_vector.size(); ++i)
{
    CLink* pLink = &(g_link_vector[i]);
    int slot_no = time_slot_no;
    float RT_capacity = pLink->TD_link_capacity[slot_no];
    if (pLink->link_type >= 0)  // do not need to be the cap reduced link
           {

                    double total_waiting_time_in_min = 0;

                    for (auto it = pLink->ExitQueue.begin(); it != pLink->ExitQueue.end(); ++it)
                    {
                        int agent_id = (*it);
                        CAgent_Simu* p_agent = g_agent_simu_vector[agent_id];


                        int current_link_seq_no = p_agent->path_link_seq_no_vector[p_agent->m_current_link_seq_no];
                        int arrival_time_in_t = p_agent->m_Veh_LinkArrivalTime_in_simu_interval[p_agent->m_current_link_seq_no];
                        int waiting_time_in_t = simu_time_interval_t - arrival_time_in_t;
                        total_waiting_time_in_min += waiting_time_in_t * number_of_seconds_per_interval / 60; // in min

                    }
                    pLink->RT_travel_time = total_waiting_time_in_min / max(1, pLink->ExitQueue.size());  // average travel time

                    if (pLink->TD_link_closure_map.find(slot_no) != pLink->TD_link_closure_map.end())
                    {
                        pLink->RT_travel_time = 9999;
                    }

                    pLink->RT_speed_vector[time_slot_no] = pLink->length / max(0.00001, pLink->RT_travel_time / 60.0);  // mph                }



                    pLink->RT_travel_time_vector[time_slot_no] = pLink->RT_travel_time;

          }
        

}

 

#pragma omp parallel for  // step 3: C++ open mp automatically create n threads., each thread has its own computing thread on a cpu core
        for (int blk = 0; blk < g_NetworkForRTSP_vector.size(); ++blk)
        {

            NetworkForSP* pNetwork = g_NetworkForRTSP_vector[blk];
            int iteration_number = 0;

            for (int o_node_index = 0; o_node_index < pNetwork->m_origin_node_vector.size(); ++o_node_index)
            {
                start_t_lc = clock();
                pNetwork->optimal_label_correcting(blk, &assignment, iteration_number, o_node_index);


                start_t_cp = clock();
                double total_origin_least_travel_time = pNetwork->backtrace_shortest_path_tree(assignment, iteration_number, o_node_index);


            }
        

        }

        return true;

}

void Assignment::UpdateRTPath(CAgent_Simu* pAgent)
{
    // updating shorest path for vehicles passing through information node

    std::map<int, CColumnPath>::iterator it, it_begin, it_end;

    CColumnVector* p_column_pool;

        int at = pAgent->at;
        int dest = pAgent->dest;
        int tau = pAgent->tau;

        int current_link_seq_no = pAgent->path_link_seq_no_vector[pAgent->m_current_link_seq_no];
        CLink* pCurrentLink = &(g_link_vector[current_link_seq_no]);
        
        if (node_seq_no_2_info_zone_seq_no_mapping.find(pCurrentLink->to_node_seq_no) == node_seq_no_2_info_zone_seq_no_mapping.end())
        {
            return;
        }

        int orig = node_seq_no_2_info_zone_seq_no_mapping[pCurrentLink->to_node_seq_no];

        p_column_pool = &(assignment.g_column_pool[orig][dest][at][tau]);
//             scan through the map with different node sum for different continuous paths
            it_begin = p_column_pool->path_node_sequence_map.begin();
            it_end = p_column_pool->path_node_sequence_map.end();

            

                for (it = it_begin; it != it_end; ++it)
                {
                    bool b_original_path_impacted_flag = false;
                    int link_sum_0 = 0;
                    int link_sum_rerouting = 0;
                    for (int ls = pAgent->m_current_link_seq_no+1; ls < pAgent->path_link_seq_no_vector.size(); ls++)
                    {
                        link_sum_0 += pAgent->path_link_seq_no_vector[ls];

                        if (g_link_vector[pAgent->path_link_seq_no_vector[ls]].capacity_reduction_flag == 1)
                            b_original_path_impacted_flag = true;
                    }

                    for (int nl = 0; nl < it->second.m_link_size-1; ++nl)  // arc a // exclude virtual link at the end;
                    {
                        link_sum_rerouting += it->second.path_link_vector[nl];
                    }

                    // detection
                    if (b_original_path_impacted_flag == false)  // this path is not impacted
                        continue; 

                    if (link_sum_0 != link_sum_rerouting)  
                    {
                        //dtalog.output() << "route is changed for agent . " << pAgent->agent_id << endl;

                        int debug_flag = 0;

                        if(debug_flag==1)
                            {
                            for (int ls = pAgent->m_current_link_seq_no + 1; ls < pAgent->path_link_seq_no_vector.size(); ls++)
                            {
                                link_sum_0 += pAgent->path_link_seq_no_vector[ls];
                                //dtalog.output() << "current no." << ls << ":" << pAgent->path_link_seq_no_vector[ls] << ";" << endl;
                            }

                            for (int nl = 0; nl < it->second.m_link_size - 1; ++nl)  // arc a // exclude virtual link at the end;
                            {
                                link_sum_rerouting += it->second.path_link_vector[nl];
                                dtalog.output() << "rerouting no." << nl << ":" << it->second.path_link_vector[nl] << ";" << endl;
                            }
                            }

                        p_column_pool->od_volume += 1;// this is used flag

                    if(it->second.m_link_size>=2 && it->second.diverted_vehicle_map.find(pAgent->agent_id) == it->second.diverted_vehicle_map.end())  // feasible rerouting and have not been informed by this sensor yet
                    {

                        if (pAgent->diversion_flag >= 1)   // a vehicle is diverted once
                            continue; 

                        if (pAgent->agent_id == 5802)
                        {
                            int idebug = 1;
                        }
                        pAgent->path_link_seq_no_vector.resize(pAgent->m_current_link_seq_no+1);  // consider virtual link
                        pAgent->m_Veh_LinkArrivalTime_in_simu_interval.resize(pAgent->m_current_link_seq_no+1);
                        pAgent->m_Veh_LinkDepartureTime_in_simu_interval.resize(pAgent->m_current_link_seq_no+1);

                        pAgent->diversion_flag += 1;
                        //expanding
                        for (int nl = 0; nl < it->second.m_link_size; ++nl)  // arc a  // we do not exclude virtual link at the end here. as the output will exclude the virtual link in trajectory.csv
                        {
                            int link_seq_no = it->second.path_link_vector[nl];
                            pAgent->path_link_seq_no_vector.push_back(link_seq_no);
                            pAgent->m_Veh_LinkArrivalTime_in_simu_interval.push_back(-1);
                            pAgent->m_Veh_LinkDepartureTime_in_simu_interval.push_back(-1);

                        }
                        it->second.path_volume += 1;

                        it->second.diverted_vehicle_map[pAgent->agent_id] = true;

                        break;
                    }

                    }
                }

  
}

void Assignment::AllocateLinkMemory4Simulation()
{
    g_number_of_simulation_intervals = (g_LoadingEndTimeInMin - g_LoadingStartTimeInMin + 60) * 60 /number_of_seconds_per_interval+2;

    g_number_of_intervals_in_sec = (g_LoadingEndTimeInMin - g_LoadingStartTimeInMin + 60) * 60;

    dtalog.output() << "LoadingStartTimeInMin = " << g_LoadingStartTimeInMin << endl;
    dtalog.output() << "g_LoadingStartTimeInMin = " << g_LoadingEndTimeInMin << endl;
    dtalog.output() << "number_of_simulation_intervals = " << g_number_of_simulation_intervals << endl;
    dtalog.output() << "number_of_simu intervals in sec = " << g_number_of_intervals_in_sec << endl;

    g_number_of_loading_intervals_in_sec = (g_LoadingEndTimeInMin - g_LoadingStartTimeInMin) * 60 ;

    g_number_of_intervals_in_min = (int)(g_number_of_simulation_intervals / number_of_simu_intervals_in_min +1);
    // add + 120 as a buffer
    g_number_of_in_memory_simulation_intervals = g_number_of_simulation_intervals;

    dtalog.output() << "allocate 2D dynamic memory LinkOutFlowCapacity..." << endl;

    m_LinkOutFlowCapacity = Allocate2DDynamicArray <float>(g_number_of_links, g_number_of_intervals_in_sec);  //1
    m_LinkOutFlowState = Allocate2DDynamicArray <int>(g_number_of_links, g_number_of_intervals_in_sec);  //1
    // discharge rate per simulation time interval
    dtalog.output() << "allocate 2D dynamic memory m_LinkCumulativeArrivalVector..." << endl;
    m_LinkCumulativeArrivalVector = Allocate2DDynamicArray <float>(g_number_of_links, g_number_of_intervals_in_min);  //2

    dtalog.output() << "allocate 2D dynamic memory m_LinkCumulativeDepartureVector..." << endl;
    m_LinkCumulativeDepartureVector = Allocate2DDynamicArray <float>(g_number_of_links, g_number_of_intervals_in_min);  //3

    m_LinkCACount = Allocate1DDynamicArray <float>(g_number_of_links);
    m_LinkCDCount = Allocate1DDynamicArray <float>(g_number_of_links);

    dtalog.output() << "allocate 2D dynamic memory m_LinkTDWaitingTime..." << endl;
    m_LinkTDWaitingTime = Allocate2DDynamicArray <float>(g_number_of_links, g_number_of_intervals_in_min); //5

    m_LinkTotalWaitingTimeVector.clear();
    for (int i = 0; i < g_number_of_links; ++i)
    {
        m_LinkTotalWaitingTimeVector.push_back(0.0);
    }
    
    dtalog.output() << "initializing time dependent capacity data..." << endl;

    unsigned int RandomSeed = 101;
    float residual;
    float random_ratio = 0;

#pragma omp parallel for
    for (int i = 0; i < g_number_of_links; ++i)
    {
        float cap_count = 0;
        float discharge_rate_per_sec = g_link_vector[i].lane_capacity * g_link_vector[i].number_of_lanes / 3600.0 ;
        float discharge_rate_after_loading = 10 * discharge_rate_per_sec;  ///to collect travel time statistics * 10 times of capacity to discharge all flow */ at the end of simulation time interval

        if (i == 24)
        {
            int debug = 1;
        }
        for (int t = 0; t < g_number_of_intervals_in_sec; ++t)
        {
            float OutFlowRate = 0;
            if (t >= g_number_of_loading_intervals_in_sec)
                OutFlowRate = discharge_rate_after_loading;
                /* 10 times of capacity to discharge all flow */
            else
                OutFlowRate = discharge_rate_per_sec;

            // free flow state
            //cell based simulation 

            m_LinkOutFlowCapacity[i][t] = discharge_rate_after_loading;

           if (g_link_vector[i].cell_type >= 0)  // DTA micro cell based simulation
               m_LinkOutFlowCapacity[i][t] = 1 * g_link_vector[i].number_of_lanes;
           else
           {
               residual = OutFlowRate - (int)(OutFlowRate);
               //RandomSeed is automatically updated.
               RandomSeed = (LCG_a * RandomSeed + LCG_c) % LCG_M;
               random_ratio = float(RandomSeed) / LCG_M;

               if (random_ratio < residual)
                   m_LinkOutFlowCapacity[i][t] = (int)(OutFlowRate) +1;
               else
                   m_LinkOutFlowCapacity[i][t] = (int)(OutFlowRate);

               if (t < g_number_of_loading_intervals_in_sec)
                   cap_count += m_LinkOutFlowCapacity[i][t];

           }

            // 1800 rule
            //if(t%2==0)
            //    m_LinkOutFlowCapacity[i][t] = 1*g_link_vector[i].number_of_lanes;
            //if (t % 2 == 1)
            //    m_LinkOutFlowCapacity[i][t] = 0;



        }

        for (int t = 0; t < g_number_of_intervals_in_min; ++t)
        {        // convert per hour capacity to per second and then to per simulation interval
        m_LinkCumulativeArrivalVector[i][t] = 0;
        m_LinkCumulativeDepartureVector[i][t] = 0;
        m_LinkTDWaitingTime[i][t] = 0;
        }


    }

    // for each link with  for link type code is 's'
    //reset the time-dependent capacity to zero
    //go through the records defined in timing_arc file
    //only enable m_LinkOutFlowCapacity[l][t] for the timestamp in the time window per time interval and reset the link TD travel time.

    for (unsigned li = 0; li < g_link_vector.size(); ++li)
    {
        if(g_link_vector[li].timing_arc_flag)
        {
            // reset for signalized links (not freeway links as type code != 'f' for the case of freeway workzones)
            // only for the loading period
            for (int t = 0; t < g_number_of_loading_intervals_in_sec; ++t)
            {
                m_LinkOutFlowCapacity[li][t] = 0;
                m_LinkOutFlowState[li][t] = 0;

            }
        }

    }

    for (unsigned l = 0; l < g_link_vector.size(); ++l)
    {
        if (g_link_vector[l].timing_arc_flag == true)  // with timing data
        {
            int number_of_cycles = (g_LoadingEndTimeInMin - g_LoadingStartTimeInMin) * 60 / max(1, g_link_vector[l].VDF_period[0].cycle_length);  // unit: seconds;

            int cycle_length = g_link_vector[l].VDF_period[0].cycle_length;
            int start_green_time = g_link_vector[l].VDF_period[0].start_green_time;
            int end_green_time = g_link_vector[l].VDF_period[0].end_green_time;

            if (end_green_time < start_green_time)
            {
                end_green_time += cycle_length;  // consider a looped end green time notation, e.g. 60-10 for cl = 100, then end green time should be 110. 
            }
            dtalog.output() << "signal timing data: link: cycle_length = " << cycle_length << 
                ",start_green_time = " << start_green_time << 
                ",end_green_time = " << end_green_time <<
                    endl;

            for (int cycle_no = 0; cycle_no < number_of_cycles; ++cycle_no)
            {
                int count = 0;

                // relative time horizon
                for (int t = cycle_no * cycle_length + start_green_time; t <= cycle_no * cycle_length + end_green_time; t += 1)
                {
                    // activate capacity for this time duration
                    m_LinkOutFlowCapacity[l][t] = g_link_vector[l].saturation_flow_rate * g_link_vector[l].number_of_lanes / 3600.0;

                    if (t % 2 == 0)
                        m_LinkOutFlowCapacity[l][t] = 1* g_link_vector[l].number_of_lanes;
                    if (t % 2 == 1)
                        m_LinkOutFlowCapacity[l][t] = 0;

                    m_LinkOutFlowState[l][t] = 1;

                    //residual = m_LinkOutFlowCapacity[l][t] - (int)(m_LinkOutFlowCapacity[l][t]);
                    ////RandomSeed is automatically updated.
                    //RandomSeed = (LCG_a * RandomSeed + LCG_c) % LCG_M;
                    //random_ratio = float(RandomSeed) / LCG_M;

                    //if (random_ratio < residual)
                    //    m_LinkOutFlowCapacity[l][t] = (int)(m_LinkOutFlowCapacity[l][t]) + 1;
                    //else
                    //    m_LinkOutFlowCapacity[l][t] = (int)(m_LinkOutFlowCapacity[l][t]);
                }
            }
        } //end f
    }

    dtalog.output() << "End of initializing time dependent capacity data." << endl;
}

void Assignment::DeallocateLinkMemory4Simulation()
{
	// g_fout << "deallocate 2D dynamic memory m_LinkOutFlowCapacity..." << endl;
    //if(m_LinkOutFlowCapacity)
    //    Deallocate2DDynamicArray(m_LinkOutFlowCapacity , g_number_of_links, g_number_of_intervals_in_sec);  //1
    
   // memory leak here: Xuesong: 11/20/2021                                                                                                         // 
// //if (m_LinkOutFlowState)
    //    Deallocate2DDynamicArray(m_LinkOutFlowState, g_number_of_links, g_number_of_intervals_in_sec);  //1
    // g_fout << "deallocate 2D dynamic memory m_LinkCumulativeArrivalVector..." << endl;
    if(m_LinkCumulativeArrivalVector)
        Deallocate2DDynamicArray(m_LinkCumulativeArrivalVector, g_number_of_links, g_number_of_intervals_in_min);  //2
	// g_fout << "deallocate 2D dynamic memory m_LinkCumulativeDepartureVector..." << endl;
    if (m_LinkCumulativeDepartureVector)
        Deallocate2DDynamicArray(m_LinkCumulativeDepartureVector, g_number_of_links, g_number_of_intervals_in_min); //3

    if(m_LinkCACount)
        Deallocate1DDynamicArray(m_LinkCACount, g_number_of_links);
    
    if (m_LinkCDCount)
        Deallocate1DDynamicArray(m_LinkCDCount, g_number_of_links);

    if (m_LinkTDWaitingTime)
        Deallocate2DDynamicArray(m_LinkTDWaitingTime, g_number_of_links, g_number_of_intervals_in_min); //4
}

void Assignment::STTrafficSimulation()
{
    //given p_agent->path_link_seq_no_vector path link sequence no for each agent
    int TotalCumulative_Arrival_Count = 0;
    int TotalCumulative_Departure_Count = 0;
    double TotalVehicleMileTraveled = 0;
    double TotalVehicleHourTraveled = 0;

    clock_t start_t;
    start_t = clock();

    AllocateLinkMemory4Simulation();

    int agent_type_size = g_AgentTypeVector.size();
    int zone_size = g_zone_vector.size();
    int demand_period_size = g_DemandPeriodVector.size();

    CColumnVector* p_column_pool;
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
                    p_column_pool = &(assignment.g_column_pool[orig][dest][at][tau]);
                    if (p_column_pool->od_volume > 0)
                    {
                        // scan through the map with different node sum for different continuous paths
                        it_begin = p_column_pool->path_node_sequence_map.begin();
                        it_end = p_column_pool->path_node_sequence_map.end();

                        for (it = it_begin; it != it_end; ++it)
                        {
                            path_toll = 0;
                            path_distance = 0;
                            path_travel_time = 0;

                            int VehicleSize = (it->second.path_volume +0.5);

                            for(int v = 0; v < VehicleSize; ++v)
                            {

                                int slot_no = assignment.g_DemandPeriodVector[tau].get_time_slot_no();

                                if (it->second.path_volume < 1)
                                    time_stamp = (slot_no - assignment.g_DemandPeriodVector[tau].starting_time_slot_no) * MIN_PER_TIMESLOT + g_agent_simu_vector.size() % 15 / 15.0 *1.0 * MIN_PER_TIMESLOT;
                                else
                                    time_stamp = (slot_no - assignment.g_DemandPeriodVector[tau].starting_time_slot_no) * MIN_PER_TIMESLOT + v *1.0 / it->second.path_volume*  MIN_PER_TIMESLOT;

                                if (it->second.m_link_size == 0)   // only load agents with physical path
                                    continue;

                                CAgent_Simu* pAgent = new CAgent_Simu();
                                // for future use of column pool
                                pAgent->at = at;
                                pAgent->dest = dest;
                                pAgent->tau = tau;

                                pAgent->PCE_unit_size = max(1, (int) (assignment.g_AgentTypeVector[at].PCE + 0.5));  // convert a possible floating point pce to an integer value for simulation
                                pAgent->time_headway = (int)(assignment.g_AgentTypeVector[at].time_headway_in_sec / number_of_seconds_per_interval + 0.5);
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
                                pAgent->m_Veh_LinkDepartureTime_in_simu_interval[0] = pAgent->m_Veh_LinkArrivalTime_in_simu_interval[0] + (int)(g_link_vector[FirstLink].free_flow_travel_time_in_min* number_of_simu_intervals_in_min);

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

    int number_of_in_min_for_RTSP = 5;

    int link_size = g_link_vector.size();

    // initialize CA and CD count for each link

    bool cell_based_simulation_mode = false;
    for (int li = 0; li < link_size; ++li)
    {
        CLink* pLink = &(g_link_vector[li]);

        if (pLink->cell_type >= 0)  // activate cell based simulation mode
            cell_based_simulation_mode = true;

            m_LinkCACount[li] = 0;
            m_LinkCDCount[li] = 0;

            g_link_vector[li].current_driving_AgentID = -1;
            g_link_vector[li].win_count = 0;
            g_link_vector[li].lose_count = 0;
            g_link_vector[li].time_to_be_released = -1;

    }


    bool bRealTimeInformationActivated = false;
    // first loop for time t
    for (int t = 0; t < g_number_of_simulation_intervals; ++t)
    {

        if(t% number_of_simu_intervals_in_min==0)  // every min
        {
            int time_in_min = t / number_of_simu_intervals_in_min;
            for (int li = 0; li < link_size; ++li)
            {
                CLink* pLink = &(g_link_vector[li]);
                {
                    m_LinkCumulativeArrivalVector[li][time_in_min] = m_LinkCACount[li];
                    m_LinkCumulativeDepartureVector[li][time_in_min] = m_LinkCDCount[li];
                }
            }
        }
        
        if (t % ( int(number_of_in_min_for_RTSP * 60/number_of_seconds_per_interval)) == 0)
        {
            int slot_no = (t * number_of_seconds_per_interval / 60 + g_LoadingStartTimeInMin) / 15;
            if (RTSP_RealTimeShortestPathFinding(slot_no, t)==true)
                bRealTimeInformationActivated = true;
        }

        int number_of_simu_interval_per_min = 60 / number_of_seconds_per_interval;
        double network_wide_speed = 60;
        double network_wide_travel_time = 0;

        if (TotalVehicleHourTraveled > 0.001)
        {
            network_wide_speed=  TotalVehicleMileTraveled / max(0.0001, TotalVehicleHourTraveled);
            network_wide_travel_time = TotalVehicleHourTraveled / max(1, TotalCumulative_Departure_Count)*60.0; // from hour to min
        }
         

        if (t % (number_of_simu_interval_per_min * 10) == 0)
            dtalog.output() << "simu time= " << t / number_of_simu_interval_per_min << " min, CA = " << TotalCumulative_Arrival_Count << " CD=" << TotalCumulative_Departure_Count
            << ", speed=" << network_wide_speed << ", travel time = " << network_wide_travel_time
            << endl;

        if (g_AgentTDListMap.find(t) != g_AgentTDListMap.end())
        {
            for (int Agent_v = 0; Agent_v < g_AgentTDListMap[t].m_AgentIDVector.size(); ++Agent_v)
            {
                int agent_id = g_AgentTDListMap[t].m_AgentIDVector[Agent_v];

                CAgent_Simu* p_agent = g_agent_simu_vector[agent_id];
                p_agent->m_bGenereated = true;
                int FirstLink = p_agent->path_link_seq_no_vector[0];
                m_LinkCACount[FirstLink] += 1;

                g_link_vector[FirstLink].EntranceQueue.push_back(p_agent->agent_id);

                g_link_vector[FirstLink].current_driving_AgentID = p_agent->agent_id;

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
                int arrival_time_in_simu_interval = p_agent->m_Veh_LinkArrivalTime_in_simu_interval[p_agent->m_current_link_seq_no];
                int link_travel_time_in_simu_interavls = max(1,g_link_vector[li].free_flow_travel_time_in_min* number_of_simu_intervals_in_min);
                p_agent->m_Veh_LinkDepartureTime_in_simu_interval[p_agent->m_current_link_seq_no] = arrival_time_in_simu_interval + g_link_vector[li].free_flow_travel_time_in_min* number_of_simu_intervals_in_min;

                //dtalog.output() << "reserve TD at time t = " << t << " on link" << g_node_vector[pLink->from_node_seq_no].node_id << " -> " << g_node_vector[pLink->to_node_seq_no].node_id <<
                //    " for agent = " << agent_id << " on link seq. no " << p_agent->m_current_link_seq_no << " with travel time in simu intervas = " << link_travel_time_in_simu_interavls << endl;

            }
        }

        /// for each link, for all vehicles in the queue 
        ///                         // back trace the resourece use vector right after the vehicle moves to the next link
        for (int li = 0; li < link_size; ++li)
        {
            CLink* pLink = &(g_link_vector[li]);
            for (auto it = pLink->ExitQueue.begin(); it != pLink->ExitQueue.end(); ++it)
            {
                int agent_id = (*it);
                CAgent_Simu* p_agent = g_agent_simu_vector[agent_id];

                if (p_agent->PCE_unit_size >= 2)
                {
                    for (int l_backtrace = 1; l_backtrace < p_agent->PCE_unit_size; l_backtrace++)
                    {
                        int local_l_index = p_agent->m_current_link_seq_no - l_backtrace;
                        if (local_l_index >= 0)
                        {
                            int current_link_seq_no = p_agent->path_link_seq_no_vector[local_l_index];
                            CLink* pCurrentLink = &(g_link_vector[current_link_seq_no]);
                            pCurrentLink->time_to_be_released = t + p_agent->time_headway;
                            // big truck/bus, 
                            //backtrace the previous l - k links, k = 1, 2, K, and set release time
                              //  K as the number of units of truck
                        }

                    }
                }
            }
        } // end of link 

                /// 

        /// 
        /// 
        /// </summary>
        int node_size = g_node_vector.size();

#pragma omp parallel for  // parallel computing for each node
        for (int node = 0; node < node_size; ++node)  // for each node
        {
            if (g_node_vector[node].node_id == 2347 && t/240>=18)
            {
                int debug = 1;
            }

            bool node_resource_competing_mode = false;
            int incoming_link_request_count = 0;
            int incoming_link_index_FIFO = -1;
            int debug_node_resource_competing_mode = 0;
            if(cell_based_simulation_mode)
            {
            std::map<int, int> next_link_for_resource_request;


            for (int i = 0; i < g_node_vector[node].m_incoming_link_seq_no_vector.size(); ++i)
            {
                int link = g_node_vector[node].m_incoming_link_seq_no_vector[i];  // we will start with different first link from the incoming link list,
                // equal change, regardless of # of lanes or main line vs. ramp, but one can use service arc, to control the effective capacity rates, e.g. through a metered ramp, to
                // allow mainline to use the remaining flow
                CLink* pLink = &(g_link_vector[link]);
                if( pLink->ExitQueue.size() >=1)
                {
                    int agent_id = pLink->ExitQueue.front();
                    CAgent_Simu* p_agent = g_agent_simu_vector[agent_id];
                    int next_link_seq_no = p_agent->path_link_seq_no_vector[p_agent->m_current_link_seq_no + 1];
                    CLink* pNextLink = &(g_link_vector[next_link_seq_no]);

                    int current_vehicle_count = m_LinkCACount[next_link_seq_no] - m_LinkCDCount[next_link_seq_no];
                    if (pNextLink->cell_type >=0 && current_vehicle_count < pNextLink->spatial_capacity_in_vehicles)  // only apply for cell mode
                    {
                        next_link_for_resource_request[next_link_seq_no] = 1;
                        incoming_link_request_count++;
                    }
                }
            }

            
            if(incoming_link_request_count>=2)
            {
                if (next_link_for_resource_request.size() == 1)
                {
//                if(g_node_vector[node].node_id == 2347 && t / 240 >= 18)
                    node_resource_competing_mode = true;

                }
            }
                // determine FIFO link with earliest departure time request
            if (node_resource_competing_mode)
            {
                int earlest_departure_time_interval = 99999999;

                for (int i = 0; i < g_node_vector[node].m_incoming_link_seq_no_vector.size(); ++i)
                {
                    int link = g_node_vector[node].m_incoming_link_seq_no_vector[i];  // we will start with different first link from the incoming link list,
                    // equal change, regardless of # of lanes or main line vs. ramp, but one can use service arc, to control the effective capacity rates, e.g. through a metered ramp, to
                    // allow mainline to use the remaining flow
                    CLink* pLink = &(g_link_vector[link]);
                    if (pLink->cell_type >=0 && pLink->ExitQueue.size() >= 1)
                    {
                        int agent_id = pLink->ExitQueue.front();
                        CAgent_Simu* p_agent = g_agent_simu_vector[agent_id];
                        if (p_agent->m_Veh_LinkDepartureTime_in_simu_interval[p_agent->m_current_link_seq_no] < earlest_departure_time_interval)
                        {
                            earlest_departure_time_interval = p_agent->m_Veh_LinkDepartureTime_in_simu_interval[p_agent->m_current_link_seq_no];
                            incoming_link_index_FIFO = i;  // i is the link index in the vector of m_incoming_link_seq_no_vector of node
                        }
                    }
                }
            }

            
            }
            // for each incoming link
            for (int i = 0; i < g_node_vector[node].m_incoming_link_seq_no_vector.size(); ++i)
            {
                int incoming_link_index = i;
                
                if (incoming_link_index_FIFO >= 0)
                {
                    incoming_link_index = (i + incoming_link_index_FIFO) % (g_node_vector[node].m_incoming_link_seq_no_vector.size());  // cycle loop
                    if(debug_node_resource_competing_mode)
                    dtalog.output() << "FIFO judgement i = " << i << endl;
                }else
                {  // random mode
                    incoming_link_index = (i + t) % (g_node_vector[node].m_incoming_link_seq_no_vector.size());
                }


                int link = g_node_vector[node].m_incoming_link_seq_no_vector[incoming_link_index];  // we will start with different first link from the incoming link list,
                // equal change, regardless of # of lanes or main line vs. ramp, but one can use service arc, to control the effective capacity rates, e.g. through a metered ramp, to
                // allow mainline to use the remaining flow
                CLink* pLink = &(g_link_vector[link]);

                // check if the current link has sufficient capacity
                // most critical and time-consuming task, check link outflow capacity

                int time_in_sec = t*number_of_seconds_per_interval;


                while (m_LinkOutFlowCapacity[link][time_in_sec] >= 1 && pLink->ExitQueue.size() >=1)
                {
                    int agent_id = pLink->ExitQueue.front();
                    CAgent_Simu* p_agent = g_agent_simu_vector[agent_id];

                    if (p_agent->m_Veh_LinkDepartureTime_in_simu_interval[p_agent->m_current_link_seq_no] > t)
                    {
                        // the future departure time on this link is later than the current time
                        break;
                    }

                    if (p_agent->m_current_link_seq_no == p_agent->path_link_seq_no_vector.size() - 2)  //-2 do not consider virtual destnation arc
                    {
                        // end of path
                        pLink->ExitQueue.pop_front();
                        pLink->current_driving_AgentID = -1;

                        p_agent->m_bCompleteTrip = true;
                        m_LinkCDCount[link] += 1;
                        p_agent->arrival_time_in_min = t*1.0/ number_of_simu_interval_per_min;
                        p_agent->path_travel_time_in_min = p_agent->arrival_time_in_min - p_agent->departure_time_in_min;

                        if (p_agent->agent_id == 0)
                        {
                            int debug = 1;
                        }
                        p_agent->path_distance = 0;

                        for (int link_s = 0; link_s < p_agent->path_link_seq_no_vector.size(); link_s++)
                        {
                            int link_seq_no = p_agent->path_link_seq_no_vector[link_s];
                            //path_toll += g_link_vector[link_seq_no].VDF_period[p_agent->tau].toll[at];
                            p_agent->path_distance += g_link_vector[link_seq_no].length;
                        }

                        #pragma omp critical
                        {
                            TotalCumulative_Departure_Count += 1;
                            TotalVehicleHourTraveled += p_agent->path_travel_time_in_min / 60.0;
                            TotalVehicleMileTraveled += p_agent->path_distance ;
                        }

                    }
                    else
                    {
                        //RT re-routing 
                        // test the one to all trees
                        // replease link sequence 
                        // contiue
                        // 
                        // 
                        // not complete the trip. move to the next link's entrance queue
                        int next_link_seq_no = p_agent->path_link_seq_no_vector[p_agent->m_current_link_seq_no + 1];
                        CLink* pNextLink = &(g_link_vector[next_link_seq_no]);

                        if (p_agent->m_current_link_seq_no >= 4)
                        {
                            int debug = 1;
                        }
                        //// spatial queue
                        ////if(pNextLink->length <=0.008)  // cell based link
                        ////{ 
                        //    pNextLink->spatial_capacity_in_vehicles = 1; // for cell based model
                        ////}

                        if (pNextLink->cell_type >= 0 || pNextLink->traffic_flow_code == 2)  // DTA micro cell based simulation
                         {
                                int current_vehicle_count = m_LinkCACount[next_link_seq_no] - m_LinkCDCount[next_link_seq_no];
                                if (current_vehicle_count >= pNextLink->spatial_capacity_in_vehicles  || (t < pNextLink->time_to_be_released ) ) 
                                    // this is an OR condition, related to reaction time tau
                            {
                                // spatical queue in the next link is blocked, break the while loop from here, as a first in first out queue.
								// g_fout << "spatical queue in the next link is blocked on link seq  " <<  g_node_vector[pNextLink->from_node_seq_no].node_id  << " -> " << g_node_vector[pNextLink->to_node_seq_no].node_id  <<endl;
                                //dtalog.output() << "spatical queue in the next link is blocked at time  = " << t << ",link" << g_node_vector[pLink->from_node_seq_no].node_id << " -> " << g_node_vector[pLink->to_node_seq_no].node_id <<
                                //    ",current_vehicle_count = " << current_vehicle_count << endl;
                                if (node_resource_competing_mode)
                                {
                                    pLink->lose_count++;
                                    if (debug_node_resource_competing_mode)
                                    { 
                                    dtalog.output() << "lose: the next link is blocked at time  = " << t << ", from link" << g_node_vector[pLink->from_node_seq_no].node_id << " -> " << g_node_vector[pLink->to_node_seq_no].node_id 
                                        << " to link" << g_node_vector[pNextLink->from_node_seq_no].node_id << " -> " << g_node_vector[pNextLink->to_node_seq_no].node_id
                                        << " time request " << p_agent->m_Veh_LinkDepartureTime_in_simu_interval[p_agent->m_current_link_seq_no] 
                                        << endl;
                                    }
                                }
                                break;  // being blocked
                            }
                            else // free to go if none of the above conditions cannot be met 
                            {
                                if (node_resource_competing_mode)
                                {
                                    pLink->win_count++;
                                    if (debug_node_resource_competing_mode)
                                    {
                                        dtalog.output() << "win:spatical queue in the next link is blocked at time  = " << t << ",link" << g_node_vector[pLink->from_node_seq_no].node_id << " -> " << g_node_vector[pLink->to_node_seq_no].node_id
                                            << " to link" << g_node_vector[pNextLink->from_node_seq_no].node_id << " -> " << g_node_vector[pNextLink->to_node_seq_no].node_id
                                            << " time request " << p_agent->m_Veh_LinkDepartureTime_in_simu_interval[p_agent->m_current_link_seq_no]
                                            << endl;
                                    }
                                }


                            }

                        }

        //                // kinematic wave
        //                if (pNextLink->traffic_flow_code == 3)  // to be discussed later with Cafer
        //                {
        //                    int lagged_time_stamp = max(0, t - 1 - pNextLink->BWTT_in_simulation_interval);
        //                    int current_vehicle_count = m_LinkCACount[next_link_seq_no] - m_LinkCDCount[next_link_seq_no];
        //                    if (current_vehicle_count > pNextLink->spatial_capacity_in_vehicles)
        //                    {
        //                        // spatical queue in the next link is blocked, break the while loop from here, as a first in first out queue.
								//// g_fout << "spatical queue in the next link is blocked on link seq  " <<  g_node_vector[pNextLink->from_node_seq_no].node_id  << " -> " << g_node_vector[pNextLink->to_node_seq_no].node_id  <<endl;
        //                        break;
        //                    }
        //                }

                        pLink->ExitQueue.pop_front();
                        pLink->current_driving_AgentID = -1;
                        pLink->time_to_be_released = t + p_agent->time_headway;

                        // back trace the resourece use vector right after the vehicle moves to the next link
                        if (p_agent->PCE_unit_size >= 2)
                        {
                            for (int l_backtrace = 1; l_backtrace < p_agent->PCE_unit_size; l_backtrace++)
                            {
                                int local_l_index = p_agent->m_current_link_seq_no - l_backtrace;
                                if(local_l_index>=0)
                                {
                                int current_link_seq_no = p_agent->path_link_seq_no_vector[local_l_index];
                                CLink* pCurrentLink = &(g_link_vector[current_link_seq_no]);
                                pCurrentLink->time_to_be_released = t + p_agent->time_headway; 
                                // big truck/bus, 
                                //backtrace the previous l - k links, k = 1, 2, K, and set release time
                                  //  K as the number of units of truck
                                }

                            }

                        }


                        pNextLink->EntranceQueue.push_back(agent_id);
                        pNextLink->current_driving_AgentID = agent_id;

                        p_agent->m_Veh_LinkDepartureTime_in_simu_interval[p_agent->m_current_link_seq_no] = t;
                        p_agent->m_Veh_LinkArrivalTime_in_simu_interval[p_agent->m_current_link_seq_no + 1] = t;

                        float travel_time_in_sec = (p_agent->m_Veh_LinkDepartureTime_in_simu_interval[p_agent->m_current_link_seq_no] - p_agent->m_Veh_LinkArrivalTime_in_simu_interval[p_agent->m_current_link_seq_no])* number_of_seconds_per_interval;
                        //for each waited vehicle
                        float waiting_time_in_min = travel_time_in_sec/60.0 - pLink->free_flow_travel_time_in_min;

                        m_LinkTDWaitingTime[link][p_agent->m_Veh_LinkArrivalTime_in_simu_interval[p_agent->m_current_link_seq_no] / number_of_simu_intervals_in_min] += waiting_time_in_min;
                        m_LinkTotalWaitingTimeVector[link] += waiting_time_in_min;
                        m_LinkCDCount[link] += 1;
                        m_LinkCACount[next_link_seq_no] += 1;

                        //test rerouting
                        if(bRealTimeInformationActivated)
                            UpdateRTPath(p_agent);
                    }

                    //move
                    p_agent->m_current_link_seq_no += 1;
                    m_LinkOutFlowCapacity[link][time_in_sec] -= 1;//deduct the capacity
                }
            }
        } // conditions
    }  // departure time events
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
            string measurement_type;
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
                    if (g_link_vector[link_seq_no].obs_count >= 1)  // data exist
                    {
                        if(upper_bound_flag==0)
                        {  // over write only if the new data are acutal counts, 
                            
                        g_link_vector[link_seq_no].obs_count = count;
                        g_link_vector[link_seq_no].upper_bound_flag = upper_bound_flag;
                        }
                        else  // if the new data are upper bound, skip it and keep the actual counts 
                        {
                            // do nothing 
                        }


                    }
                    else
                    {
                    g_link_vector[link_seq_no].obs_count = count;
                    g_link_vector[link_seq_no].upper_bound_flag = upper_bound_flag;
                    }
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

    // automatically add ultimate period capacity as the upper bound of link flow rates  

        //g_link_vector[link_seq_no].obs_count = count;
        //g_link_vector[link_seq_no].upper_bound_flag = upper_bound_flag;

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
        float total_system_travel_cost = 0;
        //step 2.1
        // we can have a recursive formulat to reupdate the current link volume by a factor of k/(k+1),
        // and use the newly generated path flow to add the additional 1/(k+1)
        double system_gap = 0;
        double gap = g_reset_and_update_link_volume_based_on_ODME_columns(g_link_vector.size(),s, system_gap);
        //step 2.2: based on newly calculated path volumn, update volume based travel time, and update volume based measurement error/deviation
                // and use the newly generated path flow to add the additional 1/(k+1)

        double gap_improvement = gap - prev_gap;

        if (s >= 1 && gap_improvement > 0.001)  // convergency criterion
            break;

        prev_gap = gap;

        int column_pool_counts = 0;
        int column_path_counts = 0;
        int column_pool_with_sensor_counts = 0;
        int column_path_with_sensor_counts = 0;

        //step 3: calculate shortest path at inner iteration of column flow updating
#pragma omp parallel for
        for (int orig = 0; orig < g_zone_vector.size(); ++orig)  // o
        {
            CColumnVector* p_column_pool;
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
                    for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); ++tau)  //tau
                    {
                        p_column_pool = &(assignment.g_column_pool[orig][dest][at][tau]);
                        if (p_column_pool->od_volume > 0)
                        {
                            column_pool_counts++;

                            column_vector_size = p_column_pool->path_node_sequence_map.size();
                            path_seq_count = 0;

                            it_begin = p_column_pool->path_node_sequence_map.begin();
                            it_end = p_column_pool->path_node_sequence_map.end();
                            int i = 0;
                            for (it = it_begin; it != it_end; ++it, ++i) // for each k
                            {
                                column_path_counts++;

                                if (s>=1 && it->second.measurement_flag == 0)  // after 1  iteration, if there are no data passing through this path column. we will skip it in the ODME process
                                    continue; 

                                path_gradient_cost = 0;
                                path_distance = 0;
                                path_travel_time = 0;
                                p_column_pool->m_passing_sensor_flag = 0;
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
                                    p_column_pool->m_passing_sensor_flag += 1;
                                    it->second.measurement_flag = 1;
                                    
                                }

                                // step 3.2 destination attraction flow gradient

                                if (g_zone_vector[dest].obs_attraction > 0)
                                {
                                    if (g_zone_vector[orig].obs_attraction_upper_bound_flag == 0)
                                        path_gradient_cost += g_zone_vector[dest].est_attraction_dev;

                                    if (g_zone_vector[orig].obs_attraction_upper_bound_flag == 1 && g_zone_vector[dest].est_attraction_dev > 0)
                                        path_gradient_cost += g_zone_vector[dest].est_attraction_dev;

                                    p_column_pool->m_passing_sensor_flag += 1;
                                    it->second.measurement_flag = 1;
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
                                        {// we only consider the over capaity value here to penalize the path flow 
                                            path_gradient_cost += g_link_vector[link_seq_no].est_count_dev;
                                            est_count_dev += g_link_vector[link_seq_no].est_count_dev;
                                        }
                                        p_column_pool->m_passing_sensor_flag += 1;
                                        it->second.measurement_flag = 1;
                                    }
                                }

                                // statistics collection 

                                if (it->second.measurement_flag >= 1)
                                    column_path_with_sensor_counts++;

                                it->second.path_gradient_cost = path_gradient_cost;

                                step_size = 0.01;
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
                            }  // end of loop for all paths in the column pools
                            
                            if (p_column_pool->m_passing_sensor_flag >= 1)
                                column_pool_with_sensor_counts++;

                        }
                    }
                }
            }
        }
        if (s == 0)
        {
            float percentage_of_OD_columns_with_sensors = column_pool_with_sensor_counts*1.0 / max(1, column_pool_counts) * 100;
            float percentage_of_paths_with_sensors = column_path_with_sensor_counts * 1.0 / max(1, column_path_counts) * 100;
            dtalog.output() << "count of all column pool vectors=" << column_pool_counts << ", "
                << "count of all paths =" << column_path_counts << ", "
                << "count of column_pools with sensors = " << column_pool_with_sensor_counts << "(" << percentage_of_OD_columns_with_sensors << "%), "
                << "count of column_paths with sensors = " << column_path_with_sensor_counts << " (" << percentage_of_paths_with_sensors << "%)" << endl;

        }

    }


    //if (assignment.g_pFileDebugLog != NULL)
    //	fprintf(assignment.g_pFileDebugLog, "CU: iteration %d: total_gap=, %f,total_relative_gap=, %f,\n", s, total_gap, total_gap / max(0.0001, total_system_travel_cost));
}