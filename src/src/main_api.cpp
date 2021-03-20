/* Portions Copyright 2019 Xuesong Zhou
 *
 * If you help write or modify the code, please also list your names here.
 * The reason of having Copyright info here is to ensure all the modified version, as a whole, under the GPL
 * and further prevent a violation of the GPL.
 *
 * More about "How to use GNU licenses for your own software"
 * http://www.gnu.org/licenses/gpl-howto.html 
 */
#ifdef _WIN32
#include "pch.h"
#endif

// Peiheng, 02/03/21, remove them later
#pragma warning(disable : 4305 4267 4018) 
// stop warning: "conversion from 'int' to 'float', possible loss of data"
#pragma warning(disable: 4244)

#include <iostream>
#include <fstream>
#include <list> 
#include <omp.h>
#include <algorithm>
#include <time.h>
#include <functional>
#include <cstdio>   
#include <cmath>
#include <stack>
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <iomanip>
#include <cstring>
#include "utils.h"
#include "config.h"

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
using std::ostringstream;

// some basic parameters setting

//Pls make sure the _MAX_K_PATH > Agentlite.cpp's g_number_of_column_generation_iterations+g_reassignment_number_of_K_paths and the _MAX_ZONE remain the same with .cpp's defination
constexpr auto _MAX_LABEL_COST = 1.0e+15;

constexpr auto _MAX_AGNETTYPES = 4; //because of the od demand store format,the MAX_demandtype must >=g_DEMANDTYPES.size()+1;
constexpr auto _MAX_TIMEPERIODS = 4; // time period set to 4: mid night, morning peak, mid-day and afternoon peak;
constexpr auto _MAX_MEMORY_BLOCKS = 20;

constexpr auto _MAX_LINK_SIZE_IN_A_PATH = 1000;
constexpr auto _MAX_LINK_SIZE_FOR_A_NODE = 200;

constexpr auto _MAX_TIMESLOT_PerPeriod = 100; // max 96 15-min slots per day
constexpr auto _default_saturation_flow_rate = 1530;

constexpr auto MIN_PER_TIMESLOT = 15;

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

FILE* g_pFileOutputLog = nullptr;
ofstream g_fout("log.txt");

int g_debug_level = 0;
int g_log_odme = 0;
int g_log_path = 0;

void fopen_ss(FILE** file, const char* fileName, const char* mode)
{
    *file = fopen(fileName, mode);
}

void g_ProgramStop()
{
    g_fout << "STALite Program stops. Press any key to terminate. Thanks!" << endl;
    getchar();
    exit(0);
}

//below shows where the functions used in Agentlite.cpp come from!
//Utility.cpp

void CCSVParser::ConvertLineStringValueToIntegers()
{
    LineIntegerVector.clear();
    for (unsigned i = 0; i < LineFieldsValue.size(); i++)
    {
        std::string si = LineFieldsValue[i];
        int value = atoi(si.c_str());

        if (value >= 1)
            LineIntegerVector.push_back(value);
    }
}

bool CCSVParser::OpenCSVFile(std::string fileName, bool b_required)
{

    mFileName = fileName;
    inFile.open(fileName.c_str());

    if (inFile.is_open())
    {
        if (IsFirstLineHeader)
        {
            std::string s;
            std::getline(inFile, s);
            std::vector<std::string> FieldNames = ParseLine(s);

            for (size_t i = 0;i < FieldNames.size();i++)
            {
                std::string tmp_str = FieldNames.at(i);
                size_t start = tmp_str.find_first_not_of(" ");

                std::string name;
                if (start == std::string::npos)
                {
                    name = "";
                }
                else
                {
                    name = tmp_str.substr(start);
                    //TRACE("%s,", name.c_str());
                }
                FieldsIndices[name] = (int)i;
            }
        }
        return true;
    }
    else
    {
        if (b_required)
        {
            g_fout << "File " << fileName << " does not exist. Please check." << std::endl;
            //g_ProgramStop();
        }
        return false;
    }
}

bool CCSVParser::ReadRecord()
{
    LineFieldsValue.clear();

    if (inFile.is_open())
    {
        std::string s;
        std::getline(inFile, s);
        if (s.length() > 0)
        {
            LineFieldsValue = ParseLine(s);
            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        return false;
    }
}


bool CCSVParser::ReadSectionHeader(std::string s)
{
    //skip // data 

    Headers.clear();
    FieldsIndices.clear();

    if (s.length() == 0)
        return true;

    std::vector<std::string> FieldNames = ParseLine(s);

    for (size_t i = 0; i < FieldNames.size(); i++)
    {
        std::string tmp_str = FieldNames.at(i);
        size_t start = tmp_str.find_first_not_of(" ");

        std::string name;
        if (start == std::string::npos)
        {
            name = "";
        }
        else
        {
            name = tmp_str.substr(start);
        }
        Headers.push_back(name);
        FieldsIndices[name] = (int)i;
    }
    return true;
}

bool CCSVParser::ReadRecord_Section()
{
    LineFieldsValue.clear();

    if (inFile.is_open())
    {
        std::string s;
        std::getline(inFile, s);
        if (s.length() > 0)
        {
            if(s.find("[") != std::string::npos)  // synchro single csv file
            {
                LineFieldsValue = ParseLine(s);

                if (LineFieldsValue.size() >= 1)
                {
                    SectionName = LineFieldsValue[0];
                }

                //re-read section header
                ReadSectionHeader(s);
                std::getline(inFile, s);

            }
            LineFieldsValue = ParseLine(s);
            return true;
        }
        else
        {
            if (m_bLastSectionRead)  // reach the last section
                return false;
            else
            {
                if (inFile.eof())
                    return false;
                else
                    return true;
            }
            
        }
    }
    else
    {
        return false;
    }
}

std::vector<std::string> CCSVParser::ParseLine(std::string line)
{
    std::vector<std::string> SeperatedStrings;
    std::string subStr;

    if (line.length() == 0)
        return SeperatedStrings;

    std::istringstream ss(line);

    if (line.find_first_of('"') == std::string::npos)
    {
        while (std::getline(ss, subStr, Delimiter))
        {
            SeperatedStrings.push_back(subStr);
        }

        if (line.at(line.length() - 1) == ',')
        {
            SeperatedStrings.push_back("");
        }
    }
    else
    {
        while (line.length() > 0)
        {
            size_t n1 = line.find_first_of(',');
            size_t n2 = line.find_first_of('"');

            if (n1 == std::string::npos && n2 == std::string::npos) //last field without double quotes
            {
                subStr = line;
                SeperatedStrings.push_back(subStr);
                break;
            }

            if (n1 == std::string::npos && n2 != std::string::npos) //last field with double quotes
            {
                size_t n3 = line.find_first_of('"', n2 + 1); // second double quote

                //extract content from double quotes
                subStr = line.substr(n2 + 1, n3 - n2 - 1);
                SeperatedStrings.push_back(subStr);

                break;
            }

            if (n1 != std::string::npos && (n1 < n2 || n2 == std::string::npos))
            {
                subStr = line.substr(0, n1);
                SeperatedStrings.push_back(subStr);
                if (n1 < line.length() - 1)
                {
                    line = line.substr(n1 + 1);
                }
                else // comma is the last char in the line string, push an empty string to the back of vector
                {
                    SeperatedStrings.push_back("");
                    break;
                }
            }

            if (n1 != std::string::npos && n2 != std::string::npos && n2 < n1)
            {
                size_t n3 = line.find_first_of('"', n2 + 1); // second double quote
                subStr = line.substr(n2 + 1, n3 - n2 - 1);
                SeperatedStrings.push_back(subStr);
                size_t idx = line.find_first_of(',', n3 + 1);

                if (idx != std::string::npos)
                {
                    line = line.substr(idx + 1);
                }
                else
                {
                    break;
                }
            }
        }

    }
    return SeperatedStrings;
}

template <class T> bool CCSVParser::GetValueByFieldName(std::string field_name, T& value, bool required_field, bool NonnegativeFlag)
{

    if (FieldsIndices.find(field_name) == FieldsIndices.end())
    {
        if (required_field)
        {
            g_fout << "Field " << field_name << " in file " << mFileName << " does not exist. Please check the file." << endl;

            g_ProgramStop();
        }
        return false;
    }
    else
    {
        if (LineFieldsValue.size() == 0)
        {
            return false;
        }

        int size = (int)(LineFieldsValue.size());
        if (FieldsIndices[field_name] >= size)
        {
            return false;
        }

        string str_value = LineFieldsValue[FieldsIndices[field_name]];

        if (str_value.length() <= 0)
        {
            return false;
        }

        istringstream ss(str_value);

        T converted_value;
        ss >> converted_value;

        if (/*!ss.eof() || */ ss.fail())
        {
            return false;
        }

        if (required_field)
        {
            if(NonnegativeFlag)
            { 
            if (converted_value < 0)
                converted_value = 0;
            }
        }

        value = converted_value;
        return true;
    }
}

bool CCSVParser::GetValueByFieldName(std::string field_name, std::string& value, bool required_field)
{
    if (FieldsIndices.find(field_name) == FieldsIndices.end())
    {
        if (required_field)
        {
            g_fout << "Field " << field_name << " in file " << mFileName << " does not exist. Please check the file." << std::endl;

            g_ProgramStop();
        }
        return false;
    }
    else
    {
        if (LineFieldsValue.size() == 0)
        {
            return false;
        }

        unsigned int index = FieldsIndices[field_name];
        if (index >= LineFieldsValue.size())
        {
            return false;
        }
        std::string str_value = LineFieldsValue[index];

        if (str_value.length() <= 0)
        {
            return false;
        }

        value = str_value;
        return true;
    }
}

template <typename T>
T** AllocateDynamicArray(int nRows, int nCols)
{
    T** dynamicArray;

    dynamicArray = new (std::nothrow) T*[nRows];

    if (!dynamicArray)
    {
        g_fout << "Error: insufficient memory.";
        g_ProgramStop();
    }

    for (int i = 0; i < nRows; ++i)
    {
        dynamicArray[i] = new (std::nothrow) T[nCols];

        if (!dynamicArray[i])
        {
            g_fout << "Error: insufficient memory.";
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
        g_fout << "Error: insufficient memory.";
        g_ProgramStop();
    }

    for (int x = 0; x < nX; ++x)
    {
        if (x % 1000 == 0)
        {
            g_fout << "allocating 3D memory for " << x << endl;
        }

        dynamicArray[x] = new (std::nothrow) T*[nY];

        if (!dynamicArray[x])
        {
            g_fout << "Error: insufficient memory.";
            g_ProgramStop();
        }

        for (int y = 0; y < nY; ++y)
        {
            dynamicArray[x][y] = new (std::nothrow) T[nZ];
            if (!dynamicArray[x][y])
            {
                g_fout << "Error: insufficient memory.";
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
        g_fout << "Error: insufficient memory.";
        g_ProgramStop();
    }

    for (int m = 0; m < nM; ++m)
    {
        if (m % 1000 == 0)
            g_fout << "allocating 4D memory for " << m << " zones" << endl;

        dynamicArray[m] = new (std::nothrow) T**[nX];

        if (!dynamicArray[m])
        {
            g_fout << "Error: insufficient memory.";
            g_ProgramStop();
        }

        for (int x = 0; x < nX; ++x)
        {
            dynamicArray[m][x] = new (std::nothrow) T*[nY];

            if (!dynamicArray[m][x])
            {
                g_fout << "Error: insufficient memory.";
                g_ProgramStop();
            }

            for (int y = 0; y < nY; ++y)
            {
                dynamicArray[m][x][y] = new (std::nothrow) T[nZ];
                if (!dynamicArray[m][x][y])
                {
                    g_fout << "Error: insufficient memory.";
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

class CDemand_Period {
public:

    CDemand_Period() : demand_period{ 0 }, starting_time_slot_no{ 0 }, ending_time_slot_no{ 0 }
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
    CAgent_type() : agent_type_no{ 1 }, value_of_time{ 1 }
    {
    }

    int agent_type_no;
    // dollar per hour
    float value_of_time;
    // link type, product consumption equivalent used, for travel time calculation
    float PCE;
    string agent_type;
};

class CLinkType
{
public:
    CLinkType() : link_type{ 1 }, number_of_links{ 0 }, traffic_flow_code{ 0 }
    {
    }

    bool AllowAgentType(string agent_type)
    {
        if (agent_type_blocklist.size() == 0)  // if the agent_type_blocklist is empty then all types are allowed.
            return true;
        else
        {
            if (agent_type_blocklist.find(agent_type) == string::npos)  // otherwise, only an agent type is listed in this "block list", then this agent is allowed to travel on this link
                return true;
            else
            {
                g_fout << "important: agent_type " << agent_type << " is prohibited " << " on link type " << link_type << endl;
                return false;
            }
        }
    }

    int link_type;
    int number_of_links;
    int traffic_flow_code;

    string link_type_name;
    string agent_type_blocklist;
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
            delete path_node_vector;
            delete path_link_vector;
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
    float path_volume;
    float path_switch_volume;
    float path_travel_time;
    float path_distance;
    float path_toll;
    // first order graident cost.
    float path_gradient_cost;
    // first order graident cost - least gradient cost
    float path_gradient_cost_difference;
    // first order graident cost - least gradient cost
    float path_gradient_cost_relative_difference;

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
    float od_volume;  
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
            delete m_Veh_LinkArrivalTime_in_simu_interval;
        if (m_Veh_LinkDepartureTime_in_simu_interval) 
            delete m_Veh_LinkDepartureTime_in_simu_interval;

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
    Assignment() : assignment_mode{ 0 }, g_number_of_memory_blocks{ 4 }, g_number_of_threads{ 1 }, g_link_type_file_loaded{ true }, g_agent_type_file_loaded{ false },
        total_demand_volume{ 0.0 }, g_origin_demand_array{ nullptr }, g_column_pool{ nullptr }, g_number_of_in_memory_simulation_intervals{ 500 },
        g_number_of_column_generation_iterations{ 20 }, g_number_of_demand_periods{ 24 },g_number_of_links{ 0 }, g_number_of_timing_arcs{ 0 },
        g_number_of_nodes{ 0 }, g_number_of_zones{ 0 }, g_number_of_agent_types{ 0 }, g_reassignment_tau0{ 999 }, debug_detail_flag{ 1 }
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
    void AllocateLinkMemory4Simulation();
    void DeallocateLinkMemory4Simulation();

    int assignment_mode;
    int g_number_of_memory_blocks;
    int g_number_of_threads;
    
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
        alpha{ 0.15f }, beta{ 4 }, rho{ 1 }, marginal_base{ 1 },
        starting_time_slot_no{ 0 }, ending_time_slot_no{ 0 },
        cycle_length{ 29 }, red_time{ 0 }, t0{ 0 }, t3{ 0 }
    {
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
        float mu = _default_saturation_flow_rate; //default saturation flow rates
        float s_bar = 1.0 / 60.0 * red * red / (2*cycle_length); // 60.0 is used to convert sec to min
        float uniform_delay = s_bar / max(1 - lambda / mu, 0.1f) ; 

        return uniform_delay;
    }

    float PerformBPR(float volume)
    {
        // take nonnegative values
        volume = max(0.0f, volume);  

        // Peiheng, 02/02/21, useless block
        if (volume > 1.0)
        {
            int debug = 1;
        }

        VOC = volume / max(0.00001f, capacity);
        avg_travel_time = FFTT + FFTT * alpha * pow(volume / max(0.00001f, capacity), beta);
        marginal_base = FFTT * alpha * beta*pow(volume / max(0.00001f, capacity), beta - 1);

        return avg_travel_time;
        // volume --> avg_traveltime
    }
    
    // input period based volume
    float PerformBPR_X(float volume)
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
            for (int tt_relative = 0; tt_relative <= L; tt_relative++)
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

                if (time_abs >= t0 && time_abs <= t3 && t3>t0)
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
                // g_fout << avg_waiting_time << endl;
                // avg_travel_time = (FFTT_in_hour + avg_waiting_time)*60.0; // avg_travel time should be per min 
            }
        }

        return avg_travel_time;
    }
    
    float m;
    // we should also pass uncongested_travel_time as length/(speed_at_capacity)
    float VOC;
    //updated BPR-X parameters
    float gamma;
    float mu;
    //peak hour factor
    float PHF;
    //standard BPR parameter 
    float alpha;
    float beta;
    float rho;
    float marginal_base;
    // in 15 min slot
    int starting_time_slot_no;
    int ending_time_slot_no;
    float cycle_length;
    float red_time;
    int t0, t3;

    bool bValidQueueData;
    string period;

    float capacity;
    float FFTT;

    float congestion_period_P;
    // inpput
    float volume;

    //output
    float avg_delay;
    float avg_travel_time = 0;
    float avg_waiting_time = 0;

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
    CLink() :main_node_id{ -1 }, NEMA_phase_number{ -1 }, obs_count{ -1 }, upper_bound_flag{ 0 }, est_count_dev{ 0 }, 
        BWTT_in_simulation_interval{ 100 }, zone_seq_no_for_outgoing_connector{ -1 }, number_of_lanes{ 1 }, lane_capacity{ 1999 },
        length{ 1 }, free_flow_travel_time_in_min{ 1 }, toll{ 0 }, route_choice_cost{ 0 }, link_spatial_capacity{ 100 },
        service_arc_flag{ false }, traffic_flow_code{ 0 }, spatial_capacity_in_vehicles{ 999999 }
    {
        for (int tau = 0; tau < _MAX_TIMEPERIODS; ++tau)
        {
            flow_volume_per_period[tau] = 0;
            queue_length_perslot[tau] = 0;
            travel_time_per_period[tau] = 0;

            for(int at = 0; at < _MAX_AGNETTYPES; ++at)
                volume_per_period_per_at[tau][at] = 0;

            TDBaseTT[tau] = 0;
            TDBaseCap[tau] = 0;
            TDBaseFlow[tau] = 0;
            TDBaseQueue[tau] = 0;
            //cost_perhour[tau] = 0;
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
        return (flow_volume_per_period[tau] + TDBaseFlow[tau]) / max(0.00001f, TDBaseCap[tau]);
    }


    float get_speed(int tau)
    {
        return length / max(travel_time_per_period[tau], 0.0001f) * 60;  // per hour
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
        float generalized_cost = travel_time_per_period[tau] + toll / assignment.g_AgentTypeVector[agent_type_no].value_of_time * 60;

        // system optimal mode or exterior panalty mode
        if (assignment.assignment_mode == 4)
            generalized_cost += travel_marginal_cost_per_period[tau][agent_type_no];

        return generalized_cost;
    }
    
    int main_node_id;

    int NEMA_phase_number;
    float obs_count;
    int upper_bound_flag;
    float est_count_dev;

    int BWTT_in_simulation_interval;
    int zone_seq_no_for_outgoing_connector;

    int number_of_lanes;
    float lane_capacity;
    float length;
    float free_flow_travel_time_in_min;
    float toll;
    float route_choice_cost;
    float link_spatial_capacity;

    bool service_arc_flag;
    int traffic_flow_code;
    int spatial_capacity_in_vehicles;

    // 1. based on BPR. 

    int link_seq_no;
    string link_id;
    string geometry;

    int from_node_seq_no;
    int to_node_seq_no;
    int link_type;

    string movement_str;

    float PCE;
    float fftt;

    CVDF_Period VDF_period[_MAX_TIMEPERIODS];

    float TDBaseTT[_MAX_TIMEPERIODS];
    float TDBaseCap[_MAX_TIMEPERIODS];
    float TDBaseFlow[_MAX_TIMEPERIODS];
    float TDBaseQueue[_MAX_TIMEPERIODS];

    int type;

    //static
    //float flow_volume;
    //float travel_time;

    float flow_volume_per_period[_MAX_TIMEPERIODS];
    float volume_per_period_per_at[_MAX_TIMEPERIODS][_MAX_AGNETTYPES];

    float queue_length_perslot[_MAX_TIMEPERIODS];  // # of vehicles in the vertical point queue
    float travel_time_per_period[_MAX_TIMEPERIODS];
    float travel_marginal_cost_per_period[_MAX_TIMEPERIODS][_MAX_AGNETTYPES];

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
    CNode() : zone_id{ -1 }, zone_org_id{ -1 }, prohibited_movement_size{ 0 }, node_seq_no{ -1 }
    {
    }

    //int accessible_node_count;

    int zone_id;
    // original zone id for non-centriod nodes
    int zone_org_id; 
    int prohibited_movement_size;
    // sequence number 
    int node_seq_no;

    //external node number 
    int node_id;

    double x;
    double y;

    std::vector<int> m_outgoing_link_seq_no_vector;
    std::vector<int> m_incoming_link_seq_no_vector;

    std::vector<int> m_to_node_seq_no_vector;
    std::map<int, int> m_to_node_2_link_seq_no_map;

    std::map<string, int> m_prohibited_movement_string_map;
};


std::vector<CNode> g_node_vector;
std::vector<CLink> g_link_vector;
std::vector<CServiceArc> g_service_arc_vector;

class COZone
{
public:
    COZone() : obs_production{ -1 }, obs_attraction{ -1 },
        est_production{ -1 }, est_attraction{ -1 },
        est_production_dev{ 0 }, est_attraction_dev{ -1 }
    {
    }

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

    int* OutgoingLinkNoArray;
    int* OutgoingNodeNoArray;
    int OutgoingLinkSize;
};

class NetworkForSP  // mainly for shortest path calculation
{
public:
    NetworkForSP() : temp_path_node_vector_size{ 1000 }, m_value_of_time{ 10 }, bBuildNetwork{ false }, m_memory_block_no{ 0 }
    {
    }

    int temp_path_node_vector_size;
    float m_value_of_time;
    bool bBuildNetwork;
    int m_memory_block_no;

    //node seq vector for each ODK
    int temp_path_node_vector[1000];
    //node seq vector for each ODK
    int temp_path_link_vector[1000];

    bool m_bSingleSP_Flag;

    // assigned nodes for computing 
    std::vector<int> m_origin_node_vector; 
    std::vector<int> m_origin_zone_seq_no_vector;

    int tau; // assigned nodes for computing 
    int m_agent_type_no; // assigned nodes for computing 

    CNodeForwardStar* NodeForwardStarArray;

    int m_threadNo;  // internal thread number 

    int m_ListFront; // used in coding SEL
    int m_ListTail;  // used in coding SEL
    int* m_SENodeList; // used in coding SEL

    // label cost for shortest path calcuating
    float* m_node_label_cost;
    // time-based cost
    float* m_label_time_array;  
    // distance-based cost
    float* m_label_distance_array; 

    // predecessor for nodes
    int* m_node_predecessor;
    // update status 
    int* m_node_status_array;
    // predecessor for this node points to the previous link that updates its label cost (as part of optimality condition) (for easy referencing)
    int* m_link_predecessor;

    float* m_link_flow_volume_array; 
    
    float* m_link_genalized_cost_array;
    int* m_link_outgoing_connector_zone_seq_no_array;

    // major function 1:  allocate memory and initialize the data 
    void AllocateMemory(int number_of_nodes, int number_of_links)
    {
        NodeForwardStarArray = new CNodeForwardStar[number_of_nodes]; 

        m_SENodeList = new int[number_of_nodes];  //1

        m_LinkBasedSEList = new int[number_of_links];  //1;  // dimension: number of links

        m_node_status_array = new int[number_of_nodes];  //2
        m_label_time_array = new float[number_of_nodes];  //3
        m_label_distance_array = new float[number_of_nodes];  //4
        m_node_predecessor = new int[number_of_nodes];  //5
        m_link_predecessor = new int[number_of_nodes];  //6
        m_node_label_cost = new float[number_of_nodes];  //7

        m_link_flow_volume_array = new float[number_of_links];  //8

        m_link_genalized_cost_array = new float[number_of_links];  //9
        m_link_outgoing_connector_zone_seq_no_array = new int[number_of_links]; //10
    }

    ~NetworkForSP()
    {
        if (m_SENodeList)  //1
            delete m_SENodeList;

        if (m_node_status_array)  //2
            delete m_node_status_array;

        if (m_label_time_array)  //3
            delete m_label_time_array;

        if (m_label_distance_array)  //4
            delete m_label_distance_array;

        if (m_node_predecessor)  //5
            delete m_node_predecessor;

        if (m_link_predecessor)  //6
            delete m_link_predecessor;

        if (m_node_label_cost)  //7
            delete m_node_label_cost;

        if (m_link_flow_volume_array)  //8
            delete m_link_flow_volume_array;

        if (m_link_genalized_cost_array) //9
            delete m_link_genalized_cost_array;
        
        if (m_link_outgoing_connector_zone_seq_no_array) //10
            delete m_link_outgoing_connector_zone_seq_no_array;

        // Peiheng, 02/02/21, useless block
        // Initialization for all non-origin nodes
        for (int i = 0; i < assignment.g_number_of_nodes; ++i)
        {
            if (NodeForwardStarArray[i].OutgoingLinkSize > 0)
            {
                //delete NodeForwardStarArray[i].OutgoingLinkNoArray;
                //delete NodeForwardStarArray[i].OutgoingNodeNoArray;
            }
        }

        if (NodeForwardStarArray)
            delete NodeForwardStarArray;
    }

    void UpdateGeneralizedLinkCost()
    {
        for (int i = 0; i < g_link_vector.size(); ++i)
        {
            CLink* pLink = &(g_link_vector[i]);
            m_link_genalized_cost_array[i] = pLink->travel_time_per_period[tau] + pLink->route_choice_cost + pLink->toll / m_value_of_time * 60;  // *60 as 60 min per hour
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
                if(p_assignment->g_LinkTypeMap[g_link_vector[link_seq_no].link_type].AllowAgentType (p_assignment->g_AgentTypeVector[m_agent_type_no].agent_type))
                { 
                    m_outgoing_link_seq_no_vector[outgoing_link_size] = link_seq_no;
                    m_to_node_seq_no_vector[outgoing_link_size] = g_node_vector[i].m_to_node_seq_no_vector[j];

                    outgoing_link_size++;

                    if (outgoing_link_size >= _MAX_LINK_SIZE_FOR_A_NODE)
                    {
                        g_fout << " Error: outgoing_link_size >= _MAX_LINK_SIZE_FOR_A_NODE" << endl;
                        g_ProgramStop();
                    }
                }
            }

            int node_seq_no = g_node_vector[i].node_seq_no;
            NodeForwardStarArray[node_seq_no].OutgoingLinkSize = outgoing_link_size;
            
            if(outgoing_link_size>=1)
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
        if (g_debug_level == 2)
        {
            g_fout << "add outgoing link data into dynamic array" << endl;

            for (int i = 0; i < g_node_vector.size(); ++i)
            {
                if (g_node_vector[i].zone_org_id > 0) // for each physical node
                { // we need to make sure we only create two way connectors between nodes and zones 
                    g_fout << "node id= " << g_node_vector[i].node_id << " with zone id " << g_node_vector[i].zone_org_id << "and "
                           << NodeForwardStarArray[i].OutgoingLinkSize << " outgoing links." << endl;

                    for (int j = 0; j < NodeForwardStarArray[i].OutgoingLinkSize; j++)
                    {
                        int link_seq_no = NodeForwardStarArray[i].OutgoingLinkNoArray[j];
                        g_fout << "  outgoing node = " << g_node_vector[g_link_vector[link_seq_no].to_node_seq_no].node_id << endl;
                    }
                }
                else
                {
                    if (g_debug_level == 3)
                    {
                        g_fout << "node id= " << g_node_vector[i].node_id << " with "
                               << NodeForwardStarArray[i].OutgoingLinkSize << " outgoing links." << endl;

                        for (int j = 0; j < NodeForwardStarArray[i].OutgoingLinkSize; ++j)
                        {
                            int link_seq_no = NodeForwardStarArray[i].OutgoingLinkNoArray[j];
                            g_fout << "  outgoing node = " << g_node_vector[g_link_vector[link_seq_no].to_node_seq_no].node_id << endl;
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

    void backtrace_shortest_path_tree(Assignment& assignment, int iteration_number, int o_node_index);

    //major function 2: // time-dependent label correcting algorithm with double queue implementation
    float optimal_label_correcting(int processor_id, Assignment* p_assignment, int iteration_k, int o_node_index, int d_node_no  = -1, bool pure_travel_time_cost = false)
    {	
        // g_debugging_flag = 1;
        int SE_loop_count = 0;

        if (iteration_k == 0)
            BuildNetwork(p_assignment);  // based on agent type and link type
        
        UpdateGeneralizedLinkCost();

        int origin_node = m_origin_node_vector[o_node_index]; // assigned nodes for computing 
        int origin_zone = m_origin_zone_seq_no_vector[o_node_index]; // assigned nodes for computing 
        int agent_type = m_agent_type_no; // assigned nodes for computing 

        if (p_assignment->g_number_of_nodes >= 1000 && origin_zone%97 == 0)
            g_fout << "label correcting for zone " << origin_zone <<  " in processor " << processor_id <<  endl;

        if (g_debug_level  >= 2)
            g_fout << "SP iteration k =  " << iteration_k << ": origin node: " << g_node_vector[origin_node].node_id << endl;

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
        float new_time = 0;
        float new_distance = 0;
        float new_to_node_cost = 0;
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

            if (g_log_path >= 2)
            { 
                g_fout << "SP:scan SE node: " << g_node_vector[from_node].node_id << " with " 
                << NodeForwardStarArray[from_node].OutgoingLinkSize  << " outgoing link(s). "<< endl;
            }
            //scan all outbound nodes of the current node
            
            int pred_link_seq_no = m_link_predecessor[from_node];

            // for each link (i,j) belong A(i)
            for (int i = 0; i < NodeForwardStarArray[from_node].OutgoingLinkSize; ++i)
            {
                to_node = NodeForwardStarArray[from_node].OutgoingNodeNoArray[i];
                link_sqe_no = NodeForwardStarArray[from_node].OutgoingLinkNoArray[i];

                if (g_log_path >= 2)
                    g_fout << "SP:  checking outgoing node " << g_node_vector[to_node].node_id << endl;

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
                            g_fout << "prohibited movement " << movement_string << " will not be used " << endl;
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
                new_to_node_cost = m_node_label_cost[from_node] + m_link_genalized_cost_array[link_sqe_no];

                if (g_log_path)
                {
                    g_fout << "SP:  checking from node " << g_node_vector[from_node].node_id  
                           << "  to node" << g_node_vector[to_node].node_id << " cost = " << new_to_node_cost << endl;
                }

                if (new_to_node_cost < m_node_label_cost[to_node]) // we only compare cost at the downstream node ToID at the new arrival time t
                {
                    if (g_log_path)
                    {
                        g_fout << "SP:  updating node: " << g_node_vector[to_node].node_id << " current cost:" << m_node_label_cost[to_node] 
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

                    if (g_log_path)
                    {
                        g_fout << "SP: add node " << g_node_vector[to_node].node_id << " new cost:" << new_to_node_cost 
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

        if (g_log_path)
        { 
            g_fout << "SPtree at iteration k = " << iteration_k <<  " origin node: " 
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
                    g_fout << "SP node: " << g_node_vector[i].node_id << " label cost " << m_node_label_cost[i] << "time "
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

NetworkForSP g_RoutingNetwork;


float g_read_float(FILE* f)
{
    /*    
        read a floating point number from the current pointer of the file,
        skip all spaces
     */
    char ch, buf[32];
    int i = 0;
    int flag = 1;

    /* returns -1 if end of file is reached */
    while (true)
    {
        ch = getc(f);
        if (ch == EOF || ch == '*' || ch == '$') return -1;
        if (isdigit(ch))
            break;

        if (ch == '-')
            flag = -1;
        else
            flag = 1;
    }
    
    if (ch == EOF) return -1;
    while (isdigit(ch) || ch == '.') {
        buf[i++] = ch;
        ch = fgetc(f);

    }
    buf[i] = 0;

    /* atof function converts a character string (char *) into a doubleing
    pointer equivalent, and if the string is not a floting point number,
    a zero will be return.
    */

    return (float)(atof(buf) * flag);
}

//split the string by "_"
vector<string> split(const string &s, const string &seperator) 
{
    vector<string> result;
    typedef string::size_type string_size;
    string_size i = 0;

    while (i != s.size()) {
        int flag = 0;
        while (i != s.size() && flag == 0) {
            flag = 1;
            for (string_size x = 0; x < seperator.size(); ++x)
                if (s[i] == seperator[x]) {
                    ++i;
                    flag = 0;
                    break;
                }
        }

        flag = 0;
        string_size j = i;
        while (j != s.size() && flag == 0) {
            for (string_size x = 0; x < seperator.size(); ++x)
                if (s[j] == seperator[x]) {
                    flag = 1;
                    break;
                }
            if (flag == 0)
                ++j;
        }
        if (i != j) {
            result.push_back(s.substr(i, j - i));
            i = j;
        }
    }

    return result;
}
//string test_str = "0300:30:120_0600:30:140";
//
//g_global_minute = g_time_parser(test_str);
//
//for (int i = 0; i < g_global_minute.size(); i++)
//{
//	g_fout << "The number of global minutes is: " << g_global_minute[i] << " minutes" << endl;
//}

vector<float> g_time_parser(string str)
{
    vector<float> output_global_minute;

    int string_lenghth = str.length();

    //ASSERT(string_lenghth < 100);

    const char* string_line = str.data(); //string to char*

    int char_length = strlen(string_line);

    char ch, buf_ddhhmm[32] = { 0 }, buf_SS[32] = { 0 }, buf_sss[32] = { 0 };
    char dd1, dd2, hh1, hh2, mm1, mm2, SS1, SS2, sss1, sss2, sss3;
    float ddf1, ddf2, hhf1, hhf2, mmf1, mmf2, SSf1, SSf2, sssf1, sssf2, sssf3;
    float global_minute = 0;
    float dd = 0, hh = 0, mm = 0, SS = 0, sss = 0;
    int i = 0;
    int buffer_i = 0, buffer_k = 0, buffer_j = 0;
    int num_of_colons = 0;

    //DDHHMM:SS:sss or HHMM:SS:sss

    while (i < char_length)
    {
        ch = string_line[i++];

        if (num_of_colons == 0 && ch != '_' && ch != ':') //input to buf_ddhhmm until we meet the colon
        {
            buf_ddhhmm[buffer_i++] = ch;
        }
        else if (num_of_colons == 1 && ch != ':') //start the Second "SS"
        {
            buf_SS[buffer_k++] = ch;
        }
        else if (num_of_colons == 2 && ch != ':') //start the Millisecond "sss"
        {
            buf_sss[buffer_j++] = ch;
        }

        if (ch == '_' || i == char_length) //start a new time string
        {
            if (buffer_i == 4) //"HHMM"
            {
                //HHMM, 0123
                hh1 = buf_ddhhmm[0]; //read each first
                hh2 = buf_ddhhmm[1];
                mm1 = buf_ddhhmm[2];
                mm2 = buf_ddhhmm[3];

                hhf1 = ((float)hh1 - 48); //convert a char to a float
                hhf2 = ((float)hh2 - 48);
                mmf1 = ((float)mm1 - 48);
                mmf2 = ((float)mm2 - 48);

                dd = 0;
                hh = hhf1 * 10 * 60 + hhf2 * 60;
                mm = mmf1 * 10 + mmf2;
            }
            else if (buffer_i == 6) //"DDHHMM"
            {
                //DDHHMM, 012345
                dd1 = buf_ddhhmm[0]; //read each first
                dd2 = buf_ddhhmm[1];
                hh1 = buf_ddhhmm[2];
                hh2 = buf_ddhhmm[3];
                mm1 = buf_ddhhmm[4];
                mm2 = buf_ddhhmm[5];

                ddf1 = ((float)dd1 - 48); //convert a char to a float
                ddf2 = ((float)dd2 - 48);
                hhf1 = ((float)hh1 - 48);
                hhf2 = ((float)hh2 - 48);
                mmf1 = ((float)mm1 - 48);
                mmf2 = ((float)mm2 - 48);

                dd = ddf1 * 10 * 24 * 60 + ddf2 * 24 * 60;
                hh = hhf1 * 10 * 60 + hhf2 * 60;
                mm = mmf1 * 10 + mmf2;
            }

            if (num_of_colons == 1 || num_of_colons == 2)
            {
                //SS, 01
                SS1 = buf_SS[0]; //read each first
                SS2 = buf_SS[1];

                SSf1 = ((float)SS1 - 48); //convert a char to a float
                SSf2 = ((float)SS2 - 48);

                SS = (SSf1 * 10 + SSf2) / 60;
            }

            if (num_of_colons == 2)
            {
                //sss, 012
                sss1 = buf_sss[0]; //read each first
                sss2 = buf_sss[1];
                sss3 = buf_sss[2];

                sssf1 = ((float)sss1 - 48); //convert a char to a float
                sssf2 = ((float)sss2 - 48);
                sssf3 = ((float)sss3 - 48);

                sss = (sssf1 * 100 + sssf2 * 10 + sssf3) / 1000;
            }

            global_minute = dd + hh + mm + SS + sss;

            output_global_minute.push_back(global_minute);

            //initialize the parameters
            buffer_i = 0;
            buffer_k = 0;
            buffer_j = 0;
            num_of_colons = 0;
        }

        if (ch == ':')
            num_of_colons += 1;
    }

    return output_global_minute;
}

//vector<float> g_time_parser(vector<string>& inputstring)
//{
//	vector<float> output_global_minute;
//
//	for (int k = 0; k < inputstring.size(); k++)
//	{
//		vector<string> sub_string = split(inputstring[k], "_");
//
//		for (int i = 0; i < sub_string.size(); i++)
//		{
//			//HHMM
//			//012345
//			char hh1 = sub_string[i].at(0);
//			char hh2 = sub_string[i].at(1);
//			char mm1 = sub_string[i].at(2);
//			char mm2 = sub_string[i].at(3);
//
//			float hhf1 = ((float)hh1 - 48);
//			float hhf2 = ((float)hh2 - 48);
//			float mmf1 = ((float)mm1 - 48);
//			float mmf2 = ((float)mm2 - 48);
//
//			float hh = hhf1 * 10 * 60 + hhf2 * 60;
//			float mm = mmf1 * 10 + mmf2;
//			float global_mm_temp = hh + mm;
//			output_global_minute.push_back(global_mm_temp);
//		}
//	}
//
//	return output_global_minute;
//} // transform hhmm to minutes 

// Peiheng, 02/02/21, inline?
inline string g_time_coding(float time_stamp)
{
    int hour = time_stamp / 60;
    int minute = time_stamp - hour * 60;
    int second = (time_stamp - hour * 60 - minute)*60;

    ostringstream strm;
    strm.fill('0');
    strm << std::setw(2) << hour << std::setw(2) << minute /*<< ":" << setw(2) << second*/;

    return strm.str();
} // transform hhmm to minutes 

//void ReadLinkTollScenarioFile(Assignment& assignment)
//{
//
//	for (unsigned li = 0; li < g_link_vector.size(); li++)
//	{
//
//		g_link_vector[li].TollMAP.erase(g_link_vector[li].TollMAP.begin(), g_link_vector[li].TollMAP.end()); // remove all previouly read records
//	}
//
//	// generate toll based on demand type code in input_link.csv file
//	int demand_flow_type_count = 0;
//
//	for (unsigned li = 0; li < g_link_vector.size(); li++)
//	{
//		if (g_link_vector[li].agent_type_code.size() >= 1)
//		{  // with data string
//
//			std::string agent_type_code = g_link_vector[li].agent_type_code;
//
//			vector<float> TollRate;
//			for (int at = 0; at < assignment.g_AgentTypeVector.size(); at++)
//			{
//				CString number;
//				number.Format(_T("%d"), at);
//
//				std::string str_number = CString2StdString(number);
//				if (agent_type_code.find(str_number) == std::string::npos)   // do not find this number
//				{
//					g_link_vector[li].TollMAP[at] = 999;
//					demand_flow_type_count++;
//				}
//				else
//				{
//					g_link_vector[li].TollMAP[at] = 0;
//				}
//
//			}  //end of pt
//		}
//	}
//}

int g_ParserIntSequence(std::string str, std::vector<int>& vect)
{
    std::stringstream ss(str);
    int i;

    while (ss >> i)
    {
        vect.push_back(i);
        if (ss.peek() == ';')
            ss.ignore();
    }

    return vect.size();
}

void g_ReadDemandFileBasedOnDemandFileList(Assignment& assignment)
{
    //	fprintf(g_pFileOutputLog, "number of zones =,%lu\n", g_zone_vector.size());

    assignment.InitializeDemandMatrix(g_zone_vector.size(), assignment.g_AgentTypeVector.size(), assignment.g_DemandPeriodVector.size());

    float total_demand_in_demand_file = 0;

    CCSVParser parser;
    g_fout << endl;
    g_fout << "Step 1.8: Reading file section [demand_file_list] in setting.csv..." << endl;
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
                    g_fout << "Please provide format_type in section [demand_file_list.]" << endl;
                    g_ProgramStop();
                }

                parser.GetValueByFieldName("agent_type", agent_type);

                int agent_type_no = 0;
                int demand_period_no = 0;

                if (assignment.demand_period_to_seqno_mapping.find(demand_period) != assignment.demand_period_to_seqno_mapping.end())
                    demand_period_no = assignment.demand_period_to_seqno_mapping[demand_period];
                else
                {
                    g_fout << "Error: demand period in section [demand_file_list]" << demand_period << "cannot be found." << endl;
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
                        g_fout << "Error: agent_type in agent_type " << agent_type << "cannot be found." << endl;
                        g_ProgramStop();
                    }
                }

                if (demand_period_no > _MAX_TIMEPERIODS)
                {
                    g_fout << "demand_period_no should be less than settings in demand_period section. Please change the parameter settings in the source code." << endl;
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
                                    g_fout << endl << "Error: Only one line has been read from file. Are there multiple columns of demand type in file " << file_name << " per line?" << endl;
                                    g_ProgramStop();
                                }
                                break;
                            }

                            if (assignment.g_zoneid_to_zone_seq_no_mapping.find(origin_zone) == assignment.g_zoneid_to_zone_seq_no_mapping.end())
                            {
                                if(error_count < 10)
                                    g_fout << endl << "Warning: origin zone " << origin_zone << "  has not been defined in node.csv" << endl;

                                error_count++;
                                 // origin zone has not been defined, skipped. 
                                continue;
                            }

                            if (assignment.g_zoneid_to_zone_seq_no_mapping.find(destination_zone) == assignment.g_zoneid_to_zone_seq_no_mapping.end())
                            {
                                if (error_count < 10)
                                    g_fout << endl << "Warning: destination zone " << destination_zone << "  has not been defined in node.csv" << endl;

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
                                g_fout << "o_zone_id:" << origin_zone << ", d_zone_id: " << destination_zone << ", value = " << demand_value << endl;

                            line_no++;
                        }  // scan lines

                        fclose(st);

                        g_fout << "total_demand_volume is " << assignment.total_demand_volume << endl << endl;
                    }
                    else
                    {
                        // open file
                        g_fout << "Error: File " << file_name << " cannot be opened.\n It might be currently used and locked by EXCEL." << endl;
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
                                g_fout << "demand_volume is " << total_demand_in_demand_file << endl;

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
                                    g_fout << "ODDemandVolume <= 0.001 for OD pair" << o_zone_id << "->" << d_zone_id  << "in routing policy file " << file_name.c_str() << ". Please check" <<  endl;
                                    g_ProgramStop();
                                }

                                agent_path_element.volume = routing_ratio * ODDemandVolume;
                                //assignment.g_origin_demand_array[from_zone_seq_no][agent_type_no][demand_period_no] should be loaded first
                            }

                            //apply for both agent csv and routing policy
                            assignment.g_column_pool[from_zone_seq_no][to_zone_seq_no][agent_type_no][demand_period_no].bfixed_route = true;

                            bool bValid = true;

                            std::string path_node_sequence;
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
                        g_fout << "Error: File " << file_name << " cannot be opened.\n It might be currently used and locked by EXCEL." << endl;
                        g_ProgramStop();
                    }
                }
                else
                {
                    g_fout << "Error: format_type = " << format_type << " is not supported. Currently STALite supports format such as column and agent_csv." << endl;
                    g_ProgramStop();
                }
            }
        }
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

void g_ReadInputData(Assignment& assignment)
{
    assignment.g_LoadingStartTimeInMin = 99999;
    assignment.g_LoadingEndTimeInMin = 0;

    //step 0:read demand period file
    CCSVParser parser_demand_period;
    g_fout << "Step 1: Reading input data" << endl;
    g_fout << "Step 1.1: Reading section [demand_period] in setting.csv..." << endl;
    
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
                    g_fout << "Error: Field demand_period in file demand_period cannot be read." << endl;
                    g_ProgramStop();
                }
                
                vector<float> global_minute_vector;

                if (!parser_demand_period.GetValueByFieldName("time_period", demand_period.time_period))
                { 
                    g_fout << "Error: Field time_period in file demand_period cannot be read." << endl;
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
            g_fout << "Error:  Section demand_period has no information." << endl;
            g_ProgramStop();
        }
    }
    else
    {
        g_fout << "Error: File settings.csv cannot be opened.\n It might be currently used and locked by EXCEL." << endl;
        g_ProgramStop();
    }

    g_fout << "number of demand periods = " << assignment.g_DemandPeriodVector.size() << endl;

    assignment.g_number_of_demand_periods = assignment.g_DemandPeriodVector.size();
    //step 1:read demand type file

    g_fout << "Step 1.2: Reading section [link_type] in setting.csv..." << endl;
    
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
                        g_fout << "Error: Field link_type cannot be found in file link_type.csv." << endl;
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
                    g_fout << "Error: Field link_type " << element.link_type << " has been defined more than once in file link_type.csv." << endl;
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

                g_fout << "important: traffic_flow_code on link type " << element.link_type  << " is " << element.traffic_flow_code  << endl;

                parser_link_type.GetValueByFieldName("agent_type_blocklist", element.agent_type_blocklist, true);

                if (element.agent_type_blocklist.size() >= 1)
                    g_fout << "important: agent type of " << element.agent_type_blocklist << " are prohibited " << " on link type " << element.link_type << endl;
            
                assignment.g_LinkTypeMap[element.link_type] = element;
                line_no++;
            }
        }

        parser_link_type.CloseCSVFile();
    }

    g_fout << "number of link types = " << assignment.g_LinkTypeMap.size() << endl;

    CCSVParser parser_agent_type;
    g_fout << "Step 1.3: Reading section [agent_type] in setting.csv..." << endl;

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
            g_fout << "Error: Section agent_type does not contain information." << endl;
    }

    if (assignment.g_AgentTypeVector.size() >= _MAX_AGNETTYPES)
    {
        g_fout << "Error: agent_type = " << assignment.g_AgentTypeVector.size() << " in section agent_type is too large. " << "_MAX_AGNETTYPES = " << _MAX_AGNETTYPES << "Please contact program developers!";
        g_ProgramStop();
    }

    g_fout << "number of agent typess = " << assignment.g_AgentTypeVector.size() << endl;

    assignment.g_number_of_nodes = 0;
    assignment.g_number_of_links = 0;  // initialize  the counter to 0

    int internal_node_seq_no = 0;
    // step 3: read node file 

    std::map<int, int> zone_id_to_centriod_node_id_mapping;  // this is an one-to-one mapping
    std::map<int, int> zone_id_mapping;  // this is used to mark if this zone_id has been identified or not

    std::map<int, int> zone_id_production; 
    std::map<int, int> zone_id_attraction; 

    CCSVParser parser;

    g_fout << "Step 1.4: Reading node data in node.csv..."<< endl;

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

            parser.GetValueByFieldName("x_coord", node.x,true, false);
            parser.GetValueByFieldName("y_coord", node.y,true, false);
            
            // this is an activity node // we do not allow zone id of zero
            if(zone_id>=1)
            { 
                // for physcial nodes because only centriod can have valid zone_id. 
                node.zone_org_id = zone_id;
                if (zone_id_mapping.find(zone_id) == zone_id_mapping.end())
                {
                    //create zone 
                    zone_id_mapping[zone_id] = node_id;
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

            /*node.x = x;
            node.y = y;*/
            internal_node_seq_no++;

            // push it to the global node vector
            g_node_vector.push_back(node);
            assignment.g_number_of_nodes++;

            if (assignment.g_number_of_nodes % 5000 == 0)
                g_fout << "reading " << assignment.g_number_of_nodes << " nodes.. " << endl;
        }

        g_fout << "number of nodes = " << assignment.g_number_of_nodes << endl;

    	// fprintf(g_pFileOutputLog, "number of nodes =,%d\n", assignment.g_number_of_nodes);
        parser.CloseCSVFile();
    }

    // initialize zone vector
    g_fout << "Step 1.5: Initializing O-D zone vector..." << endl;

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
    
    // for od calibration, I think we don't need to implement for now
    // gravity model.
    if (assignment.assignment_mode == 5)
    {
        g_fout << "writing demand.csv.." << endl;

        FILE* g_pFileODMatrix = nullptr;
        fopen_ss(&g_pFileODMatrix, "demand.csv", "w");

        if (!g_pFileODMatrix)
        {
            g_fout << "File demand.csv cannot be opened." << endl;
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
                            g_fout << "orig= " << g_zone_vector[orig].zone_id << " dest= " << g_zone_vector[dest].zone_id << ":" <<  value << endl;
                        }
                    }
                }
            }

            fclose(g_pFileODMatrix);
        }
    }

    g_fout << "number of zones = " << g_zone_vector.size() << endl;
    // step 4: read link file 

    CCSVParser parser_link;

    g_fout << "Step 1.6: Reading link data in link.csv... " << endl;
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
                g_fout << "Error: from_node_id " << from_node_id << " in file link.csv is not defined in node.csv." << endl;
                continue; //has not been defined
            }

            if (assignment.g_node_id_to_seq_no_map.find(to_node_id) == assignment.g_node_id_to_seq_no_map.end())
            {
                g_fout << "Error: to_node_id " << to_node_id << " in file link.csv is not defined in node.csv." << endl;
                continue; //has not been defined
            }

            if (assignment.g_link_id_map.find(linkID) != assignment.g_link_id_map.end())
                g_fout << "Error: link_id " << linkID.c_str() << " has been defined more than once. Please check link.csv." << endl;

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

            parser_link.GetValueByFieldName("link_type", link.link_type);

            string movement_str;
            parser_link.GetValueByFieldName("movement_str", movement_str, false);
            parser_link.GetValueByFieldName("geometry", link.geometry,false);

            // and valid
            if (movement_str.size() > 0) 
            {
                int main_node_id = -1;
                parser_link.GetValueByFieldName("main_node_id", main_node_id, true, false);

                int NEMA_phase_number = 0;
                parser_link.GetValueByFieldName("NEMA_phase_number", NEMA_phase_number, false, true);

                link.movement_str = movement_str;
                link.main_node_id = main_node_id;
                link.NEMA_phase_number = NEMA_phase_number;
            }

            if (assignment.g_LinkTypeMap.find(link.link_type) == assignment.g_LinkTypeMap.end())
            {
                g_fout << "link type " << link.link_type << " in link.csv is not defined for link " << from_node_id << "->"<< to_node_id << " in link_type.csv" <<endl;
                g_ProgramStop();
            }

            if (assignment.g_LinkTypeMap[link.link_type].type_code == "c" && g_node_vector[internal_from_node_seq_no].zone_id >=0)
            {
                if(assignment.g_zoneid_to_zone_seq_no_mapping.find(g_node_vector[internal_from_node_seq_no].zone_id) != assignment.g_zoneid_to_zone_seq_no_mapping.end())
                link.zone_seq_no_for_outgoing_connector = assignment.g_zoneid_to_zone_seq_no_mapping [g_node_vector[internal_from_node_seq_no].zone_id];
            }

            parser_link.GetValueByFieldName("toll", link.toll,false,false);
            parser_link.GetValueByFieldName("additional_cost", link.route_choice_cost, false, false);

            float length = 1.0; // km or mile
            float free_speed = 1.0;
            float k_jam = 200;
            float bwtt_speed = 12;  //miles per hour

            float lane_capacity = 1800;
            parser_link.GetValueByFieldName("length", length);
            parser_link.GetValueByFieldName("free_speed", free_speed);

            free_speed = max(0.1f, free_speed);

            int number_of_lanes = 1;
            parser_link.GetValueByFieldName("lanes", number_of_lanes);
            parser_link.GetValueByFieldName("capacity", lane_capacity);

            link.free_flow_travel_time_in_min = length / free_speed * 60;
            link.traffic_flow_code = assignment.g_LinkTypeMap[link.link_type].traffic_flow_code;

            //spatial queue and kinematic wave
            if (link.traffic_flow_code >= 2)
                link.spatial_capacity_in_vehicles = max(1.0f,length * number_of_lanes * k_jam);

            // kinematic wave
            if (link.traffic_flow_code == 3)    
                link.BWTT_in_simulation_interval = length / bwtt_speed *3600/ number_of_seconds_per_interval;

            // Peiheng, 02/03/21, useless block
            if (linkID == "10")
                int i_debug = 1;

            char VDF_field_name[20];

            for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
            {
                //setup default values
                link.VDF_period[tau].capacity = lane_capacity * number_of_lanes;
                link.VDF_period[tau].FFTT = length / free_speed * 60.0;  // 60.0 for 60 min per hour
                link.VDF_period[tau].alpha = 0.15;
                link.VDF_period[tau].beta = 4;
                link.VDF_period[tau].starting_time_slot_no = assignment.g_DemandPeriodVector[tau].starting_time_slot_no;
                link.VDF_period[tau].ending_time_slot_no = assignment.g_DemandPeriodVector[tau].ending_time_slot_no;

                int demand_period_id = assignment.g_DemandPeriodVector[tau].demand_period_id;
                sprintf (VDF_field_name, "VDF_fftt%d", demand_period_id);
                parser_link.GetValueByFieldName(VDF_field_name, link.VDF_period[tau].FFTT,false,false);  // FFTT should be per min

                sprintf (VDF_field_name, "VDF_cap%d", demand_period_id);
                parser_link.GetValueByFieldName(VDF_field_name, link.VDF_period[tau].capacity, false, false);  // capacity should be per period per link (include all lanes)

                sprintf (VDF_field_name, "VDF_alpha%d", demand_period_id);
                parser_link.GetValueByFieldName(VDF_field_name, link.VDF_period[tau].alpha, false, false);

                sprintf (VDF_field_name, "VDF_beta%d", demand_period_id);
                parser_link.GetValueByFieldName(VDF_field_name, link.VDF_period[tau].beta, false, false);

                sprintf(VDF_field_name, "VDF_PHF%d", demand_period_id);
                parser_link.GetValueByFieldName(VDF_field_name, link.VDF_period[tau].PHF, false, false); 

                sprintf (VDF_field_name, "VDF_mu%d", demand_period_id);
                parser_link.GetValueByFieldName(VDF_field_name, link.VDF_period[tau].mu, false, false);  // mu should be per hour per link, so that we can calculate congestion duration and D/mu in BPR-X

                //sprintf (VDF_field_name, "VDF_gamma%d", demand_period_id);  // remove gamma
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
            link.link_spatial_capacity = k_jam * number_of_lanes*length;

            link.length = length;
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

            assignment.g_number_of_links++;

            if (assignment.g_number_of_links % 10000 == 0)
                g_fout << "reading " << assignment.g_number_of_links << " links.. " << endl;
        }

        parser_link.CloseCSVFile();
    }
    // we now know the number of links
    g_fout << "number of links = " << assignment.g_number_of_links << endl;

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

    g_fout << "number of links =" << assignment.g_number_of_links << endl;

    if (g_debug_level == 2)
    {
        for (int i = 0; i < g_node_vector.size(); ++i)
        {
            if (g_node_vector[i].zone_org_id > 0) // for each physical node
            { 
                // we need to make sure we only create two way connectors between nodes and zones 
                g_fout << "node id= " << g_node_vector[i].node_id << " with zone id " << g_node_vector[i].zone_org_id << "and "
                    << g_node_vector[i].m_outgoing_link_seq_no_vector.size() << " outgoing links." << endl;
                for (int j = 0; j < g_node_vector[i].m_outgoing_link_seq_no_vector.size(); ++j)
                {
                    int link_seq_no = g_node_vector[i].m_outgoing_link_seq_no_vector[j];
                    g_fout << "  outgoing node = " << g_node_vector[g_link_vector[link_seq_no].to_node_seq_no].node_id << endl;
                }
            }
            else
            {
                if (g_debug_level == 3)
                {
                    g_fout << "node id= " << g_node_vector[i].node_id << " with " << g_node_vector[i].m_outgoing_link_seq_no_vector.size() << " outgoing links." << endl;
                    for (int j = 0; j < g_node_vector[i].m_outgoing_link_seq_no_vector.size(); ++j)
                    {
                        int link_seq_no = g_node_vector[i].m_outgoing_link_seq_no_vector[j];
                        g_fout << "  outgoing node = " << g_node_vector[g_link_vector[link_seq_no].to_node_seq_no].node_id << endl;
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
                g_fout << "Error: node_id " << node_id << " in file movement.csv is not defined in node.csv." << endl;
                //has not been defined
                continue;
            }

            parser_movement.GetValueByFieldName("ib_link_id", ib_link_id);
            parser_movement.GetValueByFieldName("ob_link_id", ob_link_id);

            if (assignment.g_link_id_map.find(ib_link_id) != assignment.g_link_id_map.end())
                g_fout << "Error: ib_link_id " << ib_link_id.c_str() << " has not been defined in movement.csv. Please check link.csv." << endl;

            if (assignment.g_link_id_map.find(ob_link_id) != assignment.g_link_id_map.end())
                g_fout << "Error: ob_link_id " << ob_link_id.c_str() << " has not been defined in movement.csv. Please check link.csv." << endl;

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

        g_fout << "Step XX: Reading movement.csv data with " << prohibited_count << " prohibited records." << endl;
        parser_movement.CloseCSVFile();
    }

    // for testing
    cout << "g_ReadInputData" <<endl;
};

void g_reload_service_arc_data(Assignment& assignment)
{
    g_fout << "Step 1.7: Reading service arc in timing.csv..." << endl;

    CCSVParser parser_service_arc;
    if (parser_service_arc.OpenCSVFile("timing.csv", false))
    {
        while (parser_service_arc.ReadRecord())  // if this line contains [] mark, then we will also read field headers.
        {
            int from_node_id = 0;
            if (!parser_service_arc.GetValueByFieldName("from_node_id", from_node_id))
            {
                g_fout << "Error: from_node_id in file timing.csv is not defined." << endl;
                continue;
            }
            
            int to_node_id = 0;
            if (!parser_service_arc.GetValueByFieldName("to_node_id", to_node_id))
            {
                continue;
            }

            if (assignment.g_node_id_to_seq_no_map.find(from_node_id) == assignment.g_node_id_to_seq_no_map.end())
            {
                g_fout << "Error: from_node_id " << from_node_id << " in file timing.csv is not defined in node.csv." << endl;
                //has not been defined
                continue;
            }
            if (assignment.g_node_id_to_seq_no_map.find(to_node_id) == assignment.g_node_id_to_seq_no_map.end())
            {
                g_fout << "Error: to_node_id " << to_node_id << " in file timing.csv is not defined in node.csv." << endl;
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
                g_fout << "Error: Link " << from_node_id << "->" << to_node_id << " in file timing.csv is not defined in link.csv." << endl;
                continue;
            }

            string time_period;
            if (!parser_service_arc.GetValueByFieldName("time_window", time_period))
            {
                g_fout << "Error: Field time_window in file timing.csv cannot be read." << endl;
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
                g_fout << "reading " << assignment.g_number_of_timing_arcs << " timing_arcs.. " << endl;
        }

        parser_service_arc.CloseCSVFile();
    }

    g_fout << endl;
    g_fout << "Step 1.8: Reading file section [demand_file_list] in setting.csv..." << endl;
    
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
                    g_fout << "Error: from_node_id in file timing.csv is not defined." << endl;
                    continue;
                }

                int to_node_id = 0;
                if (!parser.GetValueByFieldName("to_node_id", to_node_id))
                    continue;

                if (assignment.g_node_id_to_seq_no_map.find(from_node_id) == assignment.g_node_id_to_seq_no_map.end())
                {
                    g_fout << "Error: from_node_id " << from_node_id << " in file timing.csv is not defined in node.csv." << endl;
                    //has not been defined
                    continue;
                }
                if (assignment.g_node_id_to_seq_no_map.find(to_node_id) == assignment.g_node_id_to_seq_no_map.end())
                {
                    g_fout << "Error: to_node_id " << to_node_id << " in file timing.csv is not defined in node.csv." << endl;
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
                    g_fout << "Error: Link " << from_node_id << "->" << to_node_id << " in file timing.csv is not defined in link.csv." << endl;
                    continue;
                }

                string time_period;
                if (!parser.GetValueByFieldName("time_window", time_period))
                {
                    g_fout << "Error: Field time_window in file timing.csv cannot be read." << endl;
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

                g_fout << "reading " << assignment.g_number_of_timing_arcs << " capacity reduction scenario.. " << endl;
            }
        }
    
        parser.CloseCSVFile(); 
    }
    // we now know the number of links
    g_fout << "number of timing records = " << assignment.g_number_of_timing_arcs << endl << endl;
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

    if(iteration_index>=0)
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
                                    PCE_ratio = 1;
                                    //#pragma omp critical
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

void update_link_travel_time_and_cost()
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
}

// changes here are also for odmes, don't need to implement the changes in this function for now
void g_reset_and_update_link_volume_based_on_ODME_columns(int number_of_links, int iteration_no)
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

    // calcualte deviation for each measurement type
    for (int i = 0; i < number_of_links; ++i)
    {
        g_link_vector[i].CalculateTD_VDFunction();

        if (g_link_vector[i].obs_count >= 1)  // with data
        {
            int tau = 0;
            g_link_vector[i].est_count_dev = g_link_vector[i].flow_volume_per_period[tau] - g_link_vector[i].obs_count;
            if (g_debug_level == 2)
            { 
                g_fout << "link " << g_node_vector [g_link_vector[i].from_node_seq_no].node_id
                       << "->" << g_node_vector[g_link_vector[i].to_node_seq_no].node_id
                       << "obs:, " << g_link_vector[i].obs_count << "est:, " << g_link_vector[i].flow_volume_per_period[tau]
                       << "dev:," << g_link_vector[i].est_count_dev << endl;
            }
            total_gap += abs(g_link_vector[i].est_count_dev);
            sub_total_gap_link_count += g_link_vector[i].est_count_dev / g_link_vector[i].obs_count;
        }
    }

    for (int orig = 0; orig < g_zone_vector.size(); ++orig)  // o
    {
        if (g_zone_vector[orig].obs_attraction >= 1)  // with observation
        {
            g_zone_vector[orig].est_attraction_dev = g_zone_vector[orig].est_attraction - g_zone_vector[orig].obs_attraction;

            if (g_debug_level  == 2)
            {
                g_fout << "zone " << g_zone_vector[orig].zone_id << "A: obs:" << g_zone_vector[orig].obs_attraction
                       << ",est:," << g_zone_vector[orig].est_attraction << ",dev:," << g_zone_vector[orig].est_attraction_dev << endl;
            }

            total_gap += abs(g_zone_vector[orig].est_attraction_dev);
            sub_total_gap_A_count += g_zone_vector[orig].est_attraction_dev / g_zone_vector[orig].obs_attraction;
        }

        if (g_zone_vector[orig].obs_production >= 1)  // with observation
        {
            g_zone_vector[orig].est_production_dev = g_zone_vector[orig].est_production - g_zone_vector[orig].obs_production;

            if (g_debug_level  == 2)
            {
                g_fout << "zone " << g_zone_vector[orig].zone_id << "P: obs:" << g_zone_vector[orig].obs_production
                       << ",est:," << g_zone_vector[orig].est_production << ",dev:," << g_zone_vector[orig].est_production_dev << endl;
            }

            total_gap += abs(g_zone_vector[orig].est_production_dev);
            sub_total_gap_P_count += g_zone_vector[orig].est_production_dev / g_zone_vector[orig].obs_production;
        }
    }

    g_fout << "ODME #" << iteration_no<< " total abs gap= " << total_gap 
           << ",subg_link: " << sub_total_gap_link_count*100 
           << ",subg_P: " << sub_total_gap_P_count*100 
           << ",subg_A: " << sub_total_gap_A_count * 100 << endl;
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

        float path_toll = 0;
        float path_gradient_cost = 0;
        float path_distance = 0;
        float path_travel_time = 0;
        int link_seq_no;

        float link_travel_time;
        float total_switched_out_path_volume = 0;

        float step_size = 0;
        float previous_path_volume = 0;

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
                                path_toll += g_link_vector[link_seq_no].toll;
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
                                    it->second.path_volume = max(0.0f, it->second.path_volume - step_size * it->second.path_gradient_cost_relative_difference);

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
        g_fout << "Current iteration number: " << n << endl;
        g_update_gradient_cost_and_assigned_flow_in_column_pool(assignment, n);
        
        if(g_debug_level >=3)
        { 
            for (int i = 0; i < g_link_vector.size(); ++i) 
            {
                g_fout << "link: " << g_node_vector[g_link_vector[i].from_node_seq_no].node_id << "-->"
                       << g_node_vector[g_link_vector[i].to_node_seq_no].node_id << ", "
                       << "flow count:" << g_link_vector[i].flow_volume_per_period[0] << endl;
            }
        }
    }
}

char str_buffer[STRING_LENGTH_PER_LINE];

void g_output_simulation_result(Assignment& assignment)
{
    g_fout << "writing link_performance.csv.." << endl;

    int b_debug_detail_flag = 0;
    FILE* g_pFileLinkMOE = nullptr;

    fopen_ss(&g_pFileLinkMOE,"link_performance.csv", "w");
    if (!g_pFileLinkMOE)
    {
        g_fout << "File link_performance.csv cannot be opened." << endl;
        g_ProgramStop();
    }
    else
    {
        if (assignment.assignment_mode <= 1 || assignment.assignment_mode == 3)  //ODME
        {
            // Option 2: BPR-X function
            fprintf(g_pFileLinkMOE, "link_id,from_node_id,to_node_id,time_period,volume,travel_time,speed,VOC,queue,density,geometry,");

             //ODME
            if (assignment.assignment_mode == 3)
                fprintf(g_pFileLinkMOE, "obs_count,dev,");

            fprintf(g_pFileLinkMOE, "notes\n");

             //Initialization for all nodes
            for (int i = 0; i < g_link_vector.size(); ++i)
            {
                // virtual connectors
                if (g_link_vector[i].link_type == -1)
                    continue;

                for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
                {
                    float speed = g_link_vector[i].length / (max(0.001f, g_link_vector[i].VDF_period[tau].avg_travel_time) / 60.0);
                    fprintf(g_pFileLinkMOE, "%s,%d,%d,%s,%.3f,%.3f,%.3f,%.3f,0,0,\"%s\",",
                        g_link_vector[i].link_id.c_str(),
                        g_node_vector[g_link_vector[i].from_node_seq_no].node_id,
                        g_node_vector[g_link_vector[i].to_node_seq_no].node_id,
                        assignment.g_DemandPeriodVector[tau].time_period.c_str(),
                        g_link_vector[i].flow_volume_per_period[tau],
                        g_link_vector[i].VDF_period[tau].avg_travel_time,
                        speed,  /* 60.0 is used to convert min to hour */
                        g_link_vector[i].VDF_period[tau].VOC,
                        g_link_vector[i].geometry.c_str());

                    if (assignment.assignment_mode == 3)  //ODME
                    {
                        if (g_link_vector[i].obs_count >= 1) //ODME
                            fprintf(g_pFileLinkMOE, "%.1f,%.1f,", g_link_vector[i].obs_count, g_link_vector[i].est_count_dev);
                        else
                            fprintf(g_pFileLinkMOE, ",,");
                    }
                    fprintf(g_pFileLinkMOE, "period-based\n");

                    // print out for BPR-X 
                    bool b_print_out_for_BPR_X = false;
                    if(b_print_out_for_BPR_X)
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
        }
        else if (assignment.assignment_mode == 2)  // space time based simulation // ODME
        {
            // Option 2: BPR-X function
            fprintf(g_pFileLinkMOE, "link_id,from_node_id,to_node_id,time_period,volume,CA,CD,density,queue,travel_time,waiting_time_in_sec,speed,");
            fprintf(g_pFileLinkMOE, "notes\n");

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

                        fprintf(g_pFileLinkMOE, "%s,%d,%d,%s_%s,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,",
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

                        fprintf(g_pFileLinkMOE, "simulation-based\n");
                    }
                }  // for each time t
            }  // for each link l
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

    if (assignment.assignment_mode == 0)
    {
        FILE* g_pFileODMOE = nullptr;
        fopen_ss(&g_pFileODMOE, "agent.csv", "w");
        fclose(g_pFileODMOE);
    }
    else if(assignment.assignment_mode >= 1)
    {
        g_fout << "writing agent.csv.." << endl;

        float path_time_vector[_MAX_LINK_SIZE_IN_A_PATH];
        FILE* g_pFileODMOE = nullptr;
        fopen_ss(&g_pFileODMOE,"agent.csv", "w");

        if (!g_pFileODMOE)
        {
            g_fout << "File agent.csv cannot be opened." << endl;
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

        g_fout << "writing data for " << zone_size << "  zones " << endl;

        for (int orig = 0; orig < zone_size; ++orig)
        { 
            if(g_zone_vector[orig].zone_id %100 ==0)
                g_fout << "o zone id =  " << g_zone_vector[orig].zone_id << endl;

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
                                    g_fout << "writing " << count/1000 << "K agents with CPU time " << iteration_t / 1000.0 << " s" << endl;
                                }

                                path_toll = 0;
                                path_distance = 0;
                                path_travel_time = 0;
                                path_time_vector[0] = time_stamp;

                                for (int nl = 0; nl < it->second.m_link_size; ++nl)  // arc a
                                {
                                    int link_seq_no = it->second.path_link_vector[nl];
                                    path_toll += g_link_vector[link_seq_no].toll;
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
                                                         path_distance
                                                        );

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
                                        g_fout << "Error: buffer_len >= STRING_LENGTH_PER_LINE." << endl;
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
                                    fprintf(g_pFileODMOE, "%s", str_buffer, buffer_len);
                                    count++;
                                }
                                else
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
                                                             path_distance
                                        );

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
                                                time_in_min = assignment.g_LoadingStartTimeInMin + pAgentSimu->m_Veh_LinkArrivalTime_in_simu_interval[nt] * number_of_seconds_per_interval / 60.0 ;
                                            else
                                                time_in_min = assignment.g_LoadingStartTimeInMin + pAgentSimu->m_Veh_LinkDepartureTime_in_simu_interval[nt - 1]* number_of_seconds_per_interval / 60.0 ;  // last link in the path

                                            path_time_vector[nt] = time_in_min;
                                        }

                                        for (int nt = 0 + virtual_link_delta; nt < it->second.m_link_size + 1 - virtual_link_delta; ++nt)
                                            buffer_len += sprintf(str_buffer + buffer_len, "%s;", g_time_coding(path_time_vector[nt]).c_str());
                                        buffer_len += sprintf(str_buffer + buffer_len, ",");

                                        for (int nt = 0 + virtual_link_delta; nt < it->second.m_link_size + 1 - virtual_link_delta; ++nt)
                                            buffer_len += sprintf(str_buffer + buffer_len, "%.2f;", path_time_vector[nt]);
                                        buffer_len += sprintf(str_buffer + buffer_len, ",");

                                        for (int nt = 0 + virtual_link_delta; nt < it->second.m_link_size - virtual_link_delta; ++nt)
                                            buffer_len += sprintf(str_buffer + buffer_len, "%.2f;", path_time_vector[nt+1] - path_time_vector[nt]);
                                        buffer_len += sprintf(str_buffer + buffer_len, "\n");

                                        if (buffer_len >= STRING_LENGTH_PER_LINE - 1)
                                        {
                                            g_fout << "Error: buffer_len >= STRING_LENGTH_PER_LINE." << endl;
                                            g_ProgramStop();
                                        }

                                        fprintf(g_pFileODMOE, "%s", str_buffer, buffer_len);
                                        count++;
                                    }
                                }
                            }   
                        }    
                    }
                }
            }
        }

        fclose(g_pFileODMOE);
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

    g_fout << "writing link_performance_sig.csv.." << endl;

    int b_debug_detail_flag = 0;
    FILE* g_pFileLinkMOE = nullptr;
    
    fopen_ss(&g_pFileLinkMOE, "link_performance_sig.csv", "w");
    if (!g_pFileLinkMOE)
    {
        g_fout << "File link_performance_sig.csv cannot be opened." << endl;
        g_ProgramStop();
    }

    if (assignment.assignment_mode <= 1)
    {
        // Option 2: BPR-X function
        fprintf(g_pFileLinkMOE, "link_id,from_node_id,to_node_id,demand_period,time_period,movement_str,main_node_id,NEMA_phase_number,volume,travel_time,speed,VOC,");

        fprintf(g_pFileLinkMOE, "notes\n");

        //Initialization for all nodes
        for (int i = 0; i < g_link_vector.size(); ++i)
        {
            for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
            {
                if (g_link_vector[i].movement_str.length() >= 1)
                {
                    float speed = g_link_vector[i].length / (max(0.001f, g_link_vector[i].VDF_period[tau].avg_travel_time) / 60.0);
                    fprintf(g_pFileLinkMOE, "%s,%d,%d,%s,%s,%s,%d,%d,%.3f,%.3f,%.3f,%.3f,",
                        g_link_vector[i].link_id.c_str(),
                        g_node_vector[g_link_vector[i].from_node_seq_no].node_id,
                        g_node_vector[g_link_vector[i].to_node_seq_no].node_id,
                        assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
                        assignment.g_DemandPeriodVector[tau].time_period.c_str(),
                        g_link_vector[i].movement_str.c_str(),
                        g_link_vector[i].main_node_id,
                        g_link_vector[i].NEMA_phase_number,
                        g_link_vector[i].flow_volume_per_period[tau],
                        g_link_vector[i].VDF_period[tau].avg_travel_time,
                        speed,  /* 60.0 is used to convert min to hour */
                        g_link_vector[i].VDF_period[tau].VOC);
                    
                    fprintf(g_pFileLinkMOE, "period-based\n");
                }
            }
        }
    }

    if (assignment.assignment_mode == 2)  // space time based simulation
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
    g_fout << "Step 2: Assigning computing tasks to memory blocks..." << endl;

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
                    p_NetworkForSP->tau = tau;
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

    g_fout << "There are " << g_NetworkForSP_vector.size() << " networks in memory." << endl;
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

void g_reset_link_volume_for_all_processors()
{
#pragma omp parallel for
    for (int ProcessID = 0; ProcessID < g_NetworkForSP_vector.size(); ++ProcessID)
    {
        NetworkForSP* pNetwork = g_NetworkForSP_vector[ProcessID];
        //Initialization for all non-origin nodes
        int number_of_links = assignment.g_number_of_links;
        for (int i = 0; i < number_of_links; ++i)
            pNetwork->m_link_flow_volume_array[i] = 0;
    }
}

void g_fetch_link_volume_for_all_processors()
{
    for (int ProcessID = 0; ProcessID < g_NetworkForSP_vector.size(); ++ProcessID)
    {
        NetworkForSP* pNetwork = g_NetworkForSP_vector[ProcessID];

        for (int i = 0; i < g_link_vector.size(); ++i)
            g_link_vector[i].flow_volume_per_period[pNetwork->tau] += pNetwork->m_link_flow_volume_array[i];
    }
    // step 1: travel time based on VDF
}

//major function: update the cost for each node at each SP tree, using a stack from the origin structure 
void NetworkForSP::backtrace_shortest_path_tree(Assignment& assignment, int iteration_number_outterloop, int o_node_index)
{

    int origin_node = m_origin_node_vector[o_node_index]; // assigned no
    int m_origin_zone_seq_no = m_origin_zone_seq_no_vector[o_node_index]; // assigned no

    //if (assignment.g_number_of_nodes >= 100000 && m_origin_zone_seq_no % 100 == 0)
    //{
    //	g_fout << "backtracing for zone " << m_origin_zone_seq_no << endl;
    //}

    int departure_time = tau;
    int agent_type = m_agent_type_no;

    if (g_node_vector[origin_node].m_outgoing_link_seq_no_vector.size() == 0)
        return;

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
    float ODvolume, volume;
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

            pColumnVector = &(assignment.g_column_pool[m_origin_zone_seq_no][destination_zone_seq_no][agent_type][tau]);

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

                // backtrace the sp tree from the destination to the root (at origin) 
                while (current_node_seq_no >= 0 && current_node_seq_no < number_of_nodes)
                {
                    temp_path_node_vector[l_node_size++] = current_node_seq_no;

                    if (l_node_size >= temp_path_node_vector_size)
                    {
                        g_fout << "Error: l_node_size >= temp_path_node_vector_size" << endl;
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
                    if (pColumnVector->path_node_sequence_map.find(node_sum) == assignment.g_column_pool[m_origin_zone_seq_no][destination_zone_seq_no][agent_type][tau].path_node_sequence_map.end())
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

double network_assignment(int assignment_mode, int iteration_number, int column_updating_iterations)
{
    // for testing
    cout << "network assignment" << endl;

    int signal_updating_iterations = 0;

    // k iterations for column generation
    assignment.g_number_of_column_generation_iterations = iteration_number;
    // 0: link UE: 1: path UE, 2: Path SO, 3: path resource constraints 
    assignment.assignment_mode = assignment_mode; 
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
    g_fout << endl;
    g_fout << "Step 3: Column Generation for Traffic Assignment..." << endl;
    g_fout << "Total Column Generation iteration: " << assignment.g_number_of_column_generation_iterations << endl;
    for (int iteration_number = 0; iteration_number < assignment.g_number_of_column_generation_iterations; iteration_number++)
    {
        g_fout << endl;
        g_fout << "Current iteration number:" << iteration_number << endl;
        end_t = clock();
        iteration_t = end_t - start_t;
        g_fout << "Current CPU time: " << iteration_t / 1000.0 << " s" << endl;

        // commment out for DLL version
        // if (signal_updating_iterations >=1 && iteration_number >= signal_updating_iterations)
        // {
        //     g_fout << "use SignalAPI to recalibrate signal timing at iteration " << iteration_number << endl;
        //     SignalAPI(iteration_number, assignment_mode, 0);
        //     g_reload_service_arc_data(assignment);
        // }

        // step 3.1 update travel time and resource consumption		
        clock_t start_t_lu = clock();

        // initialization at beginning of shortest path
        update_link_travel_time_and_cost(); 

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

        if (g_debug_level  >= 3)
        {
            g_fout << "Results:" << endl;
            for (int i = 0; i < g_link_vector.size(); ++i) {
                g_fout << "link: " << g_node_vector[g_link_vector[i].from_node_seq_no].node_id << "-->"
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
        clock_t  start_t_cp = clock();

        cumulative_lc = 0;
        cumulative_cp = 0;
        cumulative_lu = 0;

#pragma omp parallel for  // step 3: C++ open mp automatically create n threads., each thread has its own computing thread on a cpu core 
        for (int ProcessID = 0; ProcessID < g_NetworkForSP_vector.size(); ++ProcessID)
        {
            int agent_type_no = g_NetworkForSP_vector[ProcessID]->m_agent_type_no;

            for (int o_node_index = 0; o_node_index < g_NetworkForSP_vector[ProcessID]->m_origin_node_vector.size(); ++o_node_index)
            {
                start_t_lc = clock();
                g_NetworkForSP_vector[ProcessID]->optimal_label_correcting(ProcessID, &assignment, iteration_number, o_node_index);
                end_t = clock();
                cumulative_lc += end_t - start_t_lc;

                start_t_cp = clock();
                g_NetworkForSP_vector[ProcessID]->backtrace_shortest_path_tree(assignment, iteration_number, o_node_index);
                end_t = clock();
                cumulative_cp += end_t - start_t_cp;
            }
            // perform one to all shortest path tree calculation
        }

        // link based computing mode, we have to collect link volume from all processors.
        if (assignment.assignment_mode == 0)
            g_fetch_link_volume_for_all_processors();

        // g_fout << "LC with CPU time " << cumulative_lc / 1000.0 << " s; " << endl;
        // g_fout << "column generation with CPU time " << cumulative_cp / 1000.0 << " s; " << endl;

        //****************************************//

        // last iteraion before performing signal timing updating
        if (signal_updating_iterations >= 1 && iteration_number >= signal_updating_iterations-1)
            g_output_simulation_result_for_signal_api(assignment);
    }
    g_fout << endl;

    // step 1.8: column updating stage: for given column pool, update volume assigned for each column
    g_fout << "Step 4: Column Pool Updating" << endl;
    g_fout << "Total Column Pool Updating iteration: " << column_updating_iterations << endl;
    start_t = clock();
    g_column_pool_optimization(assignment, column_updating_iterations);
    g_fout << endl;

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
        g_fout << "Step 5: Simulation for traffic assignment.." << endl;
        assignment.STTrafficSimulation();
        g_fout << endl;
    }

    if (assignment.assignment_mode == 3)
    {
        g_fout << "Step 6: O-D estimation for traffic assignment.." << endl;
        assignment.Demand_ODME(column_updating_iterations);
        g_fout << endl;
    }

    end_t = clock();
    total_t = (end_t - start_t);
    g_fout << "Done!" << endl;

    g_fout << "CPU Running Time for column pool updating: " << total_t / 1000.0 << " s" << endl;

    start_t = clock();

    //step 5: output simulation results of the new demand 
    g_output_simulation_result(assignment);

    end_t = clock();
    total_t = (end_t - start_t);
    g_fout << "Output for assignment with " << assignment.g_number_of_column_generation_iterations << " iterations. Traffic assignment completes!" << endl;
    g_fout << "CPU Running Time for outputting simulation results: " << total_t / 1000.0 << " s" << endl;

    g_fout << "free memory.." << endl;
    g_node_vector.clear();

    for (int i = 0; i < g_link_vector.size(); ++i)
        g_link_vector[i].free_memory();
    g_link_vector.clear();

    g_fout << "done." << endl;

    return 1;
}

void Assignment::AllocateLinkMemory4Simulation()
{
    g_number_of_simulation_intervals = (g_LoadingEndTimeInMin - g_LoadingStartTimeInMin + 60) * 60 /number_of_seconds_per_interval; 
    g_number_of_loading_intervals = (g_LoadingEndTimeInMin - g_LoadingStartTimeInMin) * 60 / number_of_seconds_per_interval;

    g_number_of_simulation_horizon_in_min = (int)(g_number_of_simulation_intervals / number_of_interval_per_min +1);
    // add + 120 as a buffer
    g_number_of_in_memory_simulation_intervals = g_number_of_simulation_intervals;

    g_fout << "allocate 2D dynamic memory LinkOutFlowCapacity..." << endl;

    m_LinkOutFlowCapacity = AllocateDynamicArray <float>(g_number_of_links, g_number_of_simulation_intervals);  //1
    // discharge rate per simulation time interval
    g_fout << "allocate 2D dynamic memory m_LinkCumulativeArrival..." << endl;
    m_LinkCumulativeArrival = AllocateDynamicArray <float>(g_number_of_links, g_number_of_simulation_intervals);  //2

    g_fout << "allocate 2D dynamic memory m_LinkCumulativeDeparture..." << endl;
    m_LinkCumulativeDeparture = AllocateDynamicArray <float>(g_number_of_links, g_number_of_simulation_intervals);  //3

    g_fout << "allocate 2D dynamic memory m_LinkTDTravelTime..." << endl;
    m_LinkTDTravelTime = AllocateDynamicArray <int>(g_number_of_links, g_number_of_simulation_horizon_in_min); //4

    g_fout << "allocate 2D dynamic memory m_LinkTDWaitingTime..." << endl;
    m_LinkTDWaitingTime = AllocateDynamicArray <float>(g_number_of_links, g_number_of_simulation_horizon_in_min); //5

    g_fout << "initializing time dependent capacity data..." << endl;

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

    g_fout << "End of initializing time dependent capacity data." << endl;
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
            g_fout << "generating " << g_agent_simu_vector.size()/1000 << " K agents for "  << orig << "  zones " << endl;

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

                                for (int nl = 0; nl < it->second.m_link_size; nl++)  // arc a
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

    g_fout << "number of simulation zones:" << zone_size << endl;

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
            g_fout << "simu time= " << t / number_of_simu_interval_per_min << " min, CA = " << TotalCumulative_Arrival_Count << " CD=" << TotalCumulative_Departure_Count << endl;
        
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
                    g_fout << "Error: from_node_id " << from_node_id << " in file measurement.csv is not defined in node.csv." << endl;
                    //has not been defined
                    continue;
                }
                if (g_node_id_to_seq_no_map.find(to_node_id) == assignment.g_node_id_to_seq_no_map.end())
                {
                    g_fout << "Error: to_node_id " << to_node_id << " in file measurement.csv is not defined in node.csv." << endl;
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
                    g_fout << "Error: Link " << from_node_id << "->" << to_node_id << " in file timing.csv is not defined in link.csv." << endl;
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

    for (int s = 0; s < OD_updating_iterations; ++s)
    {
        float total_gap = 0;
        float total_relative_gap = 0;
        float total_gap_count = 0;
        //step 2.1
        // we can have a recursive formulat to reupdate the current link volume by a factor of k/(k+1),
        // and use the newly generated path flow to add the additional 1/(k+1)
        g_reset_and_update_link_volume_based_on_ODME_columns(g_link_vector.size(),s); 
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
                            for (it = it_begin; it != it_end; ++it, i++) // for each k 
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

                                    if (g_zone_vector[orig].obs_attraction_upper_bound_flag == 1 && g_zone_vector[dest].est_attraction_dev >0)
                                        path_gradient_cost += g_zone_vector[dest].est_attraction_dev;
                                }

                                float est_count_dev = 0;
                                for (int nl = 0; nl < it->second.m_link_size; ++nl)  // arc a
                                {
                                    // step 3.3 link flow gradient 
                                    link_seq_no = it->second.path_link_vector[nl];
                                    if(g_link_vector[link_seq_no].obs_count>=1)
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

                                it->second.path_volume = max(1.0f, it->second.path_volume - change);

                                if (g_log_odme == 1)
                                { 
                                    g_fout << "OD " << orig << "-> " << dest << " path id:" << i << ", prev_vol"
                                           << prev_path_volume << ", gradient_cost = " << it->second.path_gradient_cost 
                                           << " link," << g_link_vector[link_seq_no].est_count_dev
                                           << " P," << g_zone_vector[orig].est_production_dev
                                           << " A," << g_zone_vector[orig].est_attraction_dev
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