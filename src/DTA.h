#ifndef GUARD_DTA_H
#define GUARD_DTA_H

#include <algorithm>
#include <cmath>
#include <iomanip>

// Peiheng, 02/21/22, a temporary fix (bad practice)
using std::max;

constexpr auto MAX_LABEL_COST = 1.0e+15;
constexpr auto _INFO_ZONE_ID = 100000;
constexpr auto MAX_MODETYPES = 10; //10 //because of the od demand store format,the MAX_demandtype must >=g_DEMANDTYPES.size()+1;

constexpr auto MAX_ORIGIN_DISTRICTS = 30; //origin based agreegration grids

constexpr auto MAX_MEMORY_BLOCKS = 100;

constexpr auto MAX_LINK_SIZE_IN_A_PATH = 10000;
constexpr auto MAX_LINK_SIZE_FOR_A_NODE = 10000;
constexpr auto MAX_TIMESLOT_PerPeriod = 300; // max 96 5-min slots per day
constexpr auto MAX_TIMEINTERVAL_PerDay = 300; // max 96*3 5-min slots per day
constexpr auto MAX_DAY_PerYear = 360; // max 96*3 5-min slots per day
constexpr auto _default_saturation_flow_rate = 1800;

constexpr auto MIN_PER_TIMESLOT = 5;
constexpr auto simulation_discharge_period_in_min = 10;


/* make sure we change the following two parameters together*/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

constexpr auto number_of_seconds_per_interval = 0.25;  // consistent with the cell link_distance_VDF of 7 meters
constexpr auto number_of_simu_interval_reaction_time = 4;  // reaction time as 1 second, 4 simu intervals, CAV: 0.5 seconds
constexpr auto number_of_simu_intervals_in_min = 240; // 60/0.25 number_of_seconds_per_interval

//constexpr auto number_of_seconds_per_interval = 0.05;  // consistent with the cell link_distance_VDF of 7 meters
//constexpr auto number_of_simu_interval_reaction_time = 20;  // reaction time as 1 second, 4 simu intervals, CAV: 0.5 seconds
//constexpr auto number_of_simu_intervals_in_min = 1200; // 60/0.25 number_of_seconds_per_interval


/* number_of_seconds_per_interval should satisify the ratio of 60/number_of_seconds_per_interval is an integer*/

// Linear congruential generator
constexpr auto LCG_a = 17364;
constexpr auto LCG_c = 0;
constexpr auto LCG_M = 65521;  // it should be 2^32, but we use a small 16-bit number to save memory

enum e_traffic_flow_model { point_queue = 0, spatial_queue, kinemative_wave };
enum e_assignment_mode { lue = 0, path_based_assignment= 1, simulation_dta=2};

// FILE* g_pFileOutputLog = nullptr;
extern void g_OutputModelFiles(int mode);
extern int g_related_zone_vector_size;

extern int g_number_of_active_mode_types;

extern  std::ofstream  g_DTA_log_file;

extern double g_get_random_ratio();
class CDeparture_time_Profile {
public:
    CDeparture_time_Profile() : departure_time_profile_no{ 0 }, m_RandomSeed { 101 }, cumulative_departure_time_ratio{ 0 }, starting_time_slot_no{ 0 }, ending_time_slot_no{ 0 }
    {
        for (int s = 0; s <= 96 * 3; s++)
        {
            cumulative_departure_time_ratio[s] = 0;
            departure_time_ratio[s] = 1.0/300;
        }
    }

    void compute_cumulative_profile(int starting_slot_no, int ending_slot_no, bool b_with_log = false)
    {
        for (int s = 0; s <= 96 * 3; s++)
        {
            cumulative_departure_time_ratio[s] = 0;
        }

        double total_ratio = 0.0;
        for (int s = starting_slot_no +1; s <=ending_slot_no; s++)
        {
            total_ratio += departure_time_ratio[s];
        }

        if (total_ratio < 0.000001)
            total_ratio = 0.000001;

        cumulative_departure_time_ratio[starting_slot_no] = 0;
        float cumulative_ratio = 0;
        for (int s = starting_slot_no+1; s <= ending_slot_no; s++)
        {
            cumulative_ratio += departure_time_ratio[s] / total_ratio;
            cumulative_departure_time_ratio[s] = cumulative_ratio;

            int hour = s / 12;
            int minute = s * 5 - hour * 60;

            if(s== starting_slot_no + 1 || s== ending_slot_no)
            {
                if(b_with_log)
                {
            dtalog.output() << std::setprecision(5) << "[DATA INFO] Cumulative profile no." << departure_time_profile_no << ", ratio at slot  " << s << " (" << hour << ":" << minute << ") = " <<
                departure_time_ratio[s] << '\n';
            g_DTA_log_file << std::setprecision(5) << "[DATA INFO] Cumulative profile no." << departure_time_profile_no << ", ratio at slot  " << s << " (" << hour << ":" << minute << ") = " <<
                departure_time_ratio[s] << '\n';
                }
            }
        }

        if (b_with_log)
        {
            dtalog.output() << std::setprecision(5) << "[DATA INFO] Final cumulative profile ratio = " << cumulative_departure_time_ratio[ending_slot_no - 1] << '\n';
            g_DTA_log_file << std::setprecision(5) << "[DATA INFO] Final cumulative profile ratio = " << cumulative_departure_time_ratio[ending_slot_no - 1] << '\n';
        }
    }

    int get_time_slot_no(int agent_seq_no, int agent_size)
    {

        float r = 0;
        if (agent_size >= 10)  // large number of agents, then use pure uniform sequence
            r = agent_seq_no * 1.0 / agent_size; // r is between 0 and 1
        else
            r = g_get_random_ratio();  // small sample case

        for (int s = starting_time_slot_no; s < ending_time_slot_no; s++)
        {
            if (r < cumulative_departure_time_ratio[s])
            {
                int hour = s / 12;
                int minute = s * 5 - hour * 60;
//                dtalog.output() << "s=" << s <<" (" << hour << ":" << minute << ") = "  << ending_time_slot_no << '\n';
//                g_DTA_log_file << "s=" << s <<" (" << hour << ":" << minute << ") = "  << ending_time_slot_no << '\n';

                return s;
            }
        }
        int hour = starting_time_slot_no / 12;
        int minute = starting_time_slot_no * 5 - hour * 60;
//        dtalog.output() << "s=" << starting_time_slot_no << " (" << hour << ":" << minute << ") = " << ending_time_slot_no << '\n';
//        g_DTA_log_file << "s=" << starting_time_slot_no << " (" << hour << ":" << minute << ") = " << ending_time_slot_no << '\n';
        return starting_time_slot_no;  // first time slot as the default value
    }

    double get_deparure_time_in_min(int agent_seq_no, int agent_size)
    {
        int idebug = 0;
        double r = 0;
        if (agent_size >= 10)  // large number of agents, then use pure uniform sequence
            r = agent_seq_no * 1.0 / agent_size; // r is between 0 and 1
        else
            r = g_get_random_ratio();  // small sample case

        for (int s = starting_time_slot_no+1; s <= ending_time_slot_no; s++)
        {
            int hour = s / 12;
            int minute = (int)(( s*1.0 /12 - hour) * 60 + 0.5);
            if(idebug)
            {
            dtalog.output() << "[DATA INFO] s=" << s << " (" << hour << ":" << minute << ") = " << cumulative_departure_time_ratio[s] << '\n';
            g_DTA_log_file << "[DATA INFO] s=" << s << " (" << hour << ":" << minute << ") = " << cumulative_departure_time_ratio[s] << '\n';
            }
            if (r < cumulative_departure_time_ratio[s])
            {

                double slot_fraction = cumulative_departure_time_ratio[s] - cumulative_departure_time_ratio[s-1];
                double floating_point = max(0.0, (r - cumulative_departure_time_ratio[s - 1]) / max(0.00001, slot_fraction));

                double time_in_min = (s- starting_time_slot_no + floating_point )* MIN_PER_TIMESLOT;
                if (idebug)
                {
                    dtalog.output() << "[DATA INFO]  select: s=" << s << " (" << hour << ":" << minute << ") = " << ending_time_slot_no << ", dep_time = " << time_in_min <<"," << '\n';
                    g_DTA_log_file << "[DATA INFO]  select: s=" << s << " (" << hour << ":" << minute << ") = " << ending_time_slot_no << ", dep_time = " << time_in_min <<"," << '\n';
                }
                return time_in_min;
            }
        }

        if (idebug)
        {
            int hour = starting_time_slot_no / 12;
            int minute = starting_time_slot_no * 5 - hour * 60;

            dtalog.output() << "[DATA INFO] s=" << starting_time_slot_no << " (" << hour << ":" << minute << ") = " << ending_time_slot_no << '\n';
            g_DTA_log_file << "[DATA INFO] s=" << starting_time_slot_no << " (" << hour << ":" << minute << ") = " << ending_time_slot_no << '\n';
        }
        return (r) * MIN_PER_TIMESLOT  ;  // first time slot as the default value
    }

    unsigned int m_RandomSeed;
    int departure_time_profile_no;
    int starting_time_slot_no;
    int ending_time_slot_no;
    float departure_time_ratio[MAX_TIMESLOT_PerPeriod];
    float cumulative_departure_time_ratio[MAX_TIMESLOT_PerPeriod];
};

class CModeType_Summary
{
public:
    CModeType_Summary() : count{ 0 },
        total_od_volume{ 0 }, total_agent_distance_km{ 0 }, total_agent_distance_mile{ 0 }, total_agent_travel_time{ 0 }, total_agent_co2{ 0 }, total_agent_nox{ 0 }, avg_travel_time{ 0 }, avg_travel_distance_km{ 0 }, avg_travel_distance_mile{ 0 }, avg_co2{ 0 }, avg_nox { 0}
    {}

    int count;
    double total_od_volume;
    double total_agent_distance_km;
    double total_agent_distance_mile;
    double total_agent_travel_time;
    double total_agent_delay;
    double total_agent_co2;
    double total_agent_nox;

    double avg_travel_time;
    double avg_travel_delay;
    double avg_co2;
    double avg_nox;
    double avg_travel_distance_km;
    double avg_travel_distance_mile;
};

class CSystem_Summary
{
public:
    int district_id;
    int district_name;

    std::vector<DTAGDPoint> shape_points;
    std::vector<int> district_zone_vector;


    void reset_data(int at)
    {
        if (at >= MAX_MODETYPES)
            return;

        data_by_mode_type[at].count =0;
        data_by_mode_type[at].total_od_volume = 0;
        data_by_mode_type[at].total_agent_travel_time = 0;
        data_by_mode_type[at].total_agent_distance_km = 0;
        data_by_mode_type[at].total_agent_distance_mile = 0;
        data_by_mode_type[at].total_agent_co2 = 0;
        data_by_mode_type[at].total_agent_nox = 0;

    }

    void record_mode_volume(int at, double od_volume)
    {
        if (at >= MAX_MODETYPES)
            return;

        data_by_mode_type[at].total_od_volume += od_volume;

    }

    void record_mode_od_data(CModeType_Summary element, int at)
    {
        if (at >= MAX_MODETYPES)
            return;

        data_by_mode_type[at].count += 1;
        data_by_mode_type[at].total_agent_travel_time += element.total_agent_travel_time;
        data_by_mode_type[at].total_agent_distance_km += element.total_agent_distance_km;
        data_by_mode_type[at].total_agent_distance_mile += element.total_agent_distance_mile;
        data_by_mode_type[at].total_agent_delay += element.total_agent_delay;

        data_by_mode_type[at].total_agent_co2+= element.total_agent_co2;
        data_by_mode_type[at].total_agent_nox += element.total_agent_nox;

    }

    void computer_avg_value( int at)
    {
        float count = data_by_mode_type[at].count;
        if (count >= 1)
        {
            data_by_mode_type[at].avg_travel_distance_km = data_by_mode_type[at].total_agent_distance_km / max(0.001, data_by_mode_type[at].total_od_volume);
            data_by_mode_type[at].avg_travel_distance_mile = data_by_mode_type[at].total_agent_distance_mile / max(0.001, data_by_mode_type[at].total_od_volume);
            data_by_mode_type[at].avg_travel_time = data_by_mode_type[at].total_agent_travel_time / max(0.001, data_by_mode_type[at].total_od_volume);
            data_by_mode_type[at].avg_travel_delay = data_by_mode_type[at].total_agent_delay / max(0.001, data_by_mode_type[at].total_od_volume);
            data_by_mode_type[at].avg_co2 = data_by_mode_type[at].total_agent_co2 / max(0.001, data_by_mode_type[at].total_od_volume);
            data_by_mode_type[at].avg_nox = data_by_mode_type[at].total_agent_nox / max(0.001, data_by_mode_type[at].total_od_volume);
        }
    }

    CModeType_Summary data_by_mode_type[MAX_MODETYPES];
};

class Cmode_type {
public:
    Cmode_type() : mode_type_no{ 1 }, value_of_time{ 100 }, time_headway_in_sec{ 1 }, real_time_information_type{ 0 }, access_speed{ 2 }, access_distance_lb{ 0.0001 }, access_distance_ub{ 4 }, acecss_link_k{ 4 },
        person_occupancy{ 1 }, desired_speed_ratio{ 1 }, number_of_allowed_links{ 0 }, eco_so_flag{ 0 }, eco_so_flow_switch_bound{ 0 }, pce{ 1 }
    {
    }

    int mode_type_no;
    // dollar per hour
    float value_of_time;
    int eco_so_flag;
    int eco_so_flow_switch_bound;
    // link type, product consumption equivalent used, for travel time calculation
    double person_occupancy;
    double pce; // multimodal equivalent unit: simialr to PCE: passenger car equivalent, but can be extended to bike and walk and other modes
    double desired_speed_ratio;

    float time_headway_in_sec;
    int real_time_information_type;
    std::string mode_type;
    int number_of_allowed_links;

    std::string display_code;

    std::string access_node_type;
    float access_speed;

    float access_distance_lb;
    float access_distance_ub;
    int acecss_link_k;

    std::map<int, bool> zone_id_cover_map;
};

class CLinkType
{
public:
    CLinkType() : link_type{ 1 }, number_of_links{ 0 }, traffic_flow_code{ spatial_queue }, k_jam{ 300 }
    {
        for (int at = 0; at < g_number_of_active_mode_types; at++)
        {
            for (int value_index = 0; value_index < 4; value_index++)
            {
                emissions_co2_matrix[at][value_index] = 0;
                emissions_nox_matrix[at][value_index] = 0;
            }
        }
    }

    int link_type;
    int number_of_links;
    float k_jam;
    std::string link_type_name;
    std::string type_code;

    e_traffic_flow_model traffic_flow_code;

    double emissions_co2_matrix[MAX_MODETYPES][4];
    double emissions_nox_matrix[MAX_MODETYPES][4];
};

class CColumnPath {
public:
    CColumnPath() : path_node_vector{ nullptr }, path_link_vector{ nullptr }, path_seq_no{ 0 }, route_seq_id{ 0 }, m_link_size{ 0 }, m_node_size{ 0 },
        path_switch_volume{ 0 }, path_volume{ 0 }, path_preload_volume{ 0 }, path_volume_before_ODME{ -1 }, path_volume_after_ODME{ -1 }, path_volume_before_dtm{ -1 }, path_volume_after_dtm{ -1 }, path_travel_time{ 0 }, path_distance{ 0 }, path_toll{ 0 }, UE_gap{ 0 }, UE_relative_gap{ 0 },
        path_gradient_cost{ 0 }, path_gradient_cost_difference{ 0 }, path_gradient_cost_relative_difference{ 0 }, subarea_output_flag{ 1 }, measurement_flag{ 0 }, impacted_path_flag{ 0 },
        network_design_detour_mode{ 0 }, b_RT_new_path_flag{0}, global_path_no{ -1 }, b_sensitivity_analysis_flag{ false }
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

        if (m_link_size == 0)
        {
            int i_debug = 1;
        }

        // dynamic array
    try {
        path_node_vector = new int[node_size];
        path_link_vector = new int[link_size];

    }
    catch (const std::bad_alloc& e) {
        std::cout << "Memory allocation failed: " << e.what() << std::endl;
    }

    if (backwardflag)
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
        if (m_link_size == 0)
        {
            int i_debug = 1;
            dtalog.output() << "[ERROR] m_link_size == 0 in function CColumnPath::AllocateVector()!";
            g_DTA_log_file << "[ERROR] m_link_size == 0 in function CColumnPath::AllocateVector()!";
            g_program_stop();
        }

        // dynamic array
        path_node_vector = new int[m_node_size];
        path_link_vector = new int[m_link_size];

        if (backwardflag)
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

    std::vector<int> path_link_STL_vector;
    int path_seq_no;  // path id within an OD pair
    int route_seq_id; // across all OD pairs
    int global_path_no;
    std::string path_id;
    // path volume
    double path_volume = 0;
    double path_preload_volume = 0;
    double path_volume_before_ODME = 0;
    double path_volume_after_ODME = 0;
    double path_volume_before_dtm = 0;
    double path_volume_after_dtm = 0;

    std::vector<float> departure_time_in_min;
    int subarea_output_flag;
    int measurement_flag;
    int impacted_path_flag;
    bool b_sensitivity_analysis_flag;
    int b_RT_new_path_flag;

    int network_design_detour_mode;  // network_design_mode = 1: passing through network design locations, // 2: OD pair passing through networok design, but this path is an alternative path as detour

    double path_switch_volume;
    double path_travel_time;
    double path_distance;
    double path_toll;
    // first order graident cost.
    double path_gradient_cost;
    double UE_gap;
    double UE_relative_gap;

    std::map <int, double> path_time_per_iteration_map;
    std::map <int, double> path_volume_per_iteration_map;

    std::map <int, double> path_time_per_iteration_SA_map;
    std::map <int, double> path_volume_per_iteration_SA_map;

    std::map <int, double> path_time_per_iteration_ODME_map;
    std::map <int, double> path_volume_per_iteration_ODME_map;

    // first order graident cost - least gradient cost
    double path_gradient_cost_difference;
    // first order graident cost - least gradient cost
    double path_gradient_cost_relative_difference;

    int m_node_size;
    int m_link_size;

    std::vector <int> path_sensor_link_vector;  // added with mustafa, 12/24/2022. only contains with the links with measurements
    std::vector <int> path_SA_link_vector;  // added with mustafa, 12/24/2022. only contains with the links with measurements

    int m_sensor_link_size;

    std::map <int, bool> diverted_vehicle_map;

    std::vector<int> agent_simu_id_vector;
};

class CActivityTravelPattern {
public:
    int activity_travel_pattern_index;
    std::string activity_travel_pattern_id;
    std::string mode_chain_str;
    std::string demand_period_chain_str;
    int number_of_trips;
    std::vector<std::string> demand_period_chain;
    std::vector<int> demand_period_vector;
    std::vector<std::string> mode_chain;
    std::vector<int> mode_vector;
    CActivityTravelPattern() : number_of_trips{ 0 }
    {}
};

struct SChoiceAlt {
    int mode_no;
    int demand_peroid_no;
    int o_zone_no;
    int d_zone_no;
    float volume;
};

class CChoiceAlt {
public:
    CChoiceAlt() : volume{ 0 } {}
    std::string multi_dim_choice_id;
    std::string activity_zone_chain_str;

    int choice_alternative_id;
    std::string activity_travel_pattern_id;
    std::vector<int> activity_zone_chain;

    std::vector<SChoiceAlt> activity_chain;
    float volume;
};


class CChoiceSet {
public:
    CChoiceSet() : demand_volume{ 0 }, avg_travel_time{ 0 }, avg_travel_cost{ 0 }
    {

    }
    int choice_set_index;
    std::string multi_dim_choice_id;
    std::string mode_tag, demand_period_tag, spatial_tag, travel_purpose_tag, data_tag;

    //mode_tag: This could be used to categorize models based on the type of travel mode they focus on or include.For example, the tags could be 'car', 'public_transit', 'cycling', 'walking', 'ridesharing', 'autonomous_vehicles', etc.
    //
    //demand_period_tag : This tag could indicate the temporal scope of the model.It could be 'short_term' (for daily or within - day models), 'long_term' (for models dealing with decisions that have long - term impacts like vehicle ownership or residential location), or 'peak' and 'off_peak' (for models specifically focused on peak or off - peak travel periods).
    //
    //spatial_tag : This tag would denote the geographic scale of the model.Tags could include 'neighborhood', 'city', 'metropolitan', 'regional', 'national', or 'global'.
    //
    //travel_purpose_tag : This tag would categorize models based on the travel purpose they represent.For example, tags could be 'work', 'education', 'shopping', 'leisure', 'personal_business', etc.
    //
    //data_tag : This tag would classify models based on the type of data they require.Possible tags could be 'individual_level', 'household_level', 'aggregate_level', 'revealed_preference', 'stated_preference', 'big_data', 'GPS_data', etc.


    std::vector<CChoiceAlt> choice_alt_vector;
    float demand_volume;
    float avg_travel_time;
    float avg_travel_cost;
};


class CAgentPath {
public:
    CAgentPath() : path_id{ 0 }, node_sum{ -1 }, travel_time{ 0 }, distance{ 0 }, volume{ 0 }
    {
    }

    std::string path_id;
    int node_sum;
    float travel_time;
    float distance;
    float volume;

    std::vector <int> path_link_sequence;


};

// Peiheng, 02/21/22 remove extern to make it forward declaration
class NetworkForSP;

class CColumnVector {

public:
    // this is colletion of unique paths
    CColumnVector() :  prev_od_volume{ 0 }, bfixed_route{ false }, m_passing_sensor_flag{ -1 }, information_type{ 0 }, activity_mode_type_no{ 0 },
        departure_time_profile_no{ -1 }, OD_impact_flag{ 0 }, subarea_passing_flag{ 1 }, OD_based_UE_relative_gap{ 0 }, least_travel_time{ 0 }
    {

            od_volume = 0;
            od_volume_from_route_file = 0;
            avg_travel_time = 0;
            avg_distance = 0;
            avg_distance = 0;

    }


    void reset_column_pool()
    {
        path_node_sequence_map.clear();

        OD_based_UE_relative_gap = 0;

    }
    bool subarea_passing_flag;

    std::map<int, bool> at_od_impacted_flag_map; // for each agent type

    float avg_cost;
    float avg_travel_time;
    float avg_distance;
    // od volume
    double od_volume;
    double od_volume_from_route_file;

    double least_travel_time;

    double OD_based_UE_relative_gap;
    std::map<int, double> od_volume_per_iteration_map;

    double prev_od_volume;

    bool bfixed_route;
    int information_type;
    std::vector<int> activity_zone_no_vector;

    std::string  activity_zone_sequence;
    std::string activity_mode_type_sequence;
    std::vector<int>  activity_mode_type_no_vector;

    int activity_mode_type_no;
    int departure_time_profile_no;
    int m_passing_sensor_flag;
    // first key is the sum of node id;. e.g. node 1, 3, 2, sum of those node ids is 6, 1, 4, 2 then node sum is 7.
    // Peiheng, 02/02/21, potential memory leak, fix it
    std::map <int, CColumnPath> path_node_sequence_map;
    std::map<int, CColumnPath>::iterator it, it_begin, it_end;

    int ModifyColumnPoolAfterLoadingRouteFile()
    {
        double ratio = od_volume / max(0.00001, od_volume_from_route_file);

        if(fabs(ratio - 1.0f)>0.001)
        {
            ColumnPoolDemandModification(ratio);
            return 1;
        }

        return 0;
    }

    void ColumnPoolDemandModification(double ratio)
    {
        it_begin = path_node_sequence_map.begin();
        it_end = path_node_sequence_map.end();

        for (it = it_begin; it != it_end; ++it)
        {
            it->second.path_volume *= ratio;
        }
    }

    int OD_impact_flag;  // 0: no passing network design locations; //1: all paths passing through network deign locations: //2: there are alternative detours w.r.t. network design location
};

class CAgent_Column {
public:
    CAgent_Column() : cost{ 0 }, volume{ 0 },
        travel_time{ 0 }, distance{ 0 }, agent_id { 0 }, o_zone_id { 0 }, d_zone_id { 0 }, o_node_id { 0 }, d_node_id { 0 }
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

    std::string mode_type;

    std::vector<int> path_node_vector;
    std::vector<int> path_link_vector;
    std::vector<double> path_time_vector;
};

// event structure in this "event-based" traffic simulation

class DTAVehListPerTimeInterval
{
public:
    std::vector<int> m_AgentIDVector;
};

class CAgent_Simu
{
public:
    CAgent_Simu() : agent_vector_seq_no{ -1 }, path_toll{ 0 }, departure_time_in_min{ 0 }, m_bGenereated{ false }, m_bCompleteTrip{ false },
        path_travel_time_in_min{ 0 }, path_distance{ 0 }, diverted_flag{ 0 }, time_headway{ number_of_simu_interval_reaction_time }, PCE_unit_size{ 1 }, impacted_flag{ -1 }, impacted_link_seq_no{ 99999 },
        info_receiving_flag{ 0 }, desired_free_travel_time_ratio{ 1.0 }, waiting_time_in_min{ 0 }, max_link_waiting_time_in_min{ 0 },
        max_waiting_time_link_no{ -1 }, p_RTNetwork{ NULL }, mode_type_no{ 0 },
        departure_time_in_simu_interval{ 0 }, arrival_time_in_min{ 0 }, m_current_link_seq_no{ 0 }, m_path_link_seq_no_vector_size{ 0 }
    {
    }

    ~CAgent_Simu()
    {
    }

    NetworkForSP* p_RTNetwork;
    void AllocateMemory()
    {
        m_current_link_seq_no = 0;

        m_veh_link_arrival_time_in_simu_interval.reserve(path_link_seq_no_vector.size());
        m_veh_link_departure_time_in_simu_interval.reserve(path_link_seq_no_vector.size());

        for (int i = 0; i < path_link_seq_no_vector.size(); ++i)
        {
            m_veh_link_arrival_time_in_simu_interval.push_back(-1);
            m_veh_link_departure_time_in_simu_interval.push_back(-1);
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
    double departure_time_in_min;
    double waiting_time_in_min;
    double max_link_waiting_time_in_min;
    int max_waiting_time_link_no;

    int info_receiving_flag;
//    int diverted_flag;
    int impacted_flag;
    int diverted_flag;
    std::string impacted_str;
    int impacted_link_seq_no;
    bool m_bGenereated;
    bool m_bCompleteTrip;

    int mode_type_no;
    int agent_id;

    // for column pool index
    int at;
    int tau;
    int dest;

    int m_current_link_seq_no;
    int m_path_link_seq_no_vector_size;

    // for tracking the meso link based statitics: working with Alicia to implement
    std::map<int, int> meso_link_id_map;
    std::vector<int> path_meso_link_id_vector;

    int departure_time_in_simu_interval;
    float arrival_time_in_min;
    float path_travel_time_in_min;
    float path_distance;

    unsigned int m_RandomSeed;

    // external input
    std::vector<int> path_link_seq_no_vector;
    std::vector<int>  m_veh_link_arrival_time_in_simu_interval;
    std::vector<int>  m_veh_link_departure_time_in_simu_interval;
    int time_headway;  // in terms of simulation interval
    int PCE_unit_size;  // the number of units:
    double desired_free_travel_time_ratio;
};

class DTAScenario {
public:

    DTAScenario() : scenario_index{ 0 }, year{ 2023 }, activate{ 1 }, dtm_flag{ 0 }
    {
    }

    int scenario_index;
    int year;
    std::string scenario_name;
    std::string scenario_description;
    int activate;
    int dtm_flag;
};
class Assignment {
public:
    // default is UE
    Assignment() : assignment_mode{ lue }, g_number_of_cpu_processors{ 4 }, g_number_of_threads{ 1 }, g_info_updating_freq_in_min{ 5 }, g_visual_distance_in_cells{ 5 },
        g_link_type_file_loaded{ true }, g_mode_type_file_loaded{ false }, total_route_demand_volume{ 0 }, total_real_time_demand_volume{ 0 }, g_column_pool{ nullptr }, g_number_of_in_memory_simulation_intervals{ 500 },
        g_number_of_column_generation_iterations{ 20 }, g_number_of_column_updating_iterations{ 0 }, g_number_of_ODME_iterations{ 0 }, g_number_of_sensitivity_analysis_iterations_for_dtm{ -1 }, g_number_of_links{ 0 }, g_number_of_timing_arcs{ 0 },
        g_number_of_nodes{ 0 }, g_number_of_zones{ 0 }, g_number_of_mode_types{ 0 }, debug_detail_flag{ 1 }, path_output{ 1 }, trajectory_output_count{ -1 },
        trace_output{ 0 }, major_path_volume_threshold{ 0.1 }, trajectory_sampling_rate{ 1.0 }, td_link_performance_sampling_interval_in_min{ 1 }, dynamic_link_performance_sampling_interval_hd_in_min{ 15 }, trajectory_diversion_only{ 0 }, m_GridResolution{ 0.01 },
        shortest_path_log_zone_id{ 1 }, g_number_of_analysis_districts{ 1 },
         g_length_unit_flag{ 0 }, g_speed_unit_flag{ 0 }, active_dms_count{ 0 }, active_lane_closure_count{ 0 }, g_number_of_real_time_mode_types{ 0 }, g_number_of_DMS_mode_types{ 0 }, g_first_link_type{ -1 },
        g_max_number_of_super_zones{ 100 }, b_forward_star_structure_log{ 0 }, b_sp_log{ 0 }, starting_time_slot_no{ 0 }, ending_time_slot_no{ 12 /* for 1 hour with every 5 min intervals*/}

    {
        m_LinkCumulativeArrivalVector  = NULL;
        m_LinkCumulativeDepartureVector = NULL;

        m_link_CA_count = NULL;  // CA, assign this value to m_LinkCumulativeArrivalVector at a given time in min
        m_link_CD_count = NULL;
        m_LinkOutFlowCapacity = NULL;
        m_LinkOutFlowState =  NULL;

       // sp_log_file.open("log_label_correcting.txt");
        //assignment_log_file.open("log_traffic_assignment.csv");
       /* assignment_log_file << "iteration_no,link_id,from_node_id,to_node_id,volume,travel_time" << '\n';*/

   //     log_subarea_focusing_file.open("log_subarea_focusing.txt");



        simu_log_file.open("log_simulation.txt");


//        simu_log_file << "start" << '\n';
        g_rt_network_pool = NULL;
        g_column_pool = NULL;
    }

    ~Assignment()
    {
        if (g_column_pool)
            Deallocate3DDynamicArray(g_column_pool, g_related_zone_vector_size, g_related_zone_vector_size);

        sp_log_file.close();
    //    log_subarea_focusing_file.close();
        summary_file.close();
        summary_file2.close();
        summary_corridor_file.close();
        summary_system_file.close();
        simu_log_file.close();
        assignment_log_file.close();
        DeallocateLinkMemory4Simulation();
    }

    void InitializeDemandMatrix(int number_of_signficant_zones, int number_of_zones, int number_of_mode_types)
    {
        g_number_of_zones = number_of_zones;
        g_number_of_mode_types = number_of_mode_types;

        total_demand_volume = 0;
        g_column_pool = Allocate3DDynamicArray<CColumnVector>(number_of_signficant_zones, g_related_zone_vector_size, max(1, number_of_mode_types));

        for (int i = 0; i < number_of_zones; ++i)
        {
            g_origin_demand_array[i] = 0.0;
        }

        for (int i = 0; i < number_of_mode_types; ++i)
        {
                total_demand[i] = 0.0;
                total_route_demand[i] = 0.0;
        }

        g_DemandGlobalMultiplier = 1.0f;
    }

    int get_in_memory_time(int t)
    {
        return t % g_number_of_in_memory_simulation_intervals;
    }

    std::vector<DTAGDPoint> g_subarea_shape_points;
    int g_max_number_of_super_zones;

    void STTrafficSimulation();
    void STMesoTrafficSimulation();

    void Demand_ODME(int OD_updating_iterations);
    void Sensor_Vector_based_Demand_ODME(int OD_updating_iterations);
    void AllocateLinkMemory4Simulation();
    //int update_real_time_info_path(CAgent_Simu* p_agent, int& impacted_flag_change, float updating_in_min);
    bool RTSP_real_time_travel_time_updating(int time_slot_no, int simu_interval_t);
    void DeallocateLinkMemory4Simulation();

    std::map<int, int> zone_id_to_centriod_node_no_mapping;  // this is an one-to-one mapping
    std::map<int, int> zone_id_2_node_no_mapping;  // this is used to mark if this zone_id has been identified or not
    std::map<int, __int64> zone_id_2_cell_id_mapping;  // this is used to mark if this zone_id has been identified or not
    std::map<__int64, int> cell_id_mapping;  // this is used to mark if this cell_id has been identified or not
    std::map<__int64, std::string> cell_id_2_cell_code_mapping;  // this is used to mark if this cell_id has been identified or not


    double m_GridResolution;
    e_assignment_mode assignment_mode;

    int active_dms_count;
    int active_lane_closure_count;
    int g_number_of_cpu_processors;
    int g_visual_distance_in_cells;
    float g_info_updating_freq_in_min;

    int g_number_of_threads;
    int path_output;
    int trajectory_output_count;
    int trace_output;
    float trajectory_sampling_rate;
    int trajectory_diversion_only;
    int td_link_performance_sampling_interval_in_min;
    float dynamic_link_performance_sampling_interval_hd_in_min;

    float major_path_volume_threshold;
    int shortest_path_log_zone_id;

    bool g_link_type_file_loaded;
    bool g_mode_type_file_loaded;

    float total_demand_volume;
    float total_real_time_demand_volume;

    float total_route_demand_volume;
    std::map<int, float> g_origin_demand_array;
    CColumnVector*** g_column_pool;
    NetworkForSP* ** g_rt_network_pool;

    // the data horizon in the memory
    int g_number_of_in_memory_simulation_intervals;
    int g_number_of_column_generation_iterations;
    int g_number_of_sensitivity_analysis_iterations_for_dtm;
    int g_number_of_column_updating_iterations;
    int g_max_num_significant_zones_in_subarea;
    int g_max_num_significant_zones_outside_subarea;

    int g_number_of_ODME_iterations;
    int g_length_unit_flag;
    int g_speed_unit_flag;

    int g_number_of_links;
    int g_number_of_timing_arcs;
    int g_number_of_nodes;
    int b_forward_star_structure_log;
    int b_sp_log;
    int g_number_of_zones;
    int g_number_of_mode_types;
    int g_number_of_real_time_mode_types;
    int g_number_of_DMS_mode_types;

    int starting_time_slot_no;
    int ending_time_slot_no;
    float time_period_in_hour;
    float t2_peak_in_hour;


    std::map<int, int> node_seq_no_2_zone_id_mapping;  // this is used to mark if this zone_id has been identified or not
    std::map<int, int> zone_seq_no_2_info_mapping;  // this is used to mark if this zone_id has been identified or not
    std::map<int, int> zone_seq_no_2_activity_mapping;  // this is used to mark if this zone_id has been identified or not

    std::map<int, int> zone_id_to_centriod_node_id_mapping;  // this is an one-to-one mapping
    std::map<int, int> zone_id_to_seed_zone_id_mapping;  // this is an one-to-one mapping


    int debug_detail_flag;

    // hash table, map external node number to internal node sequence no.
    std::map<int, int> g_node_id_to_seq_no_map;
    std::map<int, int> access_node_id_to_zone_id_map;

    std::map<int, int> g_zone_seq_no_to_analysis_distrct_id_mapping;

    // from integer to integer map zone_id to zone_seq_no
    std::map<int, int> g_zoneid_to_zone_seq_no_mapping;
    std::map<int, int> g_zoneid_to_zone_sindex_no_mapping;  //subarea based index

    std::map<std::string, int> g_link_id_map;

    std::map<int, double> zone_id_X_mapping;
    std::map<int, double> zone_id_Y_mapping;

  //  CDemand_Period g_DemandPeriodVector;
    std::vector<CDeparture_time_Profile> g_DepartureTimeProfileVector;

    int g_LoadingStartTimeInMin;
    int g_LoadingEndTimeInMin;

    std::vector<Cmode_type> g_ModeTypeVector;


    std::map<std::string, CActivityTravelPattern> g_ActivityTravelPatternMap;
    std::map<std::string, CChoiceSet> g_ChoiceSetMap;


    int g_number_of_analysis_districts;
    std::map<int, CLinkType> g_LinkTypeMap;
    int g_first_link_type;

    std::map<std::string, int> mode_type_2_seqno_mapping;

    std::map<int, double> o_district_id_factor_map;
    std::map<int, double> d_district_id_factor_map;
    std::map<int, double> od_district_id_factor_map;

    std::map<int, double> SA_o_district_id_factor_map;
    std::map<int, double> SA_d_district_id_factor_map;
    std::map<int, double> SA_od_district_id_factor_map;

    float total_demand[MAX_MODETYPES];
    float total_route_demand[MAX_MODETYPES];
    float g_DemandGlobalMultiplier;

    // used in ST Simulation
    float** m_LinkOutFlowCapacity;  // per second interval for simplicity
    int** m_LinkOutFlowState;  // per second interval for simplicity

    // in min
    float** m_link_TD_waiting_time;
    std::vector<float> m_link_total_waiting_time_vector;;
    // number of simulation time intervals

    float** m_LinkCumulativeArrivalVector;
    float** m_LinkCumulativeDepartureVector;

    float* m_link_CA_count;  // CA, assign this value to m_LinkCumulativeArrivalVector at a given time in min
    float* m_link_CD_count;  // CD

    int g_start_simu_interval_no;
    int g_number_of_simulation_intervals;
    // is shorter than g_number_of_simulation_intervals
    int g_number_of_loading_intervals_in_sec;
    // the data horizon in the memory in min
    int g_number_of_intervals_in_min;

    int g_number_of_intervals_in_sec;

    std::map<std::string, int> m_TMClink_map;
    std::map<std::string, int> m_TMC_corridor_map;

    std::ofstream simu_log_file;

    std::ofstream sp_log_file;
    std::ofstream assignment_log_file;

    std::ofstream log_subarea_focusing_file;

    std::ofstream summary_file;
    std::ofstream summary_file2;
    std::ofstream summary_corridor_file;
    std::ofstream summary_system_file;

};

extern Assignment assignment;

#include "VDF.h"

class CLink
{
public:
    // construction
    CLink() :main_node_id{ -1 }, free_speed{ 100 }, v_congestion_cutoff{ 100 }, v_critical{ 60 },
        length_in_meter{ 1 }, link_distance_VDF{ 0.001 }, link_DOC{ 0 },
        BWTT_in_simulation_interval{ 100 }, zone_seq_no_for_outgoing_connector{ -1 }, lane_capacity{ 1999 },
        free_flow_travel_time_in_min{ 0.01 }, link_spatial_capacity{ 100 },
        timing_arc_flag{ false }, traffic_flow_code{ 0 }, spatial_capacity_in_vehicles{ 999999 }, subarea_id{ -1 }, RT_flow_volume{ 0 },
        cell_type{ -1 }, saturation_flow_rate{ 1800 }, dynamic_link_event_start_time_in_min{ 99999 }, b_automated_generated_flag{ false }, time_to_be_released{ -1 },
        RT_waiting_time{ 0 }, FT{ 1 }, AT{ 1 }, s3_m{ 4 }, tmc_road_order{ 0 }, tmc_road_sequence{ -1 },
        tmc_corridor_id{ -1 }, from_node_id{ -1 }, to_node_id{ -1 }, kjam{ 300 }, link_distance_km{ 0 }, link_distance_mile{ 0 }, meso_link_id{ -1 }, total_simulated_delay_in_min{ 0 },
        total_simulated_meso_link_incoming_volume{ 0 }, global_minute_capacity_reduction_start{ -1 }, global_minute_capacity_reduction_end{ -1 },
         AB_flag{ 1 }, BA_link_no{ -1 }, cost{ 0 }, win_count{ 0 }, lose_count{ 0 },
        vdf_data_count{ 0 }, peak_load_factor{ 1.0 }, Q_cd{ 1.0 }, Q_n{ 1.0 }, Q_cp{ 0.28125 /*0.15*15/8*/ }, Q_s{ 4 }, vf{ 100 }, v_VDF_congestion_cutoff{ 45 }, vt2{ -1 },
        alpha{ 0.15 }, beta{ 4 },  rho{ 1 }, preload{ 0 }, penalty{ 0 }, RT_route_regeneration_penalty{ 0 }, lane_closure_final_lanes{ 0 },
        starting_time_in_hour{ 0 }, ending_time_in_hour{ 0 },
        volume_before_odme{ 0 }, volume_after_odme{ 0 },

        cycle_length{ -1 }, red_time{ 0 }, effective_green_time{ 0 }, t0{ -1 }, t3{ -1 }, start_green_time{ -1 }, end_green_time{ -1 }, L{ 1 },
        queue_length{ 0 }, obs_count{ 0 }, upper_bound_flag{ 1 }, est_count_dev{ 0 }, avg_waiting_time{ 0 }, P{ -1 }, Severe_Congestion_P{ -1 }, lane_based_D{ 0 }, lane_based_Vph{ 0 }, avg_speed_BPR{ -1 }, avg_queue_speed{ -1 }, number_of_lanes{ 1 }, sa_volume{ 0 }, t2{ 1 }, k_critical{ 45 }, link_main_volume{ 0 },
        Q_mu{ 0 }, Q_gamma{ 0 }, dynamic_traffic_management_flag{ 0 },
        volume_before_dtm{ 0 }, speed_before_dtm{ 0 }, DoC_before_dtm{ 0 }, P_before_dtm{ 0 },
        volume_after_dtm{ 0 }, speed_after_dtm{ 0 }, DoC_after_dtm{ 0 }, P_after_dtm{ 0 },
        ref_link_volume{ -1 }, free_speed_diff_link_specific{ 0 }, capacity_diff_link_specific{ 0 }
    {
        for (int at = 0; at < g_number_of_active_mode_types; at++)
        {
            {
                toll[at] = 0;
            }

            {
                allowed_uses = "";
            }
        }
    }

    float PerformSignalVDF(float hourly_per_lane_volume, float red, float cycle_length)
    {
        float lambda = hourly_per_lane_volume;
        float mu = _default_saturation_flow_rate; //default saturation flow ratesa
        float s_bar = 1.0 / 60.0 * red * red / (2 * cycle_length); // 60.0 is used to convert sec to min
        float uniform_delay = s_bar / max(1 - lambda / mu, 0.1f);

        return uniform_delay;
    }
    int vdf_data_count;

    double get_speed_from_volume(float hourly_volume, float vf, float k_critical, float s3_m)
    {
        //test data free_speed = 55.0f;
        //speed = 52;
        //k_critical = 23.14167648;

        double max_lane_capacity = k_critical * vf / pow(2, 2 / s3_m);

        if (hourly_volume > max_lane_capacity)
            hourly_volume = max_lane_capacity;
        // we should add a capacity upper bound on hourly_volume;

        double coef_a = pow(k_critical, s3_m);
        double coef_b = pow(k_critical, s3_m) * pow(vf, s3_m / 2.0);
        double coef_c = pow(hourly_volume, s3_m);  // D is hourly demand volume, which is equivalent to flow q in S3 model

        double speed = pow((coef_b + pow(coef_b * coef_b - 4.0 * coef_a * coef_c, 0.5)) / (2.0 * coef_a), 2.0 / s3_m);    //under uncongested condition
        if (speed >= vf)
            speed = vf;

        if (speed < 0)
            speed = 0;

        return speed;
    }

    double get_volume_from_speed(float speed, float vf, float k_critical, float s3_m)
    {
        //test data free_speed = 55.0f;
        //speed = 52;
        //k_critical = 23.14167648;

        if (speed < 0)
            return -1;

        double speed_ratio = vf / max(1.0f, speed);
        if (speed_ratio <= 1.00001)
            speed_ratio = 1.00001;

        /*   float volume = 0;*/
        double ratio_difference = pow(speed_ratio, s3_m / 2) - 1;

        double ratio_difference_final = max(ratio_difference, 0.00000001);

        double volume = speed * k_critical * pow(ratio_difference_final, 1 / s3_m);  // volume per hour per lane

        return volume;
    }

    double calculate_travel_time_based_on_QVDF(int at, double fftt, double volume, double mode_hourly_capacity, double peak_load_factor,
        float model_speed[MAX_TIMEINTERVAL_PerDay], float est_volume_per_hour_per_lane[MAX_TIMEINTERVAL_PerDay],
        CLinkType link_type/*, double link_avg_co2_emit_per_mode[MAX_MODETYPES], double link_avg_nox_emit_per_mode[MAX_MODETYPES]*/
    )
    {
        // step 0: setup P' initial value and time_period_in_hour, obtain Q_peak_load_factor,
        P = -0.5;

        double time_period_in_min = max(0.1, (ending_time_in_hour - starting_time_in_hour) * 60);
        double time_period_in_hour = max(0.1, (ending_time_in_hour - starting_time_in_hour));

        double avg_travel_time = 0;
        // QVDF
        double dc_transition_ratio = 1;


        // step 1: calculate lane_based D based on plf and number_of_lanes from link volume V over the analysis period  take nonnegative values
        lane_based_D = max(0.0, volume) / time_period_in_hour / max(0.000001, number_of_lanes) / peak_load_factor;
        // step 2: D_ C ratio based on lane-based D and lane-based ultimate hourly capacity,
        // uncongested states D <C
        // congested states D > C, leading to P > 1
        double DOC = lane_based_D / max(0.00001, mode_hourly_capacity);

        if (number_of_lanes < 0.6)  // dynamic lane closure scenario, we computing D, we assume number_of_lanes = 1
        {
            lane_based_D = max(0.0, volume) / time_period_in_hour / peak_load_factor;
            DOC = lane_based_D / max(0.00001, mode_hourly_capacity * number_of_lanes);

        }


        double congestion_ref_speed = v_VDF_congestion_cutoff;
        if (DOC < 1)
            congestion_ref_speed = (1 - DOC) * this->free_speed +  DOC * this->v_VDF_congestion_cutoff;



        if (volume > 1)
            int iii_debug = 1;

        //step 3.1 fetch vf and v_congestion_cutoff based on fftt, VCTT (to be compartible with transit data, such as waiting time )
        // we could have a period based fftt, so we need to convert fftt to vfree
        // if we only have one period, then we can directly use vf and v_congestion_cutoff.

        //step 3.2 calculate speed from VDF based on D/C ratio
        avg_queue_speed = congestion_ref_speed / (1.0 + alpha * pow(DOC, beta));
        // step 3.3 taking the minimum of BPR- v and Q VDF v based on log sum function

       // let us use link_length_in_km = 1 for the following calculation
        double  link_length_in_1km = 1.0;
        double RTT = 0;
        RTT = link_length_in_1km / v_VDF_congestion_cutoff;
        double Q_n_current_value = Q_n;
        link_DOC = DOC;
        //// will revisit again
        //if (DOC < dc_transition_ratio)  // free flow regime
        //{

        //    double vf_alpha = (1.0 + alpha) * vf / max(0.0001, v_VDF_congestion_cutoff) - 1.0;
        //    // fixed to pass through vcutoff point vf/ (1+vf_alpha) = vc / (1+ VDF_alpha) ->
        //    // 1+vf_alpha = vf/vc *(1+VDF_alpha)
        //    // vf_qlpha =  vf/vc *(1+VDF_alpha) - 1
        //    // revised BPR DC
        //    double vf_beta = beta; // to be calibrated
        //    double vf_avg_speed = vf / (1.0 + vf_alpha * pow(DOC, vf_beta));
        //    avg_queue_speed = vf_avg_speed; // rewrite with vf based speed
        //    Q_n_current_value = beta;
        //    RTT = link_length_in_1km / max(0.01, vf);  // in hour
        //}

        // BPR
             // step 2: D_ C ratio based on lane-based D and lane-based ultimate hourly capacity,
             // uncongested states D <C
             // congested states D > C, leading to P > 1


        //step 3.1 fetch vf and v_congestion_cutoff based on fftt, VCTT (to be compartible with transit data, such as waiting time )
        // we could have a period based fftt, so we need to convert fftt to vfree
        // if we only have one period, then we can directly use vf and v_congestion_cutoff.

        //step 3.2 calculate speed from VDF based on D/C ratio
        avg_speed_BPR = free_speed / (1.0 + alpha * pow(DOC, beta));
        //avg_travel_time = fftt * (1+ alpha * pow(DOC, beta)); // Mark: fftt should be vctt


        avg_travel_time = fftt * vf / max(0.1, avg_queue_speed); //


        avg_waiting_time = avg_travel_time - fftt;
        //step 4.4 compute vt2
    //            vt2 = avg_queue_speed * 8.0 / 15.0;  // 8/15 is a strong assumption


        P = Q_cd * pow(DOC, Q_n_current_value);  // applifed for both uncongested and congested conditions

        double base = Q_cp * pow(P, Q_s) + 1.0;
        vt2 = v_VDF_congestion_cutoff / max(0.001, base);
        //step 4.1: compute congestion duration P

       //step 4.2 t0 and t3
        t0 = t2 - 0.5 * P;
        t3 = t2 + 0.5 * P;

        if (P >= 2)
        {
            int idebug = 1;
        }
        // work on congested condition
        //step 4.3 compute mu
        Q_mu = min(mode_hourly_capacity, lane_based_D / max(0.01, P));

        //use  as the lower speed compared to 8/15 values for the congested states

        double wt2 = link_length_in_1km / vt2 - RTT; // in hour


        //step 5 compute gamma parameter is controlled by the maximum queue
        Q_gamma = wt2 * 64 * Q_mu / pow(P, 4);  // because q_tw = w*mu =1/4 * gamma (P/2)^4, => 1/vt2 * mu = 1/4 * gamma  * (P/2)^4


        double td_w = 0;
        //step scan the entire analysis period
        Severe_Congestion_P = 0;


        for (int t_in_min = starting_time_in_hour * 60; t_in_min <= ending_time_in_hour * 60; t_in_min += 5)
        {
            double t = t_in_min / 60.0;  // t in hour
            double td_queue = 0;
            double td_speed = 0;

            if (t0 <= t && t <= t3)  // within congestion duration P
            {
                //1/4*gamma*(t-t0)^2(t-t3)^2
                td_queue = 0.25 * Q_gamma * pow((t - t0), 2) * pow((t - t3), 2);
                td_w = td_queue / max(0.001, Q_mu);
                //L/[(w(t)+RTT_in_hour]
                td_speed = link_length_in_1km / (td_w + RTT);
            }
            else if (t < t0) //outside
            {
                td_queue = 0;
                double factor = (t - starting_time_in_hour) / max(0.001, t0 - starting_time_in_hour);
                td_speed = (1 - factor) * vf + factor * max(v_VDF_congestion_cutoff, avg_queue_speed);
            }
            else  // t> t3
            {
                td_queue = 0;
                double factor = (t - t3) / max(0.001, ending_time_in_hour - t3);
                td_speed = (1 - factor) * max(v_VDF_congestion_cutoff, avg_queue_speed) + (factor)*vf;
            }

            // dtalog.output() << "td_queue t" << t << " =  " << td_queue << ", speed =" << td_speed << '\n';
            // g_DTA_log_file << "td_queue t" << t << " =  " << td_queue << ", speed =" << td_speed << '\n';

            int t_interval = t_in_min / 5;

            if (t_in_min <= 410)
            {
                int idebug = 1;
            }
            double td_flow = 0; // default: get_volume_from_speed(td_speed, vf, k_critical, s3_m);
            model_speed[t_interval] = td_speed;
            est_volume_per_hour_per_lane[t_interval] = td_flow;

            if (td_speed < vf * 0.5)
                Severe_Congestion_P += 5.0 / 60;  // 5 min interval
        }


        //// peak load duration
        //double pl_t0 = t2 - max(0.5, 0.5 * P);
        //double pl_t3 = t2 + max(0.5, 0.5 * P);
        //double est_peak_load_demand = 0;
        ////est_non_peak_load_demand should not be counted, as avg non-peak rates have been determined by (V-D)/(L-P)

        //double hourly_rate_2_volume_factor = number_of_lanes / 12.0;  // /12 to convert hourly to 5 min volume;
        //// step 2
        //for (int t_in_min = starting_time_in_hour * 60; t_in_min < ending_time_in_hour * 60; t_in_min += 5)
        //{
        //    double t = t_in_min / 60.0;  // t in hour
        //    int t_interval = t_in_min / 5;
        //
        //     if (t >= pl_t0  && t <= pl_t3)
        //        {
        //         est_peak_load_demand += est_volume_per_hour_per_lane[t_interval] * hourly_rate_2_volume_factor;
        //        }
        //}
        //// step 3:
        //double peak_load_volume_scale_factor = lane_based_D / max(0.0001,est_peak_load_demand);


        ////step 4
        //for (int t_in_min = starting_time_in_hour * 60; t_in_min < ending_time_in_hour * 60; t_in_min += 5)
        //{
        //    double t = t_in_min / 60.0;  // t in hour
        //    int t_interval = t_in_min / 5;

        //    if (t < pl_t0)
        //    {
        //        est_volume_per_hour_per_lane[t_interval] = min(static_cast<float_t>(lane_based_ultimate_hourly_capacity), est_volume_per_hour_per_lane[t_interval]);
        //    }
        //    else if (t > pl_t3)
        //    {
        //        est_volume_per_hour_per_lane[t_interval] = min(static_cast<float_t>(lane_based_ultimate_hourly_capacity), est_volume_per_hour_per_lane[t_interval]);
        //    }
        //    else
        //    {
        //        est_volume_per_hour_per_lane[t_interval] = min(lane_based_ultimate_hourly_capacity, est_volume_per_hour_per_lane[t_interval] * peak_load_volume_scale_factor);
        //    }
        //}


        ////final stage: compute avg emission in peak period
        // vq: speed in miles per hour, converted from km per hour

        // apply final travel time range constraints


        if (avg_travel_time > max(15.0, time_period_in_min * 1.5))  // use 1.5 times to consider the some wide range bound
            avg_travel_time = max(15.0, time_period_in_min * 1.5);

        double vf_mph = vf / 1.609;
        double vq = vf_mph / max(0.00001, avg_travel_time / fftt) / 1.609;

        // vf_minus_vq: difference between vf and vq
        double vf_minus_vq = vf_mph - vq;

        // waiting_time_w: waiting time, computed as the difference between average travel time and Free-Flow Travel Time (fftt)
        double waiting_time_w = avg_travel_time - fftt;

        // set speed v to be equal to vq
        double v = vq;

        // lambda_emission: CO2 emission rate, computed using a quadratic function of speed v
        double lambda_emission = v * v * link_type.emissions_co2_matrix[at][1] + v * link_type.emissions_co2_matrix[at][2] + link_type.emissions_co2_matrix[at][3];
        double ratio = 0.0;

        // if the absolute difference between vf and vq is greater than 1
        if (fabs(vf_mph - vq) > 1)
        {
            // update the ratio to adjust for changes in emission rate with respect to changes in speed
            ratio = (lambda_emission * vf_mph - vq) / (vf_mph - vq);
        }

        // compute the emission rate as the product of the coefficient and the sum of fftt and waiting time, adjusted for speed changes
        double emission_rate = link_type.emissions_co2_matrix[at][0] * (fftt / 60.0 + waiting_time_w / 60.0 * ratio);

        if (emission_rate < -1)
        {
            int debug_flag = 1;
        }
        // store the computed total CO2 emissions for the mode type in the link_avg_co2_emit_per_mode matrix
        link_avg_co2_emit_per_mode[at] = emission_rate / 1000.0;  // convert to kg


        //nox emissions

        // compute the NOx emission rate using a similar process to that of CO2 emissions
        lambda_emission = v * v * link_type.emissions_nox_matrix[at][1] + v * link_type.emissions_nox_matrix[at][2] + link_type.emissions_nox_matrix[at][3];
        ratio = 0.0;

        if (fabs(vf_mph - vq) > 1)
        {
            ratio = (lambda_emission * vf_mph - vq) / (vf_mph - vq);
        }

        emission_rate = link_type.emissions_nox_matrix[at][0] * (fftt / 60.0 + waiting_time_w / 60.0 * ratio);

        // store the computed total NOx emissions for the mode type in the link_avg_nox_emit_per_mode matrix
        link_avg_nox_emit_per_mode[at] = emission_rate / 1000.0;  // convert to kg;


        // to do for Mohammad
        //for (int t_in_min = starting_time_in_hour * 60; t_in_min <= ending_time_in_hour * 60; t_in_min += 5)
        //{
        //    double t = t_in_min / 60.0;  // t in hour
        //    int t_interval = t_in_min / 5;

        //    double v = model_speed[t_interval] / 1.609;  // convert kmph to mph internally
        //    double lambda_emission = v * v * link_type.emissions_co2_matrix[at][1] + v * link_type.emissions_co2_matrix[at][2] + link_type.emissions_co2_matrix[at][3];
        //    double waiting_time_w  =
        //    double vf_minus_vq = vf-v;

        //    double emission_rate = link_type.emissions_co2_matrix[at][0] * (fftt, + );



        //    /*CLinkType link_type, double link_avg_co2_emit_per_mode[MAX_MODETYPES], double link_avg_nox_emit_per_mode[MAX_MODETYPES]*/
        //}

        return avg_travel_time;
    }

    double calculate_travel_time_based_on_BPR(int at, double fftt, double volume, double mode_hourly_capacity, double peak_load_factor,
        float model_speed[MAX_TIMEINTERVAL_PerDay], float est_volume_per_hour_per_lane[MAX_TIMEINTERVAL_PerDay],
        CLinkType link_type/*, double link_avg_co2_emit_per_mode[MAX_MODETYPES], double link_avg_nox_emit_per_mode[MAX_MODETYPES]*/
    )
    {
        double time_period_in_min = max(0.1, (ending_time_in_hour - starting_time_in_hour) * 60);
        double time_period_in_hour = max(0.1, (ending_time_in_hour - starting_time_in_hour));

        double avg_travel_time = 0;
        // QVDF
        double dc_transition_ratio = 1;


        // step 1: calculate lane_based D based on plf and number_of_lanes from link volume V over the analysis period  take nonnegative values
        lane_based_D = max(0.0, volume) / time_period_in_hour / max(0.000001, number_of_lanes) / peak_load_factor;
        // step 2: D_ C ratio based on lane-based D and lane-based ultimate hourly capacity,
        // uncongested states D <C
        // congested states D > C, leading to P > 1
        double DOC = lane_based_D / max(0.00001, mode_hourly_capacity);

        if (number_of_lanes < 0.6)  // dynamic lane closure scenario, we computing D, we assume number_of_lanes = 1
        {
            lane_based_D = max(0.0, volume) / time_period_in_hour / peak_load_factor;
            DOC = lane_based_D / max(0.00001, mode_hourly_capacity * number_of_lanes);
        }

        if (volume > 1)
            int iii_debug = 1;

        // BPR
            // step 2: D_ C ratio based on lane-based D and lane-based ultimate hourly capacity,
            // uncongested states D <C
            // congested states D > C, leading to P > 1
        link_DOC = DOC;
        //step 3.1 fetch vf and v_congestion_cutoff based on fftt, VCTT (to be compartible with transit data, such as waiting time )
        // we could have a period based fftt, so we need to convert fftt to vfree
        // if we only have one period, then we can directly use vf and v_congestion_cutoff.

        //step 3.2 calculate speed from VDF based on D/C ratio
        avg_speed_BPR = vf / (1.0 + alpha * pow(DOC, beta));
        avg_travel_time = fftt * (1 + alpha * pow(DOC, beta)); // Mark: fftt should be vctt

        //if (cycle_length >= 1)  // signal delay
        //{
        //    float s_bar = 1.0 / 60.0 * red_time * red_time / (2 * cycle_length); // 60.0 is used to convert sec to min
        //    double lambda = lane_based_D;
        //    float uniform_delay = s_bar / max(1 - lambda / saturation_flow_rate, 0.1);
        //    avg_travel_time = uniform_delay + fftt;
        //}

        if (DOC > 0.0001 && avg_travel_time > 10)
        {
            int idebug;
            idebug = 1;
        }

        avg_waiting_time = avg_travel_time - fftt;
        return avg_travel_time;
    }
    std::map<int, float> turn_link_count_map;
    std::map<int, int> restricted_turn_nodes_map;
    string restricted_turn_nodes_str;
    //double DOC;
    //double VOC;

    //updated BPR-X parameters
    double vt2;
    //peak hour factor
    double alpha;
    double beta;
    double ref_link_volume;
    //    double BPR_period_capacity_at[MAX_MODETYPES];
    double peak_load_factor;
    double Q_cd;
    double Q_cp;
    double Q_n;
    double Q_s;
    double Q_mu;
    double Q_gamma;

    double volume_before_odme;
    double volume_after_odme;
    double obs_count;
    int upper_bound_flag;
    double est_count_dev;

    string dtm_scenario_code;
    double volume_before_dtm;
    double speed_before_dtm;
    double DoC_before_dtm;
    double P_before_dtm;

    double volume_after_dtm;
    double speed_after_dtm;
    double DoC_after_dtm;
    double P_after_dtm;

    double starting_time_in_hour;
    double ending_time_in_hour;
    double t2;
    double k_critical;

    double v_VDF_congestion_cutoff;
    double vf;

    double sa_volume;
    double lane_closure_final_lanes;

    int dynamic_traffic_management_flag; // 0: normal: 1: adding lanes, -1: capacity reduction: 2: VMS: -2: induced delay
    double preload;
    double toll[MAX_MODETYPES];

    double free_speed_diff_link_specific;
    double capacity_diff_link_specific;

    double dsr[MAX_MODETYPES]; // desired speed ratio with respect to free-speed
    double penalty;
    double RT_route_regeneration_penalty;

    string allowed_uses;

    double rho;
    //    double marginal_base;
        // in 15 min slot
    float cycle_length;
    float red_time;
    float effective_green_time;
    float saturation_flow_rate;
    int start_green_time;
    int end_green_time;

    double t0, t3;

    bool bValidQueueData;
    string period;

    //    double period_capacity;  // link based period_capacity  //depreciated; will not be used.
    double lane_based_ultimate_hourly_capacity;
    double lane_based_ultimate_hourly_cap_at[MAX_MODETYPES];

    double number_of_lanes;

    // double fftt;
    double P;
    double Severe_Congestion_P;
    double L;
    double lane_based_D;
    double lane_based_Vph;
    double link_DOC;

    double avg_speed_BPR;  // normal BPR
    double avg_queue_speed;  // queue VDF speed
    // inpput
    double link_main_volume;
    //    std::map <int, double> link_volume_per_iteration_map;

        //output
    double avg_delay_0;  // 0 mean driving mode
    double avg_travel_time_0; // 0 mean driving mode

    //// t starting from starting_time_slot_no if we map back to the 24 hour horizon
    float queue_length;
    float arrival_flow_volume;
    float discharge_rate;  // period based
    float avg_waiting_time;
    float travel_time;

    void allocate_memory()
    {
        link_type = 0;
        number_of_lanes = 1;  // default all open
        free_speed = 100;
        capacity = 2000;

        total_volume_for_all_mode_types = 0;
        total_person_volume_for_all_mode_types = 0;
        queue_link_distance_VDF_perslot = 0;
        //cost_perhour = 0;
        for (int at = 0; at < g_number_of_active_mode_types; ++at)
        {
            volume_per_mode_type[at] = 0;
            converted_PCE_volume_per_at[at] = 0;
            link_avg_travel_time[at] = 0;

            link_avg_co2_emit_per_mode[at] = 0;
            link_avg_nox_emit_per_mode[at] = 0;
        }
    }

    ~CLink()
    {
        //In our open source package, we dynamically allocate instances of CLinkand store their pointers in g_link_vector for later use.Considering the fact that these instances are still being used via the pointers in g_link_vector, we don't explicitly free memory in ~CLink().

        //    Let's clarify the typical use of a destructor in C++. It is mainly to free resources that the object owns when it goes out of scope or is explicitly deleted. This can include memory, file descriptors, database connections, and other resources. However, the destructor is not meant to free resources that are still being used elsewhere, which is the case in our package.

        //    If the object is still being used via g_link_vector, deleting the object would lead to undefined behavior, as there would be dangling pointers in g_link_vector pointing to memory that's no longer valid.

        //    The proper approach that we've adopted is to keep track of all dynamically allocated CLink objects and delete them when they are no longer needed, i.e., when they are no longer accessible via g_link_vector. We handle this in a clean-up function or similar part of our code, or in the destructor of the class/structure that owns g_link_vector, making sure that by the time this destructor is called, there's no more need for the CLink objects.

        //    Even though the best practice in modern C++ is to avoid explicit newand delete whenever possible, and instead use smart pointers like std::unique_ptr or std::shared_ptr that automatically manage the object's lifetime, our requirement is to use raw pointers and dynamic allocation. That's why we've clear ownership semantics and life-cycle management to ensure no memory leaks or undefined behavior occur.

        //    Here's how we've improved our design :

        //In this design, the CLink objects will be properly deallocated in the following function free_memory() when the Network object is destroyed.We are always mindful to remove any CLink objects from g_link_vector before we delete them elsewhere in our code to prevent dangling pointers.
    }

    void calculate_dynamic_VDFunction(int inner_iteration_number, bool congestion_bottleneck_sensitivity_analysis_mode);

    double get_generalized_first_order_gradient_cost_of_second_order_loss_for_mode_type( int mode_type_no)
    {
        // *60 as 60 min per hour
        double generalized_cost = link_avg_travel_time[mode_type_no] + penalty + toll[mode_type_no] / assignment.g_ModeTypeVector[mode_type_no].value_of_time * 60;

        if (assignment.g_ModeTypeVector[mode_type_no].eco_so_flag == 1)  // eo so users;
        {
            double speed_mph = this->link_distance_VDF / 1.6 / (max(0.001, link_avg_travel_time[mode_type_no] / 60));
            generalized_cost = link_avg_travel_time[mode_type_no];
        }
        //double emission_marginal = function(link_avg_travel_time[mode_type_no], this->link_distance_VDF);

        // system optimal mode or exterior panalty mode
        //if (assignment.assignment_mode == 4)
        //    generalized_cost += travel_marginal_cost[mode_type_no];

        return generalized_cost;
    }

    int main_node_id;

    int BWTT_in_simulation_interval;
    int zone_seq_no_for_outgoing_connector;

    //double number_of_lanes;

    double lane_capacity;

    std::map <int, int> m_link_pedefined_capacity_map_in_sec;  // per sec
    std::map <int, float> m_link_pedefined_information_response_map;  // per min, absolute time

    float model_speed[MAX_TIMEINTERVAL_PerDay];

    float est_volume_per_hour_per_lane[MAX_TIMEINTERVAL_PerDay];

    float est_avg_waiting_time_in_min[MAX_TIMEINTERVAL_PerDay]; // at link level
    float est_queue_length_per_lane[MAX_TIMEINTERVAL_PerDay];

    float get_model_5_min_speed(int time_in_min)
    {
        int t = time_in_min / 5;
        float total_speed_value = 0;
        int total_speed_count = 0;

        return model_speed[t];
    }

    float get_model_15_min_speed(int time_in_min)
    {
        int t = time_in_min / 5;
        float total_speed_value = 0;
        int total_speed_count = 0;

        for (int tt = 0; tt < 3; tt++)
        {

            if (t + tt >= 0 && t + tt < MAX_TIMEINTERVAL_PerDay)
            {
                if (model_speed[t + tt] >= 1)
                {
                    total_speed_value += model_speed[t + tt];
                    total_speed_count++;
                }
            }
        }

        return total_speed_value / max(1, total_speed_count);
    }

    float get_model_hourly_speed(int time_in_min)
    {
        int t = time_in_min / 5;
        float total_speed_value = 0;
        int total_speed_count = 0;

        for (int tt = 0; tt < 12; tt++)
        {

            if (t + tt >= 0 && t + tt < MAX_TIMEINTERVAL_PerDay)
            {
                if (model_speed[t + tt] >= 1)
                {
                    total_speed_value += model_speed[t + tt];
                    total_speed_count++;
                }
            }
        }

        return total_speed_value / max(1, total_speed_count);
    }

    float get_est_hourly_volume(int time_in_min)
    {
        int t = time_in_min / 5;
        float total_volume_value = 0;
        int total_volume_count = 0;

        for (int tt = 0; tt < 12; tt++)
        {

            if (t + tt >= 0 && t + tt < MAX_TIMEINTERVAL_PerDay)
            {
                if (est_volume_per_hour_per_lane[t + tt] >= 1)
                {
                    total_volume_value += est_volume_per_hour_per_lane[t + tt];
                    total_volume_count++;
                }
            }
        }

        return total_volume_value / max(1, total_volume_count);
    }

    int dynamic_link_event_start_time_in_min;
    std::map <int, bool> dynamic_link_closure_map;
    std::map <int, std::string> dynamic_link_closure_type_map;

    double length_in_meter;
    double link_distance_VDF;
    double link_distance_km;
    double link_distance_mile;
    double free_flow_travel_time_in_min;
    double total_simulated_delay_in_min;
    int total_simulated_meso_link_incoming_volume;
    double free_speed;

    double cost;
    double link_spatial_capacity;

    bool timing_arc_flag;
    int traffic_flow_code;
    int spatial_capacity_in_vehicles;
    int time_to_be_released;

    // 1. based on BPR.

    int link_seq_no;

    std::map<int, int> capacity_reduction_map;

    int global_minute_capacity_reduction_start;
    int global_minute_capacity_reduction_end;

    std::map<int, int> vms_map;

    std::string link_id;

    std::string geometry;

    int meso_link_id;
    int FT;
    int AT;
    std::string vdf_code;
    float PCE;

    float v_congestion_cutoff; // critical speed;
    float v_critical;

    float s3_m; // m factor in s3 model

    void update_kc(float free_speed_value)
    {
        v_critical = lane_capacity / k_critical;
        s3_m = 2 * log(2) / log(free_speed_value / v_critical);

        TMC_highest_speed = free_speed_value;

    }

    double get_volume_from_speed(float speed, float free_speed_value, float lane_capacity)
    {
        //test data free_speed = 55.0f;
        //speed = 52;
        //k_critical = 23.14167648;

        if (speed > free_speed_value * 0.99)
            speed = free_speed_value * 0.99;

        if (speed < 0)
            return -1;

        k_critical = 45;  // 45 vehicles per mile per lane based on HCM
        v_critical = lane_capacity / k_critical;
        s3_m = 2 * log(2) / log(free_speed_value / v_critical);

        double speed_ratio = free_speed_value / max(1.0f, speed);
        if (speed_ratio <= 1.00001)
            speed_ratio = 1.00001;

        /*   float volume = 0;*/
        double ratio_difference = pow(speed_ratio, s3_m / 2) - 1;

        double ratio_difference_final = max(ratio_difference, 0.00000001);

        double volume = speed * k_critical * pow(ratio_difference_final, 1 / s3_m);

        if (volume > lane_capacity)
            volume = lane_capacity;
        return volume;

    }

    bool AllowModeType(std::string mode_type)
    {
        if (allowed_uses.size() == 0 || allowed_uses.empty() || allowed_uses == "null" || allowed_uses == "all")  // if the allowed_uses is empty then all types are allowed.
            return true;
        else
        {
            if (allowed_uses.find(mode_type) != std::string::npos)  // otherwise, only an agent type is listed in this "allowed_uses", then this agent type is allowed to travel on this link
                return true;
            else
            {
                return false;
            }
        }
    }

    int from_node_seq_no;
    int to_node_seq_no;

    int from_node_id;
    int to_node_id;
    int AB_flag;  // 1 and -1;
    int BA_link_no;

    int link_type;
    double capacity;

    bool b_automated_generated_flag;

    int cell_type;  // 2 lane changing
    std::string mvmt_txt_id;
    std::string link_specifical_flag_str;
    std::string tmc_corridor_name;
    std::string link_type_name;

    float kjam;
    int type;

    //static
    //float flow_volume;
    //float travel_time;

    int subarea_id;

    double total_volume_for_all_mode_types;
    double total_person_volume_for_all_mode_types;

    double RT_flow_volume;
    double background_total_volume_for_all_mode_types;

    double  volume_per_mode_type[MAX_MODETYPES];
    double  converted_PCE_volume_per_at[MAX_MODETYPES];

    double  recorded_volume_per_at[MAX_MODETYPES];
    //double **recorded_lanes_per_at;
    double  recorded_PCE_per_at[MAX_MODETYPES];
    double recorded_DOC_per_at[MAX_MODETYPES];
    double recorded_TT_per_at[MAX_MODETYPES];
    double recorded_CO2_per_at[MAX_MODETYPES];
    double recorded_NOX_per_at[MAX_MODETYPES];

    double  queue_link_distance_VDF_perslot;  // # of vehicles in the vertical point queue
    double link_avg_travel_time[MAX_MODETYPES];
    double link_avg_co2_emit_per_mode[MAX_MODETYPES];
    double link_avg_nox_emit_per_mode[MAX_MODETYPES];

    double RT_waiting_time;

//    std::map<int, float> RT_travel_time_map;
    std::map<int, float> RT_speed_vector;
    //	double  travel_marginal_cost[MAX_MODETYPES];

    //TMC
    std::string tmc_code;
    int tmc_corridor_id;

    int tmc_road_order;

    int tmc_road_sequence;
    std::string tmc_road, tmc_direction, tmc_intersection;
    float tmc_reference_speed;
    float tmc_mean_speed;

    float tmc_volume;
    DTAGDPoint TMC_from, TMC_to;
    float TMC_highest_speed;

    //end of TMC

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

    int win_count;
    int lose_count;
};



class CNode
{
public:
    CNode() : zone_id{ -1 }, zone_org_id{ -1 }, layer_no{ 0 }, MRM_gate_flag{ -1 }, prohibited_movement_size{ 0 }, node_seq_no{ -1 }, subarea_id{ -1 }, is_activity_node{ 0 }, mode_type_no{ -1 }, is_boundary{ 0 }, access_distance{ 0.04 }
    {
    }

    //int accessible_node_count;

    int zone_id;
    __int64 cell_id;
    std::string cell_str;

    std::vector<CCoordinate> zone_coordinate_vector;

    // original zone id for non-centriod nodes
    int zone_org_id;
    float access_distance;
    std::string node_type;
    std::string mode_type_str;
    int subarea_id;
    int prohibited_movement_size;
    // sequence number
    int node_seq_no;
    int layer_no;
    int MRM_gate_flag;

    //external node number
    int node_id;

    int is_activity_node;
    int is_boundary;
    int mode_type_no;

    double x;
    double y;

    std::vector<int> m_outgoing_link_seq_no_vector;
    std::vector<int> m_incoming_link_seq_no_vector;

    std::vector<int> m_to_node_seq_no_vector;
    std::map<int, int> m_to_node_2_link_seq_no_map;

    // for simulation
    std::map<int, int> next_link_for_resource_request;

    std::map<int, float> label_cost;

    std::map<std::string, int> pred_per_iteration_map;
    std::map<std::string, float> label_cost_per_iteration_map;

    std::map<std::string, int> pred_RT_map;
    std::map<std::string, float> label_cost_RT_map;
};

class CInfoCell {
public:
    __int64 cell_id;
    std::string cell_str;
    std::vector<DTAGDPoint> m_ShapePoints;

    void CreateCell(double x, double y, double grid_resolution)
    {
        __int64 xi;
        xi = floor(x / grid_resolution);

        __int64 yi;
        yi = floor(y / grid_resolution);

        double left, right, top, bottom;

        left = xi * grid_resolution;
        right = (xi + 1) * grid_resolution;

        top = (yi + 1) * grid_resolution;
        bottom = (yi)*grid_resolution;

        DTAGDPoint	pt0, pt1, pt2, pt3, pt4;

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
    CTMC_Corridor_Info() : total_PMT{ 0 }, total_PHT{ 0 },
        total_PSDT {0 } ,
        lowest_speed{ 9999 },
        highest_speed{ -1 },
        link_count{ 0 }
    {
    }

    void reset()
    {
        total_PMT = 0;
        total_PHT = 0;
        total_PSDT = 0;
        lowest_speed = 9999;
        highest_speed = -1;
        link_count = 0;
    }

    double get_avg_speed()
    {
        return total_PMT / max(0.001, total_PHT);  //miles per hour
    }
    double total_PMT;
    double total_PHT;
    double total_PSDT;

    double avg_speed;
    double lowest_speed;
    double highest_speed;
    double total_congestion_duration;
    int link_count;

    std::map<int, int> road_sequence_map;
};


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
    COZone() : cell_id{ 0 }, gravity_production{ 10 }, gravity_attraction{ 10 }, cell_x{ 0 }, cell_y{ 0 },
        gravity_est_production{ 0 }, gravity_est_attraction{ 0 }, subarea_significance_flag{ true }, preread_total_O_demand{ 0 }, sindex{ -1 }, origin_zone_impact_volume{ 0 }, subarea_inside_flag { 3 },
        superzone_index{ -100 }, bcluster_seed{ false }, b_shortest_path_computing_flag{ true }, super_seed_zone_id{ -1 }, b_distrct_cluster_seed{ false }, analysis_district_index{ -100 }, preread_total_O_related_demand {0}
    {
    }

    std::map <int, double> preread_ODdemand;
    std::map <int, double> impact_passing_ODdemand;
    double preread_total_O_demand;
    double preread_total_O_related_demand;

    int sindex;  //subarea significance index
    int superzone_index;

    bool bcluster_seed;

    bool b_shortest_path_computing_flag;  // for super zones

    int  analysis_district_index;
    bool b_distrct_cluster_seed;

    bool subarea_significance_flag;
    double origin_zone_impact_volume;
    int subarea_inside_flag;
    std::map <int, bool> subarea_impact_flag;
    __int64 cell_id;
    std::string cell_code;
    double cell_x;
    double cell_y;

    float gravity_production;
    float gravity_attraction;

    float gravity_est_production;
    float gravity_est_attraction;

    // 0, 1,
    int zone_seq_no;
    // external zone id // this is origin zone
    int zone_id;
    int super_seed_zone_id;

    int node_seq_no;

    float obs_production_upper_bound_flag;
    float obs_attraction_upper_bound_flag;

    CODMatrix m_ODMatrix;
    CODMatrix m_ODAccessibilityMatrix;

    std::vector<int> m_activity_node_vector;
};


extern std::vector<COZone> g_zone_vector;
extern int g_related_zone_vector_size;  // by default, the same size as g_zone_vector, // if there is subarea, will be non-impacted zones., g_zone_vector has org_sindex;

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


//class for vehicle scheduling states
class CODState
{
public:

    float volume;
    float value;
    int orig;
    int dest;
    int tau;
    int at;

    int subarea_inside_flag_orig;
    int subarea_inside_flag_dest;

    void setup_input(int o, int d,  int a, int subarea_inside_flag_o, int subarea_inside_flag_d)
    {
        orig = o;
        dest = d;
        at = a;
        subarea_inside_flag_orig = subarea_inside_flag_o;
        subarea_inside_flag_dest = subarea_inside_flag_d;
    }

    void input_value(float val)
    {
        value = val;
    }

    bool operator<(const CODState& other) const
    {
        return value > other.value;
    }

};

extern std::vector<CNode> g_node_vector;
extern std::vector<CLink> g_link_vector;
//extern std::map<std::string, CVDF_Type> g_vdf_type_map;

extern std::map<int, DTAVehListPerTimeInterval> g_AgentTDListMap;
extern vector<CAgent_Simu*> g_agent_simu_vector;

extern std::map<std::string, CTMC_Corridor_Info> g_tmc_corridor_vector;
extern std::map<std::string, CInfoCell> g_info_cell_map;

void g_assign_RT_computing_tasks_to_memory_blocks(Assignment& assignment);
void g_load_demand_side_scenario_file(Assignment& assignment);
extern void g_add_new_virtual_connector_link(int internal_from_node_seq_no, int internal_to_node_seq_no, string mode_type_str, int zone_seq_no);
extern double g_Find_PPP_RelativeAngle(const DTAGDPoint* p1, const DTAGDPoint* p2, const DTAGDPoint* p3, const DTAGDPoint* p4);
extern double g_GetPoint2LineDistance(const DTAGDPoint* pt, const DTAGDPoint* FromPt, const DTAGDPoint* ToPt);
extern double g_GetPoint2Point_Distance(const DTAGDPoint* p1, const DTAGDPoint* p2);

#endif
