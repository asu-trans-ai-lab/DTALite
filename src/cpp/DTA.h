#ifndef GUARD_DTA_H
#define GUARD_DTA_H
#define BUILD_EXE //self-use

#include <algorithm>
#include <iomanip>

// Peiheng, 02/21/22, a temporary fix (bad practice)
using std::max;

constexpr auto MAX_LABEL_COST = 1.0e+15;
constexpr auto _INFO_ZONE_ID = 100000;
constexpr auto MAX_SCENARIOS = 10; //because of the od demand store format,the MAX_demandtype must >=g_DEMANDTYPES.size()+1;
constexpr auto MAX_AGNETTYPES = 10; //because of the od demand store format,the MAX_demandtype must >=g_DEMANDTYPES.size()+1;
constexpr auto MAX_LINKTTYPES = 100; 
constexpr auto MAX_TIMEPERIODS = 6; // time period set to 6: AM, MD, PM, LPM, SAT_MD
//constexpr auto MAX_AGNETTYPES = 10; //because of the od demand store format,the MAX_demandtype must >=g_DEMANDTYPES.size()+1;
//constexpr auto MAX_TIMEPERIODS = 6; // time period set to 4: mid night, morning peak, mid-day and afternoon peak;
//constexpr auto MAX_ORIGIN_DISTRICTS = 30; //origin based agreegration grids
//constexpr auto MAX_ORIGIN_DISTRICTS = 30; //origin based agreegration grids
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
enum e_VDF_type { q_vdf = 0, bpr_vdf };
enum e_assignment_mode { lue = 0, path_based_assignment= 1, simulation_dta=2};

// FILE* g_pFileOutputLog = nullptr;
extern void g_OutputModelFiles(int mode);
extern int g_related_zone_vector_size;

class CDemand_Period {
public:
    CDemand_Period() : demand_period_id { 0 }, demand_period { 0 }, starting_time_slot_no{ 0 }, ending_time_slot_no{ 0 }, t2_peak_in_hour{ 0 }, time_period_in_hour{ 1 }, number_of_demand_files{ 0 }
    {
    }

    int get_time_horizon_in_min()
    {
        return (ending_time_slot_no - starting_time_slot_no) * MIN_PER_TIMESLOT;
    }

    std::string demand_period;
    int starting_time_slot_no;
    int ending_time_slot_no;
    float time_period_in_hour;
    float t2_peak_in_hour;
    std::string time_period;
    int number_of_demand_files;
    int demand_period_id;
};

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

    void compute_cumulative_profile(int starting_slot_no, int ending_slot_no)
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
         
            dtalog.output() << std::setprecision(5) << "cumulative profile no. " << departure_time_profile_no << ", ratio at slot  " << s << " (" << hour << ":" << minute << ") = " << 
                departure_time_ratio[s] << ",CR " << 
                cumulative_departure_time_ratio[s] << std::endl;
        }

       dtalog.output() << std::setprecision(5) << "final cumulative profile ratio = " << cumulative_departure_time_ratio[ending_slot_no - 1] << std::endl;
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
//                dtalog.output() << "s=" << s <<" (" << hour << ":" << minute << ") = "  << ending_time_slot_no << std::endl;

                return s;
            }
        }
        int hour = starting_time_slot_no / 12;
        int minute = starting_time_slot_no * 5 - hour * 60;
//        dtalog.output() << "s=" << starting_time_slot_no << " (" << hour << ":" << minute << ") = " << ending_time_slot_no << std::endl;
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
            dtalog.output() << "s=" << s << " (" << hour << ":" << minute << ") = " << cumulative_departure_time_ratio[s] << std::endl;
            }
            if (r < cumulative_departure_time_ratio[s])
            {

                double slot_fraction = cumulative_departure_time_ratio[s] - cumulative_departure_time_ratio[s-1];
                double floating_point = max(0.0, (r - cumulative_departure_time_ratio[s - 1]) / max(0.00001, slot_fraction));

                double time_in_min = (s- starting_time_slot_no + floating_point )* MIN_PER_TIMESLOT;
                if (idebug)
                { 
                    dtalog.output() << "select: s=" << s << " (" << hour << ":" << minute << ") = " << ending_time_slot_no << ", dep_time = " << time_in_min <<"," << std::endl;
                }
                return time_in_min;
            }
        }

        if (idebug)
        {
            int hour = starting_time_slot_no / 12;
            int minute = starting_time_slot_no * 5 - hour * 60;

            dtalog.output() << "s=" << starting_time_slot_no << " (" << hour << ":" << minute << ") = " << ending_time_slot_no << std::endl;
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
        total_od_volume{ 0 }, total_person_distance_km{ 0 }, total_person_distance_mile{ 0 }, total_person_travel_time{ 0 }, avg_travel_time {0}, avg_travel_distance_km {0}, avg_travel_distance_mile{ 0 }
    {}

    int count;
    double total_od_volume;
    double total_person_distance_km;
    double total_person_distance_mile;
    double total_person_travel_time;
    double avg_travel_time;
    double avg_travel_distance_km;
    double avg_travel_distance_mile;
};

class CSystem_Summary
{
public:

    void record_mode_volume(int tau, int at, double od_volume)
    {
        if (at >= MAX_AGNETTYPES)
            return;

        data_by_demand_period_mode_type[tau][at].total_od_volume += od_volume;

    }

    void record_mode_od_data(CModeType_Summary element, int tau, int at)
    {
        if (at >= MAX_AGNETTYPES)
            return;

        data_by_demand_period_mode_type[tau][at].count += 1;
        data_by_demand_period_mode_type[tau][at].total_person_travel_time += element.total_person_travel_time;
        data_by_demand_period_mode_type[tau][at].total_person_distance_km += element.total_person_distance_km;
        data_by_demand_period_mode_type[tau][at].total_person_distance_mile += element.total_person_distance_mile;

    }

    void computer_avg_value(int tau, int at)
    {
        float count = data_by_demand_period_mode_type[tau][at].count;
        if (count >= 1)
        {
            data_by_demand_period_mode_type[tau][at].avg_travel_distance_km = data_by_demand_period_mode_type[tau][at].total_person_distance_km / max(1, data_by_demand_period_mode_type[tau][at].total_od_volume);
            data_by_demand_period_mode_type[tau][at].avg_travel_distance_mile = data_by_demand_period_mode_type[tau][at].total_person_distance_mile / max(1, data_by_demand_period_mode_type[tau][at].total_od_volume);
            data_by_demand_period_mode_type[tau][at].avg_travel_time = data_by_demand_period_mode_type[tau][at].total_person_travel_time / max(1, data_by_demand_period_mode_type[tau][at].total_od_volume);
        }
    }

    CModeType_Summary data_by_demand_period_mode_type[MAX_TIMEPERIODS][MAX_AGNETTYPES];
};


class CAnalysisDistrict
{
public:
    int district_id;
    int district_name;
    std::vector<DTAGDPoint> shape_points;

    void record_origin_2_district_volume(int at, double od_volume)
    {
        if (at >= MAX_AGNETTYPES)
            return;

        data_by_mode_type[at].total_od_volume += od_volume;

    }

    void record_link_2_district_data(CModeType_Summary element, int at)
    {
        if (at >= MAX_AGNETTYPES)
            return;

        data_by_mode_type[at].count += 1;
        data_by_mode_type[at].total_person_travel_time += element.total_person_travel_time;
        data_by_mode_type[at].total_person_distance_km += element.total_person_distance_km;
        data_by_mode_type[at].total_person_distance_mile += element.total_person_distance_mile;

    }

    void computer_avg_value(int at)
    {
        float count = data_by_mode_type[at].count;
        if (count >= 1)
        {
            data_by_mode_type[at].avg_travel_distance_km = data_by_mode_type[at].total_person_distance_km/ max(1,data_by_mode_type[at].total_od_volume);
            data_by_mode_type[at].avg_travel_distance_mile = data_by_mode_type[at].total_person_distance_mile / max(1, data_by_mode_type[at].total_od_volume);
            data_by_mode_type[at].avg_travel_time = data_by_mode_type[at].total_person_travel_time / max(1, data_by_mode_type[at].total_od_volume);
        }
    }

    CModeType_Summary data_by_mode_type[MAX_AGNETTYPES];
};

class Cmode_type {
public:
    Cmode_type() : mode_type_no{ 1 }, value_of_time{ 100 }, time_headway_in_sec{ 1 }, real_time_information{ 0 }, access_speed{ 2 }, access_distance_lb{ 0.0001 }, access_distance_ub{ 4 }, acecss_link_k{ 4 },
         OCC{ 1 }, DSR{ 1 }, number_of_allowed_links{ 0 }, mode_specific_assignment_flag{ 0 }
    {
    }

    int mode_specific_assignment_flag; 

    int mode_type_no;
    // dollar per hour
    float value_of_time;
    // link type, product consumption equivalent used, for travel time calculation
    double OCC;
    double DSR;

    float time_headway_in_sec;
    int real_time_information;
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
    CLinkType() : link_type{ 1 }, number_of_links{ 0 }, traffic_flow_code{ spatial_queue }, k_jam{ 300 }, vdf_type{ q_vdf }
    {
        for (int at = 0; at < MAX_AGNETTYPES; at++)
        {
            free_speed_at[at] = 60;
            capacity_at[at] = 2000;
            FFTT_at[at] = 1;
            lanes_at[at] = 1;

            for (int at2 = 0; at2 < MAX_AGNETTYPES; at2++)
            {
                meu_matrix[at][at2] = 1;
            }
        }


        for (int tau = 0; tau < MAX_TIMEPERIODS; tau++)
        {
            for (int at = 0; at < MAX_AGNETTYPES; at++)
            { 
              peak_load_factor_period_at[tau][at] = 1;
            }
        }
        
    }

    int link_type;
    int number_of_links;
    float k_jam;
    std::string link_type_name;
    std::string type_code;
    e_VDF_type vdf_type;
    e_traffic_flow_model traffic_flow_code;

    double free_speed_at[MAX_AGNETTYPES];
    double lanes_at[MAX_AGNETTYPES];
    double capacity_at[MAX_AGNETTYPES];
    double FFTT_at[MAX_AGNETTYPES];
    std::string allow_uses_period[MAX_TIMEPERIODS];
    double peak_load_factor_period_at[MAX_TIMEPERIODS][MAX_AGNETTYPES];
    double meu_matrix[MAX_AGNETTYPES][MAX_AGNETTYPES];



};

class CColumnPath {
public:
    CColumnPath() : path_node_vector{ nullptr }, path_link_vector{ nullptr }, path_seq_no{ 0 }, m_link_size{ 0 }, m_node_size{ 0 },
        path_switch_volume{ 0 }, path_volume{ 0 }, path_preload_volume{ 0 }, path_volume_before_ODME{ -1 }, path_volume_after_ODME{ -1 }, path_volume_before_sa{ -1 }, path_volume_after_sa{ -1 }, path_travel_time{ 0 }, path_distance{ 0 }, path_toll{ 0 }, UE_gap{ 0 }, UE_relative_gap{ 0 },
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
        path_node_vector = new int[node_size];
        path_link_vector = new int[link_size];

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
            std::cout << "error: m_link_size == 0 in function CColumnPath::AllocateVector()!";
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
    int path_seq_no;
    int global_path_no;
    std::string path_id;
    // path volume
    double path_volume = 0;
    double path_preload_volume = 0;
    double path_volume_before_ODME = 0;
    double path_volume_after_ODME = 0;
    double path_volume_before_sa = 0;
    double path_volume_after_sa = 0;

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
public:
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
    int active_scenario_index;
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

public:


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
    CColumnVector() : avg_cost{ 0 }, avg_time{ 0 }, distance{ 0 },  prev_od_volume{ 0 }, bfixed_route{ false }, m_passing_sensor_flag{ -1 }, information_type{ 0 }, activity_mode_type_no{ 0 },
        departure_time_profile_no{ -1 }, OD_impact_flag{ 0 }, subarea_passing_flag{ 1 }, relative_OD_gap{ 0 }
    {

        for (int si = 0; si < MAX_SCENARIOS; si++)
            od_volume[si] = 0;
    }


    void reset_column_pool()
    {
        path_node_sequence_map.clear();

        avg_cost = 0;
        avg_time = 0;
        distance = 0;
        relative_OD_gap = 0;

    }
    bool subarea_passing_flag;
    
    std::map<int, bool> at_od_impacted_flag_map; // for each agent type

    float avg_cost;
    float avg_time;
    float distance;
    // od volume
    double od_volume[MAX_SCENARIOS];

    double relative_OD_gap; 
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
    std::string demand_period;

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

    DTAScenario() : scenario_index{ 0 }, year{ 2023 }, activation{ 1 }
    {
    }

    int scenario_index;
    int year;
    std::string scenario_name;
    std::string scenario_description;
    int activation;
};
class Assignment {
public:
    // default is UE
    Assignment() : assignment_mode{ lue }, g_number_of_memory_blocks{ 4 }, g_number_of_threads{ 1 }, g_info_updating_freq_in_min{ 5 }, g_visual_distance_in_cells{ 5 },
        g_link_type_file_loaded{ true }, g_mode_type_file_loaded{ false }, total_route_demand_volume{ 0 },
        total_demand_volume{ 0.0 }, g_column_pool{ nullptr }, g_number_of_in_memory_simulation_intervals{ 500 },
        g_number_of_column_generation_iterations{ 20 }, g_number_of_column_updating_iterations{ 0 }, g_number_of_ODME_iterations{ 0 }, g_number_of_sensitivity_analysis_iterations{ -1 }, g_number_of_demand_periods{ 24 }, g_number_of_links{ 0 }, g_number_of_timing_arcs{ 0 },
        g_number_of_nodes{ 0 }, g_number_of_zones{ 0 }, g_number_of_mode_types{ 0 }, debug_detail_flag{ 1 }, path_output{ 1 }, trajectory_output_count{ -1 },
        trace_output{ 0 }, major_path_volume_threshold{ 0.1 }, trajectory_sampling_rate{ 1.0 }, td_link_performance_sampling_interval_in_min{ -1 }, dynamic_link_performance_sampling_interval_hd_in_min{ 15 }, trajectory_diversion_only{ 0 }, m_GridResolution{ 0.01 },
        shortest_path_log_zone_id{ -1 }, g_number_of_analysis_districts{ 1 },
        active_scenario_index{ 0 }

    {
        m_LinkCumulativeArrivalVector  = NULL;
        m_LinkCumulativeDepartureVector = NULL;

        m_link_CA_count = NULL;  // CA, assign this value to m_LinkCumulativeArrivalVector at a given time in min
        m_link_CD_count = NULL;
        m_LinkOutFlowCapacity = NULL;
        m_LinkOutFlowState =  NULL;

        sp_log_file.open("log_label_correcting.txt");
        
        summary_file.open("final_summary.csv", std::fstream::out);
        if (summary_file &&!summary_file.is_open())
        {
            dtalog.output() << "File final_summary.csv cannot be open.";
            g_program_stop();
        }
        simu_log_file.open("log_simulation.txt");
//        simu_log_file << "start" << std::endl;
        g_rt_network_pool = NULL;
        g_column_pool = NULL;
    }

    ~Assignment()
    {
        if (g_column_pool)
            Deallocate4DDynamicArray(g_column_pool, g_related_zone_vector_size, g_related_zone_vector_size, g_number_of_mode_types);
        
        if(g_rt_network_pool)
            Deallocate3DDynamicArray(g_rt_network_pool, g_number_of_zones, g_number_of_mode_types);

        sp_log_file.close();
        
        summary_file.close();
        summary_file2.close();
        summary_corridor_file.close();
        summary_system_file.close();
        simu_log_file.close();
        DeallocateLinkMemory4Simulation();
    }

    void InitializeDemandMatrix(int number_of_signficant_zones, int number_of_zones, int number_of_mode_types, int number_of_time_periods)
    {
        total_demand_volume = 0.0;
        g_number_of_zones = number_of_zones;
        g_number_of_mode_types = number_of_mode_types;

        g_column_pool = Allocate4DDynamicArray<CColumnVector>(number_of_signficant_zones, g_related_zone_vector_size, max(1, number_of_mode_types), number_of_time_periods);

        for (int i = 0; i < number_of_zones; ++i)
        {
            g_origin_demand_array[i] = 0.0;
        }

        for (int i = 0; i < number_of_mode_types; ++i)
        {
            for (int tau = 0; tau < g_number_of_demand_periods; ++tau)
            { 
                total_demand[i][tau] = 0.0;
                total_route_demand[i][tau] = 0.0;
            }
        }

        g_DemandGlobalMultiplier = 1.0f;
    }

    int get_in_memory_time(int t)
    {
        return t % g_number_of_in_memory_simulation_intervals;
    }


    std::vector<DTAScenario> g_DTAscenario_vector;
    std::map<int, int> g_active_DTAscenario_map;
    
    std::vector<DTAGDPoint> g_subarea_shape_points;
    std::vector<DTAGDPoint> g_MRM_subarea_shape_points;
    

    std::map<int, int> g_node_id_to_MRM_subarea_mapping;  // this is an one-to-one mapping: 1: outside 1: inside 

    std::map<int, int> g_micro_node_id_to_MRM_bridge_mapping;  // this is an one-to-one mapping: for MRM bridge

    std::map<int, int> g_macro_node_id_to_MRM_incoming_bridge_mapping;  // this is an one-to-one mapping: for MRM bridge
    std::map<int, int> g_macro_node_id_to_MRM_outgoing_bridge_mapping;  // this is an one-to-one mapping: for MRM bridge





    void STTrafficSimulation();
    void STMesoTrafficSimulation();

    //OD demand estimation estimation
    void GenerateDefaultMeasurementData();
    void Demand_ODME(int OD_updating_iterations);
    void Sensor_Vector_based_Demand_ODME(int OD_updating_iterations);
    void AllocateLinkMemory4Simulation();
    int update_real_time_info_path(CAgent_Simu* p_agent, int& impacted_flag_change, float updating_in_min);
    bool RTSP_real_time_travel_time_updating(int time_slot_no, int simu_interval_t);
    void DeallocateLinkMemory4Simulation();

    std::map<int, int> zone_id_to_centriod_node_no_mapping;  // this is an one-to-one mapping
    std::map<int, int> zone_id_2_node_no_mapping;  // this is used to mark if this zone_id has been identified or not
    std::map<int, __int64> zone_id_2_cell_id_mapping;  // this is used to mark if this zone_id has been identified or not
    std::map<__int64, int> cell_id_mapping;  // this is used to mark if this cell_id has been identified or not
    std::map<__int64, std::string> cell_id_2_cell_code_mapping;  // this is used to mark if this cell_id has been identified or not


    double m_GridResolution;
    e_assignment_mode assignment_mode;

    int active_scenario_index;
    int g_number_of_memory_blocks;
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
    float total_route_demand_volume;
    std::map<int, float> g_origin_demand_array;
    CColumnVector**** g_column_pool;
    NetworkForSP* *** g_rt_network_pool;

    // the data horizon in the memory
    int g_number_of_in_memory_simulation_intervals;
    int g_number_of_column_generation_iterations;
    int g_number_of_sensitivity_analysis_iterations;
    int g_number_of_column_updating_iterations;
    int g_number_of_ODME_iterations;
    int g_number_of_demand_periods;

    int g_number_of_links;
    int g_number_of_timing_arcs;
    int g_number_of_nodes;
    int g_number_of_zones;
    int g_number_of_mode_types;

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

    std::vector<CDemand_Period> g_DemandPeriodVector;
    std::vector<CDeparture_time_Profile> g_DepartureTimeProfileVector;

    int g_LoadingStartTimeInMin;
    int g_LoadingEndTimeInMin;

    std::vector<Cmode_type> g_ModeTypeVector;


    std::map<std::string, CActivityTravelPattern> g_ActivityTravelPatternMap;
    std::map<std::string, CChoiceSet> g_ChoiceSetMap;
    

    int g_number_of_analysis_districts;
    std::map<int, CLinkType> g_LinkTypeMap;

    std::map<std::string, int> demand_period_to_seqno_mapping;
    std::map<std::string, int> mode_type_2_seqno_mapping;


    std::map<int, double> o_district_id_factor_map;
    std::map<int, double> d_district_id_factor_map;
    std::map<int, double> od_district_id_factor_map;

    std::map<int, double> SA_o_district_id_factor_map;
    std::map<int, double> SA_d_district_id_factor_map;
    std::map<int, double> SA_od_district_id_factor_map;



    float total_demand[MAX_AGNETTYPES][MAX_TIMEPERIODS];
    float total_route_demand[MAX_AGNETTYPES][MAX_TIMEPERIODS];
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
    std::ofstream summary_file;
    std::ofstream summary_file2;
    std::ofstream summary_corridor_file;
    std::ofstream summary_system_file;

};

extern Assignment assignment;

# include "VDF.h"

class CLink
{
public:
    // construction
    CLink() :main_node_id{ -1 }, free_speed{ 100 }, v_congestion_cutoff{ 100 }, v_critical { 60 },
        length_in_meter{ 1 }, link_distance_VDF {0.001}, 
        BWTT_in_simulation_interval{ 100 }, zone_seq_no_for_outgoing_connector{ -1 }, lane_capacity{ 1999 },
         free_flow_travel_time_in_min{ 0.01 }, link_spatial_capacity{ 100 }, 
        timing_arc_flag{ false }, traffic_flow_code{ 0 }, spatial_capacity_in_vehicles{ 999999 }, subarea_id{ -1 }, RT_flow_volume{ 0 },
        cell_type{ -1 }, saturation_flow_rate{ 1800 }, dynamic_link_event_start_time_in_min{ 99999 }, b_automated_generated_flag{ false }, time_to_be_released{ -1 },
        RT_waiting_time{ 0 }, FT{ 1 }, AT{ 1 }, s3_m{ 4 }, tmc_road_order{ 0 }, tmc_road_sequence{ -1 }, k_critical{ 45 }, vdf_type{ q_vdf }, 
        tmc_corridor_id{ -1 }, from_node_id{ -1 }, to_node_id{ -1 }, kjam{ 300 }, link_distance_km{ 0 }, link_distance_mile{ 0 }, meso_link_id{ -1 }, total_simulated_delay_in_min{ 0 }, 
        total_simulated_meso_link_incoming_volume{ 0 }, global_minute_capacity_reduction_start{ -1 }, global_minute_capacity_reduction_end{ -1 },
        layer_no{ 0 }, AB_flag {1}, BA_link_no {-1}
   {
        for (int si = 0; si < MAX_SCENARIOS; si++)
        {
            link_type_si[si] = 0;
            number_of_lanes_si[si] = 1;  // default all open 
            penalty_si_flag[si] = 0;

            for (int tau = 0; tau < MAX_TIMEPERIODS; ++tau)
                for (int at = 0; at < MAX_AGNETTYPES; ++at)
                {

                    penalty_si_at[si][at][tau] = 0;
                    recorded_volume_per_mode_type_per_period[tau][at][si] = 0;
                    recorded_converted_MEU_volume_per_period_per_at[tau][at][si] = 0;
                }

        }


        for (int tau = 0; tau < MAX_TIMEPERIODS; ++tau)
        {
            total_volume_for_all_mode_types_per_period[tau] = 0;
            total_person_volume_for_all_mode_types_per_period[tau] = 0;
            queue_link_distance_VDF_perslot[tau] = 0;
                       //cost_perhour[tau] = 0;
            for (int at = 0; at < MAX_AGNETTYPES; ++at)
            {
                volume_per_mode_type_per_period[tau][at] = 0; 
                converted_MEU_volume_per_period_per_at[tau][at] = 0; 
                travel_time_per_period[tau][at] = 0;

            }
           
        }

    }

    ~CLink()
    {
    }

    // Peiheng, 02/05/21, useless block
    void free_memory()
    {
    }

    void calculate_dynamic_VDFunction(int inner_iteration_number, bool congestion_bottleneck_sensitivity_analysis_mode, int vdf_type);

    void calculate_marginal_cost_for_mode_type(int tau, int mode_type_no, float PCE_mode_type)
    {
        // volume * dervative
        // BPR_term: volume * FFTT * alpha * (beta) * power(v/c, beta-1),

//        travel_marginal_cost_per_period[tau][mode_type_no] = VDF_period[tau].marginal_base * PCE_mode_type;
    }

    double get_generalized_first_order_gradient_cost_of_second_order_loss_for_mode_type(int tau, int mode_type_no)
    {
        // *60 as 60 min per hour
        double generalized_cost = travel_time_per_period[tau][mode_type_no] + VDF_period[tau].penalty + VDF_period[tau].toll[mode_type_no] / assignment.g_ModeTypeVector[mode_type_no].value_of_time * 60;

        // system optimal mode or exterior panalty mode
        //if (assignment.assignment_mode == 4)
        //    generalized_cost += travel_marginal_cost_per_period[tau][mode_type_no];

        return generalized_cost;
    }

    int main_node_id;


    int BWTT_in_simulation_interval;
    int zone_seq_no_for_outgoing_connector;

    double number_of_lanes_si[MAX_SCENARIOS];
    double penalty_si_at[MAX_SCENARIOS][MAX_AGNETTYPES][MAX_TIMEPERIODS];
    double penalty_si_flag[MAX_SCENARIOS];

    double lane_capacity;
    double saturation_flow_rate;

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
    float k_critical; // critical density;
    float s3_m; // m factor in s3 model

    void update_kc(float free_speed_value)
    {
        k_critical = 45;  // 45 vehicles per mile per lane based on HCM
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

    bool AllowModeType(std::string mode_type, int tau, int active_si)
    {
        if (VDF_period[tau].allowed_uses[active_si].size() == 0 || VDF_period[tau].allowed_uses[active_si] == "all")  // if the allowed_uses is empty then all types are allowed.
            return true;
        else
        {
            if (VDF_period[tau].allowed_uses[active_si].find(mode_type) != std::string::npos)  // otherwise, only an agent type is listed in this "allowed_uses", then this agent type is allowed to travel on this link
                return true;
            else
            {
                return false;
            }


        }
    }


    bool SA_AllowModeType(std::string mode_type, int tau)
    {
        if (VDF_period[tau].sa_allowed_uses.size() == 0 || VDF_period[tau].sa_allowed_uses == "all")  // if the allowed_uses is empty then all types are allowed.
            return true;
        else
        {
            if (VDF_period[tau].sa_allowed_uses.find(mode_type) != std::string::npos)  // otherwise, only an agent type is listed in this "sa_allowed_uses", then this agent type is allowed to travel on this link
                return true;
            else
            {
                return false;
            }


        }
    }

    int from_node_seq_no;
    int to_node_seq_no;
    int layer_no;
    int from_node_id;
    int to_node_id;
    int AB_flag;  // 1 and -1; 
    int BA_link_no;

    int link_type_si[MAX_SCENARIOS];

    bool b_automated_generated_flag;

    int cell_type;  // 2 lane changing
    std::string mvmt_txt_id;
    std::string link_code_str;
    std::string tmc_corridor_name;
    std::string link_type_name;
    std::string link_type_code;

    e_VDF_type    vdf_type;
    float kjam;

    CPeriod_VDF VDF_period[MAX_TIMEPERIODS];

    int type;

    //static
    //float flow_volume;
    //float travel_time;

    int subarea_id;
    
    double total_volume_for_all_mode_types_per_period[MAX_TIMEPERIODS];
    double total_person_volume_for_all_mode_types_per_period[MAX_TIMEPERIODS];

    double RT_flow_volume;
    double background_total_volume_for_all_mode_types_per_period[MAX_TIMEPERIODS];

    double  volume_per_mode_type_per_period[MAX_TIMEPERIODS][MAX_AGNETTYPES];
    double  converted_MEU_volume_per_period_per_at[MAX_TIMEPERIODS][MAX_AGNETTYPES];

    double  recorded_volume_per_mode_type_per_period[MAX_TIMEPERIODS][MAX_AGNETTYPES][MAX_SCENARIOS];
    double  recorded_converted_MEU_volume_per_period_per_at[MAX_TIMEPERIODS][MAX_AGNETTYPES][MAX_SCENARIOS];


    

    double  queue_link_distance_VDF_perslot[MAX_TIMEPERIODS];  // # of vehicles in the vertical point queue
    double travel_time_per_period[MAX_TIMEPERIODS][MAX_AGNETTYPES];
    double RT_waiting_time;

//    std::map<int, float> RT_travel_time_map;
    std::map<int, float> RT_speed_vector;
    //	double  travel_marginal_cost_per_period[MAX_TIMEPERIODS][MAX_AGNETTYPES];

    int number_of_periods;


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


class CVDF_Type
{
public:
    CVDF_Type() {}

    void record_qvdf_data(CPeriod_VDF element, int tau)
    {
        if (tau >= MAX_TIMEPERIODS)
            return;

        if (VDF_period_sum[tau].vdf_data_count == 0)
        {
            VDF_period_sum[tau].Q_alpha = element.Q_alpha;
            VDF_period_sum[tau].Q_beta = element.Q_beta;
            VDF_period_sum[tau].Q_cp = element.Q_cp;
            VDF_period_sum[tau].Q_n = element.Q_n;
            VDF_period_sum[tau].Q_s = element.Q_s;
            VDF_period_sum[tau].Q_cd = element.Q_cd;
        }
        else
        {
            VDF_period_sum[tau].Q_alpha +=  element.Q_alpha;
            VDF_period_sum[tau].Q_beta +=  element.Q_beta;
            VDF_period_sum[tau].Q_cp += element.Q_cp;
            VDF_period_sum[tau].Q_n +=  element.Q_n;
            VDF_period_sum[tau].Q_s += element.Q_s;
            VDF_period_sum[tau].Q_cd += element.Q_cd;

        }

        VDF_period_sum[tau].vdf_data_count++;
    }

    void computer_avg_parameter(int tau)
    {
        float count = VDF_period_sum[tau].vdf_data_count;
        if(count>=1)
        {
        VDF_period_sum[tau].Q_alpha /= count;
        VDF_period_sum[tau].Q_beta /= count;
        VDF_period_sum[tau].Q_cp /= count;
        VDF_period_sum[tau].Q_n /= count;
        VDF_period_sum[tau].Q_s /= count;
        VDF_period_sum[tau].Q_cd /= count;
        }
    }

    std::string vdf_code;
    CPeriod_VDF VDF_period_sum[MAX_TIMEPERIODS];
};

class CPeriod_Corridor
{
public:
    CPeriod_Corridor() :volume{ 0 }, count{ 0 }, speed{ 0 }, DoC{ 0 }, P{ 0 },  AvgP{ 0 }, MaxP{0}
    {}

    int count;
    double volume, speed, DoC, D, P, AvgP, MaxP;


};


class CCorridorInfo
{
public:
    CCorridorInfo() {}

    void record_link_2_corridor_data(CPeriod_Corridor element, int tau)
    {
        if (tau >= MAX_TIMEPERIODS)
            return;

        corridor_period[tau].volume += element.volume;
        corridor_period[tau].DoC += element.DoC; 
        corridor_period[tau].speed += element.speed;
        corridor_period[tau].P = max(corridor_period[tau].P, element.P);
        corridor_period[tau].count += 1;

    }

    void computer_avg_value(int tau)
    {
        float count = corridor_period[tau].count;
        if (count >= 1)
        {
            corridor_period[tau].volume /= count;
            corridor_period[tau].speed /= count;;
            corridor_period[tau].DoC /= count;;
        }
    }

    std::string tmc_corridor_name;
    CPeriod_Corridor corridor_period[MAX_TIMEPERIODS];
    CPeriod_Corridor corridor_period_before[MAX_TIMEPERIODS];
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

    std::map<std::string, int> m_prohibited_movement_string_map;
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
    COZone() : cell_id{ 0 }, obs_production { 0 }, obs_attraction{ 0 },
        est_production{ 0 }, est_attraction{ 0 },
        est_production_dev{ 0 }, est_attraction_dev{ 0 }, gravity_production{ 100 }, gravity_attraction{ 100 }, cell_x{ 0 }, cell_y{ 0 },
        gravity_est_production{ 0 }, gravity_est_attraction{ 0 }, subarea_significance_flag{ true }, preread_total_O_demand{ 0 }, sindex{ -1 }, origin_zone_impact_volume{ 0 }, subarea_inside_flag { 3 },
        superzone_index{ -100 }, bcluster_seed{ false }, b_shortest_path_computing_flag{ true }, super_seed_zone_id{ -1 }, b_distrct_cluster_seed{ false }, distrct_cluster_index{ -100 }, preread_total_O_related_demand {0}
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

    int  distrct_cluster_index;
    bool b_distrct_cluster_seed;

    bool subarea_significance_flag;
    double origin_zone_impact_volume;
    int subarea_inside_flag;
    std::map <int, bool> subarea_impact_flag;
    __int64 cell_id;
    std::string cell_code;
    double cell_x;
    double cell_y;

    float obs_production;
    float obs_attraction;

    float gravity_production;
    float gravity_attraction;

    float gravity_est_production;
    float gravity_est_attraction;

    float est_production;
    float est_attraction;

    float est_production_dev;
    float est_attraction_dev;

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


    void setup_input(int o, int d,  int a, int t, int subarea_inside_flag_o, int subarea_inside_flag_d)
    {
        orig = o;
        dest = d;
        tau = t;
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
extern std::map<std::string, CVDF_Type> g_vdf_type_map;
extern std::map<std::string, CCorridorInfo> g_corridor_info_SA_map;

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
