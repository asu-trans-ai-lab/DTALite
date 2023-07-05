/* Portions Copyright 2019-2021 Xuesong Zhou and Peiheng Li, Cafer Avci

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

#include "config.h"
#include "utils.h"

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

using std::max;
using std::min;
using std::string;
using std::vector;
using std::map;
using std::ifstream;
using std::ofstream;
using std::istringstream;


__int64 g_get_cell_ID(double x, double y, double grid_resolution)
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

string g_get_cell_code(double x, double y, double grid_resolution, double left, double top)
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

#include "DTA.h"
#include "routing.h"


// some basic parameters setting


std::vector<NetworkForSP*> g_NetworkForSP_vector;
std::vector<NetworkForSP*> g_NetworkForRTSP_vector;
NetworkForSP g_RoutingNetwork;
std::map<int, DTAVehListPerTimeInterval> g_AgentTDListMap;
vector<CAgent_Simu*> g_agent_simu_vector;
std::vector<CNode> g_node_vector;
std::vector<CLink> g_link_vector;
std::map<string, CVDF_Type> g_vdf_type_map;
std::map<int, CSystem_Summary>  g_district_summary_map;
std::map<int, CSystem_Summary>  g_scenario_summary_map;
std::vector<COZone> g_zone_vector;
int g_related_zone_vector_size;
int g_number_of_active_scenarios = 1;
int g_number_of_active_mode_types = 1;
int g_number_of_active_demand_perioids =1;


Assignment assignment;

std::map<string, CTMC_Corridor_Info> g_tmc_corridor_vector;
std::map<string, CInfoCell> g_info_cell_map;

#include "trip_generation.h"
#include "input.h"
#include "output.h"
#include "assignment.h"
#include "ODME.h"
#include "scenario_API.h"


void g_column_regeneration(Assignment& assignment, bool real_time_info_flag)  // for RT
{		//step 3: column generation stage: find shortest path and update path cost of tree using TD travel time
	clock_t start_t, end_t, iteration_t;
	start_t = clock();
	
	dtalog.output() << "[PROCESS INFO] Step 5: Column Re-Generation for Traffic Assignment..." << '\n';
	dtalog.output() << "[DATA INFO] Total Column Re-Generation iteration: " << assignment.g_number_of_sensitivity_analysis_iterations_for_dtm << '\n';
	// stage II: column re-generation at the sensitivity analysis stage

	// for K = 5 iterations, which should be sufficient


	double total_system_travel_time = 0;
	double total_least_system_travel_time = 0;
	// initialization at beginning of shortest path
	double total_distance = 0;
	total_system_travel_time = update_link_travel_time_and_cost(0, total_distance);

	int number_of_column_regeneration_iterations = 1;
	for (int iteration_number = 0; iteration_number < number_of_column_regeneration_iterations; iteration_number++)
	{


		
		dtalog.output() << "[DATA INFO] Regenerating column for real-time information provision - Iteration Number: " << iteration_number << '\n';
		end_t = clock();
		iteration_t = end_t - start_t;
		dtalog.output() << ", CPU time: " << iteration_t / 1000.0 << " s" << '\n';

		// step 3.1 update travel time and resource consumption
		clock_t start_t_lu = clock();

		//////detect induced_delay
		for (int i = 0; i < g_link_vector.size(); ++i)
		{
			for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
			{
				// used in travel time calculation
				//if (g_link_vector[i].VDF_period[tau].dynamic_traffic_management_flag == 0)
				//{
				//	if (g_link_vector[i].VDF_period[tau].DOC > 3)  // dramatic delay, 3 is a cut off ratio/threshold for us to define large delay
				//	{
				//		g_link_vector[i].VDF_period[tau].dynamic_traffic_management_flag = -2; // induced delay
				//	}

				//}

			}
		}


			//// we can have a recursive formulat to reupdate the current link volume by a factor of k/(k+1),f
			////  and use the newly generated path flow to add the additional 1/(k+1)
			//bool sensitivity_analaysis_flag = true;
			//g_reset_and_update_link_volume_based_on_columns(g_link_vector.size(), iteration_number, true, sensitivity_analaysis_flag);

		if (dtalog.debug_level() >= 3)
		{
			dtalog.output() << "[DATA INFO] Results:" << '\n';
			for (int i = 0; i < g_link_vector.size(); ++i) {
				dtalog.output() << "[DATA INFO] link: " << g_node_vector[g_link_vector[i].from_node_seq_no].node_id << "-->"
					<< g_node_vector[g_link_vector[i].to_node_seq_no].node_id << ", "
					<< "flow count:" << g_link_vector[i].total_volume_for_all_mode_types_per_period[0] << '\n';
			}
		}

		end_t = clock();
		iteration_t = end_t - start_t_lu;
		// g_fout << "Link update with CPU time " << iteration_t / 1000.0 << " s; " << (end_t - start_t) / 1000.0 << " s" << '\n';

		//****************************************//
		//step 3.2 computng block for continuous variables;

		clock_t start_t_lc = clock();
		clock_t start_t_cp = clock();
		clock_t iteration_t, cumulative_lc, cumulative_cp, cumulative_lu;
		cumulative_lc = 0;
		cumulative_cp = 0;
		cumulative_lu = 0;

		int number_of_memory_blocks = min((int)g_NetworkForSP_vector.size(), assignment.g_number_of_memory_blocks);

#pragma omp parallel for  // step 3: C++ open mp automatically create n threads., each thread has its own computing thread on a cpu core
		//for (int ProcessID = 0; ProcessID < g_NetworkForSP_vector.size(); ++ProcessID)
		//{
		//    int mode_type_no = g_NetworkForSP_vector[ProcessID]->m_mode_type_no;

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
			for (int blk = 0; blk < assignment.g_ModeTypeVector.size() * assignment.g_DemandPeriodVector.size(); ++blk)
			{
				int network_copy_no = blk * assignment.g_number_of_memory_blocks + ProcessID;
				if (network_copy_no >= g_NetworkForSP_vector.size())
					continue;

				NetworkForSP* pNetwork = g_NetworkForSP_vector[network_copy_no];


				for (int o_node_index = 0; o_node_index < pNetwork->m_origin_node_vector.size(); ++o_node_index)
				{
					start_t_lc = clock();


					int origin_zone_seq_no = pNetwork->m_origin_zone_seq_no_vector[o_node_index]; // assigned nodes for computing

				if(real_time_info_flag == true)
				{
					if (assignment.g_ModeTypeVector[pNetwork->m_mode_type_no].real_time_information_type == 2 /*DMS*/
						&& assignment.zone_seq_no_2_info_mapping.find(origin_zone_seq_no) != assignment.zone_seq_no_2_info_mapping.end())
					{
						int idebug = 1;
					}

					if (assignment.g_ModeTypeVector[pNetwork->m_mode_type_no].real_time_information_type == 0)  // for non-information users, do not generate new path "columns"
						continue;  // skip the following shortest path codes

					if (assignment.g_ModeTypeVector[pNetwork->m_mode_type_no].real_time_information_type == 2 /*DMS infomation with physcial origin*/
						&& assignment.zone_seq_no_2_info_mapping.find(origin_zone_seq_no) == assignment.zone_seq_no_2_info_mapping.end()
						)
					{
						continue;  // skip the shortest path regeneration in the SA stage
					}

				}
				else
				{  // offline case,
					// all mode types will recompute their paths based on allowed_uses and updated number of lanes
				}

				pNetwork->optimal_label_correcting(ProcessID, &assignment, iteration_number, o_node_index, false, false, true, real_time_info_flag);



					end_t = clock();
					cumulative_lc += end_t - start_t_lc;

					start_t_cp = clock();

					bool sensitivity_analaysis_flag = true;
					double total_origin_least_travel_time = pNetwork->backtrace_RT_shortest_path_tree(assignment, number_of_column_regeneration_iterations, o_node_index, sensitivity_analaysis_flag);


#pragma omp critical
					{
						total_least_system_travel_time += total_origin_least_travel_time;
					}
					end_t = clock();
					cumulative_cp += end_t - start_t_cp;
				}
			}

		}


		if (iteration_number == 0)
		{
			assignment.summary_file << "Sensitivity analysis stage" << '\n';
			assignment.summary_file << ",Iteration,Avg Travel Time(min)" << '\n';
		}

		double Avg_Travel_Time = total_system_travel_time / max(1.0f, assignment.total_demand_volume);

		assignment.summary_file << iteration_number << "," << Avg_Travel_Time << "," << '\n';
		// link based computing mode, we have to collect link volume from all processors.
		if (assignment.assignment_mode == lue)
			g_fetch_link_volume_for_all_processors();
	}

	g_reset_RT_link_penalty_in_column_pool(assignment);
}


void g_reset_link_volume_in_master_program_without_columns(int number_of_links, int iteration_index, bool b_self_reducing_path_volume)
{
	int number_of_demand_periods = assignment.g_number_of_demand_periods;

	if (iteration_index == 0)
	{
		for (int i = 0; i < number_of_links; ++i)
		{
			for (int tau = 0; tau < number_of_demand_periods; ++tau)
			{
				// used in travel time calculation
				g_link_vector[i].total_volume_for_all_mode_types_per_period[tau] = 0;
			}
		}
	}
	else
	{
		double ratio = 1;
		for (int i = 0; i < number_of_links; ++i)
		{
			if (b_self_reducing_path_volume)
			{

				for (int tau = 0; tau < number_of_demand_periods; ++tau)
				{
					ratio = double(iteration_index) / double(iteration_index + 1);
					// after link volumn "tally", self-deducting the path volume by 1/(k+1) (i.e. keep k/(k+1) ratio of previous flow)
					// so that the following shortes path will be receiving 1/(k+1) flow
					g_link_vector[i].total_volume_for_all_mode_types_per_period[tau] = g_link_vector[i].total_volume_for_all_mode_types_per_period[tau] * ratio;
					g_link_vector[i].total_person_volume_for_all_mode_types_per_period[tau] = g_link_vector[i].total_person_volume_for_all_mode_types_per_period[tau] * ratio;

					for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
					{
						g_link_vector[i].volume_per_mode_type_per_period[tau][at] *= ratio;
						g_link_vector[i].converted_MEU_volume_per_period_per_at[tau][at] *= ratio;
					}
				}

			}
		}
	}
}

//***
// major function 1:  allocate memory and initialize the data
// void AllocateMemory(int number_of_nodes)
//
//major function 2: // time-dependent label correcting algorithm with double queue implementation
//int optimal_label_correcting(int origin_node, int destination_node, int departure_time, int g_debugging_flag, FILE* g_pFileDebugLog, Assignment& assignment, int time_period_no = 0, int mode_type = 1, float VOT = 10)

//	//major function: update the cost for each node at each SP tree, using a stack from the origin structure
//int tree_cost_updating(int origin_node, int departure_time, int g_debugging_flag, FILE* g_pFileDebugLog, Assignment& assignment, int time_period_no = 0, int mode_type = 1)

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
	dtalog.output() << "[PROCESS INFO] Step 3: Allocating tasks to different memory blocks for efficient computing." << '\n';

	NetworkForSP* PointerMatrx[MAX_MODETYPES][MAX_TIMEPERIODS][MAX_MEMORY_BLOCKS] = { NULL };
	NetworkForSP* RTPointerMatrx[MAX_MODETYPES][MAX_TIMEPERIODS][MAX_MEMORY_BLOCKS] = { NULL };

	int computing_zone_count = 0;

	for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
	{
		for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); ++tau)
		{
			//assign all nodes to the corresponding thread
			for (int z = 0; z < g_zone_vector.size(); ++z)
			{

				if (z < assignment.g_number_of_memory_blocks)
				{
					NetworkForSP* p_NetworkForSP = new NetworkForSP();

					p_NetworkForSP->m_origin_node_vector.push_back(g_zone_vector[z].node_seq_no);
					p_NetworkForSP->m_origin_zone_seq_no_vector.push_back(z);

					p_NetworkForSP->m_mode_type_no = at;
					p_NetworkForSP->m_tau = tau;

					computing_zone_count++;

					p_NetworkForSP->AllocateMemory(assignment.g_number_of_nodes, assignment.g_number_of_links, assignment.g_number_of_analysis_districts);

					PointerMatrx[at][tau][z] = p_NetworkForSP;

					g_NetworkForSP_vector.push_back(p_NetworkForSP);
				}
				else  // zone seq no is greater than g_number_of_memory_blocks
				{
					if (g_zone_vector[z].subarea_significance_flag == false)  // due to subarea
					{
						continue;
					}

					if (g_zone_vector[z].b_shortest_path_computing_flag == false)  // due to super zone aggregation
					{
						continue;
					}


					if (assignment.g_origin_demand_array[z] > 0.001 ||
						assignment.zone_seq_no_2_info_mapping.find(z) != assignment.zone_seq_no_2_info_mapping.end()
						)
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


	dtalog.output() << "[DATA INFO] There are " << g_NetworkForSP_vector.size() << " shortest path (SP) networks stored in memory for processing." << '\n';
	dtalog.output() << "[DATA INFO] There are " << computing_zone_count << " different agent types across various zones. These computations will be performed by the multi CPU processors." << '\n';

}

void g_assign_RT_computing_tasks_to_memory_blocks(Assignment& assignment)
{
	//fprintf(g_pFileDebugLog, "-------g_assign_computing_tasks_to_memory_blocks-------\n");
	// step 2: assign node to thread
	dtalog.output() << "[PROCESS INFO] Step 2: Assigning real time info computing tasks to memory blocks..." << '\n';

	int computing_zone_count = 0;

	int z_size = g_zone_vector.size();

	int at_size = assignment.g_ModeTypeVector.size();

	int tau_s_size = assignment.g_DemandPeriodVector.size();
	assignment.g_rt_network_pool = Allocate3DDynamicArray<NetworkForSP*>(z_size, at_size, tau_s_size);


	for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
	{
		for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); ++tau)
		{
			//assign all nodes to the corresponding thread
			for (int z = 0; z < g_zone_vector.size(); ++z)
			{
				NetworkForSP* p_NetworkForSP = new NetworkForSP();

				p_NetworkForSP->m_RT_dest_node = g_zone_vector[z].node_seq_no;
				p_NetworkForSP->m_RT_dest_zone = z;

				p_NetworkForSP->m_mode_type_no = at;
				p_NetworkForSP->m_tau = tau;

				computing_zone_count++;


				p_NetworkForSP->AllocateMemory(assignment.g_number_of_nodes, assignment.g_number_of_links, assignment.g_number_of_analysis_districts);


				assignment.g_rt_network_pool[z][at][tau] = p_NetworkForSP;  // assign real time computing network to the online colume;


				g_NetworkForRTSP_vector.push_back(p_NetworkForSP);
			}
		}
	}

	dtalog.output() << "[DATA INFO] There are " << g_NetworkForRTSP_vector.size() << " RTSP networks in memory." << '\n';

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

	for (int n = 0; n < g_NetworkForRTSP_vector.size(); ++n)
	{
		NetworkForSP* p_NetworkForSP = g_NetworkForRTSP_vector[n];
		delete p_NetworkForSP;
	}

}



void g_reset_link_volume_for_all_processors()
{
	int number_of_memory_blocks = min((int)g_NetworkForSP_vector.size(), assignment.g_number_of_memory_blocks);

#pragma omp parallel for
	for (int ProcessID = 0; ProcessID < number_of_memory_blocks; ++ProcessID)
	{
		for (int n = 0; n < g_NetworkForSP_vector.size(); n++)
		{
			if (n % number_of_memory_blocks == ProcessID)
			{
				NetworkForSP* pNetwork = g_NetworkForSP_vector[n];
				//Initialization for all non-origin nodes
				int number_of_links = assignment.g_number_of_links;
				for (int i = 0; i < number_of_links; ++i)
				{
					pNetwork->m_link_mode_type_volume_array[i] = 0;
					pNetwork->m_link_person_volume_array[i] = 0;


				}

			}
		}

	}
}


void g_fetch_link_volume_for_all_processors()
{
	for (int ProcessID = 0; ProcessID < g_NetworkForSP_vector.size(); ++ProcessID)
	{
		NetworkForSP* pNetwork = g_NetworkForSP_vector[ProcessID];

		for (int i = 0; i < g_link_vector.size(); ++i)
		{
			g_link_vector[i].total_volume_for_all_mode_types_per_period[pNetwork->m_tau] += pNetwork->m_link_mode_type_volume_array[i];
			g_link_vector[i].total_person_volume_for_all_mode_types_per_period[pNetwork->m_tau] += pNetwork->m_link_person_volume_array[i];

			g_link_vector[i].volume_per_mode_type_per_period[pNetwork->m_tau][pNetwork->m_mode_type_no] += pNetwork->m_link_mode_type_volume_array[i];

		}
	}
	// step 1: travel time based on VDF
}

// This function is used to dynamically setup the number of lanes for a specific link in a scenario
void CLink::setup_dynamic_number_of_lanes(int scenario_index)
{
	// Fetching link type from the scenario index
	int link_type = link_type_si[scenario_index];

	// Looping over all mode types
	for (int mode_type_index = 0; mode_type_index < assignment.g_ModeTypeVector.size(); mode_type_index++)
	{
		// Looping over all demand periods
		for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
		{
			// Initializing number of lanes to 1
			int number_of_lanes = 1;

			if (mode_type_index == 0)  // If the current mode is the main mode (auto)
			{
				// Fetch the number of lanes for the current scenario
				number_of_lanes = number_of_lanes_si[assignment.active_scenario_index];
			}
			else // Non-main modes (e.g., active transportation, public transportation)
			{
				// If the current mode type requires a specific assignment, bike, walk
				if (assignment.g_ModeTypeVector[mode_type_index].mode_specific_assignment_flag  ==1)
				{
					// Fetch the link type for the active scenario
					int current_link_type = g_link_vector[link_seq_no].link_type_si[assignment.active_scenario_index];

					// Fetch the number of lanes for the current mode type from the link type
					number_of_lanes = assignment.g_LinkTypeMap[link_type].lanes_mode_type[mode_type_index];
				}
				else // Non-main modes such as hov, and truck, real time info
				{
					// Fetch the number of lanes for the current scenario
					number_of_lanes = number_of_lanes_si[assignment.active_scenario_index];
				}
			}

			// Recording the number of lanes for the current period, mode type, and scenario
			recorded_lanes_per_period_per_at[tau][mode_type_index][scenario_index] = number_of_lanes;
		}
	}
}

// This function is used to calculate dynamic Volume-Delay Function (VDF) for traffic assignment

void CLink::calculate_dynamic_VDFunction(int inner_iteration_number, bool congestion_bottleneck_sensitivity_analysis_mode, int VDF_type_no)
{
	// Resetting real time (RT) waiting time at the beginning of each simulation iteration
	RT_waiting_time = 0;

	// Fetching the link type for the active scenario
	int link_type = link_type_si[assignment.active_scenario_index];

	// Looping over all demand periods
	for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
	{
		// Looping over all mode types, starting with the base mode type
		for (int mode_type_index = 0; mode_type_index < assignment.g_ModeTypeVector.size(); mode_type_index++)
		{
			double link_volume_to_be_assigned = 0;

			// Resetting the Mobility Equivalent Unit (MEU) conversion volume
			converted_MEU_volume_per_period_per_at[tau][mode_type_index] = 0;

			if (inner_iteration_number == 1 && mode_type_index == 1)
			{
				int debug_int = 1;

			}
			// Looping over all mode types to convert volumes based on MEU
			for (int mode_type_index_secondary = 0; mode_type_index_secondary < assignment.g_ModeTypeVector.size(); mode_type_index_secondary++)
			{
				// Fetching the MEU conversion value for this mode type pair
				double MEU_conversion_value = assignment.g_LinkTypeMap[link_type].meu_matrix[mode_type_index][mode_type_index_secondary];

				if (mode_type_index == mode_type_index_secondary)
					MEU_conversion_value = 1.0;

				// Adding the converted volume to the MEU conversion volume
				converted_MEU_volume_per_period_per_at[tau][mode_type_index] += volume_per_mode_type_per_period[tau][mode_type_index_secondary] * MEU_conversion_value;
			}

			double mode_hourly_capacity = 0;
			double mode_FFFTT = 0;  // free-flow travel time
			double mode_peak_load_factor = 1;

			// Handling primary mode or non-major modes such as HOV or truck, they need to be mapped to the major auto mode for traffic assignment
			if (mode_type_index == 0 || assignment.g_ModeTypeVector[mode_type_index].mode_specific_assignment_flag == 0)
			{
				// The base mode type index
				int primary_mode_index = 0;

				mode_FFFTT = VDF_period[tau].FFTT_at[primary_mode_index];  // free flow travel time 
				mode_peak_load_factor = assignment.g_LinkTypeMap[link_type].peak_load_factor_period_at[tau][primary_mode_index];  // peak load factor 

				// Adding preloaded volume to the converted MEU volume for this period and mode type
				link_volume_to_be_assigned = converted_MEU_volume_per_period_per_at[tau][primary_mode_index] + VDF_period[tau].preload;

				// Fetching the lane based ultimate hourly capacity
				mode_hourly_capacity = VDF_period[tau].lane_based_ultimate_hourly_capacity;
			}
			else
			{
				// For active and public transportation modes, they have their own specific speed limit and capacity 

				mode_FFFTT = link_distance_VDF / assignment.g_LinkTypeMap[link_type_si[assignment.active_scenario_index]].free_speed_at[mode_type_index] * 60.0;
				mode_peak_load_factor = assignment.g_LinkTypeMap[link_type].peak_load_factor_period_at[tau][mode_type_index];

				// Fetching the converted MEU volume for this period and mode type
				link_volume_to_be_assigned = converted_MEU_volume_per_period_per_at[tau][mode_type_index];

				// Fetching the hourly capacity for this mode type
				mode_hourly_capacity = assignment.g_LinkTypeMap[link_type].capacity_at[mode_type_index];
			}

			// Fetching the number of lanes for this period and mode type
			VDF_period[tau].nlanes = recorded_lanes_per_period_per_at[tau][mode_type_index][assignment.active_scenario_index];

			if (VDF_period[tau].nlanes == 0)  // if the number of lanes is zero, skip to next iteration
			{
				int idebug;
				idebug = 1;
				link_avg_travel_time_per_period[tau][mode_type_index] = 999900;
				continue;
			}

			CLinkType link_type; 

			link_type = assignment.g_LinkTypeMap[this->link_type_si[assignment.active_scenario_index]];


			// Calculate travel time based on QVDF
			link_avg_travel_time_per_period[tau][mode_type_index] = VDF_period[tau].calculate_travel_time_based_on_QVDF(
				mode_type_index,
				VDF_period[tau].FFTT_at[mode_type_index],
				link_volume_to_be_assigned,
				mode_hourly_capacity,
				mode_peak_load_factor,
				this->model_speed,
				this->est_volume_per_hour_per_lane,
				link_type, tau, this->link_avg_co2_emit_per_mode, this->link_avg_nox_emit_per_mode);

			// Add additional penalty if exists
			if (fabs(penalty_si_at[assignment.active_scenario_index][mode_type_index][tau]) > 0.0001)
			{
				link_avg_travel_time_per_period[tau][mode_type_index] += penalty_si_at[assignment.active_scenario_index][mode_type_index][tau];
			}

			// Set link volume for the VDF period
			VDF_period[tau].link_volume = link_volume_to_be_assigned;

			if (assignment.g_ModeTypeVector[mode_type_index].mode_specific_assignment_flag == 0)
			{
				// The condition checks if the mode specific assignment flag for the given mode type index (could be truck, hov, real time in, CAV, EV mode) is zero.
				// If it is zero, we assume that the mode type does not have its own specific values and we need to use default values from the base mode type, which is 'auto' in this case.

				// The base mode type index for 'auto'.
				int primary_mode_index = 0;

				// Transferring data from the primary 'auto' mode to the associated modes when the associated mode does not have its own specific values. 

				// Transferring average travel time per period from 'auto' to the current mode type.
				link_avg_travel_time_per_period[tau][mode_type_index] = link_avg_travel_time_per_period[tau][primary_mode_index];

				// Transferring average CO2 emission per mode from 'auto' to the current mode type.
				link_avg_co2_emit_per_mode[tau][mode_type_index] = link_avg_co2_emit_per_mode[tau][primary_mode_index];

				// Transferring average NOx emission per mode from 'auto' to the current mode type.
				link_avg_nox_emit_per_mode[tau][mode_type_index] = link_avg_nox_emit_per_mode[tau][primary_mode_index];
			}
		}
	}
}
//
//void  CLink::calculate_dynamic_VDFunction(int inner_iteration_number, bool congestion_bottleneck_sensitivity_analysis_mode, int VDF_type_no)
//{
//	RT_waiting_time = 0; // reset RT_travel time for each end of simulation iteration
//
//	int link_type = link_type_si[assignment.active_scenario_index];
//	// for each time period
//	for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
//	{
//		for (int mode_type_index = 0; mode_type_index < assignment.g_ModeTypeVector.size(); mode_type_index++)  //first consider base mode index
//		{
//				double link_volume_to_be_assigned = 0;
//				// mobility equivalenet unit MEU mapping matrix
//				converted_MEU_volume_per_period_per_at[tau][mode_type_index] = 0;  // reset
//			for (int mode_type_index_secondary = 0; mode_type_index_secondary < assignment.g_ModeTypeVector.size(); mode_type_index_secondary++)  // the secondary mode through the conversion of MEU
//				{
//
//					double MEU_conversion_value = assignment.g_LinkTypeMap[link_type].meu_matrix[mode_type_index_secondary][mode_type_index];  // mode_type_index is the base mode_type_index, at_k is the external mode_type_index with impact
//					converted_MEU_volume_per_period_per_at[tau][mode_type_index]+= volume_per_mode_type_per_period[tau][mode_type_index_secondary] * MEU_conversion_value;
//
//				}
//
//				double mode_hourly_capacity = 0;
//				double mode_FFFTT = 0;
//				double mode_peak_load_factor = 1;
//
//				if (mode_type_index == 0 || assignment.g_ModeTypeVector[mode_type_index].mode_specific_assignment_flag == 0)  // primary mode or non-major mode such as HOV or truck they need to be mapped to the major auto mode to perform traffic assignment
//				{
//					int at_base = 0;  // so the base is the mode index  = 0
//					mode_FFFTT = VDF_period[tau].FFTT_at[at_base];  // free flow travel time 
//					mode_peak_load_factor = assignment.g_LinkTypeMap[link_type].peak_load_factor_period_at[tau][at_base];  // peak load factor 
//					link_volume_to_be_assigned = converted_MEU_volume_per_period_per_at[tau][at_base] + VDF_period[tau].preload; // converted link volume 
//					mode_hourly_capacity = VDF_period[tau].lane_based_ultimate_hourly_capacity;  // featch the lane based ultimate hourly capacity 
//				}
//				else
//				{ // otherwisem, it this for active transportation and public transportation assignment, they have their own mode specific dedicated speed limit and capacity 
//
//					mode_FFFTT = link_distance_VDF / assignment.g_LinkTypeMap[link_type_si[assignment.active_scenario_index]].free_speed_at[mode_type_index] * 60.0; // 60.0 for 60 min per hour
//					mode_peak_load_factor = assignment.g_LinkTypeMap[link_type].peak_load_factor_period_at[tau][mode_type_index];
//					link_volume_to_be_assigned = converted_MEU_volume_per_period_per_at[tau][mode_type_index]; // its own agent type
//					mode_hourly_capacity = assignment.g_LinkTypeMap[link_type].capacity_at[mode_type_index];
//				}
//
//
//				// direct read 
//				VDF_period[tau].nlanes = recorded_lanes_per_period_per_at[tau][mode_type_index][assignment.active_scenario_index];
//
//
//				if (VDF_period[tau].nlanes == 0)  // skip if the number of lanes is zero
//				{
//					int idebug;
//					idebug = 1;
//					link_avg_travel_time_per_period[tau][mode_type_index] = 999999;
//					continue;
//				}
//
//
//				link_avg_travel_time_per_period[tau][mode_type_index] = VDF_period[tau].calculate_travel_time_based_on_QVDF(
//					mode_type_index,
//					VDF_period[tau].FFTT_at[mode_type_index],
//					link_volume_to_be_assigned,
//					mode_hourly_capacity,
//					mode_peak_load_factor,
//					this->model_speed,
//					this->est_volume_per_hour_per_lane);
//
//				if(fabs(penalty_si_at[assignment.active_scenario_index][mode_type_index][tau]) >0.0001)
//				{
//					link_avg_travel_time_per_period[tau][mode_type_index] += penalty_si_at[assignment.active_scenario_index][mode_type_index][tau];  // add additional penalty
//				}
//
//
//				VDF_period[tau].link_volume = link_volume_to_be_assigned;
//			}
//		}
//	
//
//}
int read_route_information_to_replace_column_generation_and_ODME()
{
				int path_counts = 0;
				float sum_of_path_volume = 0;
				int line_no = 0;
				CDTACSVParser parser;
				if (parser.OpenCSVFile("route.csv", false))
				{
					int total_path_in_demand_file = 0;
					// read agent file line by line,

					int o_zone_id, d_zone_id;

					int this_departure_time_profile_no = 0;
					string mode_type, demand_period;

					std::vector <int> node_sequence;

					dtalog.output() << "[STATUS INFO] Reading preload route file" << '\n';

					while (parser.ReadRecord())
					{


						total_path_in_demand_file++;
						if (total_path_in_demand_file % 100000 == 0)
							dtalog.output() << "[DATA INFO] total_path_in_demand_file is " << total_path_in_demand_file << '\n';

						parser.GetValueByFieldName("o_zone_id", o_zone_id);
						parser.GetValueByFieldName("d_zone_id", d_zone_id);

						CAgentPath agent_path_element;

						parser.GetValueByFieldName("path_id", agent_path_element.path_id, true);


						int from_zone_seq_no = 0;
						int to_zone_seq_no = 0;
						from_zone_seq_no = assignment.g_zoneid_to_zone_seq_no_mapping[o_zone_id];
						to_zone_seq_no = assignment.g_zoneid_to_zone_seq_no_mapping[d_zone_id];
						int from_zone_sindex = g_zone_vector[from_zone_seq_no].sindex;
						if (from_zone_sindex == -1)
							continue;

						int to_zone_sindex = g_zone_vector[to_zone_seq_no].sindex;
						if (to_zone_sindex == -1)
							continue;

						double volume = 0;
						parser.GetValueByFieldName("volume", volume);
						agent_path_element.volume = volume;
						path_counts++;
						sum_of_path_volume += agent_path_element.volume;

						string mode_type, demand_period;
						int mode_type_no = 0;
						int demand_period_no = 0;

						parser.GetValueByFieldName("mode_type", mode_type);
						parser.GetValueByFieldName("demand_period", demand_period);

						//char time_interval_field_name[20];
						if (assignment.mode_type_2_seqno_mapping.find(mode_type) != assignment.mode_type_2_seqno_mapping.end())
								mode_type_no = assignment.mode_type_2_seqno_mapping[mode_type];

						if (assignment.demand_period_to_seqno_mapping.find(demand_period) != assignment.demand_period_to_seqno_mapping.end())
							demand_period_no = assignment.demand_period_to_seqno_mapping[demand_period];

						assignment.total_demand[mode_type_no][demand_period_no] += agent_path_element.volume;
						assignment.g_column_pool[from_zone_sindex][to_zone_sindex][mode_type_no][demand_period_no].od_volume[assignment.active_scenario_index] += agent_path_element.volume;
						assignment.g_column_pool[from_zone_sindex][to_zone_sindex][mode_type_no][demand_period_no].departure_time_profile_no = this_departure_time_profile_no;  // to be updated
						assignment.total_route_demand_volume += agent_path_element.volume;
						assignment.g_origin_demand_array[from_zone_seq_no] += agent_path_element.volume;

						int analysis_district_id = assignment.g_zone_seq_no_to_analysis_distrct_id_mapping[from_zone_seq_no];
						//g_district_summary_map[analysis_district_id].record_origin_2_district_volume(mode_type_no, agent_path_element.volume);

						//apply for both agent csv and routing policy
						assignment.g_column_pool[from_zone_sindex][to_zone_sindex][mode_type_no][demand_period_no].bfixed_route = true;

						bool bValid = true;

						string path_node_sequence;
						parser.GetValueByFieldName("node_sequence", path_node_sequence);

						if (path_node_sequence.size() == 0)
							continue;

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
									node_sum += internal_node_seq_no * link_seq_no;
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

							CColumnVector* pColumnVector = &(assignment.g_column_pool[from_zone_sindex][to_zone_sindex][mode_type_no][demand_period_no]);
							pColumnVector->departure_time_profile_no = this_departure_time_profile_no;
							// we cannot find a path with the same node sum, so we need to add this path into the map,
							if (pColumnVector->path_node_sequence_map.find(node_sum) == pColumnVector->path_node_sequence_map.end())
							{
								// add this unique path
								int path_count = pColumnVector->path_node_sequence_map.size();
								pColumnVector->path_node_sequence_map[node_sum].path_seq_no = path_count;
								pColumnVector->path_node_sequence_map[node_sum].path_id = agent_path_element.path_id;
								pColumnVector->path_node_sequence_map[node_sum].path_volume = 0;
								pColumnVector->path_node_sequence_map[node_sum].path_toll = 0;

								pColumnVector->path_node_sequence_map[node_sum].AllocateVector(node_no_sequence, link_no_sequence, false);
							}

							pColumnVector->path_node_sequence_map[node_sum].path_volume += agent_path_element.volume;
							pColumnVector->path_node_sequence_map[node_sum].path_preload_volume += volume;

						}

						line_no++;
						if (line_no % 100000 == 0)
							dtalog.output() << "[STATUS INFO] Reading demand file line no =  " << line_no / 1000 << "k" << '\n';

					}
					dtalog.output() << "[DATA INFO] total_demand_volume loaded from route file is " << sum_of_path_volume << " with " << path_counts << " paths." << '\n';

					return 1;
				}
				else
				{
				return 0;
				}

}

double network_assignment(int assignment_mode, int column_generation_iterations, int column_updating_iterations, int ODME_iterations, int sensitivity_analysis_iterations, int simulation_iterations, int number_of_memory_blocks, int length_unit_flag, int speed_unit_flag)
{

	assignment.summary_file << "[PROCESS INFO] Step 0: reading scenario_index_list.csv" << '\n';
	dtalog.output() << "[PROCESS INFO] Step 0.1: reading scenario_index_list.csv" << '\n';

	CDTACSVParser parser_scenario_index_file_list;
	if (parser_scenario_index_file_list.OpenCSVFile("scenario_index_list.csv", false))
	{

		while (parser_scenario_index_file_list.ReadRecord())  // if this line contains [] mark, then we will also read field headers.
		{
			DTAScenario element;

			int activate = 0;
			parser_scenario_index_file_list.GetValueByFieldName("activate", activate);
			if(activate ==1)
			{

			parser_scenario_index_file_list.GetValueByFieldName("scenario_index", element.scenario_index);
			parser_scenario_index_file_list.GetValueByFieldName("scenario_name", element.scenario_name);
			assignment.g_active_DTAscenario_map[element.scenario_index] = assignment.g_DTA_scenario_vector.size();
			assignment.g_DTA_scenario_vector.push_back(element);

			if (assignment.g_DTA_scenario_vector.size() > MAX_SCENARIOS - 1)
			{
				assignment.summary_file << "[ERROR] " << "MAX_SCENARIOS in DTALite " << MAX_SCENARIOS -1 << ". Users have too many scenarios now as " << assignment.g_DTA_scenario_vector.size() << ". Please contact the developer to obtain a version with larger scenario size if needed." << '\n';
				dtalog.output() << "[ERROR]" << "MAX_SCENARIOS in DTALite " << MAX_SCENARIOS-1 << ". Users have too many scenarios now as " << assignment.g_DTA_scenario_vector.size() << ". Please contact the developer to obtain a version with larger scenario size if needed." << '\n';
				g_program_stop();

			}
			assignment.summary_file << "[STATUS INFO] scenario_index=" << element.scenario_index << ",scenario_name=," << element.scenario_name.c_str() << '\n';

			}
		}
	}
	else
	{
		dtalog.output() << "[ERROR] File scenario_index_list.csv does not exist!" << '\n';
		g_program_stop();
	}

	g_number_of_active_scenarios = assignment.g_DTA_scenario_vector.size();
	clock_t start_t0, end_t0, total_t0;
	int signal_updating_iterations = 0;
	start_t0 = clock();
	// k iterations for column generation
	assignment.g_number_of_column_generation_iterations = column_generation_iterations;
	assignment.g_number_of_column_updating_iterations = column_updating_iterations;
	assignment.g_number_of_ODME_iterations = ODME_iterations;
	assignment.g_number_of_sensitivity_analysis_iterations_for_dtm = sensitivity_analysis_iterations;

	assignment.active_scenario_index = 0;

	// 0: link UE: 1: path UE, 2: Path SO, 3: path resource constraints

	if(assignment_mode==0)
		assignment.assignment_mode = lue;

	if (assignment_mode == 1)
		assignment.assignment_mode = path_based_assignment;

	if (assignment_mode == 2)
		assignment.assignment_mode = simulation_dta;

	assignment.g_number_of_memory_blocks = min(max(1, number_of_memory_blocks), MAX_MEMORY_BLOCKS);

	if (assignment.assignment_mode == lue)
		column_updating_iterations = 0;

	// step 1: read input data of network / demand tables / Toll
	g_read_input_data(assignment);
	// to add back: Xueosong g_ReadInformationConfiguration(assignment);
	//g_reload_timing_arc_data(assignment); // no need to load timing data from timing.csv
	// to do: g_ReadOutputFileConfiguration(assignment);

	g_ReadDemandFileBasedOnDemandFileList(assignment);

	// after we read the physical links and demand files
// we create virtual connectors
	for (int i = 0; i < g_node_vector.size(); ++i)
	{

		if (g_node_vector[i].zone_org_id >= 1) // for each physical node
		{ // we need to make sure we only create two way connectors between nodes and zones

				// for each node-zone pair: create a pair of connectors with the agent-type related acess_map
			int zone_org_id = g_node_vector[i].zone_org_id;
			int internal_from_node_seq_no, internal_to_node_seq_no, zone_seq_no;

			internal_from_node_seq_no = g_node_vector[i].node_seq_no;


			if (assignment.zone_id_to_seed_zone_id_mapping.find(zone_org_id) != assignment.zone_id_to_seed_zone_id_mapping.end())
			{
				// seed zone id has been identified.
				zone_org_id = assignment.zone_id_to_seed_zone_id_mapping[zone_org_id];  // update zone id
			}

			int node_id = assignment.zone_id_to_centriod_node_id_mapping[zone_org_id];


			internal_to_node_seq_no = assignment.g_node_id_to_seq_no_map[node_id];


			zone_seq_no = assignment.g_zoneid_to_zone_seq_no_mapping[zone_org_id];

			// we need to mark all accessble model on this access links, so we can handle that in the future for each agent type's memory block in shortest path
			// incomming virtual connector
			g_add_new_virtual_connector_link(internal_from_node_seq_no, internal_to_node_seq_no, g_node_vector[i].mode_type_str, -1);
			// outgoing virtual connector
			g_add_new_virtual_connector_link(internal_to_node_seq_no, internal_from_node_seq_no, g_node_vector[i].mode_type_str, zone_seq_no);
			// result is that: we have a unique pair of node-zone access link in the overall network, but still carry mode_type_acess_map for mode types with access on this node-zone connector

		}

	}


	//step 2: allocate memory and assign computing tasks
	g_assign_computing_tasks_to_memory_blocks(assignment); // static cost based label correcting

	// definte timestamps
	clock_t start_t, end_t, total_t;
	clock_t start_simu, end_simu, total_simu;
	start_t = clock();


	clock_t iteration_t, cumulative_lc, cumulative_cp, cumulative_lu;

	//step 3: column generation stage: find shortest path and update path cost of tree using TD travel time
	
	dtalog.output() << "[PROCESS INFO] Step 4: Optimizing traffic assignment with Column Generation, an iterative method of adding (path) variables to solve large-scale problems. The optimization process involves creating and updating a set of candidate routes, known as the column pool" << '\n';
	dtalog.output() << "[DATA INFO] Perform " << assignment.g_number_of_column_generation_iterations << " iterations for the Column Generation process" << '\n';

	assignment.summary_file << "[PROCESS INFO]Step 4: Column Generation for Traffic Assignment" << '\n';
	assignment.summary_file << ",Total number of column generation iterations =, " << assignment.g_number_of_column_generation_iterations << '\n';
	assignment.summary_file << ",Total number of column generation iterations =, " << assignment.g_number_of_column_generation_iterations << '\n';


	g_outputZonalHierarchyMapping(assignment);

	for (int sii = 0; sii < assignment.g_DTA_scenario_vector.size(); sii++)
	{



		assignment.active_scenario_index = assignment.g_DTA_scenario_vector[sii].scenario_index;
		dtalog.output() << "[DATA INFO] global scenario index = " << assignment.g_DTA_scenario_vector[sii].scenario_index << " --------------------------" << '\n';
		
		g_load_dynamic_traffic_management_file(assignment);

		//if (read_route_information_to_replace_column_generation_and_ODME() == 0)
		{  // apply column generation and ODME

			dtalog.output()  << "[DATA INFO] " 
				<< std::setw(12) << std::left << "Iter. No."
				<< std::setw(12) << std::left << "CPU time(s)"
				<< std::setw(30) << std::left << "Sys. Wide Travel Time (TT)"
				<< std::setw(20) << std::left << "Least system TT"
				<< std::setw(12) << std::left << "Gap (%)" << '\n';

			for (int iteration_number = 0; iteration_number < max(1, assignment.g_number_of_column_generation_iterations); iteration_number++)
			{
				end_t = clock();
				iteration_t = end_t - start_t;


				// step 3.1 update travel time and resource consumption
				clock_t start_t_lu = clock();

				double total_system_wide_travel_time = 0;
				double total_least_system_travel_time = 0;
				double total_travel_distance = 0;
				double avg_speed = 0;

				// initialization at beginning of shortest path
				total_system_wide_travel_time = update_link_travel_time_and_cost(iteration_number, total_travel_distance);

				if (assignment.assignment_mode == lue)
				{
					//fw link based UE
					g_reset_link_volume_in_master_program_without_columns(g_link_vector.size(), iteration_number, true);
					g_reset_link_volume_for_all_processors();
				}
				else
				{
					// we can have a recursive formula to reupdate the current link volume by a factor of k/(k+1),
					//  and use the newly generated path flow to add the additional 1/(k+1)
					g_reset_and_update_link_volume_based_on_columns(g_link_vector.size(), iteration_number, true, false);
				}

				if (dtalog.debug_level() >= 3)
				{
					dtalog.output() << "[DATA INFO] Results:" << '\n';
					for (int i = 0; i < g_link_vector.size(); ++i) {
						dtalog.output() << "[DATA INFO] link: " << g_node_vector[g_link_vector[i].from_node_seq_no].node_id << "-->"
							<< g_node_vector[g_link_vector[i].to_node_seq_no].node_id << ", "
							<< "flow count:" << g_link_vector[i].total_volume_for_all_mode_types_per_period[0] << '\n';
					}
				}

				end_t = clock();
				iteration_t = end_t - start_t_lu;
				// g_fout << "Link update with CPU time " << iteration_t / 1000.0 << " s; " << (end_t - start_t) / 1000.0 << " s" << '\n';

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
				//    int mode_type_no = g_NetworkForSP_vector[ProcessID]->m_mode_type_no;

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
					for (int blk = 0; blk < assignment.g_ModeTypeVector.size() * assignment.g_DemandPeriodVector.size(); ++blk)
					{
						int network_copy_no = blk * assignment.g_number_of_memory_blocks + ProcessID;
						if (network_copy_no >= g_NetworkForSP_vector.size())
							continue;

						NetworkForSP* pNetwork = g_NetworkForSP_vector[network_copy_no];

						for (int o_node_index = 0; o_node_index < pNetwork->m_origin_node_vector.size(); ++o_node_index)
						{
							start_t_lc = clock();
							pNetwork->optimal_label_correcting(ProcessID, &assignment, iteration_number, o_node_index, false, false, false, false);


							end_t = clock();
							cumulative_lc += end_t - start_t_lc;

							start_t_cp = clock();
							double total_origin_least_travel_time = pNetwork->backtrace_shortest_path_tree(assignment, iteration_number, o_node_index, false);


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
				if (assignment.assignment_mode == lue)
					g_fetch_link_volume_for_all_processors();

				// g_fout << "LC with CPU time " << cumulative_lc / 1000.0 << " s; " << '\n';
				// g_fout << "column generation with CPU time " << cumulative_cp / 1000.0 << " s; " << '\n';

				//****************************************//

				// last iteraion before performing signal timing updating
				double relative_gap = max((total_system_wide_travel_time - total_least_system_travel_time) / max(0.00001, total_least_system_travel_time),0.0);

				double CPU_Running_Time = cumulative_cp / 1000.0;
				int number_of_agents = 0;
				double Avg_Travel_Time = total_system_wide_travel_time / max(1.0f, assignment.total_demand_volume);

				dtalog.output() << "[DATA INFO] "
					<< std::setw(12) << std::left << iteration_number
					<< std::setw(12) << std::left << iteration_t * 1.0 / 1000.0
					<< std::scientific << std::setw(30) << std::left << total_system_wide_travel_time
					<< std::scientific << std::setw(20) << std::left << total_least_system_travel_time
					<< std::fixed << std::setw(12) << std::left << relative_gap * 100 << '\n';

				if (iteration_number == 0)
				{

					assignment.summary_file << "[DATA INFO] Iteration, CPU running time (sec), # of agents, Avg Travel Time(min),  Avg UE gap %" << '\n';
				}

				assignment.summary_file << iteration_number << "," << CPU_Running_Time << "," << assignment.total_demand_volume << "," << Avg_Travel_Time << "," << relative_gap * 100 << '\n';
			}

			if (assignment.active_scenario_index == 0)
			{
				g_OutputModelFiles(10);

			}
			// column updating stage: for given column pool, update volume assigned for each column
			dtalog.output() << "[PROCESS INFO] Step 5: Column Pool Updating" << '\n';
			dtalog.output() << "[DATA INFO] Total number of column pool updating iterations = " << column_updating_iterations << '\n';

			assignment.summary_file << "[PROCESS INFO] Step 5: column pool-based flow updating for traffic assignment " << '\n';
			assignment.summary_file << ",# of flow updating iterations=," << column_updating_iterations << '\n';

			start_t = clock();
			g_column_pool_optimization(assignment, column_updating_iterations, false);


			

			// post-processing for collecting path flow to link flow
			// post-processsing route assignment aggregation
			if (assignment.assignment_mode != lue)
			{
				// we can have a recursive formulat to reupdate the current link volume by a factor of k/(k+1),
				// and use the newly generated path flow to add the additional 1/(k+1)
				g_reset_and_update_link_volume_based_on_columns(g_link_vector.size(), column_generation_iterations, false, false);
			}
			else
				g_reset_link_volume_in_master_program_without_columns(g_link_vector.size(), column_generation_iterations, false);

			// initialization at the first iteration of shortest path
			double total_distance = 0;
			update_link_travel_time_and_cost(column_generation_iterations, total_distance);


			dtalog.output() << "[PROCESS INFO] Step 6: OD demand matrix estimation if file sensor_data.csv is provided." << '\n';
			if (assignment.g_number_of_ODME_iterations >= 1)
			{

				// record before odme volume
				for (int i = 0; i < g_link_vector.size(); ++i)
				{
					for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
					{
						g_link_vector[i].VDF_period[tau].volume_before_odme = g_link_vector[i].total_volume_for_all_mode_types_per_period[tau] + g_link_vector[i].VDF_period[tau].preload;
					}
				}



				assignment.summary_file << "[PROCESS INFO] Step 6: OD demand matrix estimation" << '\n';

				g_update_odme_volume_in_column_pool(assignment, 0);  // ODME before
				assignment.Demand_ODME(ODME_iterations);
				g_update_odme_volume_in_column_pool(assignment, 1);  // ODME after


				assignment.summary_file << ",# of ODME_iterations=," << ODME_iterations << '\n';

				// record after odme volume
				for (int i = 0; i < g_link_vector.size(); ++i)
				{
					for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
					{
						g_link_vector[i].VDF_period[tau].volume_after_odme = g_link_vector[i].total_volume_for_all_mode_types_per_period[tau] + g_link_vector[i].VDF_period[tau].preload;
					}
				}

				
			}

		}  // end of ODME
		// stage II sensitivity analysis stage


		assignment.summary_file << "[PROCESS INFO] Step 7: perform sensitivity analysis if dynamic_traffic_management.csv is provided for dtm_type = lane_closure or dms. " << '\n';
		dtalog.output() << "[PROCESS INFO] Step 7: Performing Sensitivity Analysis. Proceeds only if dynamic_traffic_management.csv is provided with dtm_type set to either 'lane_closure'. " << '\n';
		
		g_reset_link_district_performance_per_scenario(assignment);
		g_record_link_district_performance_per_scenario(assignment, 0);

		if (assignment.g_number_of_sensitivity_analysis_iterations_for_dtm >= 0)  // for real-time information, where historical path flows are kept the same
		{


			dtalog.output() << "[PROCESS INFO] Step 7.1: Applying dynamic traffic_management scenarios" << '\n';

			int count = 0; 
			// apply dynamic traffic managementscenarios
			for (int i = 0; i < g_link_vector.size(); ++i)
			{
				for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
				{
					// used in travel time calculation
					if (g_link_vector[i].VDF_period[tau].dynamic_traffic_management_flag == -1 /*sa flag*/)  // we
					{
						float old_lanes = g_link_vector[i].VDF_period[tau].nlanes;
						g_link_vector[i].VDF_period[tau].nlanes = g_link_vector[i].VDF_period[tau].lane_closure_final_lanes; // apply the absolute lane changes
						g_link_vector[i].VDF_period[tau].allowed_uses[assignment.active_scenario_index] = g_link_vector[i].VDF_period[tau].sa_allowed_uses; // apply the lane changes

						count += 1; 
						int at = 0;
						if (g_link_vector[i].VDF_period[tau].nlanes < 0.01)
							g_link_vector[i].VDF_period[tau].FFTT_at[at] = 1440; // prevent free flow travel time under zero flow due to blockage
					}

				}
			}


			dtalog.output() << "[DATA INFO]  # of DTM records = " << count << '\n';

			// apply demand side scenarios

			for (int orig = 0; orig < g_zone_vector.size(); ++orig)  // o
			{
				CColumnVector* p_column_pool;
				std::map<int, CColumnPath>::iterator it, it_begin, it_end;
				int from_zone_sindex = g_zone_vector[orig].sindex;
				if (from_zone_sindex == -1)
					continue;

				for (int dest = 0; dest < g_zone_vector.size(); ++dest) //d
				{

					int to_zone_sindex = g_zone_vector[dest].sindex;
					if (to_zone_sindex == -1)
						continue;

					double local_scale_factor = 1.0;
					int o_district = assignment.g_zone_seq_no_to_analysis_distrct_id_mapping[orig];
					int d_district = assignment.g_zone_seq_no_to_analysis_distrct_id_mapping[dest];


					if (assignment.SA_o_district_id_factor_map.find(o_district) != assignment.SA_o_district_id_factor_map.end())
					{
						local_scale_factor = assignment.SA_o_district_id_factor_map[o_district];
					}

					//d based factor
					if (assignment.SA_d_district_id_factor_map.find(d_district) != assignment.SA_d_district_id_factor_map.end())
					{
						local_scale_factor = assignment.SA_d_district_id_factor_map[d_district];
					}

					int od_district_key = o_district * 1000 + d_district;
					if (assignment.SA_od_district_id_factor_map.find(od_district_key) != assignment.SA_od_district_id_factor_map.end())
					{
						local_scale_factor = assignment.SA_od_district_id_factor_map[od_district_key];
					}


					for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)  //at
					{
						for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); ++tau)  //tau
						{
							p_column_pool = &(assignment.g_column_pool[from_zone_sindex][to_zone_sindex][at][tau]);
							if (p_column_pool->od_volume[assignment.active_scenario_index] > 0 && fabs(local_scale_factor - 1.0) > 0.001)
							{
								p_column_pool->od_volume[assignment.active_scenario_index] = p_column_pool->od_volume[assignment.active_scenario_index] * local_scale_factor;
							}
						}
					}
				}
			}

			//end of scenario

			double total_system_travel_time = 0;
			double total_travel_distance = 0;

			// initialization at beginning of shortest path
			total_system_travel_time = update_link_travel_time_and_cost(0, total_travel_distance);
			
			dtalog.output() << "[PROCESS INFO] Step 7.2: record route volume before applying dynamic traffic_management scenarios" << '\n';


			g_update_sa_volume_in_column_pool(assignment, 0);  // SA before

			dtalog.output() << "[PROCESS INFO] Step 7.3: path column regeneration for applying dynamic traffic_management scenarios" << '\n';

			g_column_regeneration(assignment, true);

			// stage II: right after column generation, we will enable the path modifications for DMS users

			if (assignment.active_dms_count >= 1)
			{
				dtalog.output() << "[PROCESS INFO] Step 7.4: path column modification for dms users" << '\n';

				g_column_pool_route_modification(assignment, column_updating_iterations);
			}
			else
			{
				dtalog.output() << "[PROCESS INFO] Step 7.4: skip: path column modification for dms users, as dms count = 0" << '\n';

			}

			dtalog.output() << "[PROCESS INFO] Step 7.5: path column pool optimization for real time information users" << '\n';
			// stage III: after DMS path modification, we enable further user equilibirum computing for real time info users by switching their routes
			g_column_pool_optimization(assignment, assignment.g_number_of_sensitivity_analysis_iterations_for_dtm, true);

			dtalog.output() << "[PROCESS INFO] Step 7.6: record route volume after applying dynamic traffic_management scenarios" << '\n';
			g_update_sa_volume_in_column_pool(assignment, 1);  // SA after
		}// SA change

		g_classification_in_column_pool(assignment);

		dtalog.output() << "[PROCESS INFO] Step 8: Executing Traffic Simulation. Proceeds only if simulation_output is set to 1 in settings.csv . " << '\n';

		assignment.summary_file << "[PROCESS INFO] Step 8: Executing Traffic Simulation. Proceeds only if simulation_output is set to 1 in settings.csv. " << '\n';

		if (simulation_iterations >= 1)
		{
			start_simu = clock();

			dtalog.output() << "[PROCESS INFO] Step 8: Simulation for traffic assignment.." << '\n';
			assignment.STTrafficSimulation();
			end_simu = clock();
			total_simu = end_simu - end_simu;

			dtalog.output() << "[DATA INFO] CPU running time for traffic simulation: " << total_simu / 1000.0 << " s" << '\n';
			
		}

		end_t = clock();
		total_t = (end_t - start_t);

		dtalog.output() << "[DATA INFO] Total CPU running time for the entire process = " << total_t / 1000.0 << " s" << '\n';

		start_t = clock();

		//step 5: output simulation results of the new demand
		dtalog.output() << "[PROCESS INFO] Step 9: Generating Zonal Hierarchy Mapping (e.g., Zone to Super-Zone and Zone to District Mapping)" << '\n';
		assignment.summary_file << "[PROCESS INFO] Step 9: Generating Zonal Hierarchy Mapping (e.g., Zone to Super-Zone and Zone to District Mapping) " << '\n';

		// zonal_hierarchy_mapping


		dtalog.output() << "[PROCESS INFO] Step 10: Outputting Traffic Assignment Results." << '\n';

		g_output_assignment_result(assignment,0);

	// g_output_demand_bin(assignment);
	if (assignment.g_subarea_shape_points.size() >= 3) // if there is a subarea defined.
	{
	g_output_assignment_result(assignment, 1);
	}

	g_output_choice_set_result(assignment);

	if (assignment.assignment_mode == simulation_dta)
	{
	g_output_simulation_result(assignment);
	}


	// end of scenario computing, clean the memory
	g_column_pool_reset(assignment);
	}
	g_output_accessibility_result(assignment);
	g_output_assignment_summary_result(assignment, 0);
	g_output_2_way_assignment_summary_result(assignment, 0);
	g_OutputSummaryFiles(assignment);
	// at the end of simulation
	// validation step if reading data are available
	bool b_sensor_reading_data_available = false;
	CDTACSVParser parser_reading;
	if (parser_reading.OpenCSVFile("Reading.csv", false))
	{
		parser_reading.CloseCSVFile();
		b_sensor_reading_data_available = true;
	}

	//    g_output_dynamic_queue_profile();
		//

	end_t = clock();
	total_t = (end_t - start_t);
	dtalog.output() << "[DATA INFO] CPU running time for outputting simulation results: " << total_t / 1000.0 << " s" << '\n';

	dtalog.output() << "[STATUS INFO] Freeing memory.." << '\n';


	// temp comment out
	//for (int i = 0; i < g_link_vector.size(); ++i)
	//	g_link_vector[i].free_memory();
	//g_link_vector.clear();

	//g_node_vector.clear();

	end_t0 = clock();
	total_t0 = (end_t0 - start_t0);
	int second = total_t0 / 1000.0;
	int min = second / 60;
	int sec = second - min * 60;
	dtalog.output() << "[DATA INFO] CPU running time for the entire process: " << min << " min " << sec << " sec" << '\n';
	dtalog.output() << "[STATUS INFO] DTALite computing process has successfully completed. Congratulations on a successful execution! Feel free to review the results and explore the generated outputs. Thank you for using DTALite and contributing to the advancement of transportation analysis and open science." << std::endl;

	return 1;
}

int Assignment::update_real_time_info_path(CAgent_Simu* p_agent, int& impacted_flag_change, float updating_in_min)
{
	// updating shorest path for vehicles passing through information node

	if (p_agent->diverted_flag >= 3)   // a vehicle is diverted once
		return 0;

	int current_link_no = p_agent->path_link_seq_no_vector[p_agent->m_current_link_seq_no];

	if (p_agent->m_current_link_seq_no > p_agent->path_link_seq_no_vector.size() - 5)  // very late step
		return 0;

	if (g_link_vector[current_link_no].capacity_reduction_map.find(p_agent->tau) != g_link_vector[current_link_no].capacity_reduction_map.end())
	{
		//simu_log_file << "trapped on incident link, return" << '\n';
		return -1;
	}

	int at = p_agent->at;
	int dest = p_agent->dest;
	int tau = p_agent->tau;

	int current_link_seq_no = p_agent->path_link_seq_no_vector[p_agent->m_current_link_seq_no];
	CLink* p_current_link = &(g_link_vector[current_link_seq_no]);


	//int link_sum0 = 0;
	int visual_distance_in_number_of_cells = 10;
	int impacted_link_no = -1;
	//for (int ls = p_agent->m_current_link_seq_no + 1; ls < min(p_agent->m_current_link_seq_no + visual_distance_in_number_of_cells, p_agent->path_link_seq_no_vector.size()); ls++)
	//{
	//    int link_no = p_agent->path_link_seq_no_vector[ls];
	//    //if (g_link_vector[link_no].capacity_reduction_map.find(tau) != g_link_vector[link_no].capacity_reduction_map.end())
	//    //{
	//    //    impacted_link_no = link_no;
	//    //    simu_log_file << " link capacity reduction " <<
	//    //        "at seq no" << ls << "," <<
	//    //        g_node_vector[g_link_vector[link_no].from_node_seq_no].node_id << "->" <<
	//    //        g_node_vector[g_link_vector[link_no].to_node_seq_no].node_id << " is detected" <<
	//    //        "for agent = " << p_agent->agent_id << '\n';

	//    //}
	//}


//    simu_log_file << "the near-sight path is impacted for agent =" << p_agent->agent_id << '\n';

	std::vector<int> path_link_vector;

	NetworkForSP* pNetwork = p_agent->p_RTNetwork;

	if (pNetwork == NULL)
		return -1;
	//simu_log_file << " rerouting planning at node " << g_node_vector[ g_link_vector[current_link_no].to_node_seq_no].node_id <<
	//    " for agent = " << p_agent->agent_id << '\n';

	/*simu_log_file << " current link status reduction= " << g_link_vector[current_link_no].capacity_reduction_map.size() << '\n';*/

	int required_updating_in_min = updating_in_min;

	//if (updating_in_min < pNetwork->updating_time_in_min + 1)
	//{
#pragma omp critical
	{
		pNetwork->optimal_backward_label_correcting_from_destination(0, this, required_updating_in_min, pNetwork->m_RT_dest_zone, pNetwork->m_RT_dest_node, impacted_link_no, 0);
	}
	//}
	//else  having shortest path data

	if (pNetwork->forwardtrace_shortest_path_tree(this, p_current_link->to_node_seq_no, path_link_vector) > 99999)
	{
		//simu_log_file << " no valid rerouting path (label cost >99999) for " <<
		//    p_agent->m_current_link_seq_no <<
		//    "at node " << g_node_vector[g_link_vector[p_agent->path_link_seq_no_vector[p_agent->m_current_link_seq_no]].to_node_seq_no].node_id <<
		//    " for agent = " << p_agent->agent_id << '\n';
	   //blocked routes
		return -1;
	}


	if (path_link_vector.size() <= 2)
	{
		//simu_log_file << " no valid rerouting path set for " <<
		//    p_agent->m_current_link_seq_no <<
		//    "for agent = " << p_agent->agent_id << '\n';
		return -1;
	}

	if (path_link_vector.size() >= 2)  // feasible rerouting and have not been informed by this sensor yet
	{

		int debug_flag = 0;
		int trace_agent_id = 609;
		////        if (p_agent->agent_id == trace_agent_id)
		//        {
		//            debug_flag = 1;
		//            simu_log_file << " current link " <<
		//                p_agent->m_current_link_seq_no <<
		//                "for agent = " << p_agent->agent_id << '\n';
		//        }
		//

			   // step 5: change the the current routes

		std::vector<int> extended_path_link_vector;
		std::vector<int> extended_path_link_arrival_time_vector;
		std::vector<int> extended_path_link_departure_time_vector;

		for (int l = 0; l <= p_agent->m_current_link_seq_no; l++)
		{
			extended_path_link_vector.push_back(p_agent->path_link_seq_no_vector[l]);
			extended_path_link_arrival_time_vector.push_back(p_agent->m_veh_link_arrival_time_in_simu_interval[l]);
			extended_path_link_departure_time_vector.push_back(p_agent->m_veh_link_departure_time_in_simu_interval[l]);

		}


		int current_time_t = p_agent->m_veh_link_arrival_time_in_simu_interval[p_agent->m_current_link_seq_no + 1];


		//expanding
		//simu_log_file << "agent id " << p_agent->agent_id << '\n';
		for (int nl = 0; nl < path_link_vector.size(); ++nl)  // arc a  // we do not exclude virtual link at the end here. as the output will exclude the virtual link in trajectory.csv
		{
			int link_seq_no = path_link_vector[nl];

			if (link_seq_no < 0)
			{
				int i_error = 1;
				break;
			}
			extended_path_link_vector.push_back(link_seq_no);
			extended_path_link_arrival_time_vector.push_back(-1);
			extended_path_link_departure_time_vector.push_back(-1);

			//    simu_log_file << g_node_vector[g_link_vector[link_seq_no].from_node_seq_no].node_id << ",";


		}
		//simu_log_file << " rerouting ends. " << '\n';

		p_agent->path_link_seq_no_vector.clear();
		p_agent->m_veh_link_arrival_time_in_simu_interval.clear();
		p_agent->m_veh_link_departure_time_in_simu_interval.clear();


		p_agent->diverted_flag = 1;

		for (int l = 0; l < extended_path_link_vector.size(); l++)
		{
			int link_no = extended_path_link_vector[l];

			if (g_link_vector[link_no].m_link_pedefined_capacity_map_in_sec.size() > 0 /*passing through incident area */ || g_link_vector[link_no].VDF_period[p_agent->tau].RT_allowed_use[p_agent->at] == false)
			{
				p_agent->diverted_flag = -1;
			}

			p_agent->path_link_seq_no_vector.push_back(extended_path_link_vector[l]);
			p_agent->m_veh_link_arrival_time_in_simu_interval.push_back(extended_path_link_arrival_time_vector[l]);
			p_agent->m_veh_link_departure_time_in_simu_interval.push_back(extended_path_link_departure_time_vector[l]);

		}

		p_agent->m_veh_link_arrival_time_in_simu_interval[p_agent->m_current_link_seq_no + 1] = current_time_t;
		// reupdating as this updating on the current link + 1 was performed before calling the RT path updating.

		int required_updating_in_min = updating_in_min;
		if (p_agent->diverted_flag == -1)
		{
			required_updating_in_min = -100; // will compute the path for sure
			pNetwork->optimal_backward_label_correcting_from_destination(0, this, required_updating_in_min, pNetwork->m_RT_dest_zone, pNetwork->m_RT_dest_node, impacted_link_no, 1);
			pNetwork->forwardtrace_shortest_path_tree(this, p_current_link->to_node_seq_no, path_link_vector);
		}

		return 1;
	}
	return 0;
}