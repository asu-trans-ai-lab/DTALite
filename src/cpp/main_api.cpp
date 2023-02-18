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
using std::cout;
using std::endl;
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
std::map<string, CCorridorInfo>  g_corridor_info_base0_map, g_corridor_info_SA_map;
std::map<int, CAnalysisDistrict>  g_district_info_base0_map, g_district_info_SA_map;
std::vector<COZone> g_zone_vector;
int g_related_zone_vector_size;

Assignment assignment;

std::map<string, CTMC_Corridor_Info> g_tmc_corridor_vector;
std::map<string, CInfoCell> g_info_cell_map;

#include "trip_generation.h"
#include "input.h"
#include "output.h"
#include "assignment.h"
#include "ODME.h"
#include "scenario_API.h"


void g_column_regeneration(Assignment& assignment)  // for RT
{		//step 3: column generation stage: find shortest path and update path cost of tree using TD travel time
	clock_t start_t, end_t, iteration_t;
	start_t = clock();
	dtalog.output() << endl;
	dtalog.output() << "Step 5: Column Re-Generation for Traffic Assignment..." << endl;
	dtalog.output() << "Total Column Re-Generation iteration: " << assignment.g_number_of_sensitivity_analysis_iterations << endl;
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


		dtalog.output() << endl;
		dtalog.output() << "SA: iteration number:" << iteration_number << endl;
		end_t = clock();
		iteration_t = end_t - start_t;
		dtalog.output() << "SA: CPU time: " << iteration_t / 1000.0 << " s" << endl;

		// step 3.1 update travel time and resource consumption
		clock_t start_t_lu = clock();

		//////detect induced_delay
		for (int i = 0; i < g_link_vector.size(); ++i)
		{
			for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
			{
				// used in travel time calculation
				//if (g_link_vector[i].VDF_period[tau].network_design_flag == 0)
				//{
				//	if (g_link_vector[i].VDF_period[tau].DOC > 3)  // dramatic delay, 3 is a cut off ratio/threshold for us to define large delay
				//	{
				//		g_link_vector[i].VDF_period[tau].network_design_flag = -2; // induced delay
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
			dtalog.output() << "Results:" << endl;
			for (int i = 0; i < g_link_vector.size(); ++i) {
				dtalog.output() << "link: " << g_node_vector[g_link_vector[i].from_node_seq_no].node_id << "-->"
					<< g_node_vector[g_link_vector[i].to_node_seq_no].node_id << ", "
					<< "flow count:" << g_link_vector[i].PCE_volume_per_period[0] << endl;
			}
		}

		end_t = clock();
		iteration_t = end_t - start_t_lu;
		// g_fout << "Link update with CPU time " << iteration_t / 1000.0 << " s; " << (end_t - start_t) / 1000.0 << " s" << endl;

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
				int network_copy_no = blk * assignment.g_number_of_memory_blocks + ProcessID;
				if (network_copy_no >= g_NetworkForSP_vector.size())
					continue;

				NetworkForSP* pNetwork = g_NetworkForSP_vector[network_copy_no];


				for (int o_node_index = 0; o_node_index < pNetwork->m_origin_node_vector.size(); ++o_node_index)
				{
					start_t_lc = clock();


					int origin_zone_seq_no = pNetwork->m_origin_zone_seq_no_vector[o_node_index]; // assigned nodes for computing


					if (assignment.g_AgentTypeVector[pNetwork->m_agent_type_no].real_time_information == 2 /*DMS*/
						&& assignment.zone_seq_no_2_info_mapping.find(origin_zone_seq_no) != assignment.zone_seq_no_2_info_mapping.end())
					{
						int idebug = 1;
					}

					if (assignment.g_AgentTypeVector[pNetwork->m_agent_type_no].real_time_information == 0)  // for non-information users, do not generate new path "columns"
						continue;  // skip the following shortest path codes

					if (assignment.g_AgentTypeVector[pNetwork->m_agent_type_no].real_time_information == 2 /*DMS infomation with physcial origin*/
						&& assignment.zone_seq_no_2_info_mapping.find(origin_zone_seq_no) == assignment.zone_seq_no_2_info_mapping.end()
						)
					{
						continue;  // skip the shortest path regeneration in the SA stage
					}

					pNetwork->optimal_label_correcting(ProcessID, &assignment, iteration_number, o_node_index, false, false, true);


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
			assignment.summary_file << "Sensitivity analysis stage" << endl;
			assignment.summary_file << ",Iteration,Avg Travel Time(min)" << endl;
		}

		double Avg_Travel_Time = total_system_travel_time / max(1, assignment.total_demand_volume);

		assignment.summary_file << iteration_number << "," << Avg_Travel_Time << "," << endl;
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
				g_link_vector[i].PCE_volume_per_period[tau] = 0;
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
					g_link_vector[i].PCE_volume_per_period[tau] = g_link_vector[i].PCE_volume_per_period[tau] * ratio;
					g_link_vector[i].person_volume_per_period[tau] = g_link_vector[i].person_volume_per_period[tau] * ratio;

					for (int at = 0; at < assignment.g_AgentTypeVector.size(); ++at)
						g_link_vector[i].person_volume_per_period_per_at[tau][at] *= ratio;
				}

				//for (int at = 0; at < assignment.g_AgentTypeVector.size(); ++at)
				//{
				//	for(int og = 0; og < assignment.g_number_of_analysis_districts; ++og)
				//		g_link_vector[i].person_volume_per_district_per_at[og][at] *= ratio;
				//}

			}
		}
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

	NetworkForSP* PointerMatrx[MAX_AGNETTYPES][MAX_TIMEPERIODS][MAX_MEMORY_BLOCKS] = { NULL };
	NetworkForSP* RTPointerMatrx[MAX_AGNETTYPES][MAX_TIMEPERIODS][MAX_MEMORY_BLOCKS] = { NULL };

	int computing_zone_count = 0;

	for (int at = 0; at < assignment.g_AgentTypeVector.size(); ++at)
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

					p_NetworkForSP->m_agent_type_no = at;
					p_NetworkForSP->m_tau = tau;
					p_NetworkForSP->PCE_ratio = assignment.g_AgentTypeVector[at].PCE;
					p_NetworkForSP->OCC_ratio = assignment.g_AgentTypeVector[at].OCC;

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


	dtalog.output() << "There are " << g_NetworkForSP_vector.size() << " SP networks in memory." << endl;
	dtalog.output() << "There are " << computing_zone_count << " agent type*zones to be computed in CPU." << endl;

}

void g_assign_RT_computing_tasks_to_memory_blocks(Assignment& assignment)
{
	//fprintf(g_pFileDebugLog, "-------g_assign_computing_tasks_to_memory_blocks-------\n");
	// step 2: assign node to thread
	dtalog.output() << "Step 2: Assigning RT info computing tasks to memory blocks..." << endl;

	int computing_zone_count = 0;

	int z_size = g_zone_vector.size();

	int at_size = assignment.g_AgentTypeVector.size();

	int tau_s_size = assignment.g_DemandPeriodVector.size();
	assignment.g_rt_network_pool = Allocate3DDynamicArray<NetworkForSP*>(z_size, at_size, tau_s_size);


	for (int at = 0; at < assignment.g_AgentTypeVector.size(); ++at)
	{
		for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); ++tau)
		{
			//assign all nodes to the corresponding thread
			for (int z = 0; z < g_zone_vector.size(); ++z)
			{
				NetworkForSP* p_NetworkForSP = new NetworkForSP();

				p_NetworkForSP->m_RT_dest_node = g_zone_vector[z].node_seq_no;
				p_NetworkForSP->m_RT_dest_zone = z;

				p_NetworkForSP->m_agent_type_no = at;
				p_NetworkForSP->m_tau = tau;
				p_NetworkForSP->PCE_ratio = assignment.g_AgentTypeVector[at].PCE;
				p_NetworkForSP->OCC_ratio = assignment.g_AgentTypeVector[at].OCC;

				computing_zone_count++;


				p_NetworkForSP->AllocateMemory(assignment.g_number_of_nodes, assignment.g_number_of_links, assignment.g_number_of_analysis_districts);


				assignment.g_rt_network_pool[z][at][tau] = p_NetworkForSP;  // assign real time computing network to the online colume;


				g_NetworkForRTSP_vector.push_back(p_NetworkForSP);
			}
		}
	}

	dtalog.output() << "There are " << g_NetworkForRTSP_vector.size() << " RTSP networks in memory." << endl;

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
					pNetwork->m_link_PCE_volume_array[i] = 0;
					pNetwork->m_link_person_volume_array[i] = 0;


					for(int og = 0; og < assignment.g_number_of_analysis_districts; og++)
					{
					pNetwork->m_link_person_volume_per_grid_array[i][og] = 0;
					}

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
			g_link_vector[i].PCE_volume_per_period[pNetwork->m_tau] += pNetwork->m_link_PCE_volume_array[i];
			g_link_vector[i].person_volume_per_period[pNetwork->m_tau] += pNetwork->m_link_person_volume_array[i];

			g_link_vector[i].person_volume_per_period_per_at[pNetwork->m_tau][pNetwork->m_agent_type_no] += pNetwork->m_link_person_volume_array[i];

			//for(int og = 0; og < assignment.g_number_of_analysis_districts; og ++)
			//{
			//g_link_vector[i].person_volume_per_district_per_at[og][pNetwork->m_agent_type_no] += pNetwork->m_link_person_volume_per_grid_array[i][og];
			//}
		}
	}
	// step 1: travel time based on VDF
}

void  CLink::calculate_dynamic_VDFunction(int inner_iteration_number, bool congestion_bottleneck_sensitivity_analysis_mode, int VDF_type_no)
{
	RT_waiting_time = 0; // reset RT_travel time for each end of simulation iteration 

	if (VDF_type_no == 0 || VDF_type_no == 1)  // BPR, QVDF
	{
		// for each time period
		for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
		{
			double link_volume_to_be_assigned = PCE_volume_per_period[tau] + VDF_period[tau].preload + VDF_period[tau].sa_volume;

			if (link_id == "10")
			{
				int idebug = 1;
			}
			if (VDF_period[tau].nlanes == 0)
			{
				int idebug;
				idebug = 1;
			}



			travel_time_per_period[tau] = VDF_period[tau].calculate_travel_time_based_on_QVDF(
				link_volume_to_be_assigned,
				this->model_speed, this->est_volume_per_hour_per_lane);


			VDF_period[tau].link_volume = link_volume_to_be_assigned;
//			VDF_period[tau].travel_time_per_iteration_map[inner_iteration_number] = travel_time_per_period[tau];
		}
	}
	else  // VDF_type_no = 2: 
	{
		// to do
		//calculate time-dependent travel time across all periods
		// and then assign the values to each periodd
		// for loop
		// for each time period
		// initialization 
		for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
		{
			VDF_period[tau].queue_length = 0;
			VDF_period[tau].arrival_flow_volume = PCE_volume_per_period[tau];
			VDF_period[tau].discharge_rate = VDF_period[tau].lane_based_ultimate_hourly_capacity * VDF_period[tau].nlanes;

			// dependend on the downstream blockagge states 

			VDF_period[tau].avg_waiting_time = 0;
		}
		//time-slot based queue evolution
		for (int tau = 1; tau < assignment.g_number_of_demand_periods; ++tau)
		{
			VDF_period[tau].queue_length = max(0.0f, VDF_period[tau - 1].queue_length + VDF_period[tau].arrival_flow_volume - VDF_period[tau].discharge_rate);

			if (inner_iteration_number == 1 && g_node_vector[this->from_node_seq_no].node_id == 1 &&
				g_node_vector[this->to_node_seq_no].node_id == 3)
			{
				int idebug = 1;

			}
		}

		// slot based total waiting time 
		for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
		{
			float prevailing_queue_length = 0;

			if (tau >= 1)
				prevailing_queue_length = VDF_period[tau - 1].queue_length;


			float total_waiting_time = (prevailing_queue_length + VDF_period[tau].queue_length) / 2.0 * VDF_period[tau].L;  // unit min
			// to do: 

			VDF_period[tau].avg_waiting_time = total_waiting_time / max(1.0f, VDF_period[tau].arrival_flow_volume);
			VDF_period[tau].avg_travel_time = VDF_period[tau].FFTT + VDF_period[tau].avg_waiting_time;
			VDF_period[tau].DOC = (prevailing_queue_length + VDF_period[tau].arrival_flow_volume) / max(0.01, VDF_period[tau].lane_based_ultimate_hourly_capacity * VDF_period[tau].nlanes);
			//to do:


			this->travel_time_per_period[tau] = VDF_period[tau].avg_waiting_time + VDF_period[tau].FFTT;

			// apply queue length at the time period to all the time slots included in this period
			for (int slot_no = assignment.g_DemandPeriodVector[tau].starting_time_slot_no; slot_no < assignment.g_DemandPeriodVector[tau].ending_time_slot_no; slot_no++)
			{
				this->est_queue_length_per_lane[slot_no] = VDF_period[tau].queue_length / max(1.0, this->number_of_lanes);
				this->est_avg_waiting_time_in_min[slot_no] = VDF_period[tau].avg_waiting_time;
				float number_of_hours = assignment.g_DemandPeriodVector[tau].time_period_in_hour;  // unit hour
				this->est_volume_per_hour_per_lane[slot_no] = VDF_period[tau].arrival_flow_volume / max(0.01f, number_of_hours) / max(1.0, this->number_of_lanes);
			}
		}
	}
}

double network_assignment(int assignment_mode, int column_generation_iterations, int column_updating_iterations, int ODME_iterations, int sensitivity_analysis_iterations, int simulation_iterations, int number_of_memory_blocks)
{
	clock_t start_t0, end_t0, total_t0;
	int signal_updating_iterations = 0;
	start_t0 = clock();
	// k iterations for column generation
	assignment.g_number_of_column_generation_iterations = column_generation_iterations;
	assignment.g_number_of_column_updating_iterations = column_updating_iterations;
	assignment.g_number_of_ODME_iterations = 100;
	assignment.g_number_of_sensitivity_analysis_iterations = sensitivity_analysis_iterations;
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
	g_ReadInformationConfiguration(assignment);

	//g_reload_timing_arc_data(assignment); // no need to load timing data from timing.csv
	g_ReadOutputFileConfiguration(assignment);
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
			g_add_new_virtual_connector_link(internal_from_node_seq_no, internal_to_node_seq_no, g_node_vector[i].agent_type_str, -1);
			// outgoing virtual connector
			g_add_new_virtual_connector_link(internal_to_node_seq_no, internal_from_node_seq_no, g_node_vector[i].agent_type_str, zone_seq_no);
			// result is that: we have a unique pair of node-zone access link in the overall network, but still carry agent_type_acess_map for agent types with access on this node-zone connector

		}

	}
	g_load_supply_side_scenario_file(assignment);

	//step 2: allocate memory and assign computing tasks
	g_assign_computing_tasks_to_memory_blocks(assignment); // static cost based label correcting

	// definte timestamps
	clock_t start_t, end_t, total_t;
	clock_t start_simu, end_simu, total_simu;
	start_t = clock();

	clock_t iteration_t, cumulative_lc, cumulative_cp, cumulative_lu;

	//step 3: column generation stage: find shortest path and update path cost of tree using TD travel time
	dtalog.output() << endl;
	dtalog.output() << "Step 4: Column Generation for Traffic Assignment..." << endl;
	dtalog.output() << "Total Column Generation iteration: " << assignment.g_number_of_column_generation_iterations << endl;


	for (int iteration_number = 0; iteration_number < max(1, assignment.g_number_of_column_generation_iterations); iteration_number++)
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
		double total_travel_distance = 0;
		double avg_speed = 0;

		// initialization at beginning of shortest path
		total_system_travel_time = update_link_travel_time_and_cost(iteration_number, total_travel_distance);

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
			dtalog.output() << "Results:" << endl;
			for (int i = 0; i < g_link_vector.size(); ++i) {
				dtalog.output() << "link: " << g_node_vector[g_link_vector[i].from_node_seq_no].node_id << "-->"
					<< g_node_vector[g_link_vector[i].to_node_seq_no].node_id << ", "
					<< "flow count:" << g_link_vector[i].PCE_volume_per_period[0] << endl;
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
				int network_copy_no = blk * assignment.g_number_of_memory_blocks + ProcessID;
				if (network_copy_no >= g_NetworkForSP_vector.size())
					continue;

				NetworkForSP* pNetwork = g_NetworkForSP_vector[network_copy_no];

				for (int o_node_index = 0; o_node_index < pNetwork->m_origin_node_vector.size(); ++o_node_index)
				{
					start_t_lc = clock();
					pNetwork->optimal_label_correcting(ProcessID, &assignment, iteration_number, o_node_index, false, false, false);


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

		// g_fout << "LC with CPU time " << cumulative_lc / 1000.0 << " s; " << endl;
		// g_fout << "column generation with CPU time " << cumulative_cp / 1000.0 << " s; " << endl;

		//****************************************//

		// last iteraion before performing signal timing updating
		double relative_gap = (total_system_travel_time - total_least_system_travel_time) / max(0.00001, total_least_system_travel_time);
		dtalog.output() << "iteration: " << iteration_number << ",systemTT: " << total_system_travel_time << endl;

		double CPU_Running_Time = 0;
		int number_of_agents = 0;
		double Avg_Travel_Time = total_system_travel_time / max(1, assignment.total_demand_volume);

		dtalog.output() << "iteration: " << iteration_number << ",systemTT: " << total_system_travel_time << ", least system TT:" <<
		    total_least_system_travel_time << ",gap=" << relative_gap << endl;

		if(iteration_number==0)
		{ 
			g_output_accessibility_result(assignment);
			assignment.summary_file << "Iteration, CPU Running Time, # of agents, Avg Travel Time(min),  Avg UE gap %" << endl;
		}

		assignment.summary_file << iteration_number << "," << CPU_Running_Time << "," << assignment.total_demand_volume << "," << Avg_Travel_Time << "," << relative_gap*100 << endl;
	}
	assignment.summary_file << "Step 4: Column Generation for Traffic Assignment" << endl;
	dtalog.output() << ",# of column generation (shortest path finding) iterations=, " << assignment.g_number_of_column_generation_iterations << endl;

	dtalog.output() << endl;

	// step 1.8: column updating stage: for given column pool, update volume assigned for each column
	dtalog.output() << "Step 4.2: Column Pool Updating" << endl;
	dtalog.output() << "Total Column Pool Updating iteration: " << column_updating_iterations << endl;

	start_t = clock();
	g_column_pool_optimization(assignment, column_updating_iterations);

	assignment.summary_file << "Step 6: column pool-based flow updating for traffic assignment " << endl;
	assignment.summary_file << ",# of flow updating iterations=," << column_updating_iterations << endl;

	dtalog.output() << endl;

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


	if (assignment.g_number_of_ODME_iterations >=1)
	{

		// record before odme volume
		for (int i = 0; i < g_link_vector.size(); ++i)
		{
			for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
			{
					g_link_vector[i].VDF_period[tau].volume_before_odme = g_link_vector[i].PCE_volume_per_period[tau] + g_link_vector[i].VDF_period[tau].preload;
			}
		}


		dtalog.output() << "Step 4.3: OD estimation for traffic assignment.." << endl;
		g_update_odme_volume_in_column_pool(assignment, 0);  // ODME before
		assignment.Demand_ODME(ODME_iterations, sensitivity_analysis_iterations);
		g_update_odme_volume_in_column_pool(assignment, 1);  // ODME after


		assignment.summary_file << "Step 4.3: OD estimation " << endl;
		assignment.summary_file << ",# of ODME_iterations=," << ODME_iterations << endl;

		// record after odme volume
		for (int i = 0; i < g_link_vector.size(); ++i)
		{
			for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
			{
				g_link_vector[i].VDF_period[tau].volume_after_odme = g_link_vector[i].PCE_volume_per_period[tau] + g_link_vector[i].VDF_period[tau].preload;
			}
		}

		dtalog.output() << endl;
	}

	// stage II sensitivity analysis stage 
	if (assignment.g_number_of_sensitivity_analysis_iterations >= 0)
	{

		int base_case_flag = 0; //base_case
		g_record_corridor_performance_summary(assignment, 0);

		// apply supply side scenarios
		for (int i = 0; i < g_link_vector.size(); ++i)
		{
			for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
			{
				// used in travel time calculation
				if (g_link_vector[i].VDF_period[tau].network_design_flag == -1 /*sa flag*/)  // we 
				{
					float old_lanes = g_link_vector[i].VDF_period[tau].nlanes;
					g_link_vector[i].VDF_period[tau].nlanes = g_link_vector[i].VDF_period[tau].sa_lanes_change; // apply the lane changes 
					g_link_vector[i].VDF_period[tau].BPR_period_capacity *= g_link_vector[i].VDF_period[tau].nlanes / max(0.00001, old_lanes);
					if (g_link_vector[i].VDF_period[tau].nlanes < 0.01)
						g_link_vector[i].VDF_period[tau].FFTT = 1440; // prevent free flow travel time under zero flow due to blockage
				}

			}
		}

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


				for (int at = 0; at < assignment.g_AgentTypeVector.size(); ++at)  //at
				{
					for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); ++tau)  //tau
					{
						p_column_pool = &(assignment.g_column_pool[from_zone_sindex][to_zone_sindex][at][tau]);
						if (p_column_pool->od_volume > 0 && fabs(local_scale_factor-1.0) > 0.001)
						{
							p_column_pool->od_volume = p_column_pool->od_volume* local_scale_factor;
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

		g_update_sa_volume_in_column_pool(assignment, 0);  // SA before

		g_column_regeneration(assignment);

		// stage II: right after column generation, we will enable the path modifications for DMS users 
		g_column_pool_route_modification(assignment, column_updating_iterations);


		// stage III: after DMS path modification, we enable further user equilibirum computing for real time info users by switching their routes 
		g_column_pool_optimization(assignment, assignment.g_number_of_sensitivity_analysis_iterations, true);
		
		g_update_sa_volume_in_column_pool(assignment, 1);  // SA after
	}

	g_classification_in_column_pool(assignment);
	

	if (simulation_iterations >= 1)
	{
		start_simu = clock();
		assignment.summary_file << "Step 5: traffic simulataion. " << endl;


		dtalog.output() << "Step 5: Simulation for traffic assignment.." << endl;
		assignment.STTrafficSimulation();
		end_simu = clock();
		total_simu = end_simu - end_simu;

		dtalog.output() << "CPU Running Time for traffic simulation: " << total_simu / 1000.0 << " s" << endl;
		dtalog.output() << endl;
	}

	end_t = clock();
	total_t = (end_t - start_t);
	dtalog.output() << "Done!" << endl;

	dtalog.output() << "CPU Running Time for the entire computing progcess: " << total_t / 1000.0 << " s" << endl;

	start_t = clock();

	//step 5: output simulation results of the new demand

	g_record_corridor_performance_summary(assignment, 1);

	g_output_assignment_result(assignment,0);

	if (assignment.g_subarea_shape_points.size() >= 3) // if there is a subarea defined.
	{
	g_output_assignment_result(assignment, 1);
	}

	if (assignment.assignment_mode == simulation_dta)
	{
	g_output_simulation_result(assignment);
	}

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
	dtalog.output() << "Output for assignment with " << assignment.g_number_of_column_generation_iterations << " iterations. Traffic assignment completes!" << endl;
	dtalog.output() << "CPU Running Time for outputting simulation results: " << total_t / 1000.0 << " s" << endl;

	dtalog.output() << "free memory.." << endl;


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
	dtalog.output() << "CPU Running Time for Entire Process: " << min << " min " << sec << " sec" << endl;
	dtalog.output() << "done." << endl;

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
		//simu_log_file << "trapped on incident link, return" << endl;
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
	//    //        "for agent = " << p_agent->agent_id << endl;

	//    //}
	//}


//    simu_log_file << "the near-sight path is impacted for agent =" << p_agent->agent_id << endl;

	std::vector<int> path_link_vector;

	NetworkForSP* pNetwork = p_agent->p_RTNetwork;

	if (pNetwork == NULL)
		return -1;
	//simu_log_file << " rerouting planning at node " << g_node_vector[ g_link_vector[current_link_no].to_node_seq_no].node_id <<
	//    " for agent = " << p_agent->agent_id << endl;

	/*simu_log_file << " current link status reduction= " << g_link_vector[current_link_no].capacity_reduction_map.size() << endl;*/

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
		//    " for agent = " << p_agent->agent_id << endl;
	   //blocked routes
		return -1;
	}


	if (path_link_vector.size() <= 2)
	{
		//simu_log_file << " no valid rerouting path set for " <<
		//    p_agent->m_current_link_seq_no <<
		//    "for agent = " << p_agent->agent_id << endl;
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
		//                "for agent = " << p_agent->agent_id << endl;
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
		//simu_log_file << "agent id " << p_agent->agent_id << endl;
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
		//simu_log_file << " rerouting ends. " << endl;

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

