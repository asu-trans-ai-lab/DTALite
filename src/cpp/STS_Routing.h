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


class VehicleScheduleNetworks {

public:
	int m_agent_type_no;
	int m_time_interval_size;  // 1440

	std::vector<CNode> m_node_vector;  // local copy of node vector, based on agent type and origin node
	std::map<int, float> g_passenger_link_profit;  // link with negative link

	void BuildNetwork(Assignment& assignment, int tau, int iteration)
	{

		int shortest_path_debugging_flag = 1;




		for (int i = 0; i < assignment.g_number_of_nodes; i++) //Initialization for all non-origin nodes
		{
			CNode node;  // create a node object 

			node.node_id = g_node_vector[i].node_id;
			node.node_seq_no = g_node_vector[i].node_seq_no;

			for (int j = 0; j < g_node_vector[i].m_outgoing_link_seq_no_vector.size(); j++)
			{

				int link_seq_no = g_node_vector[i].m_outgoing_link_seq_no_vector[j];

				if (assignment.g_LinkTypeMap[g_link_vector[link_seq_no].link_type].AllowAgentType(assignment.g_AgentTypeVector[m_agent_type_no].agent_type))  // only predefined allowed agent type can be considered
				{
					int from_node_seq_no = g_link_vector[link_seq_no].from_node_seq_no;
					node.m_outgoing_link_seq_no_vector.push_back(link_seq_no);

					g_passenger_link_profit[link_seq_no] = g_link_vector[link_seq_no].get_generalized_first_order_gradient_cost_of_second_order_loss_for_agent_type(tau, m_agent_type_no);

					if (shortest_path_debugging_flag && assignment.g_pFileDebugLog != NULL)
						fprintf(assignment.g_pFileDebugLog, "DP iteration %d: link %d->%d:  profit %.3f\n",
							iteration,
							g_node_vector[g_link_vector[link_seq_no].from_node_seq_no].node_id,
							g_node_vector[g_link_vector[link_seq_no].to_node_seq_no].node_id,
							g_passenger_link_profit[link_seq_no]);



					node.m_to_node_seq_no_vector.push_back(g_node_vector[i].m_to_node_seq_no_vector[j]);
				}

			}

			m_node_vector.push_back(node);

		}
	}


	//class for vehicle scheduling states
	class CVSState
	{
	public:
		int current_node_no;  // space dimension

		std::map<int, int> passenger_service_state;  // passenger means link with negative costs

		std::vector<int> m_visit_link_sequence;  // store link sequence

		int m_vehicle_capacity;

		float LabelCost;  // with LR price
		float LabelTime;   // sum of travel time up to now, arrival time at node


		CVSState()
		{
			LabelTime = 0;
			LabelCost = 0;
			m_vehicle_capacity = 0;
		}

		void Copy(CVSState* pSource)
		{
			current_node_no = pSource->current_node_no;
			passenger_service_state.clear();
			passenger_service_state = pSource->passenger_service_state;


			m_visit_link_sequence = pSource->m_visit_link_sequence;
			m_vehicle_capacity = pSource->m_vehicle_capacity;
			LabelCost = pSource->LabelCost;
			LabelTime = pSource->LabelTime;
		}
		int GetPassengerLinkServiceState(int link_no)
		{
			if (passenger_service_state.find(link_no) != passenger_service_state.end())
				return passenger_service_state[link_no];  // 1 or 2
			else
				return 0;
		}



		std::string generate_string_key()
		{

			stringstream s;

			s << "n";
			s << current_node_no;  // space key
			for (std::map<int, int>::iterator it = passenger_service_state.begin(); it != passenger_service_state.end(); ++it)
			{
				s << "_";

				s << it->first << "[" << it->first << "]";

			}
			string converted(s.str());
			return converted;

		}

		bool operator<(const CVSState& other) const
		{
			return LabelCost < other.LabelCost;
		}

	};

	class C_time_indexed_state_vector
	{
	public:
		int current_time;


		std::vector<CVSState> m_VSStateVector;

		std::map<std::string, int> m_state_map;

		void Reset()
		{
			current_time = 0;
			m_VSStateVector.clear();
			m_state_map.clear();
		}

		int m_find_state_index(std::string string_key)
		{

			if (m_state_map.find(string_key) != m_state_map.end())
			{
				return m_state_map[string_key];
			}
			else
				return -1;  // not found

		}

		void update_state(CVSState new_element)
		{
			std::string string_key = new_element.generate_string_key();//if it is new, string is n100, no state index
			int state_index = m_find_state_index(string_key);

			if (state_index == -1)  // no such state at this time
			{
				// add new state
				state_index = m_VSStateVector.size();
				m_VSStateVector.push_back(new_element);
				m_state_map[string_key] = state_index;
			}
			else
			{//DP
				if (new_element.LabelCost < m_VSStateVector[state_index].LabelCost)
				{
					m_VSStateVector[state_index].Copy(&new_element);
				}

			}

		}

		void Sort()
		{
			std::sort(m_VSStateVector.begin(), m_VSStateVector.end());

			m_state_map.clear(); // invalid
		}

		void SortAndCleanEndingState(int BestKValue)
		{
			if (m_VSStateVector.size() > 2 * BestKValue)
			{
				std::sort(m_VSStateVector.begin(), m_VSStateVector.end());

				m_state_map.clear(); // invalid
				m_VSStateVector.erase(m_VSStateVector.begin() + BestKValue, m_VSStateVector.end());
			}
		}

		float GetBestValue()
		{
			// LabelCost not PrimalCost when sorting
			std::sort(m_VSStateVector.begin(), m_VSStateVector.end());

			if (m_VSStateVector.size() >= 1)
			{
				std::string state_str = m_VSStateVector[0].generate_string_key();

				return m_VSStateVector[0].LabelCost;

			}
			else
				return MAX_LABEL_COST;
		}

		std::vector<int> GetBestLinkSequence()
		{
			std::vector <int> link_sequence;
			// LabelCost not PrimalCost when sorting
			std::sort(m_VSStateVector.begin(), m_VSStateVector.end());

			if (m_VSStateVector.size() >= 1)
			{
				std::string state_str = m_VSStateVector[0].generate_string_key();

				if (m_VSStateVector[0].m_visit_link_sequence.size() > 0)
					return m_VSStateVector[0].m_visit_link_sequence;
				else
					return link_sequence;

			}

			return link_sequence;
		}

	};

	//vehicle state at time t

	// for collecting the final feasible states accesible to the depot
	C_time_indexed_state_vector g_ending_state_vector;

	C_time_indexed_state_vector** g_time_dependent_state_vector;  // label cost vector [i,t,w]


	void AllocateVSNMemory(int number_of_nodes)
	{
		g_time_dependent_state_vector = AllocateDynamicArray <C_time_indexed_state_vector>(number_of_nodes, m_time_interval_size);  //1
	}

	~VehicleScheduleNetworks()
	{
		DeallocateDynamicArray(g_time_dependent_state_vector, g_node_vector.size(), m_time_interval_size);
	}


	std::vector<int> g_optimal_time_dependenet_dynamic_programming(
		int origin_node,
		int destination_node,
		int vehicle_capacity,
		int demand_time_period_no,
		//maximum choose
		int BestKSize)
		// time-dependent label correcting algorithm with double queue implementation
	{
		int arrival_time = m_time_interval_size - 1;  // restrict the search range.

		int t = 0;
		//step 2: Initialization for origin node at the preferred departure time, at departure time
		for (int i = 0; i < assignment.g_number_of_nodes; i++)
		{
			g_time_dependent_state_vector[i][0].Reset();

		}
		g_ending_state_vector.Reset();

		CVSState element;

		element.current_node_no = origin_node;
		g_time_dependent_state_vector[origin_node][0].update_state(element);


		// step 3: //dynamic programming
		for (t = 0; t < arrival_time; t++)  //first loop: time
		{

			for (int n = 0; n < m_node_vector.size(); n++)
			{
				// step 1: sort m_VSStateVector by labelCost for scan best k elements in step2
				g_time_dependent_state_vector[n][t].Sort();

				// step 2: scan the best k elements
				for (int w_index = 0; w_index < min(BestKSize, g_time_dependent_state_vector[n][t].m_VSStateVector.size()); w_index++)
				{
					CVSState* pElement = &(g_time_dependent_state_vector[n][t].m_VSStateVector[w_index]);

					int from_node = pElement->current_node_no;

					// step 2.1 link from node to toNode
					for (int i = 0; i < m_node_vector[from_node].m_outgoing_link_seq_no_vector.size(); i++)
					{
						int link_seq_no = m_node_vector[from_node].m_outgoing_link_seq_no_vector[i];
						int to_node = g_link_vector[link_seq_no].to_node_seq_no;

						float new_time = pElement->LabelTime + g_link_vector[link_seq_no].travel_time_per_period[demand_time_period_no];
						int new_time_int = max(pElement->LabelTime + 1, (int)(new_time + 0.5));  // move at least one time step further

						// step 2.2. check feasibility of node type with the current element
						if (new_time <= arrival_time)
						{

							// skip scanning when the origin/destination nodes arrival time is out of time window
							//feasible state transitions
								// do not need waiting
							CVSState new_element;
							new_element.Copy(pElement);

							new_element.current_node_no = to_node;

							new_element.LabelCost += g_link_vector[link_seq_no].travel_time_per_period[demand_time_period_no] / 60.0 * assignment.g_AgentTypeVector[m_agent_type_no].value_of_time; // 60.0 is to convert hour to 60 min as VOT is denoted as dollars per hour
							if (g_passenger_link_profit.find(link_seq_no) != g_passenger_link_profit.end())
							{
								new_element.LabelCost += g_passenger_link_profit[link_seq_no];// + negative cost
								new_element.passenger_service_state[link_seq_no] = 1;  // mark carry status
								new_element.m_vehicle_capacity -= 1;

							}


							new_element.m_visit_link_sequence.push_back(link_seq_no);
							g_time_dependent_state_vector[to_node][new_time_int].update_state(new_element);

							if (to_node == destination_node)
							{

								//time window of destination_node
								if (new_time < arrival_time)
								{
									g_ending_state_vector.update_state(new_element);
									g_ending_state_vector.SortAndCleanEndingState(BestKSize);
								}
							}
						}

					}
				}
			}  // for all nodes
		} // for all time t


	// no backf
		return g_ending_state_vector.GetBestLinkSequence();

	}

};

//
//void g_column_pool_real_time_updating(Assignment& assignment, int column_updating_iterations)
//{
//	int agent_type_size = assignment.g_AgentTypeVector.size();
//	int zone_size = g_zone_vector.size();
//	int demand_period_size = assignment.g_DemandPeriodVector.size();
//
//	CColumnVector* p_column_pool;
//
//	float path_toll = 0;
//	float path_distance = 0;
//	float path_travel_time = 0;
//	float path_delay = 0;
//	float path_FF_travel_time = 0;
//	float time_stamp = 0;
//
//	std::map<int, CColumnPath>::iterator it, it_begin, it_end;
//
//	if (assignment.major_path_volume_threshold > 0.00001)  // performing screening of path flow pattern
//	{
//
//		//initialization 
//		bool b_subarea_mode = false;
//
//		int number_of_links = g_link_vector.size();
//		for (int i = 0; i < number_of_links; ++i)
//		{
//			for (int tau = 0; tau < demand_period_size; ++tau)
//			{
//				// used in travel time calculation
//				g_link_vector[i].background_total_volume_for_all_agent_types_per_period[tau] = 0;
//			}
//
//			if (g_node_vector[g_link_vector[i].from_node_seq_no].subarea_id >= 1 && g_node_vector[g_link_vector[i].to_node_seq_no].node_id >= 1)
//			{
//
//				g_link_vector[i].subarea_id = g_node_vector[g_link_vector[i].from_node_seq_no].subarea_id;
//				b_subarea_mode = true;
//			}
//			else
//				g_link_vector[i].subarea_id = 0;
//
//		}
//
//
//		/// <summary>  screening the path flow pattern
//		for (int orig = 0; orig < zone_size; ++orig)
//		{
//
//			for (int at = 0; at < agent_type_size; ++at)
//			{
//				for (int dest = 0; dest < zone_size; ++dest)
//				{
//					for (int tau = 0; tau < demand_period_size; ++tau)
//					{
//						p_column_pool = &(assignment.g_column_pool[orig][dest][at][tau]);
//						if (p_column_pool->od_volume[assignment.active_scenario_index] > 0)
//						{
//							// scan through the map with different node sum for different continuous paths
//							it_begin = p_column_pool->path_node_sequence_map.begin();
//							it_end = p_column_pool->path_node_sequence_map.end();
//
//							for (it = it_begin; it != it_end; ++it)
//							{
//								int subarea_output_flag = 0;
//								if (b_subarea_mode == true)
//								{
//									int insubarea_flag = 0;
//
//									for (int nl = 0; nl < it->second.m_link_size; ++nl)  // arc a
//									{
//										int link_seq_no = it->second.path_link_vector[nl];
//
//										if (g_link_vector[link_seq_no].subarea_id >= 1)
//										{
//											insubarea_flag = 1;
//											break;
//										}
//
//									}
//									// 
//									if (insubarea_flag && it->second.path_volume > assignment.major_path_volume_threshold)
//									{
//										subarea_output_flag = 1;
//									}
//
//								}
//								else
//								{
//									if (it->second.path_volume > assignment.major_path_volume_threshold)
//										subarea_output_flag = 1;
//
//								}
//								if (subarea_output_flag == 0)
//								{
//									it->second.subarea_output_flag = 0;  // disable the output of this column into route_assignment.csv
//
//									for (int nl = 0; nl < it->second.m_link_size; ++nl)  // arc a
//									{
//										int link_seq_no = it->second.path_link_vector[nl];
//										g_link_vector[link_seq_no].background_total_volume_for_all_agent_types_per_period[tau] += it->second.path_volume;
//									}
//								}
//
//							}
//						}
//					}
//					/// </summary>
//					/// <param name="assignment"></param>
//				}
//
//			}
//		}
//
//		/// output background_link_volume.csv
//		dtalog.output() << "writing link_performance.csv.." << endl;
//
//		int b_debug_detail_flag = 0;
//		FILE* g_pFileLinkMOE = nullptr;
//
//		fopen_ss(&g_pFileLinkMOE, "link_background_volume.csv", "w");
//		if (!g_pFileLinkMOE)
//		{
//			dtalog.output() << "File link_background_volume.csv cannot be opened." << endl;
//			g_program_stop();
//		}
//		else
//		{
//			fprintf(g_pFileLinkMOE, "link_id,from_node_id,to_node_id,from_cell_code,time_period,volume,background_volume,major_path_volume,ratio_of_major_path_flow,geometry,");
//
//			fprintf(g_pFileLinkMOE, "notes\n");
//
//			//Initialization for all nodes
//			for (int i = 0; i < g_link_vector.size(); ++i)
//			{
//				// virtual connectors
//				if (g_link_vector[i].link_type == -1)
//					continue;
//
//				for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
//				{
//					double volume = g_link_vector[i].total_volume_for_all_agent_types_per_period[tau] + g_link_vector[i].VDF_period[tau].preload;
//					double major_path_link_volume = g_link_vector[i].total_volume_for_all_agent_types_per_period[tau] + g_link_vector[i].VDF_period[tau].preload - g_link_vector[i].background_total_volume_for_all_agent_types_per_period[tau];
//					double ratio = major_path_link_volume / max(volume, 0.000001);
//
//					if (volume < 0.0000001)
//						ratio = -1;
//					fprintf(g_pFileLinkMOE, "%s,%d,%d,%s,%s,%.3f,%.3f,%.3f,%.3f,\"%s\",",
//						g_link_vector[i].link_id.c_str(),
//						g_node_vector[g_link_vector[i].from_node_seq_no].node_id,
//						g_node_vector[g_link_vector[i].to_node_seq_no].node_id,
//						g_node_vector[g_link_vector[i].from_node_seq_no].cell_str.c_str(),
//						assignment.g_DemandPeriodVector[tau].time_period.c_str(),
//						g_link_vector[i].total_volume_for_all_agent_types_per_period[tau] + g_link_vector[i].VDF_period[tau].preload,
//						g_link_vector[i].background_total_volume_for_all_agent_types_per_period[tau],
//						major_path_link_volume,
//						ratio,
//						g_link_vector[i].geometry.c_str());
//					fprintf(g_pFileLinkMOE, "\n");
//
//				}
//
//			}
//
//			fclose(g_pFileLinkMOE);
//		}
//
//
//	} // end of path flow pattern screening 
//	dtalog.output() << "writing data for " << zone_size << "  zones " << endl;
//
//	for (int orig = 0; orig < zone_size; ++orig)
//	{
//		if (g_zone_vector[orig].zone_id % 100 == 0)
//			dtalog.output() << "o zone id =  " << g_zone_vector[orig].zone_id << endl;
//
//		for (int at = 0; at < agent_type_size; ++at)
//		{
//			for (int dest = 0; dest < zone_size; ++dest)
//			{
//				for (int tau = 0; tau < demand_period_size; ++tau)
//				{
//					p_column_pool = &(assignment.g_column_pool[orig][dest][at][tau]);
//					if (p_column_pool->od_volume[assignment.active_scenario_index] > 0 ||
//						(assignment.zone_seq_no_2_info_mapping.find(orig) != assignment.zone_seq_no_2_info_mapping.end()
//							&& assignment.g_AgentTypeVector[at].real_time_information >= 1)
//						)
//
//					{
//
//						int information_type = 0;
//
//						if (assignment.zone_seq_no_2_info_mapping.find(orig) != assignment.zone_seq_no_2_info_mapping.end()
//							&& assignment.g_AgentTypeVector[at].real_time_information >= 1)
//						{
//							p_column_pool->od_volume[assignment.active_scenario_index] = 1;  // reset the volume as 1 to enable visualization 
//
//						}
//
//
//	for (int at = 0; at < assignment.g_AgentTypeVector.size(); at++)  //m
//		if (assignment.g_AgentTypeVector[at].flow_type == 2)  // we only take into account the MAS generated agent volume into the link volume and link resource in this second optimization stage.
//		{
//			for (int o = 0; o < g_zone_vector.size(); o++)  // o
//				for (int d = 0; d < g_zone_vector.size(); d++) //d
//					for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); tau++)  //tau
//						if (assignment.g_column_pool[o][d][at][tau].od_volume[active_scenario_index] > 0)
//						{
//							for (int ai = 0; ai < assignment.g_column_pool[o][d][at][tau].discrete_agent_path_vector.size(); ai++)
//							{
//								CAgentPath agent_path = assignment.g_column_pool[o][d][at][tau].discrete_agent_path_vector[ai];
//								{
//									// internal step 1: test shortest path travel time 
//									g_RoutingNetwork.m_origin_node = agent_path.o_node_no;
//									g_RoutingNetwork.m_agent_type_no = at;
//									g_RoutingNetwork.tau = tau;
//									float route_trip_time = g_RoutingNetwork.optimal_label_correcting(assignment, 0, agent_path.d_node_no, true);
//
//									VehicleScheduleNetworks vsn;
//									vsn.m_agent_type_no = at;
//									vsn.m_time_interval_size = max(max(route_trip_time * 1.5, route_trip_time + 10), max(assignment.g_DemandPeriodVector[tau].get_time_horizon_in_min() * 1.5, assignment.g_DemandPeriodVector[tau].get_time_horizon_in_min() + 10));
//									vsn.AllocateVSNMemory(assignment.g_number_of_nodes);
//									vsn.BuildNetwork(assignment, tau, iteration_number);
//
//									assignment.g_column_pool[o][d][at][tau].discrete_agent_path_vector[ai].path_link_sequence = vsn.g_optimal_time_dependenet_dynamic_programming(agent_path.o_node_no, agent_path.d_node_no, 0, tau, 10);
//
//								}
//							}
//						}
//
//		}
//
//}
