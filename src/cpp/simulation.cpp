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
#include "DTA.h"

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

extern void g_rt_info_column_generation(Assignment* p_assignment, float current_time_in_min, int recording_flag);
extern double g_get_random_ratio();
extern void g_output_agent_csv(Assignment& assignment);

void Assignment::AllocateLinkMemory4Simulation()
{
	g_number_of_simulation_intervals = (g_LoadingEndTimeInMin - g_LoadingStartTimeInMin + simulation_discharge_period_in_min) * 60 / number_of_seconds_per_interval + 2;

	g_number_of_intervals_in_sec = (g_LoadingEndTimeInMin - g_LoadingStartTimeInMin + simulation_discharge_period_in_min) * 60;

	dtalog.output() << "[DATA INFO] LoadingStartTimeInMin = " << g_LoadingStartTimeInMin << '\n';
	dtalog.output() << "[DATA INFO] g_LoadingStartTimeInMin = " << g_LoadingEndTimeInMin << '\n';
	dtalog.output() << "[DATA INFO] number_of_simulation_intervals = " << g_number_of_simulation_intervals << '\n';
	dtalog.output() << "[DATA INFO] number_of_simu intervals in sec = " << g_number_of_intervals_in_sec << '\n';

	g_number_of_loading_intervals_in_sec = (g_LoadingEndTimeInMin - g_LoadingStartTimeInMin) * 60;

	g_number_of_intervals_in_min = (int)(g_number_of_simulation_intervals / number_of_simu_intervals_in_min + 1);
	// add + 120 as a buffer
	g_number_of_in_memory_simulation_intervals = g_number_of_simulation_intervals;

	dtalog.output() << "[STATUS INFO] allocate 2D dynamic memory for vector LinkOutFlowCapacity..." << '\n';

	m_LinkOutFlowCapacity = Allocate2DDynamicArray <float>(g_number_of_links, g_number_of_intervals_in_sec);  //1
	m_LinkOutFlowState = Allocate2DDynamicArray <int>(g_number_of_links, g_number_of_intervals_in_sec);  //1
	// discharge rate per simulation time interval
	dtalog.output() << "[STATUS INFO] allocate 2D dynamic memory for vector m_LinkCumulativeArrivalVector..." << '\n';
	m_LinkCumulativeArrivalVector = Allocate2DDynamicArray <float>(g_number_of_links, g_number_of_intervals_in_min);  //2

	dtalog.output() << "[STATUS INFO] allocate 2D dynamic memory for vector m_LinkCumulativeDepartureVector..." << '\n';
	m_LinkCumulativeDepartureVector = Allocate2DDynamicArray <float>(g_number_of_links, g_number_of_intervals_in_min);  //3

	m_link_CA_count = Allocate1DDynamicArray <float>(g_number_of_links);
	m_link_CD_count = Allocate1DDynamicArray <float>(g_number_of_links);

	dtalog.output() << "[STATUS INFO] allocate 2D dynamic memory for vector  m_link_TD_waiting_time..." << '\n';
	m_link_TD_waiting_time = Allocate2DDynamicArray <float>(g_number_of_links, g_number_of_intervals_in_min); //5

	m_link_total_waiting_time_vector.clear();
	for (int i = 0; i < g_number_of_links; ++i)
	{
		m_link_total_waiting_time_vector.push_back(0.0);
	}

	g_AgentTDListMap.clear();

	dtalog.output() << "[STATUS INFO] initializing time dependent capacity data..." << '\n';

#pragma omp parallel for
	for (int i = 0; i < g_number_of_links; ++i)
	{
		float cap_count = 0;
		float discharge_rate_per_sec = g_link_vector[i].lane_capacity * g_link_vector[i].number_of_lanes_si[assignment.active_scenario_index] / 3600.0;
		float discharge_rate_after_loading = 20 * discharge_rate_per_sec;  ///to collect travel time statistics * 10 times of capacity to discharge all flow */ at the end of simulation time interval

		if (i == 24)
		{
			int debug = 1;
		}

		unsigned int RandomSeed = 101;
		float residual;
		float random_ratio = 0;

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
				m_LinkOutFlowCapacity[i][t] = 1 * g_link_vector[i].number_of_lanes_si[assignment.active_scenario_index];
			else
			{
				residual = OutFlowRate - (int)(OutFlowRate);
				//RandomSeed is automatically updated.
				RandomSeed = (LCG_a * RandomSeed + LCG_c) % LCG_M;
				random_ratio = float(RandomSeed) / LCG_M;

				if (random_ratio < residual)
					m_LinkOutFlowCapacity[i][t] = (int)(OutFlowRate)+1;
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
			m_link_TD_waiting_time[i][t] = 0;
			m_LinkOutFlowState[i][t] = 1;
		}


	}

	// for each link with  for link type code is 's'
	//reset the time-dependent capacity to zero
	//go through the records defined in timing_arc file
	//only enable m_LinkOutFlowCapacity[l][t] for the timestamp in the time window per time interval and reset the link TD travel time.


	for (unsigned l = 0; l < g_link_vector.size(); ++l)
	{
		if (l == 40949)
		{
			int i_bebug = 1;
		}
		if (g_link_vector[l].timing_arc_flag == true)  // with timing data
		{

			if (g_link_vector[l].timing_arc_flag)
			{
				// reset for signalized links (not freeway links as type code != 'f' for the case of freeway workzones)
				// only for the loading period
				for (int t = 0; t < g_number_of_loading_intervals_in_sec; ++t)
				{
					m_LinkOutFlowCapacity[l][t] = 0;
					m_LinkOutFlowState[l][t] = 0;

				}
			}

			int number_of_cycles = (g_LoadingEndTimeInMin - g_LoadingStartTimeInMin) * 60 / max(1.0f, g_link_vector[l].VDF_period[0].cycle_length);  // unit: seconds;

			int cycle_length = g_link_vector[l].VDF_period[0].cycle_length;
			int start_green_time = g_link_vector[l].VDF_period[0].start_green_time;
			int end_green_time = g_link_vector[l].VDF_period[0].end_green_time;

			if (end_green_time < start_green_time)
			{
				end_green_time += cycle_length;  // consider a looped end green time notation, e.g. 60-10 for cl = 100, then end green time should be 110.
			}
			dtalog.output() << "[DATA INFO] signal timing data: link: cycle_length = " << cycle_length <<
				",start_green_time = " << start_green_time <<
				",end_green_time = " << end_green_time <<
				'\n';

			for (int cycle_no = 0; cycle_no < number_of_cycles; ++cycle_no)
			{
				int count = 0;

				// relative time horizon
				for (int t = cycle_no * cycle_length + start_green_time; t <= cycle_no * cycle_length + end_green_time; t += 1)
				{
					// activate capacity for this time duration
					m_LinkOutFlowCapacity[l][t] = g_link_vector[l].saturation_flow_rate * g_link_vector[l].number_of_lanes_si[assignment.active_scenario_index] / 3600.0;

					if (t % 2 == 0)
						m_LinkOutFlowCapacity[l][t] = 1 * g_link_vector[l].number_of_lanes_si[assignment.active_scenario_index];
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
		} if (g_link_vector[l].m_link_pedefined_capacity_map_in_sec.size() > 0)
		{


			for (int t = 0; t < g_number_of_loading_intervals_in_sec; ++t)
			{
				if (g_link_vector[l].m_link_pedefined_capacity_map_in_sec.find(t) != g_link_vector[l].m_link_pedefined_capacity_map_in_sec.end())
				{
					// predefined outlflow capacity
					m_LinkOutFlowCapacity[l][t] = 0;
					m_LinkOutFlowState[l][t] = 0;
				}
			}
		}

		//end f
	}

	dtalog.output() << "[STATUS INFO] End of initializing time dependent capacity data." << '\n';
}

void Assignment::DeallocateLinkMemory4Simulation()
{
	// g_fout << "deallocate 2D dynamic memory m_LinkOutFlowCapacity..." << '\n';
	if(m_LinkOutFlowCapacity)
	    Deallocate2DDynamicArray(m_LinkOutFlowCapacity , g_number_of_links, g_number_of_intervals_in_sec);  //1

	if (m_LinkOutFlowState)
	    Deallocate2DDynamicArray(m_LinkOutFlowState, g_number_of_links, g_number_of_intervals_in_sec);  //1
	// g_fout << "deallocate 2D dynamic memory m_LinkCumulativeArrivalVector..." << '\n';
	if (m_LinkCumulativeArrivalVector)
		Deallocate2DDynamicArray(m_LinkCumulativeArrivalVector, g_number_of_links, g_number_of_intervals_in_min);  //2
	// g_fout << "deallocate 2D dynamic memory m_LinkCumulativeDepartureVector..." << '\n';
	if (m_LinkCumulativeDepartureVector)
		Deallocate2DDynamicArray(m_LinkCumulativeDepartureVector, g_number_of_links, g_number_of_intervals_in_min); //3

	if (m_link_CA_count)
		Deallocate1DDynamicArray(m_link_CA_count, g_number_of_links);

	if (m_link_CD_count)
		Deallocate1DDynamicArray(m_link_CD_count, g_number_of_links);

	if (m_link_TD_waiting_time)
		Deallocate2DDynamicArray(m_link_TD_waiting_time, g_number_of_links, g_number_of_intervals_in_min); //4
}

void Assignment::STTrafficSimulation()
{
	g_agent_simu_vector.clear();

	ofstream simu_log_file;

	//given p_agent->path_link_seq_no_vector path link sequence no for each agent
	int TotalCumulative_Arrival_Count = 0;
	int TotalCumulative_rerouting_count = 0;
	int TotalCumulative_impact_count = 0;
	int TotalCumulative_Departure_Count = 0;
	double TotalVehicleMileTraveled = 0;
	double TotalVehicleHourTraveled = 0;

	clock_t start_t;
	start_t = clock();

	AllocateLinkMemory4Simulation();

	int mode_type_size = g_ModeTypeVector.size();
	int zone_size = g_zone_vector.size();
	int demand_period_size = g_DemandPeriodVector.size();
	int system_information_activate_time_in_min = 999999;
	int system_information_activate_time_in_simu_interval = 99999999;


	for (int li = 0; li < g_link_vector.size(); ++li)
	{
		CLink* p_link = &(g_link_vector[li]);
		{

			system_information_activate_time_in_min = min(g_link_vector[li].dynamic_link_event_start_time_in_min, system_information_activate_time_in_min);

			for (int at = 0; at < mode_type_size; ++at)
			{
				for (int tau = 0; tau < demand_period_size; ++tau)
				{
					if (p_link->AllowModeType(g_ModeTypeVector[at].mode_type, tau, active_scenario_index) == false)
					{
						p_link->VDF_period[tau].RT_allowed_use[at] = false;  // follow the general rule of allow_uses, such HOV, special lanes
					}
					if (p_link->SA_AllowModeType(g_ModeTypeVector[at].mode_type, tau) == false)
					{
						p_link->VDF_period[tau].SA_allowed_use[at] = false;  // follow the general rule of allow_uses, such HOV, special lanes
					}
				}
			}

			p_link->EntranceQueue.clear();
			p_link->ExitQueue.clear();
		}

	}


	dtalog.output() << "[STATUS INFO] system-wide information activate time = " << system_information_activate_time_in_min << " in min" << '\n';
	assignment.summary_file << ", system - wide information activate time = " << system_information_activate_time_in_min << " in min" << '\n';

	if (system_information_activate_time_in_min < 2000)
	{
		system_information_activate_time_in_simu_interval = (system_information_activate_time_in_min - g_LoadingStartTimeInMin) * 60 / number_of_seconds_per_interval;
		//	dtalog.output() << "system-wide information activate time in simu interval = " << system_information_activate_time_in_simu_interval << '\n';
	}

	CColumnVector* p_column_pool;
	float path_toll = 0;
	float path_distance = 0;
	float path_travel_time = 0;
	float time_stamp = 0;

	int trace_agent_id = 609;  //default can be -1
	int trace_link_no = -1;  //default can be -1
	int trace_node_no = -100;  //default can be -1

	double previous_headway_residual = 0;
	std::map<int, CColumnPath>::iterator it, it_begin, it_end;
	int global_path_no = 0;
	for (int orig = 0; orig < zone_size; ++orig)
	{
		int from_zone_sindex = g_zone_vector[orig].sindex;
		if (from_zone_sindex == -1)
			continue;


		for (int at = 0; at < mode_type_size; ++at)
		{
			for (int dest = 0; dest < zone_size; ++dest)
			{
				int to_zone_sindex = g_zone_vector[dest].sindex;
				if (to_zone_sindex == -1)
					continue;
				for (int tau = 0; tau < demand_period_size; ++tau)
				{
					p_column_pool = &(assignment.g_column_pool[from_zone_sindex][to_zone_sindex][at][tau]);



					if (p_column_pool->od_volume[assignment.active_scenario_index] > 0)
					{

						// scan through the map with different node sum for different continuous paths
						it_begin = p_column_pool->path_node_sequence_map.begin();
						it_end = p_column_pool->path_node_sequence_map.end();



						for (it = it_begin; it != it_end; ++it)
						{
							it->second.global_path_no = global_path_no++;
							if (it->second.global_path_no == 55)
							{
								int debug = 1;
							}

							path_toll = 0;
							path_distance = 0;
							path_travel_time = 0;

							int VehicleSize = (it->second.path_volume + 0.5);  // total number of vehicles
							CDeparture_time_Profile dep_time_element;
							if (p_column_pool->departure_time_profile_no < assignment.g_DepartureTimeProfileVector.size())
							{
								dep_time_element = assignment.g_DepartureTimeProfileVector[p_column_pool->departure_time_profile_no];
							}
							else
							{
								dtalog.output() << "[ERROR] error in departure_time_profile_no=  " << p_column_pool->departure_time_profile_no << '\n';
								dtalog.output() << "[DATA INFO] size of g_DepartureTimeProfileVector = " << assignment.g_DepartureTimeProfileVector.size() << '\n';
								g_program_stop();

							}


							int slot_agent_size = 5;

							if (VehicleSize > 10)
								slot_agent_size = VehicleSize / max(1, dep_time_element.ending_time_slot_no - dep_time_element.starting_time_slot_no);

							if (slot_agent_size < 0.1)
								slot_agent_size = 5;

							double headway_in_min = max(0.033, (double) MIN_PER_TIMESLOT / max(1, slot_agent_size));


//							previous_headway_residual
							for (int v = 0; v < VehicleSize; ++v)
							{

								//int slot_no = dep_time_element.get_time_slot_no(v, VehicleSize);
								//double random_first_headway_min_fraction = g_get_random_ratio();

								//double fraction = (float)(g_agent_simu_vector.size() % slot_agent_size) / slot_agent_size;

								//time_stamp = (slot_no - dep_time_element.starting_time_slot_no) * MIN_PER_TIMESLOT + fraction* MIN_PER_TIMESLOT;


								int headway_vehicle_size = 3;
								double random_headway_fraction = (float)(g_agent_simu_vector.size() % headway_vehicle_size) / headway_vehicle_size * 0.05;
								time_stamp = dep_time_element.get_deparure_time_in_min(v, VehicleSize) + random_headway_fraction;



								if (it->second.m_link_size == 0)   // only load agents with physical path
									continue;

								//if (g_agent_simu_vector.size() < 10)
								//{
								//	dtalog.output() << std::setprecision(5) << "generate vehicle: slot_no = " << slot_no << ",t = " << time_stamp << ":" << time_stamp / 60 << " min, agent id = " << g_agent_simu_vector.size() << " with path node size = "
								//		<< it->second.m_link_size << '\n';
								//}

								CAgent_Simu* p_agent = new CAgent_Simu();

								if (p_agent->agent_id == 45)
								{
									int idebug = 1;
								}
								// for future use of column pool
								p_agent->at = at;
								p_agent->dest = dest;
								p_agent->tau = tau;

								if (assignment.g_rt_network_pool != NULL)
								{
									p_agent->p_RTNetwork = assignment.g_rt_network_pool[dest][at][tau];
								}
								p_agent->PCE_unit_size = 1; //  max(1, (int)(assignment.g_ModeTypeVector[at].PCE + 0.5));  // convert a possible floating point pce to an integer value for simulation
								p_agent->desired_free_travel_time_ratio = max(1.0, 1.0 / max(0.01, assignment.g_ModeTypeVector[at].DSR));
								p_agent->time_headway = (int)(assignment.g_ModeTypeVector[at].time_headway_in_sec / number_of_seconds_per_interval + 0.5);
								p_agent->agent_id = g_agent_simu_vector.size();
								p_agent->mode_type_no = assignment.g_ModeTypeVector[at].mode_type_no;
								p_agent->departure_time_in_min = time_stamp;

								p_agent->path_travel_time_in_min = g_LoadingEndTimeInMin - p_agent->departure_time_in_min;  // by default
								it->second.agent_simu_id_vector.push_back(p_agent->agent_id);

								int simulation_time_intervalNo = (int)(p_agent->departure_time_in_min * 60 / number_of_seconds_per_interval);

								if (simulation_time_intervalNo < 0)  // reset to a zero index departure time interval
								{
									simulation_time_intervalNo = 0;
									int idebug = 1;
								}

								if (simulation_time_intervalNo > 1440 * 60 * 4)
								{
									int idebug = 1;
								}

//								 dtalog.output() << "generate vehicle: t = " << simulation_time_intervalNo << ":" << simulation_time_intervalNo / 4.0 / 60 << " min, agent id = " << p_agent->agent_id<< '\n';


								g_AgentTDListMap[simulation_time_intervalNo].m_AgentIDVector.push_back(p_agent->agent_id);

								p_agent->impacted_flag = 0;
								p_agent->impacted_link_seq_no = 99999;
								for (int nl = 0; nl < it->second.m_link_size; ++nl)  // arc a
								{
									int link_seq_no = it->second.path_link_vector[nl];

									if (g_link_vector[link_seq_no].capacity_reduction_map.size() > 0)
									{
										p_agent->impacted_flag = 1;

										p_agent->impacted_link_seq_no = min(nl, p_agent->impacted_link_seq_no);  // early impacted link seq no
									}

									p_agent->path_link_seq_no_vector.push_back(link_seq_no);
								}

								p_agent->AllocateMemory();

								int FirstLink = p_agent->path_link_seq_no_vector[0];

								p_agent->m_veh_link_arrival_time_in_simu_interval[0] = simulation_time_intervalNo;
								int fftt_simu_interval = g_link_vector[FirstLink].free_flow_travel_time_in_min * p_agent->desired_free_travel_time_ratio * number_of_simu_intervals_in_min;
								p_agent->m_veh_link_departure_time_in_simu_interval[0] = p_agent->m_veh_link_arrival_time_in_simu_interval[0] + fftt_simu_interval;

								g_agent_simu_vector.push_back(p_agent);
							}
						}
					}
				}
			}
		}
	}

	dtalog.output() << "[DATA INFO] number of simulation zones:" << zone_size << '\n';
	dtalog.output() << "[STATUS INFO] generating " << g_agent_simu_vector.size() / 1000 << " K agents" << '\n';

	//dtalog.output() << "[DATA INFO] CI = cumulative impact count, CR= cumulative rerouting count" << '\n';
	assignment.summary_file << ", # of simulated agents in agent.csv=," << g_agent_simu_vector.size() << "," <<'\n';
	dtalog.output() << std::left << std::setprecision(4)
		<< std::setw(20) << "[DATA INFO]" 
		<< std::setw(25) << "In-simulation Time (min)"
		<< std::setw(20) << "Cumu. Arrival"
		<< std::setw(20) << "Cumu. Departure"
		<< std::setw(15) << "Speed"
		<< std::setw(15) << "Travel Time" << '\n';

	int current_active_agent_id = 0;
	// the number of threads is redifined.
	int number_of_threads = omp_get_max_threads();

	int number_of_in_min_for_RTSP = 5;

	int link_size = g_link_vector.size();

	// initialize CA and CD count for each link

	bool cell_based_simulation_mode = false;
	for (int li = 0; li < link_size; ++li)
	{
		CLink* p_link = &(g_link_vector[li]);

		if (p_link->cell_type >= 0)  // activate cell based simulation mode
			cell_based_simulation_mode = true;

		m_link_CA_count[li] = 0;
		m_link_CD_count[li] = 0;

		g_link_vector[li].win_count = 0;
		g_link_vector[li].lose_count = 0;
		g_link_vector[li].time_to_be_released = -1;

	}

	//std::map<int, DTAVehListPerTimeInterval>::iterator it_agent;

	//for (it_agent = g_AgentTDListMap.begin(); it_agent != g_AgentTDListMap.end(); ++it_agent)
	//{
	//    dtalog.output() << "t = " << it_agent->first << ":" << it_agent->first/4.0/60 <<" min," << it_agent->second.m_AgentIDVector.size() << '\n';
	//}


	bool bRealTimeInformationActivated = false;
	// first loop for time t
	for (int t = 0; t < g_number_of_simulation_intervals; ++t)
	{
		float current_time_in_min = t * number_of_seconds_per_interval / 60.0 + g_LoadingStartTimeInMin;


		if (t % number_of_simu_intervals_in_min == 0)  // every min
		{
			int time_in_min = t / number_of_simu_intervals_in_min;
			for (int li = 0; li < link_size; ++li)
			{
				CLink* p_link = &(g_link_vector[li]);
				{
					m_LinkCumulativeArrivalVector[li][time_in_min] = m_link_CA_count[li];
					m_LinkCumulativeDepartureVector[li][time_in_min] = m_link_CD_count[li];

				}
			}
		}


		int number_of_simu_interval_per_min = 60 / number_of_seconds_per_interval;
		double network_wide_speed = 60;
		double network_wide_travel_time = 0;

		if (TotalVehicleHourTraveled > 0.001)
		{
			network_wide_speed = TotalVehicleMileTraveled / max(0.0001, TotalVehicleHourTraveled);
			network_wide_travel_time = TotalVehicleHourTraveled / max(1, TotalCumulative_Departure_Count) * 60.0; // from hour to min
		}

		int info_updating_freq_in_min_in_in_6sec = assignment.g_info_updating_freq_in_min * 10;
		int current_time_in_6sec = current_time_in_min * 10;  // *10: 0.1 min, 6 seconds.

		if (t % number_of_simu_intervals_in_min == 0 /*// every min*/ && current_time_in_min >= system_information_activate_time_in_min - 1 && current_time_in_6sec % info_updating_freq_in_min_in_in_6sec == 0) // RT updating every x=5 or even 0.1 min
		{ // updating real time link travel time for RT destination based shortest path
#pragma omp parallel for    // parallel computing for each link
			for (int i = 0; i < g_link_vector.size(); ++i)
			{
				CLink* p_link = &(g_link_vector[i]);
				p_link->RT_waiting_time = 0;
				if (p_link->link_type_si[0] >= 0)  // do not need to be the cap reduced link
				{

					double total_waiting_time_in_min = 0;

					for (auto it = p_link->ExitQueue.begin(); it != p_link->ExitQueue.end(); ++it)
					{
						int agent_id = (*it);
						CAgent_Simu* p_agent = g_agent_simu_vector[agent_id];

						int current_link_seq_no = p_agent->path_link_seq_no_vector[p_agent->m_current_link_seq_no];
						int arrival_time_in_t = p_agent->m_veh_link_arrival_time_in_simu_interval[p_agent->m_current_link_seq_no];
						int waiting_time_in_t = t - arrival_time_in_t;
						total_waiting_time_in_min += waiting_time_in_t * number_of_seconds_per_interval / 60; // in min

					}

					if (p_link->m_link_pedefined_capacity_map_in_sec.size() > 0)
					{
						int current_time_in_relative_sec = (current_time_in_min - assignment.g_LoadingStartTimeInMin) * 60;  // *10: 0.1 min, 6 seconds.

						if (p_link->m_link_pedefined_capacity_map_in_sec.find(current_time_in_relative_sec) != p_link->m_link_pedefined_capacity_map_in_sec.end())
						{
							p_link->RT_waiting_time = 99;
						}
					}
					else
					{
						p_link->RT_waiting_time = total_waiting_time_in_min / max((size_t)1, p_link->ExitQueue.size());  // average travel time
					}

					//int timestamp_in_min = g_LoadingStartTimeInMin + t / number_of_simu_intervals_in_min;
				//	p_link->RT_travel_time_map[timestamp_in_min] = p_link->ExitQueue.size() + p_link->RT_waiting_time;
				}
			}
			int recording_flag = 0;

			if (t / number_of_simu_interval_per_min == system_information_activate_time_in_min + 5)
				recording_flag = 0;

			g_rt_info_column_generation(this, current_time_in_min, recording_flag);  // parallel computing inside, and to provide rt information for vehicles to access

		}




		if (t % (number_of_simu_interval_per_min * 10) == 0)
		{
			// Log the information
			dtalog.output() << std::left << std::setprecision(4)
				<< std::setw(20) << "[DATA INFO]" 
				<< std::setw(25) << t / number_of_simu_interval_per_min
				<< std::setw(20) << TotalCumulative_Arrival_Count
				<< std::setw(20) << TotalCumulative_Departure_Count
				<< std::setw(15) << network_wide_speed
				<< std::setw(15) << network_wide_travel_time << '\n';


			assignment.summary_file << std::setprecision(5) << ",simu time=," << t / number_of_simu_interval_per_min << " min, CA =," << TotalCumulative_Arrival_Count << ",CD=," << TotalCumulative_Departure_Count
				<< ", speed=," << network_wide_speed << ", travel time =," << network_wide_travel_time <<
				",CI=," << TotalCumulative_impact_count <<
				",CR=," << TotalCumulative_rerouting_count
				<< '\n';
			 simu_log_file << "," << std::setprecision(5) << "simu time= " << t / number_of_simu_interval_per_min << " min, CA = " << TotalCumulative_Arrival_Count << " CD=" << TotalCumulative_Departure_Count
				<< ", speed=" << network_wide_speed << ", travel time = " << network_wide_travel_time <<
				",CI=" << TotalCumulative_impact_count <<
				",CR=" << TotalCumulative_rerouting_count
				<< '\n';
		}

		if (g_AgentTDListMap.find(t) != g_AgentTDListMap.end())
		{
			for (int Agent_v = 0; Agent_v < g_AgentTDListMap[t].m_AgentIDVector.size(); ++Agent_v)
			{
				int agent_id = g_AgentTDListMap[t].m_AgentIDVector[Agent_v];

				CAgent_Simu* p_agent = g_agent_simu_vector[agent_id];
				p_agent->m_bGenereated = true;
				int FirstLink = p_agent->path_link_seq_no_vector[0];
				m_link_CA_count[FirstLink] += 1;

#pragma omp critical  // counting impacted vehicles once t
				{
					if (p_agent->impacted_flag == 1)
					{
						TotalCumulative_impact_count++;
					}
				}

				g_link_vector[FirstLink].EntranceQueue.push_back(p_agent->agent_id);


				if (trace_agent_id == agent_id)
				{
					simu_log_file << "trace tag 1: simu time= " << t / number_of_simu_interval_per_min << " min, , traced vehicle enters the network = " << '\n';
				}


				TotalCumulative_Arrival_Count++;
			}
		}

#pragma omp parallel for    // parallel computing for each link
		for (int li = 0; li < link_size; ++li)
		{
			CLink* p_link = &(g_link_vector[li]);

			// if there are Agents in the entrance queue
			while (p_link->EntranceQueue.size() > 0)
			{
				int agent_id = p_link->EntranceQueue.front();

				p_link->EntranceQueue.pop_front();
				p_link->ExitQueue.push_back(agent_id);
				CAgent_Simu* p_agent = g_agent_simu_vector[agent_id];
				if (trace_agent_id == agent_id)
				{
					simu_log_file << "trace tag 2: simu time interval = " << t << " min, , traced vehicle moves from entrance queue to exit queue on link = "
						<< g_node_vector[p_link->from_node_seq_no].node_id << " -> " << g_node_vector[p_link->to_node_seq_no].node_id <<
						" on its link seq.no " << p_agent->m_current_link_seq_no << '\n';

					trace_link_no = li;
					trace_node_no = g_node_vector[p_link->to_node_seq_no].node_seq_no;
				}

				int arrival_time_in_simu_interval = p_agent->m_veh_link_arrival_time_in_simu_interval[p_agent->m_current_link_seq_no];
				int link_travel_time_in_simu_interavls = max(1.0, g_link_vector[li].free_flow_travel_time_in_min * number_of_simu_intervals_in_min);
				int fftt_simu_interval = (int)(g_link_vector[li].free_flow_travel_time_in_min * p_agent->desired_free_travel_time_ratio * number_of_simu_intervals_in_min + 0.5);
				p_agent->m_veh_link_departure_time_in_simu_interval[p_agent->m_current_link_seq_no] = arrival_time_in_simu_interval + fftt_simu_interval;

				//dtalog.output() << "reserve TD at time t = " << t << " on link" << g_node_vector[p_link->from_node_seq_no].node_id << " -> " << g_node_vector[p_link->to_node_seq_no].node_id <<
				//    " for agent = " << agent_id << " on link seq. no " << p_agent->m_current_link_seq_no << " with travel time in simu intervas = " << link_travel_time_in_simu_interavls << '\n';

			}
		}

		/// for each link, for all vehicles in the queue
		///                         // back trace the resourece use vector right after the vehicle moves to the next link
		for (int li = 0; li < link_size; ++li)
		{
			CLink* p_link = &(g_link_vector[li]);
			int queue_position_no = 0;
			for (auto it = p_link->ExitQueue.begin(); it != p_link->ExitQueue.end(); ++it)
			{
				int agent_id = (*it);
				CAgent_Simu* p_agent = g_agent_simu_vector[agent_id];

				if (trace_agent_id == agent_id)
				{
					simu_log_file << "trace tag 2.5: simu time= " << t / number_of_simu_interval_per_min << " min, , traced vehicle is waiting at the exit queue on link = "
						<< g_node_vector[p_link->from_node_seq_no].node_id << " -> " << g_node_vector[p_link->to_node_seq_no].node_id <<
						" on its link seq.no " << p_agent->m_current_link_seq_no <<
						" queue pos.no " << queue_position_no <<
						" of queue size = " << p_link->ExitQueue.size() << '\n';

					trace_link_no = li;
				}

				queue_position_no++;

				if (p_agent->PCE_unit_size >= 2)
				{
					for (int l_backtrace = 1; l_backtrace < p_agent->PCE_unit_size; l_backtrace++)
					{
						int local_l_index = p_agent->m_current_link_seq_no - l_backtrace;
						if (local_l_index >= 0)
						{
							int current_link_seq_no = p_agent->path_link_seq_no_vector[local_l_index];
							CLink* p_current_link = &(g_link_vector[current_link_seq_no]);
							p_current_link->time_to_be_released = t + p_agent->time_headway;
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
			if (node == trace_node_no)
			{
				simu_log_file << "trace tag 3F_1: simu time int = " << t << ", , check node incoming mode "
					<< g_node_vector[node].node_id << '\n';
			}

			bool node_resource_competing_mode = false;
			int incoming_link_request_count = 0;
			int incoming_link_index_FIFO = -1;
			int debug_node_resource_competing_mode = 0;
			if (cell_based_simulation_mode)
			{
				g_node_vector[node].next_link_for_resource_request.clear();


				CLink* p_link;
				CLink* pNextLink;
				for (int i = 0; i < g_node_vector[node].m_incoming_link_seq_no_vector.size(); ++i)
				{
					int link = g_node_vector[node].m_incoming_link_seq_no_vector[i];  // we will start with different first link from the incoming link list,
					// equal change, regardless of # of lanes or main line vs. ramp, but one can use service arc, to control the effective capacity rates, e.g. through a metered ramp, to
					// allow mainline to use the remaining flow
					p_link = &(g_link_vector[link]);

					if (node == trace_node_no)
					{
						simu_log_file << "trace tag 3F_2: simu time int = " << t << ", , check node incoming mode "
							<< g_node_vector[node].node_id << ",i= " << i << ", link id"
							<< g_node_vector[p_link->from_node_seq_no].node_id << " -> " << g_node_vector[p_link->to_node_seq_no].node_id

							<< '\n';
					}



					if (p_link->ExitQueue.size() >= 1)
					{
						int agent_id = p_link->ExitQueue.front();
						CAgent_Simu* p_agent = g_agent_simu_vector[agent_id];
						int next_link_seq_no = p_agent->path_link_seq_no_vector[p_agent->m_current_link_seq_no + 1];
						pNextLink = &(g_link_vector[next_link_seq_no]);

						if (trace_agent_id == agent_id)
						{
							simu_log_file << "trace tag 3: simu time= " << t / number_of_simu_interval_per_min << " min, , traced vehicle wants to transfer from link = "
								<< g_node_vector[p_link->from_node_seq_no].node_id << " -> " << g_node_vector[p_link->to_node_seq_no].node_id <<
								" on its link seq.no " << p_agent->m_current_link_seq_no << " to the next link. " << '\n';
						}

						int current_vehicle_count = m_link_CA_count[next_link_seq_no] - m_link_CD_count[next_link_seq_no];

						if (trace_agent_id == agent_id)
						{
							simu_log_file << "trace tag 3.5: simu time= " << t * 1.0 / number_of_simu_interval_per_min << " min, , traced vehicle checks next link = "
								<< g_node_vector[pNextLink->from_node_seq_no].node_id << " -> " << g_node_vector[pNextLink->to_node_seq_no].node_id <<
								", cell type = " << pNextLink->cell_type << ", current_vehicle_count = " << current_vehicle_count <<
								",pNextLink->spatial_capacity_in_vehicles= " << pNextLink->spatial_capacity_in_vehicles << '\n';

						}

						if (pNextLink->cell_type >= 0 && current_vehicle_count < pNextLink->spatial_capacity_in_vehicles)  // only apply for cell mode
						{
							g_node_vector[node].next_link_for_resource_request[next_link_seq_no] = 1;
							incoming_link_request_count++;
						}
					}
				}


				if (incoming_link_request_count >= 2)
				{
					if (g_node_vector[node].next_link_for_resource_request.size() >= 1)
					{
						//                if(g_node_vector[node].node_id == 2347 && t / 240 >= 18)
						node_resource_competing_mode = true;

					}
				}

				if (node == trace_node_no)
				{
					simu_log_file << "trace tag 3F_2: simu time= " << t * 1.0 / number_of_simu_interval_per_min << " min, , check node incoming mode "
						<< g_node_vector[node].node_id
						<< "next_link_for_resource_request size = " << g_node_vector[node].next_link_for_resource_request.size()
						<< "node_resource_competing_mode = " << node_resource_competing_mode
						<< '\n';
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
						p_link = &(g_link_vector[link]);

						if (trace_link_no == link)
						{
							simu_log_file << "trace tag 3.6E: simu time= " << t * 1.0 / number_of_simu_interval_per_min << " min, , traced vehicle uses FIFO rule on link "
								<< g_node_vector[p_link->from_node_seq_no].node_id << " -> " << g_node_vector[p_link->to_node_seq_no].node_id <<
								" queue size = " << p_link->ExitQueue.size() << '\n';
						}
						if (p_link->cell_type >= 0 && p_link->ExitQueue.size() >= 1)
						{
							int agent_id = p_link->ExitQueue.front();
							CAgent_Simu* p_agent = g_agent_simu_vector[agent_id];
							if (p_agent->m_veh_link_departure_time_in_simu_interval[p_agent->m_current_link_seq_no] < earlest_departure_time_interval)
							{

								if (trace_agent_id == agent_id)
								{
									simu_log_file << "trace tag 3.6A: simu time= " << t * 1.0 / number_of_simu_interval_per_min << " min, , traced vehicle uses FIFO rule on link "
										<< g_node_vector[p_link->from_node_seq_no].node_id << " -> " << g_node_vector[p_link->to_node_seq_no].node_id <<
										" expectd departure time = " << p_agent->m_veh_link_departure_time_in_simu_interval[p_agent->m_current_link_seq_no] <<
										'\n';
								}

								earlest_departure_time_interval = p_agent->m_veh_link_departure_time_in_simu_interval[p_agent->m_current_link_seq_no];
								incoming_link_index_FIFO = i;  // i is the link index in the vector of m_incoming_link_seq_no_vector of node
							}
							else
							{
								if (trace_agent_id == agent_id)
								{
									simu_log_file << "trace tag 3.6B: simu time= " << t * 1.0 / number_of_simu_interval_per_min << " min, , traced vehicle fails to use FIFO rule on link "
										<< g_node_vector[p_link->from_node_seq_no].node_id << " -> " << g_node_vector[p_link->to_node_seq_no].node_id <<
										" expectd departure time = " << p_agent->m_veh_link_departure_time_in_simu_interval[p_agent->m_current_link_seq_no] <<
										" while " <<
										" earlest_departure_time_interval = " << earlest_departure_time_interval << '\n';
								}

							}

						}
					}
				}


			}
			// for each incoming link

			if (node == trace_node_no)
			{
				simu_log_file << "trace tag 3F: simu time= " << t * 1.0 / number_of_simu_interval_per_min << " min, , check node incoming mode "
					<< g_node_vector[node].node_id
					<< "incoming_link_index_FIFO = " << incoming_link_index_FIFO
					<< '\n';
			}

			for (int i = 0; i < g_node_vector[node].m_incoming_link_seq_no_vector.size(); ++i)
			{
				int incoming_link_index = i;

				if (incoming_link_index_FIFO >= 0)
				{
					incoming_link_index = (i + incoming_link_index_FIFO) % (g_node_vector[node].m_incoming_link_seq_no_vector.size());  // cycle loop
					if (debug_node_resource_competing_mode)
						dtalog.output() << "[DATA INFO] FIFO judgement i = " << i << '\n';
				}
				else
				{  // random mode
					incoming_link_index = (i + t) % (g_node_vector[node].m_incoming_link_seq_no_vector.size());
				}


				int link = g_node_vector[node].m_incoming_link_seq_no_vector[incoming_link_index];  // we will start with different first link from the incoming link list,
				// equal change, regardless of # of lanes or main line vs. ramp, but one can use service arc, to control the effective capacity rates, e.g. through a metered ramp, to
				// allow mainline to use the remaining flow
				CLink* p_link = &(g_link_vector[link]);

				if (node == trace_node_no)
				{
					simu_log_file << "trace tag 3F: simu time= " << t * 1.0 / number_of_simu_interval_per_min << " min, , check node incoming mode "
						<< g_node_vector[node].node_id
						<< ", check link" << g_node_vector[p_link->from_node_seq_no].node_id << " -> " << g_node_vector[p_link->to_node_seq_no].node_id << '\n';

				}


				// check if the current link has sufficient capacity
				// most critical and time-consuming task, check link outflow capacity

				int time_in_sec = t * number_of_seconds_per_interval;


				while (m_LinkOutFlowCapacity[link][time_in_sec] >= 1 && p_link->ExitQueue.size() >= 1)
				{
					int agent_id = p_link->ExitQueue.front();
					CAgent_Simu* p_agent = g_agent_simu_vector[agent_id];

					if (p_agent->m_veh_link_departure_time_in_simu_interval[p_agent->m_current_link_seq_no] > t)
					{
						// the future departure time on this link is later than the current time
						break;
					}

					if (p_agent->m_current_link_seq_no == p_agent->path_link_seq_no_vector.size() - 2)  //-2 do not consider virtual destination arc
					{
						// end of path
						p_link->ExitQueue.pop_front();

						p_agent->m_bCompleteTrip = true;
						m_link_CD_count[link] += 1;
						p_agent->arrival_time_in_min = t * 1.0 / number_of_simu_interval_per_min;
						p_agent->path_travel_time_in_min = p_agent->arrival_time_in_min - p_agent->departure_time_in_min;

						if (p_agent->agent_id == 0)
						{
							int debug = 1;
						}
						p_agent->path_distance = 0;

						for (int link_s = 0; link_s < p_agent->path_link_seq_no_vector.size(); link_s++)
						{
							int link_seq_no = p_agent->path_link_seq_no_vector[link_s];
							//path_toll += g_link_vector[link_seq_no].VDF_period[p_agent->tau].toll[at][assignment.active_scenario_index];
							p_agent->path_distance += g_link_vector[link_seq_no].link_distance_VDF;
						}

#pragma omp critical
						{
							TotalCumulative_Departure_Count += 1;
							TotalVehicleHourTraveled += p_agent->path_travel_time_in_min / 60.0;
							TotalVehicleMileTraveled += p_agent->path_distance;
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

						if (next_link_seq_no < 0 && next_link_seq_no < g_link_vector.size())
						{
							int ierror = 1;
						}
						CLink* pNextLink = &(g_link_vector[next_link_seq_no]);

						if (p_agent->diverted_flag == 2)
						{
							int debug = 1;
							int arrival_time_in_t = p_agent->m_veh_link_arrival_time_in_simu_interval[p_agent->m_current_link_seq_no];
							int waiting_time_in_t = t - arrival_time_in_t;
							float total_waiting_time_in_min = waiting_time_in_t * number_of_seconds_per_interval / 60; // in min
							if (total_waiting_time_in_min > 2)
							{
								int idebug2 = 1;
							}

						}
						//// spatial queue
						if (pNextLink->link_distance_VDF <= 0.010)  // cell based link
						{
							pNextLink->spatial_capacity_in_vehicles = 1; // for cell based model
						}

						if (1/*pNextLink->cell_type >= 0 || pNextLink->traffic_flow_code == 2*/)  // DTA micro cell based simulation
						{
							int current_vehicle_count = m_link_CA_count[next_link_seq_no] - m_link_CD_count[next_link_seq_no];
							if (current_vehicle_count >= pNextLink->spatial_capacity_in_vehicles || (t < pNextLink->time_to_be_released))
								// this is an OR condition, related to reaction time tau
							{
								// spatical queue in the next link is blocked, break the while loop from here, as a first in first out queue.
								// g_fout << "spatical queue in the next link is blocked on link seq  " <<  g_node_vector[pNextLink->from_node_seq_no].node_id  << " -> " << g_node_vector[pNextLink->to_node_seq_no].node_id  <<'\n';
								//dtalog.output() << "spatical queue in the next link is blocked at time  = " << t << ",link" << g_node_vector[p_link->from_node_seq_no].node_id << " -> " << g_node_vector[p_link->to_node_seq_no].node_id <<
								//    ",current_vehicle_count = " << current_vehicle_count << '\n';
								if (node_resource_competing_mode)
								{
									p_link->lose_count++;
									if (debug_node_resource_competing_mode)
									{
										simu_log_file << "lose: the next link is blocked at time  = " << t * 1.0 << ", from link" << g_node_vector[p_link->from_node_seq_no].node_id << " -> " << g_node_vector[p_link->to_node_seq_no].node_id
											<< " to link" << g_node_vector[pNextLink->from_node_seq_no].node_id << " -> " << g_node_vector[pNextLink->to_node_seq_no].node_id
											<< " time request " << p_agent->m_veh_link_departure_time_in_simu_interval[p_agent->m_current_link_seq_no]
											<< '\n';
									}
								}
								break;  // being blocked
							}
							else // free to go if none of the above conditions cannot be met
							{
								if (node_resource_competing_mode)
								{
									p_link->win_count++;
									if (debug_node_resource_competing_mode)
									{
										dtalog.output() << "[DATA INFO] win:spatical queue in the next link is blocked at time  = " << t << ",link" << g_node_vector[p_link->from_node_seq_no].node_id << " -> " << g_node_vector[p_link->to_node_seq_no].node_id
											<< " to link" << g_node_vector[pNextLink->from_node_seq_no].node_id << " -> " << g_node_vector[pNextLink->to_node_seq_no].node_id
											<< " time request " << p_agent->m_veh_link_departure_time_in_simu_interval[p_agent->m_current_link_seq_no]
											<< '\n';
									}
								}


							}

						}

						                // spatal queue
						if (pNextLink->traffic_flow_code == spatial_queue)
						{

						    int current_vehicle_count = m_link_CA_count[next_link_seq_no] - m_link_CD_count[next_link_seq_no];
						    if (current_vehicle_count > pNextLink->spatial_capacity_in_vehicles)
						    {
						        // spatical queue in the next link is blocked, break the while loop from here, as a first in first out queue.
								// g_fout << "spatical queue in the next link is blocked on link seq  " <<  g_node_vector[pNextLink->from_node_seq_no].node_id  << " -> " << g_node_vector[pNextLink->to_node_seq_no].node_id  <<'\n';
						        break;
						    }
						}

						//                // kinematic wave
						//                if (pNextLink->traffic_flow_code == 3)  // to be discussed later with Cafer
						//                {
						//                    int lagged_time_stamp = max(0, t - 1 - pNextLink->BWTT_in_simulation_interval);
						//                    int current_vehicle_count = m_link_CA_count[next_link_seq_no] - m_link_CD_count[next_link_seq_no];
						//                    if (current_vehicle_count > pNextLink->spatial_capacity_in_vehicles)
						//                    {
						//                        // spatical queue in the next link is blocked, break the while loop from here, as a first in first out queue.
												//// g_fout << "spatical queue in the next link is blocked on link seq  " <<  g_node_vector[pNextLink->from_node_seq_no].node_id  << " -> " << g_node_vector[pNextLink->to_node_seq_no].node_id  <<'\n';
						//                        break;
						//                    }
						//                }

						p_link->ExitQueue.pop_front();

						//if(p_link->ExitQueue.size() !=0 || p_link->EntranceQueue.size() != 0)
						//{
						//    simu_log_file << "Error at time = " << t << " , link = "
						//        << g_node_vector[p_link->from_node_seq_no].node_id << " -> " << g_node_vector[p_link->to_node_seq_no].node_id <<

						//        << "ExitQueue.size() = " << p_link->ExitQueue.size()
						//        << "EntranceQueue.size() = " << p_link->EntranceQueue.size()
						//        << "spatial queue capacity =  " << p_link->spatial_capacity_in_vehicles << '\n';
						//}

						p_link->time_to_be_released = t + p_agent->time_headway;

						// back trace the resourece use vector right after the vehicle moves to the next link
						if (p_agent->PCE_unit_size >= 2)
						{
							for (int l_backtrace = 1; l_backtrace < p_agent->PCE_unit_size; l_backtrace++)
							{
								int local_l_index = p_agent->m_current_link_seq_no - l_backtrace;
								if (local_l_index >= 0)
								{
									int current_link_seq_no = p_agent->path_link_seq_no_vector[local_l_index];
									CLink* p_current_link = &(g_link_vector[current_link_seq_no]);
									p_current_link->time_to_be_released = t + p_agent->time_headway;
									// big truck/bus,
									//backtrace the previous l - k links, k = 1, 2, K, and set release time
									  //  K as the number of units of truck
								}

							}

						}

						if (trace_agent_id == agent_id)
						{
							simu_log_file << "trace tag 4: simu time interval = " << t << " , traced vehicle transfers from link = "
								<< g_node_vector[p_link->from_node_seq_no].node_id << " -> " << g_node_vector[p_link->to_node_seq_no].node_id <<
								"  to next lnk " <<
								g_node_vector[pNextLink->from_node_seq_no].node_id << " -> " << g_node_vector[pNextLink->to_node_seq_no].node_id << '\n';
						}

						pNextLink->EntranceQueue.push_back(agent_id);


						p_agent->m_veh_link_departure_time_in_simu_interval[p_agent->m_current_link_seq_no] = t;

						p_agent->m_veh_link_arrival_time_in_simu_interval[p_agent->m_current_link_seq_no + 1] = t;

						int waiting_time_in_simu_interval = p_agent->m_veh_link_departure_time_in_simu_interval[p_agent->m_current_link_seq_no] - p_agent->m_veh_link_arrival_time_in_simu_interval[p_agent->m_current_link_seq_no];
						float travel_time_in_sec = waiting_time_in_simu_interval * number_of_seconds_per_interval;
						//for each waited vehicle
						float waiting_time_in_min = travel_time_in_sec / 60.0 - p_link->free_flow_travel_time_in_min;

						m_link_TD_waiting_time[link][p_agent->m_veh_link_arrival_time_in_simu_interval[p_agent->m_current_link_seq_no] / number_of_simu_intervals_in_min] += waiting_time_in_min;
						m_link_total_waiting_time_vector[link] += waiting_time_in_min;
						p_agent->waiting_time_in_min += waiting_time_in_min;

						p_link->total_simulated_delay_in_min += waiting_time_in_min;

						int meso_link_id = p_link->meso_link_id;

						if (p_agent->meso_link_id_map.find(meso_link_id) == p_agent->meso_link_id_map.end())  // detect if this is first time visit
						{
							p_link->total_simulated_meso_link_incoming_volume += 1;
							p_agent->meso_link_id_map[meso_link_id] = 1;  // mark as visited
							p_agent->path_meso_link_id_vector.push_back(meso_link_id);
						}
						else
						{
							// do nothing
						}

						if (waiting_time_in_min < 0)
						{
							waiting_time_in_min = 0;
						}

						if (p_agent->m_current_link_seq_no >=1 && waiting_time_in_min >=0.001 && waiting_time_in_min > p_agent->max_link_waiting_time_in_min)
						{
							p_agent->max_waiting_time_link_no = link;
							p_agent->max_link_waiting_time_in_min = waiting_time_in_min;

						}
						m_link_CD_count[link] += 1;
						m_link_CA_count[next_link_seq_no] += 1;

						if (p_agent->agent_id == trace_agent_id)
						{
							simu_log_file << "Simu time interval = " << t << " , traced vehicle transfers to node = " << g_node_vector[p_link->to_node_seq_no].node_id <<
								" with travel time = " << travel_time_in_sec << " sec" << '\n';
						}

						//if (t == system_information_activate_time_in_simu_interval)
						//{
						//	int idebug = 2;
						//}
						//test rerouting
						//if (t >= system_information_activate_time_in_simu_interval)
						bool rerouting_decision_flag = false;
						int visual_distance_in_number_of_cells = max(1, g_visual_distance_in_cells);  // 5 cells ahead


						if (p_link->cell_type < 0)  // link based mode
							visual_distance_in_number_of_cells = 1;

						// case 1: visual distance based information
						if (t >= system_information_activate_time_in_simu_interval && p_agent->impacted_flag == 1 && (p_agent->info_receiving_flag != 2 && p_agent->info_receiving_flag != 3)&& p_agent->impacted_link_seq_no == p_agent->m_current_link_seq_no + visual_distance_in_number_of_cells)
						{
							rerouting_decision_flag = true;
							p_agent->info_receiving_flag = 1;
						}

						// case 2: dynamic variable message signs
						int current_time_in_min_int = current_time_in_min;
						if (t >= system_information_activate_time_in_simu_interval &&
							p_link->m_link_pedefined_information_response_map.find(current_time_in_min_int) != p_link->m_link_pedefined_information_response_map.end())
						{

							if (p_agent->agent_id % 100 <= p_link->m_link_pedefined_information_response_map[current_time_in_min_int] * 100)  // modeling response uniformly
							{
								rerouting_decision_flag = true;
								p_agent->info_receiving_flag = 2;
							}
						}

						// case 3: g_real_time_info_ratio
						if (t >= system_information_activate_time_in_simu_interval && p_agent->impacted_flag == 1 && p_agent->m_current_link_seq_no >= 1 && p_agent->info_receiving_flag == 0)
						{  // touch the physical road, and receive no information so far : info_receiving_flag =0;

							int mode_type_no = p_agent->mode_type_no;
							//requre mode_type_no be sequential
							if (mode_type_no-1 < g_ModeTypeVector.size() && g_ModeTypeVector[p_agent->mode_type_no - 1].real_time_information_type)  // modeling response uniformly
							{
								rerouting_decision_flag = true;
								p_agent->info_receiving_flag = 3;
							}
						}

						if (rerouting_decision_flag) // with DMS or visual distance info, to obtain destination based tree to use latest real time information
						{
							//simu_log_file << " start rerouting  at seq no" << p_agent->m_current_link_seq_no << " for agent" << p_agent->agent_id << '\n';

							int impacted_flag_change = 0;
							float timestamp_in_min = g_LoadingStartTimeInMin + t / number_of_simu_intervals_in_min;
							if (update_real_time_info_path(p_agent, impacted_flag_change, timestamp_in_min) == 1)
							{
#pragma omp critical
								{
									TotalCumulative_rerouting_count++;
								}
							}
						}
					}

					//move
					p_agent->m_current_link_seq_no += 1;
					m_LinkOutFlowCapacity[link][time_in_sec] -= 1;//deduct the capacity
				}
			}
		} // conditions
	}  // departure time events
	g_output_agent_csv(assignment);

}

// FILE* g_pFileOutputLog = nullptr;