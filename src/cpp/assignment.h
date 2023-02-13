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
using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::map;
using std::ifstream;
using std::ofstream;
using std::istringstream;


void g_reset_and_update_link_volume_based_on_columns(int number_of_links, int iteration_index, bool b_self_reducing_path_volume, bool b_sensitivity_analysis_flag)
{

	// reset the link volume
	for (int i = 0; i < number_of_links; ++i)
	{
		for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
		{
			// used in travel time calculation
			g_link_vector[i].PCE_volume_per_period[tau] = 0;
			g_link_vector[i].person_volume_per_period[tau] = 0;

			for (int at = 0; at < assignment.g_AgentTypeVector.size(); ++at)
				g_link_vector[i].person_volume_per_period_per_at[tau][at] = 0;
		}

		//for (int at = 0; at < assignment.g_AgentTypeVector.size(); ++at)
		//	for (int og = 0; og < assignment.g_number_of_analysis_districts; ++og)
		//	{
		//		g_link_vector[i].person_volume_per_district_per_at[og][at] = 0;
		//	}

	}

	if (iteration_index >= 0)
	{
		for (int at = 0; at < assignment.g_AgentTypeVector.size(); ++at)  //m
		{


			std::map<int, CColumnPath>::iterator it;
			int zone_size = g_zone_vector.size();
			int tau_size = assignment.g_DemandPeriodVector.size();

			float link_volume_contributed_by_path_volume;

			int link_seq_no;
			double PCE_ratio = 1;
			double OCC_ratio = 1;

			int nl;

			std::map<int, CColumnPath>::iterator it_begin;
			std::map<int, CColumnPath>::iterator it_end;

			int column_vector_size;
			CColumnVector* p_column_pool;

			for (int orig = 0; orig < zone_size; ++orig)  // o
			{   
				int from_zone_sindex = g_zone_vector[orig].sindex;
				if (from_zone_sindex == -1)
					continue;

				int analysis_district_id = assignment.g_zone_seq_no_to_analysis_distrct_id_mapping[orig];
				
				for (int dest = 0; dest < zone_size; ++dest) //d
				{
					int to_zone_sindex = g_zone_vector[dest].sindex;
					if (to_zone_sindex == -1)
						continue;

					for (int tau = 0; tau < tau_size; ++tau)  //tau
					{
						p_column_pool = &(assignment.g_column_pool[from_zone_sindex][to_zone_sindex][at][tau]);
						if (p_column_pool->od_volume > 0)
						{

							column_vector_size = p_column_pool->path_node_sequence_map.size();

							it_begin = p_column_pool->path_node_sequence_map.begin();
							it_end = p_column_pool->path_node_sequence_map.end();
							for (it = it_begin; it != it_end; ++it)
							{
								link_volume_contributed_by_path_volume = it->second.path_volume;  // assign all OD flow to this first path

								// add path volume to link volume
								for (nl = 0; nl < it->second.m_link_size; ++nl)  // arc a
								{
									link_seq_no = it->second.path_link_vector[nl];

									// MSA updating for the existing column pools
									// if iteration_index = 0; then update no flow discount is used (for the column pool case)
									PCE_ratio = g_link_vector[link_seq_no].VDF_period[tau].pce[at];  // updated on 08/16/2021 for link dependent and agent type dependent pce factor mainly for trucks 
									OCC_ratio = g_link_vector[link_seq_no].VDF_period[tau].occ[at];  // updated on 08/16/2021 for link dependent and agent type dependent pce factor mainly for trucks 
#pragma omp critical
									{
										g_link_vector[link_seq_no].PCE_volume_per_period[tau] += link_volume_contributed_by_path_volume * PCE_ratio;
										g_link_vector[link_seq_no].person_volume_per_period[tau] += link_volume_contributed_by_path_volume * OCC_ratio;
										g_link_vector[link_seq_no].person_volume_per_period_per_at[tau][at] += link_volume_contributed_by_path_volume * OCC_ratio;  // pure volume, not consider PCE
										//g_link_vector[link_seq_no].person_volume_per_district_per_at[analysis_district_id][at] += link_volume_contributed_by_path_volume * OCC_ratio;  // pure volume, not consider PCE

										

									}
								}

								if (b_sensitivity_analysis_flag == true) 
								{
									if (assignment.g_AgentTypeVector[at].real_time_information == 0)
									{
										// important condition here: working with Frank Zhang, Macro Li. we do not reset the path volume to zero for non-info user class /agent type
										if (orig == 2 && dest == 1)
										{
											int idebug = 1;
										}
										continue;
									}
									if (assignment.g_AgentTypeVector[at].real_time_information == 2 /*DMS*/)
									{

										if (assignment.zone_seq_no_2_info_mapping.find(g_zone_vector[orig].zone_seq_no) == assignment.zone_seq_no_2_info_mapping.end())
										{
											// this is DMS user with physical origin zone (other than information zone)

											continue;
										}
									}
									//only real time information == 0, or 2 as DMS with virtual information zone as the origin 
									if (!p_column_pool->bfixed_route && b_self_reducing_path_volume)
									{
										//after link volumn "tally", self-deducting the path volume by 1/(k+1) (i.e. keep k/(k+1) ratio of previous flow) so that the following shortes path will be receiving 1/(k+1) flow
										it->second.path_volume = it->second.path_volume * (double(iteration_index) / double(iteration_index + 1));
									}

								}
								else
								{  // first stage before sensitivity analysis
									// this  self-deducting action does not agents with fixed routing policies.
									if (!p_column_pool->bfixed_route && b_self_reducing_path_volume)
									{
										//after link volumn "tally", self-deducting the path volume by 1/(k+1) (i.e. keep k/(k+1) ratio of previous flow) so that the following shortes path will be receiving 1/(k+1) flow
										it->second.path_volume = it->second.path_volume * (double(iteration_index) / double(iteration_index + 1));
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

double update_link_travel_time_and_cost(int inner_iteration_number, double &total_distance)
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
		//	waiting_time_count += assignment.m_link_TD_waiting_time[l][tt/number_of_simu_intervals_in_min];   // tally total waiting cou
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
		g_link_vector[i].calculate_dynamic_VDFunction(inner_iteration_number, false, g_link_vector[i].vdf_type);

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
	total_distance = 0;
	for (int i = 0; i < g_link_vector.size(); ++i)
	{
		for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); ++tau)
		{
			total_network_travel_time += g_link_vector[i].VDF_period[tau].avg_travel_time * g_link_vector[i].VDF_period[tau].link_volume;

			total_distance += g_link_vector[i].length_in_meter * g_link_vector[i].VDF_period[tau].link_volume;
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

	double total_system_travel_cost = 0;
	double total_system_travel_time = 0;
	double total_system_demand = 0;
	double total_system_UE_gap = 0;
	// reset the link volume
	for (int i = 0; i < number_of_links; ++i)
	{
		for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
		{
			// used in travel time calculation
			g_link_vector[i].PCE_volume_per_period[tau] = 0;
			g_link_vector[i].person_volume_per_period[tau] = 0;

			for (int at = 0; at < assignment.g_AgentTypeVector.size(); ++at)
				g_link_vector[i].person_volume_per_period_per_at[tau][at] = 0;
		}

		//for (int at = 0; at < assignment.g_AgentTypeVector.size(); ++at)
		//	for (int og = 0; og < assignment.g_number_of_analysis_districts; ++og)
		//	{
		//		g_link_vector[i].person_volume_per_district_per_at[og][at] = 0;
		//	}

	}

	// reset the estimated production and attraction
	for (int orig = 0; orig < g_zone_vector.size(); ++orig)  // o
	{
		g_zone_vector[orig].est_attraction = 0;
		g_zone_vector[orig].est_production = 0;
	}

	for (int at = 0; at < assignment.g_AgentTypeVector.size(); ++at)  //m
	{

		int zone_size = g_zone_vector.size();
		int tau_size = assignment.g_DemandPeriodVector.size();
		float PCE_ratio = assignment.g_AgentTypeVector[at].PCE;
		float OCC_ratio = assignment.g_AgentTypeVector[at].OCC;


#pragma omp parallel for
		for (int orig = 0; orig < zone_size; ++orig)  // o
		{
			int from_zone_sindex = g_zone_vector[orig].sindex;
			if (from_zone_sindex == -1)
				continue;


			int analysis_district_id = assignment.g_zone_seq_no_to_analysis_distrct_id_mapping[orig];
			std::map<int, CColumnPath>::iterator it;
			float link_volume_contributed_by_path_volume;
			int nl;

			std::map<int, CColumnPath>::iterator it_begin;
			std::map<int, CColumnPath>::iterator it_end;

			int column_vector_size;
			CColumnVector* p_column_pool;
			for (int dest = 0; dest < zone_size; ++dest) //d
			{
				int to_zone_sindex = g_zone_vector[dest].sindex;

				if (to_zone_sindex == -1)
					continue;

				for (int tau = 0; tau < tau_size; ++tau)  //tau
				{
					p_column_pool = &(assignment.g_column_pool[from_zone_sindex][to_zone_sindex][at][tau]);
					if (p_column_pool->od_volume > 0)
					{


						// continuous: type 0
						column_vector_size = p_column_pool->path_node_sequence_map.size();

						it_begin = p_column_pool->path_node_sequence_map.begin();
						it_end = p_column_pool->path_node_sequence_map.end();

						double least_cost = 999999;
						int least_cost_path_seq_no = -1;
						int least_cost_path_node_sum_index = -1;
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

						least_cost = 999999;
						path_seq_count = 0;

						it_begin = p_column_pool->path_node_sequence_map.begin();
						it_end = p_column_pool->path_node_sequence_map.end();
						for (it = it_begin; it != it_end; ++it)
						{
							total_system_demand += it->second.path_volume;
							path_toll = 0;
							path_gradient_cost = 0;
							path_distance = 0;
							path_travel_time = 0;

							for (int nl = 0; nl < it->second.m_link_size; ++nl)  // arc a
							{
								link_seq_no = it->second.path_link_vector[nl];
								link_travel_time = g_link_vector[link_seq_no].travel_time_per_period[tau];
								path_travel_time += link_travel_time;

							}

							it->second.path_toll = path_toll;
							it->second.path_travel_time = path_travel_time;
							total_system_travel_time += (it->second.path_travel_time * it->second.path_volume);

							if (column_vector_size == 1)  // only one path
							{

								break;
							}

							if (path_travel_time < least_cost)
							{
								least_cost = path_travel_time;
								least_cost_path_seq_no = it->second.path_seq_no;
								least_cost_path_node_sum_index = it->first;
							}
#pragma omp critical
							{
								total_system_travel_cost += (it->second.path_travel_time * it->second.path_volume);
							}
						} // end for each path



						if (column_vector_size >= 2)
						{
							// step 2: calculate gradient cost difference for each column path
							total_switched_out_path_volume = 0;
							for (it = it_begin; it != it_end; ++it)
							{
								if (it->second.path_seq_no != least_cost_path_seq_no)  //for non-least cost path
								{
									it->second.UE_gap = it->second.path_travel_time - least_cost;
									it->second.UE_relative_gap = (it->second.path_travel_time - least_cost) / max(0.0001, least_cost);

#pragma omp critical
									{
										total_system_UE_gap += (it->second.UE_gap * it->second.path_volume);
									}




								}
							}

						}  // end for each path

						for (it = it_begin; it != it_end; ++it)  // path k
						{
							link_volume_contributed_by_path_volume = it->second.path_volume;  // assign all OD flow to this first path

#pragma omp critical
							{
								g_zone_vector[orig].est_production += it->second.path_volume;
								g_zone_vector[dest].est_attraction += it->second.path_volume;
							}
							// add path volume to link volume
							for (nl = 0; nl < it->second.m_link_size; ++nl)  // arc a
							{
								link_seq_no = it->second.path_link_vector[nl];

								// MSA updating for the existing column pools
								// if iteration_index = 0; then update no flow discount is used (for the column pool case)
#pragma omp critical
								{
									g_link_vector[link_seq_no].PCE_volume_per_period[tau] += link_volume_contributed_by_path_volume * PCE_ratio;
									g_link_vector[link_seq_no].person_volume_per_period[tau] += link_volume_contributed_by_path_volume * OCC_ratio;
									g_link_vector[link_seq_no].person_volume_per_period_per_at[tau][at] += link_volume_contributed_by_path_volume;  // pure volume, not consider PCE

									//g_link_vector[link_seq_no].person_volume_per_district_per_at[analysis_district_id][at] += link_volume_contributed_by_path_volume;  // pure volume, not consider PCE
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
		g_link_vector[i].calculate_dynamic_VDFunction(iteration_no, false, g_link_vector[i].vdf_type);
		for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); ++tau)  //tau
		{
			if (assignment.g_DemandPeriodVector[tau].number_of_demand_files == 0)
				continue;

			if (g_link_vector[i].VDF_period[tau].obs_count >= 1)  // with data
			{
				g_link_vector[i].VDF_period[tau].est_count_dev = g_link_vector[i].PCE_volume_per_period[tau] + g_link_vector[i].VDF_period[tau].preload - g_link_vector[i].VDF_period[tau].obs_count;

				if (dtalog.debug_level() == 2)
				{
					dtalog.output() << "link " << g_node_vector[g_link_vector[i].from_node_seq_no].node_id
						<< "->" << g_node_vector[g_link_vector[i].to_node_seq_no].node_id
						<< "obs:, " << g_link_vector[i].VDF_period[tau].obs_count << "est:, " << g_link_vector[i].PCE_volume_per_period[tau]
						<< "dev:," << g_link_vector[i].VDF_period[tau].est_count_dev << endl;
				}
				if (g_link_vector[i].VDF_period[tau].upper_bound_flag == 0)
				{
					total_gap += abs(g_link_vector[i].VDF_period[tau].est_count_dev);
					sub_total_gap_link_count += fabs(g_link_vector[i].VDF_period[tau].est_count_dev / g_link_vector[i].VDF_period[tau].obs_count);
					sub_total_system_gap_count += g_link_vector[i].VDF_period[tau].est_count_dev / g_link_vector[i].VDF_period[tau].obs_count;
				}
				else
				{  // upper bound constraints 
					if (g_link_vector[i].VDF_period[tau].est_count_dev > 0)
					{
						total_gap += abs(g_link_vector[i].VDF_period[tau].est_count_dev);
						sub_total_gap_link_count += fabs(g_link_vector[i].VDF_period[tau].est_count_dev / g_link_vector[i].VDF_period[tau].obs_count);
						sub_total_system_gap_count += g_link_vector[i].VDF_period[tau].est_count_dev / g_link_vector[i].VDF_period[tau].obs_count;
					}
				}
				total_link_count += 1;
			}
		}

	}
	//for (int orig = 0; orig < g_zone_vector.size(); ++orig)  // o
	//{
	//    if (g_zone_vector[orig].obs_attraction >= 1)  // with observation
	//    {
	//        g_zone_vector[orig].est_attraction_dev = g_zone_vector[orig].est_attraction - g_zone_vector[orig].obs_attraction;

	//        if (dtalog.debug_level() == 2)
	//        {
	//            dtalog.output() << "zone " << g_zone_vector[orig].zone_id << "A: obs:" << g_zone_vector[orig].obs_attraction
	//                << ",est:," << g_zone_vector[orig].est_attraction << ",dev:," << g_zone_vector[orig].est_attraction_dev << endl;
	//        }

	//        total_gap += abs(g_zone_vector[orig].est_attraction_dev);
	//        sub_total_gap_A_count += g_zone_vector[orig].est_attraction_dev / g_zone_vector[orig].obs_attraction;
	//    }

	//    if (g_zone_vector[orig].obs_production >= 1)  // with observation
	//    {
	//        g_zone_vector[orig].est_production_dev = g_zone_vector[orig].est_production - g_zone_vector[orig].obs_production;

	//        if (dtalog.debug_level() == 2)
	//        {
	//            dtalog.output() << "zone " << g_zone_vector[orig].zone_id << "P: obs:" << g_zone_vector[orig].obs_production
	//                << ",est:," << g_zone_vector[orig].est_production << ",dev:," << g_zone_vector[orig].est_production_dev << endl;
	//        }

	//        total_gap += abs(g_zone_vector[orig].est_production_dev);
	//        sub_total_gap_P_count += g_zone_vector[orig].est_production_dev / g_zone_vector[orig].obs_production;
	//    }
	//}

	assignment.summary_file << "ODME #" << iteration_no
		<< ", link MAE= " << total_gap / max(1, total_link_count)
		<< ",link_MAPE: " << (sub_total_gap_link_count) / max(1, total_link_count) * 100 <<
		"%,system_MPE: " << (sub_total_system_gap_count) / max(1, total_link_count) * 100 <<
		"%,avg_tt = " << total_system_travel_time / max(0.1, total_system_demand) << "(min) " <<
		",UE gap =" << total_system_UE_gap / max(0.00001, total_system_demand) << "(min)" <<
		" = (" << total_system_UE_gap / max(0.00001, total_system_travel_time) * 100 << " %)"
		<< endl;

	dtalog.output() << "ODME #" << iteration_no
		<< ", link MAE= " << total_gap / max(1, total_link_count)
		<< ",link_MAPE: " << (sub_total_gap_link_count) / max(1, total_link_count) * 100 <<
		"%,system_MPE: " << (sub_total_system_gap_count) / max(1, total_link_count) * 100 <<
		"%,avg_tt = " << total_system_travel_time / max(0.1, total_system_demand) << "(min) " <<
		",UE gap =" << total_system_UE_gap / max(0.00001, total_system_demand) << "(min)" <<
		" = (" << total_system_UE_gap / max(0.00001, total_system_travel_time) * 100 << " %)"
		<< endl;
	double gap = sub_total_gap_link_count / max(1, total_link_count);
	system_gap = sub_total_system_gap_count / max(1, total_link_count);

	return gap;
}


double g_reset_and_update_sensor_link_volume_based_on_ODME_columns(int number_of_links, int iteration_no, double& system_gap)
{
	float total_gap = 0;
	float sub_total_gap_link_count = 0;
	float sub_total_system_gap_count = 0;
	system_gap = 0;
	float sub_total_gap_P_count = 0;
	float sub_total_gap_A_count = 0;

	double total_system_travel_cost = 0;
	double total_system_travel_time = 0;
	double total_system_demand = 0;
	double total_system_UE_gap = 0;
	// reset the link volume
	for (int i = 0; i < number_of_links; ++i)
	{
			for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
			{

				if (g_link_vector[i].VDF_period[tau].obs_count >= 1)
				{
					// used in travel time calculation
				g_link_vector[i].PCE_volume_per_period[tau] = 0;
				g_link_vector[i].person_volume_per_period[tau] = 0;

				for (int at = 0; at < assignment.g_AgentTypeVector.size(); ++at)
					g_link_vector[i].person_volume_per_period_per_at[tau][at] = 0;
				}


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

		int zone_size = g_zone_vector.size();
		int tau_size = assignment.g_DemandPeriodVector.size();
		float PCE_ratio = assignment.g_AgentTypeVector[at].PCE;
		float OCC_ratio = assignment.g_AgentTypeVector[at].OCC;


#pragma omp parallel for
		for (int orig = 0; orig < zone_size; ++orig)  // o
		{
			int from_zone_sindex = g_zone_vector[orig].sindex;
			if (from_zone_sindex == -1)
				continue;


			int analysis_district_id = assignment.g_zone_seq_no_to_analysis_distrct_id_mapping[orig];
			std::map<int, CColumnPath>::iterator it;
			float link_volume_contributed_by_path_volume;
			int nl;

			std::map<int, CColumnPath>::iterator it_begin;
			std::map<int, CColumnPath>::iterator it_end;

			int column_vector_size;
			CColumnVector* p_column_pool;
			for (int dest = 0; dest < zone_size; ++dest) //d
			{
				int to_zone_sindex = g_zone_vector[dest].sindex;

				if (to_zone_sindex == -1)
					continue;

				for (int tau = 0; tau < tau_size; ++tau)  //tau
				{
					p_column_pool = &(assignment.g_column_pool[from_zone_sindex][to_zone_sindex][at][tau]);
					if (p_column_pool->od_volume > 0)
					{

						// continuous: type 0
						column_vector_size = p_column_pool->path_node_sequence_map.size();

						it_begin = p_column_pool->path_node_sequence_map.begin();
						it_end = p_column_pool->path_node_sequence_map.end();

						double least_cost = 999999;
						int least_cost_path_seq_no = -1;
						int least_cost_path_node_sum_index = -1;
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

						least_cost = 999999;
						path_seq_count = 0;

						it_begin = p_column_pool->path_node_sequence_map.begin();
						it_end = p_column_pool->path_node_sequence_map.end();


						for (it = it_begin; it != it_end; ++it)  // path k
						{
							link_volume_contributed_by_path_volume = it->second.path_volume;  // assign all OD flow to this first path

#pragma omp critical
							// add path volume to link volume
							for (nl = 0; nl < it->second.path_sensor_link_vector.size(); ++nl)  // arc a
							{
								link_seq_no = it->second.path_sensor_link_vector[nl];

								// MSA updating for the existing column pools
								// if iteration_index = 0; then update no flow discount is used (for the column pool case)
								{
									g_link_vector[link_seq_no].PCE_volume_per_period[tau] += link_volume_contributed_by_path_volume * PCE_ratio;
									g_link_vector[link_seq_no].person_volume_per_period[tau] += link_volume_contributed_by_path_volume * OCC_ratio;
									g_link_vector[link_seq_no].person_volume_per_period_per_at[tau][at] += link_volume_contributed_by_path_volume;  // pure volume, not consider PCE

									//g_link_vector[link_seq_no].person_volume_per_district_per_at[analysis_district_id][at] += link_volume_contributed_by_path_volume;  // pure volume, not consider PCE
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
		g_link_vector[i].calculate_dynamic_VDFunction(iteration_no, false, g_link_vector[i].vdf_type);
		for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); ++tau)  //tau
		{
			if (assignment.g_DemandPeriodVector[tau].number_of_demand_files == 0)
				continue;

			if (g_link_vector[i].VDF_period[tau].obs_count >= 1)  // with data
			{
				g_link_vector[i].VDF_period[tau].est_count_dev = g_link_vector[i].PCE_volume_per_period[tau] + g_link_vector[i].VDF_period[tau].preload - g_link_vector[i].VDF_period[tau].obs_count;

				if (dtalog.debug_level() == 2)
				{
					dtalog.output() << "link " << g_node_vector[g_link_vector[i].from_node_seq_no].node_id
						<< "->" << g_node_vector[g_link_vector[i].to_node_seq_no].node_id
						<< "obs:, " << g_link_vector[i].VDF_period[tau].obs_count << "est:, " << g_link_vector[i].PCE_volume_per_period[tau]
						<< "dev:," << g_link_vector[i].VDF_period[tau].est_count_dev << endl;
				}
				if (g_link_vector[i].VDF_period[tau].upper_bound_flag == 0)
				{
					total_gap += abs(g_link_vector[i].VDF_period[tau].est_count_dev);
					sub_total_gap_link_count += fabs(g_link_vector[i].VDF_period[tau].est_count_dev / g_link_vector[i].VDF_period[tau].obs_count);
					sub_total_system_gap_count += g_link_vector[i].VDF_period[tau].est_count_dev / g_link_vector[i].VDF_period[tau].obs_count;
				}
				else
				{  // upper bound constraints 
					if (g_link_vector[i].VDF_period[tau].est_count_dev > 0)
					{
						total_gap += abs(g_link_vector[i].VDF_period[tau].est_count_dev);
						sub_total_gap_link_count += fabs(g_link_vector[i].VDF_period[tau].est_count_dev / g_link_vector[i].VDF_period[tau].obs_count);
						sub_total_system_gap_count += g_link_vector[i].VDF_period[tau].est_count_dev / g_link_vector[i].VDF_period[tau].obs_count;
					}
				}
				total_link_count += 1;
			}
		}

	}
	assignment.summary_file << "ODME #" << iteration_no
		<< ", link MAE= " << total_gap / max(1, total_link_count)
		<< ",link_MAPE: " << (sub_total_gap_link_count) / max(1, total_link_count) * 100 <<
		"%,system_MPE: " << (sub_total_system_gap_count) / max(1, total_link_count) * 100 <<
		//"%,avg_tt = " << total_system_travel_time / max(0.1, total_system_demand) << "(min) " <<
		//",UE gap =" << total_system_UE_gap / max(0.00001, total_system_demand) << "(min)" <<
		//" = (" << total_system_UE_gap / max(0.00001, total_system_travel_time) * 100 << " %)"
		endl;

	dtalog.output() << "ODME #" << iteration_no
		<< ", link MAE= " << total_gap / max(1, total_link_count)
		<< ",link_MAPE: " << (sub_total_gap_link_count) / max(1, total_link_count) * 100 <<
		"%,system_MPE: " << (sub_total_system_gap_count) / max(1, total_link_count) * 100 <<
		"%,avg_tt = " << total_system_travel_time / max(0.1, total_system_demand) << "(min) " <<
		",UE gap =" << total_system_UE_gap / max(0.00001, total_system_demand) << "(min)" <<
		" = (" << total_system_UE_gap / max(0.00001, total_system_travel_time) * 100 << " %)"
		<< endl;
	double gap = sub_total_gap_link_count / max(1, total_link_count);
	system_gap = sub_total_system_gap_count / max(1, total_link_count);

	return gap;
}

double g_reset_and_update_sensor_link_volume_based_on_ODME_columns_complete_implementation(int number_of_links, int iteration_no, double& system_gap)
{
	float total_gap = 0;
	float sub_total_gap_link_count = 0;
	float sub_total_system_gap_count = 0;
	system_gap = 0;
	float sub_total_gap_P_count = 0;
	float sub_total_gap_A_count = 0;

	double total_system_travel_cost = 0;
	double total_system_travel_time = 0;
	double total_system_demand = 0;
	double total_system_UE_gap = 0;
	// reset the link volume
	for (int i = 0; i < number_of_links; ++i)
	{
		for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
		{

			if (g_link_vector[i].VDF_period[tau].obs_count >= 1)
			{
				// used in travel time calculation
				g_link_vector[i].PCE_volume_per_period[tau] = 0;
				g_link_vector[i].person_volume_per_period[tau] = 0;

				for (int at = 0; at < assignment.g_AgentTypeVector.size(); ++at)
					g_link_vector[i].person_volume_per_period_per_at[tau][at] = 0;
			}


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

		int zone_size = g_zone_vector.size();
		int tau_size = assignment.g_DemandPeriodVector.size();
		float PCE_ratio = assignment.g_AgentTypeVector[at].PCE;
		float OCC_ratio = assignment.g_AgentTypeVector[at].OCC;


#pragma omp parallel for
		for (int orig = 0; orig < zone_size; ++orig)  // o
		{
			int from_zone_sindex = g_zone_vector[orig].sindex;
			if (from_zone_sindex == -1)
				continue;


			int analysis_district_id = assignment.g_zone_seq_no_to_analysis_distrct_id_mapping[orig];
			std::map<int, CColumnPath>::iterator it;
			float link_volume_contributed_by_path_volume;
			int nl;

			std::map<int, CColumnPath>::iterator it_begin;
			std::map<int, CColumnPath>::iterator it_end;

			int column_vector_size;
			CColumnVector* p_column_pool;
			for (int dest = 0; dest < zone_size; ++dest) //d
			{
				int to_zone_sindex = g_zone_vector[dest].sindex;

				if (to_zone_sindex == -1)
					continue;

				for (int tau = 0; tau < tau_size; ++tau)  //tau
				{
					p_column_pool = &(assignment.g_column_pool[from_zone_sindex][to_zone_sindex][at][tau]);
					if (p_column_pool->od_volume > 0)
					{

						// continuous: type 0
						column_vector_size = p_column_pool->path_node_sequence_map.size();

						it_begin = p_column_pool->path_node_sequence_map.begin();
						it_end = p_column_pool->path_node_sequence_map.end();

						double least_cost = 999999;
						int least_cost_path_seq_no = -1;
						int least_cost_path_node_sum_index = -1;
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

						least_cost = 999999;
						path_seq_count = 0;

						it_begin = p_column_pool->path_node_sequence_map.begin();
						it_end = p_column_pool->path_node_sequence_map.end();
						for (it = it_begin; it != it_end; ++it)
						{
							total_system_demand += it->second.path_volume;
							path_toll = 0;
							path_gradient_cost = 0;
							path_distance = 0;
							path_travel_time = 0;

							for (int nl = 0; nl < it->second.m_link_size; ++nl)  // arc a
							{
								link_seq_no = it->second.path_link_vector[nl];
								link_travel_time = g_link_vector[link_seq_no].travel_time_per_period[tau];
								path_travel_time += link_travel_time;

							}

							it->second.path_toll = path_toll;
							it->second.path_travel_time = path_travel_time;
							total_system_travel_time += (it->second.path_travel_time * it->second.path_volume);

							if (column_vector_size == 1)  // only one path
							{
								break;
							}

							if (path_travel_time < least_cost)
							{
								least_cost = path_travel_time;
								least_cost_path_seq_no = it->second.path_seq_no;
								least_cost_path_node_sum_index = it->first;
							}
#pragma omp critical
							{
								total_system_travel_cost += (it->second.path_travel_time * it->second.path_volume);
							}
						} // end for each path




						for (it = it_begin; it != it_end; ++it)  // path k
						{
							link_volume_contributed_by_path_volume = it->second.path_volume;  // assign all OD flow to this first path

#pragma omp critical
							// add path volume to link volume
							for (nl = 0; nl < it->second.path_sensor_link_vector.size(); ++nl)  // arc a
							{
								link_seq_no = it->second.path_sensor_link_vector[nl];

								// MSA updating for the existing column pools
								// if iteration_index = 0; then update no flow discount is used (for the column pool case)
								{
									g_link_vector[link_seq_no].PCE_volume_per_period[tau] += link_volume_contributed_by_path_volume * PCE_ratio;
									g_link_vector[link_seq_no].person_volume_per_period[tau] += link_volume_contributed_by_path_volume * OCC_ratio;
									g_link_vector[link_seq_no].person_volume_per_period_per_at[tau][at] += link_volume_contributed_by_path_volume;  // pure volume, not consider PCE

									//g_link_vector[link_seq_no].person_volume_per_district_per_at[analysis_district_id][at] += link_volume_contributed_by_path_volume;  // pure volume, not consider PCE
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
		g_link_vector[i].calculate_dynamic_VDFunction(iteration_no, false, g_link_vector[i].vdf_type);
		for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); ++tau)  //tau
		{
			if (assignment.g_DemandPeriodVector[tau].number_of_demand_files == 0)
				continue;

			if (g_link_vector[i].VDF_period[tau].obs_count >= 1)  // with data
			{
				g_link_vector[i].VDF_period[tau].est_count_dev = g_link_vector[i].PCE_volume_per_period[tau] + g_link_vector[i].VDF_period[tau].preload - g_link_vector[i].VDF_period[tau].obs_count;

				if (dtalog.debug_level() == 2)
				{
					dtalog.output() << "link " << g_node_vector[g_link_vector[i].from_node_seq_no].node_id
						<< "->" << g_node_vector[g_link_vector[i].to_node_seq_no].node_id
						<< "obs:, " << g_link_vector[i].VDF_period[tau].obs_count << "est:, " << g_link_vector[i].PCE_volume_per_period[tau]
						<< "dev:," << g_link_vector[i].VDF_period[tau].est_count_dev << endl;
				}
				if (g_link_vector[i].VDF_period[tau].upper_bound_flag == 0)
				{
					total_gap += abs(g_link_vector[i].VDF_period[tau].est_count_dev);
					sub_total_gap_link_count += fabs(g_link_vector[i].VDF_period[tau].est_count_dev / g_link_vector[i].VDF_period[tau].obs_count);
					sub_total_system_gap_count += g_link_vector[i].VDF_period[tau].est_count_dev / g_link_vector[i].VDF_period[tau].obs_count;
				}
				else
				{  // upper bound constraints 
					if (g_link_vector[i].VDF_period[tau].est_count_dev > 0)
					{
						total_gap += abs(g_link_vector[i].VDF_period[tau].est_count_dev);
						sub_total_gap_link_count += fabs(g_link_vector[i].VDF_period[tau].est_count_dev / g_link_vector[i].VDF_period[tau].obs_count);
						sub_total_system_gap_count += g_link_vector[i].VDF_period[tau].est_count_dev / g_link_vector[i].VDF_period[tau].obs_count;
					}
				}
				total_link_count += 1;
			}
		}

	}
	//for (int orig = 0; orig < g_zone_vector.size(); ++orig)  // o
	//{
	//    if (g_zone_vector[orig].obs_attraction >= 1)  // with observation
	//    {
	//        g_zone_vector[orig].est_attraction_dev = g_zone_vector[orig].est_attraction - g_zone_vector[orig].obs_attraction;

	//        if (dtalog.debug_level() == 2)
	//        {
	//            dtalog.output() << "zone " << g_zone_vector[orig].zone_id << "A: obs:" << g_zone_vector[orig].obs_attraction
	//                << ",est:," << g_zone_vector[orig].est_attraction << ",dev:," << g_zone_vector[orig].est_attraction_dev << endl;
	//        }

	//        total_gap += abs(g_zone_vector[orig].est_attraction_dev);
	//        sub_total_gap_A_count += g_zone_vector[orig].est_attraction_dev / g_zone_vector[orig].obs_attraction;
	//    }

	//    if (g_zone_vector[orig].obs_production >= 1)  // with observation
	//    {
	//        g_zone_vector[orig].est_production_dev = g_zone_vector[orig].est_production - g_zone_vector[orig].obs_production;

	//        if (dtalog.debug_level() == 2)
	//        {
	//            dtalog.output() << "zone " << g_zone_vector[orig].zone_id << "P: obs:" << g_zone_vector[orig].obs_production
	//                << ",est:," << g_zone_vector[orig].est_production << ",dev:," << g_zone_vector[orig].est_production_dev << endl;
	//        }

	//        total_gap += abs(g_zone_vector[orig].est_production_dev);
	//        sub_total_gap_P_count += g_zone_vector[orig].est_production_dev / g_zone_vector[orig].obs_production;
	//    }
	//}

	assignment.summary_file << "ODME #" << iteration_no
		<< ", link MAE= " << total_gap / max(1, total_link_count)
		<< ",link_MAPE: " << (sub_total_gap_link_count) / max(1, total_link_count) * 100 <<
		"%,system_MPE: " << (sub_total_system_gap_count) / max(1, total_link_count) * 100 <<
		"%,avg_tt = " << total_system_travel_time / max(0.1, total_system_demand) << "(min) " <<
		",UE gap =" << total_system_UE_gap / max(0.00001, total_system_demand) << "(min)" <<
		" = (" << total_system_UE_gap / max(0.00001, total_system_travel_time) * 100 << " %)"
		<< endl;

	dtalog.output() << "ODME #" << iteration_no
		<< ", link MAE= " << total_gap / max(1, total_link_count)
		<< ",link_MAPE: " << (sub_total_gap_link_count) / max(1, total_link_count) * 100 <<
		"%,system_MPE: " << (sub_total_system_gap_count) / max(1, total_link_count) * 100 <<
		"%,avg_tt = " << total_system_travel_time / max(0.1, total_system_demand) << "(min) " <<
		",UE gap =" << total_system_UE_gap / max(0.00001, total_system_demand) << "(min)" <<
		" = (" << total_system_UE_gap / max(0.00001, total_system_travel_time) * 100 << " %)"
		<< endl;
	double gap = sub_total_gap_link_count / max(1, total_link_count);
	system_gap = sub_total_system_gap_count / max(1, total_link_count);

	return gap;
}
void g_update_gradient_cost_and_assigned_flow_in_column_pool(Assignment& assignment, int inner_iteration_number, bool b_sensitivity_analysis_flag)
{
	double total_system_cost_gap = 0;
	float total_relative_gap = 0;
	double total_system_travel_cost = 0;
	double total_system_travel_time = 0;
	double total_system_demand = 0;

	// we can have a recursive formulat to reupdate the current link volume by a factor of k/(k+1),
	// and use the newly generated path flow to add the additional 1/(k+1)
	g_reset_and_update_link_volume_based_on_columns(g_link_vector.size(), inner_iteration_number, false, b_sensitivity_analysis_flag);

	if (b_sensitivity_analysis_flag == true)  // check estimation counts
	{
		for (int i = 0; i < g_link_vector.size(); ++i)
		{
			for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
			{

				if (g_link_vector[i].VDF_period[tau].obs_count >= 1)  // with data
				{
					g_link_vector[i].VDF_period[tau].est_count_dev = g_link_vector[i].PCE_volume_per_period[tau] + g_link_vector[i].VDF_period[tau].preload - g_link_vector[i].VDF_period[tau].obs_count;
				}
			}
		}
	}
	// step 4: based on newly calculated path volumn, update volume based travel time, and update volume based resource balance, update gradie
	double total_distance = 0;
	update_link_travel_time_and_cost(inner_iteration_number, total_distance);
	// step 0


//   assignment.summary_file << ",iteration,key,o,d,at,tau,volume,"<< endl;
	//step 1: calculate shortest path at inner iteration of column flow updating
//#pragma omp parallel for
	for (int orig = 0; orig < g_zone_vector.size(); ++orig)  // o
	{
		CColumnVector* p_column_pool;
		std::map<int, CColumnPath>::iterator it, it_begin, it_end;
		int column_vector_size;

		double least_gradient_cost = 999999;
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

		int from_zone_sindex = g_zone_vector[orig].sindex;
		if (from_zone_sindex == -1)
			continue;

		for (int dest = 0; dest < g_zone_vector.size(); ++dest) //d
		{
			int to_zone_sindex = g_zone_vector[dest].sindex;
			if (to_zone_sindex == -1)
				continue;

			for (int at = 0; at < assignment.g_AgentTypeVector.size(); ++at)  //m
			{
				for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); ++tau)  //tau
				{
					p_column_pool = &(assignment.g_column_pool[from_zone_sindex][to_zone_sindex][at][tau]);
					if (p_column_pool->od_volume > 0)
					{

						double diff = p_column_pool->od_volume - p_column_pool->prev_od_volume;

						if (b_sensitivity_analysis_flag && inner_iteration_number >= 1)
						{
							if (diff < -0.0001 || diff > 0.0001)
							{
								int idebug = 1;
							}

							if (inner_iteration_number >= 1)
								diff = p_column_pool->od_volume - p_column_pool->od_volume_per_iteration_map[inner_iteration_number - 1];

							if (diff < -0.0001 || diff > 0.0001)
							{
								int idebug = 1;
							}

						}

						if (b_sensitivity_analysis_flag)
						{
							if (g_zone_vector[orig].zone_id == 6 && g_zone_vector[dest].zone_id == 2)
							{
								int idebug = 1;
							}
						}

						p_column_pool->prev_od_volume = p_column_pool->od_volume;

						column_vector_size = p_column_pool->path_node_sequence_map.size();


						if (b_sensitivity_analysis_flag)
						{
							p_column_pool->od_volume_per_iteration_map[inner_iteration_number] = p_column_pool->od_volume;
						}

						// scan through the map with different node sum for different paths
						/// step 1: update gradient cost for each column path

						least_gradient_cost = 999999;
						least_gradient_cost_path_seq_no = -1;
						least_gradient_cost_path_node_sum_index = -1;
						path_seq_count = 0;

						it_begin = p_column_pool->path_node_sequence_map.begin();
						it_end = p_column_pool->path_node_sequence_map.end();

						bool least_path_passing_improvement_flag = false;
						for (it = it_begin; it != it_end; ++it)
						{
							path_toll = 0;
							path_gradient_cost = 0;
							path_distance = 0;
							path_travel_time = 0;

							/// <summary>
							/// 
							if(b_sensitivity_analysis_flag && assignment.g_AgentTypeVector[at].real_time_information ==1 )  //real time information
							{
								for (int nl = 0; nl < it->second.m_link_size; ++nl)  // arc a
								{
									// step 3.3 link flow gradient
									link_seq_no = it->second.path_link_vector[nl];

									if (g_link_vector[link_seq_no].VDF_period[tau].network_design_flag == -1)  // affected by capacity reduction supply side scenario
									{
										
										p_column_pool->OD_impact_flag = 1;
										p_column_pool->at_od_impacted_flag_map[at] = true;

										break;

									}
								}
							}

							/// </summary>
							/// <param name="assignment"></param>
							/// <param name="inner_iteration_number"></param>
							/// <param name="b_sensitivity_analysis_flag"></param>

							for (int nl = 0; nl < it->second.m_link_size; ++nl)  // arc a
							{
								link_seq_no = it->second.path_link_vector[nl];
								path_toll += g_link_vector[link_seq_no].VDF_period[tau].toll[at];
								path_distance += g_link_vector[link_seq_no].link_distance_VDF;
								link_travel_time = g_link_vector[link_seq_no].travel_time_per_period[tau];
								path_travel_time += link_travel_time;

								path_gradient_cost += g_link_vector[link_seq_no].get_generalized_first_order_gradient_cost_of_second_order_loss_for_agent_type(tau, at);
							}

							it->second.path_toll = path_toll;
							it->second.path_travel_time = path_travel_time;
							it->second.path_gradient_cost = path_gradient_cost;


							if (b_sensitivity_analysis_flag == false)
								it->second.path_time_per_iteration_map[inner_iteration_number] = path_travel_time;
							else  // SA mode
								it->second.path_time_per_iteration_SA_map[inner_iteration_number] = path_travel_time;

#pragma omp critical
							{
								total_system_travel_time += (it->second.path_travel_time * it->second.path_volume);
								total_system_demand += it->second.path_volume;


								if (column_vector_size == 1)  // only one path
								{
									total_system_travel_cost += (it->second.path_gradient_cost * it->second.path_volume);

								}
							}

							if (path_gradient_cost < least_gradient_cost)
							{
								least_gradient_cost = path_gradient_cost;
								least_gradient_cost_path_seq_no = it->second.path_seq_no;
								least_gradient_cost_path_node_sum_index = it->first;
								if (it->second.impacted_path_flag)
								{
									least_path_passing_improvement_flag = 1;
								}
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

									//if(it->second.path_gradient_cost_difference >0.0001f)
									{
										it->second.path_gradient_cost_relative_difference = it->second.path_gradient_cost_difference / max(0.0001, least_gradient_cost);
									}

#pragma omp critical
									{
										total_system_cost_gap += (it->second.path_gradient_cost_difference * it->second.path_volume);
										total_system_travel_cost += (it->second.path_gradient_cost * it->second.path_volume);
									}

									if (b_sensitivity_analysis_flag == true)  // SA stages
									{

										// during SA stage,
										// regular drivers do not switch at all (except for those encounter the dramatic induced delay)
										// real time and DMS users only switch when they see the DMS or pass through the incident on the original routes 

										step_size = 1.0 / (inner_iteration_number + 2) * p_column_pool->od_volume; // small changes 

									}
									else
									{  // column updating step size
										step_size = 1.0 / (inner_iteration_number + 2) * p_column_pool->od_volume;
									}


									previous_path_volume = it->second.path_volume; //b

									double flow_shift = step_size * max(0.0000, it->second.path_gradient_cost_relative_difference);  //c, must be positive


									if (it->second.path_gradient_cost_relative_difference > 3 * 60)  // difference more than 3 hours
									{
										flow_shift = it->second.path_volume;  //switch out
									}

									if (flow_shift > it->second.path_volume)
									{
										flow_shift = it->second.path_volume;
									}

									if (flow_shift >= 0.000001)
									{
										int idebug = 1;
									}

									// special condition for sensitivity iteration 0: reset the route flow switing size  as 0: no route flow switching is alloed for this base case: 
									// Discussion with Chongnan. 
									if (b_sensitivity_analysis_flag == true)
									{
										if (p_column_pool->OD_impact_flag == 0 )
										// by default no switching
										{
											step_size = 0;
											flow_shift = 0;
										}
										else
										{  // p_column_pool->OD_impact_flag ==1 
											// one more layer of checking
											if (p_column_pool->at_od_impacted_flag_map.find(at) == p_column_pool->at_od_impacted_flag_map.end())
											{ // for RT users, only one condition is used to mark at_od_impacted_flag_map
											  // for DMS users, two conditions are used to mark at_od_impacted_flag_map
												step_size = 0;
												flow_shift = 0;
											}

										}
									}
									
									//recall that it->second.path_gradient_cost_difference >=0
									// step 3.1: shift flow from nonshortest path to shortest path
									it->second.path_volume = max(0.0, it->second.path_volume - flow_shift);
									//d
									// 
									//we use min(step_size to ensure a path is not switching more than 1/n proportion of flow
									it->second.path_switch_volume = (previous_path_volume - it->second.path_volume);

									// d-b
									// should be nonnegative
									total_switched_out_path_volume += (previous_path_volume - it->second.path_volume);

									if (fabs(total_switched_out_path_volume) > 0.00001)
									{
										int idebug = 1;
									}
								}
							}

							//step 3.2 consider least cost path, receive all volume shifted from non-shortest path
							if (least_gradient_cost_path_seq_no != -1 && p_column_pool->path_node_sequence_map.find(least_gradient_cost_path_node_sum_index) != p_column_pool->path_node_sequence_map.end())
							{
								if (least_gradient_cost_path_node_sum_index < 100)
								{
									int i_debug = 1;
								}
								p_column_pool->path_node_sequence_map[least_gradient_cost_path_node_sum_index].path_volume += total_switched_out_path_volume;

								if (b_sensitivity_analysis_flag == false)
									p_column_pool->path_node_sequence_map[least_gradient_cost_path_node_sum_index].path_volume_per_iteration_map[inner_iteration_number] = p_column_pool->path_node_sequence_map[least_gradient_cost_path_node_sum_index].path_volume;
								else
									p_column_pool->path_node_sequence_map[least_gradient_cost_path_node_sum_index].path_volume_per_iteration_SA_map[inner_iteration_number] = p_column_pool->path_node_sequence_map[least_gradient_cost_path_node_sum_index].path_volume;
#pragma omp critical
								{

									total_system_travel_cost += (p_column_pool->path_node_sequence_map[least_gradient_cost_path_node_sum_index].path_gradient_cost *
										p_column_pool->path_node_sequence_map[least_gradient_cost_path_node_sum_index].path_volume);
								}
							}

						}
						// record path flow for all paths( including shortst path and non_shortest path)
						for (it = it_begin; it != it_end; ++it)
						{

							if (b_sensitivity_analysis_flag == false)
								it->second.path_volume_per_iteration_map[inner_iteration_number] = it->second.path_volume;
							else  //SA mode
								it->second.path_volume_per_iteration_SA_map[inner_iteration_number] = it->second.path_volume;


							if (b_sensitivity_analysis_flag == true)
							{
								if(inner_iteration_number == 0)
								it->second.path_volume_before_sa = it->second.path_volume;

								if (inner_iteration_number >= 1)
								it->second.path_volume_after_sa = it->second.path_volume;
							}
						}


					}
				}
			}
		}
	}

	double avg_travel_time = total_system_travel_time / max(0.001, total_system_demand);

	assignment.summary_file << "column updating: iteration= " << inner_iteration_number << ", avg travel time = " <<
		avg_travel_time << "(min), optimization obj = " << total_system_cost_gap
		<< ",Relative_gap=" << total_system_cost_gap * 100.0 / max(0.00001, total_system_travel_cost) << " %" << endl;

	dtalog.output() << "column updating: iteration= " << inner_iteration_number << ", avg travel time = " <<
		avg_travel_time << "(min), optimization obj = " << total_system_cost_gap
		<< ",Relative_gap=" << total_system_cost_gap * 100.0 / max(0.00001, total_system_travel_cost) << " %" << endl;

	string stage_str;
	stage_str = "column updating";

	if (b_sensitivity_analysis_flag)
		stage_str = "sensitivity analaysis";

	assignment.summary_file2 << stage_str.c_str() << ",iteration," << inner_iteration_number <<
		",total_system_demand," << total_system_demand <<
		",avg travel time," << avg_travel_time <<
		",optimization obj," << total_system_cost_gap <<
		",relative_gap," << total_system_cost_gap * 100.0 / max(0.00001, total_system_travel_cost) << "," << endl;
}

void g_classification_in_column_pool(Assignment& assignment)
{
	int impact_OD_counts = 0;
	int impact_OD_counts_detour = 0;
	//#pragma omp parallel for
	for (int orig = 0; orig < g_zone_vector.size(); ++orig)  // o
	{
		CColumnVector* p_column_pool;
		std::map<int, CColumnPath>::iterator it, it_begin, it_end;
		int column_vector_size;
		int from_zone_sindex = g_zone_vector[orig].sindex;
		if (from_zone_sindex == -1)
			continue;

		int link_seq_no;

		for (int dest = 0; dest < g_zone_vector.size(); ++dest) //d
		{
			int to_zone_sindex = g_zone_vector[dest].sindex;
			if (to_zone_sindex == -1)
				continue;

			for (int at = 0; at < assignment.g_AgentTypeVector.size(); ++at)  //m
			{
				//if (assignment.g_AgentTypeVector[at].real_time_information != 0)  // users with information, continue;
				//	continue; 

				for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); ++tau)  //tau
				{

					p_column_pool = &(assignment.g_column_pool[from_zone_sindex][to_zone_sindex][at][tau]);

					if (p_column_pool->od_volume > 0)
					{

						column_vector_size = p_column_pool->path_node_sequence_map.size();


						// scan through the map with different node sum for different paths
						/// step 1: update gradient cost for each column path

						it_begin = p_column_pool->path_node_sequence_map.begin();
						it_end = p_column_pool->path_node_sequence_map.end();

						bool least_path_passing_improvement_flag = false;
						// scan all paths in this OD pair
						int path_count = 0;
						int network_design_path_count = 0;

						for (it = it_begin; it != it_end; ++it)  // for each path
						{
							for (int nl = 0; nl < it->second.m_link_size; ++nl)  // arc a
							{
								link_seq_no = it->second.path_link_vector[nl];

								if (g_link_vector[link_seq_no].VDF_period[tau].network_design_flag == -1 /* sa impacted, but not DMS location*/)  // screening condition 1: passing through the network design location
								{
									it->second.impacted_path_flag = 1;

									// to be revised: passing through work zone, and with signal timing enhancemnets
								}

							}

							if (it->second.impacted_path_flag)
								network_design_path_count++;

							path_count++;
						}


						if (network_design_path_count >= 1)
						{
							p_column_pool->OD_impact_flag = 1;
						}

					}
				} // for each tau
			}// for each agent type mode
		} // for each d
	}
	string stage_str;
	stage_str = "classification";

//	assignment.summary_file2 << stage_str.c_str() << ",impact_OD_counts," << impact_OD_counts <<
//		",impact_OD_counts_with_detour," << impact_OD_counts_detour << endl;



}
CColumnPath* g_add_new_column(CColumnVector* pColumnVector, std::vector <int> link_seq_vector, float volume)
{
	int l_node_size = link_seq_vector.size() + 1;
	int l_link_size = link_seq_vector.size();

	//node seq vector for each ODK
	int temp_path_node_vector[MAX_LINK_SIZE_IN_A_PATH];
	//node seq vector for each ODK
	int temp_path_link_vector[MAX_LINK_SIZE_IN_A_PATH];

	temp_path_node_vector[0] = g_link_vector[link_seq_vector[0]].from_node_seq_no;
	int node_sum = 0;
	for (int l = 0; l < link_seq_vector.size(); l++)
	{
		temp_path_link_vector[l] = link_seq_vector[l];
		temp_path_node_vector[l + 1] = g_link_vector[link_seq_vector[l]].to_node_seq_no;
		node_sum += temp_path_link_vector[l] + temp_path_node_vector[l + 1];
	}


	// add this unique path

	pColumnVector->path_node_sequence_map[node_sum].path_seq_no = pColumnVector->path_node_sequence_map.size();;

	pColumnVector->path_node_sequence_map[node_sum].AllocateVector(
		l_node_size,
		temp_path_node_vector,
		l_link_size,
		temp_path_link_vector, false);

	pColumnVector->path_node_sequence_map[node_sum].path_volume = volume;  // we add 1/K * OD volume to a new path or an existing path with same node sum.

	return &(pColumnVector->path_node_sequence_map[node_sum]);  // a pointer to be used by other programs
}

void g_column_pool_optimization(Assignment& assignment, int column_updating_iterations, bool sensitivity_analysis_flag = false)
{
	assignment.summary_file << "column updating" << endl;

	// column_updating_iterations is internal numbers of column updating
	for (int n = 0; n < column_updating_iterations; ++n)
	{
		g_update_gradient_cost_and_assigned_flow_in_column_pool(assignment, n, sensitivity_analysis_flag);

		if (dtalog.debug_level() >= 3)
		{
			for (int i = 0; i < g_link_vector.size(); ++i)
			{
				dtalog.output() << "link: " << g_node_vector[g_link_vector[i].from_node_seq_no].node_id << "-->"
					<< g_node_vector[g_link_vector[i].to_node_seq_no].node_id << ", "
					<< "flow count:" << g_link_vector[i].PCE_volume_per_period[0] << endl;
			}
		}
	}
}

void g_column_regeneration(Assignment& assignment);
void g_reset_link_volume_in_master_program_without_columns(int number_of_links, int iteration_index, bool b_self_reducing_path_volume);
void g_reset_link_volume_for_all_processors();
void g_fetch_link_volume_for_all_processors();

void g_column_pool_route_modification(Assignment& assignment, int inner_iteration_number)  // for DMS users
{

	//step 1: calculate shortest path at inner iteration of column flow updating
#pragma omp parallel for
	for (int orig = 0; orig < g_zone_vector.size(); ++orig)  // o
	{
		CColumnVector* p_column_pool;

		std::map<int, CColumnPath>::iterator it, it_begin, it_end;
		int column_vector_size;

		int path_seq_count = 0;

		double path_toll = 0;
		double path_gradient_cost = 0;
		double path_distance = 0;
		double path_travel_time = 0;
		int link_seq_no;

		int from_zone_sindex = g_zone_vector[orig].sindex;
		if (from_zone_sindex == -1)
			continue;

		for (int dest = 0; dest < g_zone_vector.size(); ++dest) //d
		{
			int to_zone_sindex = g_zone_vector[dest].sindex;

			if (to_zone_sindex == -1)
				continue;

			for (int at = 0; at < assignment.g_AgentTypeVector.size(); ++at)  //m
			{
				if (assignment.g_AgentTypeVector[at].real_time_information == 2)   // case of DMS 
				{
					for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); ++tau)  //tau
				{
					p_column_pool = &(assignment.g_column_pool[from_zone_sindex][to_zone_sindex][at][tau]);
					if (p_column_pool->od_volume > 0)
					{
							column_vector_size = p_column_pool->path_node_sequence_map.size();

							// scan through the map with different node sum for different paths


							it_begin = p_column_pool->path_node_sequence_map.begin();
							it_end = p_column_pool->path_node_sequence_map.end();

							//test condition 1: passing through information zone
							for (it = it_begin; it != it_end; ++it)  // scan each first-stage original path
							{
								path_seq_count = 0;
								bool b_passing_information_zone = false;
								int new_orig_zone_id = 0;


								//test condition 2: passing through capacity impact area
								bool b_passing_capacity_impact_area = false;

								if (it->second.path_volume < 0.00001)  // positive flow
									continue;

								std::vector <int> link_seq_vector;  // reset local variable for each path/column

								for (int nl = 0; nl < it->second.m_link_size; ++nl)  // arc a  // first nl
								{
									link_seq_no = it->second.path_link_vector[nl];
									CLink* p_current_link = &(g_link_vector[link_seq_no]);


									if (b_passing_information_zone == false &&
										p_current_link->VDF_period[tau].network_design_flag ==2 /*DMS link*/ &&
										assignment.node_seq_no_2_zone_id_mapping.find(p_current_link->to_node_seq_no) != assignment.node_seq_no_2_zone_id_mapping.end())
										  // this node been defined as zone
									{

										cout << "nl= " << nl << " link sequence = " << link_seq_no << endl;
										p_column_pool->OD_impact_flag = 1;
										p_column_pool->at_od_impacted_flag_map[at];

										int zone_id = assignment.node_seq_no_2_zone_id_mapping[p_current_link->to_node_seq_no];

										if (assignment.g_zoneid_to_zone_seq_no_mapping.find(zone_id) == assignment.g_zoneid_to_zone_seq_no_mapping.end())  // not found
											continue; 

										int zone_no = assignment.g_zoneid_to_zone_seq_no_mapping[zone_id];

										b_passing_information_zone = true;
										new_orig_zone_id = zone_id;  // zone id to zone no.

											for (int nl2 = 0; nl2 <= nl; ++nl2)  // arc a
											{  // copy the existing link sequence up to the downstream node id corresponding to the info zone

												// part A subpath
												int l_link_seq_no = it->second.path_link_vector[nl2];
												link_seq_vector.push_back(l_link_seq_no);
											}
									}

									if (g_link_vector[link_seq_no].VDF_period[tau].network_design_flag == -1)  // affected by capacity reduction
									{
										b_passing_capacity_impact_area = true;
										it->second.impacted_path_flag = 1; // impacted
									}

								}

								// two conditions are satisfied

								if (b_passing_capacity_impact_area == true && b_passing_information_zone == true)
								{

									CColumnVector* p_2_stage_column_pool;

									// first step: fetch another column from the newly defined origin_zone id (i.e. information zone)
									int info_orig = assignment.g_zoneid_to_zone_seq_no_mapping[new_orig_zone_id];

									int from_zone_sindex = g_zone_vector[info_orig].sindex;
									if (from_zone_sindex == -1)
										continue;

									int to_zone_sindex = g_zone_vector[dest].sindex;
									if (to_zone_sindex == -1)
										continue;

									//step 2: fetch the related column pool from the information node/zone
									p_2_stage_column_pool = &(assignment.g_column_pool[from_zone_sindex][to_zone_sindex][at][tau]);  // we come from info_orig but going to  the same destination with same at, and assignment period tau 
									//             scan through the map with different node sum for different continuous paths
									// part C subpath

									if (p_2_stage_column_pool->path_node_sequence_map.size() == 0)  // if there is no path available int the precomputed column pool
									{
										cout << "if there is no path available int the precomputed column pool";
										g_program_stop();
										continue;  // 

									}


									std::map<int, CColumnPath>::iterator it2, it_begin2, it_end2;

									it_begin2 = p_2_stage_column_pool->path_node_sequence_map.begin();
									it_end2 = p_2_stage_column_pool->path_node_sequence_map.end();

									bool b_sa_column_found_flag = false;
									// check non-diverted path travel time 
									float non_diverted_path_travel_time = 999999; // default value
									for (it2 = it_begin2; it2 != it_end2; ++it2)
									{
										for (int nl2 = 1; nl2 < it2->second.m_link_size; ++nl2)  // arc a // starting from 1 to exclude virtual link at the beginning
										{
											// compute it2->second.path_travel_time
											it2->second.path_travel_time += g_link_vector[it2->second.path_link_vector[nl2]].travel_time_per_period[tau];

										}

										// diversion status determination
										bool b_diversion_flag = true;


										for (int nl2 = 1; nl2 < it2->second.m_link_size; ++nl2)  // arc a // starting from 1 to exclude virtual link at the beginning
											if (g_link_vector[it2->second.path_link_vector[nl2]].VDF_period[tau].network_design_flag <= -1)  // affected by capacity reduction
											{
												b_diversion_flag = false;
											}

										if (!b_diversion_flag)
										{
											float abs_BR_travel_time = 3;
											float ratio_BR_travel_time = 0.15;
											if (it2->second.path_travel_time < non_diverted_path_travel_time - abs_BR_travel_time  && it2->second.path_travel_time < (non_diverted_path_travel_time*(1- ratio_BR_travel_time)))
												non_diverted_path_travel_time = it2->second.path_travel_time;
										}

									}
								

									// step 3: we can still have k-path from the info zone to to final destination so we need to select alternative diverted path with lower travel time compared to the historical route
									for (it2 = it_begin2; it2 != it_end2; ++it2)
									{
										bool b_diversion_flag = true;

																		
										for (int nl2 = 1; nl2 < it2->second.m_link_size; ++nl2)  // arc a // starting from 1 to exclude virtual link at the beginning
											{


												if (g_link_vector[it2->second.path_link_vector[nl2]].VDF_period[tau].network_design_flag <= -1)  // affected by capacity reduction
												{
													b_diversion_flag = false;
												}

											}

											
										if (b_diversion_flag == true && b_sa_column_found_flag == false) // second sub path: 3 conditions: diverted. and first alternative path
										{
											if(it2->second.path_travel_time < non_diverted_path_travel_time)
											{ //select alternative diverted path with lower travel time compared to the historical route
												for (int nl2 = 1; nl2 < it2->second.m_link_size; ++nl2)  // arc a // starting from 1 to exclude virtual link at the beginning
												{
													link_seq_vector.push_back(it2->second.path_link_vector[nl2]);
													// construct sub path C
													cout << "reroute B: l=" << nl2 << "," << g_node_vector[g_link_vector[it2->second.path_link_vector[nl2]].to_node_seq_no].node_id << endl;

												}

											b_sa_column_found_flag = true;
											}
										}

									}

									// re-assmble the path vector of path_link_vector only if a feasible sa column is found.
									if (b_sa_column_found_flag  == true &&  it->second.path_link_vector != NULL)
									{
										// copy the updated path (stage1 + stage 2) back to the path link vector 
										bool b_optional_DMS_flag = true;

										//required diversion									
										if(b_optional_DMS_flag== false)
										{
										delete it->second.path_link_vector;
										it->second.path_link_vector = new int[link_seq_vector.size()];
										for (int l = 0; l < link_seq_vector.size(); l++)
										{
											it->second.path_link_vector[l] = link_seq_vector[l];
										}

										it->second.m_link_size = link_seq_vector.size();

										// copy the updated path (stage1 + stage 2) back to the path node vector 
										delete it->second.path_node_vector;
										it->second.path_node_vector = new int[link_seq_vector.size() + 1];


										// first node
										it->second.path_node_vector[0] = g_link_vector[link_seq_vector[0]].from_node_seq_no;


										// remaining nodes to the end of path

										for (int l = 0; l < link_seq_vector.size(); l++)
										{
											it->second.path_node_vector[l + 1] = g_link_vector[link_seq_vector[l]].to_node_seq_no;
											cout << "reroute: l=" << l << "," << g_node_vector[g_link_vector[link_seq_vector[l]].to_node_seq_no].node_id << endl; 

										}
										it->second.m_node_size = link_seq_vector.size() + 1;

										it->second.impacted_path_flag = 2; // diverted 
										}
										else
										{  //optional detour

											// add a new path to p_column_pool->path_node_sequence_map
											// simple logit model 

											float base_path_travel_time = 0;
											for (int nl0 = 0; nl0 < it->second.m_link_size; ++nl0)  // arc a // starting from 1 to exclude virtual link at the beginning
											{
												// compute it2->second.path_travel_time
												base_path_travel_time += g_link_vector[it->second.path_link_vector[nl0]].travel_time_per_period[tau];

											}
											float alt_path_travel_time = 0;
											for (int nl3 = 1; nl3 < link_seq_vector.size(); ++nl3)  // arc a // starting from 1 to exclude virtual link at the beginning
											{
												// compute it2->second.path_travel_time
												alt_path_travel_time += g_link_vector[link_seq_vector[nl3]].travel_time_per_period[tau];

											}

											float beta = -0.1; 
											float prob_base = exp(beta * base_path_travel_time) / (exp(beta * base_path_travel_time) + exp(beta * alt_path_travel_time));
											float base_volume = it->second.path_volume;
											it->second.path_volume = base_volume * prob_base;
											CColumnPath* pColumnPath = g_add_new_column(p_column_pool,link_seq_vector, base_volume * (1- prob_base));

											pColumnPath->impacted_path_flag = 3; // diverted and new colume

										}
									}
									p_2_stage_column_pool->information_type = 1;
								
									// path D as A + C


								}  // two conditions satisified

							}  //end of scanning for the first stage path in the column pool

						} // agent type is real time agent type

					}  // with positve OD volume


				} // tau
			}  //agent type
		} //dest
	}  // orig


	dtalog.output() << " updating";


}


void g_rt_info_column_generation(Assignment* p_assignment, float current_time_in_min, int recording_flag = 0)
{
	//dtalog.output() << "Begin the computing of " << g_NetworkForRTSP_vector.size() << " RTSP networks in CPU." << endl;

	clock_t start_t0, end_t0, total_t0;

	start_t0 = clock();
#pragma omp parallel for  // step 3: C++ open mp automatically create n threads., each thread has its own computing thread on a cpu core
	for (int blk = 0; blk < g_NetworkForRTSP_vector.size(); ++blk)
	{

	NetworkForSP* pNetwork = g_NetworkForRTSP_vector[blk];

		if (assignment.g_DemandPeriodVector[pNetwork->m_tau].starting_time_slot_no * MIN_PER_TIMESLOT > current_time_in_min)  // RT network is for a later time interval
			continue;

		pNetwork->optimal_backward_label_correcting_from_destination(blk, p_assignment, current_time_in_min, pNetwork->m_RT_dest_zone, pNetwork->m_RT_dest_node, -1, recording_flag);

	}
	end_t0 = clock();

	total_t0 = (end_t0 - start_t0);
	int second = total_t0 / 1000.0;
	int min = second / 60;
	int sec = second - min * 60;
	//dtalog.output() << "CPU Running Time for RT shortest path: " << min << " min " << sec << " sec" << endl;
	assignment.summary_file << ", RT shortest path at time =," << current_time_in_min << "min" << endl;

}


void g_column_pool_activity_scheduling(Assignment& assignment, int inner_iteration_number)
{

	//step 1: calculate shortest path at inner iteration of column flow updating

	for (int orig = 0; orig < g_zone_vector.size(); ++orig)  // o
	{
		int from_zone_sindex = g_zone_vector[orig].sindex;
		if (from_zone_sindex == -1)
			continue;

		CColumnVector* p_column_pool;
		int path_seq_count = 0;

		double path_toll = 0;
		double path_gradient_cost = 0;
		double path_distance = 0;
		double path_travel_time = 0;

		for (int dest = 0; dest < g_zone_vector.size(); ++dest) //d
		{
			int to_zone_sindex = g_zone_vector[dest].sindex;

			if (to_zone_sindex == -1)
				continue;

			for (int at = 0; at < assignment.g_AgentTypeVector.size(); ++at)  //m
			{
				for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); ++tau)  //tau
				{
					p_column_pool = &(assignment.g_column_pool[from_zone_sindex][to_zone_sindex][at][tau]);
					if (p_column_pool->od_volume > 0)
					{
						if (p_column_pool->activity_zone_no_vector.size())   // case of activity zones
						{
							p_column_pool->path_node_sequence_map.clear();  // remove existing single OD pair based routes


							std::vector <int> link_seq_vector;
							// for each origin and detination pair in activity zone no to perform routing continuously

							for (int az = 0; az < p_column_pool->activity_zone_no_vector.size() - 1; az++) // key step: go through each activty OD pair
							{ // 0 will the origin
								// last one will destination
								int aat = p_column_pool->activity_agent_type_no_vector[az];

								CColumnVector* p_2_stage_column_pool;

								int activity_orig = p_column_pool->activity_zone_no_vector[az];
								int activity_dest = p_column_pool->activity_zone_no_vector[az + 1];

								int from_zone_sindex = g_zone_vector[activity_orig].sindex;

								if (from_zone_sindex == -1)
									continue;

								int to_zone_sindex = g_zone_vector[activity_dest].sindex;
								if (to_zone_sindex == -1)
									continue;

								//step 2: fetch the related column pool from the information node/zone
								p_2_stage_column_pool = &(assignment.g_column_pool[from_zone_sindex][to_zone_sindex][aat][tau]);  // we come from info_orig but going to  the same destination with same at, and assignment period tau 
								//             scan through the map with different node sum for different continuous paths

								std::map<int, CColumnPath>::iterator it2, it_begin2, it_end2;

								it_begin2 = p_2_stage_column_pool->path_node_sequence_map.begin();
								it_end2 = p_2_stage_column_pool->path_node_sequence_map.end();


								for (it2 = it_begin2; it2 != it_end2; ++it2)  // we can still have k-path from the info zone to to final destination so we need to random select one 
								{

									for (int nl = 1; nl < it2->second.m_link_size - 1; ++nl)  // arc a // exclude virtual link in the beginning and at the end;
									{
										link_seq_vector.push_back(it2->second.path_link_vector[nl]);
									}
									break; // only connect with the first available second stage path
								}


							}

							if (link_seq_vector.size() == 0)
							{
								int i_debug = 1;
								continue;
							}

							int node_sum = 0;
							for (int l = 0; l < link_seq_vector.size(); l++)
							{
								node_sum += link_seq_vector[l];
							}


							// add this unique path  // later we can add k activity paths
							int path_count = p_column_pool->path_node_sequence_map.size();
							p_column_pool->path_node_sequence_map[node_sum].path_seq_no = path_count;
							p_column_pool->path_node_sequence_map[node_sum].path_volume = p_column_pool->od_volume;
							p_column_pool->path_node_sequence_map[node_sum].path_toll = 0;


							p_column_pool->path_node_sequence_map[node_sum].path_link_vector = new int[link_seq_vector.size()];
							p_column_pool->path_node_sequence_map[node_sum].path_node_vector = new int[link_seq_vector.size() + 1];

							for (int l = 0; l < link_seq_vector.size(); l++)
							{
								p_column_pool->path_node_sequence_map[node_sum].path_link_vector[l] = link_seq_vector[l];
								p_column_pool->path_node_sequence_map[node_sum].path_link_STL_vector.push_back(link_seq_vector[l]);
							}
							p_column_pool->path_node_sequence_map[node_sum].m_link_size = link_seq_vector.size();

							// copy the updated path (stage1 + stage 2) back to the path node vector 

							// first node
							p_column_pool->path_node_sequence_map[node_sum].path_node_vector[0] = g_link_vector[link_seq_vector[0]].from_node_seq_no;

							// remaining nodes to the end of path

							for (int l = 0; l < link_seq_vector.size(); l++)
							{
								p_column_pool->path_node_sequence_map[node_sum].path_node_vector[l + 1] = g_link_vector[link_seq_vector[l]].to_node_seq_no;
							}
							p_column_pool->path_node_sequence_map[node_sum].m_node_size = link_seq_vector.size() + 1;
						} //end of conditions for activity chain
					}  // with positve OD volume


				} // tau
			}  //agent type
		} //dest
	}  // orig


	dtalog.output() << " updating";
}
