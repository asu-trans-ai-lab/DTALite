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

// updates for OD re-generations
void Assignment::Demand_ODME(int OD_updating_iterations)
{
	int sensor_count = 0;
	if (OD_updating_iterations >= 1)
	{
		// step 1: read sensor_data.csv
		CDTACSVParser parser_measurement;

		if (parser_measurement.OpenCSVFile("sensor_data.csv", true))
		{
			while (parser_measurement.ReadRecord())  // if this line contains [] mark, then we will also read field headers.
			{

					int from_node_id;
					if (!parser_measurement.GetValueByFieldName("from_node_id", from_node_id))
						continue;

					int to_node_id;
					if (!parser_measurement.GetValueByFieldName("to_node_id", to_node_id))
						continue;

					
					int scenario_index = 0;
					parser_measurement.GetValueByFieldName("scenario_index", scenario_index);

					if (g_scenario_summary_map.find(scenario_index) == g_scenario_summary_map.end())
					{
						dtalog.output() << "[ERROR] scenario_index =  " << scenario_index << " in file sensor_data.csv is not defined in scenario_index_list.csv." << '\n';
						//has not been defined
						continue;

					}

					int activate_flag = 0; 
					parser_measurement.GetValueByFieldName("activate", activate_flag);

					if (activate_flag == 0)
						continue; 

					// add the to node id into the outbound (adjacent) node list
					if (g_node_id_to_seq_no_map.find(from_node_id) == assignment.g_node_id_to_seq_no_map.end())
					{
						dtalog.output() << "[ERROR] from_node_id " << from_node_id << " in file sensor_data.csv is not defined in node.csv." << '\n';
						//has not been defined
						continue;
					}
					if (g_node_id_to_seq_no_map.find(to_node_id) == assignment.g_node_id_to_seq_no_map.end())
					{
						dtalog.output() << "[ERROR] to_node_id " << to_node_id << " in file sensor_data.csv is not defined in node.csv." << '\n';
						//has not been defined
						continue;
					}

					float count = -1;

					//demand_period
					//mode_type 
					string demand_period_str;
					parser_measurement.GetValueByFieldName("demand_period", demand_period_str);
					int tau = 0;

					if (assignment.demand_period_to_seqno_mapping.find(demand_period_str) != assignment.demand_period_to_seqno_mapping.end())
						tau = assignment.demand_period_to_seqno_mapping[demand_period_str];

					{
						int upper_bound_flag = 0;
						{
							int demand_period_id = assignment.g_DemandPeriodVector[tau].demand_period_id;

							if (assignment.g_DemandPeriodVector[tau].number_of_demand_files == 0)
								continue;
							parser_measurement.GetValueByFieldName("count", count, true);
							parser_measurement.GetValueByFieldName("upper_bound_flag", upper_bound_flag, true);

						}
						// map external node number to internal node seq no.
						int internal_from_node_seq_no = assignment.g_node_id_to_seq_no_map[from_node_id];
						int internal_to_node_seq_no = assignment.g_node_id_to_seq_no_map[to_node_id];

						if (g_node_vector[internal_from_node_seq_no].m_to_node_2_link_seq_no_map.find(internal_to_node_seq_no) != g_node_vector[internal_from_node_seq_no].m_to_node_2_link_seq_no_map.end())
						{
							int link_seq_no = g_node_vector[internal_from_node_seq_no].m_to_node_2_link_seq_no_map[internal_to_node_seq_no];
							if (g_link_vector[link_seq_no].VDF_period[tau].obs_count[scenario_index] >= 1)  // data exist
							{
								if (upper_bound_flag == 0)
								{  // over write only if the new data are acutal counts,

									g_link_vector[link_seq_no].VDF_period[tau].obs_count[scenario_index] = count;
									g_link_vector[link_seq_no].VDF_period[tau].upper_bound_flag[scenario_index] = upper_bound_flag;
									sensor_count++;
								}
								else  // if the new data are upper bound, skip it and keep the actual counts
								{
									// do nothing
								}


							}
							else
							{
								g_link_vector[link_seq_no].VDF_period[tau].obs_count[scenario_index] = count;
								g_link_vector[link_seq_no].VDF_period[tau].upper_bound_flag[scenario_index] = upper_bound_flag;
									sensor_count++;
							}
						}
						else
						{
							dtalog.output() << "[ERROR] Link " << from_node_id << "->" << to_node_id << " in file sensor_data.csv is not defined in link.csv." << '\n';
							continue;
						}
					}

		

			}

			parser_measurement.CloseCSVFile();
		}

		assignment.summary_file << "ODME stage: # of sensors =," << sensor_count << '\n';

		// step 1: input the measurements of
		// Pi
		// Dj
		// link l


		// step 2: loop for adjusting OD demand
		double prev_gap = 9999999;
		for (int s = 0; s < OD_updating_iterations; ++s)
		{
			double total_system_demand = 0;
			float total_gap = 0;
			float total_relative_gap = 0;
			float total_system_travel_cost = 0;
			//step 2.1
			// we can have a recursive formulat to reupdate the current link volume by a factor of k/(k+1),
			// and use the newly generated path flow to add the additional 1/(k+1)
			double system_gap = 0;

			double gap = g_reset_and_update_link_volume_based_on_ODME_columns(g_link_vector.size(), s, system_gap, false);
			//step 2.2: based on newly calculated path volumn, update volume based travel time, and update volume based measurement error/deviation
					// and use the newly generated path flow to add the additional 1/(k+1)

			double gap_increase = gap - prev_gap;

			if (s >= 5 && gap_increase > 0.01 )  // convergency criterion  // comment out to maintain consistency
			{
				dtalog.output() << "[PROCESS INFO] ODME stage terminates with gap increase = " << gap_increase * 100 << "% at iteration = " << s << '\n';
				assignment.summary_file << "ODME stage terminates with gap increase = " << gap_increase*100 << "% at iteration = " << s <<  '\n';
			    break;

			}


			if (gap < 0.01)  // convergency criterion  // comment out to maintain consistency
			{
				dtalog.output() << "[PROCESS INFO] ODME stage terminates with gap < 1% as " << gap * 100 << "% at iteration = " << s << '\n';
				assignment.summary_file << "ODME stage terminates with gap < 1% as " << gap * 100 << "% at iteration = " << s << '\n';
				break;

			}

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
				int from_zone_sindex = g_zone_vector[orig].sindex;
				if (from_zone_sindex == -1)
					continue;

				for (int dest = 0; dest < g_zone_vector.size(); ++dest) //d
				{

					int to_zone_sindex = g_zone_vector[dest].sindex;
					if (to_zone_sindex == -1)
						continue;

					for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)  //at agent type
					{
						for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); ++tau)  //tau, assginment period
						{
							p_column_pool = &(assignment.g_column_pool[from_zone_sindex][to_zone_sindex][at][tau]);
							if (p_column_pool->od_volume[assignment.active_scenario_index] > 0)
							{
								column_pool_counts++;

								column_vector_size = p_column_pool->path_node_sequence_map.size();
								path_seq_count = 0;

								it_begin = p_column_pool->path_node_sequence_map.begin();
								it_end = p_column_pool->path_node_sequence_map.end();

								//stage 1: least cost
								double least_cost = 999999;
								int least_cost_path_seq_no = -1;
								int least_cost_path_node_sum_index = -1;
								path_seq_count = 0;

								it_begin = p_column_pool->path_node_sequence_map.begin();
								it_end = p_column_pool->path_node_sequence_map.end();
								for (it = it_begin; it != it_end; ++it)
								{
									path_toll = 0;
									path_distance = 0;
									path_travel_time = 0;

									if (s == 0) // perform computing of path travel time for the first iteration
									{

										for (int nl = 0; nl < it->second.m_link_size; ++nl)  // arc a along the path
										{
											link_seq_no = it->second.path_link_vector[nl];
											path_toll += g_link_vector[link_seq_no].VDF_period[tau].toll[at];
											path_distance += g_link_vector[link_seq_no].link_distance_VDF;
											double link_travel_time = g_link_vector[link_seq_no].link_avg_travel_time_per_period[tau][at];
											path_travel_time += link_travel_time;

											if (g_link_vector[link_seq_no].VDF_period[tau].obs_count[assignment.active_scenario_index] >= 1)  // added with mustafa 12/24/2022, verified
											{
												it->second.measurement_flag = 1;  // this path column has measurement
												it->second.path_sensor_link_vector.push_back(link_seq_no);
											}
										}


									it->second.path_toll = path_toll;
									it->second.path_travel_time = path_travel_time;


									if (path_travel_time < least_cost)
									{
										least_cost = path_travel_time;
										least_cost_path_seq_no = it->second.path_seq_no;
										least_cost_path_node_sum_index = it->first;
									}
									}

								}


								//stage 2: deviation based on observation
								int i = 0;
								for (it = it_begin; it != it_end; ++it, ++i) // for each k
								{

									column_path_counts++;

									if (s >= 1 && it->second.measurement_flag == 0)  // after 1  iteration, if there are no data passing through this path column. we will skip it in the ODME process
									{
										total_system_demand += it->second.path_volume;
										continue;
									}



									it->second.UE_gap = (it->second.path_travel_time - least_cost);
									path_gradient_cost = 0;
									path_distance = 0;
									path_travel_time = 0;
									p_column_pool->m_passing_sensor_flag = -1;
									// step 3.1 origin production flow gradient

									// est_production_dev = est_production - obs_production;
									// requirement: when equality flag is 1,

								
									float est_count_dev = 0;
									for (int nl = 0; nl < it->second.path_sensor_link_vector.size(); ++nl)  // arc a  // modified with mustafa, 12/24/2022, verified
									{
										// step 3.3 link flow gradient
										link_seq_no = it->second.path_sensor_link_vector[nl];
										if (g_link_vector[link_seq_no].VDF_period[tau].obs_count[assignment.active_scenario_index] >= 1)
										{
											if (g_link_vector[link_seq_no].VDF_period[tau].upper_bound_flag[0] == 0)
											{
												path_gradient_cost += g_link_vector[link_seq_no].VDF_period[tau].est_count_dev[assignment.active_scenario_index];
												est_count_dev += g_link_vector[link_seq_no].VDF_period[tau].est_count_dev[assignment.active_scenario_index];
											}

											if (g_link_vector[link_seq_no].VDF_period[tau].upper_bound_flag[0] == 1 && g_link_vector[link_seq_no].VDF_period[tau].est_count_dev[assignment.active_scenario_index] > 0)
											{// we only consider the over capaity value here to penalize the path flow
												path_gradient_cost += g_link_vector[link_seq_no].VDF_period[tau].est_count_dev[assignment.active_scenario_index];
												est_count_dev += g_link_vector[link_seq_no].VDF_period[tau].est_count_dev[assignment.active_scenario_index];
											}
											p_column_pool->m_passing_sensor_flag += 1;
											it->second.measurement_flag = 1;
										}
									}

									// statistics collection

									if (it->second.measurement_flag >= 1)
										column_path_with_sensor_counts++;

									it->second.path_gradient_cost = path_gradient_cost;

									step_size = 1.0 / (OD_updating_iterations + 2.0);   // memo: with alicia, 0.05
									// memo: with Peiheng, use 1.0/(OD_updating_iterations+2.0) for stability
									double prev_path_volume = it->second.path_volume;

									double weight_of_measurements = 1;  // ad hoc weight on the measurements with respect to the UE gap// because unit of UE gap  is around 1-5 mins, measurement error is around 100 vehicles per hour per lane

									double change = step_size * (weight_of_measurements * it->second.path_gradient_cost + (1 - weight_of_measurements) * it->second.UE_gap);

									//dtalog.output() <<" path =" << i << ", gradient cost of measurements =" << it->second.path_gradient_cost << ", UE gap=" << it->second.UE_gap << '\n';

									float bound = 0.1;
									float change_lower_bound = it->second.path_volume * bound * (-1);
									float change_upper_bound = it->second.path_volume * bound;

									// reset
									if (change < change_lower_bound)
										change = change_lower_bound;

									// reset
									if (change > change_upper_bound)
										change = change_upper_bound;

									it->second.path_volume = max(0.000, it->second.path_volume - change);

									total_system_demand += it->second.path_volume;

									if (dtalog.log_odme() == 1)
									{
										dtalog.output() << "[DATA INFO] OD " << orig << "-> " << dest << " path id:" << i << ", prev_vol"
											<< prev_path_volume << ", gradient_cost = " << it->second.path_gradient_cost
											<< " UE gap," << it->second.UE_gap
											<< " link," << g_link_vector[link_seq_no].VDF_period[tau].est_count_dev
											<< "proposed change = " << step_size * it->second.path_gradient_cost
											<< "actual change = " << change
											<< "new vol = " << it->second.path_volume << '\n';
									}
								}  // end of loop for all paths in the column pools


							 //// record adjustment results
								//for (it = it_begin; it != it_end; ++it) // for each k
								//{
								//	it->second.path_time_per_iteration_ODME_map[s] = path_travel_time;
								//	it->second.path_volume_per_iteration_ODME_map[s] = it->second.path_volume;
								//}



								if (p_column_pool->m_passing_sensor_flag >= 1)
								    column_pool_with_sensor_counts++;

							}
						}
					}
				}
			}
			if (s == 0)
			{
				float percentage_of_OD_columns_with_sensors = column_pool_with_sensor_counts * 1.0 / max(1, column_pool_counts) * 100;
				float percentage_of_paths_with_sensors = column_path_with_sensor_counts * 1.0 / max(1, column_path_counts) * 100;
				dtalog.output() << "[DATA INFO] count of all column pool vectors=" << column_pool_counts << ", "
					<< "count of all paths =" << column_path_counts << ", "
					<< "count of column_pools with sensors = " << column_pool_with_sensor_counts << "(" << percentage_of_OD_columns_with_sensors << "%), "
					<< "count of column_paths with sensors = " << column_path_with_sensor_counts << " (" << percentage_of_paths_with_sensors << "%)" << '\n';

			}

			dtalog.output() << "[DATA INFO] total_system_demand =" << total_system_demand << ", ";

		}


		// post-procese link volume based on OD volumns
		// very import: noted by Peiheng and Xuesong on 01/30/2022
		double system_gap = 0;
		g_reset_and_update_link_volume_based_on_ODME_columns(g_link_vector.size(), OD_updating_iterations, system_gap, true);
		// we now have a consistent link-to-path volumne in g_link_vector[link_seq_no].total_volume_for_all_mode_types_per_period[tau]
	}


}


