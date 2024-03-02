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

void g_reset_link_district_performance_per_scenario(Assignment& assignment)
{

	for (int i = 0; i < g_link_vector.size(); ++i)
	{

		// virtual connectors
		if (g_link_vector[i].link_type_si[0] == -1)
			continue;
		for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
		{


			g_link_vector[i].VDF_period[tau].volume_before_dtm = 0;
			g_link_vector[i].VDF_period[tau].volume_after_dtm = 0;

			g_link_vector[i].VDF_period[tau].speed_before_dtm = 0;
			g_link_vector[i].VDF_period[tau].DoC_before_dtm = 0;
			g_link_vector[i].VDF_period[tau].P_before_dtm = 0;
			g_link_vector[i].VDF_period[tau].speed_after_dtm = 0;
			g_link_vector[i].VDF_period[tau].DoC_after_dtm = 0;
			g_link_vector[i].VDF_period[tau].P_after_dtm = 0;


		}

	}


}

void g_record_link_district_performance_per_scenario(Assignment& assignment, int base_case_flag)
{


	int total_volume = 0; 

	for (int i = 0; i < g_link_vector.size(); ++i)
	{

		// virtual connectors
		if (g_link_vector[i].link_type_si[0] == -1)
			continue;
		int mode_type_size = assignment.g_ModeTypeVector.size();
		for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
		{

			if (assignment.g_DemandPeriodVector[tau].number_of_demand_files == 0)
				continue;

			float speed = g_link_vector[i].free_speed;  // default speed
			if (g_link_vector[i].VDF_period[tau].avg_travel_time_0 > 0.001f)
				speed = g_link_vector[i].link_distance_VDF / max(0.000001, g_link_vector[i].VDF_period[tau].avg_travel_time_0) *60;

			float speed_ratio = speed / max(0.001, g_link_vector[i].free_speed);  // default speed

			for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
			{

			float mode_type_volume = g_link_vector[i].volume_per_mode_type_per_period[tau][at] + g_link_vector[i].VDF_period[tau].preload;
			//VMT,VHT,PMT,PHT,PDT
			float preload = g_link_vector[i].VDF_period[tau].preload;
			float VMT = mode_type_volume * g_link_vector[i].link_distance_VDF;

			float VHT = mode_type_volume * g_link_vector[i].VDF_period[tau].avg_travel_time_0 / 60.0;
			float PMT = mode_type_volume * g_link_vector[i].link_distance_VDF;
			float PHT = mode_type_volume * g_link_vector[i].VDF_period[tau].avg_travel_time_0 / 60.0;
			float PDT_vf = mode_type_volume * (g_link_vector[i].VDF_period[tau].avg_travel_time_0 - g_link_vector[i].VDF_period[tau].FFTT_at[0]) / 60.0;
			float PDT_vc = max(0.0, mode_type_volume * (g_link_vector[i].VDF_period[tau].avg_travel_time_0 - g_link_vector[i].VDF_period[tau].FFTT_at[0] * g_link_vector[i].free_speed / max(0.001f, g_link_vector[i].v_congestion_cutoff)) / 60.0);

			int scenario_no = assignment.g_active_DTAscenario_map[ assignment.active_scenario_index];
				g_link_vector[i].recorded_volume_per_period_per_at[tau][at][scenario_no] = g_link_vector[i].volume_per_mode_type_per_period[tau][at];
				total_volume += g_link_vector[i].recorded_volume_per_period_per_at[tau][at][scenario_no];

				g_link_vector[i].recorded_MEU_per_period_per_at[tau][at][scenario_no] = g_link_vector[i].converted_MEU_volume_per_period_per_at[tau][at];
				g_link_vector[i].recorded_capacity_per_period_per_at[tau][at][scenario_no] =
					g_link_vector[i].VDF_period[tau].capacity_at[at] * g_link_vector[i].number_of_lanes_si[assignment.active_scenario_index];
				g_link_vector[i].recorded_DOC_per_period_per_at[tau][at][scenario_no] = g_link_vector[i].VDF_period[tau].DOC_mode[at];
				g_link_vector[i].recorded_TT_per_period_per_at[tau][at][scenario_no] = g_link_vector[i].link_avg_travel_time_per_period[tau][at];

				g_link_vector[i].recorded_CO2_per_period_per_at[tau][at][scenario_no] = g_link_vector[i].link_avg_co2_emit_per_mode[tau][at]* mode_type_volume;
				g_link_vector[i].recorded_NOX_per_period_per_at[tau][at][scenario_no] = g_link_vector[i].link_avg_nox_emit_per_mode[tau][at]*mode_type_volume;

				if (base_case_flag == 0)
				{  // base case
					g_link_vector[i].VDF_period[tau].volume_before_dtm += mode_type_volume;
				}
				else
				{
					g_link_vector[i].VDF_period[tau].volume_after_dtm += mode_type_volume;
				}

				if (at == 0)
				{
					g_link_vector[i].VDF_period[tau].speed_before_dtm = speed;
					g_link_vector[i].VDF_period[tau].DoC_before_dtm = g_link_vector[i].VDF_period[tau].DOC_mode[0];
					g_link_vector[i].VDF_period[tau].P_before_dtm = g_link_vector[i].VDF_period[tau].P;
				}
				else  /// SA case
				{
					g_link_vector[i].VDF_period[tau].speed_after_dtm = speed;
					g_link_vector[i].VDF_period[tau].DoC_after_dtm = g_link_vector[i].VDF_period[tau].DOC_mode[0];
					g_link_vector[i].VDF_period[tau].P_after_dtm = g_link_vector[i].VDF_period[tau].P;
				}

			}


		}
		}

	//dtalog.output() << "[DATA INFO] Record link based and district_based performance for scenario_index =" << assignment.active_scenario_index << " with a total volume = " << total_volume <<  '\n';
	//g_DTA_log_file << "[DATA INFO] Record link based and district_based performance for scenario_index =" << assignment.active_scenario_index << " with a total volume = " << total_volume <<  '\n';

}

void g_outputZonalHierarchyMapping(Assignment& assignment)
{

	for (int district = 0; district < assignment.g_number_of_analysis_districts; district++)
	{
		std::vector<DTAGDPoint> points, points_in_polygon;

		//DTAGDPoint pt;
		//pt.x = 0; pt.y = 4; points.push_back(pt);
		//pt.x = 1; pt.y = 1; points.push_back(pt);
		//pt.x = 2; pt.y = 2; points.push_back(pt);
		//pt.x = 4; pt.y = 4; points.push_back(pt);
		//pt.x = 3; pt.y = 1; points.push_back(pt);
		//pt.x = 3; pt.y = 3; points.push_back(pt);
		//g_find_convex_hull(points, points_in_polygon);

		//test

		//clustering
		for (int orig = 0; orig < g_zone_vector.size(); orig++)  // o
		{
			if (g_zone_vector[orig].analysis_district_index == district)
			{
				DTAGDPoint pt;
				pt.x = g_zone_vector[orig].cell_x;
				pt.y = g_zone_vector[orig].cell_y;
				g_district_summary_map[district].district_zone_vector.push_back(g_zone_vector[orig].zone_id);
				points.push_back(pt);
			}
		}


		if (points.size() >= 4)  // 4 as restrictive conditions
		{
			g_find_convex_hull(points, points_in_polygon);
		}
		else
		{
			for (int i = 0; i < points.size(); i++)
			{
				points_in_polygon.push_back(points[i]);

			}
		}

		for (int i = 0; i < points_in_polygon.size(); i++)
		{
			g_district_summary_map[district].shape_points.push_back(points_in_polygon[i]);

		}


	}

	dtalog.output() << "[PROCESS INFO] Creating a file named 'zonal_hierarchy_mapping.csv', which contains data about zones related to each subarea. This file will help us focus our approach on specific areas. " << '\n';
	g_DTA_log_file << "[PROCESS INFO] Creating a file named 'zonal_hierarchy_mapping.csv', which contains data about zones related to each subarea. This file will help us focus our approach on specific areas. " << '\n';
	FILE* g_pFileZone = nullptr;
	g_pFileZone = fopen("zonal_hierarchy_mapping.csv", "w");

	if (g_pFileZone == NULL)
	{
		dtalog.output() << "[ERROR] File zonal_hierarchy_mapping.csv cannot be opened." << '\n';
		g_DTA_log_file << "[ERROR] File zonal_hierarchy_mapping.csv cannot be opened." << '\n';
		return; 
	}

	fprintf(g_pFileZone, "first_column,zone_id,x_coord,y_coord,super_zone_id,analysis_district_id,demand,subarea_significance_flag,inside_flag\n");
	for (int orig = 0; orig < g_zone_vector.size(); orig++)  // o
	{

		int analysis_district_id = assignment.g_zone_seq_no_to_analysis_distrct_id_mapping[g_zone_vector[orig].zone_seq_no];

		fprintf(g_pFileZone, "0,%d,%f,%f,%d,%d,%f,%d,%d\n",
			g_zone_vector[orig].zone_id, g_zone_vector[orig].cell_x, g_zone_vector[orig].cell_y,
			g_zone_vector[orig].sindex,
			analysis_district_id,
			g_zone_vector[orig].origin_zone_impact_volume,
			g_zone_vector[orig].subarea_significance_flag,
			g_zone_vector[orig].subarea_inside_flag);
	}

	fclose(g_pFileZone);
}

void g_OutputSummaryFiles(Assignment& assignment)
{


	assignment.summary_system_file.open("system_performance_summary.csv", std::fstream::out);

	assignment.summary_system_file << "first_column,scenario_index,scenario_name,demand_period,mode_type,od_volume,number_of_routes,total_distance_km,total_distance_mile,total_travel_time_min,total_travel_delay_min,total_co2_kg,total_nox_kg,avg_distance_km,avg_distance_mile,avg_travel_time_in_min,avg_travel_delay_in_min,avg_speed_kmph,avg_speed_mph,avg_co2_kg,avg_nox_kg," << '\n';
	int mode_type_size = assignment.g_ModeTypeVector.size();
	std::map<int, CSystem_Summary>::iterator it_s;
	for (it_s = g_scenario_summary_map.begin(); it_s != g_scenario_summary_map.end(); ++it_s)
	{
		for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
		{
			for (int at = 0; at < mode_type_size; ++at)
			{
				int scenario_no = assignment.g_active_DTAscenario_map[it_s->first];
				assignment.summary_system_file << "0," << it_s->first << "," << assignment.g_DTA_scenario_vector[scenario_no].scenario_name.c_str() << ",";
				assignment.summary_system_file << assignment.g_DemandPeriodVector[tau].demand_period.c_str() << ",";

				assignment.summary_system_file << assignment.g_ModeTypeVector[at].mode_type.c_str() << ",";
				assignment.summary_system_file << it_s->second.data_by_demand_period_mode_type[tau][at].total_od_volume << ",";

				it_s->second.computer_avg_value(tau, at);

				assignment.summary_system_file <<
					it_s->second.data_by_demand_period_mode_type[tau][at].count << "," <<
					it_s->second.data_by_demand_period_mode_type[tau][at].total_agent_distance_km << "," <<
					it_s->second.data_by_demand_period_mode_type[tau][at].total_agent_distance_mile << "," <<
					it_s->second.data_by_demand_period_mode_type[tau][at].total_agent_travel_time << "," <<
					it_s->second.data_by_demand_period_mode_type[tau][at].total_agent_delay << "," <<
					it_s->second.data_by_demand_period_mode_type[tau][at].total_agent_co2 << "," <<
					it_s->second.data_by_demand_period_mode_type[tau][at].total_agent_nox << "," <<

					it_s->second.data_by_demand_period_mode_type[tau][at].avg_travel_distance_km << "," <<
					it_s->second.data_by_demand_period_mode_type[tau][at].avg_travel_distance_mile << "," <<
					it_s->second.data_by_demand_period_mode_type[tau][at].avg_travel_time << "," << 
					it_s->second.data_by_demand_period_mode_type[tau][at].avg_travel_delay << "," <<
					it_s->second.data_by_demand_period_mode_type[tau][at].avg_travel_distance_km / max(0.0001, it_s->second.data_by_demand_period_mode_type[tau][at].avg_travel_time / 60.0) << "," <<
					it_s->second.data_by_demand_period_mode_type[tau][at].avg_travel_distance_mile / max(0.0001, it_s->second.data_by_demand_period_mode_type[tau][at].avg_travel_time / 60.0) << "," << it_s->second.data_by_demand_period_mode_type[tau][at].avg_co2 << "," <<
					it_s->second.data_by_demand_period_mode_type[tau][at].avg_nox << "," << '\n';

			}
		}
		// end of scenario
	}
	assignment.summary_system_file.close();


}

// FILE* g_pFileOutputLog = nullptr;

void g_output_dynamic_queue_profile()  // generated from VDF, from numerical queue evolution calculation
{
	dtalog.output() << "[STATUS INFO] writing link_queue_profile.csv.." << '\n';
	g_DTA_log_file << "[STATUS INFO] writing link_queue_profile.csv.." << '\n';

	int b_debug_detail_flag = 0;
	FILE* g_pFileLinkMOE = nullptr;

	string file_name = "link_queue_profile.csv";

	fopen_ss(&g_pFileLinkMOE, file_name.c_str(), "w");

	if (!g_pFileLinkMOE)
	{
		dtalog.output() << "[ERROR] File " << file_name.c_str() << " cannot be opened." << '\n';
		g_DTA_log_file << "[ERROR] File " << file_name.c_str() << " cannot be opened." << '\n';
		return; 
	}
	else
	{

		// Option 2: BPR-X function
		fprintf(g_pFileLinkMOE, "link_id,tmc_corridor_name,tmc_road_sequence,tmc,link_type_name,from_node_id,to_node_id,geometry,");

		fprintf(g_pFileLinkMOE, "link_type_code,FT,AT,nlanes,link_distance_meter,free_speed_kmph,capacity,k_critical,v_congestion_cutoff,");
		for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); tau++)
		{
			fprintf(g_pFileLinkMOE, "%s_Volume,%s_speed_BPR,%s_speed_QVDF, %s_t0,%s_t3,%s_P, %s_D,%s_DC_ratio,%s_mu,%s_gamma,",
				assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
				assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
				assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
				assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
				assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
				assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
				assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
				assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
				assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
				assignment.g_DemandPeriodVector[tau].demand_period.c_str()

			);
		}

		fprintf(g_pFileLinkMOE, "");

		// hourly data
		for (int t = 6 * 60; t < 20 * 60; t += 60)
		{
			int hour = t / 60;
			int minute = t - hour * 60;

			fprintf(g_pFileLinkMOE, "vh%02d,", hour);
		}

		for (int t = 6 * 60; t < 20 * 60; t += 60)
		{
			int hour = t / 60;
			int minute = t - hour * 60;

			fprintf(g_pFileLinkMOE, "vrh%02d,", hour);
		}

		for (int t = 6 * 60; t < 20 * 60; t += 15)
		{
			int hour = t / 60;
			int minute = t - hour * 60;
			fprintf(g_pFileLinkMOE, "v%02d:%02d,", hour, minute);
		}

		for (int t = 6 * 60; t < 20 * 60; t += 15)
		{
			int hour = t / 60;
			int minute = t - hour * 60;
			fprintf(g_pFileLinkMOE, "vr%02d:%02d,", hour, minute);
		}

		//for (int t = 6 * 60; t < 20 * 60; t += 5)
		//{
		//    int hour = t / 60;
		//    int minute = t - hour * 60;
		//    fprintf(g_pFileLinkMOE, "v%02d:%02d,", hour, minute);
		//}


		fprintf(g_pFileLinkMOE, "\n");

		//Initialization for all nodes
		for (int i = 0; i < g_link_vector.size(); ++i)
		{
			// virtual connectors
			if (g_link_vector[i].link_type_si[0] == -1)
				continue;

			double vehicle_volume_0 = 0;
			for (int tau = 0; tau < min((size_t)3, assignment.g_DemandPeriodVector.size()); tau++)
			{

				vehicle_volume_0 += g_link_vector[i].total_volume_for_all_mode_types_per_period[tau] + g_link_vector[i].VDF_period[tau].preload + g_link_vector[i].VDF_period[tau].sa_volume;
			}

			if (vehicle_volume_0 < 1)  // no volume links, skip
				continue;


			fprintf(g_pFileLinkMOE, "%s,%s,%d,%s,%s,%d,%d,\"%s\",",
				g_link_vector[i].link_id.c_str(),
				g_link_vector[i].tmc_corridor_name.c_str(),
				g_link_vector[i].tmc_road_sequence,
				g_link_vector[i].tmc_code.c_str(),

				g_link_vector[i].link_type_name.c_str(),

				g_node_vector[g_link_vector[i].from_node_seq_no].node_id,
				g_node_vector[g_link_vector[i].to_node_seq_no].node_id,
				g_link_vector[i].geometry.c_str());


			fprintf(g_pFileLinkMOE, "%s,%d,%d,%d,%f,%f,%f,%f,%f,",
				g_link_vector[i].link_type_code.c_str(),
				g_link_vector[i].FT,
				g_link_vector[i].AT,
				g_link_vector[i].number_of_lanes_si[0],
				g_link_vector[i].link_distance_VDF,
				g_link_vector[i].free_speed,
				g_link_vector[i].lane_capacity,
				g_link_vector[i].k_critical,
				g_link_vector[i].v_congestion_cutoff);

			//AM g_link_vector[i].VDF_period[0].

			int period_index = 0;

			double assignment_PMT;
			double assignment_PHT;
			double assignment_PSDT;
			double assignment_VCDT = 0;

			double assignment_VMT;
			double assignment_VHT;

			for (int tau = 0; tau < min((size_t)3, assignment.g_DemandPeriodVector.size()); tau++)
			{
				double vehicle_volume = g_link_vector[i].total_volume_for_all_mode_types_per_period[period_index] + g_link_vector[i].VDF_period[period_index].preload + g_link_vector[i].VDF_period[period_index].sa_volume;
				double mode_type_volume = g_link_vector[i].total_agent_volume_for_all_mode_types_per_period[period_index] + g_link_vector[i].VDF_period[period_index].preload + g_link_vector[i].VDF_period[period_index].sa_volume;

				assignment_VMT = g_link_vector[i].link_distance_VDF * vehicle_volume;
				assignment_VHT = g_link_vector[i].link_avg_travel_time_per_period[period_index][0] * vehicle_volume / 60.0;  // 60.0 converts min to hour

				assignment_PMT = g_link_vector[i].link_distance_VDF * mode_type_volume;
				assignment_PHT = g_link_vector[i].link_avg_travel_time_per_period[period_index][0] * mode_type_volume / 60.0;  // 60.0 converts min to hour

				assignment_PSDT = (g_link_vector[i].link_avg_travel_time_per_period[period_index][0] - g_link_vector[i].free_flow_travel_time_in_min) * mode_type_volume / 60.0;  // 60.0 converts min to hour

				double VCTT = g_link_vector[i].link_distance_VDF / max(1.0f, g_link_vector[i].v_congestion_cutoff) * 60;
				assignment_VCDT = max(0.0, g_link_vector[i].link_avg_travel_time_per_period[period_index][0] - VCTT) * g_link_vector[i].total_volume_for_all_mode_types_per_period[period_index] / 60.0;  // 60.0 converts min to hour


				fprintf(g_pFileLinkMOE, "%f,%f,%f, %f,%f,%f, %f,%f,%f,%f,",
					g_link_vector[i].total_volume_for_all_mode_types_per_period[period_index] + g_link_vector[i].VDF_period[period_index].preload + g_link_vector[i].VDF_period[period_index].sa_volume,
					g_link_vector[i].VDF_period[period_index].avg_speed_BPR,
					g_link_vector[i].VDF_period[period_index].avg_queue_speed,

					g_link_vector[i].VDF_period[period_index].t0,
					g_link_vector[i].VDF_period[period_index].t3,
					g_link_vector[i].VDF_period[period_index].P,

					g_link_vector[i].VDF_period[period_index].lane_based_D,
					g_link_vector[i].VDF_period[period_index].DOC_mode[0],
					g_link_vector[i].VDF_period[period_index].Q_mu,
					g_link_vector[i].VDF_period[period_index].Q_gamma
				);
			}

			for (int t = 6 * 60; t < 20 * 60; t += 60)
			{
				float speed = g_link_vector[i].get_model_hourly_speed(t);
				fprintf(g_pFileLinkMOE, "%.3f,", speed);
			}

			for (int t = 6 * 60; t < 20 * 60; t += 60)
			{
				float speed_ratio = g_link_vector[i].get_model_hourly_speed(t) / max(1.0f, g_link_vector[i].v_congestion_cutoff);
				if (speed_ratio > 1)
					speed_ratio = 1;

				fprintf(g_pFileLinkMOE, "%.3f,", speed_ratio);
			}


			for (int t = 6 * 60; t < 20 * 60; t += 15)
			{
				int time_interval = t / 5;
				float speed = g_link_vector[i].model_speed[time_interval];
				fprintf(g_pFileLinkMOE, "%.3f,", speed);
			}

			for (int t = 6 * 60; t < 20 * 60; t += 15)
			{
				int time_interval = t / 5;
				float speed_ratio = g_link_vector[i].model_speed[time_interval] / max(1.0f, g_link_vector[i].v_congestion_cutoff);
				if (speed_ratio > 1)
					speed_ratio = 1;

				fprintf(g_pFileLinkMOE, "%.3f,", speed_ratio);
			}

			fprintf(g_pFileLinkMOE, "\n");

		}  // for each link l
		fclose(g_pFileLinkMOE);
	}//assignment mode 2 as simulation

}

void g_output_district_performance_result(Assignment& assignment)
{


	char district_performance_file_name[50];
	int scenario_no = assignment.g_active_DTAscenario_map[assignment.active_scenario_index];
	sprintf(district_performance_file_name, "district_performance_s%d_%s.csv", assignment.active_scenario_index, assignment.g_DTA_scenario_vector[scenario_no].scenario_name.c_str());
	
	
	dtalog.output() << "[STATUS INFO] writing district_performance_s.csv.." << '\n';
	g_DTA_log_file << "[STATUS INFO] writing district_performance_s.csv.." << '\n';


	assignment.summary_system_file.open(district_performance_file_name, std::fstream::out);
	assignment.summary_system_file << "first_column,district_seq_id,scenario_name,demand_period,mode_type,od_volume,total_distance_km,total_distance_mile,total_travel_time_min,avg_distance_km,avg_distance_mile,avg_travel_time_in_min,zone_vector,geometry" << '\n'; 

	int mode_type_size = assignment.g_ModeTypeVector.size();
	std::map<int, CSystem_Summary>::iterator it_s;
	for (it_s = g_district_summary_map.begin(); it_s != g_district_summary_map.end(); ++it_s)
	{
		for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
		{
			for (int at = 0; at < mode_type_size; ++at)
			{

				assignment.summary_system_file << "0," << it_s->first << "," << assignment.g_DTA_scenario_vector[scenario_no].scenario_name.c_str() << ",";
				assignment.summary_system_file << assignment.g_DemandPeriodVector[tau].demand_period.c_str() << ",";

				assignment.summary_system_file << assignment.g_ModeTypeVector[at].mode_type.c_str() << ",";
				assignment.summary_system_file << it_s->second.data_by_demand_period_mode_type[tau][at].total_od_volume << ",";

				it_s->second.computer_avg_value(tau, at);

				assignment.summary_system_file <<
					it_s->second.data_by_demand_period_mode_type[tau][at].total_agent_distance_km << "," <<
					it_s->second.data_by_demand_period_mode_type[tau][at].total_agent_distance_mile << "," <<
					it_s->second.data_by_demand_period_mode_type[tau][at].total_agent_travel_time << "," <<
					it_s->second.data_by_demand_period_mode_type[tau][at].avg_travel_distance_km << "," <<
					it_s->second.data_by_demand_period_mode_type[tau][at].avg_travel_distance_mile << "," <<
					it_s->second.data_by_demand_period_mode_type[tau][at].avg_travel_time << ",";

				for (int i = 0; i < it_s->second.district_zone_vector.size(); i++)
				{
					assignment.summary_system_file << it_s->second.district_zone_vector[i] << ";";
						
				}
				assignment.summary_system_file << ",";

				assignment.summary_system_file << "\"POLYGON ((";

				for (int i = 0; i < it_s->second.shape_points.size(); i++)
				{
					assignment.summary_system_file << it_s->second.shape_points[i].x << " " << it_s->second.shape_points[i].y << ",";

				}
				assignment.summary_system_file << "))\"";
				assignment.summary_system_file << "\n";

				it_s->second.reset_data(tau, at);  // reset the statistics 
			}
		}
		// end of scenario

	}
	assignment.summary_system_file.close();

}

void g_output_route_assignment_results(Assignment& assignment, int subarea_id)
{

	double path_time_vector[MAX_LINK_SIZE_IN_A_PATH];
	char route_assignment_file_name[50];

	int scenario_no = assignment.g_active_DTAscenario_map[assignment.active_scenario_index];
	sprintf(route_assignment_file_name, "route_assignment_s%d_%s.csv", assignment.active_scenario_index, assignment.g_DTA_scenario_vector[scenario_no].scenario_name.c_str());



	FILE* g_pFilePathMOE = nullptr;
	fopen_ss(&g_pFilePathMOE, route_assignment_file_name, "w");
	dtalog.output() << "[STATUS INFO] writing route_assignment.csv for each scenario" << '\n';
	g_DTA_log_file << "[STATUS INFO] writing route_assignment.csv for each scenario" << '\n';

	if (!g_pFilePathMOE)
	{
		dtalog.output() << "[ERROR] File " << route_assignment_file_name << " cannot be opened." << '\n';
		g_DTA_log_file << "[ERROR] File " << route_assignment_file_name << " cannot be opened." << '\n';
		return;
	}

	fprintf(g_pFilePathMOE, "first_column,route_seq_id,o_zone_id,d_zone_id,o_super_zone_index,d_super_zone_index,od_pair_key,information_type,mode_type,demand_period,volume,");
	fprintf(g_pFilePathMOE, "distance_km,distance_mile,travel_time,FFTT,co2,nox,UE_OD_based_relative_gap_compared_to_least_time_path_travel_time,UE_path_based_gap,preload_volume_from_route_input_file,ODME_volume_before,ODME_volume_after,ODME_volume_diff,SIMU_volume,SIMU_travel_time,");

	fprintf(g_pFilePathMOE, "toll,number_of_nodes,node_sequence,link_id_sequence,");
	//, ODME_#_of_sensor_links,
	// subarea_flag,
	fprintf(g_pFilePathMOE, "geometry,");

	fprintf(g_pFilePathMOE, "link_special_flag_sequence,sequential_link_delay,sequential_link_FFTT,");

	fprintf(g_pFilePathMOE, "DTM_OD_impact,DTM_path_impact,DTM_#_of_lane_closure_links,lane_closure_link_vector,DTM_new_path_generated,DTM_volume_before,DTM_volume_after,DTM_volume_diff, ");



	//// stage 1: column updating
	//for (int iteration_number = 0; iteration_number < min(20, assignment.g_number_of_column_updating_iterations); iteration_number++)
	//{
	//    fprintf(g_pFilePathMOE, "TT_%d,", iteration_number);
	//}

	//for (int iteration_number = 0; iteration_number < min(20, assignment.g_number_of_column_updating_iterations); iteration_number++)
	//{
	//    fprintf(g_pFilePathMOE, "VOL_%d,", iteration_number);
	//}

	//stage 2: ODME
	//for (int iteration_number = 0; iteration_number < min(20, assignment.g_number_of_ODME_iterations); iteration_number++)
	//{
	//	fprintf(g_pFilePathMOE, "ODME_TT_%d,", iteration_number);
	//}

	//for (int iteration_number = max(0, assignment.g_number_of_ODME_iterations - 10); iteration_number < assignment.g_number_of_ODME_iterations; iteration_number++)
	//{
	//	fprintf(g_pFilePathMOE, "ODME_Vol_%d,", iteration_number);
	//}

	//stage 3: sensitivity analysis
	for (int iteration_number = 0; iteration_number <= assignment.g_number_of_sensitivity_analysis_iterations_for_dtm; iteration_number++)
	{
		fprintf(g_pFilePathMOE, "DTM_TT_%d,", iteration_number);
	}
	for (int iteration_number = 0; iteration_number <= assignment.g_number_of_sensitivity_analysis_iterations_for_dtm; iteration_number++)
	{
		if (iteration_number == 0)
			fprintf(g_pFilePathMOE, "DTM_Vol_%d,", iteration_number);
		else
			fprintf(g_pFilePathMOE, "DTM_Delta_Vol_%d,", iteration_number);
	}
	fprintf(g_pFilePathMOE, "\n");

	int count = 1;

	clock_t start_t, end_t;
	start_t = clock();
	clock_t iteration_t;


	int mode_type_size = assignment.g_ModeTypeVector.size();
	int zone_size = g_zone_vector.size();
	int demand_period_size = assignment.g_DemandPeriodVector.size();

	CColumnVector* p_column_pool;

	float path_toll = 0;
	double path_distance = 0;
	double path_FFTT = 0;
	double path_distance_km = 0;
	double path_distance_ml = 0;
	double path_travel_time = 0;
	double path_travel_delay = 0;



	float path_travel_time_without_access_link = 0;
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
			//for (int tau = 0; tau < demand_period_size; ++tau)
			//{
			//	// used in travel time calculation
			//	g_link_vector[i].background_total_volume_for_all_mode_types_per_period[tau] = 0;
			//}

			if (g_node_vector[g_link_vector[i].from_node_seq_no].subarea_id >= 1 && g_node_vector[g_link_vector[i].to_node_seq_no].node_id >= 1)
			{

				b_subarea_mode = true;
			}

			for (int tau = 0; tau < demand_period_size; ++tau)
			{
			g_link_vector[i].VDF_period[tau].turn_link_count_map.clear();


			}

		}


		/// <summary>  screening the path flow pattern
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
					if (to_zone_sindex == -1 || to_zone_sindex == from_zone_sindex)
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
								if (subarea_output_flag == 0)
								{
									for (int nl = 0; nl < it->second.m_link_size; ++nl)  // arc a
									{
										int link_seq_no = it->second.path_link_vector[nl];
										g_link_vector[link_seq_no].background_total_volume_for_all_mode_types_per_period[tau] += it->second.path_volume;
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

		///// output background_link_volume.csv
		//dtalog.output() << "[STATUS INFO] writing link_background_volume.csv.." << '\n';
		//g_DTA_log_file << "[STATUS INFO] writing link_background_volume.csv.." << '\n';

		//int b_debug_detail_flag = 0;
		//FILE* g_pFileLinkMOE = nullptr;
		//int b_background_link_volume_file = 0;

		//if (b_background_link_volume_file)
		//{
		//	fopen_ss(&g_pFileLinkMOE, "link_background_volume.csv", "w");
		//	if (!g_pFileLinkMOE)
		//	{
		//		dtalog.output() << "[ERROR] File link_background_volume.csv cannot be opened." << '\n';
		//		g_DTA_log_file << "[ERROR] File link_background_volume.csv cannot be opened." << '\n';
		//		g_program_stop();
		//	}
		//	else
		//	{
		//		fprintf(g_pFileLinkMOE, "link_id,from_node_id,to_node_id,subarea_id,time_period,volume,background_volume,major_path_volume,ratio_of_major_path_flow,geometry,");

		//		fprintf(g_pFileLinkMOE, "notes\n");

		//		//Initialization for all nodes
		//		for (int i = 0; i < g_link_vector.size(); ++i)
		//		{
		//			//   virtual connectors
		//			if (g_link_vector[i].link_type_si[0] == -1)
		//				continue;

		//			for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
		//			{
		//				double volume = g_link_vector[i].total_volume_for_all_mode_types_per_period[tau] + g_link_vector[i].VDF_period[tau].preload;
		//				double major_path_link_volume = g_link_vector[i].total_volume_for_all_mode_types_per_period[tau] + g_link_vector[i].VDF_period[tau].preload - g_link_vector[i].background_total_volume_for_all_mode_types_per_period[tau];
		//				double ratio = major_path_link_volume / max(volume, 0.000001);

		//				if (volume < 0.0000001)
		//					ratio = -1;
		//				fprintf(g_pFileLinkMOE, "%s,%d,%d,%d,%s,%.3f,%.3f,%.3f,%.3f,\"%s\",",
		//					g_link_vector[i].link_id.c_str(),
		//					g_node_vector[g_link_vector[i].from_node_seq_no].node_id,
		//					g_node_vector[g_link_vector[i].to_node_seq_no].node_id,
		//					g_link_vector[i].subarea_id,
		//					assignment.g_DemandPeriodVector[tau].time_period.c_str(),
		//					g_link_vector[i].total_volume_for_all_mode_types_per_period[tau] + g_link_vector[i].VDF_period[tau].preload,
		//					g_link_vector[i].background_total_volume_for_all_mode_types_per_period[tau],
		//					major_path_link_volume,
		//					ratio,
		//					g_link_vector[i].geometry.c_str());
		//				fprintf(g_pFileLinkMOE, "\n");

		//			}

		//		}

		//		fclose(g_pFileLinkMOE);
			//}

		//}


	} // end of path flow pattern screening

	for (int orig = 0; orig < zone_size; ++orig)
	{
		if (g_zone_vector[orig].zone_id % 100 == 0)
		{
			dtalog.output() << "[DATA INFO] o zone id =  " << g_zone_vector[orig].zone_id << '\n';
			g_DTA_log_file << "[DATA INFO] o zone id =  " << g_zone_vector[orig].zone_id << '\n';
		}
		int from_zone_sindex = g_zone_vector[orig].sindex;
		if (from_zone_sindex == -1)
			continue;

		if (g_zone_vector[orig].b_shortest_path_computing_flag == false)
			continue;

		for (int dest = 0; dest < zone_size; ++dest)
		{
			int to_zone_sindex = g_zone_vector[dest].sindex;
			if (to_zone_sindex == -1)
				continue;

			if (g_zone_vector[dest].b_shortest_path_computing_flag == false)
				continue;

			for (int tau = 0; tau < demand_period_size; ++tau)
			{

				int global_od_impact_flag_across_all_mode_types = 0;
				for (int at = 0; at < mode_type_size; ++at)
				{
					p_column_pool = &(assignment.g_column_pool[from_zone_sindex][to_zone_sindex][at][tau]);

					if (p_column_pool->OD_impact_flag == 1 && assignment.g_ModeTypeVector[at].real_time_information_type == 0)  // as long as regular users of this OD pair is impacted.
					{
						global_od_impact_flag_across_all_mode_types = 1;
					}
				}

				for (int at = 0; at < mode_type_size; ++at)
				{
					p_column_pool = &(assignment.g_column_pool[from_zone_sindex][to_zone_sindex][at][tau]);



					if (p_column_pool->od_volume[assignment.active_scenario_index] > 0)
					{
						double od_toll = 0;
						double od_trave_time = 0;
						double od_distance_km = 0;

						int at_OD_impact_flag;
						if (p_column_pool->at_od_impacted_flag_map.find(at) != p_column_pool->at_od_impacted_flag_map.end())
						{
							at_OD_impact_flag = 1;  // marked
						}
						else
						{
							at_OD_impact_flag = 0;
						}



						int information_type = p_column_pool->information_type;


						time_stamp = (assignment.g_DemandPeriodVector[tau].starting_time_slot_no + assignment.g_DemandPeriodVector[tau].ending_time_slot_no) / 2.0 * MIN_PER_TIMESLOT;

						// scan through the map with different node sum for different continuous paths
						it_begin = p_column_pool->path_node_sequence_map.begin();
						it_end = p_column_pool->path_node_sequence_map.end();

						for (it = it_begin; it != it_end; ++it)
						{
							if (at == 1)
							{
								int idebug = 1;
							}

							if (count % 100000 == 0)
							{
								end_t = clock();
								iteration_t = end_t - start_t;
								dtalog.output() << "[STATUS INFO] writing " << count / 1000 << "K paths with CPU time " << iteration_t / 1000.0 << " s" << '\n';
								g_DTA_log_file << "[STATUS INFO] writing " << count / 1000 << "K paths with CPU time " << iteration_t / 1000.0 << " s" << '\n';
							}

							path_toll = 0;
							path_FFTT = 0;
							path_distance = 0;
							path_distance_km = 0;
							path_distance_ml = 0;
							path_travel_time = 0;
							path_travel_delay = 0;
							path_travel_time_without_access_link = 0;
							path_FF_travel_time = 0;

							double path_co2 = 0;
							double path_nox = 0;

							path_time_vector[0] = time_stamp;
							path_time_vector[1] = time_stamp;

							string link_specifical_flag_str;


							if (it->second.m_link_size >= MAX_LINK_SIZE_IN_A_PATH - 2)
							{
								dtalog.output() << "[ERROR] it->second.m_link_size < MAX_LINK_SIZE_IN_A_PATH" << '\n';
								g_DTA_log_file << "[ERROR] it->second.m_link_size < MAX_LINK_SIZE_IN_A_PATH" << '\n';
								dtalog.output() << "[DATA INFO] o= " << g_zone_vector[orig].zone_id << '\n';
								g_DTA_log_file << "[DATA INFO] o= " << g_zone_vector[orig].zone_id << '\n';
								dtalog.output() << "[DATA INFO] d= " << g_zone_vector[dest].zone_id << '\n';
								g_DTA_log_file << "[DATA INFO] d= " << g_zone_vector[dest].zone_id << '\n';
								dtalog.output() << "[DATA INFO] agent type= " << assignment.g_ModeTypeVector[at].mode_type.c_str() << '\n';
								g_DTA_log_file << "[DATA INFO] agent type= " << assignment.g_ModeTypeVector[at].mode_type.c_str() << '\n';

								for (int nl = 0; nl < it->second.m_link_size; ++nl)  // arc a
								{
									dtalog.output() << "[DATA INFO] node no." << nl << " =" << g_node_vector[it->second.path_node_vector[nl]].node_id << '\n';
									g_DTA_log_file << "[DATA INFO] node no." << nl << " =" << g_node_vector[it->second.path_node_vector[nl]].node_id << '\n';
								}

								continue; 
							}

							for (int nl = 0; nl < it->second.m_link_size; ++nl)  // arc a
							{

								int link_seq_no = it->second.path_link_vector[nl];
								if (g_link_vector[link_seq_no].link_type_si[assignment.active_scenario_index] >= 0)
								{
									path_toll += g_link_vector[link_seq_no].VDF_period[tau].toll[at][assignment.active_scenario_index];
									path_FFTT += g_link_vector[link_seq_no].VDF_period[tau].FFTT_at[at];
									path_distance += g_link_vector[link_seq_no].link_distance_VDF;
									path_distance_km += g_link_vector[link_seq_no].link_distance_km;
									path_distance_ml += g_link_vector[link_seq_no].link_distance_mile;
									double link_travel_time = g_link_vector[link_seq_no].link_avg_travel_time_per_period[tau][at];
									double link_CO2 = g_link_vector[link_seq_no].link_avg_co2_emit_per_mode[tau][at];
									double link_NOX = g_link_vector[link_seq_no].link_avg_nox_emit_per_mode[tau][at];

									path_travel_time += link_travel_time;
									path_travel_delay += max(0.0,link_travel_time - g_link_vector[link_seq_no].VDF_period[tau].FFTT_at[at]);
									path_co2 += link_CO2;
									path_nox += link_NOX;


									od_toll += g_link_vector[link_seq_no].VDF_period[tau].toll[at][assignment.active_scenario_index] * it->second.path_volume;
									od_trave_time += link_travel_time * it->second.path_volume;
									od_distance_km += g_link_vector[link_seq_no].link_distance_VDF * it->second.path_volume;

									if (g_link_vector[link_seq_no].link_type_si[assignment.active_scenario_index] != 1000)   // skip access link
									{
										path_travel_time_without_access_link += link_travel_time;
									}

									path_FF_travel_time += g_link_vector[link_seq_no].VDF_period[tau].FFTT_at[at];
									time_stamp += link_travel_time;
									path_time_vector[nl + 1] = time_stamp;
									link_specifical_flag_str += g_link_vector[link_seq_no].link_specifical_flag_str;
								}
								else
								{
									int virtual_link = 0;
								}
							}

							CModeType_Summary l_element;

							l_element.count = 1;
							l_element.total_agent_travel_time = path_travel_time * it->second.path_volume;
							l_element.total_agent_delay = path_travel_delay * it->second.path_volume;
							
							l_element.total_agent_distance_km = path_distance_km * it->second.path_volume;
							l_element.total_agent_distance_mile = path_distance_ml * it->second.path_volume;
							l_element.total_agent_co2 = path_co2 * it->second.path_volume;
							l_element.total_agent_nox = path_nox * it->second.path_volume;

							int analysis_district_id = assignment.g_zone_seq_no_to_analysis_distrct_id_mapping[orig];

							g_district_summary_map[analysis_district_id].record_mode_volume(tau, at, it->second.path_volume);  // based on analysis district id, we count the total volume 
							g_district_summary_map[analysis_district_id].record_mode_od_data(l_element, tau, at);  // based on analysis district id, we count the OD performance 

							g_scenario_summary_map[assignment.active_scenario_index].record_mode_od_data(l_element, tau, at);  // we count system wide performance. the total volume has been counted in demand reading part 

							double total_agent_path_travel_time = 0;

							for (int vi = 0; vi < it->second.agent_simu_id_vector.size(); ++vi)
							{
								int agent_simu_id = it->second.agent_simu_id_vector[vi];
								CAgent_Simu* pAgentSimu = g_agent_simu_vector[agent_simu_id];
								total_agent_path_travel_time += pAgentSimu->path_travel_time_in_min;
							}

							double final_simu_path_travel_time = path_travel_time;  // by default

							if (it->second.agent_simu_id_vector.size() > 1)  // with simulated agents
							{
								final_simu_path_travel_time = total_agent_path_travel_time / it->second.agent_simu_id_vector.size();
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

							if (information_type == 1)  // information diversion start from physical nodes
							{
								virtual_first_link_delta = 0;
								virtual_last_link_delta = 1;
							}

							if (p_column_pool->activity_zone_no_vector.size() > 0)  // with activity zones
							{
								virtual_first_link_delta = 0;
								virtual_last_link_delta = 0;
							}



							// assignment_mode = 1, path flow mode
							if (assignment.assignment_mode != lue)
							{
								if (it->second.m_node_size - virtual_first_link_delta - virtual_last_link_delta <= 1)
									continue;

								int impacted_path_flag = it->second.impacted_path_flag;  //NA by default


								double volume_before_ODME = max(0.0, it->second.path_volume_before_ODME);
								double volume_after_ODME = max(0.0, it->second.path_volume_after_ODME);
								double volume = it->second.path_volume;

								double volume_diff_ODME = 0;

								if (volume_before_ODME >= -0.000001)
								{
									volume_diff_ODME = volume_after_ODME - volume_before_ODME;
								}

								double volume_before_dtm = max(0.0, it->second.path_volume_before_dtm);
								double volume_after_dtm = max(0.0, it->second.path_volume_after_dtm);
								double volume_diff_dtm = 0;

								if (volume_before_dtm >= -0.000001)
								{
									volume_diff_dtm = volume_after_dtm - volume_before_dtm;
								}

								// keep this record
								it->second.route_seq_id= count;
								fprintf(g_pFilePathMOE, "0,%d,%d,%d,%d,%d,%d->%d,%d,%s,%s,",
									count,
									g_zone_vector[orig].zone_id,
									g_zone_vector[dest].zone_id,
									g_zone_vector[orig].sindex,
									g_zone_vector[dest].sindex,
									g_zone_vector[orig].zone_id,
									g_zone_vector[dest].zone_id,
									information_type,
									assignment.g_ModeTypeVector[at].mode_type.c_str(),
									assignment.g_DemandPeriodVector[tau].demand_period.c_str());

								fprintf(g_pFilePathMOE, "%.4f,%.4f,%.4f,%.4f,%.4f,%f,%f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%d,%.4f,%.4f,%d,",
									volume,
									path_distance_km,
									path_distance_ml,
									path_travel_time,
									path_FFTT,
									path_co2,
									path_nox,
									p_column_pool->OD_based_UE_relative_gap,
									it->second.path_gradient_cost_relative_difference,
									it->second.path_preload_volume,
									volume_before_ODME,
									volume_after_ODME,
									volume_diff_ODME,
									it->second.agent_simu_id_vector.size(),
									final_simu_path_travel_time,
									path_toll,
									it->second.m_node_size - virtual_first_link_delta - virtual_last_link_delta

								);

								// link code sequenece
								if(at ==0)
								{
									for (int nl = 0 + virtual_first_link_delta; nl < it->second.m_link_size - virtual_last_link_delta - 1; ++nl)  // check all movements
									{
										int link_seq_no = it->second.path_link_vector[nl];
										int next_link_seq_no = it->second.path_link_vector[nl + 1];

										if (g_link_vector[link_seq_no].VDF_period[tau].turn_link_count_map.find(next_link_seq_no) == g_link_vector[link_seq_no].VDF_period[tau].turn_link_count_map.end())
											g_link_vector[link_seq_no].VDF_period[tau].turn_link_count_map[next_link_seq_no] = volume;
										else 
											g_link_vector[link_seq_no].VDF_period[tau].turn_link_count_map[next_link_seq_no] += volume;

									}
								}

								/* Format and print various data */
								for (int ni = 0 + virtual_first_link_delta; ni < it->second.m_node_size - virtual_last_link_delta; ++ni)
								{
									fprintf(g_pFilePathMOE, "%d;", g_node_vector[it->second.path_node_vector[ni]].node_id);
									//if (g_node_vector[it->second.path_node_vector[ni]].node_id == 87)
									//{
									//	int debug_i = 1;
									//}
								}

								fprintf(g_pFilePathMOE, ",");
								int link_seq_no;
								// link id sequence
								for (int nl = 0 + virtual_first_link_delta; nl < it->second.m_link_size - virtual_last_link_delta; ++nl)
								{
									link_seq_no = it->second.path_link_vector[nl];
									fprintf(g_pFilePathMOE, "%s;", g_link_vector[link_seq_no].link_id.c_str());
								}
								fprintf(g_pFilePathMOE, ",");

								if (it->second.m_node_size - virtual_first_link_delta - virtual_last_link_delta >= 2)
								{
									fprintf(g_pFilePathMOE, "\"LINESTRING (");
									for (int ni = 0 + virtual_first_link_delta; ni < it->second.m_node_size - virtual_last_link_delta; ++ni)
									{
										fprintf(g_pFilePathMOE, "%f %f", g_node_vector[it->second.path_node_vector[ni]].x,
											g_node_vector[it->second.path_node_vector[ni]].y);

										if (ni != it->second.m_node_size - virtual_last_link_delta - 1)
											fprintf(g_pFilePathMOE, ", ");
									}

									fprintf(g_pFilePathMOE, ")\"");
								}

								fprintf(g_pFilePathMOE, ",");
								if (assignment.g_number_of_nodes < 5000 || (assignment.g_number_of_nodes >= 5000 && it->second.path_volume >= 5))  // critcal path volume
								{
									//fprintf(g_pFilePathMOE, ",");
									// link type name sequenece
									//for (int nl = 0 + virtual_first_link_delta; nl < it->second.m_link_size - virtual_last_link_delta; ++nl)
									//{
									//    link_seq_no = it->second.path_link_vector[nl];
									//    fprintf(g_pFilePathMOE, "%s;", g_link_vector[link_seq_no].link_type_name.c_str());
									//}


									// link code sequenece
									for (int nl = 0 + virtual_first_link_delta; nl < it->second.m_link_size - virtual_last_link_delta; ++nl)
									{
										link_seq_no = it->second.path_link_vector[nl];

										if (g_link_vector[link_seq_no].link_specifical_flag_str.size() >= 1)
											fprintf(g_pFilePathMOE, "%s;", g_link_vector[link_seq_no].link_specifical_flag_str.c_str());
									}
									fprintf(g_pFilePathMOE, ",");

									// link link_distance_VDF sequenece
									for (int nl = 0 + virtual_first_link_delta; nl < it->second.m_link_size - virtual_last_link_delta; ++nl)
									{
										link_seq_no = it->second.path_link_vector[nl];
										fprintf(g_pFilePathMOE, "%.3f;", g_link_vector[link_seq_no].link_avg_travel_time_per_period[tau][at] - g_link_vector[link_seq_no].free_flow_travel_time_in_min);
									}
									fprintf(g_pFilePathMOE, ",");

									// link FFTT sequenece
									for (int nl = 0 + virtual_first_link_delta; nl < it->second.m_link_size - virtual_last_link_delta; ++nl)
									{
										link_seq_no = it->second.path_link_vector[nl];
										fprintf(g_pFilePathMOE, "%.3f;", g_link_vector[link_seq_no].free_flow_travel_time_in_min);
									}
									fprintf(g_pFilePathMOE, ",");
								}
								else
								{
									fprintf(g_pFilePathMOE, ",,,");
								}

								int scenaro_no = assignment.g_active_DTAscenario_map[assignment.active_scenario_index];
								if (assignment.g_DTA_scenario_vector[scenaro_no].dtm_flag && global_od_impact_flag_across_all_mode_types)
								{
									fprintf(g_pFilePathMOE, "%d,%d,%d,",
										//									//	it->second.path_sensor_link_vector.size()*1.0f, 

										// "DTM_OD_impact_all_mode_types,DTM_OD_impact_mode_type,DTM_path_impact,DTM_#_of_lane_closure_links,DTM_new_path_generated,DTM_volume_before,DTM_volume_after,DTM_volume_diff, ");
										global_od_impact_flag_across_all_mode_types,
										impacted_path_flag,
										it->second.path_SA_link_vector.size()
									);

									for (int nl = 0; nl < it->second.path_SA_link_vector.size(); ++nl)
									{
										link_seq_no = it->second.path_SA_link_vector[nl];

										fprintf(g_pFilePathMOE, "%s;", g_link_vector[link_seq_no].link_id.c_str());
									}
									

										fprintf(g_pFilePathMOE, ",%d,%4f, %.4f, %4f",
										it->second.b_RT_new_path_flag,
										volume_before_dtm,
										volume_after_dtm,
										volume_diff_dtm
									);
								}
								else
								{
									fprintf(g_pFilePathMOE, "0,,,,,,,,");

								}

								fprintf(g_pFilePathMOE, ",");
								//for (int nt = 0 + virtual_first_link_delta; nt < min(5000-1, it->second.m_node_size - virtual_last_link_delta); ++nt)
								//{
								//    if(path_time_vector[nt] < assignment.g_LoadingEndTimeInMin)
								//    {
								//    fprintf(g_pFilePathMOE, "%s;", g_time_coding(path_time_vector[nt]).c_str());
								//    }
								//}
								//fprintf(g_pFilePathMOE, ",");

								//for (int nt = 0 + virtual_first_link_delta; nt < it->second.m_link_size+1 - virtual_last_link_delta; ++nt)
								//    fprintf(g_pFilePathMOE, "%.2f;", path_time_vector[nt]);

								//fprintf(g_pFilePathMOE, ",");

								//for (int nt = 0 + virtual_first_link_delta; nt < it->second.m_link_size - virtual_last_link_delta; ++nt)
								//    fprintf(g_pFilePathMOE, "%.2f;", path_time_vector[nt+1]- path_time_vector[nt]);

								//fprintf(g_pFilePathMOE, ",");
								// output the TT and vol of column updating

								// stage 1:
								//for (int iteration_number = 0; iteration_number < min(20, assignment.g_number_of_column_updating_iterations); iteration_number++)
								//{
								//    double TT = -1;
								//    if (it->second.path_time_per_iteration_map.find(iteration_number) != it->second.path_time_per_iteration_map.end())
								//    {
								//        TT = it->second.path_time_per_iteration_map[iteration_number];
								//    }


								//    fprintf(g_pFilePathMOE, "%f,", TT);

								//}
								//for (int iteration_number = 0; iteration_number < min(20, assignment.g_number_of_column_updating_iterations); iteration_number++)
								//{
								//    double Vol = 0;
								//    if (it->second.path_volume_per_iteration_map.find(iteration_number) != it->second.path_volume_per_iteration_map.end())
								//    {
								//        Vol = it->second.path_volume_per_iteration_map[iteration_number];
								//    }

								//    fprintf(g_pFilePathMOE, "%f,", Vol);

								//}
								// stage II: ODME
								//for (int iteration_number = 0; iteration_number < min(20, assignment.g_number_of_ODME_iterations); iteration_number++)
								//{
								//	double TT = -1;
								//	if (it->second.path_time_per_iteration_ODME_map.find(iteration_number) != it->second.path_time_per_iteration_ODME_map.end())
								//	{
								//		TT = it->second.path_time_per_iteration_ODME_map[iteration_number];
								//	}


								//	fprintf(g_pFilePathMOE, "%f,", TT);

								//}
								//for (int iteration_number = 0; iteration_number < min(20, assignment.g_number_of_ODME_iterations); iteration_number++)
								//{
								//	double Vol = 0;
								//	if (it->second.path_volume_per_iteration_ODME_map.find(iteration_number) != it->second.path_volume_per_iteration_ODME_map.end())
								//	{
								//		Vol = it->second.path_volume_per_iteration_ODME_map[iteration_number];
								//	}

								//	fprintf(g_pFilePathMOE, "%f,", Vol);

								//}

								//stage III:

								// output the TT and vol of sensitivity analysis
								for (int iteration_number = 0; iteration_number <= assignment.g_number_of_sensitivity_analysis_iterations_for_dtm; iteration_number++)
								{
									double TT = -1;
									if (it->second.path_time_per_iteration_SA_map.find(iteration_number) != it->second.path_time_per_iteration_SA_map.end())
									{
										TT = it->second.path_time_per_iteration_SA_map[iteration_number];
									}

									fprintf(g_pFilePathMOE, "%f,", TT);

								}
								double prev_value = 0;
								for (int iteration_number = 0; iteration_number <= assignment.g_number_of_sensitivity_analysis_iterations_for_dtm; iteration_number++)
								{
									double Vol = 0;
									if (it->second.path_volume_per_iteration_SA_map.find(iteration_number) != it->second.path_volume_per_iteration_SA_map.end())
									{
										Vol = it->second.path_volume_per_iteration_SA_map[iteration_number];
									}

									fprintf(g_pFilePathMOE, "%f,", Vol - prev_value);
									prev_value = Vol;

								}



								fprintf(g_pFilePathMOE, "\n");
								count++;
							}

						}  // for each route
						p_column_pool->avg_travel_time[assignment.active_scenario_index] = od_trave_time / p_column_pool->od_volume[assignment.active_scenario_index];
						p_column_pool->avg_distance[assignment.active_scenario_index] = od_distance_km / p_column_pool->od_volume[assignment.active_scenario_index];
						p_column_pool->avg_cost[assignment.active_scenario_index] = od_toll / p_column_pool->od_volume[assignment.active_scenario_index];

					} // od level
				}
			}
		}
	}
	fclose(g_pFilePathMOE);

}



void g_output_assignment_result(Assignment& assignment, int subarea_id)
{
	g_record_link_district_performance_per_scenario(assignment, 1);
	g_output_route_assignment_results(assignment, 1);

	int b_debug_detail_flag = 0;
	FILE* g_pFileLinkMOE = nullptr;

	double total_VMT = 0;
	double total_VHT = 0;
	double total_PMT = 0;
	double total_PHT = 0;
	double total_link_volume = 0;
	double total_link_speed = 0;
	double total_link_speed_ratio = 0;
	int MOE_count = 0;

	if (subarea_id == 0)
	{
		assignment.summary_file << "Output Link Performance:" << '\n';

		char link_performance_file_name[50];
		int scenario_no = assignment.g_active_DTAscenario_map[assignment.active_scenario_index];
		sprintf(link_performance_file_name, "link_performance_s%d_%s.csv", assignment.active_scenario_index, assignment.g_DTA_scenario_vector[scenario_no].scenario_name.c_str());

		dtalog.output() << "[STATUS INFO] writing file " << link_performance_file_name << '\n';
		g_DTA_log_file << "[STATUS INFO] writing file " << link_performance_file_name << '\n';
		fopen_ss(&g_pFileLinkMOE, link_performance_file_name, "w");
		if (!g_pFileLinkMOE)
		{
			dtalog.output() << "[ERROR] File " << link_performance_file_name << " cannot be opened." << '\n';
			g_DTA_log_file << "[ERROR] File " << link_performance_file_name << " cannot be opened." << '\n';
			return; 
		}


	}
	else
	{
		assignment.summary_file << "Output subarea Link Performance:" << '\n';

		dtalog.output() << "[STATUS INFO] writing subarea_link_performance.csv.." << '\n';
		g_DTA_log_file << "[STATUS INFO] writing subarea_link_performance.csv.." << '\n';
		fopen_ss(&g_pFileLinkMOE, "subarea_link_performance.csv", "w");
		if (!g_pFileLinkMOE)
		{
			dtalog.output() << "[ERROR] File subarea_link_performance.csv cannot be opened." << '\n';
			g_DTA_log_file << "[ERROR] File subarea_link_performance.csv cannot be opened." << '\n';
			return;
		}
	}


	int ref_volume_count = 0;
	float total_ref_volume_dev_abs_percentage = 0;

	{

		fprintf(g_pFileLinkMOE, "link_seq_id,link_id,vdf_type,from_node_id,to_node_id,lanes,distance_km,distance_mile,fftt,meso_link_id,meso_link_incoming_volume,tmc,tmc_corridor_name,tmc_corridor_id,tmc_road_order,tmc_road_sequence,subarea_id,link_type,link_type_code,vdf_code,time_period,volume,ref_volume,ref_diff,restricted_turn_nodes,auto_turn_volume,ODME_volume_before,ODME_volume_after,ODME_volume_diff,ODME_obs_count,ODME_upperbound_capacity_type,ODME_obs_count_dev,");

		fprintf(g_pFileLinkMOE, "preload_volume,mode_type_volume,travel_time,speed_kmph,speed_mph,speed_ratio,VOC,DOC,lane_capacity,link_capacity,queue,total_simu_waiting_time_in_min,avg_simu_waiting_time_in_min,plf,lanes,D_per_hour_per_lane,QVDF_cd,QVDF_n,P,severe_congestion_duration_in_h,vf,v_congestion_cutoff,QVDF_cp,QVDF_s,QVDF_v,vt2,VMT,VHT,PMT,PHT,PDT_vf,PDT_vc,geometry,");

		for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
			fprintf(g_pFileLinkMOE, "moving_agent_vol_%s,", assignment.g_ModeTypeVector[at].mode_type.c_str());

		for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
			fprintf(g_pFileLinkMOE, "MEU_vol_%s,", assignment.g_ModeTypeVector[at].mode_type.c_str());

		for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
			fprintf(g_pFileLinkMOE, "mode_link_cap_%s,", assignment.g_ModeTypeVector[at].mode_type.c_str());

		for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
			fprintf(g_pFileLinkMOE, "mode_DOC_%s,", assignment.g_ModeTypeVector[at].mode_type.c_str());

		for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
			fprintf(g_pFileLinkMOE, "mode_TT_%s,", assignment.g_ModeTypeVector[at].mode_type.c_str());

		//for (int og = 0; og < assignment.g_number_of_analysis_districts; ++og)
		//	for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
		//	{
		//		fprintf(g_pFileLinkMOE, "moving_agent_vol_district_%d_%s,", og, assignment.g_ModeTypeVector[at].mode_type.c_str());

		//	}


//		fprintf(g_pFileLinkMOE, "obs_count,upper_bound_type,dev,");

		for (int t = 6 * 60; t < 20 * 60; t += 15)
		{
			int hour = t / 60;
			int minute = t - hour * 60;

			fprintf(g_pFileLinkMOE, "v%02d:%02d,", hour, minute);
		}

		fprintf(g_pFileLinkMOE, "DTM_type,DTM_volume_before,DTM_volume_after,DTM_volume_diff,DTM_speed_before,DTM_speed_after,DTM_speed_diff,");
		fprintf(g_pFileLinkMOE, "DTM_DoC_before,DTM_DoC_after,DTM_Doc_diff,DTM_P_before,DTM_P_after,DTM_P_diff,");

		fprintf(g_pFileLinkMOE, "notes\n");

		//Initialization for all nodes
		int physical_link_count = 0; 
		for (int i = 0; i < g_link_vector.size(); ++i)
		{
			// virtual connectors
			if (g_link_vector[i].link_type_si[0] == -1)
				continue;

			physical_link_count++; 
			if (subarea_id == 1 && g_link_vector[i].subarea_id <= 0)
			{
				continue;
			}

			string vdf_type_str;

			if (g_link_vector[i].vdf_type == bpr_vdf)
			{
				vdf_type_str = "bpr";
			}
			else
			{
				vdf_type_str = "qvdf";
			}

			float total_vehicle_volume = 0;
			int scenario_code_count = 0;
			for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
			{
				total_vehicle_volume += g_link_vector[i].VDF_period[tau].volume_before_dtm;
				total_vehicle_volume += g_link_vector[i].total_volume_for_all_mode_types_per_period[tau] + g_link_vector[i].VDF_period[tau].preload;
				scenario_code_count += g_link_vector[i].VDF_period[tau].dtm_scenario_code.size();

			}

			//if (scenario_code_count == 0 && total_vehicle_volume < 0.001 && g_link_vector.size() > 20000)
			//	continue;

			for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
			{

				if (assignment.g_DemandPeriodVector[tau].number_of_demand_files == 0)
					continue;



				float speed = g_link_vector[i].free_speed;  // default speed


				if (g_link_vector[i].VDF_period[tau].avg_travel_time_0 > 0.001f)
					speed = g_link_vector[i].link_distance_VDF / max(0.1, g_link_vector[i].VDF_period[tau].avg_travel_time_0) * 60.0  ;

				float speed_ratio = speed / max(0.001, g_link_vector[i].free_speed);  // default speed
				float vehicle_volume = g_link_vector[i].total_volume_for_all_mode_types_per_period[tau] + g_link_vector[i].VDF_period[tau].preload;
				float ref_volume_diff = 0;

				if (g_link_vector[i].VDF_period[tau].ref_link_volume > 1)
				{
					ref_volume_diff = vehicle_volume - g_link_vector[i].VDF_period[tau].ref_link_volume;

					ref_volume_count++;
					total_ref_volume_dev_abs_percentage += fabs(ref_volume_diff / g_link_vector[i].VDF_period[tau].ref_link_volume * 100);
				}


				float mode_type_volume = g_link_vector[i].total_agent_volume_for_all_mode_types_per_period[tau] + g_link_vector[i].VDF_period[tau].preload;
				//VMT,VHT,PMT,PHT,PDT
				float preload = g_link_vector[i].VDF_period[tau].preload;
				float VMT = vehicle_volume * g_link_vector[i].link_distance_mile;

				float VHT = vehicle_volume * g_link_vector[i].VDF_period[tau].avg_travel_time_0 / 60.0;
				float PMT = mode_type_volume * g_link_vector[i].link_distance_mile;
				float PHT = mode_type_volume * g_link_vector[i].VDF_period[tau].avg_travel_time_0 / 60.0;
				float PDT_vf = mode_type_volume * (g_link_vector[i].VDF_period[tau].avg_travel_time_0 - g_link_vector[i].VDF_period[tau].FFTT_at[0]) / 60.0;
				float PDT_vc = max(0.0, mode_type_volume * (g_link_vector[i].VDF_period[tau].avg_travel_time_0 - g_link_vector[i].VDF_period[tau].FFTT_at[0] * g_link_vector[i].free_speed / max(0.001f, g_link_vector[i].v_congestion_cutoff)) / 60.0);



				total_VMT += VMT;
				total_VHT += VHT;
				total_PMT += PMT;
				total_PHT += PHT;
				total_link_volume += vehicle_volume;
				total_link_speed += speed;
				total_link_speed_ratio += speed_ratio;

				MOE_count++;

				fprintf(g_pFileLinkMOE, "%d,%s,%s,%d,%d,%f,%f,%f,%f,%d,%d,",

					//end with obs_volume_diff
					physical_link_count,
					g_link_vector[i].link_id.c_str(),
					vdf_type_str.c_str(),
					g_node_vector[g_link_vector[i].from_node_seq_no].node_id,
					g_node_vector[g_link_vector[i].to_node_seq_no].node_id,
					g_link_vector[i].VDF_period[tau].nlanes,
					g_link_vector[i].link_distance_km,
					g_link_vector[i].link_distance_mile,
					g_link_vector[i].free_flow_travel_time_in_min,
					g_link_vector[i].meso_link_id,
					g_link_vector[i].total_simulated_meso_link_incoming_volume);

				fprintf(g_pFileLinkMOE, "%s,%s,%d,%d,%d,%d,%d,%s,%s,%s,",
					g_link_vector[i].tmc_code.c_str(),
					g_link_vector[i].tmc_corridor_name.c_str(),
					g_link_vector[i].tmc_corridor_id,
					g_link_vector[i].tmc_road_order,
					g_link_vector[i].tmc_road_sequence,
					g_link_vector[i].subarea_id,
					g_link_vector[i].link_type_si[assignment.active_scenario_index],
					g_link_vector[i].link_type_code.c_str(),
					g_link_vector[i].vdf_code.c_str(),
					assignment.g_DemandPeriodVector[tau].time_period.c_str());

				fprintf(g_pFileLinkMOE, "%.3f,",
					vehicle_volume);

				fprintf(g_pFileLinkMOE, "%.3f,%.3f,",
					g_link_vector[i].VDF_period[tau].ref_link_volume,
					ref_volume_diff);
				fprintf(g_pFileLinkMOE, "%s,", g_link_vector[i].VDF_period[tau].restricted_turn_nodes_str.c_str());

				std::map<int, float>::iterator it, it_begin, it_end;
				
				it_begin = g_link_vector[i].VDF_period[tau].turn_link_count_map.begin();
				it_end = g_link_vector[i].VDF_period[tau].turn_link_count_map.end();

				for (it = it_begin; it != it_end; ++it)  //
				{
					int link_no = it->first;  // the first key is link_no
					int to_node_id = g_node_vector[g_link_vector[link_no].to_node_seq_no].node_id;  // outgoing node id of next link 
					fprintf(g_pFileLinkMOE, "%d:%d;", to_node_id, (int)(it->second));
				}


				fprintf(g_pFileLinkMOE, ",");



				fprintf(g_pFileLinkMOE, "%.3f,%.3f,%.3f,",
					g_link_vector[i].VDF_period[tau].volume_before_odme,
					g_link_vector[i].VDF_period[tau].volume_after_odme,
					g_link_vector[i].VDF_period[tau].volume_after_odme - g_link_vector[i].VDF_period[tau].volume_before_odme);

				if (g_link_vector[i].VDF_period[tau].obs_count[assignment.active_scenario_index] >= 1)
				{
					fprintf(g_pFileLinkMOE, "%.3f,%.3f,%.3f,",
						g_link_vector[i].VDF_period[tau].obs_count[assignment.active_scenario_index],
						g_link_vector[i].VDF_period[tau].upper_bound_flag,
						g_link_vector[i].VDF_period[tau].obs_count[assignment.active_scenario_index] - g_link_vector[i].VDF_period[tau].volume_after_odme);
				}
				else
				{
					fprintf(g_pFileLinkMOE, ",,,");
				}

				int scenario_no = assignment.g_active_DTAscenario_map[assignment.active_scenario_index];
				int mode_type_index = 0; 
				fprintf(g_pFileLinkMOE, "%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,",
					preload,
					mode_type_volume,
					g_link_vector[i].link_avg_travel_time_per_period[tau][0],
					speed,  /* 60.0 is used to convert min to hour */
					speed / 1.609,  /* 60.0 is used to convert min to hour */
					speed_ratio,
					g_link_vector[i].VDF_period[tau].DOC_mode[0],
					g_link_vector[i].VDF_period[tau].DOC_mode[0],
					g_link_vector[i].VDF_period[tau].lane_based_ultimate_hourly_capacity,
					g_link_vector[i].VDF_period[tau].lane_based_ultimate_hourly_capacity * g_link_vector[i].recorded_lanes_per_period_per_at[tau][mode_type_index][scenario_no],
					g_link_vector[i].VDF_period[tau].queue_length,
					max(0.0, g_link_vector[i].total_simulated_delay_in_min));


				fprintf(g_pFileLinkMOE, "%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,",

					max(0.0, g_link_vector[i].total_simulated_delay_in_min) / max(1.0f, vehicle_volume),
					g_link_vector[i].VDF_period[tau].Q_peak_load_factor,
					g_link_vector[i].VDF_period[tau].nlanes,
					g_link_vector[i].VDF_period[tau].lane_based_D,
					g_link_vector[i].VDF_period[tau].Q_cd,
					g_link_vector[i].VDF_period[tau].Q_n,
					g_link_vector[i].VDF_period[tau].P,
					g_link_vector[i].VDF_period[tau].Severe_Congestion_P,
					g_link_vector[i].VDF_period[tau].vf,
					g_link_vector[i].VDF_period[tau].v_congestion_cutoff);


				fprintf(g_pFileLinkMOE, "%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,",
					g_link_vector[i].VDF_period[tau].Q_cp,
					g_link_vector[i].VDF_period[tau].Q_s,
					g_link_vector[i].VDF_period[tau].avg_queue_speed,
					g_link_vector[i].VDF_period[tau].vt2,
					VMT,
					VHT,
					PMT,
					PHT,
					PDT_vf,
					PDT_vc);

				fprintf(g_pFileLinkMOE, "\"%s\",",
					g_link_vector[i].geometry.c_str());


				for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
					fprintf(g_pFileLinkMOE, "%.3f,", g_link_vector[i].volume_per_mode_type_per_period[tau][at]);

				for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
					fprintf(g_pFileLinkMOE, "%.3f,", g_link_vector[i].converted_MEU_volume_per_period_per_at[tau][at]);

				for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
					fprintf(g_pFileLinkMOE, "%.3f,", g_link_vector[i].VDF_period[tau].capacity_at[at] * g_link_vector[i].number_of_lanes_si[assignment.active_scenario_index]); //



				for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
				{
					fprintf(g_pFileLinkMOE, "%f,", g_link_vector[i].VDF_period[tau].DOC_mode[at]);

					if (g_link_vector[i].from_node_id == 1 && g_link_vector[i].to_node_id == 3)
					{

						int idebug = 1;

					}


				}

				for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
					fprintf(g_pFileLinkMOE, "%.3f,", g_link_vector[i].link_avg_travel_time_per_period[tau][at]);


				for (int t = 6 * 60; t < 20 * 60; t += 15)
				{
					float speed = g_link_vector[i].get_model_15_min_speed(t);
					fprintf(g_pFileLinkMOE, "%.3f,", speed);
				}

				if(g_link_vector[i].VDF_period[tau].dtm_scenario_code.size()>0)
				{
				fprintf(g_pFileLinkMOE, "%s,", g_link_vector[i].VDF_period[tau].dtm_scenario_code.c_str());
				}
				else
				{
					fprintf(g_pFileLinkMOE, ",");

				}

				fprintf(g_pFileLinkMOE, "%.3f,%.3f,%.3f,", g_link_vector[i].VDF_period[tau].volume_before_dtm,
					g_link_vector[i].VDF_period[tau].volume_after_dtm,
					g_link_vector[i].VDF_period[tau].volume_after_dtm - g_link_vector[i].VDF_period[tau].volume_before_dtm);

				fprintf(g_pFileLinkMOE, "%.3f,%.3f,%.3f,", g_link_vector[i].VDF_period[tau].speed_before_dtm,
					g_link_vector[i].VDF_period[tau].speed_after_dtm,
					g_link_vector[i].VDF_period[tau].speed_after_dtm - g_link_vector[i].VDF_period[tau].speed_before_dtm);

				fprintf(g_pFileLinkMOE, "%.3f,%.3f,%.3f,", g_link_vector[i].VDF_period[tau].DoC_before_dtm,
					g_link_vector[i].VDF_period[tau].DoC_after_dtm,
					g_link_vector[i].VDF_period[tau].DoC_after_dtm - g_link_vector[i].VDF_period[tau].DoC_before_dtm);


				fprintf(g_pFileLinkMOE, "%.3f,%.3f,%.3f,", g_link_vector[i].VDF_period[tau].P_before_dtm,
					g_link_vector[i].VDF_period[tau].P_after_dtm,
					g_link_vector[i].VDF_period[tau].P_after_dtm - g_link_vector[i].VDF_period[tau].P_before_dtm);


				fprintf(g_pFileLinkMOE, "period-based\n");
			}

		}
		fclose(g_pFileLinkMOE);
	}

	assignment.summary_file << ",ref_link_vol_count=," << ref_volume_count << "," << "MAPE=," << total_ref_volume_dev_abs_percentage / max(1, ref_volume_count) << "%" << '\n';

	assignment.summary_file << ",VMT=," << total_VMT << "," << "VKT=," << total_VMT * 1.6090 << '\n';
	assignment.summary_file << ",VHT=," << total_VHT << '\n';
	assignment.summary_file << ",network vehicle speed (MPH) =," << total_VMT / max(0.0001, total_VHT) << ",";
	assignment.summary_file << ",network vehicle speed (KPH) =," << total_VMT * 1.609 / max(0.0001, total_VHT) << '\n';

	assignment.summary_file << ",PMT=," << total_PMT << ",";
	assignment.summary_file << "PKT=," << total_PMT * 1.609 << '\n';
	assignment.summary_file << ",PHT=," << total_PHT << '\n';
	assignment.summary_file << ",network person speed (MPH) =," << total_PMT / max(0.0001, total_PHT) << ",";
	assignment.summary_file << ",network person speed (KPH) =," << total_PMT * 1.609 / max(0.0001, total_PHT) << '\n';


	assignment.summary_file << ",simple avg link volume=," << total_link_volume / max(1, MOE_count) << '\n';
	assignment.summary_file << ",simple avg link speed=," << total_link_speed / max(1, MOE_count) << '\n';
	assignment.summary_file << ",simple avg link speed ratio=," << total_link_speed_ratio / max(1, MOE_count) << '\n';



	if (subarea_id == 1)
		return;
	g_output_district_performance_result(assignment);
	//	g_output_dynamic_queue_profile();
}

void g_output_choice_set_result(Assignment& assignment)
{


	double path_time_vector[MAX_LINK_SIZE_IN_A_PATH];
	char route_assignment_file_name[50];

	int scenario_no = assignment.g_active_DTAscenario_map[assignment.active_scenario_index];
	sprintf(route_assignment_file_name, "choice_set_output_%d_%s.csv", assignment.active_scenario_index, assignment.g_DTA_scenario_vector[scenario_no].scenario_name.c_str());

	FILE* g_pFilePathMOE = nullptr;
	fopen_ss(&g_pFilePathMOE, route_assignment_file_name, "w");
	//dtalog.output() << "[STATUS INFO] writing choice_set_output.csv.." << '\n';
	//g_DTA_log_file << "[STATUS INFO] writing choice_set_output.csv.." << '\n';

	if (!g_pFilePathMOE)
	{
		dtalog.output() << "[ERROR] File " << route_assignment_file_name << " cannot be opened." << '\n';
		g_DTA_log_file << "[ERROR] File " << route_assignment_file_name << " cannot be opened." << '\n';
		return;
	}

	fprintf(g_pFilePathMOE, "first_column,multi_dim_choice_id,mode_tag,demand_period_tag,spatial_tag,travel_purpose_tag,data_tag,choice_alternative_id,activity_travel_pattern_id,mode_chain,demand_period_chain,zone_chain,volume,travel_time,travel_distance,");
	fprintf(g_pFilePathMOE, "geometry,\n");

	std::map<std::string, CChoiceSet>::iterator it, it_begin, it_end;
	CColumnVector* p_column_pool;
	it_begin = assignment.g_ChoiceSetMap.begin();
	it_end = assignment.g_ChoiceSetMap.end();

	for (it = it_begin; it != it_end; ++it)  // for each choice set, 1st loop
	{
		for (int alt_no = 0; alt_no < it->second.choice_alt_vector.size(); alt_no++)  // for each alternative, 2nd loop
		{
			fprintf(g_pFilePathMOE, "0,");  //multi_dim_choice_id,mode_tag,demand_period_tag,spatial_tag,travel_purpose_tag,data_tag,
			fprintf(g_pFilePathMOE, "%s,%s,%s,%s,%s,%s,", it->first.c_str(),
				it->second.mode_tag.c_str(), it->second.demand_period_tag.c_str(), it->second.spatial_tag.c_str(), it->second.travel_purpose_tag.c_str(), it->second.data_tag.c_str());

			int alt_travel_time = 0;
			int alt_travel_cost = 0;


			for (int trip_no = 0; trip_no < it->second.choice_alt_vector[alt_no].activity_chain.size(); trip_no++)  // for each trip leg, 3rd loop
			{
				int from_zone_sindex = it->second.choice_alt_vector[alt_no].activity_chain[trip_no].o_zone_no;
				int to_zone_sindex = it->second.choice_alt_vector[alt_no].activity_chain[trip_no].d_zone_no;
				int at = it->second.choice_alt_vector[alt_no].activity_chain[trip_no].mode_no;
				int tau = it->second.choice_alt_vector[alt_no].activity_chain[trip_no].demand_peroid_no;

				p_column_pool = &(assignment.g_column_pool[from_zone_sindex][to_zone_sindex][at][tau]);


				std::map<int, CColumnPath>::iterator it_p, it_begin_p, it_end_p;
				{
					// scan through the map with different node sum for different continuous paths
					it_begin_p = p_column_pool->path_node_sequence_map.begin();  // for each route leg, 3rd loop
					it_end_p = p_column_pool->path_node_sequence_map.end();

					float path_toll = 0;
					double path_distance = 0;
					double path_distance_km = 0;
					double path_distance_ml = 0;
					float path_travel_time = 0;
					float path_travel_time_without_access_link = 0;
					float path_delay = 0;
					float path_FF_travel_time = 0;
					float time_stamp = 0;

					for (it_p = it_begin_p; it_p != it_end_p; ++it_p)
					{  // first path in the column pool

						for (int nl = 0; nl < it_p->second.m_link_size; ++nl)  // for each link in the route
						{

							int link_seq_no = it_p->second.path_link_vector[nl];
							if (g_link_vector[link_seq_no].link_type_si[assignment.active_scenario_index] >= 0)
							{
								path_toll += g_link_vector[link_seq_no].VDF_period[tau].toll[at][assignment.active_scenario_index];
								path_distance += g_link_vector[link_seq_no].link_distance_VDF;
								path_distance_km += g_link_vector[link_seq_no].link_distance_km;
								path_distance_ml += g_link_vector[link_seq_no].link_distance_mile;
								float link_travel_time = g_link_vector[link_seq_no].link_avg_travel_time_per_period[tau][at];
								path_travel_time += link_travel_time;

								if (g_link_vector[link_seq_no].link_type_si[assignment.active_scenario_index] != 1000)   // skip access link
								{
									path_travel_time_without_access_link += link_travel_time;
								}

								path_FF_travel_time += g_link_vector[link_seq_no].VDF_period[tau].FFTT_at[at];
								time_stamp += link_travel_time;
								path_time_vector[nl + 1] = time_stamp;

							}
						}

						alt_travel_time += path_travel_time;
						alt_travel_cost += path_distance;
						break;  // at the first path
					}



				}

			}// end of trip lags
			it->second.avg_travel_time = alt_travel_time;
			it->second.avg_travel_cost = alt_travel_cost;

			//choice_alternative_id,activity_travel_pattern_id,zone_chain,volume,
			fprintf(g_pFilePathMOE, "%d,%s,%s,%s,%s,%.1f,%.1f,%.1f,",
				it->second.choice_alt_vector[alt_no].choice_alternative_id,
				it->second.choice_alt_vector[alt_no].activity_travel_pattern_id.c_str(),
				assignment.g_ActivityTravelPatternMap[it->second.choice_alt_vector[alt_no].activity_travel_pattern_id].mode_chain_str.c_str(),
				assignment.g_ActivityTravelPatternMap[it->second.choice_alt_vector[alt_no].activity_travel_pattern_id].demand_period_chain_str.c_str(),
				it->second.choice_alt_vector[alt_no].activity_zone_chain_str.c_str(),
				it->second.choice_alt_vector[alt_no].volume,
				it->second.avg_travel_time,
				it->second.avg_travel_cost
			);

			fprintf(g_pFilePathMOE, "\"LINESTRING (");

			for (int trip_no = 0; trip_no < it->second.choice_alt_vector[alt_no].activity_chain.size(); trip_no++)  // for each trip leg, 3rd loop
			{
				int from_zone_sindex = it->second.choice_alt_vector[alt_no].activity_chain[trip_no].o_zone_no;
				int to_zone_sindex = it->second.choice_alt_vector[alt_no].activity_chain[trip_no].d_zone_no;
				int at = it->second.choice_alt_vector[alt_no].activity_chain[trip_no].mode_no;
				int tau = it->second.choice_alt_vector[alt_no].activity_chain[trip_no].demand_peroid_no;

				p_column_pool = &(assignment.g_column_pool[from_zone_sindex][to_zone_sindex][at][tau]);

				std::map<int, CColumnPath>::iterator it_p, it_begin_p, it_end_p;
				{
					// scan through the map with different node sum for different continuous paths
					it_begin_p = p_column_pool->path_node_sequence_map.begin();
					it_end_p = p_column_pool->path_node_sequence_map.end();

					for (it_p = it_begin_p; it_p != it_end_p; ++it_p)
					{
						for (int ni = 1; ni < it_p->second.m_node_size - 1; ++ni)					// for links in the, 4th loop
						{
							fprintf(g_pFilePathMOE, "%f %f,", g_node_vector[it_p->second.path_node_vector[ni]].x,
								g_node_vector[it_p->second.path_node_vector[ni]].y);

						}
						break;
					}
				}

			}  // for each alternative's trip leg
			fprintf(g_pFilePathMOE, ")\"");

			fprintf(g_pFilePathMOE, "\n");

		}
	}


	fclose(g_pFilePathMOE);

	//	g_output_dynamic_queue_profile();
}
void g_output_assignment_summary_result(Assignment& assignment, int subarea_id)
{

	int b_debug_detail_flag = 0;
	FILE* g_pFileLinkMOE = nullptr;

	double total_VMT = 0;
	double total_VHT = 0;
	double total_PMT = 0;
	double total_PHT = 0;
	double total_link_volume = 0;
	double total_link_speed = 0;
	double total_link_speed_ratio = 0;
	int MOE_count = 0;

	if (subarea_id == 0)
	{
		assignment.summary_file << "Output Link Performance Summary" << '\n';

		dtalog.output() << "[STATUS INFO] writing link_performances_summary.csv.." << '\n';
		g_DTA_log_file << "[STATUS INFO] writing link_performances_summary.csv.." << '\n';
		fopen_ss(&g_pFileLinkMOE, "link_performance_summary.csv", "w");
		if (!g_pFileLinkMOE)
		{
			dtalog.output() << "[ERROR] File link_performance_summary.csv cannot be opened." << '\n';
			g_DTA_log_file << "[ERROR] File link_performance_summary.csv cannot be opened." << '\n';
			return; 
		}


	}
	else
	{
		assignment.summary_file << "Output subarea Link Performance:" << '\n';

		dtalog.output() << "[STATUS INFO] writing subarea_link_performance.csv.." << '\n';
		g_DTA_log_file << "[STATUS INFO] writing subarea_link_performance.csv.." << '\n';
		fopen_ss(&g_pFileLinkMOE, "subarea_link_performance.csv", "w");
		if (!g_pFileLinkMOE)
		{
			dtalog.output() << "[ERROR] File subarea_link_performance.csv cannot be opened." << '\n';
			g_DTA_log_file << "[ERROR] File subarea_link_performance.csv cannot be opened." << '\n';
			return; 
		}
	}

	if (g_pFileLinkMOE != NULL)
	{

		int ref_volume_count = 0;
		float total_ref_volume_dev_abs_percentage = 0;


		fprintf(g_pFileLinkMOE, "link_seq_id,link_id,base_s0_link_type,from_node_id,to_node_id,geometry,distance_km,distance_mile,fftt,free_speed_base_mode,link_specific_speed_diff,capacity_base_mode,link_specific_capacity_diff,meso_link_id,tmc,tmc_corridor_name,tmc_corridor_id,tmc_road_order,tmc_road_sequence,subarea_id,");

		for (int sii = 0; sii < assignment.g_DTA_scenario_vector.size(); sii++)
		{
			int scenario_index = assignment.g_DTA_scenario_vector[sii].scenario_index;

			fprintf(g_pFileLinkMOE, "lanes_s%d,", scenario_index);
		}

		for (int sii = 0; sii < assignment.g_DTA_scenario_vector.size(); sii++)
		{
			int scenario_index = assignment.g_DTA_scenario_vector[sii].scenario_index;
			fprintf(g_pFileLinkMOE, "link_type_s%d,", scenario_index);
		}
		for (int sii = 0; sii < assignment.g_DTA_scenario_vector.size(); sii++)
		{
			int scenario_index = assignment.g_DTA_scenario_vector[sii].scenario_index;
			fprintf(g_pFileLinkMOE, "free_speed_s%d,", scenario_index);
		}
		for (int sii = 0; sii < assignment.g_DTA_scenario_vector.size(); sii++)
		{
			int scenario_index = assignment.g_DTA_scenario_vector[sii].scenario_index;
			fprintf(g_pFileLinkMOE, "capacity_s%d,", scenario_index);
		}

		for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
		{
			for (int sii = 0; sii < assignment.g_DTA_scenario_vector.size(); sii++)
			{
				int scenario_index = assignment.g_DTA_scenario_vector[sii].scenario_index;
				{

					for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
					{

						fprintf(g_pFileLinkMOE, "lanes_p%d_s%d_%s,", tau + 1, scenario_index, assignment.g_ModeTypeVector[at].mode_type.c_str());
					}

				}
			}

		}

		for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
		{
			for (int sii = 0; sii < assignment.g_DTA_scenario_vector.size(); sii++)
			{
				int scenario_index = assignment.g_DTA_scenario_vector[sii].scenario_index;
				{

					for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
					{

						fprintf(g_pFileLinkMOE, "VOL_p%d_s%d_%s,", tau+1, scenario_index, assignment.g_ModeTypeVector[at].mode_type.c_str());
					}

				}
			}

		}
		for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
		{
			for (int sii = 0; sii < assignment.g_DTA_scenario_vector.size(); sii++)
			{
				int scenario_index = assignment.g_DTA_scenario_vector[sii].scenario_index;
				{

					for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
					{

						fprintf(g_pFileLinkMOE, "MEU_p%d_s%d_%s,", tau+1, scenario_index, assignment.g_ModeTypeVector[at].mode_type.c_str());
					}

				}
			}

		}


		for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
		{
			for (int sii = 0; sii < assignment.g_DTA_scenario_vector.size(); sii++)
			{
				int scenario_index = assignment.g_DTA_scenario_vector[sii].scenario_index;
				{

					for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
					{

						fprintf(g_pFileLinkMOE, "CAP_p%d_s%d_%s,", tau+1, scenario_index, assignment.g_ModeTypeVector[at].mode_type.c_str());
					}

				}
			}

		}

		for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
		{
			for (int sii = 0; sii < assignment.g_DTA_scenario_vector.size(); sii++)
			{
				int scenario_index = assignment.g_DTA_scenario_vector[sii].scenario_index;
				{

					for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
					{

						fprintf(g_pFileLinkMOE, "PEN_p%d_s%d_%s,", tau + 1, scenario_index, assignment.g_ModeTypeVector[at].mode_type.c_str());
					}

				}
			}

		}

		for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
		{
			for (int sii = 0; sii < assignment.g_DTA_scenario_vector.size(); sii++)
			{
				int scenario_index = assignment.g_DTA_scenario_vector[sii].scenario_index;
				{

					for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
					{

						fprintf(g_pFileLinkMOE, "DOC_p%d_s%d_%s,", tau+1, scenario_index, assignment.g_ModeTypeVector[at].mode_type.c_str());
					}

				}
			}

		}

		for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
		{
			for (int sii = 0; sii < assignment.g_DTA_scenario_vector.size(); sii++)
			{
				int scenario_index = assignment.g_DTA_scenario_vector[sii].scenario_index;
				{

					for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
					{

						fprintf(g_pFileLinkMOE, "TT_p%d_s%d_%s,", tau+1, scenario_index, assignment.g_ModeTypeVector[at].mode_type.c_str());
					}

				}
			}

		}

		for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
		{
			for (int sii = 0; sii < assignment.g_DTA_scenario_vector.size(); sii++)
			{
				int scenario_index = assignment.g_DTA_scenario_vector[sii].scenario_index;
				{

					for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
					{

						fprintf(g_pFileLinkMOE, "CO2_p%d_s%d_%s,", tau + 1, scenario_index, assignment.g_ModeTypeVector[at].mode_type.c_str());
					}

				}
			}

		}


		for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
		{
			for (int sii = 0; sii < assignment.g_DTA_scenario_vector.size(); sii++)
			{
				int scenario_index = assignment.g_DTA_scenario_vector[sii].scenario_index;
				{

					for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
					{

						fprintf(g_pFileLinkMOE, "NOX_p%d_s%d_%s,", tau + 1, scenario_index, assignment.g_ModeTypeVector[at].mode_type.c_str());
					}

				}
			}

		}

		for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); ++tau)
		{
			for (int sii = 0; sii < assignment.g_DTA_scenario_vector.size(); sii++)
			{
				int scenario_index = assignment.g_DTA_scenario_vector[sii].scenario_index;

				int scenario_no = assignment.g_active_DTAscenario_map[scenario_index];
				for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
				{
					fprintf(g_pFileLinkMOE, "VOL_%s_%s_%s,", assignment.g_DemandPeriodVector[tau].demand_period.c_str(), assignment.g_DTA_scenario_vector[scenario_no].scenario_name.c_str(), assignment.g_ModeTypeVector[at].mode_type.c_str());
					fprintf(g_pFileLinkMOE, "MEU_%s_%s_%s,", assignment.g_DemandPeriodVector[tau].demand_period.c_str(), assignment.g_DTA_scenario_vector[scenario_no].scenario_name.c_str(), assignment.g_ModeTypeVector[at].mode_type.c_str());
					fprintf(g_pFileLinkMOE, "CAP_%s_%s_%s,", assignment.g_DemandPeriodVector[tau].demand_period.c_str(), assignment.g_DTA_scenario_vector[scenario_no].scenario_name.c_str(), assignment.g_ModeTypeVector[at].mode_type.c_str());
					fprintf(g_pFileLinkMOE, "PEN_%s_%s_%s,", assignment.g_DemandPeriodVector[tau].demand_period.c_str(), assignment.g_DTA_scenario_vector[scenario_no].scenario_name.c_str(), assignment.g_ModeTypeVector[at].mode_type.c_str());
					fprintf(g_pFileLinkMOE, "DOC_%s_%s_%s,", assignment.g_DemandPeriodVector[tau].demand_period.c_str(), assignment.g_DTA_scenario_vector[scenario_no].scenario_name.c_str(), assignment.g_ModeTypeVector[at].mode_type.c_str());
					fprintf(g_pFileLinkMOE, "TT_%s_%s_%s,", assignment.g_DemandPeriodVector[tau].demand_period.c_str(), assignment.g_DTA_scenario_vector[scenario_no].scenario_name.c_str(), assignment.g_ModeTypeVector[at].mode_type.c_str());
					fprintf(g_pFileLinkMOE, "CO2_%s_%s_%s,", assignment.g_DemandPeriodVector[tau].demand_period.c_str(), assignment.g_DTA_scenario_vector[scenario_no].scenario_name.c_str(), assignment.g_ModeTypeVector[at].mode_type.c_str());
					fprintf(g_pFileLinkMOE, "NOX_%s_%s_%s,", assignment.g_DemandPeriodVector[tau].demand_period.c_str(), assignment.g_DTA_scenario_vector[scenario_no].scenario_name.c_str(), assignment.g_ModeTypeVector[at].mode_type.c_str());
				}


			}
		}
		fprintf(g_pFileLinkMOE, "\n");
		//		fprintf(g_pFileLinkMOE, "obs_count,upper_bound_type,dev,");
				// header

				//Initialization for all nodes
		for (int i = 0; i < g_link_vector.size(); ++i)
		{
			// virtual connectors
			if (g_link_vector[i].link_type_si[0] == -1)
				continue;

			string vdf_type_str;


			float total_vehicle_volume = 0;
			int scenario_code_count = 0;
			for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
			{
				total_vehicle_volume += g_link_vector[i].VDF_period[tau].volume_before_dtm;
				total_vehicle_volume += g_link_vector[i].total_volume_for_all_mode_types_per_period[tau] + g_link_vector[i].VDF_period[tau].preload;
				scenario_code_count += g_link_vector[i].VDF_period[tau].dtm_scenario_code.size();

			}

			//if (scenario_code_count == 0 && total_vehicle_volume < 0.001 && g_link_vector.size() > 20000)
			//	continue;




			float speed = g_link_vector[i].free_speed;  // default speed


			fprintf(g_pFileLinkMOE, "%d,%s,%s,%d,%d,",
				i+1,
				//end with obs_volume_diff
				g_link_vector[i].link_id.c_str(),
				assignment.g_LinkTypeMap[g_link_vector[i].link_type_si[0]].link_type_name.c_str(),
				g_node_vector[g_link_vector[i].from_node_seq_no].node_id,
				g_node_vector[g_link_vector[i].to_node_seq_no].node_id);

			if (g_link_vector[i].geometry.size() > 1)
			{
				fprintf(g_pFileLinkMOE, "\"%s\",",
					g_link_vector[i].geometry.c_str());
			}
			else
			{

				fprintf(g_pFileLinkMOE, "\"LINESTRING (");
				fprintf(g_pFileLinkMOE, "%f %f,", g_node_vector[g_link_vector[i].from_node_seq_no].x, g_node_vector[g_link_vector[i].from_node_seq_no].y);

				fprintf(g_pFileLinkMOE, "%f %f", g_node_vector[g_link_vector[i].to_node_seq_no].x, g_node_vector[g_link_vector[i].to_node_seq_no].y);

				fprintf(g_pFileLinkMOE, ")\",");
			}

			int at_base = 0; 
			int tau_base = 0; 
			fprintf(g_pFileLinkMOE, "%f,%f,%f,%f,%f,%f,%f,%d,",

				g_link_vector[i].link_distance_km,
				g_link_vector[i].link_distance_mile,
				g_link_vector[i].free_flow_travel_time_in_min,

				g_link_vector[i].VDF_period[tau_base].free_speed_at[at_base],
				g_link_vector[i].VDF_period[tau_base].free_speed_diff_link_specific,
				g_link_vector[i].VDF_period[tau_base].capacity_at[at_base],
				g_link_vector[i].VDF_period[tau_base].capacity_diff_link_specific,

				g_link_vector[i].meso_link_id);

			fprintf(g_pFileLinkMOE, "%s,%s,%d,%d,%d,%d,",
				g_link_vector[i].tmc_code.c_str(),
				g_link_vector[i].tmc_corridor_name.c_str(),
				g_link_vector[i].tmc_corridor_id,
				g_link_vector[i].tmc_road_order,
				g_link_vector[i].tmc_road_sequence,
				g_link_vector[i].subarea_id);

			for (int sii = 0; sii < assignment.g_DTA_scenario_vector.size(); sii++)
			{
				int scenario_index = assignment.g_DTA_scenario_vector[sii].scenario_index;
				fprintf(g_pFileLinkMOE, "%.1f,", g_link_vector[i].number_of_lanes_si[scenario_index]);

			}

			for (int sii = 0; sii < assignment.g_DTA_scenario_vector.size(); sii++)
			{
				int scenario_index = assignment.g_DTA_scenario_vector[sii].scenario_index;
				fprintf(g_pFileLinkMOE, "%s,", assignment.g_LinkTypeMap[g_link_vector[i].link_type_si[scenario_index]].link_type_name.c_str());
			}


			for (int sii = 0; sii < assignment.g_DTA_scenario_vector.size(); sii++)
			{
				int scenario_index = assignment.g_DTA_scenario_vector[sii].scenario_index;
				fprintf(g_pFileLinkMOE, "%.1f,", g_link_vector[i].free_speed_si[scenario_index]);

			}


			for (int sii = 0; sii < assignment.g_DTA_scenario_vector.size(); sii++)
			{
				int scenario_index = assignment.g_DTA_scenario_vector[sii].scenario_index;
				fprintf(g_pFileLinkMOE, "%.1f,", g_link_vector[i].capacity_si[scenario_index]);

			}

			for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
			{
				for (int sii = 0; sii < assignment.g_DTA_scenario_vector.size(); sii++)
				{
					int scenario_index = assignment.g_DTA_scenario_vector[sii].scenario_index;
					int scenario_no = assignment.g_active_DTAscenario_map[scenario_index];
					{

						for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
						{

							fprintf(g_pFileLinkMOE, "%.1f,", g_link_vector[i].recorded_lanes_per_period_per_at[tau][at][scenario_no]);  // pk format

						}
					}
				}

			}

			float total_recorded_volume = 0; 
			for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
			{
				for (int sii = 0; sii < assignment.g_DTA_scenario_vector.size(); sii++)
				{
					int scenario_index = assignment.g_DTA_scenario_vector[sii].scenario_index;
					{
						int scenario_no = assignment.g_active_DTAscenario_map[scenario_index];

						for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
						{

							total_recorded_volume += g_link_vector[i].recorded_volume_per_period_per_at[tau][at][scenario_no];
							fprintf(g_pFileLinkMOE, "%.1f,", g_link_vector[i].recorded_volume_per_period_per_at[tau][at][scenario_no]);  // pk format

						}
					}
				}

			}

			int debug_flag;
			debug_flag = 1;

			for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
			{
				for (int sii = 0; sii < assignment.g_DTA_scenario_vector.size(); sii++)
				{
					int scenario_index = assignment.g_DTA_scenario_vector[sii].scenario_index;
					{
						int scenario_no = assignment.g_active_DTAscenario_map[scenario_index];

						for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
						{

							fprintf(g_pFileLinkMOE, "%.1f,", g_link_vector[i].recorded_MEU_per_period_per_at[tau][at][scenario_no]);  // MEU

						}
					}
				}

			}
			for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
			{
				for (int sii = 0; sii < assignment.g_DTA_scenario_vector.size(); sii++)
				{
					int scenario_index = assignment.g_DTA_scenario_vector[sii].scenario_index;
					{
						int scenario_no = assignment.g_active_DTAscenario_map[scenario_index];

						for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
						{

							fprintf(g_pFileLinkMOE, "%.1f,", g_link_vector[i].VDF_period[tau].capacity_at[at] * g_link_vector[i].number_of_lanes_si[scenario_no]);  // CAP

						}
					}
				}

			}

			for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
			{
				for (int sii = 0; sii < assignment.g_DTA_scenario_vector.size(); sii++)
				{
					int scenario_index = assignment.g_DTA_scenario_vector[sii].scenario_index;
					{
						int scenario_no = assignment.g_active_DTAscenario_map[scenario_index];

						for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
						{

							fprintf(g_pFileLinkMOE, "%.1f,", g_link_vector[i].penalty_si_at[tau][at][scenario_no]);  // PEN

						}
					}
				}

			}			


			for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
			{
				for (int sii = 0; sii < assignment.g_DTA_scenario_vector.size(); sii++)
				{
					int scenario_index = assignment.g_DTA_scenario_vector[sii].scenario_index;
					{
						int scenario_no = assignment.g_active_DTAscenario_map[scenario_index];

						for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
						{

							fprintf(g_pFileLinkMOE, "%.4f,", g_link_vector[i].recorded_DOC_per_period_per_at[tau][at][scenario_no]);  // DOC

						}
					}
				}

			}
			for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
			{
				for (int sii = 0; sii < assignment.g_DTA_scenario_vector.size(); sii++)
				{
					int scenario_index = assignment.g_DTA_scenario_vector[sii].scenario_index;
					{
						int scenario_no = assignment.g_active_DTAscenario_map[scenario_index];

						for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
						{

							fprintf(g_pFileLinkMOE, "%.1f,", g_link_vector[i].recorded_TT_per_period_per_at[tau][at][scenario_no]);  // TT

						}
					}
				}

			}

			for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
			{
				for (int sii = 0; sii < assignment.g_DTA_scenario_vector.size(); sii++)
				{
					int scenario_index = assignment.g_active_DTAscenario_map[assignment.g_DTA_scenario_vector[sii].scenario_index];
					{	int scenario_no = assignment.g_active_DTAscenario_map[scenario_index];

						for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
						{
						
							fprintf(g_pFileLinkMOE, "%.1f,", g_link_vector[i].recorded_CO2_per_period_per_at[tau][at][scenario_no]);  // TT

						}
					}
				}

			}

			for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
			{
				for (int sii = 0; sii < assignment.g_DTA_scenario_vector.size(); sii++)
				{
					int scenario_index = assignment.g_DTA_scenario_vector[sii].scenario_index;
					{
						int scenario_no = assignment.g_active_DTAscenario_map[scenario_index];
						for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
						{

							fprintf(g_pFileLinkMOE, "%.1f,", g_link_vector[i].recorded_NOX_per_period_per_at[tau][at][scenario_no]);  // TT

						}
					}
				}

			}
			//
			for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
			{
				for (int sii = 0; sii < assignment.g_DTA_scenario_vector.size(); sii++)
				{
					int scenario_index = assignment.g_DTA_scenario_vector[sii].scenario_index;
					{

						for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
						{
							int scenario_no = assignment.g_active_DTAscenario_map[scenario_index];
							fprintf(g_pFileLinkMOE, "%.1f,", g_link_vector[i].recorded_volume_per_period_per_at[tau][at][scenario_no]);
							fprintf(g_pFileLinkMOE, "%.1f,", g_link_vector[i].recorded_MEU_per_period_per_at[tau][at][scenario_no]);
							fprintf(g_pFileLinkMOE, "%.3f,", g_link_vector[i].recorded_capacity_per_period_per_at[tau][at][scenario_no]);
							fprintf(g_pFileLinkMOE, "%.3f,", g_link_vector[i].penalty_si_at[tau][at][scenario_no]);
							fprintf(g_pFileLinkMOE, "%.4f,", g_link_vector[i].recorded_DOC_per_period_per_at[tau][at][scenario_no]);
							fprintf(g_pFileLinkMOE, "%.3f,", g_link_vector[i].recorded_TT_per_period_per_at[tau][at][scenario_no]);
							fprintf(g_pFileLinkMOE, "%.3f,", g_link_vector[i].recorded_CO2_per_period_per_at[tau][at][scenario_no]);
							fprintf(g_pFileLinkMOE, "%.3f,", g_link_vector[i].recorded_NOX_per_period_per_at[tau][at][scenario_no]);

						}
					}
				}

			}



			fprintf(g_pFileLinkMOE, "\n");
		}


		fclose(g_pFileLinkMOE);
	}

}
void g_output_2_way_assignment_summary_result(Assignment& assignment, int subarea_id)
{

	int b_debug_detail_flag = 0;
	FILE* g_pFileLinkMOE = nullptr;

	double total_VMT = 0;
	double total_VHT = 0;
	double total_PMT = 0;
	double total_PHT = 0;
	double total_link_volume = 0;
	double total_link_speed = 0;
	double total_link_speed_ratio = 0;
	int MOE_count = 0;

	if (subarea_id == 0)
	{
		assignment.summary_file << "Output 2 Way Link Performance Summary" << '\n';

		dtalog.output() << "[STATUS INFO] writing link_performances_summary_2way.csv.." << '\n';
		g_DTA_log_file << "[STATUS INFO] writing link_performances_summary_2way.csv.." << '\n';
		fopen_ss(&g_pFileLinkMOE, "link_performance_summary_2way.csv", "w");
		if (!g_pFileLinkMOE)
		{
			dtalog.output() << "[ERROR] File link_performance_summary2_way.csv cannot be opened." << '\n';
			g_DTA_log_file << "[ERROR] File link_performance_summary2_way.csv cannot be opened." << '\n';
			return; 
		}


	}
	else
	{
		assignment.summary_file << "Output subarea Link Performance:" << '\n';

		dtalog.output() << "[STATUS INFO] writing subarea_link_performance.csv.." << '\n';
		g_DTA_log_file << "[STATUS INFO] writing subarea_link_performance.csv.." << '\n';
		fopen_ss(&g_pFileLinkMOE, "subarea_link_performance.csv", "w");
		if (!g_pFileLinkMOE)
		{
			dtalog.output() << "[ERROR] File subarea_link_performance.csv cannot be opened." << '\n';
			g_DTA_log_file << "[ERROR] File subarea_link_performance.csv cannot be opened." << '\n';
			return;
		}
	}

	if (g_pFileLinkMOE != NULL)
	{

		int ref_volume_count = 0;
		float total_ref_volume_dev_abs_percentage = 0;


		fprintf(g_pFileLinkMOE, "link_seq_id,link_id,link_type,from_node_id,to_node_id,geometry,distance_km,distance_mile,fftt,meso_link_id,tmc,tmc_corridor_name,tmc_corridor_id,tmc_road_order,tmc_road_sequence,subarea_id,");

		for (int sii = 0; sii < assignment.g_DTA_scenario_vector.size(); sii++)
		{
			int scenario_index = assignment.g_DTA_scenario_vector[sii].scenario_index;
			fprintf(g_pFileLinkMOE, "lanes_s%d,", scenario_index);
		}

		for (int sii = 0; sii < assignment.g_DTA_scenario_vector.size(); sii++)
		{
			int scenario_index = assignment.g_DTA_scenario_vector[sii].scenario_index;
			fprintf(g_pFileLinkMOE, "link_type_s%d,", scenario_index);
		}



		for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); ++tau)
		{
			for (int sii = 0; sii < assignment.g_DTA_scenario_vector.size(); sii++)
			{
				int scenario_index = assignment.g_DTA_scenario_vector[sii].scenario_index;
				{
					int scenario_no = assignment.g_active_DTAscenario_map[scenario_index];

					for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
					{

						fprintf(g_pFileLinkMOE, "V%s_%s_%s,", assignment.g_DemandPeriodVector[tau].demand_period.c_str(), assignment.g_DTA_scenario_vector[scenario_no].scenario_name.c_str(), assignment.g_ModeTypeVector[at].mode_type.c_str());
					}

				}
			}

		}
		for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); ++tau)
		{
			for (int sii = 0; sii < assignment.g_DTA_scenario_vector.size(); sii++)
			{
				int scenario_index = assignment.g_DTA_scenario_vector[sii].scenario_index;
				{
					int scenario_no = assignment.g_active_DTAscenario_map[scenario_index];

					for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
					{

						fprintf(g_pFileLinkMOE, "VAB_%s_%s_%s,", assignment.g_DemandPeriodVector[tau].demand_period.c_str(), assignment.g_DTA_scenario_vector[scenario_no].scenario_name.c_str(), assignment.g_ModeTypeVector[at].mode_type.c_str());
					}

				}
			}

		}
		for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); ++tau)
		{
			for (int sii = 0; sii < assignment.g_DTA_scenario_vector.size(); sii++)
			{
				int scenario_index = assignment.g_DTA_scenario_vector[sii].scenario_index;
				{
					int scenario_no = assignment.g_active_DTAscenario_map[scenario_index];
					for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
					{

						fprintf(g_pFileLinkMOE, "VBA_%s_%s_%s,", assignment.g_DemandPeriodVector[tau].demand_period.c_str(), assignment.g_DTA_scenario_vector[scenario_no].scenario_name.c_str(), assignment.g_ModeTypeVector[at].mode_type.c_str());
					}

				}
			}

		}

		for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); ++tau)
		{
			for (int sii = 0; sii < assignment.g_DTA_scenario_vector.size(); sii++)
			{
				int scenario_index = assignment.g_DTA_scenario_vector[sii].scenario_index;

				int scenario_no = assignment.g_active_DTAscenario_map[scenario_index];
				for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
				{
					fprintf(g_pFileLinkMOE, "MEUVAB%s%s%s,", assignment.g_DemandPeriodVector[tau].demand_period.c_str(), assignment.g_DTA_scenario_vector[scenario_no].scenario_name.c_str(), assignment.g_ModeTypeVector[at].mode_type.c_str());
					fprintf(g_pFileLinkMOE, "capAB%s%s%s,", assignment.g_DemandPeriodVector[tau].demand_period.c_str(), assignment.g_DTA_scenario_vector[scenario_no].scenario_name.c_str(), assignment.g_ModeTypeVector[at].mode_type.c_str());
					fprintf(g_pFileLinkMOE, "PENAB%s%s%s,", assignment.g_DemandPeriodVector[tau].demand_period.c_str(), assignment.g_DTA_scenario_vector[scenario_no].scenario_name.c_str(), assignment.g_ModeTypeVector[at].mode_type.c_str());
					fprintf(g_pFileLinkMOE, "DOCAB%s%s%s,", assignment.g_DemandPeriodVector[tau].demand_period.c_str(), assignment.g_DTA_scenario_vector[scenario_no].scenario_name.c_str(), assignment.g_ModeTypeVector[at].mode_type.c_str());
					fprintf(g_pFileLinkMOE, "TTAB%s_%s_%s,", assignment.g_DemandPeriodVector[tau].demand_period.c_str(), assignment.g_DTA_scenario_vector[scenario_no].scenario_name.c_str(), assignment.g_ModeTypeVector[at].mode_type.c_str());
				}


			}
		}

		for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); ++tau)
		{
			for (int sii = 0; sii < assignment.g_DTA_scenario_vector.size(); sii++)
			{
				int scenario_index = assignment.g_DTA_scenario_vector[sii].scenario_index;
				int scenario_no = assignment.g_active_DTAscenario_map[scenario_index];
				for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
				{
					fprintf(g_pFileLinkMOE, "MEUVBA%s%s%s,", assignment.g_DemandPeriodVector[tau].demand_period.c_str(), assignment.g_DTA_scenario_vector[scenario_no].scenario_name.c_str(), assignment.g_ModeTypeVector[at].mode_type.c_str());
					fprintf(g_pFileLinkMOE, "capBA%s%s%s,", assignment.g_DemandPeriodVector[tau].demand_period.c_str(), assignment.g_DTA_scenario_vector[scenario_no].scenario_name.c_str(), assignment.g_ModeTypeVector[at].mode_type.c_str());
					fprintf(g_pFileLinkMOE, "PENBA%s%s%s,", assignment.g_DemandPeriodVector[tau].demand_period.c_str(), assignment.g_DTA_scenario_vector[scenario_no].scenario_name.c_str(), assignment.g_ModeTypeVector[at].mode_type.c_str());
					fprintf(g_pFileLinkMOE, "DOCBA%s%s%s,", assignment.g_DemandPeriodVector[tau].demand_period.c_str(), assignment.g_DTA_scenario_vector[scenario_no].scenario_name.c_str(), assignment.g_ModeTypeVector[at].mode_type.c_str());
					fprintf(g_pFileLinkMOE, "TTBA%s_%s_%s,", assignment.g_DemandPeriodVector[tau].demand_period.c_str(), assignment.g_DTA_scenario_vector[scenario_no].scenario_name.c_str(), assignment.g_ModeTypeVector[at].mode_type.c_str());
				}


			}
		}


	}

	fprintf(g_pFileLinkMOE, "\n");
	//		fprintf(g_pFileLinkMOE, "obs_count,upper_bound_type,dev,");
			// header

			//Initialization for all nodes
	for (int i = 0; i < g_link_vector.size(); ++i)
	{
		// virtual connectors
		if (g_link_vector[i].link_type_si[0] == -1)
			continue;

		if (g_link_vector[i].AB_flag == -1)  // skip BA links
			continue;

		string vdf_type_str;


		float total_vehicle_volume = 0;
		int scenario_code_count = 0;
		for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
		{
			total_vehicle_volume += g_link_vector[i].VDF_period[tau].volume_before_dtm;
			total_vehicle_volume += g_link_vector[i].total_volume_for_all_mode_types_per_period[tau] + g_link_vector[i].VDF_period[tau].preload;
			scenario_code_count += g_link_vector[i].VDF_period[tau].dtm_scenario_code.size();

		}

		//if (scenario_code_count == 0 && total_vehicle_volume < 0.001 && g_link_vector.size() > 20000)
		//	continue;




		float speed = g_link_vector[i].free_speed;  // default speed


		fprintf(g_pFileLinkMOE, "%d,%s,%s,%d,%d,",
			i+1,
			//end with obs_volume_diff
			g_link_vector[i].link_id.c_str(),
			assignment.g_LinkTypeMap[g_link_vector[i].link_type_si[0]].link_type_name.c_str(),
			g_node_vector[g_link_vector[i].from_node_seq_no].node_id,
			g_node_vector[g_link_vector[i].to_node_seq_no].node_id);

		if (g_link_vector[i].geometry.size() > 1)
		{
			fprintf(g_pFileLinkMOE, "\"%s\",",
				g_link_vector[i].geometry.c_str());
		}
		else
		{

			fprintf(g_pFileLinkMOE, "\"LINESTRING (");
			fprintf(g_pFileLinkMOE, "%f %f,", g_node_vector[g_link_vector[i].from_node_seq_no].x, g_node_vector[g_link_vector[i].from_node_seq_no].y);

			fprintf(g_pFileLinkMOE, "%f %f", g_node_vector[g_link_vector[i].to_node_seq_no].x, g_node_vector[g_link_vector[i].to_node_seq_no].y);

			fprintf(g_pFileLinkMOE, ")\",");
		}

		fprintf(g_pFileLinkMOE, "%f,%f,%f,%d,",

			g_link_vector[i].link_distance_km,
			g_link_vector[i].link_distance_mile,
			g_link_vector[i].free_flow_travel_time_in_min,
			g_link_vector[i].meso_link_id);

		fprintf(g_pFileLinkMOE, "%s,%s,%d,%d,%d,%d,",
			g_link_vector[i].tmc_code.c_str(),
			g_link_vector[i].tmc_corridor_name.c_str(),
			g_link_vector[i].tmc_corridor_id,
			g_link_vector[i].tmc_road_order,
			g_link_vector[i].tmc_road_sequence,
			g_link_vector[i].subarea_id);


		for (int sii = 0; sii < assignment.g_DTA_scenario_vector.size(); sii++)
		{
			int scenario_index = assignment.g_DTA_scenario_vector[sii].scenario_index;
			fprintf(g_pFileLinkMOE, "%.1f,", g_link_vector[i].number_of_lanes_si[scenario_index]);

		}

		for (int sii = 0; sii < assignment.g_DTA_scenario_vector.size(); sii++)
		{
			int scenario_index = assignment.g_DTA_scenario_vector[sii].scenario_index;
			fprintf(g_pFileLinkMOE, "%s,", assignment.g_LinkTypeMap[g_link_vector[i].link_type_si[scenario_index]].link_type_name.c_str());
		}

		for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
		{
			for (int sii = 0; sii < assignment.g_DTA_scenario_vector.size(); sii++)
			{
				int scenario_index = assignment.g_DTA_scenario_vector[sii].scenario_index;
				{
					int scenario_no = assignment.g_active_DTAscenario_map[scenario_index];
					for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
					{

						double AB_volume = g_link_vector[i].recorded_volume_per_period_per_at[tau][at][scenario_no];
						double BA_volume = 0;

						if (g_link_vector[i].BA_link_no >= 0)
						{
							BA_volume = g_link_vector[g_link_vector[i].BA_link_no].recorded_volume_per_period_per_at[tau][at][scenario_no];
						}

						double total_volume = AB_volume + BA_volume;

						fprintf(g_pFileLinkMOE, "%.1f,", total_volume);

					}
				}
			}

		}

		for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
		{
			for (int sii = 0; sii < assignment.g_DTA_scenario_vector.size(); sii++)
			{
				int scenario_index = assignment.g_DTA_scenario_vector[sii].scenario_index;
				{
					int scenario_no = assignment.g_active_DTAscenario_map[scenario_index];
					for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
					{

						fprintf(g_pFileLinkMOE, "%.1f,", g_link_vector[i].recorded_volume_per_period_per_at[tau][at][scenario_no]);

					}
				}
			}

		}

		for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
		{
			for (int sii = 0; sii < assignment.g_DTA_scenario_vector.size(); sii++)
			{
				int scenario_index = assignment.g_DTA_scenario_vector[sii].scenario_index;
				{
					int scenario_no = assignment.g_active_DTAscenario_map[scenario_index];

					for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
					{
						if (g_link_vector[i].BA_link_no >= 0)
							fprintf(g_pFileLinkMOE, "%.1f,", g_link_vector[g_link_vector[i].BA_link_no].recorded_volume_per_period_per_at[tau][at][scenario_no]);
						else
							fprintf(g_pFileLinkMOE, ",");

					}
				}
			}

		}

		for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
		{
			for (int sii = 0; sii < assignment.g_DTA_scenario_vector.size(); sii++)
			{
				int scenario_index = assignment.g_DTA_scenario_vector[sii].scenario_index;
				{
					int scenario_no = assignment.g_active_DTAscenario_map[scenario_index];
					for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
					{

						fprintf(g_pFileLinkMOE, "%.1f,", g_link_vector[i].recorded_MEU_per_period_per_at[tau][at][scenario_no]);
						fprintf(g_pFileLinkMOE, "%.1f,", g_link_vector[i].recorded_capacity_per_period_per_at[tau][at][scenario_no]);
						fprintf(g_pFileLinkMOE, "%.1f,", g_link_vector[i].penalty_si_at[tau][at][scenario_no]);

						
						fprintf(g_pFileLinkMOE, "%.4f,", g_link_vector[i].recorded_DOC_per_period_per_at[tau][at][scenario_no]);
						fprintf(g_pFileLinkMOE, "%.3f,", g_link_vector[i].recorded_TT_per_period_per_at[tau][at][scenario_no]);

					}
				}
			}

		}
		//

		{

			for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); ++tau)
			{
				for (int sii = 0; sii < assignment.g_DTA_scenario_vector.size(); sii++)
				{
					int scenario_index = assignment.g_DTA_scenario_vector[sii].scenario_index;
					{
						int scenario_no = assignment.g_active_DTAscenario_map[scenario_index];
						for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
						{
							if (g_link_vector[i].BA_link_no >= 0)
							{
								fprintf(g_pFileLinkMOE, "%.1f,", g_link_vector[g_link_vector[i].BA_link_no].recorded_MEU_per_period_per_at[tau][at][scenario_no]);
								fprintf(g_pFileLinkMOE, "%.3f,", g_link_vector[g_link_vector[i].BA_link_no].recorded_capacity_per_period_per_at[tau][at][scenario_no]);
								fprintf(g_pFileLinkMOE, "%.3f,", g_link_vector[g_link_vector[i].BA_link_no].penalty_si_at[tau][at][scenario_no]);
							
								fprintf(g_pFileLinkMOE, "%.4f,", g_link_vector[g_link_vector[i].BA_link_no].recorded_DOC_per_period_per_at[tau][at][scenario_no]);
								fprintf(g_pFileLinkMOE, "%.3f,", g_link_vector[g_link_vector[i].BA_link_no].recorded_TT_per_period_per_at[tau][at][scenario_no]);
							}
							else
								fprintf(g_pFileLinkMOE, ",,,,");


						}
					}
				}

			}
		}
		fprintf(g_pFileLinkMOE, "\n");
	}

	fclose(g_pFileLinkMOE);


}
void g_output_accessibility_result(Assignment& assignment)
{
	if (assignment.assignment_mode == lue)
	{
		dtalog.output() << "[WARNING] link based user equilibrum mode: no od_performance_summary.csv output.." << '\n';
		g_DTA_log_file << "[WARNING] link based user equilibrum mode: no od_performance_summary.csv output.." << '\n';
		return;
	}

	FILE* g_pFileODMOE = nullptr;
	fopen_ss(&g_pFileODMOE, "od_performance_summary.csv", "w");

	dtalog.output() << "[STATUS INFO] writing od_performance_summary.csv.." << '\n';
	g_DTA_log_file << "[STATUS INFO] writing od_performance_summary.csv.." << '\n';

	double path_time_vector[MAX_LINK_SIZE_IN_A_PATH];

	fopen_ss(&g_pFileODMOE, "od_performance_summary.csv", "w");

	if (!g_pFileODMOE)
	{
		dtalog.output() << "[ERROR] File od_performance_summary.csv cannot be opened." << '\n';
		g_DTA_log_file << "[ERROR] File od_performance_summary.csv cannot be opened." << '\n';
		return;
	}

	fprintf(g_pFileODMOE, "od_seq_id,o_zone_id,d_zone_id,o_sindex,d_sindex,o_district_id,d_district_id,mode_type,demand_period,connectivity_flag,s_x_coord,s_y_coord,t_x_coord,t_y_coord,");
	for (int sii = 0; sii < assignment.g_DTA_scenario_vector.size(); sii++)
	{
		int scenario_index = assignment.g_DTA_scenario_vector[sii].scenario_index;
		fprintf(g_pFileODMOE, "volume_s%d,", scenario_index);
	}
	for (int sii = 0; sii < assignment.g_DTA_scenario_vector.size(); sii++)
	{
		int scenario_index = assignment.g_DTA_scenario_vector[sii].scenario_index;
		fprintf(g_pFileODMOE, "avg_tt_s%d,", scenario_index);
	}
	for (int sii = 0; sii < assignment.g_DTA_scenario_vector.size(); sii++)
	{
		int scenario_index = assignment.g_DTA_scenario_vector[sii].scenario_index;
		fprintf(g_pFileODMOE, "avg_distance_km_s%d,", scenario_index);
	}	for (int sii = 0; sii < assignment.g_DTA_scenario_vector.size(); sii++)
	{
		int scenario_index = assignment.g_DTA_scenario_vector[sii].scenario_index;
		fprintf(g_pFileODMOE, "avg_cost_s%d,", scenario_index);
	}



	fprintf(g_pFileODMOE, "\n");
	int count = 1;

	clock_t start_t, end_t;
	start_t = clock();
	clock_t iteration_t;


	int mode_type_size = assignment.g_ModeTypeVector.size();
	int zone_size = g_zone_vector.size();
	int demand_period_size = assignment.g_DemandPeriodVector.size();

	CColumnVector* p_column_pool;

	float path_toll = 0;
	float path_distance = 0;
	float path_travel_time = 0;
	float path_travel_time_without_access_link = 0;
	float path_delay = 0;
	float path_FF_travel_time = 0;
	float time_stamp = 0;

	bool column_pool_ready_flag = false;

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
				g_link_vector[i].background_total_volume_for_all_mode_types_per_period[tau] = 0;
			}

			if (g_node_vector[g_link_vector[i].from_node_seq_no].subarea_id >= 1 && g_node_vector[g_link_vector[i].to_node_seq_no].node_id >= 1)
			{

				b_subarea_mode = true;
			}

		}


		/// <summary>  screening the path flow pattern
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
							column_pool_ready_flag = true;

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
								if (subarea_output_flag == 0)
								{
									it->second.subarea_output_flag = 0;  // disable the output of this column into route_assignment.csv

									for (int nl = 0; nl < it->second.m_link_size; ++nl)  // arc a
									{
										int link_seq_no = it->second.path_link_vector[nl];
										g_link_vector[link_seq_no].background_total_volume_for_all_mode_types_per_period[tau] += it->second.path_volume;
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


	} // end of path flow pattern screening
	dtalog.output() << "[STATUS INFO] writing data for " << zone_size << "  zones " << '\n';
	g_DTA_log_file << "[STATUS INFO] writing data for " << zone_size << "  zones " << '\n';

	int number_of_connected_OD_pairs = 0;

	std::map<string, bool> l_origin_destination_map;
	std::map<string, bool> l_origin_destination_disconnected_map;

	for (int orig = 0; orig < zone_size; ++orig)
	{
		int from_zone_sindex = g_zone_vector[orig].sindex;
		if (from_zone_sindex == -1)
			continue;

		if (g_zone_vector[orig].b_shortest_path_computing_flag == false)
			continue;

		if (g_zone_vector[orig].zone_id % 100 == 0)
		{
			dtalog.output() << "[DATA INFO] o zone id =  " << g_zone_vector[orig].zone_id << '\n';
			g_DTA_log_file << "[DATA INFO] o zone id =  " << g_zone_vector[orig].zone_id << '\n';
		}
		for (int dest = 0; dest < zone_size; ++dest)
		{
			int to_zone_sindex = g_zone_vector[dest].sindex;

			if (to_zone_sindex == -1)
				continue;

			if (g_zone_vector[dest].b_shortest_path_computing_flag == false)
				continue;

			for (int at = 0; at < mode_type_size; ++at)
			{
				for (int tau = 0; tau < demand_period_size; ++tau)
				{
					p_column_pool = &(assignment.g_column_pool[from_zone_sindex][to_zone_sindex][at][tau]);
					if (p_column_pool->od_volume[assignment.active_scenario_index] > 0)
					{

						int information_type = 0;

						if (assignment.zone_seq_no_2_info_mapping.find(orig) != assignment.zone_seq_no_2_info_mapping.end())
						{
							information_type = 1;
						}

						time_stamp = (assignment.g_DemandPeriodVector[tau].starting_time_slot_no + assignment.g_DemandPeriodVector[tau].ending_time_slot_no) / 2.0 * MIN_PER_TIMESLOT;

						string od_string;

						od_string = std::to_string(g_zone_vector[orig].zone_id) + "->" + std::to_string(g_zone_vector[dest].zone_id);

						if (p_column_pool->path_node_sequence_map.size() >= 1)
						{
							l_origin_destination_map[od_string] = true;
						}
						else
						{
							if (orig != dest)
							{
								od_string = od_string + ":at = " + assignment.g_ModeTypeVector[at].mode_type.c_str();
								od_string = od_string + ":dp = " + assignment.g_DemandPeriodVector[tau].demand_period.c_str();
								l_origin_destination_disconnected_map[od_string] = false;
							}
						}



							fprintf(g_pFileODMOE, "%d,%d,%d,%d,%d,%d,%d,",
								count,
								g_zone_vector[orig].zone_id,
								g_zone_vector[dest].zone_id,
								g_zone_vector[orig].sindex,
								g_zone_vector[dest].sindex,
								assignment.g_zone_seq_no_to_analysis_distrct_id_mapping[orig],
								assignment.g_zone_seq_no_to_analysis_distrct_id_mapping[dest]

							);
							fprintf(g_pFileODMOE, "%s,%s,0,",
								assignment.g_ModeTypeVector[at].mode_type.c_str(),
								assignment.g_DemandPeriodVector[tau].demand_period.c_str());

							fprintf(g_pFileODMOE, "%f,%f,%f,%f,",
								g_zone_vector[orig].cell_x,
								g_zone_vector[orig].cell_y,
								g_zone_vector[dest].cell_x,
								g_zone_vector[dest].cell_y);

							for (int sii = 0; sii < assignment.g_DTA_scenario_vector.size(); sii++)
							{
								int scenario_index = assignment.g_DTA_scenario_vector[sii].scenario_index;
								fprintf(g_pFileODMOE, "%4.2f,", p_column_pool->od_volume[scenario_index]);
							}
							for (int sii = 0; sii < assignment.g_DTA_scenario_vector.size(); sii++)
							{
								int scenario_index = assignment.g_DTA_scenario_vector[sii].scenario_index;
								fprintf(g_pFileODMOE, "%4.2f,", p_column_pool->avg_travel_time[scenario_index]);
							}
							for (int sii = 0; sii < assignment.g_DTA_scenario_vector.size(); sii++)
							{
								int scenario_index = assignment.g_DTA_scenario_vector[sii].scenario_index;
								fprintf(g_pFileODMOE, "%4.2f,", p_column_pool->avg_distance[scenario_index]);
							}
							for (int sii = 0; sii < assignment.g_DTA_scenario_vector.size(); sii++)
							{
								int scenario_index = assignment.g_DTA_scenario_vector[sii].scenario_index;
								fprintf(g_pFileODMOE, "%4.2f,", p_column_pool->avg_cost[scenario_index]);
							}

							fprintf(g_pFileODMOE, "\n");

							continue;
						}

						
						}
					}
				}
	
	
	}
	fclose(g_pFileODMOE);

	if (column_pool_ready_flag)
	{
		assignment.summary_file << "     Check OD connectivity and accessibility in od_performance_summary.csv" << '\n';
		assignment.summary_file << ", # of connected OD pairs =, " << l_origin_destination_map.size() << '\n';
		assignment.summary_file << ", # of OD/mode_type/demand_type columns without paths =, " << l_origin_destination_disconnected_map.size() << '\n';
//		dtalog.output() << ", # of connected OD pairs = " << l_origin_destination_map.size() << '\n';  // to be checked, Xuesong, 2023
//		g_DTA_log_file << ", # of connected OD pairs = " << l_origin_destination_map.size() << '\n';  // to be checked, Xuesong, 2023

	}
	//if (l_origin_destination_map.size() == 0)
	//{
	//	g_OutputModelFiles(10); // label cost tree
	//	dtalog.output() << "Please check the connectivity of OD pairs and in network and field allow_uses in link.csv." << '\n';
	//	g_DTA_log_file << "Please check the connectivity of OD pairs and in network and field allow_uses in link.csv." << '\n';
	//	dtalog.output() << "Please check the model_shortest_path_tree.csv file." << '\n';
	//	g_DTA_log_file << "Please check the model_shortest_path_tree.csv file." << '\n';
	//	//			g_program_stop();
	//}


}

void g_output_dynamic_link_performance(Assignment& assignment, int output_mode = 1)
{
	//dtalog.output() << "writing dynamic_waiting_performance_profile.csv.." << '\n';
	//g_DTA_log_file << "writing dynamic_waiting_performance_profile.csv.." << '\n';

	//int b_debug_detail_flag = 0;
	//FILE* g_pFileLinkMOE = nullptr;

	//string file_name = "dynamic_link_waiting_time_profile.csv";

	// fopen_ss(&g_pFileLinkMOE, file_name.c_str(), "w");

	//if (!g_pFileLinkMOE)
	//{
	//    dtalog.output() << "File " << file_name.c_str() << " cannot be opened." << '\n';
	//    g_DTA_log_file << "File " << file_name.c_str() << " cannot be opened." << '\n';
	//    g_program_stop();
	//}

	//    // Option 2: BPR-X function
	//    fprintf(g_pFileLinkMOE, "link_id,scenario_flag,from_node_id,to_node_id,geometry,");
	//    for (int t = 1; t < assignment.g_number_of_intervals_in_min; ++t)
	//        fprintf(g_pFileLinkMOE, "WT%d,",t);

	//    for (int t = 1; t < assignment.g_number_of_intervals_in_min; ++t)
	//        fprintf(g_pFileLinkMOE, "QL%d,", t);

	//    fprintf(g_pFileLinkMOE, "notes\n");


	//    //Initialization for all nodes
	//    for (int i = 0; i < g_link_vector.size(); ++i)
	//    {
	//        // virtual connectors
	//        if (g_link_vector[i].link_type == -1)
	//            continue;

	//                fprintf(g_pFileLinkMOE, "%s,%d,%d,%d,\"%s\",",
	//                    g_link_vector[i].link_id.c_str(),
	//                    g_link_vector[i].capacity_reduction_map.size(),
	//                    g_node_vector[g_link_vector[i].from_node_seq_no].node_id,
	//                    g_node_vector[g_link_vector[i].to_node_seq_no].node_id,
	//                    g_link_vector[i].geometry.c_str());

	//                // first loop for time t
	//                for (int t = 1; t < assignment.g_number_of_intervals_in_min; ++t)
	//                {
	//                    float waiting_time_in_min = 0;

	//                       int timestamp_in_min = assignment.g_LoadingStartTimeInMin + t;

	//                        //if (g_link_vector[i].RT_travel_time_map.find(timestamp_in_min) != g_link_vector[i].RT_travel_time_map.end())
	//                        //{
	//                        //    waiting_time_in_min = g_link_vector[i].RT_travel_time_map[timestamp_in_min];
	//                        //}

	//                        if(timestamp_in_min>0.05)
	//                            fprintf(g_pFileLinkMOE, "%.2f,", waiting_time_in_min);
	//                        else
	//                            fprintf(g_pFileLinkMOE, ",");

	//                }

	//                // first loop for time t
	//                for (int t = 1; t < assignment.g_number_of_intervals_in_min; ++t)
	//                {
	//                    float QL = (assignment.m_LinkCumulativeArrivalVector[i][t] - assignment.m_LinkCumulativeDepartureVector[i][t]);

	//                        fprintf(g_pFileLinkMOE, "%.1f,", QL);
	//                }
	//                fprintf(g_pFileLinkMOE, "\n");
	//    }  // for each link l
	//    fclose(g_pFileLinkMOE);
}

void g_output_TD_link_performance(Assignment& assignment, int output_mode = 1)
{
	dtalog.output() << "[STATUS INFO] writing td_link_performance.csv.." << '\n';
	g_DTA_log_file << "[STATUS INFO] writing td_link_performance.csv.." << '\n';
	dtalog.output() << "[STATUS INFO] writing td_link_performance.csv.." << '\n';
	g_DTA_log_file << "[STATUS INFO] writing td_link_performance.csv.." << '\n';

	int b_debug_detail_flag = 0;
	FILE* g_pFileLinkMOE = nullptr;


	char TD_link_performance_file_name[50];

	int scenario_no = assignment.g_active_DTAscenario_map[assignment.active_scenario_index];
	sprintf(TD_link_performance_file_name, "td_link_performance%d_%s.csv", assignment.active_scenario_index, assignment.g_DTA_scenario_vector[scenario_no].scenario_name.c_str());

	fopen_ss(&g_pFileLinkMOE, TD_link_performance_file_name, "w");

	if (!g_pFileLinkMOE)
	{
		dtalog.output() << "[ERROR] File " << TD_link_performance_file_name << " cannot be opened." << '\n';
		g_DTA_log_file << "[ERROR] File " << TD_link_performance_file_name << " cannot be opened." << '\n';
		return;
	}
	else
	{

		// Option 2: BPR-X function
		fprintf(g_pFileLinkMOE, "link_id,tmc_corridor_name,link_type_name,from_node_id,to_node_id,meso_link_id,from_cell_code,lanes,length,free_speed_kmph,free_speed_mph,FFTT,time_period,inflow_volume,volume,CA,CD,density,queue,queue_ratio,discharge_cap,TD_free_flow_travel_time,waiting_time_in_sec,speed,geometry,");
		fprintf(g_pFileLinkMOE, "notes\n");

		int sampling_time_interval = assignment.td_link_performance_sampling_interval_in_min; // min by min

		if (sampling_time_interval <= 0)
		{ // apply automated configuration
			if (g_link_vector.size() > 5000)
				sampling_time_interval = 15;

			if (g_link_vector.size() > 10000)
				sampling_time_interval = 30;

			if (g_link_vector.size() > 50000)
				sampling_time_interval = 60;

			if (g_link_vector.size() > 500000)
				sampling_time_interval = 120;
		}
		else
		{
			//with user input

		}

		//Initialization for all nodes
		for (int i = 0; i < g_link_vector.size(); ++i)
		{
			// virtual connectors
			if (g_link_vector[i].link_type_si[0] == -1)
				continue;

			// first loop for time t
			for (int t = 1; t < assignment.g_number_of_intervals_in_min; ++t)
			{


				if (t % (sampling_time_interval) == 0)
				{
					int time_in_min = t;  //relative time

					float inflow_volume = 0;
					float outflow_volume = 0;
					float queue = 0;
					float waiting_time_in_sec = 0;
					int arrival_rate = 0;
					float avg_waiting_time_in_sec = 0;

					float travel_time = (float)(g_link_vector[i].free_flow_travel_time_in_min + avg_waiting_time_in_sec / 60.0);
					float speed = g_link_vector[i].link_distance_VDF / (g_link_vector[i].free_flow_travel_time_in_min / 60.0);
					float virtual_arrival = 0;

					float discharge_rate = g_link_vector[i].lane_capacity * g_link_vector[i].number_of_lanes_si[0] * sampling_time_interval / 60.0;

					if (time_in_min >= 1)
					{
						inflow_volume = assignment.m_LinkCumulativeArrivalVector[i][t] - assignment.m_LinkCumulativeArrivalVector[i][t - sampling_time_interval];
						outflow_volume = assignment.m_LinkCumulativeDepartureVector[i][t] - assignment.m_LinkCumulativeDepartureVector[i][t - sampling_time_interval];

						queue = assignment.m_LinkCumulativeArrivalVector[i][t] - assignment.m_LinkCumulativeDepartureVector[i][t];
						//							waiting_time_in_min = queue / (max(1, volume));

						float waiting_time_count = 0;

						waiting_time_in_sec = assignment.m_link_TD_waiting_time[i][t] * number_of_seconds_per_interval;

						arrival_rate = assignment.m_LinkCumulativeArrivalVector[i][t] - assignment.m_LinkCumulativeArrivalVector[i][t - sampling_time_interval];
						avg_waiting_time_in_sec = waiting_time_in_sec / max(1, arrival_rate);

						travel_time = (float)(g_link_vector[i].free_flow_travel_time_in_min + avg_waiting_time_in_sec / 60.0);
						speed = g_link_vector[i].link_distance_VDF / (max(0.00001, travel_time / 60.0));
					}

					if (speed >= 1000)
					{
						int i_debug = 1;
					}
					float density = (assignment.m_LinkCumulativeArrivalVector[i][t] - assignment.m_LinkCumulativeDepartureVector[i][t]) / (g_link_vector[i].link_distance_VDF * g_link_vector[i].number_of_lanes_si[assignment.active_scenario_index]);
					float queue_ratio = (assignment.m_LinkCumulativeArrivalVector[i][t] - assignment.m_LinkCumulativeDepartureVector[i][t]) / (g_link_vector[i].link_distance_VDF * g_link_vector[i].number_of_lanes_si[assignment.active_scenario_index] * g_link_vector[i].kjam);
					if (queue_ratio > 1)
						queue_ratio = 1;
					if (output_mode == 0)
					{
						if (assignment.m_LinkCumulativeArrivalVector[i][t] < 1000)
							continue; //skip
					}

					fprintf(g_pFileLinkMOE, "%s,%s,%s,%d,%d,%d,%s,%.1f,%.3f,%.3f,%.3f,%.3f,%s_%s,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,\"%s\",",
						g_link_vector[i].link_id.c_str(),
						g_link_vector[i].tmc_corridor_name.c_str(),
						g_link_vector[i].link_type_name.c_str(),

						g_node_vector[g_link_vector[i].from_node_seq_no].node_id,
						g_node_vector[g_link_vector[i].to_node_seq_no].node_id,
						g_link_vector[i].meso_link_id,
						g_node_vector[g_link_vector[i].from_node_seq_no].cell_str.c_str(),
						g_link_vector[i].number_of_lanes_si[0],
						g_link_vector[i].length_in_meter,
						g_link_vector[i].free_speed,
						g_link_vector[i].free_speed / 1.609,
						g_link_vector[i].free_flow_travel_time_in_min,

						g_time_coding(assignment.g_LoadingStartTimeInMin + time_in_min - sampling_time_interval).c_str(),
						g_time_coding(assignment.g_LoadingStartTimeInMin + time_in_min).c_str(),
						inflow_volume,
						outflow_volume,
						assignment.m_LinkCumulativeArrivalVector[i][t],
						assignment.m_LinkCumulativeDepartureVector[i][t],
						density,
						queue,
						queue_ratio,
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

void g_output_dynamic_link_state(Assignment& assignment, int output_mode = 1)
{
	dtalog.output() << "[STATUS INFO] writing log_dynamic_link_state.txt.." << '\n';
	g_DTA_log_file << "[STATUS INFO] writing log_dynamic_link_state.txt.." << '\n';

	int b_debug_detail_flag = 0;
	FILE* g_pFileLinkMOE = nullptr;

	string file_name = "log_dynamic_link_state.txt";

	fopen_ss(&g_pFileLinkMOE, file_name.c_str(), "w");

	if (!g_pFileLinkMOE)
	{
		dtalog.output() << "[ERROR] File " << file_name.c_str() << " cannot be opened." << '\n';
		g_DTA_log_file << "[ERROR] File " << file_name.c_str() << " cannot be opened." << '\n';
		return;
	}
	else
	{

		// Option 2: BPR-X function
		fprintf(g_pFileLinkMOE, "link_id,from_node_id,to_node_id,time_period,duration_in_sec,state,state_code\n");


		//Initialization for all nodes
		for (unsigned li = 0; li < g_link_vector.size(); ++li)
		{

			if (g_link_vector[li].timing_arc_flag)  // only output the capaicty for signalized data
			{
				// reset for signalized links (not freeway links as type code != 'f' for the case of freeway workzones)
				// only for the loading period

				int t = 0;
				int last_t = t;
				if (assignment.m_LinkOutFlowState == NULL)
					break;

				int current_state = assignment.m_LinkOutFlowState[li][t];
				if (current_state == 0)  //close in the beginning
					current_state = -1; // change the initial state to be able to be record the change below

				while (t < assignment.g_number_of_loading_intervals_in_sec - 1)
				{
					int next_state = assignment.m_LinkOutFlowState[li][t + 1];

					bool print_out_flag = false;

					if (print_out_flag == false && next_state == current_state && t < assignment.g_number_of_loading_intervals_in_sec - 2)
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

						if (state_str.size() > 0)  // with data to output
						{
							fprintf(g_pFileLinkMOE, "%s,%d,%d,%s_%s,%d,%d,%s\n",
								g_link_vector[li].link_id.c_str(),
								g_node_vector[g_link_vector[li].from_node_seq_no].node_id,
								g_node_vector[g_link_vector[li].to_node_seq_no].node_id,
								g_time_coding(assignment.g_LoadingStartTimeInMin + last_t / 60.0).c_str(),
								g_time_coding(assignment.g_LoadingStartTimeInMin + (t + 1) / 60.0).c_str(),
								t + 1 - last_t,
								current_state,
								state_str.c_str());
						}

						last_t = t + 1;
						current_state = assignment.m_LinkOutFlowState[li][t + 1];

						if (t + 1 == assignment.g_number_of_loading_intervals_in_sec - 2)
						{
							//boundary condition anyway
							current_state = -1;
						}

					}
					t++;
				}


			}
			if (g_link_vector[li].m_link_pedefined_capacity_map_in_sec.size() > 0)  // only output the capaicty for signalized data
			{
				// reset for signalized links (not freeway links as type code != 'f' for the case of freeway workzones)
				// only for the loading period

				fprintf(g_pFileLinkMOE, "%s,%d,%d,%s_%s,%d,0,w\n",
					g_link_vector[li].link_id.c_str(),
					g_node_vector[g_link_vector[li].from_node_seq_no].node_id,
					g_node_vector[g_link_vector[li].to_node_seq_no].node_id,
					g_time_coding(g_link_vector[li].global_minute_capacity_reduction_start).c_str(),
					g_time_coding(g_link_vector[li].global_minute_capacity_reduction_end).c_str(),
					g_link_vector[li].global_minute_capacity_reduction_end - g_link_vector[li].global_minute_capacity_reduction_start);

			}


		}



		fclose(g_pFileLinkMOE);

	}
}

//void g_output_simulation_agents(Assignment& assignment)
//{
//    if (assignment.assignment_mode == lue || assignment.trajectory_output_count == 0)  //LUE
//    {
//        FILE* g_pFilePathMOE = nullptr;
//        fopen_ss(&g_pFilePathMOE, "agent.csv", "w");
//        fclose(g_pFilePathMOE);
//    }
//    else if (assignment.assignment_mode >= 1)  //UE mode, or ODME, DTA
//    {
//        dtalog.output() << "writing agent.csv.." << '\n';
//        g_DTA_log_file << "writing agent.csv.." << '\n';
//
//        double path_time_vector[MAX_LINK_SIZE_IN_A_PATH];
//        FILE* g_pFileAgent = nullptr;
//        fopen_ss(&g_pFileAgent, "agent.csv", "w");
//
//        if (!g_pFileAgent)
//        {
//            dtalog.output() << "File agent.csv cannot be opened." << '\n';
//            g_DTA_log_file << "File agent.csv cannot be opened." << '\n';
//            g_program_stop();
//        }
//
//        fprintf(g_pFileAgent, "agent_id,o_zone_id,d_zone_id,OD_index,path_id,#_of_links,diverted_flag,mode_type,demand_period,volume,toll,travel_time,distance_km,speed,departure_time_in_min,arrival_time_in_min,departure_time_slot_no,\n");
//
//        int count = 1;
//
//        clock_t start_t, end_t;
//        start_t = clock();
//        clock_t iteration_t;
//
//        int mode_type_size = assignment.g_ModeTypeVector.size();
//        int zone_size = g_zone_vector.size();
//        int demand_period_size = assignment.g_DemandPeriodVector.size();
//
//        CColumnVector* p_column_pool;
//
//        float path_toll = 0;
//        float path_distance = 0;
//        float path_travel_time = 0;
//        float path_travel_time_without_access_link = 0;
//        float time_stamp = 0;
//
//        if (assignment.trajectory_sampling_rate < 0.01)
//            assignment.trajectory_sampling_rate = 0.01;
//        int sampling_step = 100 / int(100 * assignment.trajectory_sampling_rate + 0.5);
//
//        std::map<int, CColumnPath>::iterator it, it_begin, it_end;
//
//        dtalog.output() << "writing data for " << zone_size << "  zones " << '\n';
//        g_DTA_log_file << "writing data for " << zone_size << "  zones " << '\n';
//
//        for (int orig = 0; orig < zone_size; ++orig)
//        {
//            if (g_zone_vector[orig].zone_id % 100 == 0)
//                dtalog.output() << "o zone id =  " << g_zone_vector[orig].zone_id << '\n';
//                g_DTA_log_file << "o zone id =  " << g_zone_vector[orig].zone_id << '\n';
//
//            for (int at = 0; at < mode_type_size; ++at)
//            {
//                for (int dest = 0; dest < zone_size; ++dest)
//                {
//                    for (int tau = 0; tau < demand_period_size; ++tau)
//                    {
//                        p_column_pool = &(assignment.g_column_pool[orig][dest][at][tau]);
//                        if (p_column_pool->od_volume[assignment.active_scenario_index] > 0)
//                        {
//
//                            // scan through the map with different node sum for different continuous paths
//                            it_begin = p_column_pool->path_node_sequence_map.begin();
//                            it_end = p_column_pool->path_node_sequence_map.end();
//
//                            for (it = it_begin; it != it_end; ++it)
//                            {
//                                if (count % 100000 == 0)
//                                {
//                                    end_t = clock();
//                                    iteration_t = end_t - start_t;
//                                    dtalog.output() << "writing " << count / 1000 << "K agents with CPU time " << iteration_t / 1000.0 << " s" << '\n';
//                                    g_DTA_log_file << "writing " << count / 1000 << "K agents with CPU time " << iteration_t / 1000.0 << " s" << '\n';
//                                }
//
//                                if (count % sampling_step != 0)
//                                    continue;
//
//
//                                path_toll = 0;
//                                path_distance = 0;
//                                path_travel_time = 0;
//                                path_travel_time_without_access_link = 0;
//                                path_time_vector[0] = time_stamp;
//
//
//
//                                // assignment_mode = 1, path flow mode
//                                {
//                                    // assignment_mode = 2, simulated agent flow mode //DTA simulation
//
//                                    for (int vi = 0; vi < it->second.agent_simu_id_vector.size(); ++vi)
//                                    {
//                                        int agent_simu_id = it->second.agent_simu_id_vector[vi];
//                                        CAgent_Simu* pAgentSimu = g_agent_simu_vector[agent_simu_id];
//                                        time_stamp = assignment.g_LoadingStartTimeInMin + pAgentSimu->departure_time_in_min;
//
//                                        float departure_time_in_min = time_stamp;
//                                        float arrival_time_in_min = assignment.g_LoadingStartTimeInMin + pAgentSimu->arrival_time_in_min;
//                                        int departure_time_in_slot_no = time_stamp / MIN_PER_TIMESLOT;
//                                        float speed = pAgentSimu->path_distance / max(0.001, pAgentSimu->path_travel_time_in_min) * 60;
//
//                                        if (it->second.m_link_size >= MAX_LINK_SIZE_IN_A_PATH -1)
//                                        {
//                                            dtalog.output() << "error: it->second.m_link_size >= MAX_LINK_SIZE_IN_A_PATH-1" << '\n';
//                                            g_DTA_log_file << "error: it->second.m_link_size >= MAX_LINK_SIZE_IN_A_PATH-1" << '\n';
//                                            g_program_stop();
//                                        }
//
//                                        fprintf(g_pFileAgent, "%d,%d,%d,%d->%d,%d,%d,%d,%s,%s,1,%.1f,%.4f,%.4f,%.4f,%.4f,%.4f,%d",
//                                            pAgentSimu->agent_id,
//                                            g_zone_vector[orig].zone_id,
//                                            g_zone_vector[dest].zone_id,
//                                            g_zone_vector[orig].zone_id,
//                                            g_zone_vector[dest].zone_id,
//                                            it->second.path_seq_no,
//                                            it->second.m_link_size,
//                                            pAgentSimu->diverted_flag,
//                                            assignment.g_ModeTypeVector[at].mode_type.c_str(),
//                                            assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
//                                            path_toll,
//                                            pAgentSimu->path_travel_time_in_min,
//                                            pAgentSimu->path_distance, speed,
//                                            departure_time_in_min,
//                                            arrival_time_in_min,
//                                            departure_time_in_slot_no);
//
//                                        count++;
//
//                                        fprintf(g_pFileAgent, "\n");
//                                    }
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//        }
//        fclose(g_pFileAgent);
//    }
//}

void g_output_agent_csv(Assignment& assignment)
{
	int total_impacted_vehicle_count = 0;
	int total_diverted_vehicle_count = 0;
	int total_non_diverted_vehicle_count = 0;

	double total_diverted_vehicle_travel_time = 0;
	double total_diverted_vehicle_travel_distance = 0;

	int total_diverted_DMS_vehicle_count = 0;
	double total_diverted_DMS_vehicle_travel_time = 0;
	double total_diverted_DMS_vehicle_travel_distance = 0;

	double total_non_diverted_vehicle_travel_time = 0;
	double total_non_diverted_vehicle_travel_distance = 0;

	if (assignment.assignment_mode == lue)  //LUE
	{
		FILE* g_pFilePathMOE = nullptr;
		fopen_ss(&g_pFilePathMOE, "agent.csv", "w");
		fclose(g_pFilePathMOE);
	}
	else if (assignment.assignment_mode >= 1)  //UE mode, or ODME, DTA
	{

		char agent_file_name[50];

		int scenario_no = assignment.g_active_DTAscenario_map[assignment.active_scenario_index];
		sprintf(agent_file_name, "agent_%d_%s.csv", assignment.active_scenario_index, assignment.g_DTA_scenario_vector[scenario_no].scenario_name.c_str());

		dtalog.output() << "[STATUS INFO] writing " << agent_file_name << '\n';
		g_DTA_log_file << "[STATUS INFO] writing " << agent_file_name << '\n';
		assignment.summary_file << ",summary by multi-modal agents,mode_type,#_agents,avg_distance,avg_travel_time_in_min,avg_free_speed," << '\n';

		double path_time_vector[MAX_LINK_SIZE_IN_A_PATH];
		FILE* g_pFilePathMOE = nullptr;
		fopen_ss(&g_pFilePathMOE, agent_file_name, "w");

		if (!g_pFilePathMOE)
		{
			dtalog.output() << "[ERROR] File " << agent_file_name  << "cannot be opened." << '\n';
			g_DTA_log_file << "[ERROR] File " << agent_file_name  << "cannot be opened." << '\n';
			return;
		}

		fprintf(g_pFilePathMOE, "first_column,agent_id,o_zone_id,d_zone_id,OD_key,route_seq_id,path_no,impacted_flag,info_receiving_flag,diverted_flag,mode_type_no,mode_type,demand_period,volume,toll,departure_time,dt_hhmm,departure_time_in_1min,departure_time_in_5_min,travel_time,distance_mile,speed_mph,waiting_time,max_link_waiting_time,max_wait_link,\n");

		int count = 1;

		clock_t start_t, end_t;
		start_t = clock();
		clock_t iteration_t;

		int mode_type_size = assignment.g_ModeTypeVector.size();
		int zone_size = g_zone_vector.size();
		int demand_period_size = assignment.g_DemandPeriodVector.size();

		CColumnVector* p_column_pool;

		float path_toll = 0;
		float path_distance = 0;
		float path_travel_time = 0;
		float speed_mph = 0;
		float time_stamp = 0;

		if (assignment.trajectory_sampling_rate < 0.01)
			assignment.trajectory_sampling_rate = 0.01;
		int sampling_step = 100 / int(100 * assignment.trajectory_sampling_rate + 0.5);

		std::map<int, CColumnPath>::iterator it, it_begin, it_end;

		dtalog.output() << "[STATUS INFO] writing data for " << zone_size << "  zones " << '\n';
		g_DTA_log_file << "[STATUS INFO] writing data for " << zone_size << "  zones " << '\n';



		for (int at = 0; at < mode_type_size; ++at)
		{
			int agent_count = 0;
			double total_agent_distance = 0;
			double total_agent_travel_time = 0;

			for (int orig = 0; orig < zone_size; ++orig)
			{
				if (g_zone_vector[orig].zone_id % 100 == 0)
				{
					dtalog.output() << "[DATA INFO] o zone id =  " << g_zone_vector[orig].zone_id << '\n';
					g_DTA_log_file << "[DATA INFO] o zone id =  " << g_zone_vector[orig].zone_id << '\n';
				}
				int from_zone_sindex = g_zone_vector[orig].sindex;
				if (from_zone_sindex == -1)
					continue;


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

								bool route_impacted_flag = false;

								{

									for (int vi = 0; vi < it->second.agent_simu_id_vector.size(); ++vi)
									{

										int agent_simu_id = it->second.agent_simu_id_vector[vi];
										CAgent_Simu* pAgentSimu = g_agent_simu_vector[agent_simu_id];

										if (pAgentSimu->impacted_flag == 1)
										{
											route_impacted_flag = true;
										}
									}

								}

								if (count % 100000 == 0)
								{
									end_t = clock();
									iteration_t = end_t - start_t;
									dtalog.output() << "[STATUS INFO] writing " << count / 1000 << "K agents with CPU time " << iteration_t / 1000.0 << " s" << '\n';
									g_DTA_log_file << "[STATUS INFO] writing " << count / 1000 << "K agents with CPU time " << iteration_t / 1000.0 << " s" << '\n';
								}

								if (count % sampling_step != 0)
									continue;

								if (count > assignment.trajectory_output_count && assignment.trajectory_output_count >= 1)  // if trajectory_output_count ==-1, this constraint does not apply
									continue;

								path_toll = 0;
								path_distance = 0;
								path_travel_time = 0;
								path_time_vector[0] = time_stamp;

								for (int nl = 0; nl < it->second.m_link_size; ++nl)  // arc a
								{
									int link_seq_no = it->second.path_link_vector[nl];
									path_toll += g_link_vector[link_seq_no].VDF_period[tau].toll[at][assignment.active_scenario_index];
									path_distance += g_link_vector[link_seq_no].link_distance_VDF;
									float link_travel_time = g_link_vector[link_seq_no].link_avg_travel_time_per_period[tau][at];

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

									agent_count += it->second.agent_simu_id_vector.size();
									for (int vi = 0; vi < it->second.agent_simu_id_vector.size(); ++vi)
									{


										int agent_simu_id = it->second.agent_simu_id_vector[vi];
										CAgent_Simu* pAgentSimu = g_agent_simu_vector[agent_simu_id];
										if (pAgentSimu->agent_id == 81)
										{
											int idebug = 1;
										}
										if (assignment.trajectory_diversion_only == 1 && pAgentSimu->diverted_flag == 0)  // diversion flag only, then we skip the non-diversion path
											continue;

										float vehicle_departure_time = assignment.g_LoadingStartTimeInMin + pAgentSimu->departure_time_in_min;

										float waiting_time_in_min = pAgentSimu->waiting_time_in_min;
										float max_link_waiting_time_in_min = pAgentSimu->max_link_waiting_time_in_min;

										int from_node_id_max_wait = -1;
										int to_node_id_max_wait = -1;

										if (pAgentSimu->max_waiting_time_link_no >= 0)
										{
											from_node_id_max_wait = g_node_vector[g_link_vector[pAgentSimu->max_waiting_time_link_no].from_node_seq_no].node_id;
											to_node_id_max_wait = g_node_vector[g_link_vector[pAgentSimu->max_waiting_time_link_no].to_node_seq_no].node_id;
										}

										time_stamp = assignment.g_LoadingStartTimeInMin + pAgentSimu->departure_time_in_min;
										for (int nt = 0 + virtual_link_delta; nt < pAgentSimu->path_link_seq_no_vector.size() + 1 - virtual_link_delta; ++nt)
										{
											float time_in_min = 0;

											if (nt < pAgentSimu->path_link_seq_no_vector.size() - virtual_link_delta)
												time_in_min = assignment.g_LoadingStartTimeInMin + pAgentSimu->m_veh_link_arrival_time_in_simu_interval[nt] * number_of_seconds_per_interval / 60.0;
											else
												time_in_min = assignment.g_LoadingStartTimeInMin + pAgentSimu->m_veh_link_departure_time_in_simu_interval[nt - 1] * number_of_seconds_per_interval / 60.0;  // last link in the path

											path_time_vector[nt] = time_in_min;
										}

										float vehicle_travel_time = pAgentSimu->path_travel_time_in_min;
										speed_mph = path_distance / max(0.01, vehicle_travel_time / 60.0);
										//if (count >= 2000)  // only output up to 2000
										//    break;

										int hour = vehicle_departure_time / 60;
										int minute = (int)((vehicle_departure_time / 60.0 - hour) * 60);

										if (route_impacted_flag == true && pAgentSimu->m_bCompleteTrip == true)
										{
											total_impacted_vehicle_count++;

											if (pAgentSimu->diverted_flag >= 1)
											{
												total_diverted_vehicle_count++;
												total_diverted_vehicle_travel_time += vehicle_travel_time;
												total_diverted_vehicle_travel_distance += path_distance;

												if (pAgentSimu->info_receiving_flag == 2) // DMS
												{
													total_diverted_DMS_vehicle_count++;
													total_diverted_DMS_vehicle_travel_time += vehicle_travel_time;
													total_diverted_DMS_vehicle_travel_distance += path_distance;

												}
											}
											else
											{
												total_non_diverted_vehicle_count++;
												total_non_diverted_vehicle_travel_time += vehicle_travel_time;
												total_non_diverted_vehicle_travel_distance += path_distance;

											}
										}


										// some bugs for output link performances before
										fprintf(g_pFilePathMOE, "0,%d,%d,%d,%d->%d,%d,%d,%d,%d,%d,%d,%s,",
											pAgentSimu->agent_id,
											g_zone_vector[orig].zone_id,
											g_zone_vector[dest].zone_id,
											g_zone_vector[orig].zone_id,
											g_zone_vector[dest].zone_id,
											it->second.route_seq_id,
											it->second.path_seq_no,
											pAgentSimu->impacted_flag,
											pAgentSimu->info_receiving_flag,
											pAgentSimu->diverted_flag,
											pAgentSimu->mode_type_no,
											assignment.g_ModeTypeVector[at].mode_type.c_str());

										fprintf(g_pFilePathMOE, "%s,1,0,%.4f,T%02d%02d,%d,%d,%.4f,%.4f,",
											assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
											vehicle_departure_time,
											hour, minute,
											int(vehicle_departure_time),
											int(vehicle_departure_time / 5) * 5,
											vehicle_travel_time,
											path_distance, speed_mph);

										total_agent_distance += path_distance;

										total_agent_travel_time += vehicle_travel_time;

										fprintf(g_pFilePathMOE, "%.1f,%.1f,", waiting_time_in_min, max_link_waiting_time_in_min);

										if (max_link_waiting_time_in_min >= 0.1)
											fprintf(g_pFilePathMOE, "%d->%d,", from_node_id_max_wait, to_node_id_max_wait);
										else
											fprintf(g_pFilePathMOE, ",");


										//fprintf(g_pFilePathMOE, ")\"");
										fprintf(g_pFilePathMOE, "\n");

										count++;
									}
								}
							}
						}
					}
				}
			}

			assignment.summary_file << ",," << assignment.g_ModeTypeVector[at].mode_type.c_str() << "," << agent_count << "," << total_agent_distance / max(1, agent_count)
				<< "," << total_agent_travel_time / max(1, agent_count) << "," << total_agent_distance / max(0.0001, total_agent_travel_time / 60.0) << '\n';


		}
		fclose(g_pFilePathMOE);
	}
	assignment.summary_file << ",total_impacted_vehicles=," << total_impacted_vehicle_count << "," << '\n';
	assignment.summary_file << ",total_diverted_vehicle_count=," << total_diverted_vehicle_count << "," <<
		",precentage=," << total_diverted_vehicle_count * 100.0 / max(total_impacted_vehicle_count, 1) <<
		"%, avg travel time =, " << total_diverted_vehicle_travel_time / max(1, total_diverted_vehicle_count) <<
		", avg distance =, " << total_diverted_vehicle_travel_distance / max(1, total_diverted_vehicle_count) <<
		'\n';

	assignment.summary_file << ",total_diverted_DMS_veh_count=," << total_diverted_DMS_vehicle_count << "," <<
		",precentage=," << total_diverted_DMS_vehicle_count * 100.0 / max(total_impacted_vehicle_count, 1) <<
		"%, avg travel time =, " << total_diverted_DMS_vehicle_travel_time / max(1, total_diverted_DMS_vehicle_count) <<
		", avg distance =, " << total_diverted_DMS_vehicle_travel_distance / max(1, total_diverted_DMS_vehicle_count) <<
		'\n';

	assignment.summary_file << ",total_non_diverted_vehicle_count=," << total_non_diverted_vehicle_count << "," <<
		",precentage=," << total_non_diverted_vehicle_count * 100.0 / max(total_impacted_vehicle_count, 1) <<
		"%, avg travel time =, " << total_non_diverted_vehicle_travel_time / max(1, total_non_diverted_vehicle_count) <<
		", avg distance =, " << total_non_diverted_vehicle_travel_distance / max(1, total_non_diverted_vehicle_count) <<
		'\n';
}

void g_output_trajectory_csv(Assignment& assignment)
{

	if (assignment.assignment_mode == lue || assignment.trajectory_output_count == 0)  //LUE
	{
		FILE* g_pFilePathMOE = nullptr;
		fopen_ss(&g_pFilePathMOE, "trajectory.csv", "w");
		fclose(g_pFilePathMOE);
	}
	else if (assignment.assignment_mode >= 1)  //UE mode, or ODME, DTA
	{
		int total_number_of_nodes = 0;
		dtalog.output() << "[STATUS INFO] writing trajectory.csv.." << '\n';
		g_DTA_log_file << "[STATUS INFO] writing trajectory.csv.." << '\n';

		double path_time_vector[MAX_LINK_SIZE_IN_A_PATH];

		char trajectory_file_name[50];

		int scenario_no = assignment.g_active_DTAscenario_map[assignment.active_scenario_index];
		sprintf(trajectory_file_name, "trajectory_%d_%s.csv", assignment.active_scenario_index, assignment.g_DTA_scenario_vector[scenario_no].scenario_name.c_str());

		FILE* g_pFilePathMOE = nullptr;
		fopen_ss(&g_pFilePathMOE, trajectory_file_name, "w");

		if (!g_pFilePathMOE)
		{
			dtalog.output() << "[ERROR] File trajectory.csv cannot be opened." << '\n';
			g_DTA_log_file << "[ERROR] File trajectory.csv cannot be opened." << '\n';
			return;
		}

		fprintf(g_pFilePathMOE, "first_column,agent_id,o_zone_id,d_zone_id,route_seq_id,path_no,impacted_flag,info_receiving_flag,diverted_flag,mode_type_no,mode_type,demand_period,volume,toll,departure_time,travel_time,distance_km,distance_mile,speed_kmph,speed_mph,waiting_time,max_link_waiting_time,max_wait_link,node_sequence,time_sequence,geometry\n");

		int count = 1;

		clock_t start_t, end_t;
		start_t = clock();
		clock_t iteration_t;

		int mode_type_size = assignment.g_ModeTypeVector.size();
		int zone_size = g_zone_vector.size();
		int demand_period_size = assignment.g_DemandPeriodVector.size();

		CColumnVector* p_column_pool;

		float path_toll = 0;
		double path_distance_km = 0;
		double path_distance_mile = 0;
		float path_travel_time = 0;
		float speed_mph = 0;
		float speed_kmph = 0;
		float time_stamp = 0;

		if (assignment.trajectory_sampling_rate < 0.01)
			assignment.trajectory_sampling_rate = 0.01;
		int sampling_step = 100 / int(100 * assignment.trajectory_sampling_rate + 0.5);

		std::map<int, CColumnPath>::iterator it, it_begin, it_end;

		dtalog.output() << "[STATUS INFO] writing data for " << zone_size << "  zones " << '\n';
		g_DTA_log_file << "[STATUS INFO] writing data for " << zone_size << "  zones " << '\n';

		for (int orig = 0; orig < zone_size; ++orig)
		{
			int from_zone_sindex = g_zone_vector[orig].sindex;
			if (from_zone_sindex == -1)
				continue;

			if (g_zone_vector[orig].zone_id % 100 == 0)
			{ 
				dtalog.output() << "[DATA INFO] o zone id =  " << g_zone_vector[orig].zone_id << '\n';
				g_DTA_log_file << "[DATA INFO] o zone id =  " << g_zone_vector[orig].zone_id << '\n';
			}
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
								if (count % 100000 == 0)
								{
									end_t = clock();
									iteration_t = end_t - start_t;
									dtalog.output() << "[STATUS INFO] writing " << count / 1000 << "K agents with CPU time " << iteration_t / 1000.0 << " s" << '\n';
									g_DTA_log_file << "[STATUS INFO] writing " << count / 1000 << "K agents with CPU time " << iteration_t / 1000.0 << " s" << '\n';
								}

								if (count % sampling_step != 0)
									continue;

								if (count > assignment.trajectory_output_count && assignment.trajectory_output_count >= 1)  // if trajectory_output_count ==-1, this constraint does not apply
									continue;

								path_toll = 0;
								path_distance_km = 0;
								path_distance_mile = 0;
								path_travel_time = 0;
								path_time_vector[0] = time_stamp;

								for (int nl = 0; nl < it->second.m_link_size; ++nl)  // arc a
								{
									int link_seq_no = it->second.path_link_vector[nl];
									path_toll += g_link_vector[link_seq_no].VDF_period[tau].toll[at][assignment.active_scenario_index];
									path_distance_km += g_link_vector[link_seq_no].link_distance_km;
									path_distance_mile += g_link_vector[link_seq_no].link_distance_mile;
									float link_travel_time = g_link_vector[link_seq_no].link_avg_travel_time_per_period[tau][at];

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
										if (assignment.trajectory_diversion_only == 1 && pAgentSimu->diverted_flag == 0)  // diversion flag only, then we skip the non-diversion path
											continue;

										float vehicle_departure_time = assignment.g_LoadingStartTimeInMin + pAgentSimu->departure_time_in_min;

										float waiting_time_in_min = pAgentSimu->waiting_time_in_min;
										float max_link_waiting_time_in_min = pAgentSimu->max_link_waiting_time_in_min;

										int from_node_id_max_wait = -1;
										int to_node_id_max_wait = -1;

										if (pAgentSimu->max_waiting_time_link_no >= 0)
										{
											from_node_id_max_wait = g_node_vector[g_link_vector[pAgentSimu->max_waiting_time_link_no].from_node_seq_no].node_id;
											to_node_id_max_wait = g_node_vector[g_link_vector[pAgentSimu->max_waiting_time_link_no].to_node_seq_no].node_id;
										}

										time_stamp = assignment.g_LoadingStartTimeInMin + pAgentSimu->departure_time_in_min;
										for (int nt = 0 + virtual_link_delta; nt < pAgentSimu->path_link_seq_no_vector.size() + 1 - virtual_link_delta; ++nt)
										{
											float time_in_min = 0;

											if (nt < pAgentSimu->path_link_seq_no_vector.size() - virtual_link_delta)
												time_in_min = assignment.g_LoadingStartTimeInMin + pAgentSimu->m_veh_link_arrival_time_in_simu_interval[nt] * number_of_seconds_per_interval / 60.0;
											else
												time_in_min = assignment.g_LoadingStartTimeInMin + pAgentSimu->m_veh_link_departure_time_in_simu_interval[nt - 1] * number_of_seconds_per_interval / 60.0;  // last link in the path

											path_time_vector[nt] = time_in_min;
										}

										float vehicle_travel_time = pAgentSimu->path_travel_time_in_min;
										speed_mph = path_distance_mile / max(0.01, vehicle_travel_time / 60.0);
										speed_kmph = path_distance_km / max(0.01, vehicle_travel_time / 60.0);
										//if (count >= 2000)  // only output up to 2000
										//    break;

										// some bugs for output link performances before
										fprintf(g_pFilePathMOE, "0,%d,%d,%d,%d,%d,%d,%d,%d,%d,%s,",
											pAgentSimu->agent_id,
											g_zone_vector[orig].zone_id,
											g_zone_vector[dest].zone_id,
											it->second.route_seq_id,
											it->second.path_seq_no,
											pAgentSimu->impacted_flag,
											pAgentSimu->info_receiving_flag,
											pAgentSimu->diverted_flag,
											pAgentSimu->mode_type_no,
											assignment.g_ModeTypeVector[at].mode_type.c_str());

										fprintf(g_pFilePathMOE, "%s,1,0,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,",
											assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
											vehicle_departure_time,
											vehicle_travel_time,
											path_distance_km, path_distance_mile, speed_kmph, speed_mph);

										fprintf(g_pFilePathMOE, "%.1f,%.1f,", waiting_time_in_min, max_link_waiting_time_in_min);

										if (max_link_waiting_time_in_min >= 0.1)
											fprintf(g_pFilePathMOE, "%d->%d,", from_node_id_max_wait, to_node_id_max_wait);
										else
											fprintf(g_pFilePathMOE, ",");


										total_number_of_nodes += (pAgentSimu->path_link_seq_no_vector.size() + 1);
										/* Format and print various data */
										for (int ni = 0 + virtual_link_delta; ni < pAgentSimu->path_link_seq_no_vector.size(); ++ni)
										{
											int node_id = g_node_vector[g_link_vector[pAgentSimu->path_link_seq_no_vector[ni]].from_node_seq_no].node_id;
											fprintf(g_pFilePathMOE, "%d;", node_id);
										}

										fprintf(g_pFilePathMOE, ",");

										//for (int nt = 0 + virtual_link_delta; nt < pAgentSimu->path_link_seq_no_vector.size() + 1 - virtual_link_delta; ++nt)
										//	fprintf(g_pFilePathMOE, "%f;", path_time_vector[nt]); ! we do not use this format from nowL 0308 2023

										// time coded
										for (int nt = 0 + virtual_link_delta; nt < pAgentSimu->path_link_seq_no_vector.size() + 1 - virtual_link_delta; ++nt)
											fprintf(g_pFilePathMOE, "%s;", g_time_coding(path_time_vector[nt]).c_str());

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


										fprintf(g_pFilePathMOE, ")\"");
										fprintf(g_pFilePathMOE, "\n");

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
		assignment.summary_file << ", # of simulated agents in trajectory.csv=," << count << ",avg # of nodes per agent=" << total_number_of_nodes * 1.0 / max(1, count) << '\n';
	}

	int b_trace_file = false;

	// output trace file
	if (assignment.assignment_mode == simulation_dta)  //LUE
	{
		FILE* g_pFilePathMOE = nullptr;
		fopen_ss(&g_pFilePathMOE, "trace.csv", "w");
		dtalog.output() << "[STATUS INFO] writing trace.csv.." << '\n';
		g_DTA_log_file << "[STATUS INFO] writing trace.csv.." << '\n';

		double path_time_vector[MAX_LINK_SIZE_IN_A_PATH];

		fopen_ss(&g_pFilePathMOE, "trace.csv", "w");

		if (!g_pFilePathMOE)
		{
			dtalog.output() << "[ERROR] File trace.csv cannot be opened." << '\n';
			g_DTA_log_file << "[ERROR] File trace.csv cannot be opened." << '\n';
			return;
		}

		fprintf(g_pFilePathMOE, "agent_id,seq_no,node_id,timestamp,timestamp_in_min,trip_time_in_min,travel_time_in_sec,waiting_time_in_simu_interval,x_coord,y_coord\n");

		int count = 1;

		clock_t start_t, end_t;
		start_t = clock();
		clock_t iteration_t;

		int mode_type_size = assignment.g_ModeTypeVector.size();
		int zone_size = g_zone_vector.size();
		int demand_period_size = assignment.g_DemandPeriodVector.size();

		CColumnVector* p_column_pool;

		float path_toll = 0;
		float path_distance = 0;
		float path_travel_time = 0;
		float time_stamp = 0;

		if (assignment.trajectory_sampling_rate < 0.01)
			assignment.trajectory_sampling_rate = 0.01;
		int sampling_step = 1;

		std::map<int, CColumnPath>::iterator it, it_begin, it_end;

		dtalog.output() << "[STATUS INFO] writing data for " << zone_size << "  zones " << '\n';
		g_DTA_log_file << "[STATUS INFO] writing data for " << zone_size << "  zones " << '\n';

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
								if (count % 100000 == 0)
								{
									end_t = clock();
									iteration_t = end_t - start_t;
									dtalog.output() << "[STATUS INFO] writing " << count / 1000 << "K agents with CPU time " << iteration_t / 1000.0 << " s" << '\n';
									g_DTA_log_file << "[STATUS INFO] writing " << count / 1000 << "K agents with CPU time " << iteration_t / 1000.0 << " s" << '\n';
								}

								if (count % sampling_step != 0)
									continue;

								if (count >= 2000)  // only output up to 2000
									break;

								path_toll = 0;
								path_distance = 0;
								path_travel_time = 0;
								path_time_vector[0] = time_stamp;

								for (int nl = 0; nl < it->second.m_link_size; ++nl)  // arc a
								{
									int link_seq_no = it->second.path_link_vector[nl];
									path_toll += g_link_vector[link_seq_no].VDF_period[tau].toll[at][assignment.active_scenario_index];
									path_distance += g_link_vector[link_seq_no].link_distance_VDF;
									float link_travel_time = g_link_vector[link_seq_no].link_avg_travel_time_per_period[tau][at];
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
										//if (assignment.trajectory_diversion_only == 1 && pAgentSimu->diverted_flag == 0)  // diversion flag only, then we skip the non-diversion path
										//     continue;

										for (int nt = 0 + virtual_link_delta; nt < pAgentSimu->path_link_seq_no_vector.size() + 1 - virtual_link_delta; ++nt)
										{
											float time_in_min = 0;

											if (nt < pAgentSimu->path_link_seq_no_vector.size() - virtual_link_delta)
												time_in_min = assignment.g_LoadingStartTimeInMin + pAgentSimu->m_veh_link_arrival_time_in_simu_interval[nt] * number_of_seconds_per_interval / 60.0;
											else
												time_in_min = assignment.g_LoadingStartTimeInMin + pAgentSimu->m_veh_link_departure_time_in_simu_interval[nt - 1] * number_of_seconds_per_interval / 60.0;  // last link in the path

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
											int link_seq_no = pAgentSimu->path_link_seq_no_vector[nt];
											if (link_seq_no < 0)
											{
												int i_error = 1;
												continue;
											}
											float travel_time_in_sec = (path_time_vector[nt + 1] - path_time_vector[nt]) * 60;
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

void g_output_trajectory_bin(Assignment& assignment)
{
	if (assignment.assignment_mode == 0 || g_agent_simu_vector.size() <= 200000)
		return;


	if (assignment.assignment_mode == 0 || assignment.trajectory_output_count == 0)  //LUE
	{
		FILE* g_pFilePathMOE = nullptr;
		fopen_ss(&g_pFilePathMOE, "trajectory.bin", "wb");
		fclose(g_pFilePathMOE);
	}
	else if (assignment.assignment_mode >= 1)  //UE mode, or ODME, DTA
	{
		dtalog.output() << "[STATUS INFO] writing trajectory.bin.." << '\n';
		g_DTA_log_file << "[STATUS INFO] writing trajectory.bin.." << '\n';

		int path_node_vector[MAX_LINK_SIZE_IN_A_PATH];
		double path_time_vector[MAX_LINK_SIZE_IN_A_PATH];
		FILE* g_pFilePathMOE = nullptr;
		fopen_ss(&g_pFilePathMOE, "trajectory.bin", "wb");

		if (!g_pFilePathMOE)
		{
			dtalog.output() << "[ERROR] File trajectory.bin cannot be opened." << '\n';
			g_DTA_log_file << "[ERROR] File trajectory.bin cannot be opened." << '\n';
			return; 
		}


		struct STrajectoryHeader
		{
			int agent_id, o_zone_id, d_zone_id, path_id;
			int display_code;
			int impacted_flag, info_receiving_flag, diverted_flag;
			int mode_type_no, PCE_unit, demand_period;
			int node_size, link_size;
			float volume, toll, travel_time, distance_km;
			int reserved1;
			int reserved2;
			float distance_mile;
			float reserved4;

		};

		STrajectoryHeader header;

		int count = 1;

		clock_t start_t, end_t;
		start_t = clock();
		clock_t iteration_t;

		int mode_type_size = assignment.g_ModeTypeVector.size();
		int zone_size = g_zone_vector.size();
		int demand_period_size = assignment.g_DemandPeriodVector.size();

		CColumnVector* p_column_pool;

		float path_toll = 0;
		float path_distance_km = 0;
		float path_distance_mile = 0;
		float path_travel_time = 0;
		float time_stamp = 0;

		if (assignment.trajectory_sampling_rate < 0.01)
			assignment.trajectory_sampling_rate = 0.01;
		int sampling_step = 100 / int(100 * assignment.trajectory_sampling_rate + 0.5);

		std::map<int, CColumnPath>::iterator it, it_begin, it_end;

		dtalog.output() << "[STATUS INFO] writing data for " << zone_size << "  zones " << '\n';
		g_DTA_log_file << "[STATUS INFO] writing data for " << zone_size << "  zones " << '\n';

		for (int orig = 0; orig < zone_size; ++orig)
		{
			int from_zone_sindex = g_zone_vector[orig].sindex;
			if (from_zone_sindex == -1)
				continue;

			if (g_zone_vector[orig].zone_id % 100 == 0)
			{
				dtalog.output() << "[DATA INFO] o zone id =  " << g_zone_vector[orig].zone_id << '\n';
				g_DTA_log_file << "[DATA INFO] o zone id =  " << g_zone_vector[orig].zone_id << '\n';
			}
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
								if (count % 100000 == 0)
								{
									end_t = clock();
									iteration_t = end_t - start_t;
									dtalog.output() << "[STATUS INFO] writing " << count / 1000 << "K binary agents with CPU time " << iteration_t / 1000.0 << " s" << '\n';
									g_DTA_log_file << "[STATUS INFO] writing " << count / 1000 << "K binary agents with CPU time " << iteration_t / 1000.0 << " s" << '\n';
								}

								if (count % sampling_step != 0)
									continue;


								path_toll = 0;
								path_distance_km = 0;
								path_distance_mile = 0;
								path_travel_time = 0;
								path_time_vector[0] = time_stamp;

								for (int nl = 0; nl < it->second.m_link_size; ++nl)  // arc a
								{
									int link_seq_no = it->second.path_link_vector[nl];
									path_toll += g_link_vector[link_seq_no].VDF_period[tau].toll[at][assignment.active_scenario_index];
									path_distance_km += g_link_vector[link_seq_no].link_distance_km;
									path_distance_mile += g_link_vector[link_seq_no].link_distance_mile;
									float link_travel_time = g_link_vector[link_seq_no].link_avg_travel_time_per_period[tau][at];

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
										if (assignment.trajectory_diversion_only == 1 && pAgentSimu->diverted_flag == 0)  // diversion flag only, then we skip the non-diversion path
											continue;

										time_stamp = assignment.g_LoadingStartTimeInMin + pAgentSimu->departure_time_in_min;
										for (int nt = 0 + virtual_link_delta; nt < pAgentSimu->path_link_seq_no_vector.size() + 1 - virtual_link_delta; ++nt)
										{
											double time_in_min = 0;

											if (nt < pAgentSimu->path_link_seq_no_vector.size() - virtual_link_delta)
												time_in_min = assignment.g_LoadingStartTimeInMin + pAgentSimu->m_veh_link_arrival_time_in_simu_interval[nt] * number_of_seconds_per_interval / 60.0;
											else
												time_in_min = assignment.g_LoadingStartTimeInMin + pAgentSimu->m_veh_link_departure_time_in_simu_interval[nt - 1] * number_of_seconds_per_interval / 60.0;  // last link in the path

											path_time_vector[nt - virtual_link_delta] = time_in_min;
										}

										float vehicle_travel_time = pAgentSimu->path_travel_time_in_min;


										header.agent_id = pAgentSimu->agent_id;
										header.o_zone_id = g_zone_vector[orig].zone_id;
										header.d_zone_id = g_zone_vector[dest].zone_id;
										header.path_id = it->second.path_seq_no + 1;

										header.impacted_flag = pAgentSimu->impacted_flag; // assignment.g_ModeTypeVector[at].display_code_no;
										header.info_receiving_flag = pAgentSimu->info_receiving_flag; // assignment.g_ModeTypeVector[at].display_code_no;
										header.diverted_flag = pAgentSimu->diverted_flag; // assignment.g_ModeTypeVector[at].display_code_no;

										header.mode_type_no = pAgentSimu->mode_type_no; // assignment.g_ModeTypeVector[at].mode_type_no;
										header.PCE_unit = pAgentSimu->PCE_unit_size; // pAgentSimu->PCE_unit_size;
										header.demand_period = tau;
										header.toll = path_toll;
										header.travel_time = vehicle_travel_time;
										header.distance_km = path_distance_km;
										header.distance_mile = path_distance_mile;
										header.node_size = pAgentSimu->path_link_seq_no_vector.size() - virtual_link_delta;
										header.link_size = pAgentSimu->path_link_seq_no_vector.size() + 1 - virtual_link_delta;

										fwrite(&header, sizeof(struct STrajectoryHeader), 1, g_pFilePathMOE);

										/* Format and print various data */

										for (int ni = 0 + virtual_link_delta; ni < pAgentSimu->path_link_seq_no_vector.size(); ++ni)
										{
											int link_seq_no = pAgentSimu->path_link_seq_no_vector[ni];
											if (link_seq_no < 0)
											{
												int i_error = 1;
												continue;
											}

											path_node_vector[ni - virtual_link_delta] = g_node_vector[g_link_vector[pAgentSimu->path_link_seq_no_vector[ni]].from_node_seq_no].node_id;
										}

										fwrite(&path_node_vector, sizeof(int), header.node_size, g_pFilePathMOE);
										fwrite(&path_time_vector, sizeof(double), header.link_size, g_pFilePathMOE);


										count++;
									}
								}
							}
						}
					}
				}
			}
		}

		end_t = clock();
		iteration_t = end_t - start_t;
		dtalog.output() << "[STATUS INFO] Comlete writing " << count / 1000 << "K binary agents with CPU time " << iteration_t / 1000.0 << " s." << '\n';
		g_DTA_log_file << "[STATUS INFO] Comlete writing " << count / 1000 << "K binary agents with CPU time " << iteration_t / 1000.0 << " s." << '\n';
		fclose(g_pFilePathMOE);
	}

}

void g_output_demand_bin(Assignment& assignment)
{

	dtalog.output() << "[STATUS INFO] writing demand.bin.." << '\n';
	g_DTA_log_file << "[STATUS INFO] writing demand.bin.." << '\n';

	FILE* g_pFilePathMOE = nullptr;
	fopen_ss(&g_pFilePathMOE, "output_demand.bin", "wb");

	if (!g_pFilePathMOE)
	{
		dtalog.output() << "[ERROR] File demand.bin cannot be opened." << '\n';
		g_DTA_log_file << "[ERROR] File demand.bin cannot be opened." << '\n';
		return;
	}


	struct SDemandHeader
	{
		int o_zone_id, d_zone_id, mode_type_no, demand_period;
		double volume;
	};

	SDemandHeader header;

	int count = 1;

	clock_t start_t, end_t;
	start_t = clock();
	clock_t iteration_t;

	int mode_type_size = assignment.g_ModeTypeVector.size();
	int zone_size = g_zone_vector.size();
	int demand_period_size = assignment.g_DemandPeriodVector.size();

	CColumnVector* p_column_pool;

	std::map<int, CColumnPath>::iterator it, it_begin, it_end;

	dtalog.output() << "[STATUS INFO] writing data for " << zone_size << "  zones " << '\n';
	g_DTA_log_file << "[STATUS INFO] writing data for " << zone_size << "  zones " << '\n';

	for (int orig = 0; orig < zone_size; ++orig)
	{
		int from_zone_sindex = g_zone_vector[orig].sindex;
		if (from_zone_sindex == -1)
			continue;

		if (g_zone_vector[orig].zone_id % 100 == 0)
		{
			dtalog.output() << "[DATA INFO] o zone id =  " << g_zone_vector[orig].zone_id << '\n';
			g_DTA_log_file << "[DATA INFO] o zone id =  " << g_zone_vector[orig].zone_id << '\n';
		}
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

						header.o_zone_id = g_zone_vector[orig].zone_id;
						header.d_zone_id = g_zone_vector[dest].zone_id;
						header.mode_type_no = at; // assignment.g_ModeTypeVector[at].mode_type_no;
						header.demand_period = tau;
						header.volume = p_column_pool->od_volume[assignment.active_scenario_index];


						fwrite(&header, sizeof(struct SDemandHeader), 1, g_pFilePathMOE);
						count++;

					}
				}
			}
		}
	}

	end_t = clock();
	iteration_t = end_t - start_t;
	dtalog.output() << "[STATUS INFO] Complete writing " << count / 1000 << "K binary demand pairs with CPU time " << iteration_t / 1000.0 << " s." << '\n';
	g_DTA_log_file << "[STATUS INFO] Complete writing " << count / 1000 << "K binary demand pairs with CPU time " << iteration_t / 1000.0 << " s." << '\n';
	fclose(g_pFilePathMOE);
}


void g_output_simulation_result(Assignment& assignment)
{

	//    g_output_dynamic_link_performance(assignment, 2);
	g_output_TD_link_performance(assignment, 1);
	g_output_dynamic_link_state(assignment, 1);
	g_output_trajectory_bin(assignment);
	//g_output_simulation_agents(assignment);
	g_output_trajectory_csv(assignment);


}

void g_OutputModelFiles(int mode)
{
	if (mode == 1)
	{
		FILE* g_pFileModelNode = fopen("model_node.csv", "w");

		if (g_pFileModelNode != NULL)
		{
			fprintf(g_pFileModelNode, "node_id,node_no,layer_no,MRM_gate_flag,node_type,is_boundary,#_of_outgoing_nodes,activity_node_flag,mode_type,zone_id,cell_code,info_zone_flag,x_coord,y_coord\n");
			for (int i = 0; i < g_node_vector.size(); i++)
			{


				if (g_node_vector[i].node_id >= 0)  //for all physical links
				{

					fprintf(g_pFileModelNode, "%d,%d,%d,%d,%s,%d,%d,%d,%s, %ld,%s,%d,%f,%f\n",
						g_node_vector[i].node_id,
						g_node_vector[i].node_seq_no,
						g_node_vector[i].layer_no,
						g_node_vector[i].MRM_gate_flag,
						g_node_vector[i].node_type.c_str(),
						g_node_vector[i].is_boundary,
						g_node_vector[i].m_outgoing_link_seq_no_vector.size(),
						g_node_vector[i].is_activity_node,
						g_node_vector[i].mode_type_str.c_str(),

						g_node_vector[i].zone_org_id,
						g_node_vector[i].cell_str.c_str(),
						0,
						g_node_vector[i].x,
						g_node_vector[i].y
					);

				}

			}

			fclose(g_pFileModelNode);
		}
		else
		{
			dtalog.output() << "[ERROR] File model_node.csv cannot be opened.\n It might be currently used and locked by EXCEL." << '\n';
			g_DTA_log_file << "[ERROR] File model_node.csv cannot be opened.\n It might be currently used and locked by EXCEL." << '\n';
			return;


		}

	}

	if (mode == 2)
	{
		FILE* g_pFileModelLink = fopen("model_link.csv", "w");

		if (g_pFileModelLink != NULL)
		{
			fprintf(g_pFileModelLink, "link_id,link_no,layer_no,from_node_id,to_node_id,from_gate_flag,to_gate_flag,link_type,link_type_name,link_type_code,lanes,link_distance_VDF,free_speed,cutoff_speed,fftt,capacity,allow_uses,");
			fprintf(g_pFileModelLink, "BPR_plf,BPR_alpha,BPR_beta,QVDF_plf,QVDF_alpha,QVDF_beta,QVDF_cd,QVDF_n,");
			fprintf(g_pFileModelLink, "geometry\n");

			//VDF_fftt1,VDF_cap1,VDF_alpha1,VDF_beta1
			for (int i = 0; i < g_link_vector.size(); i++)
			{
				//if (g_link_vector[i].link_type_si[0] <= -100)  // invisible link
				//{
				//	continue;
				//}

				fprintf(g_pFileModelLink, "%s,%d,%d,%d,%d,%d,%d,%d,%s,%s,%f,%f,%f,%f,%f,%f,%s,%f,%f,%f,%f,%f,%f,%f,%f,",
					g_link_vector[i].link_id.c_str(),
					g_link_vector[i].link_seq_no,
					g_link_vector[i].layer_no,
					g_node_vector[g_link_vector[i].from_node_seq_no].node_id,
					g_node_vector[g_link_vector[i].to_node_seq_no].node_id,
					g_node_vector[g_link_vector[i].from_node_seq_no].MRM_gate_flag,
					g_node_vector[g_link_vector[i].to_node_seq_no].MRM_gate_flag,
					g_link_vector[i].link_type_si[0],
					g_link_vector[i].link_type_name.c_str(),
					g_link_vector[i].link_type_name.c_str(),
					g_link_vector[i].number_of_lanes_si[0],
					g_link_vector[i].link_distance_VDF,
					g_link_vector[i].free_speed,
					g_link_vector[i].v_congestion_cutoff,
					g_link_vector[i].free_flow_travel_time_in_min,
					g_link_vector[i].lane_capacity,
					g_link_vector[i].VDF_period[0].allowed_uses[0].c_str(),
					g_link_vector[i].VDF_period[0].Q_peak_load_factor,
					g_link_vector[i].VDF_period[0].alpha,
					g_link_vector[i].VDF_period[0].beta,
					g_link_vector[i].VDF_period[0].Q_peak_load_factor,
					g_link_vector[i].VDF_period[0].Q_alpha,
					g_link_vector[i].VDF_period[0].Q_beta,
					g_link_vector[i].VDF_period[0].Q_cd,
					g_link_vector[i].VDF_period[0].Q_n
				);

				if (g_link_vector[i].geometry.size() > 0)
				{
					fprintf(g_pFileModelLink, "\"%s\",", g_link_vector[i].geometry.c_str());
				}
				else
				{
					fprintf(g_pFileModelLink, "\"LINESTRING (");

					fprintf(g_pFileModelLink, "%f %f,", g_node_vector[g_link_vector[i].from_node_seq_no].x, g_node_vector[g_link_vector[i].from_node_seq_no].y);
					fprintf(g_pFileModelLink, "%f %f", g_node_vector[g_link_vector[i].to_node_seq_no].x, g_node_vector[g_link_vector[i].to_node_seq_no].y);

					fprintf(g_pFileModelLink, ")\"");
				}

				fprintf(g_pFileModelLink, "\n");

			}


			fclose(g_pFileModelLink);
		}
		else
		{
			dtalog.output() << "[ERROR] File model_link.csv cannot be opened.\n It might be currently used and locked by EXCEL." << '\n';
			g_DTA_log_file << "[ERROR] File model_link.csv cannot be opened.\n It might be currently used and locked by EXCEL." << '\n';
			return;

		}

	}
	if (mode == 3)
	{
		int connector_count = 0;
		for (int i = 0; i < g_link_vector.size(); i++)
		{
			if (g_link_vector[i].link_type_si[0] == 1000)  // connector
			{
				connector_count += 1;
			}
		}

		if (connector_count >= 1)
		{
			FILE* g_pFileModelLink = fopen("access_link.csv", "w");

			if (g_pFileModelLink != NULL)
			{
				fprintf(g_pFileModelLink, "link_id,link_no,from_node_id,to_node_id,link_type,link_type_name,lanes,link_distance_VDF,free_speed,fftt,capacity,allow_uses,geometry\n");

				//VDF_fftt1,VDF_cap1,VDF_alpha1,VDF_beta1
				for (int i = 0; i < g_link_vector.size(); i++)
				{
					if (g_link_vector[i].link_type_si[0] == 1000)  // connector
					{

						fprintf(g_pFileModelLink, "%s,%d,%d,%d,%d,%s,%d,%f,%f,%f,%f,%s,",
							g_link_vector[i].link_id.c_str(),
							g_link_vector[i].link_seq_no,
							g_node_vector[g_link_vector[i].from_node_seq_no].node_id,
							g_node_vector[g_link_vector[i].to_node_seq_no].node_id,
							g_link_vector[i].link_type_si[0],
							g_link_vector[i].link_type_name.c_str(),
							g_link_vector[i].number_of_lanes_si[0],
							g_link_vector[i].link_distance_VDF,
							g_link_vector[i].free_speed,
							g_link_vector[i].free_flow_travel_time_in_min,
							g_link_vector[i].lane_capacity,
							g_link_vector[i].VDF_period[0].allowed_uses[0].c_str()
							//g_link_vector[i].VDF_period[0].FFTT,
						   //g_link_vector[i].VDF_period[0].period_capacity,
						   //g_link_vector[i].VDF_period[0].alpha,
						   //g_link_vector[i].VDF_period[0].beta,
						);

						if (g_link_vector[i].geometry.size() > 0)
						{
							fprintf(g_pFileModelLink, "\"%s\",", g_link_vector[i].geometry.c_str());
						}
						else
						{
							fprintf(g_pFileModelLink, "\"LINESTRING (");

							fprintf(g_pFileModelLink, "%f %f,", g_node_vector[g_link_vector[i].from_node_seq_no].x, g_node_vector[g_link_vector[i].from_node_seq_no].y);
							fprintf(g_pFileModelLink, "%f %f", g_node_vector[g_link_vector[i].to_node_seq_no].x, g_node_vector[g_link_vector[i].to_node_seq_no].y);

							fprintf(g_pFileModelLink, ")\"");
						}

						fprintf(g_pFileModelLink, "\n");

					}


				}

				fclose(g_pFileModelLink);
			}
			else
			{
				dtalog.output() << "[ERROR] File access_link.csv cannot be opened.\n It might be currently used and locked by EXCEL." << '\n';
				g_DTA_log_file << "[ERROR] File access_link.csv cannot be opened.\n It might be currently used and locked by EXCEL." << '\n';
				return;

			}

		}

	}

	if (mode == 3)  // cell
	{
		//FILE* g_pFileZone = nullptr;
		//g_pFileZone = fopen("model_cell.csv", "w");

		//if (g_pFileZone == NULL)
		//{
		//    dtalog.output() << "File model_cell.csv cannot be opened." << '\n';
		//    g_DTA_log_file << "File model_cell.csv cannot be opened." << '\n';
		//    g_program_stop();
		//}
		//else
		//{


		//    fprintf(g_pFileZone, "cell_code,geometry\n");

		//    std::map<string, CInfoCell>::iterator it;

		//    for (it = g_info_cell_map.begin(); it != g_info_cell_map.end(); ++it)
		//    {

		//        fprintf(g_pFileZone, "%s,", it->first.c_str());
		//        fprintf(g_pFileZone, "\"LINESTRING (");

		//        for (int s = 0; s < it->second.m_ShapePoints.size(); s++)
		//        {
		//            fprintf(g_pFileZone, "%f %f,", it->second.m_ShapePoints[s].x, it->second.m_ShapePoints[s].y);
		//        }

		//        fprintf(g_pFileZone, ")\"");
		//        fprintf(g_pFileZone, "\n");
		//    }
		//    fclose(g_pFileZone);
		//}
	}

	if (mode == 10)
	{
		FILE* g_pFileModel_LC = fopen("log_shortest_path_tree.csv", "w");

		if (g_pFileModel_LC != NULL)
		{
			fprintf(g_pFileModel_LC, "iteration,mode_type,zone_id,node_id,d_zone_id,connected_flag,pred,label_cost,pred_link_cost,x_coord,y_coord,geometry\n");
			for (int i = 0; i < g_node_vector.size(); i++)
			{


				//                if (g_node_vector[i].node_id >= 0)  //for all physical links
				{
					std::map<string, float> ::iterator it;

					for (it = g_node_vector[i].label_cost_per_iteration_map.begin(); it != g_node_vector[i].label_cost_per_iteration_map.end(); ++it)
					{
						int node_pred_id = -1;
						int pred_no = g_node_vector[i].pred_per_iteration_map[it->first];
						float pred_link_cost = 0;
						if (pred_no >= 0)
						{
							std::map<string, float> ::iterator it2;

							float pred_node_label_cost = 0;
							for (it2 = g_node_vector[pred_no].label_cost_per_iteration_map.begin(); it2 != g_node_vector[pred_no].label_cost_per_iteration_map.end(); ++it2)
							{
								pred_node_label_cost = it2->second;
							}

							pred_link_cost = it->second - pred_node_label_cost;
							node_pred_id = g_node_vector[pred_no].node_id;

						}
						int d_zone_id = g_node_vector[i].zone_id;
						int connected_flag = 0;

						if (it->second < 100000)
							connected_flag = 1;

						//if (connected_flag == 1)
						{
							fprintf(g_pFileModel_LC, "%s,%d,%d,%d,%f,%f,%f,%f,", it->first.c_str(), d_zone_id, connected_flag, node_pred_id, it->second,
								pred_link_cost,
								g_node_vector[i].x, g_node_vector[i].y);

							int pred_no_2 = i;
							if (pred_no >= 0)
								pred_no_2 = pred_no;

							fprintf(g_pFileModel_LC, "\"LINESTRING (");
							fprintf(g_pFileModel_LC, "%f %f,", g_node_vector[i].x, g_node_vector[i].y);

							fprintf(g_pFileModel_LC, "%f %f", g_node_vector[pred_no_2].x, g_node_vector[pred_no_2].y);

							fprintf(g_pFileModel_LC, ")\"");

							fprintf(g_pFileModel_LC, "\n");
						}
					}

				}

			}

			fclose(g_pFileModel_LC);
		}
		else
		{
			dtalog.output() << "[ERROR] File model_label_cost_tree.csv cannot be opened.\n It might be currently used and locked by EXCEL." << '\n';
			g_DTA_log_file << "[ERROR] File model_label_cost_tree.csv cannot be opened.\n It might be currently used and locked by EXCEL." << '\n';
			return;


		}

	}

	if (mode == 20)
	{
		FILE* g_pFileModel_LC = fopen("model_RT_shortest_tree.csv", "w");

		if (g_pFileModel_LC != NULL)
		{
			fprintf(g_pFileModel_LC, "demand_period,mode_type,d_zone_id,node_id,o_zone_id,connected_flag,pred,label_cost,geometry\n");
			for (int i = 0; i < g_node_vector.size(); i++)
			{
				//                if (g_node_vector[i].node_id >= 0)  //for all physical links
				{
					std::map<string, float> ::iterator it;

					for (it = g_node_vector[i].label_cost_RT_map.begin(); it != g_node_vector[i].label_cost_RT_map.end(); ++it)
					{
						int node_pred_id = -1;
						int pred_no = g_node_vector[i].pred_RT_map[it->first];
						if (pred_no >= 0)
							node_pred_id = g_node_vector[pred_no].node_id;
						int o_zone_id = g_node_vector[i].zone_id;
						int connected_flag = 0;

						if (it->second < 1000)  // feasible link only
						{
							connected_flag = 1;
							fprintf(g_pFileModel_LC, "%s,%d,%d,%d,%f,", it->first.c_str(), o_zone_id, connected_flag, node_pred_id, it->second);

							if (node_pred_id >= 0)
							{
								fprintf(g_pFileModel_LC, "\"LINESTRING (");
								fprintf(g_pFileModel_LC, "%f %f,", g_node_vector[i].x, g_node_vector[i].y);
								fprintf(g_pFileModel_LC, "%f %f,", g_node_vector[pred_no].x, g_node_vector[pred_no].y);
								fprintf(g_pFileModel_LC, ")\"");
								fprintf(g_pFileModel_LC, ",");
							}

							fprintf(g_pFileModel_LC, "\n");
						}
					}

				}

			}

			fclose(g_pFileModel_LC);
		}
	}
}
