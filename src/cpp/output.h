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

#include "DTA.h"

void g_record_corridor_performance_summary(Assignment& assignment, int base_case_flag)
{
	for (int i = 0; i < g_link_vector.size(); ++i)
	{
		// virtual connectors
		if (g_link_vector[i].link_type_si[0] == -1)
			continue;


		for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
		{

			if (assignment.g_DemandPeriodVector[tau].number_of_demand_files == 0)
				continue;

			float speed = g_link_vector[i].free_speed;  // default speed 
			if (g_link_vector[i].VDF_period[tau].avg_travel_time > 0.001f)
				speed = g_link_vector[i].free_speed / max(0.000001, g_link_vector[i].VDF_period[tau].avg_travel_time) * g_link_vector[i].VDF_period[tau].FFTT_at[0];

			float speed_ratio = speed / max(1.0, g_link_vector[i].free_speed);  // default speed 
			float vehicle_volume = g_link_vector[i].total_volume_for_all_mode_types_per_period[tau] + g_link_vector[i].VDF_period[tau].preload;
			float person_volume = g_link_vector[i].total_person_volume_for_all_mode_types_per_period[tau] + g_link_vector[i].VDF_period[tau].preload;
			//VMT,VHT,PMT,PHT,PDT
			float preload = g_link_vector[i].VDF_period[tau].preload;
			float VMT = vehicle_volume * g_link_vector[i].link_distance_VDF;

			float VHT = vehicle_volume * g_link_vector[i].VDF_period[tau].avg_travel_time / 60.0;
			float PMT = person_volume * g_link_vector[i].link_distance_VDF;
			float PHT = person_volume * g_link_vector[i].VDF_period[tau].avg_travel_time / 60.0;
			float PDT_vf = person_volume * (g_link_vector[i].VDF_period[tau].avg_travel_time - g_link_vector[i].VDF_period[tau].FFTT_at[0]) / 60.0;
			float PDT_vc = max(0.0, person_volume * (g_link_vector[i].VDF_period[tau].avg_travel_time - g_link_vector[i].VDF_period[tau].FFTT_at[0] * g_link_vector[i].free_speed / max(0.001f, g_link_vector[i].v_congestion_cutoff)) / 60.0);

			for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
			{
				g_link_vector[i].recorded_volume_per_mode_type_per_period[tau][at][assignment.active_scenario_index] = g_link_vector[i].volume_per_mode_type_per_period[tau][at];
				g_link_vector[i].recorded_converted_MEU_volume_per_period_per_at[tau][at][assignment.active_scenario_index] = g_link_vector[i].converted_MEU_volume_per_period_per_at[tau][at];
			}


			if (base_case_flag == 0)
			{  // base case
				g_link_vector[i].VDF_period[tau].volume_before_sa = vehicle_volume;
				g_link_vector[i].VDF_period[tau].speed_before_sa = speed;
				g_link_vector[i].VDF_period[tau].DoC_before_sa = g_link_vector[i].VDF_period[tau].DOC_mode[0];
				g_link_vector[i].VDF_period[tau].P_before_sa = g_link_vector[i].VDF_period[tau].P;

			}
			else  /// SA case
			{
				g_link_vector[i].VDF_period[tau].volume_after_sa = vehicle_volume;
				g_link_vector[i].VDF_period[tau].speed_after_sa = speed;
				g_link_vector[i].VDF_period[tau].DoC_after_sa = g_link_vector[i].VDF_period[tau].DOC_mode[0];
				g_link_vector[i].VDF_period[tau].P_after_sa = g_link_vector[i].VDF_period[tau].P;

			}

			if (g_link_vector[i].tmc_corridor_name.size() > 0 || g_link_vector[i].VDF_period[tau].network_design_flag == 1)
			{  // with corridor name
				CPeriod_Corridor l_element;
				l_element.volume = vehicle_volume;
				l_element.speed = speed;
				l_element.DoC = g_link_vector[i].VDF_period[tau].DOC_mode[0];
				l_element.P = g_link_vector[i].VDF_period[tau].P;

				string key = g_link_vector[i].tmc_corridor_name;
				if (g_link_vector[i].VDF_period[tau].network_design_flag != 0)
				{
					int idebug = 1;
				}
				if (base_case_flag == 0)
				{  // base case
					g_corridor_info_base0_map[key].record_link_2_corridor_data(l_element, tau);
					if (g_link_vector[i].VDF_period[tau].network_design_flag != 0)
					{
						key = g_link_vector[i].VDF_period[tau].scenario_code + " link_id:" + g_link_vector[i].link_id;
						g_corridor_info_base0_map[key].record_link_2_corridor_data(l_element, tau);
					}
				}
				else  /// SA case
				{
					g_corridor_info_SA_map[key].record_link_2_corridor_data(l_element, tau);

					if (g_link_vector[i].VDF_period[tau].network_design_flag != 0)
					{
						key = g_link_vector[i].VDF_period[tau].scenario_code + " link_id:" + g_link_vector[i].link_id;
						g_corridor_info_SA_map[key].record_link_2_corridor_data(l_element, tau);
					}
				}

			}

			int mode_type_size = assignment.g_ModeTypeVector.size();
			for (int ad = 0; ad < assignment.g_number_of_analysis_districts; ++ad)
				for (int at = 0; at < mode_type_size; ++at)
				{

					CModeType_District l_element;


					//l_element.total_person_travel_time = g_link_vector[i].person_volume_per_district_per_at[ad][at] * g_link_vector[i].VDF_period[at].avg_travel_time;
					//l_element.total_person_distance_km = g_link_vector[i].person_volume_per_district_per_at[ad][at] * g_link_vector[i].link_distance_km;
					//l_element.total_person_distance_mile = g_link_vector[i].person_volume_per_district_per_at[ad][at] * g_link_vector[i].link_distance_mile;

					// district
					if (base_case_flag == 0)
						g_district_info_base0_map[ad].record_link_2_district_data(l_element, at);
					else
						g_district_info_SA_map[ad].record_link_2_district_data(l_element, at);

				}


		}
	}
}

void g_OutputSummaryFiles(Assignment& assignment)
{
	// if 
	if (g_corridor_info_base0_map.size() > 0)
	{

		assignment.summary_corridor_file.open("corridor_performance.csv", std::fstream::out);

		assignment.summary_corridor_file << "tmc_corridor_name,demand_period,links,vol0,speed0,DoC0,max_P0,vol,speed,Doc,P,diff_vol,diff_spd,diff_Doc,diff_P,d%_vol,d%_spd," << endl;

		std::map<string, CCorridorInfo>::iterator it;
		for (it = g_corridor_info_base0_map.begin(); it != g_corridor_info_base0_map.end(); ++it)
		{
			assignment.summary_corridor_file << it->first.c_str() << ",";
			for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); tau++)
			{
				if (assignment.g_DemandPeriodVector[tau].number_of_demand_files >= 1)
				{
					it->second.computer_avg_value(tau);

					assignment.summary_corridor_file <<
						assignment.g_DemandPeriodVector[tau].demand_period.c_str() << "," <<
						it->second.corridor_period[tau].count << "," <<
						it->second.corridor_period[tau].volume << "," <<
						it->second.corridor_period[tau].speed << "," <<
						it->second.corridor_period[tau].DoC << "," <<
						it->second.corridor_period[tau].P << ",";


					if (g_corridor_info_SA_map.find(it->first) != g_corridor_info_SA_map.end())
					{
						g_corridor_info_SA_map[it->first].computer_avg_value(tau);

						CCorridorInfo corridor1 = g_corridor_info_SA_map[it->first];

						assignment.summary_corridor_file <<
							corridor1.corridor_period[tau].volume << "," <<
							corridor1.corridor_period[tau].speed << "," <<
							corridor1.corridor_period[tau].DoC << "," <<
							corridor1.corridor_period[tau].P << "," <<
							corridor1.corridor_period[tau].volume - it->second.corridor_period[tau].volume << "," <<
							corridor1.corridor_period[tau].speed - it->second.corridor_period[tau].speed << "," <<
							corridor1.corridor_period[tau].DoC - it->second.corridor_period[tau].DoC << "," <<
							corridor1.corridor_period[tau].P - it->second.corridor_period[tau].P << "," <<
							(corridor1.corridor_period[tau].volume - it->second.corridor_period[tau].volume) / max(1.0, it->second.corridor_period[tau].volume) * 100.0 << "," <<
							(corridor1.corridor_period[tau].speed - it->second.corridor_period[tau].speed) / max(1.0, it->second.corridor_period[tau].speed) * 100.0 << ",";

					}

					assignment.summary_corridor_file << endl;
				}
			}
		}
	}
	///  //output analysis_district performance file

	if (g_district_info_base0_map.size() >= 2)
	{
		assignment.summary_district_file.open("analysis_district_performance.csv", std::fstream::out);


		assignment.summary_district_file << "district_id,mode_type,od_volume,number_of_links,total_distance_km,total_distance_mile,total_travel_time_min,avg_distance_km,avg_distance_mile,avg_time,SA_total_travel_time_min,SA_avg_distance_km_,SA_avg_distance_mile,SA_avg_time,SA_distance_diff_perc,SA_time_diff_perc,geometry," << endl;
		int mode_type_size = assignment.g_ModeTypeVector.size();
		std::map<int, CAnalysisDistrict>::iterator it;
		for (it = g_district_info_base0_map.begin(); it != g_district_info_base0_map.end(); ++it)
			for (int at = 0; at < mode_type_size; ++at)
			{
				assignment.summary_district_file << it->first << ",";
				assignment.summary_district_file << assignment.g_ModeTypeVector[at].mode_type.c_str() << ",";
				assignment.summary_district_file << it->second.data_by_mode_type[at].total_od_volume << ",";

				{
					it->second.computer_avg_value(at);

					assignment.summary_district_file <<
						it->second.data_by_mode_type[at].count_of_links << "," <<
						it->second.data_by_mode_type[at].total_person_distance_km << "," <<
						it->second.data_by_mode_type[at].total_person_distance_mile << "," <<
						it->second.data_by_mode_type[at].total_person_travel_time << "," <<
						it->second.data_by_mode_type[at].avg_travel_distance_km << "," <<
						it->second.data_by_mode_type[at].avg_travel_distance_mile << "," <<
						it->second.data_by_mode_type[at].avg_travel_time << ",";


					if (g_district_info_SA_map.find(it->first) != g_district_info_SA_map.end())
					{
						g_district_info_SA_map[it->first].data_by_mode_type[at].total_od_volume = it->second.data_by_mode_type[at].total_od_volume;
						g_district_info_SA_map[it->first].computer_avg_value(at);

						CAnalysisDistrict dist_sa = g_district_info_SA_map[it->first];
						assignment.summary_district_file <<
							dist_sa.data_by_mode_type[at].total_person_travel_time << "," <<
							dist_sa.data_by_mode_type[at].avg_travel_distance_km << "," <<
							dist_sa.data_by_mode_type[at].avg_travel_distance_mile << "," <<
							dist_sa.data_by_mode_type[at].avg_travel_time << "," <<
							(dist_sa.data_by_mode_type[at].avg_travel_distance_mile - it->second.data_by_mode_type[at].avg_travel_distance_mile) / max(0.01, it->second.data_by_mode_type[at].avg_travel_distance_mile) * 100.0 << "," <<
							(dist_sa.data_by_mode_type[at].avg_travel_time - it->second.data_by_mode_type[at].avg_travel_time) / max(0.01, it->second.data_by_mode_type[at].avg_travel_time) * 100.0 << ",";
					}
					else
					{
						assignment.summary_district_file << ",,,,,,";

					}

					assignment.summary_district_file << "\"POLYGON ((";
					for (int i = 0; i < it->second.shape_points.size(); i++)
					{
						assignment.summary_district_file << it->second.shape_points[i].x << " " <<
							it->second.shape_points[i].y << ",";

					}

					assignment.summary_district_file << "))\"";

					assignment.summary_district_file << endl;
				}
			}
	}  // end of district

}

// FILE* g_pFileOutputLog = nullptr;

void g_output_dynamic_queue_profile()  // generated from VDF, from numerical queue evolution calculation
{
	dtalog.output() << "writing link_queue_profile.csv.." << endl;

	int b_debug_detail_flag = 0;
	FILE* g_pFileLinkMOE = nullptr;

	string file_name = "link_queue_profile.csv";

	fopen_ss(&g_pFileLinkMOE, file_name.c_str(), "w");

	if (!g_pFileLinkMOE)
	{
		dtalog.output() << "File " << file_name.c_str() << " cannot be opened." << endl;
		g_program_stop();
	}
	else
	{

		// Option 2: BPR-X function
		fprintf(g_pFileLinkMOE, "link_id,tmc_corridor_name,tmc_road_sequence,tmc,link_type_name,from_node_id,to_node_id,geometry,");

		fprintf(g_pFileLinkMOE, "link_type_code,FT,AT,nlanes,link_distance_VDF,free_speed,capacity,k_critical,v_congestion_cutoff,");
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
				double person_volume = g_link_vector[i].total_person_volume_for_all_mode_types_per_period[period_index] + g_link_vector[i].VDF_period[period_index].preload + g_link_vector[i].VDF_period[period_index].sa_volume;

				assignment_VMT = g_link_vector[i].link_distance_VDF * vehicle_volume;
				assignment_VHT = g_link_vector[i].travel_time_per_period[period_index][0] * vehicle_volume / 60.0;  // 60.0 converts min to hour

				assignment_PMT = g_link_vector[i].link_distance_VDF * person_volume;
				assignment_PHT = g_link_vector[i].travel_time_per_period[period_index][0] * person_volume / 60.0;  // 60.0 converts min to hour

				assignment_PSDT = (g_link_vector[i].travel_time_per_period[period_index][0] - g_link_vector[i].free_flow_travel_time_in_min) * person_volume / 60.0;  // 60.0 converts min to hour

				double VCTT = g_link_vector[i].link_distance_VDF / max(1.0f, g_link_vector[i].v_congestion_cutoff) * 60;
				assignment_VCDT = max(0.0, g_link_vector[i].travel_time_per_period[period_index][0] - VCTT) * g_link_vector[i].total_volume_for_all_mode_types_per_period[period_index] / 60.0;  // 60.0 converts min to hour


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

void g_output_assignment_result(Assignment& assignment, int subarea_id)
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
		assignment.summary_file << "Output Link Performance:" << endl;

		char link_performance_file_name[50];
		int scenario_index_no = assignment.g_active_DTAscenario_map[assignment.active_scenario_index];
		sprintf(link_performance_file_name, "link_performance_s%d_%s.csv", assignment.active_scenario_index, assignment.g_DTAscenario_vector[scenario_index_no].scenario_name.c_str());

		dtalog.output() << "writing link_performance_s" << assignment.active_scenario_index << ".csv" << endl;
		fopen_ss(&g_pFileLinkMOE, link_performance_file_name, "w");
		if (!g_pFileLinkMOE)
		{
			dtalog.output() << "File link_performance_s" << assignment.active_scenario_index << ".csv" << "cannot be opened." << endl;
			g_program_stop();
		}


	}
	else
	{
		assignment.summary_file << "Output subarea Link Performance:" << endl;

		dtalog.output() << "writing subarea_link_performance.csv.." << endl;
		fopen_ss(&g_pFileLinkMOE, "subarea_link_performance.csv", "w");
		if (!g_pFileLinkMOE)
		{
			dtalog.output() << "File subarea_link_performance.csv cannot be opened." << endl;
			g_program_stop();
		}
	}


	int ref_volume_count = 0;
	float total_ref_volume_dev_abs_percentage = 0;

	{

		fprintf(g_pFileLinkMOE, "link_id,vdf_type,from_node_id,to_node_id,lanes,distance_km,distance_mile,fftt,meso_link_id,meso_link_incoming_volume,tmc,tmc_corridor_name,tmc_corridor_id,tmc_road_order,tmc_road_sequence,subarea_id,link_type,link_type_code,vdf_code,time_period,volume,ref_volume,ref_diff,volume_before_odme,volume_after_odme,volume_diff_odme,obs_count,upper_bound_type,obs_count_dev_odme,");



		fprintf(g_pFileLinkMOE, "preload_volume,person_volume,travel_time,speed_kmph,speed_mph,speed_ratio,VOC,DOC,capacity,queue,total_simu_waiting_time_in_min,avg_simu_waiting_time_in_min,plf,lanes,D_per_hour_per_lane,QVDF_cd,QVDF_n,P,severe_congestion_duration_in_h,vf,v_congestion_cutoff,QVDF_cp,QVDF_s,QVDF_v,vt2,VMT,VHT,PMT,PHT,PDT_vf,PDT_vc,geometry,");

		for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
			fprintf(g_pFileLinkMOE, "person_vol_%s,", assignment.g_ModeTypeVector[at].mode_type.c_str());

		for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
			fprintf(g_pFileLinkMOE, "MEU_vol_%s,", assignment.g_ModeTypeVector[at].mode_type.c_str());

		for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
			fprintf(g_pFileLinkMOE, "mode_cap_%s,", assignment.g_ModeTypeVector[at].mode_type.c_str());

		for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
			fprintf(g_pFileLinkMOE, "mode_VOC_%s,", assignment.g_ModeTypeVector[at].mode_type.c_str());

		//for (int og = 0; og < assignment.g_number_of_analysis_districts; ++og)
		//	for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
		//	{
		//		fprintf(g_pFileLinkMOE, "person_vol_district_%d_%s,", og, assignment.g_ModeTypeVector[at].mode_type.c_str());

		//	}


//		fprintf(g_pFileLinkMOE, "obs_count,upper_bound_type,dev,");

		for (int t = 6 * 60; t < 20 * 60; t += 15)
		{
			int hour = t / 60;
			int minute = t - hour * 60;

			fprintf(g_pFileLinkMOE, "v%02d:%02d,", hour, minute);
		}

		fprintf(g_pFileLinkMOE, "scenario_code,volume_before_sa,volume_after_sa,volume_diff_sa,speed_before_sa,speed_after_sa,speed_diff_sa,DoC_before_sa,DoC_after_sa,Doc_diff_sa,P_before_sa,P_after_sa,P_diff_sa,");


		fprintf(g_pFileLinkMOE, "notes\n");

		//Initialization for all nodes
		for (int i = 0; i < g_link_vector.size(); ++i)
		{
			// virtual connectors
			if (g_link_vector[i].link_type_si[0] == -1)
				continue;

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
				total_vehicle_volume += g_link_vector[i].VDF_period[tau].volume_before_sa;
				total_vehicle_volume += g_link_vector[i].total_volume_for_all_mode_types_per_period[tau] + g_link_vector[i].VDF_period[tau].preload;
				scenario_code_count += g_link_vector[i].VDF_period[tau].scenario_code.size();

			}

			//if (scenario_code_count == 0 && total_vehicle_volume < 0.001 && g_link_vector.size() > 20000)
			//	continue;

			for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
			{

				if (assignment.g_DemandPeriodVector[tau].number_of_demand_files == 0)
					continue;



				float speed = g_link_vector[i].free_speed;  // default speed 


				if (g_link_vector[i].VDF_period[tau].avg_travel_time > 0.001f)
					speed = g_link_vector[i].free_speed / max(0.000001, g_link_vector[i].VDF_period[tau].avg_travel_time) * g_link_vector[i].VDF_period[tau].FFTT_at[0];

				float speed_ratio = speed / max(1.0, g_link_vector[i].free_speed);  // default speed 
				float vehicle_volume = g_link_vector[i].total_volume_for_all_mode_types_per_period[tau] + g_link_vector[i].VDF_period[tau].preload;
				float ref_volume_diff = 0;

				if (g_link_vector[i].VDF_period[tau].ref_link_volume > 1)
				{
					ref_volume_diff = vehicle_volume - g_link_vector[i].VDF_period[tau].ref_link_volume;

					ref_volume_count++;
					total_ref_volume_dev_abs_percentage += fabs(ref_volume_diff / g_link_vector[i].VDF_period[tau].ref_link_volume * 100);
				}


				float person_volume = g_link_vector[i].total_person_volume_for_all_mode_types_per_period[tau] + g_link_vector[i].VDF_period[tau].preload;
				//VMT,VHT,PMT,PHT,PDT
				float preload = g_link_vector[i].VDF_period[tau].preload;
				float VMT = vehicle_volume * g_link_vector[i].link_distance_mile;

				float VHT = vehicle_volume * g_link_vector[i].VDF_period[tau].avg_travel_time / 60.0;
				float PMT = person_volume * g_link_vector[i].link_distance_mile;
				float PHT = person_volume * g_link_vector[i].VDF_period[tau].avg_travel_time / 60.0;
				float PDT_vf = person_volume * (g_link_vector[i].VDF_period[tau].avg_travel_time - g_link_vector[i].VDF_period[tau].FFTT_at[0]) / 60.0;
				float PDT_vc = max(0.0, person_volume * (g_link_vector[i].VDF_period[tau].avg_travel_time - g_link_vector[i].VDF_period[tau].FFTT_at[0] * g_link_vector[i].free_speed / max(0.001f, g_link_vector[i].v_congestion_cutoff)) / 60.0);



				total_VMT += VMT;
				total_VHT += VHT;
				total_PMT += PMT;
				total_PHT += PHT;
				total_link_volume += vehicle_volume;
				total_link_speed += speed;
				total_link_speed_ratio += speed_ratio;

				MOE_count++;

				fprintf(g_pFileLinkMOE, "%s,%s,%d,%d,%f,%f,%f,%f,%d,%d,",

					//end with obs_volume_diff
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

				fprintf(g_pFileLinkMOE, "%.3f,%.3f,%.3f,",
					vehicle_volume,
					g_link_vector[i].VDF_period[tau].ref_link_volume,
					ref_volume_diff);

				fprintf(g_pFileLinkMOE, "%.3f,%.3f,%.3f,",
					g_link_vector[i].VDF_period[tau].volume_before_odme,
					g_link_vector[i].VDF_period[tau].volume_after_odme,
					g_link_vector[i].VDF_period[tau].volume_after_odme - g_link_vector[i].VDF_period[tau].volume_before_odme);

				if (g_link_vector[i].VDF_period[tau].obs_count >= 1)
				{
					fprintf(g_pFileLinkMOE, "%.3f,%.3f,%.3f,",
						g_link_vector[i].VDF_period[tau].obs_count,
						g_link_vector[i].VDF_period[tau].upper_bound_flag,
						g_link_vector[i].VDF_period[tau].obs_count - g_link_vector[i].VDF_period[tau].volume_after_odme);
				}
				else
				{
					fprintf(g_pFileLinkMOE, ",,,");
				}


				fprintf(g_pFileLinkMOE, "%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,",
					preload,
					person_volume,
					g_link_vector[i].VDF_period[tau].avg_travel_time,
					speed,  /* 60.0 is used to convert min to hour */
					speed / 1.609,  /* 60.0 is used to convert min to hour */
					speed_ratio,
					g_link_vector[i].VDF_period[tau].DOC_mode[0],
					g_link_vector[i].VDF_period[tau].DOC_mode[0],
					g_link_vector[i].VDF_period[tau].lane_based_ultimate_hourly_capacity,
					g_link_vector[i].VDF_period[tau].queue_length,
					max(0, g_link_vector[i].total_simulated_delay_in_min));


				fprintf(g_pFileLinkMOE, "%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,",

					max(0, g_link_vector[i].total_simulated_delay_in_min) / max(1, vehicle_volume),
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
					fprintf(g_pFileLinkMOE, "%.3f,", g_link_vector[i].VDF_period[tau].DOC_mode[at]); // no complete

				for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
					fprintf(g_pFileLinkMOE, "%.3f,", g_link_vector[i].VDF_period[tau].DOC_mode[at]);

				//for (int og = 0; og < assignment.g_number_of_analysis_districts; ++og)
				//	for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
				//	{
				//		fprintf(g_pFileLinkMOE, "%.3f,", g_link_vector[i].person_volume_per_district_per_at[og][at]);
				//	}


				for (int t = 6 * 60; t < 20 * 60; t += 15)
				{
					float speed = g_link_vector[i].get_model_15_min_speed(t);
					fprintf(g_pFileLinkMOE, "%.3f,", speed);
				}

				fprintf(g_pFileLinkMOE, "%s,", g_link_vector[i].VDF_period[tau].scenario_code.c_str());

				fprintf(g_pFileLinkMOE, "%.3f,%.3f,%.3f,", g_link_vector[i].VDF_period[tau].volume_before_sa,
					g_link_vector[i].VDF_period[tau].volume_after_sa,
					g_link_vector[i].VDF_period[tau].volume_after_sa - g_link_vector[i].VDF_period[tau].volume_before_sa);

				fprintf(g_pFileLinkMOE, "%.3f,%.3f,%.3f,", g_link_vector[i].VDF_period[tau].speed_before_sa,
					g_link_vector[i].VDF_period[tau].speed_after_sa,
					g_link_vector[i].VDF_period[tau].speed_after_sa - g_link_vector[i].VDF_period[tau].speed_before_sa);

				fprintf(g_pFileLinkMOE, "%.3f,%.3f,%.3f,", g_link_vector[i].VDF_period[tau].DoC_before_sa,
					g_link_vector[i].VDF_period[tau].DoC_after_sa,
					g_link_vector[i].VDF_period[tau].DoC_after_sa - g_link_vector[i].VDF_period[tau].DoC_before_sa);


				fprintf(g_pFileLinkMOE, "%.3f,%.3f,%.3f,", g_link_vector[i].VDF_period[tau].P_before_sa,
					g_link_vector[i].VDF_period[tau].P_after_sa,
					g_link_vector[i].VDF_period[tau].P_after_sa - g_link_vector[i].VDF_period[tau].P_before_sa);


				fprintf(g_pFileLinkMOE, "period-based\n");
			}

		}
		fclose(g_pFileLinkMOE);
	}

	assignment.summary_file << ",ref_link_vol_count=," << ref_volume_count << "," << "MAPE=," << total_ref_volume_dev_abs_percentage / max(1, ref_volume_count) << "%" << endl;

	assignment.summary_file << ",VMT=," << total_VMT << "," << "VKT=," << total_VMT * 1.6090 << endl;
	assignment.summary_file << ",VHT=," << total_VHT << endl;
	assignment.summary_file << ",network vehicle speed (MPH) =," << total_VMT / max(0.0001, total_VHT) << ",";
	assignment.summary_file << ",network vehicle speed (KPH) =," << total_VMT * 1.609 / max(0.0001, total_VHT) << endl;

	assignment.summary_file << ",PMT=," << total_PMT << ",";
	assignment.summary_file << "PKT=," << total_PMT * 1.609 << endl;
	assignment.summary_file << ",PHT=," << total_PHT << endl;
	assignment.summary_file << ",network person speed (MPH) =," << total_PMT / max(0.0001, total_PHT) << ",";
	assignment.summary_file << ",network person speed (KPH) =," << total_PMT * 1.609 / max(0.0001, total_PHT) << endl;


	assignment.summary_file << ",simple avg link volume=," << total_link_volume / max(1, MOE_count) << endl;
	assignment.summary_file << ",simple avg link speed=," << total_link_speed / max(1, MOE_count) << endl;
	assignment.summary_file << ",simple avg link speed ratio=," << total_link_speed_ratio / max(1, MOE_count) << endl;




	if (subarea_id == 1)
		return;

	{

		double path_time_vector[MAX_LINK_SIZE_IN_A_PATH];
		char route_assignment_file_name[50];

		int scenario_index_no = assignment.g_active_DTAscenario_map[assignment.active_scenario_index];
		sprintf(route_assignment_file_name, "route_assignment_s%d_%s.csv", assignment.active_scenario_index, assignment.g_DTAscenario_vector[scenario_index_no].scenario_name.c_str());



		FILE* g_pFilePathMOE = nullptr;
		fopen_ss(&g_pFilePathMOE, route_assignment_file_name, "w");
		dtalog.output() << "writing route_assignment.csv.." << endl;

		if (!g_pFilePathMOE)
		{
			dtalog.output() << "File " << route_assignment_file_name << " cannot be opened." << endl;
			g_program_stop();
		}

		fprintf(g_pFilePathMOE, "first_column,path_no,o_zone_id,d_zone_id,od_pair,o_sindex,d_sindex,within_OD_path_no,path_id,activity_zone_sequence,activity_mode_type_sequence,information_type,mode_type,demand_period,volume,OD_relative_gap,travel_time,path_gap,preload_volume,volume_before_odme,volume_after_odme,volume_diff_odme,rt_new_path_flag,volume_before_sa,volume_after_sa,volume_diff_sa,simu_volume,subarea_flag,OD_impact_flag,at_OD_impact_flag,");
		fprintf(g_pFilePathMOE, "path_impact_flag,toll,#_of_nodes,#_of_sensor_links,#_of_SA_links,travel_time,VDF_travel_time,VDF_travel_time_without_access_link,distance_km,distance_mile,node_sequence,link_sequence, ");

		//// stage 1: column updating
		//for (int iteration_number = 0; iteration_number < min(20, assignment.g_number_of_column_updating_iterations); iteration_number++)
		//{ 
		//    fprintf(g_pFilePathMOE, "TT_%d,", iteration_number);
		//}

		//for (int iteration_number = 0; iteration_number < min(20, assignment.g_number_of_column_updating_iterations); iteration_number++)
		//{
		//    fprintf(g_pFilePathMOE, "Vol_%d,", iteration_number);
		//}

		//stage 2: ODME
		for (int iteration_number = 0; iteration_number < min(20, assignment.g_number_of_ODME_iterations); iteration_number++)
		{
			fprintf(g_pFilePathMOE, "ODME_TT_%d,", iteration_number);
		}

		for (int iteration_number = max(0, assignment.g_number_of_ODME_iterations - 10); iteration_number < assignment.g_number_of_ODME_iterations; iteration_number++)
		{
			fprintf(g_pFilePathMOE, "ODME_Vol_%d,", iteration_number);
		}

		//stage 3: sensitivity analysi
		for (int iteration_number = 0; iteration_number <= assignment.g_number_of_sensitivity_analysis_iterations; iteration_number++)
		{
			fprintf(g_pFilePathMOE, "SA_TT_%d,", iteration_number);
		}
		for (int iteration_number = 0; iteration_number <= assignment.g_number_of_sensitivity_analysis_iterations; iteration_number++)
		{
			if (iteration_number == 0)
				fprintf(g_pFilePathMOE, "SA_Vol_%d,", iteration_number);
			else
				fprintf(g_pFilePathMOE, "SA_Delta_Vol_%d,", iteration_number);
		}
		fprintf(g_pFilePathMOE, "geometry,");

		fprintf(g_pFilePathMOE, "link_type_name_sequence,link_code_sequence,link_link_distance_VDF_sequence,link_FFTT_sequence,");
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
		double path_distance_km = 0;
		double path_distance_ml = 0;
		float path_travel_time = 0;
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

			/// output background_link_volume.csv
			dtalog.output() << "writing link_background_volume.csv.." << endl;

			int b_debug_detail_flag = 0;
			FILE* g_pFileLinkMOE = nullptr;
			int b_background_link_volume_file = 0;

			if (b_background_link_volume_file)
			{
				fopen_ss(&g_pFileLinkMOE, "link_background_volume.csv", "w");
				if (!g_pFileLinkMOE)
				{
					dtalog.output() << "File link_background_volume.csv cannot be opened." << endl;
					g_program_stop();
				}
				else
				{
					fprintf(g_pFileLinkMOE, "link_id,from_node_id,to_node_id,subarea_id,time_period,volume,background_volume,major_path_volume,ratio_of_major_path_flow,geometry,");

					fprintf(g_pFileLinkMOE, "notes\n");

					//Initialization for all nodes
					for (int i = 0; i < g_link_vector.size(); ++i)
					{
						//   virtual connectors
						if (g_link_vector[i].link_type_si[0] == -1)
							continue;

						for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
						{
							double volume = g_link_vector[i].total_volume_for_all_mode_types_per_period[tau] + g_link_vector[i].VDF_period[tau].preload;
							double major_path_link_volume = g_link_vector[i].total_volume_for_all_mode_types_per_period[tau] + g_link_vector[i].VDF_period[tau].preload - g_link_vector[i].background_total_volume_for_all_mode_types_per_period[tau];
							double ratio = major_path_link_volume / max(volume, 0.000001);

							if (volume < 0.0000001)
								ratio = -1;
							fprintf(g_pFileLinkMOE, "%s,%d,%d,%d,%s,%.3f,%.3f,%.3f,%.3f,\"%s\",",
								g_link_vector[i].link_id.c_str(),
								g_node_vector[g_link_vector[i].from_node_seq_no].node_id,
								g_node_vector[g_link_vector[i].to_node_seq_no].node_id,
								g_link_vector[i].subarea_id,
								assignment.g_DemandPeriodVector[tau].time_period.c_str(),
								g_link_vector[i].total_volume_for_all_mode_types_per_period[tau] + g_link_vector[i].VDF_period[tau].preload,
								g_link_vector[i].background_total_volume_for_all_mode_types_per_period[tau],
								major_path_link_volume,
								ratio,
								g_link_vector[i].geometry.c_str());
							fprintf(g_pFileLinkMOE, "\n");

						}

					}

					fclose(g_pFileLinkMOE);
				}

			}


		} // end of path flow pattern screening 
		dtalog.output() << "writing data for " << zone_size << "  zones " << endl;





		for (int orig = 0; orig < zone_size; ++orig)
		{
			if (g_zone_vector[orig].zone_id % 100 == 0)
				dtalog.output() << "o zone id =  " << g_zone_vector[orig].zone_id << endl;

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

						if (p_column_pool->OD_impact_flag == 1 && assignment.g_ModeTypeVector[at].real_time_information == 0)  // as long as regular users of this OD pair is impacted. 
						{
							global_od_impact_flag_across_all_mode_types = 1;
						}
					}

					for (int at = 0; at < mode_type_size; ++at)
					{
						p_column_pool = &(assignment.g_column_pool[from_zone_sindex][to_zone_sindex][at][tau]);
						if (p_column_pool->od_volume[assignment.active_scenario_index] > 0)
						{

							int at_OD_impact_flag;
							if (p_column_pool->at_od_impacted_flag_map.find(at) != p_column_pool->at_od_impacted_flag_map.end())
							{
								at_OD_impact_flag = 1;  // marked
							}
							else
							{
								at_OD_impact_flag = 0;
							}


							if (g_zone_vector[orig].zone_id == 4 && g_zone_vector[dest].zone_id == 2 && at == 0)
							{
								int idebug = 1;
							}

							int information_type = p_column_pool->information_type;


							time_stamp = (assignment.g_DemandPeriodVector[tau].starting_time_slot_no + assignment.g_DemandPeriodVector[tau].ending_time_slot_no) / 2.0 * MIN_PER_TIMESLOT;

							// scan through the map with different node sum for different continuous paths
							it_begin = p_column_pool->path_node_sequence_map.begin();
							it_end = p_column_pool->path_node_sequence_map.end();

							for (it = it_begin; it != it_end; ++it)
							{

								if (count % 100000 == 0)
								{
									end_t = clock();
									iteration_t = end_t - start_t;
									dtalog.output() << "writing " << count / 1000 << "K agents with CPU time " << iteration_t / 1000.0 << " s" << endl;
								}

								path_toll = 0;
								path_distance = 0;
								path_distance_km = 0;
								path_distance_ml = 0;
								path_travel_time = 0;
								path_travel_time_without_access_link = 0;
								path_FF_travel_time = 0;

								path_time_vector[0] = time_stamp;
								path_time_vector[1] = time_stamp;

								string link_code_str;


								if (it->second.m_link_size >= MAX_LINK_SIZE_IN_A_PATH - 2)
								{
									dtalog.output() << "error: it->second.m_link_size < MAX_LINK_SIZE_IN_A_PATH" << endl;
									dtalog.output() << "o= " << g_zone_vector[orig].zone_id << endl;
									dtalog.output() << "d= " << g_zone_vector[dest].zone_id << endl;
									dtalog.output() << "agent type= " << assignment.g_ModeTypeVector[at].mode_type.c_str() << endl;

									for (int nl = 0; nl < it->second.m_link_size; ++nl)  // arc a
									{
										dtalog.output() << "node no." << nl << " =" << g_node_vector[it->second.path_node_vector[nl]].node_id << endl;
									}

									g_program_stop();
								}

								for (int nl = 0; nl < it->second.m_link_size; ++nl)  // arc a
								{

									int link_seq_no = it->second.path_link_vector[nl];
									if (g_link_vector[link_seq_no].link_type_si[assignment.active_scenario_index] >= 0)
									{
										path_toll += g_link_vector[link_seq_no].VDF_period[tau].toll[at];
										path_distance += g_link_vector[link_seq_no].link_distance_VDF;
										path_distance_km += g_link_vector[link_seq_no].link_distance_km;
										path_distance_ml += g_link_vector[link_seq_no].link_distance_mile;
										float link_travel_time = g_link_vector[link_seq_no].travel_time_per_period[tau][at];
										path_travel_time += link_travel_time;

										if (g_link_vector[link_seq_no].link_type_si[assignment.active_scenario_index] != 1000)   // skip access link
										{
											path_travel_time_without_access_link += link_travel_time;
										}

										path_FF_travel_time += g_link_vector[link_seq_no].VDF_period[tau].FFTT_at[at];
										time_stamp += link_travel_time;
										path_time_vector[nl + 1] = time_stamp;
										link_code_str += g_link_vector[link_seq_no].link_code_str;
									}
									else
									{
										int virtual_link = 0;
									}
								}

								double total_agent_path_travel_time = 0;

								for (int vi = 0; vi < it->second.agent_simu_id_vector.size(); ++vi)
								{
									int agent_simu_id = it->second.agent_simu_id_vector[vi];
									CAgent_Simu* pAgentSimu = g_agent_simu_vector[agent_simu_id];
									total_agent_path_travel_time += pAgentSimu->path_travel_time_in_min;
								}

								double final_path_travel_time = path_travel_time;  // by default

								if (it->second.agent_simu_id_vector.size() > 1)  // with simulated agents
								{
									final_path_travel_time = total_agent_path_travel_time / it->second.agent_simu_id_vector.size();
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
									if (it->second.m_node_size - virtual_first_link_delta - virtual_last_link_delta <= 2)
										continue;

									int impacted_path_flag = it->second.impacted_path_flag;  //NA by default


									double volume_before_ODME = max(0, it->second.path_volume_before_ODME);
									double volume_after_ODME = max(0, it->second.path_volume_after_ODME);
									double volume = it->second.path_volume;

									double volume_diff_ODME = 0;

									if (volume_before_ODME >= -0.000001)
									{
										volume_diff_ODME = volume_after_ODME - volume_before_ODME;
									}

									double volume_before_sa = max(0, it->second.path_volume_before_sa);
									double volume_after_sa = max(0, it->second.path_volume_after_sa);
									double volume_diff_sa = 0;

									if (volume_before_sa >= -0.000001)
									{
										volume_diff_sa = volume_after_sa - volume_before_sa;
									}


									fprintf(g_pFilePathMOE, ",%d,%d,%d,%d,%d,%d->%d,%d,%d,%s,%s,%d,%s,%s,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%d,%d,%d,%d,%d,%.1f,%d,%d,%d,%.1f,%.4f,%.4f,%.4f,%.4f,",
										count,
										g_zone_vector[orig].zone_id,
										g_zone_vector[dest].zone_id,
										g_zone_vector[orig].sindex,
										g_zone_vector[dest].sindex,
										g_zone_vector[orig].zone_id,
										g_zone_vector[dest].zone_id,
										it->second.path_seq_no,
										it->second.path_seq_no + 1,
										p_column_pool->activity_zone_sequence.c_str(),
										p_column_pool->activity_mode_type_sequence.c_str(),
										information_type,
										assignment.g_ModeTypeVector[at].mode_type.c_str(),
										assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
										volume,
										p_column_pool->relative_OD_gap,
										path_travel_time,
										it->second.path_gradient_cost_relative_difference,
										it->second.path_preload_volume,
										volume_before_ODME,
										volume_after_ODME,
										volume_diff_ODME,
										it->second.b_RT_new_path_flag,
										volume_before_sa,
										volume_after_sa,
										volume_diff_sa,
										it->second.agent_simu_id_vector.size(),
										it->second.subarea_output_flag,
										global_od_impact_flag_across_all_mode_types,
										at_OD_impact_flag,
										impacted_path_flag,
										path_toll,
										it->second.m_node_size - virtual_first_link_delta - virtual_last_link_delta,
										it->second.path_sensor_link_vector.size(),
										it->second.path_SA_link_vector.size(),
										final_path_travel_time,
										path_travel_time,
										path_travel_time_without_access_link,
										//path_FF_travel_time,
										//final_path_travel_time- path_FF_travel_time,
										path_distance_km,
										path_distance_ml);

									/* Format and print various data */
									for (int ni = 0 + virtual_first_link_delta; ni < it->second.m_node_size - virtual_last_link_delta; ++ni)
										fprintf(g_pFilePathMOE, "%d;", g_node_vector[it->second.path_node_vector[ni]].node_id);

									fprintf(g_pFilePathMOE, ",");
									int link_seq_no;
									// link id sequence
									for (int nl = 0 + virtual_first_link_delta; nl < it->second.m_link_size - virtual_last_link_delta; ++nl)
									{
										link_seq_no = it->second.path_link_vector[nl];
										fprintf(g_pFilePathMOE, "%s;", g_link_vector[link_seq_no].link_id.c_str());
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
									for (int iteration_number = 0; iteration_number < min(20, assignment.g_number_of_ODME_iterations); iteration_number++)
									{
										double TT = -1;
										if (it->second.path_time_per_iteration_ODME_map.find(iteration_number) != it->second.path_time_per_iteration_ODME_map.end())
										{
											TT = it->second.path_time_per_iteration_ODME_map[iteration_number];
										}


										fprintf(g_pFilePathMOE, "%f,", TT);

									}
									for (int iteration_number = 0; iteration_number < min(20, assignment.g_number_of_ODME_iterations); iteration_number++)
									{
										double Vol = 0;
										if (it->second.path_volume_per_iteration_ODME_map.find(iteration_number) != it->second.path_volume_per_iteration_ODME_map.end())
										{
											Vol = it->second.path_volume_per_iteration_ODME_map[iteration_number];
										}

										fprintf(g_pFilePathMOE, "%f,", Vol);

									}

									//stage III: 

									// output the TT and vol of sensitivity analysis
									for (int iteration_number = 0; iteration_number <= assignment.g_number_of_sensitivity_analysis_iterations; iteration_number++)
									{
										double TT = -1;
										if (it->second.path_time_per_iteration_SA_map.find(iteration_number) != it->second.path_time_per_iteration_SA_map.end())
										{
											TT = it->second.path_time_per_iteration_SA_map[iteration_number];
										}

										fprintf(g_pFilePathMOE, "%f,", TT);

									}
									double prev_value = 0;
									for (int iteration_number = 0; iteration_number <= assignment.g_number_of_sensitivity_analysis_iterations; iteration_number++)
									{
										double Vol = 0;
										if (it->second.path_volume_per_iteration_SA_map.find(iteration_number) != it->second.path_volume_per_iteration_SA_map.end())
										{
											Vol = it->second.path_volume_per_iteration_SA_map[iteration_number];
										}

										fprintf(g_pFilePathMOE, "%f,", Vol - prev_value);
										prev_value = Vol;

									}


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

									if (it->second.path_volume >= 5)  // critcal path volume
									{
										fprintf(g_pFilePathMOE, ",");
										// link type name sequenece
										//for (int nl = 0 + virtual_first_link_delta; nl < it->second.m_link_size - virtual_last_link_delta; ++nl)
										//{
										//    link_seq_no = it->second.path_link_vector[nl];
										//    fprintf(g_pFilePathMOE, "%s;", g_link_vector[link_seq_no].link_type_name.c_str());
										//}
										fprintf(g_pFilePathMOE, ",");

										// link code sequenece
										for (int nl = 0 + virtual_first_link_delta; nl < it->second.m_link_size - virtual_last_link_delta; ++nl)
										{
											link_seq_no = it->second.path_link_vector[nl];

											if (g_link_vector[link_seq_no].link_code_str.size() > 1)
												fprintf(g_pFilePathMOE, "%s;", g_link_vector[link_seq_no].link_code_str.c_str());
										}
										fprintf(g_pFilePathMOE, ",");

										// link link_distance_VDF sequenece
										for (int nl = 0 + virtual_first_link_delta; nl < it->second.m_link_size - virtual_last_link_delta; ++nl)
										{
											link_seq_no = it->second.path_link_vector[nl];
											fprintf(g_pFilePathMOE, "%.3f;", g_link_vector[link_seq_no].link_distance_VDF);
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
									fprintf(g_pFilePathMOE, "\n");
									count++;
								}

							}
						}
					}
				}
			}
		}
		fclose(g_pFilePathMOE);

	}
	//	g_output_dynamic_queue_profile();
}

void g_output_choice_set_result(Assignment& assignment)
{


	double path_time_vector[MAX_LINK_SIZE_IN_A_PATH];
	char route_assignment_file_name[50];

	int scenario_index_no = assignment.g_active_DTAscenario_map[assignment.active_scenario_index];
	sprintf(route_assignment_file_name, "choice_set_output_%d_%s.csv", assignment.active_scenario_index, assignment.g_DTAscenario_vector[scenario_index_no].scenario_name.c_str());

	FILE* g_pFilePathMOE = nullptr;
	fopen_ss(&g_pFilePathMOE, route_assignment_file_name, "w");
	dtalog.output() << "writing choice_set_output.csv.." << endl;

	if (!g_pFilePathMOE)
	{
		dtalog.output() << "File " << route_assignment_file_name << " cannot be opened." << endl;
		g_program_stop();
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
			fprintf(g_pFilePathMOE, ",");  //multi_dim_choice_id,mode_tag,demand_period_tag,spatial_tag,travel_purpose_tag,data_tag,
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
								path_toll += g_link_vector[link_seq_no].VDF_period[tau].toll[at];
								path_distance += g_link_vector[link_seq_no].link_distance_VDF;
								path_distance_km += g_link_vector[link_seq_no].link_distance_km;
								path_distance_ml += g_link_vector[link_seq_no].link_distance_mile;
								float link_travel_time = g_link_vector[link_seq_no].travel_time_per_period[tau][at];
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
		assignment.summary_file << "Output Link Performance Summary" << endl;

		dtalog.output() << "writing link_performances_summary.csv.." << endl;
		fopen_ss(&g_pFileLinkMOE, "link_performance_summary.csv", "w");
		if (!g_pFileLinkMOE)
		{
			dtalog.output() << "File link_performance_summary.csv cannot be opened." << endl;
			g_program_stop();
		}


	}
	else
	{
		assignment.summary_file << "Output subarea Link Performance:" << endl;

		dtalog.output() << "writing subarea_link_performance.csv.." << endl;
		fopen_ss(&g_pFileLinkMOE, "subarea_link_performance.csv", "w");
		if (!g_pFileLinkMOE)
		{
			dtalog.output() << "File subarea_link_performance.csv cannot be opened." << endl;
			g_program_stop();
		}
	}

	if (g_pFileLinkMOE != NULL)
	{

		int ref_volume_count = 0;
		float total_ref_volume_dev_abs_percentage = 0;


		fprintf(g_pFileLinkMOE, "link_id,based_link_type,from_node_id,to_node_id,geometry,distance_km,distance_mile,fftt,meso_link_id,tmc,tmc_corridor_name,tmc_corridor_id,tmc_road_order,tmc_road_sequence,subarea_id,");

		for (int sii = 0; sii < assignment.g_DTAscenario_vector.size(); sii++)
		{
			int scenario_index = assignment.g_DTAscenario_vector[sii].scenario_index;
			fprintf(g_pFileLinkMOE, "lanes_s%d,", scenario_index);
		}

		for (int sii = 0; sii < assignment.g_DTAscenario_vector.size(); sii++)
		{
			int scenario_index = assignment.g_DTAscenario_vector[sii].scenario_index;
			fprintf(g_pFileLinkMOE, "type_s%d,", scenario_index);
		}

		for (int sii = 0; sii < assignment.g_DTAscenario_vector.size(); sii++)
		{
			int scenario_index = assignment.g_DTAscenario_vector[sii].scenario_index;
			fprintf(g_pFileLinkMOE, "penaltys%d,", scenario_index);
		}


		for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
		{
			for (int sii = 0; sii < assignment.g_DTAscenario_vector.size(); sii++)
			{
				int scenario_index = assignment.g_DTAscenario_vector[sii].scenario_index;
				{

					for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
					{

						fprintf(g_pFileLinkMOE, "p%d_s%d_%s_vol,", tau, scenario_index, assignment.g_ModeTypeVector[at].mode_type.c_str());
					}

				}
			}

		}

		for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); ++tau)
		{
			for (int sii = 0; sii < assignment.g_DTAscenario_vector.size(); sii++)
			{
				int scenario_index = assignment.g_DTAscenario_vector[sii].scenario_index;
				{
					int scenario_index_no = assignment.g_active_DTAscenario_map[scenario_index];

					for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
					{

						fprintf(g_pFileLinkMOE, "%s_s_%s_%s_vol,", assignment.g_DemandPeriodVector[tau].demand_period.c_str(), assignment.g_DTAscenario_vector[scenario_index_no].scenario_name.c_str(), assignment.g_ModeTypeVector[at].mode_type.c_str());
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

			string vdf_type_str;


			float total_vehicle_volume = 0;
			int scenario_code_count = 0;
			for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
			{
				total_vehicle_volume += g_link_vector[i].VDF_period[tau].volume_before_sa;
				total_vehicle_volume += g_link_vector[i].total_volume_for_all_mode_types_per_period[tau] + g_link_vector[i].VDF_period[tau].preload;
				scenario_code_count += g_link_vector[i].VDF_period[tau].scenario_code.size();

			}

			//if (scenario_code_count == 0 && total_vehicle_volume < 0.001 && g_link_vector.size() > 20000)
			//	continue;




			float speed = g_link_vector[i].free_speed;  // default speed 


			fprintf(g_pFileLinkMOE, "%s,%s,%d,%d,",

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

			for (int sii = 0; sii < assignment.g_DTAscenario_vector.size(); sii++)
			{
				int scenario_index = assignment.g_DTAscenario_vector[sii].scenario_index;
				fprintf(g_pFileLinkMOE, "%.1f,", g_link_vector[i].number_of_lanes_si[scenario_index]);

			}

			for (int sii = 0; sii < assignment.g_DTAscenario_vector.size(); sii++)
			{
				int scenario_index = assignment.g_DTAscenario_vector[sii].scenario_index;
				fprintf(g_pFileLinkMOE, "%s,", assignment.g_LinkTypeMap[g_link_vector[i].link_type_si[scenario_index]].link_type_name.c_str());
			}
			for (int sii = 0; sii < assignment.g_DTAscenario_vector.size(); sii++)
			{
				int scenario_index = assignment.g_DTAscenario_vector[sii].scenario_index;
				fprintf(g_pFileLinkMOE, "%.1f,", g_link_vector[i].penalty_si_flag[scenario_index]);
			}

			for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
			{
				for (int sii = 0; sii < assignment.g_DTAscenario_vector.size(); sii++)
				{
					int scenario_index = assignment.g_DTAscenario_vector[sii].scenario_index;
					{

						for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
						{

							fprintf(g_pFileLinkMOE, "%.1f,", g_link_vector[i].recorded_volume_per_mode_type_per_period[tau][at][scenario_index]);

						}
					}
				}

			}

			//

			for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); ++tau)
			{
				for (int sii = 0; sii < assignment.g_DTAscenario_vector.size(); sii++)
				{
					int scenario_index = assignment.g_DTAscenario_vector[sii].scenario_index;
					{

						for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
						{

							fprintf(g_pFileLinkMOE, "%.1f,", g_link_vector[i].recorded_volume_per_mode_type_per_period[tau][at][scenario_index]);

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
		assignment.summary_file << "Output 2 Way Link Performance Summary" << endl;

		dtalog.output() << "writing link_performances_summary_2way.csv.." << endl;
		fopen_ss(&g_pFileLinkMOE, "link_performance_summary_2way.csv", "w");
		if (!g_pFileLinkMOE)
		{
			dtalog.output() << "File link_performance_summary2_way.csv cannot be opened." << endl;
			g_program_stop();
		}


	}
	else
	{
		assignment.summary_file << "Output subarea Link Performance:" << endl;

		dtalog.output() << "writing subarea_link_performance.csv.." << endl;
		fopen_ss(&g_pFileLinkMOE, "subarea_link_performance.csv", "w");
		if (!g_pFileLinkMOE)
		{
			dtalog.output() << "File subarea_link_performance.csv cannot be opened." << endl;
			g_program_stop();
		}
	}

	if (g_pFileLinkMOE != NULL)
	{

		int ref_volume_count = 0;
		float total_ref_volume_dev_abs_percentage = 0;


		fprintf(g_pFileLinkMOE, "link_id,link_type,from_node_id,to_node_id,geometry,distance_km,distance_mile,fftt,meso_link_id,tmc,tmc_corridor_name,tmc_corridor_id,tmc_road_order,tmc_road_sequence,subarea_id,");

		for (int sii = 0; sii < assignment.g_DTAscenario_vector.size(); sii++)
		{
			int scenario_index = assignment.g_DTAscenario_vector[sii].scenario_index;
			fprintf(g_pFileLinkMOE, "lanes_s%d,", scenario_index);
		}

		for (int sii = 0; sii < assignment.g_DTAscenario_vector.size(); sii++)
		{
			int scenario_index = assignment.g_DTAscenario_vector[sii].scenario_index;
			fprintf(g_pFileLinkMOE, "type_s%d,", scenario_index);
		}

		for (int sii = 0; sii < assignment.g_DTAscenario_vector.size(); sii++)
		{
			int scenario_index = assignment.g_DTAscenario_vector[sii].scenario_index;
			fprintf(g_pFileLinkMOE, "penaltys%d,", scenario_index);
		}

		for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); ++tau)
		{
			for (int sii = 0; sii < assignment.g_DTAscenario_vector.size(); sii++)
			{
				int scenario_index = assignment.g_DTAscenario_vector[sii].scenario_index;
				{
					int scenario_index_no = assignment.g_active_DTAscenario_map[scenario_index];

					for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
					{

						fprintf(g_pFileLinkMOE, "V%s%s%s,", assignment.g_DemandPeriodVector[tau].demand_period.c_str(), assignment.g_DTAscenario_vector[scenario_index_no].scenario_name.c_str(), assignment.g_ModeTypeVector[at].mode_type.c_str());
					}

				}
			}

		}
		for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); ++tau)
		{
			for (int sii = 0; sii < assignment.g_DTAscenario_vector.size(); sii++)
			{
				int scenario_index = assignment.g_DTAscenario_vector[sii].scenario_index;
				{
					int scenario_index_no = assignment.g_active_DTAscenario_map[scenario_index];

					for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
					{

						fprintf(g_pFileLinkMOE, "VAB%s%s%s,", assignment.g_DemandPeriodVector[tau].demand_period.c_str(), assignment.g_DTAscenario_vector[scenario_index_no].scenario_name.c_str(), assignment.g_ModeTypeVector[at].mode_type.c_str());
					}

				}
			}

		}
		for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); ++tau)
		{
			for (int sii = 0; sii < assignment.g_DTAscenario_vector.size(); sii++)
			{
				int scenario_index = assignment.g_DTAscenario_vector[sii].scenario_index;
				{
					int scenario_index_no = assignment.g_active_DTAscenario_map[scenario_index];
					for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
					{

						fprintf(g_pFileLinkMOE, "VBA%s%s%s,", assignment.g_DemandPeriodVector[tau].demand_period.c_str(), assignment.g_DTAscenario_vector[scenario_index_no].scenario_name.c_str(), assignment.g_ModeTypeVector[at].mode_type.c_str());
					}

				}
			}

		}
		fprintf(g_pFileLinkMOE, "block,");

		for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); ++tau)
		{
			for (int sii = 0; sii < assignment.g_DTAscenario_vector.size(); sii++)
			{
				int scenario_index = assignment.g_DTAscenario_vector[sii].scenario_index;

				int scenario_index_no = assignment.g_active_DTAscenario_map[scenario_index];
				for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
				{
					fprintf(g_pFileLinkMOE, "%s%s%sMEUVAB,", assignment.g_DemandPeriodVector[tau].demand_period.c_str(), assignment.g_DTAscenario_vector[scenario_index_no].scenario_name.c_str(), assignment.g_ModeTypeVector[at].mode_type.c_str());
					fprintf(g_pFileLinkMOE, "%s%s%scapAB,", assignment.g_DemandPeriodVector[tau].demand_period.c_str(), assignment.g_DTAscenario_vector[scenario_index_no].scenario_name.c_str(), assignment.g_ModeTypeVector[at].mode_type.c_str());
					fprintf(g_pFileLinkMOE, "%s%s%sDOCAB,", assignment.g_DemandPeriodVector[tau].demand_period.c_str(), assignment.g_DTAscenario_vector[scenario_index_no].scenario_name.c_str(), assignment.g_ModeTypeVector[at].mode_type.c_str());
					fprintf(g_pFileLinkMOE, "%s%s%sTTAB,", assignment.g_DemandPeriodVector[tau].demand_period.c_str(), assignment.g_DTAscenario_vector[scenario_index_no].scenario_name.c_str(), assignment.g_ModeTypeVector[at].mode_type.c_str());
				}


			}
		}

		for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); ++tau)
		{
			for (int sii = 0; sii < assignment.g_DTAscenario_vector.size(); sii++)
			{
				int scenario_index = assignment.g_DTAscenario_vector[sii].scenario_index;
				int scenario_index_no = assignment.g_active_DTAscenario_map[scenario_index];
				for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
				{
					fprintf(g_pFileLinkMOE, "%s%s%sMEUVBA,", assignment.g_DemandPeriodVector[tau].demand_period.c_str(), assignment.g_DTAscenario_vector[scenario_index_no].scenario_name.c_str(), assignment.g_ModeTypeVector[at].mode_type.c_str());
					fprintf(g_pFileLinkMOE, "%s%s%scapBA,", assignment.g_DemandPeriodVector[tau].demand_period.c_str(), assignment.g_DTAscenario_vector[scenario_index_no].scenario_name.c_str(), assignment.g_ModeTypeVector[at].mode_type.c_str());
					fprintf(g_pFileLinkMOE, "%s%s%sDOCBA,", assignment.g_DemandPeriodVector[tau].demand_period.c_str(), assignment.g_DTAscenario_vector[scenario_index_no].scenario_name.c_str(), assignment.g_ModeTypeVector[at].mode_type.c_str());
					fprintf(g_pFileLinkMOE, "%s%s%sTTBA,", assignment.g_DemandPeriodVector[tau].demand_period.c_str(), assignment.g_DTAscenario_vector[scenario_index_no].scenario_name.c_str(), assignment.g_ModeTypeVector[at].mode_type.c_str());
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
			total_vehicle_volume += g_link_vector[i].VDF_period[tau].volume_before_sa;
			total_vehicle_volume += g_link_vector[i].total_volume_for_all_mode_types_per_period[tau] + g_link_vector[i].VDF_period[tau].preload;
			scenario_code_count += g_link_vector[i].VDF_period[tau].scenario_code.size();

		}

		//if (scenario_code_count == 0 && total_vehicle_volume < 0.001 && g_link_vector.size() > 20000)
		//	continue;




		float speed = g_link_vector[i].free_speed;  // default speed 


		fprintf(g_pFileLinkMOE, "%s,%s,%d,%d,",

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


		for (int sii = 0; sii < assignment.g_DTAscenario_vector.size(); sii++)
		{
			int scenario_index = assignment.g_DTAscenario_vector[sii].scenario_index;
			fprintf(g_pFileLinkMOE, "%.1f,", g_link_vector[i].number_of_lanes_si[scenario_index]);

		}

		for (int sii = 0; sii < assignment.g_DTAscenario_vector.size(); sii++)
		{
			int scenario_index = assignment.g_DTAscenario_vector[sii].scenario_index;
			fprintf(g_pFileLinkMOE, "%s,", assignment.g_LinkTypeMap[g_link_vector[i].link_type_si[scenario_index]].link_type_name.c_str());
		}
		for (int sii = 0; sii < assignment.g_DTAscenario_vector.size(); sii++)
		{
			int scenario_index = assignment.g_DTAscenario_vector[sii].scenario_index;
			fprintf(g_pFileLinkMOE, "%.1f,", g_link_vector[i].penalty_si_flag[scenario_index]);
		}

		for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
		{
			for (int sii = 0; sii < assignment.g_DTAscenario_vector.size(); sii++)
			{
				int scenario_index = assignment.g_DTAscenario_vector[sii].scenario_index;
				{

					for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
					{

						double AB_volume = g_link_vector[i].recorded_volume_per_mode_type_per_period[tau][at][scenario_index];
						double BA_volume = 0;

						if (g_link_vector[i].BA_link_no >= 0)
						{
							BA_volume = g_link_vector[g_link_vector[i].BA_link_no].recorded_volume_per_mode_type_per_period[tau][at][scenario_index];
						}

						double total_volume = AB_volume + BA_volume;

						fprintf(g_pFileLinkMOE, "%.1f,", total_volume);

					}
				}
			}

		}

		for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
		{
			for (int sii = 0; sii < assignment.g_DTAscenario_vector.size(); sii++)
			{
				int scenario_index = assignment.g_DTAscenario_vector[sii].scenario_index;
				{

					for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
					{

						fprintf(g_pFileLinkMOE, "%.1f,", g_link_vector[i].recorded_volume_per_mode_type_per_period[tau][at][scenario_index]);

					}
				}
			}

		}

		for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
		{
			for (int sii = 0; sii < assignment.g_DTAscenario_vector.size(); sii++)
			{
				int scenario_index = assignment.g_DTAscenario_vector[sii].scenario_index;
				{

					for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
					{
						if (g_link_vector[i].BA_link_no >= 0)
							fprintf(g_pFileLinkMOE, "%.1f,", g_link_vector[g_link_vector[i].BA_link_no].recorded_volume_per_mode_type_per_period[tau][at][scenario_index]);
						else
							fprintf(g_pFileLinkMOE, ",");

					}
				}
			}

		}

		fprintf(g_pFileLinkMOE, ",");

		for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
		{
			for (int sii = 0; sii < assignment.g_DTAscenario_vector.size(); sii++)
			{
				int scenario_index = assignment.g_DTAscenario_vector[sii].scenario_index;
				{

					for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
					{

						fprintf(g_pFileLinkMOE, "%.1f,", g_link_vector[i].converted_MEU_volume_per_period_per_at[tau][at]);
						fprintf(g_pFileLinkMOE, "%.1f,", g_link_vector[i].VDF_period[tau].capacity_at[at]);
						fprintf(g_pFileLinkMOE, "%.3f,", g_link_vector[i].VDF_period[tau].DOC_mode[at]);
						fprintf(g_pFileLinkMOE, "%.3f,", g_link_vector[i].travel_time_per_period[tau][at]);

					}
				}
			}

		}
		//

		{

			for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); ++tau)
			{
				for (int sii = 0; sii < assignment.g_DTAscenario_vector.size(); sii++)
				{
					int scenario_index = assignment.g_DTAscenario_vector[sii].scenario_index;
					{

						for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
						{
							if (g_link_vector[i].BA_link_no >= 0)
							{
								fprintf(g_pFileLinkMOE, "%.1f,", g_link_vector[g_link_vector[i].BA_link_no].converted_MEU_volume_per_period_per_at[tau][at]);
								fprintf(g_pFileLinkMOE, "%.3f,", g_link_vector[g_link_vector[i].BA_link_no].VDF_period[tau].capacity_at[at]);
								fprintf(g_pFileLinkMOE, "%.3f,", g_link_vector[g_link_vector[i].BA_link_no].VDF_period[tau].DOC_mode[at]);
								fprintf(g_pFileLinkMOE, "%.3f,", g_link_vector[g_link_vector[i].BA_link_no].travel_time_per_period[tau][at]);
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
		dtalog.output() << "link based user equilibrum mode: no od_performance.csv output.." << endl;
		return;
	}

	FILE* g_pFilePathMOE = nullptr;
	fopen_ss(&g_pFilePathMOE, "od_performance.csv", "w");

	dtalog.output() << "writing od_performance.csv.." << endl;

	double path_time_vector[MAX_LINK_SIZE_IN_A_PATH];

	fopen_ss(&g_pFilePathMOE, "od_performance.csv", "w");

	if (!g_pFilePathMOE)
	{
		dtalog.output() << "File od_performance.csv cannot be opened." << endl;
		g_program_stop();
	}

	fprintf(g_pFilePathMOE, "od_no,o_zone_id,d_zone_id,o_sindex,d_sindex,o_district_id,d_district_id,mode_type,demand_period,volume,connectivity_flag,s_x_coord,s_y_coord,t_x_coord,t_y_coord,path_FF_travel_time_min,distance_km,distance_mile,");

	for (int d = 0; d < assignment.g_number_of_sensitivity_analysis_iterations; d++)
	{
		fprintf(g_pFilePathMOE, "OD%d,", d);
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
	dtalog.output() << "writing data for " << zone_size << "  zones " << endl;

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
			dtalog.output() << "o zone id =  " << g_zone_vector[orig].zone_id << endl;
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


						if (p_column_pool->path_node_sequence_map.size() == 0)
						{

							fprintf(g_pFilePathMOE, "%d,%d,%d,%d,%d,%d,%d,",
								count,
								g_zone_vector[orig].zone_id,
								g_zone_vector[dest].zone_id,
								g_zone_vector[orig].sindex,
								g_zone_vector[dest].sindex,
								assignment.g_zone_seq_no_to_analysis_distrct_id_mapping[orig],
								assignment.g_zone_seq_no_to_analysis_distrct_id_mapping[dest]

							);
							fprintf(g_pFilePathMOE, "%s,%s,%.2f,0,",
								assignment.g_ModeTypeVector[at].mode_type.c_str(),
								assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
								p_column_pool->od_volume[assignment.active_scenario_index]);

							fprintf(g_pFilePathMOE, "%f,%f,%f,%f\n",
								g_zone_vector[orig].cell_x,
								g_zone_vector[orig].cell_y,
								g_zone_vector[dest].cell_x,
								g_zone_vector[dest].cell_y);
							continue;
						}
						// scan through the map with different node sum for different continuous paths
						it_begin = p_column_pool->path_node_sequence_map.begin();
						it_end = p_column_pool->path_node_sequence_map.end();


						for (it = it_begin; it != it_end; ++it)
						{
							if (it->second.subarea_output_flag == 0)
								continue;

							if (count % 100000 == 0)
							{
								end_t = clock();
								iteration_t = end_t - start_t;
								dtalog.output() << "writing " << count / 1000 << "K agents with CPU time " << iteration_t / 1000.0 << " s" << endl;
							}

							path_toll = 0;
							path_distance = 0;
							path_travel_time = 0;
							path_travel_time_without_access_link = 0;
							path_FF_travel_time = 0;

							path_time_vector[0] = time_stamp;
							string link_code_str;

							for (int nl = 0; nl < it->second.m_link_size; ++nl)  // arc a
							{

								int link_seq_no = it->second.path_link_vector[nl];
								if (g_link_vector[link_seq_no].link_type_si[assignment.active_scenario_index] >= 0)
								{
									path_toll += g_link_vector[link_seq_no].VDF_period[tau].toll[at];
									path_distance += g_link_vector[link_seq_no].link_distance_km;
									float link_travel_time = g_link_vector[link_seq_no].travel_time_per_period[tau][at];
									path_travel_time += link_travel_time;

									if (g_link_vector[link_seq_no].link_type_si[assignment.active_scenario_index] != 1000)   // skip access link
									{
										path_travel_time_without_access_link += link_travel_time;
									}
									path_FF_travel_time += g_link_vector[link_seq_no].VDF_period[tau].FFTT_at[at];
									time_stamp += link_travel_time;
									path_time_vector[nl + 1] = time_stamp;

									link_code_str += g_link_vector[link_seq_no].link_code_str;
								}
								else
								{
									int virtual_link = 0;
								}
							}

							double total_agent_path_travel_time = 0;

							for (int vi = 0; vi < it->second.agent_simu_id_vector.size(); ++vi)
							{
								int agent_simu_id = it->second.agent_simu_id_vector[vi];
								CAgent_Simu* pAgentSimu = g_agent_simu_vector[agent_simu_id];
								total_agent_path_travel_time += pAgentSimu->path_travel_time_in_min;
							}

							double final_path_travel_time = path_travel_time;  // by default

							if (it->second.agent_simu_id_vector.size() > 1)  // with simulated agents
							{
								final_path_travel_time = total_agent_path_travel_time / it->second.agent_simu_id_vector.size();
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


							// assignment_mode = 1, path flow mode
							if (assignment.assignment_mode >= 1)
							{

								int o_zone_mode_type_cover_flag = 0;
								int d_zone_mode_type_cover_flag = 0;

								if (assignment.g_ModeTypeVector[at].zone_id_cover_map.find(g_zone_vector[orig].zone_id) != assignment.g_ModeTypeVector[at].zone_id_cover_map.end())
								{
									o_zone_mode_type_cover_flag = 1;
								}

								if (assignment.g_ModeTypeVector[at].zone_id_cover_map.find(g_zone_vector[dest].zone_id) != assignment.g_ModeTypeVector[at].zone_id_cover_map.end())
								{
									d_zone_mode_type_cover_flag = 1;
								}

								int number_of_nodes = it->second.m_link_size + 1;

								if (number_of_nodes <= 1)
								{
									int i_debug = 1;
								}
								int route_no = it->second.global_path_no;
								if (route_no == -1)
									route_no = count;

								fprintf(g_pFilePathMOE, "%d,%d,%d,%d,%d,%d,%d,",
									count,
									g_zone_vector[orig].zone_id,
									g_zone_vector[dest].zone_id,
									g_zone_vector[orig].sindex,
									g_zone_vector[dest].sindex,
									assignment.g_zone_seq_no_to_analysis_distrct_id_mapping[orig],
									assignment.g_zone_seq_no_to_analysis_distrct_id_mapping[dest]);

								fprintf(g_pFilePathMOE, "%s,%s,%.2f,1,%f,%f,%f,%f,%.4f,%.4f,%.4f,",
									assignment.g_ModeTypeVector[at].mode_type.c_str(),
									assignment.g_DemandPeriodVector[tau].demand_period.c_str(),
									it->second.path_volume,
									g_zone_vector[orig].cell_x,
									g_zone_vector[orig].cell_y,
									g_zone_vector[dest].cell_x,
									g_zone_vector[dest].cell_y,

									path_FF_travel_time,
									path_distance,
									path_distance / 1.609
								);
								//it->second.path_seq_no,
								//information_type,
								//it->second.agent_simu_id_vector.size(),
								//it->second.subarea_output_flag,
								//it->second.measurement_flag,
								//path_toll,
								//final_path_travel_time,
								//path_travel_time,
								//path_travel_time_without_access_link,
								////path_FF_travel_time,
								////final_path_travel_time- path_FF_travel_time,
								//number_of_nodes,
								//path_distance);

								double prev_value = 0;

								for (int d = 0; d < assignment.g_number_of_sensitivity_analysis_iterations; d++)
								{
									double value = 0;

									if (p_column_pool->od_volume_per_iteration_map.find(d) != p_column_pool->od_volume_per_iteration_map.end())
									{
										value = p_column_pool->od_volume_per_iteration_map[d];
									}

									fprintf(g_pFilePathMOE, "%f,", value - prev_value);
									prev_value = value;
								}
								fprintf(g_pFilePathMOE, "\n");


								count++;
							}

						}
					}
				}
			}
		}
	}
	fclose(g_pFilePathMOE);

	if (column_pool_ready_flag)
	{
		assignment.summary_file << "     Check OD connectivity and accessibility in od_performance.csv" << endl;
		assignment.summary_file << ", # of connected OD pairs =, " << l_origin_destination_map.size() << endl;
		assignment.summary_file << ", # of OD/mode_type/demand_type columns without paths =, " << l_origin_destination_disconnected_map.size() << endl;
		dtalog.output() << ", # of connected OD pairs = " << l_origin_destination_map.size() << endl;

	}
	if (l_origin_destination_map.size() == 0)
	{
		g_OutputModelFiles(10); // label cost tree
		dtalog.output() << "Please check the connectivity of OD pairs and in network and field allow_uses in link.csv." << endl;
		cout << "Please check the model_shortest_path_tree.csv file." << endl;
		//			g_program_stop();
	}


}

void g_output_dynamic_link_performance(Assignment& assignment, int output_mode = 1)
{
	//dtalog.output() << "writing dynamic_waiting_performance_profile.csv.." << endl;

	//int b_debug_detail_flag = 0;
	//FILE* g_pFileLinkMOE = nullptr;

	//string file_name = "dynamic_link_waiting_time_profile.csv";

	// fopen_ss(&g_pFileLinkMOE, file_name.c_str(), "w");

	//if (!g_pFileLinkMOE)
	//{
	//    dtalog.output() << "File " << file_name.c_str() << " cannot be opened." << endl;
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
	dtalog.output() << "writing td_link_performance.csv.." << endl;
	cout << "writing td_link_performance.csv.." << endl;

	int b_debug_detail_flag = 0;
	FILE* g_pFileLinkMOE = nullptr;

	string file_name = "td_link_performance.csv";

	fopen_ss(&g_pFileLinkMOE, file_name.c_str(), "w");

	if (!g_pFileLinkMOE)
	{
		dtalog.output() << "File " << file_name.c_str() << " cannot be opened." << endl;
		g_program_stop();
	}
	else
	{

		// Option 2: BPR-X function
		fprintf(g_pFileLinkMOE, "link_id,tmc_corridor_name,link_type_name,from_node_id,to_node_id,meso_link_id,from_cell_code,lanes,length,free_speed_kmph,free_speed_mph,FFTT,time_period,inflow_volume,outflow_volume,CA,CD,density,queue_length_in_process,queue_ratio,discharge_cap,TD_free_flow_travel_time,waiting_time_in_sec,speed_kmph,speed_mph,geometry,");
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

					fprintf(g_pFileLinkMOE, "%s,%s,%s,%d,%d,%d,%s,%d,%.3f,%.3f,%.3f,%.3f,%s_%s,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,\"%s\",",
						g_link_vector[i].link_id.c_str(),
						g_link_vector[i].tmc_corridor_name.c_str(),
						g_link_vector[i].link_type_name.c_str(),

						g_node_vector[g_link_vector[i].from_node_seq_no].node_id,
						g_node_vector[g_link_vector[i].to_node_seq_no].node_id,
						g_link_vector[i].meso_link_id,
						g_node_vector[g_link_vector[i].from_node_seq_no].cell_str.c_str(),
						g_link_vector[i].number_of_lanes_si[0],
						g_link_vector[i].link_distance_VDF,
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
	dtalog.output() << "writing log_dynamic_link_state.txt.." << endl;

	int b_debug_detail_flag = 0;
	FILE* g_pFileLinkMOE = nullptr;

	string file_name = "log_dynamic_link_state.txt";

	fopen_ss(&g_pFileLinkMOE, file_name.c_str(), "w");

	if (!g_pFileLinkMOE)
	{
		dtalog.output() << "File " << file_name.c_str() << " cannot be opened." << endl;
		g_program_stop();
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
//        dtalog.output() << "writing agent.csv.." << endl;
//
//        double path_time_vector[MAX_LINK_SIZE_IN_A_PATH];
//        FILE* g_pFileAgent = nullptr;
//        fopen_ss(&g_pFileAgent, "agent.csv", "w");
//
//        if (!g_pFileAgent)
//        {
//            dtalog.output() << "File agent.csv cannot be opened." << endl;
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
//        dtalog.output() << "writing data for " << zone_size << "  zones " << endl;
//
//        for (int orig = 0; orig < zone_size; ++orig)
//        {
//            if (g_zone_vector[orig].zone_id % 100 == 0)
//                dtalog.output() << "o zone id =  " << g_zone_vector[orig].zone_id << endl;
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
//                                    dtalog.output() << "writing " << count / 1000 << "K agents with CPU time " << iteration_t / 1000.0 << " s" << endl;
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
//                                            dtalog.output() << "error: it->second.m_link_size >= MAX_LINK_SIZE_IN_A_PATH-1" << endl;
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
		dtalog.output() << "writing agent.csv.." << endl;
		assignment.summary_file << ",summary by multi-modal agents,mode_type,#_agents,avg_distance,avg_travel_time_in_min,avg_free_speed," << endl;

		double path_time_vector[MAX_LINK_SIZE_IN_A_PATH];
		FILE* g_pFilePathMOE = nullptr;
		fopen_ss(&g_pFilePathMOE, "agent.csv", "w");

		if (!g_pFilePathMOE)
		{
			dtalog.output() << "File agent.csv cannot be opened." << endl;
			g_program_stop();
		}

		fprintf(g_pFilePathMOE, "first_column,agent_id,o_zone_id,d_zone_id,OD_key,path_id,activity_zone_sequence,activity_mode_type_sequence,display_code,impacted_flag,info_receiving_flag,diverted_flag,mode_type_no,mode_type,PCE_unit,demand_period,volume,toll,departure_time,dt_hhmm,departure_time_in_1min,departure_time_in_5_min,travel_time,distance_mile,speed_mph,waiting_time,max_link_waiting_time,max_wait_link,\n");

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

		dtalog.output() << "writing data for " << zone_size << "  zones " << endl;



		for (int at = 0; at < mode_type_size; ++at)
		{
			int agent_count = 0;
			double total_agent_distance = 0;
			double total_agent_travel_time = 0;

			for (int orig = 0; orig < zone_size; ++orig)
			{
				if (g_zone_vector[orig].zone_id % 100 == 0)
					dtalog.output() << "o zone id =  " << g_zone_vector[orig].zone_id << endl;

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
									dtalog.output() << "writing " << count / 1000 << "K agents with CPU time " << iteration_t / 1000.0 << " s" << endl;
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
									path_toll += g_link_vector[link_seq_no].VDF_period[tau].toll[at];
									path_distance += g_link_vector[link_seq_no].link_distance_VDF;
									float link_travel_time = g_link_vector[link_seq_no].travel_time_per_period[tau][at];

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
										fprintf(g_pFilePathMOE, ",%d,%d,%d,%d->%d,%d,%s,%s,%d,%d,%d,%d,%d,%s,",
											pAgentSimu->agent_id,
											g_zone_vector[orig].zone_id,
											g_zone_vector[dest].zone_id,
											g_zone_vector[orig].zone_id,
											g_zone_vector[dest].zone_id,
											it->second.path_seq_no,
											p_column_pool->activity_zone_sequence.c_str(),
											p_column_pool->activity_mode_type_sequence.c_str(),
											0,
											pAgentSimu->impacted_flag,
											pAgentSimu->info_receiving_flag,
											pAgentSimu->diverted_flag,
											pAgentSimu->mode_type_no,
											assignment.g_ModeTypeVector[at].mode_type.c_str());

										fprintf(g_pFilePathMOE, "%d,%s,1,0,%.4f,T%02d%02d,%d,%d,%.4f,%.4f,",
											pAgentSimu->PCE_unit_size,
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
				<< "," << total_agent_travel_time / max(1, agent_count) << "," << total_agent_distance / max(0.0001, total_agent_travel_time / 60.0) << endl;


		}
		fclose(g_pFilePathMOE);
	}
	assignment.summary_file << ",total_impacted_vehicles=," << total_impacted_vehicle_count << "," << endl;
	assignment.summary_file << ",total_diverted_vehicle_count=," << total_diverted_vehicle_count << "," <<
		",precentage=," << total_diverted_vehicle_count * 100.0 / max(total_impacted_vehicle_count, 1) <<
		"%, avg travel time =, " << total_diverted_vehicle_travel_time / max(1, total_diverted_vehicle_count) <<
		", avg distance =, " << total_diverted_vehicle_travel_distance / max(1, total_diverted_vehicle_count) <<
		endl;

	assignment.summary_file << ",total_diverted_DMS_veh_count=," << total_diverted_DMS_vehicle_count << "," <<
		",precentage=," << total_diverted_DMS_vehicle_count * 100.0 / max(total_impacted_vehicle_count, 1) <<
		"%, avg travel time =, " << total_diverted_DMS_vehicle_travel_time / max(1, total_diverted_DMS_vehicle_count) <<
		", avg distance =, " << total_diverted_DMS_vehicle_travel_distance / max(1, total_diverted_DMS_vehicle_count) <<
		endl;

	assignment.summary_file << ",total_non_diverted_vehicle_count=," << total_non_diverted_vehicle_count << "," <<
		",precentage=," << total_non_diverted_vehicle_count * 100.0 / max(total_impacted_vehicle_count, 1) <<
		"%, avg travel time =, " << total_non_diverted_vehicle_travel_time / max(1, total_non_diverted_vehicle_count) <<
		", avg distance =, " << total_non_diverted_vehicle_travel_distance / max(1, total_non_diverted_vehicle_count) <<
		endl;
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
		dtalog.output() << "writing trajectory.csv.." << endl;

		double path_time_vector[MAX_LINK_SIZE_IN_A_PATH];
		FILE* g_pFilePathMOE = nullptr;
		fopen_ss(&g_pFilePathMOE, "trajectory.csv", "w");

		if (!g_pFilePathMOE)
		{
			dtalog.output() << "File trajectory.csv cannot be opened." << endl;
			g_program_stop();
		}

		fprintf(g_pFilePathMOE, "first_column,agent_id,o_zone_id,d_zone_id,path_id,display_code,impacted_flag,info_receiving_flag,diverted_flag,mode_type_no,mode_type,PCE_unit,demand_period,volume,toll,departure_time,travel_time,distance_km,distance_mile,speed_kmph,speed_mph,waiting_time,max_link_waiting_time,max_wait_link,node_sequence,time_sequence,geometry\n");

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

		dtalog.output() << "writing data for " << zone_size << "  zones " << endl;

		for (int orig = 0; orig < zone_size; ++orig)
		{
			int from_zone_sindex = g_zone_vector[orig].sindex;
			if (from_zone_sindex == -1)
				continue;

			if (g_zone_vector[orig].zone_id % 100 == 0)
				dtalog.output() << "o zone id =  " << g_zone_vector[orig].zone_id << endl;

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
									dtalog.output() << "writing " << count / 1000 << "K agents with CPU time " << iteration_t / 1000.0 << " s" << endl;
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
									path_toll += g_link_vector[link_seq_no].VDF_period[tau].toll[at];
									path_distance_km += g_link_vector[link_seq_no].link_distance_km;
									path_distance_mile += g_link_vector[link_seq_no].link_distance_mile;
									float link_travel_time = g_link_vector[link_seq_no].travel_time_per_period[tau][at];

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
										fprintf(g_pFilePathMOE, ",%d,%d,%d,%d,%d,%d,%d,%d,%d,%s,",
											pAgentSimu->agent_id,
											g_zone_vector[orig].zone_id,
											g_zone_vector[dest].zone_id,
											it->second.path_seq_no,
											0,
											pAgentSimu->impacted_flag,
											pAgentSimu->info_receiving_flag,
											pAgentSimu->diverted_flag,
											pAgentSimu->mode_type_no,
											assignment.g_ModeTypeVector[at].mode_type.c_str());

										fprintf(g_pFilePathMOE, "%d,%s,1,0,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,",
											pAgentSimu->PCE_unit_size,
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
		assignment.summary_file << ", # of simulated agents in trajectory.csv=," << count << ",avg # of nodes per agent=" << total_number_of_nodes * 1.0 / max(1, count) << endl;
	}

	int b_trace_file = false;

	// output trace file
	if (assignment.assignment_mode == simulation_dta)  //LUE
	{
		FILE* g_pFilePathMOE = nullptr;
		fopen_ss(&g_pFilePathMOE, "trace.csv", "w");
		dtalog.output() << "writing trace.csv.." << endl;

		double path_time_vector[MAX_LINK_SIZE_IN_A_PATH];

		fopen_ss(&g_pFilePathMOE, "trace.csv", "w");

		if (!g_pFilePathMOE)
		{
			dtalog.output() << "File trace.csv cannot be opened." << endl;
			g_program_stop();
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

		dtalog.output() << "writing data for " << zone_size << "  zones " << endl;

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
									dtalog.output() << "writing " << count / 1000 << "K agents with CPU time " << iteration_t / 1000.0 << " s" << endl;
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
									path_toll += g_link_vector[link_seq_no].VDF_period[tau].toll[at];
									path_distance += g_link_vector[link_seq_no].link_distance_VDF;
									float link_travel_time = g_link_vector[link_seq_no].travel_time_per_period[tau][at];
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
		dtalog.output() << "writing trajectory.bin.." << endl;

		int path_node_vector[MAX_LINK_SIZE_IN_A_PATH];
		double path_time_vector[MAX_LINK_SIZE_IN_A_PATH];
		FILE* g_pFilePathMOE = nullptr;
		fopen_ss(&g_pFilePathMOE, "trajectory.bin", "wb");

		if (!g_pFilePathMOE)
		{
			dtalog.output() << "File trajectory.bin cannot be opened." << endl;
			g_program_stop();
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

		dtalog.output() << "writing data for " << zone_size << "  zones " << endl;

		for (int orig = 0; orig < zone_size; ++orig)
		{
			int from_zone_sindex = g_zone_vector[orig].sindex;
			if (from_zone_sindex == -1)
				continue;

			if (g_zone_vector[orig].zone_id % 100 == 0)
				dtalog.output() << "o zone id =  " << g_zone_vector[orig].zone_id << endl;

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
									dtalog.output() << "writing " << count / 1000 << "K binary agents with CPU time " << iteration_t / 1000.0 << " s" << endl;
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
									path_toll += g_link_vector[link_seq_no].VDF_period[tau].toll[at];
									path_distance_km += g_link_vector[link_seq_no].link_distance_km;
									path_distance_mile += g_link_vector[link_seq_no].link_distance_mile;
									float link_travel_time = g_link_vector[link_seq_no].travel_time_per_period[tau][at];

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
		dtalog.output() << "Comlete writing " << count / 1000 << "K binary agents with CPU time " << iteration_t / 1000.0 << " s." << endl;
		fclose(g_pFilePathMOE);
	}

}

void g_output_demand_bin(Assignment& assignment)
{

	dtalog.output() << "writing demand.bin.." << endl;

	FILE* g_pFilePathMOE = nullptr;
	fopen_ss(&g_pFilePathMOE, "output_demand.bin", "wb");

	if (!g_pFilePathMOE)
	{
		dtalog.output() << "File demand.bin cannot be opened." << endl;
		g_program_stop();
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

	dtalog.output() << "writing data for " << zone_size << "  zones " << endl;

	for (int orig = 0; orig < zone_size; ++orig)
	{
		int from_zone_sindex = g_zone_vector[orig].sindex;
		if (from_zone_sindex == -1)
			continue;

		if (g_zone_vector[orig].zone_id % 100 == 0)
			dtalog.output() << "o zone id =  " << g_zone_vector[orig].zone_id << endl;

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
	dtalog.output() << "Complete writing " << count / 1000 << "K binary demand pairs with CPU time " << iteration_t / 1000.0 << " s." << endl;
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
			dtalog.output() << "Error: File model_node.csv cannot be opened.\n It might be currently used and locked by EXCEL." << endl;
			g_program_stop();


		}

	}

	if (mode == 2)
	{
		FILE* g_pFileModelLink = fopen("model_link.csv", "w");

		if (g_pFileModelLink != NULL)
		{
			fprintf(g_pFileModelLink, "link_id,link_no,layer_no,from_node_id,to_node_id,from_gate_flag,to_gate_flag,link_type,link_type_name,lanes,link_distance_VDF,free_speed,cutoff_speed,fftt,capacity,allow_uses,");
			fprintf(g_pFileModelLink, "BPR_plf,BPR_alpha,BPR_beta,QVDF_plf,QVDF_alpha,QVDF_beta,QVDF_cd,QVDF_n,");
			fprintf(g_pFileModelLink, "geometry\n");

			//VDF_fftt1,VDF_cap1,VDF_alpha1,VDF_beta1
			for (int i = 0; i < g_link_vector.size(); i++)
			{
				if (g_link_vector[i].link_type_si[0] <= -100)  // invisible link
				{
					continue;
				}

				fprintf(g_pFileModelLink, "%s,%d,%d,%d,%d,%d,%d,%d,%s,%f,%f,%f,%f,%f,%f,%s,%f,%f,%f,%f,%f,%f,%f,%f,",
					g_link_vector[i].link_id.c_str(),
					g_link_vector[i].link_seq_no,
					g_link_vector[i].layer_no,
					g_node_vector[g_link_vector[i].from_node_seq_no].node_id,
					g_node_vector[g_link_vector[i].to_node_seq_no].node_id,
					g_node_vector[g_link_vector[i].from_node_seq_no].MRM_gate_flag,
					g_node_vector[g_link_vector[i].to_node_seq_no].MRM_gate_flag,
					g_link_vector[i].link_type_si[0],
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
			dtalog.output() << "Error: File model_link.csv cannot be opened.\n It might be currently used and locked by EXCEL." << endl;
			g_program_stop();

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
				dtalog.output() << "Error: File access_link.csv cannot be opened.\n It might be currently used and locked by EXCEL." << endl;
				g_program_stop();

			}

		}

	}

	if (mode == 3)  // cell
	{
		//FILE* g_pFileZone = nullptr;
		//g_pFileZone = fopen("model_cell.csv", "w");

		//if (g_pFileZone == NULL)
		//{
		//    cout << "File model_cell.csv cannot be opened." << endl;
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

						if (connected_flag == 1)
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
			dtalog.output() << "Error: File model_label_cost_tree.csv cannot be opened.\n It might be currently used and locked by EXCEL." << endl;
			g_program_stop();


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
