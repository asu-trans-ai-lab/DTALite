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
#include "DTA.h"
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

void g_load_dynamic_traffic_management_file(Assignment& assignment)
{
	dtalog.output() << "[PROCESS INFO] Step 4.1: Reading dynamic_traffic_management data..." << '\n';


	// we setup the initial number of lanes per demand period, per mode and for each scenario, then we will load the dynamic_traffic_management file to overwrite the # of lanes for dynamic lane use in some special cases
	for (int i = 0; i < g_link_vector.size(); ++i)
	{
		g_link_vector[i].setup_dynamic_number_of_lanes(assignment.active_scenario_index);
	}

	CDTACSVParser parser;

	int capacity_count = 0;
	int dtm_dynamic_lane_use_count = 0;
	int dtm_lane_closure_count = 0;

	int dtm_dms_count = 0;

	int record_no = 0;

	FILE* g_pFileModel_LC = fopen("log_scenario.txt", "w");
	if (g_pFileModel_LC == NULL)
		return;

	fprintf(g_pFileModel_LC, "dtm_type,demand_period,from_node,to_node,geometry,values\n");

	assignment.active_dms_count = 0;
	assignment.active_lane_closure_count = 0;

	bool bFoundFlag = false;

	if (parser.OpenCSVFile("dynamic_traffic_management.csv", false))
	{
		fprintf(g_pFileModel_LC, "reading record %d\n", record_no+1);
		record_no = record_no+1;

		while (parser.ReadRecord())
		{

			int acitcivate_flag = 0;

			parser.GetValueByFieldName("activate", acitcivate_flag);


			if (acitcivate_flag == 0)  // continue to skip this record
				continue;

			int DTM_id = 0;
			if (!parser.GetValueByFieldName("dtm_id", DTM_id))
			{
				dtalog.output() << "[ERROR] Field dtm_id in file dynamic_traffic_management.csv does not have a valid value." << '\n';
				fprintf(g_pFileModel_LC, "[ERROR] field dtm_id %d in file dynamic_traffic_management.csv does not have a valid value\n");
				g_program_stop();
			}
			int from_node_id = 0;
			if (!parser.GetValueByFieldName("from_node_id", from_node_id))
			{
				dtalog.output() << "[ERROR] from_node_id in file dynamic_traffic_management.csv is not defined." << '\n';
				fprintf(g_pFileModel_LC, "[ERROR] from_node_id %d in file dynamic_traffic_management.csv is not defined\n", from_node_id);
				g_program_stop();
			}

			int to_node_id = 0;
			if (!parser.GetValueByFieldName("to_node_id", to_node_id))
			{
				dtalog.output() << "[ERROR] to_node_id in file dynamic_traffic_management.csv is not defined." << '\n';
				fprintf(g_pFileModel_LC, "[ERROR] to_node_id %d in file dynamic_traffic_management.csv is not defined\n", to_node_id);
				g_program_stop();
			}


			if (from_node_id == 0 && to_node_id == 0)
			{
				continue;
			}

			if (assignment.g_node_id_to_seq_no_map.find(from_node_id) == assignment.g_node_id_to_seq_no_map.end())
			{
				dtalog.output() << "[ERROR] from_node_id " << from_node_id << " in file scenario.csv is not defined in node.csv." << '\n';
				fprintf(g_pFileModel_LC, "[ERROR] from_node_id %d in file scenario.csv is not defined in node.csv\n", from_node_id);
				continue;
			}
			if (assignment.g_node_id_to_seq_no_map.find(to_node_id) == assignment.g_node_id_to_seq_no_map.end())
			{
				dtalog.output() << "[ERROR] to_node_id " << to_node_id << " in file scenario.csv is not defined in node.csv." << '\n';
				fprintf(g_pFileModel_LC, "[ERROR] to_node_id %d in file scenario.csv is not defined in node.csv\n", to_node_id);
				//has not been defined
				continue;
			}

			// create a link object


			// map external node number to internal node seq no.
			int internal_from_node_seq_no = assignment.g_node_id_to_seq_no_map[from_node_id];
			int internal_to_node_seq_no = assignment.g_node_id_to_seq_no_map[to_node_id];

			int link_seq_no = 0;
			if (g_node_vector[internal_from_node_seq_no].m_to_node_2_link_seq_no_map.find(internal_to_node_seq_no) != g_node_vector[internal_from_node_seq_no].m_to_node_2_link_seq_no_map.end())
			{
				link_seq_no = g_node_vector[internal_from_node_seq_no].m_to_node_2_link_seq_no_map[internal_to_node_seq_no];
			}
			else
			{
				dtalog.output() << "[ERROR] Link " << from_node_id << "->" << to_node_id << " in file dynamic_traffic_management.csv is not defined in link.csv." << '\n';
				fprintf(g_pFileModel_LC, "[ERROR] link %d-> %d in file scenario.csv is not defined in link.csv\n", from_node_id, to_node_id);
				continue;
			}

			vector<float> global_minute_vector;
			string demand_period;
			parser.GetValueByFieldName("demand_period", demand_period);
			int tau = 0;
			if (assignment.demand_period_to_seqno_mapping.find(demand_period) != assignment.demand_period_to_seqno_mapping.end())
			{
				tau = assignment.demand_period_to_seqno_mapping[demand_period];
			}

			fprintf(g_pFileModel_LC, "%s,%d,%d,", demand_period.c_str(), from_node_id, to_node_id);
			fprintf(g_pFileModel_LC, "\"LINESTRING (");
			fprintf(g_pFileModel_LC, "%f %f,", g_node_vector[internal_from_node_seq_no].x, g_node_vector[internal_from_node_seq_no].y);
			fprintf(g_pFileModel_LC, "%f %f,", g_node_vector[internal_to_node_seq_no].x, g_node_vector[internal_to_node_seq_no].y);
			fprintf(g_pFileModel_LC, ")\"");

			string str_mode_type;
			parser.GetValueByFieldName("mode_type", str_mode_type); 


			int mode_type_no = 0;
			if (str_mode_type.size() > 0 && assignment.mode_type_2_seqno_mapping.find(str_mode_type) != assignment.mode_type_2_seqno_mapping.end())
			{
				mode_type_no = assignment.mode_type_2_seqno_mapping[str_mode_type];

			}
			else
			{
				dtalog.output() << "[ERROR] Please ensure that the mode type '" << str_mode_type << "' specified in dynamic_traffic_management.csv for record number " << DTM_id << " is present in the mode_type.csv file." << '\n';
			}


			string DTM_type, scenario_index_vector_str;
			parser.GetValueByFieldName("dtm_type", DTM_type);
			
			parser.GetValueByFieldName("scenario_index_vector", scenario_index_vector_str);
			std::vector<int> scenario_index_vector;

			g_ParserIntSequence(scenario_index_vector_str, scenario_index_vector);

			g_link_vector[link_seq_no].link_specifical_flag_str.clear();
			g_link_vector[link_seq_no].VDF_period[tau].dynamic_traffic_management_flag = 0;
			g_link_vector[link_seq_no].VDF_period[tau].lane_closure_final_lanes = 0;  // apply the change
			g_link_vector[link_seq_no].VDF_period[tau].dtm_scenario_code.clear();





			for (int scenario_index_i = 0; scenario_index_i < scenario_index_vector.size(); scenario_index_i++)
			{ 
				if (scenario_index_vector[scenario_index_i] == assignment.active_scenario_index)
				{
					assignment.g_DTA_scenario_vector[assignment.active_scenario_index].dtm_flag = 1;
					bFoundFlag = true;
				}
				
			}

			if (bFoundFlag == false)  // this active scenario index 
				continue; 
			if (DTM_type == "dynamic_lane_use")
			{
				// capacity in the space time arcs
				float final_lanes = 0;
				parser.GetValueByFieldName("final_lanes", final_lanes);

				g_link_vector[link_seq_no].recorded_lanes_per_period_per_at[tau][mode_type_no][assignment.active_scenario_index] = final_lanes;  // apply the change
				fprintf(g_pFileModel_LC, "dynamic lane use,final_number_of_lanes=%f,", final_lanes);
				dtm_dynamic_lane_use_count++;
			}
			else if (DTM_type == "lane_closure")
			{
				// capacity in the space time arcs
				float final_lanes = 0;
				parser.GetValueByFieldName("final_lanes", final_lanes);

				g_link_vector[link_seq_no].VDF_period[tau].lane_closure_final_lanes = final_lanes;  // apply the change
				g_link_vector[link_seq_no].VDF_period[tau].dynamic_traffic_management_flag = -1;
				g_link_vector[link_seq_no].link_specifical_flag_str = "lane_closure";
				//
				g_link_vector[link_seq_no].VDF_period[tau].dtm_scenario_code = "lane_closure";
				dtm_lane_closure_count++;

				if (assignment.g_ModeTypeVector[mode_type_no].real_time_information_type == 0) {
					dtalog.output() << "[ERROR] Lane closure record " << DTM_id << " in the file dynamic_traffic_management.csv is being checked. Please ensure that the mode type = " << str_mode_type << " specified in mode_type.csv has DTM_real_time_info_type = 1 to enable real-time information modeling." << '\n';
					g_program_stop();
				}



				fprintf(g_pFileModel_LC, "lane closure,final_number_of_lanes=%f,", final_lanes);
			}
			else if (DTM_type == "dms")
			{
				g_link_vector[link_seq_no].VDF_period[tau].dynamic_traffic_management_flag = 2;
				g_link_vector[link_seq_no].VDF_period[tau].dtm_scenario_code = "dms";
				g_link_vector[link_seq_no].link_specifical_flag_str = "dms";

				if (assignment.node_seq_no_2_zone_id_mapping.find(g_link_vector[link_seq_no].to_node_seq_no) == assignment.node_seq_no_2_zone_id_mapping.end())
				{
					dtalog.output() << "[ERROR]  The upstream node ('to_node_id' = " << to_node_id << ") of a DMS link must be associated with a 'zone_id' in 'node.csv' for the simulation of DMS information provision strategies." << '\n';
					g_program_stop();

				}

				int zone_id = assignment.node_seq_no_2_zone_id_mapping[g_link_vector[link_seq_no].to_node_seq_no];

				if (assignment.g_zoneid_to_zone_seq_no_mapping.find(zone_id) == assignment.g_zoneid_to_zone_seq_no_mapping.end())  // not found
					continue;

				int zone_no = assignment.g_zoneid_to_zone_seq_no_mapping[zone_id];



				assignment.zone_seq_no_2_info_mapping[zone_no] = 1;  // set information zone flag


				assignment.active_dms_count +=1;
				dtm_dms_count++;
			}
			else
			{
				dtalog.output() << "[ERROR] DTALite does not support DTM_type = " << DTM_type.c_str() << " in dynamic_traffic_management.csv " << '\n';
				dtalog.output() << "[INFO] DTALite only support DTM_type = lane_closure or dms in dynamic_traffic_management.csv " << '\n';
				g_program_stop();
			}



			//if (DTM_type == "incident")
			//{

			//	string time_period;
			//	if (!parser.GetValueByFieldName("time_period", time_period))
			//	{
			//		dtalog.output() << "Error: Field time_window in file scenario.csv cannot be read." << '\n';
			//		g_program_stop();
			//		break;
			//	}



			//	//input_string includes the start and end time of a time period with

			//	global_minute_vector = g_time_parser(time_period); //global_minute_vector incldue the starting and ending time
			//	if (global_minute_vector.size() == 2)
			//	{
			//		if (global_minute_vector[0] < assignment.g_LoadingStartTimeInMin)
			//			global_minute_vector[0] = assignment.g_LoadingStartTimeInMin;

			//		if (global_minute_vector[0] > assignment.g_LoadingEndTimeInMin)
			//			global_minute_vector[0] = assignment.g_LoadingEndTimeInMin;

			//		if (global_minute_vector[1] < assignment.g_LoadingStartTimeInMin)
			//			global_minute_vector[1] = assignment.g_LoadingStartTimeInMin;

			//		if (global_minute_vector[1] > assignment.g_LoadingEndTimeInMin)
			//			global_minute_vector[1] = assignment.g_LoadingEndTimeInMin;

			//		if (global_minute_vector[1] < global_minute_vector[0])
			//			global_minute_vector[1] = global_minute_vector[0];

			//		g_link_vector[link_seq_no].dynamic_link_event_start_time_in_min = min((float)g_link_vector[link_seq_no].dynamic_link_event_start_time_in_min, global_minute_vector[0]);
			//	}
			//	else
			//	{
			//		continue;
			//	}
			//}
			//if (DTM_type == "incident")
			//{
			//	// capacity in the space time arcs
			//	float capacity = 1;
			//	parser.GetValueByFieldName("capacity", capacity);

			//	unsigned int RandomSeed = 101;
			//	float residual;
			//	float random_ratio = 0;

			//	double per_sec_capacity = capacity / 3600;
			//	for (int s = global_minute_vector[0] * 60; s <= global_minute_vector[1] * 60; s++)
			//	{
			//		int t_simu_second = (s - assignment.g_LoadingStartTimeInMin * 60);

			//		if (t_simu_second < 0)
			//			t_simu_second = 0;


			//		if (capacity < 1)
			//		{
			//			g_link_vector[link_seq_no].m_link_pedefined_capacity_map_in_sec[t_simu_second] = 0;
			//		}
			//		else
			//		{
			//			residual = per_sec_capacity - (int)(per_sec_capacity);
			//			RandomSeed = (LCG_a * RandomSeed + LCG_c) % LCG_M;
			//			random_ratio = float(RandomSeed) / LCG_M;

			//			if (random_ratio < residual)
			//				g_link_vector[link_seq_no].m_link_pedefined_capacity_map_in_sec[t_simu_second] = 1;
			//			else
			//				g_link_vector[link_seq_no].m_link_pedefined_capacity_map_in_sec[t_simu_second] = 0;
			//		}


			////
			//	if (capacity < g_link_vector[link_seq_no].lane_capacity * g_link_vector[link_seq_no].number_of_lanes * 0.8)
			//	{
			//		g_link_vector[link_seq_no].capacity_reduction_map[tau] = 1;
			//		g_link_vector[link_seq_no].global_minute_capacity_reduction_start = global_minute_vector[0];
			//		g_link_vector[link_seq_no].global_minute_capacity_reduction_end = global_minute_vector[1];
			//	}
			//	fprintf(g_pFileModel_LC, "incident,cap=%f,", capacity);
			//	incident_count++;
			//
			//else if (DTM_type == "dms")
			//{
			//	// capacity in the space time arcs
			//	float response_rate = 1;
			//	parser.GetValueByFieldName("response_rate", response_rate);

			//	for (int t = global_minute_vector[0]; t <= global_minute_vector[1]; t++)
			//	{
			//		g_link_vector[link_seq_no].m_link_pedefined_information_response_map[t] = response_rate;
			//	}

			//	g_link_vector[link_seq_no].vms_map[tau] = 1;
			//	fprintf(g_pFileModel_LC, "dms,response_rate=%f,", response_rate);

			//	dtm_dms_count++;


			//}
			//else
			//{
			//	if (DTM_type.size() >= 1)
			//	{
			//		dtalog.output() << "Error: DTM_type = " << DTM_type << " is not supported. Currently DTALite supports DTM_type such as sa, incident, dms" << '\n';
			//		g_program_stop();
			//	}
			//}
			fprintf(g_pFileModel_LC, "\n");

		}


	}
	// allocate

	if(dtm_dynamic_lane_use_count = 0)
		dtalog.output() << "No dynamic lane use scenarios found in dynamic_traffic_management.csv." << '\n';
	else 
		dtalog.output() << "[STATUS INFO] loading " << dtm_dynamic_lane_use_count << " dynamic lane use scenarios in dynamic_traffic_management.csv  " << '\n';

	if (dtm_lane_closure_count = 0)
		dtalog.output() << "No lane closure scenarios found in dynamic_traffic_management.csv." << '\n';
	else
		dtalog.output() << "[STATUS INFO] loading " << dtm_dynamic_lane_use_count << " dynamic lane use scenarios in dynamic_traffic_management.csv  " << '\n';



	if (dtm_lane_closure_count  > 0  && assignment.g_number_of_real_time_mode_types ==0)
	{
		dtalog.output() << "[ERROR] No mode type with 'DTM_real_time_info_type' = 1 is defined in 'mode_type.csv'. If a lane closure scenario is activated, a corresponding real-time information user class (i.e., mode type) must be defined in mode_type.csv." << '\n';
		g_program_stop();
	}


	if (bFoundFlag && dtm_lane_closure_count == 0 && assignment.total_real_time_demand_volume > 0.1) {
		dtalog.output() << "[ERROR] There is positive demand volume specified for real-time information users, for scenario index = " << assignment.active_scenario_index <<  " but no lane closure scenarios are activated in the dynamic_traffic_management.csv file." << '\n';

		for (int i = 0; i < assignment.g_ModeTypeVector.size(); i++) {
			if (assignment.g_ModeTypeVector[i].real_time_information_type != 0) {
				dtalog.output() << "[INFO] The real-time information mode '" << assignment.g_ModeTypeVector[i].mode_type.c_str() << "' is defined in the mode_type.csv file." << '\n';
			}
		}

		g_program_stop();
	}


	//dtalog.output() << "[STATUS INFO] loading " << dtm_dms_count << " dms scenarios dynamic_traffic_management.csv " << '\n';
	if (dtm_dms_count > 0 && assignment.g_number_of_DMS_mode_types == 0)
	{
		dtalog.output() << "[ERROR] No mode type with 'DTM_real_time_info_type' = 2  is defined in 'mode_type.csv'. If a DMS scenario is activated, a corresponding real-time information user class (i.e., mode type) must be defined in mode_type.csv.." << '\n';
		g_program_stop();
	}
	//if (dtm_dms_count > 0 && dtm_lane_closure_count ==0 )
	//{
	//	dtalog.output() << "[ERROR] If a DMS scenario is activated, a corresponding lane closure must be defined in dynamic_traffic_management.csv.." << '\n';
	//	g_program_stop();
	//}


	//if (incident_count >= 1)
	//{
	//	g_assign_RT_computing_tasks_to_memory_blocks(assignment);
	//}
	fclose(g_pFileModel_LC);
	parser.CloseCSVFile();

	// we now know the number of links

	assignment.summary_file << ",read dynamic traffic managementscenario" << '\n';
	//assignment.summary_file << ", # of records in dynamic_traffic_management.csv=," << dtm_lane_closure_count + incident_count + dtm_dms_count << "," << '\n';
	assignment.summary_file << ", # of lane closure records in dynamic_traffic_management.csv=," << dtm_lane_closure_count << "," << '\n';
	//assignment.summary_file << ", # of incident records in dynamic_traffic_management.csv=," << incident_count << "," << '\n';
	assignment.summary_file << ", # of dms records in dynamic_traffic_management.csv=," << dtm_dms_count << "," << '\n';

	//if (capacity_count > 0 || sa_capacity_count > 0)
	//{
	//	assignment.g_number_of_sensitivity_analysis_iterations_for_dtm = max(0, assignment.g_number_of_sensitivity_analysis_iterations_for_dtm);
	//}

}


void g_load_demand_side_scenario_file(Assignment& assignment)
{
	dtalog.output() << "[PROCESS INFO] Step 2.0: Reading demand side scenario data..." << '\n';

	CDTACSVParser parser;

	int demand_scenario_count = 0;


	if (parser.OpenCSVFile("demand_side_scenario.csv", false))
	{
		while (parser.ReadRecord())
		{

			int o_district_id = 0;
			parser.GetValueByFieldName("o_district_id", o_district_id);

			int d_district_id = 0;
			parser.GetValueByFieldName("d_district_id", d_district_id);

			if (o_district_id == 0 && d_district_id == 0)
			{
				continue;
			}

			string demand_period;

			parser.GetValueByFieldName("demand_period", demand_period);
			int tau = 0;
			if (assignment.demand_period_to_seqno_mapping.find(demand_period) != assignment.demand_period_to_seqno_mapping.end())
			{
				tau = assignment.demand_period_to_seqno_mapping[demand_period];
			}
			else
			{
				dtalog.output() <<"[ERROR] demand_period" << "is not defined in settings.csv." << '\n';
				g_program_stop();
			}

			double scale_factor = 1;

			parser.GetValueByFieldName("scale_factor", scale_factor, false);

			string DTM_type;
			parser.GetValueByFieldName("DTM_type", DTM_type);
			// o_distrct_id_factor_map,
			// d_district_id_factor_map,
			// od_district_id_factor_map

			if (DTM_type == "o_based")
			{
				assignment.o_district_id_factor_map[o_district_id] = scale_factor;

			}else if (DTM_type == "d_based")
			{
				assignment.d_district_id_factor_map[d_district_id] = scale_factor;
			}else if (DTM_type == "od_based")
			{
				int od_district_key = o_district_id * 1000 + d_district_id;

				assignment.od_district_id_factor_map[od_district_key] = scale_factor;
			}else if (DTM_type == "sa_o_based")
			{
				assignment.SA_o_district_id_factor_map[o_district_id] = scale_factor;

			}
			else if (DTM_type == "sa_d_based")
			{
				assignment.SA_d_district_id_factor_map[d_district_id] = scale_factor;
			}
			else if (DTM_type == "sa_od_based")
			{
				int od_district_key = o_district_id * 1000 + d_district_id;

				assignment.SA_od_district_id_factor_map[od_district_key] = scale_factor;
			}
			else
			{
				dtalog.output() << "[ERROR] scenario type in demand_side_scenario.csv supports only o_based, d_based, od_based, SA_o_based, SA_d_based, SA_od_based." << '\n';
				dtalog.output() << DTM_type.c_str() << "is not supported." << '\n';
				g_program_stop();
			}

			demand_scenario_count++;

			dtalog.output() << "[STATUS INFO] reading " << demand_scenario_count << " demand side scenario(s).. " << '\n';

		}
	}
	parser.CloseCSVFile();

	// we now know the number of links
	assignment.summary_file << ", # of demand side records in demand_side_scenario.csv=," << demand_scenario_count << "," << '\n';

}


//char lr_price_field_name[50];

//for (int at = 0; at < assignment.g_ModeTypeVector.size(); at++)
//{

//	double LR_ini_price = 0;
//	double LR_RT_price = 0;

//	sprintf(lr_price_field_name, "lr_price_%s", assignment.g_ModeTypeVector[at].mode_type.c_str());
//	parser.GetValueByFieldName(lr_price_field_name, LR_ini_price, true, false);
//	g_link_vector[link_seq_no].VDF_period[tau].LR_price[at] = LR_ini_price;

//	//if (capacity < 1)
//	g_link_vector[link_seq_no].VDF_period[tau].RT_allowed_use[at] = false;

//	if (assignment.g_ModeTypeVector[at].real_time_information_type >= 1)
//	{
//		sprintf(lr_price_field_name, "lr_rt_price_%s", assignment.g_ModeTypeVector[at].mode_type.c_str());
//		parser.GetValueByFieldName(lr_price_field_name, LR_RT_price, true, false);
//		g_link_vector[link_seq_no].VDF_period[tau].LR_RT_price[at] = LR_RT_price;

//		if (fabs(LR_RT_price) > 0.001)
//		{
//			dtalog.output() << "link " << from_node_id << "->" << to_node_id << " has a lr RT price of " << g_link_vector[link_seq_no].VDF_period[tau].LR_RT_price[at] << " for agent type "
//				<< assignment.g_ModeTypeVector[at].mode_type.c_str() << " at demand period " << demand_period.c_str() << '\n';
//		}

//	}

//	if (fabs(LR_ini_price) > 0.001)
//	{
//		dtalog.output() << "link " << from_node_id << "->" << to_node_id << " has a lr price of " << g_link_vector[link_seq_no].VDF_period[tau].LR_price[at] << " for agent type "
//			<< assignment.g_ModeTypeVector[at].mode_type.c_str() << " at demand period " << demand_period.c_str() << '\n';
//	}

//}

