/* Portions Copyright 2019-2021 Xuesong Zhou and Peiheng Li, Cafer Avci

 * If you help write or modify the code, please also list your names here.
 * The reason of having Copyright info here is to ensure all the modified version, as a whole, under the GPL
 * and further prevent a violation of the GPL.
 *
 * More about "How to use GNU licenses for your own software"
 * http://www.gnu.org/licenses/gpl-howto.html
 */


 // Peiheng, 02/03/21, remove them later after adopting better casting
#pragma warning(disable : 4305 4267 4018 )
// stop warning: "conversion from 'int' to 'float', possible loss of data"
#pragma warning(disable: 4244)
#pragma warning(disable: 4267)
#pragma warning(disable: 4477)


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
#include <iomanip> // Include the <iomanip> library for formatting output

using std::max;
using std::min;
using std::string;
using std::vector;
using std::map;
using std::ifstream;
using std::ofstream;
using std::istringstream;


bool isGlobalPolygon(std::string polygon) {
	std::string globalPolygon = "POLYGON ((-180 -90, 180 -90, 180 90, -180 90, -180 -90))";

	// Trim leading and trailing whitespace from the input string
	polygon.erase(polygon.begin(), std::find_if(polygon.begin(), polygon.end(), [](unsigned char ch) { return !std::isspace(ch); }));
	polygon.erase(std::find_if(polygon.rbegin(), polygon.rend(), [](unsigned char ch) { return !std::isspace(ch); }).base(), polygon.end());

	return polygon == globalPolygon;
}

bool g_read_subarea_CSV_file(Assignment& assignment)
{
	CDTACSVParser parser;
	if (parser.OpenCSVFile("subarea.csv", false))
	{
		while (parser.ReadRecord())
		{
			string geo_string;
			if (parser.GetValueByFieldName("geometry", geo_string))
			{


				if (isGlobalPolygon(geo_string) == true && g_zone_vector.size() > 1000)
				{
					dtalog.output() << "[PROCESS INFO] The file 'subarea.csv' currently contains the default global coordinate system, bypassing the need for subarea processing, for a system with more than 1,000 zones." << '\n';
					g_DTA_log_file << "[PROCESS INFO] The file 'subarea.csv' currently contains the default global coordinate system, bypassing the need for subarea processing, for a system with more than 1,000 zones." << '\n';
					break;
				}
					



				// overwrite when the field "geometry" exists
				CDTAGeometry geometry(geo_string);

				std::vector<CCoordinate> CoordinateVector = geometry.GetCoordinateList();

				if (CoordinateVector.size() > 0)
				{
					for (int p = 0; p < CoordinateVector.size(); p++)
					{
						DTAGDPoint pt;
						pt.x = CoordinateVector[p].X;
						pt.y = CoordinateVector[p].Y;

						assignment.g_subarea_shape_points.push_back(pt);
					}
				}

				dtalog.output() << "[PROCESS INFO] Validated subarea.csv. It contains " << assignment.g_subarea_shape_points.size() << " geometric points. Proceeding to FOCUSING analysis step." << '\n';
				g_DTA_log_file << "[PROCESS INFO] Validated subarea.csv. It contains " << assignment.g_subarea_shape_points.size() << " geometric points. Proceeding to FOCUSING analysis step." << '\n';

				break;
			}

		}
		parser.CloseCSVFile();
	}


	return true;
}


void g_read_departure_time_profile(Assignment& assignment)
{
	CDTACSVParser parser;
	
	dtalog.output() << "[PROCESS INFO] Step 2.0: Reading file departure_time_profile.csv" << '\n';
	g_DTA_log_file << "[PROCESS INFO] Step 2.0: Reading file departure_time_profile.csv" << '\n';


	if (parser.OpenCSVFile("departure_time_profile.csv", false))
	{
		while (parser.ReadRecord())
		{
			int departure_time_profile_no = 0;
			string time_period;

			CDeparture_time_Profile dep_time;

			if (!parser.GetValueByFieldName("departure_time_profile_no", departure_time_profile_no))
				break;

			dep_time.departure_time_profile_no = departure_time_profile_no;
			if (departure_time_profile_no != assignment.g_DepartureTimeProfileVector.size())
			{
				dtalog.output() << "[ERROR] Field departure_time_profile_no in field departure_time_profile should be sequential as a value of ." << assignment.g_DepartureTimeProfileVector.size() << '\n';
				g_DTA_log_file << "[ERROR] Field departure_time_profile_no in field departure_time_profile should be sequential as a value of ." << assignment.g_DepartureTimeProfileVector.size() << '\n';
				departure_time_profile_no = assignment.g_DepartureTimeProfileVector.size();

			}

			vector<float> global_minute_vector;

			if (!parser.GetValueByFieldName("time_period", time_period))
			{
				dtalog.output() << "[ERROR] Field time_period in field departure_time_profile cannot be read." << '\n';
				g_DTA_log_file << "[ERROR] Field time_period in field departure_time_profile cannot be read." << '\n';
				time_period = "0700_0800"; 
			}


			//input_string includes the start and end time of a time period with hhmm format
			global_minute_vector = g_time_parser(time_period); //global_minute_vector incldue the starting and ending time

			if (global_minute_vector.size() == 2)
			{
				dep_time.starting_time_slot_no = global_minute_vector[0] / MIN_PER_TIMESLOT;  // read the data
				dep_time.ending_time_slot_no = global_minute_vector[1] / MIN_PER_TIMESLOT;    // read the data from settings.csv

			}

			char time_interval_field_name[20];
			char time_interval_field_name2[20];


			for (int s = dep_time.starting_time_slot_no; s <= dep_time.ending_time_slot_no; s += 1)
			{
				int hour = s / 12;
				int minute = (int)((s / 12.0 - hour) * 60 + 0.5);

				double value = 0;


				sprintf(time_interval_field_name, "T%02d%02d", hour, minute);

				if (parser.GetValueByFieldName(time_interval_field_name, value, false) == true)
				{
					dep_time.departure_time_ratio[s] = value;
				}

				if (s == dep_time.starting_time_slot_no || s == dep_time.ending_time_slot_no)
				{

					dtalog.output() << "[DATA INFO] At time " << hour << ":" << minute << " (" << time_interval_field_name <<") , the demand ratio is " << value << "." << '\n';
					g_DTA_log_file << "[DATA INFO] At time " << hour << ":" << minute << " (" << time_interval_field_name <<") , the demand ratio is " << value << "." << '\n';

				}

			}

			dep_time.compute_cumulative_profile(dep_time.starting_time_slot_no, dep_time.ending_time_slot_no, true);
			assignment.g_DepartureTimeProfileVector.push_back(dep_time);


		}

	}
}

double g_pre_read_demand_file(Assignment& assignment)
{

	float total_demand_in_demand_file = 0;

	CDTACSVParser parser;
	

	assignment.summary_file << "pre-read demand_file_list.csv." << '\n';

	//The aim of preloading the demand matrix under the subarea scenario is to gain insights into the structure of the Origin - Destination(OD) demand
	//	across the entire network.By doing so, we can identify critical OD pairsand zones based on their relevanceand importance.This understanding aids in the subsequent loading 
	//	of the demand matrix, enabling a more focusedand efficient approach.Consequently, this process allows for a significant reduction in the number of zones, 
	//	enhancing the efficiencyand performance of our traffic assignment tasks*/


	if (parser.OpenCSVFile("demand_file_list.csv", false))
	{
		while (parser.ReadRecord())
		{
			int this_departure_time_profile_no = 0;

			int file_sequence_no = 1;

			string format_type = "null";

			int demand_format_flag = 0;

			if (!parser.GetValueByFieldName("file_sequence_no", file_sequence_no))
				break;

			// skip negative sequence no
			if (file_sequence_no <= -1)
				continue;

			double loading_scale_factor = 1.0;
			string file_name, demand_period_str, mode_type;
			parser.GetValueByFieldName("file_name", file_name);
			parser.GetValueByFieldName("demand_period", demand_period_str);
			parser.GetValueByFieldName("format_type", format_type);


			if (format_type.find("null") != string::npos)  // skip negative sequence no
			{
				dtalog.output() << "[ERROR] Please provide format_type in section [demand_file_list.]" << '\n';
				g_DTA_log_file << "[ERROR] Please provide format_type in section [demand_file_list.]" << '\n';
				format_type = "csv"; 
			}


			if (format_type.find("column") != string::npos || format_type.find("bin") != string::npos)  // or muliti-column
			{

				// try to detect if we have a route.csv as preload file, if yes, we skip the following reading .
				CDTACSVParser parser_route;

				struct SDemandHeader
				{
					int o_zone_id, d_zone_id, mode_type_no, demand_period;
					double volume;
				};

				SDemandHeader header;


				CDTACSVParser parser;
				FILE* pFile;
				int file_exists = 0;
				int line_no = 0;
				int file_format = 0;

				fopen_ss(&pFile, "demand.bin", "rb");
				if (pFile != NULL)
				{
					file_format = 2;
					dtalog.output() << "[STATUS INFO] pre-reading demand.bin in fast binary file reading mode." << '\n';
					g_DTA_log_file << "[STATUS INFO] pre-reading demand.bin in fast binary file reading mode." << '\n';
				}


				if (file_format == 0 && format_type.find("column") != string::npos)
				{
					if (parser.OpenCSVFile(file_name, false) == false)
					{
						dtalog.output() << "[ERROR] column file " << file_name.c_str() << " does not exist." << '\n';
						g_DTA_log_file << "[ERROR] column file " << file_name.c_str() << " does not exist." << '\n';
						break;

					}
					file_format = 1;
				}

				int error_count = 0;
				int critical_OD_count = 0;
				double critical_OD_volume = 0;
				// read the file formaly after the test.
				while (file_format >= 1)
				{
					int o_zone_id, d_zone_id;
					double demand_value = 0;

					if (file_format == 1)  // column
					{
						int flag = parser.ReadRecord();
						if (flag == 0)
							break;

						parser.GetValueByFieldName("o_zone_id", o_zone_id);
						parser.GetValueByFieldName("d_zone_id", d_zone_id);
						parser.GetValueByFieldName("volume", demand_value);

					}
					if (file_format == 2)  // binary
					{
						if (feof(pFile))
							break;

						size_t result;

						fread(&header, sizeof(header), 1, pFile);
						o_zone_id = header.o_zone_id;
						d_zone_id = header.d_zone_id;
						demand_value = header.volume;

					}

					// end of pre-reading  each record

					if (assignment.g_zoneid_to_zone_seq_no_mapping.find(o_zone_id) == assignment.g_zoneid_to_zone_seq_no_mapping.end())
					{
						// origin zone has not been defined, skipped.
						continue;
					}

					if (assignment.g_zoneid_to_zone_seq_no_mapping.find(d_zone_id) == assignment.g_zoneid_to_zone_seq_no_mapping.end())
					{
						continue;
					}

					int from_zone_seq_no = 0;
					int to_zone_seq_no = 0;
					from_zone_seq_no = assignment.g_zoneid_to_zone_seq_no_mapping[o_zone_id];
					to_zone_seq_no = assignment.g_zoneid_to_zone_seq_no_mapping[d_zone_id];

					g_zone_vector[from_zone_seq_no].preread_total_O_demand += demand_value;
					g_zone_vector[from_zone_seq_no].preread_ODdemand[to_zone_seq_no] += demand_value;

					total_demand_in_demand_file += demand_value;
					// encounter return
					if (demand_value < -99)
						break;


					line_no++;

					if (line_no % 1000000 == 0)
					{
						int line_no_output = line_no / 1000000;
						dtalog.output() << "[STATUS INFO] Pre-Reading demand file line no =  " << line_no_output << " million lines" << '\n';
						g_DTA_log_file << "[STATUS INFO] Pre-Reading demand file line no =  " << line_no_output << " million lines" << '\n';
					}
				}  // scan lines


			// pre-read complete.
				break;

			}
		}
	}
	return total_demand_in_demand_file;
}

int g_load_demand_from_route_file_in_settings()
{
	CDTACSVParser parser;

	if (parser.OpenCSVFile("demand_file_list.csv", false))
	{
		while (parser.ReadRecord())
		{
			int this_departure_time_profile_no = 0;

			int file_sequence_no = 1;

			string format_type = "null";

			int demand_format_flag = 0;

			if (!parser.GetValueByFieldName("file_sequence_no", file_sequence_no))
				break;

			// skip negative sequence no
			if (file_sequence_no <= -1)
				continue;

			double loading_scale_factor = 1.0;
			string file_name, demand_period_str, mode_type;
			parser.GetValueByFieldName("file_name", file_name);
			parser.GetValueByFieldName("format_type", format_type);


			if (format_type.compare("route") == 0)
			{
				parser.CloseCSVFile();
				return 1;
			}
		}
		parser.CloseCSVFile();

	}
	return 0;
}
void g_check_demand_volume_with_mode_type_log(Assignment& assignment)
{
	// Set the column widths and precision for formatting
	const int fieldWidth = 15;
	dtalog.output() << "[PROCESS INFO]"; 
	g_DTA_log_file << "[PROCESS INFO]"; 
	dtalog.output() << '\n';
	g_DTA_log_file << '\n';

	dtalog.output() << "";
	g_DTA_log_file << "";
	dtalog.output() << std::setw(fieldWidth) << "demand_period and mode_type";
	g_DTA_log_file << std::setw(fieldWidth) << "demand_period and mode_type";
	dtalog.output() << std::setw(fieldWidth) << "total_demand";
	g_DTA_log_file << std::setw(fieldWidth) << "total_demand";
	dtalog.output() << std::setw(fieldWidth) << "#_of_links";
	g_DTA_log_file << std::setw(fieldWidth) << "#_of_links";
	dtalog.output() << std::setw(fieldWidth) << "speed_mph";
	g_DTA_log_file << std::setw(fieldWidth) << "speed_mph";
	dtalog.output() << std::setw(fieldWidth) << "speed_kmph";
	g_DTA_log_file << std::setw(fieldWidth) << "speed_kmph";
	dtalog.output() << std::setw(fieldWidth) << "length km";
	g_DTA_log_file << std::setw(fieldWidth) << "length km";
	dtalog.output() << std::setw(fieldWidth) << "avg_lane_cap";
	g_DTA_log_file << std::setw(fieldWidth) << "avg_lane_cap";
	dtalog.output() << std::setw(fieldWidth) << "avg_length_m";
	g_DTA_log_file << std::setw(fieldWidth) << "avg_length_m";
	dtalog.output() << '\n';
	g_DTA_log_file << '\n';

	for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); ++tau)
	{
		for (int at = 0; at < assignment.g_ModeTypeVector.size(); at++)
		{
			dtalog.output() << std::setw(fieldWidth) << assignment.g_DemandPeriodVector[tau].demand_period.c_str() << "," << assignment.g_ModeTypeVector[at].mode_type.c_str();
			g_DTA_log_file << std::setw(fieldWidth) << assignment.g_DemandPeriodVector[tau].demand_period.c_str() << "," << assignment.g_ModeTypeVector[at].mode_type.c_str();
			dtalog.output() << std::setw(fieldWidth) << assignment.total_demand[at][tau] << "";
			g_DTA_log_file << std::setw(fieldWidth) << assignment.total_demand[at][tau] << "";
		

			int link_count = 0;
			double total_speed = 0;
			double total_length = 0;
			double total_lane_capacity = 0;
			double total_link_capacity = 0;

			for (int i = 0; i < g_link_vector.size(); i++)
			{
				if (g_link_vector[i].link_type_si[0] >= 0 && g_link_vector[i].AllowModeType(assignment.g_ModeTypeVector[at].mode_type, tau, assignment.active_scenario_index))
				{
					link_count++;
					total_speed += g_link_vector[i].free_speed;
					total_length += g_link_vector[i].length_in_meter * g_link_vector[i].number_of_lanes_si[assignment.active_scenario_index];
					total_lane_capacity += g_link_vector[i].lane_capacity;
					total_link_capacity += g_link_vector[i].lane_capacity * g_link_vector[i].number_of_lanes_si[assignment.active_scenario_index];
				}
			}


			const int precision = 2;

			dtalog.output() << std::setw(fieldWidth) << link_count;
			g_DTA_log_file << std::setw(fieldWidth) << link_count;
			dtalog.output() << std::setw(fieldWidth) << std::fixed << std::setprecision(precision) << (total_speed / max(1, link_count));
			g_DTA_log_file << std::setw(fieldWidth) << std::fixed << std::setprecision(precision) << (total_speed / max(1, link_count));
			dtalog.output() << std::setw(fieldWidth) << std::fixed << std::setprecision(precision) << (total_speed / max(1, link_count) * 1.60934);
			g_DTA_log_file << std::setw(fieldWidth) << std::fixed << std::setprecision(precision) << (total_speed / max(1, link_count) * 1.60934);
			dtalog.output() << std::setw(fieldWidth) << std::fixed << std::setprecision(precision) << (total_length / 1000.0);
			g_DTA_log_file << std::setw(fieldWidth) << std::fixed << std::setprecision(precision) << (total_length / 1000.0);
			dtalog.output() << std::setw(fieldWidth) << std::fixed << std::setprecision(precision) << (total_lane_capacity / max(1, link_count));
			g_DTA_log_file << std::setw(fieldWidth) << std::fixed << std::setprecision(precision) << (total_lane_capacity / max(1, link_count));
			dtalog.output() << std::setw(fieldWidth) << std::fixed << std::setprecision(precision) << (total_length / max(1, link_count)) << '\n';
			g_DTA_log_file << std::setw(fieldWidth) << std::fixed << std::setprecision(precision) << (total_length / max(1, link_count)) << '\n';

			if (assignment.total_demand[at][tau] > 0 && link_count == 0)
			{
				dtalog.output() << "[ERROR] Demand period = " << assignment.g_DemandPeriodVector[tau].demand_period.c_str() << ", mode_type = " << assignment.g_ModeTypeVector[at].mode_type.c_str() << " has a total demand of " << assignment.total_demand[at][tau] << ", but this mode type is not allowed on any link in the network." << '\n';
				g_DTA_log_file << "[ERROR] Demand period = " << assignment.g_DemandPeriodVector[tau].demand_period.c_str() << ", mode_type = " << assignment.g_ModeTypeVector[at].mode_type.c_str() << " has a total demand of " << assignment.total_demand[at][tau] << ", but this mode type is not allowed on any link in the network." << '\n';
				dtalog.output() << "Please check the link_type.csv file to see if the 'allowed_uses_p1' field allows this type of travel demand for major facilities, such as real-time information or HOV for highway facilities." << '\n';
				g_DTA_log_file << "Please check the link_type.csv file to see if the 'allowed_uses_p1' field allows this type of travel demand for major facilities, such as real-time information or HOV for highway facilities." << '\n';

				g_program_stop();
			}
		}
	}


}
void g_check_demand_volume_with_mode_type(Assignment& assignment)
{
	assignment.summary_file << ",summary by multi-modal and demand types,demand_period,mode_type,total demmand volume in scenario 0, # of links,avg_free_speed_mph,avg_free_speed_kmph,total_length_in_km,total_capacity,avg_lane_capacity,avg_length_in_meter," << '\n';

	for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); ++tau)
		for (int at = 0; at < assignment.g_ModeTypeVector.size(); at++)
		{
			assignment.summary_file << ",," << assignment.g_DemandPeriodVector[tau].demand_period.c_str() << ",";
			assignment.summary_file << assignment.g_ModeTypeVector[at].mode_type.c_str() << "," << assignment.total_demand[at][tau] << ",";

			int link_count = 0;
			double total_speed = 0;
			double total_length = 0;
			double total_lane_capacity = 0;
			double total_link_capacity = 0;

			for (int i = 0; i < g_link_vector.size(); i++)
			{
				if (g_link_vector[i].link_type_si[0] >= 0 && g_link_vector[i].AllowModeType(assignment.g_ModeTypeVector[at].mode_type, tau, assignment.active_scenario_index))
				{
					link_count++;
					total_speed += g_link_vector[i].free_speed;
					total_length += g_link_vector[i].length_in_meter * g_link_vector[i].number_of_lanes_si[assignment.active_scenario_index];
					total_lane_capacity += g_link_vector[i].lane_capacity;
					total_link_capacity += g_link_vector[i].lane_capacity * g_link_vector[i].number_of_lanes_si[assignment.active_scenario_index];
				}
			}
			assignment.summary_file << link_count << "," <<
				total_speed / max(1, link_count) << "," <<
				total_speed / max(1, link_count) * 1.60934 << "," <<
				total_length / 1000.0 << "," <<
				total_link_capacity << "," <<
				total_lane_capacity / max(1, link_count) << "," << total_length / max(1, link_count) << "," << '\n';

		}

	g_check_demand_volume_with_mode_type_log(assignment);
}


void g_create_subarea_related_zone_structure(int max_num_significant_zones_in_subarea, int max_num_significant_zones_outside_subarea)
{

	//The code is analyzing various zones for demand, determining which zones are "inside" or "outside" a certain area, and then adjusting the zone list based on certain criteria.
	//	Here's a summary of what the code does:
	//	1.	First, the code initializes a number of variables representing various kinds of demandand zone counts.
	//	2.	The code then preloads an Origin - Destination(OD) demand for the simulation.If the preloaded demand is negligible, the program stops.
	//	3.	If the preloaded demand is not negligible, the code goes through all origin zones.
	//	4.	For each origin zone, it checks whether the zone is within a specified subarea.If it is, the total demand variables are updated, the zone's status is updated, and the inside zone count is increased.
	//	5.	If the zone is not within the subarea, the code checks if a line from this origin to every other destination intersects the subarea.If it does, and if there is demand from the origin to that destination, the total demandand zone impact volume variables are updated, and the zone's status is updated as a related external zone.
	//	6.	After processing all destinations for a particular origin, if the zone impact volume is below a certain threshold, the zone is classified as a non - impact zone.Otherwise, it is marked as a significant related external zone.
	//	7.	After all the zones have been processed, the code logs the countsand total demand associated with various zone classifications.
	//	8.	The code then calculates a cutoff for the number of origin zones to keep, based on the number of zones with significant external volume.
	//	9.	If the number of zones is greater than the cutoff, the code sorts the zones based on their impact volume, determines a volume cutoff value, and resets the status of zones that fall below this volume cutoff value.
	//	10.	The code then logs the cutoff zone size, the total demand from the zones that made the cutoff, the volume cutoff value used to determine significance, and the ratio of remaining external related demand after the cutoff.
	//	This code effectively focuses the simulation on zonesand routes that have a significant impact on the demand within the designated subarea, reducing the computational load by ignoring zonesand routes that do not meet certain significance criteria.


	//Here's an attempt to simplify the process with a sort of mathematical notation.
	//	Let's denote the following:
	//	•	Z be the set of all zones
	//	•	I be the set of "inside" zones, such that I ⊆ Z
	//	•	E be the set of related external zones, such that E ⊆ Zand I ∩ E = ∅
	//	•	S be the set of significant external zones, such that S ⊆ E
	//	•	D(o, d) be the demand from zone o to zone d
	//	•	V(z) be the volume impact of zone z
	//	•	R be the ratio of origin zone impact volume to total demand
	//	The key steps in the process can then be represented as follows :
	//1.	Preload total demand :
	//D_total = Σ[D(o, d) for all o, d ∈ Z]
	//	2.	Check each zone to see if it's "inside" or "external":
	//	I = { o ∈ Z : o is inside the subarea }
	//	E = { o ∈ Z : o is not in I and there exists d ∈ Z such that a line from o to d intersects the subarea }
	//	3.	Calculate the total demand and impact volume associated with each zone :
	//D_I = Σ[D(o, d) for o ∈ Iand for all d ∈ Z]
	//	D_E = Σ[D(o, d) for o ∈ Eand for all d ∈ Z]
	//	V(o) = Σ[D(o, d) for o ∈ Zand for all d ∈ Z if a line from o to d intersects the subarea]
	//	4.	Classify each external zone as significant or non - impact based on its impact volume :
	//S = { o ∈ E : V(o) ≥ threshold_value }
	//	non_impact_zones = E - S
	//	5.	Adjust the set of zones based on a volume cutoff value :
	//If size(S) > cutoff_value, sort S by V(o), find a V(cutoff), and reclassify zones in S :
	//S = { o ∈ S : V(o) ≥ V(cutoff) }
	//	non_impact_zones = non_impact_zones ∪{ o ∈ S : V(o) < V(cutoff) }



	double total_demand = 0;
	double total_related_demand = 0;
	double total_related_demand_inside_signifant = 0;
	double total_related_demand_inside_insignifant = 0;

	double total_related_demand_from_internal = 0;
	double total_related_demand_from_external = 0;
	double total_related_demand_from_external_cutoff = 0;
	double volume_cut_off_value = 0;

	int non_impact_zone_count = 0;
	int inside_significant_zone_count = 0;
	int inside_insignificant_zone_count = 0;

	int related_external_zone_count = 0;
	int related_external_zone_count_with_significant_volume = 0;
	// pre read demand

	dtalog.output() << "[PROCESS INFO] Step 1.8: Preloading Origin-Destination (OD) demand data as part of the focusing approach. " <<  '\n';
	g_DTA_log_file << "[PROCESS INFO] Step 1.8: Preloading Origin-Destination (OD) demand data as part of the focusing approach. " <<  '\n';

	double total_preload_demand = g_pre_read_demand_file(assignment);



	if (total_preload_demand <= 0.001)
	{
		dtalog.output() << "[ERROR] demand is missing " <<'\n';
		g_DTA_log_file << "[ERROR] demand is missing " <<'\n';
		g_program_stop();
	}
	else
	{
		// first stage: scan all origin zones
		// output: g_zone_vector[orig].subarea_significance_flag: true for inside zone, and related zone
		// total_related_demand: inside to all, outside passing to all
		//assignment.log_subarea_focusing_file << "[PROCESS INFO] Starting 0 stage: Scanning all origin zones, determine the cut off for the zones inside...\n";

		std::vector <float> significant_subarea_inside_origin_zone_volumes;

		//// Loop through all origin zones


		for (int orig = 0; orig < g_zone_vector.size(); orig++)  // o
		{
			DTAGDPoint Pt;
			Pt.x = g_zone_vector[orig].cell_x;
			Pt.y = g_zone_vector[orig].cell_y;

			total_demand += g_zone_vector[orig].preread_total_O_demand;
			// Check if the zone lies inside the subarea

			if (g_test_point_in_polygon(Pt, assignment.g_subarea_shape_points) == 1)
			{
				significant_subarea_inside_origin_zone_volumes.push_back(g_zone_vector[orig].preread_total_O_demand);
			}


		}

		 //Check if the vector has elements

			// Sort the vector in ascending order, this is a simple bubble sorting method 
			std::sort(significant_subarea_inside_origin_zone_volumes.begin(), significant_subarea_inside_origin_zone_volumes.end());

			int subarea_cutoff_zone_size = max_num_significant_zones_in_subarea;
			// Determine the cut-off index for the inside zones 
			int cut_off_subarea_inside_zone_index = (significant_subarea_inside_origin_zone_volumes.size() - subarea_cutoff_zone_size <= 0) ? 0 : significant_subarea_inside_origin_zone_volumes.size() - subarea_cutoff_zone_size;
			// Get the volume at the cut-off index
			float volume_cut_off_value_subarea_inside_zone = 0;
			if (cut_off_subarea_inside_zone_index >= 0 && cut_off_subarea_inside_zone_index < significant_subarea_inside_origin_zone_volumes.size() - 1)  // to prevent out of ranage error
			{
				volume_cut_off_value_subarea_inside_zone = significant_subarea_inside_origin_zone_volumes[cut_off_subarea_inside_zone_index];
			}


		assignment.log_subarea_focusing_file << "[PROCESS INFO] Starting 1st stage: Scanning all origin zones...\n";

		for (int orig = 0; orig < g_zone_vector.size(); orig++)  // for zone 
		{
			DTAGDPoint Pt;
			Pt.x = g_zone_vector[orig].cell_x;
			Pt.y = g_zone_vector[orig].cell_y;

			total_demand += g_zone_vector[orig].preread_total_O_demand;
			// Check if the zone lies inside the subarea and the 

			if (g_test_point_in_polygon(Pt, assignment.g_subarea_shape_points) == 1 )
			{
				if(g_zone_vector[orig].preread_total_O_demand >= volume_cut_off_value_subarea_inside_zone)  // volume is greater than the cutoff volume  for inside significant
				{ 
				total_related_demand_inside_signifant += g_zone_vector[orig].preread_total_O_demand;
				total_related_demand_from_internal += g_zone_vector[orig].preread_total_O_demand;
				total_related_demand += g_zone_vector[orig].preread_total_O_demand;

				g_zone_vector[orig].origin_zone_impact_volume = g_zone_vector[orig].preread_total_O_demand;
				g_zone_vector[orig].subarea_inside_flag = 3;  // zone is inside the subarea
				inside_significant_zone_count++;
				}
				else  // inside but not significant 
				{
					total_related_demand_inside_insignifant += g_zone_vector[orig].preread_total_O_demand;
					total_related_demand += g_zone_vector[orig].preread_total_O_demand;
					g_zone_vector[orig].subarea_inside_flag = 4;  // zone is inside the subarea but not significant 
					inside_insignificant_zone_count++;
				}

				continue; //inside the subarea, keep the zone anyway

			}

			// Continue processing if the zone is outside the subarea
			assignment.log_subarea_focusing_file << "[PROCESS INFO] inside_significant_zone_count =  " << inside_significant_zone_count << "n";
			assignment.log_subarea_focusing_file << "[PROCESS INFO] inside_insignificant_zone_count =  " << inside_insignificant_zone_count << "n";



			double origin_zone_impact_volume = 0;
			// determine impact OD 
			for (int dest = 0; dest < g_zone_vector.size(); dest++) //d
			{
				if (orig == dest)
					continue;

				double Ax, Ay;
				double Bx, By;

				Ax = g_zone_vector[orig].cell_x;
				Ay = g_zone_vector[orig].cell_y;

				Bx = g_zone_vector[dest].cell_x;
				By = g_zone_vector[dest].cell_y;

				bool intersect_flag = g_get_line_polygon_intersection(Ax, Ay, Bx, By, assignment.g_subarea_shape_points);  //test if OD pair passes through the subarea

				if (intersect_flag == true)
				{
					if (g_zone_vector[orig].preread_ODdemand.find(dest) != g_zone_vector[orig].preread_ODdemand.end())
					{
						total_related_demand += g_zone_vector[orig].preread_ODdemand[dest];
						total_related_demand_from_external += g_zone_vector[orig].preread_ODdemand[dest];
						origin_zone_impact_volume += g_zone_vector[orig].preread_ODdemand[dest];
						g_zone_vector[orig].preread_total_O_related_demand += g_zone_vector[orig].preread_ODdemand[dest];
						g_zone_vector[orig].impact_passing_ODdemand[dest] = g_zone_vector[orig].preread_ODdemand[dest];;
						if (g_zone_vector[orig].subarea_inside_flag == 0)
						{
							g_zone_vector[orig].subarea_inside_flag = 2;
							g_zone_vector[orig].subarea_significance_flag = true;
						}
					}
				}
			}  // end of destination zone loop

			related_external_zone_count++;
			assignment.log_subarea_focusing_file << "[PROCESS INFO] Zone " << orig << " processed. Related external zone count: " << related_external_zone_count << '\n';

			double impact_demand_ratio = origin_zone_impact_volume / max(0.0001, g_zone_vector[orig].preread_total_O_demand);
			//			if (origin_zone_impact_volume < 5 && impact_demand_ratio <0.1)

			g_zone_vector[orig].origin_zone_impact_volume = origin_zone_impact_volume;

			int threadshold_cut_off = 5;

			if (origin_zone_impact_volume < threadshold_cut_off)
			{
				g_zone_vector[orig].subarea_inside_flag = 0;
				g_zone_vector[orig].subarea_significance_flag = false;  //reset this as no signalized zone
				non_impact_zone_count++;
				// if there are more than y = 5 vehicles from a zone passing through the subarea(with the mark from step 2),
				//then we keep them in shortest path computing.
			}
			else
			{
				related_external_zone_count_with_significant_volume++;
			}

		} // end of origin zone loop
		assignment.log_subarea_focusing_file << "[PROCESS INFO] Finished scanning all origin zones\n";

		assignment.summary_file << ",first stage scanning" << '\n';
		assignment.summary_file << ",# of inside zones = " << inside_significant_zone_count+ inside_insignificant_zone_count << ", total_related_demand = " << total_related_demand << '\n';
		assignment.summary_file << ",# inside_significant_zone_count = " << inside_significant_zone_count << ", total_related_demand_from_internal zones = " << total_related_demand_inside_signifant << '\n';
		assignment.summary_file << ",# inside_insignificant_zone_count = " << inside_insignificant_zone_count << ", total_related_demand_from_internal zones = " << total_related_demand_inside_insignifant << '\n';
		assignment.summary_file << ",# related_external_zone_count = " << related_external_zone_count << ", external related demand = " << total_related_demand_from_external << '\n';
		assignment.summary_file << ",# related_external_zone_count_with_significant_volume = " << related_external_zone_count_with_significant_volume << '\n';

		dtalog.output() << "[DATA INFO] Initiating the first stage of the subarea scanning process, where we determine cut-off zones.\n";
		g_DTA_log_file << "[DATA INFO] Initiating the first stage of the subarea scanning process, where we determine cut-off zones.\n";

		// Information about the number of inside zones and total related demand.
		dtalog.output() << "[DATA INFO] Total demand related to these zones : " << total_related_demand << ".\n";
		g_DTA_log_file << "[DATA INFO] Total demand related to these zones : " << total_related_demand << ".\n";

		dtalog.output() << "[DATA INFO] Number of significant zones inside the subarea: " << inside_significant_zone_count << ". Total related demand inside significant zones: " << total_related_demand_inside_signifant << ".\n";
		g_DTA_log_file << "[DATA INFO] Number of significant zones inside the subarea: " << inside_significant_zone_count << ". Total related demand inside significant zones: " << total_related_demand_inside_signifant << ".\n";
		float ratio_significant = inside_significant_zone_count / max(1, inside_significant_zone_count + inside_insignificant_zone_count); 
		dtalog.output() << "[DATA INFO] Number of insignificant zones inside the subarea: " << inside_insignificant_zone_count << ". Total related demand inside insignificant zones: " << total_related_demand_inside_insignifant << ".\n";
		g_DTA_log_file << "[DATA INFO] Number of insignificant zones inside the subarea: " << inside_insignificant_zone_count << ". Total related demand inside insignificant zones: " << total_related_demand_inside_insignifant << ".\n";
		dtalog.output() << "[DATA INFO] Ratio of total significant zone to  zones inside the subarea: " << ratio_significant << ".\n";
		g_DTA_log_file << "[DATA INFO] Ratio of total significant zone to  zones inside the subarea: " << ratio_significant << ".\n";
		// Information about the number of inside zones and total related demand from internal zones.
		// Information about the number of external zones related to the subarea and total related demand from external zones.
		dtalog.output() << "[DATA INFO] Number of external zones related to the subarea: " << related_external_zone_count << ". Total demand related to the subarea from these external zones: " << total_related_demand_from_external << ".\n";
		g_DTA_log_file << "[DATA INFO] Number of external zones related to the subarea: " << related_external_zone_count << ". Total demand related to the subarea from these external zones: " << total_related_demand_from_external << ".\n";

		// Information about the number of external zones related to the subarea with significant volume.
		dtalog.output() << "[DATA INFO] Number of external zones related to the subarea with significant volume: " << related_external_zone_count_with_significant_volume << ".\n";
		g_DTA_log_file << "[DATA INFO] Number of external zones related to the subarea with significant volume: " << related_external_zone_count_with_significant_volume << ".\n";

		// use the following rule to determine how many zones to keep as origin zones, we still keep all destination zones
		int cutoff_zone_size = min(g_zone_vector.size(), static_cast<size_t>(max(related_external_zone_count_with_significant_volume * 0.25, 100.0)));

		if (cutoff_zone_size > max_num_significant_zones_outside_subarea)
			cutoff_zone_size = max_num_significant_zones_outside_subarea;

		// Check if the number of zones remaining after the cut-off and inside zones are still greater than total zones
		if ((cutoff_zone_size + inside_significant_zone_count) < g_zone_vector.size()) {

			dtalog.output() << "[PROCESS INFO] Handling additional zones...\n";
			g_DTA_log_file << "[PROCESS INFO] Handling additional zones...\n";

			// Vector to hold the volume of significant origin zones
			std::vector <float> significant_origin_zone_volumes;

			// Loop through all origin zones
			for (int orig = 0; orig < g_zone_vector.size(); orig++) {
				// Check if the zone is significant and its inside flag is set to 2
				if (g_zone_vector[orig].subarea_significance_flag == true && g_zone_vector[orig].subarea_inside_flag == 2) {
					// Add the volume of the zone to the vector
					significant_origin_zone_volumes.push_back(g_zone_vector[orig].origin_zone_impact_volume);
				}
			}

			// Check if the vector has elements
			if (significant_origin_zone_volumes.size() > 0) {
				dtalog.output() << "[PROCESS INFO] Sorting zone volumes...\n";
				g_DTA_log_file << "[PROCESS INFO] Sorting zone volumes...\n";

				// Sort the vector in ascending order
				std::sort(significant_origin_zone_volumes.begin(), significant_origin_zone_volumes.end());

				// Determine the cut-off index
				int cut_off_zone_index = (significant_origin_zone_volumes.size() - cutoff_zone_size <= 0) ? 0 : significant_origin_zone_volumes.size() - cutoff_zone_size;

				// Get the volume at the cut-off index
				float volume_cut_off_value = significant_origin_zone_volumes[cut_off_zone_index];

				// Loop through all origin zones again
				int log_non_significant_count = 0; 
				int log_significantly_related_count = 0; 
				for (int orig = 0; orig < g_zone_vector.size(); orig++) {
					// Check if the zone is significant and its inside flag is set to 2
					if (g_zone_vector[orig].subarea_inside_flag == 2 && g_zone_vector[orig].subarea_significance_flag == true) {
						// If the volume of the zone is less than the cut-off value, set the zone to non-significant
						if (g_zone_vector[orig].origin_zone_impact_volume < volume_cut_off_value && log_non_significant_count <3) {
							dtalog.output() << "[PROCESS INFO] Zone " << orig << " is non-significant.\n";
							g_DTA_log_file << "[PROCESS INFO] Zone " << orig << " is non-significant.\n";
							log_non_significant_count++;
							g_zone_vector[orig].subarea_significance_flag = false;
							g_zone_vector[orig].subarea_inside_flag = 1;  // zone is set to boundary
						}
						else {
							if(log_significantly_related_count<3)
							{
							dtalog.output() << "[PROCESS INFO] Zone " << orig << " is significantly related.\n";
							g_DTA_log_file << "[PROCESS INFO] Zone " << orig << " is significantly related.\n";
							log_significantly_related_count++; 
							}

							g_zone_vector[orig].subarea_inside_flag = 2;  // zone is set to significantly related
							total_related_demand_from_external_cutoff += g_zone_vector[orig].preread_total_O_related_demand;
						}
					}
				}
			}
		} // End of handling additional zones





		assignment.summary_file << " During the second stage of cut off, the cut off zone includes all inside zones along with significant external related zones." << '\n';
		assignment.summary_file << ",# cut off zone size = " << cutoff_zone_size << ", total_related_demand_from_external_cutoff  = " << total_related_demand_from_external_cutoff << '\n';
		assignment.summary_file << ",cut off volume threshold to determine significance = " << volume_cut_off_value << ".External origin zones with a volume larger than this threshold will be retained." << '\n';
		double cut_off_demand_ratio = total_related_demand_from_external_cutoff / max(0.001, total_related_demand_from_external);
		assignment.summary_file << ", remaining external related demand percentage after the cut off  = " << cut_off_demand_ratio * 100 << " %" << '\n';


		dtalog.output() << "[DATA INFO]" << " During the second stage of cut off, the cut off zone includes all inside zones along with significant external related zones " << '\n';
		g_DTA_log_file << "[DATA INFO]" << " During the second stage of cut off, the cut off zone includes all inside zones along with significant external related zones " << '\n';
		dtalog.output() << "[DATA INFO]" << " # cut-off zone size = " << cutoff_zone_size << ",  total related demand from external cut-off  = " << total_related_demand_from_external_cutoff << '\n';
		g_DTA_log_file << "[DATA INFO]" << " # cut-off zone size = " << cutoff_zone_size << ",  total related demand from external cut-off  = " << total_related_demand_from_external_cutoff << '\n';
		dtalog.output() << "[DATA INFO]" << " Cut-off volume threshold to determine significance = " << volume_cut_off_value << ". External origin zones with a volume larger than this threshold will be retained." << '\n';
		g_DTA_log_file << "[DATA INFO]" << " Cut-off volume threshold to determine significance = " << volume_cut_off_value << ". External origin zones with a volume larger than this threshold will be retained." << '\n';
		cut_off_demand_ratio = total_related_demand_from_external_cutoff / max(0.001, total_related_demand_from_external);
		dtalog.output() << "[DATA INFO]" << " Remaining external related demand percentage after the cut-off:  = " << cut_off_demand_ratio * 100 << " %" << '\n';
		g_DTA_log_file << "[DATA INFO]" << " Remaining external related demand percentage after the cut-off:  = " << cut_off_demand_ratio * 100 << " %" << '\n';

	}

	// we need additional handlingh here to select most important zones inside the subarea 
	// 
	// 		number of super zones
	int sindex = 0;
	for (int orig = 0; orig < g_zone_vector.size(); orig++)  // o
	{
		if (g_zone_vector[orig].subarea_significance_flag == false)  // no significant: skip
		{
			g_zone_vector[orig].sindex = -1;

		}
		else
		{
			g_zone_vector[orig].sindex = sindex++;  // resort the index
		}

	}

	if (sindex > 0)
	{
		g_related_zone_vector_size = sindex;
	}
	else
	{

		// no zone intersects with the given subarea area, return back to the normal mode 
		int sindex = 0;
		for (int orig = 0; orig < g_zone_vector.size(); orig++)  // o
		{
			g_zone_vector[orig].subarea_significance_flag = true; 
			g_zone_vector[orig].sindex = sindex++;  // resort the index

		}
	}


	////
	if (assignment.g_subarea_shape_points.size() >= 3) // if there is a subarea defined.
	{

		if (assignment.g_number_of_analysis_districts <= 1)  // no external distriction is given
		{// generate analaysis district
			// output:
			//assignment.g_zone_seq_no_to_analysis_distrct_id_mapping
			//assignment.g_number_of_analysis_districts = zone_id_to_analysis_district_id_mapping[ozone.zone_id] + 1;

			double distrct_ratio = 20.0 / g_related_zone_vector_size;

			int distrct_index = 1;

			for (int orig = 0; orig < g_zone_vector.size(); orig++)  // o
			{
				if (g_zone_vector[orig].sindex >= 0)  // no significant: skip
				{
					if (g_get_random_ratio() < distrct_ratio)
					{
						assignment.g_zone_seq_no_to_analysis_distrct_id_mapping[orig] = distrct_index;
						g_zone_vector[orig].b_distrct_cluster_seed = true;
						g_zone_vector[orig].analysis_district_index = distrct_index;

						distrct_index++;
					}
					else
					{
						g_zone_vector[orig].analysis_district_index = -1;
					}

				}

			}
			assignment.g_number_of_analysis_districts = distrct_index;

			//clustering
			for (int orig = 0; orig < g_zone_vector.size(); orig++)  // o
			{
				if (g_zone_vector[orig].analysis_district_index == -1)
				{
					double min_distance = 99999;

					for (int d = 0; d < g_zone_vector.size(); d++)  // o
					{
						if (g_zone_vector[d].b_distrct_cluster_seed == true)
						{
							double delta_x = g_zone_vector[orig].cell_x - g_zone_vector[d].cell_x;
							double delta_y = g_zone_vector[orig].cell_y - g_zone_vector[d].cell_y;

							double local_distance = pow(delta_x * delta_x + delta_y * delta_y, 0.5);

							if (local_distance < min_distance)
							{
								min_distance = local_distance;
								g_zone_vector[orig].analysis_district_index = g_zone_vector[d].analysis_district_index;
								assignment.g_zone_seq_no_to_analysis_distrct_id_mapping[orig] = g_zone_vector[d].analysis_district_index;

							}

						}
					}

				}

			}


		}

		// end of generating analaysis district


		if (g_related_zone_vector_size >= 30000)  // do not use this function for now.
		{
			// k mean method
			//random select 300 zones then find the closest zone to aggregate

			double super_zone_ratio = 300.0 / g_related_zone_vector_size;

			int super_zone_index = assignment.g_number_of_analysis_districts + 1;  // automated district  id start from the externally given district id

			for (int orig = 0; orig < g_zone_vector.size(); orig++)  // o
			{
				if (g_zone_vector[orig].sindex >= 0)  // no significant: skip
				{
					if (g_get_random_ratio() < super_zone_ratio)
					{
						g_zone_vector[orig].superzone_index = super_zone_index;
						g_zone_vector[orig].bcluster_seed = true;
						super_zone_index++;
					}
					else
					{
						g_zone_vector[orig].superzone_index = -1;
						g_zone_vector[orig].b_shortest_path_computing_flag = false;
					}

				}
				else
				{
					g_zone_vector[orig].b_shortest_path_computing_flag = false;
				}

			}
			g_related_zone_vector_size = super_zone_index;
			//clustering
			for (int orig = 0; orig < g_zone_vector.size(); orig++)  // o
			{
				if (g_zone_vector[orig].superzone_index == -1)
				{
					double min_distance = 99999;

					for (int d = 0; d < g_zone_vector.size(); d++)  // o
					{
						if (g_zone_vector[d].bcluster_seed == true)
						{
							double delta_x = g_zone_vector[orig].cell_x - g_zone_vector[d].cell_x;
							double delta_y = g_zone_vector[orig].cell_y - g_zone_vector[d].cell_y;

							double local_distance = pow(delta_x * delta_x + delta_y * delta_y, 0.5);

							if (local_distance < min_distance)
							{
								min_distance = local_distance;
								g_zone_vector[orig].superzone_index = g_zone_vector[d].superzone_index;
								assignment.zone_id_to_seed_zone_id_mapping[g_zone_vector[orig].zone_id] = g_zone_vector[d].zone_id;

							}

						}
					}

				}

			}

			// assign the value back
			for (int orig = 0; orig < g_zone_vector.size(); orig++)  // o
			{
				if (g_zone_vector[orig].sindex >= 0) // related zones
				{
					g_zone_vector[orig].sindex = g_zone_vector[orig].superzone_index;  //rewrite using super zones
				}

			}


		}
	}




	
	double non_related_zone_ratio = non_impact_zone_count * 1.0 / g_zone_vector.size();
	dtalog.output() <<"[DATA INFO] "<< non_impact_zone_count << " zones are not significantly related to the subarea, accounting for " << non_related_zone_ratio * 100 << " % of the total." << '\n';
	g_DTA_log_file <<"[DATA INFO] "<< non_impact_zone_count << " zones are not significantly related to the subarea, accounting for " << non_related_zone_ratio * 100 << " % of the total." << '\n';
	dtalog.output() << "[DATA INFO] "<<inside_significant_zone_count << " zones are inside the subarea. " << '\n';
	g_DTA_log_file << "[DATA INFO] "<<inside_significant_zone_count << " zones are inside the subarea. " << '\n';



	double related_ratio = total_related_demand / max(0.001, total_demand);
	dtalog.output() << "[DATA INFO] Total demand = " << total_demand << ". Subarea-related demand = " << total_related_demand << " ("
		<< related_ratio * 100 << "%), accounting for " << (related_ratio) * 100 << "% inside the subarea, i.e. " << (1-related_ratio) * 100 << " % outside the subarea." << '\n';

	g_DTA_log_file << "[DATA INFO] Total demand = " << total_demand << ". Subarea-related demand = " << total_related_demand << " ("
		<< related_ratio * 100 << "%), accounting for " << (related_ratio) * 100 << "% inside the subarea, i.e. " << (1 - related_ratio) * 100 << " % outside the subarea." << '\n';

		//assignment.summary_file << ",# of subarea related zones =," << g_related_zone_vector_size << '\n';
	//assignment.summary_file << "," <<" # of inside zones = " << inside_significant_zone_count <<  '\n';
	//assignment.summary_file << "," << " # of related outside zones = " << related_external_zone_count << '\n';
}
int g_detector_route_file()
{
	CDTACSVParser parser;
	if (parser.OpenCSVFile("route.csv", false))
	{
		parser.CloseCSVFile();
		return 1;
	}

	return 0;
}

void g_ReadDemandFileBasedOnDemandFileList(Assignment& assignment)
{
//	g_load_demand_side_scenario_file(assignment);

	g_related_zone_vector_size = g_zone_vector.size();
	for (int orig = 0; orig < g_zone_vector.size(); orig++)  // o
	{
		g_zone_vector[orig].sindex = orig;  //copy this index
	}

	// subarea handling step 1: reading

	dtalog.output() << "[PROCESS INFO] Step 1.7: reading input subarea.csv" << '\n';
	g_DTA_log_file << "[PROCESS INFO] Step 1.7: reading input subarea.csv" << '\n';
	g_read_subarea_CSV_file(assignment);
	// subarea handling step 2:
	// for each OD

	if (assignment.g_subarea_shape_points.size() >= 3) // if there is a subarea defined and not using preloaded route file
	{

		for (int orig = 0; orig < g_zone_vector.size(); orig++)  // o
		{
			g_zone_vector[orig].subarea_inside_flag = 0;  //reset the subarea inside flag to 0 once there is a subarea.csv
		}


		// for node layer
		for (int i = 0; i < g_node_vector.size(); i++)
		{
			DTAGDPoint pt;
			pt.x = g_node_vector[i].x;
			pt.y = g_node_vector[i].y;

			if (g_test_point_in_polygon(pt, assignment.g_subarea_shape_points) == 1)
			{
				g_node_vector[i].subarea_id = 1;

			}
		}

		// for link layer
		for (int l = 0; l < g_link_vector.size(); l++)
		{
			if (g_node_vector[g_link_vector[l].from_node_seq_no].subarea_id >= 1 || g_node_vector[g_link_vector[l].to_node_seq_no].subarea_id >= 1)
			{
				g_link_vector[l].subarea_id = g_node_vector[g_link_vector[l].from_node_seq_no].subarea_id;
			}
		}

		// create subarea_related_zone information


			g_create_subarea_related_zone_structure(assignment.g_max_num_significant_zones_in_subarea, assignment.g_max_num_significant_zones_outside_subarea);

	}

	//	fprintf(g_pFileOutputLog, "number of zones =,%lu\n", g_zone_vector.size());
	g_read_departure_time_profile(assignment);

	assignment.InitializeDemandMatrix(g_related_zone_vector_size, g_zone_vector.size(), assignment.g_ModeTypeVector.size(), assignment.g_DemandPeriodVector.size());

	float related_flow_matrix[4][4] = { 0 };

	float total_demand_in_demand_file = 0;

	CDTACSVParser parser;
	
	dtalog.output() << "[PROCESS INFO] Step 2.1: Reading file demand_file_list.csv..." << '\n';
	g_DTA_log_file << "[PROCESS INFO] Step 2.1: Reading file demand_file_list.csv..." << '\n';

	assignment.summary_file << "[PROCESS INFO] Step 2.1: read demand, defined in demand_file_list.csv." << '\n';
	int scenario_index_vector_error_count = 0; 
	int reading_demand_file_log_count = 0;

	int count = 0; 

	if (parser.OpenCSVFile("demand_file_list.csv", false))
	{
		while (parser.ReadRecord())
		{
			int this_departure_time_profile_no = 0;

			int file_sequence_no = 1;

			string format_type = "null";

			int demand_format_flag = 0;

			if (!parser.GetValueByFieldName("file_sequence_no", file_sequence_no))
				break;

			// skip negative sequence no
			if (file_sequence_no <= -1)
				continue;

			double loading_scale_factor = 1.0;
			string file_name, demand_period_str, mode_type;
			parser.GetValueByFieldName("file_name", file_name);
			parser.GetValueByFieldName("demand_period", demand_period_str);
			parser.GetValueByFieldName("format_type", format_type);
			parser.GetValueByFieldName("scale_factor", loading_scale_factor, false);
			parser.GetValueByFieldName("departure_time_profile_no", this_departure_time_profile_no, false);

			string scenario_index_vector_str;
			if (parser.GetValueByFieldName("scenario_index_vector", scenario_index_vector_str, false, false) == false)
			{
				if(count ==0 )
				dtalog.output() << "[WARNING] Field scenario_index_vector is missing in file demand_file_list.csv." << '\n';
				g_DTA_log_file << "[WARNING] Field scenario_index_vector is missing in file demand_file_list.csv." << '\n';

			}
			std::vector<int> scenario_index_vector;

			if (scenario_index_vector_str.size() == 0)  // default for scenario 0
				scenario_index_vector_str = "0";

			g_ParserIntSequence(scenario_index_vector_str, scenario_index_vector);

			for (int scenario_index_i = 0; scenario_index_i < scenario_index_vector.size(); scenario_index_i++)
			{

				if(reading_demand_file_log_count <5)
				{
				dtalog.output() << "[STATUS INFO] reading demand file " << file_name.c_str () << " for scenario index = " << scenario_index_i << '\n';
				g_DTA_log_file << "[STATUS INFO] reading demand file " << file_name.c_str () << " for scenario index = " << scenario_index_i << '\n';
				reading_demand_file_log_count++; 
				}


				assignment.summary_file << "[STATUS INFO] reading demand file for scenario index = " << scenario_index_i << '\n';
				if (loading_scale_factor < 0.0001)
					continue;

				int si = scenario_index_vector[scenario_index_i];
				if (assignment.g_active_DTAscenario_map.find(si) == assignment.g_active_DTAscenario_map.end())
				{
					if(scenario_index_vector_error_count <3)
					{
					dtalog.output() << "[WARNING] scenario_index = " << si << " in  the field of scenario_index_vector in file demand_file_list.csv  has not been defined in file scenario_file_list.csv." << '\n';
					g_DTA_log_file << "[WARNING] scenario_index = " << si << " in  the field of scenario_index_vector in file demand_file_list.csv  has not been defined in file scenario_file_list.csv." << '\n';
					scenario_index_vector_error_count++; 
					}
					continue;
				}


				if (this_departure_time_profile_no >= assignment.g_DepartureTimeProfileVector.size())
				{
					dtalog.output() << "[ERROR] departure_time_profile_no = " << this_departure_time_profile_no << " in  file demand_file_list.csv  has not been defined in file departure_time_profile.csv." << '\n';
					g_DTA_log_file << "[ERROR] departure_time_profile_no = " << this_departure_time_profile_no << " in  file demand_file_list.csv  has not been defined in file departure_time_profile.csv." << '\n';
					this_departure_time_profile_no = 0;

				}

				if (parser.GetValueByFieldName("mode_type", mode_type, false, false) == false)
					mode_type = assignment.g_ModeTypeVector[0].mode_type; 

				int mode_type_no = 0;
				int demand_period_no = 0;

				std::string str  = demand_period_str;

				std::transform(str.begin(), str.end(), str.begin(),
					[](unsigned char c) { return std::tolower(c); });

				demand_period_str = str;

				if (assignment.demand_period_to_seqno_mapping.find(demand_period_str) != assignment.demand_period_to_seqno_mapping.end())
					demand_period_no = assignment.demand_period_to_seqno_mapping[demand_period_str];
				else
				{
					if (demand_period_str.size() == 0)
					{
						dtalog.output() << "[WARNING] file_sequence_no =" << file_sequence_no << " demand period = is empty in demand_file_list.csv.  " << '\n';
						g_DTA_log_file << "[WARNING] file_sequence_no =" << file_sequence_no << " demand period = is empty in demand_file_list.csv.  " << '\n';
						demand_period_no = 0; 
					}
					else
					{
						dtalog.output() << "[ERROR]  demand_period= " << demand_period_str.c_str() << " in  demand_file_list.csv has not been defined in demand_period.csv." << '\n';
						g_DTA_log_file << "[ERROR]  demand_period= " << demand_period_str.c_str() << " in  demand_file_list.csv has not been defined in demand_period.csv." << '\n';
						demand_period_no = 0;
					}

				}

				//char time_interval_field_name[20];
				CDemand_Period  demand_period = assignment.g_DemandPeriodVector[demand_period_no];
				assignment.g_DemandPeriodVector[demand_period_no].number_of_demand_files++;


				if (demand_period.starting_time_slot_no * MIN_PER_TIMESLOT < assignment.g_LoadingStartTimeInMin)
					assignment.g_LoadingStartTimeInMin = demand_period.starting_time_slot_no * MIN_PER_TIMESLOT;

				if (demand_period.ending_time_slot_no * MIN_PER_TIMESLOT > assignment.g_LoadingEndTimeInMin)
					assignment.g_LoadingEndTimeInMin = demand_period.ending_time_slot_no * MIN_PER_TIMESLOT;

				if (assignment.g_LoadingEndTimeInMin < assignment.g_LoadingStartTimeInMin)
				{
					assignment.g_LoadingEndTimeInMin = assignment.g_LoadingStartTimeInMin + 1; // in case user errror
				}

				if (format_type.find("null") != string::npos)  // skip negative sequence no
				{
					dtalog.output() << "[ERROR] Please provide format_type in file demand_file_list." << '\n';
					g_DTA_log_file << "[ERROR] Please provide format_type in file demand_file_list." << '\n';
					format_type = "cvs";
				}



					if (assignment.mode_type_2_seqno_mapping.find(mode_type) != assignment.mode_type_2_seqno_mapping.end())
						mode_type_no = assignment.mode_type_2_seqno_mapping[mode_type];
					else
					{
						dtalog.output() << "[ERROR] mode_type = " << mode_type.c_str() << " in field mode_type of file demand_file_list.csv is not defined in the file mode_type.csv yet." << '\n';
						g_DTA_log_file << "[ERROR] mode_type = " << mode_type.c_str() << " in field mode_type of file demand_file_list.csv is not defined in the file mode_type.csv yet." << '\n';
						mode_type = "auto";
					}

				if (demand_period_no > MAX_TIMEPERIODS)
				{
					dtalog.output() << "[ERROR] demand_period_no should be less than settings in demand_period section. Please change the parameter settings in the source code." << '\n';
					g_DTA_log_file << "[ERROR] demand_period_no should be less than settings in demand_period section. Please change the parameter settings in the source code." << '\n';
					demand_period_no = 0;
				}



				if (format_type.find("column") != string::npos || format_type.find("bin") != string::npos)  // or muliti-column
				{

					// try to detect if we have a route.csv as preload file, if yes, we skip the following reading .
					CDTACSVParser parser_route;

					struct SDemandHeader
					{
						int o_zone_id, d_zone_id, mode_type_no, demand_period;
						double volume;
					};

					SDemandHeader header;


					CDTACSVParser parser;
					FILE* pFile;
					int file_exists = 0;
					int line_no = 0;
					int file_format = 0;

					if (parser_route.OpenCSVFile("route.csv", false))
					{
						// a route file exists
						dtalog.output() << "[STATUS INFO] route.csv exists as preload file so we skip the reading for the column based demand file." << '\n';
						g_DTA_log_file << "[STATUS INFO] route.csv exists as preload file so we skip the reading for the column based demand file." << '\n';
					}
					else
					{
						fopen_ss(&pFile, "demand.bin", "rb");
						if (pFile != NULL)
						{
							file_format = 2;
							dtalog.output() << "[STATUS INFO] reading demand.bin in fast binary file reading mode." << '\n';
							g_DTA_log_file << "[STATUS INFO] reading demand.bin in fast binary file reading mode." << '\n';
						}


						if (file_format == 0 && format_type.find("column") != string::npos)
						{
							if (parser.OpenCSVFile(file_name, false) == false)
							{
								dtalog.output() << "[ERROR] column file " << file_name.c_str() << ".csv does not exist." << '\n';
								g_DTA_log_file << "[ERROR] column file " << file_name.c_str() << ".csv does not exist." << '\n';
								assignment.summary_file << "[ERROR] File Missing!,file_sequence_no=," << file_sequence_no << ",file_name =, " << file_name.c_str() << '\n';
								continue; 

							}
							file_format = 1;
						}
					}

					int error_count = 0;
					int critical_OD_count = 0;
					double critical_OD_volume = 0;
					// read the file formaly after the test.
					std::map<int, int> missing_zone_map; 
					while (file_format >= 1)
					{
						int o_zone_id, d_zone_id;
						double demand_value = 0;

						if (file_format == 1)  // column
						{
							int flag = parser.ReadRecord();
							if (flag == 0)  // this is end of file.
								break;

							parser.GetValueByFieldName("o_zone_id", o_zone_id);
							parser.GetValueByFieldName("d_zone_id", d_zone_id);
							parser.GetValueByFieldName("volume", demand_value);

						}
						if (file_format == 2)  // binary
						{
							if (feof(pFile))
								break;

							size_t result; fread(&header, sizeof(header), 1, pFile);
							o_zone_id = header.o_zone_id;
							d_zone_id = header.d_zone_id;
							demand_value = header.volume;

							mode_type_no = header.mode_type_no;
							demand_period_no = header.demand_period;

						}

						// end of reading read each record


						if (o_zone_id == d_zone_id)
						{
							continue;
						}

						if (assignment.g_zoneid_to_zone_seq_no_mapping.find(o_zone_id) == assignment.g_zoneid_to_zone_seq_no_mapping.end())
						{
							if (error_count < 10 && missing_zone_map.find(o_zone_id)== missing_zone_map.end())
							{
								dtalog.output() << '\n' << "[WARNING] origin zone " << o_zone_id << "  has no activity nodes defined in node.csv or zone.csv." << '\n';
								g_DTA_log_file << '\n' << "[WARNING] origin zone " << o_zone_id << "  has no activity nodes defined in node.csv or zone.csv." << '\n';
								missing_zone_map[o_zone_id] = 1;
							} 

							error_count++;
							// origin zone has not been defined, skipped.
							continue;
						}

						if (assignment.g_zoneid_to_zone_seq_no_mapping.find(d_zone_id) == assignment.g_zoneid_to_zone_seq_no_mapping.end())
						{
							if (error_count < 10 && missing_zone_map.find(d_zone_id) == missing_zone_map.end())
							{
								dtalog.output() << '\n' << "[WARNING] destination zone " << d_zone_id << "  has no activity nodes defined in node.csv or zone.csv." << '\n';
								g_DTA_log_file << '\n' << "[WARNING] destination zone " << d_zone_id << "  has no activity nodes defined in node.csv or zone.csv." << '\n';
								missing_zone_map[d_zone_id] = 1;
							}

							error_count++;
							// destination zone has not been defined, skipped.
							continue;
						}

						int from_zone_seq_no = 0;
						int to_zone_seq_no = 0;
						from_zone_seq_no = assignment.g_zoneid_to_zone_seq_no_mapping[o_zone_id];
						to_zone_seq_no = assignment.g_zoneid_to_zone_seq_no_mapping[d_zone_id];

						// test if bypass subaera
						int inside_flag_from = g_zone_vector[from_zone_seq_no].subarea_inside_flag;
						int inside_flag_to = g_zone_vector[to_zone_seq_no].subarea_inside_flag;



						if (o_zone_id == 1394 && d_zone_id == 39)
						{
							int idebug = 1;
						}

						if (inside_flag_from == 2)
						{
							if (inside_flag_to == 2)
							{
								int idebug = 1;
							}

						}


						if (inside_flag_from == 0 || inside_flag_to == 0)  // one of trip ends are not related the subarea
						{
							if (g_zone_vector[from_zone_seq_no].impact_passing_ODdemand.find(to_zone_seq_no) == g_zone_vector[from_zone_seq_no].impact_passing_ODdemand.end())
								continue;  // skip all non passing OD pairs
						}

						//I = I: 2-2
						//I - E1: 2-1
						//I - E2: 2-0

						//E1 - I: 1-2
						//E2 - I: 0-2

						///E1-E1: 1-1

						//E2-E2
						//E1-E2
						//E1-E1

						related_flow_matrix[inside_flag_from][inside_flag_to] += demand_value;


						//if ((inside_flag_from == 2 || inside_flag_from == 3) && (inside_flag_to == 2 || inside_flag_to == 3))
						//{
						//	total_matrix_demand_value += demand_value;
						//	assignment.summary_file << "," << o_zone_id << "," << d_zone_id << "," << inside_flag_from << "," << inside_flag_to << "," << demand_value << ",cd= " << total_matrix_demand_value << '\n';
						//} // print out all OD pairs in the matrix


						int from_zone_sindex = g_zone_vector[from_zone_seq_no].sindex;
						if (from_zone_sindex == -1)
							continue;

						int to_zone_sindex = g_zone_vector[to_zone_seq_no].sindex;
						if (to_zone_sindex == -1)
							continue;


						double local_scale_factor = 1.0;

						//redundancy analysis
						//o based factor

						int o_district = assignment.g_zone_seq_no_to_analysis_distrct_id_mapping[from_zone_seq_no];
						int d_district = assignment.g_zone_seq_no_to_analysis_distrct_id_mapping[to_zone_seq_no];


						if (assignment.o_district_id_factor_map.find(o_district) != assignment.o_district_id_factor_map.end())
						{
							local_scale_factor = assignment.o_district_id_factor_map[o_district];
						}

						//d based factor
						if (assignment.d_district_id_factor_map.find(d_district) != assignment.d_district_id_factor_map.end())
						{
							local_scale_factor = assignment.d_district_id_factor_map[d_district];
						}

						int od_district_key = o_district * 1000 + d_district;
						if (assignment.od_district_id_factor_map.find(od_district_key) != assignment.od_district_id_factor_map.end())
						{
							local_scale_factor = assignment.od_district_id_factor_map[od_district_key];
						}

						if (assignment.shortest_path_log_zone_id == -1)  // set this to the first zone
							assignment.shortest_path_log_zone_id = o_zone_id;

						// encounter return
						if (demand_value < -99)
							break;

						demand_value *= (loading_scale_factor * local_scale_factor);


						if (demand_value >= 5)
						{
							critical_OD_volume += demand_value;
							critical_OD_count += 1;
							//dtalog.output() << origin_zone << "," << destination_zone << "," << demand_value << "," << "\"LINESTRING( " <<
							//g_DTA_log_file << origin_zone << "," << destination_zone << "," << demand_value << "," << "\"LINESTRING( " <<
							//    assignment.zone_id_X_mapping[origin_zone] << " " << assignment.zone_id_Y_mapping[origin_zone] << "," <<
							//    assignment.zone_id_X_mapping[destination_zone] << " " << assignment.zone_id_Y_mapping[destination_zone] << ")\" " << '\n';

						}

						g_scenario_summary_map[si].record_mode_volume(demand_period_no,mode_type_no, demand_value);

						assignment.total_demand[mode_type_no][demand_period_no] += demand_value;
						assignment.g_column_pool[from_zone_sindex][to_zone_sindex][mode_type_no][demand_period_no].od_volume[si] += demand_value;
						assignment.g_column_pool[from_zone_sindex][to_zone_sindex][mode_type_no][demand_period_no].departure_time_profile_no = this_departure_time_profile_no;
						assignment.total_demand_volume[si] += demand_value;

						if(assignment.g_ModeTypeVector[mode_type_no].real_time_information_type>=1)
						{
							assignment.total_real_time_demand_volume += demand_value;
	
						}

						assignment.g_origin_demand_array[from_zone_seq_no] += demand_value;



						// we generate vehicles here for each OD data line
						if (line_no <= 1 && demand_value >= 0.001)
						{// read only one line, but has not reached the end of the line
							dtalog.output() << "[DATA INFO] o_zone_id:" << o_zone_id << ", d_zone_id: " << d_zone_id << ", value = " << demand_value << '\n';
							g_DTA_log_file << "[DATA INFO] o_zone_id:" << o_zone_id << ", d_zone_id: " << d_zone_id << ", value = " << demand_value << '\n';
						}
						line_no++;

					}

					dtalog.output() << "[DATA INFO] reading file " << file_name.c_str() << ", scenario no : " << si << ", cumulative total demand volume is " << assignment.total_demand_volume[si] << '\n';
					g_DTA_log_file << "[DATA INFO] reading file " << file_name.c_str() << ", scenario no: " << si << ", cumulative total demand volume is " << assignment.total_demand_volume[si] << '\n';
					//dtalog.output() << "crtical demand volume has " << critical_OD_count << " OD pairs in size," << critical_OD_volume << ", " << ", account for " << critical_OD_volume / max(0.1f, assignment.total_demand_volume) * 100 << "%%" << '\n';
					//g_DTA_log_file << "crtical demand volume has " << critical_OD_count << " OD pairs in size," << critical_OD_volume << ", " << ", account for " << critical_OD_volume / max(0.1f, assignment.total_demand_volume) * 100 << "%%" << '\n';

					//dtalog.output() << "crtical OD zones volume has " << critical_OD_count << " OD pairs in size," << critical_OD_volume << ", " << ", account for " << critical_OD_volume / max(0.1f, assignment.total_demand_volume) * 100 << "%%" << '\n';
					//g_DTA_log_file << "crtical OD zones volume has " << critical_OD_count << " OD pairs in size," << critical_OD_volume << ", " << ", account for " << critical_OD_volume / max(0.1f, assignment.total_demand_volume) * 100 << "%%" << '\n';

				}
				else if (format_type.compare("route") == 0)
				{

					int path_counts = 0;
					float sum_of_path_volume = 0;
					CDTACSVParser parser;
					if (parser.OpenCSVFile(file_name, false))
					{
						int total_path_in_demand_file = 0;
						// read agent file line by line,

						int o_zone_id, d_zone_id;
						string mode_type, demand_period;

						std::vector <int> node_sequence;

						while (parser.ReadRecord())
						{
							total_path_in_demand_file++;
							if (total_path_in_demand_file % 1000 == 0)
								dtalog.output() << "[DATA INFO] total_path_in_demand_file is " << total_path_in_demand_file << '\n';
								g_DTA_log_file << "[DATA INFO] total_path_in_demand_file is " << total_path_in_demand_file << '\n';

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
							volume *= loading_scale_factor;
							agent_path_element.volume = volume;
							path_counts++;
							sum_of_path_volume += agent_path_element.volume;

							assignment.total_demand[mode_type_no][demand_period_no] += agent_path_element.volume;
							assignment.g_column_pool[from_zone_sindex][to_zone_sindex][mode_type_no][demand_period_no].od_volume[si] += agent_path_element.volume;
							assignment.g_column_pool[from_zone_sindex][to_zone_sindex][mode_type_no][demand_period_no].departure_time_profile_no = this_departure_time_profile_no;
							assignment.total_demand_volume[si] += agent_path_element.volume;
							assignment.g_origin_demand_array[from_zone_seq_no] += agent_path_element.volume;

							g_scenario_summary_map[si].record_mode_volume(demand_period_no,mode_type_no, agent_path_element.volume);
							//apply for both agent csv and routing policy
							assignment.g_column_pool[from_zone_seq_no][to_zone_seq_no][mode_type_no][demand_period_no].bfixed_route = true;

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
						}
						dtalog.output() << "[DATA INFO] total_demand_volume loaded from path file is " << sum_of_path_volume << " with " << path_counts << "paths." << '\n';
						g_DTA_log_file << "[DATA INFO] total_demand_volume loaded from path file is " << sum_of_path_volume << " with " << path_counts << "paths." << '\n';

					}
					else
					{
						//open file
						dtalog.output() << "[ERROR] File " << file_name << " cannot be opened.\n It might be currently used and locked by EXCEL." << '\n';
						g_DTA_log_file << "[ERROR] File " << file_name << " cannot be opened.\n It might be currently used and locked by EXCEL." << '\n';
						continue; 
					}
				}
				else if (format_type.compare("activity_chain") == 0)
				{
					////////////////////////
	// Step 1.3: Reading choice_set.csv...
					dtalog.output() << "[PROCESS INFO] Step 1.3: Reading choice_set.csv..." << '\n';
					g_DTA_log_file << "[PROCESS INFO] Step 1.3: Reading choice_set.csv..." << '\n';
					assignment.summary_file << "[PROCESS INFO] Step 1.3: Reading choice_set.csv..." << '\n';

					CDTACSVParser parser_choice_set;
					int Field_volume_log_count = 0; 
					if (parser_choice_set.OpenCSVFile("choice_set.csv", false))
					{
						while (parser_choice_set.ReadRecord())
						{
							int index = 0;
							std::string multi_id;
							int alt_id;
							std::string pattern_id;


							CChoiceAlt choice_alt;

							if (!parser_choice_set.GetValueByFieldName("choice_set_index", index))
								break;

							if (!parser_choice_set.GetValueByFieldName("multi_dim_choice_id", multi_id))
							{
								dtalog.output() << "[ERROR] Field multi_dim_choice_id in file choice_set.csv cannot be read." << '\n';
								g_DTA_log_file << "[ERROR] Field multi_dim_choice_id in file choice_set.csv cannot be read." << '\n';
								continue; 
							}

							std::string mode_tag, demand_period_tag, spatial_tag, travel_purpose_tag, data_tag;

							parser_choice_set.GetValueByFieldName("mode_tag", mode_tag);
							parser_choice_set.GetValueByFieldName("demand_period_tag", demand_period_tag);
							parser_choice_set.GetValueByFieldName("spatial_tag", spatial_tag);
							parser_choice_set.GetValueByFieldName("travel_purpose_tag", travel_purpose_tag);
							parser_choice_set.GetValueByFieldName("data_tag", data_tag);

							assignment.g_ChoiceSetMap[multi_id].active_scenario_index = si;
							assignment.g_ChoiceSetMap[multi_id].multi_dim_choice_id = multi_id;
							assignment.g_ChoiceSetMap[multi_id].choice_set_index = index;

							if (mode_tag.size() > 0)
								assignment.g_ChoiceSetMap[multi_id].mode_tag = mode_tag;

							if (demand_period_tag.size() > 0)
								assignment.g_ChoiceSetMap[multi_id].demand_period_tag = demand_period_tag;

							if (spatial_tag.size() > 0)
								assignment.g_ChoiceSetMap[multi_id].spatial_tag = spatial_tag;

							if (travel_purpose_tag.size() > 0)
								assignment.g_ChoiceSetMap[multi_id].travel_purpose_tag = travel_purpose_tag;

							if (data_tag.size() > 0)
								assignment.g_ChoiceSetMap[multi_id].data_tag = data_tag;


							if (!parser_choice_set.GetValueByFieldName("choice_alternative_id", choice_alt.choice_alternative_id))
							{
								dtalog.output() << "[ERROR] Field choice_alternative_id in file choice_set.csv cannot be read." << '\n';
								g_DTA_log_file << "[ERROR] Field choice_alternative_id in file choice_set.csv cannot be read." << '\n';
								continue; 
							}

							std::string activity_travel_pattern_id;
							if (!parser_choice_set.GetValueByFieldName("activity_travel_pattern_id", activity_travel_pattern_id))
							{
								dtalog.output() << "[ERROR] Field activity_travel_pattern_id in file choice_set.csv cannot be read." << '\n';
								g_DTA_log_file << "[ERROR] Field activity_travel_pattern_id in file choice_set.csv cannot be read." << '\n';
								continue;
							}

							std::string activity_zone_vector_str;
							if (!parser_choice_set.GetValueByFieldName("zone_chain", activity_zone_vector_str))
							{
								dtalog.output() << "[ERROR] Field zone_chain in file activity_travel_pattern.csv cannot be read." << '\n';
								g_DTA_log_file << "[ERROR] Field zone_chain in file activity_travel_pattern.csv cannot be read." << '\n';
								continue;
							}

							choice_alt.activity_zone_chain_str = activity_zone_vector_str;

							int number_of_odpairs = g_ParserIntSequence(activity_zone_vector_str, choice_alt.activity_zone_chain);

							if (assignment.g_ActivityTravelPatternMap.find(activity_travel_pattern_id) != assignment.g_ActivityTravelPatternMap.end())
							{
								if (number_of_odpairs != assignment.g_ActivityTravelPatternMap[activity_travel_pattern_id].number_of_trips + 1)
								{
									dtalog.output() << "[ERROR] number_of_trips != assignment.g_ActivityTravelPatternMap[activity_travel_pattern_id].number_of_trips in file activity_travel_pattern.csv and choice_set.csv." << '\n';
									g_DTA_log_file << "[ERROR] number_of_trips != assignment.g_ActivityTravelPatternMap[activity_travel_pattern_id].number_of_trips in file activity_travel_pattern.csv and choice_set.csv." << '\n';
									continue; 
								}
							}


							float volume = 0;
							if (!parser_choice_set.GetValueByFieldName("volume", volume) && Field_volume_log_count<=2)
							{
								dtalog.output() << "[ERROR] Field volume in file activity_travel_pattern.csv cannot be read." << '\n';
								g_DTA_log_file << "[ERROR] Field volume in file activity_travel_pattern.csv cannot be read." << '\n';
								Field_volume_log_count++; 
								continue; 
							}

							choice_alt.multi_dim_choice_id = multi_id;
							choice_alt.volume = volume;
							choice_alt.activity_travel_pattern_id = activity_travel_pattern_id;

							int mode_type_no;
							int demand_period_no;
							int from_zone_sindex;
							int to_zone_sindex;

							if (assignment.g_ActivityTravelPatternMap.find(activity_travel_pattern_id) != assignment.g_ActivityTravelPatternMap.end())
							{
								// for each lag
								for (int trip_no = 0; trip_no < assignment.g_ActivityTravelPatternMap[activity_travel_pattern_id].number_of_trips; trip_no++)
								{
									mode_type_no = assignment.g_ActivityTravelPatternMap[activity_travel_pattern_id].mode_vector[trip_no];
									demand_period_no = assignment.g_ActivityTravelPatternMap[activity_travel_pattern_id].demand_period_vector[trip_no];
									int o_zone_id = choice_alt.activity_zone_chain[trip_no];
									int d_zone_id = choice_alt.activity_zone_chain[trip_no + 1];
									int from_zone_seq_no = 0;
									int to_zone_seq_no = 0;
									from_zone_seq_no = assignment.g_zoneid_to_zone_seq_no_mapping[o_zone_id];
									to_zone_seq_no = assignment.g_zoneid_to_zone_seq_no_mapping[d_zone_id];

									// test if bypass subaera
									int inside_flag_from = g_zone_vector[from_zone_seq_no].subarea_inside_flag;
									int inside_flag_to = g_zone_vector[to_zone_seq_no].subarea_inside_flag;



									if (o_zone_id == 1394 && d_zone_id == 39)
									{
										int idebug = 1;
									}

									if (inside_flag_from == 2)
									{
										if (inside_flag_to == 2)
										{
											int idebug = 1;
										}

									}


									if (inside_flag_from == 0 || inside_flag_to == 0)  // one of trip ends are not related the subarea
									{
										if (g_zone_vector[from_zone_seq_no].impact_passing_ODdemand.find(to_zone_seq_no) == g_zone_vector[from_zone_seq_no].impact_passing_ODdemand.end())
											continue;  // skip all non passing OD pairs
									}

									SChoiceAlt s_alt;
									s_alt.mode_no = mode_type_no;
									s_alt.demand_peroid_no = demand_period_no;
									s_alt.o_zone_no = from_zone_seq_no;
									s_alt.d_zone_no = to_zone_seq_no;

									choice_alt.activity_chain.push_back(s_alt);

									assignment.total_demand[mode_type_no][demand_period_no] += volume;
									assignment.g_column_pool[from_zone_seq_no][to_zone_seq_no][mode_type_no][demand_period_no].od_volume[si] += volume;
									assignment.total_demand_volume[si] += volume;
									assignment.g_origin_demand_array[from_zone_seq_no] += volume;
									g_scenario_summary_map[si].record_mode_volume(demand_period_no,mode_type_no, volume);
								}
								//end of trip leg loop

								assignment.g_ChoiceSetMap[multi_id].choice_alt_vector.push_back(choice_alt);
								assignment.g_ChoiceSetMap[multi_id].choice_set_index = index;



								if (assignment.g_ChoiceSetMap.size() == 0)
								{
									dtalog.output() << "[ERROR] file choice_set has no information." << '\n';
									g_DTA_log_file << "[ERROR] file choice_set has no information." << '\n';
									continue; 
								}
							}
						}
						parser_choice_set.CloseCSVFile();
					}
					else
					{
						dtalog.output() << "[ERROR] File choice_set.csv cannot be opened.\n It might be currently used and locked by EXCEL." << '\n';
						g_DTA_log_file << "[ERROR] File choice_set.csv cannot be opened.\n It might be currently used and locked by EXCEL." << '\n';
						continue; 
					}


				}
			}
		}
	}
	/////
	/// summary
	//////


	/// <summary>
	///  subarea
	/// </summary>
	/// <param name="assignment"></param>
	/*assignment.summary_file << ",total demand =, " << assignment.total_demand_volume[si] << '\n';*/

	g_check_demand_volume_with_mode_type(assignment);

	std::vector<CODState> ODStateVector;
	for (int orig = 0; orig < g_zone_vector.size(); orig++)  // o
	{
		CColumnVector* p_column_pool;
		int path_seq_count = 0;

		int from_zone_sindex = g_zone_vector[orig].sindex;
		if (from_zone_sindex == -1)
			continue;

		for (int dest = 0; dest < g_zone_vector.size(); dest++) //d
		{
			int to_zone_sindex = g_zone_vector[dest].sindex;
			if (to_zone_sindex == -1)
				continue;

			int inside_flag_from = g_zone_vector[orig].subarea_inside_flag;
			int inside_flag_to = g_zone_vector[dest].subarea_inside_flag;

			for (int at = 0; at < assignment.g_ModeTypeVector.size(); at++)  //m
			{
				for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); tau++)  //tau
				{
					p_column_pool = &(assignment.g_column_pool[from_zone_sindex][to_zone_sindex][at][tau]);
					if (p_column_pool->od_volume[assignment.active_scenario_index] > 0)
					{
						CODState ods;
						ods.setup_input(orig, dest, at, tau, inside_flag_from, inside_flag_to);

						float passing_demand = p_column_pool->od_volume[assignment.active_scenario_index];
						if (g_zone_vector[orig].impact_passing_ODdemand.find(dest) != g_zone_vector[orig].impact_passing_ODdemand.end())
						{
							passing_demand = g_zone_vector[orig].impact_passing_ODdemand[dest];
						}

						ods.input_value(passing_demand);
						ODStateVector.push_back(ods);
					}
				}
			}
		}
	}

	assignment.summary_file << ",from_flag,to_flag,volume" << '\n';
	assignment.summary_file << ",from_flag,to_flag,volume" << '\n';
	assignment.summary_file << "FOCUSING internal step 1: Focus-subarea Approach,0=not related, 1=significantly related external but not in cutoff, 2= cut-off external, 3= inside" << '\n';
	assignment.summary_file << ",";
	for (int subare_inside_flag_to = 0; subare_inside_flag_to <= 3; subare_inside_flag_to++)
		assignment.summary_file << subare_inside_flag_to << ",";

	assignment.summary_file << '\n';

	for (int subare_inside_flag_from = 0; subare_inside_flag_from <= 3; subare_inside_flag_from++)
	{
		assignment.summary_file << subare_inside_flag_from;
		for (int subare_inside_flag_to = 0; subare_inside_flag_to <= 3; subare_inside_flag_to++)
		{
			assignment.summary_file << "," << related_flow_matrix[subare_inside_flag_from][subare_inside_flag_to];
		}
		assignment.summary_file << '\n';
	}



	if(ODStateVector.size()>0)
	{ 

		std::sort(ODStateVector.begin(), ODStateVector.end());

	assignment.summary_file << "FOCUSING internal step 2: Origin-based flow extraction,top 10 OD,rank,o,d,inside_flag_o,inside_flag_d,mode_type,departure_time,volume" << '\n';

	for (int k = 0; k < min(size_t(100), ODStateVector.size()); k++)
	{
		int o = ODStateVector[k].orig;
		int d = ODStateVector[k].dest;
		int at = ODStateVector[k].at;
		int tau = ODStateVector[k].tau;

		assignment.summary_file << ",," << k + 1 << "," << g_zone_vector[o].zone_id << "," << g_zone_vector[d].zone_id << "," <<
			ODStateVector[k].subarea_inside_flag_orig << "," << ODStateVector[k].subarea_inside_flag_dest << "," <<
			assignment.g_ModeTypeVector[at].mode_type.c_str() << "," << assignment.g_DemandPeriodVector[tau].demand_period.c_str()
			<< "," << ODStateVector[k].value << '\n';
	}

	}
}

void g_ReadOutputFileConfiguration(Assignment& assignment)
{
	dtalog.output() << "[PROCESS INFO] Step 1.9: Reading file section [output_file_configuration] in settings.csv..." << '\n';
	g_DTA_log_file << "[PROCESS INFO] Step 1.9: Reading file section [output_file_configuration] in settings.csv..." << '\n';

	CDTACSVParser parser;
	if (parser.OpenCSVFile("settings.csv", false))
	{
		while (parser.ReadRecord())
		{
			parser.GetValueByKeyName("path_output", assignment.path_output, false, false);
			parser.GetValueByKeyName("major_path_volume_threshold", assignment.major_path_volume_threshold, false, false);
			parser.GetValueByKeyName("shortest_path_log_zone_id", assignment.shortest_path_log_zone_id, false, false);
			parser.GetValueByKeyName("trajectory_output_count", assignment.trajectory_output_count, false, false);
			parser.GetValueByKeyName("trace_output", assignment.trace_output, false, false);

			parser.GetValueByKeyName("trajectory_sampling_rate", assignment.trajectory_sampling_rate, false, false);
			parser.GetValueByKeyName("trajectory_diversion_only", assignment.trajectory_diversion_only, false, false);
			parser.GetValueByKeyName("td_link_performance_sampling_interval_in_min", assignment.td_link_performance_sampling_interval_in_min, false, false);


			dtalog.output() << "[DATA INFO] dynamic_link_performance_sampling_interval_in_min = " << assignment.td_link_performance_sampling_interval_in_min << " min" << '\n';
			g_DTA_log_file << "[DATA INFO] dynamic_link_performance_sampling_interval_in_min = " << assignment.td_link_performance_sampling_interval_in_min << " min" << '\n';
			dtalog.output() << "[DATA INFO] dynamic_link_performance_sampling_interval_hd_in_min = " << assignment.dynamic_link_performance_sampling_interval_hd_in_min << " min" << '\n';
			g_DTA_log_file << "[DATA INFO] dynamic_link_performance_sampling_interval_hd_in_min = " << assignment.dynamic_link_performance_sampling_interval_hd_in_min << " min" << '\n';

		}

		parser.CloseCSVFile();
	}
}

void g_ReadInformationConfiguration(Assignment& assignment)
{
	dtalog.output() << ",Reading file section [real_time_info] in settings.csv..." << '\n';
	g_DTA_log_file << ",Reading file section [real_time_info] in settings.csv..." << '\n';

	dtalog.output() << "[STATUS INFO] Reading file section [real_time_info] in settings.csv..." << '\n';
	g_DTA_log_file << "[STATUS INFO] Reading file section [real_time_info] in settings.csv..." << '\n';

	CDTACSVParser parser;
	if (parser.OpenCSVFile("settings.csv", false))
	{

		while (parser.ReadRecord())
		{
			string str_key;
			parser.GetValueByFieldName("key", str_key);

			if (str_key == "info_updating_freq_in_min")
			{
				parser.GetValueByFieldName("value", assignment.g_info_updating_freq_in_min, false, false);
				dtalog.output() << "[DATA INFO] info_updating_freq_in_min = " << assignment.g_info_updating_freq_in_min << " min" << '\n';
				g_DTA_log_file << "[DATA INFO] info_updating_freq_in_min = " << assignment.g_info_updating_freq_in_min << " min" << '\n';

			}

			if (str_key == "info_updating_freq_in_min")
			{
				parser.GetValueByFieldName("value", assignment.g_visual_distance_in_cells, false, false);
				dtalog.output() << "[DATA INFO] visual_distance_in_cells = " << assignment.g_visual_distance_in_cells << " cells" << '\n';
				g_DTA_log_file << "[DATA INFO] visual_distance_in_cells = " << assignment.g_visual_distance_in_cells << " cells" << '\n';
			}

		}

		parser.CloseCSVFile();
	}
	assignment.summary_file << ",info_updating_freq_in_min= " << assignment.g_info_updating_freq_in_min << " min" << '\n';
	assignment.summary_file << ",visual_distance_in_cells= " << assignment.g_visual_distance_in_cells << " cells" << '\n';
}

void g_add_new_virtual_connector_link(int internal_from_node_seq_no, int internal_to_node_seq_no, string mode_type_str, int zone_seq_no = -1)
{
	// create a link object
	CLink link;
	link.link_id = "connector";
	link.link_type_code = "connector";

	link.from_node_seq_no = internal_from_node_seq_no;
	link.to_node_seq_no = internal_to_node_seq_no;
	link.link_seq_no = assignment.g_number_of_links;
	link.to_node_seq_no = internal_to_node_seq_no;
	//virtual connector

	for (int si = 0; si < g_number_of_active_scenarios; si++)
	{
		link.link_type_si[si] = -1;
		link.number_of_lanes_si[si] = 20;  // default all open
	}
	//only for outgoing connectors
	link.zone_seq_no_for_outgoing_connector = zone_seq_no;

	//BPR
	link.traffic_flow_code = 0;

	link.spatial_capacity_in_vehicles = 99999;
	link.lane_capacity = 999999;
	link.link_spatial_capacity = 99999;
	link.link_distance_VDF = 0.00001;
	link.free_flow_travel_time_in_min = 0.1;

	for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
	{
		//setup default values
		link.VDF_period[tau].lane_based_ultimate_hourly_capacity = 99999;
		// 60.0 for 60 min per hour

		link.VDF_period[tau].alpha = 0;
		link.VDF_period[tau].beta = 0;

		for (int at = 0; at < assignment.g_ModeTypeVector.size(); at++)
		{
			link.VDF_period[tau].FFTT_at[at] = 0.0001;
			link.link_avg_travel_time_per_period[tau][at] = 0;
		}




	}

	// add this link to the corresponding node as part of outgoing node/link
	g_node_vector[internal_from_node_seq_no].m_outgoing_link_seq_no_vector.push_back(link.link_seq_no);
	// add this link to the corresponding node as part of outgoing node/link
	g_node_vector[internal_to_node_seq_no].m_incoming_link_seq_no_vector.push_back(link.link_seq_no);
	// add this link to the corresponding node as part of outgoing node/link
	g_node_vector[internal_from_node_seq_no].m_to_node_seq_no_vector.push_back(link.to_node_seq_no);
	// add this link to the corresponding node as part of outgoing node/link
	g_node_vector[internal_from_node_seq_no].m_to_node_2_link_seq_no_map[link.to_node_seq_no] = link.link_seq_no;

	g_link_vector.push_back(link);

	assignment.g_number_of_links++;
}



double g_CheckActivityNodes(Assignment& assignment)
{

	int activity_node_count = 0;
	for (int i = 0; i < g_node_vector.size(); i++)
	{

		if (g_node_vector[i].is_activity_node >= 1)
		{
			activity_node_count++;
		}
	}


	if (activity_node_count <= 1)
	{
		activity_node_count = 0;
		int sampling_rate = 10;

		for (int i = 0; i < g_node_vector.size(); i++)
		{

			if (i % sampling_rate == 0)
			{
				g_node_vector[i].is_activity_node = 10;//random generation
				activity_node_count++;
			}
		}

		//if (activity_node_count <= 1)
		//{
		//    activity_node_count = 0;
		//    sampling_rate = 2;

		//    for (int i = 0; i < g_node_vector.size(); i++)
		//    {

		//        if (i % sampling_rate == 0)
		//        {
		//            g_node_vector[i].is_activity_node = 10;//random generation
		//            activity_node_count++;
		//        }
		//    }
		//     still no activity nodes, define all nodes as activity nodes
		//    if (activity_node_count <= 1)
		//    {
		//        activity_node_count = 0;

		//        for (int i = 0; i < g_node_vector.size(); i++)
		//        {

		//            g_node_vector[i].is_activity_node = 10;//random generation
		//            activity_node_count++;
		//        }
		//    }
		//}


	}


	// calculate avg near by distance;
	double total_near_by_distance = 0;
	activity_node_count = 0;
	for (int i = 0; i < g_node_vector.size(); i++)
	{
		double min_near_by_distance = 100;
		if (g_node_vector[i].is_activity_node)
		{
			activity_node_count++;
			for (int j = 0; j < g_node_vector.size(); j++)
			{
				if (i != j && g_node_vector[j].is_activity_node)
				{



					double near_by_distance = g_calculate_p2p_distance_in_meter_from_latitude_longitude(g_node_vector[i].x, g_node_vector[i].y, g_node_vector[j].x, g_node_vector[j].y);

					if (near_by_distance < min_near_by_distance)
						min_near_by_distance = near_by_distance;

				}

			}

			total_near_by_distance += min_near_by_distance;
			activity_node_count++;
		}
	}

	double nearby_distance = total_near_by_distance / max(1, activity_node_count);
	return nearby_distance;

}

int g_detect_if_demand_data_provided(Assignment& assignment)
{

	CDTACSVParser parser;
	
	dtalog.output() << "[STATUS INFO] Reading file demand_file_list.csv..." << '\n';
	g_DTA_log_file << "[STATUS INFO] Reading file demand_file_list.csv..." << '\n';
	if (parser.OpenCSVFile("demand_file_list.csv", false))
	{
		while (parser.ReadRecord())
		{

			int file_sequence_no = 1;

			string format_type = "null";

			int demand_format_flag = 0;

			if (!parser.GetValueByFieldName("file_sequence_no", file_sequence_no))
				break;

			// skip negative sequence no
			if (file_sequence_no <= -1)
				continue;

			string file_name, demand_period_str, mode_type;
			parser.GetValueByFieldName("file_name", file_name);
			parser.GetValueByFieldName("format_type", format_type);

			if (format_type.find("column") != string::npos)  // or muliti-column
			{
				// read the file formaly after the test.
				CDTACSVParser parser;

				if (parser.OpenCSVFile(file_name, false)  == false)
				{
					if (file_name == "demand.csv")
						return 0;
					else
					{
						dtalog.output() << "[ERROR INFO] demand file " << file_name << " does not exist. Users are expected to handle this e.g., by checking if the file names are provided correctly in the demand_file_list.csv field.." << '\n';
						g_DTA_log_file << "[ERROR INFO] demand file " << file_name << " does not exist. Users are expected to handle this e.g., by checking if the file names are provided correctly in the demand_file_list.csv field.." << '\n';
						continue; 
					}
				}
				else
				{
					return 1; // colulmn format demand file is needed.
				}
			}
			if (format_type.find("matrix") != string::npos)
			{
				// read the file formaly after the test.
				CDTACSVParser parser;

				if (parser.OpenCSVFile(file_name, false))
				{
					return 0;
				}
				else
				{
					return 2; // matrix format demand file is needed.
				}
			}
			if (format_type.find("activity_plan") != string::npos)
			{
				// read the file formaly after the test.
				CDTACSVParser parser;

				if (parser.OpenCSVFile(file_name, false))
				{
					return 0;
				}
				else
				{
					return 2; // matrix format demand file is needed.
				}
			}
			if (format_type.find("route") != string::npos)
			{
				// read the file formaly after the test.
				CDTACSVParser parser;

				if (parser.OpenCSVFile(file_name, false))
				{
					return 0;
				}
				else
				{
					return 2; // matrix format demand file is needed.
				}
			}
		}
		parser.CloseCSVFile();
	}

	return 100;  //default
}

int g_detect_if_zones_defined_in_node_csv(Assignment& assignment)
{
	CDTACSVParser parser;

	int number_of_activity_nodes = 0;
	int number_of_boundary_nodes = 0;


	if (parser.OpenCSVFile("node.csv", true))
	{
		while (parser.ReadRecord())  // if this line contains [] mark, then we will also read field headers.
		{
			int node_id;
			if (!parser.GetValueByFieldName("node_id", node_id))
				continue;

			int zone_id = 0;
			int is_boundary = 0;
			parser.GetValueByFieldName("zone_id", zone_id);
			parser.GetValueByFieldName("is_boundary", is_boundary, false);

			if (zone_id >= 1)
			{
				number_of_activity_nodes++;
				if (is_boundary != 0)
				{
					number_of_boundary_nodes++;
				}
			}

		}

		parser.CloseCSVFile();

		if (number_of_activity_nodes >= 2)  // if node.csv or zone.csv have 2 more zones;
		{
			assignment.summary_file << ", number of activity nodes defined in node.csv=, " << number_of_activity_nodes << '\n';
			assignment.summary_file << ", number of boundary activity nodes defined in node.csv=, " << number_of_boundary_nodes << '\n';



			return number_of_activity_nodes;
		}

	}
	else
	{
		dtalog.output() << "[ERROR] The critical input GMNS file 'node.csv' is missing. Please make sure the file is included in the appropriate directory." << '\n';
		g_DTA_log_file << "[ERROR] The critical input GMNS file 'node.csv' is missing. Please make sure the file is included in the appropriate directory." << '\n';
		g_program_stop();
	}

	return 0;
}


void g_read_link_qvdf_data(Assignment& assignment)
{
	CDTACSVParser parser;

	if (parser.OpenCSVFile("link_qvdf.csv", true))
	{
		while (parser.ReadRecord())  // if this line contains [] mark, then we will also read field headers.
		{
			string data_type;
			parser.GetValueByFieldName("data_type", data_type);

			if (data_type == "vdf_code")
			{
				string vdf_code;
				parser.GetValueByFieldName("vdf_code", vdf_code);


				for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
				{
					int demand_period_id = assignment.g_DemandPeriodVector[tau].demand_period_id;
					CLink this_link;
					char CSV_field_name[50];
					bool VDF_required_field_flag = true;
					//					sprintf(CSV_field_name, "QVDF_plf%d", demand_period_id);
					//					parser.GetValueByFieldName(CSV_field_name, this_link.VDF_period[tau].peak_load_factor, VDF_required_field_flag, false);
					sprintf(CSV_field_name, "QVDF_alpha%d", demand_period_id);
					parser.GetValueByFieldName(CSV_field_name, this_link.VDF_period[tau].Q_alpha, VDF_required_field_flag, false);
					sprintf(CSV_field_name, "QVDF_beta%d", demand_period_id);
					parser.GetValueByFieldName(CSV_field_name, this_link.VDF_period[tau].Q_beta, VDF_required_field_flag, false);
					sprintf(CSV_field_name, "QVDF_cd%d", demand_period_id);
					parser.GetValueByFieldName(CSV_field_name, this_link.VDF_period[tau].Q_cd, VDF_required_field_flag, false);
					sprintf(CSV_field_name, "QVDF_cp%d", demand_period_id);
					parser.GetValueByFieldName(CSV_field_name, this_link.VDF_period[tau].Q_cp, VDF_required_field_flag, false);
					sprintf(CSV_field_name, "QVDF_n%d", demand_period_id);
					parser.GetValueByFieldName(CSV_field_name, this_link.VDF_period[tau].Q_n, VDF_required_field_flag, false);
					sprintf(CSV_field_name, "QVDF_s%d", demand_period_id);
					parser.GetValueByFieldName(CSV_field_name, this_link.VDF_period[tau].Q_s, VDF_required_field_flag, false);
					g_vdf_type_map[vdf_code].record_qvdf_data(this_link.VDF_period[tau], tau);
				}

			}
			else
			{

				int from_node_id;
				if (!parser.GetValueByFieldName("from_node_id", from_node_id))
					continue;

				int to_node_id;
				if (!parser.GetValueByFieldName("to_node_id", to_node_id))
					continue;

				// add the to node id into the outbound (adjacent) node list
				if (assignment.g_node_id_to_seq_no_map.find(from_node_id) == assignment.g_node_id_to_seq_no_map.end())
				{
					dtalog.output() << "[WARNING] from_node_id " << from_node_id << " in file link_qvdf.csv is not defined in node.csv." << '\n';
					g_DTA_log_file << "[WARNING] from_node_id " << from_node_id << " in file link_qvdf.csv is not defined in node.csv." << '\n';
					//has not been defined
					continue;
				}
				if (assignment.g_node_id_to_seq_no_map.find(to_node_id) == assignment.g_node_id_to_seq_no_map.end())
				{
					dtalog.output() << "[WARNING] to_node_id " << to_node_id << " in file link_qvdf.csv is not defined in node.csv." << '\n';
					g_DTA_log_file << "[WARNING] to_node_id " << to_node_id << " in file link_qvdf.csv is not defined in node.csv." << '\n';
					//has not been defined
					continue;
				}

				for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
				{
					int demand_period_id = assignment.g_DemandPeriodVector[tau].demand_period_id;
					// map external node number to internal node seq no.
					int internal_from_node_seq_no = assignment.g_node_id_to_seq_no_map[from_node_id];
					int internal_to_node_seq_no = assignment.g_node_id_to_seq_no_map[to_node_id];

					if (g_node_vector[internal_from_node_seq_no].m_to_node_2_link_seq_no_map.find(internal_to_node_seq_no) != g_node_vector[internal_from_node_seq_no].m_to_node_2_link_seq_no_map.end())
					{
						int link_seq_no = g_node_vector[internal_from_node_seq_no].m_to_node_2_link_seq_no_map[internal_to_node_seq_no];
						if (link_seq_no >= 0 && g_link_vector[link_seq_no].vdf_type == q_vdf  /*QVDF*/)  // data exist
						{
							CLink* p_link = &(g_link_vector[link_seq_no]);
							char CSV_field_name[50];
							bool VDF_required_field_flag = true;
							sprintf(CSV_field_name, "QVDF_plf%d", demand_period_id);
							parser.GetValueByFieldName(CSV_field_name, p_link->VDF_period[tau].Q_peak_load_factor, VDF_required_field_flag, false);
							sprintf(CSV_field_name, "QVDF_alpha%d", demand_period_id);
							parser.GetValueByFieldName(CSV_field_name, p_link->VDF_period[tau].Q_alpha, VDF_required_field_flag, false);
							sprintf(CSV_field_name, "QVDF_beta%d", demand_period_id);
							parser.GetValueByFieldName(CSV_field_name, p_link->VDF_period[tau].Q_beta, VDF_required_field_flag, false);
							sprintf(CSV_field_name, "QVDF_cd%d", demand_period_id);
							parser.GetValueByFieldName(CSV_field_name, p_link->VDF_period[tau].Q_cd, VDF_required_field_flag, false);
							sprintf(CSV_field_name, "QVDF_n%d", demand_period_id);
							parser.GetValueByFieldName(CSV_field_name, p_link->VDF_period[tau].Q_n, VDF_required_field_flag, false);
							sprintf(CSV_field_name, "QVDF_cp%d", demand_period_id);
							parser.GetValueByFieldName(CSV_field_name, p_link->VDF_period[tau].Q_cp, VDF_required_field_flag, false);
							sprintf(CSV_field_name, "QVDF_s%d", demand_period_id);
							parser.GetValueByFieldName(CSV_field_name, p_link->VDF_period[tau].Q_s, VDF_required_field_flag, false);

						}
					}
				}
			}
		}
		parser.CloseCSVFile();
	}

}

extern unsigned int g_RandomSeed;
extern void InitWELLRNG512a(unsigned int* init);

void g_detector_file_open_status(Assignment& assignment)
{
	FILE* g_pFilePathMOE = nullptr;

	fopen_ss(&g_pFilePathMOE, "final_summary.csv", "w");
	if (!g_pFilePathMOE)
	{
		dtalog.output() << "[ERROR] File final_summary.csv cannot be opened." << '\n';
		g_DTA_log_file << "[ERROR] File final_summary.csv cannot be opened." << '\n';
		return; 
	}
	else
	{
		fclose(g_pFilePathMOE);
	}
	fopen_ss(&g_pFilePathMOE, "od_performance_summary.csv", "w");
	if (!g_pFilePathMOE)
	{
		dtalog.output() << "[ERROR] File od_performance_summary.csv cannot be opened." << '\n';
		g_DTA_log_file << "[ERROR] File od_performance_summary.csv cannot be opened." << '\n';
		return; 
	}
	else
	{
		fclose(g_pFilePathMOE);
	}

	if (assignment.g_subarea_shape_points.size() >= 3)
	{

		fopen_ss(&g_pFilePathMOE, "subarea_link_performance.csv", "w");
		if (!g_pFilePathMOE)
		{
			dtalog.output() << "[ERROR] File subarea_link_performance.csv cannot be opened." << '\n';
			g_DTA_log_file << "[ERROR] File subarea_link_performance.csv cannot be opened." << '\n';
			return; 
		}
		else
		{
			fclose(g_pFilePathMOE);
		}
	}

	if (assignment.assignment_mode == simulation_dta)
	{
		fopen_ss(&g_pFilePathMOE, "agent.csv", "w");

		if (!g_pFilePathMOE)
		{
			dtalog.output() << "[ERROR] File agent.csv cannot be opened." << '\n';
			g_DTA_log_file << "[ERROR] File agent.csv cannot be opened." << '\n';
			return; 
		}
		else
		{
			fclose(g_pFilePathMOE);
		}

		fopen_ss(&g_pFilePathMOE, "td_link_performance.csv", "w");
		if (!g_pFilePathMOE)
		{
			dtalog.output() << "[ERROR] File td_link_performance.csv cannot be opened." << '\n';
			g_DTA_log_file << "[ERROR] File td_link_performance.csv cannot be opened." << '\n';
			return;
		}
		else
		{
			fclose(g_pFilePathMOE);
		}

	}



}
void g_read_input_data(Assignment& assignment)
{
	g_detector_file_open_status(assignment);

	unsigned int state[16];

	for (int k = 0; k < 16; ++k)
	{
		state[k] = k + g_RandomSeed;
	}

	InitWELLRNG512a(state);

	assignment.g_LoadingStartTimeInMin = 99999;
	assignment.g_LoadingEndTimeInMin = 0;

	//step 0:read demand period file
	CDTACSVParser parser_demand_period;

	dtalog.output() << "[PROCESS INFO] Step 1: Reading input data" << '\n';
	g_DTA_log_file << "[PROCESS INFO] Step 1: Reading input data" << '\n';
	dtalog.output() << "[PROCESS INFO] Step 1.1: Reading demand_period.csv..." << '\n';
	g_DTA_log_file << "[PROCESS INFO] Step 1.1: Reading demand_period.csv..." << '\n';
	assignment.summary_file << "[PROCESS INFO] Step 1.1: Reading demand_period.csv..." << '\n';


	if (parser_demand_period.OpenCSVFile("demand_period.csv", false))
	{
		while (parser_demand_period.ReadRecord())
		{
			CDemand_Period demand_period;

			if (!parser_demand_period.GetValueByFieldName("demand_period_id", demand_period.demand_period_id))
				break;

			if (!parser_demand_period.GetValueByFieldName("demand_period", demand_period.demand_period))
			{
				dtalog.output() << "[ERROR] Field demand_period in file demand_period cannot be read." << '\n';
				g_DTA_log_file << "[ERROR] Field demand_period in file demand_period cannot be read." << '\n';
				continue; 
			}

			vector<float> global_minute_vector;

			if (!parser_demand_period.GetValueByFieldName("time_period", demand_period.time_period))
			{
				dtalog.output() << "[ERROR] Field time_period in file demand_period cannot be read." << '\n';
				g_DTA_log_file << "[ERROR] Field time_period in file demand_period cannot be read." << '\n';
				continue; 
			}



			//input_string includes the start and end time of a time period with hhmm format
			global_minute_vector = g_time_parser(demand_period.time_period); //global_minute_vector incldue the starting and ending time

			if (global_minute_vector.size() == 2)
			{
				demand_period.starting_time_slot_no = global_minute_vector[0] / MIN_PER_TIMESLOT;  // read the data
				demand_period.ending_time_slot_no = global_minute_vector[1] / MIN_PER_TIMESLOT;    // read the data from settings.csv
				demand_period.time_period_in_hour = (global_minute_vector[1] - global_minute_vector[0]) / 60.0;
				demand_period.t2_peak_in_hour = (global_minute_vector[0] + global_minute_vector[1]) / 2 / 60;

				//g_fout << global_minute_vector[0] << '\n';
				//g_fout << global_minute_vector[1] << '\n';


				string peak_time_str;
				if (parser_demand_period.GetValueByFieldName("peak_time", peak_time_str, false))
				{
					demand_period.t2_peak_in_hour = g_timestamp_parser(peak_time_str) / 60.0;

				}

			}

			assignment.demand_period_to_seqno_mapping[demand_period.demand_period] = assignment.g_DemandPeriodVector.size();
			assignment.g_DemandPeriodVector.push_back(demand_period);

			CDeparture_time_Profile dep_time;
			dep_time.starting_time_slot_no = demand_period.starting_time_slot_no;
			dep_time.ending_time_slot_no = demand_period.ending_time_slot_no;

			for (int s = 0; s <= 96 * 3; s++)
			{
				dep_time.departure_time_ratio[s] = 1.0 / 300.0;
			}

//			dtalog.output() << "[DATA INFO] A default flat departure time profile is used..." << '\n';
//			g_DTA_log_file << "[DATA INFO] A default flat departure time profile is used..." << '\n';
			dep_time.compute_cumulative_profile(demand_period.starting_time_slot_no, demand_period.ending_time_slot_no, false);

			if (assignment.g_DepartureTimeProfileVector.size() == 0)
			{
				//default profile
				assignment.g_DepartureTimeProfileVector.push_back(dep_time);
			}

			assignment.summary_file << ",demand_period= " << demand_period.demand_period_id <<
				", " << demand_period.demand_period.c_str() << ",time_period=" << demand_period.time_period.c_str() << '\n';

		}

		parser_demand_period.CloseCSVFile();

		if (assignment.g_DemandPeriodVector.size() == 0)
		{
			dtalog.output() << "[ERROR] File demand_period has no information." << '\n';
			g_DTA_log_file << "[ERROR] File demand_period has no information." << '\n';
			return; 
		}
	}
	else
	{
		dtalog.output() << "[ERROR] File demand_period.csv cannot be opened.\n It might be currently used and locked by EXCEL." << '\n';
		g_DTA_log_file << "[ERROR] File demand_period.csv cannot be opened.\n It might be currently used and locked by EXCEL." << '\n';
		return;
	}

	dtalog.output() << "[DATA INFO] number of demand periods = " << assignment.g_DemandPeriodVector.size() << '\n';
	g_DTA_log_file << "[DATA INFO] number of demand periods = " << assignment.g_DemandPeriodVector.size() << '\n';

	assignment.g_number_of_demand_periods = assignment.g_DemandPeriodVector.size();
	g_number_of_active_demand_perioids = assignment.g_DemandPeriodVector.size();;

	if (assignment.g_number_of_demand_periods > MAX_TIMEPERIODS)
	{
		dtalog.output() << "[ERROR] the number of demand periods in settings.csv os greater than the internal size of MAX_TIMEPERIODS.\nPlease contact developers" << '\n';
		g_DTA_log_file << "[ERROR] the number of demand periods in settings.csv os greater than the internal size of MAX_TIMEPERIODS.\nPlease contact developers" << '\n';
		g_program_stop();
	}
	//step 1:read demand type file

	CDTACSVParser parser_mode_type;
	dtalog.output() << "[PROCESS INFO] Step 1.2: Reading mode_type.csv..." << '\n';
	g_DTA_log_file << "[PROCESS INFO] Step 1.2: Reading mode_type.csv..." << '\n';
	assignment.summary_file << "[PROCESS INFO] Step 1.2: Reading file mode_type.csv..." << '\n';
	if (parser_mode_type.OpenCSVFile("mode_type.csv", false))
	{

		assignment.g_ModeTypeVector.clear();
		while (parser_mode_type.ReadRecord())
		{
			Cmode_type mode_type;

			if (!parser_mode_type.GetValueByFieldName("mode_type", mode_type.mode_type))
				break;

			int activate_flag = 1;
			parser_mode_type.GetValueByFieldName("activate", activate_flag,false,false);

			if (activate_flag == 0)
				continue;

			mode_type.mode_type_no = assignment.g_ModeTypeVector.size() + 1;

			int mode_specific_assignment_flag = 1;
			if (parser_mode_type.GetValueByFieldName("multimodal_dedicated_assignment_flag", mode_specific_assignment_flag, false, false) == false)
			{
				dtalog.output() << "[WARNING] Field multimodal_dedicated_assignment_flag is missing from mode_type.csv, a default value of 1 is used.  Please add flags to clearly identify assignment cost structures." << '\n';
				g_DTA_log_file << "[WARNING] Field multimodal_dedicated_assignment_flag is missing from mode_type.csv, a default value of 1 is used.  Please add flags to clearly identify assignment cost structures." << '\n';

			}



			//substring overlapping checking

			{
				for (int at = 0; at < assignment.g_ModeTypeVector.size(); at++)
				{
					if (assignment.g_ModeTypeVector[at].mode_type.find(mode_type.mode_type) != string::npos)
					{
						dtalog.output() << "[ERROR] Error substring duplication checking : mode_type = " << assignment.g_ModeTypeVector[at].mode_type.c_str() <<
							" in section mode_type is overlapping with " << mode_type.mode_type.c_str() << ". Please add flags such as to avoid overlapping in the use of allowe_uses field.";
				
						g_DTA_log_file << "[ERROR] Error substring duplication checking : mode_type = " << assignment.g_ModeTypeVector[at].mode_type.c_str() <<
							" in section mode_type is overlapping with " << mode_type.mode_type.c_str() << ". Please add flags such as to avoid overlapping in the use of allowe_uses field.";

					}

				}

			}

			parser_mode_type.GetValueByFieldName("vot", mode_type.value_of_time, false, false);

			// scan through the map with different node sum for different paths

			parser_mode_type.GetValueByFieldName("person_occupancy", mode_type.OCC,false,false);
			parser_mode_type.GetValueByFieldName("desired_speed_ratio", mode_type.DSR, false, false);

			parser_mode_type.GetValueByFieldName("headway_in_sec", mode_type.time_headway_in_sec, false, false);
			parser_mode_type.GetValueByFieldName("display_code", mode_type.display_code, false);
			parser_mode_type.GetValueByFieldName("DTM_real_time_info_type", mode_type.real_time_information_type, false,false);

			if (mode_type.real_time_information_type == 1  && mode_specific_assignment_flag ==1)
			{
					dtalog.output() << "[WARNING] The mode type '"
					<< mode_type.mode_type.c_str()
					<< "' specified in the 'mode_type.csv' file is not intended to have a dedicated travel time function,  because 'DTM_real_time_info_type' is set to '1' in 'link_type.csv'. The 'multimodal_dedicated_assignment_flag' has been reset to '0'.";
					mode_specific_assignment_flag = 0;

					g_DTA_log_file << "[WARNING] The mode type '"
						<< mode_type.mode_type.c_str()
						<< "' specified in the 'mode_type.csv' file is not intended to have a dedicated travel time function,  because 'DTM_real_time_info_type' is set to '1' in 'link_type.csv'. The 'multimodal_dedicated_assignment_flag' has been reset to '0'.";
					mode_specific_assignment_flag = 0;
			}
				

			if (mode_specific_assignment_flag == 1 || assignment.g_ModeTypeVector.size() == 0)
			{
				mode_type.mode_specific_assignment_flag= mode_specific_assignment_flag;

				if(mode_type.real_time_information_type!=0)
				{	
					dtalog.output() << "[DATA INFO] The mode type "
						<< mode_type.mode_type.c_str()
						<< " defined in the mode_type.csv file will utilize dedicated multimodal assignment features, complete with its own capacity and free_speed parameters in link_type.csv";

					g_DTA_log_file << "[DATA INFO] The mode type "
						<< mode_type.mode_type.c_str()
						<< " defined in the mode_type.csv file will utilize dedicated multimodal assignment features, complete with its own capacity and free_speed parameters in link_type.csv";
				}

			}

			if (mode_type.real_time_information_type == 1)  // real time info
			{
				assignment.g_number_of_real_time_mode_types++;
			}

			if (mode_type.real_time_information_type == 2)  //dms
			{
				assignment.g_number_of_DMS_mode_types++;
			}	
			
			parser_mode_type.GetValueByFieldName("access_node_type", mode_type.access_node_type, false);

			if (mode_type.access_node_type.size() > 0)
			{
				parser_mode_type.GetValueByFieldName("access_speed", mode_type.access_speed);
				parser_mode_type.GetValueByFieldName("access_distance_lb", mode_type.access_distance_lb);
				parser_mode_type.GetValueByFieldName("access_distance_ub", mode_type.access_distance_ub);

				if (mode_type.access_distance_ub < 100)
				{
					dtalog.output() << "[ERROR] access_distance_ub = " << mode_type.access_distance_ub << "< 100. Please ensure the unit is meter." << '\n';
					g_DTA_log_file << "[ERROR] access_distance_ub = " << mode_type.access_distance_ub << "< 100. Please ensure the unit is meter." << '\n';
					g_program_stop();
				}
				parser_mode_type.GetValueByFieldName("acecss_link_k", mode_type.acecss_link_k);
			}

			assignment.mode_type_2_seqno_mapping[mode_type.mode_type] = assignment.g_ModeTypeVector.size();

			assignment.g_ModeTypeVector.push_back(mode_type);
			assignment.summary_file << "mode_type =, " << mode_type.mode_type.c_str() << ", real time info flag = " << mode_type.real_time_information_type << '\n';

		}

		assignment.g_number_of_mode_types = assignment.g_ModeTypeVector.size();
		g_number_of_active_mode_types = assignment.g_ModeTypeVector.size();

		parser_mode_type.CloseCSVFile();
	}

	if (assignment.g_ModeTypeVector.size() == 0)
	{
		dtalog.output() << "[ERROR] File mode_type does not have information." << '\n';
		g_DTA_log_file << "[ERROR] File mode_type does not have information." << '\n';
		g_program_stop();
	}

	if (assignment.g_ModeTypeVector.size() >= MAX_MODETYPES)
	{
		dtalog.output() << "[ERROR] mode_type = " << assignment.g_ModeTypeVector.size() << " in section mode_type is too large. " << '\n' << "MAX_MODETYPES = " << MAX_MODETYPES << "Please contact program developers!";
		g_DTA_log_file << "[ERROR] mode_type = " << assignment.g_ModeTypeVector.size() << " in section mode_type is too large. " << '\n' << "MAX_MODETYPES = " << MAX_MODETYPES << "Please contact program developers!";
		g_program_stop();
	}

	dtalog.output() << "[DATA INFO] number of mode types = " << assignment.g_ModeTypeVector.size() << '\n';
	g_DTA_log_file << "[DATA INFO] number of mode types = " << assignment.g_ModeTypeVector.size() << '\n';


	////// Step 1.2: Reading activity_travel_pattern.csv...
	////dtalog.output() << "[PROCESS INFO] Step 1.25: Reading optional activity_travel_pattern.csv..." << '\n';
	////g_DTA_log_file << "[PROCESS INFO] Step 1.25: Reading optional activity_travel_pattern.csv..." << '\n';
	////assignment.summary_file << "[PROCESS INFO] Step 1.25: Reading optional activity_travel_pattern.csv..." << '\n';

	//CDTACSVParser parser_activity_travel_pattern;

	//if (parser_activity_travel_pattern.OpenCSVFile("activity_travel_pattern.csv", false))
	//{
	//	while (parser_activity_travel_pattern.ReadRecord())
	//	{
	//		CActivityTravelPattern travel_pattern;

	//		if (!parser_activity_travel_pattern.GetValueByFieldName("activity_travel_pattern_index", travel_pattern.activity_travel_pattern_index))
	//			break;

	//		if (!parser_activity_travel_pattern.GetValueByFieldName("activity_travel_pattern_id", travel_pattern.activity_travel_pattern_id))
	//		{
	//			dtalog.output() << "[ERROR] Field activity_travel_pattern_id in file activity_travel_pattern.csv cannot be read." << '\n';
	//			g_DTA_log_file << "[ERROR] Field activity_travel_pattern_id in file activity_travel_pattern.csv cannot be read." << '\n';
	//			g_program_stop();
	//		}

	//		if (assignment.g_ActivityTravelPatternMap.find(travel_pattern.activity_travel_pattern_id) != assignment.g_ActivityTravelPatternMap.end())
	//		{
	//			dtalog.output() << "[ERROR] Field activity_travel_pattern_id = " << travel_pattern.activity_travel_pattern_id.c_str() << " has been defined twice in file activity_travel_pattern.csv." << '\n';
	//			g_DTA_log_file << "[ERROR] Field activity_travel_pattern_id = " << travel_pattern.activity_travel_pattern_id.c_str() << " has been defined twice in file activity_travel_pattern.csv." << '\n';
	//			g_program_stop();
	//		}


	//		if (!parser_activity_travel_pattern.GetValueByFieldName("number_of_trips", travel_pattern.number_of_trips))
	//		{
	//			dtalog.output() << "[ERROR] Field number_of_trips in file activity_travel_pattern.csv cannot be read." << '\n';
	//			g_DTA_log_file << "[ERROR] Field number_of_trips in file activity_travel_pattern.csv cannot be read." << '\n';
	//			g_program_stop();
	//		}

	//		std::string demand_period_vector_str;
	//		if (!parser_activity_travel_pattern.GetValueByFieldName("demand_period_chain", demand_period_vector_str))
	//		{
	//			dtalog.output() << "[ERROR] Field demand_period_vector in file activity_travel_pattern.csv cannot be read." << '\n';
	//			g_DTA_log_file << "[ERROR] Field demand_period_vector in file activity_travel_pattern.csv cannot be read." << '\n';
	//			g_program_stop();
	//		}

	//		std::string mode_vector_str;
	//		if (!parser_activity_travel_pattern.GetValueByFieldName("mode_chain", mode_vector_str))
	//		{
	//			dtalog.output() << "[ERROR] Field mode_chain in file activity_travel_pattern.csv cannot be read." << '\n';
	//			g_DTA_log_file << "[ERROR] Field mode_chain in file activity_travel_pattern.csv cannot be read." << '\n';
	//			g_program_stop();
	//		}

	//		travel_pattern.mode_chain_str = mode_vector_str;
	//		travel_pattern.demand_period_chain_str = demand_period_vector_str;

	//		// split the string into a vector by ';' as the delimiter and set the demand period vector

	//		int demand_period_size = g_ParserStringSequence(demand_period_vector_str, travel_pattern.demand_period_chain);

	//		int mode_size = g_ParserStringSequence(mode_vector_str, travel_pattern.mode_chain);

	//		if (demand_period_size != travel_pattern.number_of_trips)
	//		{
	//			dtalog.output() << "[ERROR] demand_period_size != number_of_trips in file activity_travel_pattern.csv." << '\n';
	//			g_DTA_log_file << "[ERROR] demand_period_size != number_of_trips in file activity_travel_pattern.csv." << '\n';
	//			g_program_stop();
	//		}

	//		if (mode_size != travel_pattern.number_of_trips)
	//		{
	//			dtalog.output() << "[ERROR] mode_size != number_of_trips in file activity_travel_pattern.csv. for index = " <<
	//			g_DTA_log_file << "[ERROR] mode_size != number_of_trips in file activity_travel_pattern.csv. for index = " <<
	//				travel_pattern.activity_travel_pattern_index << " id = " << travel_pattern.activity_travel_pattern_id.c_str() << '\n';
	//			g_program_stop();
	//		}

	//		for (int trip_no = 0; trip_no < travel_pattern.mode_chain.size(); trip_no++)
	//		{
	//			std::string mode_type = travel_pattern.mode_chain[trip_no];

	//			if (assignment.mode_type_2_seqno_mapping.find(mode_type) != assignment.mode_type_2_seqno_mapping.end())
	//			{
	//				int mode_type_no = assignment.mode_type_2_seqno_mapping[mode_type];
	//				travel_pattern.mode_vector.push_back(mode_type_no);
	//			}
	//			else
	//			{
	//				dtalog.output() << "[ERROR] Field mode_chain  in file activity_travel_pattern.csv has a value " << mode_type.c_str() << "not defined yet in mode_type.csv." << '\n';
	//				g_DTA_log_file << "[ERROR] Field mode_chain  in file activity_travel_pattern.csv has a value " << mode_type.c_str() << "not defined yet in mode_type.csv." << '\n';
	//				g_program_stop();
	//			}

	//			std::string assignment_period = travel_pattern.demand_period_chain[trip_no];

	//			if (assignment.demand_period_to_seqno_mapping.find(assignment_period) != assignment.demand_period_to_seqno_mapping.end())
	//			{
	//				int assignment_period_no = assignment.demand_period_to_seqno_mapping[mode_type];
	//				travel_pattern.demand_period_vector.push_back(assignment_period_no);
	//			}
	//			else
	//			{
	//				dtalog.output() << "[ERROR] Field mode_chain  in file activity_travel_pattern.csv has a value " << mode_type << "not defined yet in mode_type.csv." << '\n';
	//				g_DTA_log_file << "[ERROR] Field mode_chain  in file activity_travel_pattern.csv has a value " << mode_type << "not defined yet in mode_type.csv." << '\n';
	//				g_program_stop();
	//			}

	//		}



	//		if (mode_size != travel_pattern.number_of_trips)
	//		{
	//			dtalog.output() << "[ERROR] mode_size != travel_pattern.number_of_trips in file activity_travel_pattern.csv for record" << travel_pattern.activity_travel_pattern_id.c_str() << '\n';
	//			g_DTA_log_file << "[ERROR] mode_size != travel_pattern.number_of_trips in file activity_travel_pattern.csv for record" << travel_pattern.activity_travel_pattern_id.c_str() << '\n';
	//			g_program_stop();
	//		}

	//		if (demand_period_size != travel_pattern.number_of_trips)
	//		{
	//			dtalog.output() << "[ERROR] demand_period_size != travel_pattern.number_of_trips in file activity_travel_pattern.csv for record" << travel_pattern.activity_travel_pattern_id.c_str() << '\n';
	//			g_DTA_log_file << "[ERROR] demand_period_size != travel_pattern.number_of_trips in file activity_travel_pattern.csv for record" << travel_pattern.activity_travel_pattern_id.c_str() << '\n';
	//			g_program_stop();
	//		}

	//		// split the string into a vector by ';' as the delimiter and set the mode vector
	//	//	travel_pattern.setModes(g_string_to_vector(mode_vector_str, ';'));

	//		assignment.g_ActivityTravelPatternMap[travel_pattern.activity_travel_pattern_id] = travel_pattern;

	//		assignment.summary_file << "[DATA INFO] activity_travel_pattern_index= " << travel_pattern.activity_travel_pattern_id <<
	//			", activity_travel_pattern_id=" << travel_pattern.activity_travel_pattern_id.c_str() << ", number_of_trips=" << travel_pattern.number_of_trips << '\n';
	//	}

	//	parser_activity_travel_pattern.CloseCSVFile();

	//	if (assignment.g_ActivityTravelPatternMap.size() == 0)
	//	{
	//		dtalog.output() << "[ERROR] File activity_travel_pattern has no information." << '\n';
	//		g_DTA_log_file << "[ERROR] File activity_travel_pattern has no information." << '\n';
	//		//			g_program_stop();
	//	}
	//}
	//else
	//{
	//	//dtalog.output() << "[WARNING] File activity_travel_pattern.csv cannot be opened. Note that, file activity_travel_pattern.csv is optional." << '\n';
	//	//g_DTA_log_file << "[WARNING] File activity_travel_pattern.csv cannot be opened. Note that, file activity_travel_pattern.csv is optional." << '\n';
	//	//	g_program_stop();
	//}



	dtalog.output() << "[PROCESS INFO] Step 1.3: Reading link_type.csv" << '\n';
	g_DTA_log_file << "[PROCESS INFO] Step 1.3: Reading link_type.csv" << '\n';

	CDTACSVParser parser_link_type;

	int emission_log_count = 0; 
	int meu_log_count = 0; 
	int peak_load_factor_log_count = 0; 
	int allowed_use_log_count = 0; 
	if (parser_link_type.OpenCSVFile("link_type.csv", false))
	{
		// create a special link type as virtual connector
		CLinkType element_vc;
		// -1 is for virutal connector
		element_vc.link_type = -1;
		element_vc.type_code = "c";
		element_vc.traffic_flow_code = spatial_queue;
		assignment.g_LinkTypeMap[element_vc.link_type] = element_vc;
		//end of create special link type for virtual connectors

		int line_no = 0;

		std::vector <int> overlapping_link_type_vector; 
		while (parser_link_type.ReadRecord())
		{
			CLinkType element;

			if (!parser_link_type.GetValueByFieldName("link_type", element.link_type))
			{
				if (line_no == 0)
				{
					dtalog.output() << "[ERROR] Field link_type cannot be found in file link_type.csv." << '\n';
					g_DTA_log_file << "[ERROR] Field link_type cannot be found in file link_type.csv." << '\n';
					continue; 
				}
				else
				{
					// read empty line
					break;
				}
			}

			parser_link_type.GetValueByFieldName("link_type_name", element.link_type_name);

			if (assignment.g_LinkTypeMap.find(element.link_type) != assignment.g_LinkTypeMap.end())
			{
				overlapping_link_type_vector.push_back(element.link_type);

				continue;
			}

			for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
			{
				char CSV_field_name[50];
				sprintf(CSV_field_name, "allowed_uses_p%d", tau + 1);
				
				if (parser_link_type.GetValueByFieldName(CSV_field_name, element.allow_uses_period[tau], false, false))
				{
					if (line_no == 0 && allowed_use_log_count < 2) {
						allowed_use_log_count++;
						dtalog.output() << "[WARNING] Field '" << CSV_field_name << "' not found in 'link_type.csv'. All modes will be allowed for link type:" << element.link_type_name << ". Consider adding '" << CSV_field_name << "' to the 'link_type.csv' for more accurate results." << '\n';
						g_DTA_log_file << "[WARNING] Field '" << CSV_field_name << "' not found in 'link_type.csv'. All modes will be allowed for link type:" << element.link_type_name << ". Consider adding '" << CSV_field_name << "' to the 'link_type.csv' for more accurate results." << '\n';

					}

				}
				for (int at = 0; at < assignment.g_ModeTypeVector.size(); at++)
				{
					sprintf(CSV_field_name, "peak_load_factor_p%d_%s", tau + 1, assignment.g_ModeTypeVector[at].mode_type.c_str());
					if (parser_link_type.GetValueByFieldName(CSV_field_name, element.peak_load_factor_period_at[tau][at], false) == false)
					{
						if (line_no == 0 && peak_load_factor_log_count < 2) {
							peak_load_factor_log_count++;
							dtalog.output() << "[WARNING] Field '" << CSV_field_name << "' not found in 'link_type.csv'. The default peak load factor 1.0 was used. Consider adding '" << CSV_field_name << "' to the 'link_type.csv' for more accurate results." << '\n';
							g_DTA_log_file << "[WARNING] Field '" << CSV_field_name << "' not found in 'link_type.csv'. The default peak load factor 1.0 was used. Consider adding '" << CSV_field_name << "' to the 'link_type.csv' for more accurate results." << '\n';

						}

					}
				}
			}
			//=----

			char CSV_field_name[50];
			int at_base = 0;
			double lanes_mode_type = -1; // default
			sprintf(CSV_field_name, "lanes_%s", assignment.g_ModeTypeVector[at_base].mode_type.c_str());
			if (parser_link_type.GetValueByFieldName(CSV_field_name, lanes_mode_type, false, false) == true)
			{

				dtalog.output() << "[WARNING] Field '" << CSV_field_name << "' is found in 'link_type.csv'. ";
				g_DTA_log_file << "[WARNING] Field '" << CSV_field_name << "' is found in 'link_type.csv'. ";
				dtalog.output() << " To ensure data consistency, we advise users to specify the number of lanes for the main base mode " << assignment.g_ModeTypeVector[at_base].mode_type.c_str() << " in both the 'lanes' and 'lanes_s0' fields in the master 'link.csv' file. To avoid confusion, please remove this field from the 'link_type.csv'." << '\n';
				g_DTA_log_file << " To ensure data consistency, we advise users to specify the number of lanes for the main base mode " << assignment.g_ModeTypeVector[at_base].mode_type.c_str() << " in both the 'lanes' and 'lanes_s0' fields in the master 'link.csv' file. To avoid confusion, please remove this field from the 'link_type.csv'." << '\n';

			}



			for (int at = 0; at < assignment.g_ModeTypeVector.size(); at++)
			{
				char CSV_field_name[50];

				if (assignment.g_ModeTypeVector[at].mode_specific_assignment_flag == 1)
				{
					double capacity_at = 2000; // default
					sprintf(CSV_field_name, "capacity_%s", assignment.g_ModeTypeVector[at].mode_type.c_str());
					if (parser_link_type.GetValueByFieldName(CSV_field_name, capacity_at,false,false)==false)
					{
						if(line_no == 0){
							dtalog.output() << "[WARNING] Field '" << CSV_field_name << "' not found in 'link_type.csv'. The default capacity of 2000 was used. Consider adding '" << CSV_field_name << "' to the 'link_type.csv' for more accurate results." << '\n';
							g_DTA_log_file << "[WARNING] Field '" << CSV_field_name << "' not found in 'link_type.csv'. The default capacity of 2000 was used. Consider adding '" << CSV_field_name << "' to the 'link_type.csv' for more accurate results." << '\n';
							
						}
						capacity_at = 2000;
					}


					if (capacity_at > 0.1)  // log
					{
						element.capacity_at[at] = capacity_at;
					}
				}
				//----

				double free_speed_at = -1; // default

				sprintf(CSV_field_name, "free_speed_%s", assignment.g_ModeTypeVector[at].mode_type.c_str());
				if (parser_link_type.GetValueByFieldName(CSV_field_name, free_speed_at, false, false) == false)
				{
					if (line_no == 0) {
						dtalog.output() << "[WARNING] Field '" << CSV_field_name << "' not found in 'link_type.csv'. The default free speed 60 was used. Consider adding '" << CSV_field_name << "' to the 'link_type.csv' for more accurate results." << '\n';
						g_DTA_log_file << "[WARNING] Field '" << CSV_field_name << "' not found in 'link_type.csv'. The default free speed 60 was used. Consider adding '" << CSV_field_name << "' to the 'link_type.csv' for more accurate results." << '\n';
					}
					free_speed_at = 60;
				}


				if (free_speed_at > 0.1)  // log
				{
					element.free_speed_at[at] = free_speed_at;
				}


				if (at >= 1 && assignment.g_ModeTypeVector[at].mode_specific_assignment_flag == 1)
				{

					if (assignment.g_speed_unit_flag == 1)  // mph;
						free_speed_at = free_speed_at / 1.609; // convert from mile per hour to km per hour


					double lanes_mode_type = -1; // default

					sprintf(CSV_field_name, "lanes_%s", assignment.g_ModeTypeVector[at].mode_type.c_str());
					if (parser_link_type.GetValueByFieldName(CSV_field_name, lanes_mode_type, false, false) == false)
					{
						if (line_no == 0) {
							dtalog.output() << "[WARNING] Field '" << CSV_field_name << "' not found in 'link_type.csv'. The default number of lanes 0 was used. Consider adding '" << CSV_field_name << "' to the 'link_type.csv' for more accurate results." << '\n';
							g_DTA_log_file << "[WARNING] Field '" << CSV_field_name << "' not found in 'link_type.csv'. The default number of lanes 0 was used. Consider adding '" << CSV_field_name << "' to the 'link_type.csv' for more accurate results." << '\n';
						}
						lanes_mode_type = 0;
					}

					element.lanes_mode_type[at] = lanes_mode_type;
				}

				//				element.lanes_mode_type[at] = lanes_mode_type;

				for (int at2 = 0; at2 < assignment.g_ModeTypeVector.size(); at2++)
				{
					double meu_value = 0.0;

					if (at2 != at)  // let us use the MEU as 0 for simplicity 
						meu_value = 0.0;

					if (at2 == at)
					{
						element.meu_matrix[at][at2] = 1.0;
						continue; 
					}

					sprintf(CSV_field_name, "meu_%s_%s", assignment.g_ModeTypeVector[at].mode_type.c_str(), assignment.g_ModeTypeVector[at].mode_type.c_str());
					if(parser_link_type.GetValueByFieldName(CSV_field_name, meu_value, false, false)== false)
					{
						if (line_no == 0 && meu_log_count < 3) {
							meu_log_count++; 

							dtalog.output() << "[WARNING] Field '" << CSV_field_name << "' not found in 'link_type.csv'. The MEU = 0.0 was used. Consider adding '" << CSV_field_name << "' to the 'link_type.csv' for more accurate results." << '\n';
							g_DTA_log_file << "[WARNING] Field '" << CSV_field_name << "' not found in 'link_type.csv'. The MEU = 0.0 was used. Consider adding '" << CSV_field_name << "' to the 'link_type.csv' for more accurate results." << '\n';

						}
						element.meu_matrix[at][at2] = meu_value;

					}
				}

			}



			string traffic_flow_code_str;
			parser_link_type.GetValueByFieldName("type_code", element.type_code, true);

			string vdf_type_str;

			element.vdf_type = bpr_vdf;
			parser_link_type.GetValueByFieldName("vdf_type", vdf_type_str, false);


			if (vdf_type_str == "bpr")
				element.vdf_type = bpr_vdf;
			if (vdf_type_str == "qvdf")
				element.vdf_type = q_vdf;


			element.traffic_flow_code = spatial_queue;

			parser_link_type.GetValueByFieldName("traffic_flow_model", traffic_flow_code_str, false);
			parser_link_type.GetValueByFieldName("k_jam_km", element.k_jam, false);

			// by default bpr


			if (traffic_flow_code_str == "point_queue")
				element.traffic_flow_code = point_queue;

			if (traffic_flow_code_str == "spatial_queue")
				element.traffic_flow_code = spatial_queue;

			if (traffic_flow_code_str == "kw")
				element.traffic_flow_code = kinemative_wave;

			//				dtalog.output() << "important: traffic_flow_code on link type " << element.link_type << " is " << element.traffic_flow_code << '\n';
			//				g_DTA_log_file << "important: traffic_flow_code on link type " << element.link_type << " is " << element.traffic_flow_code << '\n';

			assignment.summary_file << "link_type =, " << element.link_type << ", link_type_name = " << element.link_type_name.c_str() << '\n';

			for (int at = 0; at < assignment.g_ModeTypeVector.size(); at++)
			{
				assignment.summary_file << ",mode_type= " << assignment.g_ModeTypeVector[at].mode_type.c_str() <<
					",cap= " << element.capacity_at[at] << ",free_speed=" << element.free_speed_at[at] << '\n';
				//",lanes = " << element.lanes_mode_type[at] << '\n';
			}


			// Loop over all mode types (e.g. auto, bike, bus, etc.)
			for (int at = 0; at < assignment.g_ModeTypeVector.size(); at++)
			{
				// Generate the field name in the CSV file for CO2 emissions for this mode type
				char CSV_field_name[50];
				sprintf(CSV_field_name, "emissions_%s_co2", assignment.g_ModeTypeVector[at].mode_type.c_str());

				string emissions_co2_str;

				// Check if the field could not be found in the CSV file
				if (parser_link_type.GetValueByFieldName(CSV_field_name, emissions_co2_str, false, false) == false)
				{
					// If the field cannot be found, output a warning and use a default value of 0.0
					if (line_no == 0 && emission_log_count < 4) {
						emission_log_count++;
						dtalog.output() << "[WARNING] Field '" << CSV_field_name << "' not found in 'link_type.csv'. The default value of 0 was used. Consider adding '" << CSV_field_name << "' to the 'link_type.csv' for more accurate results." << '\n';
						g_DTA_log_file << "[WARNING] Field '" << CSV_field_name << "' not found in 'link_type.csv'. The default value of 0 was used. Consider adding '" << CSV_field_name << "' to the 'link_type.csv' for more accurate results." << '\n';
					}
				}

				// Convert the string containing CO2 emissions coefficients into a vector of doubles
				std::vector<double> emissions_co2_coeff_vector;
				g_ParserDoubleSequence(emissions_co2_str, emissions_co2_coeff_vector);

				// Store up to the first 4 emissions coefficients in the element's matrix
				for (int i = 0; i < min(4, emissions_co2_coeff_vector.size()); i++)
				{
					element.emissions_co2_matrix[at][i] = emissions_co2_coeff_vector[i];
				}

				// Repeat the same process for NOx emissions:
				sprintf(CSV_field_name, "emissions_%s_nox", assignment.g_ModeTypeVector[at].mode_type.c_str());

				string emissions_nox_str;
				if (parser_link_type.GetValueByFieldName(CSV_field_name, emissions_nox_str, false, false) == false)
				{
					if (line_no == 0 && emission_log_count < 4) {
						emission_log_count++; 
						dtalog.output() << "[WARNING] Field '" << CSV_field_name << "' not found in 'link_type.csv'. The default value of 0 was used. Consider adding '" << CSV_field_name << "' to the 'link_type.csv' for more accurate results." << '\n';
						g_DTA_log_file << "[WARNING] Field '" << CSV_field_name << "' not found in 'link_type.csv'. The default value of 0 was used. Consider adding '" << CSV_field_name << "' to the 'link_type.csv' for more accurate results." << '\n';
					}
				}

				std::vector<double> emissions_nox_coeff_vector;
				g_ParserDoubleSequence(emissions_nox_str, emissions_nox_coeff_vector);

				for (int i = 0; i < min(4, emissions_nox_coeff_vector.size()); i++)
				{
					element.emissions_nox_matrix[at][i] = emissions_nox_coeff_vector[i];
				}
			}


			
			if (assignment.g_first_link_type < 0)
				assignment.g_first_link_type = element.link_type;

			assignment.g_LinkTypeMap[element.link_type] = element;
			line_no++;
		}

		if (overlapping_link_type_vector.size() > 0)
		{
			dtalog.output() << "[WARNING] the following link type(s) have been defined more than once in file link_type.csv:  ";
			g_DTA_log_file << "[WARNING] the following link type(s) have been defined more than once in file link_type.csv:  ";

			for (int i = 0; i < overlapping_link_type_vector.size(); i++)
			{
				dtalog.output() << overlapping_link_type_vector[i] << ", "; 
				g_DTA_log_file << overlapping_link_type_vector[i] << ", "; 
			}

			
		}


		

		parser_link_type.CloseCSVFile();
	}
	else
	{
		dtalog.output() << "[ERROR] File link_type.csv cannot be openned" << '\n';
		g_DTA_log_file << "[ERROR] File link_type.csv cannot be openned" << '\n';
		return; 
	}




	dtalog.output() << "[DATA INFO] number of link types = " << assignment.g_LinkTypeMap.size() << '\n';
	g_DTA_log_file << "[DATA INFO] number of link types = " << assignment.g_LinkTypeMap.size() << '\n';



	assignment.g_number_of_nodes = 0;
	assignment.g_number_of_links = 0;  // initialize  the counter to 0


	int number_of_zones = g_detect_if_zones_defined_in_node_csv(assignment);
	// = 1: normal
	//= 0, there are is boundary
	//=-1, no information



	if (number_of_zones <= 1)
	{
		CDTACSVParser parser_z;
		if (parser_z.OpenCSVFile("zone.csv", true))
		{
			parser_z.CloseCSVFile();
		}
		else
		{   // without zone.csv file

			if (g_TAZ_2_GMNS_zone_generation(assignment) == false)
			{
				g_grid_zone_generation(assignment);
			}

		}

	}

	int internal_node_seq_no = 0;
	// step 3: read node file


	std::map<int, int> zone_id_mapping;  // this is used to mark if this zone_id has been identified or not
	std::map<int, double> zone_id_x;
	std::map<int, double> zone_id_y;



	std::map<int, float> zone_id_production;
	std::map<int, float> zone_id_attraction;

	CDTACSVParser parser;

	int multmodal_activity_node_count = 0;

	dtalog.output() << "[PROCESS INFO] Step 1.35: Reading optional zone data in zone.csv..." << '\n';
	g_DTA_log_file << "[PROCESS INFO] Step 1.35: Reading optional zone data in zone.csv..." << '\n';


	if (parser.OpenCSVFile("zone.csv", true))
	{
		while (parser.ReadRecord())  // if this line contains [] mark, then we will also read field headers.
		{
			int zone_id = 0;
			if (!parser.GetValueByFieldName("zone_id", zone_id))
				continue;

			if (zone_id <= 0)
			{
				continue;
			}

			string access_node_vector_str;
			parser.GetValueByFieldName("access_node_vector", access_node_vector_str,false,false);

			std::vector<int> access_node_vector;

			g_ParserIntSequence(access_node_vector_str, access_node_vector);

			for (int i = 0; i < access_node_vector.size(); i++)
			{
				assignment.access_node_id_to_zone_id_map[access_node_vector[i]] = zone_id;
				zone_id_mapping[zone_id] = access_node_vector[i];

			}

			float production = 0;
			float attraction = 0;
			parser.GetValueByFieldName("production", production, false);
			parser.GetValueByFieldName("attraction", attraction, false);


			zone_id_production[zone_id] = production;
			zone_id_attraction[zone_id] = attraction;
			// push it to the global node vector
		}

		dtalog.output() << "[STATUS INFO] reading " << assignment.access_node_id_to_zone_id_map.size() << " access nodes from zone.csv.. " << '\n';
		g_DTA_log_file << "[STATUS INFO] reading " << assignment.access_node_id_to_zone_id_map.size() << " access nodes from zone.csv.. " << '\n';
		parser.CloseCSVFile();
	}



	dtalog.output() << "[PROCESS INFO] Step 1.4: Reading node data in node.csv..." << '\n';
	g_DTA_log_file << "[PROCESS INFO] Step 1.4: Reading node data in node.csv..." << '\n';
	std::map<int, int> zone_id_to_analysis_district_id_mapping;

	// MRM: stage 1: filtering, output: gate nodes, bridge nodes for micro and macro networks

	// MRM: stage 2: read node layer
	int max_layer_flag = 0;


	//	assignment.summary_file << ", # of node layers to be read," << max_layer_flag+1 << "," << '\n';

	for (int layer_no = 0; layer_no <= max_layer_flag; layer_no++)
	{
		dtalog.output() << "[STATUS INFO] Reading node data " << '\n';
		g_DTA_log_file << "[STATUS INFO] Reading node data " << '\n';

		// master file: reading nodes

		string file_name = "node.csv";


		if (parser.OpenCSVFile(file_name.c_str(), true))
		{

			while (parser.ReadRecord())  // if this line contains [] mark, then we will also read field headers.
			{
				// create a node object
				CNode node;

				int node_id;
				if (!parser.GetValueByFieldName("node_id", node_id))
					continue;


				if (node_id == 9439)
				{
					int trace_flag = 1;
				}


				if (max_layer_flag == 1) // MRM mode
				{
					if (layer_no == 0)  // macro node layer
					{  // skip all inside non-gate node

						if (assignment.g_node_id_to_MRM_subarea_mapping.find(node_id) != assignment.g_node_id_to_MRM_subarea_mapping.end())
						{
							if (assignment.g_node_id_to_MRM_subarea_mapping[node_id] == 1
								// inside and not defined as gate node
								&&
								assignment.g_macro_node_id_to_MRM_incoming_bridge_mapping.find(node_id) == assignment.g_macro_node_id_to_MRM_incoming_bridge_mapping.end()
								&&
								assignment.g_macro_node_id_to_MRM_outgoing_bridge_mapping.find(node_id) == assignment.g_macro_node_id_to_MRM_outgoing_bridge_mapping.end()


								)
							{
								continue;

							}

						}



					}


				}

				if (assignment.g_node_id_to_seq_no_map.find(node_id) != assignment.g_node_id_to_seq_no_map.end())
				{
					//has been defined
					continue;
				}




				double x_coord;
				double y_coord;

				parser.GetValueByFieldName("x_coord", x_coord, true, false);
				parser.GetValueByFieldName("y_coord", y_coord, true, false);


				// micro network filter:
				DTAGDPoint pt;
				pt.x = x_coord;
				pt.y = y_coord;

				assignment.g_node_id_to_seq_no_map[node_id] = internal_node_seq_no;

				node.node_id = node_id;
				node.node_seq_no = internal_node_seq_no;
				node.layer_no = layer_no;
				node.x = x_coord;
				node.y = y_coord;

				//if (assignment.g_micro_node_id_to_MRM_bridge_mapping.find(node_id) != assignment.g_micro_node_id_to_MRM_bridge_mapping.end())
				//	node.MRM_gate_flag = assignment.g_micro_node_id_to_MRM_bridge_mapping[node_id];

				//if (assignment.g_macro_node_id_to_MRM_incoming_bridge_mapping.find(node_id) != assignment.g_macro_node_id_to_MRM_incoming_bridge_mapping.end())
				//	node.MRM_gate_flag = assignment.g_macro_node_id_to_MRM_incoming_bridge_mapping[node_id];

				//if (assignment.g_macro_node_id_to_MRM_outgoing_bridge_mapping.find(node_id) != assignment.g_macro_node_id_to_MRM_outgoing_bridge_mapping.end())
				//	node.MRM_gate_flag = assignment.g_macro_node_id_to_MRM_outgoing_bridge_mapping[node_id];


				int zone_id = -1;

				parser.GetValueByFieldName("is_boundary", node.is_boundary, false, false);
				parser.GetValueByFieldName("node_type", node.node_type, false);// step 1 for adding access links: read node type
				parser.GetValueByFieldName("zone_id", zone_id);

				int analysis_district_id = 0;  // default

				parser.GetValueByFieldName("district_id", analysis_district_id, false);

				if (analysis_district_id >= MAX_ORIGIN_DISTRICTS)
				{
				dtalog.output() << "[ERROR] district_id in node.csv is larger than predefined value of "  << MAX_ORIGIN_DISTRICTS  << "." << '\n';
				g_DTA_log_file << "[ERROR] district_id in node.csv is larger than predefined value of "  << MAX_ORIGIN_DISTRICTS  << "." << '\n';
				g_program_stop();
				}

				if (node_id == 2235)
				{
					int idebug = 1;
				}
				//read from mapping created in zone file
				if (zone_id == -1 && assignment.access_node_id_to_zone_id_map.find(node_id) != assignment.access_node_id_to_zone_id_map.end())
				{
					zone_id = assignment.access_node_id_to_zone_id_map[node_id];
				}

				if (zone_id >= 1)
				{
					zone_id_to_analysis_district_id_mapping[zone_id] = analysis_district_id;

					node.zone_org_id = zone_id;  // this note here, we use zone_org_id to ensure we will only have super centriods with zone id positive.
					if (zone_id >= 1)
						node.is_activity_node = 1;  // from zone

					string str_mode_type;
					parser.GetValueByFieldName("mode_type", str_mode_type, false); //step 2 for adding access links: read mode_type for adding access links

					if (str_mode_type.size() > 0 && assignment.mode_type_2_seqno_mapping.find(str_mode_type) != assignment.mode_type_2_seqno_mapping.end())
					{
						node.mode_type_str = str_mode_type;
						node.mode_type_no = assignment.mode_type_2_seqno_mapping[str_mode_type];
						multmodal_activity_node_count++;
					}
				}


				int subarea_id = -1;
				parser.GetValueByFieldName("subarea_id", subarea_id, false);
				node.subarea_id = subarea_id;
				// this is an activity node // we do not allow zone id of zero
				if (zone_id >= 1)
				{
					// for physcial nodes because only centriod can have valid zone_id.
					node.zone_org_id = zone_id;
					if (zone_id_mapping.find(zone_id) == zone_id_mapping.end())
					{
						//create zone
						zone_id_mapping[zone_id] = node_id;
					}


				}
				if (zone_id >= 1)
				{

					assignment.node_seq_no_2_zone_id_mapping[internal_node_seq_no] = zone_id;
				}

				/*node.x = x;
				node.y = y;*/
				internal_node_seq_no++;

				// push it to the global node vector
				g_node_vector.push_back(node);


				assignment.g_number_of_nodes++;

				if (assignment.g_number_of_nodes % 5000 == 0)
				{
					dtalog.output() << "[STATUS INFO] reading " << assignment.g_number_of_nodes << " nodes.. " << '\n';
					g_DTA_log_file << "[STATUS INFO] reading " << assignment.g_number_of_nodes << " nodes.. " << '\n';
				}
			}



			dtalog.output() << "[DATA INFO] number of nodes = " << assignment.g_number_of_nodes << '\n';
			g_DTA_log_file << "[DATA INFO] number of nodes = " << assignment.g_number_of_nodes << '\n';
			//dtalog.output() << "[DATA INFO] number of multimodal activity nodes = " << multmodal_activity_node_count << '\n';
			//g_DTA_log_file << "[DATA INFO] number of multimodal activity nodes = " << multmodal_activity_node_count << '\n';
			//dtalog.output() << "[NOTE] One can add mode_type in node.csv to denote transit stations as part of efforts for modeling multmodal activities" << '\n';
			//g_DTA_log_file << "[NOTE] One can add mode_type in node.csv to denote transit stations as part of efforts for modeling multmodal activities" << '\n';


			// fprintf(g_pFileOutputLog, "number of nodes =,%d\n", assignment.g_number_of_nodes);
			parser.CloseCSVFile();
		}
	}



	int debug_line_count = 0;


	/// <summary>  mappping node to zone
	// hanlding multimodal access link: stage 1
	//step 3 for adding access links: there is node type restriction defined in agent type section of settings.csv
	for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at) // first loop for each agent type
	{
		if (assignment.g_ModeTypeVector[at].access_node_type.size() > 0)  // for each multmodal agent type
		{
			// find the closest zone id

			if (debug_line_count <= 20)
			{

				dtalog.output() << "[DATA INFO] multimodal access link generation condition 1: agent type " << assignment.g_ModeTypeVector[at].mode_type.c_str() << " has access node type" << assignment.g_ModeTypeVector[at].access_node_type.size() << '\n';
				g_DTA_log_file << "[DATA INFO] multimodal access link generation condition 1: agent type " << assignment.g_ModeTypeVector[at].mode_type.c_str() << " has access node type" << assignment.g_ModeTypeVector[at].access_node_type.size() << '\n';
				// zone without multimodal access
				debug_line_count++;
			}

			for (int a_k = 0; a_k < g_node_vector.size(); a_k++)
			{
				if (g_node_vector[a_k].is_activity_node == 1 && g_node_vector[a_k].mode_type_no == at) //second loop for mode_specific activity node
				{

					int zone_id = g_node_vector[a_k].zone_org_id;
					int zone_seq_no = zone_seq_no = assignment.g_zoneid_to_zone_seq_no_mapping[zone_id];

					if (debug_line_count <= 20)
					{

						dtalog.output() << "[DATA INFO] multimodal access link generation condition 2: agent type no = " << at << " for node no. " << a_k << "as activity node with zone_id >=1" << '\n';
						g_DTA_log_file << "[DATA INFO] multimodal access link generation condition 2: agent type no = " << at << " for node no. " << a_k << "as activity node with zone_id >=1" << '\n';
						// zone without multimodal access
						debug_line_count++;
					}


					// stage 2:  // min_distance
					double min_distance = 9999999;
					int min_distance_node_seq_no = -1;

					// stage 1:  // preferreed distance range
					double min_distance_within_range = 9999999;
					int min_distance_node_id_within_range = -1;

					std::vector<int> access_node_seq_vector;
					std::vector<float> access_node_distance_vector;

					for (int i = 0; i < g_node_vector.size(); i++)
					{
						if (g_node_vector[i].node_type.size() > 0)  // stop or station  //third loop for each stop or station node
						{

							if (assignment.g_ModeTypeVector[at].access_node_type.find(g_node_vector[i].node_type) != string::npos)  // check allowed access code
							{

								double zone_x = g_node_vector[a_k].x;
								double zone_y = g_node_vector[a_k].y;

								//test                                double near_by_distance_1 = g_calculate_p2p_distance_in_meter_from_latitude_longitude(-77.429293, 39.697895, -77.339847, 38.947676);

								double distance = g_calculate_p2p_distance_in_meter_from_latitude_longitude(zone_x, zone_y, g_node_vector[i].x, g_node_vector[i].y);
								// calculate the distance

								if (distance < min_distance)
								{
									min_distance = distance;
									min_distance_node_seq_no = i;
								}

								if (distance >= assignment.g_ModeTypeVector[at].access_distance_lb && distance <= assignment.g_ModeTypeVector[at].access_distance_ub)  // check the range
								{
									min_distance_within_range = distance;
									min_distance_node_id_within_range = i;
									access_node_seq_vector.push_back(i);
									access_node_distance_vector.push_back(distance);
								}
							}

						}
					}  // scan for all nodes


					// check access node vector for each pair of zone and agent type
					//
					if (access_node_seq_vector.size() > 0)  // preferred: access link within the range
					{
						float distance_k_cut_off_value = 99999;

						if (access_node_distance_vector.size() > assignment.g_ModeTypeVector[at].acecss_link_k)
						{

							std::vector<float> access_node_distance_vector_temp;
							access_node_distance_vector_temp = access_node_distance_vector;
							std::sort(access_node_distance_vector_temp.begin(), access_node_distance_vector_temp.end());

							distance_k_cut_off_value = access_node_distance_vector_temp[max(0, assignment.g_ModeTypeVector[at].acecss_link_k - 1)];
							//distance_k can be dynamically determined based on the density of stops and stations at different areas, e.g.CBM vs. rual area
						}

						for (int an = 0; an < access_node_seq_vector.size(); an++)
						{
							if (access_node_distance_vector[an] < distance_k_cut_off_value)  // within the shortest k ranage
							{
								g_add_new_access_link(a_k, access_node_seq_vector[an], access_node_distance_vector[an], at, -1);
								//incoming connector from station to activity centers
								g_add_new_access_link(access_node_seq_vector[an], a_k, access_node_distance_vector[an], at, -1);
								assignment.g_ModeTypeVector[at].zone_id_cover_map[zone_id] = true;
							}
						}

					}
					else if (min_distance_node_seq_no >= 0 && min_distance < assignment.g_ModeTypeVector[at].access_distance_ub)  // no node in the  preferred range, just use any feasible node with minimum distance by default
					{
						g_add_new_access_link(a_k, min_distance_node_seq_no, min_distance, at, -1);
						g_add_new_access_link(min_distance_node_seq_no, a_k, min_distance, at, -1);
						assignment.g_ModeTypeVector[at].zone_id_cover_map[zone_id] = true;

					}
					else {

						//                        dtalog.output() << " zone" << g_node_vector[a_k].zone_org_id << " with agent type = " << assignment.g_ModeTypeVector[at].mode_type.c_str() << " has no access to stop or station" << '\n';
						//                        g_DTA_log_file << " zone" << g_node_vector[a_k].zone_org_id << " with agent type = " << assignment.g_ModeTypeVector[at].mode_type.c_str() << " has no access to stop or station" << '\n';
												// zone without multimodal access
					}

				}  // for each zone

			}
		}
	}// for each agent type

	//g_InfoZoneMapping(assignment);

	// initialize zone vector
	dtalog.output() << "[PROCESS INFO] Step 1.5: Initializing O-D zone vector..." << '\n';
	g_DTA_log_file << "[PROCESS INFO] Step 1.5: Initializing O-D zone vector..." << '\n';

	int connector_link_warning_log_count = 0; 
	std::map<int, int>::iterator it;
	// creating zone centriod
	for (it = zone_id_mapping.begin(); it != zone_id_mapping.end(); ++it)
	{
		COZone ozone;

		// for each zone, we have to also create centriod
		ozone.zone_id = it->first;  // zone_id


		ozone.zone_seq_no = g_zone_vector.size();
		//ozone.obs_production = zone_id_production[it->first];
		//ozone.obs_attraction = zone_id_attraction[it->first];

		int zone_node_id = it->second; 

		if (assignment.g_node_id_to_seq_no_map.find(zone_node_id) != assignment.g_node_id_to_seq_no_map.end())
		{
			CNode node;
			node = g_node_vector[assignment.g_node_id_to_seq_no_map[zone_node_id]];
			ozone.cell_x = node.x;
			ozone.cell_y = node.y;

			assignment.zone_id_X_mapping[it->first] = node.x;
			assignment.zone_id_Y_mapping[it->first] = node.y;

		}

		if(zone_id_production[it->first] > 0)
		{
		ozone.gravity_production = zone_id_production[it->first];
		}
		if(zone_id_attraction[it->first] > 0)
		{
		ozone.gravity_attraction = zone_id_attraction[it->first];
		}

		assignment.g_zoneid_to_zone_seq_no_mapping[ozone.zone_id] = ozone.zone_seq_no;  // create the zone id to zone se&zq no mapping
		assignment.g_zoneid_to_zone_sindex_no_mapping[ozone.zone_id] = ozone.zone_seq_no;  // create the zone id to zone se&zq no mapping

		if (zone_id_to_analysis_district_id_mapping[ozone.zone_id] + 1 > assignment.g_number_of_analysis_districts)
			assignment.g_number_of_analysis_districts = zone_id_to_analysis_district_id_mapping[ozone.zone_id] + 1;

		assignment.g_zone_seq_no_to_analysis_distrct_id_mapping[ozone.zone_seq_no] = zone_id_to_analysis_district_id_mapping[ozone.zone_id];
		ozone.analysis_district_index = zone_id_to_analysis_district_id_mapping[ozone.zone_id];

		// create a centriod
		CNode node;
		// very large number as a special id
		node.node_id = -1 * ozone.zone_id;
		node.node_seq_no = g_node_vector.size();
		assignment.g_node_id_to_seq_no_map[node.node_id] = node.node_seq_no;
		node.zone_id = ozone.zone_id;
		node.x = ozone.cell_x;
		node.y = ozone.cell_y;
		// push it to the global node vector
		g_node_vector.push_back(node);
		assignment.g_number_of_nodes++;

		ozone.node_seq_no = node.node_seq_no;
		// this should be the only one place that defines this mapping
		assignment.zone_id_to_centriod_node_id_mapping[ozone.zone_id] = node.node_id;
		// add element into vector
		g_zone_vector.push_back(ozone);
	}


	dtalog.output() << "[DATA INFO] number of zones = " << g_zone_vector.size() << '\n';
	g_DTA_log_file << "[DATA INFO] number of zones = " << g_zone_vector.size() << '\n';

	g_zone_to_access(assignment);  // only under zone2connector mode
//	g_OutputModelFiles(3);  // node


	int demand_data_mode = g_detect_if_demand_data_provided(assignment);


	if (demand_data_mode == 0)  // if the demand files are not provided
	{
		dtalog.output() << "[PROCESS INFO] The demand files are not provided, proceeds to trip generation mode." << '\n';
		g_DTA_log_file << "[PROCESS INFO] The demand files are not provided, proceeds to trip generation mode." << '\n';

		g_demand_file_generation(assignment);
	}


	// MRM: stage 3: read link layers
	// step 4: read link file

	CDTACSVParser parser_link;

	int link_type_warning_count = 0;
	int link_penalty_info_count = 0;

	int length_missing_error = 0; 
	int free_speed_missing_error = 0;
	int capacity_missing_error = 0; 
	bool length_in_km_waring = false;
	dtalog.output() << "[PROCESS INFO] Step 1.6: Reading link data in link.csv... " << '\n';
	g_DTA_log_file << "[PROCESS INFO] Step 1.6: Reading link data in link.csv... " << '\n';

	int toll_message_count = 0;
	std::map<int, int> missing_link_type_mapping;

	for (int layer_no = 0; layer_no <= max_layer_flag; layer_no++)
	{
//		dtalog.output() << "link layer= " << layer_no << '\n';
//		g_DTA_log_file << "link layer= " << layer_no << '\n';


		string file_name = "link.csv";


		if (parser_link.OpenCSVFile(file_name.c_str(), true))
		{

			int line_no = 0; 

			while (parser_link.ReadRecord())  // if this line contains [] mark, then we will also read field headers.
			{

				// create a link object  // this should be place inside
				CLink link;

				string link_type_name_str;
				parser_link.GetValueByFieldName("link_type_name", link_type_name_str, false);

				long from_node_id = -1;
				if (!parser_link.GetValueByFieldName("from_node_id", from_node_id))
					continue;

				long to_node_id = -1;
				if (!parser_link.GetValueByFieldName("to_node_id", to_node_id))
					continue;

				string linkID;
				parser_link.GetValueByFieldName("link_id", linkID, false);
				// add the to node id into the outbound (adjacent) node list

				if (max_layer_flag == 1) // MRM mode
				{
					if (layer_no == 0)  // macro link layer
					{  // skip all links inside non-gate node

						if (assignment.g_node_id_to_MRM_subarea_mapping.find(from_node_id) != assignment.g_node_id_to_MRM_subarea_mapping.end()
							|| assignment.g_node_id_to_MRM_subarea_mapping.find(to_node_id) != assignment.g_node_id_to_MRM_subarea_mapping.end())
						{
							if (assignment.g_node_id_to_MRM_subarea_mapping[from_node_id] == 1
								&& assignment.g_node_id_to_MRM_subarea_mapping[to_node_id] == 1)  // both upstream and downstream inside
							{
								continue;

							}
							// keep all outside link or bridge link with a gate node inside the suabrea MRM

						}

					}


				}

				if (assignment.g_node_id_to_seq_no_map.find(from_node_id) == assignment.g_node_id_to_seq_no_map.end())
				{

					int MRM_subarea_mapping_flag = assignment.g_node_id_to_MRM_subarea_mapping[from_node_id];
					int MRM_subarea_mapping_flag_to = assignment.g_node_id_to_MRM_subarea_mapping[to_node_id];

					dtalog.output() << "[ERROR] from_node_id " << from_node_id << " in file link.csv is not defined in node.csv." << '\n';
					g_DTA_log_file << "[ERROR] from_node_id " << from_node_id << " in file link.csv is not defined in node.csv." << '\n';
					dtalog.output() << "   from_node_id mapping =" << MRM_subarea_mapping_flag << '\n';
					g_DTA_log_file << "   from_node_id mapping =" << MRM_subarea_mapping_flag << '\n';
					continue; //has not been defined
				}

				if (assignment.g_node_id_to_seq_no_map.find(to_node_id) == assignment.g_node_id_to_seq_no_map.end())
				{

					int MRM_subarea_mapping_flag = assignment.g_node_id_to_MRM_subarea_mapping[from_node_id];
					int MRM_subarea_mapping_flag_to = assignment.g_node_id_to_MRM_subarea_mapping[to_node_id];

					dtalog.output() << "[ERROR] to_node_id " << to_node_id << " in file link.csv is not defined in node.csv." << '\n';
					g_DTA_log_file << "[ERROR] to_node_id " << to_node_id << " in file link.csv is not defined in node.csv." << '\n';
					continue; //has not been defined
				}

				//if (assignment.g_link_id_map.find(linkID) != assignment.g_link_id_map.end())
				//    dtalog.output() << "Error: link_id " << linkID.c_str() << " has been defined more than once. Please check link.csv." << '\n';
				//    g_DTA_log_file << "Error: link_id " << linkID.c_str() << " has been defined more than once. Please check link.csv." << '\n';

				int internal_from_node_seq_no = assignment.g_node_id_to_seq_no_map[from_node_id];  // map external node number to internal node seq no.
				int internal_to_node_seq_no = assignment.g_node_id_to_seq_no_map[to_node_id];



				link.from_node_seq_no = internal_from_node_seq_no;
				link.to_node_seq_no = internal_to_node_seq_no;
				link.link_seq_no = assignment.g_number_of_links;
				link.to_node_seq_no = internal_to_node_seq_no;
				link.layer_no = layer_no;
				link.link_id = linkID;
				link.subarea_id = -1;

				if (g_node_vector[link.from_node_seq_no].subarea_id >= 1 || g_node_vector[link.to_node_seq_no].subarea_id >= 1)
				{
					link.subarea_id = g_node_vector[link.from_node_seq_no].subarea_id;
				}

				assignment.g_link_id_map[link.link_id] = 1;

				string movement_str;
				parser_link.GetValueByFieldName("mvmt_txt_id", movement_str, false);
				int cell_type = -1;
				if (parser_link.GetValueByFieldName("cell_type", cell_type, false) == true)
					link.cell_type = cell_type;

				int meso_link_id = -1;
				parser_link.GetValueByFieldName("meso_link_id", meso_link_id, false);

				if (meso_link_id >= 0)
					link.meso_link_id = meso_link_id;

				parser_link.GetValueByFieldName("geometry", link.geometry, false);

				if (link.geometry.size() == 0)
				{
					std::ostringstream out;

					out << "LINESTRING (";

					out << g_node_vector[link.from_node_seq_no].x  << " " << g_node_vector[link.from_node_seq_no].y  << ", ";
					out << g_node_vector[link.to_node_seq_no].x << " " << g_node_vector[link.to_node_seq_no].y;


					out << ")";

					link.geometry = out.str();
				}

				parser_link.GetValueByFieldName("link_special_flag", link.link_specifical_flag_str, false);

				link.tmc_corridor_name = "network_wide";
				parser_link.GetValueByFieldName("tmc_corridor_name", link.tmc_corridor_name, false);
				parser_link.GetValueByFieldName("link_type_name", link.link_type_name, false);

				parser_link.GetValueByFieldName("link_type_code", link.link_type_code, false);

				// and valid
				if (movement_str.size() > 0)
				{
					int main_node_id = -1;


					link.mvmt_txt_id = movement_str;
					link.main_node_id = main_node_id;
				}


				int default_link_type = assignment.g_first_link_type;
				char link_type_field_name[50];

				parser_link.GetValueByFieldName("link_type", default_link_type, false);
				for (int sii = 0; sii < assignment.g_DTA_scenario_vector.size(); sii++)
				{
					int scenario_index = assignment.g_DTA_scenario_vector[sii].scenario_index;
					sprintf(link_type_field_name, "link_type_s%d", scenario_index);

					link.link_type_si[scenario_index] = default_link_type;
					if (parser_link.GetValueByFieldName(link_type_field_name, link.link_type_si[scenario_index], false) == false && line_no == 0)
					{
						dtalog.output() << "[WARNING] Field " << link_type_field_name << "  in link.csv is not defined. The default value in the field link_type in link.csv is used." << '\n';
						g_DTA_log_file << "[WARNING] Field " << link_type_field_name << "  in link.csv is not defined. The default value in the field link_type in link.csv is used." << '\n';
						continue;

					}

					if (assignment.g_LinkTypeMap.find(link.link_type_si[scenario_index]) == assignment.g_LinkTypeMap.end())
					{
						missing_link_type_mapping[link.link_type_si[scenario_index]] = 1;

						link.link_type_si[scenario_index] = assignment.g_first_link_type;
					}


				}




				if (assignment.g_LinkTypeMap[link.link_type_si[0]].type_code == "c")  // suggestion: we can move "c" as "connector" in allowed_uses
				{
					if (g_node_vector[internal_from_node_seq_no].zone_org_id >= 1)
					{
						int zone_org_id = g_node_vector[internal_from_node_seq_no].zone_org_id;
						if (assignment.g_zoneid_to_zone_seq_no_mapping.find(zone_org_id) != assignment.g_zoneid_to_zone_seq_no_mapping.end())
							link.zone_seq_no_for_outgoing_connector = assignment.g_zoneid_to_zone_seq_no_mapping[zone_org_id];
					}
					else
					{
						if (g_node_vector[internal_to_node_seq_no].zone_org_id < 0 && connector_link_warning_log_count < 2) // this connector's upstream and downstream nodes are not associated with any zone id
						{
							//						link.zone_seq_no_for_outgoing_connector = 999999; // only for the purpose of being skipped in the routing
							dtalog.output() << "[WARNING] The upstream node of the connector link from node " << from_node_id << "to node " << to_node_id << " in link.csv does not have a defined zone ID. Please ensure that the zone ID is defined for the upstream node in order to accurately assign the traffic from the respective zone through this connector link." << '\n';
							g_DTA_log_file << "[WARNING] The upstream node of the connector link from node " << from_node_id << "to node " << to_node_id << " in link.csv does not have a defined zone ID. Please ensure that the zone ID is defined for the upstream node in order to accurately assign the traffic from the respective zone through this connector link." << '\n';
							connector_link_warning_log_count++;
						}
					}

				}

				char lanes_scenario_field_name[50];

				int default_number_of_lanes = 1;
				parser_link.GetValueByFieldName("lanes", default_number_of_lanes, false);
				link.number_of_lanes_si[0] = default_number_of_lanes;

				for (int sii = 0; sii < assignment.g_DTA_scenario_vector.size(); sii++)
				{
					int scenario_index = assignment.g_DTA_scenario_vector[sii].scenario_index;
					sprintf(lanes_scenario_field_name, "lanes_s%d", scenario_index);
					double number_of_lanes_value = default_number_of_lanes;  // default value
					parser_link.GetValueByFieldName(lanes_scenario_field_name, number_of_lanes_value, false);

					link.number_of_lanes_si[scenario_index] = number_of_lanes_value;
				}

				double length_in_meter = 1000.0; // km or mile
				double free_speed = 60.0;
				double lane_capacity = 2000;
				// first step, get the initial value based on link type for speed and capacity
				for (int at = 0; at < assignment.g_ModeTypeVector.size(); at++)
				{
					for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
					{
						link.VDF_period[tau].capacity_at[at] = assignment.g_LinkTypeMap[link.link_type_si[0]].capacity_at[at];
						link.VDF_period[tau].free_speed_at[at] = assignment.g_LinkTypeMap[link.link_type_si[0]].free_speed_at[at];
						link.VDF_period[tau].lanes_mode_type[at] = link.number_of_lanes_si[0];
					}

				}

				double cutoff_speed = 1.0;
				double k_jam = assignment.g_LinkTypeMap[link.link_type_si[0]].k_jam;
				double bwtt_speed = 12;  //miles per hour

				double length_original_input = 0; 

				if (parser_link.GetValueByFieldName("length", length_original_input, false, false) == false)
				{
					if(length_missing_error==0)
					{
					dtalog.output() << "[ERORR] Field length in link.csv is missing." << '\n';
					g_DTA_log_file << "[ERORR] Field length in link.csv is missing." << '\n';
					length_in_meter = 999.0;
					length_missing_error++; 
					}
				}

				if (assignment.g_length_unit_flag == 1)  // mile;
					length_in_meter = length_original_input * 1609.0; // convert from mile to meter
				if (assignment.g_length_unit_flag == 2)  // km;
					length_in_meter = length_original_input *1000.0; // convert from km to meter
				if (assignment.g_length_unit_flag == 0)  // meter;
					length_in_meter = length_original_input; // no conversion

				parser_link.GetValueByFieldName("FT", link.FT, false, true);
				parser_link.GetValueByFieldName("AT", link.AT, false, true);
				parser_link.GetValueByFieldName("vdf_code", link.vdf_code, false);

				if (length_in_km_waring == false && length_in_meter < 0.1)
				{
					dtalog.output() << "[WARNING] link link_length =" << length_in_meter << " in link.csv for link " << from_node_id << "->" << to_node_id << ". Please ensure the unit of the link link_distance_VDF is meter." << '\n';
					g_DTA_log_file << "[WARNING] link link_length =" << length_in_meter << " in link.csv for link " << from_node_id << "->" << to_node_id << ". Please ensure the unit of the link link_distance_VDF is meter." << '\n';
					// link.link_type has been taken care by its default constructor
					length_in_km_waring = true;
				}

				if (length_in_meter < 1)
				{
					length_in_meter = 1;  // minimum link_distance_VDF
				}

				// second step, we read the link-specific value (only for based mode)
				if (parser_link.GetValueByFieldName("free_speed", free_speed, false, false) == false)
				{

					if (free_speed_missing_error == 0)
					{
						dtalog.output() << "[ERORR] Field free_speed in link.csv is missing." << '\n';
						g_DTA_log_file << "[ERORR] Field free_speed in link.csv is missing." << '\n';
						free_speed = 60;
						free_speed_missing_error++; 
					}
				}
				if (assignment.g_speed_unit_flag == 1)  // mph;
					free_speed = free_speed / 1.609; // convert from mile per hour to km per hour

				parser_link.GetValueByFieldName("capacity", lane_capacity, false);  // optional for capacity

								// second step, we read the link-specific value (only for based mode)
				if (parser_link.GetValueByFieldName("capacity", lane_capacity, false, false) == false)
				{

					if (capacity_missing_error == 0)
					{
						dtalog.output() << "[ERORR] Field capacity in link.csv is missing." << '\n';
						g_DTA_log_file << "[ERORR] Field capacity in link.csv is missing." << '\n';
						lane_capacity = 1800;
						capacity_missing_error++;
					}
				}
				double free_speed_value;
				//speed_unit: km/ph us_customary  si1
				// mile per h our ->



				cutoff_speed = free_speed * 0.85; //default;
				parser_link.GetValueByFieldName("cutoff_speed", cutoff_speed, false);


				free_speed = max(0.1, free_speed);

				link.free_speed = free_speed;
				link.lane_capacity = lane_capacity;
				int at_base = 0;

				// third step, we write the link-specific value (only for based mode)
				for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
				{
					link.VDF_period[tau].capacity_at[at_base] = lane_capacity;
					link.VDF_period[tau].free_speed_at[at_base] = free_speed;
				}

				link.v_congestion_cutoff = cutoff_speed;



				parser_link.GetValueByFieldName("saturation_flow_rate", link.saturation_flow_rate, false);


				link.free_flow_travel_time_in_min = length_in_meter / 1000.0 / free_speed * 60;  // link_distance_VDF in meter
				float fftt_in_sec = link.free_flow_travel_time_in_min * 60;  // link_distance_VDF in meter

				link.length_in_meter = length_in_meter;
				link.link_distance_VDF = length_in_meter / 1000.0;
				// to do: settings.csv
				// meter -> value /1000
				// km -> value
				// mile -> /1.609

				link.link_distance_km = length_in_meter / 1000.0;
				link.link_distance_mile = length_in_meter / 1609.0;

				link.traffic_flow_code = assignment.g_LinkTypeMap[link.link_type_si[0]].traffic_flow_code;

				//spatial queue and kinematic wave
				link.spatial_capacity_in_vehicles = max(1.0, link.link_distance_VDF * link.number_of_lanes_si[0] * k_jam);

				// kinematic wave
				if (link.traffic_flow_code == 3)
					link.BWTT_in_simulation_interval = link.link_distance_VDF / bwtt_speed * 3600 / number_of_seconds_per_interval;

				link.vdf_type = assignment.g_LinkTypeMap[link.link_type_si[0]].vdf_type;
				link.kjam = assignment.g_LinkTypeMap[link.link_type_si[0]].k_jam;
				char CSV_field_name[50];

				// based on link type
				for (int at = 0; at < assignment.g_ModeTypeVector.size(); at++)
				{
					for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
					{
						link.VDF_period[tau].capacity_at[at] = assignment.g_LinkTypeMap[link.link_type_si[0]].capacity_at[at];
						link.VDF_period[tau].free_speed_at[at] = assignment.g_LinkTypeMap[link.link_type_si[0]].free_speed_at[at];
						link.VDF_period[tau].FFTT_at[at] = link.link_distance_VDF / max(0.0001, link.VDF_period[tau].free_speed_at[at]) * 60.0;  // 60.0 for 60 min per hour
					}

				}

				// if we have link_specific values in the master link.csv, we will use it as the final value for the driving mode 
				for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
				{
					if (lane_capacity > 0)
					{
						link.VDF_period[tau].capacity_diff_link_specific = lane_capacity - link.VDF_period[tau].capacity_at[at_base];
						link.VDF_period[tau].capacity_at[at_base] = lane_capacity;
					}

					if (free_speed > 0)
					{
						link.VDF_period[tau].free_speed_diff_link_specific = free_speed - link.VDF_period[tau].free_speed_at[at_base];
						link.VDF_period[tau].free_speed_at[at_base] = free_speed;
						link.VDF_period[tau].FFTT_at[at_base] = link.link_distance_VDF / max(0.0001, free_speed) * 60.0;  // 60.0 for 60 min per hour

					}
				}


				// reading for VDF related functions
				// step 1 read type


					//data initialization

				for (int time_index = 0; time_index < MAX_TIMEINTERVAL_PerDay; time_index++)
				{
					link.model_speed[time_index] = free_speed;
					link.est_volume_per_hour_per_lane[time_index] = 0;
					link.est_avg_waiting_time_in_min[time_index] = 0;
					link.est_queue_length_per_lane[time_index] = 0;
				}


				for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
				{
					//setup default values
					link.VDF_period[tau].vdf_type = assignment.g_LinkTypeMap[link.link_type_si[0]].vdf_type;
					link.VDF_period[tau].lane_based_ultimate_hourly_capacity = lane_capacity;
					link.VDF_period[tau].nlanes = link.number_of_lanes_si[0];

					int at_base = 0;
					link.VDF_period[tau].vf = link.free_speed;
					link.VDF_period[tau].v_congestion_cutoff = link.v_congestion_cutoff;
					link.VDF_period[tau].alpha = 0.15;
					link.VDF_period[tau].beta = 4;
					link.VDF_period[tau].preload = 0;
					parser_link.GetValueByFieldName("VDF_alpha", link.VDF_period[tau].alpha, false);
					parser_link.GetValueByFieldName("VDF_beta", link.VDF_period[tau].beta, false);

					for (int at = 0; at < assignment.g_ModeTypeVector.size(); at++)
					{
						link.VDF_period[tau].toll[at][assignment.active_scenario_index] = 0;
						link.VDF_period[tau].LR_price[at] = 0;
						link.VDF_period[tau].LR_RT_price[at] = 0;

						for (int si = 0; si < g_number_of_active_scenarios; si++)
						{
							link.VDF_period[tau].allowed_uses[si] = assignment.g_LinkTypeMap[link.link_type_si[si]].allow_uses_period[tau];
						}
					}

					link.VDF_period[tau].starting_time_in_hour = assignment.g_DemandPeriodVector[tau].starting_time_slot_no * MIN_PER_TIMESLOT / 60.0;
					link.VDF_period[tau].ending_time_in_hour = assignment.g_DemandPeriodVector[tau].ending_time_slot_no * MIN_PER_TIMESLOT / 60.0;
					link.VDF_period[tau].L = assignment.g_DemandPeriodVector[tau].time_period_in_hour;
					link.VDF_period[tau].t2 = assignment.g_DemandPeriodVector[tau].t2_peak_in_hour;

					sprintf(CSV_field_name, "restricted_turn_nodes_p%d", tau + 1);

					string restricted_turn_nodes_str;
					parser_link.GetValueByFieldName(CSV_field_name, restricted_turn_nodes_str, false);

					//"restricted_turn_nodes: This field represents the nodes within the network where left turns and possible U-turns are prevented, to facilitate better traffic flow and manage congestion."
					std::vector<int> restricted_turn_nodes;

					if (restricted_turn_nodes_str.size() > 0)
					{
						g_ParserIntSequence(restricted_turn_nodes_str, restricted_turn_nodes);

						for (int node_i = 0; node_i < restricted_turn_nodes.size(); node_i++)
						{

							if (assignment.g_node_id_to_seq_no_map.find(restricted_turn_nodes[node_i]) == assignment.g_node_id_to_seq_no_map.end())
							{

								dtalog.output() << "[ERROR] restricted_turn_nodes: " << restricted_turn_nodes[node_i] << " is not defined in node.csv, for link from node " << from_node_id << " to node " << to_node_id << " in the 'link.csv' file " << ".\n";
								g_DTA_log_file << "[ERROR] restricted_turn_nodes: " << restricted_turn_nodes[node_i] << " is not defined in node.csv, for link from node " << from_node_id << " to node " << to_node_id << " in the 'link.csv' file " << ".\n";
								continue; 
							}
							int node_no = assignment.g_node_id_to_seq_no_map[restricted_turn_nodes[node_i]];
							g_node_vector[link.to_node_seq_no].prohibited_movement_size++;  // we consider the movement mark as link.to_node_seq_no 
							link.VDF_period[tau].restricted_turn_nodes_str = restricted_turn_nodes_str;
							link.VDF_period[tau].restricted_turn_nodes_map[node_no] = 1;

						}


					}

					double value_ref = -1;
					sprintf(CSV_field_name, "ref_volume_p%d", tau + 1);
					parser_link.GetValueByFieldName(CSV_field_name, value_ref, false, false);
					if (value_ref > 1)
						link.VDF_period[tau].ref_link_volume = value_ref;

					sprintf(CSV_field_name, "VDF_preload_p%d", tau + 1);
					parser_link.GetValueByFieldName(CSV_field_name, link.VDF_period[tau].preload, false, false);

					for (int at = 0; at < assignment.g_ModeTypeVector.size(); at++)
					{

						//sprintf(CSV_field_name, "VDF_toll%s%d", assignment.g_ModeTypeVector[at].mode_type.c_str(), tau);
						//parser_link.GetValueByFieldName(CSV_field_name, link.VDF_period[tau].toll[at], false, false);

						//if (link.VDF_period[tau].toll[at] > 0.001 && toll_message_count < 10)
						//{
						//	dtalog.output() << "[DATA INFO] link " << from_node_id << "->" << to_node_id << " has a toll of " << link.VDF_period[tau].toll[at] << " for agent type "
						//	g_DTA_log_file << "[DATA INFO] link " << from_node_id << "->" << to_node_id << " has a toll of " << link.VDF_period[tau].toll[at] << " for agent type "
						//		<< assignment.g_ModeTypeVector[at].mode_type.c_str() << " at demand period " << at << '\n';
						//	toll_message_count++;
						//}
						//sprintf(CSV_field_name, "VDF_penalty%d", at);
						/*					parser_link.GetValueByFieldName(CSV_field_name, link.VDF_period[tau].penalty, false, false);*/

					}

					if (link.cell_type >= 2) // micro lane-changing arc
					{
						// additinonal min: 1 second 1/60.0 min
						link.VDF_period[tau].penalty += 1.0 / 60.0;
					}

					parser_link.GetValueByFieldName("cycle_length", link.VDF_period[tau].cycle_length, false, false);

					if (link.VDF_period[tau].cycle_length >= 1)
					{
						link.timing_arc_flag = true;

						parser_link.GetValueByFieldName("start_green_time", link.VDF_period[tau].start_green_time);
						parser_link.GetValueByFieldName("end_green_time", link.VDF_period[tau].end_green_time);

						link.VDF_period[tau].effective_green_time = link.VDF_period[tau].end_green_time - link.VDF_period[tau].start_green_time;

						if (link.VDF_period[tau].effective_green_time < 0)
							link.VDF_period[tau].effective_green_time = link.VDF_period[tau].cycle_length;

						link.VDF_period[tau].red_time = max(1.0f, link.VDF_period[tau].cycle_length - link.VDF_period[tau].effective_green_time);
						parser_link.GetValueByFieldName("red_time", link.VDF_period[tau].red_time, false);
						parser_link.GetValueByFieldName("green_time", link.VDF_period[tau].effective_green_time, false);

						if (link.saturation_flow_rate > 1000)  // protect the data attributes to be reasonable
						{
							link.VDF_period[tau].saturation_flow_rate = link.saturation_flow_rate;
						}
					}

				}

				// for each period

				float default_cap = 1000;
				float default_BaseTT = 1;


				//link.m_OutflowNumLanes = number_of_lanes;//visum lane_cap is actually link_cap

				link.update_kc(free_speed);
				link.link_spatial_capacity = k_jam * link.number_of_lanes_si[0] * link.link_distance_VDF;
				// min // calculate link cost based link_distance_VDF and speed limit // later we should also read link_capacity, calculate delay

								//// TMC reading
				string tmc_code;
				parser_link.GetValueByFieldName("tmc", link.tmc_code, false);
				if (link.tmc_code.size() > 0)
				{
					parser_link.GetValueByFieldName("tmc_corridor_name", link.tmc_corridor_name, false);
					link.tmc_corridor_id = 1;
					link.tmc_road_sequence = 1;
					parser_link.GetValueByFieldName("tmc_corridor_id", link.tmc_corridor_id, false);
					parser_link.GetValueByFieldName("tmc_road_sequence", link.tmc_road_sequence, false);
				}


				

				//int sequential_copying = 0;
				//
				//parser_link.GetValueByFieldName("sequential_copying", sequential_copying);

					// read optional penalty term
				for (int sii = 0; sii < assignment.g_DTA_scenario_vector.size(); sii++)
				{
					int scenario_index = assignment.g_DTA_scenario_vector[sii].scenario_index;

					for (int at = 0; at < assignment.g_ModeTypeVector.size(); at++)
					{
						char penalty_field_name[50];
						sprintf(penalty_field_name, "penalty_%s_s%d", assignment.g_ModeTypeVector[at].mode_type.c_str(), scenario_index);
						double penalty = 0;
						parser_link.GetValueByFieldName(penalty_field_name, penalty, false);


						if (fabs(penalty) > 0.0001)
						{
							if(link_penalty_info_count<3)
							{
								dtalog.output() << "[DATA INFO] Link penalty detected! For the link from node " << from_node_id << " to node " << to_node_id << " in the 'link.csv' file, the penalty field named '" << penalty_field_name << "' has a value of " << penalty << ".\n";
								g_DTA_log_file << "[DATA INFO] Link penalty detected! For the link from node " << from_node_id << " to node " << to_node_id << " in the 'link.csv' file, the penalty field named '" << penalty_field_name << "' has a value of " << penalty << ".\n";
								link_penalty_info_count++;
							}
							assignment.summary_file << "[DATA INFO] link penalty detected! ! link  " << from_node_id << " ->" << to_node_id << "  in link.csv has penalty_field_name = " << penalty_field_name << "of " << penalty << '\n';

							

							for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
							{

								link.penalty_si_at[tau][at][scenario_index] = penalty;
							}


						}
					}


				}
				int dir_flag = 1;
				parser_link.GetValueByFieldName("dir_flag", dir_flag, false);

				int link_code_start = 1;
				int link_code_end = 1;

				if (dir_flag == -1) // reversed
				{
					link_code_start = 2; link_code_end = 2;
				}

				if (dir_flag == 0) // two-directional link
				{
					link_code_start = 1; link_code_end = 2;
				}

				if (dir_flag == 2) // bi-directional link  for single rail track
				{
					link_code_start = 1; link_code_end = 2;
				}

				int layer_no = 0;

				for (int link_code = link_code_start; link_code <= link_code_end; link_code++)
				{
					int internal_from_node_seq_no_dir = -1;
					int internal_to_node_seq_no_dir = -1;
					int link_seq_no = link.link_seq_no;
					if (link_code == 1)
					{
						internal_from_node_seq_no_dir = internal_from_node_seq_no;
						internal_to_node_seq_no_dir = internal_to_node_seq_no;
						link_seq_no = link.link_seq_no;
						link.AB_flag = 1;

						if (dir_flag == 0 || dir_flag == 2)
						{
							link.BA_link_no = link.link_seq_no + 1;
						}
						else
							link.BA_link_no = -1;


					}
					else //link_code_start ==2
					{
						// switch the node sequence no
						int temp_node_no = internal_to_node_seq_no;
						internal_to_node_seq_no = internal_from_node_seq_no;
						internal_from_node_seq_no = temp_node_no;

						int temp_node_id = link.to_node_id;
						link.to_node_id = link.from_node_id;
						link.from_node_id = temp_node_id;


						link.AB_flag = -1;
						link.BA_link_no = -1;


						link.link_seq_no += 1;  // update link seq no by + 1 for BA link
						link_seq_no = link.link_seq_no;


						link.from_node_seq_no = internal_from_node_seq_no;
						link.to_node_seq_no = internal_to_node_seq_no;


					}

					g_node_vector[internal_from_node_seq_no].m_outgoing_link_seq_no_vector.push_back(link_seq_no);  // add this link to the corresponding node as part of outgoing node/link
					g_node_vector[internal_to_node_seq_no].m_incoming_link_seq_no_vector.push_back(link_seq_no);  // add this link to the corresponding node as part of outgoing node/link

					g_node_vector[internal_from_node_seq_no].m_to_node_seq_no_vector.push_back(link.to_node_seq_no);  // add this link to the corresponding node as part of outgoing node/link
					g_node_vector[internal_from_node_seq_no].m_to_node_2_link_seq_no_map[link.to_node_seq_no] = link.link_seq_no;  // add this link to the corresponding node as part of outgoing node/link
					g_link_vector.push_back(link);


					assignment.g_number_of_links++;
				}

				if (assignment.g_number_of_links % 10000 == 0)
				{
					dtalog.output() << "[STATUS INFO] reading " << assignment.g_number_of_links << " links.. " << '\n';
					g_DTA_log_file << "[STATUS INFO] reading " << assignment.g_number_of_links << " links.. " << '\n';
				}
				line_no++;
			}
			parser_link.CloseCSVFile();

		}
		else
		{
		dtalog.output() << "[ERROR] The critical input GMNS file 'link.csv' is missing. Please make sure the file is included in the appropriate directory." << '\n';
		g_DTA_log_file << "[ERROR] The critical input GMNS file 'link.csv' is missing. Please make sure the file is included in the appropriate directory." << '\n';
		g_program_stop();
		}

	}


	

	if (missing_link_type_mapping.size() > 0)
	{
		dtalog.output() << "[WARNING] The following link types in link.csv are not defined in link_type.csv. ";
		g_DTA_log_file << "[WARNING] The following link types in link.csv are not defined in link_type.csv. ";
			// Iterate over the map and print all values
			for (const auto& pair : missing_link_type_mapping) 
			{
				int value = pair.first;
				dtalog.output() << value << "; ";
				g_DTA_log_file << value << "; ";
			}

			dtalog.output() <<  '\n';
			g_DTA_log_file <<  '\n';
			dtalog.output() << "The default value of link type = " << assignment.g_first_link_type << " is used." << '\n';
			g_DTA_log_file << "The default value of link type = " << assignment.g_first_link_type << " is used." << '\n';
	}



	assignment.summary_file << "[PROCESS INFO] Step 1: read network node.csv, link.csv, zone.csv " << '\n';
	assignment.summary_file << ",# of nodes = ," << g_node_vector.size() << '\n';
	assignment.summary_file << ",# of links =," << g_link_vector.size() << '\n';
	assignment.summary_file << ",# of zones =," << g_zone_vector.size() << '\n';
	const int fieldWidth = 12;
	dtalog.output() << "[PROCESS INFO] Step 1: read network node.csv, link.csv, zone.csv " << '\n';
	g_DTA_log_file << "[PROCESS INFO] Step 1: read network node.csv, link.csv, zone.csv " << '\n';
	dtalog.output() << "[DATA INFO] " << std::setw(fieldWidth) << "# of nodes = " << g_node_vector.size() << '\n';
	g_DTA_log_file << "[DATA INFO] " << std::setw(fieldWidth) << "# of nodes = " << g_node_vector.size() << '\n';
	dtalog.output() << "[DATA INFO] " << std::setw(fieldWidth) << "# of links = " << g_link_vector.size() << '\n';
	g_DTA_log_file << "[DATA INFO] " << std::setw(fieldWidth) << "# of links = " << g_link_vector.size() << '\n';
	dtalog.output() << "[DATA INFO] " << std::setw(fieldWidth) << "# of zones = " << g_zone_vector.size() << '\n';
	g_DTA_log_file << "[DATA INFO] " << std::setw(fieldWidth) << "# of zones = " << g_zone_vector.size() << '\n';

	assignment.summary_file << "[PROCESS INFO] Step 1: read network node.csv, link.csv, zone.csv " << '\n';
	assignment.summary_file << ",# of nodes = ," << g_node_vector.size() << '\n';
	assignment.summary_file << ",# of links =," << g_link_vector.size() << '\n';
	assignment.summary_file << ",# of zones =," << g_zone_vector.size() << '\n';
	//assignment.summary_file << ",# of subarea links =," << bridge_link_count << '\n';
	//assignment.summary_file << ",# of subarea nodes =," << number_of_micro_gate_nodes << '\n';


	/// <summary>
	/// ///////////////////////////
	/// </summary>
	/// <param name="assignment"></param>
	assignment.summary_file << ",summary by road link type,link_type,link_type_name,# of links,avg_free_speed_mph,avg_free_speed_kmph,,total_length,total_capacity,avg_lane_capacity,avg_length_in_meter," << '\n';
	std::map<int, CLinkType>::iterator it_link_type;
	int count_zone_demand = 0;
	for (it_link_type = assignment.g_LinkTypeMap.begin(); it_link_type != assignment.g_LinkTypeMap.end(); ++it_link_type)
	{
		assignment.summary_file << ",," << it_link_type->first << "," << it_link_type->second.link_type_name.c_str() << ",";

		int link_count = 0;
		double total_speed = 0;
		double total_length = 0;
		double total_lane_capacity = 0;
		double total_link_capacity = 0;

		for (int i = 0; i < g_link_vector.size(); i++)
		{
			if (g_link_vector[i].link_type_si[assignment.active_scenario_index] >= 0 && g_link_vector[i].link_type_si[assignment.active_scenario_index] == it_link_type->first)
			{
				link_count++;
				total_speed += g_link_vector[i].free_speed;
				total_length += g_link_vector[i].length_in_meter * g_link_vector[i].number_of_lanes_si[assignment.active_scenario_index];
				total_lane_capacity += g_link_vector[i].lane_capacity;
				total_link_capacity += g_link_vector[i].lane_capacity * g_link_vector[i].number_of_lanes_si[assignment.active_scenario_index];
			}
		}
		assignment.summary_file << link_count << "," <<
			total_speed / max(1, link_count) << "," <<
			total_speed / max(1, link_count) * 1.60934 << "," <<
			total_length / 1000.0 << "," <<
			total_link_capacity << "," <<
			total_lane_capacity / max(1, link_count) << "," << total_length / max(1, link_count) << "," << '\n';


	}

	{
	
		const int fieldWidth = 12;

		dtalog.output() << "[PROCESS INFO] ";
		g_DTA_log_file << "[PROCESS INFO] ";
		dtalog.output() << std::setw(fieldWidth) << "link_type,";
		g_DTA_log_file << std::setw(fieldWidth) << "link_type,";
		dtalog.output() << std::setw(fieldWidth) << "link_name,";
		g_DTA_log_file << std::setw(fieldWidth) << "link_name,";
		dtalog.output() << std::setw(fieldWidth) << "# links,";
		g_DTA_log_file << std::setw(fieldWidth) << "# links,";
		dtalog.output() << std::setw(fieldWidth) << "avg_mph,";
		g_DTA_log_file << std::setw(fieldWidth) << "avg_mph,";
		dtalog.output() << std::setw(fieldWidth) << "avg_kmph,";
		g_DTA_log_file << std::setw(fieldWidth) << "avg_kmph,";
		dtalog.output() << std::setw(fieldWidth) << "total_len_km,";
		g_DTA_log_file << std::setw(fieldWidth) << "total_len_km,";
		dtalog.output() << std::setw(fieldWidth) << "total_cap,";
		g_DTA_log_file << std::setw(fieldWidth) << "total_cap,";
		dtalog.output() << std::setw(fieldWidth) << "avg_lane,";
		g_DTA_log_file << std::setw(fieldWidth) << "avg_lane,";
		dtalog.output() << std::setw(fieldWidth) << "avg_len," << '\n';
		g_DTA_log_file << std::setw(fieldWidth) << "avg_len," << '\n';

		std::map<int, CLinkType>::iterator it_link_type;
		int count_zone_demand = 0;

		for (it_link_type = assignment.g_LinkTypeMap.begin(); it_link_type != assignment.g_LinkTypeMap.end(); ++it_link_type)
		{
			assignment.summary_file << ",," << it_link_type->first << "," << it_link_type->second.link_type_name.c_str() << ",";

			int link_count = 0;
			double total_speed = 0;
			double total_length = 0;
			double total_lane_capacity = 0;
			double total_link_capacity = 0;

			for (int i = 0; i < g_link_vector.size(); i++)
			{
				if (g_link_vector[i].link_type_si[assignment.active_scenario_index] >= 0 && g_link_vector[i].link_type_si[assignment.active_scenario_index] == it_link_type->first)
				{
					link_count++;
					total_speed += g_link_vector[i].free_speed;
					total_length += g_link_vector[i].length_in_meter * g_link_vector[i].number_of_lanes_si[assignment.active_scenario_index];
					total_lane_capacity += g_link_vector[i].lane_capacity;
					total_link_capacity += g_link_vector[i].lane_capacity * g_link_vector[i].number_of_lanes_si[assignment.active_scenario_index];
				}
			}

			if (link_count > 0)
			{
				dtalog.output() << "[PROCESS INFO]";
				g_DTA_log_file << "[PROCESS INFO]";

				dtalog.output() << std::setw(fieldWidth) << it_link_type->first;
				g_DTA_log_file << std::setw(fieldWidth) << it_link_type->first;
				dtalog.output() << std::setw(fieldWidth) << it_link_type->second.link_type_name.c_str();
				g_DTA_log_file << std::setw(fieldWidth) << it_link_type->second.link_type_name.c_str();
				dtalog.output() << std::setw(fieldWidth) << link_count;
				g_DTA_log_file << std::setw(fieldWidth) << link_count;

				dtalog.output() << std::setw(fieldWidth) << (total_speed / max(1, link_count));
				g_DTA_log_file << std::setw(fieldWidth) << (total_speed / max(1, link_count));
				dtalog.output() << std::setw(fieldWidth) << (total_speed / max(1, link_count) * 1.60934);
				g_DTA_log_file << std::setw(fieldWidth) << (total_speed / max(1, link_count) * 1.60934);
				dtalog.output() << std::setw(fieldWidth) << (total_length / 1000.0);
				g_DTA_log_file << std::setw(fieldWidth) << (total_length / 1000.0);
				dtalog.output() << std::setw(fieldWidth) << total_link_capacity;
				g_DTA_log_file << std::setw(fieldWidth) << total_link_capacity;
				dtalog.output() << std::setw(fieldWidth) << (total_lane_capacity / max(1, link_count));
				g_DTA_log_file << std::setw(fieldWidth) << (total_lane_capacity / max(1, link_count));
				dtalog.output() << std::setw(fieldWidth) << (total_length / max(1, link_count)) << '\n';
				g_DTA_log_file << std::setw(fieldWidth) << (total_length / max(1, link_count)) << '\n';
			}
		}
	}


	if (dtalog.debug_level() == 2)
	{
		for (int i = 0; i < g_node_vector.size(); ++i)
		{
			if (g_node_vector[i].zone_org_id >= 1) // for each physical node
			{
				// we need to make sure we only create two way connectors between nodes and zones
				dtalog.output() << "[DATA INFO] node id= " << g_node_vector[i].node_id << " with zone id " << g_node_vector[i].zone_org_id << "and "
					<< g_node_vector[i].m_outgoing_link_seq_no_vector.size() << " outgoing links." << '\n';

				g_DTA_log_file << "[DATA INFO] node id= " << g_node_vector[i].node_id << " with zone id " << g_node_vector[i].zone_org_id << "and "
					<< g_node_vector[i].m_outgoing_link_seq_no_vector.size() << " outgoing links." << '\n';

				for (int j = 0; j < g_node_vector[i].m_outgoing_link_seq_no_vector.size(); ++j)
				{
					int link_seq_no = g_node_vector[i].m_outgoing_link_seq_no_vector[j];
					dtalog.output() << "[DATA INFO] outgoing node = " << g_node_vector[g_link_vector[link_seq_no].to_node_seq_no].node_id << '\n';
					g_DTA_log_file << "[DATA INFO] outgoing node = " << g_node_vector[g_link_vector[link_seq_no].to_node_seq_no].node_id << '\n';
				}
			}
			else
			{
				if (dtalog.debug_level() == 3)
				{
					dtalog.output() << "[DATA INFO] node id= " << g_node_vector[i].node_id << " with " << g_node_vector[i].m_outgoing_link_seq_no_vector.size() << " outgoing links." << '\n';
					g_DTA_log_file << "[DATA INFO] node id= " << g_node_vector[i].node_id << " with " << g_node_vector[i].m_outgoing_link_seq_no_vector.size() << " outgoing links." << '\n';
					for (int j = 0; j < g_node_vector[i].m_outgoing_link_seq_no_vector.size(); ++j)
					{
						int link_seq_no = g_node_vector[i].m_outgoing_link_seq_no_vector[j];
						dtalog.output() << "[DATA INFO] outgoing node = " << g_node_vector[g_link_vector[link_seq_no].to_node_seq_no].node_id << '\n';
						g_DTA_log_file << "[DATA INFO] outgoing node = " << g_node_vector[g_link_vector[link_seq_no].to_node_seq_no].node_id << '\n';
					}
				}
			}
		}
	}
	//if (assignment.assignment_mode != 11)   // not tmc mode
	//	g_read_link_qvdf_data(assignment);


}






//
//
//void g_reload_timing_arc_data(Assignment& assignment)
//{
//    dtalog.output() << "Step 1.7: Reading service arc in timing.csv..." << '\n';
//    g_DTA_log_file << "Step 1.7: Reading service arc in timing.csv..." << '\n';
//
//    CDTACSVParser parser_timing_arc;
//    if (parser_timing_arc.OpenCSVFile("timing.csv", false))
//    {
//        while (parser_timing_arc.ReadRecord())  // if this line contains [] mark, then we will also read field headers.
//        {
//            string mvmt_key;
//            if (!parser_timing_arc.GetValueByFieldName("mvmt_key", mvmt_key))
//            {
//                dtalog.output() << "Error: mvmt_key in file timing.csv is not defined." << '\n';
//                g_DTA_log_file << "Error: mvmt_key in file timing.csv is not defined." << '\n';
//                continue;
//            }
//            // create a link object
//            CSignalTiming timing_arc;
//
//            if (assignment.g_mvmt_key_to_link_no_map.find(mvmt_key) == assignment.g_mvmt_key_to_link_no_map.end())
//            {
//                dtalog.output() << "Error: mvmt_key " << mvmt_key << " in file timing.csv is not defined in link.csv." << '\n';
//                g_DTA_log_file << "Error: mvmt_key " << mvmt_key << " in file timing.csv is not defined in link.csv." << '\n';
//                //has not been defined
//                continue;
//            }
//            else
//            {
//                timing_arc.link_seq_no = assignment.g_mvmt_key_to_link_no_map[mvmt_key];
//                g_link_vector[timing_arc.link_seq_no].timing_arc_flag = true;
//            }
//
//            string time_period;
//            if (!parser_timing_arc.GetValueByFieldName("time_window", time_period))
//            {
//                dtalog.output() << "Error: Field time_window in file timing.csv cannot be read." << '\n';
//                g_DTA_log_file << "Error: Field time_window in file timing.csv cannot be read." << '\n';
//                g_program_stop();
//                break;
//            }
//
//            vector<float> global_minute_vector;
//
//            //input_string includes the start and end time of a time period with hhmm format
//            global_minute_vector = g_time_parser(time_period); //global_minute_vector incldue the starting and ending time
//            if (global_minute_vector.size() == 2)
//            {
//                if (global_minute_vector[0] < assignment.g_LoadingStartTimeInMin)
//                    global_minute_vector[0] = assignment.g_LoadingStartTimeInMin;
//
//                if (global_minute_vector[0] > assignment.g_LoadingEndTimeInMin)
//                    global_minute_vector[0] = assignment.g_LoadingEndTimeInMin;
//
//                if (global_minute_vector[1] < assignment.g_LoadingStartTimeInMin)
//                    global_minute_vector[1] = assignment.g_LoadingStartTimeInMin;
//
//                if (global_minute_vector[1] > assignment.g_LoadingEndTimeInMin)
//                    global_minute_vector[1] = assignment.g_LoadingEndTimeInMin;
//
//                if (global_minute_vector[1] < global_minute_vector[0])
//                    global_minute_vector[1] = global_minute_vector[0];
//
//            }
//            else
//                continue;
//
//            float time_interval = 0;
//
//
//            // capacity in the space time arcs
//            float capacity = 1;
//            parser_timing_arc.GetValueByFieldName("capacity", capacity);
//            timing_arc.VDF_capacity = max(0.0f, capacity);
//
//            // capacity in the space time arcs
//            parser_timing_arc.GetValueByFieldName("cycle_length", timing_arc.cycle_length);
//
//            // capacity in the space time arcs
//            parser_timing_arc.GetValueByFieldName("red_time", timing_arc.red_time);
//            parser_timing_arc.GetValueByFieldName("start_green_time", timing_arc.start_green_time);
//            parser_timing_arc.GetValueByFieldName("end_green_time", timing_arc.end_green_time);
//
//            for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
//            {
//                    // to do: we need to consider multiple periods in the future, Xuesong Zhou, August 20, 2020.
//                g_link_vector[timing_arc.link_seq_no].VDF_period[tau].red_time = timing_arc.red_time;
//                g_link_vector[timing_arc.link_seq_no].VDF_period[tau].cycle_length = timing_arc.cycle_length;
//            }
//
//
//            g_signal_timing_arc_vector.push_back(timing_arc);
//            assignment.g_number_of_timing_arcs++;
//
//            if (assignment.g_number_of_timing_arcs % 10000 == 0)
//                dtalog.output() << "reading " << assignment.g_number_of_timing_arcs << " timing_arcs.. " << '\n';
//                g_DTA_log_file << "reading " << assignment.g_number_of_timing_arcs << " timing_arcs.. " << '\n';
//        }
//
//        parser_timing_arc.CloseCSVFile();
//    }
//
//    
//    dtalog.output() << "Step 1.8: Reading file section [demand_file_list] in settings.csv..." << '\n';
//    g_DTA_log_file << "Step 1.8: Reading file section [demand_file_list] in settings.csv..." << '\n';
//    // we now know the number of links
//    dtalog.output() << "number of timing records = " << assignment.g_number_of_timing_arcs << '\n' << '\n';
//    g_DTA_log_file << "number of timing records = " << assignment.g_number_of_timing_arcs << '\n' << '\n';
//}
//


