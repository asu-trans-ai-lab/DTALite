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

#ifdef _WIN32
#define YAML_CPP_STATIC_DEFINE
#endif
#include <yaml-cpp/yaml.h>

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

	YAML::Node config = YAML::LoadFile("settings.yml");

	// Create a vector to hold all the demand file configurations
	// Check if 'demand_files' is a sequence and iterate over it
	if (config["subarea_analysis"].IsSequence()) {
		for (const YAML::Node& node : config["subarea_analysis"]) {

		int activate  = node["activate"].as<int>(0);
		if (activate == 0)
			return false;

			std::string geo_string = node["subarea_geometry"].as<std::string>();
			
				if (isGlobalPolygon(geo_string) == true && g_zone_vector.size() > 1000)
				{
					dtalog.output() << "[PROCESS INFO] The section subarea_geometry currently contains the default global coordinate system, bypassing the need for subarea processing, for a system with more than 1,000 zones." << '\n';
					g_DTA_log_file << "[PROCESS INFO] The section subarea_geometry currently contains the default global coordinate system, bypassing the need for subarea processing, for a system with more than 1,000 zones." << '\n';
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

				dtalog.output() << "[PROCESS INFO] Validated subarea_geometry. It contains " << assignment.g_subarea_shape_points.size() << " geometric points. Proceeding to FOCUSING analysis step." << '\n';
				g_DTA_log_file << "[PROCESS INFO] Validated subarea_geometry. It contains " << assignment.g_subarea_shape_points.size() << " geometric points. Proceeding to FOCUSING analysis step." << '\n';


				// the start interation of generating signals, if there is no signals set this number larger than the itertion number
				assignment.g_max_number_of_super_zones =  node["max_number_of_super_zones"].as<int>(100);


				break;
			}

		}
		
	


	return true;
}


void g_read_departure_time_profile(Assignment& assignment)
{

	dtalog.output() << "[PROCESS INFO] Step 2.0: Reading section departure_time_profile" << '\n';
	g_DTA_log_file << "[PROCESS INFO] Step 2.0: Reading section departure_time_profile" << '\n';

	YAML::Node config = YAML::LoadFile("settings.yml");

	if (config["departure_time_profile"].IsSequence()) {
		for (const YAML::Node& node : config["departure_time_profile"]) {
			{
				int departure_time_profile_no = 0;
				string time_period;

				CDeparture_time_Profile dep_time;
				departure_time_profile_no = node["departure_time_profile_no"].as<int>(1);

				dep_time.departure_time_profile_no = departure_time_profile_no;
				if (departure_time_profile_no != assignment.g_DepartureTimeProfileVector.size())
				{
					dtalog.output() << "[ERROR] Field departure_time_profile_no in field departure_time_profile should be sequential as a value of ." << assignment.g_DepartureTimeProfileVector.size() << '\n';
					g_DTA_log_file << "[ERROR] Field departure_time_profile_no in field departure_time_profile should be sequential as a value of ." << assignment.g_DepartureTimeProfileVector.size() << '\n';
					departure_time_profile_no = assignment.g_DepartureTimeProfileVector.size();

				}

			

				for (int s = assignment.starting_time_slot_no; s <= assignment.ending_time_slot_no; s += 1)
				{
					int hour = s / 12;
					int minute = (int)((s / 12.0 - hour) * 60 + 0.5);

					double value = 0;

					char time_interval_field_name[100];
					sprintf(time_interval_field_name, "T%02d%02d", hour, minute);

					dep_time.departure_time_ratio[s]  = node[time_interval_field_name].as<float>(1);


					if (s == dep_time.starting_time_slot_no || s == dep_time.ending_time_slot_no)
					{

						dtalog.output() << "[DATA INFO] At time " << hour << ":" << minute << " (" << time_interval_field_name << ") , the demand ratio is " << value << "." << '\n';
						g_DTA_log_file << "[DATA INFO] At time " << hour << ":" << minute << " (" << time_interval_field_name << ") , the demand ratio is " << value << "." << '\n';

					}

				}

				dep_time.compute_cumulative_profile(dep_time.starting_time_slot_no, dep_time.ending_time_slot_no, true);
				assignment.g_DepartureTimeProfileVector.push_back(dep_time);


			}

		}
	}
}

double g_pre_read_demand_file(Assignment& assignment)
{

	float total_demand_in_demand_file = 0;

	CDTACSVParser parser;
	
	// Define a struct to hold the demand file information
	struct DemandFile {
		std::string file_name;
		std::string format_type;
	};

	assignment.summary_file << "pre-read demand_file_list section." << '\n';

	//The aim of preloading the demand matrix under the subarea scenario is to gain insights into the structure of the Origin - Destination(OD) demand
	//	across the entire network.By doing so, we can identify critical OD pairsand zones based on their relevanceand importance.This understanding aids in the subsequent loading 
	//	of the demand matrix, enabling a more focusedand efficient approach.Consequently, this process allows for a significant reduction in the number of zones, 
	//	enhancing the efficiencyand performance of our traffic assignment tasks*/

	YAML::Node config = YAML::LoadFile("settings.yml");

	// Create a vector to hold all the demand file configurations
	std::vector<DemandFile> demandFiles;

	// Check if 'demand_files' is a sequence and iterate over it
	if (config["demand_files"].IsSequence()) {
		for (const YAML::Node& node : config["demand_files"]) {

			double loading_scale_factor = 1.0;
			string mode_type;

			DemandFile df;
			df.file_name = node["file_name"].as<std::string>("demand.csv");
			df.format_type = node["format_type"].as<std::string>("column");

			demandFiles.push_back(df);

			if (df.format_type.find("column") != string::npos || df.format_type.find("bin") != string::npos)  // or muliti-column
			{

				// try to detect if we have a route.csv as preload file, if yes, we skip the following reading .
				CDTACSVParser parser_route;

				struct SDemandHeader
				{
					int o_zone_id, d_zone_id, mode_type_no;
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


				if (file_format == 0 && df.format_type.find("column") != string::npos)
				{
					if (parser.OpenCSVFile(df.file_name, false) == false)
					{
						dtalog.output() << "[ERROR] column file " << df.file_name.c_str() << " does not exist." << '\n';
						g_DTA_log_file << "[ERROR] column file " << df.file_name.c_str() << " does not exist." << '\n';
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
			string file_name,  mode_type;
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
	const int fieldWidth = 12;
	dtalog.output() << "[PROCESS INFO]"; 
	g_DTA_log_file << "[PROCESS INFO]"; 
	dtalog.output() << '\n';
	g_DTA_log_file << '\n';

	dtalog.output() << "";
	g_DTA_log_file << "";
	dtalog.output() << std::setw(fieldWidth) << " mode_type";
	g_DTA_log_file << std::setw(fieldWidth) << "mode_type";
	dtalog.output() << std::setw(fieldWidth) << "total_demand";
	g_DTA_log_file << std::setw(fieldWidth) << "total_demand";
	dtalog.output() << std::setw(fieldWidth) << "#_of_links";
	g_DTA_log_file << std::setw(fieldWidth) << "#_of_links";
	dtalog.output() << std::setw(fieldWidth) << "speed_mph";
	g_DTA_log_file << std::setw(fieldWidth) << "speed_mph";
	dtalog.output() << std::setw(fieldWidth) << "speed_kmph";
	g_DTA_log_file << std::setw(fieldWidth) << "speed_kmph";
	dtalog.output() << std::setw(fieldWidth) << "length_km";
	g_DTA_log_file << std::setw(fieldWidth) << "length_km";
	dtalog.output() << std::setw(fieldWidth) << "avg_lane_cap";
	g_DTA_log_file << std::setw(fieldWidth) << "avg_lane_cap";
	dtalog.output() << std::setw(fieldWidth) << "avg_length_m";
	g_DTA_log_file << std::setw(fieldWidth) << "avg_length_m";
	dtalog.output() << '\n';
	g_DTA_log_file << '\n';

		for (int at = 0; at < assignment.g_ModeTypeVector.size(); at++)
		{
			dtalog.output() << std::setw(fieldWidth) <<  "," << assignment.g_ModeTypeVector[at].mode_type.c_str();
			g_DTA_log_file << std::setw(fieldWidth) << assignment.g_ModeTypeVector[at].mode_type.c_str();
			dtalog.output() << std::setw(fieldWidth) << assignment.total_demand[at] << "";
			g_DTA_log_file << std::setw(fieldWidth) << assignment.total_demand[at] << "";
		

			int link_count = 0;
			double total_speed = 0;
			double total_length = 0;
			double total_lane_capacity = 0;
			double total_link_capacity = 0;

			for (int i = 0; i < g_link_vector.size(); i++)
			{
				if (g_link_vector[i].link_type >= 0 && g_link_vector[i].AllowModeType(assignment.g_ModeTypeVector[at].mode_type))
				{
					link_count++;
					total_speed += g_link_vector[i].free_speed;
					total_length += g_link_vector[i].length_in_meter * g_link_vector[i].number_of_lanes;
					total_lane_capacity += g_link_vector[i].lane_capacity;
					total_link_capacity += g_link_vector[i].lane_capacity * g_link_vector[i].number_of_lanes;
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

			if (assignment.total_demand[at] > 0 && link_count == 0)
			{
				dtalog.output() << ", mode_type = " << assignment.g_ModeTypeVector[at].mode_type.c_str() << " has a total demand of " << assignment.total_demand[at] << ", but this mode type is not allowed on any link in the network." << '\n';
				g_DTA_log_file << ", mode_type = " << assignment.g_ModeTypeVector[at].mode_type.c_str() << " has a total demand of " << assignment.total_demand[at] << ", but this mode type is not allowed on any link in the network." << '\n';
 				dtalog.output() << "Please check the section link_type file to see if the 'allowed_uses' field allows this type of travel demand for major facilities, such as real-time information or HOV for highway facilities." << '\n';
				g_DTA_log_file << "Please check the section link_type file to see if the 'allowed_uses' field allows this type of travel demand for major facilities, such as real-time information or HOV for highway facilities." << '\n';

				g_program_stop();
			}
		}



}
void g_check_demand_volume_with_mode_type(Assignment& assignment)
{
	assignment.summary_file << ",summary by multi-modal and demand types,mode_type,total demmand volume in scenario 0, # of links,avg_free_speed_mph,avg_free_speed_kmph,total_length_in_km,total_capacity,avg_lane_capacity,avg_length_in_meter," << '\n';


		for (int at = 0; at < assignment.g_ModeTypeVector.size(); at++)
		{
			assignment.summary_file << ",," << ",";
			assignment.summary_file << assignment.g_ModeTypeVector[at].mode_type.c_str() << "," << assignment.total_demand[at] << ",";

			int link_count = 0;
			double total_speed = 0;
			double total_length = 0;
			double total_lane_capacity = 0;
			double total_link_capacity = 0;

			for (int i = 0; i < g_link_vector.size(); i++)
			{
				if (g_link_vector[i].link_type >= 0 && g_link_vector[i].AllowModeType(assignment.g_ModeTypeVector[at].mode_type))
				{
					link_count++;
					total_speed += g_link_vector[i].free_speed;
					total_length += g_link_vector[i].length_in_meter * g_link_vector[i].number_of_lanes;
					total_lane_capacity += g_link_vector[i].capacity;
					total_link_capacity += g_link_vector[i].capacity * g_link_vector[i].number_of_lanes;
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


void g_create_subarea_related_zone_structure(int max_number_of_super_zones)
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

			int subarea_cutoff_zone_size = 1000;
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


		if (g_related_zone_vector_size >= max_number_of_super_zones)  // 
		{
			// k mean method
			//random select 300 zones then find the closest zone to aggregate

			double super_zone_ratio = max_number_of_super_zones *1.0 / g_related_zone_vector_size;

			int super_zone_index = assignment.g_number_of_analysis_districts + 1;  // automated district  id start from the externally given district id

			for (int orig = 0; orig < g_zone_vector.size(); orig++)  // o
			{
				if (g_zone_vector[orig].sindex >= 0)  // no significant: skip
				{
					double selection_ratio = 1; 

					if (g_zone_vector[orig].subarea_inside_flag == 3) //inside 
					{
						selection_ratio = super_zone_ratio * 2;
					}
					if (g_zone_vector[orig].subarea_inside_flag == 2) //significant outside 
					{
						selection_ratio = super_zone_ratio * 0.7;
					}
					if (g_zone_vector[orig].subarea_inside_flag == 1) //less significant outside 
					{
						selection_ratio = super_zone_ratio * 0.2;
					}

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

void g_UpdateColumnPoolAfterLoadingRouteFile()
{
  int path_seq_count_modified = 0;
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

			for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)  //m
			{
					p_column_pool = &(assignment.g_column_pool[from_zone_sindex][to_zone_sindex][at]);
					if (p_column_pool->od_volume > 0)
					{
						path_seq_count_modified +=p_column_pool->ModifyColumnPoolAfterLoadingRouteFile();


					}

			}
		}
	}
	g_DTA_log_file << "[DATA INFO] Total path modified = " << path_seq_count_modified << '\n';
		dtalog.output() << "[DATA INFO] Total path modified = " << path_seq_count_modified << '\n';

}

void g_ReadRouteFile(Assignment& assignment)
{

			int path_counts = 0;
			float sum_of_path_volume = 0;

			CDTACSVParser parser;
			if (parser.OpenCSVFile("route.csv", false))
			{
				int total_path_in_demand_file = 0;
				// read agent file line by line,

				int o_zone_id, d_zone_id;
				int  mode_type_no = 0;


				std::vector <int> node_sequence;

				while (parser.ReadRecord())
				{
					total_path_in_demand_file++;
					if (total_path_in_demand_file % 10000 == 0)
					{
					
					dtalog.output() << "[DATA INFO] total routes in route file is " << total_path_in_demand_file << '\n';
					g_DTA_log_file << "[DATA INFO] total routes in route file is is " << total_path_in_demand_file << '\n';
					}

					parser.GetValueByFieldName("o_zone_id", o_zone_id);
					parser.GetValueByFieldName("d_zone_id", d_zone_id);
					parser.GetValueByFieldName("mode_type_id", mode_type_no);
					
					if (mode_type_no >= assignment.g_ModeTypeVector.size())
					{
						// error output 
						break; 
					}

					CAgentPath agent_path_element;

					parser.GetValueByFieldName("route_seq_id", agent_path_element.path_id, true);


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
					//volume *= loading_scale_factor;
					agent_path_element.volume = volume;
					path_counts++;
					sum_of_path_volume += agent_path_element.volume;

				 assignment.g_column_pool[from_zone_sindex][to_zone_sindex][mode_type_no].od_volume_from_route_file += agent_path_element.volume;

	
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

						CColumnVector* pColumnVector = &(assignment.g_column_pool[from_zone_sindex][to_zone_sindex][mode_type_no]);
						pColumnVector->departure_time_profile_no = 0;
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
				dtalog.output() << "[DATA INFO] total_demand_volume loaded from route file is " << sum_of_path_volume << " with " << path_counts << " paths." << '\n';
				g_DTA_log_file << "[DATA INFO] total_demand_volume loaded from route file is " << sum_of_path_volume << " with " << path_counts << " paths." << '\n';
			}
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

	dtalog.output() << "[PROCESS INFO] Step 1.7: reading section subarea" << '\n';
	g_DTA_log_file << "[PROCESS INFO] Step 1.7: reading section subarea" << '\n';
	g_read_subarea_CSV_file(assignment);
	// subarea handling step 2:
	// for each OD

	if (assignment.g_subarea_shape_points.size() >= 3) // if there is a subarea defined and not using preloaded route file
	{

		for (int orig = 0; orig < g_zone_vector.size(); orig++)  // o
		{
			g_zone_vector[orig].subarea_inside_flag = 0;  //reset the subarea inside flag to 0 once there is a subarea
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


		g_create_subarea_related_zone_structure(assignment.g_max_number_of_super_zones);

	}

	//	fprintf(g_pFileOutputLog, "number of zones =,%lu\n", g_zone_vector.size());
	g_read_departure_time_profile(assignment);

	assignment.InitializeDemandMatrix(g_related_zone_vector_size, g_zone_vector.size(), assignment.g_ModeTypeVector.size());

	float related_flow_matrix[4][4] = { 0 };

	float total_demand_in_demand_file = 0;

	CDTACSVParser parser;

	dtalog.output() << "[PROCESS INFO] Step 2.1: Reading file demand_file_list section..." << '\n';
	g_DTA_log_file << "[PROCESS INFO] Step 2.1: Reading file demand_file_list section..." << '\n';


	// Define a struct to hold the demand file information
	struct DemandFile {
		int file_sequence_no;
		std::string file_name;
		std::string mode_type;
		std::string format_type;
		float scale_factor;
		int departure_time_profile_no;
	};

	assignment.summary_file << "[PROCESS INFO] Step 2.1: read demand, defined in demand_files section." << '\n';
	int scenario_index_vector_error_count = 0;
	int reading_demand_file_log_count = 0;

	int count = 0;
	YAML::Node config = YAML::LoadFile("settings.yml");

	// Create a vector to hold all the demand file configurations
	std::vector<DemandFile> demandFiles;

	// Check if 'demand_files' is a sequence and iterate over it
	if (config["demand_files"].IsSequence()) {
		for (const YAML::Node& node : config["demand_files"]) {

			double loading_scale_factor = 1.0;
			string file_name,  mode_type, format_type;

			DemandFile df;
			df.file_sequence_no = node["file_sequence_no"].as<int>(1);

			df.file_name = node["file_name"].as<std::string>("demand.csv");
			df.mode_type = node["mode_type"].as<std::string>("auto");
			df.format_type = node["format_type"].as<std::string>("column");
			df.scale_factor = node["scale_factor"].as<float>(1);
			df.departure_time_profile_no = node["departure_time_profile_no"].as<int>(1);

			demandFiles.push_back(df);

			int file_sequence_no = df.file_sequence_no;
			file_name = df.file_name;
			mode_type = df.mode_type;
			format_type = df.format_type;
			loading_scale_factor = df.scale_factor;
			int this_departure_time_profile_no = df.departure_time_profile_no;



			if (reading_demand_file_log_count < 5)
			{
				dtalog.output() << "[STATUS INFO] reading demand file " << file_name.c_str() << '\n';
				g_DTA_log_file << "[STATUS INFO] reading demand file " << file_name.c_str() << '\n';
				reading_demand_file_log_count++;
			}


			if (loading_scale_factor < 0.0001)
				continue;



			if (this_departure_time_profile_no >= assignment.g_DepartureTimeProfileVector.size())
			{
				dtalog.output() << "[ERROR] departure_time_profile_no = " << this_departure_time_profile_no << " in  demand_files section has not been defined in section departure_time_profile." << '\n';
				g_DTA_log_file << "[ERROR] departure_time_profile_no = " << this_departure_time_profile_no << " in  demand_files section has not been defined in section departure_time_profile." << '\n';
				this_departure_time_profile_no = 0;

			}

			int mode_type_no = 0;


			//char time_interval_field_name[20];



			if (assignment.starting_time_slot_no * MIN_PER_TIMESLOT < assignment.g_LoadingStartTimeInMin)
				assignment.g_LoadingStartTimeInMin = assignment.starting_time_slot_no * MIN_PER_TIMESLOT;

			if (assignment.ending_time_slot_no * MIN_PER_TIMESLOT > assignment.g_LoadingEndTimeInMin)
				assignment.g_LoadingEndTimeInMin = assignment.ending_time_slot_no * MIN_PER_TIMESLOT;

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
				dtalog.output() << "[ERROR] mode_type = " << mode_type.c_str() << " in field mode_type of demand_files section is not defined in the section mode_type yet." << '\n';
				g_DTA_log_file << "[ERROR] mode_type = " << mode_type.c_str() << " in field mode_type of demand_files section is not defined in the section mode_type yet." << '\n';
				mode_type = "auto";
			}



			if (format_type.find("column") != string::npos || format_type.find("bin") != string::npos)  // or muliti-column
			{

				// try to detect if we have a route.csv as preload file, if yes, we skip the following reading .
				CDTACSVParser parser_route;

				struct SDemandHeader
				{
					int o_zone_id, d_zone_id, mode_type_no;
					double volume;
				};

				SDemandHeader header;


				CDTACSVParser parser;
				FILE* pFile;
				int file_exists = 0;
				int line_no = 0;
				int file_format = 0;

				//if (parser_route.OpenCSVFile("route.csv", false))
				//{
				//	// a route file exists
				//	dtalog.output() << "[STATUS INFO] route.csv exists as preload file so we skip the reading for the column based demand file." << '\n';
				//	g_DTA_log_file << "[STATUS INFO] route.csv exists as preload file so we skip the reading for the column based demand file." << '\n';
				//}
				//else
				//{
				//	fopen_ss(&pFile, "demand.bin", "rb");
				//	if (pFile != NULL)
				//	{
				//		file_format = 2;
				//		dtalog.output() << "[STATUS INFO] reading demand.bin in fast binary file reading mode." << '\n';
				//		g_DTA_log_file << "[STATUS INFO] reading demand.bin in fast binary file reading mode." << '\n';
				//	}


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

					}

					// end of reading read each record


					if (o_zone_id == d_zone_id)
					{
						continue;
					}

					if (assignment.g_zoneid_to_zone_seq_no_mapping.find(o_zone_id) == assignment.g_zoneid_to_zone_seq_no_mapping.end())
					{
						if (error_count < 10 && missing_zone_map.find(o_zone_id) == missing_zone_map.end())
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

					g_scenario_summary.record_mode_volume(mode_type_no, demand_value);

					assignment.total_demand[mode_type_no] += demand_value;
					assignment.g_column_pool[from_zone_sindex][to_zone_sindex][mode_type_no].od_volume += demand_value;
					assignment.g_column_pool[from_zone_sindex][to_zone_sindex][mode_type_no].departure_time_profile_no = this_departure_time_profile_no;
					assignment.total_demand_volume += demand_value;

					if (assignment.g_ModeTypeVector[mode_type_no].real_time_information_type >= 1)
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

				dtalog.output() << "[DATA INFO] reading file " << file_name.c_str() << ", cumulative total demand volume is " << assignment.total_demand_volume << '\n';
				g_DTA_log_file << "[DATA INFO] reading file " << file_name.c_str() << ", cumulative total demand volume is " << assignment.total_demand_volume << '\n';
				//dtalog.output() << "crtical demand volume has " << critical_OD_count << " OD pairs in size," << critical_OD_volume << ", " << ", account for " << critical_OD_volume / max(0.1f, assignment.total_demand_volume) * 100 << "%%" << '\n';
				//g_DTA_log_file << "crtical demand volume has " << critical_OD_count << " OD pairs in size," << critical_OD_volume << ", " << ", account for " << critical_OD_volume / max(0.1f, assignment.total_demand_volume) * 100 << "%%" << '\n';

				//dtalog.output() << "crtical OD zones volume has " << critical_OD_count << " OD pairs in size," << critical_OD_volume << ", " << ", account for " << critical_OD_volume / max(0.1f, assignment.total_demand_volume) * 100 << "%%" << '\n';
				//g_DTA_log_file << "crtical OD zones volume has " << critical_OD_count << " OD pairs in size," << critical_OD_volume << ", " << ", account for " << critical_OD_volume / max(0.1f, assignment.total_demand_volume) * 100 << "%%" << '\n';
			}
			/////
			/// summary
			//////


			g_ReadRouteFile(assignment);
			g_UpdateColumnPoolAfterLoadingRouteFile();

			/// <summary>
			///  subarea
			/// </summary>
			/// <param name="assignment"></param>
			/*assignment.summary_file << ",total demand =, " << assignment.total_demand_volume << '\n';*/

			g_check_demand_volume_with_mode_type(assignment);
		
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
	link.allocate_memory();
	link.link_id = "connector";

	link.from_node_seq_no = internal_from_node_seq_no;
	link.to_node_seq_no = internal_to_node_seq_no;
	link.link_seq_no = assignment.g_number_of_links;
	link.to_node_seq_no = internal_to_node_seq_no;
	//virtual connector
		link.link_type = -1;
		link.free_speed = 100;
		link.capacity = 2000;
		link.number_of_lanes = 20;  // default all open

	//only for outgoing connectors
	link.zone_seq_no_for_outgoing_connector = zone_seq_no;

	//BPR
	link.traffic_flow_code = 0;

	link.spatial_capacity_in_vehicles = 99999;
	link.lane_capacity = 999999;
	link.link_spatial_capacity = 99999;
	link.link_distance_VDF = 0.00001;
	link.free_flow_travel_time_in_min = 0.1;

		//setup default values
		link.lane_based_ultimate_hourly_capacity = 99999;
		// 60.0 for 60 min per hour

		link.alpha = 0;
		link.beta = 0;

		for (int at = 0; at < assignment.g_ModeTypeVector.size(); at++)
		{
			link.link_avg_travel_time[at] = 0;
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
	
	dtalog.output() << "[STATUS INFO] Reading demand_files section..." << '\n';
	g_DTA_log_file << "[STATUS INFO] Reading demand_files section..." << '\n';
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

			string file_name, mode_type;
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

extern unsigned int g_RandomSeed;
extern void InitWELLRNG512a(unsigned int* init);

void g_detector_file_open_status(Assignment& assignment)
{
	FILE* g_pFilePathMOE = nullptr;

	
	fopen_ss(&g_pFilePathMOE, "od_performance.csv", "w");
	if (!g_pFilePathMOE)
	{
		dtalog.output() << "[ERROR] File od_performance.csv cannot be opened." << '\n';
		g_DTA_log_file << "[ERROR] File od_performance.csv cannot be opened." << '\n';
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

template<typename T>
bool ReadValueFromNode(const YAML::Node& node, const std::string& key, T& value, bool required) {
	if (node[key]) {
		value = node[key].as<T>();
		return true;
	}
	else if (required) {
		throw std::runtime_error("Required field '" + key + "' is missing.");
	}
	return false;
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

	std::string fileName = "settings.yml";
	std::ifstream inFile(fileName.c_str());
	YAML::Node settings;
	if (!inFile.is_open()) {
		dtalog.output() << "Error opening file: " << fileName << std::endl;
		return;
	}

	try {
		settings = YAML::Load(inFile);
		// Now you can work with the 'settings' as a YAML::Node object
	}
	catch (const YAML::ParserException& e) {
		dtalog.output() << "Error parsing the file: " << e.what() << std::endl;
		return;
	}
	dtalog.output() << "[PROCESS INFO] Step 1: Reading input data" << '\n';
	g_DTA_log_file << "[PROCESS INFO] Step 1: Reading input data" << '\n';
	dtalog.output() << "[PROCESS INFO] Step 1.1: Reading section demand_period ..." << '\n';
	g_DTA_log_file << "[PROCESS INFO] Step 1.1: Reading section  demand_period..." << '\n';
	assignment.summary_file << "[PROCESS INFO] Step 1.1: Reading section demand_period ..." << '\n';


	try {
	const auto& demand_periods = settings["demand_period"];
	
	try
	{
		for (const YAML::Node& demand_period_node : demand_periods)
		{


			std::string time_period;

				time_period = demand_period_node["time_period"].as<std::string>("0700-0800");

				vector<float> global_minute_vector;

				//input_string includes the start and end time of a time period with hhmm format
				global_minute_vector = g_time_parser(time_period); //global_minute_vector incldue the starting and ending time

			if (global_minute_vector.size() == 2)
			{
				assignment.starting_time_slot_no = global_minute_vector[0] / MIN_PER_TIMESLOT;  // read the data
				assignment.ending_time_slot_no = global_minute_vector[1] / MIN_PER_TIMESLOT;    // read the data from settings.csv
				assignment.time_period_in_hour = (global_minute_vector[1] - global_minute_vector[0]) / 60.0;
				assignment.t2_peak_in_hour = (global_minute_vector[0] + global_minute_vector[1]) / 2 / 60;

				//g_fout << global_minute_vector[0] << '\n';
				//g_fout << global_minute_vector[1] << '\n';


			//	string peak_time_str = demand_period_node["peak_time"].as<std::string>("0000");

			//	demand_period.t2_peak_in_hour = g_timestamp_parser(peak_time_str) / 60.0;


			}

			CDeparture_time_Profile dep_time;
			dep_time.starting_time_slot_no = assignment.starting_time_slot_no;
			dep_time.ending_time_slot_no = assignment.ending_time_slot_no;

			for (int s = 0; s <= 96 * 3; s++)
			{
				dep_time.departure_time_ratio[s] = 1.0 / 300.0;
			}

			//			dtalog.output() << "[DATA INFO] A default flat departure time profile is used..." << '\n';
			//			g_DTA_log_file << "[DATA INFO] A default flat departure time profile is used..." << '\n';
			dep_time.compute_cumulative_profile(assignment.starting_time_slot_no, assignment.ending_time_slot_no, false);

			if (assignment.g_DepartureTimeProfileVector.size() == 0)
			{
				//default profile
				assignment.g_DepartureTimeProfileVector.push_back(dep_time);
			}

			assignment.summary_file << ",demand_period= " << global_minute_vector[0] << "-> " << global_minute_vector[1] << '\n';

		}
	}

	catch (const std::exception& e)
	{
		std::cerr << "Exception occurred: " << e.what() << std::endl;
	}

	}
	catch (const YAML::ParserException& e) {
		dtalog.output() << "Error parsing the section demand_period in settings.yml: " << e.what() << std::endl;
		return;
	}
	

//step 1:read demand type file

	dtalog.output() << "[PROCESS INFO] Step 1.2: Reading section mode_type..." << '\n';
	g_DTA_log_file << "[PROCESS INFO] Step 1.2: Reading section mode_type..." << '\n';
	assignment.summary_file << "[PROCESS INFO] Step 1.2: Reading section mode_type..." << '\n';


	unsigned short i = 0;
	const auto& mode_types = settings["mode_types"];

	assignment.g_ModeTypeVector.clear();


		try
		{
			for (const YAML::Node& mode_node : mode_types)
			{
			Cmode_type mode_type;

			// auto type_ = a["type"];
			auto mode_type_name = mode_node["mode_type"].as<std::string>("auto");
			////if (contains_agent_name(name))
			//{
			//	std::cerr << "duplicate agent type found: " << name << '\n';
			//	continue;
			//}
			mode_type.mode_type = mode_type_name;
			mode_type.mode_type_no = assignment.g_ModeTypeVector.size() + 1;

			auto mode_type_index = mode_node["mode_type_index"].as<unsigned short>(1);

			mode_type.value_of_time = mode_node["vot"].as<float>(20);
/*			mode_type.eco_so_flag = a["eco_so_flag"].as<float>(20)*/;



		// Assuming the YAML content is under a node called 'sensor_data'

			mode_type.pce = mode_node["pce"].as<int>(1);

			mode_type.person_occupancy = mode_node["person_occupancy"].as<float>(1.0);
			mode_type.desired_speed_ratio = mode_node["desired_speed_ratio"].as<float>(1.0);
			mode_type.time_headway_in_sec = mode_node["time_headway_in_sec"].as<float>(1);
			mode_type.display_code = mode_node["display_code"].as<int>(1);
			mode_type.real_time_information_type = mode_node["DTM_real_time_info_type"].as<float>(0);

			assignment.mode_type_2_seqno_mapping[mode_type.mode_type] = assignment.g_ModeTypeVector.size();

			assignment.g_ModeTypeVector.push_back(mode_type);
			assignment.summary_file << "mode_type =, " << mode_type.mode_type.c_str() << ", real time info flag = " << mode_type.real_time_information_type << '\n';


			//substring overlapping checking

			//{
			//	for (int at = 0; at < assignment.g_ModeTypeVector.size(); at++)
			//	{
			//		if (assignment.g_ModeTypeVector[at].mode_type.find(mode_type.mode_type) != string::npos)
			//		{
			//			dtalog.output() << "[ERROR] Error substring duplication checking : mode_type = " << assignment.g_ModeTypeVector[at].mode_type.c_str() <<
			//				" in section mode_type is overlapping with " << mode_type.mode_type.c_str() << ". Please add flags such as to avoid overlapping in the use of allowe_uses field.";
			//	
			//			g_DTA_log_file << "[ERROR] Error substring duplication checking : mode_type = " << assignment.g_ModeTypeVector[at].mode_type.c_str() <<
			//				" in section mode_type is overlapping with " << mode_type.mode_type.c_str() << ". Please add flags such as to avoid overlapping in the use of allowe_uses field.";

			//		}

			//	}

			//}



			if (mode_type.real_time_information_type == 1)  // real time info
			{
				assignment.g_number_of_real_time_mode_types++;
			}

			if (mode_type.real_time_information_type == 2)  //dms
			{
				assignment.g_number_of_DMS_mode_types++;
			}


		}
	

		}
		catch (const std::exception& e)
		{
			std::cerr << "Exception occurred: " << e.what() << std::endl;
		}
		assignment.g_number_of_mode_types = assignment.g_ModeTypeVector.size();
		g_number_of_active_mode_types = assignment.g_ModeTypeVector.size();



	if (assignment.g_ModeTypeVector.size() == 0)
	{
		dtalog.output() << "[ERROR] File mode_type does not have information." << '\n';
		g_DTA_log_file << "[ERROR] File mode_type does not have information." << '\n';
		g_program_stop();
	}

	if (assignment.g_ModeTypeVector.size() >= MAX_MODETYPES)
	{
		dtalog.output() << "[ERROR] mode_type = " << assignment.g_ModeTypeVector.size() << " in section mode_types is too large. " << '\n' << "MAX_MODETYPES = " << MAX_MODETYPES << "Please contact program developers!";
		g_DTA_log_file << "[ERROR] mode_type = " << assignment.g_ModeTypeVector.size() << " in section mode_types is too large. " << '\n' << "MAX_MODETYPES = " << MAX_MODETYPES << "Please contact program developers!";
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
	//				dtalog.output() << "[ERROR] Field mode_chain  in file activity_travel_pattern.csv has a value " << mode_type.c_str() << "not defined yet in section mode_type." << '\n';
	//				g_DTA_log_file << "[ERROR] Field mode_chain  in file activity_travel_pattern.csv has a value " << mode_type.c_str() << "not defined yet in section mode_type." << '\n';
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
	//				dtalog.output() << "[ERROR] Field mode_chain  in file activity_travel_pattern.csv has a value " << mode_type << "not defined yet in section mode_type." << '\n';
	//				g_DTA_log_file << "[ERROR] Field mode_chain  in file activity_travel_pattern.csv has a value " << mode_type << "not defined yet in section mode_type." << '\n';
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



	dtalog.output() << "[PROCESS INFO] Step 1.3: Reading section link_type" << '\n';
	g_DTA_log_file << "[PROCESS INFO] Step 1.3: Reading section link_type" << '\n';



	int emission_log_count = 0; 
	int meu_log_count = 0; 
	int peak_load_factor_log_count = 0; 
	int allowed_use_log_count = 0; 
	
	// Create a vector to hold all the scenario configurations


	// Check if 'scenarios' is a sequence and iterate over it
	if (settings["link_types"].IsSequence())
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
		for (const YAML::Node& node : settings["link_types"]) {
			CLinkType element;

			element.link_type = node["link_type"].as<int>(1);
			element.link_type_name = node["link_type_name"].as<std::string>("auto");

			if (assignment.g_LinkTypeMap.find(element.link_type) != assignment.g_LinkTypeMap.end())
			{
				overlapping_link_type_vector.push_back(element.link_type);

				continue;
			}

				string traffic_flow_code_str;
				element.type_code = node["type_code"].as<std::string>("a");



				element.traffic_flow_code = spatial_queue;

				traffic_flow_code_str = node["traffic_flow_model"].as<std::string>("spatial_queue");
				element.k_jam = node["k_jam_km"].as<float>(200);

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


				// Loop over all mode types (e.g. auto, bike, bus, etc.)
				for (int at = 0; at < assignment.g_ModeTypeVector.size(); at++)
				{
					// Generate the field name in the CSV file for CO2 emissions for this mode type
					char CSV_field_name[50];
					sprintf(CSV_field_name, "emissions_%s_co2", assignment.g_ModeTypeVector[at].mode_type.c_str());

					string emissions_co2_str;

					// Check if the field could not be found in the CSV file


					try {
						emissions_co2_str = node[CSV_field_name].as<std::string>("");

					}
					catch (const YAML::ParserException& e) {
					
						// If the field cannot be found, output a warning and use a default value of 0.0
						if (line_no == 0 && emission_log_count < 4) {
							emission_log_count++;
							dtalog.output() << "[WARNING] Field '" << CSV_field_name << "' not found in section link_type. The default value of 0 was used. Consider adding '" << CSV_field_name << "' to the section link_type for more accurate results." << '\n';
							g_DTA_log_file << "[WARNING] Field '" << CSV_field_name << "' not found in section link_type. The default value of 0 was used. Consider adding '" << CSV_field_name << "' to the section link_type for more accurate results." << '\n';
						}
					}

					// Convert the string containing CO2 emissions coefficients into a vector of doubles
					std::vector<double> emissions_co2_coeff_vector;
					g_ParserDoubleSequence(emissions_co2_str, emissions_co2_coeff_vector);

					// Store up to the first 4 emissions coefficients in the element's matrix
					for (int i = 0; i < min(static_cast<size_t>(4), emissions_co2_coeff_vector.size()); i++)
					{
						element.emissions_co2_matrix[at][i] = emissions_co2_coeff_vector[i];
					}

					// Repeat the same process for NOx emissions:
					sprintf(CSV_field_name, "emissions_%s_nox", assignment.g_ModeTypeVector[at].mode_type.c_str());

					string emissions_nox_str;

					try {
						emissions_nox_str = node[CSV_field_name].as<std::string>("");
					}
					catch (const YAML::ParserException& e) 	{
						if (line_no == 0 && emission_log_count < 4) {
							emission_log_count++;
							dtalog.output() << "[WARNING] Field '" << CSV_field_name << "' not found in section link_type. The default value of 0 was used. Consider adding '" << CSV_field_name << "' to the section link_type for more accurate results." << '\n';
							g_DTA_log_file << "[WARNING] Field '" << CSV_field_name << "' not found in section link_type. The default value of 0 was used. Consider adding '" << CSV_field_name << "' to the section link_type for more accurate results." << '\n';
						}
					}

					std::vector<double> emissions_nox_coeff_vector;
					g_ParserDoubleSequence(emissions_nox_str, emissions_nox_coeff_vector);

					for (int i = 0; i < min(static_cast<size_t>(4), emissions_nox_coeff_vector.size()); i++)
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
				dtalog.output() << "[WARNING] the following link type(s) have been defined more than once in file section link_type:  ";
				g_DTA_log_file << "[WARNING] the following link type(s) have been defined more than once in file section link_type:  ";

				for (int i = 0; i < overlapping_link_type_vector.size(); i++)
				{
					dtalog.output() << overlapping_link_type_vector[i] << ", ";
					g_DTA_log_file << overlapping_link_type_vector[i] << ", ";
				}


			}





		

	}




	dtalog.output() << "[DATA INFO] number of link types = " << assignment.g_LinkTypeMap.size() << '\n';
	g_DTA_log_file << "[DATA INFO] number of link types = " << assignment.g_LinkTypeMap.size() << '\n';



	assignment.g_number_of_nodes = 0;
	assignment.g_number_of_links = 0;  // initialize  the counter to 0


	int number_of_zones = g_detect_if_zones_defined_in_node_csv(assignment);
	// = 1: normal
	//= 0, there are is boundary
	//=-1, no information



	

	int internal_node_seq_no = 0;
	// step 3: read node file


	std::map<int, int> zone_id_mapping;  // this is used to mark if this zone_id has been identified or not
	std::map<int, double> zone_id_x;
	std::map<int, double> zone_id_y;



	std::map<int, float> zone_id_production;
	std::map<int, float> zone_id_attraction;

	CDTACSVParser parser;

	int multmodal_activity_node_count = 0;


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
				link.allocate_memory();

				string link_type_name_str;
				parser_link.GetValueByFieldName("link_type_name", link_type_name_str, false);

				long from_node_id = -1;
				if (!parser_link.GetValueByFieldName("from_node_id", from_node_id))
					continue;

				long to_node_id = -1;
				if (!parser_link.GetValueByFieldName("to_node_id", to_node_id))
					continue;

				link.from_node_id = from_node_id;
				link.to_node_id = to_node_id;

				string linkID;
				parser_link.GetValueByFieldName("link_id", linkID, false);
				// add the to node id into the outbound (adjacent) node list


				if (assignment.g_node_id_to_seq_no_map.find(from_node_id) == assignment.g_node_id_to_seq_no_map.end())
				{


					dtalog.output() << "[ERROR] from_node_id " << from_node_id << " in file link.csv is not defined in node.csv." << '\n';
					g_DTA_log_file << "[ERROR] from_node_id " << from_node_id << " in file link.csv is not defined in node.csv." << '\n';

					continue; //has not been defined
				}

				if (assignment.g_node_id_to_seq_no_map.find(to_node_id) == assignment.g_node_id_to_seq_no_map.end())
				{


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

				if (parser_link.GetValueByFieldName("link_special_flag", link.link_specifical_flag_str, false) && link.link_specifical_flag_str.size() > 0)
				{
					dtalog.output() << "[INFO] link_special_flag =  " << link.link_specifical_flag_str.c_str() << " is defined for link " << from_node_id << "->" << to_node_id << " in link.csv." << '\n';
					g_DTA_log_file << "[INFO] link_special_flag =  " << link.link_specifical_flag_str.c_str() << " is defined for link " << from_node_id << "->" << to_node_id << " in link.csv." << '\n';

				}

				link.tmc_corridor_name = "network_wide";
				parser_link.GetValueByFieldName("tmc_corridor_name", link.tmc_corridor_name, false);
				parser_link.GetValueByFieldName("link_type_name", link.link_type_name, false);


				// and valid
				if (movement_str.size() > 0)
				{
					int main_node_id = -1;


					link.mvmt_txt_id = movement_str;
					link.main_node_id = main_node_id;
				}


				int default_link_type = assignment.g_first_link_type;
				char link_type_field_name[50];

					link.link_type = default_link_type;
					parser_link.GetValueByFieldName("link_type", link.link_type);
					

					if (assignment.g_LinkTypeMap.find(link.link_type) == assignment.g_LinkTypeMap.end())
					{
						missing_link_type_mapping[link.link_type] = 1;

						link.link_type = assignment.g_first_link_type;
						dtalog.output() << "[WARNING] Field link_type " << link.link_type << "  in link.csv is not defined. The default value in the field link_type in link.csv is used." << '\n';
						g_DTA_log_file << "[WARNING] Field link_type " << link.link_type << "  in link.csv is not defined. The default value in the field link_type in link.csv is used." << '\n';
						line_no++;
						continue;
					}




				if (assignment.g_LinkTypeMap[link.link_type].type_code == "c")  // suggestion: we can move "c" as "connector" in allowed_uses
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
				char free_speed_scenario_field_name[50];
				char capacity_scenario_field_name[50];

				int default_number_of_lanes = 1;
				parser_link.GetValueByFieldName("lanes", default_number_of_lanes, false);
				link.number_of_lanes = default_number_of_lanes;

				
	
				double length_in_meter = 1000.0; // km or mile
				double free_speed_original = 100.0;
				double free_speed = 60.0;

				double lane_capacity = 2000;
				// first step, get the initial value based on link type for speed and capacity


				double cutoff_speed = 1.0;
				double k_jam = assignment.g_LinkTypeMap[link.link_type].k_jam;
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
				//-----------------
				// second step, we read the link-specific value (only for based mode)
				if (parser_link.GetValueByFieldName("free_speed", free_speed_original, false, false) == false)
				{

					if (free_speed_missing_error == 0)
					{
						dtalog.output() << "[ERORR] Field free_speed in link.csv is missing." << '\n';
						g_DTA_log_file << "[ERORR] Field free_speed in link.csv is missing." << '\n';
						free_speed_original = 100;
						free_speed_missing_error++; 
					}
				}
				if (assignment.g_speed_unit_flag == 1)  // mph;
					free_speed = free_speed_original * 1.609; // convert from mile per hour to km per hour ( by multipling the original speed values by 1.609)

				double default_free_speed = free_speed; 
				
				{
					
					sprintf(free_speed_scenario_field_name, "free_speed");

					if (parser_link.GetValueByFieldName(free_speed_scenario_field_name, free_speed_original, false))
					{
						if (assignment.g_speed_unit_flag == 1)  // mph;
							free_speed = free_speed_original * 1.609; // convert from mile per hour to km per hour ( by multipling the original speed values by 1.609)
						
						link.free_speed = free_speed;
					}
					else
					{
						link.free_speed = default_free_speed;  // use default value
					}

				}

				//-----------------

				// second step, we read the link-specific value (only for based mode)
				if (parser_link.GetValueByFieldName("capacity", lane_capacity, false, false) == true)
				{
					
					{
						

						link.capacity = lane_capacity;
					}
	
				}


				double free_speed_value;
				//speed_unit: km/ph us_customary  si1
				// mile per h our ->

				if (link.link_id == "5")
					int iiidebug = 1;

				cutoff_speed = free_speed * 0.75; //default;

				//cutoff_speed = free_speed * 1.0; //default;

				if(parser_link.GetValueByFieldName("cutoff_speed", cutoff_speed, false) == true)  //read data 
				{
				
				if (assignment.g_speed_unit_flag == 1)  // mph;
					cutoff_speed = cutoff_speed * 1.609; // convert from mile per hour to km per hour ( by multipling the original speed values by 1.609)

				}


				link.v_congestion_cutoff = cutoff_speed;;

				free_speed = max(0.1, free_speed);

				link.free_speed = free_speed;
				link.lane_capacity = lane_capacity;
				int at_base = 0;

				// third step, we write the link-specific value (only for based mode)


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

				link.traffic_flow_code = assignment.g_LinkTypeMap[link.link_type].traffic_flow_code;

				//spatial queue and kinematic wave
				link.spatial_capacity_in_vehicles = max(1.0, link.link_distance_VDF * link.number_of_lanes * k_jam);

				// kinematic wave
				if (link.traffic_flow_code == 3)
					link.BWTT_in_simulation_interval = link.link_distance_VDF / bwtt_speed * 3600 / number_of_seconds_per_interval;

				link.kjam = assignment.g_LinkTypeMap[link.link_type].k_jam;
				char CSV_field_name[50];



				// reading for VDF related functions
				// step 1 read type


					////data initialization

					//link.est_volume_per_hour_per_lane[time_index] = 0;
					//link.est_avg_waiting_time_in_min[time_index] = 0;
					//link.est_queue_length_per_lane[time_index] = 0;
	


					//setup default values
					link.lane_based_ultimate_hourly_capacity = lane_capacity;
					link.number_of_lanes = link.number_of_lanes;

					link.vf = link.free_speed;
					link.v_VDF_congestion_cutoff = link.v_congestion_cutoff;   // mwe need to assign the value to the VDF data type
					link.alpha = 0.15;
					link.beta = 4;
					link.preload = 0;
					parser_link.GetValueByFieldName("VDF_alpha", link.alpha, false);
					parser_link.GetValueByFieldName("VDF_beta", link.beta, false);
					parser_link.GetValueByFieldName("VDF_plf", link.peak_load_factor, false);

					parser_link.GetValueByFieldName("VDF_cd", link.Q_cd, false);
					parser_link.GetValueByFieldName("VDF_cp", link.Q_cp, false);
					parser_link.GetValueByFieldName("VDF_n", link.Q_n, false);
					parser_link.GetValueByFieldName("VDF_s", link.Q_s, false);

					if (link.peak_load_factor < 0.1)
						link.peak_load_factor = 1;

					if (link.peak_load_factor > 1)
						link.peak_load_factor = 1;

					//double VDF_fftt = 0; 
					//if (parser_link.GetValueByFieldName("VDF_fftt", VDF_fftt, false, false) == true)
					//{
					//		link.free_flow_travel_time_in_min = VDF_fftt;
					//}

					parser_link.GetValueByFieldName("lanes", link.number_of_lanes);
					parser_link.GetValueByFieldName("allowed_uses", link.allowed_uses);

					

					for (int at = 0; at < assignment.g_ModeTypeVector.size(); at++)
					{
						link.toll[at] = 0;

					}

					link.starting_time_in_hour = assignment.starting_time_slot_no * MIN_PER_TIMESLOT / 60.0;
					link.ending_time_in_hour = assignment.ending_time_slot_no * MIN_PER_TIMESLOT / 60.0;
					link.L = assignment.time_period_in_hour;
					link.t2 = assignment.t2_peak_in_hour;

					string restricted_turn_nodes_str;
					parser_link.GetValueByFieldName("restricted_turn_nodes", restricted_turn_nodes_str, false);

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
							link.restricted_turn_nodes_str = restricted_turn_nodes_str;
							link.restricted_turn_nodes_map[node_no] = 1;

						}


					}

					double value_ref = -1;
					parser_link.GetValueByFieldName("ref_volume", value_ref, false, false);
					if (value_ref > 1)
						link.ref_link_volume = value_ref;

					parser_link.GetValueByFieldName("VDF_preload", link.preload, false, false);

					for (int at = 0; at < assignment.g_ModeTypeVector.size(); at++)
					{

						sprintf(CSV_field_name, "toll_%s", assignment.g_ModeTypeVector[at].mode_type.c_str());
						parser_link.GetValueByFieldName(CSV_field_name, link.toll[at], false, false);

						if (link.toll[at] > 0.001 && toll_message_count < 10)
						{
							dtalog.output() << "[DATA INFO] link " << from_node_id << "->" << to_node_id << " has a toll of " << link.toll[at] << " for agent type " << '\n';
							g_DTA_log_file << "[DATA INFO] link " << from_node_id << "->" << to_node_id << " has a toll of " << link.toll[at] << " for agent type " << 
								assignment.g_ModeTypeVector[at].mode_type.c_str() << " at demand period " << at << '\n';
							toll_message_count++;
						}
						//sprintf(CSV_field_name, "VDF_penalty_p%d", at);
						/*					parser_link.GetValueByFieldName(CSV_field_name, link.penalty, false, false);*/

					}

					if (link.cell_type >= 2) // micro lane-changing arc
					{
						// additinonal min: 1 second 1/60.0 min
						link.penalty += 1.0 / 60.0;
					}

					parser_link.GetValueByFieldName("cycle_length", link.cycle_length, false, false);

					if (link.cycle_length >= 1)
					{
						link.timing_arc_flag = true;

						parser_link.GetValueByFieldName("start_green_time", link.start_green_time);
						parser_link.GetValueByFieldName("end_green_time", link.end_green_time);

						link.effective_green_time = link.end_green_time - link.start_green_time;

						if (link.effective_green_time < 0)
							link.effective_green_time = link.cycle_length;

						link.red_time = max(1.0f, link.cycle_length - link.effective_green_time);
						parser_link.GetValueByFieldName("red_time", link.red_time, false);
						parser_link.GetValueByFieldName("green_time", link.effective_green_time, false);

						if (link.saturation_flow_rate > 1000)  // protect the data attributes to be reasonable
						{
							link.saturation_flow_rate = link.saturation_flow_rate;
						}
					}



				float default_cap = 1000;
				float default_BaseTT = 1;


				//link.m_OutflowNumLanes = number_of_lanes;//visum lane_cap is actually link_cap

				link.update_kc(free_speed);
				link.link_spatial_capacity = k_jam * link.number_of_lanes * link.link_distance_VDF;
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
				
				{
					
					
					for (int at = 0; at < assignment.g_ModeTypeVector.size(); at++)
					{
						char penalty_field_name[50];
						sprintf(penalty_field_name, "penalty_%s", assignment.g_ModeTypeVector[at].mode_type.c_str());
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

						if (dir_flag == 0 || dir_flag == 2)//
						{
							link.BA_link_no = link.link_seq_no + 1;
						}
						else
							link.BA_link_no = -1;

						g_node_vector[internal_from_node_seq_no].m_outgoing_link_seq_no_vector.push_back(link_seq_no);  // add this link to the corresponding node as part of outgoing node/link
						g_node_vector[internal_to_node_seq_no].m_incoming_link_seq_no_vector.push_back(link_seq_no);  // add this link to the corresponding node as part of outgoing node/link

						g_node_vector[internal_from_node_seq_no].m_to_node_seq_no_vector.push_back(link.to_node_seq_no);  // add this link to the corresponding node as part of outgoing node/link
						g_node_vector[internal_from_node_seq_no].m_to_node_2_link_seq_no_map[link.to_node_seq_no] = link.link_seq_no;  // add this link to the corresponding node as part of outgoing node/link
						g_link_vector.push_back(link);  // copy the memory to the global vector 
											
						assignment.g_number_of_links++;




					}
					else //link_code_start ==2
					{
						CLink link2; 
						link2 = link;  // first copy the content from link to link 2

						link2.allocate_memory();  // then allocate memory 

						// the following copy the content for fixed and dynamic array


								link2.allowed_uses = link.allowed_uses;


								link2.link_type = link.link_type;
								link2.number_of_lanes = link.number_of_lanes;  
								link2.free_speed = link.free_speed;
								link2.capacity = link.capacity;



						link2.to_node_id = link.from_node_id;
						link2.from_node_id = link.to_node_id;


						link2.AB_flag = -1;
						link2.BA_link_no = -1;


						link2.link_seq_no = link.link_seq_no + 1;  // update link seq no by + 1 for BA link
						link_seq_no = link2.link_seq_no;


						link2.from_node_seq_no = internal_to_node_seq_no;
						link2.to_node_seq_no = internal_from_node_seq_no;

						g_node_vector[internal_to_node_seq_no].m_outgoing_link_seq_no_vector.push_back(link_seq_no);  // add this link to the corresponding node as part of outgoing node/link
						g_node_vector[internal_from_node_seq_no].m_incoming_link_seq_no_vector.push_back(link_seq_no);  // add this link to the corresponding node as part of outgoing node/link

						g_node_vector[internal_to_node_seq_no].m_to_node_seq_no_vector.push_back(link2.to_node_seq_no);  // add this link to the corresponding node as part of outgoing node/link
						g_node_vector[internal_to_node_seq_no].m_to_node_2_link_seq_no_map[link2.to_node_seq_no] = link2.link_seq_no;  // add this link to the corresponding node as part of outgoing node/link
						g_link_vector.push_back(link2);  // copy the memory to the global vector 


						assignment.g_number_of_links++;

					}

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

	


	

	if (missing_link_type_mapping.size() > 0)
	{
		dtalog.output() << "[WARNING] The following link types in link.csv are not defined in section link_type. ";
		g_DTA_log_file << "[WARNING] The following link types in link.csv are not defined in section link_type. ";
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



	assignment.summary_file << "[PROCESS INFO] Step 1: read network node.csv, link.csv " << '\n';
	assignment.summary_file << ",# of nodes = ," << g_node_vector.size() << '\n';
	assignment.summary_file << ",# of links =," << g_link_vector.size() << '\n';
	assignment.summary_file << ",# of zones =," << g_zone_vector.size() << '\n';
	const int fieldWidth = 12;
	dtalog.output() << "[PROCESS INFO] Step 1: read network node.csv, link.csv" << '\n';
	g_DTA_log_file << "[PROCESS INFO] Step 1: read network node.csv, link.csv" << '\n';
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
			if (g_link_vector[i].link_type >= 0 && g_link_vector[i].link_type == it_link_type->first)
			{
				link_count++;
				total_speed += g_link_vector[i].free_speed;
				total_length += g_link_vector[i].length_in_meter * g_link_vector[i].number_of_lanes;
				total_lane_capacity += g_link_vector[i].lane_capacity;
				total_link_capacity += g_link_vector[i].lane_capacity * g_link_vector[i].number_of_lanes;
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
	
		const int fieldWidth = 10;

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
		dtalog.output() << std::setw(fieldWidth) << "avg_lane_cap,";
		g_DTA_log_file << std::setw(fieldWidth) << "avg_lane_cap,";
		dtalog.output() << std::setw(fieldWidth) << "avg_len_meter," << '\n';
		g_DTA_log_file << std::setw(fieldWidth) << "avg_len_meter," << '\n';

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
				if (g_link_vector[i].link_type >= 0 && g_link_vector[i].link_type == it_link_type->first)
				{
					link_count++;
					total_speed += g_link_vector[i].free_speed;
					total_length += g_link_vector[i].length_in_meter * g_link_vector[i].number_of_lanes;
					total_lane_capacity += g_link_vector[i].lane_capacity;
					total_link_capacity += g_link_vector[i].lane_capacity * g_link_vector[i].number_of_lanes;
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
//                g_link_vector[timing_arc.link_seq_no].red_time = timing_arc.red_time;
//                g_link_vector[timing_arc.link_seq_no].cycle_length = timing_arc.cycle_length;
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

//// new Yaml file input

//////////////////////////////////////////
void read_settings_yml(const std::string& file_path)
{


	//// it is possible that no AgentType is set up
	//if (this->ats.empty())
	//	this->ats.push_back(new AgentType());

	
}

void read_settings(const std::string& dir)
{

	read_settings_yml("settings.yml");

	//path file_path = dir + '/' + "settings.yml";
	//if (exists(file_path))
	//	read_settings_yml(file_path.string());
	//else
	//	auto_setup();
}



