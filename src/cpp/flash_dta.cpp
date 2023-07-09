/* Portions Copyright 2019 Xuesong Zhou and Peiheng Li
 *
 * If you help write or modify the code, please also list your names here.
 * The reason of having Copyright info here is to ensure all the modified version, as a whole, under the GPL
 * and further prevent a violation of the GPL.
 *
 * More about "How to use GNU licenses for your own software"
 * http://www.gnu.org/licenses/gpl-howto.html
 */

#ifdef _WIN32
#include "pch.h"
#endif

#include "config.h"
#include "utils.h"
#include <tuple>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <filesystem>

// to do:
// dynamic fluid based analytical based DTA and, agent based simulation, consistency between agent based simulation and dynamic fluid based ADTA
// dynamic fluid based DTA is an extention of static traffic assignment based on steady state
// default district_id
// loading vehicle trajectory file in demand through map matching

using namespace std;

void write_default_setting_file_if_not_exist()
{

	std::string filename = "settings.csv";

	// Check if the file exists
	std::ifstream inFile(filename);
	if (inFile.good()) 
	{
		return;  // If the file exists, we'll stop the program here.
	}

	// Define the data to write to the file.
		std::vector<std::tuple<std::string, std::string, std::string>> data = {
			{"assignment", "number_of_iterations", "20"},
			{"assignment", "route_output", "1"},
			{"assignment", "simulation_output", "0"},
			{"assignment", "UE_convergence_percentage", "0.1"},
			{"cpu", "number_of_memory_blocks", "4"},
			{"unit", "length_unit", "meter"},
			{"unit", "speed_unit", "kmph"},
			{"subarea", "max_num_significant_zones_in_subarea","50000"},
			{"subarea", "max_num_significant_zones_outside_subarea","50000"}
		};

	// Open the file for output.
	std::ofstream outFile("settings.csv");

	// Write the header to the file.
	outFile << "section,key,value\n";

	// Write the data to the file.
	for (const auto& item : data) {
		outFile << std::get<0>(item) << ","
			<< std::get<1>(item) << ","
			<< std::get<2>(item) << "\n";
	}

	// Close the file.
	outFile.close();
	return;
}
bool write_default_scenario_index_file_if_not_exist()
{

	std::string filename = "scenario_index_list.csv";

	// Check if the file exists
	std::ifstream inFile(filename);
	if (inFile.good())
	{
		return 0;  // If the file exists, we'll stop the program here.
	}

	// Define the data to write to the file.
	std::vector<std::tuple<std::string, std::string, std::string, std::string, std::string, std::string>> data = {
		{"0", "0", "2025", "25nb", "2025 no built", "1"},
		{"0", "1","2040", "2040", "2040 future year", "0"},
	};

	// Open the file for output.
	std::ofstream outFile(filename);

	// Write the header to the file.
	outFile << "first_column,scenario_index,year,scenario_name,scenario_description,activate\n";

	// Write the data to the file.
	for (const auto& item : data) {
		outFile << std::get<0>(item) << ","
			<< std::get<1>(item) << ","
			<< std::get<2>(item) << ","
			<< std::get<3>(item) << ","
			<< std::get<4>(item) << ","
			<< std::get<5>(item) << "\n";
	}


	// Close the file.
	outFile.close();
	return 1;
}


int write_default_demand_period_file_if_not_exist() {
	// Define the output file name
	std::string filename = "demand_period.csv";

	// Check if the file exists
	{
		std::ifstream inFile(filename);
		if (inFile.good()) {
			std::cout << "The file " << filename << " already exists.\n";
			return 0;  // If the file exists, we'll stop the program here.
		}
	}  // inFile goes out of scope here and is closed.

 // Define the data to write to the file.
	std::vector<std::tuple<std::string, std::string, std::string, std::string, std::string, std::string>> data = {
		{"0", "1", "am", "weekday", "0700_0800", "730"},
	};

	// Open the file for output.
	std::ofstream outFile(filename);

	// Write the header to the file.
	outFile << "first_column,demand_period_id,demand_period,notes,time_period,peak_time\n";

	// Write the data to the file.
	for (const auto& item : data) {
		outFile << std::get<0>(item) << ","
			<< std::get<1>(item) << ","
			<< std::get<2>(item) << ","
			<< std::get<3>(item) << ","
			<< std::get<4>(item) << ","
			<< std::get<5>(item) << "\n";
	}

	// Close the file.
	outFile.close();

	return 0;
}

void write_default_mode_type_file_if_not_exist() {
	// Define the output file name
	std::string filename = "mode_type.csv";

	// Check if the file exists
	{
		std::ifstream inFile(filename);
		if (inFile.good()) {
			std::cout << "The file " << filename << " already exists.\n";
			return;  // If the file exists, we'll return from the function here.
		}
	}  // inFile goes out of scope here and is closed.

	// Define the data to write to the file.
	std::vector<std::tuple<std::string, std::string, std::string, std::string, std::string, std::string, std::string, std::string, std::string, std::string>> data = {
		{"0", "auto", "0", "auto", "10", "1", "1", "1.5", "0", "1"},
		{"0", "walk", "1", "walk", "10", "1", "1", "1", "0", "0"},
		{"0", "bike", "2", "bike", "10", "1", "1", "1", "0", "0"},
		{"0", "bus", "3", "bus", "10", "1", "10", "0.1", "0", "0"},
		{"0", "truck", "4", "truck", "10", "0", "1", "0.2", "0", "0"},
		{"0", "cav", "5", "cav", "10", "0", "1", "0.75", "1", "0"},
		{"0", "ev", "6", "ev", "10", "0", "1", "1.5", "0", "0"},
		{"0", "hov", "7", "hov", "10", "0", "2", "1.5", "0", "0"},
	};

	// Open the file for output.
	std::ofstream outFile(filename);

	// Write the header to the file.
	outFile << "first_column,mode_type,mode_type_index,name,vot,multimodal_dedicated_assignment_flag,person_occupancy,headway_in_sec,DTM_real_time_info_type,activate\n";

	// Write the data to the file.
	for (const auto& item : data) {
		outFile << std::get<0>(item) << ","
			<< std::get<1>(item) << ","
			<< std::get<2>(item) << ","
			<< std::get<3>(item) << ","
			<< std::get<4>(item) << ","
			<< std::get<5>(item) << ","
			<< std::get<6>(item) << ","
			<< std::get<7>(item) << ","
			<< std::get<8>(item) << ","
			<< std::get<9>(item) << "\n";
	}

	// Close the file.
	outFile.close();
}


struct LinkTypeData {
	std::string first_column;
	std::string link_type;
	std::string link_type_name;
	std::string name_description;
	std::string type_code;
	std::string traffic_flow_model;
	std::string allowed_uses_p1;
	std::string allowed_uses_p2;
	std::string allowed_uses_p3;
	std::string peak_load_factor_p1_auto;
	std::string peak_load_factor_p1_bike;
	std::string peak_load_factor_p1_walk;
	std::string free_speed_auto;
	std::string free_speed_bike;
	std::string free_speed_walk;
	std::string capacity_auto;
	std::string capacity_bike;
	std::string capacity_walk;
	std::string lanes_bike;
	std::string lanes_walk;
	std::string k_jam_km;
	std::string meu_auto_bike;
	std::string meu_auto_walk;
	std::string meu_bike_walk;
	std::string meu_bike_auto;
	std::string meu_walk_bike;
	std::string meu_walk_auto;
	std::string emissions_auto_co2;
	std::string emissions_auto_nox;
	std::string emissions_bike_co2;
	std::string emissions_bike_nox;
	std::string emissions_walk_co2;
	std::string emissions_walk_nox;
	std::string emissions_ev_co2;
	std::string emissions_ev_nox;
	std::string emissions_truck_co2;
	std::string emissions_truck_nox;
	std::string emissions_bus_co2;
	std::string emissions_bus_nox;
	std::string emissions_hov_co2;
	std::string emissions_hov_nox;
};

void write_default_link_type_file_if_not_exist() {
	std::string filename = "link_type.csv";

	std::ifstream inFile(filename.c_str());
	if (inFile.good()) {
		std::cout << "The file " << filename << " already exists.\n";
		return;
	}

	std::vector<LinkTypeData> data;

	// Add each line as a dataItem in the vector data
	// First record
	LinkTypeData dataItem1 = {
		"0", "1", "motorway", "motorway", "f", "kw", "auto", "", "",
		"1", "1", "1", "120", "25", "5", "2000", "300", "200",
		"0", "0", "300", "0", "0", "0", "0", "0", "0",
		"20785.99541;0.0002;0.0042;0.3412", "5.53516;0.0003;0.0043;0.0959",
		"0;0;0;0", "0;0;0;0", "0;0;0;0", "0;0;0;0",
		"10392.99771;0.0002;0.0042;0.3412", "2;0.0003;0.0043;0.0959",
		"23816.14583;0.0002;0.0042;0.3412", "6.342370833;0.0003;0.0043;0.0959",
		"25115.20833;0.0002;0.0042;0.3412", "6.688318333;0.0003;0.0043;0.0959",
		"10392.99771;0.0002;0.0042;0.3412", "2;0.0003;0.0043;0.0959"
	};
	data.push_back(dataItem1);

	// Second record
	LinkTypeData dataItem2 = {
		"0", "2", "trunk", "trunk", "a", "spatial_queue", "auto", "", "",
		"1", "1", "1", "100", "25", "5", "1800", "300", "200",
		"0", "0", "300", "0", "0", "0", "0", "0", "0",
		"20785.99541;0.0002;0.0042;0.3413", "5.53516;0.0003;0.0043;0.0960",
		"0;0;0;0", "0;0;0;0", "0;0;0;0", "0;0;0;0",
		"10392.99771;0.0002;0.0042;0.3412", "2;0.0003;0.0043;0.0960",
		"23816.14583;0.0002;0.0042;0.3413", "6.342370833;0.0003;0.0043;0.0960",
		"25115.20833;0.0002;0.0042;0.3413", "6.688318333;0.0003;0.0043;0.0960",
		"10392.99771;0.0002;0.0042;0.3413", "2;0.0003;0.0043;0.0960"
	};
	data.push_back(dataItem2);
	LinkTypeData dataItem3 = { "0", "3", "primary", "primary", "a", "spatial_queue", "", "", "", "1", "1", "1", "80", "25", "5", "1500", "300", "200", "1", "1", "300", "0", "0", "0", "0", "0", "0", "20785.99541;0.0002;0.0042;0.3414", "5.53516;0.0003;0.0043;0.0961", "0;0;0;0", "0;0;0;0", "0;0;0;0", "0;0;0;0", "10392.99771;0.0002;0.0042;0.3413", "2;0.0003;0.0043;0.0961", "23816.14583;0.0002;0.0042;0.3414", "6.342370833;0.0003;0.0043;0.0961", "25115.20833;0.0002;0.0042;0.3414", "6.688318333;0.0003;0.0043;0.0961", "10392.99771;0.0002;0.0042;0.3414", "2;0.0003;0.0043;0.0961" };
	data.push_back(dataItem3);

	LinkTypeData dataItem4 = { "0", "4", "residential", "residential", "a", "point_queue", "", "", "", "1", "1", "1", "30", "25", "5", "1000", "300", "200", "1", "1", "300", "0", "0", "0", "0", "0", "0", "20785.99541;0.0002;0.0042;0.3415", "5.53516;0.0003;0.0043;0.0962", "0;0;0;0", "0;0;0;0", "0;0;0;0", "0;0;0;0", "10392.99771;0.0002;0.0042;0.3414", "2;0.0003;0.0043;0.0962", "23816.14583;0.0002;0.0042;0.3415", "6.342370833;0.0003;0.0043;0.0962", "25115.20833;0.0002;0.0042;0.3415", "6.688318333;0.0003;0.0043;0.0962", "10392.99771;0.0002;0.0042;0.3415", "2;0.0003;0.0043;0.0962" };
	data.push_back(dataItem4);

	LinkTypeData dataItem5 = { "0", "5", "secondary", "secondary", "a", "point_queue", "", "", "", "1", "1", "1", "60", "25", "5", "1400", "300", "200", "1", "1", "300", "0", "0", "0", "0", "0", "0", "20785.99541;0.0002;0.0042;0.3416", "5.53516;0.0003;0.0043;0.0963", "0;0;0;0", "0;0;0;0", "0;0;0;0", "0;0;0;0", "10392.99771;0.0002;0.0042;0.3415", "2;0.0003;0.0043;0.0963", "23816.14583;0.0002;0.0042;0.3416", "6.342370833;0.0003;0.0043;0.0963", "25115.20833;0.0002;0.0042;0.3416", "6.688318333;0.0003;0.0043;0.0963", "10392.99771;0.0002;0.0042;0.3416", "2;0.0003;0.0043;0.0963" };
	data.push_back(dataItem5);

	LinkTypeData dataItem6 = { "0", "6", "tertiary", "tertiary", "a", "point_queue", "", "", "", "1", "1", "1", "50", "25", "5", "1300", "300", "200", "1", "1", "300", "0", "0", "0", "0", "0", "0", "20785.99541;0.0002;0.0042;0.3417", "5.53516;0.0003;0.0043;0.0964", "0;0;0;0", "0;0;0;0", "0;0;0;0", "0;0;0;0", "10392.99771;0.0002;0.0042;0.3416", "2;0.0003;0.0043;0.0964", "23816.14583;0.0002;0.0042;0.3417", "6.342370833;0.0003;0.0043;0.0964", "25115.20833;0.0002;0.0042;0.3417", "6.688318333;0.0003;0.0043;0.0964", "10392.99771;0.0002;0.0042;0.3417", "2;0.0003;0.0043;0.0964" };
	data.push_back(dataItem6);

	LinkTypeData dataItem7 = { "0", "7", "unclassified", "unclassified", "a", "point_queue", "", "", "", "1", "1", "1", "40", "25", "5", "1200", "300", "200", "1", "1", "300", "0", "0", "0", "0", "0", "0", "20785.99541;0.0002;0.0042;0.3418", "5.53516;0.0003;0.0043;0.0965", "0;0;0;0", "0;0;0;0", "0;0;0;0", "0;0;0;0", "10392.99771;0.0002;0.0042;0.3417", "2;0.0003;0.0043;0.0965", "23816.14583;0.0002;0.0042;0.3418", "6.342370833;0.0003;0.0043;0.0965", "25115.20833;0.0002;0.0042;0.3418", "6.688318333;0.0003;0.0043;0.0965", "10392.99771;0.0002;0.0042;0.3418", "2;0.0003;0.0043;0.0965" };
	data.push_back(dataItem7);

	LinkTypeData dataItem9 = { "0", "9", "collector", "collector", "a", "point_queue", "", "", "", "1", "1", "1", "50", "25", "5", "1200", "300", "200", "1", "1", "300", "0", "0", "0", "0", "0", "0", "20785.99541;0.0002;0.0042;0.3419", "5.53516;0.0003;0.0043;0.0966", "0;0;0;0", "0;0;0;0", "0;0;0;0", "0;0;0;0", "10392.99771;0.0002;0.0042;0.3418", "2;0.0003;0.0043;0.0966", "23816.14583;0.0002;0.0042;0.3419", "6.342370833;0.0003;0.0043;0.0966", "25115.20833;0.0002;0.0042;0.3419", "6.688318333;0.0003;0.0043;0.0966", "10392.99771;0.0002;0.0042;0.3419", "2;0.0003;0.0043;0.0966" };
	data.push_back(dataItem9);


	LinkTypeData dataItem10 = { "0", "10", "connector", "connector", "c", "point_queue", "", "", "", "1", "1", "1", "60", "25", "5", "15000", "300", "200", "1", "1", "300", "0", "0", "0", "0", "0", "0" };
	data.push_back(dataItem10);

	LinkTypeData dataItem100 = { "0", "100", "shared_use", "", "a", "point_queue", "", "", "", "1", "1", "1", "20", "25", "5", "200", "300", "200", "1", "1", "300", "0.5", "0", "0", "2", "0", "0", "20785.99541;0.0002;0.0042;0.3419", "5.53516;0.0003;0.0043;0.0966", "0;0;0;0", "0;0;0;0", "0;0;0;0", "0;0;0;0", "10392.99771;0.0002;0.0042;0.3418", "2;0.0003;0.0043;0.0966", "23816.14583;0.0002;0.0042;0.3419", "6.342370833;0.0003;0.0043;0.0966", "25115.20833;0.0002;0.0042;0.3419", "6.688318333;0.0003;0.0043;0.0966", "10392.99771;0.0002;0.0042;0.3419", "2;0.0003;0.0043;0.0966" };
	data.push_back(dataItem100);

	LinkTypeData dataItem200 = { "0", "200", "bikeonly", "Bike only", "a", "point_queue", "bike", "bike", "bike", "1", "1", "1", "5", "25", "5", "200", "300", "200", "1", "1", "300", "0", "0", "0", "0", "0", "0"};
	data.push_back(dataItem200);

	LinkTypeData dataItem300 = { "0", "300", "walkonly", "walk only", "a", "point_queue", "walk", "walk", "walk", "1", "1", "1", "25", "25", "5", "300", "300", "200", "1", "1", "300", "0", "0", "0", "0", "0", "0"};
	data.push_back(dataItem300);


	std::ofstream outFile(filename.c_str());

	outFile << "first_column,link_type,link_type_name,name_description,type_code,traffic_flow_model,allowed_uses_p1,allowed_uses_p2,allowed_uses_p3,"
		"peak_load_factor_p1_auto,peak_load_factor_p1_bike,peak_load_factor_p1_walk,free_speed_auto,free_speed_bike,free_speed_walk,capacity_auto,capacity_bike,"
		"capacity_walk,lanes_bike,lanes_walk,k_jam_km,meu_auto_bike,meu_auto_walk,meu_bike_walk,meu_bike_auto,meu_walk_bike,meu_walk_auto,"
		"emissions_auto_co2,emissions_auto_nox,emissions_bike_co2,emissions_bike_nox,emissions_walk_co2,emissions_walk_nox,emissions_ev_co2,"
		"emissions_ev_nox,emissions_truck_co2,emissions_truck_nox,emissions_bus_co2,emissions_bus_nox,emissions_hov_co2,emissions_hov_nox\n";

	for (const auto& item : data) {
		outFile << item.first_column << ","
			<< item.link_type << ","
			<< item.link_type_name << ","
			<< item.name_description << ","
			<< item.type_code << ","
			<< item.traffic_flow_model << ","
			<< item.allowed_uses_p1 << ","
			<< item.allowed_uses_p2 << ","
			<< item.allowed_uses_p3 << ","
			<< item.peak_load_factor_p1_auto << ","
			<< item.peak_load_factor_p1_bike << ","
			<< item.peak_load_factor_p1_walk << ","
			<< item.free_speed_auto << ","
			<< item.free_speed_bike << ","
			<< item.free_speed_walk << ","
			<< item.capacity_auto << ","
			<< item.capacity_bike << ","
			<< item.capacity_walk << ","
			<< item.lanes_bike << ","
			<< item.lanes_walk << ","
			<< item.k_jam_km << ","
			<< item.meu_auto_bike << ","
			<< item.meu_auto_walk << ","
			<< item.meu_bike_walk << ","
			<< item.meu_bike_auto << ","
			<< item.meu_walk_bike << ","
			<< item.meu_walk_auto << ","
			<< item.emissions_auto_co2 << ","
			<< item.emissions_auto_nox << ","
			<< item.emissions_bike_co2 << ","
			<< item.emissions_bike_nox << ","
			<< item.emissions_walk_co2 << ","
			<< item.emissions_walk_nox << ","
			<< item.emissions_ev_co2 << ","
			<< item.emissions_ev_nox << ","
			<< item.emissions_truck_co2 << ","
			<< item.emissions_truck_nox << ","
			<< item.emissions_bus_co2 << ","
			<< item.emissions_bus_nox << ","
			<< item.emissions_hov_co2 << ","
			<< item.emissions_hov_nox << "\n";


	}

	outFile.close();
}


struct DemandFileListData {
	std::string first_column;
	std::string file_sequence_no;
	std::string scenario_index_vector;
	std::string file_name;
	std::string demand_period;
	std::string mode_type;
	std::string format_type;
	std::string scale_factor;
	std::string departure_time_profile_no;
	std::string comment;
};
void write_default_demand_file_list_if_not_exist() {
	std::string filename = "demand_file_list.csv";

	std::ifstream inFile(filename.c_str());
	if (inFile.good()) {
		std::cout << "The file " << filename << " already exists.\n";
		return;
	}

	std::vector<DemandFileListData> data;

	DemandFileListData dataItem = { "0", "1", "0", "demand.csv", "am", "auto", "column", "1", "1", "" };
	DemandFileListData dataItem_1 = { "0", "2", "1", "demand.csv", "am", "auto", "column", "2", "1", "" };

	data.push_back(dataItem);
	//data.push_back(dataItem_1);

	// Add more entries in the same way...

	std::ofstream outFile(filename.c_str());

	outFile << "first_column,file_sequence_no,scenario_index_vector,file_name,demand_period,mode_type,format_type,"
		"scale_factor,departure_time_profile_no,comment\n";

	for (const auto& item : data) {
		outFile << item.first_column << ","
			<< item.file_sequence_no << ","
			<< item.scenario_index_vector << ","
			<< item.file_name << ","
			<< item.demand_period << ","
			<< item.mode_type << ","
			<< item.format_type << ","
			<< item.scale_factor << ","
			<< item.departure_time_profile_no << ","
			<< item.comment << "\n";
	}

	outFile.close();
}
bool CheckSensorFileExist()
{
	CDTACSVParser parser_measurement;
	if (parser_measurement.OpenCSVFile("sensor_data.csv", false))
	{

		int count = 0;
		while (parser_measurement.ReadRecord())  // if this line contains [] mark, then we will also read field headers.
		{
			string sensor_type;


			int activate_flag = 0;
			parser_measurement.GetValueByFieldName("activate", activate_flag);
					
			if (activate_flag)
			{
				count++; break;
			}
			

		}

		parser_measurement.CloseCSVFile();

		if (count >= 1)
			return true;
		else
			return false;
	}
	return false;
}

#include <fstream>
#include <vector>

struct DepartureTimeProfileData {
	std::string first_column;
	std::string departure_time_profile_no;
	std::string time_period; // Add time_period
	std::vector<double> time_points;
};

void write_default_departure_time_profile_if_not_exist() {
	std::string filename = "departure_time_profile.csv";

	std::ifstream inFile(filename.c_str());
	if (inFile.good()) {
		std::cout << "The file " << filename << " already exists.\n";
		return;
	}

	std::vector<DepartureTimeProfileData> data;

	DepartureTimeProfileData dataItem = { "0", "1", "0600_0900", { // Include time_period here
0.000571,0.000571,0.000571,0.000571,0.000571,0.000571,0.000506,0.000506,0.000506,0.000445,0.000445,0.000445,0.000391,0.000391,0.000391,0.000357,0.000357,0.000357,0.000328,0.000328,0.000328,0.000319,0.000319,0.000319,0.000302,0.000302,0.000302,0.000292,0.000292,0.000292,0.000296,0.000296,0.000296,0.00031,0.00031,0.00031,0.00031,0.00031,0.00031,0.000319,0.000319,0.000319,0.000383,0.000383,0.000383,0.000496,0.000496,0.000496,0.000568,0.000568,0.000568,0.000656,0.000656,0.000656,0.00095,0.00095,0.00095,0.001368,0.001368,0.001368,0.001587,0.001587,0.001587,0.00175,0.00175,0.00175,0.002288,0.002288,0.002288,0.002921,0.002921,0.002921,0.003242,0.003242,0.003242,0.003218,0.003218,0.003218,0.003803,0.003803,0.003803,0.004459,0.004459,0.004459,0.005002,0.005002,0.005002,0.005207,0.005207,0.005207,0.005677,0.005677,0.005677,0.005994,0.005994,0.005994,0.006018,0.006018,0.006018,0.005508,0.005508,0.005508,0.00529,0.00529,0.00529,0.005058,0.005058,0.005058,0.004833,0.004833,0.004833,0.004421,0.004421,0.004421,0.004327,0.004327,0.004327,0.004364,0.004364,0.004364,0.004343,0.004343,0.004343,0.004139,0.004139,0.004139,0.004201,0.004201,0.004201,0.004291,0.004291,0.004291,0.00435,0.00435,0.00435,0.004409,0.004409,0.004409,0.004566,0.004566,0.004566,0.004674,0.004674,0.004674,0.004761,0.004761,0.004761,0.004827,0.004827,0.004827,0.004882,0.004882,0.004882,0.0049,0.0049,0.0049,0.004887,0.004887,0.004887,0.004835,0.004835,0.004835,0.004899,0.004899,0.004899,0.005023,0.005023,0.005023,0.005065,0.005065,0.005065,0.005162,0.005162,0.005162,0.005436,0.005436,0.005436,0.005772,0.005772,0.005772,0.005907,0.005907,0.005907,0.005877,0.005877,0.005877,0.00605,0.00605,0.00605,0.006196,0.006196,0.006196,0.006248,0.006248,0.006248,0.006308,0.006308,0.006308,0.006404,0.006404,0.006404,0.006391,0.006391,0.006391,0.006401,0.006401,0.006401,0.006526,0.006526,0.006526,0.006574,0.006574,0.006574,0.006271,0.006271,0.006271,0.005937,0.005937,0.005937,0.005578,0.005578,0.005578,0.005293,0.005293,0.005293,0.004834,0.004834,0.004834,0.004387,0.004387,0.004387,0.00403,0.00403,0.00403,0.003748,0.003748,0.003748,0.003382,0.003382,0.003382,0.003121,0.003121,0.003121,0.002963,0.002963,0.002963,0.00289,0.00289,0.00289,0.002671,0.002671,0.002671,0.002468,0.002468,0.002468,0.002365,0.002365,0.002365,0.002249,0.002249,0.002249,0.002015,0.002015,0.002015,0.001784,0.001784,0.001784,0.00164,0.00164,0.00164,0.001474,0.001474,0.001474,0.001312,0.001312,0.001312,0.001132,0.001132,0.001132,0.001005,0.001005,0.001005,0.000889,0.000889,0.000889,0.000778,0.000778,0.000778,0.000676

		//...continue for all time points
	} };
	data.push_back(dataItem);

	// Add more entries in the same way...

	std::ofstream outFile(filename.c_str());

	outFile << "first_column,departure_time_profile_no,time_period,"; // Include time_period in the header
	for (int i = 0; i < 288; ++i) {
		char buffer[10];
		sprintf(buffer, "T%04d", i * 5);
		outFile << buffer;
		if (i != 287)
			outFile << ",";
	}
	outFile << "\n";

	for (const auto& item : data) {
		outFile << item.first_column << ","
			<< item.departure_time_profile_no << ","
			<< item.time_period << ","; // Include time_period in the data writing
		for (size_t i = 0; i < item.time_points.size(); ++i) {
			outFile << item.time_points[i];
			if (i != item.time_points.size() - 1)
				outFile << ",";
		}
		outFile << "\n";
	}

	outFile.close();
}

void write_default_subarea_file_if_not_exist() {
	std::string filename = "subarea.csv";
	std::ifstream inFile(filename.c_str());

	// If the file does not exist, create it
	if (!inFile.good()) {
		std::ofstream outFile;
		outFile.open(filename.c_str());

		// Write the headers
		outFile << "notes,geometry\n";

		// Write the data
		outFile << "subarea_polygon,\"POLYGON ((-180 -90, 180 -90, 180 90, -180 90, -180 -90))\"\n";

		// Close the file
		outFile.close();

		std::cout << "The file " << filename << " was generated successfully!\n";
	}
	else {
		std::cout << "The file " << filename << " already exists.\n";
	}

	inFile.close();
}

bool CheckSupplySideScenarioFileExist()
{

	CDTACSVParser parser;

	int sa_count = 0;

	if (parser.OpenCSVFile("dynamic_traffic_management.csv", false))
	{
		while (parser.ReadRecord())
		{

			int activate_flag = 0;
			parser.GetValueByFieldName("activate", activate_flag);

			if(activate_flag==1)
			{

			int from_node_id = 0;
			parser.GetValueByFieldName("from_node_id", from_node_id);

			int to_node_id = 0;
			parser.GetValueByFieldName("to_node_id", to_node_id);

			if (from_node_id == 0 && to_node_id == 0)
			{
				continue;
			}

			sa_count++;
			}
		}
		parser.CloseCSVFile();
	}

	if (sa_count == 0)
		return false;
	else
		return true;
}

void write_default_dynamic_traffic_management_file_if_not_exist() {
	std::string filename = "dynamic_traffic_management.csv";

	std::ifstream inFile(filename.c_str());
	if (inFile.good()) {
		std::cout << "The file " << filename << " already exists.\n";
		return;
	}

	std::ofstream outFile;
	outFile.open(filename.c_str());

	outFile << "dtm_id,dtm_type,from_node_id,to_node_id,final_lanes,demand_period,mode_type,scenario_index_vector,activate\n";
	outFile << "1,lane_closure,1,3,0.5,am,info,0,0\n";

	outFile.close();
	std::cout << "The file " << filename << " has been created with default values.\n";
}
void write_default_sensor_data_file_if_not_exist() {
	std::string filename = "sensor_data.csv";
	std::ifstream inFile(filename.c_str());

	// If the file does not exist, create it
	if (!inFile.good()) {
		std::ofstream outFile;
		outFile.open(filename.c_str());

		// Write the headers
		outFile << "sensor_id,from_node_id,to_node_id,demand_period,count,scenario_index,activate\n";

		// Write the data
		outFile << "1,483,481,am,3000.975,0,0\n";

		// Close the file
		outFile.close();

		std::cout << "The file " << filename << " was generated successfully!\n";
	}
	else {
		std::cout << "The file " << filename << " already exists.\n";
	}

	inFile.close();
}

int main()
{
	std::ofstream::sync_with_stdio(false);
	// reset all the log files to defult 0: not output; if want to output these logs set to 1



	dtalog.debug_level() = 0;
	dtalog.log_sig() = 0;
	dtalog.log_odme() = 0;
	dtalog.log_path() = 0;
	dtalog.log_dta() = 0;
	dtalog.log_ue() = 0;

	int column_generation_iterations = 0;
	int column_updating_iterations = 40;
	int ODME_iterations = 0;
	int sensitivity_analysis_iterations = 0;
	int number_of_memory_blocks = 4;
	float info_updating_freq_in_min = 5;
	int simulation_output = 0;
	int max_num_significant_zones_in_subarea = 50000;
	int max_num_significant_zones_outside_subarea = 50000;



	double UE_convergence_percentage = 1.0; 

	int signal_updating_output = 0;
	// generate link performance and agent file
	int assignment_mode = 1;
	bool flag_default = false;
	int default_volume = 1;
	int length_unit_flag = 0; //0: meter, 1: mile,
	int speed_unit_flag = 0;  //0: kmph, 1: mph

	dtalog.output() << "Logbook for DTALite: The Open-Source, Lightweight Dynamic Traffic Assignment Solution" << std::fixed << std::setw(12) << '\n';
	dtalog.output() << " Overview of files and process\n";
	dtalog.output() << "1. Input Files:\n";
	dtalog.output() << "   |--- Physical Layer (node.csv, link.csv)\n";
	dtalog.output() << "   |--- Demand Layer (demand.csv, mode_type.csv, demand_period.csv, departure_time_profile.csv, demand_file_list.csv, sensor_data.csv, subarea.csv)\n";
	dtalog.output() << "   |--- Configuration Files (settings.csv, scenario_index_list.csv)\n";
	dtalog.output() << "   |--- Supply Layer (link_type.csv, dynamic_traffic_management.csv)\n";

	dtalog.output() << "\n2. Traffic Assignment and Simulation Process:\n";
	dtalog.output() << "   |--- Demand estimation based on sensor data\n";
	dtalog.output() << "   |--- Traffic assignment based on network and demand data\n";
	dtalog.output() << "   |--- Simulation of traffic based on assignment results and scenario configurations\n";
	dtalog.output() << "   |--- Performance evaluation based on simulation results and performance criteria\n";

	dtalog.output() << "\n3. Output Files:\n";
	dtalog.output() << "   |--- Link performance (link_performance_s.csv, link_performance_summary.csv)\n";
	dtalog.output() << "   |--- Route assignment (route_assignment_s.csv)\n";
	dtalog.output() << "   |--- OD pair and district performance (od_performance_summary.csv, district_performance_s.csv)\n";
	dtalog.output() << "   |--- Trajectory performance (agent_s.csv, trajectory.csv)\n";
	dtalog.output() << "   |--- System performance (system_performance_summary.csv, final_summary.csv)\n";
	dtalog.output() << "   |--- Logs and subarea mapping summary(log_main.txt, log_label_correcting, zonal_hierarchy_mapping.csv)\n";
	dtalog.output() << "--------------------------" << '\n';
	dtalog.output() << "Please provide feedback or report any issues you encounter on our GitHub site: "
		<< "https://github.com/asu-trans-ai-lab/DTALite/issues. Your input helps us enhance the software, address any concerns, and contribute to the open-source transportation ecosystem." << '\n';

	//dtalog.output() << "Input Files:" << '\n';
	//dtalog.output() << "  Physical layer:" << '\n';
	//dtalog.output() << "    node.csv: Defines nodes in the network." << '\n';
	//dtalog.output() << "    link.csv: Defines links in the network with essential attributes for assignment." << '\n';
	//dtalog.output() << "    zone.csv: Optional, as zone_id can be defined in node.csv." << '\n';
	//dtalog.output() << "  Demand layer:" << '\n';
	//dtalog.output() << "    demand.csv: Defines the demand of passengers on each OD pair. This information could be extracted by demand_file_list.csv." << '\n';
	//dtalog.output() << "    demand_period.csv: Defines demand period, which could be extracted by demand_file_list.csv." << '\n';
	//dtalog.output() << "    departure_time_profile.csv: Defines departure time in the agent-based simulation." << '\n';
	//dtalog.output() << "    demand_file_list.csv: Defines demand type, period, and format type." << '\n';
	//dtalog.output() << "    sensor_data.csv: Contains observed link volume for OD demand estimation." << '\n';
	//dtalog.output() << "    choice_set.csv: Contains choice set data for agent-based modeling." << '\n';
	//dtalog.output() << "    activity_travel_pattern.csv: (Optional) Defines activity and travel patterns of agents in the simulation." << '\n';
	//dtalog.output() << "  Supply layer:" << '\n';
	//dtalog.output() << "    dynamic_traffic_management.csv: Defines different dynamic traffic managementscenarios." << '\n';
	//dtalog.output() << "    signal_timing.csv: Contains information about signal timings at intersections." << '\n';
	//dtalog.output() << "  Configuration files:" << '\n';
	//dtalog.output() << "    settings.csv: Defines basic setting for the network, the number of iterations, etc." << '\n';
	//dtalog.output() << "    mode_type.csv: Defines attributes of each type of agent, including value of time (vot in dollars per hour) and passenger car equivalent (pce)." << '\n';
	//dtalog.output() << "    link_type.csv: Defines types of links in the network." << '\n';
	//dtalog.output() << "    link_vdf.csv: Contains analytical volume demand function parameters." << '\n';
	//dtalog.output() << "  Scenarios settings:" << '\n';
	//dtalog.output() << "    scenario_index_list.csv: Defines scenario name, scenario description and activate state." << '\n';
	//dtalog.output() << "    subarea.csv: optional, extracts the subarea polygon information using NeXTA tool." << '\n';
	//dtalog.output() << "--------------------------" << '\n';

	//dtalog.output() << "Output Files:" << '\n';
	//dtalog.output() << "  link_performance_s(scenario_index)_(scenario_name).csv: Shows the performance of each link under different scenarios, including the travel time, volume, and resource balance." << '\n';
	//dtalog.output() << "  route_assignment_s(scenario_index)_(scenario_name).csv: Shows the results of the assignment under different scenarios, including the volume, toll, travel time and distance of each path of each agent, as well as the link sequence and time sequence." << '\n';
	//dtalog.output() << "  choice_set_output_(scenario_index)_(scenario_name).csv: Shows the results of activity travel and mode choice." << '\n';
	//dtalog.output() << "  od_performance_summary.csv: Shows the performance of the OD pairs, including the o_zone_id, d_zone_id and volume." << '\n';
	//dtalog.output() << "  link_performance_summary.csv: Shows the summary of the performance of each link." << '\n';
	//dtalog.output() << "  system_performance_summary.csv: Shows the performance of the whole transportation system, including total travel time, average distance, and total distance." << '\n';
	//dtalog.output() << "  final_summary.csv: Shows a comprehensive summary of the output." << '\n';
	//dtalog.output() << "  subarea_related_zone.csv: Shows the subarea internal zones and impacted zones." << '\n';
	//dtalog.output() << "--------------------------" << '\n';
	write_default_setting_file_if_not_exist();
	write_default_scenario_index_file_if_not_exist();
	write_default_demand_period_file_if_not_exist();
	write_default_mode_type_file_if_not_exist();
	write_default_link_type_file_if_not_exist();
	write_default_demand_file_list_if_not_exist();
	write_default_departure_time_profile_if_not_exist();
	write_default_subarea_file_if_not_exist();
	write_default_sensor_data_file_if_not_exist();
	write_default_dynamic_traffic_management_file_if_not_exist();
	CDTACSVParser parser_settings;

	parser_settings.IsFirstLineHeader = true;
	if (parser_settings.OpenCSVFile("settings.csv", true))
	{

		dtalog.output() << "[PROCESS INFO] Step 0.0: Reading settings.csv." << '\n';
		
		std::string assignment_mode_str;

		int number_of_iterations=  1;
		if (parser_settings.GetValueByKeyName("number_of_iterations", number_of_iterations, false) == true)
			column_generation_iterations = number_of_iterations; // update

		dtalog.output() << "[DATA INFO] number_of_iterations = " << number_of_iterations << " in settings.csv." << '\n';
		// these are the assignment modes
		// two usually methods are ue (user equilibrium) and dta (dynamic traffic assignment)
		// the main difference of these two methods are different output in link_performance.csv
		// for basic uses set assignment mode to 'ue'
		// for more detailed link performances (one minute) set 'dta'1
		assignment_mode = 1;

		if (CheckSensorFileExist() == true)
		{
			column_updating_iterations = 5;
			ODME_iterations = 50;
		}
		else
			ODME_iterations = 0;


		if (CheckSupplySideScenarioFileExist() == true)
			sensitivity_analysis_iterations = 10;
		else
			sensitivity_analysis_iterations = -1;

		int route_output = 1;
		int route_output_value = -1;

		int s = -1;

	
			if (parser_settings.GetValueByKeyName("UE_convergence_percentage", UE_convergence_percentage, false, false) == true)
			{
				if (UE_convergence_percentage < -0.1)
					UE_convergence_percentage = 1;
			}

		dtalog.output() << "[DATA INFO] UE_convergence_percentage = " << UE_convergence_percentage << " (%) in settings.csv." << '\n';


		simulation_output = 0;  //default
		int simulation_output_value = -1;
		if (parser_settings.GetValueByKeyName("simulation_output", simulation_output_value, false, false))
		{
			simulation_output = simulation_output_value;
			if (simulation_output == 1)
				assignment_mode = 2;
		}

		dtalog.output() << "[DATA INFO] simulation_output = " << simulation_output << " in settings.csv." << '\n';


		// the start interation of generating signals, if there is no signals set this number larger than the itertion number
		int number_of_memory_blocks_values = 1;

		if(parser_settings.GetValueByKeyName("number_of_memory_blocks", number_of_memory_blocks_values, false, false))
		{
			number_of_memory_blocks = number_of_memory_blocks_values;
			dtalog.output() << "[DATA INFO] number_of_memory_blocks = " << number_of_memory_blocks << " in settings.csv." << '\n';
		}

		if (parser_settings.GetValueByKeyName("max_num_significant_zones_in_subarea", max_num_significant_zones_in_subarea, false, false))
		{
			dtalog.output() << "[DATA INFO] max_num_significant_zones_in_subarea = " << max_num_significant_zones_in_subarea << " in settings.csv." << '\n';
		}
		

		if (parser_settings.GetValueByKeyName("max_num_significant_zones_outside_subarea", max_num_significant_zones_outside_subarea, false, false))
		{
			dtalog.output() << "[DATA INFO] max_num_significant_zones_outside_subarea = " << max_num_significant_zones_outside_subarea << " in settings.csv." << '\n';
		}	



		string length_unit_str; 
		if (parser_settings.GetValueByKeyName("length_unit", length_unit_str, false, false))
		{
			dtalog.output() << "length_unit = " << length_unit_str.c_str() << " in settings.csv." << '\n';

			if (length_unit_str == "mile")
				length_unit_flag = 1;
			else
				length_unit_flag = 0; // always as default 0

		}


		string speed_unit_str;
		if (parser_settings.GetValueByKeyName("speed_unit", speed_unit_str, false, false))
		{
			dtalog.output() << "speed_unit = " << speed_unit_str.c_str() << " in settings.csv." << '\n';

			if (speed_unit_str == "mph")
				speed_unit_flag = 1;
			else
				speed_unit_flag = 0; // always as default 0 

		}


		

		double number_of_seconds_per_interval_input = 0.25;
		//parser_settings.GetValueByFieldName("number_of_seconds_per_interval", number_of_seconds_per_interval_input, false, false);

		//if (number_of_seconds_per_interval_input > 0.0001)
		//	number_of_seconds_per_interval = number_of_seconds_per_interval_input;

		if (parser_settings.SectionName == "[log]")
		{
			parser_settings.GetValueByFieldName("sig", dtalog.log_sig(), false);
			parser_settings.GetValueByFieldName("odme", dtalog.log_odme(), false);
			parser_settings.GetValueByFieldName("path", dtalog.log_path(), false);
			parser_settings.GetValueByFieldName("ue", dtalog.log_ue(), false);
		}
	}

	// scenario
	//
	// obtain initial flow values
	network_assignment(assignment_mode, column_generation_iterations, column_updating_iterations, ODME_iterations, sensitivity_analysis_iterations, simulation_output, number_of_memory_blocks, length_unit_flag, speed_unit_flag, UE_convergence_percentage, max_num_significant_zones_in_subarea, max_num_significant_zones_outside_subarea);

	return 0;
}

void DTALiteAPI()
{
	main();
}