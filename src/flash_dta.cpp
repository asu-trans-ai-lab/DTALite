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


#ifdef _WIN32
#define YAML_CPP_STATIC_DEFINE
#endif
#include <yaml-cpp/yaml.h>

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
std::ofstream  g_DTA_log_file;
void write_default_setting_file_if_not_exist()
{

	std::string filename = "settings.yml";

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
			{"cpu", "number_of_cpu_processors", "4"},
			{"unit", "length_unit", "meter"},
			{"unit", "speed_unit", "kmph"},
			{"assignment", "odme_activate","1"},
		};

	// Open the file for output.
	std::ofstream outFile("settings.yml");

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



#include <fstream>
#include <vector>

struct DepartureTimeProfileData {
	std::string first_column;
	std::string departure_time_profile_no;
	std::string time_period; // Add time_period
	std::vector<double> time_points;
};




bool CheckSupplySideScenarioFileExist(YAML::Node config)
{

		// Access the "sensor_data" section
		const YAML::Node& dynamic_traffic_management_data_node = config["dynamic_traffic_management_data"];

		// Check if "dynamic_traffic_management_data" is a sequence
		if (dynamic_traffic_management_data_node.IsSequence())
		{
			// Iterate over the sequence and read each entry
			for (const auto& data : dynamic_traffic_management_data_node)
			{
				int activate = data["activate"].as<int>(0);
				if (activate == 1)
					return true;

			}

		}
		return false;


}


int main()
{
	std::ofstream::sync_with_stdio(false);
	// reset all the log files to defult 0: not output; if want to output these logs set to 1

	g_DTA_log_file.open("log_DTA.txt");

	dtalog.debug_level() = 0;
	dtalog.log_sig() = 0;
	dtalog.log_odme() = 0;
	dtalog.log_path() = 0;
	dtalog.log_dta() = 0;
	dtalog.log_ue() = 0;

	int column_generation_iterations = 20;
	int column_updating_iterations = 40;
	int ODME_iterations = 0;
	int sensitivity_analysis_iterations = 0;
	int number_of_cpu_processors = 4;
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
	dtalog.output() << "   |--- Demand Layer (demand.csv, mode_type, departure_time_profile, subarea)\n";
	dtalog.output() << "   |--- Supply Layer (link_type)\n";
	dtalog.output() << "   |--- Configuration Files (settings.yml, scenario)\n";

	dtalog.output() << "\n2. Traffic Assignment and Simulation Process:\n";
	dtalog.output() << "   |--- Traffic assignment based on network and demand data\n";
	dtalog.output() << "   |--- OD demand matrix estimation based on sensor data\n";
	dtalog.output() << "   |--- Simulation of traffic based on assignment results and scenario configurations\n";
	dtalog.output() << "   |--- Performance evaluation based on simulation results and performance criteria\n";

	dtalog.output() << "\n3. Output Files:\n";
	dtalog.output() << "   |--- Link performance (link_performance.csv, link_performance_summary.csv)\n";
	dtalog.output() << "   |--- Route assignment (route_assignment.csv)\n";
	dtalog.output() << "   |--- OD pair and district performance (od_performance_summary.csv, district_performance_s.csv)\n";
	dtalog.output() << "   |--- Trajectory performance (agent.csv, trajectory.csv)\n";
	dtalog.output() << "   |--- System performance (system_performance_summary.csv)\n";
	dtalog.output() << "   |--- Logs and subarea mapping summary(log_main.txt, log_label_correcting, zone_mapping.csv)\n";
	dtalog.output() << "--------------------------" << '\n';
	dtalog.output() << "Please provide feedback or report any issues you encounter on our GitHub site: "
		<< "https://github.com/asu-trans-ai-lab/DTALite/issues. Your input helps us enhance the software, address any concerns, and contribute to the open-source transportation ecosystem." << '\n';

	g_DTA_log_file << "Logbook for DTALite: The Open-Source, Lightweight Dynamic Traffic Assignment Solution" << std::fixed << std::setw(12) << '\n';
	g_DTA_log_file << " Overview of files and process\n";
	g_DTA_log_file << "1. Input Files:\n";
	g_DTA_log_file << "   |--- Physical Layer (node.csv, link.csv)\n";
	g_DTA_log_file << "   |--- Demand Layer (demand.csv, mode_type, departure_time_profile, subarea)\n";
	g_DTA_log_file << "   |--- Supply Layer (link_type)\n";
	g_DTA_log_file << "   |--- Configuration Files (settings.yml, scenario)\n";


	g_DTA_log_file << "\n2. Traffic Assignment and Simulation Process:\n";
	g_DTA_log_file << "   |--- Traffic assignment based on network and demand data\n";
	g_DTA_log_file << "   |--- OD demand Demand estimation based on sensor data\n";
	g_DTA_log_file << "   |--- Simulation of traffic based on assignment results and scenario configurations\n";
	g_DTA_log_file << "   |--- Performance evaluation based on simulation results and performance criteria\n";

	g_DTA_log_file << "\n3. Output Files:\n";
	g_DTA_log_file << "   |--- Link performance (link_performance.csv, link_performance_summary.csv)\n";
	g_DTA_log_file << "   |--- Route assignment (route_assignment.csv)\n";
	g_DTA_log_file << "   |--- OD pair and district performance (od_performance_summary.csv, district_performance_s.csv)\n";
	g_DTA_log_file << "   |--- Trajectory performance (agent.csv, trajectory.csv)\n";
	g_DTA_log_file << "   |--- System performance (system_performance_summary.csv, final_summary.csv)\n";
	g_DTA_log_file << "   |--- Logs and subarea mapping summary(log_main.txt, log_label_correcting, zone_mapping.csv)\n";
	g_DTA_log_file << "--------------------------" << '\n';
	g_DTA_log_file << "Please provide feedback or report any issues you encounter on our GitHub site: "
		<< "https://github.com/asu-trans-ai-lab/DTALite/issues. Your input helps us enhance the software, address any concerns, and contribute to the open-source transportation ecosystem." << '\n';


	g_DTA_log_file << "Input Files:" << '\n';
	g_DTA_log_file << "  Physical layer:" << '\n';
	g_DTA_log_file << "    node.csv: Defines nodes in the network." << '\n';
	g_DTA_log_file << "    link.csv: Defines links in the network with essential attributes for assignment." << '\n';
	g_DTA_log_file << "  Configuration files:" << '\n';
	g_DTA_log_file << "    settings.yml: Defines basic setting for the network, the number of iterations, etc." << '\n';

	g_DTA_log_file << "  Demand layer:" << '\n';
	g_DTA_log_file << "    demand.csv: Defines the demand of passengers on each OD pair. This information could be extracted by demand_file_list.csv." << '\n';
	g_DTA_log_file << "    demand_period : Defines demand period, which could be extracted by demand_file_list.csv." << '\n';
	g_DTA_log_file << "    departure_time_profile: Defines departure time in the agent-based simulation." << '\n';
	g_DTA_log_file << "    demand_file_list: Defines demand type, period, and format type." << '\n';
	g_DTA_log_file << "    sensor_data : Contains observed link volume for OD demand estimation." << '\n';
	//g_DTA_log_file << "    choice_set.csv: Contains choice set data for agent-based modeling." << '\n';
	//g_DTA_log_file << "    activity_travel_pattern.csv: (Optional) Defines activity and travel patterns of agents in the simulation." << '\n';
	g_DTA_log_file << "  Supply layer:" << '\n';
	g_DTA_log_file << "    dynamic_traffic_management: Defines different dynamic traffic management scenarios." << '\n';
	g_DTA_log_file << "    signal_timing  which contains information about signal timings at intersections, coded in link.csv." << '\n';
	g_DTA_log_file << "    mode_type: Defines attributes of each type of agent, including value of time (vot in dollars per hour) and passenger car equivalent (pce)." << '\n';
	g_DTA_log_file << "    link_type: Defines types of links in the network." << '\n';
	g_DTA_log_file << "    link_qvdf.csv: Contains analytical volume demand function parameters." << '\n';
	g_DTA_log_file << "  Scenarios settings:" << '\n';
	g_DTA_log_file << "    scenario_index_list: Defines scenario name, scenario description and activate state." << '\n';
	g_DTA_log_file << "    subarea: extracts the subarea polygon information using NeXTA tool." << '\n';
	g_DTA_log_file << "--------------------------" << '\n';

	g_DTA_log_file << "Output Files:" << '\n';
	g_DTA_log_file << "  link_performance_s(scenario_index)_(scenario_name).csv: Shows the performance of each link under different scenarios, including the travel time, volume, and resource balance." << '\n';
	g_DTA_log_file << "  route_assignment_s(scenario_index)_(scenario_name).csv: Shows the results of the assignment under different scenarios, including the volume, toll, travel time and distance of each path of each agent, as well as the link sequence and time sequence." << '\n';
	//g_DTA_log_file << "  choice_set_output_(scenario_index)_(scenario_name).csv: Shows the results of activity travel and mode choice." << '\n';
	g_DTA_log_file << "  od_performance_summary.csv: Shows the performance of the OD pairs, including the o_zone_id, d_zone_id and volume." << '\n';
	g_DTA_log_file << "  link_performance_summary.csv: Shows the summary of the performance of each link." << '\n';
	g_DTA_log_file << "  system_performance_summary.csv: Shows the performance of the whole transportation system, including total travel time, average distance, and total distance." << '\n';
	g_DTA_log_file << "  final_summary.csv: Shows a comprehensive summary of the output." << '\n';
	g_DTA_log_file << "  internal_zone_mapping.csv: Shows the subarea internal zones and impacted zones." << '\n';
	g_DTA_log_file << "--------------------------" << '\n';
	write_default_setting_file_if_not_exist();

	bool bCorrectSettingFormat = false; 

	std::string fileName = "settings.yml";
	std::ifstream inFile(fileName.c_str());
	YAML::Node settings;
	if (!inFile.is_open()) {
		dtalog.output() << "Error opening file: " << fileName << std::endl;
		return 0;
	}

	try {
		settings = YAML::Load(inFile);
		// Now you can work with the 'settings' as a YAML::Node object
	}
	catch (const YAML::ParserException& e) {
		dtalog.output() << "Error parsing the file: " << e.what() << std::endl;
		return 0;
	}
		dtalog.output() << "[PROCESS INFO] Step 0.0: Reading settings.yml." << '\n';
		g_DTA_log_file << "[PROCESS INFO] Step 0.0: Reading settings.yml." << '\n';

		int number_of_iterations = settings["assignment"]["number_of_iterations"].as<int>(1);
		UE_convergence_percentage =  settings["assignment"]["UE_convergence_percentage"].as<float>(0.1);
		if (UE_convergence_percentage < -0.1)
			UE_convergence_percentage = 1;

		dtalog.output() << "[DATA INFO] UE_convergence_percentage = " << UE_convergence_percentage << " (%) in settings.yml." << '\n';
		g_DTA_log_file << "[DATA INFO] UE_convergence_percentage = " << UE_convergence_percentage << " (%) in settings.yml." << '\n';

		int number_of_column_updating_iterations = settings["assignment"]["number_of_column_updating_iterations"].as<int>(20);

		assignment_mode = 0;

		int route_output = 1;

		route_output = settings["assignment"]["route_output"].as<int>(0);
		if (route_output == 1)
			assignment_mode = 1;

		if (column_updating_iterations < 1)
			column_updating_iterations = 1;

		dtalog.output() << "[DATA INFO] number_of_column_updating_iterations = " << column_updating_iterations << " in settings.yml." << '\n';
		g_DTA_log_file << "[DATA INFO] number_of_column_updating_iterations = " << column_updating_iterations << "  in settings.yml." << '\n';

		simulation_output = settings["assignment"]["simulation_output"].as<int>(0);
			if (simulation_output == 1)
				assignment_mode = 2;
		

		dtalog.output() << "[DATA INFO] simulation_output = " << simulation_output << " in settings.yml." << '\n';
		g_DTA_log_file << "[DATA INFO] simulation_output = " << simulation_output << " in settings.yml." << '\n';

		// Accessing nested elements
		number_of_cpu_processors = settings["cpu"]["number_of_cpu_processors"].as<int>(4);
		dtalog.output() << "[DATA INFO] number_of_cpu_processors = " << number_of_cpu_processors << " in settings.yml." << '\n';
		g_DTA_log_file << "[DATA INFO] number_of_cpu_processors = " << number_of_cpu_processors << " in settings.yml." << '\n';


		column_generation_iterations = number_of_iterations; // update

		dtalog.output() << "[DATA INFO] number_of_iterations = " << number_of_iterations << " in settings.yml." << '\n';
		g_DTA_log_file << "[DATA INFO] number_of_iterations = " << number_of_iterations << " in settings.yml." << '\n';
		// these are the assignment modes
		// two usually methods are ue (user equilibrium) and dta (dynamic traffic assignment)
		// the main difference of these two methods are different output in link_performance.csv
		// for basic uses set assignment mode to 'ue'
		// for more detailed link performances (one minute) set 'dta'1

		int odme_activate = settings["assignment"]["odme_activate"].as<int>(0);

		if (odme_activate==1)
		{
			route_output = 1;
			assignment_mode = 1;
			column_updating_iterations = 5;
			ODME_iterations = 50;
		}
		else
			ODME_iterations = 0;


		if (CheckSupplySideScenarioFileExist(settings) == true)
		{

			sensitivity_analysis_iterations = 20;
			//int sensitivity_analysis_iterations_value = -1;
			//if (parser_settings.GetValueByKeyName("sensitivity_analysis_iterations", sensitivity_analysis_iterations_value, false, false))
			//{
			//	sensitivity_analysis_iterations = sensitivity_analysis_iterations_value;
			//}
		}





			std::string length_unit_str = settings["unit"]["length_unit"].as<std::string>("meter");
	
			dtalog.output() << "length_unit = " << length_unit_str.c_str() << " in settings.yml." << '\n';
			g_DTA_log_file << "length_unit = " << length_unit_str.c_str() << " in settings.yml." << '\n';

			if (length_unit_str == "mile")
				length_unit_flag = 1;
			else if (length_unit_str == "km")
				length_unit_flag = 2;
			else if(length_unit_str == "meter")
				length_unit_flag = 0; // always as default 0
			else
			{
				if (length_unit_str.size() > 0)
				{
					dtalog.output() << "[ERROR] length_unit = " << length_unit_str.c_str() << " in settings.yml  is not supported. The supported unit of length is meter, km or mile.." << '\n';
					g_DTA_log_file << "[ERROR] length_unit = " << length_unit_str.c_str() << " in settings.yml  is not supported. The supported unit of length is meter, km or mile.." << '\n';
					length_unit_flag = 0;
				}
			}

			
			

			std::string speed_unit_str = settings["unit"]["speed_unit"].as<std::string>("kmph");



			dtalog.output() << "speed_unit = " << speed_unit_str.c_str() << " in settings.yml." << '\n';
			g_DTA_log_file << "speed_unit = " << speed_unit_str.c_str() << " in settings.yml." << '\n';

			if (speed_unit_str == "mph")
				speed_unit_flag = 1;
			else
				speed_unit_flag = 0; // always as default 0 



		

		double number_of_seconds_per_interval_input = 0.25;
		//parser_settings.GetValueByFieldName("number_of_seconds_per_interval", number_of_seconds_per_interval_input, false, false);

		//if (number_of_seconds_per_interval_input > 0.0001)
		//	number_of_seconds_per_interval = number_of_seconds_per_interval_input;

		//if (parser_settings.SectionName == "[log]")
		//{
		//	parser_settings.GetValueByFieldName("sig", dtalog.log_sig(), false);
		//	parser_settings.GetValueByFieldName("odme", dtalog.log_odme(), false);
		//	parser_settings.GetValueByFieldName("path", dtalog.log_path(), false);
		//	parser_settings.GetValueByFieldName("ue", dtalog.log_ue(), false);
		//}
	

	// scenario
	//
	// obtain initial flow values
	network_assignment(assignment_mode, column_generation_iterations, column_updating_iterations, ODME_iterations, sensitivity_analysis_iterations, simulation_output, number_of_cpu_processors, length_unit_flag, speed_unit_flag, UE_convergence_percentage, max_num_significant_zones_in_subarea, max_num_significant_zones_outside_subarea);

	if (g_DTA_log_file.is_open())
		g_DTA_log_file.close(); 
	return 0;
}

void DTALiteAPI()
{
	main();
}