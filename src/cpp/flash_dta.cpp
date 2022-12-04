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

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

void write_default_setting_file_if_not_exist()
{

	CDTACSVParser parser_settings;
	parser_settings.IsFirstLineHeader = false;

	if (parser_settings.OpenCSVFile("settings.csv", false))
	{
		parser_settings.CloseCSVFile();
		return;
	}

	ofstream myfile;
	myfile.open("settings.csv");
	myfile << "[assignment],number_of_iterations,route_output,simulation_output" << endl;
	myfile << ",1,1,0," << endl;
	myfile << "[agent_type],agent_type_no,agent_type,name,display_code,vot,flow_type,pce,person_occupancy,desired_speed_ratio,headway,real_time_info," << endl;
	myfile << ",1,auto,auto,auto,10,0,1,1,1,1,0," << endl;
	//myfile << ",2,walk,walk,walk,10,0,1,1,0.1,1,0," << endl;
	//myfile << ",3,bike,bike,bike,10,0,,1,0.2,1,0" << endl;
	//myfile << ",4,bus,bus,bus,10,0,1,10,0.75,5,0" << endl;
	//myfile << ",5,truck,truck,truck,10,0,1,1,1,3,0" << endl;
	//myfile << ",6,cav,auto,auto,10,0,1,1,1,1,0," << endl;
	//myfile << ",7,ev,ev,ev,10,0,1,1,1,1,0," << endl;
	myfile << "" << endl;
	myfile << "[link_type],link_type,link_type_name,type_code,traffic_flow_code,vdf_type," << endl;
	myfile << ",1,motorway,,f,0,qvdf," << endl;
	myfile << ",2,trunk,,a,0,qvdf," << endl;
	myfile << ",3,primary,,a,0,bpr," << endl;
	myfile << ",4,residential,,a,0,bpr," << endl;
	myfile << ",5,secondary,,a,0,bpr," << endl;
	myfile << ",6,tertiary,,a,0,bpr," << endl;
	myfile << ",7,unclassified,,a,0,bpr," << endl;
	myfile << "[demand_period],demand_period_id,demand_period,,time_period" << endl;
	myfile << ", 1,AM, ,0800_0900,, , , " << endl;
	myfile << "[departure_time_profile], departure_time_profile_no, , ,time_period, , , ,T0000,T0005,T0010,T0015,T0020,T0025,T0030,T0035,T0040,T0045,T0050,T0055,T0100,T0105,T0110,T0115,T0120,T0125,T0130,T0135,T0140,T0145,T0150,T0155,T0200,T0205,T0210,T0215,T0220,T0224,T0230,T0235,T0239,T0245,T0250,T0255,T0300,T0305,T0310,T0315,T0320,T0325,T0330,T0335,T0340,T0345,T0350,T0355,T0400,T0404,T0410,T0415,T0420,T0425,T0430,T0435,T0440,T0445,T0450,T0455,T0500,T0505,T0510,T0515,T0520,T0525,T0530,T0535,T0540,T0545,T0550,T0555,T0600,T0605,T0610,T0615,T0620,T0625,T0630,T0635,T0640,T0645,T0650,T0655,T0700,T0705,T0710,T0715,T0720,T0725,T0730,T0735,T0740,T0745,T0750,T0755,T0800,T0805,T0810,T0815,T0820,T0825,T0830,T0835,T0840,T0845,T0850,T0855,T0900,T0905,T0910,T0915,T0920,T0925,T0930,T0935,T0940,T0945,T0950,T0955,T1000,T1005,T1010,T1015,T1020,T1025,T1030,T1035,T1040,T1045,T1050,T1055,T1100,T1105,T1110,T1115,T1120,T1125,T1130,T1135,T1140,T1145,T1150,T1155,T1200,T1205,T1210,T1215,T1220,T1225,T1230,T1235,T1240,T1245,T1250,T1255,T1300,T1305,T1310,T1315,T1320,T1325,T1330,T1335,T1340,T1345,T1350,T1355,T1400,T1405,T1410,T1415,T1420,T1425,T1430,T1435,T1440,T1445,T1450,T1455,T1500,T1505,T1510,T1515,T1520,T1525,T1530,T1535,T1540,T1545,T1550,T1555,T1600,T1605,T1610,T1615,T1620,T1625,T1630,T1635,T1640,T1645,T1650,T1655,T1700,T1705,T1710,T1715,T1720,T1725,T1730,T1735,T1740,T1745,T1750,T1755,T1800,T1805,T1810,T1815,T1820,T1825,T1830,T1835,T1840,T1845,T1850,T1855,T1900,T1905,T1910,T1915,T2020,T1925,T1930,T1935,T1940,T1945,T1950,T1955,T2000,T2005,T2010,T2015,T2020,T2025,T2030,T2035,T2040,T2045,T2050,T2055,T2100,T2105,T2110,T2115,T2120,T2125,T2130,T2135,T2140,T2145,T2150,T2155,T2200,T2205,T2210,T2215,T2220,T2225,T2230,T2235,T2240,T2245,T2250,T2255,T2300,T2305,T2310,T2315,T2320,T2325,T2330,T2335,T2340,T2345,T2350,T2355,T2400" << endl;
	myfile << ", 1, , ,0900_0930, , , ,0.000571,0.000571,0.000571,0.000571,0.000571,0.000571,0.000506,0.000506,0.000506,0.000445,0.000445,0.000445,0.000391,0.000391,0.000391,0.000357,0.000357,0.000357,0.000328,0.000328,0.000328,0.000319,0.000319,0.000319,0.000302,0.000302,0.000302,0.000292,0.000292,0.000292,0.000296,0.000296,0.000296,0.00031,0.00031,0.00031,0.00031,0.00031,0.00031,0.000319,0.000319,0.000319,0.000383,0.000383,0.000383,0.000496,0.000496,0.000496,0.000568,0.000568,0.000568,0.000656,0.000656,0.000656,0.00095,0.00095,0.00095,0.001368,0.001368,0.001368,0.001587,0.001587,0.001587,0.00175,0.00175,0.00175,0.002288,0.002288,0.002288,0.002921,0.002921,0.002921,0.003242,0.003242,0.003242,0.003218,0.003218,0.003218,0.003803,0.003803,0.003803,0.004459,0.004459,0.004459,0.005002,0.005002,0.005002,0.005207,0.005207,0.005207,0.005677,0.005677,0.005677,0.005994,0.005994,0.005994,0.006018,0.006018,0.006018,0.005508,0.005508,0.005508,0.00529,0.00529,0.00529,0.005058,0.005058,0.005058,0.004833,0.004833,0.004833,0.004421,0.004421,0.004421,0.004327,0.004327,0.004327,0.004364,0.004364,0.004364,0.004343,0.004343,0.004343,0.004139,0.004139,0.004139,0.004201,0.004201,0.004201,0.004291,0.004291,0.004291,0.00435,0.00435,0.00435,0.004409,0.004409,0.004409,0.004566,0.004566,0.004566,0.004674,0.004674,0.004674,0.004761,0.004761,0.004761,0.004827,0.004827,0.004827,0.004882,0.004882,0.004882,0.0049,0.0049,0.0049,0.004887,0.004887,0.004887,0.004835,0.004835,0.004835,0.004899,0.004899,0.004899,0.005023,0.005023,0.005023,0.005065,0.005065,0.005065,0.005162,0.005162,0.005162,0.005436,0.005436,0.005436,0.005772,0.005772,0.005772,0.005907,0.005907,0.005907,0.005877,0.005877,0.005877,0.00605,0.00605,0.00605,0.006196,0.006196,0.006196,0.006248,0.006248,0.006248,0.006308,0.006308,0.006308,0.006404,0.006404,0.006404,0.006391,0.006391,0.006391,0.006401,0.006401,0.006401,0.006526,0.006526,0.006526,0.006574,0.006574,0.006574,0.006271,0.006271,0.006271,0.005937,0.005937,0.005937,0.005578,0.005578,0.005578,0.005293,0.005293,0.005293,0.004834,0.004834,0.004834,0.004387,0.004387,0.004387,0.00403,0.00403,0.00403,0.003748,0.003748,0.003748,0.003382,0.003382,0.003382,0.003121,0.003121,0.003121,0.002963,0.002963,0.002963,0.00289,0.00289,0.00289,0.002671,0.002671,0.002671,0.002468,0.002468,0.002468,0.002365,0.002365,0.002365,0.002249,0.002249,0.002249,0.002015,0.002015,0.002015,0.001784,0.001784,0.001784,0.00164,0.00164,0.00164,0.001474,0.001474,0.001474,0.001312,0.001312,0.001312,0.001132,0.001132,0.001132,0.001005,0.001005,0.001005,0.000889,0.000889,0.000889,0.000778,0.000778,0.000778,0.000676" << endl;
	myfile << "[demand_file_list],file_sequence_no,file_name,,format_type,demand_period,agent_type,scale_factor,departure_time_profile_no," << endl;
	myfile << ",0,demand.csv,,column,AM,auto,0.1,1," << endl;
	myfile.close();
	return;
}

bool CheckMeasurementFileExist()
{
	CDTACSVParser parser_measurement;
	if (parser_measurement.OpenCSVFile("measurement.csv", false))
	{

		int count = 0;
		while (parser_measurement.ReadRecord())  // if this line contains [] mark, then we will also read field headers.
		{
			string measurement_type;
			parser_measurement.GetValueByFieldName("measurement_type", measurement_type);

			if (measurement_type == "link")
			{
				int from_node_id;
				if (!parser_measurement.GetValueByFieldName("from_node_id", from_node_id))
					continue;

				int to_node_id;
				if (!parser_measurement.GetValueByFieldName("to_node_id", to_node_id))
					continue;

				count++;

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


bool CheckSupplySideScenarioFileExist()
{

	CDTACSVParser parser;

	int sa_count = 0;

	if (parser.OpenCSVFile("supply_side_scenario.csv", false))
	{
		while (parser.ReadRecord())
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
		parser.CloseCSVFile();
	}

	if (sa_count == 0)
		return false;
	else
		return true;
}


int main()
{
	// reset all the log files to defult 0: not output; if want to output these logs set to 1
	dtalog.output() << "DTALite Log" << std::fixed << std::setw(12) << '\n';
	dtalog.debug_level() = 0;
	dtalog.log_sig() = 0;
	dtalog.log_odme() = 0;
	dtalog.log_path() = 0;
	dtalog.log_dta() = 0;
	dtalog.log_ue() = 0;

	int column_generation_iterations = 0;
	int column_updating_iterations = 2;
	int ODME_iterations = 0;
	int sensitivity_analysis_iterations = 0;
	int number_of_memory_blocks = 4;
	float info_updating_freq_in_min = 5;
	int simulation_output = 0;

	int signal_updating_output = 0;
	// generate link performance and agent file
	int assignment_mode = 1;
	bool flag_default = false;
	int default_volume = 1;
	int link_length_in_meter_flag = 0;

	write_default_setting_file_if_not_exist();
	CDTACSVParser parser_settings;
	parser_settings.IsFirstLineHeader = false;


	if (parser_settings.OpenCSVFile("settings.csv", false))
	{
		while (parser_settings.ReadRecord_Section())
		{
			if (parser_settings.SectionName.compare("assignment") >= 0);
			{
				std::string assignment_mode_str;

				parser_settings.GetValueByFieldName("number_of_iterations", column_generation_iterations,false);
				// these are the assignment modes
				// two usually methods are ue (user equilibrium) and dta (dynamic traffic assignment)
				// the main difference of these two methods are different output in link_performance.csv
				// for basic uses set assignment mode to 'ue'
				// for more detailed link performances (one minute) set 'dta'1
				assignment_mode = 1;

				if (CheckMeasurementFileExist() == true)
				{
					column_updating_iterations = 5;
					ODME_iterations = 100;
				}
				else
					ODME_iterations = 0;


				if (CheckSupplySideScenarioFileExist() == true)
					sensitivity_analysis_iterations = 20;
				else
					sensitivity_analysis_iterations = -1;

				int route_output = 1;
				parser_settings.GetValueByFieldName("route_output", route_output, false, false);
				if (route_output == 0)   // reset back to lue mode
					assignment_mode = 0;

				simulation_output = 0;
				parser_settings.GetValueByFieldName("simulation_output", simulation_output, false, false);
				if (simulation_output ==1)
					assignment_mode = 2;

				// the start interation of generating signals, if there is no signals set this number larger than the iteration number
				number_of_memory_blocks = 1;
				parser_settings.GetValueByFieldName("number_of_memory_blocks", number_of_memory_blocks, false, false);
				double number_of_seconds_per_interval_input = 0.25;
				//parser_settings.GetValueByFieldName("number_of_seconds_per_interval", number_of_seconds_per_interval_input, false, false);

				//if (number_of_seconds_per_interval_input > 0.0001)
				//	number_of_seconds_per_interval = number_of_seconds_per_interval_input;

				dtalog.output() << "number_of_memory_blocks = " << number_of_memory_blocks << " in settings.csv." << std::endl;

				// just one record
				break;
			}

			if (parser_settings.SectionName == "[log]")
			{
				parser_settings.GetValueByFieldName("sig", dtalog.log_sig(), false);
				parser_settings.GetValueByFieldName("odme", dtalog.log_odme(), false);
				parser_settings.GetValueByFieldName("path", dtalog.log_path(), false);
				parser_settings.GetValueByFieldName("ue", dtalog.log_ue(), false);

				// just one record
				break;
			}
		}
	}
	// obtain initial flow values
	network_assignment(assignment_mode, column_generation_iterations, column_updating_iterations, ODME_iterations, sensitivity_analysis_iterations, simulation_output, number_of_memory_blocks);

	return 0;
}