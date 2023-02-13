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

bool g_read_subarea_CSV_file(Assignment& assignment)
{
	CDTACSVParser parser;
	if (parser.OpenCSVFile("subarea.csv",false) )
	{
		while (parser.ReadRecord())
		{
			string geo_string;
			if (parser.GetValueByFieldName("geometry", geo_string))
			{
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
	dtalog.output() << endl;
	dtalog.output() << "Step 1.8: Reading file section [departure_time_profile] in setting.csv..." << endl;
	parser.IsFirstLineHeader = false;


	if (parser.OpenCSVFile("settings.csv", false))
	{
		while (parser.ReadRecord_Section())
		{
			int departure_time_profile_no = 0;
			string time_period;

			CDeparture_time_Profile dep_time;

			if (parser.SectionName == "[departure_time_profile]")
			{

				if (!parser.GetValueByFieldName("departure_time_profile_no", departure_time_profile_no))
					break;

				dep_time.departure_time_profile_no = departure_time_profile_no;
				if (departure_time_profile_no != assignment.g_DepartureTimeProfileVector.size())
				{
					dtalog.output() << "Error: Field departure_time_profile_no in field departure_time_profile should be sequential as a value of ." << assignment.g_DepartureTimeProfileVector.size() << endl;
					g_program_stop();

				}

				vector<float> global_minute_vector;

				if (!parser.GetValueByFieldName("time_period", time_period))
				{
					dtalog.output() << "Error: Field time_period in field departure_time_profile cannot be read." << endl;
					g_program_stop();
				}


				//input_string includes the start and end time of a time period with hhmm format
				global_minute_vector = g_time_parser(time_period); //global_minute_vector incldue the starting and ending time

				if (global_minute_vector.size() == 2)
				{
					dep_time.starting_time_slot_no = global_minute_vector[0] / MIN_PER_TIMESLOT;  // read the data
					dep_time.ending_time_slot_no = global_minute_vector[1] / MIN_PER_TIMESLOT;    // read the data from setting.csv

				}

				char time_interval_field_name[20];
				char time_interval_field_name2[20];


				for (int s = dep_time.starting_time_slot_no; s <= dep_time.ending_time_slot_no; s += 1)
				{
					int hour = s / 12;
					int minute = (int)( (s / 12.0 - hour) * 60 + 0.5);

					double value = 0;
					sprintf(time_interval_field_name, "T%02d%02d", hour, minute);

					if (parser.GetValueByFieldName(time_interval_field_name, value, false) == true)
					{
						dep_time.departure_time_ratio[s] = value;
					}
					dtalog.output() << "T" << hour << "h" << minute << "min" << "=" << value << endl;

				}

				dep_time.compute_cumulative_profile(dep_time.starting_time_slot_no, dep_time.ending_time_slot_no);
				assignment.g_DepartureTimeProfileVector.push_back(dep_time);

				dtalog.output() << "compute_cumulative_profile!" << endl;

			}



		}

	}
}

double g_pre_read_demand_file(Assignment& assignment)
{

	float total_demand_in_demand_file = 0;

	CDTACSVParser parser;
	dtalog.output() << endl;
	dtalog.output() << "Step 1.8: Reading file section [demand_file_list] in setting.csv..." << endl;
	parser.IsFirstLineHeader = false;

	assignment.summary_file << "pre-read demand, defined in [demand_file_list] in settings.csv." << endl;

	if (parser.OpenCSVFile("settings.csv", false))
	{
		while (parser.ReadRecord_Section())
		{
			int this_departure_time_profile_no = 0;

			if (parser.SectionName == "[demand_file_list]")
			{
				int file_sequence_no = 1;

				string format_type = "null";

				int demand_format_flag = 0;

				if (!parser.GetValueByFieldName("file_sequence_no", file_sequence_no))
					break;

				// skip negative sequence no
				if (file_sequence_no <= -1)
					continue;

				double loading_scale_factor = 1.0;
				string file_name, demand_period_str, agent_type;
				parser.GetValueByFieldName("file_name", file_name);
				parser.GetValueByFieldName("demand_period", demand_period_str);
				parser.GetValueByFieldName("format_type", format_type);


				if (format_type.find("null") != string::npos)  // skip negative sequence no
				{
					dtalog.output() << "Please provide format_type in section [demand_file_list.]" << endl;
					g_program_stop();
				}


				if (format_type.find("column") != string::npos)  // or muliti-column
				{
					bool bFileReady = false;
					int error_count = 0;
					int critical_OD_count = 0;
					double critical_OD_volume = 0;

					// read the file formaly after the test.
					CDTACSVParser parser;
					int line_no = 0;
					if (parser.OpenCSVFile(file_name, false))
					{
						// read agent file line by line,

						int o_zone_id, d_zone_id;
						string agent_type, demand_period;

						std::vector <int> node_sequence;

						while (parser.ReadRecord())
						{
							float demand_value = 0;
							parser.GetValueByFieldName("o_zone_id", o_zone_id);
							parser.GetValueByFieldName("d_zone_id", d_zone_id);
							parser.GetValueByFieldName("volume", demand_value);

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
						}  // scan lines
					}

					// pre-read complete.
					break;
				}
				else if (format_type.compare("matrix") == 0)
				{
					bool bFileReady = false;
					int error_count = 0;
					int critical_OD_count = 0;
					double critical_OD_volume = 0;

					vector<int> LineIntegerVector;

					CDTACSVParser parser;
					parser.IsFirstLineHeader = false;
					if (parser.OpenCSVFile(file_name, true))
					{
						int i = 0;
						if (parser.ReadRecord())
						{
							parser.ConvertLineStringValueToIntegers();
							LineIntegerVector = parser.LineIntegerVector;
						}
						parser.CloseCSVFile();
					}

					int number_of_zones = LineIntegerVector.size();


					bFileReady = false;

					FILE* st;
					fopen_ss(&st, file_name.c_str(), "r");
					if (st != NULL)
					{
						// read the first line
						g_read_a_line(st);

						cout << "number of zones to be read = " << number_of_zones << endl;

						//test if a zone has been defined. 
						for (int destination_zone_index = 0; destination_zone_index < number_of_zones; destination_zone_index++)
						{
							int zone = LineIntegerVector[destination_zone_index];
							if (assignment.g_zoneid_to_zone_seq_no_mapping.find(zone) == assignment.g_zoneid_to_zone_seq_no_mapping.end())
							{
								if (error_count < 10)
									dtalog.output() << endl << "Warning: destination zone " << zone << "  has not been defined in node.csv" << endl;

								error_count++;
								// destination zone has not been defined, skipped.
								continue;
							}

						}


						int line_no = 0;
						for (int origin_zone_index = 0; origin_zone_index < number_of_zones; origin_zone_index++)
						{
							int origin_zone = (int)(g_read_float(st)); // read the origin zone number

							if (assignment.g_zoneid_to_zone_seq_no_mapping.find(origin_zone) == assignment.g_zoneid_to_zone_seq_no_mapping.end())
							{
								if (error_count < 10)
									dtalog.output() << endl << "Warning: destination zone " << origin_zone << "  has not been defined in node.csv" << endl;

								for (int destination_zone_index = 0; destination_zone_index < number_of_zones; destination_zone_index++)
								{
									float demand_value = g_read_float(st);
								}

								error_count++;
								// destination zone has not been defined, skipped.
								continue;
							}

							cout << "Reading file no." << file_sequence_no << " " << file_name << " at zone " << origin_zone << " ... " << endl;

							for (int destination_zone_index = 0; destination_zone_index < number_of_zones; destination_zone_index++)
							{
								int destination_zone = LineIntegerVector[destination_zone_index];

								float demand_value = g_read_float(st);

								int from_zone_seq_no = 0;
								int to_zone_seq_no = 0;
								from_zone_seq_no = assignment.g_zoneid_to_zone_seq_no_mapping[origin_zone];
								to_zone_seq_no = assignment.g_zoneid_to_zone_seq_no_mapping[destination_zone];

								g_zone_vector[from_zone_seq_no].preread_total_O_demand += demand_value;
								g_zone_vector[from_zone_seq_no].preread_ODdemand[to_zone_seq_no] += demand_value;

								total_demand_in_demand_file += demand_value;
								line_no++;
							}  // scan lines

						}

						fclose(st);
					}

					//preread complete
					break;
				}
			}
		}
	}
	return total_demand_in_demand_file;
}

void g_ReadDemandFileBasedOnDemandFileList(Assignment& assignment)
{
	g_load_demand_side_scenario_file(assignment);

	g_related_zone_vector_size = g_zone_vector.size();
	for (int orig = 0; orig < g_zone_vector.size(); orig++)  // o
	{
		g_zone_vector[orig].sindex = orig;  //copy this index 
	}

	// subarea handling step 1: reading
	g_read_subarea_CSV_file(assignment);
	// subarea handling step 2:
	// for each OD 

	if (assignment.g_subarea_shape_points.size() >= 3) // if there is a subarea defined.
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
		double total_demand = 0;
		double total_related_demand = 0;
		double total_related_demand_from_internal = 0;
		double total_related_demand_from_external = 0;
		double total_related_demand_from_external_cutoff = 0;
		double volume_cut_off_value = 0;

		int non_impact_zone_count = 0;
		int inside_zone_count = 0;
		int related_external_zone_count = 0;
		int related_external_zone_count_with_significant_volume = 0;
		// pre read demand 
		double total_preload_demand = g_pre_read_demand_file(assignment);

		cout << "total pre-load demand = " << total_preload_demand << endl;
		if (total_preload_demand <= 1)
		{
			g_program_stop();
		}
		else
		{
			// first stage: scan all origin zones
			// output: g_zone_vector[orig].subarea_significance_flag: true for inside zone, and related zone
			// total_related_demand: inside to all, outside passing to all
		for (int orig = 0; orig < g_zone_vector.size(); orig++)  // o
		{
			DTAGDPoint Pt;
			Pt.x = g_zone_vector[orig].cell_x;
			Pt.y = g_zone_vector[orig].cell_y;

			total_demand += g_zone_vector[orig].preread_total_O_demand;

			if (g_test_point_in_polygon(Pt, assignment.g_subarea_shape_points) == 1)
			{
				total_related_demand += g_zone_vector[orig].preread_total_O_demand;

				total_related_demand_from_internal += g_zone_vector[orig].preread_total_O_demand;

				g_zone_vector[orig].origin_zone_impact_volume = g_zone_vector[orig].preread_total_O_demand;
				g_zone_vector[orig].subarea_inside_flag = 3;  // zone is inside the subarea
				inside_zone_count++;
				continue; //inside the subarea, keep the zone anyway
			}



			double origin_zone_impact_volume = 0;

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

			double impact_demand_ratio = origin_zone_impact_volume / max(0.0001, g_zone_vector[orig].preread_total_O_demand);
			//			if (origin_zone_impact_volume < 5 && impact_demand_ratio <0.1)

			g_zone_vector[orig].origin_zone_impact_volume = origin_zone_impact_volume;

			if (origin_zone_impact_volume < 5)
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

		}

		assignment.summary_file << ",first stage scanning" << endl;
		assignment.summary_file << ",# of zones = " << inside_zone_count << ", total = " << total_related_demand << endl;
		assignment.summary_file << ",# inside_zone_count = " << inside_zone_count << ", inside demand = " << total_related_demand_from_internal << endl;
		assignment.summary_file << ",# related_external_zone_count = " << related_external_zone_count << ", external related demand = " << total_related_demand_from_external << endl;
		assignment.summary_file << ",# related_external_zone_count_with_significant_volume = " << related_external_zone_count_with_significant_volume << endl;


		// use the following rule to determine how many zones to keep as origin zones, we still keep all destination zones
		int cutoff_zone_size = min ( g_zone_vector.size(), max( related_external_zone_count_with_significant_volume*0.25, 100));


		if (cutoff_zone_size + inside_zone_count < g_zone_vector.size())  // additional handling  remaining zones are still greater than 1000
		// if (cutoff_zone_size > inside_zone_count)  // additional handling  remaining zones are still greater than 1000
		{

			std::vector <float> origin_zone_volume_vector;

			for (int orig = 0; orig < g_zone_vector.size(); orig++)  // o
			{
				if (g_zone_vector[orig].subarea_significance_flag == true)
				{
					if(g_zone_vector[orig].subarea_inside_flag == 2)
						origin_zone_volume_vector.push_back(g_zone_vector[orig].origin_zone_impact_volume);
					else
					{
						int idebug = 1;
					}
				}

			}

			// sort 

			if(origin_zone_volume_vector.size() > 0)
			{
			std::sort(origin_zone_volume_vector.begin(), origin_zone_volume_vector.end());

			// determine the cut off value;
			int cut_off_zone_index = max(0, origin_zone_volume_vector.size() - cutoff_zone_size);

			volume_cut_off_value = origin_zone_volume_vector[cut_off_zone_index];



			for (int orig = 0; orig < g_zone_vector.size(); orig++)  // o
			{
				if (g_zone_vector[orig].subarea_inside_flag == 2 && g_zone_vector[orig].subarea_significance_flag == true)
				{
					if (g_zone_vector[orig].origin_zone_impact_volume < volume_cut_off_value)  // volume_cut_off_value is based on the value just established
					{
						g_zone_vector[orig].subarea_significance_flag = false; // set non significant
						g_zone_vector[orig].subarea_inside_flag = 1;  // zone is set to boundary
//						related_zone_count++;
					}
					else
					{
						g_zone_vector[orig].subarea_inside_flag = 2;  // zone is set to significantly related
						total_related_demand_from_external_cutoff += g_zone_vector[orig].preread_total_O_related_demand;
						
					}

				}

			}

			}

			//re select the impact origin zones 
		}
		
		assignment.summary_file << ",second stage cut off" << endl;
		assignment.summary_file << ",# cut off zone size = " << cutoff_zone_size << ", total_related_demand_from_external_cutoff  = " << total_related_demand_from_external_cutoff << endl;
		assignment.summary_file << ",cut off threshold = " << volume_cut_off_value << endl;
		double cut_off_demand_ratio = total_related_demand_from_external_cutoff / max(1, total_related_demand_from_external);
		assignment.summary_file << ",cut off demand ratio = " << cut_off_demand_ratio << endl;

		
		}

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
		g_related_zone_vector_size = sindex;

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
							g_zone_vector[orig].distrct_cluster_index = distrct_index;

							distrct_index++;
						}
						else
						{
							g_zone_vector[orig].distrct_cluster_index = -1;
						}

					}

				}
				assignment.g_number_of_analysis_districts = distrct_index;

				//clustering 
				for (int orig = 0; orig < g_zone_vector.size(); orig++)  // o
				{
					if (g_zone_vector[orig].distrct_cluster_index == -1)
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
									g_zone_vector[orig].distrct_cluster_index = g_zone_vector[d].distrct_cluster_index;
									assignment.g_zone_seq_no_to_analysis_distrct_id_mapping[orig] = g_zone_vector[d].distrct_cluster_index;

								}

							}
						}

					}

				}



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
						if (g_zone_vector[orig].distrct_cluster_index == district)
						{
							DTAGDPoint pt;
							pt.x = g_zone_vector[orig].cell_x;
							pt.y = g_zone_vector[orig].cell_y;

							points.push_back(pt);
						}
					}


					if(points.size()>=4)  // 4 as restrictive conditions
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
						g_district_info_base0_map[district].shape_points.push_back(points_in_polygon[i]);

					}


				}

			}

			// end of generating analaysis district


			if (g_related_zone_vector_size >= 30000)  // do not use this function for now.
			{
				// k mean method
				//random select 300 zones then find the closest zone to aggregate

				double super_zone_ratio = 300.0/g_related_zone_vector_size;

				int super_zone_index = assignment.g_number_of_analysis_districts+1;  // automated district  id start from the externally given district id 

				for (int orig = 0; orig < g_zone_vector.size(); orig++)  // o
				{
					if (g_zone_vector[orig].sindex>=0)  // no significant: skip
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
						if(g_zone_vector[orig].sindex >=0) // related zones
						{
							g_zone_vector[orig].sindex = g_zone_vector[orig].superzone_index;  //rewrite using super zones
						}

				}


			}
		}



		FILE* g_pFileZone = nullptr;
		g_pFileZone = fopen("subarea_related_zone.csv", "w");

		if (g_pFileZone == NULL)
		{
			cout << "File subarea_related_zone.csv cannot be opened." << endl;
			g_program_stop();
		}

		fprintf(g_pFileZone, "first_column,zone_id,x_coord,y_coord,super_zone_id,analysis_district_id,demand,subarea_significance_flag,inside_flag\n");
		for (int orig = 0; orig < g_zone_vector.size(); orig++)  // o
		{
				fprintf(g_pFileZone, ",%d,%f,%f,%d,%d,%f,%d,%d\n", 
					g_zone_vector[orig].zone_id, g_zone_vector[orig].cell_x, g_zone_vector[orig].cell_y, 
					g_zone_vector[orig].sindex,
					g_zone_vector[orig].distrct_cluster_index,
					g_zone_vector[orig].origin_zone_impact_volume,
					g_zone_vector[orig].subarea_significance_flag, 
					g_zone_vector[orig].subarea_inside_flag);
		}

		fclose(g_pFileZone);

		dtalog.output() << endl;
		double non_related_zone_ratio = non_impact_zone_count * 1.0 / g_zone_vector.size();
		dtalog.output() << non_impact_zone_count << " zones are not signifantly related by the subarea. = (" << non_related_zone_ratio << ")" <<endl;
		dtalog.output() << inside_zone_count << " zones are inside the subarea. " << endl;



		double related_ratio = total_related_demand / max(1, total_demand);
		dtalog.output() << "total demand = " << total_demand << " total subarea-related demand = " << total_related_demand << " ("
			<< related_ratio << " )" << endl;
		
		//assignment.summary_file << ",# of subarea related zones =," << g_related_zone_vector_size << endl;
		//assignment.summary_file << "," <<" # of inside zones = " << inside_zone_count <<  endl;
		//assignment.summary_file << "," << " # of related outside zones = " << related_external_zone_count << endl;
	}

	//	fprintf(g_pFileOutputLog, "number of zones =,%lu\n", g_zone_vector.size());
	g_read_departure_time_profile(assignment);

	assignment.InitializeDemandMatrix(g_related_zone_vector_size, g_zone_vector.size(), assignment.g_AgentTypeVector.size(), assignment.g_DemandPeriodVector.size());
	
	float related_flow_matrix[4][4] = { 0 };

	float total_demand_in_demand_file = 0;

	CDTACSVParser parser;
	dtalog.output() << endl;
	dtalog.output() << "Step 1.8: Reading file section [demand_file_list] in setting.csv..." << endl;
	parser.IsFirstLineHeader = false;

	assignment.summary_file << "Step 2.1: read demand, defined in [demand_file_list] in settings.csv." << endl;

	if (parser.OpenCSVFile("settings.csv", false))
	{
		while (parser.ReadRecord_Section())
		{
			int this_departure_time_profile_no = 0;

			if (parser.SectionName == "[demand_file_list]")
			{
				int file_sequence_no = 1;

				string format_type = "null";

				int demand_format_flag = 0;

				if (!parser.GetValueByFieldName("file_sequence_no", file_sequence_no))
					break;

				// skip negative sequence no
				if (file_sequence_no <= -1)
					continue;

				double loading_scale_factor = 1.0;
				string file_name, demand_period_str, agent_type;
				parser.GetValueByFieldName("file_name", file_name);
				parser.GetValueByFieldName("demand_period", demand_period_str);
				parser.GetValueByFieldName("format_type", format_type);
				parser.GetValueByFieldName("scale_factor", loading_scale_factor,false);
				parser.GetValueByFieldName("departure_time_profile_no", this_departure_time_profile_no, false);

				if (this_departure_time_profile_no >= assignment.g_DepartureTimeProfileVector.size())
				{
					dtalog.output() << "Error: departure_time_profile_no = " << this_departure_time_profile_no << " in  section [demand_file_list]  has not been defined in section [departure_time_profile]." << endl;
					g_program_stop();

				}

				parser.GetValueByFieldName("agent_type", agent_type);

				int agent_type_no = 0;
				int demand_period_no = 0;

				if (assignment.demand_period_to_seqno_mapping.find(demand_period_str) != assignment.demand_period_to_seqno_mapping.end())
					demand_period_no = assignment.demand_period_to_seqno_mapping[demand_period_str];
				else
				{
					dtalog.output() << "Error: demand period = " << demand_period_str.c_str() << " in  section [demand_file_list]  has not been defined." << endl;
					g_program_stop();
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
					dtalog.output() << "Please provide format_type in section [demand_file_list.]" << endl;
					g_program_stop();
				}



				bool b_multi_agent_list = false;

				if (agent_type == "multi_agent_list")
					b_multi_agent_list = true;
				else
				{
					if (assignment.agent_type_2_seqno_mapping.find(agent_type) != assignment.agent_type_2_seqno_mapping.end())
						agent_type_no = assignment.agent_type_2_seqno_mapping[agent_type];
					else
					{
						dtalog.output() << "Error: agent_type = " << agent_type.c_str() << " in field agent_type of section [demand_file_list] in file setting.csv cannot be found." << endl;
						g_program_stop();
					}
				}

				if (demand_period_no > MAX_TIMEPERIODS)
				{
					dtalog.output() << "demand_period_no should be less than settings in demand_period section. Please change the parameter settings in the source code." << endl;
					g_program_stop();
				}

				if (format_type.find("column") != string::npos)  // or muliti-column
				{
					bool bFileReady = false;
					int error_count = 0;
					int critical_OD_count = 0;
					double critical_OD_volume = 0;


					// read the file formaly after the test.
					CDTACSVParser parser;
					int line_no = 0;
					if (parser.OpenCSVFile(file_name, false))
					{
						// read agent file line by line,

						int o_zone_id, d_zone_id;
						string agent_type, demand_period;

						std::vector <int> node_sequence;

						assignment.summary_file << ",o_zone_id,d_zone_id,inside_flag_from,inside_flag_to,volume" << endl;

						float total_matrix_demand_value = 0;

						while (parser.ReadRecord())
						{
							float demand_value = 0;							
							parser.GetValueByFieldName("o_zone_id", o_zone_id);
							parser.GetValueByFieldName("d_zone_id", d_zone_id);
							parser.GetValueByFieldName("volume", demand_value);

							if (o_zone_id == d_zone_id)
							{
								continue;
							}

							if (assignment.g_zoneid_to_zone_seq_no_mapping.find(o_zone_id) == assignment.g_zoneid_to_zone_seq_no_mapping.end())
							{
								if (error_count < 10)
									dtalog.output() << endl << "Warning: origin zone " << o_zone_id << "  has not been defined in node.csv" << endl;

								error_count++;
								// origin zone has not been defined, skipped.
								continue;
							}

							if (assignment.g_zoneid_to_zone_seq_no_mapping.find(d_zone_id) == assignment.g_zoneid_to_zone_seq_no_mapping.end())
							{
								if (error_count < 10)
									dtalog.output() << endl << "Warning: destination zone " << d_zone_id << "  has not been defined in node.csv" << endl;

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
								if(inside_flag_to == 2)
								{ 
								int idebug = 1;
								}

							}


							if (inside_flag_from == 0 || inside_flag_to == 0)  // one of trip ends are not related the subarea
							{
								if (g_zone_vector[from_zone_seq_no].impact_passing_ODdemand.find(to_zone_seq_no)== g_zone_vector[from_zone_seq_no].impact_passing_ODdemand.end())
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
							//	assignment.summary_file << "," << o_zone_id << "," << d_zone_id << "," << inside_flag_from << "," << inside_flag_to << "," << demand_value << ",cd= " << total_matrix_demand_value << endl;
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

							demand_value *= (loading_scale_factor* local_scale_factor);


							if (demand_value >= 5)
							{
								critical_OD_volume += demand_value;
								critical_OD_count += 1;
								//dtalog.output() << origin_zone << "," << destination_zone << "," << demand_value << "," << "\"LINESTRING( " <<
								//    assignment.zone_id_X_mapping[origin_zone] << " " << assignment.zone_id_Y_mapping[origin_zone] << "," <<
								//    assignment.zone_id_X_mapping[destination_zone] << " " << assignment.zone_id_Y_mapping[destination_zone] << ")\" " << endl;

							}

							int analysis_district_id = assignment.g_zone_seq_no_to_analysis_distrct_id_mapping[from_zone_seq_no];
							g_district_info_base0_map[analysis_district_id].record_origin_2_district_volume(agent_type_no, demand_value);

							assignment.total_demand[agent_type_no][demand_period_no] += demand_value;
							assignment.g_column_pool[from_zone_sindex][to_zone_sindex][agent_type_no][demand_period_no].od_volume += demand_value;
							assignment.g_column_pool[from_zone_sindex][to_zone_sindex][agent_type_no][demand_period_no].departure_time_profile_no = this_departure_time_profile_no;
							assignment.total_demand_volume += demand_value;
							assignment.g_origin_demand_array[from_zone_seq_no] += demand_value;
							
						
							//assignment.summary_file << ",,,,,,,,," << o_zone_id << "," << d_zone_id << "," << inside_flag_from << "," << inside_flag_to << "," << demand_value << ",cd= " << assignment.total_demand_volume << endl;
						
							//// print out all OD pairs in the matrix

							if (total_matrix_demand_value > assignment.total_demand_volume)
							{
								g_program_stop();
							}

							// we generate vehicles here for each OD data line
							if (line_no <= 5)  // read only one line, but has not reached the end of the line
								dtalog.output() << "o_zone_id:" << o_zone_id << ", d_zone_id: " << d_zone_id << ", value = " << demand_value << endl;

							line_no++;
						}  // scan lines
					}
						else
					{
						// open file
						dtalog.output() << "Error: File " << file_name << " cannot be opened.\n It might be currently used and locked by EXCEL." << endl;
						g_program_stop();
					}



					dtalog.output() << "total demand volume is " << assignment.total_demand_volume << endl;
					dtalog.output() << "crtical demand volume has " << critical_OD_count << " OD pairs in size," << critical_OD_volume << ", " << ", account for " << critical_OD_volume / max(0.1f, assignment.total_demand_volume) * 100 << "%%" << endl;

					dtalog.output() << "crtical OD zones volume has " << critical_OD_count << " OD pairs in size," << critical_OD_volume << ", " << ", account for " << critical_OD_volume / max(0.1f, assignment.total_demand_volume) * 100 << "%%" << endl;



					std::map<int, float>::iterator it;

				}
				else if (format_type.compare("path") == 0)
				{

					int path_counts = 0;
					float sum_of_path_volume = 0;
					CDTACSVParser parser;
					if (parser.OpenCSVFile(file_name, false))
					{
						int total_path_in_demand_file = 0;
						// read agent file line by line,

						int o_zone_id, d_zone_id;
						string agent_type, demand_period;

						std::vector <int> node_sequence;

						while (parser.ReadRecord())
						{
							total_path_in_demand_file++;
							if (total_path_in_demand_file % 1000 == 0)
								dtalog.output() << "total_path_in_demand_file is " << total_path_in_demand_file << endl;

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

								assignment.total_demand[agent_type_no][demand_period_no] += agent_path_element.volume;
								assignment.g_column_pool[from_zone_sindex][to_zone_sindex][agent_type_no][demand_period_no].od_volume += agent_path_element.volume;
								assignment.g_column_pool[from_zone_sindex][to_zone_sindex][agent_type_no][demand_period_no].departure_time_profile_no = this_departure_time_profile_no;
								assignment.total_demand_volume += agent_path_element.volume;
								assignment.g_origin_demand_array[from_zone_seq_no] += agent_path_element.volume;

								int analysis_district_id = assignment.g_zone_seq_no_to_analysis_distrct_id_mapping[from_zone_seq_no];
								g_district_info_base0_map[analysis_district_id].record_origin_2_district_volume(agent_type_no, agent_path_element.volume);

							//apply for both agent csv and routing policy
							assignment.g_column_pool[from_zone_seq_no][to_zone_seq_no][agent_type_no][demand_period_no].bfixed_route = true;

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

								CColumnVector* pColumnVector = &(assignment.g_column_pool[from_zone_sindex][to_zone_sindex][agent_type_no][demand_period_no]);
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
							}
						}
						dtalog.output() << "total_demand_volume loaded from path file is " << sum_of_path_volume << " with " << path_counts << "paths." << endl;

					}
					else
					{
						//open file
						dtalog.output() << "Error: File " << file_name << " cannot be opened.\n It might be currently used and locked by EXCEL." << endl;
						g_program_stop();
					}
				}
				else if (format_type.compare("activity_chain") == 0)
				{
					/////////////////////
					int path_counts = 0;
					float sum_of_path_volume = 0;
					CDTACSVParser parser;
					if (parser.OpenCSVFile(file_name, false))
					{
						int total_path_in_demand_file = 0;
						// read agent file line by line,

						int agent_id, o_zone_id, d_zone_id;
						string agent_type, demand_period;

						std::vector <int> node_sequence;

						while (parser.ReadRecord())
						{
							total_path_in_demand_file++;
							if (total_path_in_demand_file % 1000 == 0)
								dtalog.output() << "total_path_in_demand_file is " << total_path_in_demand_file << endl;

							parser.GetValueByFieldName("agent_id", agent_id);
							parser.GetValueByFieldName("o_zone_id", o_zone_id);
							parser.GetValueByFieldName("d_zone_id", d_zone_id);

							CAgentPath agent_path_element;


							int from_zone_seq_no = 0;
							int to_zone_seq_no = 0;

							// add protection later for activity and path data types
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

							assignment.total_demand[agent_type_no][demand_period_no] += agent_path_element.volume;
							assignment.g_column_pool[from_zone_sindex][to_zone_sindex][agent_type_no][demand_period_no].od_volume += agent_path_element.volume;
							assignment.g_column_pool[from_zone_sindex][to_zone_sindex][agent_type_no][demand_period_no].departure_time_profile_no = this_departure_time_profile_no;
							assignment.total_demand_volume += agent_path_element.volume;
							assignment.g_origin_demand_array[from_zone_seq_no] += agent_path_element.volume;

							int analysis_district_id = assignment.g_zone_seq_no_to_analysis_distrct_id_mapping[from_zone_seq_no];
							g_district_info_base0_map[analysis_district_id].record_origin_2_district_volume(agent_type_no, agent_path_element.volume);
							//apply for both agent csv and routing policy

							bool bValid = true;

							string activity_zone_sequence;
							parser.GetValueByFieldName("activity_zone_sequence", activity_zone_sequence);

							string activity_agent_type_sequence;
							parser.GetValueByFieldName("activity_agent_type_sequence", activity_agent_type_sequence);

							std::vector<int> zone_id_sequence;
							std::vector<string> agent_type_string_sequence;

							int activity_zone_size = g_ParserIntSequence(activity_zone_sequence, zone_id_sequence);
							int agent_type_size = g_ParserStringSequence(activity_agent_type_sequence, agent_type_string_sequence);

							if (activity_zone_size > 0 && agent_type_size != activity_zone_size + 1)
							{
							    dtalog.output() << "Error: agent_type_size != activity_zone_size + 1 in activity plan csv" << endl;
								dtalog.output() << "agent_id = " << agent_id << endl;
								dtalog.output() << "zone_sequence = " << activity_zone_sequence.c_str() << endl;
								dtalog.output() << "agent_type_sequence = " << activity_agent_type_sequence.c_str() << endl;
								g_program_stop();
							}
							std::vector<int> agent_type_no_sequence;

							for (int i = 0; i < agent_type_string_sequence.size(); ++i)
							{
								if (assignment.agent_type_2_seqno_mapping.find(agent_type_string_sequence[i]) != assignment.agent_type_2_seqno_mapping.end())
								{
									int activty_type_no = assignment.agent_type_2_seqno_mapping[agent_type_string_sequence[i]];  // normal agent type
									agent_type_no_sequence.push_back(activty_type_no);
								}
								else
								{
									dtalog.output() << "Error: agent_type in activity_agent_type_sequence cannot be found in activity plan csv" << endl;
									dtalog.output() << "agent_type = " << agent_type_string_sequence[i].c_str() << endl;
									continue;
								}

							}
							assignment.agent_type_2_seqno_mapping[activity_agent_type_sequence];

							std::vector<int> zone_no_sequence;

							for (int i = 0; i < zone_id_sequence.size(); ++i)
							{
								zone_no_sequence.push_back(from_zone_seq_no);  // first segment

								if (assignment.g_zoneid_to_zone_seq_no_mapping.find(zone_id_sequence[i]) != assignment.g_zoneid_to_zone_seq_no_mapping.end())
								{
									int internal_zone_seq_no = assignment.g_zoneid_to_zone_seq_no_mapping[zone_id_sequence[i]];  // map external node number to internal node seq no.
									assignment.zone_seq_no_2_activity_mapping[internal_zone_seq_no] = 1; // notify the column pool to prepare the paths for this activity zones, even no direct OD volume 
									zone_no_sequence.push_back(internal_zone_seq_no);
								}
								else
								{
									dtalog.output() << "Error: zone in activity_zone_sequence cannot be found in activity plan csv" << endl;
									dtalog.output() << "zone = " << zone_id_sequence[i]<< endl;
									continue;
								}

							}

							int prev_zone_seq_no = from_zone_seq_no;


							int prev_activty_type_no = -1;
							
							if(zone_id_sequence.size() >=1)
							{
								prev_activty_type_no = agent_type_no_sequence[0];
								zone_no_sequence.push_back(to_zone_seq_no);  // last segment;
							}
							

							for (int i = 0; i < zone_no_sequence.size()-1; ++i)
							{  // loop through activity plan with multiple intermediate destinations
								int current_activty_type_no = agent_type_no_sequence[i];  // normal agent type which is not activity plan type
								assignment.g_column_pool[zone_no_sequence[i]][zone_no_sequence[i+1]][current_activty_type_no][demand_period_no].od_volume += 0.01;  // small volume to ensure column pool exist 
								
							}

							if (bValid)
							{
								// first origin to final destination OD pair 
								CColumnVector* pColumnVector = &(assignment.g_column_pool[from_zone_seq_no][to_zone_seq_no][agent_type_no][demand_period_no]);

								pColumnVector->od_volume = volume;
								pColumnVector->departure_time_profile_no = this_departure_time_profile_no;
									// to be protected later by checking

								if (zone_no_sequence.size() >= 1)
								{

								pColumnVector->activity_zone_no_vector = zone_no_sequence;
								pColumnVector->activity_zone_sequence = activity_zone_sequence;
								pColumnVector->activity_agent_type_sequence = activity_agent_type_sequence;
								pColumnVector->activity_agent_type_no_vector = agent_type_no_sequence;


								}
								// in total, we have # of activity zones + 2 as zone no sequence vector

							}

						}//end of file
					} //open file
					else
					{
						//open file
						dtalog.output() << "Error: File " << file_name << " cannot be opened.\n It might be currently used and locked by EXCEL." << endl;
						g_program_stop();
					}

				}// activity
				else if (format_type.compare("matrix") == 0)
				{
					bool bFileReady = false;
					int error_count = 0;
					int critical_OD_count = 0;
					double critical_OD_volume = 0;

					vector<int> LineIntegerVector;

					CDTACSVParser parser;
					parser.IsFirstLineHeader = false;
					if (parser.OpenCSVFile(file_name, true))
					{
						int i = 0;
						if (parser.ReadRecord())
						{
							parser.ConvertLineStringValueToIntegers();
							LineIntegerVector = parser.LineIntegerVector;
						}
						parser.CloseCSVFile();
					}

					int number_of_zones = LineIntegerVector.size();


					bFileReady = false;

					FILE* st;
					fopen_ss(&st, file_name.c_str(), "r");
					if (st != NULL)
					{
						// read the first line
						g_read_a_line(st);

						cout << "number of zones to be read = " << number_of_zones << endl;

						//test if a zone has been defined. 
						for (int destination_zone_index = 0; destination_zone_index < number_of_zones; destination_zone_index++)
						{
							int zone = LineIntegerVector[destination_zone_index];
							if (assignment.g_zoneid_to_zone_seq_no_mapping.find(zone) == assignment.g_zoneid_to_zone_seq_no_mapping.end())
							{
								if (error_count < 10)
									dtalog.output() << endl << "Warning: destination zone " << zone << "  has not been defined in node.csv" << endl;

								error_count++;
								// destination zone has not been defined, skipped.
								continue;
							}

						}


						int line_no = 0;
						for (int origin_zone_index = 0; origin_zone_index < number_of_zones; origin_zone_index++)
						{
							int origin_zone = (int)(g_read_float(st)); // read the origin zone number

							if (assignment.g_zoneid_to_zone_seq_no_mapping.find(origin_zone) == assignment.g_zoneid_to_zone_seq_no_mapping.end())
							{
								if (error_count < 10)
									dtalog.output() << endl << "Warning: destination zone " << origin_zone << "  has not been defined in node.csv" << endl;

								for (int destination_zone_index = 0; destination_zone_index < number_of_zones; destination_zone_index++)
								{
									float demand_value = g_read_float(st);
								}

								error_count++;
								// destination zone has not been defined, skipped.
								continue;
							}

							cout << "Reading file no." << file_sequence_no << " " << file_name << " at zone " << origin_zone << " ... " << endl;

							for (int destination_zone_index = 0; destination_zone_index < number_of_zones; destination_zone_index++)
							{
								int destination_zone = LineIntegerVector[destination_zone_index];

								float demand_value = g_read_float(st);

								int from_zone_seq_no = 0;
								int to_zone_seq_no = 0;
								from_zone_seq_no = assignment.g_zoneid_to_zone_seq_no_mapping[origin_zone];
								to_zone_seq_no = assignment.g_zoneid_to_zone_seq_no_mapping[destination_zone];

								int from_zone_sindex = g_zone_vector[from_zone_seq_no].sindex;
								if (from_zone_sindex == -1)
									continue;

								int to_zone_sindex = g_zone_vector[to_zone_seq_no].sindex;
								if (to_zone_sindex == -1)
									continue;

								// encounter return
								if (demand_value < -99)
									break;


								demand_value *= loading_scale_factor;
								if (demand_value >= 1)
								{
									if (assignment.shortest_path_log_zone_id == -1)  // set this to the first zone
										assignment.shortest_path_log_zone_id = origin_zone;

									critical_OD_volume += demand_value;
									critical_OD_count += 1;
									//dtalog.output() << origin_zone << "," << destination_zone << "," << demand_value << "," << "\"LINESTRING( " <<
									//    assignment.zone_id_X_mapping[origin_zone] << " " << assignment.zone_id_Y_mapping[origin_zone] << "," <<
									//    assignment.zone_id_X_mapping[destination_zone] << " " << assignment.zone_id_Y_mapping[destination_zone] << ")\" " << endl;

								}

								assignment.total_demand[agent_type_no][demand_period_no] += demand_value;
								assignment.g_column_pool[from_zone_sindex][to_zone_sindex][agent_type_no][demand_period_no].od_volume += demand_value;
								assignment.g_column_pool[from_zone_sindex][to_zone_sindex][agent_type_no][demand_period_no].departure_time_profile_no = this_departure_time_profile_no;
								assignment.total_demand_volume += demand_value;
								assignment.g_origin_demand_array[from_zone_seq_no] += demand_value;

								int analysis_district_id = assignment.g_zone_seq_no_to_analysis_distrct_id_mapping[from_zone_seq_no];
								g_district_info_base0_map[analysis_district_id].record_origin_2_district_volume(agent_type_no, demand_value);

								// we generate vehicles here for each OD data line
								if (line_no <= 5)  // read only one line, but has not reached the end of the line
									dtalog.output() << "o_zone_id:" << origin_zone << ", d_zone_id: " << destination_zone << ", value = " << demand_value << endl;

								line_no++;
							}  // scan lines

						}

						fclose(st);

						dtalog.output() << "total demand volume is " << assignment.total_demand_volume << endl;
						dtalog.output() << "crtical demand volume has " << critical_OD_count << " OD pairs in size," << critical_OD_volume << ", " << ", account for " << critical_OD_volume / max(0.1f, assignment.total_demand_volume) * 100 << "%%" << endl;

						dtalog.output() << "crtical OD zones volume has " << critical_OD_count << " OD pairs in size," << critical_OD_volume << ", " << ", account for " << critical_OD_volume / max(0.1f, assignment.total_demand_volume) * 100 << "%%" << endl;

						//assignment.summary_file << "crtical demand volume has " << critical_OD_count << " OD pairs in size," << critical_OD_volume << ", " << ", account for " << critical_OD_volume / max(0.1, assignment.total_demand_volume) * 100 << "%%" << endl;
						//assignment.summary_file << "crtical OD zones volume has " << critical_OD_count << " OD pairs in size," << critical_OD_volume << ", " << ", account for " << critical_OD_volume / max(0.1, assignment.total_demand_volume) * 100 << "%%" << endl;

						std::map<int, float>::iterator it;
						int count_zone_demand = 0;
						for (it = assignment.g_origin_demand_array.begin(); it != assignment.g_origin_demand_array.end(); ++it)
						{
							if (it->second > 0.001)
							{
								dtalog.output() << "o_zone " << it->first << ", demand=," << it->second << endl;
								count_zone_demand++;
							}
						}
						dtalog.output() << "There are  " << count_zone_demand << " zones with positive demand" << endl;


					} //end reading file
					else
					{
						// open file
						dtalog.output() << "Error: File " << file_name << " cannot be opened.\n It might be currently used and locked by EXCEL." << endl;
						g_program_stop();
					}

				}
				else {
					dtalog.output() << "Error: format_type = " << format_type << " is not supported. Currently DTALite supports format such as column, matrix, activity_plan, path." << endl;
					g_program_stop();
				}

				
				assignment.summary_file << ",file_sequence_no=," << file_sequence_no << ",file_name =, "<< file_name.c_str() << ",demand_period =, "
					<< demand_period_str.c_str() << ",departure_time_profile_no=,"<< this_departure_time_profile_no << ",cumulative demand =, " << assignment.total_demand_volume << endl;

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
	assignment.summary_file << ",total demand =, " << assignment.total_demand_volume << endl;
	
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

			for (int at = 0; at < assignment.g_AgentTypeVector.size(); at++)  //m
			{
				for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); tau++)  //tau
				{
					p_column_pool = &(assignment.g_column_pool[from_zone_sindex][to_zone_sindex][at][tau]);
					if (p_column_pool->od_volume > 0)
					{
						CODState ods;
						ods.setup_input(orig, dest, at, tau, inside_flag_from, inside_flag_to);

						float passing_demand = p_column_pool->od_volume;
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

	assignment.summary_file << ",from_flag,to_flag,volume" << endl;
	assignment.summary_file << ",from_flag,to_flag,volume" << endl;
	assignment.summary_file << ",0=not related, 1=signficiantly related external but not in cutoff, 2= cut-off external, 3= inside" << endl;
	assignment.summary_file << ",";
	for (int subare_inside_flag_to = 0; subare_inside_flag_to <= 3; subare_inside_flag_to++)
		assignment.summary_file << subare_inside_flag_to << ",";

	assignment.summary_file << endl;

	for(int subare_inside_flag_from = 0; subare_inside_flag_from <=3; subare_inside_flag_from++)
	{
		assignment.summary_file <<  subare_inside_flag_from ;
		for (int subare_inside_flag_to = 0; subare_inside_flag_to <=3; subare_inside_flag_to++)
		{ 
			assignment.summary_file << "," << related_flow_matrix[subare_inside_flag_from][subare_inside_flag_to];
		}
		assignment.summary_file << endl;
	}


	std::sort(ODStateVector.begin(), ODStateVector.end());

	assignment.summary_file << ",top 10 OD,rank,o,d,inside_flag_o,inside_flag_d,agent_type,departure_time,volume" << endl;



	for (int k = 0; k < min(size_t(100), ODStateVector.size()); k++)
	{
		int o = ODStateVector[k].orig;
		int d = ODStateVector[k].dest;
		int at = ODStateVector[k].at;
		int tau = ODStateVector[k].tau;

		assignment.summary_file << ",," << k+1 << "," << g_zone_vector[o].zone_id << "," << g_zone_vector[d].zone_id << "," <<
			ODStateVector[k].subarea_inside_flag_orig << ","  << ODStateVector[k].subarea_inside_flag_dest << "," <<
			assignment.g_AgentTypeVector[at].agent_type.c_str() << "," << assignment.g_DemandPeriodVector[tau].demand_period.c_str()
			<< "," << ODStateVector[k].value << endl;
	}
}

void g_ReadOutputFileConfiguration(Assignment& assignment)
{
	dtalog.output() << "Step 1.9: Reading file section [output_file_configuration] in setting.csv..." << endl;

	cout << "Step 1.8: Reading file section [output_file_configuration] in setting.csv..." << endl;

	CDTACSVParser parser;
	parser.IsFirstLineHeader = false;
	if (parser.OpenCSVFile("settings.csv", false))
	{
		while (parser.ReadRecord_Section())
		{
			if (parser.SectionName == "[output_file_configuration]")
			{
				parser.GetValueByFieldName("path_output", assignment.path_output, false, false);
				parser.GetValueByFieldName("major_path_volume_threshold", assignment.major_path_volume_threshold, false, false);
				parser.GetValueByFieldName("shortest_path_log_zone_id", assignment.shortest_path_log_zone_id, false, false);
				parser.GetValueByFieldName("trajectory_output_count", assignment.trajectory_output_count, false, false);
				parser.GetValueByFieldName("trace_output", assignment.trace_output, false, false);

				parser.GetValueByFieldName("trajectory_sampling_rate", assignment.trajectory_sampling_rate, false, false);
				parser.GetValueByFieldName("trajectory_diversion_only", assignment.trajectory_diversion_only, false, false);
				parser.GetValueByFieldName("td_link_performance_sampling_interval_in_min", assignment.td_link_performance_sampling_interval_in_min, false, false);


				dtalog.output() << "dynamic_link_performance_sampling_interval_in_min= " << assignment.td_link_performance_sampling_interval_in_min << " min" << endl;
				dtalog.output() << "dynamic_link_performance_sampling_interval_hd_in_min= " << assignment.dynamic_link_performance_sampling_interval_hd_in_min << " min" << endl;

			}
		}

		parser.CloseCSVFile();
	}
}

void g_ReadInformationConfiguration(Assignment& assignment)
{
	dtalog.output() << "Step 1.91: Reading file section [real_time_info] in setting.csv..." << endl;

	cout << "Step 1.91: Reading file section [real_time_info] in setting.csv..." << endl;

	CDTACSVParser parser;
	parser.IsFirstLineHeader = false;
	if (parser.OpenCSVFile("settings.csv", false))
	{

		while (parser.ReadRecord_Section())
		{
			if (parser.SectionName == "[real_time_info]")
			{

				parser.GetValueByFieldName("info_updating_freq_in_min", assignment.g_info_updating_freq_in_min, false, false);
				dtalog.output() << "info_updating_freq_in_min= " << assignment.g_info_updating_freq_in_min << " min" << endl;

				parser.GetValueByFieldName("visual_distance_in_cells", assignment.g_visual_distance_in_cells, false, false);
				dtalog.output() << "visual_distance_in_cells= " << assignment.g_visual_distance_in_cells << " cells" << endl;

			}
		}

		parser.CloseCSVFile();
	}
	assignment.summary_file << "[real_time_info] " << endl;

	assignment.summary_file << ",info_updating_freq_in_min= " << assignment.g_info_updating_freq_in_min << " min" << endl;
	
	assignment.summary_file << ",visual_distance_in_cells= " << assignment.g_visual_distance_in_cells << " cells" << endl;
	

}
void g_add_new_virtual_connector_link(int internal_from_node_seq_no, int internal_to_node_seq_no, string agent_type_str, int zone_seq_no = -1)
{
	// create a link object
	CLink link;
	link.link_id = "connector";
	link.from_node_seq_no = internal_from_node_seq_no;
	link.to_node_seq_no = internal_to_node_seq_no;
	link.link_seq_no = assignment.g_number_of_links;
	link.to_node_seq_no = internal_to_node_seq_no;
	//virtual connector
	link.link_type = -1;
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
		link.VDF_period[tau].FFTT = 0.0001;
		link.VDF_period[tau].alpha = 0;
		link.VDF_period[tau].beta = 0;
		link.VDF_period[tau].allowed_uses = agent_type_str;

		link.travel_time_per_period[tau] = 0;

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
	dtalog.output() << endl;
	dtalog.output() << "Step 1.8: Reading file section [demand_file_list] in setting.csv..." << endl;
	parser.IsFirstLineHeader = false;

	if (parser.OpenCSVFile("settings.csv", false))
	{
		while (parser.ReadRecord_Section())
		{

			if (parser.SectionName == "[demand_file_list]")
			{
				int file_sequence_no = 1;

				string format_type = "null";

				int demand_format_flag = 0;

				if (!parser.GetValueByFieldName("file_sequence_no", file_sequence_no))
					break;

				// skip negative sequence no
				if (file_sequence_no <= -1)
					continue;

				string file_name, demand_period_str, agent_type;
				parser.GetValueByFieldName("file_name", file_name);
				parser.GetValueByFieldName("format_type", format_type);

				if (format_type.find("column") != string::npos)  // or muliti-column
				{
					// read the file formaly after the test.
					CDTACSVParser parser;

					if (parser.OpenCSVFile(file_name, false))
					{
						return 0;
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
				if (format_type.find("path") != string::npos)
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
		}
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
			assignment.summary_file << ", number of activity nodes defined in node.csv=, " << number_of_activity_nodes << endl;
			assignment.summary_file << ", number of boundary activity nodes defined in node.csv=, " << number_of_boundary_nodes << endl;

			

			return number_of_activity_nodes;
		}

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
					char VDF_field_name[50];
					bool VDF_required_field_flag = true;
//					sprintf(VDF_field_name, "QVDF_plf%d", demand_period_id);
//					parser.GetValueByFieldName(VDF_field_name, this_link.VDF_period[tau].queue_demand_factor, VDF_required_field_flag, false);
					sprintf(VDF_field_name, "QVDF_alpha%d", demand_period_id);
					parser.GetValueByFieldName(VDF_field_name, this_link.VDF_period[tau].Q_alpha, VDF_required_field_flag, false);
					sprintf(VDF_field_name, "QVDF_beta%d", demand_period_id);
					parser.GetValueByFieldName(VDF_field_name, this_link.VDF_period[tau].Q_beta, VDF_required_field_flag, false);
					sprintf(VDF_field_name, "QVDF_cd%d", demand_period_id);
					parser.GetValueByFieldName(VDF_field_name, this_link.VDF_period[tau].Q_cd, VDF_required_field_flag, false);
					sprintf(VDF_field_name, "QVDF_cp%d", demand_period_id);
					parser.GetValueByFieldName(VDF_field_name, this_link.VDF_period[tau].Q_cp, VDF_required_field_flag, false);
					sprintf(VDF_field_name, "QVDF_n%d", demand_period_id);
					parser.GetValueByFieldName(VDF_field_name, this_link.VDF_period[tau].Q_n, VDF_required_field_flag, false);
					sprintf(VDF_field_name, "QVDF_s%d", demand_period_id);
					parser.GetValueByFieldName(VDF_field_name, this_link.VDF_period[tau].Q_s, VDF_required_field_flag, false);
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
					dtalog.output() << "Error: from_node_id " << from_node_id << " in file measurement.csv is not defined in node.csv." << endl;
					//has not been defined
					continue;
				}
				if (assignment.g_node_id_to_seq_no_map.find(to_node_id) == assignment.g_node_id_to_seq_no_map.end())
				{
					dtalog.output() << "Error: to_node_id " << to_node_id << " in file measurement.csv is not defined in node.csv." << endl;
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
							char VDF_field_name[50];
							bool VDF_required_field_flag = true;
							sprintf(VDF_field_name, "QVDF_qdf%d", demand_period_id);
							parser.GetValueByFieldName(VDF_field_name, p_link->VDF_period[tau].queue_demand_factor, VDF_required_field_flag, false);
							sprintf(VDF_field_name, "QVDF_alpha%d", demand_period_id);
							parser.GetValueByFieldName(VDF_field_name, p_link->VDF_period[tau].Q_alpha, VDF_required_field_flag, false);
							sprintf(VDF_field_name, "QVDF_beta%d", demand_period_id);
							parser.GetValueByFieldName(VDF_field_name, p_link->VDF_period[tau].Q_beta, VDF_required_field_flag, false);
							sprintf(VDF_field_name, "QVDF_cd%d", demand_period_id);
							parser.GetValueByFieldName(VDF_field_name, p_link->VDF_period[tau].Q_cd, VDF_required_field_flag, false);
							sprintf(VDF_field_name, "QVDF_n%d", demand_period_id);
							parser.GetValueByFieldName(VDF_field_name, p_link->VDF_period[tau].Q_n, VDF_required_field_flag, false);
							sprintf(VDF_field_name, "QVDF_cp%d", demand_period_id);
							parser.GetValueByFieldName(VDF_field_name, p_link->VDF_period[tau].Q_cp, VDF_required_field_flag, false);
							sprintf(VDF_field_name, "QVDF_s%d", demand_period_id);
							parser.GetValueByFieldName(VDF_field_name, p_link->VDF_period[tau].Q_s, VDF_required_field_flag, false);

						}
					}
				}
			}
		}
		parser.CloseCSVFile();
	}


	for (int i = 0; i < g_link_vector.size(); i++)
	{
		for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
		{
			CLink* p_link = &(g_link_vector[i]);

			if (p_link->VDF_period[tau].queue_demand_factor < 0.0) // no data
			{
				string vdf_code = "all";
				if (p_link->vdf_code.size() > 0)
					vdf_code = p_link->vdf_code;

				if (g_vdf_type_map.find(vdf_code) == g_vdf_type_map.end())
				{
					vdf_code = "all";  // default if vdf_code  has no data
				}
				// apply default data 
				if (g_vdf_type_map.find(vdf_code) != g_vdf_type_map.end())
				{

					p_link->VDF_period[tau].queue_demand_factor = g_vdf_type_map[vdf_code].VDF_period_sum[tau].queue_demand_factor;
					p_link->VDF_period[tau].Q_alpha = g_vdf_type_map[vdf_code].VDF_period_sum[tau].Q_alpha;
					p_link->VDF_period[tau].Q_beta = g_vdf_type_map[vdf_code].VDF_period_sum[tau].Q_beta;
					p_link->VDF_period[tau].Q_cd = g_vdf_type_map[vdf_code].VDF_period_sum[tau].Q_cd;
					p_link->VDF_period[tau].Q_n = g_vdf_type_map[vdf_code].VDF_period_sum[tau].Q_n;
					p_link->VDF_period[tau].Q_s = g_vdf_type_map[vdf_code].VDF_period_sum[tau].Q_s;
					p_link->VDF_period[tau].Q_cp = g_vdf_type_map[vdf_code].VDF_period_sum[tau].Q_cp;
				}
				else
				{
					p_link->VDF_period[tau].queue_demand_factor = 1.0 / max(0.01, p_link->VDF_period[tau].L);  // final default, even not link_qvd.csv is provided

				}



			} // we should have QVDF data for sure


		}

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
		dtalog.output() << "File final_summary.csv cannot be opened." << endl;
		g_program_stop();
	}
	else
	{
		fclose(g_pFilePathMOE);
	}

	fopen_ss(&g_pFilePathMOE, "link_performance.csv", "w");
	if (!g_pFilePathMOE)
	{
		dtalog.output() << "File link_performance.csv cannot be opened." << endl;
		g_program_stop();
	}
	else
	{
		fclose(g_pFilePathMOE);
	}

	fopen_ss(&g_pFilePathMOE, "od_performance.csv", "w");
	if (!g_pFilePathMOE)
	{
		dtalog.output() << "File od_performance.csv cannot be opened." << endl;
		g_program_stop();
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
		dtalog.output() << "File subarea_link_performance.csv cannot be opened." << endl;
		g_program_stop();
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
			dtalog.output() << "File agent.csv cannot be opened." << endl;
			g_program_stop();
		}
		else
		{
			fclose(g_pFilePathMOE);
		}

		fopen_ss(&g_pFilePathMOE, "td_link_performance.csv", "w");
		if (!g_pFilePathMOE)
		{
			dtalog.output() << "File td_link_performance.csv cannot be opened." << endl;
			g_program_stop();
		}
		else
		{
			fclose(g_pFilePathMOE);
		}

		fopen_ss(&g_pFilePathMOE, "trajectory.csv", "w");
		if (!g_pFilePathMOE)
		{
			dtalog.output() << "File trajectory.csv cannot be opened." << endl;
			g_program_stop();
		}
		else
		{
			fclose(g_pFilePathMOE);
		}

	}


	fopen_ss(&g_pFilePathMOE, "route_assignment.csv", "w");
	if (!g_pFilePathMOE)
	{
		dtalog.output() << "File route_assignment.csv cannot be opened." << endl;
		g_program_stop();
	}
	else
	{
		fclose(g_pFilePathMOE);
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
	dtalog.output() << "_____________" << endl;
	dtalog.output() << "Step 1: Reading input data" << endl;
	dtalog.output() << "_____________" << endl;

	dtalog.output() << "Step 1.1: Reading section [demand_period] in setting.csv..." << endl;

	parser_demand_period.IsFirstLineHeader = false;
	if (parser_demand_period.OpenCSVFile("settings.csv", false))
	{
		while (parser_demand_period.ReadRecord_Section())
		{
			if (parser_demand_period.SectionName == "[demand_period]")
			{
				CDemand_Period demand_period;

				if (!parser_demand_period.GetValueByFieldName("demand_period_id", demand_period.demand_period_id))
					break;

				if (!parser_demand_period.GetValueByFieldName("demand_period", demand_period.demand_period))
				{
					dtalog.output() << "Error: Field demand_period in file demand_period cannot be read." << endl;
					g_program_stop();
				}

				vector<float> global_minute_vector;

				if (!parser_demand_period.GetValueByFieldName("time_period", demand_period.time_period))
				{
					dtalog.output() << "Error: Field time_period in file demand_period cannot be read." << endl;
					g_program_stop();
				}



				//input_string includes the start and end time of a time period with hhmm format
				global_minute_vector = g_time_parser(demand_period.time_period); //global_minute_vector incldue the starting and ending time

				if (global_minute_vector.size() == 2)
				{
					demand_period.starting_time_slot_no = global_minute_vector[0] / MIN_PER_TIMESLOT;  // read the data
					demand_period.ending_time_slot_no = global_minute_vector[1] / MIN_PER_TIMESLOT;    // read the data from setting.csv
					demand_period.time_period_in_hour = (global_minute_vector[1] - global_minute_vector[0]) / 60.0;
					demand_period.t2_peak_in_hour = (global_minute_vector[0] + global_minute_vector[1]) / 2 / 60;

					//g_fout << global_minute_vector[0] << endl;
					//g_fout << global_minute_vector[1] << endl;


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
					dep_time.departure_time_ratio[s] = 1.0/300.0;
				}

				dep_time.compute_cumulative_profile(demand_period.starting_time_slot_no, demand_period.ending_time_slot_no);

				if (assignment.g_DepartureTimeProfileVector.size() == 0)
				{
					//default profile 
					assignment.g_DepartureTimeProfileVector.push_back(dep_time);
				}

			}
		}

		parser_demand_period.CloseCSVFile();

		if (assignment.g_DemandPeriodVector.size() == 0)
		{
			dtalog.output() << "Error:  Section demand_period has no information." << endl;
			g_program_stop();
		}
	}
	else
	{
		dtalog.output() << "Error: File settings.csv cannot be opened.\n It might be currently used and locked by EXCEL." << endl;
		g_program_stop();
	}

	dtalog.output() << "number of demand periods = " << assignment.g_DemandPeriodVector.size() << endl;

	assignment.g_number_of_demand_periods = assignment.g_DemandPeriodVector.size();

	if (assignment.g_number_of_demand_periods >= MAX_TIMEPERIODS)
	{
		dtalog.output() << "Error: the number of demand periods in settings.csv os greater than the internal size of MAX_TIMEPERIODS.\nPlease contact developers" << endl;
		g_program_stop();
	}
	//step 1:read demand type file

	dtalog.output() << "Step 1.2: Reading section [link_type] in setting.csv..." << endl;

	CDTACSVParser parser_link_type;
	parser_link_type.IsFirstLineHeader = false;
	if (parser_link_type.OpenCSVFile("settings.csv", false))
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

		while (parser_link_type.ReadRecord_Section())
		{
			if (parser_link_type.SectionName == "[link_type]")
			{
				CLinkType element;

				if (!parser_link_type.GetValueByFieldName("link_type", element.link_type))
				{
					if (line_no == 0)
					{
						dtalog.output() << "Error: Field link_type cannot be found in file link_type section." << endl;
						g_program_stop();
					}
					else
					{
						// read empty line
						break;
					}
				}

				if (assignment.g_LinkTypeMap.find(element.link_type) != assignment.g_LinkTypeMap.end())
				{
					dtalog.output() << "Error: Field link_type " << element.link_type << " has been defined more than once in file link_type section" << endl;
					g_program_stop();
				}

				string traffic_flow_code_str;
				parser_link_type.GetValueByFieldName("type_code", element.type_code, true);

				string vdf_type_str;

				element.vdf_type = bpr_vdf;
				parser_link_type.GetValueByFieldName("vdf_type", vdf_type_str,false);


				if (vdf_type_str == "bpr")
					element.vdf_type = bpr_vdf;
				if (vdf_type_str == "qvdf")
					element.vdf_type = q_vdf;


				element.traffic_flow_code = spatial_queue;

				parser_link_type.GetValueByFieldName("traffic_flow_model", traffic_flow_code_str,false);
				parser_link_type.GetValueByFieldName("k_jam", element.k_jam, false);

				// by default bpr


				if (traffic_flow_code_str == "point_queue")
					element.traffic_flow_code = point_queue;

				if (traffic_flow_code_str == "spatial_queue")
					element.traffic_flow_code = spatial_queue;

				if (traffic_flow_code_str == "kw")
					element.traffic_flow_code = kinemative_wave;

				dtalog.output() << "important: traffic_flow_code on link type " << element.link_type << " is " << element.traffic_flow_code << endl;


				assignment.g_LinkTypeMap[element.link_type] = element;
				line_no++;
			}
		}

		parser_link_type.CloseCSVFile();
	}

	dtalog.output() << "number of link types = " << assignment.g_LinkTypeMap.size() << endl;

	CDTACSVParser parser_agent_type;
	dtalog.output() << "Step 1.3: Reading section [agent_type] in setting.csv..." << endl;

	parser_agent_type.IsFirstLineHeader = false;
	if (parser_agent_type.OpenCSVFile("settings.csv", false))
	{

		assignment.g_AgentTypeVector.clear();
		while (parser_agent_type.ReadRecord_Section())
		{

			if (parser_agent_type.SectionName == "[agent_type]")
			{
				CAgent_type agent_type;

				if (!parser_agent_type.GetValueByFieldName("agent_type", agent_type.agent_type))
					break;

				agent_type.agent_type_no = -1;
				parser_agent_type.GetValueByFieldName("agent_type_no", agent_type.agent_type_no, false);

				if (agent_type.agent_type_no == -1)
				{
					agent_type.agent_type_no = assignment.g_AgentTypeVector.size() + 1;
				}


				//substring overlapping checking 

				{
					for (int at = 0; at < assignment.g_AgentTypeVector.size(); at++)
					{
						if (assignment.g_AgentTypeVector[at].agent_type.find(agent_type.agent_type) != string::npos)
						{
							dtalog.output() << "Error substring duplication checking : agent_type = " << assignment.g_AgentTypeVector[at].agent_type.c_str() <<
								" in section agent_type is overlapping with " << agent_type.agent_type.c_str() << ". Please add flags such as _only to avoid overlapping in the use of allowe_uses field.";
							g_program_stop();

						}

					}

				}

				parser_agent_type.GetValueByFieldName("vot", agent_type.value_of_time, false, false);

				// scan through the map with different node sum for different paths
				parser_agent_type.GetValueByFieldName("pce", agent_type.PCE, false, false);
				parser_agent_type.GetValueByFieldName("person_occupancy", agent_type.OCC, false, false);
				parser_agent_type.GetValueByFieldName("desired_speed_ratio", agent_type.DSR, false, false);

				parser_agent_type.GetValueByFieldName("headway", agent_type.time_headway_in_sec, false, false);
				parser_agent_type.GetValueByFieldName("display_code", agent_type.display_code, false);
				parser_agent_type.GetValueByFieldName("real_time_info", agent_type.real_time_information, false);

				if (agent_type.agent_type == "rt")  // set the real time information type = 1 for rt class by default
					agent_type.real_time_information = 1;

				if (agent_type.agent_type == "dms")  // set the real time information type = 1 for dms class by default
					agent_type.real_time_information = 2;


				parser_agent_type.GetValueByFieldName("access_node_type", agent_type.access_node_type, false);

				if (agent_type.access_node_type.size() > 0)
				{
					parser_agent_type.GetValueByFieldName("access_speed", agent_type.access_speed);
					parser_agent_type.GetValueByFieldName("access_distance_lb", agent_type.access_distance_lb);
					parser_agent_type.GetValueByFieldName("access_distance_ub", agent_type.access_distance_ub);

					if (agent_type.access_distance_ub < 100)
					{
						dtalog.output() << "Error: access_distance_ub = " << agent_type.access_distance_ub << "< 100. Please ensure the unit is meter." << endl;
						g_program_stop();
					}
					parser_agent_type.GetValueByFieldName("acecss_link_k", agent_type.acecss_link_k);
				}

				assignment.agent_type_2_seqno_mapping[agent_type.agent_type] = assignment.g_AgentTypeVector.size();

				assignment.g_AgentTypeVector.push_back(agent_type);

				assignment.summary_file << "agent_type =, " << agent_type.agent_type.c_str() << ", real time info flag = " << agent_type.real_time_information << endl;

				assignment.g_number_of_agent_types = assignment.g_AgentTypeVector.size();
			}
		}
		parser_agent_type.CloseCSVFile();

		if (assignment.g_AgentTypeVector.size() == 0)
			dtalog.output() << "Error: Section agent_type does not contain information." << endl;
	}

	if (assignment.g_AgentTypeVector.size() >= MAX_AGNETTYPES)
	{
		dtalog.output() << "Error: agent_type = " << assignment.g_AgentTypeVector.size() << " in section agent_type is too large. " << endl <<"MAX_AGNETTYPES = " << MAX_AGNETTYPES << "Please contact program developers!";
		g_program_stop();
	}

	dtalog.output() << "number of agent typess = " << assignment.g_AgentTypeVector.size() << endl;


	assignment.g_number_of_nodes = 0;
	assignment.g_number_of_links = 0;  // initialize  the counter to 0


	int number_of_zones = g_detect_if_zones_defined_in_node_csv(assignment);
	// = 1: normal
	//= 0, there are is boundary 
	//=-1, no information



	if (number_of_zones <=1)
	{
		CDTACSVParser parser_z;
		if (parser_z.OpenCSVFile("zone.csv", true))
		{
			parser_z.CloseCSVFile();
		}else
		{   // without zone.csv file

				if(g_TAZ_2_GMNS_zone_generation(assignment)==false)
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

	dtalog.output() << "Step 1.3: Reading zone data in zone.csv..." << endl;


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
			parser.GetValueByFieldName("access_node_vector", access_node_vector_str);

			std::vector<int> access_node_vector;

			g_ParserIntSequence(access_node_vector_str, access_node_vector);

			for (int i = 0; i < access_node_vector.size(); i++)
			{
				assignment.access_node_id_to_zone_id_map[access_node_vector[i]] = zone_id;
			}

			float production = 0;
			float attraction = 0;
			parser.GetValueByFieldName("production", production, false);
			parser.GetValueByFieldName("attraction", attraction, false);


			zone_id_production[zone_id] = production;
			zone_id_attraction[zone_id] = attraction;
			// push it to the global node vector
			dtalog.output() << "reading " << assignment.access_node_id_to_zone_id_map.size() << " access nodes from zone.csv.. " << endl;
		}

		parser.CloseCSVFile();
	}



	dtalog.output() << "Step 1.4: Reading node data in node.csv..." << endl;
	std::map<int, int> zone_id_to_analysis_district_id_mapping;

	// MRM: stage 1: filtering, output: gate nodes, bridge nodes for micro and macro networks

	// MRM: stage 2: read node layer
	int max_layer_flag = 0;
	

//	assignment.summary_file << ", # of node layers to be read," << max_layer_flag+1 << "," << endl;

	for(int layer_no = 0; layer_no <= max_layer_flag; layer_no ++)
	{
		dtalog.output() << "Reading node data at layer" << layer_no << endl;

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


				if(max_layer_flag == 1) // MRM mode
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
			
	/*			if (analysis_district_id >= MAX_ORIGIN_DISTRICTS)
				{
					dtalog.output() << "Error: grid_id in node.csv is larger than predefined value of 20." << endl;
					g_program_stop();
				}*/

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

					string str_agent_type;
					parser.GetValueByFieldName("agent_type", str_agent_type, false); //step 2 for adding access links: read agent_type for adding access links

					if (str_agent_type.size() > 0 && assignment.agent_type_2_seqno_mapping.find(str_agent_type) != assignment.agent_type_2_seqno_mapping.end())
					{
						node.agent_type_str = str_agent_type;
						node.agent_type_no = assignment.agent_type_2_seqno_mapping[str_agent_type];
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

						assignment.zone_id_X_mapping[zone_id] = node.x;
						assignment.zone_id_Y_mapping[zone_id] = node.y;
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
					dtalog.output() << "reading " << assignment.g_number_of_nodes << " nodes.. " << endl;
			}



			dtalog.output() << "number of nodes = " << assignment.g_number_of_nodes << endl;
			dtalog.output() << "number of multimodal activity nodes = " << multmodal_activity_node_count << endl;
			dtalog.output() << "number of zones = " << zone_id_mapping.size() << endl;

			// fprintf(g_pFileOutputLog, "number of nodes =,%d\n", assignment.g_number_of_nodes);
			parser.CloseCSVFile();
		}
	}



	int debug_line_count = 0;


	/// <summary>  mappping node to zone
	// hanlding multimodal access link: stage 1
	//step 3 for adding access links: there is node type restriction defined in agent type section of settings.csv
	for (int at = 0; at < assignment.g_AgentTypeVector.size(); ++at) // first loop for each agent type
	{
		if (assignment.g_AgentTypeVector[at].access_node_type.size() > 0)  // for each multmodal agent type
		{
			// find the closest zone id

			if (debug_line_count <= 20)
			{

				dtalog.output() << " multimodal access link generation condition 1: agent type " << assignment.g_AgentTypeVector[at].agent_type.c_str() << " has access node type" << assignment.g_AgentTypeVector[at].access_node_type.size() << endl;
				// zone without multimodal access
				debug_line_count++;
			}

			for (int a_k = 0; a_k < g_node_vector.size(); a_k++)
			{
				if (g_node_vector[a_k].is_activity_node == 1 && g_node_vector[a_k].agent_type_no == at) //second loop for mode_specific activity node
				{

					int zone_id = g_node_vector[a_k].zone_org_id;
					int zone_seq_no = zone_seq_no = assignment.g_zoneid_to_zone_seq_no_mapping[zone_id];

					if (debug_line_count <= 20)
					{

						dtalog.output() << " multimodal access link generation condition 2: agent type no = " << at << " for node no. " << a_k << "as activity node with zone_id >=1" << endl;
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

							if (assignment.g_AgentTypeVector[at].access_node_type.find(g_node_vector[i].node_type) != string::npos)  // check allowed access code
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

								if (distance >= assignment.g_AgentTypeVector[at].access_distance_lb && distance <= assignment.g_AgentTypeVector[at].access_distance_ub)  // check the range 
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

						if (access_node_distance_vector.size() > assignment.g_AgentTypeVector[at].acecss_link_k)
						{

							std::vector<float> access_node_distance_vector_temp;
							access_node_distance_vector_temp = access_node_distance_vector;
							std::sort(access_node_distance_vector_temp.begin(), access_node_distance_vector_temp.end());

							distance_k_cut_off_value = access_node_distance_vector_temp[max(0, assignment.g_AgentTypeVector[at].acecss_link_k - 1)];
							//distance_k can be dynamically determined based on the density of stops and stations at different areas, e.g.CBM vs. rual area
						}

						for (int an = 0; an < access_node_seq_vector.size(); an++)
						{
							if (access_node_distance_vector[an] < distance_k_cut_off_value)  // within the shortest k ranage 
							{
								g_add_new_access_link(a_k, access_node_seq_vector[an], access_node_distance_vector[an], at, -1);
								//incoming connector from station to activity centers
								g_add_new_access_link(access_node_seq_vector[an], a_k, access_node_distance_vector[an], at, -1);
								assignment.g_AgentTypeVector[at].zone_id_cover_map[zone_id] = true;
							}
						}

					}
					else if (min_distance_node_seq_no >= 0 && min_distance < assignment.g_AgentTypeVector[at].access_distance_ub)  // no node in the  preferred range, just use any feasible node with minimum distance by default
					{
						g_add_new_access_link(a_k, min_distance_node_seq_no, min_distance, at, -1);
						g_add_new_access_link(min_distance_node_seq_no, a_k, min_distance, at, -1);
						assignment.g_AgentTypeVector[at].zone_id_cover_map[zone_id] = true;

					}
					else {

						//                        dtalog.output() << " zone" << g_node_vector[a_k].zone_org_id << " with agent type = " << assignment.g_AgentTypeVector[at].agent_type.c_str() << " has no access to stop or station" << endl;
												// zone without multimodal access
					}

				}  // for each zone

			}
		}
	}// for each agent type 

	//g_InfoZoneMapping(assignment);
	//g_OutputModelFiles(1);  // node


	// initialize zone vector
	dtalog.output() << "Step 1.5: Initializing O-D zone vector..." << endl;

	std::map<int, int>::iterator it;
	// creating zone centriod
	for (it = zone_id_mapping.begin(); it != zone_id_mapping.end(); ++it)
	{
		COZone ozone;

		// for each zone, we have to also create centriod
		ozone.zone_id = it->first;  // zone_id


		ozone.zone_seq_no = g_zone_vector.size();
		ozone.obs_production = zone_id_production[it->first];
		ozone.obs_attraction = zone_id_attraction[it->first];
		ozone.cell_x = assignment.zone_id_X_mapping[it->first];
		ozone.cell_y = assignment.zone_id_Y_mapping[it->first];
		ozone.gravity_production = zone_id_production[it->first];
		ozone.gravity_attraction = zone_id_attraction[it->first];

		assignment.g_zoneid_to_zone_seq_no_mapping[ozone.zone_id] = ozone.zone_seq_no;  // create the zone id to zone se&zq no mapping
		assignment.g_zoneid_to_zone_sindex_no_mapping[ozone.zone_id] = ozone.zone_seq_no;  // create the zone id to zone se&zq no mapping

		if (zone_id_to_analysis_district_id_mapping[ozone.zone_id]+1 > assignment.g_number_of_analysis_districts )
			assignment.g_number_of_analysis_districts = zone_id_to_analysis_district_id_mapping[ozone.zone_id]+1;

		assignment.g_zone_seq_no_to_analysis_distrct_id_mapping[ozone.zone_seq_no] = zone_id_to_analysis_district_id_mapping[ozone.zone_id];

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


	dtalog.output() << "number of zones = " << g_zone_vector.size() << endl;

	g_zone_to_access(assignment);  // only under zone2connector mode
//	g_OutputModelFiles(3);  // node


	int demand_data_mode = g_detect_if_demand_data_provided(assignment);


	if (assignment.assignment_mode== simulation_dta && demand_data_mode >= 1)
	{
		g_demand_file_generation(assignment);
	}


	// MRM: stage 3: read link layers
	// step 4: read link file

	CDTACSVParser parser_link;

	int link_type_warning_count = 0;
	bool length_in_km_waring = false;
	dtalog.output() << "Step 1.6: Reading link data in link.csv... " << endl;

	for (int layer_no = 0; layer_no <= max_layer_flag; layer_no++)
	{
		dtalog.output() << "link layer= " << layer_no << endl;


		string file_name = "link.csv";


		if (parser_link.OpenCSVFile(file_name.c_str(), true))
		{



			while (parser_link.ReadRecord())  // if this line contains [] mark, then we will also read field headers.
			{

				// create a link object  // this should be place inside
				CLink link;

				string link_type_name_str;
				parser_link.GetValueByFieldName("link_type_name", link_type_name_str, false);

				long from_node_id= -1;
				if (!parser_link.GetValueByFieldName("from_node_id", from_node_id))
					continue;

				long to_node_id =-1;
				if (!parser_link.GetValueByFieldName("to_node_id", to_node_id))
					continue;

				string linkID;
				parser_link.GetValueByFieldName("link_id", linkID,false);
				// add the to node id into the outbound (adjacent) node list

				if (from_node_id == 9439)
				{
					int trace_id = 1;
				}


				if (max_layer_flag == 1) // MRM mode
				{
					if (layer_no == 0)  // macro link layer
					{  // skip all links inside non-gate node

						if (assignment.g_node_id_to_MRM_subarea_mapping.find(from_node_id) != assignment.g_node_id_to_MRM_subarea_mapping.end()
							|| assignment.g_node_id_to_MRM_subarea_mapping.find(to_node_id) != assignment.g_node_id_to_MRM_subarea_mapping.end())
						{
							if (assignment.g_node_id_to_MRM_subarea_mapping[from_node_id] == 1
						&&		assignment.g_node_id_to_MRM_subarea_mapping[to_node_id] == 1)  // both upstream and downstream inside
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

					dtalog.output() << "Error: from_node_id " << from_node_id << " in file link.csv is not defined in node.csv." << endl;
					dtalog.output() << " from_node_id mapping =" << MRM_subarea_mapping_flag << endl;
					continue; //has not been defined
				}

				if (assignment.g_node_id_to_seq_no_map.find(to_node_id) == assignment.g_node_id_to_seq_no_map.end())
				{

					int MRM_subarea_mapping_flag = assignment.g_node_id_to_MRM_subarea_mapping[from_node_id];
					int MRM_subarea_mapping_flag_to = assignment.g_node_id_to_MRM_subarea_mapping[to_node_id];

					dtalog.output() << "Error: to_node_id " << to_node_id << " in file link.csv is not defined in node.csv." << endl;
					continue; //has not been defined
				}

				//if (assignment.g_link_id_map.find(linkID) != assignment.g_link_id_map.end())
				//    dtalog.output() << "Error: link_id " << linkID.c_str() << " has been defined more than once. Please check link.csv." << endl;

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

				int meso_link_id=-1;
				parser_link.GetValueByFieldName("meso_link_id", meso_link_id, false);

				if(meso_link_id>=0)
					link.meso_link_id = meso_link_id;

				parser_link.GetValueByFieldName("geometry", link.geometry, false);
				parser_link.GetValueByFieldName("link_code", link.link_code_str, false);

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

				// Peiheng, 05/13/21, if setting.csv does not have corresponding link type or the whole section is missing, set it as 2 (i.e., Major arterial)
				int link_type = 2;
				parser_link.GetValueByFieldName("link_type", link_type, false);

				if (assignment.g_LinkTypeMap.find(link_type) == assignment.g_LinkTypeMap.end())
				{
					if (link_type_warning_count < 10)
					{
						dtalog.output() << "link type " << link_type << " in link.csv is not defined for link " << from_node_id << "->" << to_node_id << " in link_type section in setting.csv" << endl;
						link_type_warning_count++;

						CLinkType element_vc;
						// -1 is for virutal connector
						element_vc.link_type = link_type;
						element_vc.link_type_name = link.link_type_name;
						assignment.g_LinkTypeMap[element_vc.link_type] = element_vc;

					}
					// link.link_type has been taken care by its default constructor
					//g_program_stop();
				}
				else
				{
					// link type should be defined in settings.csv
					link.link_type = link_type;
				}

				if (assignment.g_LinkTypeMap[link.link_type].type_code == "c")  // suggestion: we can move "c" as "connector" in allowed_uses
				{
					if (g_node_vector[internal_from_node_seq_no].zone_org_id >= 1)
					{
						int zone_org_id = g_node_vector[internal_from_node_seq_no].zone_org_id;
						if (assignment.g_zoneid_to_zone_seq_no_mapping.find(zone_org_id) != assignment.g_zoneid_to_zone_seq_no_mapping.end())
							link.zone_seq_no_for_outgoing_connector = assignment.g_zoneid_to_zone_seq_no_mapping[zone_org_id];
					}
				}

				double length_in_meter = 1.0; // km or mile
				double free_speed = 60.0;
				double cutoff_speed = 1.0;
				double k_jam = assignment.g_LinkTypeMap[link.link_type].k_jam;
				double bwtt_speed = 12;  //miles per hour

				double lane_capacity = 1800;
				parser_link.GetValueByFieldName("length", length_in_meter);  // in meter
				parser_link.GetValueByFieldName("FT", link.FT, false, true);
				parser_link.GetValueByFieldName("AT", link.AT, false, true);
				parser_link.GetValueByFieldName("vdf_code", link.vdf_code, false);

				if (length_in_km_waring == false && length_in_meter < 0.1)
				{
					dtalog.output() << "warning: link link_length =" << length_in_meter << " in link.csv for link " << from_node_id << "->" << to_node_id << ". Please ensure the unit of the link link_distance_VDF is meter." << endl;
					// link.link_type has been taken care by its default constructor
					length_in_km_waring = true;
				}

				if (length_in_meter < 1)
				{
					length_in_meter = 1;  // minimum link_distance_VDF
				}
				parser_link.GetValueByFieldName("free_speed", free_speed);  // free speed as km

				if (link.link_id == "201065AB")
				{
					int idebug = 1;
				}

				cutoff_speed = free_speed * 0.85; //default; 
				parser_link.GetValueByFieldName("cutoff_speed", cutoff_speed, false);


				if (free_speed <= 0.1)
					free_speed = 60;

				free_speed = max(0.1, free_speed);

				link.free_speed = free_speed;
				link.v_congestion_cutoff = cutoff_speed;



				double number_of_lanes = 1;
				parser_link.GetValueByFieldName("lanes", number_of_lanes);
				parser_link.GetValueByFieldName("capacity", lane_capacity);
				parser_link.GetValueByFieldName("saturation_flow_rate", link.saturation_flow_rate,false);

			
				link.free_flow_travel_time_in_min = length_in_meter/1000.0 / free_speed * 60;  // link_distance_VDF in meter 
				float fftt_in_sec = link.free_flow_travel_time_in_min * 60;  // link_distance_VDF in meter 

				link.length_in_meter = length_in_meter;
				link.link_distance_VDF = length_in_meter / 1000.0;
				link.link_distance_km = length_in_meter / 1000.0;
				link.link_distance_mile = length_in_meter / 1609.0;

				link.traffic_flow_code = assignment.g_LinkTypeMap[link.link_type].traffic_flow_code;
				link.number_of_lanes = number_of_lanes;
				link.lane_capacity = lane_capacity;

				//spatial queue and kinematic wave
				link.spatial_capacity_in_vehicles = max(1.0, link.link_distance_VDF * number_of_lanes * k_jam);

				// kinematic wave
				if (link.traffic_flow_code == 3)
					link.BWTT_in_simulation_interval = link.link_distance_VDF / bwtt_speed * 3600 / number_of_seconds_per_interval;

				link.vdf_type = assignment.g_LinkTypeMap[link.link_type].vdf_type;
				link.kjam = assignment.g_LinkTypeMap[link.link_type].k_jam;
				char VDF_field_name[50];

				for (int at = 0; at < assignment.g_AgentTypeVector.size(); at++)
				{
					double pce_at = -1; // default
					sprintf(VDF_field_name, "VDF_pce%s", assignment.g_AgentTypeVector[at].agent_type.c_str());

					parser_link.GetValueByFieldName(VDF_field_name, pce_at, false, true);

					if (pce_at > 1.001)  // log
					{
						//dtalog.output() << "link " << from_node_id << "->" << to_node_id << " has a pce of " << pce_at << " for agent type "
						//    << assignment.g_AgentTypeVector[at].agent_type.c_str() << endl;
					}


					for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
					{
						if (pce_at > 0)
						{   // if there are link based PCE values
							link.VDF_period[tau].pce[at] = pce_at;
						}
						else
						{
							link.VDF_period[tau].pce[at] = assignment.g_AgentTypeVector[at].PCE;

						}
						link.VDF_period[tau].occ[at] = assignment.g_AgentTypeVector[at].OCC;  //occ;
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
					link.VDF_period[tau].vdf_type = assignment.g_LinkTypeMap[link.link_type].vdf_type;
					link.VDF_period[tau].lane_based_ultimate_hourly_capacity = lane_capacity;
					link.VDF_period[tau].nlanes = number_of_lanes;

					link.VDF_period[tau].FFTT = link.link_distance_VDF / max(0.0001, link.free_speed) * 60.0;  // 60.0 for 60 min per hour
					link.VDF_period[tau].BPR_period_capacity = link.lane_capacity * link.number_of_lanes * assignment.g_DemandPeriodVector[tau].time_period_in_hour;
					link.VDF_period[tau].vf = link.free_speed;
					link.VDF_period[tau].v_congestion_cutoff = link.v_congestion_cutoff;
					link.VDF_period[tau].alpha = 0.15;
					link.VDF_period[tau].beta = 4;
					link.VDF_period[tau].preload = 0;

					for (int at = 0; at < assignment.g_AgentTypeVector.size(); at++)
					{
						link.VDF_period[tau].toll[at] = 0;
						link.VDF_period[tau].LR_price[at] = 0;
						link.VDF_period[tau].LR_RT_price[at] = 0;
					}

					link.VDF_period[tau].starting_time_in_hour = assignment.g_DemandPeriodVector[tau].starting_time_slot_no * MIN_PER_TIMESLOT / 60.0;
					link.VDF_period[tau].ending_time_in_hour = assignment.g_DemandPeriodVector[tau].ending_time_slot_no * MIN_PER_TIMESLOT / 60.0;
					link.VDF_period[tau].L = assignment.g_DemandPeriodVector[tau].time_period_in_hour;
					link.VDF_period[tau].t2 = assignment.g_DemandPeriodVector[tau].t2_peak_in_hour;
					link.VDF_period[tau].peak_load_factor = 1;

					int demand_period_id = assignment.g_DemandPeriodVector[tau].demand_period_id;
					sprintf(VDF_field_name, "VDF_fftt%d", demand_period_id);

					double FFTT = -1;
					parser_link.GetValueByFieldName(VDF_field_name, FFTT, false, false);  // FFTT should be per min

					bool VDF_required_field_flag = false;
					if (FFTT >= 0)
					{
						link.VDF_period[tau].FFTT = FFTT;
						VDF_required_field_flag = true;
					}

					if (link.VDF_period[tau].FFTT > 100)

					{
						dtalog.output() << "link " << from_node_id << "->" << to_node_id << " has a FFTT of " << link.VDF_period[tau].FFTT << " min at demand period " << demand_period_id
							<< " " << assignment.g_DemandPeriodVector[tau].demand_period.c_str() << endl;
					}

					sprintf(VDF_field_name, "VDF_cap%d", demand_period_id);
					parser_link.GetValueByFieldName(VDF_field_name, link.VDF_period[tau].BPR_period_capacity, VDF_required_field_flag, false);
					sprintf(VDF_field_name, "VDF_alpha%d", demand_period_id);
					parser_link.GetValueByFieldName(VDF_field_name, link.VDF_period[tau].alpha, VDF_required_field_flag, false);
					sprintf(VDF_field_name, "VDF_beta%d", demand_period_id);
					parser_link.GetValueByFieldName(VDF_field_name, link.VDF_period[tau].beta, VDF_required_field_flag, false);

					sprintf(VDF_field_name, "VDF_ref_link_volume%d", demand_period_id);
					parser_link.GetValueByFieldName(VDF_field_name, link.VDF_period[tau].ref_link_volume, false, false);

					

					sprintf(VDF_field_name, "VDF_allowed_uses%d", demand_period_id);
					if (parser_link.GetValueByFieldName(VDF_field_name, link.VDF_period[tau].allowed_uses, false) == false)
					{
						parser_link.GetValueByFieldName("allowed_uses", link.VDF_period[tau].allowed_uses, false);
					}

					sprintf(VDF_field_name, "VDF_preload%d", demand_period_id);
					parser_link.GetValueByFieldName(VDF_field_name, link.VDF_period[tau].preload, false, false);

					for (int at = 0; at < assignment.g_AgentTypeVector.size(); at++)
					{

						sprintf(VDF_field_name, "VDF_toll%s%d", assignment.g_AgentTypeVector[at].agent_type.c_str(), demand_period_id);
						parser_link.GetValueByFieldName(VDF_field_name, link.VDF_period[tau].toll[at], false, false);

if (link.VDF_period[tau].toll[at] > 0.001)
{
	dtalog.output() << "link " << from_node_id << "->" << to_node_id << " has a toll of " << link.VDF_period[tau].toll[at] << " for agent type "
		<< assignment.g_AgentTypeVector[at].agent_type.c_str() << " at demand period " << demand_period_id << endl;
}

					}
					sprintf(VDF_field_name, "VDF_penalty%d", demand_period_id);
					parser_link.GetValueByFieldName(VDF_field_name, link.VDF_period[tau].penalty, false, false);

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

						link.VDF_period[tau].red_time = max(1, link.VDF_period[tau].cycle_length - link.VDF_period[tau].effective_green_time);
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
				link.link_spatial_capacity = k_jam * number_of_lanes * link.link_distance_VDF;

				link.link_distance_VDF = max(0.00001, link.link_distance_VDF);
				for (int tau = 0; tau < assignment.g_number_of_demand_periods; ++tau)
					link.travel_time_per_period[tau] = link.link_distance_VDF / free_speed * 60;

				// min // calculate link cost based link_distance_VDF and speed limit // later we should also read link_capacity, calculate delay

				//int sequential_copying = 0;
				//
				//parser_link.GetValueByFieldName("sequential_copying", sequential_copying);

				g_node_vector[internal_from_node_seq_no].m_outgoing_link_seq_no_vector.push_back(link.link_seq_no);  // add this link to the corresponding node as part of outgoing node/link
				g_node_vector[internal_to_node_seq_no].m_incoming_link_seq_no_vector.push_back(link.link_seq_no);  // add this link to the corresponding node as part of outgoing node/link

				g_node_vector[internal_from_node_seq_no].m_to_node_seq_no_vector.push_back(link.to_node_seq_no);  // add this link to the corresponding node as part of outgoing node/link
				g_node_vector[internal_from_node_seq_no].m_to_node_2_link_seq_no_map[link.to_node_seq_no] = link.link_seq_no;  // add this link to the corresponding node as part of outgoing node/link


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

				g_link_vector.push_back(link);


				assignment.g_number_of_links++;

				if (assignment.g_number_of_links % 10000 == 0)
					dtalog.output() << "reading " << assignment.g_number_of_links << " links.. " << endl;
			

			}



			parser_link.CloseCSVFile();

			dtalog.output() << "number of links =" << g_link_vector.size() << endl;
		}

	 }

	dtalog.output() << "number of links =" << assignment.g_number_of_links << endl;
	
	assignment.summary_file << "Step 1: read network node.csv, link.csv, zone.csv "<< endl;
	assignment.summary_file << ",# of nodes = ," << g_node_vector.size() << endl;
	assignment.summary_file << ",# of links =," << g_link_vector.size() << endl;
	assignment.summary_file << ",# of zones =," << g_zone_vector.size() << endl;
	//assignment.summary_file << ",# of subarea links =," << bridge_link_count << endl;
	//assignment.summary_file << ",# of subarea nodes =," << number_of_micro_gate_nodes << endl;


	


	assignment.summary_file << ",summary by multi-modal and demand types,demand_period,agent_type,# of links,avg_free_speed_mph,avg_free_speed_kmph,total_length_in_km,total_capacity,avg_lane_capacity,avg_length_in_meter," << endl;
	for (int tau = 0; tau < assignment.g_DemandPeriodVector.size(); ++tau)
		for (int at = 0; at < assignment.g_AgentTypeVector.size(); at++)
		{
			assignment.summary_file << ",," << assignment.g_DemandPeriodVector[tau].demand_period.c_str() << ",";
			assignment.summary_file << assignment.g_AgentTypeVector[at].agent_type.c_str() << ",";

			int link_count = 0;
			double total_speed = 0;
			double total_length = 0;
			double total_lane_capacity = 0;
			double total_link_capacity = 0;

			for (int i = 0; i < g_link_vector.size(); i++)
			{
				if (g_link_vector[i].link_type >= 0 && g_link_vector[i].AllowAgentType(assignment.g_AgentTypeVector[at].agent_type, tau))
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
				total_speed / max(1, link_count)* 1.60934 << "," <<
				total_length/1000.0 << "," <<
				total_link_capacity << "," <<
				total_lane_capacity / max(1, link_count) << "," << total_length / max(1, link_count) << "," << endl;
		}
	/// <summary>
	/// ///////////////////////////
	/// </summary>
	/// <param name="assignment"></param>
	assignment.summary_file << ",summary by road link type,link_type,link_type_name,# of links,avg_free_speed_mph,avg_free_speed_kmph,,total_length,total_capacity,avg_lane_capacity,avg_length_in_meter," << endl;
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
			total_length/1000.0 << "," <<
			total_link_capacity << "," <<
			total_lane_capacity / max(1, link_count) << "," << total_length / max(1, link_count) << "," << endl;
	}


	//g_OutputModelFiles(2);
	//g_OutputModelFiles(3);

	if (dtalog.debug_level() == 2)
	{
		for (int i = 0; i < g_node_vector.size(); ++i)
		{
			if (g_node_vector[i].zone_org_id >= 1) // for each physical node
			{
				// we need to make sure we only create two way connectors between nodes and zones
				dtalog.output() << "node id= " << g_node_vector[i].node_id << " with zone id " << g_node_vector[i].zone_org_id << "and "
					<< g_node_vector[i].m_outgoing_link_seq_no_vector.size() << " outgoing links." << endl;
				for (int j = 0; j < g_node_vector[i].m_outgoing_link_seq_no_vector.size(); ++j)
				{
					int link_seq_no = g_node_vector[i].m_outgoing_link_seq_no_vector[j];
					dtalog.output() << "  outgoing node = " << g_node_vector[g_link_vector[link_seq_no].to_node_seq_no].node_id << endl;
				}
			}
			else
			{
				if (dtalog.debug_level() == 3)
				{
					dtalog.output() << "node id= " << g_node_vector[i].node_id << " with " << g_node_vector[i].m_outgoing_link_seq_no_vector.size() << " outgoing links." << endl;
					for (int j = 0; j < g_node_vector[i].m_outgoing_link_seq_no_vector.size(); ++j)
					{
						int link_seq_no = g_node_vector[i].m_outgoing_link_seq_no_vector[j];
						dtalog.output() << "  outgoing node = " << g_node_vector[g_link_vector[link_seq_no].to_node_seq_no].node_id << endl;
					}
				}
			}
		}
	}
	if (assignment.assignment_mode != 11)   // not tmc mode
		g_read_link_qvdf_data(assignment);
}
//CDTACSVParser parser_movement;
//int prohibited_count = 0;

//if (parser_movement.OpenCSVFile("movement.csv", false))  // not required
//{
//    while (parser_movement.ReadRecord())
//    {
//        string ib_link_id;
//        int node_id = 0;
//        string ob_link_id;
//        int prohibited_flag = 0;

//        if (!parser_movement.GetValueByFieldName("node_id", node_id))
//            break;

//        if (assignment.g_node_id_to_seq_no_map.find(node_id) == assignment.g_node_id_to_seq_no_map.end())
//        {
//            dtalog.output() << "Error: node_id " << node_id << " in file movement.csv is not defined in node.csv." << endl;
//            //has not been defined
//            continue;
//        }

//        parser_movement.GetValueByFieldName("ib_link_id", ib_link_id);
//        parser_movement.GetValueByFieldName("ob_link_id", ob_link_id);

//        if (assignment.g_link_id_map.find(ib_link_id) != assignment.g_link_id_map.end())
//            dtalog.output() << "Error: ib_link_id " << ib_link_id.c_str() << " has not been defined in movement.csv. Please check link.csv." << endl;

//        if (assignment.g_link_id_map.find(ob_link_id) != assignment.g_link_id_map.end())
//            dtalog.output() << "Error: ob_link_id " << ob_link_id.c_str() << " has not been defined in movement.csv. Please check link.csv." << endl;

//        float penalty = 0;
//        parser_movement.GetValueByFieldName("penalty", penalty);

//        if (penalty >= 99)
//        {
//            string	movement_string;
//            movement_string = ib_link_id + "->" + ob_link_id;

//            int node_no = assignment.g_node_id_to_seq_no_map[node_id];
//            g_node_vector[node_no].prohibited_movement_size++;
//            g_node_vector[node_no].m_prohibited_movement_string_map[movement_string] = 1;

//            prohibited_count++;
//        }
//    }

//    dtalog.output() << "Step XX: Reading movement.csv data with " << prohibited_count << " prohibited records." << endl;
//    parser_movement.CloseCSVFile();
//}



//
//
//void g_reload_timing_arc_data(Assignment& assignment)
//{
//    dtalog.output() << "Step 1.7: Reading service arc in timing.csv..." << endl;
//
//    CDTACSVParser parser_timing_arc;
//    if (parser_timing_arc.OpenCSVFile("timing.csv", false))
//    {
//        while (parser_timing_arc.ReadRecord())  // if this line contains [] mark, then we will also read field headers.
//        {
//            string mvmt_key;
//            if (!parser_timing_arc.GetValueByFieldName("mvmt_key", mvmt_key))
//            {
//                dtalog.output() << "Error: mvmt_key in file timing.csv is not defined." << endl;
//                continue;
//            }
//            // create a link object
//            CSignalTiming timing_arc;
//
//            if (assignment.g_mvmt_key_to_link_no_map.find(mvmt_key) == assignment.g_mvmt_key_to_link_no_map.end())
//            {
//                dtalog.output() << "Error: mvmt_key " << mvmt_key << " in file timing.csv is not defined in link.csv." << endl;
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
//                dtalog.output() << "Error: Field time_window in file timing.csv cannot be read." << endl;
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
//                dtalog.output() << "reading " << assignment.g_number_of_timing_arcs << " timing_arcs.. " << endl;
//        }
//
//        parser_timing_arc.CloseCSVFile();
//    }
//
//    dtalog.output() << endl;
//    dtalog.output() << "Step 1.8: Reading file section [demand_file_list] in setting.csv..." << endl;
//    // we now know the number of links
//    dtalog.output() << "number of timing records = " << assignment.g_number_of_timing_arcs << endl << endl;
//}
//


