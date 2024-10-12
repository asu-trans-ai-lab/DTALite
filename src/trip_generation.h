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
#include "DTA_geometry.h"
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

void g_add_new_access_link(int internal_from_node_seq_no, int internal_to_node_seq_no, float link_distance_VDF, int mode_type_no, int zone_seq_no = -1)
{
	// create a link object
	CLink link;
	link.allocate_memory();
	link.b_automated_generated_flag = true;
	link.from_node_seq_no = internal_from_node_seq_no;
	link.to_node_seq_no = internal_to_node_seq_no;
	link.link_seq_no = assignment.g_number_of_links;
	link.to_node_seq_no = internal_to_node_seq_no;


	link.link_type = -1;  // access_link

	//only for outgoing connectors
	link.zone_seq_no_for_outgoing_connector = zone_seq_no;

	//BPR
	link.traffic_flow_code = 0;

	link.spatial_capacity_in_vehicles = 99999;
	link.lane_capacity = 999999;
	link.link_spatial_capacity = 99999;
	link.link_distance_VDF = link_distance_VDF;
	link.free_speed = assignment.g_ModeTypeVector[mode_type_no].access_speed;


		//setup default values
		link.lane_based_ultimate_hourly_capacity = 99999;
		// 60.0 for 60 min per hour
		link.free_flow_travel_time_in_min = link_distance_VDF / max(0.001, link.free_speed) * 60;
		link.penalty = 99;
		link.alpha = 0;
		link.beta = 0;

		for (int at = 0; at < assignment.g_ModeTypeVector.size(); at++)
		{
		link.link_avg_travel_time[at] = link.free_flow_travel_time_in_min;
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

double get_activity_node_distance(Assignment& assignment)
{


	// calculate avg near by distance;
	double total_near_by_distance = 0;
	int activity_node_count = 0;
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

void g_grid_zone_generation(Assignment& assignment)
{
	int number_of_is_boundary = 0;
	int number_of_nodes = 0;

	CDTACSVParser parser;
	if (parser.OpenCSVFile("node.csv", true))
	{
		while (parser.ReadRecord())  // if this line contains [] mark, then we will also read field headers.
		{

			int node_id;
			if (!parser.GetValueByFieldName("node_id", node_id))
				continue;

			int is_boundary = 0;


			parser.GetValueByFieldName("is_boundary", is_boundary, false, false);


			if (is_boundary != 0)
				number_of_is_boundary++;

			number_of_nodes++;

		}

		// fprintf(g_pFileOutputLog, "number of nodes =,%d\n", assignment.g_number_of_nodes);
		parser.CloseCSVFile();
	}

	int sampling_rate = 10;
	if (number_of_is_boundary <= 1)       // random generation of activity locations // no boundary nodes
	{
		if (number_of_nodes > 1000)
			sampling_rate = number_of_nodes / 100;
	}


	std::vector<CNode> l_node_vector; // l as local

	if (parser.OpenCSVFile("node.csv", true))
	{
		while (parser.ReadRecord())  // if this line contains [] mark, then we will also read field headers.
		{

			int node_id;
			if (!parser.GetValueByFieldName("node_id", node_id))
				continue;

			CNode node;
			node.node_id = node_id;
			int is_boundary = 0;


			parser.GetValueByFieldName("is_boundary", is_boundary, false, false);

			if (number_of_is_boundary <= 1)       // random generation of activity locations // no boundary nodes
			{
				if (l_node_vector.size() % sampling_rate == 0)
				{
					is_boundary = 2;
				}
			}

			if (is_boundary != 0)
				node.is_activity_node = 1;

			parser.GetValueByFieldName("x_coord", node.x, true, false);
			parser.GetValueByFieldName("y_coord", node.y, true, false);

			l_node_vector.push_back(node);

		}

		// fprintf(g_pFileOutputLog, "number of nodes =,%d\n", assignment.g_number_of_nodes);
		parser.CloseCSVFile();
	}

	dtalog.output() << "[PROCESS INFO] Step 1.4.1: mode for creating node 2 zone mapping" << '\n';
	g_DTA_log_file << "[PROCESS INFO] Step 1.4.1: mode for creating node 2 zone mapping" << '\n';

	FILE* g_pFileZone = nullptr;
	g_pFileZone = fopen("zone.csv", "w");

	if (g_pFileZone == NULL)
	{
		dtalog.output() << "[WARNING] File zone.csv cannot be opened." << '\n';
		g_DTA_log_file << "[WARNING] File zone.csv cannot be opened." << '\n';
		return;
	}

	fprintf(g_pFileZone, "first_column,zone_id,access_node_vector,cell_code,cell_id,access_distance,x_coord,y_coord,");
	fprintf(g_pFileZone, "geometry,");

	fprintf(g_pFileZone, "production,attraction,");

	fprintf(g_pFileZone, "\n");



	double activity_nearbydistance = get_activity_node_distance(assignment);
	// initialization of grid rectangle boundary
	double left = 100000000;
	double right = -100000000;
	double top = -1000000000;
	double  bottom = 1000000000;

	for (int i = 0; i < l_node_vector.size(); i++)
	{
		// exapnd the grid boundary according to the nodes
		left = min(left, l_node_vector[i].x);
		right = max(right, l_node_vector[i].x);
		top = max(top, l_node_vector[i].y);
		bottom = min(bottom, l_node_vector[i].y);

	}

	int grid_size = 8;

	if (l_node_vector.size() > 3000)
		grid_size = 10;


	double temp_resolution = (((right - left) / grid_size + (top - bottom) / grid_size)) / 2.0;

	//if (activity_nearbydistance * 4 < temp_resolution)
	//{
	//    temp_resolution = activity_nearbydistance * 40;

	//}


	vector<double> ResolutionVector;

	ResolutionVector.push_back(0.00005);
	ResolutionVector.push_back(0.0001);
	ResolutionVector.push_back(0.0002);
	ResolutionVector.push_back(0.00025);
	ResolutionVector.push_back(0.0005);
	ResolutionVector.push_back(0.00075);
	ResolutionVector.push_back(0.001);
	ResolutionVector.push_back(0.002);
	ResolutionVector.push_back(0.0025);
	ResolutionVector.push_back(0.005);
	ResolutionVector.push_back(0.0075);
	ResolutionVector.push_back(0.01);
	ResolutionVector.push_back(0.02);
	ResolutionVector.push_back(0.025);
	ResolutionVector.push_back(0.05);
	ResolutionVector.push_back(0.075);
	ResolutionVector.push_back(0.1);
	ResolutionVector.push_back(0.2);
	ResolutionVector.push_back(0.25);
	ResolutionVector.push_back(0.5);
	ResolutionVector.push_back(0.75);
	ResolutionVector.push_back(1);
	ResolutionVector.push_back(2);
	ResolutionVector.push_back(2.5);
	ResolutionVector.push_back(5);
	ResolutionVector.push_back(7.5);
	ResolutionVector.push_back(10);
	ResolutionVector.push_back(20);
	ResolutionVector.push_back(25);
	ResolutionVector.push_back(50);
	ResolutionVector.push_back(75);

	double ClosestResolution = 1;

	if (temp_resolution < ResolutionVector[0])
		temp_resolution = ResolutionVector[0];

	for (unsigned int i = 0; i < ResolutionVector.size() - 1; i++)
	{
		if ((temp_resolution > ResolutionVector[i] + 0.000001) && temp_resolution < ResolutionVector[i + 1])
		{
			temp_resolution = ResolutionVector[i + 1]; // round up
			break;

		}
	}

	assignment.m_GridResolution = temp_resolution;

	dtalog.output() << "[PROCESS INFO] Step 1.4.2: Grid Resolution " << assignment.m_GridResolution << '\n';
	g_DTA_log_file << "[PROCESS INFO] Step 1.4.2: Grid Resolution " << assignment.m_GridResolution << '\n';

	int activity_node_count = 0;
	// check # of access nodes
	// check # of nodes

	int total_number_of_zones = 0;
	float total_number_of_trips_expected_lower_bound = 5000;

	if (l_node_vector.size() > 100000) //  regional model
	{
		total_number_of_trips_expected_lower_bound = 10000;
	}


	std::map<__int64, int> local_cell_id_2_zone_mapping;

	for (int i = 0; i < l_node_vector.size(); i++)
	{

		if (l_node_vector[i].is_activity_node != 0)
		{

			__int64 cell_id = g_get_cell_ID(l_node_vector[i].x, l_node_vector[i].y, assignment.m_GridResolution);
			if (local_cell_id_2_zone_mapping.find(cell_id) == local_cell_id_2_zone_mapping.end())  // create a cell
			{
				//create zone
				local_cell_id_2_zone_mapping[cell_id] = l_node_vector[i].node_id;
			}
		}
	}

	total_number_of_zones = local_cell_id_2_zone_mapping.size();
	int number_of_activity_nodes = 0;
	for (int i = 0; i < l_node_vector.size(); i++)
	{

		if (l_node_vector[i].is_activity_node != 0)
		{

			number_of_activity_nodes++;
		}
	}

	float production_rate_per_activity_node = total_number_of_trips_expected_lower_bound / max(1, number_of_activity_nodes);


	for (int i = 0; i < l_node_vector.size(); i++)
	{

		if (l_node_vector[i].is_activity_node != 0)
		{

			if (l_node_vector[i].node_id == 966)
			{
				int itest = 1;
			}
			__int64 cell_id = g_get_cell_ID(l_node_vector[i].x, l_node_vector[i].y, assignment.m_GridResolution);
			if (assignment.cell_id_mapping.find(cell_id) == assignment.cell_id_mapping.end())  // create a cell
			{
				//create zone
				assignment.cell_id_mapping[cell_id] = l_node_vector[i].node_id;
				string cell_code = g_get_cell_code(l_node_vector[i].x, l_node_vector[i].y, assignment.m_GridResolution, left, top);


				int x_i = floor(l_node_vector[i].x / assignment.m_GridResolution);
				int y_i = floor(l_node_vector[i].y / assignment.m_GridResolution);

				double x_coord_left = x_i * assignment.m_GridResolution;
				double y_coord_bottom = y_i * assignment.m_GridResolution;
				double x_coord_right = x_coord_left + assignment.m_GridResolution;
				double y_coord_top = y_coord_bottom + assignment.m_GridResolution;

				fprintf(g_pFileZone, "0,");
				fprintf(g_pFileZone, "%ld,", assignment.cell_id_mapping.size() + 1);
				// generate access nodes
				std::vector <int> access_node_vector;
				float max_distance = 0;

				float zone_x = (l_node_vector[i].x * 2 + x_coord_left + x_coord_right) / 4;
				float zone_y = (l_node_vector[i].y * 2 + y_coord_top + y_coord_bottom) / 4;
				for (int j = 0; j < l_node_vector.size(); j++)
				{

					if (l_node_vector[j].is_activity_node != 0)  // is_boundary flag
					{

						__int64 cell_id_j = g_get_cell_ID(l_node_vector[j].x, l_node_vector[j].y, assignment.m_GridResolution);
						if (cell_id == cell_id_j)
						{
							double distance = g_calculate_p2p_distance_in_meter_from_latitude_longitude(zone_x, zone_y, l_node_vector[j].x, l_node_vector[j].y);
							if (distance > max_distance)
							{
								max_distance = distance;
							}

							access_node_vector.push_back(l_node_vector[j].node_id);
							fprintf(g_pFileZone, "%d;", l_node_vector[j].node_id);
						}
					}
				}

				fprintf(g_pFileZone, ",%s,%d,", cell_code.c_str(), assignment.cell_id_mapping[cell_id]);
				fprintf(g_pFileZone, "%f,", max_distance);


				fprintf(g_pFileZone, "%f,%f,", (l_node_vector[i].x * 2 + x_coord_left + x_coord_right) / 4, (l_node_vector[i].y * 2 + y_coord_top + y_coord_bottom) / 4);

				fprintf(g_pFileZone, "\"LINESTRING (");

				fprintf(g_pFileZone, "%f %f,", x_coord_left, y_coord_top);
				fprintf(g_pFileZone, "%f %f,", x_coord_right, y_coord_top);
				fprintf(g_pFileZone, "%f %f,", x_coord_right, y_coord_bottom);
				fprintf(g_pFileZone, "%f %f,", x_coord_left, y_coord_bottom);
				fprintf(g_pFileZone, "%f %f,", x_coord_left, y_coord_top);
				fprintf(g_pFileZone, ")\",");

				for (int at = 0; at < assignment.g_ModeTypeVector.size(); ++at)
				{
					fprintf(g_pFileZone, "%.2f,%.2f,", access_node_vector.size() * production_rate_per_activity_node, access_node_vector.size() * production_rate_per_activity_node);
					break;  // one agent tyu
				}


				fprintf(g_pFileZone, "\n");

			}




		}
	}

	dtalog.output() << "[PROCESS INFO] Step 1.4.3: creating " << assignment.cell_id_mapping.size() << " zones." << '\n';
	g_DTA_log_file << "[PROCESS INFO] Step 1.4.3: creating " << assignment.cell_id_mapping.size() << " zones." << '\n';
	fclose(g_pFileZone);

}

int InsidePolygon(std::vector<CCoordinate> polygon, CCoordinate p)
{
	int N = polygon.size();
	int counter = 0;
	int i;
	double xinters;
	CCoordinate p1, p2;

	p1 = polygon[0];
	for (i = 1; i <= N; i++) {
		p2 = polygon[i % N];
		if (p.Y > min(p1.Y, p2.Y)) {
			if (p.Y <= max(p1.Y, p2.Y)) {
				if (p.X <= max(p1.X, p2.X)) {
					if (p1.Y != p2.Y) {
						xinters = (p.Y - p1.Y) * (p2.X - p1.X) / (p2.Y - p1.Y) + p1.X;
						if (p1.X == p2.X || p.X <= xinters)
							counter++;
					}
				}
			}
		}
		p1 = p2;
	}

	if (counter % 2 == 0)
		return 0;
	else
		return 1;
}


bool g_TAZ_2_GMNS_zone_generation(Assignment& assignment)
{


	// step 1: read TAZ.csv. zone_id x, y
	int number_of_nodes = 0;

	std::vector<CNode> l_TAZ_vector; // l as local

	CDTACSVParser parser;
	if (parser.OpenCSVFile("TAZ.csv", true))
	{
		// fprintf(g_pFileOutputLog, "number of nodes =,%d\n", assignment.g_number_of_nodes);
		while (parser.ReadRecord())  // if this line contains [] mark, then we will also read field headers.
		{

			int zone_id;
			if (!parser.GetValueByFieldName("zone_id", zone_id))
				continue;

			CNode node;
			node.zone_org_id = zone_id;

			parser.GetValueByFieldName("x_coord", node.x, true, false);
			parser.GetValueByFieldName("y_coord", node.y, true, false);

			string geo_string;
			parser.GetValueByFieldName("geometry", geo_string);

			CDTAGeometry geometry(geo_string);

			std::vector<CCoordinate> CoordinateVector = geometry.GetCoordinateList();
			node.zone_coordinate_vector = CoordinateVector;

			l_TAZ_vector.push_back(node);
		}
		parser.CloseCSVFile();
	}
	else
	{
		assignment.summary_file << "[NOTE] If zone.csv has not been prepared, please prepare TAZ.csv from the planning model so that DTALite can help with generating access nodes that connect zone to nodes in the network." << '\n';
		//dtalog.output() << "If zone.csv has not been prepared, please prepare TAZ.csv from the planning model so that DTALite can help with generating access nodes that connect zone to nodes in the network." << '\n';
		//g_DTA_log_file << "If zone.csv has not been prepared, please prepare TAZ.csv from the planning model so that DTALite can help with generating access nodes that connect zone to nodes in the network." << '\n';

		//g_program_stop();
		return false;
	}

	assignment.summary_file << "[DATA INFO] # of zones defined in TAZ.csv=,"<< l_TAZ_vector.size() << '\n';
	dtalog.output() << "[DATA INFO] # of zones defined in TAZ.csv=," << l_TAZ_vector.size() << '\n';
	g_DTA_log_file << "[DATA INFO] # of zones defined in TAZ.csv=," << l_TAZ_vector.size() << '\n';

	std::vector<CNode> l_node_vector; // l as local

	if (parser.OpenCSVFile("node.csv", true))
	{
		while (parser.ReadRecord())  // if this line contains [] mark, then we will also read field headers.
		{

			int node_id;
			if (!parser.GetValueByFieldName("node_id", node_id))
				continue;

			CNode node;

			node.node_id = node_id;
			int is_boundary = 0;

			parser.GetValueByFieldName("is_boundary", node.is_boundary, false, false);

			parser.GetValueByFieldName("x_coord", node.x, true, false);
			parser.GetValueByFieldName("y_coord", node.y, true, false);

			l_node_vector.push_back(node);

		}

		// fprintf(g_pFileOutputLog, "number of nodes =,%d\n", assignment.g_number_of_nodes);
		parser.CloseCSVFile();
	}



	dtalog.output() << "[PROCESS INFO] Step 1.4.0: QEM mode for creating TAZ 2 zone mapping with " << l_TAZ_vector.size() << " TAZ and " << l_node_vector.size() << " nodes." << '\n';
	g_DTA_log_file << "[PROCESS INFO] Step 1.4.0: QEM mode for creating TAZ 2 zone mapping with " << l_TAZ_vector.size() << " TAZ and " << l_node_vector.size() << " nodes." << '\n';

	// step 3:
	FILE* g_pFileZone = nullptr;
	g_pFileZone = fopen("zone.csv", "w");

	if (g_pFileZone == NULL)
	{
		dtalog.output() << "[WARNING] File zone.csv cannot be opened." << '\n';
		g_DTA_log_file << "[WARNING] File zone.csv cannot be opened." << '\n';
		return false;
	}


	double min_distance_without_range = 9999999;
	int min_distance_node_id_without_range = -1;


	fprintf(g_pFileZone, "first_column,zone_id,access_node_vector,access_distance_vector,x_coord,y_coord,access_link_geometry,\n");

	int access_node_size = 0;
	for (int taz = 0; taz < l_TAZ_vector.size(); taz++)
	{

		std::vector<int> access_node_seq_vector;
		std::vector<float> access_node_distance_vector;

		double zone_x = l_TAZ_vector[taz].x;
		double zone_y = l_TAZ_vector[taz].y;

		double min_distance_cutoff = 2000; //1KM

				//stage 1: is boundary nodes
		for (int i = 0; i < l_node_vector.size(); i++)
		{
			if (l_node_vector[i].is_boundary != 0)
			{
				//test                                double near_by_distance_1 = g_calculate_p2p_distance_in_meter_from_latitude_longitude(-77.429293, 39.697895, -77.339847, 38.947676);

				double distance = g_calculate_p2p_distance_in_meter_from_latitude_longitude(zone_x, zone_y, l_node_vector[i].x, l_node_vector[i].y);
				// calculate the distance

				if (distance < min_distance_cutoff)
				{
					access_node_distance_vector.push_back(distance);
				}

			}

		}
		// stage 2: default to all nearby nodes
		if (access_node_distance_vector.size() > 0)
		{
			std::sort(access_node_distance_vector.begin(), access_node_distance_vector.end());
			int acecss_link_k = 4;
			double distance_k_cut_off_value = access_node_distance_vector[max(0, acecss_link_k - 1)];

			access_node_distance_vector.clear();

			for (int i = 0; i < l_node_vector.size(); i++)
			{
				if (l_node_vector[i].is_boundary != 0)
				{
					//test                                double near_by_distance_1 = g_calculate_p2p_distance_in_meter_from_latitude_longitude(-77.429293, 39.697895, -77.339847, 38.947676);

					double distance = g_calculate_p2p_distance_in_meter_from_latitude_longitude(zone_x, zone_y, l_node_vector[i].x, l_node_vector[i].y);
					// calculate the distance

					if (distance < distance_k_cut_off_value)
					{
						access_node_seq_vector.push_back(i);
						access_node_distance_vector.push_back(distance);
					}

				}

			}


		}
		else
		{
			min_distance_cutoff = 2000; //1KM

			for (int i = 0; i < l_node_vector.size(); i++)
			{
				double distance = g_calculate_p2p_distance_in_meter_from_latitude_longitude(zone_x, zone_y, l_node_vector[i].x, l_node_vector[i].y);
				// calculate the distance
				if (distance < min_distance_cutoff)
				{
					access_node_distance_vector.push_back(distance);
				}

			}

			if (access_node_distance_vector.size() > 1)
			{

				std::sort(access_node_distance_vector.begin(), access_node_distance_vector.end());
				int acecss_link_k = 4;
				double distance_k_cut_off_value = access_node_distance_vector[max(0, acecss_link_k - 1)];
				access_node_distance_vector.clear();

				for (int i = 0; i < l_node_vector.size(); i++)
				{

					double distance = g_calculate_p2p_distance_in_meter_from_latitude_longitude(zone_x, zone_y, l_node_vector[i].x, l_node_vector[i].y);
					// calculate the distance

					if (distance < distance_k_cut_off_value)
					{
						access_node_seq_vector.push_back(i);
						access_node_distance_vector.push_back(distance);
					}


				}
			}


		}


		fprintf(g_pFileZone, "0,%d,", l_TAZ_vector[taz].zone_org_id);

		for (int n = 0; n < access_node_seq_vector.size(); n++)
		{
			int node_seq_no = access_node_seq_vector[n];
			fprintf(g_pFileZone, "%d;", l_node_vector[node_seq_no].node_id);

		}

		fprintf(g_pFileZone, ",");

		for (int n = 0; n < access_node_seq_vector.size(); n++)
		{
			fprintf(g_pFileZone, "%.3f;", access_node_distance_vector[n]);

		}

		fprintf(g_pFileZone, ",%f,%f,", zone_x, zone_y);

		fprintf(g_pFileZone, "\"MULTILINESTRING  (");

		access_node_size += access_node_seq_vector.size();
		for (int n = 0; n < access_node_seq_vector.size(); n++)
		{
			fprintf(g_pFileZone, "(");

			int node_seq_no = access_node_seq_vector[n];

			fprintf(g_pFileZone, "%f %f,", zone_x, zone_y);
			fprintf(g_pFileZone, "%f %f,", l_node_vector[node_seq_no].x, l_node_vector[node_seq_no].y);
			fprintf(g_pFileZone, ") ,");
		}
		fprintf(g_pFileZone, ")\"");

		fprintf(g_pFileZone, "\n");

	}

	dtalog.output() << "[PROCESS INFO] Step 1.4.3: creating " << l_TAZ_vector.size() << " zones." << '\n';
	g_DTA_log_file << "[PROCESS INFO] Step 1.4.3: creating " << l_TAZ_vector.size() << " zones." << '\n';
	fclose(g_pFileZone);

	assignment.summary_file << "[DATA INFO] # of zones created based on zone definition in TAZ.csv=," << l_TAZ_vector.size() << '\n';
	assignment.summary_file << "[DATA INFO] # of access nodes created based on node.csv and zone definition in TAZ.csv=," << access_node_size << '\n';


	return true;
}




void g_create_zone_vector(Assignment& assignment)
{
	std::map<int, int> waring_message_link_type_map;
	// initialize zone vector
	dtalog.output() << "[PROCESS INFO] Step 1.5: Initializing O-D zone vector..." << '\n';
	g_DTA_log_file << "[PROCESS INFO] Step 1.5: Initializing O-D zone vector..." << '\n';

	std::map<int, int>::iterator it;

	for (it = assignment.zone_id_2_node_no_mapping.begin(); it != assignment.zone_id_2_node_no_mapping.end(); ++it)
	{
		COZone ozone;

		if (it->first == 966)
		{
			int itest = 1;
		}
		// for each zone, we have to also create centriod
		ozone.zone_id = it->first;  // zone_id
		ozone.cell_id = assignment.zone_id_2_cell_id_mapping[it->first];
		ozone.cell_code = assignment.cell_id_2_cell_code_mapping[ozone.cell_id];
		ozone.zone_seq_no = g_zone_vector.size();
		ozone.cell_x = g_node_vector[it->second].x;
		ozone.cell_y = g_node_vector[it->second].y;

		dtalog.output() << "[DATA INFO] create zone id = " << ozone.zone_id << " with representive node id " << it->second << ",x = " << g_node_vector[it->second].x << ",y=" <<
			ozone.cell_y << '\n';

		g_DTA_log_file << "[DATA INFO] create zone id = " << ozone.zone_id << " with representive node id " << it->second << ",x = " << g_node_vector[it->second].x << ",y=" <<
			ozone.cell_y << '\n';

		assignment.g_zoneid_to_zone_seq_no_mapping[ozone.zone_id] = ozone.zone_seq_no;  // create the zone id to zone seq no mapping

		// create a centriod
		CNode node;
		// very large number as a special id
		node.node_id = -1 * (ozone.zone_id) - 1000000;
		node.node_seq_no = g_node_vector.size();
		assignment.g_node_id_to_seq_no_map[node.node_id] = node.node_seq_no;
		node.zone_id = ozone.zone_id;
		// push it to the global node vector
		g_node_vector.push_back(node);
		assignment.g_number_of_nodes++;

		ozone.node_seq_no = node.node_seq_no;
		// this should be the only one place that defines this mapping
		assignment.zone_id_to_centriod_node_no_mapping[ozone.zone_id] = node.node_seq_no;
		// add element into vector
		g_zone_vector.push_back(ozone);
	}

}

// This function calculates the distance between all pair of zones, and based on this,
// it calculates and assigns trips between these zones
void g_trip_generation(Assignment& assignment)
{

	// Accessibility Calculation
	// Iterate through all pairs of origin and destination zones
	for (int orig = 0; orig < g_zone_vector.size(); ++orig)  // o
	{
		for (int dest = 0; dest < g_zone_vector.size(); ++dest)  // d
		{
			// Avoid calculating distance from the zone to itself
			if (orig != dest)
			{
				// Calculate distance in miles between origin and destination zones
				float distance_in_mile = g_calculate_p2p_distance_in_meter_from_latitude_longitude(g_zone_vector[orig].cell_x, g_zone_vector[orig].cell_y, g_zone_vector[dest].cell_x, g_zone_vector[dest].cell_y) / 1609.0;

				// Store the calculated distance in the corresponding cell of the distance matrix
				g_zone_vector[orig].m_ODAccessibilityMatrix.distance_map[dest] = distance_in_mile;

				// Calculate travel time in minutes assuming the default speed is 30 miles per hour
				float travel_time_in_min = distance_in_mile / 30 * 60;

				// Store the calculated travel time in the corresponding cell of the travel time matrix
				g_zone_vector[orig].m_ODAccessibilityMatrix.value_map[dest] = travel_time_in_min;
			}
		}
	}


	FILE* g_pFileODMatrix = nullptr;
	fopen_ss(&g_pFileODMatrix, "gc_distance.csv", "w");

	if (!g_pFileODMatrix)
	{
		dtalog.output() << "File gc_distance.csv cannot be opened." << '\n';
		g_DTA_log_file << "File gc_distance.csv cannot be opened." << '\n';
		return; 
	}
	else
	{

		fprintf(g_pFileODMatrix, "o_zone_id,d_zone_id,distance,geometry\n");
		int demand_writing_log_count = 0;
		// reset the estimated production and attraction
		for (int orig = 0; orig < g_zone_vector.size(); ++orig)  // o
		{
			for (int dest = 0; dest < g_zone_vector.size(); ++dest)  // d
			{
				if (g_zone_vector[orig].gravity_production >= 0)
				{
					if (g_zone_vector[dest].gravity_attraction > 0)
					{
						float value = -1;
						if (g_zone_vector[orig].m_ODAccessibilityMatrix.distance_map.find(dest) != g_zone_vector[orig].m_ODAccessibilityMatrix.distance_map.end())
						{
							value = g_zone_vector[orig].m_ODAccessibilityMatrix.distance_map[dest];
						}

						if (value > 0.000001)
						{
							fprintf(g_pFileODMatrix, "%d,%d,%.4f,", g_zone_vector[orig].zone_id, g_zone_vector[dest].zone_id, value);

							fprintf(g_pFileODMatrix, "\"LINESTRING (");

							fprintf(g_pFileODMatrix, "%f %f,", g_zone_vector[orig].cell_x, g_zone_vector[orig].cell_y);
							fprintf(g_pFileODMatrix, "%f %f,", g_zone_vector[dest].cell_x, g_zone_vector[dest].cell_y);
							fprintf(g_pFileODMatrix, ")\"");
							fprintf(g_pFileODMatrix, "\n");
						}

					}
				}
			}

		}

		fclose(g_pFileODMatrix);
	}

	int out_of_bound_log_count = 0;
	int trip_accessibility_log_count = 0;
	int trip_distribution_log_count = 0;

	for (int orig = 0; orig < g_zone_vector.size(); ++orig)  // o
	 // gravity model;
	{
		if (g_zone_vector[orig].gravity_production > 0.00001)
		{
			float total_attraction_utility = 0;
			int count = 0;

			for (int d = 0; d < g_zone_vector.size(); ++d)
			{
				if (orig != d)
				{

					g_zone_vector[d].gravity_est_attraction = 0;

					//                            float cut_off = assignment.g_ModeTypeVector[at].trip_time_budget_in_min;
					float cut_off = 40;
					if (g_zone_vector[d].gravity_attraction > 0)
					{
						//double disutility = g_zone_vector[orig].m_ODAccessibilityMatrix.value_map[d] * beta;
						double exp_disutility = g_zone_vector[d].gravity_attraction;
						g_zone_vector[orig].m_ODAccessibilityMatrix.disutility_map[d] = exp_disutility;

						if (g_zone_vector[orig].m_ODAccessibilityMatrix.value_map[d] < cut_off)
						{
							if (trip_accessibility_log_count <= 100)
							{
								dtalog.output() << ", o: " << orig << ",d:" << d <<
									", gc distance = " << g_zone_vector[orig].m_ODAccessibilityMatrix.distance_map[d] <<
									", travel time =" << g_zone_vector[orig].m_ODAccessibilityMatrix.value_map[d] <<
									",value = " << exp_disutility << '\n';

								g_DTA_log_file << ", o: " << orig << ",d:" << d <<
									", gc distance = " << g_zone_vector[orig].m_ODAccessibilityMatrix.distance_map[d] <<
									", travel time =" << g_zone_vector[orig].m_ODAccessibilityMatrix.value_map[d] <<
									",value = " << exp_disutility << '\n';
							}
							total_attraction_utility += exp_disutility;
							trip_accessibility_log_count++;
							count++;
						}
						else
						{
							if (out_of_bound_log_count < 10)
							{
								dtalog.output() << "[DATA INFO] out of bound: " << ",o:" << orig << ",d:" << d <<
									", gc distance = " << g_zone_vector[orig].m_ODAccessibilityMatrix.distance_map[d] <<
									", travel time =" << g_zone_vector[orig].m_ODAccessibilityMatrix.value_map[d] <<
									",value = " << exp_disutility << '\n';

								g_DTA_log_file << "[DATA INFO] out of bound: " << ",o:" << orig << ",d:" << d <<
									", gc distance = " << g_zone_vector[orig].m_ODAccessibilityMatrix.distance_map[d] <<
									", travel time =" << g_zone_vector[orig].m_ODAccessibilityMatrix.value_map[d] <<
									",value = " << exp_disutility << '\n';

							}

							out_of_bound_log_count++;

						}
					}
				}
			}

			//dtalog.output() << "[DATA INFO] o: " << orig << ", total_attraction_utility =" << total_attraction_utility << '\n';
			//g_DTA_log_file << "[DATA INFO] o: " << orig << ", total_attraction_utility =" << total_attraction_utility << '\n';

			if (count > 0)
			{
				for (int dest = 0; dest < g_zone_vector.size(); ++dest)  // d
				{
					if (orig != dest)
					{
						if (g_zone_vector[dest].gravity_attraction > 0)
						{
							//                                    float cut_off = assignment.g_ModeTypeVector[at].trip_time_budget_in_min;

																//if (g_zone_vector[orig].m_ODAccessibilityMatrix.value_map[dest] < cut_off)
							{

								double ratio = g_zone_vector[orig].m_ODAccessibilityMatrix.disutility_map[dest] / total_attraction_utility;


								g_zone_vector[orig].m_ODMatrix.value_map[dest] = g_zone_vector[orig].gravity_production * ratio;

								if (trip_distribution_log_count < 100)
								{
									//dtalog.output() << ", o: " << orig << ",d:" << dest << ", ratio =" << ratio <<
									//g_DTA_log_file << ", o: " << orig << ",d:" << dest << ", ratio =" << ratio <<
									//	",trip = " << g_zone_vector[orig].m_ODMatrix.value_map[dest] << '\n';
								}
								trip_distribution_log_count++;

								g_zone_vector[dest].gravity_est_attraction += g_zone_vector[orig].m_ODMatrix.value_map[dest];
							}
						}

					}
				}
			}
		}
	}

}


void g_writing_demand_files(Assignment& assignment)
{
	dtalog.output() << "[STATUS INFO] writing demand.csv.." << '\n';
	g_DTA_log_file << "[STATUS INFO] writing demand.csv.." << '\n';

	FILE* g_pFileODMatrix = nullptr;
	fopen_ss(&g_pFileODMatrix, "demand.csv", "w");

	if (!g_pFileODMatrix)
	{
		dtalog.output() << "[ERROR] File demand.csv cannot be opened." << '\n';
		g_DTA_log_file << "[ERROR] File demand.csv cannot be opened." << '\n';
		return; 
	}
	else
	{

		fprintf(g_pFileODMatrix, "o_zone_id,d_zone_id,volume,geometry\n");
		int demand_writing_log_count = 0;
		// reset the estimated production and attraction
		for (int orig = 0; orig < g_zone_vector.size(); ++orig)  // o
		{
			for (int dest = 0; dest < g_zone_vector.size(); ++dest)  // d
			{
						if (g_zone_vector[orig].gravity_production >= 0)
						{
							if (g_zone_vector[dest].gravity_attraction > 0)
							{
								float value = 1;
								if (g_zone_vector[orig].m_ODMatrix.value_map.find(dest) != g_zone_vector[orig].m_ODMatrix.value_map.end())
								{
									value = g_zone_vector[orig].m_ODMatrix.value_map[dest];
								}

								if (value > 0.000001)
								{
									fprintf(g_pFileODMatrix, "%d,%d,%.4f,", g_zone_vector[orig].zone_id, g_zone_vector[dest].zone_id, value);

									demand_writing_log_count++;

									fprintf(g_pFileODMatrix, "\"LINESTRING (");

									fprintf(g_pFileODMatrix, "%f %f,", g_zone_vector[orig].cell_x, g_zone_vector[orig].cell_y);
									fprintf(g_pFileODMatrix, "%f %f,", g_zone_vector[dest].cell_x, g_zone_vector[dest].cell_y);
									fprintf(g_pFileODMatrix, ")\"");
									fprintf(g_pFileODMatrix, "\n");
								}

							}
					}

			}
		}

		fclose(g_pFileODMatrix);
	}


	return;

}

void g_demand_file_generation(Assignment& assignment)
{

	g_trip_generation(assignment);
	g_writing_demand_files(assignment);

}


void g_zone_to_access(Assignment& assignment)
{
	int debug_line_count = 0;
	if (assignment.assignment_mode == 21 && assignment.g_ModeTypeVector.size() > 0)  //zone2connector
	{
		int at = 0;

		if (assignment.g_ModeTypeVector[at].access_node_type.size() == 0)  // for each simple agent type
		{
			// find the closest zone id

			if (debug_line_count <= 20)
			{

				dtalog.output() << "[DATA INFO] connector generation condition 1: agent type " << assignment.g_ModeTypeVector[at].mode_type.c_str() << " has access node type" << assignment.g_ModeTypeVector[at].access_node_type.size() << '\n';
				g_DTA_log_file << "[DATA INFO] connector generation condition 1: agent type " << assignment.g_ModeTypeVector[at].mode_type.c_str() << " has access node type" << assignment.g_ModeTypeVector[at].access_node_type.size() << '\n';
				// zone without multimodal access
				debug_line_count++;
			}

			for (int a_k = 0; a_k < g_node_vector.size(); a_k++)
			{
				if (g_node_vector[a_k].is_activity_node == 100) //first loop for mode_specific activity node with zone id
				{

					int zone_id = g_node_vector[a_k].zone_org_id;
					int zone_seq_no = zone_seq_no = assignment.g_zoneid_to_zone_seq_no_mapping[zone_id];
					float access_distance = g_node_vector[a_k].access_distance;

					if (debug_line_count <= 20)
					{

						dtalog.output() << "[DATA INFO] connector generation generation condition 2: agent type no = " << at << " for node no. " << a_k << "as activity node with zone_id >=1" << '\n';
						g_DTA_log_file << "[DATA INFO] connector generation generation condition 2: agent type no = " << at << " for node no. " << a_k << "as activity node with zone_id >=1" << '\n';
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
						if (g_node_vector[i].is_activity_node == 2)  // is boundary
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

								if (distance <= access_distance * 1.01)  // check the range
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

						int acecss_link_k = 4;

						if (access_node_distance_vector.size() > acecss_link_k)
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
								if (g_node_vector[access_node_seq_vector[an]].is_boundary == 1)
									g_add_new_access_link(a_k, access_node_seq_vector[an], access_node_distance_vector[an], at, -1);

								if (g_node_vector[access_node_seq_vector[an]].is_boundary == -1)
									g_add_new_access_link(access_node_seq_vector[an], a_k, access_node_distance_vector[an], at, -1);

								assignment.g_ModeTypeVector[at].zone_id_cover_map[zone_id] = true;
							}
						}

					}

				}  // for each zone

			}

		}// for each agent type

		g_OutputModelFiles(3);  // node
		g_program_exit();

	}
}

