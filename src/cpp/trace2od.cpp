// trace2route.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

//Your code uses a function, class member, variable, or typedef that's marked deprecated. Symbols are deprecated by using a __declspec(deprecated) modifier, or the C++14 [[deprecated]] attribute. The actual C4996 warning message is specified
//by the deprecated modifier or attribute of the declaration.
//#include "MapMatching4GMNS.h"
#pragma warning(disable : 4996)

//#ifdef _WIN32
//#include "pch.h"
//#endif

#include "config.h"
#include "utils.h"


#include <algorithm>
#include <fstream>
#include <iostream>
#include <list>
#include <vector>
#include <map>
#include <cstdio>
//#include <omp.h>
//#include <time.h>
#include <ctime>
//#include <math.h>
#include <cmath>
using std::cout;
using std::endl;
using std::max;



#define MAX_LABEL_COST_ 999999
#define MAX_GRID_SIZE_ 1000


#define PI_ 3.1415
float g_NonHitDistanceRatio = 10;

float g_GridResolution = 0.05;  // in terms of long/lat 
float g_TimeResolution_inMin = 0.05;  //min --> 3 seconds
float g_SampleTimeResolution_inMin = 0.05; //min --> 3 seconds
float g_StartTimeinMin = 999999;
float g_EndTimeinMin = 0;
int g_TimeRangeInterval = 10;

int g_time_dependent_computing_mode = 0;

int g_max_number_of_threads = 4;
int g_number_of_nodes = 0;
int g_number_of_links = 0;
int g_number_of_agents = 0;
int g_grid_size = 1;

/*
constexpr auto g_max_number_of_threads = 4;
constexpr auto g_number_of_nodes = 0;
constexpr auto g_number_of_links = 0;
constexpr auto g_number_of_agents = 0;
constexpr auto g_grid_size = 1;
*/
using std::max;
using std::min;
using std::string;
using std::vector;
using std::map;
using std::ifstream;
using std::ofstream;
using std::istringstream;

//std::map<int, int> g_internal_node_seq_no_map;
//std::map<int, int> g_internal_link_no_map;
//std::map<string, int> g_internal_agent_no_map;

std::map<int, int> g_internal_node_seq_no_map;
std::map<string, int> g_internal_link_no_map;
std::map<string, int> g_internal_agent_no_map;
std::map<int, string> g_internal_trace_no_2_trace_id_map;

std::map<__int64, int> g_cell_id_2_zone_id_map;  // cell 2 zone mapping
std::map<__int64, int> g_cell_id_2_node_map;  // cell 2 node mapping

FILE* g_pFileLog = nullptr;


extern void g_Program_stop();
extern void g_OutputInputAgentCSVFile();

double g_left = 100000000;
double g_right = -100000000;
double g_top = -1000000000;
double g_bottom = 1000000000;

void g_Program_stop()
{

    cout << "Program stops. Press any key to terminate. Thanks!" << endl;
    getchar();
    exit(0);
};


__int64 g_GetCellID(double x, double y)
{
    __int64 xi;
    xi = floor(x / g_GridResolution);

    __int64 yi;
    yi = floor(y / g_GridResolution);

    return xi * 1000000 + yi;
};

int g_GetCellXID(double x, double x_min)
{
    int xi;

    xi = floor(x / g_GridResolution);

    int xmin;
    xmin = floor(x_min / g_GridResolution);

    return xi - xmin;
};

int g_GetCellYID(double y, double y_min)
{
    int yi;

    yi = floor(y / g_GridResolution);

    int ymin = floor(y_min / g_GridResolution);

    return yi - ymin;
};

int g_GetTimeInterval(float time_in_min)
{
    int time_interval = 0;
    time_interval = int((time_in_min - g_StartTimeinMin) / g_TimeResolution_inMin + 0.5);

    if (time_interval < 0)
        time_interval = 0;

    if (time_interval > g_TimeRangeInterval - 1)
        time_interval = g_TimeRangeInterval - 1;

    return time_interval;

}


float g_GetTimeInMinFromInterval(int time_interval)
{
    float time_in_min = g_StartTimeinMin + time_interval * g_TimeResolution_inMin;

    return time_in_min;

}

double g_findMedian(vector<int> a, int n)
{
    // source: https://www.geeksforgeeks.org/finding-median-of-unsorted-array-in-linear-time-using-c-stl/ 
     // If size of the arr[] is even
    if (n % 2 == 0) {

        // Applying nth_element
        // on n/2th index
        nth_element(a.begin(),
            a.begin() + n / 2,
            a.end());

        // Applying nth_element
        // on (n-1)/2 th index
        nth_element(a.begin(),
            a.begin() + (n - 1) / 2,
            a.end());

        // Find the average of value at
        // index N/2 and (N-1)/2
        return (double)(a[(n - 1) / 2]
            + a[n / 2])
            / 2.0;
    }

    // If size of the arr[] is odd
    else {

        // Applying nth_element
        // on n/2
        nth_element(a.begin(),
            a.begin() + n / 2,
            a.end());

        // Value at index (N/2)th
        // is the median
        return (double)a[n / 2];
    }
}
struct GDPoint //geometry data
{
    double x;
    double y;
};


double g_Find_P2P_Angle(const GDPoint* p1, const GDPoint* p2)
{
    double delta_x = p2->x - p1->x;
    double delta_y = p2->y - p1->y;

    if (fabs(delta_x) < 0.00001)
        delta_x = 0;

    if (fabs(delta_y) < 0.00001)
        delta_y = 0;

    int angle = atan2(delta_y, delta_x) * 180 / PI_ + 0.5;
    // angle = 90 - angle;

    while (angle < 0)
        angle += 360;

    while (angle > 360)
        angle -= 360;

    return angle;
}

double g_Find_PPP_RelativeAngle(const GDPoint* p1, const GDPoint* p2, const GDPoint* p3, const GDPoint* p4)
{
    int relative_angle;

    int angle1 = g_Find_P2P_Angle(p1, p2);
    int angle2 = g_Find_P2P_Angle(p3, p4);
    relative_angle = angle2 - angle1;

    while (relative_angle > 180)
        relative_angle -= 360;

    while (relative_angle < -180)
        relative_angle += 360;

    return relative_angle;
}

class CGeometry
{
public:
    enum GeometryType
    {
        POINT,
        LINE,
        POLYGON,
        UNKNOWN
    };

private:
    GeometryType m_Type;
    int m_NumOfCoordinates;
    std::vector<CCoordinate> v_Coordinates;
    bool ReadPointCoordinate(string s);
    bool ReadLineStringCoordinates(string s);
    bool ReadPolygonCoordinates(string s);

public:
    CGeometry(string s);
    //~CGeometry(void);
    ~CGeometry();

    //GeometryType GetGeometryType(void);
    //std::vector<CCoordinate> GetCoordinateList(void);
    //int GetNumberOfCoordinates(void);
    GeometryType GetGeometryType();
    std::vector<CCoordinate> GetCoordinateList();
    int GetNumberOfCoordinates();
};

CGeometry::CGeometry(string s)
{
    m_NumOfCoordinates = 0;

    string tmp;
    if (s.find("POINT") != std::string::npos)
    {
        tmp = s.substr(s.find_first_not_of(' '));
        size_t start_idx = tmp.find_first_of('(');
        size_t end_idx = tmp.find_first_of(')');

        if (start_idx == std::string::npos || end_idx == std::string::npos)
            return;

        string type_str = tmp.substr(0, start_idx);
        type_str.erase(type_str.find_last_not_of(" ") + 1); // works for 'LINESTRING (....' and 'LINESTRING(....'

        string start_tag = "(";
        string end_tag = ")";

        start_idx = tmp.find(start_tag);
        start_idx += start_tag.length();
        end_idx = tmp.find(end_tag);

        tmp = tmp.substr(start_idx, end_idx - start_idx);

        m_Type = POINT;
    }
    else if (s.find("LINESTRING") != std::string::npos)
    {
        tmp = s.substr(s.find_first_not_of(' '));
        size_t start_idx = tmp.find_first_of('(');
        size_t end_idx = tmp.find_first_of(')');

        if (start_idx == std::string::npos || end_idx == std::string::npos)
            return;

        string type_str = tmp.substr(0, start_idx);
        type_str.erase(type_str.find_last_not_of(" ") + 1); // works for 'LINESTRING (....' and 'LINESTRING(....'

        string start_tag = "(";
        string end_tag = ")";

        start_idx = tmp.find(start_tag);
        start_idx += start_tag.length();
        end_idx = tmp.find(end_tag);

        tmp = tmp.substr(start_idx, end_idx - start_idx);

        m_Type = LINE;
    }
    else if (s.find("POLYGON") != std::string::npos)
    {
        tmp = s.substr(s.find_first_not_of(' '));
        size_t start_idx = tmp.find('((');
        size_t end_idx = tmp.find('))');

        if (start_idx == std::string::npos || end_idx == std::string::npos)
            return;

        string type_str = tmp.substr(0, start_idx);
        type_str.erase(type_str.find_last_not_of(" ") + 1); // works for 'LINESTRING (....' and 'LINESTRING(....'

        string start_tag = "((";
        string end_tag = "))";

        start_idx = tmp.find(start_tag);
        start_idx += start_tag.length();
        end_idx = tmp.find(end_tag);

        tmp = tmp.substr(start_idx, end_idx - start_idx);

        m_Type = POLYGON;
    }
    else
    {
        m_Type = UNKNOWN;
    }

    switch (m_Type)
    {
    case POINT:
        ReadPointCoordinate(tmp);
        break;
    case LINE:
        ReadLineStringCoordinates(tmp);
        break;
    case POLYGON:
        ReadPolygonCoordinates(tmp);
        break;
    default:
        break;
    }
}

CGeometry::~CGeometry(void)
{
}

CGeometry::GeometryType CGeometry::GetGeometryType(void)
{
    return m_Type;
}

int CGeometry::GetNumberOfCoordinates(void)
{
    return m_NumOfCoordinates;
}

std::vector<CCoordinate> CGeometry::GetCoordinateList(void)
{
    return v_Coordinates;
}

bool CGeometry::ReadLineStringCoordinates(string s)
{
    istringstream ss(s);
    string sub_str;

    if (std::string::npos == s.find_first_of("0123456789"))
    {
        // "digit not found!, empty string//
        return false;
    }

    while (std::getline(ss, sub_str, ','))
    {
        sub_str = sub_str.substr(sub_str.find_first_not_of(' '));

        CCoordinate coordinate;
        std::istringstream sub_ss(sub_str);
        string tmp;

        std::getline(sub_ss, tmp, ' ');
        istringstream x_ss(tmp);
        x_ss >> coordinate.X;

        std::getline(sub_ss, tmp, ' ');
        istringstream y_ss(tmp);
        y_ss >> coordinate.Y;

        v_Coordinates.push_back(coordinate);
        m_NumOfCoordinates += 1;
    }
    return true;
}

bool CGeometry::ReadPolygonCoordinates(string s)
{
    istringstream ss(s);
    string sub_str;
    if (std::string::npos == s.find_first_of("0123456789"))
    {
        // "digit not found!, empty string//
        return false;
    }

    while (std::getline(ss, sub_str, ','))
    {
        sub_str = sub_str.substr(sub_str.find_first_not_of(' '));

        CCoordinate coordinate;
        istringstream sub_ss(sub_str);
        string tmp;

        std::getline(sub_ss, tmp, ' ');
        istringstream x_ss(tmp);
        x_ss >> coordinate.X;

        std::getline(sub_ss, tmp, ' ');
        istringstream y_ss(tmp);
        y_ss >> coordinate.Y;

        v_Coordinates.push_back(coordinate);
        m_NumOfCoordinates += 1;
    }
    return true;
}
bool CGeometry::ReadPointCoordinate(string s)
{
    CCoordinate coordinate;
    istringstream ss(s);

    string sub_str;
    std::getline(ss, sub_str, ' ');
    istringstream x_ss(sub_str);

    std::getline(ss, sub_str, ' ');
    istringstream y_ss(sub_str);
    x_ss >> coordinate.X;
    y_ss >> coordinate.Y;
    coordinate.Z = 0.0;

    v_Coordinates.push_back(coordinate);
    m_NumOfCoordinates = 1;

    return true;
}


class CMapmatchingNode
{
public:
    CMapmatchingNode()
    {
        bInsideFlag = false;

    }
    bool  bInsideFlag;

    int node_seq_no; // sequence number
    int node_id;     //external node number
    int zone_id;
    __int64 cell_id;
    string name;
    std::vector<int> m_outgoing_link_seq_no_vector;

    std::map<int, int> m_outgoing_link_seq_no_map;
    GDPoint pt;


};

class CMapmatchingLink
{
public:
    CMapmatchingLink()
    {
        lanes = 1;
        bInsideFlag = false;
        length = 1;
        FFTT_in_min = 0;
        FFTT_in_sec = 0;
        free_speed;
        x_key = 0;
        y_key = 0;
        hit_count = 0;
        likelihood_distance = 999999;

        o_distance = 999999;
        d_distance = 999999;

        AccessibilityTime = 999999;
        likely_trace_no = -1;
        use_count = 0;
        balance = 0;
        Possible_Dwell_time_in_min = 0;
        dwell_start_time_in_min = 0;
        lane_capacity = -1;
    }

    string link_id;
    __int64 cell_id;
    string name;
    string geometry;

    std::vector<GDPoint> m_PointVector;

    int from_node_id;
    int to_node_id;
    double length;
    double free_speed;
    double lane_capacity;
    int lanes;
    double FFTT_in_min;
    string link_type_code;
    string link_type_name;
    double Possible_Dwell_time_in_min;
    int dwell_start_time_in_min;
    int FFTT_in_sec;
    int x_key;
    int y_key;

    int link_seq_no;
    int from_node_seq_no;
    int to_node_seq_no;
    double likelihood_distance;
    int likely_trace_no;
    double o_distance;
    double d_distance;
    double link_distance;
    double seg_distance;
    int hit_count;
    int use_count;
    double balance;
    bool bInsideFlag;

    int AccessibilityTime;

};

std::vector<CMapmatchingNode> g_mm_node_vector;
std::vector<CMapmatchingLink> g_mm_link_vector;

double g_GetPoint2Point_Distance(const GDPoint* p1, const GDPoint* p2)
{
    return pow(((p1->x - p2->x) * (p1->x - p2->x) + (p1->y - p2->y) * (p1->y - p2->y)), 0.5);
}

double g_GetPoint2LineDistance(const GDPoint* pt, const GDPoint* FromPt, const GDPoint* ToPt, double UnitGridResolution, bool no_intersection_requirement)
{
    double U;
    GDPoint Intersection;

    double LineLength = g_GetPoint2Point_Distance(FromPt, ToPt);

    U = ((pt->x - ToPt->x) * (FromPt->x - ToPt->x) + (pt->y - ToPt->y) * (FromPt->y - ToPt->y)) / (LineLength * LineLength);

    if (no_intersection_requirement == false)
    {

        if (U < 0.0 || U > 1.0)
            return g_GridResolution; // intersection does not fall within the segment
    }
    Intersection.x = ToPt->x + U * (FromPt->x - ToPt->x);
    Intersection.y = ToPt->y + U * (FromPt->y - ToPt->y);

    double distance_1 = g_GetPoint2Point_Distance(pt, &Intersection);
    double distance_0 = g_GetPoint2Point_Distance(pt, FromPt);
    double distance_2 = g_GetPoint2Point_Distance(pt, ToPt);

    if (no_intersection_requirement)
    {
        return min(min(distance_1, distance_0), distance_2);
    }
    else
        return distance_1;
}

bool g_GetPoint2LineIntersectionFlag(const GDPoint* pt, const GDPoint* FromPt, const GDPoint* ToPt)
{
    double U;
    GDPoint Intersection;

    double LineLength = g_GetPoint2Point_Distance(FromPt, ToPt);

    U = ((pt->x - ToPt->x) * (FromPt->x - ToPt->x) + (pt->y - ToPt->y) * (FromPt->y - ToPt->y)) / (LineLength * LineLength);

    if (U < 0.0 || U > 1.0)
        return false;

    Intersection.x = ToPt->x + U * (FromPt->x - ToPt->x);
    Intersection.y = ToPt->y + U * (FromPt->y - ToPt->y);

    double distance_1 = g_GetPoint2Point_Distance(pt, &Intersection);
    double distance_0 = g_GetPoint2Point_Distance(pt, FromPt);
    double distance_2 = g_GetPoint2Point_Distance(pt, ToPt);


    if (distance_1 < 0.5 * distance_0 & distance_1 < 0.5 * distance_2 && distance_1 < 0.3 * LineLength)
        return true;
    else
        return false;
}

bool g_GetTwoPoints2LineIntersectionFlag(const GDPoint* pt0, const GDPoint* pt1, const GDPoint* FromPt, const GDPoint* ToPt)
{
    double U;
    GDPoint Intersection;
    GDPoint pt;
    pt.x = (pt0->x + pt1->x) / 2;
    pt.y = (pt0->y + pt1->y) / 2;

    double LineLength = g_GetPoint2Point_Distance(FromPt, ToPt);

    U = ((pt.x - ToPt->x) * (FromPt->x - ToPt->x) + (pt.y - ToPt->y) * (FromPt->y - ToPt->y)) / (LineLength * LineLength);

    if (U < -0.3 || U > 1.3)  // larger range
        return false;

    Intersection.x = ToPt->x + U * (FromPt->x - ToPt->x);
    Intersection.y = ToPt->y + U * (FromPt->y - ToPt->y);

    double distance_1 = g_GetPoint2Point_Distance(&pt, &Intersection);
    double distance_0 = g_GetPoint2Point_Distance(&pt, FromPt);
    double distance_2 = g_GetPoint2Point_Distance(&pt, ToPt);

    double relative_angle = fabs(g_Find_PPP_RelativeAngle(pt0, pt1, FromPt, ToPt));
    if (relative_angle > 45)
        return false;

    if (distance_1 < 1 * distance_0 & distance_1 < 1 * distance_2 && distance_1 < 1 * LineLength)
        //    if (distance_1 < 0.5 * distance_0 & distance_1 < 0.5 * distance_2 && distance_1 < 0.3 * LineLength)
        return true;
    else
        return false;
}

//ill conditioning detection
bool g_ill_conditioning_detection(double link_distance, double GPS_segment_distance)
{
    double cutoff_ratio = 3;
    double ratio = link_distance / max(0.00000001, GPS_segment_distance);
    if ((1.0 / cutoff_ratio) < ratio && ratio < cutoff_ratio)
        return false; //this indicates: good_conditioning, so ill_conditioning = false
    else
        return true;
}

class CGPSPoint
{
public:
    CGPSPoint()
    {
        bInsideGrid = false;
        Inside_index = -1;
        interval_in_second = 0;
        relative_time_in_second = 0;
        trace_no = -1;
    }

public:
    GDPoint pt;
    __int64 cell_id;
    int trace_no;

    //double time_interval_no;
    int dd;
    int global_time_in_second;
    int relative_time_in_second;
    float distance;
    float speed;
    int interval_in_second;

    bool bInsideGrid;
    int Inside_index;
};

class GridNodeSet
{
public:
    GridNodeSet()
    {
        origin_cell_flag = -1;
        destination_cell_flag = -1;
        possible_dwell_cell_flag = -1;
        zone_id = -1;
    }
    double x;
    double y;
    __int64 cell_id;

    int zone_id;
    std::vector<int> m_NodeVector;
    std::vector<int> m_LinkNoVector;
    std::vector<CGPSPoint> m_GPSPointVector;

    int origin_cell_flag;
    int destination_cell_flag;
    int possible_dwell_cell_flag;
    float start_min;
    float end_min;
    float dwell_time_min;

};
class CAgent
{
public:
    CAgent()
    {
        o_node_no = -1;
        d_node_no = -1;
        o_node_id = -1;
        d_node_id = -1;
        speed = -1;

        matching_link_no = -1;
        avg_GPS_segment_distance = 0;
        first_segment_distance = 0;
        last_segment_distance = 0;

        head_gps_index = -1;
        tail_gps_index = -1;

        o_cell_id = -1;
        d_cell_id = -1;

        origin_zone_id = -1;
        destination_zone_id = -1;
        distance = 0;
        travel_time = 0;
        sampling_rate_in_min = 0.5;
    }

    ~CAgent()
    {

    }


    int head_gps_index;
    int tail_gps_index;

    string agent_id;
    string agent_name;


    std::vector <int> likely_trace_no_vector;
    float volume;
    float speed;
    float distance;
    float travel_time;
    int agent_no;

    int o_node_no;
    int d_node_no;
    int o_node_id;
    int d_node_id;

    string allowed_link_type_code;
    string blocked_link_type_code;



    int matching_link_no;
    int origin_node_seq_no;
    int destination_node_seq_no;
    int origin_zone_id;
    int destination_zone_id;
    __int64 o_cell_id;
    __int64 d_cell_id;

    float start_time_in_min;
    float end_time_in_min;
    float duration_in_min;
    float sampling_rate_in_min;

    std::vector<CGPSPoint> m_GPSPointVector;
    std::map<int, int> m_ExcessiveDwellMap;  // time interval seq. no, time interval

    double avg_GPS_segment_distance;
    double first_segment_distance;
    double last_segment_distance;

    int upstream_node_matched_time_in_sec;

    int m_node_size;

    std::vector<int> path_link_vector, path_node_vector, path_time_vector;
    std::vector<float> path_cost_vector;

    int* path_link_matched_trace_id;    // for each link
    void AllocatePathNodeVector(int node_size, const int* node_vector, bool backwardflag = false)
    {
        m_node_size = node_size;

        if (backwardflag)
        {
            //copy backward
            for (int i = 0; i < m_node_size; ++i)
            {
                path_node_vector.push_back(node_vector[m_node_size - 1 - i]);
            }
        }

    }

    void AllocatePathNodeVector(int node_size, std::vector <int> node_vector, std::vector <int> time_index_vector, std::vector <float> cost_vector, std::vector <int> link_vector)
    {
        m_node_size = node_size;
        path_link_vector.clear();
        path_time_vector.clear();
        path_cost_vector.clear();
        path_node_vector.clear();

        //copy backward
        path_node_vector = node_vector;
        path_time_vector = time_index_vector;
        path_cost_vector = cost_vector;
        path_link_vector = link_vector;

        std::reverse(path_node_vector.begin(), path_node_vector.end());
        std::reverse(path_time_vector.begin(), path_time_vector.end());
        std::reverse(path_cost_vector.begin(), path_cost_vector.end());
        std::reverse(path_link_vector.begin(), path_link_vector.end());
    }

    GDPoint m_o_boundary_point[2]; // from point and to point
    GDPoint m_d_boundary_point[2];
};
vector<CAgent> g_agent_vector;

class NetworkForSP // mainly for shortest path calculation
{
public:
    NetworkForSP()
    {
    }

    GridNodeSet** m_GridMatrix; //important data structure for creating grid matrix
    double m_left;              // boundary of grid matrix
    double m_right;
    double m_top;
    double m_bottom;

    double m_GridXStep; // x resolution of grid cell
    double m_GridYStep;

    void BuildGridSystem()
    {

        m_GridMatrix = Allocate2DDynamicArray<GridNodeSet>(MAX_GRID_SIZE_, MAX_GRID_SIZE_);

        // initialization of grid rectangle boundary
        m_left = 100000000;
        m_right = -100000000;
        m_top = -1000000000;
        m_bottom = 1000000000;

        // exapnd the grid boundary according to the nodes
        for (int i = 0; i < g_mm_node_vector.size(); i++)
        {
            m_left = min(m_left, g_mm_node_vector[i].pt.x);
            m_right = max(m_right, g_mm_node_vector[i].pt.x);
            m_top = max(m_top, g_mm_node_vector[i].pt.y);
            m_bottom = min(m_bottom, g_mm_node_vector[i].pt.y);
        }


        m_GridXStep = g_GridResolution;
        m_GridYStep = g_GridResolution;

        g_grid_size = max((m_right - m_left) / g_GridResolution + 2, (m_top - m_bottom) / g_GridResolution + 2);

         cout << "g_GridResolution= " << g_GridResolution << ", grid size = " << g_grid_size << endl;


        // put nodes into grid cell

        for (int i = 0; i < g_mm_node_vector.size(); i++)
        {
            int x_key = g_GetCellXID(g_mm_node_vector[i].pt.x, m_left);
            int y_key = g_GetCellYID(g_mm_node_vector[i].pt.y, m_bottom);

            m_GridMatrix[x_key][y_key].m_NodeVector.push_back(i);

            if (m_GridMatrix[x_key][y_key].zone_id <= 0 && g_mm_node_vector[i].zone_id > 0)
            {
                m_GridMatrix[x_key][y_key].zone_id = g_mm_node_vector[i].zone_id;  // use the zone id for the first node in the cell

            }
        }

        // assign zone id for cell with nodes
        for (int x_i = 0; x_i <= g_grid_size; x_i++)
            for (int y_i = 0; y_i <= g_grid_size; y_i++)
            {

                if (m_GridMatrix[x_i][y_i].m_NodeVector.size() > 0 && m_GridMatrix[x_i][y_i].zone_id <= 0)
                {
                    // find near by zone id.

                    double min_distance = 10000000;
                    int min_distance_x_j = -1;
                    int min_distance_y_j = -1;

                    for (int x_j = 0; x_j <= g_grid_size; x_j++)
                        for (int y_j = 0; y_j <= g_grid_size; y_j++)
                        {

                            if (m_GridMatrix[x_j][y_j].m_NodeVector.size() > 0 && m_GridMatrix[x_j][y_j].zone_id >= 1)
                            {
                                double distance = pow(pow(x_i - x_j, 2) + pow(y_i - y_j, 2), 0.5);
                                if (distance < min_distance)
                                {
                                    min_distance = distance;
                                    min_distance_x_j = x_j;
                                    min_distance_y_j = y_j;

                                }

                            }
                        }

                    if (min_distance_x_j >= 0 && min_distance_y_j >= 0)
                    {
                        m_GridMatrix[x_i][y_i].zone_id = m_GridMatrix[min_distance_x_j][min_distance_y_j].zone_id;
                        fprintf(g_pFileLog, "cell [%d,%d] uses zone id=%d from cell [%d,%d]\n", x_i, y_i, m_GridMatrix[x_i][y_i].zone_id,
                            min_distance_x_j, min_distance_y_j);
                    }
                }
            }


        // if the node id still has no assigned zone id, use the corresponding zone id from its grid cell
        for (int i = 0; i < g_mm_node_vector.size(); i++)  // assign a zone id 
        {
            if (g_mm_node_vector[i].zone_id <= 0)
            {
                int x_key = g_GetCellXID(g_mm_node_vector[i].pt.x, m_left);
                int y_key = g_GetCellYID(g_mm_node_vector[i].pt.y, m_bottom);

                g_mm_node_vector[i].zone_id = m_GridMatrix[x_key][y_key].zone_id;


            }
        }
        // put links into grid cell

        for (int l = 0; l < g_mm_link_vector.size(); l++)
        {
            int x_key = g_GetCellXID(g_mm_node_vector[g_mm_link_vector[l].from_node_seq_no].pt.x, m_left);
            int y_key = g_GetCellYID(g_mm_node_vector[g_mm_link_vector[l].from_node_seq_no].pt.y, m_bottom);

            m_GridMatrix[x_key][y_key].m_LinkNoVector.push_back(l);
            g_mm_link_vector[l].cell_id = g_GetCellID(g_mm_node_vector[g_mm_link_vector[l].from_node_seq_no].pt.x, g_mm_node_vector[g_mm_link_vector[l].from_node_seq_no].pt.y);


            int from_x_key = x_key;
            int from_y_key = y_key;

            x_key = g_GetCellXID(g_mm_node_vector[g_mm_link_vector[l].to_node_seq_no].pt.x, m_left);
            y_key = g_GetCellYID(g_mm_node_vector[g_mm_link_vector[l].to_node_seq_no].pt.y, m_bottom);

            g_mm_link_vector[l].x_key = x_key;
            g_mm_link_vector[l].y_key = y_key;

            if (from_x_key != x_key || from_y_key != y_key) // when the from node and to node of a link belong to  different cells.
            {
                m_GridMatrix[x_key][y_key].m_LinkNoVector.push_back(l);
            }
            /// put this link to the next cells.
        }
    }

    void IdentifyGPSODPoints(int agent_no)
    {
        // first step, determine the inside flag

        // default settings if all GPS points inside
        g_agent_vector[agent_no].head_gps_index = -1;
        g_agent_vector[agent_no].tail_gps_index = g_agent_vector[agent_no].m_GPSPointVector.size() - 1;

        int g_Inside_index = -1;
        for (int g = 0; g < g_agent_vector[agent_no].m_GPSPointVector.size(); g++) // for each GPS point
        {                                                                      // x_key and y_key are relative index of grid

            int x_key = g_GetCellXID(g_agent_vector[agent_no].m_GPSPointVector[g].pt.x, m_left);
            int y_key = g_GetCellYID(g_agent_vector[agent_no].m_GPSPointVector[g].pt.y, m_bottom);

            __int64 cell_id = g_GetCellID(g_agent_vector[agent_no].m_GPSPointVector[g].pt.x, g_agent_vector[agent_no].m_GPSPointVector[g].pt.y);

            fprintf(g_pFileLog, "trace index %d, trace no: %d, cell %d, %d, %jd \n", g, g_agent_vector[agent_no].m_GPSPointVector[g].trace_no, x_key, y_key);

            m_GridMatrix[x_key][y_key].m_GPSPointVector.push_back(g_agent_vector[agent_no].m_GPSPointVector[g]);
            m_GridMatrix[x_key][y_key].cell_id = cell_id;


            if (g_agent_vector[agent_no].m_GPSPointVector[g].global_time_in_second / 60 > m_GridMatrix[x_key][y_key].end_min)
                m_GridMatrix[x_key][y_key].end_min = g_agent_vector[agent_no].m_GPSPointVector[g].global_time_in_second / 60;

            if (g_agent_vector[agent_no].m_GPSPointVector[g].global_time_in_second / 60 < m_GridMatrix[x_key][y_key].start_min)
                m_GridMatrix[x_key][y_key].start_min = g_agent_vector[agent_no].m_GPSPointVector[g].global_time_in_second / 60;

            if (g_cell_id_2_node_map.find(cell_id) == g_cell_id_2_node_map.end())  // the subarea network not defined yet
                continue;

            if (g_agent_vector[agent_no].head_gps_index == -1)
                g_agent_vector[agent_no].head_gps_index = g;

            g_agent_vector[agent_no].tail_gps_index = g;

        }

        if (g_agent_vector[agent_no].m_GPSPointVector.size() == 0)
            return;


        // second step for origin GPS index
        int  g = g_agent_vector[agent_no].head_gps_index;


        g_agent_vector[agent_no].m_o_boundary_point[0] = g_agent_vector[agent_no].m_GPSPointVector[g].pt;  // setup the default first and second origin GPS points from the head GPS point
        g_agent_vector[agent_no].m_o_boundary_point[1] = g_agent_vector[agent_no].m_GPSPointVector[g].pt;

        if (g + 1 < g_agent_vector[agent_no].m_GPSPointVector.size())
        {
            g_agent_vector[agent_no].m_o_boundary_point[1] = g_agent_vector[agent_no].m_GPSPointVector[g + 1].pt;  // second to last GPS point
        }

        // third step for destination GPS index

        g = g_agent_vector[agent_no].tail_gps_index;
        g_agent_vector[agent_no].m_d_boundary_point[0] = g_agent_vector[agent_no].m_GPSPointVector[g].pt; // setup the default last and second-to-last destination GPS points from the head GPS point
        g_agent_vector[agent_no].m_d_boundary_point[1] = g_agent_vector[agent_no].m_GPSPointVector[g].pt;

        if (g - 1 >= 0)
        {
            g_agent_vector[agent_no].m_d_boundary_point[0] = g_agent_vector[agent_no].m_GPSPointVector[g - 1].pt;  // second to last GPS point
        }


    }

    void IdentifyNetworkONode(int agent_no)
    {
        double min_distance_to_boundary_point = MAX_LABEL_COST_;

        if (g_agent_vector[agent_no].o_node_no >= 0)
        {
            origin_node_no = g_agent_vector[agent_no].o_node_no;
            return;
        }


        for (int l = 0; l < g_mm_link_vector.size(); l++) // for all links in this cell
        {

            string s1 = g_agent_vector[agent_no].allowed_link_type_code;
            string s2 = g_mm_link_vector[l].link_type_code;

            if (s1.size() > 0 && s2.size() > 0) // both agent and link with link type code
            {
                if (s1.find(s2) == std::string::npos)
                {
                    continue;  // s2 is not contained in s1, skip this type of links in identifying origins
                }
            }

            s1 = g_agent_vector[agent_no].blocked_link_type_code;
            if (s1.size() > 0 && s2.size() > 0) // both agent and link with link type code
            {
                if (s1.find(s2) != std::string::npos) // s2 is contained in s1
                {
                    continue;  // skip this type of links in identifying origins
                }
            }
            double distance = MAX_LABEL_COST_;

            //TRACE("%d->%d\n", g_mm_link_vector[l].from_node_id, g_mm_link_vector[l].to_node_id);

            //if (g_mm_link_vector[l].from_node_id == 402 && g_mm_link_vector[l].to_node_id == 131)
            //    TRACE("%d->%d\n", g_mm_link_vector[l].from_node_id, g_mm_link_vector[l].to_node_id);

            //if (g_mm_link_vector[l].from_node_id == 433 && g_mm_link_vector[l].to_node_id == 341)
            //    TRACE("%d->%d\n", g_mm_link_vector[l].from_node_id, g_mm_link_vector[l].to_node_id);

            //ill conditioning detection

            double p2l_distance = 999;
            if (g_ill_conditioning_detection(g_mm_link_vector[l].link_distance, g_agent_vector[agent_no].first_segment_distance) == false)
            { // case of good conditioning
                double distance_from = g_GetPoint2LineDistance(&g_agent_vector[agent_no].m_o_boundary_point[0], &g_mm_node_vector[g_mm_link_vector[l].from_node_seq_no].pt, &g_mm_node_vector[g_mm_link_vector[l].to_node_seq_no].pt,
                    1, false);

                double distance_to = 0;

                if (g_agent_vector[agent_no].m_GPSPointVector.size() >= 2)
                {
                    distance_to = g_GetPoint2LineDistance(&g_agent_vector[agent_no].m_o_boundary_point[1],
                        &g_mm_node_vector[g_mm_link_vector[l].from_node_seq_no].pt, &g_mm_node_vector[g_mm_link_vector[l].to_node_seq_no].pt,
                        1, false);
                }

                p2l_distance = (distance_from + distance_to) / 2;
                double distance_from_p2p = 0;
                double distance_to_p2p = 0;

                distance_from_p2p = g_GetPoint2Point_Distance(&g_agent_vector[agent_no].m_o_boundary_point[0], &g_mm_node_vector[g_mm_link_vector[l].from_node_seq_no].pt);
                distance_to_p2p = g_GetPoint2Point_Distance(&g_agent_vector[agent_no].m_o_boundary_point[1], &g_mm_node_vector[g_mm_link_vector[l].to_node_seq_no].pt);
                distance = (distance_from + distance_to + distance_from_p2p + distance_to_p2p) / 4;
            }
            else
            {
                //case of ill conditioning, we have no help, we can only take the minimum of point to point distance
                double distance_from_p2p = 0;
                double distance_to_p2p = 0;

                distance_from_p2p = g_GetPoint2Point_Distance(&g_agent_vector[agent_no].m_o_boundary_point[0], &g_mm_node_vector[g_mm_link_vector[l].from_node_seq_no].pt);
                distance_to_p2p = g_GetPoint2Point_Distance(&g_agent_vector[agent_no].m_o_boundary_point[1], &g_mm_node_vector[g_mm_link_vector[l].to_node_seq_no].pt);

                distance = min(distance_from_p2p, distance_to_p2p);

                // consider the minimal distance of any point, and avg distance of cross-section distance
            }

            // we check this relative angle condition for both ill and good conditions,
            double relative_angle = fabs(g_Find_PPP_RelativeAngle(
                &g_agent_vector[agent_no].m_o_boundary_point[0],
                &g_agent_vector[agent_no].m_o_boundary_point[1],
                &g_mm_node_vector[g_mm_link_vector[l].from_node_seq_no].pt,
                &g_mm_node_vector[g_mm_link_vector[l].to_node_seq_no].pt));

            if (relative_angle > 45)
            {
                // add penalty for opposite direction
                distance = distance * 10; /// 10 times as panalty
            }

            int i_trace = 0;

            if (distance < g_mm_link_vector[l].o_distance)
                g_mm_link_vector[l].o_distance = distance;

            if (distance < min_distance_to_boundary_point)
            {

                min_distance_to_boundary_point = distance;
                g_mm_link_vector[l].o_distance = distance;
                origin_node_no = g_mm_link_vector[l].from_node_seq_no;
                fprintf(g_pFileLog, "finding origin_node: %d -> %d, %f \n", g_mm_link_vector[l].from_node_id, g_mm_link_vector[l].to_node_id, distance);
                g_agent_vector[agent_no].matching_link_no = l;
                //							g_agent_vector[agent_no].upstream_node_matched_time_in_sec = t;

            }
        }


    }

    void IdentifyNetworkDNode(int agent_no)
    {
        if (g_agent_vector[agent_no].d_node_no >= 0)
        {
            destination_node_no = g_agent_vector[agent_no].d_node_no;
            return;
        }

        double min_distance_to_boundary_point = MAX_LABEL_COST_;
        for (int l = 0; l < g_mm_link_vector.size(); l++) // for all links in this cell
        {
            //TRACE("%d->%d\n", g_mm_link_vector[l].from_node_id, g_mm_link_vector[l].to_node_id);

            string s1 = g_agent_vector[agent_no].allowed_link_type_code;
            string s2 = g_mm_link_vector[l].link_type_code;

            if (s1.size() > 0 && s2.size() > 0)  // both agent and link with link type code
            {
                if (s1.find(s2) == std::string::npos)
                {
                    continue;  // s2 is not contained in s1, skip this type of links in identifying origins
                }
            }

            s1 = g_agent_vector[agent_no].blocked_link_type_code;
            if (s1.size() > 0 && s2.size() > 0) // both agent and link with link type code
            {
                if (s1.find(s2) != std::string::npos) // s2 is contained in s1
                {
                    continue;  // skip this type of links in identifying origins
                }
            }

            double distance = MAX_LABEL_COST_;

            if (g_ill_conditioning_detection(g_mm_link_vector[l].link_distance, g_agent_vector[agent_no].last_segment_distance) == false)
            { // case of good conditioning
                double distance_from = g_GetPoint2LineDistance(&g_agent_vector[agent_no].m_d_boundary_point[0],
                    &g_mm_node_vector[g_mm_link_vector[l].from_node_seq_no].pt, &g_mm_node_vector[g_mm_link_vector[l].to_node_seq_no].pt,
                    1, false);
                double distance_to = 0;

                distance_to = g_GetPoint2LineDistance(&g_agent_vector[agent_no].m_d_boundary_point[1],
                    &g_mm_node_vector[g_mm_link_vector[l].from_node_seq_no].pt, &g_mm_node_vector[g_mm_link_vector[l].to_node_seq_no].pt,
                    1, false);

                //ill conditioning detection

                double distance_from_p2p = 0;
                double distance_to_p2p = 0;

                distance_from_p2p = g_GetPoint2Point_Distance(&g_agent_vector[agent_no].m_d_boundary_point[0], &g_mm_node_vector[g_mm_link_vector[l].from_node_seq_no].pt);
                distance_to_p2p = g_GetPoint2Point_Distance(&g_agent_vector[agent_no].m_d_boundary_point[1], &g_mm_node_vector[g_mm_link_vector[l].to_node_seq_no].pt);
                distance = (distance_from + distance_to + distance_from_p2p + distance_to_p2p) / 4;
            }
            else
            { // case of ill conditionning
                double distance_from_p2p = 0;
                double distance_to_p2p = 0;

                distance_from_p2p = g_GetPoint2Point_Distance(&g_agent_vector[agent_no].m_d_boundary_point[0], &g_mm_node_vector[g_mm_link_vector[l].from_node_seq_no].pt);
                //[1] ending point of GPS point segment
                distance_to_p2p = g_GetPoint2Point_Distance(&g_agent_vector[agent_no].m_d_boundary_point[1], &g_mm_node_vector[g_mm_link_vector[l].to_node_seq_no].pt);
                distance = min(distance_from_p2p, distance_to_p2p);
            }

            // we check this relative angle condition for both ill and good conditions,
            double relative_angle = fabs(g_Find_PPP_RelativeAngle(
                &g_agent_vector[agent_no].m_d_boundary_point[0],
                &g_agent_vector[agent_no].m_d_boundary_point[1],
                &g_mm_node_vector[g_mm_link_vector[l].from_node_seq_no].pt,
                &g_mm_node_vector[g_mm_link_vector[l].to_node_seq_no].pt));

            if (relative_angle > 45)
            {
                // add penalty for opposite direction
                distance = distance * 10; /// 10 times as panalty
            }

            if (distance < g_mm_link_vector[l].d_distance)
                g_mm_link_vector[l].d_distance = distance;

            if (distance < min_distance_to_boundary_point)
            {
                min_distance_to_boundary_point = distance;
                destination_node_no = g_mm_link_vector[l].to_node_seq_no;
                fprintf(g_pFileLog, "finding destination_node: %d -> %d, %f \n", g_mm_link_vector[l].from_node_id, g_mm_link_vector[l].to_node_id, distance);
                //TRACE("%d -> %d, %f \n", g_mm_link_vector[l].from_node_id, g_mm_link_vector[l].to_node_id, distance);
            }
        }


    }

    void output_TD_label_cost()
    {
        for (int l = 0; l < g_mm_link_vector.size(); l++) // 
        {
            if (g_mm_link_vector[l].likelihood_distance < 998)
            {
                for (int t = 0; t < g_TimeRangeInterval; t++)
                {
                    if (m_TD_link_generalised_cost_array[l][t] < 998)
                    {

                        float timestamp_in_min = t * g_TimeResolution_inMin + g_StartTimeinMin;
                        int t_in_sec = (timestamp_in_min - g_StartTimeinMin) * 60;
                        fprintf(g_pFileLog, "arc cost %d -> %d at time index %d and sec %d = %f, hit flag: %d: GPS point ID: %d \n", g_mm_link_vector[l].from_node_id, g_mm_link_vector[l].to_node_id,
                            t, t_in_sec, m_TD_link_generalised_cost_array[l][t], m_TD_link_GPS_hit_array[l][t], m_TD_link_GPS_point_array[l][t]
                        );

                    }
                }

            }

        }
    }


    bool AddGPSPointsIntoGridSystem(int agent_no)
    { // for every agent

        for (int i = 0; i < g_mm_link_vector.size(); i++) // reset the cost for all links
        {
            g_mm_link_vector[i].likely_trace_no = -1;
            m_link_matching_trace_no_array[i] = -1;
        }

        for (int x_i = 0; x_i <= g_grid_size; x_i++)
            for (int y_i = 0; y_i <= g_grid_size; y_i++)
            {
                m_GridMatrix[x_i][y_i].possible_dwell_cell_flag = 0;
                m_GridMatrix[x_i][y_i].dwell_time_min = 0;
                m_GridMatrix[x_i][y_i].end_min = 0;
                m_GridMatrix[x_i][y_i].start_min = 9999999;
                m_GridMatrix[x_i][y_i].m_GPSPointVector.clear(); //reset the existing GPS point records in this grid
            }

        for (int i = 0; i < g_mm_link_vector.size(); i++) // reset the cost for all links
        {
            m_link_generalised_cost_array[i] = MAX_LABEL_COST_ / 1000; // feasible range

            g_mm_link_vector[i].bInsideFlag = false;
            g_mm_link_vector[i].hit_count = 0;
            for (int t = 0; t < g_TimeRangeInterval; t++)
            {
                m_TD_link_generalised_cost_array[i][t] = MAX_LABEL_COST_ / 1000; // feasible range
                m_TD_link_GPS_point_array[i][t] = -1;
                m_TD_link_GPS_hit_array[i][t] = -1;
            }

        }

        // put GPS points into grid cell
        IdentifyGPSODPoints(agent_no);

        // calculate avg distance
        double total_GPS_distance = 0;
        int g;

        if (g_agent_vector[agent_no].m_GPSPointVector.size() >= 2)
        {
            std::vector<int> time_interval_vector;


            // find median 
            for (g = 1; g < g_agent_vector[agent_no].m_GPSPointVector.size(); g++) // for each GPS point
            {
                g_agent_vector[agent_no].m_GPSPointVector[g].interval_in_second = g_agent_vector[agent_no].m_GPSPointVector[g].global_time_in_second - g_agent_vector[agent_no].m_GPSPointVector[g - 1].global_time_in_second;
                g_agent_vector[agent_no].m_GPSPointVector[g].relative_time_in_second =
                    g_agent_vector[agent_no].m_GPSPointVector[g].global_time_in_second - (int)(g_agent_vector[agent_no].start_time_in_min * 60 + 0.5);
                g_agent_vector[agent_no].m_GPSPointVector[g].distance = g_GetPoint2Point_Distance(&g_agent_vector[agent_no].m_GPSPointVector[g].pt,
                    &g_agent_vector[agent_no].m_GPSPointVector[g - 1].pt);


                time_interval_vector.push_back(g_agent_vector[agent_no].m_GPSPointVector[g].interval_in_second);
            }

            if (time_interval_vector.size() >= 4)
                g_agent_vector[agent_no].sampling_rate_in_min = g_findMedian(time_interval_vector, time_interval_vector.size()) / 60.0;
            else
                g_agent_vector[agent_no].sampling_rate_in_min = g_SampleTimeResolution_inMin;


            // find excessive dwell time 

            // find median 
            float ratio_excessive = 5;
            for (g = 1; g < g_agent_vector[agent_no].m_GPSPointVector.size(); g++) // for each GPS point
            {

                if (g_agent_vector[agent_no].m_GPSPointVector[g].interval_in_second > g_agent_vector[agent_no].sampling_rate_in_min * 60 * ratio_excessive)
                {
                    float time_stamp_in_min = g_agent_vector[agent_no].m_GPSPointVector[g - 1].global_time_in_second / 60.0;  // previous GPS point 
                    int timestamp_in_interval = g_GetTimeInterval(time_stamp_in_min);
                    int dwell_time_interval = g_agent_vector[agent_no].m_GPSPointVector[g].interval_in_second / 60.0 / g_TimeResolution_inMin;

                    g_agent_vector[agent_no].m_ExcessiveDwellMap[timestamp_in_interval] = dwell_time_interval;
                }
            }


            for (g = 0; g < g_agent_vector[agent_no].m_GPSPointVector.size() - 1; g++) // for each GPS point
            {
                double segment_distance;
                segment_distance = g_GetPoint2Point_Distance(&g_agent_vector[agent_no].m_GPSPointVector[g].pt, &g_agent_vector[agent_no].m_GPSPointVector[g + 1].pt);

                total_GPS_distance += segment_distance;
                if (g == 0)
                {
                    g_agent_vector[agent_no].first_segment_distance = segment_distance;
                }

                if (g == g_agent_vector[agent_no].m_GPSPointVector.size() - 2)
                {
                    g_agent_vector[agent_no].last_segment_distance = segment_distance;
                }
            }
        }

        g_agent_vector[agent_no].avg_GPS_segment_distance = total_GPS_distance / max(1, (int)(g_agent_vector[agent_no].m_GPSPointVector.size() - 1));

        IdentifyNetworkONode(agent_no);
        IdentifyNetworkDNode(agent_no);


        // fourth step
        // for each grid matrix cell
        //scan x and y index in the grid
        // m_link_generalised_cost_array is the link cost used in shortest path

        for (int x_i = 0; x_i <= g_grid_size; x_i++)
            for (int y_i = 0; y_i <= g_grid_size; y_i++)
            {

                if (m_GridMatrix[x_i][y_i].m_GPSPointVector.size() > 0)
                {

                    if (m_GridMatrix[x_i][y_i].end_min - m_GridMatrix[x_i][y_i].start_min > 120) // 2 hours
                    {
                        m_GridMatrix[x_i][y_i].possible_dwell_cell_flag = 1;
                        m_GridMatrix[x_i][y_i].dwell_time_min = m_GridMatrix[x_i][y_i].end_min - m_GridMatrix[x_i][y_i].start_min;
                    }
                }
            }


        for (int x_i = 0; x_i <= g_grid_size; x_i++)
            for (int y_i = 0; y_i <= g_grid_size; y_i++)
            {

                if (m_GridMatrix[x_i][y_i].m_GPSPointVector.size() > 0)
                {
                    // compute the average distance from the GPS points (g, g+1) to the ending points of a link
                    for (int local_l = 0; local_l < m_GridMatrix[x_i][y_i].m_LinkNoVector.size(); local_l++) // for all links in this cell
                    {
                        int l = m_GridMatrix[x_i][y_i].m_LinkNoVector[local_l];
                        m_link_generalised_cost_array[l] = g_GridResolution; // feasible range
                        g_mm_link_vector[l].bInsideFlag = true;

                        if (m_GridMatrix[x_i][y_i].possible_dwell_cell_flag == 1)
                        {
                            g_mm_link_vector[l].Possible_Dwell_time_in_min = m_GridMatrix[x_i][y_i].dwell_time_min;
                            g_mm_link_vector[l].dwell_start_time_in_min = m_GridMatrix[x_i][y_i].start_min;
                        }
                        g_mm_node_vector[g_mm_link_vector[l].from_node_seq_no].bInsideFlag = true;
                        g_mm_node_vector[g_mm_link_vector[l].to_node_seq_no].bInsideFlag = true;

                    }


                }

                // for each grid cell
                // second, we now select the mininum of GPS point (in the same cell) to link distance to set the link cost
                for (int g = 0; g < m_GridMatrix[x_i][y_i].m_GPSPointVector.size(); g++) // for each GPS point of an agent in the cell
                {

                    // compute the average distance from the GPS points (g, g+1) to the ending points of a link
                    for (int local_l = 0; local_l < m_GridMatrix[x_i][y_i].m_LinkNoVector.size(); local_l++) // for all links in this cell
                    {
                        int l = m_GridMatrix[x_i][y_i].m_LinkNoVector[local_l];


                        //if (g_mm_link_vector[l].from_node_id == 207 && g_mm_link_vector[l].to_node_id == 208)
                        //    TRACE("%d->%d\n", g_mm_link_vector[l].from_node_id, g_mm_link_vector[l].to_node_id);


                        g_mm_link_vector[l].o_distance = MAX_LABEL_COST_;
                        g_mm_link_vector[l].d_distance = MAX_LABEL_COST_;

                        //if (m_GridMatrix[x_i][y_i].m_GPSPointVector[g].trace_no == 4)
                        //{
                        //    //if (g_mm_link_vector[l].from_node_id == 10386)
                        //    //{
                        //    //    TRACE("");
                        //    //}
                        //}

                        //if (g_mm_link_vector[l].from_node_id == 201 && g_mm_link_vector[l].to_node_id == 202 && m_GridMatrix[x_i][y_i].m_GPSPointVector[g].trace_no == 11)
                        //{
                        //    TRACE("GPS point %d", g);
                        //}



                        bool bHitCount = false;


                        double p2l_distance_hit = g_GetPoint2LineDistance(&m_GridMatrix[x_i][y_i].m_GPSPointVector[g].pt, &g_mm_node_vector[g_mm_link_vector[l].from_node_seq_no].pt, &g_mm_node_vector[g_mm_link_vector[l].to_node_seq_no].pt,
                            1, false);

                        for (int ls = 0; ls < g_mm_link_vector[l].m_PointVector.size() - 1; ls++)  //for a pair of shape points along the link
                        {

                            if (m_GridMatrix[x_i][y_i].m_GPSPointVector.size() >= 2 && g <= m_GridMatrix[x_i][y_i].m_GPSPointVector.size() - 2
                                && g_GetTwoPoints2LineIntersectionFlag(&m_GridMatrix[x_i][y_i].m_GPSPointVector[g].pt, &m_GridMatrix[x_i][y_i].m_GPSPointVector[g + 1].pt,
                                    &g_mm_link_vector[l].m_PointVector[ls], &g_mm_link_vector[l].m_PointVector[ls + 1]))
                            {
                                g_mm_link_vector[l].hit_count += 1;
                                bHitCount = true;
                                break;
                            }
                        }

                        double p2l_distance = g_GetPoint2LineDistance(&m_GridMatrix[x_i][y_i].m_GPSPointVector[g].pt, &g_mm_node_vector[g_mm_link_vector[l].from_node_seq_no].pt, &g_mm_node_vector[g_mm_link_vector[l].to_node_seq_no].pt,
                            1, true);  // recacluate the point 2 line distance, without intersection requirement 

                        // we consider GPS segment to the link shape point segment distance
                        for (int p = 0; p < g_mm_link_vector[l].m_PointVector.size(); p++)
                        {

                            double distance_from = g_GetPoint2Point_Distance(&m_GridMatrix[x_i][y_i].m_GPSPointVector[g].pt, &g_mm_link_vector[l].m_PointVector[p]);
                            double distance_to = 0;


                            if (g_agent_vector[agent_no].m_GPSPointVector.size() >= 2 && p < g_mm_link_vector[l].m_PointVector.size() - 1 &&
                                (g != m_GridMatrix[x_i][y_i].m_GPSPointVector.size() - 1)) // boundary points)
                            {
                                distance_to = g_GetPoint2Point_Distance(&m_GridMatrix[x_i][y_i].m_GPSPointVector[g + 1].pt, &g_mm_link_vector[l].m_PointVector[p + 1]);
                            }

                            // we do not need to detect ill conditionning here, as the link in the cell , and the GPS trace information is all we have.
                            // distance from is the distance from GPS point g to from node of link
                            // distance to is the distance from GPS point g to to node of link
                            double nonhit_distance = (p2l_distance * 10 + distance_from + distance_to) / 10;


                            if (nonhit_distance > g_mm_link_vector[l].link_distance * 2)  // filter out long deviation
                                continue;

                            if (nonhit_distance < m_link_generalised_cost_array[l])  // for static case
                            {
                                m_link_generalised_cost_array[l] = nonhit_distance; // use this distance as the likelihood cost
                                m_link_matching_trace_no_array[l] = m_GridMatrix[x_i][y_i].m_GPSPointVector[g].trace_no;

                            }

                            int time_interval = g_GetTimeInterval(m_GridMatrix[x_i][y_i].m_GPSPointVector[g].global_time_in_second / 60.0);

                            double distance;
                            if (bHitCount)
                                distance = p2l_distance_hit;
                            else
                                distance = nonhit_distance * g_NonHitDistanceRatio;

                            if (distance < m_TD_link_generalised_cost_array[l][time_interval])  // TD case
                            {
                                //if (g_mm_link_vector[l].from_node_id == 201 && g_mm_link_vector[l].to_node_id == 202)
                                //{
                                //    TRACE("GPS point %d", g);
                                //}
                                int offset_in_interval = max(1.0f, g_agent_vector[agent_no].sampling_rate_in_min / g_TimeResolution_inMin);
                                int t0 = max(0, time_interval - offset_in_interval);
                                int t1 = min(time_interval + offset_in_interval, g_TimeRangeInterval);

                                for (int t = t0; t < t1; t++)
                                {
                                    m_TD_link_generalised_cost_array[l][t] = distance;
                                    m_TD_link_GPS_point_array[l][t] = m_GridMatrix[x_i][y_i].m_GPSPointVector[g].trace_no;
                                    m_TD_link_GPS_hit_array[l][t] = bHitCount;
                                }

                            }

                        }


                    }
                }



                // stage 2 handling for no matching trace_no for the link


                for (int local_l = 0; local_l < m_GridMatrix[x_i][y_i].m_LinkNoVector.size(); local_l++) // for all links in this cell
                {
                    int l = m_GridMatrix[x_i][y_i].m_LinkNoVector[local_l];

                    if (m_link_matching_trace_no_array[l] < 0)  // no matching GPS points
                    {
                        double min_distance = 99999999;
                        int trace_no = -1;
                        for (int g = 0; g < m_GridMatrix[x_i][y_i].m_GPSPointVector.size(); g++) // for each GPS point of an agent in the cell
                        {
                            double p2l_distance = g_GetPoint2LineDistance(&m_GridMatrix[x_i][y_i].m_GPSPointVector[g].pt, &g_mm_node_vector[g_mm_link_vector[l].from_node_seq_no].pt, &g_mm_node_vector[g_mm_link_vector[l].to_node_seq_no].pt,
                                1, false);
                            if (p2l_distance < min_distance)
                            {
                                min_distance = p2l_distance;
                                trace_no = m_GridMatrix[x_i][y_i].m_GPSPointVector[g].trace_no;
                            }

                        }

                        m_link_matching_trace_no_array[l] = trace_no;

                    }

                }
            }

        for (int i = 0; i < g_mm_link_vector.size(); i++) // reset the cost for all links
        {
            g_mm_link_vector[i].likelihood_distance = m_link_generalised_cost_array[i];

            g_mm_link_vector[i].likely_trace_no = m_link_matching_trace_no_array[i];

        }



        output_TD_label_cost();
        return true;
    }

    void output_grid_file()
    {

        FILE* g_pFileGrid = nullptr;
        g_pFileGrid = fopen("grid.csv", "w");

        if (g_pFileGrid == NULL)
        {
            cout << "File grid.csv cannot be opened." << endl;
            g_Program_stop();
        }
        else
        {
            fprintf(g_pFileGrid, "cell_id,x_index,y_index,zone_id,origin_flag,destination_flag,dwell_flag,n_GPS_points,n_nodes,x_coord,y_coord,geometry,\n"); //hhmmss,trace_id,travel_time,delay,geometry
            for (int x_i = 0; x_i <= g_grid_size; x_i++)
                for (int y_i = 0; y_i <= g_grid_size; y_i++)
                {

                    if (m_GridMatrix[x_i][y_i].m_NodeVector.size() > 0 || m_GridMatrix[x_i][y_i].m_GPSPointVector.size())
                    {
                        int left_int = floor(m_left / g_GridResolution);
                        int bottom_int = floor(m_bottom / g_GridResolution);

                        double x_coord_left = (left_int * g_GridResolution) + x_i * g_GridResolution;
                        double y_coord_bottom = (bottom_int * g_GridResolution) + y_i * g_GridResolution;
                        double x_coord_right = x_coord_left + g_GridResolution;
                        double y_coord_top = y_coord_bottom + g_GridResolution;
                        int origin_cell_flag = m_GridMatrix[x_i][y_i].origin_cell_flag;
                        int destination_cell_flag = m_GridMatrix[x_i][y_i].destination_cell_flag;
                        int possible_dwell_flag = m_GridMatrix[x_i][y_i].possible_dwell_cell_flag;

                        fprintf(g_pFileGrid, "%jd,%d,%d,%d,%d,%d,%d,%d,%d,", m_GridMatrix[x_i][y_i].cell_id, x_i, y_i, m_GridMatrix[x_i][y_i].zone_id, origin_cell_flag, destination_cell_flag, possible_dwell_flag,
                            m_GridMatrix[x_i][y_i].m_GPSPointVector.size(),
                            m_GridMatrix[x_i][y_i].m_NodeVector.size()
                        ); //hhmmss,trace_id,travel_time,delay,geometry
                        fprintf(g_pFileGrid, "%f,%f,", x_coord_left, y_coord_top);
                        fprintf(g_pFileGrid, "\"LINESTRING (");

                        fprintf(g_pFileGrid, "%f %f,", x_coord_left, y_coord_top);
                        fprintf(g_pFileGrid, "%f %f,", x_coord_right, y_coord_top);
                        fprintf(g_pFileGrid, "%f %f,", x_coord_right, y_coord_bottom);
                        fprintf(g_pFileGrid, "%f %f,", x_coord_left, y_coord_bottom);
                        fprintf(g_pFileGrid, "%f %f,", x_coord_left, y_coord_top);
                        fprintf(g_pFileGrid, ")\"");
                        fprintf(g_pFileGrid, "\n");
                    }

                }


            fclose(g_pFileGrid);
        }
    }

    int UpdateLinkLRPriceGridSystem(int agent_no)
    { // for every agent

        int balance_count = 0;
        float benefit = 100;
        for (int l = 0; l < g_mm_link_vector.size(); l++)
        {

            if (g_mm_link_vector[l].hit_count >= 1 && g_mm_link_vector[l].use_count == 0)// subgradient algorithm for adjusting price 
            {
                g_mm_link_vector[l].balance = g_mm_link_vector[l].use_count - g_mm_link_vector[l].hit_count;

                for (int t = 0; t < g_TimeRangeInterval; t++)
                {
                    if (m_TD_link_GPS_hit_array[l][t] >= 1)  // TD case
                        m_TD_link_generalised_cost_array[l][t] += benefit * g_mm_link_vector[l].balance;
                }

                balance_count += 1;
            }

        }

        fprintf(g_pFileLog, "LR balance count = %d \n", balance_count);
        cout << "balance_count = " << balance_count << endl;
        return balance_count;



    }




    std::vector<int> m_agent_vector;
    int m_memory_block_no;

    std::vector<int> m_origin_node_vector; // assigned nodes for computing
    std::vector<int> m_origin_zone_seq_no_vector;

    int tau;             // assigned nodes for computing
    int m_mode_type_no; // assigned nodes for computing
    double m_value_of_time;

    int m_threadNo; // internal thread number

    int m_ListFront; // used in coding SEL
    int m_ListTail;  // used in coding SEL

    int* m_SENodeList; // used in coding SEL

    double* m_node_label_cost;      // label cost // for shortest path calcuating
    double* m_label_time_array;     // time-based cost
    double* m_label_distance_array; // distance-based cost

    int* m_node_predecessor;  // predecessor for nodes
    int* m_node_status_array; // update status
    int* m_link_predecessor;  // predecessor for this node points to the previous link that updates its label cost (as part of optimality condition) (for easy referencing)

    double* m_link_flow_volume_array;

    double* m_link_generalised_cost_array;
    int* m_link_matching_trace_no_array;

    double** m_TD_link_generalised_cost_array;
    int** m_TD_link_GPS_point_array;
    int** m_TD_link_GPS_hit_array;

    double** m_TD_node_label_cost;
    int** m_TD_node_predecessor;  // predecessor for nodes
    int** m_TD_link_predecessor;  // predecessor for this node points to the previous link that updates its label cost (as part of optimality condition) (for easy referencing)
    int** m_TD_time_predecessor;  // predecessor for this node points to the previous link that updates its label cost (as part of optimality condition) (for easy referencing)

    int* temp_path_node_vector;
    // major function 1:  allocate memory and initialize the data
    void AllocateMemory(int number_of_nodes, int number_of_links)
    {

        m_SENodeList = new int[number_of_nodes]; //1

        m_node_status_array = new int[number_of_nodes];       //2
        m_label_time_array = new double[number_of_nodes];     //3
        m_label_distance_array = new double[number_of_nodes]; //4
        m_node_predecessor = new int[number_of_nodes];        //5
        m_link_predecessor = new int[number_of_nodes];        //6
        m_node_label_cost = new double[number_of_nodes];      //7

        m_link_flow_volume_array = new double[number_of_links]; //8

        m_link_generalised_cost_array = new double[number_of_links]; //9

        m_link_matching_trace_no_array = new int[number_of_links]; //9
        temp_path_node_vector = new int[number_of_nodes];

        m_TD_link_generalised_cost_array = Allocate2DDynamicArray<double>(number_of_links, g_TimeRangeInterval);

        m_TD_link_GPS_point_array = Allocate2DDynamicArray<int>(number_of_links, g_TimeRangeInterval);

        m_TD_link_GPS_hit_array = Allocate2DDynamicArray<int>(number_of_links, g_TimeRangeInterval);
        m_TD_node_label_cost = Allocate2DDynamicArray<double>(number_of_nodes, g_TimeRangeInterval);
        m_TD_node_predecessor = Allocate2DDynamicArray<int>(number_of_nodes, g_TimeRangeInterval);
        m_TD_link_predecessor = Allocate2DDynamicArray<int>(number_of_nodes, g_TimeRangeInterval);
        m_TD_time_predecessor = Allocate2DDynamicArray<int>(number_of_nodes, g_TimeRangeInterval);

    }

    ~NetworkForSP()
    {

        if (m_SENodeList != NULL) //1
            delete[] m_SENodeList;

        if (m_node_status_array != NULL) //2
            delete[] m_node_status_array;

        if (m_label_time_array != NULL) //3
            delete[] m_label_time_array;

        if (m_label_distance_array != NULL) //4
            delete[] m_label_distance_array;

        if (m_node_predecessor != NULL) //5
            delete[] m_node_predecessor;

        if (m_link_predecessor != NULL) //6
            delete[] m_link_predecessor;

        if (m_node_label_cost != NULL) //7
            delete[] m_node_label_cost;

        if (m_link_flow_volume_array != NULL) //8
            delete[] m_link_flow_volume_array;

        if (m_link_generalised_cost_array != NULL) //9
            delete[] m_link_generalised_cost_array;

        if (m_link_matching_trace_no_array != NULL) //9
            delete[] m_link_matching_trace_no_array;

        if (temp_path_node_vector != NULL) //9
            delete[] temp_path_node_vector;

        if (m_GridMatrix)
            Deallocate2DDynamicArray<GridNodeSet>(m_GridMatrix, MAX_GRID_SIZE_);

        if (m_TD_link_generalised_cost_array)
            Deallocate2DDynamicArray<double>(m_TD_link_generalised_cost_array, g_TimeRangeInterval);

        if (m_TD_link_GPS_point_array)
            Deallocate2DDynamicArray<int>(m_TD_link_GPS_point_array, g_TimeRangeInterval);

        if (m_TD_link_GPS_hit_array)
            Deallocate2DDynamicArray<int>(m_TD_link_GPS_hit_array, g_TimeRangeInterval);


        if (m_TD_node_label_cost)
            Deallocate2DDynamicArray<double>(m_TD_node_label_cost, g_TimeRangeInterval);

        if (m_TD_node_predecessor)
            Deallocate2DDynamicArray<int>(m_TD_node_predecessor, g_TimeRangeInterval);

        if (m_TD_link_predecessor)
            Deallocate2DDynamicArray<int>(m_TD_link_predecessor, g_TimeRangeInterval);

        if (m_TD_time_predecessor)
            Deallocate2DDynamicArray<int>(m_TD_time_predecessor, g_TimeRangeInterval);
    }

    // SEList: scan eligible List implementation: the reason for not using STL-like template is to avoid overhead associated pointer allocation/deallocation
    inline void SEList_clear()
    {
        m_ListFront = -1;
        m_ListTail = -1;
    }

    inline void SEList_push_front(int node)
    {
        if (m_ListFront == -1) // start from empty
        {
            m_SENodeList[node] = -1;
            m_ListFront = node;
            m_ListTail = node;
        }
        else
        {
            m_SENodeList[node] = m_ListFront;
            m_ListFront = node;
        }
    }
    inline void SEList_push_back(int node)
    {
        if (m_ListFront == -1) // start from empty
        {
            m_ListFront = node;
            m_ListTail = node;
            m_SENodeList[node] = -1;
        }
        else
        {
            m_SENodeList[m_ListTail] = node;
            m_SENodeList[node] = -1;
            m_ListTail = node;
        }
    }

    inline bool SEList_empty()
    {
        return (m_ListFront == -1);
    }

    //major function: update the cost for each node at each SP tree, using a stack from the origin structure

    int origin_node_no;
    int destination_node_no;

    //major function 2: // abel correcting algorithm with double queue implementation
    double optimal_label_correcting(int agent_no, string allowed_link_type_code, string blocked_link_type_code)
    {
        origin_node_no = -1;
        destination_node_no = -1;

        AddGPSPointsIntoGridSystem(agent_no); //find the origin and destination nodes
        int number_of_nodes = g_mm_node_vector.size();
        int i;
        for (i = 0; i < number_of_nodes; i++) //Initialization for all non-origin nodes
        {
            m_node_status_array[i] = 0; // not scanned
            m_node_label_cost[i] = MAX_LABEL_COST_;
            m_link_predecessor[i] = -1; // pointer to previous NODE INDEX from the current label at current node and time
            m_node_predecessor[i] = -1; // pointer to previous NODE INDEX from the current label at current node and time
            m_label_time_array[i] = 0;
            // comment out to speed up comuting
            ////m_label_distance_array[i] = 0;
        }

        int internal_debug_flag = 0;
        if (origin_node_no == -1 || destination_node_no == -1)
        {
            return 0;
        }

        for (int l = 0; l < g_mm_link_vector.size(); l++)
        {
            if (allowed_link_type_code.size() > 0 && g_mm_link_vector[l].link_type_code.size() > 0)
            {
                string s1 = allowed_link_type_code;
                string s2 = g_mm_link_vector[l].link_type_code;
                if (s1.find(s2) == std::string::npos)  // s2 is not contained in s1
                {
                    m_link_generalised_cost_array[l] = MAX_LABEL_COST_ / 100.0;
                }
            }

            if (blocked_link_type_code.size() > 0 && g_mm_link_vector[l].link_type_code.size() > 0)
            {
                string s1 = blocked_link_type_code;
                string s2 = g_mm_link_vector[l].link_type_code;
                if (s1.find(s2) != std::string::npos)  // s2 is contained in s1
                {
                    m_link_generalised_cost_array[l] = MAX_LABEL_COST_ / 100.0;
                }
            }

        }


        cout << "origin_node_id =" << g_mm_node_vector[origin_node_no].node_id << "; "
            << "destination_node_id =" << g_mm_node_vector[destination_node_no].node_id << endl;

        //Initialization for origin node at the preferred departure time, at departure time, cost = 0, otherwise, the delay at origin node
        m_label_time_array[origin_node_no] = g_agent_vector[agent_no].start_time_in_min;
        m_node_label_cost[origin_node_no] = 0.0;
        //Mark:	m_label_distance_array[origin_node_no] = 0.0;

        m_link_predecessor[origin_node_no] = -1; // pointer to previous NODE INDEX from the current label at current node and time
        m_node_predecessor[origin_node_no] = -1; // pointer to previous NODE INDEX from the current label at current node and time

        SEList_clear();
        SEList_push_back(origin_node_no);

        int from_node, to_node;
        int link_sqe_no;
        double new_time = 0;
        double new_distance = 0;
        double new_to_node_cost = 0;
        int tempFront;
        while (!(m_ListFront == -1)) //SEList_empty()
        {
            //          from_node = SEList_front();
            //			SEList_pop_front();  // remove current node FromID from the SE list

            from_node = m_ListFront; //pop a node FromID for scanning
            tempFront = m_ListFront;
            m_ListFront = m_SENodeList[m_ListFront];
            m_SENodeList[tempFront] = -1;

            m_node_status_array[from_node] = 2;

            int pred_link_seq_no = m_link_predecessor[from_node];

            for (i = 0; i < g_mm_node_vector[from_node].m_outgoing_link_seq_no_vector.size(); i++) // for each link (i,j) belong A(i)
            {

                link_sqe_no = g_mm_node_vector[from_node].m_outgoing_link_seq_no_vector[i];
                to_node = g_mm_link_vector[link_sqe_no].to_node_seq_no;

                //remark: the more complicated implementation can be found in paper Shortest Path Algorithms In Transportation Models: Classical and Innovative Aspects
                //	A note on least time path computation considering delays and prohibitions for intersection movements
                //very important: only origin zone can access the outbound connectors,
                //the other zones do not have access to the outbound connectors

                new_to_node_cost = m_node_label_cost[from_node] + m_link_generalised_cost_array[link_sqe_no]; // m_link_generalised_cost_array is the likelihood cost determined externally by the GPS points to link distance

                if (new_to_node_cost < m_node_label_cost[to_node]) // we only compare cost at the downstream node ToID at the new arrival time t
                {

                    m_label_time_array[origin_node_no] = 0;
                    m_node_label_cost[to_node] = new_to_node_cost;


                    m_node_predecessor[to_node] = from_node;   // pointer to previous physical NODE INDEX from the current label at current node and time
                    m_link_predecessor[to_node] = link_sqe_no; // pointer to previous physical NODE INDEX from the current label at current node and time

                    // deque updating rule for m_node_status_array
                    if (m_node_status_array[to_node] == 0)
                    {
                        ///// SEList_push_back(to_node);
                        ///// begin of inline block
                        if (m_ListFront == -1) // start from empty
                        {
                            m_ListFront = to_node;
                            m_ListTail = to_node;
                            m_SENodeList[to_node] = -1;
                        }
                        else
                        {
                            m_SENodeList[m_ListTail] = to_node;
                            m_SENodeList[to_node] = -1;
                            m_ListTail = to_node;
                        }
                        ///// end of inline block

                        m_node_status_array[to_node] = 1;
                    }
                    if (m_node_status_array[to_node] == 2)
                    {
                        /////SEList_push_front(to_node);
                        ///// begin of inline block
                        if (m_ListFront == -1) // start from empty
                        {
                            m_SENodeList[to_node] = -1;
                            m_ListFront = to_node;
                            m_ListTail = to_node;
                        }
                        else
                        {
                            m_SENodeList[to_node] = m_ListFront;
                            m_ListFront = to_node;
                        }
                        ///// end of inline block

                        m_node_status_array[to_node] = 1;
                    }
                }
            }
        }

        return m_node_label_cost[destination_node_no];
    }

    double time_dependent_label_correcting(int agent_no)
    {
        origin_node_no = -1;
        destination_node_no = -1;
        for (int link = 0; link < g_mm_link_vector.size(); link++)
        {
            g_mm_link_vector[link].hit_count = 0;
            g_mm_link_vector[link].use_count = 0;
            g_mm_link_vector[link].Possible_Dwell_time_in_min = 0;

        }
        if (AddGPSPointsIntoGridSystem(agent_no) == false)
            return 0; //step 1: find the origin and destination nodes


        int LR_iteration_size = 1;  // step 2: LR 

        int balance_count = 0;

        for (int lr_i = 0; lr_i < LR_iteration_size; lr_i++)
        {
            if (lr_i >= 1)
            {
                int balance_count_current = UpdateLinkLRPriceGridSystem(agent_no);

                if (lr_i >= 2)
                {
                    int balance_count_change = balance_count_current - balance_count;

                    if (balance_count_current == 0)  // no improvement
                        break;


                }
                balance_count = balance_count_current;


            }


            int number_of_nodes = g_mm_node_vector.size();
            int i;
            for (i = 0; i < number_of_nodes; i++) //Initialization for all non-origin nodes
            {
                if (g_mm_node_vector[i].bInsideFlag == false)
                    continue;

                for (int t = 0; t < g_TimeRangeInterval; t++)
                {
                    m_TD_node_label_cost[i][t] = MAX_LABEL_COST_;
                    m_TD_link_predecessor[i][t] = -1; // pointer to previous NODE INDEX from the current label at current node and time
                    m_TD_node_predecessor[i][t] = -1; // pointer to previous NODE INDEX from the current label at current node and time
                    m_TD_time_predecessor[i][t] = -1; // pointer to previous NODE INDEX from the current label at current node and time
                }
            }
            fprintf(g_pFileLog, "OD node index: %d -> %d\n", origin_node_no, destination_node_no);

            int internal_debug_flag = 0;
            if (origin_node_no == -1 || destination_node_no == -1)
            {

                return 0;
            }

            cout << "origin_node_no =" << g_mm_node_vector[origin_node_no].node_id << "; "
                << "destination_node_no =" << g_mm_node_vector[destination_node_no].node_id << endl;

            fprintf(g_pFileLog, "OD node id: %d -> %d\n", g_mm_node_vector[origin_node_no].node_id, g_mm_node_vector[destination_node_no].node_id);

            //Initialization for origin node at the preferred departure time, at departure time, cost = 0, otherwise, the delay at origin node
            m_label_time_array[origin_node_no] = g_agent_vector[agent_no].start_time_in_min;

            int start_timestamp_in_interval = g_GetTimeInterval(g_agent_vector[agent_no].start_time_in_min);

            g_agent_vector[agent_no].end_time_in_min = g_agent_vector[agent_no].m_GPSPointVector[g_agent_vector[agent_no].m_GPSPointVector.size() - 1].global_time_in_second / 60.0;

            int end_timestamp_in_interval = g_GetTimeInterval(g_agent_vector[agent_no].end_time_in_min);
            m_TD_node_label_cost[origin_node_no][start_timestamp_in_interval] = 0.0;

            int from_node, to_node;
            int link_sqe_no;
            double new_time = 0;
            double new_distance = 0;
            double new_to_node_cost = 0;

            fprintf(g_pFileLog, "finding time index: %d -> %d\n", start_timestamp_in_interval, end_timestamp_in_interval);

            for (int tr = start_timestamp_in_interval; tr <= end_timestamp_in_interval; tr++)  // tr as current time
            {
                int update_count = 0;

                for (int link = 0; link < g_mm_link_vector.size(); link++)
                {
                    if (g_mm_link_vector[link].bInsideFlag == false)
                        continue;

                    from_node = g_mm_link_vector[link].from_node_seq_no;
                    to_node = g_mm_link_vector[link].to_node_seq_no;

                    //if (tr == 60 && g_mm_link_vector[link].from_node_id == 9175 && g_mm_link_vector[link].to_node_id == 9212)
                    //{
                    //    TRACE("");



                    //}

                    float timestamp_in_min = tr * g_TimeResolution_inMin + g_StartTimeinMin;
                    int timestamp_in_interval = max(0, g_GetTimeInterval(timestamp_in_min));

                    int max_TT_Dwell_time_in_int = 0;

                    if (g_agent_vector[agent_no].m_ExcessiveDwellMap.find(tr) != g_agent_vector[agent_no].m_ExcessiveDwellMap.end())
                    {
                        max_TT_Dwell_time_in_int = g_agent_vector[agent_no].m_ExcessiveDwellMap[timestamp_in_interval];
                    }

                    float over_speed_limit_ratio = 0.7;
                    int min_FFTT_in_interval = max(1.0, g_mm_link_vector[link].FFTT_in_min / g_TimeResolution_inMin * over_speed_limit_ratio + 0.5);
                    int FFTT_in_interval = max(1.0, g_mm_link_vector[link].FFTT_in_min / g_TimeResolution_inMin + 0.5);
                    int stage_size = 10;
                    int max_TT_in_interval = max(1, max(max_TT_Dwell_time_in_int, FFTT_in_interval * stage_size));
                    int step_size = 1;

                    for (int travel_time_in_interval = min_FFTT_in_interval; travel_time_in_interval < max_TT_in_interval; travel_time_in_interval += step_size) //++1 might give more precise travel time
                    {
                        int final_tt_in_interval = travel_time_in_interval;
                        //dynamic step size 
                        //float ratio = travel_time_in_interval * 1.0 / FFTT_in_interval;
                        //if (ratio > 5 && travel_time_in_interval* g_TimeResolution_inMin>=2)
                        //{
                        //    step_size = max(1.0, ratio / 3); // to speed up calculation
                        //}

                         //last stage
                        //if (travel_time_in_interval == max_TT_in_interval-1 && g_mm_link_vector[link].Possible_Dwell_time_in_min>=1 && 
                        //    (timestamp_in_min>g_mm_link_vector[link].dwell_start_time_in_min && timestamp_in_min < g_mm_link_vector[link].dwell_start_time_in_min+1))
                        //{
                        //    final_tt_in_interval = g_mm_link_vector[link].Possible_Dwell_time_in_min/ g_TimeResolution_inMin;  // hour

                        //}


                        int tnext = min(g_TimeRangeInterval - 1, tr + final_tt_in_interval);

                        float benefit = 0;

                        float speed_deviation_penalty = 0;  // only apply for normal cases

                        //if (max_TT_Dwell_time_in_int == 0)
                        //{
                        //    speed_deviation_penalty = abs(final_tt_in_interval - FFTT_in_interval) * 0.000001;
                        //}
                        //else
                        //{

                        //}

                        new_to_node_cost = m_TD_node_label_cost[from_node][tr] + m_TD_link_generalised_cost_array[link][timestamp_in_interval] + benefit + speed_deviation_penalty; // m_link_generalised_cost_array is the likelihood cost determined externally by the GPS points to link distance

                        if (new_to_node_cost < m_TD_node_label_cost[to_node][tnext]) // we only compare cost at the downstream node ToID at the new arrival time t
                        {

                            m_TD_node_label_cost[to_node][tnext] = new_to_node_cost;
                            m_TD_node_predecessor[to_node][tnext] = from_node;
                            m_TD_link_predecessor[to_node][tnext] = link;
                            m_TD_time_predecessor[to_node][tnext] = tr;

                            if (g_mm_link_vector[link].AccessibilityTime > new_to_node_cost)
                            {
                                g_mm_link_vector[link].AccessibilityTime = tnext;
                            }


                            //if(agent_no==0)
                            //{
                            //fprintf(g_pFileLog, "%d: %d, link %d, time %d-> %d: cost %f \n", update_count, link, g_mm_link_vector[link].to_node_id, tr, tnext, new_to_node_cost);
                            //}

                            update_count++;
                        }

                    }
                }


            }

            int	node_size = 0;
            std::vector<int>  reversed_path_node_sequence, reversed_path_link_sequence, reversed_path_time_sequence;
            std::vector<float>  reversed_path_cost_sequence;

            int min_end_timestamp_in_interval = end_timestamp_in_interval;
            float min_label_cost = m_TD_node_label_cost[destination_node_no][end_timestamp_in_interval];

            for (int tt = max(0, end_timestamp_in_interval - 5); tt < end_timestamp_in_interval; tt++)
            {
                if (m_TD_node_label_cost[destination_node_no][tt] < min_label_cost)
                {
                    min_label_cost = m_TD_node_label_cost[destination_node_no][tt];
                    min_end_timestamp_in_interval = tt;
                }

            }

            reversed_path_node_sequence.push_back(destination_node_no); //record the first node backward, destination node
            reversed_path_time_sequence.push_back(min_end_timestamp_in_interval); //record the first node backward, destination node
            reversed_path_cost_sequence.push_back(min_label_cost);

            fprintf(g_pFileLog, "backtracing time index: %d\n", min_end_timestamp_in_interval);

            int pred_node = m_TD_node_predecessor[destination_node_no][end_timestamp_in_interval];
            int pred_time_r = m_TD_time_predecessor[destination_node_no][end_timestamp_in_interval];
            int pred_link = -1;
            int pred_node_record, pred_time_record_r;

            while (pred_node != -1 && pred_time_r >= start_timestamp_in_interval) // scan backward in the predessor array of the shortest path calculation results
            {

                //record current values of node and time predecessors, and update PredNode and PredTime
                pred_node_record = pred_node;
                pred_time_record_r = pred_time_r;

                pred_node = m_TD_node_predecessor[pred_node_record][pred_time_record_r];
                pred_time_r = m_TD_time_predecessor[pred_node_record][pred_time_record_r];
                pred_link = m_TD_link_predecessor[pred_node_record][pred_time_record_r];



                if (pred_node != -1)
                {

                    fprintf(g_pFileLog, "backtracing time index: %d", pred_time_r);
                    fprintf(g_pFileLog, "backtracing link %d: %d->%d\n", pred_link, g_mm_link_vector[pred_link].from_node_id, g_mm_link_vector[pred_link].to_node_id);

                    reversed_path_link_sequence.push_back(pred_link);

                    reversed_path_node_sequence.push_back(pred_node);
                    reversed_path_time_sequence.push_back(pred_time_r);
                    reversed_path_cost_sequence.push_back(m_TD_node_label_cost[pred_node_record][pred_time_record_r]);

                }

            }

            g_agent_vector[agent_no].AllocatePathNodeVector(reversed_path_node_sequence.size(), reversed_path_node_sequence, reversed_path_time_sequence, reversed_path_cost_sequence, reversed_path_link_sequence);

            fprintf(g_pFileLog, "AllocatePathNodeVector\n");

            if (origin_node_no >= 0 && destination_node_no >= 0) //feasible origin and destination nodes
            {
                g_agent_vector[agent_no].o_node_id = g_mm_node_vector[origin_node_no].node_id;
                g_agent_vector[agent_no].d_node_id = g_mm_node_vector[destination_node_no].node_id;
            }

            for (int i = 0; i < reversed_path_link_sequence.size(); i++)
            {
                int link_no = reversed_path_link_sequence[i];
                g_mm_link_vector[link_no].use_count += 1;

            }

        }
        return m_node_label_cost[destination_node_no];
    }

    void find_path_for_agents_assigned_for_this_thread()
    {
        int return_value;

        for (int i = 0; i < m_agent_vector.size(); i++)
        {
            CAgent* p_agent = &(g_agent_vector[m_agent_vector[i]]);

            cout << "agent_id =" << p_agent->agent_id.c_str() << endl;

            return_value = optimal_label_correcting(p_agent->agent_no, p_agent->allowed_link_type_code, p_agent->blocked_link_type_code);

            if (return_value == -1)
            {
                continue;
            }
            // step 2: backtrack to the origin (based on node and time predecessors)
            int current_node_seq_no = destination_node_no; // destination node
            int current_link_seq_no = -1;

            int l_node_size = 0;
            // backtrace the shortest path tree from the destination to the root (at origin)
            while (current_node_seq_no >= 0 && current_node_seq_no < g_number_of_nodes)
            {

                temp_path_node_vector[l_node_size++] = current_node_seq_no;

                if (l_node_size >= 10000)
                {
                    cout << "Error: l_node_size >= temp_path_node_vector_size" << endl;
                    g_Program_stop();
                }

                current_node_seq_no = m_node_predecessor[current_node_seq_no]; // update node seq no
            }

            p_agent->AllocatePathNodeVector(l_node_size, temp_path_node_vector, true);
            if (origin_node_no >= 0 && destination_node_no >= 0) //feasible origin and destination nodes
            {
                p_agent->o_node_id = g_mm_node_vector[origin_node_no].node_id;
                p_agent->d_node_id = g_mm_node_vector[destination_node_no].node_id;
            }

            for (int i = 0; i < p_agent->m_node_size - 2; i++)
            {
                int link_no = g_mm_node_vector[p_agent->path_node_vector[i]].m_outgoing_link_seq_no_map[p_agent->path_node_vector[i + 1]];

                if (link_no < g_mm_link_vector.size())
                {
                    p_agent->likely_trace_no_vector.push_back(g_mm_link_vector[link_no].likely_trace_no);
                }
            }

            for (int link = 0; link < g_mm_link_vector.size(); link++)
            {
                g_mm_link_vector[link].likely_trace_no = -1;
            }
        }
    }
    void find_TD_path_for_agents_assigned_for_this_thread()
    {
        int return_value;

        for (int i = 0; i < m_agent_vector.size(); i++)
        {
            CAgent* p_agent = &(g_agent_vector[m_agent_vector[i]]);

            cout << "agent_id =" << p_agent->agent_id.c_str() << endl;


            return_value = time_dependent_label_correcting(p_agent->agent_no);

            if (return_value == -1)
            {
                continue;
            }
        }
    }
};
NetworkForSP* g_pNetworkVector = nullptr;

vector<double> g_timestr2second(string str)
{
    vector<double> output_global_minute;
    vector<double> output_global_second;

    int string_lenghth = str.length();

    const char* string_line = str.data(); //string to char*

    int char_length = strlen(string_line);

    char ch, buf_ddhhmm[32] = { 0 }, buf_SS[32] = { 0 }, buf_sss[32] = { 0 };
    char dd1, dd2, hh1, hh2, mm1, mm2, SS1, SS2, sss1, sss2, sss3;
    double ddf1, ddf2, hhf1, hhf2, mmf1, mmf2, SSf1, SSf2, sssf1, sssf2, sssf3;
    double global_minute = 0;
    double dd = 0, hh = 0, mm = 0, SS = 0, sss = 0;
    int i = 0;
    int buffer_i = 0, buffer_k = 0, buffer_j = 0;
    int num_of_colons = 0;

    //DDHHMM:SS:sss or HHMM:SS:sss

    while (i < char_length)
    {
        ch = string_line[i++];

        if (num_of_colons == 0 && ch != '_' && ch != ':') //input to buf_ddhhmm until we meet the colon
        {
            buf_ddhhmm[buffer_i++] = ch;
        }
        else if (num_of_colons == 1 && ch != ':') //start the Second "SS"
        {
            buf_SS[buffer_k++] = ch;
        }
        else if (num_of_colons == 2 && ch != ':') //start the Millisecond "sss"
        {
            buf_sss[buffer_j++] = ch;
        }

        if (ch == '_' || ch == ';' || i == char_length) //start a new time string
        {
            if (buffer_i == 4) //"HHMM"
            {
                //HHMM, 0123
                hh1 = buf_ddhhmm[0]; //read each first
                hh2 = buf_ddhhmm[1];
                mm1 = buf_ddhhmm[2];
                mm2 = buf_ddhhmm[3];

                hhf1 = ((double)hh1 - 48); //convert a char to a double
                hhf2 = ((double)hh2 - 48);
                mmf1 = ((double)mm1 - 48);
                mmf2 = ((double)mm2 - 48);

                dd = 0;
                hh = hhf1 * 10 * 60 + hhf2 * 60;
                mm = mmf1 * 10 + mmf2;
            }
            else if (buffer_i == 6) //"DDHHMM"
            {
                //DDHHMM, 012345
                dd1 = buf_ddhhmm[0]; //read each first
                dd2 = buf_ddhhmm[1];
                hh1 = buf_ddhhmm[2];
                hh2 = buf_ddhhmm[3];
                mm1 = buf_ddhhmm[4];
                mm2 = buf_ddhhmm[5];

                ddf1 = ((double)dd1 - 48); //convert a char to a double
                ddf2 = ((double)dd2 - 48);
                hhf1 = ((double)hh1 - 48);
                hhf2 = ((double)hh2 - 48);
                mmf1 = ((double)mm1 - 48);
                mmf2 = ((double)mm2 - 48);

                dd = ddf1 * 10 * 24 * 60 + ddf2 * 24 * 60;
                hh = hhf1 * 10 * 60 + hhf2 * 60;
                mm = mmf1 * 10 + mmf2;
            }

            if (num_of_colons == 1 || num_of_colons == 2)
            {
                //SS, 01
                SS1 = buf_SS[0]; //read each first
                SS2 = buf_SS[1];

                SSf1 = ((double)SS1 - 48); //convert a char to a double
                SSf2 = ((double)SS2 - 48);

                SS = (SSf1 * 10 + SSf2) / 60;
            }

            if (num_of_colons == 2)
            {
                //sss, 012
                sss1 = buf_sss[0]; //read each first
                sss2 = buf_sss[1];
                sss3 = buf_sss[2];

                sssf1 = ((double)sss1 - 48); //convert a char to a double
                sssf2 = ((double)sss2 - 48);
                sssf3 = ((double)sss3 - 48);

                sss = (sssf1 * 100 + sssf2 * 10 + sssf3) / 1000;
            }

            global_minute = dd + hh + mm + SS + sss;
            double global_second = (dd + hh + mm + SS + sss) * 60.0;

            output_global_second.push_back(global_second);

            //initialize the parameters
            buffer_i = 0;
            buffer_k = 0;
            buffer_j = 0;
            num_of_colons = 0;
        }

        if (ch == ':')
        {
            num_of_colons += 1;
        }
    }

    return output_global_second;
}
int timestr2second(string time_str)
{ //hhmmss
    string hh = time_str.substr(0, 2);
    string mm = time_str.substr(2, 2);
    string ss = time_str.substr(5, 2);
    int hhi = stoi(hh);
    int mmi = stoi(mm);
    int ssi = stoi(ss);
    return hhi * 3600 + mmi * 60 + ssi;
}

string second2timestr(int time_int)
{
    int hhi = time_int / 3600;
    int mmi = (time_int - 3600 * hhi) / 60;
    int ssi = time_int - 3600 * hhi - 60 * mmi;
    string hh = hhi < 10 ? "0" + std::to_string(hhi) : std::to_string(hhi);
    string mm = mmi < 10 ? "0" + std::to_string(mmi) : std::to_string(mmi);
    string ss = ssi < 10 ? "0" + std::to_string(ssi) : std::to_string(ssi);
    return hh + mm + ":" + ss;
}

void g_DetermineResolution()
{

    CDTACSVParser parser;
    if (parser.OpenCSVFile("node.csv", true))
    {
        // initialization of grid rectangle boundary
        int node_count = 0;
        while (parser.ReadRecord()) // if this line contains [] mark, then we will also read field headers.
        {
            int node_id;

            if (parser.GetValueByFieldName("node_id", node_id) == false)
                continue;

            CMapmatchingNode node;
            node.node_id = node_id;

            parser.GetValueByFieldName("x_coord", node.pt.x, false);
            parser.GetValueByFieldName("y_coord", node.pt.y, false);


            // exapnd the grid boundary according to the nodes
            g_left = min(g_left, node.pt.x);
            g_right = max(g_right, node.pt.x);
            g_top = max(g_top, node.pt.y);
            g_bottom = min(g_bottom, node.pt.y);
            node_count++;
        }

        int grid_size = 5;

        if (node_count > 10000)
            grid_size = 8;

        if (node_count > 40000)
            grid_size = 10;


        double temp_resolution = (((g_right - g_left) / grid_size + (g_top - g_bottom) / grid_size)) / 2.0;

        vector<double> ResolutionVector;

        ResolutionVector.push_back(0.00005);
        ResolutionVector.push_back(0.0001);
        ResolutionVector.push_back(0.0002);
        ResolutionVector.push_back(0.0005);
        ResolutionVector.push_back(0.001);
        ResolutionVector.push_back(0.002);
        ResolutionVector.push_back(0.005);
        ResolutionVector.push_back(0.01);
        ResolutionVector.push_back(0.02);
        ResolutionVector.push_back(0.05);
        ResolutionVector.push_back(0.1);
        ResolutionVector.push_back(0.2);
        ResolutionVector.push_back(0.5);
        ResolutionVector.push_back(1);
        ResolutionVector.push_back(2);
        ResolutionVector.push_back(5);
        ResolutionVector.push_back(10);
        ResolutionVector.push_back(20);
        ResolutionVector.push_back(50);



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

        g_GridResolution = temp_resolution;

        cout << "g_GridResolution = " << temp_resolution << endl;

        parser.CloseCSVFile();
    }

}
void g_ReadInputData()
{
    //g_DetermineResolution();

    CDTACSVParser parser;
    if (parser.OpenCSVFile("node.csv", true))
    {

        while (parser.ReadRecord()) // if this line contains [] mark, then we will also read field headers.
        {
            int node_id;

            if (parser.GetValueByFieldName("node_id", node_id) == false)
                continue;

            if (g_internal_node_seq_no_map.find(node_id) != g_internal_node_seq_no_map.end())
            {
                cout << "warning: duplicate definition of node " << node_id << " was detected\n";
                continue;
            }


            int zone_id = 0;

            CMapmatchingNode node;
            node.node_id = node_id;
            node.node_seq_no = g_number_of_nodes++;
            parser.GetValueByFieldName("x_coord", node.pt.x, false);
            parser.GetValueByFieldName("y_coord", node.pt.y, false);
            parser.GetValueByFieldName("zone_id", zone_id, true);

            // exapnd the grid boundary according to the nodes
            g_left = min(g_left, node.pt.x);
            g_right = max(g_right, node.pt.x);
            g_top = max(g_top, node.pt.y);
            g_bottom = min(g_bottom, node.pt.y);

            if (g_time_dependent_computing_mode == 1)
            {
                zone_id = node_id;
            }

            __int64 cell_id = g_GetCellID(node.pt.x, node.pt.y);
            node.cell_id = cell_id;
            g_cell_id_2_node_map[cell_id] = 1;

            if (zone_id > 0) // zone id is feasible
            {
                if (g_cell_id_2_zone_id_map.find(cell_id) == g_cell_id_2_zone_id_map.end())  // not defined yet
                    g_cell_id_2_zone_id_map[cell_id] = zone_id;

            }

            g_mm_node_vector.push_back(node);
            g_internal_node_seq_no_map[node_id] = node.node_seq_no;
            if (g_number_of_nodes % 1000 == 0)
                cout << "reading " << g_number_of_nodes << " nodes.. " << endl;
        }

        //in this map matching process for planning model, for each node, we need to assign a zone id internally.so that for the map matched path, we can always find the origin zoneand destination zone along the path.

        for (int i = 0; i < g_mm_node_vector.size(); i++)  // assign a zone id 
        {
            if (g_mm_node_vector[i].zone_id <= 0)
            {
                if (g_cell_id_2_zone_id_map.find(g_mm_node_vector[i].cell_id) != g_cell_id_2_zone_id_map.end())
                {
                    //we can use the zone id in the same cell of the grid, if the zone id is not defined in node.csv
                    g_mm_node_vector[i].zone_id = g_cell_id_2_zone_id_map[g_mm_node_vector[i].cell_id];

                }
            }
        }

        cout << "number of nodes = " << g_number_of_nodes << endl;
        parser.CloseCSVFile();
    }
    else
    {
        cout << "Cannot open file node.csv" << endl;
        g_Program_stop();
    }

    CDTACSVParser parser_link;
    if (parser_link.OpenCSVFile("link.csv", true))
    {
        while (parser_link.ReadRecord())
        {
            CMapmatchingLink link;

            if (parser_link.GetValueByFieldName("link_id", link.link_id) == false)
                continue;
            if (parser_link.GetValueByFieldName("from_node_id", link.from_node_id) == false)
                continue;
            if (parser_link.GetValueByFieldName("to_node_id", link.to_node_id) == false)
                continue;

            parser_link.GetValueByFieldName("length", link.length);

            parser_link.GetValueByFieldName("free_speed", link.free_speed);

            parser_link.GetValueByFieldName("capacity", link.lane_capacity);

            parser_link.GetValueByFieldName("lanes", link.lanes);

            parser_link.GetValueByFieldName("VDF_fftt1", link.FFTT_in_min);
            parser_link.GetValueByFieldName("link_type_code", link.link_type_code);
            parser_link.GetValueByFieldName("link_type_name", link.link_type_name);


            link.FFTT_in_min = max(0.1, link.FFTT_in_min);
            link.FFTT_in_sec = link.FFTT_in_min * 60.0;


            if (g_internal_node_seq_no_map.find(link.from_node_id) == g_internal_node_seq_no_map.end())
            {
                cout << "warning: from_node_id " << link.from_node_id << " of link " << link.link_id << " has not been defined in node.csv\n";
                continue;
            }
            if (g_internal_node_seq_no_map.find(link.to_node_id) == g_internal_node_seq_no_map.end())
            {
                cout << "warning: to_node_id " << link.to_node_id << " of link " << link.link_id << " has not been defined in node.csv\n";
                continue;
            }

            string geometry_str;
            parser_link.GetValueByFieldName("geometry", geometry_str);

            link.from_node_seq_no = g_internal_node_seq_no_map[link.from_node_id];
            link.to_node_seq_no = g_internal_node_seq_no_map[link.to_node_id];
            link.geometry = geometry_str;

            // overwrite when the field "geometry" exists
            CGeometry geometry(geometry_str);
            std::vector<CCoordinate> CoordinateVector;
            CoordinateVector = geometry.GetCoordinateList();

            link.link_distance = 0;
            GDPoint Point;
            GDPoint Point_next;

            if (CoordinateVector.size() >= 2)
            {
                for (int l = 0; l < CoordinateVector.size(); l++)
                {

                    Point.x = CoordinateVector[l].X;
                    Point.y = CoordinateVector[l].Y;
                    //				GPSPoint.time_str = time_stamp;

                    link.m_PointVector.push_back(Point);

                    if (l < CoordinateVector.size() - 1) // consider shape points in each segment of a link
                    {
                        Point_next.x = CoordinateVector[l + 1].X;
                        Point_next.y = CoordinateVector[l + 1].Y;

                        link.link_distance += g_GetPoint2Point_Distance(&Point, &Point_next);
                    }
                }
                link.seg_distance = link.link_distance / (CoordinateVector.size() - 1); // take the average of segment distance
            }
            else  // no geometry
            {
                link.m_PointVector.push_back(g_mm_node_vector[link.from_node_seq_no].pt);
                link.m_PointVector.push_back(g_mm_node_vector[link.to_node_seq_no].pt);

                link.link_distance += g_GetPoint2Point_Distance(&(g_mm_node_vector[link.from_node_seq_no].pt), &(g_mm_node_vector[link.to_node_seq_no].pt));
                link.seg_distance = link.link_distance;

            }


            link.link_seq_no = g_number_of_links++;

            g_internal_link_no_map[link.link_id] = link.link_seq_no;
            g_mm_node_vector[link.from_node_seq_no].m_outgoing_link_seq_no_vector.push_back(link.link_seq_no);
            g_mm_node_vector[link.from_node_seq_no].m_outgoing_link_seq_no_map[link.to_node_seq_no] = link.link_seq_no;
            g_mm_link_vector.push_back(link);

            if (g_number_of_links % 1000 == 0)
                cout << "reading " << g_number_of_links << " links.. " << endl;
        }

        cout << "number of links = " << g_number_of_links << endl;
        parser_link.CloseCSVFile();
    }
    else
    {
        cout << "Cannot open file link.csv" << endl;
        g_Program_stop();
    }

}

bool g_ReadInputTrajectoryCSVFile()
{
    CDTACSVParser gps_parser;
    int gps_point_count = 0;
    if (gps_parser.OpenCSVFile("input_traejctory.csv", true))
    {
        double x, y;
        string time_stamp;
        int global_time_in_second;
        while (gps_parser.ReadRecord())
        {

            string agent_id;
            if (gps_parser.GetValueByFieldName("agent_id", agent_id) == false)
                continue;

            string time_sequence_str;
            gps_parser.GetValueByFieldName("time_sequence", time_sequence_str);

            float volume = 0;
            gps_parser.GetValueByFieldName("volume", volume);


            std::vector<int> time_sequence;
            g_ParserIntSequence(time_sequence_str, time_sequence);

            string geometry_str;
            gps_parser.GetValueByFieldName("geometry", geometry_str);

            // overwrite when the field "geometry" exists
            CGeometry geometry(geometry_str);
            std::vector<CCoordinate> CoordinateVector;
            CoordinateVector = geometry.GetCoordinateList();

            std::vector<CGPSPoint> l_GPSPointVector;
            if (CoordinateVector.size() >= 2)
            {
                for (int l = 0; l < CoordinateVector.size(); l++)
                {

                    CGPSPoint GPSPoint;
                    GPSPoint.pt.x = CoordinateVector[l].X;
                    GPSPoint.pt.y = CoordinateVector[l].Y;

                    if (l < time_sequence.size())
                        GPSPoint.global_time_in_second = time_sequence[l];

                    __int64 cell_id = g_GetCellID(GPSPoint.pt.x, GPSPoint.pt.y);

                    if (g_cell_id_2_node_map.find(cell_id) != g_cell_id_2_node_map.end())  //only consider the GPS points passing through the subarea
                    {
                        GPSPoint.cell_id = cell_id;
                        l_GPSPointVector.push_back(GPSPoint);
                        gps_point_count++;

                    }

                }
            }


            if (l_GPSPointVector.size() <= 1) // if the trace is outside the current network
                continue;

            if (g_internal_agent_no_map.find(agent_id) == g_internal_agent_no_map.end())
            {

                CAgent agent;
                agent.agent_id = agent_id;
                agent.volume = volume;
                agent.agent_no = g_agent_vector.size();
                g_internal_agent_no_map[agent_id] = agent.agent_no;  // assign the internal agent no as the current size of the map.
                agent.m_GPSPointVector = l_GPSPointVector;

                g_agent_vector.push_back(agent);
            }

        }

        cout << "number of agents = " << g_agent_vector.size() << endl;
        cout << "number of GPS points = " << gps_point_count << endl;

        gps_parser.CloseCSVFile();

        return true;
    }

    return false;
}




void g_OutputCell2ZoneCSVFile()
{
    FILE* g_pFileCell2Zone = nullptr;
    g_pFileCell2Zone = fopen("cell2zone.csv", "w");

    if (g_pFileCell2Zone == NULL)
    {
        cout << "File cell2zone.csv cannot be opened." << endl;
        g_Program_stop();
    }
    else
    {
        fprintf(g_pFileCell2Zone, "cell_id,zone_id\n");

        map<__int64, int>::iterator it;

        for (it = g_cell_id_2_zone_id_map.begin(); it != g_cell_id_2_zone_id_map.end(); it++)
        {
            fprintf(g_pFileCell2Zone, "%jd,%d\n",
                it->first,
                it->second);
        }



        fclose(g_pFileCell2Zone);
    }
}

void g_OutputDemandCSVFile()
{
    FILE* g_pFileAgent = nullptr;
    g_pFileAgent = fopen("demand.csv", "w");

    if (g_pFileAgent == NULL)
    {
        cout << "File agent.csv cannot be opened." << endl;
        //  g_Program_stop();
    }
    else
    {
        fprintf(g_pFileAgent, "agent_id,o_node_id,d_node_id,o_zone_id,d_zone_id,volume,geometry\n");

        for (int a = 0; a < g_agent_vector.size(); a++)
        {
            CAgent* p_agent = &(g_agent_vector[a]);
            int matching_link_from_node_id = -1;
            int matching_link_to_node_id = -1;
            string matching_link_id;

            if (p_agent->path_node_vector.size() > 0)
            {
                p_agent->o_cell_id = g_mm_node_vector[p_agent->path_node_vector[0]].cell_id;
                p_agent->d_cell_id = g_mm_node_vector[p_agent->path_node_vector[p_agent->path_node_vector.size() - 1]].cell_id;

                p_agent->origin_zone_id = g_mm_node_vector[p_agent->path_node_vector[0]].zone_id;
                p_agent->destination_zone_id = g_mm_node_vector[p_agent->path_node_vector[p_agent->path_node_vector.size() - 1]].zone_id;
            }


            if (p_agent->matching_link_no >= 0)
            {
                matching_link_from_node_id = g_mm_link_vector[p_agent->matching_link_no].from_node_id;
                matching_link_to_node_id = g_mm_link_vector[p_agent->matching_link_no].to_node_id;
                matching_link_id = g_mm_link_vector[p_agent->matching_link_no].link_id;
            }

            p_agent->distance = 0;
            for (int i = 0; i < p_agent->m_node_size - 2; i++)
            {
                int link_no = g_mm_node_vector[p_agent->path_node_vector[i]].m_outgoing_link_seq_no_map[p_agent->path_node_vector[i + 1]];

                p_agent->distance += g_mm_link_vector[link_no].length;
            }

            p_agent->travel_time = p_agent->end_time_in_min - p_agent->start_time_in_min;
            fprintf(g_pFileAgent, "%s,%d,%d,%d,%d,%f",
                p_agent->agent_id.c_str(),
                p_agent->o_node_id,
                p_agent->d_node_id,
                p_agent->origin_zone_id,
                p_agent->destination_zone_id,
                p_agent->volume
            );

            fprintf(g_pFileAgent, "\"LINESTRING (");

            int shape_point_count = 0;
            for (int i = 0; i < p_agent->path_link_vector.size(); i++)
            {
                int link_no = p_agent->path_link_vector[i];


                for (int gl = 0; gl < g_mm_link_vector[link_no].m_PointVector.size(); gl++)
                {
                    fprintf(g_pFileAgent, "%f %f,", g_mm_link_vector[link_no].m_PointVector[gl].x,
                        g_mm_link_vector[link_no].m_PointVector[gl].y);
                    shape_point_count++;
                }
            }

            if (shape_point_count == 0)  // link shape points do not exist
            {
                for (int i = 0; i < p_agent->path_node_vector.size(); i++)
                {
                    fprintf(g_pFileAgent, "%f %f,", g_mm_node_vector[p_agent->path_node_vector[i]].pt.x,
                        g_mm_node_vector[p_agent->path_node_vector[i]].pt.y);
                }

            }

            fprintf(g_pFileAgent, ")\"");
            fprintf(g_pFileAgent, "\n");
        }

        fclose(g_pFileAgent);
    }
}

bool g_LikelyRouteFinding()
{
    //int number_of_threads = g_number_of_CPU_threads();
    int number_of_threads = 1;
    g_pNetworkVector = new NetworkForSP[number_of_threads]; // create n copies of network, each for a subset of agents to use

    cout << "number of CPU threads = " << number_of_threads << endl;

    NetworkForSP* p_Network;

    for (int i = 0; i < number_of_threads; i++)
    {
        g_pNetworkVector[i].AllocateMemory(g_number_of_nodes, g_number_of_links);
        g_pNetworkVector[i].BuildGridSystem(); // called once
    }

    for (int a = 0; a < g_agent_vector.size(); a++) //
    {

        p_Network = &g_pNetworkVector[a % number_of_threads];

        p_Network->m_agent_vector.push_back(a);
    }

    //#pragma omp parallel for
    for (int thread_no = 0; thread_no < number_of_threads; thread_no++)
    {
            g_pNetworkVector[thread_no].find_path_for_agents_assigned_for_this_thread();
            g_OutputDemandCSVFile();

    }

    cout << "End of Sequential Optimization Process. " << endl;
    fprintf(g_pFileLog, "end of optimization process\n");
    return true;
}
int trace2od()
{
    clock_t start_t, end_t, total_t;

    if (g_ReadInputTrajectoryCSVFile() == false)
    {
        return 0;

    }
    g_pFileLog = fopen("log.txt", "w");

    if (g_pFileLog == NULL)
    {
        cout << "File log.txt cannot be opened." << endl;
        g_Program_stop();

    }
    g_ReadInputData();

    start_t = clock();
    g_LikelyRouteFinding();

    //  g_OutputLinklikelihoodCSVFile();
    end_t = clock();
    total_t = (end_t - start_t);
    cout << "CPU Running Time = " << total_t / 1000.0 << " seconds" << endl;
    cout << "free memory.." << endl;
    cout << "done." << endl;

    //g_mm_node_vector.clear();
    //g_mm_link_vector.clear();
    //g_agent_vector.clear();

    fclose(g_pFileLog);
    return 1;
}

