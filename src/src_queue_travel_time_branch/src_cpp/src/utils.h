/* Portions Copyright 2021 Xuesong Zhou and Peiheng Li
 *
 * If you help write or modify the code, please also list your names here.
 * The reason of having Copyright info here is to ensure all the modified version, as a whole, under the GPL
 * and further prevent a violation of the GPL.
 *
 * More about "How to use GNU licenses for your own software"
 * http://www.gnu.org/licenses/gpl-howto.html
 */

#ifndef GUARD_UTILS_H
#define GUARD_UTILS_H

#define BUILD_EXE  // define only for Windows Executable

#ifndef NOMINMAX

#ifndef max
#define max(a,b)            (((a) > (b)) ? (a) : (b))
#endif

#ifndef min
#define min(a,b)            (((a) < (b)) ? (a) : (b))
#endif

#endif  /* NOMINMAX */
#include <vector>
#include <map>
#include <iostream>
#include <sstream>
#include <string>
#include "teestream.h"
#include <cmath>
#include <ctime>
#include <stdio.h>
#include <string>
#include <vector>
using std::istringstream;

using std::string;

// if you are using cmake, please #include <build_config.h>
#ifndef _WIN32
#include <build_config.h>
#endif

// utilities functions
void g_ProgramStop();

void fopen_ss(FILE** file, const char* fileName, const char* mode);
float g_read_float(FILE* f);

std::vector<std::string> split(const std::string& s, const std::string& seperator);
int g_ParserIntSequence(std::string str, std::vector<int>& vect);

std::vector<float> g_time_parser(std::string str);
float g_measurement_tstamp_parser(std::string str, int& weekday_flag, int &day_of_year);
std::string g_time_coding(float time_stamp);
std::string g_time_coding_HHMM(float time_stamp);
int g_dayofweek(int y, int m, int d);
int g_dayofyear(int y, int m, int d);

// Peiheng, 04/01/21, this is just a temporary fix on logging in DTALite
// it creates another global variable (i.e. dtalog right after class DTALog)
// shared by all translation units, which is really bad. This is a common issue
// across the current implementation. It will be addressed properly in the refactoring.
class DTALog{
    std::ofstream logfile;
    teestream ts;
public:
    int db;
    int sig;
    int odme;
    int path;
    int dta;
    int ue;


DTALog(): logfile {"log.txt"}, ts {std::cout, logfile}
{
}

~DTALog() = default;

#ifdef BUILD_EXE
    teestream& output() {return ts;}
#else
    std::ofstream& output() {return logfile;}
#endif

    int& debug_level() {return db;}
    int debug_level() const {return db;}

    int& log_sig() {return sig;}
    int log_sig() const {return sig;}

    int& log_odme() {return odme;}
    int log_odme() const {return odme;}

    int& log_path() {return path;}
    int log_path() const {return 2;}

    int& log_dta() {return dta;}
    int log_dta() const {return dta;}

    int& log_ue() {return ue;}
    int log_ue() const {return ue;}
};

static DTALog dtalog;

struct GDPoint //geometry data
{
	double x;
	double y;
};

typedef struct
{
	double X, Y, Z;
}CCoordinate;


class CGeometry
{
public:
	enum GeometryType {
		POINT,
		LINE,
		POLYGON,
		UNKNOWN
	};
private:
	GeometryType m_Type;
	int m_NumOfCoordinates;
	std::vector<CCoordinate> v_Coordinates;

public:
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
			type_str.erase(type_str.find_last_not_of(" ") + 1);		// works for 'LINESTRING (....' and 'LINESTRING(....'

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
			type_str.erase(type_str.find_last_not_of(" ") + 1);		// works for 'LINESTRING (....' and 'LINESTRING(....'

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
			type_str.erase(type_str.find_last_not_of(" ") + 1);		// works for 'LINESTRING (....' and 'LINESTRING(....'

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
				return UnitGridResolution; // intersection does not fall within the segment
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

	double g_Find_P2P_Angle(const GDPoint* p1, const GDPoint* p2)
	{
		double PI_ = 3.14159265359;
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
	bool g_GetTwoPoint2LineIntersectionFlag(const GDPoint* pt0, const GDPoint* pt1, const GDPoint* FromPt, const GDPoint* ToPt)
	{
		int test_code = g_GetTwoPointLineIntersectionFlag_LargeRange(pt0, pt1, FromPt, ToPt);

		if (test_code >= 100)
			return false;
		double L2Ldistance = -1;

		test_code = g_GetTwoPointSegment2LineIntersectionFlag(pt0, pt1, FromPt, ToPt,L2Ldistance);
		if (test_code <= 1)
			return true;

		// continue, with at least intersection 

		for (int segNo = 10; segNo >= 1; segNo--)
		{

			float ratio = segNo * 1.0 / 10;
			float ratio_next = (segNo - 1) * 1.0 / 10;

			GDPoint pt_segment_0, pt_segment_1;
			pt_segment_0.x = ratio * pt0->x + (1 - ratio) * pt1->x;
			pt_segment_0.y = ratio * pt0->y + (1 - ratio) * pt1->y;

			pt_segment_1.x = ratio_next * pt0->x + (1 - ratio_next) * pt1->x;
			pt_segment_1.y = ratio_next * pt0->y + (1 - ratio_next) * pt1->y;
			double L2Ldistance = -1;
			if (g_GetTwoPointSegment2LineIntersectionFlag(&pt_segment_0, &pt_segment_1, FromPt, ToPt,L2Ldistance) <= 1)
				return true;

		}
		return false;
	}

	int g_GetTwoPointLineIntersectionFlag_LargeRange(const GDPoint* pt0, const GDPoint* pt1, const GDPoint* FromPt, const GDPoint* ToPt)
	{
		double U;
		GDPoint Intersection;
		GDPoint pt;
		pt.x = (pt0->x + pt1->x) / 2;
		pt.y = (pt0->y + pt1->y) / 2;

		double LineLength = g_GetPoint2Point_Distance(FromPt, ToPt);

		U = ((pt.x - ToPt->x) * (FromPt->x - ToPt->x) + (pt.y - ToPt->y) * (FromPt->y - ToPt->y)) / (LineLength * LineLength);

		if (U < -5 || U > 5)  // larger range
			return 100;

		return 50;
	}



	int g_GetTwoPointSegment2LineIntersectionFlag(const GDPoint* pt0, const GDPoint* pt1, const GDPoint* FromPt, const GDPoint* ToPt, double& L2LDistance)
	{
		double U;
		GDPoint Intersection;
		GDPoint pt;

		L2LDistance = 9999999;

		pt.x = (pt0->x + pt1->x) / 2;
		pt.y = (pt0->y + pt1->y) / 2;

		double LineLength = g_GetPoint2Point_Distance(FromPt, ToPt);

		U = ((pt.x - ToPt->x) * (FromPt->x - ToPt->x) + (pt.y - ToPt->y) * (FromPt->y - ToPt->y)) / (LineLength * LineLength);

		if (U < 0 || U > 1)  // larger range
			return 100;

		Intersection.x = ToPt->x + U * (FromPt->x - ToPt->x);
		Intersection.y = ToPt->y + U * (FromPt->y - ToPt->y);

		double distance_1 = g_GetPoint2Point_Distance(&pt, &Intersection);
		double distance_0 = g_GetPoint2Point_Distance(&pt, FromPt);
		double distance_2 = g_GetPoint2Point_Distance(&pt, ToPt);

		L2LDistance = (distance_0 + distance_1 + distance_2) / 3;

		if (distance_1 > LineLength * 3)  // longer than 3 times of link length, use an exit condition 
		{
			L2LDistance = 100;
		}

		double relative_angle = fabs(g_Find_PPP_RelativeAngle(pt0, pt1, FromPt, ToPt));
		if (relative_angle > 45)
			return 10;

		float matching_ratio = 0.3;

		float ratio_0 = distance_1 / distance_0;
		float ratio_1 = distance_1 / distance_2;
		float ratio_line_length = distance_1 / LineLength;

		if ((ratio_0 < matching_ratio || ratio_1 < matching_ratio)&& ratio_line_length < matching_ratio)
			return 1;  // being matched;
		else
			return 5;  // close but not matched. 

	}

};


class CCSVParser{
public:
    char Delimiter;
    bool IsFirstLineHeader;
    // for DataHub CSV files
    bool m_bSkipFirstLine;
    bool m_bDataHubSingleCSVFile;
    bool m_bLastSectionRead;

    std::ifstream inFile;
    std::string mFileName;
    std::string m_DataHubSectionName;
    std::string SectionName;

    std::vector<std::string> LineFieldsValue;
    std::vector<int> LineIntegerVector;
    std::vector<std::string> Headers;
    std::map<std::string, int> FieldsIndices;

    CCSVParser() : Delimiter{ ',' }, IsFirstLineHeader{ true }, m_bSkipFirstLine{ false }, m_bDataHubSingleCSVFile{ false }, m_bLastSectionRead{ false }
    {
    }

    ~CCSVParser()
    {
        if (inFile.is_open())
            inFile.close();
    }

    // inline member functions
    std::vector<std::string> GetHeaderVector()
    {
        return Headers;
    }
    void CloseCSVFile()
    {
        inFile.close();
    }

    void ConvertLineStringValueToIntegers();
    bool OpenCSVFile(std::string fileName, bool b_required);
    bool ReadRecord();
    bool ReadSectionHeader(std::string s);
    bool ReadRecord_Section();
    std::vector<std::string> ParseLine(std::string line);
    bool GetValueByFieldName(std::string field_name, std::string& value, bool required_field = true);
    template <class T> bool GetValueByFieldName(std::string field_name, T& value, bool required_field = true, bool NonnegativeFlag = true);
};

// Peiheng, 03/22/21, to avoid implicit instantiations in flash_dta.cpp and main_api.cpp for this template function only
// all the other non-inline functions are implemented in utils.cpp
template <class T>
bool CCSVParser::GetValueByFieldName(std::string field_name, T& value, bool required_field, bool NonnegativeFlag)
{
    if (FieldsIndices.find(field_name) == FieldsIndices.end())
    {
        if (required_field)
        {
            dtalog.output() << "Field " << field_name << " in file " << mFileName.c_str() << " does not exist. Please check the file." << std::endl;
            g_ProgramStop();
        }
        return false;
    }
    else
    {
        if (LineFieldsValue.size() == 0)
        {
            return false;
        }

        int size = (int)(LineFieldsValue.size());
        if (FieldsIndices[field_name] >= size)
        {
            return false;
        }

        std::string str_value = LineFieldsValue[FieldsIndices[field_name]];

        if (str_value.length() <= 0)
        {
            return false;
        }

        std::istringstream ss(str_value);

        T converted_value;
        ss >> converted_value;

        if (/*!ss.eof() || */ ss.fail())
        {
            return false;
        }

        if (required_field)
        {
            if(NonnegativeFlag)
            {
                if (converted_value < 0)
                    converted_value = 0;
            }
        }

        value = converted_value;
        return true;
    }
}

#endif