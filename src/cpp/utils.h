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
#define BUILD_EXE //self-use
// if you are using cmake, please #include <build_config.h>
#ifndef _WIN32
#include <build_config.h>
using __int64 = long long;
#endif

#include "teestream.h"

#include <vector>
#include <map>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

constexpr auto _PI = 3.1415926;
// utilities functions

struct DTAGDPoint //geometry data
{
    double x;
    double y;
};
void g_program_stop();
void g_program_exit();
bool g_get_line_polygon_intersection(
    double Ax, double Ay,
    double Bx, double By,
    std::vector<DTAGDPoint> subarea_shape_points);
void g_find_convex_hull(std::vector<DTAGDPoint> points, std::vector<DTAGDPoint> &points_in_polygon);
int g_test_point_in_polygon(DTAGDPoint Pt, std::vector<DTAGDPoint> V);

double g_calculate_p2p_distance_in_meter_from_latitude_longitude(double p1_x, double p1_y, double p2_x, double p2_y);

void fopen_ss(FILE** file, const char* fileName, const char* mode);
float g_read_float(FILE* f);

std::vector<std::string> split(const std::string& s, const std::string& seperator);
int g_ParserIntSequence(std::string str, std::vector<int>& vect);
int g_ParserStringSequence(std::string str, std::vector<std::string>& vect);

std::vector<float> g_time_parser(std::string str);
float g_timestamp_parser(std::string str);
std::string g_time_coding(float time_stamp);
bool g_read_a_line(FILE* f);
// Peiheng, 04/01/21, this is just a temporary fix on logging in DTALite
// it creates another global variable (i.e. dtalog right after class DTALog)
// shared by all translation units, which is really bad. This is a common issue
// across the current implementation. It will be addressed properly in the refactoring.
class DTALog{
    std::ofstream logfile;
    teestream ts;

    int db;
    int sig;
    int odme;
    int path;
    int dta;
    int ue;
public:

    DTALog() : logfile{ "log_main.txt" }, ts{ std::cout, logfile }, db{ 0 }, sig{ 0 }, odme{ 0 }, path{ 0 }, dta{ 0 }, ue{ 0 }
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
    int log_path() const {return path;}

    int& log_dta() {return dta;}
    int log_dta() const {return dta;}

    int& log_ue() {return ue;}
    int log_ue() const {return ue;}
};

static DTALog dtalog;



typedef struct
{
    double X, Y, Z;
}CCoordinate;

class CDTACSVParser{
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

    CDTACSVParser() : Delimiter{ ',' }, IsFirstLineHeader{ true }, m_bSkipFirstLine{ false }, m_bDataHubSingleCSVFile{ false }, m_bLastSectionRead{ false }
    {
    }

    ~CDTACSVParser()
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
bool CDTACSVParser::GetValueByFieldName(std::string field_name, T& value, bool required_field, bool NonnegativeFlag)
{
    if (FieldsIndices.find(field_name) == FieldsIndices.end())
    {
        if (required_field)
        {
            dtalog.output() << "Field " << field_name << " in file " << mFileName.c_str() << " does not exist. Please check the file." << std::endl;
            g_program_stop();
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
