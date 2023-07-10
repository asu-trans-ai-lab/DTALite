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
extern  std::ofstream  g_DTA_log_file;
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
int g_ParserDoubleSequence(std::string str, std::vector<double>& vect);
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


template <typename T>
T* Allocate1DDynamicArray(int nRows)
{
    T* dynamicVector;

    dynamicVector = new (std::nothrow) T[nRows]();

    if (dynamicVector == NULL)
    {
        exit(1);

    }
    return dynamicVector;
}

template <typename T>
void Deallocate1DDynamicArray(T* dVector, int nRows)
{
    if (!dVector)
        return;
    delete[] dVector;
}

template <typename T>
T** Allocate2DDynamicArray(int nRows, int nCols)
{
    T** dynamicArray;

    dynamicArray = new (std::nothrow) T * [nRows];

    if (!dynamicArray)
    {
        dtalog.output() << "[ERROR] insufficient memory.";
        g_DTA_log_file << "[ERROR] insufficient memory.";
        g_program_stop();
    }

    for (int i = 0; i < nRows; ++i)
    {
        dynamicArray[i] = new (std::nothrow) T[nCols];

        if (!dynamicArray[i])
        {
            dtalog.output() << "[ERROR] insufficient memory.";
            g_DTA_log_file << "[ERROR] insufficient memory.";
            g_program_stop();
        }
    }

    return dynamicArray;
}

template <typename T>
void Deallocate2DDynamicArray(T** dArray, int nRows, int nCols)
{
    if (!dArray)
        return;

    for (int x = 0; x < nRows; ++x)
        delete[] dArray[x];

    delete[] dArray;
}

template <typename T>
void Deallocate2DDynamicArray(T** dArray, int nRows)
{
    if (!dArray)
        return;

    for (int x = 0; x < nRows; ++x)
        delete[] dArray[x];

    delete[] dArray;
}
template <typename T>
T*** Allocate3DDynamicArray(int nX, int nY, int nZ)
{
    T*** dynamicArray = new (std::nothrow) T * *[nX];

    if (!dynamicArray)
    {
        dtalog.output() << "[ERROR] insufficient memory.";
        g_DTA_log_file << "[ERROR] insufficient memory.";
        g_program_stop();
    }

    for (int x = 0; x < nX; ++x)
    {
        if (x % 1000 == 0)
        {
            //dtalog.output() << "[DATA INFO] allocating 3D memory for " << x << '\n';
            //g_DTA_log_file << "[DATA INFO] allocating 3D memory for " << x << '\n';
        }

        dynamicArray[x] = new (std::nothrow) T * [nY];

        if (!dynamicArray[x])
        {
            dtalog.output() << "[ERROR] insufficient memory.";
            g_DTA_log_file << "[ERROR] insufficient memory.";
            g_program_stop();
        }

        for (int y = 0; y < nY; ++y)
        {
            dynamicArray[x][y] = new (std::nothrow) T[nZ];
            if (!dynamicArray[x][y])
            {
                dtalog.output() << "[ERROR] insufficient memory.";
                g_DTA_log_file << "[ERROR] insufficient memory.";
                g_program_stop();
            }
        }
    }

    for (int x = 0; x < nX; ++x)
        for (int y = 0; y < nY; ++y)
            for (int z = 0; z < nZ; ++z)
                dynamicArray[x][y][z] = 0;

    return dynamicArray;
}

template <typename T>
void Deallocate3DDynamicArray(T*** dArray, int nX, int nY)
{
    if (!dArray)
        return;

    for (int x = 0; x < nX; ++x)
    {
        for (int y = 0; y < nY; ++y)
            delete[] dArray[x][y];

        delete[] dArray[x];
    }

    delete[] dArray;
}

template <typename T>
T**** Allocate4DDynamicArray(int nM, int nX, int nY, int nZ)
{
    T**** dynamicArray = new (std::nothrow) T * **[nX];

    if (!dynamicArray)
    {
        dtalog.output() << "[ERROR] Insufficient memory.";
        g_DTA_log_file << "[ERROR] Insufficient memory.";
        g_program_stop();
    }

    if (nM == 0 || nX == 0 || nY == 0 || nZ == 0)
    {
        dtalog.output() << "[ERROR] Allocating 4D memory but size = 0 in 1 dimension.";
        g_DTA_log_file << "[ERROR] Allocating 4D memory but size = 0 in 1 dimension.";
        g_program_stop();
    }

    for (int m = 0; m < nM; ++m)
    {
        if (m % 1000 == 0)
        {
            dtalog.output() << "[DATA INFO] Allocating 4D memory for zone index (start from 0) " << m << " with the following dimensions: "
                << "nM = " << nM << ", "
                << "nX = " << nX << ", "
                << "nY = " << nY << ", "
                << "nZ = " << nZ << '\n';
            g_DTA_log_file << "[DATA INFO] Allocating 4D memory for zone index (start from 0) " << m << " with the following dimensions: "
                << "nM = " << nM << ", "
                << "nX = " << nX << ", "
                << "nY = " << nY << ", "
                << "nZ = " << nZ << '\n';
        }

        dynamicArray[m] = new (std::nothrow) T * *[nX];

        if (!dynamicArray[m])
        {
            dtalog.output() << "[ERROR] insufficient memory.";
            g_DTA_log_file << "[ERROR] insufficient memory.";
            g_program_stop();
        }

        for (int x = 0; x < nX; ++x)
        {
            dynamicArray[m][x] = new (std::nothrow) T * [nY];

            if (!dynamicArray[m][x])
            {
                dtalog.output() << "[ERROR] insufficient memory.";
                g_DTA_log_file << "[ERROR] insufficient memory.";
                g_program_stop();
            }

            for (int y = 0; y < nY; ++y)
            {
                dynamicArray[m][x][y] = new (std::nothrow) T[nZ];
                if (!dynamicArray[m][x][y])
                {
                    dtalog.output() << "[ERROR] insufficient memory.";
                    g_DTA_log_file << "[ERROR] insufficient memory.";
                    g_program_stop();
                }
            }
        }
    }

    return dynamicArray;
}

template <typename T>
void Deallocate4DDynamicArray(T**** dArray, int nM, int nX, int nY)
{
    if (!dArray)
        return;

    for (int m = 0; m < nM; ++m)
    {
        for (int x = 0; x < nX; ++x)
        {
            for (int y = 0; y < nY; ++y)
                delete[] dArray[m][x][y];

            delete[] dArray[m][x];
        }

        delete[] dArray[m];
    }

    delete[] dArray;
}



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
    template <class T> bool GetValueByKeyName(std::string field_name, T& value, bool required_field = true, bool NonnegativeFlag = true);
    bool CDTACSVParser::CheckingSettingFormat()
    {
        ReadRecord();
        std::string  field_name = "key";

        if (FieldsIndices.find("key") == FieldsIndices.end())
        {
            dtalog.output() << "[CRITICAL ERROR] The 'key' column in the file '" << mFileName.c_str() << "' cannot be found. Please ensure the settings file adheres to the 'section;key;value' format." << '\n';
            g_DTA_log_file << "[CRITICAL ERROR] The 'key' column in the file '" << mFileName.c_str() << "' cannot be found. Please ensure the settings file adheres to the 'section;key;value' format." << '\n';
            return false;
        }

        if (FieldsIndices.find("value") == FieldsIndices.end())
        {
            dtalog.output() << "[CRITICAL ERROR] The 'value' column in the file '" << mFileName.c_str() << "' cannot be found. Please ensure the settings file adheres to the 'section;key;value' format." << '\n';
            g_DTA_log_file << "[CRITICAL ERROR] The 'value' column in the file '" << mFileName.c_str() << "' cannot be found. Please ensure the settings file adheres to the 'section;key;value' format." << '\n';
            return false;
        }
        else
        {
            return true;
        }
        return false;
    }
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
            dtalog.output() << "[ERROR] Field " << field_name << " in file " << mFileName.c_str() << " does not exist. Please check the file." << '\n';
            g_DTA_log_file << "[ERROR] Field " << field_name << " in file " << mFileName.c_str() << " does not exist. Please check the file." << '\n';
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

        //if (required_field)
        //{
        //    if(NonnegativeFlag)
        //    {
        //        if (converted_value < 0)
        //            converted_value = 0;
        //    }
        //}

        value = converted_value;
        return true;
    }
}

template <class T>
bool CDTACSVParser::GetValueByKeyName(std::string key_name, T& value, bool required_field, bool NonnegativeFlag)
{
    if (inFile.is_open())
        inFile.close();

    OpenCSVFile(mFileName,false);
    ReadRecord();
 
        do
        {
            std::string  key_name_record = "key";
            GetValueByFieldName("key", key_name_record);

            if (key_name_record == key_name)
            {
                GetValueByFieldName("value", value);
                return true;
            }
        } while (ReadRecord());
    return false;
}


#endif