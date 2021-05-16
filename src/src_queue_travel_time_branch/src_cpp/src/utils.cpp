/* Portions Copyright 2021 Xuesong Zhou and Peiheng Li
 *
 * If you help write or modify the code, please also list your names here.
 * The reason of having Copyright info here is to ensure all the modified version, as a whole, under the GPL
 * and further prevent a violation of the GPL.
 *
 * More about "How to use GNU licenses for your own software"
 * http://www.gnu.org/licenses/gpl-howto.html
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cstring>
#include <vector>
#include "teestream.h"
#include "utils.h"

using std::endl;
using std::string;
using std::vector;
using std::ofstream;
using std::istringstream;
using std::ostringstream;

void g_ProgramStop()
{
    dtalog.output() << "STALite Program stops. Press any key to terminate. Thanks!" << endl;
    getchar();
    exit(0);
}

void fopen_ss(FILE** file, const char* fileName, const char* mode)
{
    *file = fopen(fileName, mode);
}

float g_read_float(FILE* f)
{
    /*
        read a floating point number from the current pointer of the file,
        skip all spaces
     */
    char ch, buf[32];
    int i = 0;
    int flag = 1;

    /* returns -1 if end of file is reached */
    while (true)
    {
        ch = getc(f);
        if (ch == EOF || ch == '*' || ch == '$') return -1;
        if (isdigit(ch))
            break;

        if (ch == '-')
            flag = -1;
        else
            flag = 1;
    }

    if (ch == EOF) return -1;
    while (isdigit(ch) || ch == '.') {
        buf[i++] = ch;
        ch = fgetc(f);

    }
    buf[i] = 0;

    /* atof function converts a character string (char *) into a doubleing
    pointer equivalent, and if the string is not a floting point number,
    a zero will be return.
    */

    return (float)(atof(buf) * flag);
}

//split the string by "_"
vector<string> split(const string &s, const string &seperator)
{
    vector<string> result;
    typedef string::size_type string_size;
    string_size i = 0;

    while (i != s.size()) {
        int flag = 0;
        while (i != s.size() && flag == 0) {
            flag = 1;
            for (string_size x = 0; x < seperator.size(); ++x)
                if (s[i] == seperator[x]) {
                    ++i;
                    flag = 0;
                    break;
                }
        }

        flag = 0;
        string_size j = i;
        while (j != s.size() && flag == 0) {
            for (string_size x = 0; x < seperator.size(); ++x)
                if (s[j] == seperator[x]) {
                    flag = 1;
                    break;
                }
            if (flag == 0)
                ++j;
        }
        if (i != j) {
            result.push_back(s.substr(i, j - i));
            i = j;
        }
    }

    return result;
}

vector<float> g_time_parser(string str)
{
    vector<float> output_global_minute;

    int string_lenghth = str.length();

    //ASSERT(string_lenghth < 100);

    const char* string_line = str.data(); //string to char*

    int char_length = strlen(string_line);

    char ch, buf_ddhhmm[32] = { 0 }, buf_SS[32] = { 0 }, buf_sss[32] = { 0 };
    char dd1, dd2, hh1, hh2, mm1, mm2, SS1, SS2, sss1, sss2, sss3;
    float ddf1, ddf2, hhf1, hhf2, mmf1, mmf2, SSf1, SSf2, sssf1, sssf2, sssf3;
    float global_minute = 0;
    float dd = 0, hh = 0, mm = 0, SS = 0, sss = 0;
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

        if (ch == '_' || i == char_length) //start a new time string
        {
            if (buffer_i == 4) //"HHMM"
            {
                //HHMM, 0123
                hh1 = buf_ddhhmm[0]; //read each first
                hh2 = buf_ddhhmm[1];
                mm1 = buf_ddhhmm[2];
                mm2 = buf_ddhhmm[3];

                hhf1 = ((float)hh1 - 48); //convert a char to a float
                hhf2 = ((float)hh2 - 48);
                mmf1 = ((float)mm1 - 48);
                mmf2 = ((float)mm2 - 48);

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

                ddf1 = ((float)dd1 - 48); //convert a char to a float
                ddf2 = ((float)dd2 - 48);
                hhf1 = ((float)hh1 - 48);
                hhf2 = ((float)hh2 - 48);
                mmf1 = ((float)mm1 - 48);
                mmf2 = ((float)mm2 - 48);

                dd = ddf1 * 10 * 24 * 60 + ddf2 * 24 * 60;
                hh = hhf1 * 10 * 60 + hhf2 * 60;
                mm = mmf1 * 10 + mmf2;
            }

            if (num_of_colons == 1 || num_of_colons == 2)
            {
                //SS, 01
                SS1 = buf_SS[0]; //read each first
                SS2 = buf_SS[1];

                SSf1 = ((float)SS1 - 48); //convert a char to a float
                SSf2 = ((float)SS2 - 48);

                SS = (SSf1 * 10 + SSf2) / 60;
            }

            if (num_of_colons == 2)
            {
                //sss, 012
                sss1 = buf_sss[0]; //read each first
                sss2 = buf_sss[1];
                sss3 = buf_sss[2];

                sssf1 = ((float)sss1 - 48); //convert a char to a float
                sssf2 = ((float)sss2 - 48);
                sssf3 = ((float)sss3 - 48);

                sss = (sssf1 * 100 + sssf2 * 10 + sssf3) / 1000;
            }

            global_minute = dd + hh + mm + SS + sss;

            output_global_minute.push_back(global_minute);

            //initialize the parameters
            buffer_i = 0;
            buffer_k = 0;
            buffer_j = 0;
            num_of_colons = 0;
        }

        if (ch == ':')
            num_of_colons += 1;
    }

    return output_global_minute;
}

string g_time_coding(float time_stamp)
{
    int hour = static_cast<int>(time_stamp / 60);
    int minute = static_cast<int>(time_stamp - hour * 60);
    int second = static_cast<int>((time_stamp - hour * 60 - minute) * 60);

    ostringstream strm;
    strm.fill('0');
    strm << std::setw(2) << hour << std::setw(2) << minute /*<< ":" << setw(2) << second*/;

    return strm.str();
}

int g_ParserIntSequence(std::string str, vector<int>& vect)
{
    std::stringstream ss(str);
    int i;

    while (ss >> i)
    {
        vect.push_back(i);
        if (ss.peek() == ';')
            ss.ignore();
    }

    return vect.size();
}

// definitions of CCSVParser member functions
void CCSVParser::ConvertLineStringValueToIntegers()
{
    LineIntegerVector.clear();
    for (unsigned i = 0; i < LineFieldsValue.size(); ++i)
    {
        string si = LineFieldsValue[i];
        int value = atoi(si.c_str());

        if (value >= 1)
            LineIntegerVector.push_back(value);
    }
}

bool CCSVParser::OpenCSVFile(string fileName, bool b_required)
{
    mFileName = fileName;
    inFile.open(fileName.c_str());

    if (inFile.is_open())
    {
        if (IsFirstLineHeader)
        {
            string s;
            std::getline(inFile, s);
            vector<string> FieldNames = ParseLine(s);

            for (size_t i = 0;i < FieldNames.size();i++)
            {
                string tmp_str = FieldNames.at(i);
                size_t start = tmp_str.find_first_not_of(" ");

                string name;
                if (start == string::npos)
                {
                    name = "";
                }
                else
                {
                    name = tmp_str.substr(start);
                    //TRACE("%s,", name.c_str());
                }
                FieldsIndices[name] = (int)i;
            }
        }
        return true;
    }
    else
    {
        if (b_required)
        {
            dtalog.output() << "File " << fileName << " does not exist. Please check." << std::endl;
            //g_ProgramStop();
        }
        return false;
    }
}

bool CCSVParser::ReadRecord()
{
    LineFieldsValue.clear();

    if (inFile.is_open())
    {
        string s;
        std::getline(inFile, s);
        if (s.length() > 0)
        {
            LineFieldsValue = ParseLine(s);
            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        return false;
    }
}

bool CCSVParser::ReadSectionHeader(string s)
{
    //skip // data
    Headers.clear();
    FieldsIndices.clear();

    if (s.length() == 0)
        return true;

    vector<string> FieldNames = ParseLine(s);

    for (size_t i = 0; i < FieldNames.size(); i++)
    {
        string tmp_str = FieldNames.at(i);
        size_t start = tmp_str.find_first_not_of(" ");

        string name;
        if (start == string::npos)
        {
            name = "";
        }
        else
        {
            name = tmp_str.substr(start);
        }
        Headers.push_back(name);
        FieldsIndices[name] = (int)i;
    }
    return true;
}

bool CCSVParser::ReadRecord_Section()
{
    LineFieldsValue.clear();

    if (inFile.is_open())
    {
        string s;
        std::getline(inFile, s);
        if (s.length() > 0)
        {
            if(s.find("[") != string::npos)  // synchro single csv file
            {
                LineFieldsValue = ParseLine(s);

                if (LineFieldsValue.size() >= 1)
                {
                    SectionName = LineFieldsValue[0];
                }

                //re-read section header
                ReadSectionHeader(s);
                std::getline(inFile, s);
            }
            LineFieldsValue = ParseLine(s);
            return true;
        }
        else
        {
            if (m_bLastSectionRead)  // reach the last section
                return false;
            else
            {
                if (inFile.eof())
                    return false;
                else
                    return true;
            }
        }
    }
    else
    {
        return false;
    }
}

vector<string> CCSVParser::ParseLine(string line)
{
    vector<string> SeperatedStrings;
    string subStr;

    if (line.length() == 0)
        return SeperatedStrings;

    std::istringstream ss(line);

    if (line.find_first_of('"') == string::npos)
    {
        while (std::getline(ss, subStr, Delimiter))
        {
            SeperatedStrings.push_back(subStr);
        }

        if (line.at(line.length() - 1) == ',')
        {
            SeperatedStrings.push_back("");
        }
    }
    else
    {
        while (line.length() > 0)
        {
            size_t n1 = line.find_first_of(',');
            size_t n2 = line.find_first_of('"');

            if (n1 == string::npos && n2 == string::npos) //last field without double quotes
            {
                subStr = line;
                SeperatedStrings.push_back(subStr);
                break;
            }

            if (n1 == string::npos && n2 != string::npos) //last field with double quotes
            {
                size_t n3 = line.find_first_of('"', n2 + 1); // second double quote

                //extract content from double quotes
                subStr = line.substr(n2 + 1, n3 - n2 - 1);
                SeperatedStrings.push_back(subStr);

                break;
            }

            if (n1 != string::npos && (n1 < n2 || n2 == string::npos))
            {
                subStr = line.substr(0, n1);
                SeperatedStrings.push_back(subStr);
                if (n1 < line.length() - 1)
                {
                    line = line.substr(n1 + 1);
                }
                else // comma is the last char in the line string, push an empty string to the back of vector
                {
                    SeperatedStrings.push_back("");
                    break;
                }
            }

            if (n1 != string::npos && n2 != string::npos && n2 < n1)
            {
                size_t n3 = line.find_first_of('"', n2 + 1); // second double quote
                subStr = line.substr(n2 + 1, n3 - n2 - 1);
                SeperatedStrings.push_back(subStr);
                size_t idx = line.find_first_of(',', n3 + 1);

                if (idx != string::npos)
                {
                    line = line.substr(idx + 1);
                }
                else
                {
                    break;
                }
            }
        }
    }
    return SeperatedStrings;
}

bool CCSVParser::GetValueByFieldName(string field_name, string& value, bool required_field)
{
    if (FieldsIndices.find(field_name) == FieldsIndices.end())
    {
        if (required_field)
        {
            dtalog.output() << "Field " << field_name << " in file " << mFileName << " does not exist. Please check the file." << std::endl;
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

        unsigned int index = FieldsIndices[field_name];
        if (index >= LineFieldsValue.size())
        {
            return false;
        }
        string str_value = LineFieldsValue[index];

        if (str_value.length() <= 0)
        {
            return false;
        }

        value = str_value;
        return true;
    }
}