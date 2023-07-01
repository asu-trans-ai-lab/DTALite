/* Portions Copyright 2021 Xuesong Zhou, Peiheng Li, Cafer Avci
 *
 * If you help write or modify the code, please also list your names here.
 * The reason of having Copyright info here is to ensure all the modified version, as a whole, under the GPL
 * and further prevent a violation of the GPL.
 *
 * More about "How to use GNU licenses for your own software"
 * http://www.gnu.org/licenses/gpl-howto.html
 */

#include "teestream.h"
#include "utils.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cstring>
#include <vector>
#include <cmath>
#include <stack>

using std::string;
using std::vector;
using std::ofstream;
using std::istringstream;
using std::ostringstream;

using std::asin;
using std::sin;
using std::cos;
using std::pow;
using std::sqrt;
using std::min;
using std::fmin;
#include "shared_code.h"

void g_program_stop()
{
    dtalog.output() << "[ERROR] DTALite Program stops!" << std::endl;
    exit(1);
}

void g_program_exit()
{
    dtalog.output() << "[STATUS INFO] DTALite Program completes!" << std::endl;
    exit(0);
}

void fopen_ss(FILE** file, const char* fileName, const char* mode)
{
    *file = fopen(fileName, mode);
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


float g_timestamp_parser(string str)
{
   float output_global_minute;

    int string_lenghth = str.length();

    //ASSERT(string_lenghth < 100);

    const char* string_line = str.data(); //string to char*

    int char_length = strlen(string_line);

    char ch, buf_ddhhmm[32] = { 0 }, buf_SS[32] = { 0 }, buf_sss[32] = { 0 };
    char dd1, dd2, hh1, hh2, mm1, mm2, SS1, SS2, sss1, sss2, sss3;
    float ddf1, ddf2, hhf1, hhf2, mmf1, mmf2, SSf1, SSf2, sssf1, sssf2, sssf3;
    float global_minute = 0;
    float dd = 0, hh = 0, mm = 0, SS = 0, sss = 0;
    int i = 1;  // skip T as the first letter
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

        if (i == char_length) //start a new time string
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

            output_global_minute = global_minute;

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






// definitions of CDTACSVParser member functions
void CDTACSVParser::ConvertLineStringValueToIntegers()
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

bool CDTACSVParser::OpenCSVFile(string fileName, bool b_required)
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
            dtalog.output() << "[WARNING] File " << fileName << " does not exist. Please check." << '\n';
            //g_program_stop();
        }
        return false;
    }
}

bool CDTACSVParser::ReadRecord()
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

bool CDTACSVParser::ReadSectionHeader(string s)
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

bool CDTACSVParser::ReadRecord_Section()
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

vector<string> CDTACSVParser::ParseLine(string line)
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

bool CDTACSVParser::GetValueByFieldName(string field_name, string& value, bool required_field)
{
    if (FieldsIndices.find(field_name) == FieldsIndices.end())
    {
        if (required_field)
        {
            dtalog.output() << "[ERROR] Field " << field_name << " in file " << mFileName << " does not exist. Please check the file." << '\n';
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

bool g_read_a_line(FILE* f)
/* read a line from the current line from the file */
{

    char ch;

    while (1) {
        ch = getc(f);
        if (ch != 13 && ch != 10 && ch != EOF)
        {
            // do nothing
        }
        else { /* terminate if it's end of line or end of file */
            {
                // do nothing
            }
            if (ch == EOF)
                return false;

            return true;
        }
    }
}

double g_calculate_p2p_distance_in_meter_from_latitude_longitude(double longitud1, double latitud1, double longitud2, double latitud2)
{
    double PI = 3.14159265358979323846;
    double RADIO_TERRESTRE = 6372797.56085;
    double GRADOS_RADIANES = PI / 180;

    double haversine;
    double temp;
    double distancia_puntos;

    latitud1 = latitud1 * GRADOS_RADIANES;
    longitud1 = longitud1 * GRADOS_RADIANES;
    latitud2 = latitud2 * GRADOS_RADIANES;
    longitud2 = longitud2 * GRADOS_RADIANES;

    haversine = (pow(sin((1.0 / 2) * (latitud2 - latitud1)), 2)) + ((cos(latitud1)) * (cos(latitud2)) * (pow(sin((1.0 / 2) * (longitud2 - longitud1)), 2)));
    temp = 2 * asin(fmin(1.0, sqrt(haversine)));
    distancia_puntos = RADIO_TERRESTRE * temp;

    return distancia_puntos;
}




/* ***************************************************************************** */
/* Copyright:      Francois Panneton and Pierre L'Ecuyer, University of Montreal */
/*                 Makoto Matsumoto, Hiroshima University                        */
/* Notice:         This code can be used freely for personal, academic,          */
/*                 or non-commercial purposes. For commercial purposes,          */
/*                 please contact P. L'Ecuyer at: lecuyer@iro.UMontreal.ca       */
/* ***************************************************************************** */
#define W 32
#define R 16
#define P 0
#define M1 13
#define M2 9
#define M3 5

#define MAT0POS(t,v) (v^(v>>t))
#define MAT0NEG(t,v) (v^(v<<(-(t))))
#define MAT3NEG(t,v) (v<<(-(t)))
#define MAT4NEG(t,b,v) (v ^ ((v<<(-(t))) & b))

#define V0            STATE[state_i                   ]
#define VM1           STATE[(state_i+M1) & 0x0000000fU]
#define VM2           STATE[(state_i+M2) & 0x0000000fU]
#define VM3           STATE[(state_i+M3) & 0x0000000fU]
#define VRm1          STATE[(state_i+15) & 0x0000000fU]
#define VRm2          STATE[(state_i+14) & 0x0000000fU]
#define newV0         STATE[(state_i+15) & 0x0000000fU]
#define newV1         STATE[state_i                 ]
#define newVRm1       STATE[(state_i+14) & 0x0000000fU]

#define FACT 2.32830643653869628906e-10

static unsigned int state_i = 0;
static unsigned int STATE[R];
static unsigned int z0, z1, z2;
unsigned int g_RandomSeed = 100;
void InitWELLRNG512a(unsigned int* init) {
    int j;
    state_i = 0;
    for (j = 0; j < R; j++)
        STATE[j] = init[j];
}

double WELLRNG512a(void) {
    z0 = VRm1;
    z1 = MAT0NEG(-16, V0) ^ MAT0NEG(-15, VM1);
    z2 = MAT0POS(11, VM2);
    newV1 = z1 ^ z2;
    newV0 = MAT0NEG(-2, z0) ^ MAT0NEG(-18, z1) ^ MAT3NEG(-28, z2) ^ MAT4NEG(-5, 0xda442d24U, newV1);
    state_i = (state_i + 15) & 0x0000000fU;
    return ((double)STATE[state_i]) * FACT;
}



double g_get_random_ratio()
{
    //	g_RandomSeed = (g_LCG_a * g_RandomSeed + g_LCG_c) % g_LCG_M;  //m_RandomSeed is automatically updated.
    //	return float(g_RandomSeed)/g_LCG_M;

    return WELLRNG512a();
}


//  public domain function by Darel Rex Finley, 2006
//  Determines the intersection point of the line defined by points A and B with the
//  line defined by points C and D.
//
//  Returns YES if the intersection point was found, and stores that point in X,Y.
//  Returns NO if there is no determinable intersection point, in which case X,Y will
//  be unmodified.

bool g_get_line_intersection(
    double Ax, double Ay,
    double Bx, double By,
    double Cx, double Cy,
    double Dx, double Dy)
{
    float X = 0;
    float Y = 0;

    double  distAB, theCos, theSin, newX, ABpos;

    //  Fail if either line segment is zero-length.
  //  if (Ax==Bx && Ay==By || Cx==Dx && Cy==Dy) return false;
    if (Ax == Bx && Ay == By) return false;  // comment: C and D can be the same point from a vehile with the same timestamp

    //  Fail if the segments share an end-point.
    if (Ax == Cx && Ay == Cy || Bx == Cx && By == Cy
        || Ax == Dx && Ay == Dy || Bx == Dx && By == Dy) {
        return false;
    }

    //  (1) Translate the system so that point A is on the origin.
    Bx -= Ax; By -= Ay;
    Cx -= Ax; Cy -= Ay;
    Dx -= Ax; Dy -= Ay;

    //  Discover the length of segment A-B.
    distAB = sqrt(Bx * Bx + By * By);

    //  (2) Rotate the system so that point B is on the positive X axis.
    theCos = Bx / distAB;
    theSin = By / distAB;
    newX = Cx * theCos + Cy * theSin;
    Cy = Cy * theCos - Cx * theSin; Cx = newX;
    newX = Dx * theCos + Dy * theSin;
    Dy = Dy * theCos - Dx * theSin; Dx = newX;

    //  Fail if segment C-D doesn't cross line A-B.
    if (Cy < 0. && Dy < 0. || Cy >= 0. && Dy >= 0.) return false;

    //  (3) Discover the position of the intersection point along line A-B.
    ABpos = Dx + (Cx - Dx) * Dy / (Dy - Cy);

    //  Fail if segment C-D crosses line A-B outside of segment A-B.
    if (ABpos<0. || ABpos>distAB) return false;

    //  (4) Apply the discovered position to line A-B in the original coordinate system.
    X = Ax + ABpos * theCos;
    Y = Ay + ABpos * theSin;

    //  Success.
    return true;
}

bool g_get_line_polygon_intersection(
    double Ax, double Ay,
    double Bx, double By,
    std::vector<DTAGDPoint> subarea_shape_points)
{

    double Cx, Cy;
    double Dx, Dy;
    for (int i = 0; i < subarea_shape_points.size()-1; i++)
    {
        Cx = subarea_shape_points[i].x;
        Cy = subarea_shape_points[i].y;

        Dx = subarea_shape_points[i+1].x;
        Dy = subarea_shape_points[i+1].y;

        if (g_get_line_intersection(Ax, Ay, Bx, By, Cx, Cy, Dx, Dy) == true)
            return true;
    }

    return false;
}

int g_test_point_in_polygon(DTAGDPoint Pt, std::vector<DTAGDPoint> V)
{
    int n = V.size()-1;
    int    cn = 0;    // the  crossing number counter

    // loop through all edges of the polygon
    for (int i = 0; i < n; i++)
    {    // edge from V[i]  to V[i+1]
        if (((V[i].y <= Pt.y) && (V[i+1].y > Pt.y))     // an upward crossing
            || ((V[i].y > Pt.y) && (V[i+1].y <= Pt.y)))
        { // a downward crossing
                // compute  the actual edge-ray intersect x-coordinate
            float vt = (float)(Pt.y - V[i].y) / (V[i+1].y - V[i].y);
            if (Pt.x < V[i].x + vt * (V[i+1].x - V[i].x)) // P.x < intersect
                ++cn;   // a valid crossing of y=P.y right of P.x
        }
    }
    return (cn & 1);    // 0 if even (out), and 1 if  odd (in)
}


// A globle point needed for  sorting points with reference
// to  the first point Used in compare function of qsort()
// reference: https://iq.opengenus.org/graham-scan-convex-hull/

DTAGDPoint p0;
// A utility function to find next to top in a stack

DTAGDPoint nextToTop(std::stack<DTAGDPoint>& S)
{
    DTAGDPoint p = S.top();
    S.pop();
    DTAGDPoint res = S.top();
    S.push(p);
    return res;
}
// A utility function to swap two points
int swap(DTAGDPoint& p1, DTAGDPoint& p2)
{
    DTAGDPoint temp = p1;
    p1 = p2;
    p2 = temp;
    return 1;
}
// A utility function to return square of distance
// between p1 and p2
double distSq(DTAGDPoint p1, DTAGDPoint p2)
{
    return (p1.x - p2.x) * (p1.x - p2.x) +
        (p1.y - p2.y) * (p1.y - p2.y);
}
// To find orientation of ordered triplet (p, q, r).
// The function returns following values
// 0 --> p, q and r are colinear
// 1 --> Clockwise
// 2 --> Counterclockwise
int orientation(DTAGDPoint p, DTAGDPoint q, DTAGDPoint r)
{
    double val = (q.y - p.y) * (r.x - q.x) -
        (q.x - p.x) * (r.y - q.y);
    if (fabs(val) <= 0.0001) return 0;  // colinear
    return (val > 0) ? 1 : 2; // clock or counterclock wise
}
// A function used by library function qsort() to sort an array of
// points with respect to the first DTAGDPoint
int compare(const void* vp1, const void* vp2)
{
    DTAGDPoint* p1 = (DTAGDPoint*)vp1;
    DTAGDPoint* p2 = (DTAGDPoint*)vp2;
    // Find orientation
    int o = orientation(p0, *p1, *p2);
    if (o == 0)
        return (distSq(p0, *p2) >= distSq(p0, *p1)) ? -1 : 1;
    return (o == 2) ? -1 : 1;
}
// Prints convex hull of a set of n points.
void g_find_convex_hull(std::vector<DTAGDPoint> points, std::vector<DTAGDPoint> &points_in_polygon)
{
    int n = points.size();
    // Find the bottommost DTAGDPoint
    int ymin = points[0].y, min = 0;
    for (int i = 1; i < n; i++)
    {
        int y = points[i].y;
        // Pick the bottom-most or chose the left
        // most DTAGDPoint in case of tie
        if ((y < ymin) || (ymin == y &&
            points[i].x < points[min].x))
            ymin = points[i].y, min = i;
    }
    // Place the bottom-most DTAGDPoint at first position
    swap(points[0], points[min]);
    // Sort n-1 points with respect to the first DTAGDPoint.
    // A DTAGDPoint p1 comes before p2 in sorted ouput if p2
    // has larger polar angle (in counterclockwise
    // direction) than p1
    p0 = points[0];
    qsort(&points[1], n - 1, sizeof(DTAGDPoint), compare);
    // If two or more points make same angle with p0,
    // Remove all but the one that is farthest from p0
    // Remember that, in above sorting, our criteria was
    // to keep the farthest DTAGDPoint at the end when more than
    // one points have same angle.
    int m = 1; // Initialize size of modified array
    for (int i = 1; i < n; i++)
    {
        // Keep removing i while angle of i and i+1 is same
        // with respect to p0
        while (i < n - 1 && orientation(p0, points[i],
            points[i + 1]) == 0)
            i++;
        points[m] = points[i];
        m++;  // Update size of modified array
    }
    // If modified array of points has less than 3 points,
    // convex hull is not possible
    if (m < 3) return;
    // Create an empty stack and push first three points
    // to it.
    std::stack<DTAGDPoint> S;
    S.push(points[0]);
    S.push(points[1]);
    S.push(points[2]);
    // Process remaining n-3 points
    for (int i = 3; i < m; i++)
    {
        // Keep removing top while the angle formed by
        // points next-to-top, top, and points[i] makes
        // a non-left turn
        while (S.size() >= 2 && orientation(nextToTop(S), S.top(), points[i]) != 2)
        {
               S.pop();
        }
        S.push(points[i]);
    }
    // Now stack has the output points, print contents of stack
    while (!S.empty())
    {
        DTAGDPoint p = S.top();
//        cout << "(" << p.x << ", " << p.y << ")" << '\n';
        points_in_polygon.push_back(p);
        S.pop();
    }

    if(points_in_polygon.size() > 0)
    {
    points_in_polygon.push_back(points_in_polygon[0]);
    }
}

double g_Find_P2P_Angle(const DTAGDPoint* p1, const DTAGDPoint* p2)
{
    double PI = 3.14159265358979323846;
    double delta_x = p2->x - p1->x;
    double delta_y = p2->y - p1->y;

    if (fabs(delta_x) < 0.00001)
        delta_x = 0;

    if (fabs(delta_y) < 0.00001)
        delta_y = 0;

    int angle = atan2(delta_y, delta_x) * 180 / PI + 0.5;
    // angle = 90 - angle;

    while (angle < 0)
        angle += 360;

    while (angle > 360)
        angle -= 360;

    return angle;
}

double g_Find_PPP_RelativeAngle(const DTAGDPoint* p1, const DTAGDPoint* p2, const DTAGDPoint* p3, const DTAGDPoint* p4)
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

double g_GetPoint2Point_Distance(const DTAGDPoint* p1, const DTAGDPoint* p2)
{
    return pow(((p1->x - p2->x) * (p1->x - p2->x) + (p1->y - p2->y) * (p1->y - p2->y)), 0.5);
}
double g_GetPoint2LineDistance(const DTAGDPoint* pt, const DTAGDPoint* FromPt, const DTAGDPoint* ToPt)
{
    double U;
    DTAGDPoint Intersection;

    double LineLength = g_GetPoint2Point_Distance(FromPt, ToPt);

    U = ((pt->x - ToPt->x) * (FromPt->x - ToPt->x) + (pt->y - ToPt->y) * (FromPt->y - ToPt->y)) / (LineLength * LineLength);


    Intersection.x = ToPt->x + U * (FromPt->x - ToPt->x);
    Intersection.y = ToPt->y + U * (FromPt->y - ToPt->y);

    double distance_1 = g_GetPoint2Point_Distance(pt, &Intersection);
    double distance_0 = g_GetPoint2Point_Distance(pt, FromPt);
    double distance_2 = g_GetPoint2Point_Distance(pt, ToPt);

    if (U < 0.0 || U > 1.0)
        return min(distance_0, distance_2); // intersection does not fall within the segment
    else  // intersect
        return distance_1;
}