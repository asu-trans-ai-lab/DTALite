#include "utils.h"
#include "DTA_geometry.h"

#include <sstream>
using std::istringstream;
using std::string;

CDTAGeometry::CDTAGeometry(string s)
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
	else if (s.find("MULTILINESTRING ") != std::string::npos)
	{
		tmp = s.substr(s.find_first_not_of(' '));
		size_t start_idx = tmp.find_first_of('((');
		size_t end_idx = tmp.find_first_of('))');

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


		m_Type = LINE;
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
	}else if (s.find("MULTIPOLYGON") != std::string::npos)
	{
		tmp = s.substr(s.find_first_not_of(' '));
		size_t start_idx = tmp.find('((');
		size_t end_idx = tmp.find('))');

		if (start_idx == std::string::npos || end_idx == std::string::npos)
			return;

		string type_str = tmp.substr(0, start_idx);
		type_str.erase(type_str.find_last_not_of(" ") + 1);		// works for 'LINESTRING (....' and 'LINESTRING(....'

		string start_tag = "(((";
		string end_tag = ")))";

		start_idx = tmp.find(start_tag);
		start_idx += start_tag.length();
		end_idx = tmp.find(end_tag);

		tmp = tmp.substr(start_idx, end_idx - start_idx);


		m_Type = POLYGON;
	}else if (s.find("POLYGON") != std::string::npos)
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
	}else
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

CDTAGeometry::~CDTAGeometry(void)
{
}

CDTAGeometry::GeometryType CDTAGeometry::GetGeometryType(void)
{
	return m_Type;
}

int CDTAGeometry::GetNumberOfCoordinates(void)
{
	return m_NumOfCoordinates;
}

std::vector<CCoordinate> CDTAGeometry::GetCoordinateList(void)
{
	return v_Coordinates;
}

bool CDTAGeometry::ReadLineStringCoordinates(string s)
{
	istringstream ss(s);
	string sub_str;

	if(std::string::npos == s.find_first_of("0123456789"))
	{
			// "digit not found!, empty string//
		return false;
	}

	while(std::getline(ss,sub_str, ','))
	{
		sub_str = sub_str.substr(sub_str.find_first_not_of(' '));

		CCoordinate coordinate;
		istringstream sub_ss(sub_str);
		string tmp;

		std::getline(sub_ss,tmp,' ');
		istringstream x_ss(tmp);
		x_ss >> coordinate.X;

		std::getline(sub_ss,tmp,' ');
		istringstream y_ss(tmp);
		y_ss >> coordinate.Y;

		v_Coordinates.push_back(coordinate);
		m_NumOfCoordinates += 1;
	}
	return true;
}

bool CDTAGeometry::ReadPolygonCoordinates(string s)
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

bool CDTAGeometry::ReadPointCoordinate(string s)
{
	CCoordinate coordinate;
	istringstream ss(s);

	string sub_str;
	std::getline(ss,sub_str,' ');
	istringstream x_ss(sub_str);

	std::getline(ss,sub_str,' ');
	istringstream y_ss(sub_str);
	x_ss >> coordinate.X;
	y_ss >> coordinate.Y;
	coordinate.Z = 0.0;

	v_Coordinates.push_back(coordinate);
	m_NumOfCoordinates = 1;

	return true;
}