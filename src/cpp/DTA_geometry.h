#pragma once
#include <string>
#include <vector>

using std::string;

struct DTAGDPoint;


class CDTAGeometry
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
	bool ReadPointCoordinate(string s);
	bool ReadLineStringCoordinates(string s);
	bool ReadPolygonCoordinates(string s);

public:
	CDTAGeometry(string s);
	~CDTAGeometry(void);

	GeometryType GetGeometryType(void);
	std::vector<CCoordinate> GetCoordinateList(void);
	int GetNumberOfCoordinates(void);
};




 