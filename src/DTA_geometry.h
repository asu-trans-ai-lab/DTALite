#pragma once
#include <string>
#include <vector>

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
	bool ReadPointCoordinate(std::string s);
	bool ReadLineStringCoordinates(std::string s);
	bool ReadPolygonCoordinates(std::string s);

public:
	CDTAGeometry(std::string s);
	~CDTAGeometry(void);

	GeometryType GetGeometryType(void);
	std::vector<CCoordinate> GetCoordinateList(void);
	int GetNumberOfCoordinates(void);
};