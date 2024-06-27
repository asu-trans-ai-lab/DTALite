#pragma once

#include <string>
#include <vector>

struct DTAGDPoint;

class CDTAGeometry {
public:
    enum GeometryType { POINT, LINE, POLYGON, UNKNOWN };

public:
    CDTAGeometry(std::string s);
    ~CDTAGeometry();

    GeometryType GetGeometryType();
    std::vector<CCoordinate> GetCoordinateList();
    int GetNumberOfCoordinates();

private:
    GeometryType m_Type;
    int m_NumOfCoordinates;
    std::vector<CCoordinate> v_Coordinates;
    bool ReadPointCoordinate(std::string s);
    bool ReadLineStringCoordinates(std::string s);
    bool ReadPolygonCoordinates(std::string s);
};