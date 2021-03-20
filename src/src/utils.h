#ifndef GUARD_UTILS_H
#define GUARD_UTILS_H

#include <vector>
#include <map>
#include <string>
#include <iostream>

using std::vector;
using std::map;
using std::string;
using std::ifstream;

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

    CCSVParser(void)
    {
        Delimiter = ',';
        IsFirstLineHeader = true;
        m_bSkipFirstLine = false;
        m_bDataHubSingleCSVFile = false;
        m_bLastSectionRead = false;
    }

    ~CCSVParser(void)
    {
        if (inFile.is_open()) inFile.close();
    }

    // inline member functions
    std::vector<std::string> GetHeaderVector() 
    {
        return Headers;
    }

    void CloseCSVFile()
    {
        inFile.close();
    };

    void ConvertLineStringValueToIntegers();
    bool OpenCSVFile(std::string fileName, bool b_required);
    bool ReadRecord();
    bool ReadSectionHeader(std::string s);
    bool ReadRecord_Section();
    std::vector<std::string> ParseLine(std::string line);
    template <class T> bool GetValueByFieldName(std::string field_name, T& value, bool required_field = true, bool NonnegativeFlag = true);
    bool GetValueByFieldName(std::string field_name, std::string& value, bool required_field =true);
};

void g_ProgramStop();

#endif