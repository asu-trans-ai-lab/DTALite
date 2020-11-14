//  Portions Copyright 2019
// Xuesong (Simon) Zhou
//   If you help write or modify the code, please also list your names here.
//   The reason of having Copyright info here is to ensure all the modified version, as a whole, under the GPL 
//   and further prevent a violation of the GPL.

// More about "How to use GNU licenses for your own software"
// http://www.gnu.org/licenses/gpl-howto.html

#pragma warning( disable : 4305 4267 4018) 
#include <iostream>
#include <fstream>
#include <list> 
#include <omp.h>
#include <algorithm>
#include <time.h>
#include <functional>
#include <stdio.h>   
#include <math.h>


#include <stack>
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <iostream>
#include <iomanip>
using namespace std;
using std::string;
using std::ifstream;
using std::vector;
using std::map;
using std::istringstream;
using std::max;

string g_info_String = "";

class CMainSigModual {
public:
	CMainSigModual()
	{
		b_with_loaded_data = false;
		g_number_of_links = 0;
		g_number_of_timing_arcs = 0;
		g_number_of_nodes = 0;

		b_debug_detail_flag = 1;

		g_informationCount = 0;


		SYSTEMTIME st;
		GetLocalTime(&st);
		messageType_int_to_string[0] = "MSG";//message
		messageType_int_to_string[1] = "OUT";//output
		messageType_int_to_string[2] = "WAR";//Warning

	}

	int signal_updating_output;

	int b_debug_detail_flag;
	std::map<int, int> g_internal_node_to_seq_no_map;  // hash table, map external node number to internal node sequence no. 


	std::map<string, int> g_link_id_map;


	int g_number_of_links;
	int g_number_of_timing_arcs;
	int g_number_of_nodes;
	int g_number_of_zones;

	int g_LoadingStartTimeInMin;
	int g_LoadingEndTimeInMin;
	int g_informationCount;

	bool b_with_loaded_data;


	std::map<int, char*> messageType_int_to_string;

	void WriteLog(int messageType, string info, int level)
	{
		if (log_out.log_sig)
		{
			string padding = "";
		if (level > 1)
		{
			for (size_t l = 2; l <= level; l++)
			{
				padding.append("--");
			}
		}
		if (level == -1)
		{
			padding = "++++++++++++++++++++++++++++++";
		}
		if (level == -2)
		{
			padding = "++++++++++++++++++++++++++++++";
		}

		log_out << "Info_ID: " << g_informationCount << "| Info_Type: " << messageType_int_to_string[messageType] << " \t| Info: "<< padding.c_str()  << info.c_str() << endl;
		
		//OutputDebugStringA(szLog);
		g_informationCount++;
		}
	}

	string Output_Combine_Message_And_Value_Int(string info,int value)
	{
		string value_String = to_string(value);

		info.append(": ");

		info.append(value_String);

		return info;
	}
	string Output_Combine_Message_And_Value_Double(string info, double value)
	{
		string value_String = to_string(value);

		info.append(": ");

		info.append(value_String);

		return info;
	}

	std::vector<CNode> m_sig_node_vector;
	std::vector<CLink> m_sig_link_vector;


};

CMainSigModual MainSigModual;
struct SMovementData
{
	//SMovementData()
	//{

	//}
	bool Enable = false;//no value->false
	int LinkSeqNo;
	int PhaseNo;
	//enum stage StageNo;
	vector<enum stage> StageNo_in_Order; //
	int Volume;
	int Lanes;
	int SharedLanes;
	enum group GroupNo;//L->1;T/R->2
	enum direction DirectionNo;



	//by Acquisition through processing
	int Assignment_Order;// not in use, Differentiate primary and secondary movements, then determine stages  according to volumes
	enum left_Turn_Treatment Left_Turn_Treatment;

	string linkID;

};

#pragma region Fields

//enum
enum movement_Index { EBL = 1, EBT = 2, EBR = 3, WBL = 4, WBT = 5, WBR = 6, NBL = 7, NBT = 8, NBR = 9, SBL = 10, SBT = 11, SBR = 12 };
enum stage { no_stage = -1, stage1 = 1, stage2 = 2, stage3 = 3, stage4 = 4 };
enum direction { E = 1, W = 2, N = 3, S = 4 };
enum group { L = 1, T_AND_R = 2 };
enum left_Turn_Treatment { perm = 0, prot = 1, no_business = -1 };

//array size, for constructing matirx or array
const int laneColumnSize = 32;
const int movementSize = 32;
const int NEMA_PhaseSize = 32;//temp enable=false
const int stageSize = 5;
const int ringSize = 3;
const int directionSize = 5;
const int groupSize = 5;

//index range, for indexing, from 1 to 0+range (including the 0+range th)

//parameters
double l = 12;
double x_c_Input = 0.9;
double PHF = 1;

double f_1 = 1;
double f_2 = 1;
double t_L = 4;
double t_Yellow = 4;
double t_AR = 2;

double minGreenTime = 5;
#pragma endregion

class CSignalNode
{
public:
	CSignalNode()//Constructor
	{
		y_StageMax = 1;
		x_c_output = 0.9;
		c_Min = 60;

		movement_str_to_index_map["EBL"] = EBL;
		movement_str_to_index_map["EBT"] = EBT;
		movement_str_to_index_map["EBR"] = EBR;

		movement_str_to_index_map["WBL"] = WBL;
		movement_str_to_index_map["WBT"] = WBT;
		movement_str_to_index_map["WBR"] = WBR;

		movement_str_to_index_map["NBL"] = NBL;
		movement_str_to_index_map["NBT"] = NBT;
		movement_str_to_index_map["NBR"] = NBR;

		movement_str_to_index_map["SBL"] = SBL;
		movement_str_to_index_map["SBT"] = SBT;
		movement_str_to_index_map["SBR"] = SBR;

		movement_str_array[EBL] = "EBL";
		movement_str_array[EBT] = "EBT";
		movement_str_array[EBR] = "EBR";
		movement_str_array[WBL] = "WBL";
		movement_str_array[WBT] = "WBT";
		movement_str_array[WBR] = "WBR";
		movement_str_array[NBL] = "NBL";
		movement_str_array[NBT] = "NBT";
		movement_str_array[NBR] = "NBR";
		movement_str_array[SBL] = "SBL";
		movement_str_array[SBT] = "SBT";
		movement_str_array[SBR] = "SBR";


		movement_str_to_direction_map["EBL"] = E;
		movement_str_to_direction_map["EBT"] = E;
		movement_str_to_direction_map["EBR"] = E;

		movement_str_to_direction_map["WBL"] = W;
		movement_str_to_direction_map["WBT"] = W;
		movement_str_to_direction_map["WBR"] = W;

		movement_str_to_direction_map["NBL"] = N;
		movement_str_to_direction_map["NBT"] = N;
		movement_str_to_direction_map["NBR"] = N;

		movement_str_to_direction_map["SBL"] = S;
		movement_str_to_direction_map["SBT"] = S;
		movement_str_to_direction_map["SBR"] = S;


		left_Movement_Opposing_Index_Map[EBL] = WBT;
		left_Movement_Opposing_Index_Map[WBL] = EBT;
		left_Movement_Opposing_Index_Map[NBL] = SBT;
		left_Movement_Opposing_Index_Map[SBL] = NBT;


		left_Movement_Counterpart_Index_Map[EBL] = EBT;
		left_Movement_Counterpart_Index_Map[WBL] = WBT;
		left_Movement_Counterpart_Index_Map[NBL] = NBT;
		left_Movement_Counterpart_Index_Map[SBL] = SBT;

		left_Movement_Counterpart_Right_Trun_Index_Map[EBL] = EBR;
		left_Movement_Counterpart_Right_Trun_Index_Map[WBL] = WBR;
		left_Movement_Counterpart_Right_Trun_Index_Map[NBL] = NBR;
		left_Movement_Counterpart_Right_Trun_Index_Map[SBL] = SBR;

		movement_Index_to_Group_Map[EBL] = group(1);
		movement_Index_to_Group_Map[EBT] = group(2);
		movement_Index_to_Group_Map[EBR] = group(2);
		movement_Index_to_Group_Map[WBL] = group(1);
		movement_Index_to_Group_Map[WBT] = group(2);
		movement_Index_to_Group_Map[WBR] = group(2);
		movement_Index_to_Group_Map[NBL] = group(1);
		movement_Index_to_Group_Map[NBT] = group(2);
		movement_Index_to_Group_Map[NBR] = group(2);
		movement_Index_to_Group_Map[SBL] = group(1);
		movement_Index_to_Group_Map[SBT] = group(2);
		movement_Index_to_Group_Map[SBR] = group(2);

		direction_index_to_str_map[E] = "E";
		direction_index_to_str_map[W] = "W";
		direction_index_to_str_map[N] = "N";
		direction_index_to_str_map[S] = "S";

		intersection_Average_Delay = 0;
		intersection_Total_Delay = 0;
		intersection_Total_Volume = 0;

		left_Turn_Treatment_index_to_str_map[prot] = "Protected";
		left_Turn_Treatment_index_to_str_map[perm] = "Permissive";
		left_Turn_Treatment_index_to_str_map[no_business] = "Null";


		for (int s = 0; s <= stageSize; s++)
		{
			y_Max_Stage_Array[s] = 0;

			for (int m = 0; m <= movementSize; m++)
			{
				saturation_Flow_Rate_Matrix[s][m] = 0;
				capacity_by_Stage_and_Movement_Matrix[s][m] = 0;
			}

			for (int d = 0; d <= directionSize; d++)
			{
				for (int g = 0; g <= groupSize; g++)
				{
					stage_Direction_Candidates_Matrix[s][d][g] = 0;

				}
				approach_Average_Delay_Array[d] = 0;
				approach_Total_Delay_Array[d] = 0;
				approach_Total_Volume_Array[d] = 0;
			}

		}

		for (int m = 0; m <= movementSize; m++)
		{
			movement_Array[m].Enable = false;
			movement_Array[m].Volume = 0;
		}
	}

	int movement_Range = 12;
	int direction_Range = 4;

	std::map<string, enum movement_Index> movement_str_to_index_map;
	std::map<enum direction, string > direction_index_to_str_map;
	std::map<enum left_Turn_Treatment, string > left_Turn_Treatment_index_to_str_map;


	string movement_str_array[movementSize + 1];
		
	std::map<string, enum direction> movement_str_to_direction_map;

	std::map<enum movement_Index, enum movement_Index> left_Movement_Opposing_Index_Map;

	std::map<enum movement_Index, enum movement_Index> left_Movement_Counterpart_Index_Map;//L -> T

	std::map<enum movement_Index, enum movement_Index> left_Movement_Counterpart_Right_Trun_Index_Map;//L -> R


	std::map<enum movement_Index, enum group> movement_Index_to_Group_Map;


	//array

	SMovementData movement_Array[movementSize + 1]; //indexed by enum movement_Direction   --Step 2
	int green_Start_Stage_Array[stageSize + 1];
	int green_End_Stage_Array[stageSize+1];
	double y_Max_Stage_Array[stageSize+1];
	double green_Time_Stage_Array[stageSize + 1];
	double cumulative_Green_Start_Time_Stage_Array[stageSize + 1];
	double cumulative_Green_End_Time_Stage_Array[stageSize + 1];
	double effective_Green_Time_Stage_Array[stageSize + 1];
	double cumulative_Effective_Green_Start_Time_Stage_Array[stageSize + 1];
	double cumulative_Effective_Green_End_Time_Stage_Array[stageSize + 1];
	double ratio_of_Effective_Green_Time_to_Cycle_Length_Array[stageSize + 1];

	double approach_Average_Delay_Array[directionSize+1];
	double approach_Total_Delay_Array[directionSize + 1];
	double approach_Total_Volume_Array[directionSize + 1];



	//matrix
	double saturation_Flow_Rate_Matrix[stageSize+1][movementSize + 1];//s
	double y_Stage_Movement_Matrix[stageSize + 1][movementSize + 1];
	double stage_Direction_Candidates_Matrix[stageSize + 1][directionSize + 1][groupSize + 1];

	double capacity_by_Stage_and_Movement_Matrix[stageSize + 1][movementSize + 1];
	double v_over_C_by_Stage_and_Movement_Matrix[stageSize + 1][movementSize + 1];

	double average_Uniform_Delay_Matrix[stageSize + 1][movementSize + 1];


	int NEMA_Phase_Matrix[5][3] = { 0,0,0,0,1,5,0,2,6,0,3,7,0,4,8 };//row->NEMA_phases; col->rings
	int green_Start_NEMA_Phase[5][3];

	SMovementData Stage_Ring_Movement_Matrix[ringSize + 1][stageSize + 1];

	//variables
	int stage_Range;
	double y_StageMax;
	double x_c_output;
	double c_Min;
	double c_Optimal;
	double c_Final;

	double intersection_Average_Delay;
	double intersection_Total_Delay;
	double intersection_Total_Volume;

	string LOS;

	double default_volume_T = 30;

	void PerformQEM(int nodeID)
	{
		MainSigModual.WriteLog(0, "Data Loading : Movement Volume of Main Node", 1);
		AddMovementVolume_from_LinkPerformance();

		MainSigModual.WriteLog(0, "Step 2: Perform QEM", 1);

		MainSigModual.WriteLog(0, "Step 2.1: Set Left Turn Treatments",2);
		Set_Left_Turn_Treatment();

		MainSigModual.WriteLog(0, "Step 2.2: Set StageNos",2);
		Set_StageNo_for_Movements(); 

		MainSigModual.WriteLog(0, "Step 2.3: Set Saturation Flow Rate Matrix",2);
		Set_Saturation_Flow_Rate_Matrix(); 

		MainSigModual.WriteLog(0, "Step 2.4: Calculate Flow Ratio",2);
		Calculate_Flow_of_Ratio_Max();
		MainSigModual.WriteLog(1, MainSigModual.Output_Combine_Message_And_Value_Double("Max y Value", y_StageMax),3);

		MainSigModual.WriteLog(0, "Step 2.5: Calculate Total Cycle Lost Time", 2);
		Calculate_Total_Cycle_Lost_Time();
		MainSigModual.WriteLog(1, MainSigModual.Output_Combine_Message_And_Value_Double("Total Cycle Lost Time", l), 3);

		MainSigModual.WriteLog(0, "Step 2.6: Calculate the Minimum and Optimal Cycle Length",2);
		Calculate_the_Minimum_And_Optimal_Cycle_Length();
		MainSigModual.WriteLog(1, MainSigModual.Output_Combine_Message_And_Value_Double("Min Cycle Length", c_Min),3);
		MainSigModual.WriteLog(1, MainSigModual.Output_Combine_Message_And_Value_Double("Optimal Cycle Length", c_Optimal), 3);


		MainSigModual.WriteLog(0, "Step 2.7: Recalculate xc",2);
		Calculate_the_x_c_Output();
		MainSigModual.WriteLog(1, MainSigModual.Output_Combine_Message_And_Value_Double("Recalculated xc", x_c_output),3);

		MainSigModual.WriteLog(0, "Step 2.8: Timing for Stages",2);
		Calculate_Green_Time_for_Stages();
		Printing_Green_Time_for_Stages();

		MainSigModual.WriteLog(0, "Step 2.9: Calculate capacity and ratio V/C", 2);
		Calculate_Capacity_And_Ratio_V_over_C();

		MainSigModual.WriteLog(0, "Step 2.10: Calculate Signal Delay", 2);
		Calculate_Signal_Delay(nodeID);

		MainSigModual.WriteLog(0, "Step 2.11: evaluate LOS", 2);
		Judge_Signal_LOS(nodeID);
	}


	void AddMovementStructure(int link_seq_no, string str, int lanes, int sharedLanes, string linkID)
	{
		enum movement_Index mi = movement_str_to_index_map[str];
		enum direction di = movement_str_to_direction_map[str];

		movement_Array[mi].Enable = true;
		movement_Array[mi].LinkSeqNo = link_seq_no;
		//movement_Array[mi].Volume = volume;//0
		//movement_Array[mi].StageNo = stage(movement_to_Stage_Array[mi]);
		movement_Array[mi].GroupNo = movement_Index_to_Group_Map[mi];
		movement_Array[mi].DirectionNo = di;
		movement_Array[mi].Left_Turn_Treatment = no_business;


		movement_Array[mi].Lanes = lanes;
		movement_Array[mi].SharedLanes = sharedLanes;
		movement_Array[mi].linkID = linkID;
	}

	void AddMovementVolume_from_LinkPerformance()
	{
		CCSVParser parser_link_Performance;
		int record = 0;
		if (parser_link_Performance.OpenCSVFile("link_performance_sig.csv", true))
		{

			float demand_duration_in_min = max(15, MainSigModual.g_LoadingEndTimeInMin - MainSigModual.g_LoadingStartTimeInMin);


			bool flag_reading = parser_link_Performance.ReadRecord();
			if (record == 0 && flag_reading == false)
			{
				for (size_t m = 1; m <= movement_Range; m++)
				{
					if (m % 3 == 2&& movement_Array[m].Enable==true)
					{
						float hourly_volume;
						hourly_volume = default_volume_T / (demand_duration_in_min / 60);
						movement_Array[m].Volume = hourly_volume;
					}

				}

				return;
			}
			while (flag_reading)  // if this line contains [] mark, then we will also read field headers.
			{
				record++;
				string linkID_Link_Performance;
				parser_link_Performance.GetValueByFieldName("link_id", linkID_Link_Performance);

				string demand_period;
				parser_link_Performance.GetValueByFieldName("time_period", demand_period);

				if (demand_period.length() > 0)
				{
					int flag_vacancy = 0;
					float volume_demand_period_array[13];
					for (size_t m = 1; m <= movement_Range; m++)
					{
						if (movement_Array[m].Enable == true && linkID_Link_Performance == movement_Array[m].linkID)
						{
							float volume_demand_period;
							parser_link_Performance.GetValueByFieldName("volume", volume_demand_period);
							if (volume_demand_period == 0)
							{
								flag_vacancy++;
							}
							volume_demand_period_array[m] = volume_demand_period;
						}
					}

					for (size_t m = 1; m <= movement_Range; m++)
					{
						float hourly_volume;
						if (flag_vacancy == 12 && m % 3 == 2&& movement_Array[m].Enable == true)
						{
							hourly_volume = default_volume_T / (demand_duration_in_min / 60);
						}
						else
						{
							hourly_volume = volume_demand_period_array[m] / (demand_duration_in_min / 60);

						}
						movement_Array[m].Volume = hourly_volume;

					}

				}
			}
		}
		parser_link_Performance.CloseCSVFile();
	}

	void Set_Left_Turn_Treatment()
	{
		//20200808 Add a completeness check. Assign high privileges if only left.
		for (size_t m = 1; m <= movement_Range; m++)
		{
			int final_decision;
			if (movement_Array[m].GroupNo == 1)
			{
				final_decision = 0;
				//(1)	Left-turn Lane Check
				if (movement_Array[m].Lanes > 1)
				{
					final_decision = 1;
				}
				//(2)	Minimum Volume Check
				if (movement_Array[m].Volume >= 240)
				{
					final_decision = 1;
				}
				//(3)	Opposing Through Lanes Check
				int op_Movement_Index = left_Movement_Opposing_Index_Map[movement_Index(m)];
				if (movement_Array[op_Movement_Index].Lanes >= 4)
				{
					final_decision = 1;
				}
				//(4)	Opposing Traffic Speed Check

				//(5)	Minimum Cross-Product Check
				int co_Movement_Index = left_Movement_Opposing_Index_Map[movement_Index(m)];
				if (movement_Array[co_Movement_Index].Lanes > 1)
				{
					if (movement_Array[co_Movement_Index].Volume * movement_Array[m].Volume >= 100000)
					{
						final_decision = 1;
					}
				}
				else
				{

					if (movement_Array[co_Movement_Index].Volume * movement_Array[m].Volume >= 50000)
					{
						final_decision = 1;
					}
				}
				//(6) if there is no T movement, then the left movement should be protected.
				if (movement_Array[left_Movement_Counterpart_Index_Map[movement_Index(m)]].Enable==false)
				{
					final_decision = 1;
				}
			}
			else
			{
				final_decision = -1;
			}
			movement_Array[m].Left_Turn_Treatment = left_Turn_Treatment(final_decision);

			if (final_decision != -1)
			{
				g_info_String = "Left Turn Treatment of Movement_";
				g_info_String.append(movement_str_array[m]);
				g_info_String.append(": ");
				g_info_String.append(left_Turn_Treatment_index_to_str_map[left_Turn_Treatment(final_decision)]);
				MainSigModual.WriteLog(1, g_info_String, 3);
			}

		}
	}

	void Set_StageNo_for_Movements()
	{
		//Determining the main direction
		int east_And_West_Volume = movement_Array[EBL].Volume + movement_Array[EBT].Volume + movement_Array[EBR].Volume +
			movement_Array[WBL].Volume + movement_Array[WBT].Volume + movement_Array[WBR].Volume;

		int north_And_South_Volume = movement_Array[NBL].Volume + movement_Array[NBT].Volume + movement_Array[NBR].Volume +
			movement_Array[SBL].Volume + movement_Array[SBT].Volume + movement_Array[SBR].Volume;

		stage_Range = 2;
		bool east_And_West_Flag = false;
		bool north_And_South_Flag = false;

		if (movement_Array[EBL].Left_Turn_Treatment == prot)
		{
			stage_Range++;
			east_And_West_Flag = true;
		}
		else if (movement_Array[WBL].Left_Turn_Treatment == prot)
		{
			stage_Range++;
			east_And_West_Flag = true;

		}

		if (movement_Array[NBL].Left_Turn_Treatment == prot)
		{
			stage_Range++;
			north_And_South_Flag = true;
		}
		else if (movement_Array[SBL].Left_Turn_Treatment == prot)
		{
			stage_Range++;
			north_And_South_Flag = true;

		}

		enum movement_Index firstL;
		enum movement_Index firstT;
		enum movement_Index firstR;
		enum movement_Index secondL;
		enum movement_Index secondT;
		enum movement_Index secondR;
		enum movement_Index thridL;
		enum movement_Index thridT;
		enum movement_Index thridR;
		enum movement_Index fouthL;
		enum movement_Index fouthT;
		enum movement_Index fouthR;
		if (east_And_West_Volume >= north_And_South_Volume)
		{
			//east and west first
			firstL = EBL;
			firstT = EBT;
			firstR = EBR;
			secondL = WBL;
			secondT = WBT;
			secondR = WBR;
			thridL = NBL;
			thridT = NBT;
			thridR = NBR;
			fouthL = SBL;
			fouthT = SBT;
			fouthR = SBR;
			MainSigModual.WriteLog(0, "Main Approaches: E & W", 3);
		}
		else
		{

			firstL = NBL;
			firstT = NBT;
			firstR = NBR;
			secondL = SBL;
			secondT = SBT;
			secondR = SBR;
			thridL = EBL;
			thridT = EBT;
			thridR = EBR;
			fouthL = WBL;
			fouthT = WBT;
			fouthR = WBR;
			MainSigModual.WriteLog(0, "Main Approaches: N & S", 3);

		}


		if (east_And_West_Flag)
		{
			movement_Array[firstL].StageNo_in_Order.push_back(stage1);
			movement_Array[firstT].StageNo_in_Order.push_back(stage2);
			movement_Array[firstR].StageNo_in_Order.push_back(stage2);
			movement_Array[secondL].StageNo_in_Order.push_back(stage1);
			movement_Array[secondT].StageNo_in_Order.push_back(stage2);
			movement_Array[secondR].StageNo_in_Order.push_back(stage2);
			if (north_And_South_Flag)
			{
				movement_Array[thridL].StageNo_in_Order.push_back(stage3);
				movement_Array[thridT].StageNo_in_Order.push_back(stage4);
				movement_Array[thridR].StageNo_in_Order.push_back(stage4);
				movement_Array[fouthL].StageNo_in_Order.push_back(stage3);
				movement_Array[fouthT].StageNo_in_Order.push_back(stage4);
				movement_Array[fouthR].StageNo_in_Order.push_back(stage4);
			}
			else
			{
				movement_Array[thridL].StageNo_in_Order.push_back(stage3);
				movement_Array[thridT].StageNo_in_Order.push_back(stage3);
				movement_Array[thridR].StageNo_in_Order.push_back(stage3);
				movement_Array[fouthL].StageNo_in_Order.push_back(stage3);
				movement_Array[fouthT].StageNo_in_Order.push_back(stage3);
				movement_Array[fouthR].StageNo_in_Order.push_back(stage3);
			}

		}
		else
		{
			movement_Array[firstL].StageNo_in_Order.push_back(stage1);
			movement_Array[firstT].StageNo_in_Order.push_back(stage1);
			movement_Array[firstR].StageNo_in_Order.push_back(stage1);
			movement_Array[secondL].StageNo_in_Order.push_back(stage1);
			movement_Array[secondT].StageNo_in_Order.push_back(stage1);
			movement_Array[secondR].StageNo_in_Order.push_back(stage1);
			if (north_And_South_Flag)
			{
				movement_Array[thridL].StageNo_in_Order.push_back(stage2);
				movement_Array[thridT].StageNo_in_Order.push_back(stage3);
				movement_Array[thridR].StageNo_in_Order.push_back(stage3);
				movement_Array[fouthL].StageNo_in_Order.push_back(stage2);
				movement_Array[fouthT].StageNo_in_Order.push_back(stage3);
				movement_Array[fouthR].StageNo_in_Order.push_back(stage3);
			}
			else
			{
				movement_Array[thridL].StageNo_in_Order.push_back(stage2);
				movement_Array[thridT].StageNo_in_Order.push_back(stage2);
				movement_Array[thridR].StageNo_in_Order.push_back(stage2);
				movement_Array[fouthL].StageNo_in_Order.push_back(stage2);
				movement_Array[fouthT].StageNo_in_Order.push_back(stage2);
				movement_Array[fouthR].StageNo_in_Order.push_back(stage2);
			}
		}

		//20200808 check the enable property, delete stage from enable=false movement
		//20200812 modified this part by using vector StageNo_in_Order
		int initialStage_Range = stage_Range;
		int criticalStage = -1;

		for (size_t s = initialStage_Range; s >= 1; s--)
		{
			int checkNumber = 0;
			if (criticalStage != -1)
			{
				for (size_t m = 1; m <= movement_Range; m++)
				{
					for (size_t so = 0; so < movement_Array[m].StageNo_in_Order.size(); so++)
					{
						if (movement_Array[m].StageNo_in_Order[so] > criticalStage && movement_Array[m].Enable == true)
						{
							movement_Array[m].StageNo_in_Order[so] = stage(movement_Array[m].StageNo_in_Order[so] - 1);
						}
					}

				}
			}


			for (size_t m = 1; m <= movement_Range; m++)
			{
				for (size_t so = 0; so < movement_Array[m].StageNo_in_Order.size(); so++)
				{
					if (movement_Array[m].Enable == true && movement_Array[m].StageNo_in_Order[so] == s)
					{
						checkNumber++;
					}
				}
			}

			if (checkNumber == 0)
			{
				stage_Range--;
				criticalStage = s;
			}
			else
			{
				criticalStage = -1;
			}
		}
		for (size_t m = 1; m <= movement_Range; m++)
		{
			if (movement_Array[m].Enable == false)
			{
				movement_Array[m].StageNo_in_Order[0] = stage(- 1);
			}
		}

		//TODO: we can scan stages row-wise and mark the property of each stage, left-protected for example.

		//for (size_t s = 0; s < stage_Range; s++)
		//{

		//}

		//20200812 add right-turn treatment for movements 
		for (size_t m = 1; m <= movement_Range; m++)
		{
			if (movement_Array[m].Enable==true && movement_Array[m].GroupNo==1)
			{
				if (movement_Array[left_Movement_Counterpart_Right_Trun_Index_Map[movement_Index(m)]].Enable == true)
				{
					movement_Array[left_Movement_Counterpart_Right_Trun_Index_Map[movement_Index(m)]].StageNo_in_Order.insert(movement_Array[left_Movement_Counterpart_Right_Trun_Index_Map[movement_Index(m)]].StageNo_in_Order.begin(), movement_Array[m].StageNo_in_Order[0]);
				}
			}

		}



		//movement_Array[EBL].StageNo = stage1;
		//movement_Array[EBT].StageNo = stage2;
		//movement_Array[EBR].StageNo = stage2;
		//movement_Array[WBL].StageNo = stage1;
		//movement_Array[WBT].StageNo = stage2;
		//movement_Array[WBR].StageNo = stage2;
		//movement_Array[NBL].StageNo = stage3;
		//movement_Array[NBT].StageNo = stage3;
		//movement_Array[NBR].StageNo = stage3;
		//movement_Array[SBL].StageNo = stage3;
		//movement_Array[SBT].StageNo = stage3;
		//movement_Array[SBR].StageNo = stage3;

		g_info_String = "Number of Stages: ";
		g_info_String.append(to_string(stage_Range));
		MainSigModual.WriteLog(1, g_info_String, 3);


		for (size_t m = 1; m <= movement_Range; m++)
		{
			if (movement_Array[m].Enable == true)
			{
				for (size_t so = 0; so < movement_Array[m].StageNo_in_Order.size(); so++)
				{
					if (movement_Array[m].StageNo_in_Order[so] > stage_Range)
					{
						movement_Array[m].StageNo_in_Order[so] = stage(stage_Range);
					}
					g_info_String = "StageNo of Movement_";
					g_info_String.append(movement_str_array[m]);
					g_info_String.append(": ");
					g_info_String.append("Stage_");
					g_info_String.append(to_string(movement_Array[m].StageNo_in_Order[so]));
					MainSigModual.WriteLog(1, g_info_String, 4);
				}
			}

		}
	}

	void Set_Saturation_Flow_Rate_Matrix()
	{
		//movement_Array left turn movement
		// we need to use the saturation flow rate values based on protected and permitted

		for (size_t m = 1; m <= movement_Range; m++)
		{
			if (movement_Array[m].Enable == false)
			{
				continue;
			}
			for (size_t so = 0; so < movement_Array[m].StageNo_in_Order.size(); so++)
			{
				if (movement_Array[m].Left_Turn_Treatment == prot)
				{
					saturation_Flow_Rate_Matrix[movement_Array[m].StageNo_in_Order[so]][m] = 1530 * movement_Array[m].Lanes * PHF;
				}
				else if (movement_Array[m].Left_Turn_Treatment == perm)
				{
					int op_Movement_Index = left_Movement_Opposing_Index_Map[movement_Index(m)];
					int op_volume = movement_Array[op_Movement_Index].Volume;
					saturation_Flow_Rate_Matrix[movement_Array[m].StageNo_in_Order[so]][m] = f_1 * f_2 * op_volume * (exp(-op_volume * 4.5 / 3600)) / (1 - exp(op_volume * 2.5 / 3600));
				}
				else
				{
					saturation_Flow_Rate_Matrix[movement_Array[m].StageNo_in_Order[so]][m] = 1530 * movement_Array[m].Lanes;//temp!!
				}
				g_info_String = "Saturation Flow Rate of Movement_";
				g_info_String.append(movement_str_array[m]);
				g_info_String.append(" and ");
				g_info_String.append("Stage_");
				g_info_String.append(to_string(movement_Array[m].StageNo_in_Order[so]));
				g_info_String.append(": ");
				g_info_String.append(to_string(saturation_Flow_Rate_Matrix[movement_Array[m].StageNo_in_Order[so]][m]));
				MainSigModual.WriteLog(1, g_info_String, 3);
			}
		}

		//saturation_Flow_Rate_Matrix[1][EBL] = 1750;
		//saturation_Flow_Rate_Matrix[1][WBL] = 1750;
		//saturation_Flow_Rate_Matrix[2][EBT] = 3400;
		//saturation_Flow_Rate_Matrix[2][EBR] = 3400;
		//saturation_Flow_Rate_Matrix[2][WBT] = 3400;
		//saturation_Flow_Rate_Matrix[2][WBR] = 3400;
		//saturation_Flow_Rate_Matrix[3][NBL] = 475;
		//saturation_Flow_Rate_Matrix[3][NBT] = 1800;
		//saturation_Flow_Rate_Matrix[3][NBR] = 1800;
		//saturation_Flow_Rate_Matrix[3][SBL] = 450;
		//saturation_Flow_Rate_Matrix[3][SBT] = 1800;
		//saturation_Flow_Rate_Matrix[3][SBR] = 1800;
	}

	void Calculate_Flow_of_Ratio_Max()
	{
		//y_Stage_Movement_Matrix
		//y_StageMax
		for (size_t s = 1; s <= stage_Range; s++)
		{
			y_Max_Stage_Array[s] = 0;

			for (size_t m = 1; m <= movement_Range; m++)
			{
				for (size_t so = 0; so < movement_Array[m].StageNo_in_Order.size(); so++)
				{
					if (saturation_Flow_Rate_Matrix[s][m] != 0 && movement_Array[m].Enable && movement_Array[m].StageNo_in_Order[so] == s)
					{
						y_Stage_Movement_Matrix[s][m] = double(movement_Array[m].Volume) / double(saturation_Flow_Rate_Matrix[s][m]);

						//double stage_Direction_Candidates_Matrix[stageSize][directionSize][groupSize]
						stage_Direction_Candidates_Matrix[s][movement_Array[m].DirectionNo][movement_Array[m].GroupNo] += y_Stage_Movement_Matrix[s][m];

						// we tally the movement matrix from this direction and this group number, so we can distingush movements belonging to different directions 
						if (stage_Direction_Candidates_Matrix[s][movement_Array[m].DirectionNo][movement_Array[m].GroupNo] >= y_Max_Stage_Array[s])
						{
							y_Max_Stage_Array[s] = stage_Direction_Candidates_Matrix[s][movement_Array[m].DirectionNo][movement_Array[m].GroupNo];
						}
					}
				}
			}
		}

		y_StageMax = 0;
		for (size_t i = 1; i <= stage_Range; i++)
		{
			y_StageMax += y_Max_Stage_Array[i];
		}

	}

	void Calculate_Total_Cycle_Lost_Time()
	{
		l = t_L * stage_Range;
	}

	void Calculate_the_Minimum_And_Optimal_Cycle_Length()
	{
		c_Min = max(60, (l - x_c_Input) / (x_c_Input - y_StageMax));

		c_Min = int( (c_Min + 5) / 10) * 10;  // take a number near a multiplier of 10

		c_Optimal= max(60, (1.5*l +5) / (1 - y_StageMax));
	}

	void Calculate_the_x_c_Output()
	{
		x_c_output = (y_StageMax * c_Min) / (c_Min - l);
	}

	void Calculate_Green_Time_for_Stages()
	{
		for (size_t s = 1; s <= stage_Range; s++)
		{
			//			green_Time_Stage_Array[i] = y_Max_Stage_Array[i] * c_Min / x_c_output;

			if (y_StageMax == 0)
			{
				green_Time_Stage_Array[s] = minGreenTime;
			}
			else
			{
				green_Time_Stage_Array[s] = max(minGreenTime, y_Max_Stage_Array[s] * c_Min / max(0.1, y_StageMax));

			}

			green_Time_Stage_Array[s] = (int)(green_Time_Stage_Array[s] + 0.5);
			effective_Green_Time_Stage_Array[s] = green_Time_Stage_Array[s] - t_L + t_Yellow + t_AR;
			ratio_of_Effective_Green_Time_to_Cycle_Length_Array[s] = effective_Green_Time_Stage_Array[s] / c_Min;
		}

	}

	void Printing_Green_Time_for_Stages()//output Effective green time 
	{
		c_Final = 0;

		cumulative_Green_Start_Time_Stage_Array[1] = 0;
		cumulative_Green_End_Time_Stage_Array[1] = green_Time_Stage_Array[1];

		c_Final = cumulative_Green_End_Time_Stage_Array[1];

		MainSigModual.WriteLog(1, MainSigModual.Output_Combine_Message_And_Value_Double("Green Time of Stage 1", green_Time_Stage_Array[1]),3);
		MainSigModual.WriteLog(1, MainSigModual.Output_Combine_Message_And_Value_Double("Start Green Time of Stage 1", cumulative_Green_Start_Time_Stage_Array[1]), 4);
		MainSigModual.WriteLog(1, MainSigModual.Output_Combine_Message_And_Value_Double("End Green Time of Stage 1", cumulative_Green_End_Time_Stage_Array[1]), 4);


		for (size_t i = 2; i <= stage_Range; i++)
		{
			g_info_String = "Green Time of Stage ";
			g_info_String.append(to_string(i));
			MainSigModual.WriteLog(1, MainSigModual.Output_Combine_Message_And_Value_Double(g_info_String, green_Time_Stage_Array[i]),3);
			cumulative_Green_Start_Time_Stage_Array[i] = cumulative_Green_End_Time_Stage_Array[i - 1];
			cumulative_Green_End_Time_Stage_Array[i] = cumulative_Green_Start_Time_Stage_Array[i] + green_Time_Stage_Array[i];

			if (cumulative_Green_End_Time_Stage_Array[i] > c_Final)
				c_Final = cumulative_Green_End_Time_Stage_Array[i];  // calculating the final cycle as a result of loss and green time

			g_info_String = "Start Green Time of Stage ";
			g_info_String.append(to_string(i));
			MainSigModual.WriteLog(1, MainSigModual.Output_Combine_Message_And_Value_Double(g_info_String, cumulative_Green_Start_Time_Stage_Array[i]), 4);

			g_info_String = "End Green Time of Stage ";
			g_info_String.append(to_string(i));
			MainSigModual.WriteLog(1, MainSigModual.Output_Combine_Message_And_Value_Double(g_info_String, cumulative_Green_End_Time_Stage_Array[i]), 4);
		}



		cumulative_Effective_Green_Start_Time_Stage_Array[1] = 0;
		cumulative_Effective_Green_End_Time_Stage_Array[1] = effective_Green_Time_Stage_Array[1];

		MainSigModual.WriteLog(1, MainSigModual.Output_Combine_Message_And_Value_Double("Effective Green Time of Stage 1", effective_Green_Time_Stage_Array[1]), 3);
		MainSigModual.WriteLog(1, MainSigModual.Output_Combine_Message_And_Value_Double("Start Effective Green Time of Stage 1", cumulative_Effective_Green_Start_Time_Stage_Array[1]), 4);
		MainSigModual.WriteLog(1, MainSigModual.Output_Combine_Message_And_Value_Double("End Effective Green Time of Stage 1", cumulative_Effective_Green_End_Time_Stage_Array[1]), 4);


		for (size_t i = 2; i <= stage_Range; i++)
		{
			g_info_String = "Effective Green Time of Stage ";
			g_info_String.append(to_string(i));
			MainSigModual.WriteLog(1, MainSigModual.Output_Combine_Message_And_Value_Double(g_info_String, effective_Green_Time_Stage_Array[i]), 3);
			cumulative_Effective_Green_Start_Time_Stage_Array[i] = cumulative_Effective_Green_End_Time_Stage_Array[i - 1];
			cumulative_Effective_Green_End_Time_Stage_Array[i] = cumulative_Effective_Green_Start_Time_Stage_Array[i] + effective_Green_Time_Stage_Array[i];

			g_info_String = "Start Effective Green Time of Stage ";
			g_info_String.append(to_string(i));
			MainSigModual.WriteLog(1, MainSigModual.Output_Combine_Message_And_Value_Double(g_info_String, cumulative_Effective_Green_Start_Time_Stage_Array[i]), 4);

			g_info_String = "End Effective Green Time of Stage ";
			g_info_String.append(to_string(i));
			MainSigModual.WriteLog(1, MainSigModual.Output_Combine_Message_And_Value_Double(g_info_String, cumulative_Effective_Green_End_Time_Stage_Array[i]), 4);
		}
	}

	void Calculate_Capacity_And_Ratio_V_over_C()
	{
		for (size_t s = 1; s <= stage_Range; s++)
		{
			for (size_t m = 1; m <= movement_Range; m++)
			{
				for (size_t so = 0; so < movement_Array[m].StageNo_in_Order.size(); so++)
				{
					if (saturation_Flow_Rate_Matrix[s][m] != 0 && movement_Array[m].Enable && movement_Array[m].StageNo_in_Order[so] == s)
					{
						capacity_by_Stage_and_Movement_Matrix[s][m] = saturation_Flow_Rate_Matrix[s][m] * ratio_of_Effective_Green_Time_to_Cycle_Length_Array[s];
						g_info_String = "a. Capacity of Stage_";
						g_info_String.append(to_string(s));
						g_info_String.append(" and Movement_");
						g_info_String.append(movement_str_array[m]);
						g_info_String.append(": ");
						g_info_String.append(to_string(capacity_by_Stage_and_Movement_Matrix[s][m]));
						MainSigModual.WriteLog(1, g_info_String, 3);
						v_over_C_by_Stage_and_Movement_Matrix[s][m] = movement_Array[m].Volume / capacity_by_Stage_and_Movement_Matrix[s][m];
						g_info_String = "b. V/C of Stage_";
						g_info_String.append(to_string(s));
						g_info_String.append(" and Movement_");
						g_info_String.append(movement_str_array[m]);
						g_info_String.append(": ");
						g_info_String.append(to_string(v_over_C_by_Stage_and_Movement_Matrix[s][m]));
						MainSigModual.WriteLog(1, g_info_String, 3);						
					}
				}
			}
		}
	}

	void Calculate_Signal_Delay(int nodeID)
	{
		for (size_t s = 1; s <= stage_Range; s++)
		{
			for (size_t m = 1; m <= movement_Range; m++)
			{
				for (size_t so = 0; so < movement_Array[m].StageNo_in_Order.size(); so++)
				{
					if (saturation_Flow_Rate_Matrix[s][m] != 0 && movement_Array[m].Enable && movement_Array[m].StageNo_in_Order[so] == s)
					{
						average_Uniform_Delay_Matrix[s][m] = 
							(0.5 * capacity_by_Stage_and_Movement_Matrix[s][m] * pow((1 - ratio_of_Effective_Green_Time_to_Cycle_Length_Array[s]), 2)) 
							/ (1 - v_over_C_by_Stage_and_Movement_Matrix[s][m] * ratio_of_Effective_Green_Time_to_Cycle_Length_Array[s]);

						//(2)	Average Incremental Delay          next time 
						g_info_String = "Average Uniform Delay of Stage_";
						g_info_String.append(to_string(s));
						g_info_String.append(" and Movement_");
						g_info_String.append(movement_str_array[m]);
						g_info_String.append(": ");
						g_info_String.append(to_string(average_Uniform_Delay_Matrix[s][m]));
						MainSigModual.WriteLog(1, g_info_String, 3);


						approach_Total_Delay_Array[movement_Array[m].DirectionNo] += movement_Array[m].Volume * average_Uniform_Delay_Matrix[s][m];
						approach_Total_Volume_Array[movement_Array[m].DirectionNo] += movement_Array[m].Volume;
					}
				}
			}
		}
		for (size_t d = 1; d <= direction_Range; d++)
		{
			if (approach_Total_Volume_Array[d] == 0)
			{
				continue;
			}
			approach_Average_Delay_Array[d] = approach_Total_Delay_Array[d] / approach_Total_Volume_Array[d];
			g_info_String = "Average Delay of Approach_";
			g_info_String.append(direction_index_to_str_map[direction(d)]);
			g_info_String.append(": ");
			g_info_String.append(to_string(approach_Average_Delay_Array[d]));
			MainSigModual.WriteLog(1, g_info_String, 3);



			intersection_Total_Delay += approach_Average_Delay_Array[d] * approach_Total_Volume_Array[d];
			intersection_Total_Volume += approach_Total_Volume_Array[d];
		}
		intersection_Average_Delay = intersection_Total_Delay / intersection_Total_Volume;

		g_info_String = "Total Delay of Intersection_NodeID_";
		g_info_String.append(to_string(nodeID));
		g_info_String.append(": ");
		g_info_String.append(to_string(intersection_Total_Delay));
		MainSigModual.WriteLog(1, g_info_String, 3);



		g_info_String = "Average Delay of Intersection_NodeID_";
		g_info_String.append(to_string(nodeID));
		g_info_String.append(": ");
		g_info_String.append(to_string(intersection_Average_Delay));
		MainSigModual.WriteLog(1, g_info_String, 3);
	}

	void Judge_Signal_LOS(int nodeID)
	{
		if (intersection_Total_Delay <= 10)
		{
			LOS = "A";
		}
		else if (intersection_Total_Delay <= 20)
		{
			LOS = "B";
		}
		else if (intersection_Total_Delay <= 35)
		{
			LOS = "C";
		}
		else if (intersection_Total_Delay <= 55)
		{
			LOS = "D";
		}
		else if (intersection_Total_Delay <= 80)
		{
			LOS = "E";
		}
		else
		{
			LOS = "F";
		}
		g_info_String = "LOS of Intersection_NodeID_";
		g_info_String.append(to_string(nodeID));
		g_info_String.append(": ");
		g_info_String.append(LOS);
		MainSigModual.WriteLog(1, g_info_String, 3);
	}

	int signal_node_seq_no;  // sequence number 
	int main_node_id;      //external node number 

	std::vector<int> m_movement_link_seq_no_vector;

	//void Set_Movement_Array()
//{
//	//movement_Array 
//	movement_Array[1].Enable = true;
//	movement_Array[1].Volume = 300;
//	movement_Array[1].StageNo = stage(movement_to_Stage_Array[1]);
//	movement_Array[1].GroupNo = group(1);
//	movement_Array[1].DirectionNo = E;

//	movement_Array[2].Enable = true;
//	movement_Array[2].Volume = 900;
//	movement_Array[2].StageNo = stage(movement_to_Stage_Array[2]);
//	movement_Array[2].GroupNo = group(2);
//	movement_Array[2].DirectionNo = E;

//	movement_Array[3].Enable = true;
//	movement_Array[3].Volume = 200;
//	movement_Array[3].StageNo = stage(movement_to_Stage_Array[3]);
//	movement_Array[3].GroupNo = group(2);
//	movement_Array[3].DirectionNo = E;

//	movement_Array[4].Enable = true;
//	movement_Array[4].Volume = 250;
//	movement_Array[4].StageNo = stage(movement_to_Stage_Array[4]);
//	movement_Array[4].GroupNo = group(1);
//	movement_Array[4].DirectionNo = W;

//	movement_Array[5].Enable = true;
//	movement_Array[5].Volume = 1000;
//	movement_Array[5].StageNo = stage(movement_to_Stage_Array[5]);
//	movement_Array[5].GroupNo = group(2);
//	movement_Array[5].DirectionNo = W;

//	movement_Array[6].Enable = true;
//	movement_Array[6].Volume = 150;
//	movement_Array[6].StageNo = stage(movement_to_Stage_Array[6]);
//	movement_Array[6].GroupNo = group(2);
//	movement_Array[6].DirectionNo = W;

//	movement_Array[7].Enable = true;
//	movement_Array[7].Volume = 70;
//	movement_Array[7].StageNo = stage(movement_to_Stage_Array[7]);
//	movement_Array[7].GroupNo = group(1);
//	movement_Array[7].DirectionNo = N;

//	movement_Array[8].Enable = true;
//	movement_Array[8].Volume = 310;
//	movement_Array[8].StageNo = stage(movement_to_Stage_Array[8]);
//	movement_Array[8].GroupNo = group(2);
//	movement_Array[8].DirectionNo = N;

//	movement_Array[9].Enable = true;
//	movement_Array[9].Volume = 60;
//	movement_Array[9].StageNo = stage(movement_to_Stage_Array[9]);
//	movement_Array[9].GroupNo = group(2);
//	movement_Array[9].DirectionNo = N;

//	movement_Array[10].Enable = true;
//	movement_Array[10].Volume = 90;
//	movement_Array[10].StageNo = stage(movement_to_Stage_Array[10]);
//	movement_Array[10].GroupNo = group(1);
//	movement_Array[10].DirectionNo = S;

//	movement_Array[11].Enable = true;
//	movement_Array[11].Volume = 340;
//	movement_Array[11].StageNo = stage(movement_to_Stage_Array[11]);
//	movement_Array[11].GroupNo = group(2);
//	movement_Array[11].DirectionNo = S;

//	movement_Array[12].Enable = true;
//	movement_Array[12].Volume = 50;
//	movement_Array[12].StageNo = stage(movement_to_Stage_Array[12]);
//	movement_Array[12].GroupNo = group(2);
//	movement_Array[12].DirectionNo = S;

//}

};

std::map<int, CSignalNode> g_signal_node_map;  // first key is signal node id

vector<float> g_time_parser(vector<string>& inputstring)
{
	vector<float> output_global_minute;

	for (int k = 0; k < inputstring.size(); k++)
	{
		vector<string> sub_string = split(inputstring[k], "_");

		for (int i = 0; i < sub_string.size(); i++)
		{
			//HHMM
			//012345
			char hh1 = sub_string[i].at(0);
			char hh2 = sub_string[i].at(1);
			char mm1 = sub_string[i].at(2);
			char mm2 = sub_string[i].at(3);

			float hhf1 = ((float)hh1 - 48);
			float hhf2 = ((float)hh2 - 48);
			float mmf1 = ((float)mm1 - 48);
			float mmf2 = ((float)mm2 - 48);

			float hh = hhf1 * 10 * 60 + hhf2 * 60;
			float mm = mmf1 * 10 + mmf2;
			float global_mm_temp = hh + mm;
			output_global_minute.push_back(global_mm_temp);
		}
	}

	return output_global_minute;
} // transform hhmm to minutes 


void g_sig_ReadInputData(CMainSigModual& MainSigModual)
{
	if (MainSigModual.b_with_loaded_data == true)
		return;

	//step 0:read demand period file
	MainSigModual.g_LoadingStartTimeInMin = assignment.g_LoadingStartTimeInMin;
	MainSigModual.g_LoadingEndTimeInMin = assignment.g_LoadingEndTimeInMin;


	MainSigModual.g_number_of_nodes = 0;
	MainSigModual.g_number_of_links = 0;  // initialize  the counter to 0


	int internal_node_seq_no = 0;
	// step 3: read node file 

	CCSVParser parser;
	if (parser.OpenCSVFile("node.csv", true))
	{
		while (parser.ReadRecord())  // if this line contains [] mark, then we will also read field headers.
		{

			int node_id;

			if (parser.GetValueByFieldName("node_id", node_id) == false)
				continue;

			if (MainSigModual.g_internal_node_to_seq_no_map.find(node_id) != MainSigModual.g_internal_node_to_seq_no_map.end())
			{
				continue; //has been defined
			}
			MainSigModual.g_internal_node_to_seq_no_map[node_id] = internal_node_seq_no;



			CNode node;  // create a node object 

			node.node_id = node_id;
			node.node_seq_no = internal_node_seq_no;

			internal_node_seq_no++;

			MainSigModual.m_sig_node_vector.push_back(node);  // push it to the global node vector

			MainSigModual.g_number_of_nodes++;
			if (MainSigModual.g_number_of_nodes % 5000 == 0)
			{
				log_out << "reading " << MainSigModual.g_number_of_nodes << " nodes.. " << endl;
			}
		}
		g_info_String = "Number of Nodes = ";
		g_info_String.append(to_string(MainSigModual.g_number_of_nodes));
		MainSigModual.WriteLog(0, g_info_String, 2);

	//	fprintf(g_pFileOutputLog, "number of nodes =,%d\n", MainSigModual.g_number_of_nodes);

		parser.CloseCSVFile();
	}

	// step 4: read link file 

	CCSVParser parser_link;
	if (parser_link.OpenCSVFile("link.csv", true))
	{
		while (parser_link.ReadRecord())  // if this line contains [] mark, then we will also read field headers.
		{
			int from_node_id;
			int to_node_id;
			if (parser_link.GetValueByFieldName("from_node_id", from_node_id) == false)
				continue;
			if (parser_link.GetValueByFieldName("to_node_id", to_node_id) == false)
				continue;

			string linkID;
			parser_link.GetValueByFieldName("link_id", linkID);


			// add the to node id into the outbound (adjacent) node list

			if (MainSigModual.g_internal_node_to_seq_no_map.find(from_node_id) == MainSigModual.g_internal_node_to_seq_no_map.end())
			{
				log_out << "Error: from_node_id " << from_node_id << " in file link.csv is not defined in node.csv." << endl;

				continue; //has not been defined
			}
			if (MainSigModual.g_internal_node_to_seq_no_map.find(to_node_id) == MainSigModual.g_internal_node_to_seq_no_map.end())
			{
				log_out << "Error: to_node_id " << to_node_id << " in file link.csv is not defined in node.csv." << endl;
				continue; //has not been defined
			}

			if (MainSigModual.g_link_id_map.find(linkID) != MainSigModual.g_link_id_map.end())
			{
				log_out << "Error: link_id " << linkID.c_str() << " has been defined more than once. Please check link.csv." << endl;
				continue; //has not been defined
			}


			int internal_from_node_seq_no = MainSigModual.g_internal_node_to_seq_no_map[from_node_id];  // map external node number to internal node seq no. 
			int internal_to_node_seq_no = MainSigModual.g_internal_node_to_seq_no_map[to_node_id];


			CLink link;  // create a link object 

			link.from_node_seq_no = internal_from_node_seq_no;
			link.to_node_seq_no = internal_to_node_seq_no;
			link.link_seq_no = MainSigModual.g_number_of_links;
			link.to_node_seq_no = internal_to_node_seq_no;
			link.link_id = linkID;

			MainSigModual.g_link_id_map[link.link_id] = 1;



			parser_link.GetValueByFieldName("link_type", link.link_type);
			parser_link.GetValueByFieldName("geometry", link.geometry, false);

			float length = 1.0; // km or mile
			float free_speed = 1.0;

			float lane_capacity = 1800;
			parser_link.GetValueByFieldName("length", length);
			parser_link.GetValueByFieldName("free_speed", free_speed);
			free_speed = max(0.1, free_speed);

			int number_of_lanes = 1;
			parser_link.GetValueByFieldName("lanes", number_of_lanes);
			parser_link.GetValueByFieldName("capacity", lane_capacity);

			float default_cap = 1000;
			float default_BaseTT = 1;

			link.free_flow_travel_time_in_min = length / free_speed * 60;
			link.number_of_lanes = number_of_lanes;
			link.lane_capacity = lane_capacity;


			link.length = length;
			link.free_flow_travel_time_in_min = length / free_speed * 60;


			MainSigModual.m_sig_node_vector[internal_from_node_seq_no].m_outgoing_link_seq_no_vector.push_back(link.link_seq_no);  // add this link to the corresponding node as part of outgoing node/link
			MainSigModual.m_sig_node_vector[internal_to_node_seq_no].m_incoming_link_seq_no_vector.push_back(link.link_seq_no);  // add this link to the corresponding node as part of outgoing node/link



			MainSigModual.m_sig_node_vector[internal_from_node_seq_no].m_to_node_seq_no_vector.push_back(link.to_node_seq_no);  // add this link to the corresponding node as part of outgoing node/link
			MainSigModual.m_sig_node_vector[internal_from_node_seq_no].m_to_node_2_link_seq_no_map[link.to_node_seq_no] = link.link_seq_no;  // add this link to the corresponding node as part of outgoing node/link

			MainSigModual.m_sig_link_vector.push_back(link);

			MainSigModual.g_number_of_links++;

			// map link data to signal node map.

			string movement_str;
			parser_link.GetValueByFieldName("movement_str", movement_str,false);


			if (movement_str.size() > 0)  // and valid
			{
				int main_node_id = -1;

				parser_link.GetValueByFieldName("main_node_id", main_node_id, true, false);


				int NEMA_phase_number = 0;
				parser_link.GetValueByFieldName("NEMA_phase_number", NEMA_phase_number, true, false);

				// ini_stage_number
				// ini_ min green time and max green time and extention flag

				int lanes = 0;
				parser_link.GetValueByFieldName("lanes", lanes);

				int sharedLanes = 0;
				parser_link.GetValueByFieldName("sharedLanes", sharedLanes, false,false);


				if (main_node_id >= 1)
				{
					g_signal_node_map[main_node_id].AddMovementStructure(link.link_seq_no, movement_str, lanes, sharedLanes, linkID);
				}
			}
		}
	}
		parser_link.CloseCSVFile();
	// we now know the number of links
		g_info_String = "Number of Links = ";
		g_info_String.append(to_string(MainSigModual.g_number_of_links));
		MainSigModual.WriteLog(0, g_info_String, 2);

		MainSigModual.b_with_loaded_data = true;
};


double SignalAPI(int iteration_number, int MainSigModual_mode, int signal_updating_output)
{
	// step 1: read input data of network 

	MainSigModual.signal_updating_output = signal_updating_output;

	g_sig_ReadInputData(MainSigModual);

	for (std::map<int, CSignalNode>::iterator it = g_signal_node_map.begin(); it != g_signal_node_map.end(); ++it)
	{
		g_info_String = "Start Signal Main Node ID ";
		g_info_String.append(to_string(it->first));
		g_info_String.append(":");
		MainSigModual.WriteLog(0, g_info_String, -1);

		it->second.PerformQEM(it->first);

		g_info_String = "End Signal Main Node ID ";
		g_info_String.append(to_string(it->first));
		MainSigModual.WriteLog(0, g_info_String, -2);

	}

	//output timing.csv

	FILE* g_pFileTimingArc = NULL;

	fopen_ss(&g_pFileTimingArc, "timing.csv", "w");

	if (g_pFileTimingArc == NULL)
	{
		log_out << "File timing.csv cannot be opened." << endl;
		return -1;
	}
	else
	{
		fprintf(g_pFileTimingArc, "link_id,from_node_id,to_node_id,time_window,time_interval,travel_time_delta,capacity,cycle_no,cycle_length,green_time,red_time,main_node_id,stage,movement_str,notes,geometry\n");

		for (std::map<int, CSignalNode>::iterator it = g_signal_node_map.begin(); it != g_signal_node_map.end(); ++it)
		{
			CSignalNode sn = it->second;

			log_out << "Signal timing updating for main node id " << it->first << ":" << endl;

			int cycle_time_in_sec = max(10, sn.c_Final + 0.5);
			int number_of_cycles = (MainSigModual.g_LoadingEndTimeInMin - MainSigModual.g_LoadingStartTimeInMin) * 60 / cycle_time_in_sec;  // unit: seconds;
			int offset_in_sec = 0;
			int g_loading_start_time_in_sec = MainSigModual.g_LoadingStartTimeInMin * 60 + offset_in_sec;
			int ci = 0;
			//for (int ci = 0; ci <= number_of_cycles; ci++)
			{

				for (int m = 1; m < movementSize; m++)
				{
					if (sn.movement_Array[m].Enable)
					{
						for (size_t so = 0; so < sn.movement_Array[m].StageNo_in_Order.size(); so++)
						{
							int StageNo = sn.movement_Array[m].StageNo_in_Order[so];
							// we should also consider offset.
							int global_start_time_in_sec = sn.cumulative_Green_Start_Time_Stage_Array[StageNo] + cycle_time_in_sec * ci + g_loading_start_time_in_sec;
							int global_end_time_in_sec = sn.cumulative_Green_End_Time_Stage_Array[StageNo] + cycle_time_in_sec * ci + g_loading_start_time_in_sec;

							//0300:30
							int start_hour = global_start_time_in_sec / 3600;
							int start_min = global_start_time_in_sec / 60 - start_hour * 60;
							int start_sec = global_start_time_in_sec % 60;

							int end_hour = global_end_time_in_sec / 3600;
							int end_min = global_end_time_in_sec / 60 - end_hour * 60;
							int end_sec = global_end_time_in_sec % 60;

							int from_node_id = MainSigModual.m_sig_node_vector[MainSigModual.m_sig_link_vector[sn.movement_Array[m].LinkSeqNo].from_node_seq_no].node_id;
							int to_node_id = MainSigModual.m_sig_node_vector[MainSigModual.m_sig_link_vector[sn.movement_Array[m].LinkSeqNo].to_node_seq_no].node_id;
							//						float capacity = sn.green_Time_Stage_Array[StageNo] * sn.saturation_Flow_Rate_Matrix[StageNo][m] / 3600.0;
							float capacity = sn.green_Time_Stage_Array[StageNo] * 1800.0 / 3600.0;
							//float capacity = sn.capacity_by_Stage_and_Movement_Matrix[StageNo][m]/60;

							float greenTime = sn.green_Time_Stage_Array[StageNo];

							float redTime = cycle_time_in_sec - sn.green_Time_Stage_Array[StageNo];

							fprintf(g_pFileTimingArc, "%s,%d,%d,%02d%02d:%02d_%02d%02d:%02d,-1,-1,%f,%d,%d,%f,%f,%d,%d,%s,%d;%d,\"%s\"\n",
								MainSigModual.m_sig_link_vector[sn.movement_Array[m].LinkSeqNo].link_id.c_str(),
								from_node_id,
								to_node_id,
								start_hour,
								start_min,
								start_sec,
								end_hour,
								end_min,
								end_sec,
								capacity,
								ci,
								cycle_time_in_sec,
								greenTime,
								redTime,
								it->first,
								StageNo,
								sn.movement_str_array[m].c_str(),
								m, so,
								MainSigModual.m_sig_link_vector[sn.movement_Array[m].LinkSeqNo].geometry.c_str()
							);
						}

					}  // per movement

				}
			}  // per cycle
		}  // per signal node

		fclose(g_pFileTimingArc);
	}

	MainSigModual.WriteLog(0, "-----------------------Finished----------------------- ", 1);

	return 1;
}