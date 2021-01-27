//  Portions Copyright 2019
// Xuesong Zhou
//   If you help write or modify the code, please also list your names here.
//   The reason of having Copyright info here is to ensure all the modified version, as a whole, under the GPL
//   and further prevent a violation of the GPL.

// More about "How to use GNU licenses for your own software"
// http://www.gnu.org/licenses/gpl-howto.html

#include "stdafx.h"

// beginning of log file
//a class that simultaneously outputs to ostream and ofstream objects 
//// http://www.cplusplus.com/forum/general/64174/#msg347154 
#include <iostream>

template < typename C, typename T = std::char_traits<C> >

struct basic_teebuf : public std::basic_streambuf<C, T>
{
	typedef std::basic_streambuf<C, T> streambuf_type;
	typedef typename T::int_type int_type;

	basic_teebuf(streambuf_type* buff_a, streambuf_type* buff_b)
		: first(buff_a), second(buff_b) {}

protected:
	virtual int_type overflow(int_type c)
	{
		const int_type eof = T::eof();
		if (T::eq_int_type(c, eof)) return T::not_eof(c);
		else
		{
			const C ch = T::to_char_type(c);
			if (T::eq_int_type(first->sputc(ch), eof) ||
				T::eq_int_type(second->sputc(ch), eof))
				return eof;
			else return c;
		}
	}

	virtual int sync()
	{
		return !first->pubsync() && !second->pubsync() ? 0 : -1;
	}

private:
	streambuf_type* first;
	streambuf_type* second;
};

template < typename C, typename T = std::char_traits<C> >
struct basic_teestream : public std::basic_ostream<C, T>
{
	// add more controls of the debug logs
	//
	int debug_level;
	int log_sig;
	int log_odme;
	int log_path;
	int log_dta;
	int log_ue;


	typedef std::basic_ostream<C, T> stream_type;
	typedef basic_teebuf<C, T> streambuff_type;

	basic_teestream(stream_type& first, stream_type& second)
		: stream_type(&stmbuf), stmbuf(first.rdbuf(), second.rdbuf()) {}

	basic_teestream(streambuff_type* first, streambuff_type* second)
		: stream_type(&stmbuf), stmbuf(first, second) {}

	~basic_teestream() { stmbuf.pubsync(); }

private: streambuff_type stmbuf;
};

typedef basic_teebuf<char> teebuf;
typedef basic_teestream<char> teestream;

#include <fstream>
#include <iomanip>
#include <algorithm>
#include <iterator>

// end of log file

extern teestream log_out;
extern double SignalAPI(int iteration_number, int MainSigModual_mode, int signal_updating_output);
#include "..\AgentLite\main_api.cpp";
#include "..\AgentLite\signal_api.cpp";

std::ofstream file("log.csv");
teestream log_out(std::cout, file);

int main(int argc, TCHAR* argv[], TCHAR* envp[])
{

	// reset all the log files to defult 0: not output; if want to output these logs set to 1
	log_out << "STALite Log" << std::fixed << std::setw(12) << '\n';
	log_out.debug_level = 0;
	log_out.log_sig = 0;
	log_out.log_odme = 0;
	log_out.log_path = 0;
	log_out.log_dta = 0;
	log_out.log_ue = 0;

	int iteration_number = 20;
	int column_updating_iterations = 40;
	int signal_updating_iterations = -1;
	int signal_updating_output = 0;
	int assignment_mode = 1;  // generate link performance and agent file
	bool flag_default = false;
	int default_volume = 1;



	CCSVParser parser_settings;
	parser_settings.IsFirstLineHeader = false;

	if (parser_settings.OpenCSVFile("settings.csv", false))
	{
		while (parser_settings.ReadRecord_Section())
		{
			if (parser_settings.SectionName == "[assignment]")
			{

				int value_int;
				string assignment_mode_str;
				parser_settings.GetValueByFieldName("number_of_iterations", iteration_number, true, true);
				parser_settings.GetValueByFieldName("assignment_mode", assignment_mode_str);
				// these are the assignment modes
				// two usually methods are ue (user equilibrium) and dta (dynamic traffic assignment)
				// the main difference of these two methods are different output in link_performance.csv
				// for basic uses set assignment mode to 'ue'
				// for more detailed link performances (one minute) set 'dta'
				if (assignment_mode_str == "lue")
					assignment_mode = 0;
				else if (assignment_mode_str == "ue")
					assignment_mode = 1;
				else if (assignment_mode_str == "dta")
					assignment_mode = 2;
				else if (assignment_mode_str == "odme")
					assignment_mode = 3;
				else if (assignment_mode_str == "sig")
					assignment_mode = 4;
				else if (assignment_mode_str == "gravity")
					assignment_mode = 5;
				else
				{
					log_out << "assignment_mode " << assignment_mode_str.c_str() << " in settings.csv is invalid." << endl;
					g_ProgramStop();
				}

				// iteration number of reassignment
				parser_settings.GetValueByFieldName("column_updating_iterations", column_updating_iterations, true, true);

				// the start interation of generating signals, if there is no signals set this number larger than the iteration number 
				parser_settings.GetValueByFieldName("signal_updating_iterations", signal_updating_iterations, true, false);
				

				break;  // just one record
			}

			if (parser_settings.SectionName == "[log]")
			{


				parser_settings.GetValueByFieldName("sig", log_out.log_sig,false);
				parser_settings.GetValueByFieldName("odme", log_out.log_odme, false);
				parser_settings.GetValueByFieldName("path", log_out.log_path, false);
				parser_settings.GetValueByFieldName("ue", log_out.log_ue, false);

				break;  // just one record
			}
		}
	}

	// only when assignement equals to 4, will generate the signal timing based on the link volume and capacity
	if (assignment_mode == 4)
	{
		SignalAPI(0, assignment_mode, 1);
		return 0;
	}


	network_assignment(iteration_number, assignment_mode, column_updating_iterations, signal_updating_iterations);  // obtain initial flow values

}
