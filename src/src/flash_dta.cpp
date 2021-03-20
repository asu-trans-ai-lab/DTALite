/* Portions Copyright 2019 Xuesong Zhou
 *
 * If you help write or modify the code, please also list your names here.
 * The reason of having Copyright info here is to ensure all the modified version, as a whole, under the GPL
 * and further prevent a violation of the GPL.
 *
 * More about "How to use GNU licenses for your own software"
 * http://www.gnu.org/licenses/gpl-howto.html 
 */
#ifdef _WIN32
#include "pch.h"
#endif
#include <iostream>
#include <fstream>
#include <iomanip>
#include "teestream.h"
#include "utils.h"
#include "config.h"

int main()
{
    std::ofstream file("log.csv");
    teestream log_out(std::cout, file);
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
    // generate link performance and agent file
    int assignment_mode = 1;
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
                std::string assignment_mode_str;
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
                else
                {
                    log_out << "assignment_mode " << assignment_mode_str.c_str() << " in settings.csv is invalid." << std::endl;
                    g_ProgramStop();
                }

                // iteration number of reassignment
                parser_settings.GetValueByFieldName("column_updating_iterations", column_updating_iterations, true, true);

                // the start interation of generating signals, if there is no signals set this number larger than the iteration number 
                parser_settings.GetValueByFieldName("signal_updating_iterations", signal_updating_iterations, true, false);
                
                // just one record
                break;
            }

            if (parser_settings.SectionName == "[log]")
            {
                parser_settings.GetValueByFieldName("sig", log_out.log_sig,false);
                parser_settings.GetValueByFieldName("odme", log_out.log_odme, false);
                parser_settings.GetValueByFieldName("path", log_out.log_path, false);
                parser_settings.GetValueByFieldName("ue", log_out.log_ue, false);
                
                // just one record
                break; 
            }
        }
    }
    // obtain initial flow values
    network_assignment(assignment_mode, iteration_number, column_updating_iterations); 

    return 0;
}