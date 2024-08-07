/* Portions Copyright 2021 Xuesong Zhou and Peiheng Li
 *
 * If you help write or modify the code, please also list your names here.
 * The reason of having Copyright info here is to ensure all the modified version, as a whole, under
 * the GPL and further prevent a violation of the GPL.
 *
 * More about "How to use GNU licenses for your own software"
 * http://www.gnu.org/licenses/gpl-howto.html
 */

#ifndef GUARD_CONFIG_H
#define GUARD_CONFIG_H

#include <build_config.h>

#ifdef BUILD_EXE
double network_assignment(int assignment_mode,
                          int iteration_number,
                          int column_updating_iterations,
                          int ODME_iterations,
                          int sensitivity_analysis_iterations,
                          int simulation_output,
                          int number_of_memory_blocks,
                          int length_unit_flag,
                          int speed_unit_flag,
                          double UE_convergence_percentage,
                          int max_num_significant_zones_in_subarea,
                          int max_num_significant_zones_outside_subarea);
#else
#ifdef _WIN32
#define DTALIBRARY_API __declspec(dllexport)
#else
#define DTALIBRARY_API
#endif

extern "C" DTALIBRARY_API double network_assignment(int assignment_mode,
                                                    int iteration_number,
                                                    int column_updating_iterations,
                                                    int ODME_iterations,
                                                    int sensitivity_analysis_iterations,
                                                    int simulation_output,
                                                    int number_of_memory_blocks,
                                                    int length_unit_flag,
                                                    int speed_unit_flag,
                                                    double UE_convergence_percentage,
                                                    int max_num_significant_zones_in_subarea,
                                                    int max_num_significant_zones_outside_subarea);
extern "C" DTALIBRARY_API void DTALiteAPI();
#endif

#endif