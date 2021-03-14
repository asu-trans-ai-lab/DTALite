#ifndef GUARD_MATHLIBRARY_H

#ifdef _WIN32
#define MATHLIBRARY_API __declspec(dllexport)
#else
#define MATHLIBRARY_API
#endif

extern "C" MATHLIBRARY_API double network_assignment(int iteration_number, int assignment_mode, int column_updating_iterations);

#endif