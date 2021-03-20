#ifndef GUARD_CONFIG_H

#ifdef BUILD_EXE
    double network_assignment(int assignment_mode, int iteration_number, int column_updating_iterations);
#else
    #ifdef _WIN32
        #define DTALIBRARY_API __declspec(dllexport)
    #else
        #define DTALIBRARY_API
    #endif

    extern "C" DTALIBRARY_API double network_assignment(int assignment_mode, int iteration_number, int column_updating_iterations);
#endif

#endif