//  Portions Copyright 2019
// Xuesong Zhou
//   If you help write or modify the code, please also list your names here.
//   The reason of having Copyright info here is to ensure all the modified version, as a whole, under the GPL
//   and further prevent a violation of the GPL.

// More about "How to use GNU licenses for your own software"
// http://www.gnu.org/licenses/gpl-howto.html

#include "stdafx.h"


#include "C:\SourceCode\DTALite_DLL\DTALite\src\Exe_src\AgentLite\main_api.cpp"


int main(int argc, TCHAR* argv[], TCHAR* envp[])
{
	int assignment_mode = 0;
	int iteration_number = 20;
	int column_updating_iterations = 20;
	network_assignment(iteration_number, assignment_mode, column_updating_iterations);

}
