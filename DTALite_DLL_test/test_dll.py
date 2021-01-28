import time
import os.path
import ctypes
from sys import platform

iteration_number = 20
assignment_mode = 1
# 0: Link based UE, only produces link performance file without agent path file 
# 1: Path based UE, produces link performance file and agent path file 
# 2: UE + dynamic traffic assignment and simulation , produces link performance file and agent path file 
# 3: ODME

column_generation_number = 20
# with time-dependent PointQueueX approximation

def perform_DTA():

    cwd_ = os.path.abspath(__file__)

    if platform.startswith('win32'):
        dll_file = os.path.join(os.path.dirname(cwd_), 'DTALite.dll')
    elif platform.startswith('linux'):
        dll_file = os.path.join(os.path.dirname(cwd_), 'DTALite.so')
    else:
        raise Exception('Please build the shared library compatible to your OS\
                        using source files')
                    
    cdll = ctypes.cdll.LoadLibrary(dll_file)

    time_start = time.time()

    network_compu = cdll.network_assignment
    # add_fun.argtypes=[ctypes.c_float,ctypes.c_float]
    network_compu.restype = ctypes.c_double
    d = network_compu(iteration_number,assignment_mode, column_generation_number)

    time_end = time.time()
    print('DTALite')
    print('total time: {}'.format(time_end-time_start))


if __name__ == "__main__":
    perform_DTA()