#pragma once

#ifdef _WIN32
# define MATHLIBRARY_API __declspec(dllexport)
#else
#define MATHLIBRARY_API
#endif

// extern "C" MATHLIBRARY_API double network_assignment(int a, int b);

extern "C" MATHLIBRARY_API void shortest_path(int node_size, int link_size,
											  int* from_node_id_vector, int* to_node_id_vector, double* cost_vector,
											  int* FirstLinkFrom, int* LastLinkFrom, int* sorted_link_no_vector, 
											  int o_node_id, int d_node_id, 
											  int* node_pred, int* link_pred, 
											  int* QueueNext, double* CostTo);