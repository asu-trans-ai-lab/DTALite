//#ifdef _WIN64
#include "pch.h"
//#endif

#ifndef _WIN64  //linux
#define NULL 0
#endif

#include "MathLibrary.h"

double add_c(int a, int b)  // for testing
{
	double c = 0.0;

	c += a * b;

	for (int i = 0; i < 100000000; i++)
	{
		c += 0.001;
	}
	return c;
}



double shortest_path(int node_size, int link_size,
	int* from_node_no_vector, int* to_node_no_vector, double* cost_vector,
	int* FirstLinkFrom, int* LastLinkFrom, int* sorted_link_no_vector,
	int o_node_no, int d_node_no, 
	int* node_pred, int* link_pred, int* SEList_queue_next_vector, double* label_cost_vector)
{
// The following deque based implementation is motivated and adpated by the efficient implementiation by Dr. Hillel Bar-Gera
	// similar implementation can be also found at DYNASMART system design by Dr. Hani Mahmassani and original code in DTALite by Dr. Xuesong Zhou
	//http://www.bgu.ac.il/~bargera/tntp/
	//http://www.bgu.ac.il/~bargera/tntp/FW.zip

	double cost = 0;
	int node_indiex, current_node, NewNode, k, SEList_Returnwith_Q_Count = 0;
	float NewCost;
	int INVALID = -1;
	int SEList_WAS_IN_QUEUE_Flag = -7;
	int SEList_queue_first, SEList_queue_last;
	int FirstThruNode = 0;  // used t filter out the TAZ based centriods

	//SEList_queue_next_vector is the implementation of scan eligible list for active nodes in label correcting 

	for (node_indiex = 0; node_indiex < node_size; node_indiex++) {
		SEList_queue_next_vector[node_indiex] = INVALID;  // scan eligible list
		label_cost_vector[node_indiex] = 1.0e+15;  // label cost
		link_pred[node_indiex] = INVALID;
		node_pred[node_indiex] = INVALID;
	}

	//initialization 
	current_node = o_node_no;
	SEList_queue_next_vector[current_node] = SEList_WAS_IN_QUEUE_Flag;
	link_pred[current_node] = INVALID;
	node_pred[node_indiex] = INVALID;
	label_cost_vector[current_node] = 0.0;

	//SEList initialization
	SEList_queue_first = SEList_queue_last = INVALID;

	//label correction scanning
	while ((current_node != INVALID) && (current_node != SEList_WAS_IN_QUEUE_Flag))
	{
		if (current_node >= FirstThruNode || current_node == o_node_no) {

			//cout << "Scan node_indiex " << current_node << " with " << " node_indiex cost " << label_cost_vector[current_node]  <<  endl;
//scan all outgoing link
//		FirstLinkFrom, int* LastLinkFrom, int* sorted_link_no_vector,
//			for (int i = 0; i < NodeForwardStarArray[current_node].OutgoingLinkSize; i++)  // for each link (i,j) belong A(i)
			for (k = FirstLinkFrom[current_node]; k < LastLinkFrom[current_node]; k++) 
		{
				int link_seq_no = sorted_link_no_vector[k];
				NewNode = to_node_no_vector[link_seq_no];

				NewCost = label_cost_vector[current_node] + cost_vector[link_seq_no];

				if (label_cost_vector[NewNode] > NewCost)
				{
					label_cost_vector[NewNode] = NewCost;
					link_pred[NewNode] = link_seq_no;
					node_pred[NewNode] = from_node_no_vector[link_seq_no];

/* If the new node_indiex was in the queue before, add it as the first in the queue. */
					if (SEList_queue_next_vector[NewNode] == SEList_WAS_IN_QUEUE_Flag) {
						SEList_queue_next_vector[NewNode] = SEList_queue_first;
						SEList_queue_first = NewNode;

						if (SEList_queue_last == INVALID)
							SEList_queue_last = NewNode;
						SEList_Returnwith_Q_Count++;

					}

					/* If the new node_indiex is not in the queue, and wasn't there before, add it at the end of the queue. */
					else if (SEList_queue_next_vector[NewNode] == INVALID && NewNode != SEList_queue_last) 
					{
						if (SEList_queue_last != INVALID) { 					/*Usually*/
							SEList_queue_next_vector[SEList_queue_last] = NewNode;
							SEList_queue_last = NewNode;
						}
						else {			  /* If the queue is empty, initialize it. */
							SEList_queue_first = SEList_queue_last = NewNode;
							SEList_queue_next_vector[SEList_queue_last] = INVALID;
						}

					}

					/* If the new node_indiex is in the queue, just leave it there. (Do nothing) */
				}
			}
		}

		/* Get the first node_indiex out of the queue, and use it as the current node_indiex. */
		current_node = SEList_queue_first;
		if ((current_node == INVALID) || (current_node == SEList_WAS_IN_QUEUE_Flag))
			break;

		SEList_queue_first = SEList_queue_next_vector[current_node];
		SEList_queue_next_vector[current_node] = SEList_WAS_IN_QUEUE_Flag;
		if (SEList_queue_last == current_node)
			SEList_queue_last = INVALID;
	}
	return cost;
}

