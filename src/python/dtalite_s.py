import csv
import time
import numpy 
import operator
from random import choice
import ctypes
import collections
import heapq


# =========================
# -*- data block define -*- 
# =========================
MAX_LABEL_COST_IN_SHORTEST_PATH = 10000

# ==========================================
# -*- data block simulation time horizon -*-
# ==========================================
# the length of simulation time (min)
LENGTH_OF_SIMULATION_TIME_HORIZON_IN_MIN = 90 
# the number of seconds per simulation interval  
NUMBER_OF_SECONDS_PER_SIMU_INTERVAL = 6  
LENGTH_OF_SIMULATION_TIME_HORIZON_IN_INTERVAL = int(
    LENGTH_OF_SIMULATION_TIME_HORIZON_IN_MIN * 60 / NUMBER_OF_SECONDS_PER_SIMU_INTERVAL
)
NUMBER_OF_SIMU_INTERVALS_PER_MIN = int(60 / NUMBER_OF_SECONDS_PER_SIMU_INTERVAL)

# ===========================
# -*- simulation settings -*-
# ===========================
# the number of assignment iterations: 1 as default and we can set it's value
NUMBER_OF_ASSIGNMENT_ITERATIONS = 1
# 1: assignment only, 2: assignment + simulation
g_modeling_method = 1 

# ================================================
# -*- data block start and end time of horizon -*-
# ================================================
# start time of the simulation: 9999 as default value to be updated in reading 
# function for agent file
g_simulation_start_time_in_min = 9999
# end time of the simulation:0 as default value to be updated in reading 
# function for agent file
g_simulation_end_time_in_min = 0
# the number of start simualation interval 
g_start_simu_interval_no = 0
# the number of end simualation interval 
g_end_simu_interval_no = 0

# ========================================
# -*- data block input data statistics -*-
# ========================================
g_number_of_nodes = 0
g_number_of_links = 0
g_number_of_agents = 0

# ===============================================
# -*- data block simulation result statistics -*-
# ===============================================
g_cumulative_arrival_count = 0
g_cumulative_departure_count = 0


class Node:         
    """ external_node_id: the id of node
    node_seq_no: the index of the node and we call the node by its index
    
    we use g_internal_node_seq_no_dict(id to index) and 
    g_external_node_id_dict(index to id) to map them
    """

    def __init__(self, node_seq_no, external_node_id, zone_id): 
        """ the attribute of node  """ 
        self.node_seq_no = node_seq_no
        self.external_node_id = int(external_node_id)
        self.outgoing_link_list = list()
        self.incoming_link_list = list()
        if len(zone_id) == 0:
            self.zone_id = -1
        else:    
            self.zone_id = int(zone_id)
        

class Link:

    def __init__(self, link_seq_no, from_node_no, to_node_no, 
                 from_node_id, to_node_id, length, lanes,
                 free_speed, capacity, link_type, VDF_alpha, VDF_beta):   
        """ the attribute of link """
        self.link_seq_no = link_seq_no
        self.from_node_seq_no = from_node_no
        self.to_node_seq_no = to_node_no
        self.external_from_node = int(from_node_id)
        self.external_to_node = int(to_node_id)
        # 1:one direction 2:two way
        self.type = int(link_type)
        self.lanes = int(lanes)
        self.BPR_alpha = float(VDF_alpha)
        self.BPR_beta = float(VDF_beta)
        self.flow_volume = 0
        # capacity is lane capacity per hour
        self.link_capacity = float(capacity) * int(lanes)
        # length is mile or km
        self.length = float(length) 
        # length:km, free_speed: km/h
        self.free_flow_travel_time_in_min = self.length / max(0.001,int(free_speed)) * 60  
        self.cost = self.free_flow_travel_time_in_min
        # the capacity of discharging for each link,  
        # td is time-dependent to be used in traffic signal control 
        # or work zone scheduling applications
        # self.link_capacity/g_number_of_simulation_time) should be per interval capacity 
        self.td_link_out_flow_capacity = [int(self.link_capacity/(60*NUMBER_OF_SIMU_INTERVALS_PER_MIN))]*(LENGTH_OF_SIMULATION_TIME_HORIZON_IN_INTERVAL+1) 
        # count the cumulative arrival and departure of links 
        # each simulation interval, td is time-dependent
        # [0] is setting up zero as the default 
        self.td_link_cumulative_arrival = [0]*(LENGTH_OF_SIMULATION_TIME_HORIZON_IN_INTERVAL+1) 
        # [0] is setting up zero as the default
        self.td_link_cumulative_departure = [0]*(LENGTH_OF_SIMULATION_TIME_HORIZON_IN_INTERVAL+1)   
        # the time-dependent waiting time, td is time-dependent
        # [0] is setting up zero as the default
        self.td_link_waiting_time = [0]*(LENGTH_OF_SIMULATION_TIME_HORIZON_IN_MIN+1)  
        # link-in queue  of each link
        self.entrance_queue = list()
        # link-out queue  of each link
        self.exit_queue = list()

    def ResetMOE(self):
        """ reset the measures of effectiveness """
        self.cumulative_arrival_count = 0
        self.cumulative_departure_count = 0
        self.cumulative_virtual_delay_count = 0

    def CalculateBPRFunction(self):
        """ update the cost of link """
        # travel time in min
        self.general_travel_time_in_min = self.free_flow_travel_time_in_min*(1 + self.BPR_alpha*((self.flow_volume / max(0.00001, self.link_capacity))**self.BPR_beta))  
        # travel time in simulation interval
        self.general_travel_time_in_simu_interval = self.general_travel_time_in_min * NUMBER_OF_SIMU_INTERVALS_PER_MIN  
        self.cost = self.general_travel_time_in_min
        

class Agent:
    """ 
    comments: agent_id: the id of agent
    agent_seq_no: the index of the agent and we call the agent by its index
    """
    
    def __init__(self, agent_id, agent_type, o_zone_id, d_zone_id):
        """ the attribute of agent """ 
        self.agent_id = agent_id
        # vehicle 
        self.agent_type = agent_type  
        self.o_zone_id = int(o_zone_id) 
        self.d_zone_id = int(d_zone_id)
        self.o_node_id = 0
        self.d_node_id = 0
        self.agent_seq_no = 0
        self.path_node_seq_no_list = list()
        self.path_link_seq_no_list = list()
        self.current_link_seq_no_in_path = 0 
        self.departure_time_in_min = 0
        # Passenger Car Equivalent (PCE) of the agent
        self.PCE_factor = 1  
        self.path_cost = 0
        self.b_generated = False
        self.b_complete_trip = False
        self.departure_time_in_simu_interval = int(self.departure_time_in_min*60/NUMBER_OF_SECONDS_PER_SIMU_INTERVAL + 0.5)
        self.feasible_path_exist_flag = False
        
    # def Initialization(self): 
    #     """ update the number of agents """
    #     global g_number_of_agents
    #     self.agent_seq_no = g_number_of_agents
    #     g_number_of_agents += 1
        
    def Initialize_for_simulation(self): 
        if self.path_node_seq_no_list:
            # if the agent can find its path
            # the arrival time on each link in the agent's path  
            self.veh_link_arrival_time_in_simu_interval = [-1]*(len(self.path_node_seq_no_list)-1)
            # the depature time on each link in the agent's path
            self.veh_link_departure_time_in_simu_interval = [-1]*(len(self.path_node_seq_no_list)-1)
            self.veh_link_arrival_time_in_simu_interval[0] = self.departure_time_in_simu_interval 


class Network:
    
    def __init__(self):
        self.node_list = []
        self.link_list = []
        self.agent_list = []
        self.node_size = 0
        self.link_size = 0
        self.agenet_size = 0
        # key: external node id, value:internal node id
        self.internal_node_seq_no_dict = {}
        # key: internal node id, value:external node id
        self.external_node_id_dict = {}
        # td:time-dependent, key:simulation time interval, 
        # value:agents(list) need to be activated
        self.agent_td_list_dict = {}
        # key: zone id, value: node id list
        self.zone_to_nodes_dict = {}
        self.CAPI_allocated_data = 0
        self.cdll = None

    def allocate(self):
        """ initialization for traffic assignment """ 
        # set up the return type and argument types for the shortest path 
        # function in dll.
        self.cdll = ctypes.cdll.LoadLibrary(r"../../lib/libstalite.dll")
        self.cdll.shortest_path.argtypes = [
                ctypes.c_int, ctypes.c_int, 
                numpy.ctypeslib.ndpointer(dtype=numpy.int32),
                numpy.ctypeslib.ndpointer(dtype=numpy.int32),
                numpy.ctypeslib.ndpointer(dtype=numpy.float64),
                numpy.ctypeslib.ndpointer(dtype=numpy.int32),
                numpy.ctypeslib.ndpointer(dtype=numpy.int32),
                numpy.ctypeslib.ndpointer(dtype=numpy.int32),                                        
                ctypes.c_int, ctypes.c_int, 
                numpy.ctypeslib.ndpointer(dtype=numpy.int32),
                numpy.ctypeslib.ndpointer(dtype=numpy.int32),
                numpy.ctypeslib.ndpointer(dtype=numpy.int32),
                numpy.ctypeslib.ndpointer(dtype=numpy.float64)
        ]
    
        # allocate 
        self.node_size = len(self.node_list)
        self.link_size = len(self.link_list)
        self.agenet_size = len(self.agent_list)
        self.link_cost_array = numpy.array(
            [-1]*self.link_size, 
            dtype=numpy.float64
        )
        self.link_volume_array = [0]*self.link_size
        self.node_predecessor = numpy.array(
            [-1]*self.node_size, 
            dtype=numpy.int32
        )
        self.node_label_cost = numpy.array(
            [MAX_LABEL_COST_IN_SHORTEST_PATH]*self.node_size,
            dtype=numpy.float64
        )
        self.link_predecessor = numpy.array(
            [-1]*self.node_size,
            dtype=numpy.int32
        )
        self.queue_next = numpy.array(
            [0]*self.node_size,
            dtype=numpy.int32
        )
        self.from_node_no_array = numpy.array(
            [-1]*self.link_size,
            dtype=numpy.int32
        )
        self.to_node_no_array = numpy.array(
            [-1]*self.link_size,
            dtype=numpy.int32
        )
        self.FirstLinkFrom =numpy.array(
            [-1]*self.node_size,
            dtype=numpy.int32
        )
        self.LastLinkFrom = numpy.array(
            [-1]*self.node_size,
            dtype=numpy.int32
        )
        self.sorted_link_no_vector = numpy.array(
            [-1]*self.link_size,
            dtype=numpy.int32
        )
  
        for j in range(self.link_size):
            self.from_node_no_array[j] = self.link_list[j].from_node_seq_no
            self.to_node_no_array [j] = self.link_list[j].to_node_seq_no
            self.link_cost_array[j] = self.link_list[j].cost
                
        node_OutgoingLinkSize = numpy.array(
            [0]*self.node_size,
            dtype=numpy.int32
        )

        # count the size of outgoing links for each node
        for j in range(self.link_size):
            node_OutgoingLinkSize[self.link_list[j].from_node_seq_no] += 1

        cumulative_count = 0
        for i in range(self.node_size):
            self.FirstLinkFrom[i] = cumulative_count
            self.LastLinkFrom[i] = self.FirstLinkFrom[i] + node_OutgoingLinkSize[i]
            cumulative_count +=  node_OutgoingLinkSize[i]

        # reset the counter # need to construct sorted_link_no_vector
        # we are converting a 2 dimensional dynamic array to a fixed size 
        # one-dimisonal array, with the link size 
        for i in range(self.node_size):
            node_OutgoingLinkSize[i] = 0

        # count again the current size of outgoing links for each node
        for j in range(self.link_size):
            # fetch the curent from node seq no of this link
            from_node_seq_no = self.link_list[j].from_node_seq_no
            # j is the link sequence no in the original link block
            self.sorted_link_no_vector[self.FirstLinkFrom[from_node_seq_no]+node_OutgoingLinkSize[from_node_seq_no]] = j
            # continue to count, increase by 1
            node_OutgoingLinkSize[ self.link_list[j].from_node_seq_no] += 1  

               
    def optimal_label_correcting(self, origin_node, destination_node, departure_time, sp_algm='fifo'):
        """ input : origin_node, destination_node, departure_time
            output : the shortest path
        """
        origin_node = self.internal_node_seq_no_dict[origin_node]
        destination_node = self.internal_node_seq_no_dict[destination_node]
        if not self.node_list[origin_node].outgoing_link_list:
            return 0
        
        # Initialization for all nodes
        self.node_label_cost = [MAX_LABEL_COST_IN_SHORTEST_PATH] * self.node_size
        # pointer to previous node index from the current label at current node
        self.node_predecessor = [-1] * self.node_size
        # pointer to previous node index from the current label at current node
        self.link_predecessor = [-1] * self.node_size
        
        self.node_label_cost[origin_node] = departure_time
        status = [0] * self.node_size

        if sp_algm == 'fifo':
            # scan eligible list
            SEList = []  
            SEList.append(origin_node)

            while SEList:
                from_node = SEList.pop(0)
                status[from_node] = 0
                for k in range(len(self.node_list[from_node].outgoing_link_list)):
                    to_node = self.node_list[from_node].outgoing_link_list[k].to_node_seq_no 
                    new_to_node_cost = self.node_label_cost[from_node] + self.link_cost_array[self.node_list[from_node].outgoing_link_list[k].link_seq_no]
                    # we only compare cost at the downstream node ToID at the new arrival time t
                    if new_to_node_cost < self.node_label_cost[to_node]:
                        # update cost label and node/time predecessor
                        self.node_label_cost[to_node] = new_to_node_cost
                        # pointer to previous physical node index from the current label at current node and time
                        self.node_predecessor[to_node] = from_node 
                        # pointer to previous physical node index from the current label at current node and time
                        self.link_predecessor[to_node] = self.node_list[from_node].outgoing_link_list[k].link_seq_no  
                        if not status[to_node]:
                            SEList.append(to_node)
                            status[to_node] = 1

        elif sp_algm == 'deque':
            SEList = collections.deque()
            SEList.append(origin_node)

            while SEList:
                from_node = SEList.popleft()
                status[from_node] = 2
                for k in range(len(self.node_list[from_node].outgoing_link_list)):
                    to_node = self.node_list[from_node].outgoing_link_list[k].to_node_seq_no 
                    new_to_node_cost = self.node_label_cost[from_node] + self.link_cost_array[self.node_list[from_node].outgoing_link_list[k].link_seq_no]
                    # we only compare cost at the downstream node ToID at the new arrival time t
                    if new_to_node_cost < self.node_label_cost[to_node]:
                        # update cost label and node/time predecessor
                        self.node_label_cost[to_node] = new_to_node_cost
                        # pointer to previous physical node index from the current label at current node and time
                        self.node_predecessor[to_node] = from_node 
                        # pointer to previous physical node index from the current label at current node and time
                        self.link_predecessor[to_node] = self.node_list[from_node].outgoing_link_list[k].link_seq_no  
                        if status[to_node] != 1:
                            if status[to_node] == 2:
                                SEList.appendleft(to_node)
                            else:
                                SEList.append(to_node)
                            status[to_node] = 1

        elif sp_algm == 'dijkstra':
            # scan eligible list
            SEList = []
            heapq.heapify(SEList)
            heapq.heappush(SEList, (self.node_label_cost[origin_node], origin_node))

            while SEList:
                (label_cost, from_node) = heapq.heappop(SEList)
                for k in range(len(self.node_list[from_node].outgoing_link_list)):
                    to_node = self.node_list[from_node].outgoing_link_list[k].to_node_seq_no 
                    new_to_node_cost = label_cost + self.link_cost_array[self.node_list[from_node].outgoing_link_list[k].link_seq_no]
                    # we only compare cost at the downstream node ToID at the new arrival time t
                    if new_to_node_cost < self.node_label_cost[to_node]:
                        # update cost label and node/time predecessor
                        self.node_label_cost[to_node] = new_to_node_cost
                        # pointer to previous physical node index from the current label at current node and time
                        self.node_predecessor[to_node] = from_node 
                        # pointer to previous physical node index from the current label at current node and time
                        self.link_predecessor[to_node] = self.node_list[from_node].outgoing_link_list[k].link_seq_no  
                        heapq.heappush(SEList, (self.node_label_cost[to_node], to_node))
        # end of sp_algm == 'fifo':
        
        if (destination_node >= 0 and self.node_label_cost[destination_node] < MAX_LABEL_COST_IN_SHORTEST_PATH):
            return 1
        else: 
            return -1

    def optimal_label_correcting_CAPI(self, origin_node, destination_node, departure_time):
        """ input : origin_node,destination_node,departure_time
            output : the shortest path
        """
        o_node_no = self.internal_node_seq_no_dict[origin_node]
        d_node_no = self.internal_node_seq_no_dict[destination_node]

        if not self.node_list[o_node_no].outgoing_link_list:
            return 0
        
        self.cdll.shortest_path(self.node_size, 
                                self.link_size, 
                                self.from_node_no_array,
                                self.to_node_no_array,
                                self.link_cost_array,
                                self.FirstLinkFrom,
                                self.LastLinkFrom,
                                self.sorted_link_no_vector, 
                                o_node_no,
                                d_node_no, 
                                self.node_predecessor,
                                self.link_predecessor,
                                self.queue_next,
                                self.node_label_cost)
    
        #print('shortest path')
        #print("path cost is {}".format(node_label_cost[d_node_no]))
        
        if (o_node_no >= 0 and self.node_label_cost[o_node_no] < MAX_LABEL_COST_IN_SHORTEST_PATH):
            return {
                "path_flag": 1,
                "node_label_cost": self.node_label_cost,
                "node_predecessor": self.node_predecessor,
                "link_predecessor": self.link_predecessor
            }
        else: 
            return {"path_flag": -1}

    def find_path_for_agents_CAPI(self, iteration_no):
        """ input : the number of iteration, output : updated link volume """
        for s in range(self.link_size):
            self.link_volume_array[s] = 0

        # step 1: find shortest path if needed 
        for i in range(self.agenet_size):
            residual = i % (iteration_no + 1)
            # no need to compute a new path at this iteration
            # that is, it will reuse the path from the previous iteration, 
            # stored at p_agent->path_link_seq_no_list.
            if (residual != 0):  
                continue         
            # else move to the next line for finding the shortest path 
            self.agent_list[i].path_link_seq_no_list = []
            self.agent_list[i].path_node_seq_no_list = []
            
            # step 2 build SP tree
            return_value = self.optimal_label_correcting_CAPI(
                self.agent_list[i].o_node_id, 
                self.agent_list[i].d_node_id, 
                self.agent_list[i].departure_time_in_min
            )
            
            # step 3 update path
            if (return_value["path_flag"] == -1):
                #print('agent ',i,'can not find destination node')
                continue

            current_node_seq_no = self.internal_node_seq_no_dict[self.agent_list[i].d_node_id]
            self.agent_list[i].path_cost = return_value["node_label_cost"][current_node_seq_no]

            while current_node_seq_no >= 0:
                current_link_seq_no = return_value["link_predecessor"][current_node_seq_no]
                if current_link_seq_no >= 0:
                    self.agent_list[i].path_link_seq_no_list.insert(0,current_link_seq_no)      
                self.agent_list[i].path_node_seq_no_list.insert(0,current_node_seq_no)
                current_node_seq_no = return_value["node_predecessor"][current_node_seq_no]
            
            if self.agent_list[i].path_node_seq_no_list:
                self.agent_list[i].feasible_path_exist_flag = True
        
        # step 2:  scan the shortest path to compute the link volume, 
        for i in range(self.agenet_size):
            #for each link in the path of this agent
            for j in range(len(self.agent_list[i].path_link_seq_no_list)):
                link_seq_no = self.agent_list[i].path_link_seq_no_list[j]
                self.link_volume_array[link_seq_no] += self.agent_list[i].PCE_factor

    def find_path_for_agents_withoutCAPI(self, iteration_no):
        for s in range(self.link_size):
            self.link_volume_array[s] = 0

        # step 1: find shortest path if needed 
        for i in range(self.agenet_size):
            residual = i % (iteration_no + 1)
            # no need to compute a new path at this iteration
            # that is, it will reuse the path from the previous iteration, 
            # stored at p_agent->path_link_seq_no_list.
            if (residual != 0):  
                continue         
            # else move to the next line for finding the shortest path 
            self.agent_list[i].path_link_seq_no_list = []
            self.agent_list[i].path_node_seq_no_list = []
            
            # step 2 buil SP tree
            return_value = self.optimal_label_correcting(
                self.agent_list[i].o_node_id, 
                self.agent_list[i].d_node_id,
                self.agent_list[i].departure_time_in_min,
                'fifo'
            )
           
            # step 3 update path
            if (return_value == -1):
                #print('agent ',i,'can not find destination node')
                continue

            current_node_seq_no = self.internal_node_seq_no_dict[self.agent_list[i].d_node_id]
            self.agent_list[i].path_cost = self.node_label_cost[current_node_seq_no]
            
            while current_node_seq_no >= 0:
                current_link_seq_no = self.link_predecessor[current_node_seq_no]
                if current_link_seq_no >= 0:
                    self.agent_list[i].path_link_seq_no_list.insert(0, current_link_seq_no)      
                self.agent_list[i].path_node_seq_no_list.insert(0, current_node_seq_no)
                current_node_seq_no = self.node_predecessor[current_node_seq_no]

            if self.agent_list[i].path_node_seq_no_list:
                self.agent_list[i].feasible_path_exist_flag = True

        # step 2:  scan the shortest path to compute the link volume, 
        for i in range(self.agenet_size):
            # for each link in the path of this agent
            for j in range(len(self.agent_list[i].path_link_seq_no_list)):
                link_seq_no = self.agent_list[i].path_link_seq_no_list[j]
                self.link_volume_array[link_seq_no] += self.agent_list[i].PCE_factor


def g_A2R_simu_interval(abslute_time_in_simu_intetrval_no): 
    """ convert absolute simulation interval to relative simulation interval
    
    absolute 24_hour_world clock time in simulation interval to relative time 
    in simulation interval from 0
    """
    return abslute_time_in_simu_intetrval_no - g_start_simu_interval_no


def g_R2A_simu_interval(relative_time_in_simu_intetrval_no):
    """ convert relative simulation interval to absolute simulation interval
    
    relative time in simulation interval from 0 to absolute 24_hour_world clock
    time in simulation interval
    """
    return relative_time_in_simu_intetrval_no + g_start_simu_interval_no


def time_stamp_to_HHMMSS(time_in_minutes):
    hours = int(time_in_minutes/60)
    minutes = int(time_in_minutes - hours*60)
    seconds = int((time_in_minutes - hours*60 - minutes)*60)
    return (time_int_to_str(hours) 
            + time_int_to_str(minutes) 
            + ":" 
            + time_int_to_str(seconds))


def time_int_to_str(time_int):
    if time_int < 10:
        return "0" + str(time_int)
    else:
        return str(time_int)


def g_TrafficAssignment(network):    
    print('Finding shortest path for all agents......')

    for i in range(NUMBER_OF_ASSIGNMENT_ITERATIONS):
        print('iteration_no', i, '......')
        for j in range(g_number_of_links):
            network.link_list[j].CalculateBPRFunction()
            network.link_cost_array[j] = network.link_list[j].cost

        # network.find_path_for_agents_withoutCAPI(i)     
        network.find_path_for_agents_CAPI(i)     
                        
        for k in range(g_number_of_links):
            network.link_list[k].flow_volume = 0
            network.link_list[k].flow_volume += network.link_volume_array[k]


def g_TrafficSimulation(node_list, link_list, agent_list, agent_td_list_dict):
    global g_cumulative_arrival_count
    global g_cumulative_departure_count

    link_list_size = len(link_list)
    #initialization for each link
    for li in range(link_list_size):
        link_list[li].ResetMOE()
    #initialization for each agent
    for agent_no in range(len(agent_list)):
        agent_list[agent_no].Initialize_for_simulation()

    # g_active_agent_list=list()
    # current_active_agent_id = 0
    
    for t in range(g_start_simu_interval_no, g_end_simu_interval_no, 1):
        number_of_simu_interval_per_min = 60 / NUMBER_OF_SECONDS_PER_SIMU_INTERVAL
        if t % number_of_simu_interval_per_min == 0:
            print(
                "simu time= ", int(t/number_of_simu_interval_per_min), "min, ", 
                "CA=", g_cumulative_arrival_count, 
                ", CD=", g_cumulative_departure_count
            )

        relative_t = g_A2R_simu_interval(t)
        for li in range(link_list_size):
            link = link_list[li]
            if relative_t >= 1:
                link_list[link.link_seq_no].td_link_cumulative_departure[relative_t] = link_list[link.link_seq_no].td_link_cumulative_departure[relative_t-1]
                link_list[link.link_seq_no].td_link_cumulative_arrival[relative_t] = link_list[link.link_seq_no].td_link_cumulative_arrival[relative_t-1]
  
        if t in agent_td_list_dict.keys():
            # activate the agent td: time dependent
            for j in range(len(agent_td_list_dict[t])):
                agent_no = agent_td_list_dict[t][j]
                agent = agent_list[agent_no]
                if agent.feasible_path_exist_flag:
                    agent.b_generated = True
                    first_link_seq = agent.path_link_seq_no_list[0]
                    link_list[first_link_seq].td_link_cumulative_arrival[relative_t] += 1
                    link_list[first_link_seq].entrance_queue.append(agent.agent_seq_no)
                    g_cumulative_arrival_count += 1
                
        for li in range(link_list_size):
            link = link_list[li]
            # there are agents in the entrance_queue
            while link.entrance_queue:  
                agent_seq = link.entrance_queue[0]
                link.entrance_queue.pop(0)
                link.exit_queue.append(agent_seq)
                agent_list[agent_seq].veh_link_departure_time_in_simu_interval[agent_list[agent_seq].current_link_seq_no_in_path] = (
                    agent_list[agent_seq].veh_link_arrival_time_in_simu_interval[agent_list[agent_seq].current_link_seq_no_in_path] 
                    + link.general_travel_time_in_simu_interval
                )

        for no in range(len(node_list)):
            node = node_list[no]
            incoming_link_list_size = len(node.incoming_link_list)
            for i in range(incoming_link_list_size):
                incoming_link_index = (i + t) % (incoming_link_list_size) 
                # we will start with different first link from the incoming 
                # link list, equal change, regardless of # of lanes 
                # or main line vs. ramp, but one can use service arc, 
                # to control the effective capacity rates, e.g. through a 
                # metered ramp, to allow mainline to use the remaining flow

                link_seq_no = node.incoming_link_list[incoming_link_index].link_seq_no

                #check if the current link has sufficient capacity
                #check if there are agents in exit queue
             
                while link_list[link_seq_no].td_link_out_flow_capacity[relative_t] >= 1 and len(link_list[link_seq_no].exit_queue) >= 1:
                    agent_no = node.incoming_link_list[incoming_link_index].exit_queue[0]
                   
                    if agent_list[agent_no].veh_link_departure_time_in_simu_interval[agent_list[agent_no].current_link_seq_no_in_path] > t:
                        # the future departure time on this link is later than
                        # the current time
                        break        
               
                    if agent_list[agent_no].current_link_seq_no_in_path == len(agent_list[agent_no].path_link_seq_no_list)-1: 
                        # end of the path
                        node.incoming_link_list[incoming_link_index].exit_queue.pop(0)
                        agent_list[agent_no].b_complete_trip = True
                        link_list[link_seq_no].td_link_cumulative_departure[relative_t] += 1
                        g_cumulative_departure_count += 1
                    else:  
                        # not complete the trip. 
                        # move to the next link's entrance queue
                        next_link_seq = agent_list[agent_no].path_link_seq_no_list[agent_list[agent_no].current_link_seq_no_in_path+1]
                        node.incoming_link_list[incoming_link_index].exit_queue.pop(0)
                        link_list[next_link_seq].entrance_queue.append(agent_no)
                        agent_list[agent_no].veh_link_departure_time_in_simu_interval[agent_list[agent_no].current_link_seq_no_in_path] = t
                        agent_list[agent_no].veh_link_arrival_time_in_simu_interval[agent_list[agent_no].current_link_seq_no_in_path+1] = t
                        
                        actual_travel_time = t - agent_list[agent_no].veh_link_arrival_time_in_simu_interval[agent_list[agent_no].current_link_seq_no_in_path]
			            #for each waited vehicle
                        waiting_time = actual_travel_time - link_list[link_seq_no].general_travel_time_in_min
                        temp_relative_time = g_A2R_simu_interval(agent_list[agent_no].veh_link_arrival_time_in_simu_interval[agent_list[agent_no].current_link_seq_no_in_path])
                        time_in_min = int(temp_relative_time/NUMBER_OF_SIMU_INTERVALS_PER_MIN)
                  
                        link_list[link_seq_no].td_link_waiting_time[time_in_min] += waiting_time
                        link_list[link_seq_no].td_link_cumulative_departure[relative_t] += 1
                        link_list[next_link_seq].td_link_cumulative_arrival[relative_t] += 1
                    
                    agent_list[agent_no].current_link_seq_no_in_path += 1
                    link_list[link_seq_no].td_link_out_flow_capacity[relative_t] -= 1 


def g_ReadInputData(node_list, 
                    link_list, 
                    agent_list,
                    internal_node_seq_no_dict,
                    external_node_id_dict,
                    agent_td_list_dict,
                    zone_to_nodes_dict):
    """
    input: node.csv, link.csv and demand.csv
    output: craet the network and agent 
    """
    global g_simulation_start_time_in_min
    global g_simulation_end_time_in_min
    global g_start_simu_interval_no
    global g_end_simu_interval_no
    global g_number_of_nodes
    global g_number_of_links
    global g_number_of_agents

    #step 1: read input_node
    with open('../../test/node.csv', 'r', encoding='utf-8') as fp:
        reader = csv.DictReader(fp)
        node_seq_no = 0
        for line in reader:
            node = Node(node_seq_no, line['node_id'], line['zone_id'])
            node_list.append(node)
            internal_node_seq_no_dict[node.external_node_id] = node_seq_no
            external_node_id_dict[node_seq_no] = node.external_node_id
            if node.zone_id not in zone_to_nodes_dict.keys():
                zone_to_nodes_dict[int(node.zone_id)] = list()
                zone_to_nodes_dict[int(node.zone_id)].append(node.external_node_id)
            else:
                zone_to_nodes_dict[int(node.zone_id)].append(node.external_node_id)
            node_seq_no += 1
        g_number_of_nodes = node_seq_no
    print('the number of nodes is', g_number_of_nodes)

    #step 2: read input_link
    with open('../../test/link.csv', 'r', encoding='utf-8') as fp:
        reader = csv.DictReader(fp)
        link_seq_no = 0
        for line in reader:
            from_node_no = internal_node_seq_no_dict[int(line['from_node_id'])]
            to_node_no = internal_node_seq_no_dict[int(line['to_node_id'])]
            link = Link(link_seq_no, 
                        from_node_no, 
                        to_node_no,
                        int(line['from_node_id']),
                        int(line['to_node_id']),
                        line['length'],
                        line['lanes'],
                        line['free_speed'],
                        line['capacity'],
                        line['link_type'],
                        line['VDF_alpha1'],
                        line['VDF_beta1'])
            node_list[link.from_node_seq_no].outgoing_link_list.append(link)
            node_list[link.to_node_seq_no].incoming_link_list.append(link)
            link_list.append(link)
            link_seq_no += 1
        g_number_of_links = link_seq_no
    print('the number of links is', g_number_of_links)

    #step 3:read input_agent
    with open('../../test/demand.csv', 'r', encoding='utf-8') as fp:
        reader = csv.DictReader(fp)
        agent_id = 1
        agent_type = 'v'
        for line in reader:
            volume = line['volume']
            volume_agent_size = int(float(volume) + 1)
    
            # only test up to 10k
            if agent_id >= 10000 :
                break 
    
            for i in range(volume_agent_size):
                agent = Agent(agent_id,agent_type,line['o_zone_id'], line['d_zone_id'])

                # step 3.1 generate o_node_id and d_node_id randomly according 
                # to o_zone_id and d_zone_id 
                if zone_to_nodes_dict.get(agent.o_zone_id, -1) == -1 : 
                     continue
                if zone_to_nodes_dict.get(agent.d_zone_id, -1) == -1 : 
                     continue 

                agent_id = agent_id + 1
                agent.o_node_id = choice(zone_to_nodes_dict[agent.o_zone_id])
                agent.d_node_id = choice(zone_to_nodes_dict[agent.d_zone_id])
                
                # step 3.2 the initialization of the agent
                # agent.Initialization()

                # step 3.3: update the g_simulation_start_time_in_min and 
                # g_simulation_end_time_in_min 
                if agent.departure_time_in_min < g_simulation_start_time_in_min:
                    g_simulation_start_time_in_min = agent.departure_time_in_min
                if agent.departure_time_in_min > g_simulation_end_time_in_min:
                    g_simulation_end_time_in_min = agent.departure_time_in_min

                #step 3.4: add the agent to the time dependent agent list     
                if agent.departure_time_in_simu_interval not in agent_td_list_dict.keys():
                    agent_td_list_dict[agent.departure_time_in_simu_interval] = list()
                    agent_td_list_dict[agent.departure_time_in_simu_interval].append(agent.agent_seq_no)
                else:
                    agent_td_list_dict[agent.departure_time_in_simu_interval].append(agent.agent_seq_no)
                agent_list.append(agent)

    #step 3.6:sort agents by the departure time
    sort_fun = operator.attrgetter("departure_time_in_min")
    agent_list.sort(key=sort_fun)
    g_number_of_agents = len(agent_list)
    
    print('the number of agents is', g_number_of_agents)

    for agent_no in range(g_number_of_agents):
        agent_list[agent_no].agent_seq_no = agent_no
    
        
    #step 3.7:start simulation interval and end simulation interval
    g_start_simu_interval_no = int(g_simulation_start_time_in_min * 60 / NUMBER_OF_SECONDS_PER_SIMU_INTERVAL)
    g_end_simu_interval_no = g_start_simu_interval_no + LENGTH_OF_SIMULATION_TIME_HORIZON_IN_INTERVAL

        
def g_OutputFiles(link_list, agent_list, external_node_id_dict):
    
    if g_modeling_method == 2:
        with open ('../../test/link_performance.csv', 'w', newline='') as fp:
            writer = csv.writer(fp)
            # write header
            line = ["link_id", "from_node_id", "to_node_id", "time_period",
                    "volume", "CA", "CD", "density", "queue", "travel_time",
                    "waiting_time_in_sec", "speed"]
            writer.writerow(line)
            
            # write each line
            for link in link_list:
                for relative_t in range(LENGTH_OF_SIMULATION_TIME_HORIZON_IN_INTERVAL):
                    abs_t = g_R2A_simu_interval(relative_t)
                    # if(relative_t%(60/g_number_of_seconds_per_simu_interval)==0):
                    time_in_min =  int(relative_t / (60/NUMBER_OF_SECONDS_PER_SIMU_INTERVAL))
                    abs_t_in_min = int(abs_t / (60/NUMBER_OF_SECONDS_PER_SIMU_INTERVAL))
                
                    volume = 0
                    queue = 0
                    waiting_time_in_sec = 0
                    arrival_rate = 0
                    avg_waiting_time_in_sec = 0
                    avg_travel_time_in_min = float(link.general_travel_time_in_min + avg_waiting_time_in_sec/60)
                    speed = link.length / (max(0.00001, avg_travel_time_in_min/60))
                    virtual_arrival = 0

                    if time_in_min >= 1:
                        volume = link.td_link_cumulative_departure[relative_t] - link.td_link_cumulative_departure[relative_t - int(60/NUMBER_OF_SECONDS_PER_SIMU_INTERVAL)]
    
                        if relative_t - link.general_travel_time_in_simu_interval >= 0:  
                            virtual_arrival = link.td_link_cumulative_arrival[int(relative_t - link.general_travel_time_in_simu_interval)]
                                                                                
                        queue = virtual_arrival - link.td_link_cumulative_departure[relative_t]        
                        # waiting_time_count = 0
                        waiting_time_in_sec = link.td_link_waiting_time[time_in_min] * NUMBER_OF_SECONDS_PER_SIMU_INTERVAL 
                        if (relative_t + int(60 / NUMBER_OF_SECONDS_PER_SIMU_INTERVAL) < LENGTH_OF_SIMULATION_TIME_HORIZON_IN_INTERVAL):
                            arrival_rate = link.td_link_cumulative_arrival[relative_t + int(60 / NUMBER_OF_SECONDS_PER_SIMU_INTERVAL)] - link.td_link_cumulative_arrival[relative_t]
                        avg_waiting_time_in_sec = waiting_time_in_sec / max(1, arrival_rate)
                        avg_travel_time_in_min = link.general_travel_time_in_min + avg_waiting_time_in_sec/60
                        speed = link.length / (max(0.00001, avg_travel_time_in_min) / 60)
                    # end of if ime_in_min >= 1:

                    density = (link.td_link_cumulative_arrival[relative_t] - link.td_link_cumulative_departure[relative_t]) / (link.length * link.lanes)
    
                    line = [
                        link.link_seq_no+1, 
                        link.external_from_node, 
                        link.external_to_node,
                        time_stamp_to_HHMMSS(abs_t_in_min)+"_"+time_stamp_to_HHMMSS(abs_t_in_min+1),
                        volume,
                        link.td_link_cumulative_arrival[relative_t],
                        link.td_link_cumulative_departure[relative_t],
                        density,
                        queue,
                        avg_travel_time_in_min,
                        waiting_time_in_sec,
                        speed
                    ]
                # end of for relative_t in ...
                writer.writerow(line)
            # end of link in link_list:
        
    with open("../../test/agent.csv", "w", newline='') as fp:
        writer = csv.writer(fp)
        # write header
        line = ["agent_id", "agent_type", "o_node_id",
                "d_node_id", "o_zone_id", "d_zone_id",
                "travel_time", "node_sequence", "time_sequence"]
        writer.writerow(line)
        # abs_satrt_time_in_min=int(g_start_simu_interval_no/g_number_of_simu_intervals_per_min )
        
        # write each line
        for agent in agent_list:
            if not agent.feasible_path_exist_flag:
                continue

            # define the default values
            cost = 0 
            path_node_str = ';'.join([str(external_node_id_dict[i]) for i in agent.path_node_seq_no_list])
            path_time_sequence = ''
            
            # if simulation has been performed
            if g_modeling_method == 2:
                path_time_list = list()
                cost = (agent.veh_link_departure_time_in_simu_interval[-1]-agent.veh_link_arrival_time_in_simu_interval[0])\
                        * NUMBER_OF_SECONDS_PER_SIMU_INTERVAL / 60
                for i in range(len(agent.path_link_seq_no_list)):
                    TA_in_min = agent.veh_link_arrival_time_in_simu_interval[i] * NUMBER_OF_SECONDS_PER_SIMU_INTERVAL / 60
                    TD_in_min = agent.veh_link_departure_time_in_simu_interval[i] * NUMBER_OF_SECONDS_PER_SIMU_INTERVAL / 60
                    if i == 0:
                        path_time_list.extend([time_stamp_to_HHMMSS(TA_in_min), 
                                               time_stamp_to_HHMMSS(TD_in_min)])
                    else:
                        path_time_list.append(time_stamp_to_HHMMSS(TD_in_min))
                    path_time_sequence = ';'.join(path_time_list)
                    
            line = [agent.agent_id,
                    agent.agent_type,
                    agent.o_node_id,
                    agent.d_node_id,
                    agent.o_zone_id,
                    agent.d_zone_id,cost,
                    path_node_str,
                    path_time_sequence]
            writer.writerow(line)

    with open("../../test/solution.csv", 'w', newline='') as fp:
        writer = csv.writer(fp)
        line = ["number_of_nodes", 
                "number_of_links", 
                "number_of_agents",
                "CPU running time"]
        writer.writerow(line)
        line_value = [g_number_of_nodes, g_number_of_links, 
                      g_number_of_agents, end_time-begin_time]
        writer.writerow(line_value)

         
if __name__=="__main__":

    network = Network()
    g_ReadInputData(network.node_list, 
                    network.link_list, 
                    network.agent_list,
                    network.internal_node_seq_no_dict,
                    network.external_node_id_dict,
                    network.agent_td_list_dict,
                    network.zone_to_nodes_dict)

    network.allocate()
    
    begin_time = time.time()

    g_TrafficAssignment(network)
    #g_TrafficSimulation(network.node_list, network.link_list, 
    #                    network.agent_list, network.agent_td_list_dict)

    end_time = time.time()
    g_OutputFiles(network.link_list, 
                  network.agent_list, 
                  network.external_node_id_dict)
        
    # print('end time: {0: .2f}'.format(end_time))
    print('shortest path processing time: {0: .2f}'.format(end_time-begin_time) 
          + 's')