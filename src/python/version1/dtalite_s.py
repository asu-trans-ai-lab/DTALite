import csv
from time import *
import operator

#data block define
_max_number_of_nodes=10000   
_max_label_cost_in_shortest_path = 10000

#data block simulation time horizon 
g_length_of_simulation_time_horizon_in_min=90    #the length of simulation time(min)
g_number_of_seconds_per_simu_interval = 6;  #the number of seconds per simulation interval
g_length_of_simulation_time_horizon_in_interval =int(g_length_of_simulation_time_horizon_in_min*60/g_number_of_seconds_per_simu_interval)
g_number_of_simu_intervals_per_min=int(60/g_number_of_seconds_per_simu_interval)  

#data block start and end time of horizon 
g_simulation_start_time_in_min = 9999  # start time of the simulation: 9999 as default value to be updated in reading function for agent file
g_simulation_end_time_in_min = 0 # end time of the simulation:0 as default value to be updated in reading function for agent file
g_start_simu_interval_no=0  # the number of start simualation interval 
g_end_simu_interval_no=0 # the number of end simualation interval 

#data block input data statistics 
g_number_of_nodes=0
g_number_of_links=0
g_number_of_agents=0

#data block simulation result statistics
g_cumulative_arrival_count=0
g_cumulative_departure_count = 0

g_number_of_assignment_iterations=1   #the number of assignment iterations: 1 as default and we can set it's value 

#list definition 
g_node_list=list()
g_link_list=list()
g_agent_list=list()

#dictionary definition 
g_internal_node_seq_no_dict=dict() #key: external node id   value:internal node id
g_external_node_id_dict=dict()  #key: internal node id   value:external node id
g_agent_td_list_dict=dict()     #td:time-dependent    key:simulation time interval   value:agents(list) need to be activated

# data conversion utility functions
def time_stamp_to_HHMMSS(time_in_minutes):
    hours=int(time_in_minutes/60)
    minutes=int(time_in_minutes-hours*60)
    seconds=int((time_in_minutes-hours*60-minutes)*60)
    return time_int_to_str(hours)+time_int_to_str(minutes)+":"+time_int_to_str(seconds)

def time_int_to_str(time_int):
    if time_int<10:
        return "0"+str(time_int)
    else:
        return str(time_int)
    
# convert absolute simulation interval to relative simulation interval
# A2R:absolute 24_hour_world clock time in simulation interval to relative time in simulation interval from 0
def g_A2R_simu_interval(abslute_time_in_simu_intetrval_no):
    return (abslute_time_in_simu_intetrval_no-g_start_simu_interval_no)
# convert relative simulation interval to absolute simulation interval
# R2A:relative time in simulation interval from 0 to absolute 24_hour_world clock time in simulation interval
def g_R2A_simu_interval(relative_time_in_simu_intetrval_no):
    return (relative_time_in_simu_intetrval_no+g_start_simu_interval_no)

# comments: external_node_id: the id of node
#           node_seq_no: the index of the node and we call the node by its index
#           we use g_internal_node_seq_no_dict(id to index) and g_external_node_id_dict(index to id) to map them
class Node:
    def __init__(self,external_node_id):  # the attribute of node 
        self.node_seq_no=0
        self.external_node_id=external_node_id
        self.outgoing_link_list=list()
        self.incoming_link_list=list()
        self.Initialization()

    def Initialization(self):  # establish mapping relation between external node id and internal node id
        global g_number_of_nodes
        self.node_seq_no=g_number_of_nodes
        g_internal_node_seq_no_dict[int(self.external_node_id)]=self.node_seq_no
        g_external_node_id_dict[int(self.node_seq_no)]=self.external_node_id
        g_number_of_nodes+=1
        

class Link:
    def __init__(self,from_node_id,to_node_id,length,lanes,free_speed,capacity,   # the attribute of link
                 link_type,VDF_alpha1,VDF_beta1):   
        self.from_node_seq_no=g_internal_node_seq_no_dict[int(from_node_id)]
        self.to_node_seq_no=g_internal_node_seq_no_dict[int(to_node_id)]
        self.external_from_node=int(from_node_id)
        self.external_to_node=int(to_node_id)
        self.type=int(link_type) # 1:one direction 2:two way
        self.lanes=int(lanes)
        self.BPR_alpha=float(VDF_alpha1)
        self.BPR_beta=float(VDF_beta1)
        self.flow_volume=0
        self.link_capacity = float(capacity)* int(lanes) # capacity is lane capacity per hour
        self.length=float(length) # length is mile or km 
        self.free_flow_travel_time_in_min=self.length / max(0.001,int(free_speed)) * 60  # length:km   free_speed: km/h
        self.cost=self.free_flow_travel_time_in_min
        #the capacity of discharging for each link,  td is time-dependent, to used in traffic signal control or work zone scheduling applications 
        self.td_link_outflow_capacity=[int(self.link_capacity/(60*g_number_of_simu_intervals_per_min))]*(g_length_of_simulation_time_horizon_in_interval+1) # self.link_capacity/g_number_of_simulation_time) should be per interval capacity 

        #count the cumulative arrival and departure of links each simulation interval, td is time-dependent
        self.td_link_cumulative_arrival=[0]*(g_length_of_simulation_time_horizon_in_interval+1) # [0] is setting up zero as the default 
        self.td_link_cumulative_departure=[0]*(g_length_of_simulation_time_horizon_in_interval+1)  # [0] is setting up zero as the default 

        #the time-dependent waiting time, td is time-dependent
        self.td_link_waiting_time=[0]*(g_length_of_simulation_time_horizon_in_min+1) # [0] is setting up zero as the default 
        
        self.entrance_queue=list() # link-in queue  of each link
        self.exit_queue=list() # link-out queue  of each link
        

        self.Initialization()

    def Initialization(self): # update the number of links 
        global g_number_of_links
        self.link_seq_no=g_number_of_links
        g_number_of_links +=1
    def ResetMOE(self): # reset the measures of effectiveness
        self.cumulative_arrival_count = 0
        self.cumulative_departure_count = 0
        self.cumulative_virtual_delay_count = 0
    def CalculateBPRFunction(self): # update the cost of link
        self.general_travel_time_in_min = self.free_flow_travel_time_in_min*(1 + self.BPR_alpha*((self.flow_volume / max(0.00001, self.link_capacity))**self.BPR_beta))  #travel time in min
        self.general_travel_time_in_simu_interval=self.general_travel_time_in_min*g_number_of_simu_intervals_per_min  #travel time in simulation interval
        self.cost = self.general_travel_time_in_min
        
# comments: agent_id: the id of agent
#           agent_seq_no: the index of the agent and we call the agent by its index
class Agent():
    def __init__(self,agent_id,agent_type,o_node_id,    # the attribute of agent
                 d_node_id,departure_time_in_min,PCE,
                 ):
        self.agent_id=agent_id
        self.agent_type=agent_type  #vehicle 
        self.o_node_id=int(o_node_id) 
        self.d_node_id=int(d_node_id) 
        self.path_node_seq_no_list=list()
        self.path_link_seq_no_list=list()
        self.current_link_seq_no_in_path=0 
        self.departure_time_in_min=float(departure_time_in_min)
        self.PCE_factor=int(PCE)  #Passenger Car Equivalent (PCE) of the agent
        self.path_cost=0
        self.b_generated=False
        self.b_complete_trip=False

        self.departure_time_in_simu_interval = int(self.departure_time_in_min*60.0 / g_number_of_seconds_per_simu_interval + 0.5)
        
    def Initialization(self): # update the number of agents 
        global g_number_of_agents
        self.agent_seq_no=g_number_of_agents
        g_number_of_agents+=1
        
    def Initialize_for_simulation(self): 
        if(len(self.path_node_seq_no_list)!=0):      # if the agent can find its path 
            self.veh_link_arrival_time_in_simu_interval=[-1]*(len(self.path_node_seq_no_list)-1)  # the arrival time on each link in the agent's path
            self.veh_link_departure_time_in_simu_interval=[-1]*(len(self.path_node_seq_no_list)-1) # the depature time on each link in the agent's path
            self.veh_link_arrival_time_in_simu_interval[0]=self.departure_time_in_simu_interval 
            self.feasible_path_exist_flag=True # the path exist
        else:
            self.feasible_path_exist_flag=False


def g_read_input_data():
    global g_node_list
    global g_link_list
    global g_agent_list
    global g_simulation_start_time_in_min
    global g_simulation_end_time_in_min
    global g_start_simu_interval_no
    global g_end_simu_interval_no

    # input: node.csv, link.csv and demand.csv
    # output: craet the network and agent 

    #step 1: read input_node
    with open('node.csv','r',encoding='utf-8') as fp:
        reader=csv.DictReader(fp)
        for line in reader:
            node=Node(line['node_id'])
            g_node_list.append(node)
    print('the number of nodes is',g_number_of_nodes)

    #step 2: read input_link
    with open('link.csv','r',encoding='utf-8') as fp:
        reader=csv.DictReader(fp)
        for line in reader:
            link=Link(line['from_node_id'],line['to_node_id'],
                      line['length'],line['lanes'],line['free_speed'],
                      line['capacity'],line['link_type'],line['VDF_alpha1'],
                      line['VDF_beta1'])
            g_node_list[link.from_node_seq_no].outgoing_link_list.append(link)
            g_node_list[link.to_node_seq_no].incoming_link_list.append(link)

            g_link_list.append(link)
    print('the number of links is',g_number_of_links)

    #step 3:read input_agent
    with open('demand.csv','r',encoding='utf-8') as fp:
        reader=csv.DictReader(fp)
        for line in reader:
            agent=Agent(line['agent_id'],line["agent_type"],line['o_node_id'],
                        line['d_node_id'],line['departure_time_in_min'],
                        line['PCE']
                       )
            # step 3.1 check if the o_node_id or d_node_id is not defined in nodes
            if (agent.o_node_id not in g_internal_node_seq_no_dict.keys() or
                agent.d_node_id not in g_internal_node_seq_no_dict.keys()):
                continue
            # step 3.2 check if the o_node_id equals d_node_id 
            if(agent.o_node_id==agent.d_node_id ):
                continue
            # step 3.3 the initialization of the agent
            agent.Initialization()

            # step 3.4: update the g_simulation_start_time_in_min and g_simulation_end_time_in_min 
            if agent.departure_time_in_min<g_simulation_start_time_in_min:
                g_simulation_start_time_in_min=agent.departure_time_in_min
            if agent.departure_time_in_min>g_simulation_end_time_in_min:
                g_simulation_end_time_in_min=agent.departure_time_in_min

            #step 3.5: add the agent to the time dependent agent list     
            if(agent.departure_time_in_simu_interval not in g_agent_td_list_dict.keys()):
                g_agent_td_list_dict[agent.departure_time_in_simu_interval]=list()
                g_agent_td_list_dict[agent.departure_time_in_simu_interval].append(agent.agent_seq_no)
            else:
                g_agent_td_list_dict[agent.departure_time_in_simu_interval].append(agent.agent_seq_no)
            g_agent_list.append(agent)
            
    print('the number of agents is',g_number_of_agents)

    #step 3.6:sort agents by the departure time
    sort_fun = operator.attrgetter("departure_time_in_min")
    g_agent_list.sort(key=sort_fun)

    for agent_no in range(len(g_agent_list)):
        g_agent_list[agent_no].agent_seq_no=agent_no
        
    #step 3.7:start simulation interval and end simulation interval
    g_start_simu_interval_no =int(g_simulation_start_time_in_min * 60 / g_number_of_seconds_per_simu_interval)
    g_end_simu_interval_no =g_start_simu_interval_no + g_length_of_simulation_time_horizon_in_interval
    

class Network:
    
    def allocate(self,number_of_nodes,number_of_links):
        self.node_predecessor = [-1]*number_of_nodes
        self.node_label_cost = [_max_label_cost_in_shortest_path]*number_of_nodes
        self.link_predecessor = [-1]*number_of_nodes
        self.link_cost_array = [1.0]*number_of_links
        self.link_volume_array = [0.0]*number_of_links  # initialization for traffic assignment 

        for j in range(number_of_links):
            self.link_volume_array[j] = 0.0
            self.link_cost_array[j] = g_link_list[j].cost
        
    def optimal_label_correcting(self,origin_node,destination_node,departure_time):
    # input : origin_node,destination_node,departure_time
    # output : the shortest path
        origin_node=g_internal_node_seq_no_dict[origin_node]
        destination_node=g_internal_node_seq_no_dict[destination_node]
        if len(g_node_list[origin_node].outgoing_link_list) == 0:
            return 0
        
        for i in range(g_number_of_nodes): #Initialization for all nodes
            self.node_label_cost[i] = _max_label_cost_in_shortest_path
            self.node_predecessor[i] = -1 #pointer to previous node index from the current label at current node 
            self.link_predecessor[i] = -1 #pointer to previous node index from the current label at current node
        
        self.node_label_cost[origin_node] = departure_time
        SEList = []  # scan eligible list
        
        SEList.append(origin_node)

        while len(SEList)>0:
            from_node = SEList[0]
            del SEList[0]
            for k in range(len(g_node_list[from_node].outgoing_link_list)):
                to_node = g_node_list[from_node].outgoing_link_list[k].to_node_seq_no
             
                new_to_node_cost = self.node_label_cost[from_node] + self.link_cost_array[g_node_list[from_node].outgoing_link_list[k].link_seq_no]
                
                if (new_to_node_cost < self.node_label_cost[to_node]):  #we only compare cost at the downstream node ToID at the new arrival time t
                    # update cost label and node/time predecessor
                    self.node_label_cost[to_node] = new_to_node_cost
                    self.node_predecessor[to_node] = from_node  #pointer to previous physical node index from the current label at current node and time
                    self.link_predecessor[to_node] = g_node_list[from_node].outgoing_link_list[k].link_seq_no  #pointer to previous physical node index from the current label at current node and time
                    SEList.append(to_node)
                                        

        if (destination_node >= 0 and self.node_label_cost[destination_node] < _max_label_cost_in_shortest_path):
            return 1
        else: 
            return -1



    def find_path_for_agents(self,iteration_no):
    # input : the number of iterations
    # output : updated link volume
            
        for s in range(g_number_of_links):
            self.link_volume_array[s] = 0

        #step 1: find shortest path if needed 
        # MSA: only 1/k agents will swtich to the shortest path; (k-1)/K will not have changes to select shortest path, stay on the current path
        for i in range(len(g_agent_list)):
            residual = i % (iteration_no + 1)
            if (residual != 0):  #no need to compute a new path at this iteration
                continue         #that is, it will reuse the path from the previous iteration, stored at p_agent->path_link_seq_no_list.
                
            #else move to the next line for finding the shortest path 
            g_agent_list[i].path_link_seq_no_list = []
            g_agent_list[i].path_node_seq_no_list = []
            #step 2 buil SP tree
            
            return_value = self.optimal_label_correcting(g_agent_list[i].o_node_id, g_agent_list[i].d_node_id, g_agent_list[i].departure_time_in_min)
            #step 3 update path

            if (return_value == -1):
                #print('agent ',i,'cannot find destination node')
                continue

            current_node_seq_no = g_internal_node_seq_no_dict[g_agent_list[i].d_node_id]
            g_agent_list[i].path_cost = self.node_label_cost[current_node_seq_no]
            
            while (current_node_seq_no>=0):
                if (current_node_seq_no >= 0):  
                    current_link_seq_no = self.link_predecessor[current_node_seq_no]
                
                    if(current_link_seq_no>=0):
                        g_agent_list[i].path_link_seq_no_list.append(current_link_seq_no)      

                
                    g_agent_list[i].path_node_seq_no_list.append(current_node_seq_no)

                current_node_seq_no = self.node_predecessor[current_node_seq_no]
            g_agent_list[i].path_node_seq_no_list.reverse()
            g_agent_list[i].path_link_seq_no_list.reverse()
         
        #step 2:  scan the shortest path to compute the link volume, 
        
        for i in range(len(g_agent_list)):

            for j in range(len(g_agent_list[i].path_link_seq_no_list)):#for each link in the path of this agent
                link_seq_no = g_agent_list[i].path_link_seq_no_list[j]
                self.link_volume_array[link_seq_no] += g_agent_list[i].PCE_factor 

def g_traffic_assignment():
    global g_link_list
    
    network = Network()
    network.allocate(g_number_of_nodes,g_number_of_links)

    print('Finding shortest path for all agents......')

    for i in range(g_number_of_assignment_iterations): # perform traffic assignment to reach a UE or SO criterion with given number of iterations 
        print('iteration_no',i,'......')
        for l in range(g_number_of_links):
            g_link_list[l].CalculateBPRFunction()  # use the link volume from previous iteration to calculate the link travel time for each link 
            network.link_cost_array[l] = g_link_list[l].cost # use link travel time as cost for shortest path calculation 

            network.find_path_for_agents(i) # find the paths for all agents, we do have a MSA rule being applied to select which agents will update the shortest path     
                        
            for k in range(g_number_of_links): # tally link volume 
                g_link_list[k].flow_volume = network.link_volume_array[k]
        
        
        
def g_traffic_simulation():
    global g_cumulative_arrival_count
    global g_cumulative_departure_count
     
    #step 1 initialization 
    #initialization for each link
    for li in range(len(g_link_list)):
        g_link_list[li].ResetMOE()
    #initialization for each agent
    for agent_no in range(len(g_agent_list)):
        g_agent_list[agent_no].Initialize_for_simulation()

    g_active_agent_list=list()
    current_active_agent_id = 0
    
    #step 2 time-based loop for simulation 
    for t in range(g_start_simu_interval_no,g_end_simu_interval_no,1):
        number_of_simu_interval_per_min = 60 / g_number_of_seconds_per_simu_interval;
        if(t%number_of_simu_interval_per_min==0):
            print("simu time= ",int(t / number_of_simu_interval_per_min),"min, ", "CA=",g_cumulative_arrival_count,", CD=",g_cumulative_departure_count)

        relative_t = g_A2R_simu_interval(t) # obtain relative time index from global absolute time stamp
        for li in range(0,len(g_link_list)):# step 2.1. carry the cumulative count statistics from the previous time step
            link=g_link_list[li]
            if (relative_t >= 1):
                g_link_list[link.link_seq_no].td_link_cumulative_departure[relative_t]=g_link_list[link.link_seq_no].td_link_cumulative_departure[relative_t-1]
                g_link_list[link.link_seq_no].td_link_cumulative_arrival[relative_t] =g_link_list[link.link_seq_no].td_link_cumulative_arrival[relative_t-1]
  
        #step 2.2. transfer active agents to the first link of their paths
        if(t in g_agent_td_list_dict.keys()):  #activate the agent td: time dependent
            for j in range(0,len(g_agent_td_list_dict[t])):
                agent_no=g_agent_td_list_dict[t][j]
                agent=g_agent_list[agent_no]
                if(agent.feasible_path_exist_flag==True):  # we only consider the agents with feasible paths
                    agent.b_generated=True
                    first_link_seq=agent.path_link_seq_no_list[0]
                    g_link_list[first_link_seq].td_link_cumulative_arrival[relative_t]+=1
                    g_link_list[first_link_seq].entrance_queue.append(agent.agent_seq_no) # put active agent to the entrance queue of the first link 
                    g_cumulative_arrival_count+=1
        #step 2.3 link transfer process of moving vehicles from entrance queue to the exit queue of a link         
        for li in range(0,len(g_link_list)):
            link=g_link_list[li]
            while(len(link.entrance_queue)>0):  #there are agents in the entrance_queue
                agent_seq=link.entrance_queue[0]
                link.entrance_queue.pop(0)
                link.exit_queue.append(agent_seq)
                g_agent_list[agent_seq].veh_link_departure_time_in_simu_interval[g_agent_list[agent_seq].current_link_seq_no_in_path]=\
                    g_agent_list[agent_seq].veh_link_arrival_time_in_simu_interval[g_agent_list[agent_seq].current_link_seq_no_in_path]+link.general_travel_time_in_simu_interval
        #step 2.4 node transfer process of moving vehicles from exit queue to the entrance queue of a link  
        for no in range(0,len(g_node_list)):
            node=g_node_list[no]
            for l in range(0,len(node.incoming_link_list)):
                incoming_link_index=(l+t)%(len(node.incoming_link_list)) 
                #we will start with different first link from the incoming link list, 
                #equal change, regardless of # of lanes or main line vs. ramp, but one can use service arc, to control the effective capacity rates, e.g. through a metered ramp, to 
		        #allow mainline to use the remaining flow

                link_seq_no=node.incoming_link_list[incoming_link_index].link_seq_no

                #check if the current link has sufficient capacity
                #check if there are agents in exit queue
                #step 2.4.1 check condition of outflow capacity
                while (g_link_list[link_seq_no].td_link_outflow_capacity[relative_t]>=1 and len(g_link_list[link_seq_no].exit_queue)>=1):
                    agent_no=node.incoming_link_list[incoming_link_index].exit_queue[0]
                   
                    if(g_agent_list[agent_no].veh_link_departure_time_in_simu_interval[g_agent_list[agent_no].current_link_seq_no_in_path]>t):
                        break        #the future departure time on this link is later than the current time
               
                    #step 2.4.2 end of path checking 
                    if(g_agent_list[agent_no].current_link_seq_no_in_path==len(g_agent_list[agent_no].path_link_seq_no_list)-1): # end of the path
                        node.incoming_link_list[incoming_link_index].exit_queue.pop(0)
                        g_agent_list[agent_no].b_complete_trip = True
                        g_link_list[link_seq_no].td_link_cumulative_departure[relative_t] += 1
                        g_cumulative_departure_count += 1
                    #step 2.4.3 move from exit queue to the entrace queue
                    else:  #  not complete the trip. move to the next link's entrance queue
                        next_link_seq=g_agent_list[agent_no].path_link_seq_no_list[g_agent_list[agent_no].current_link_seq_no_in_path+1]
                        node.incoming_link_list[incoming_link_index].exit_queue.pop(0)
                        g_link_list[next_link_seq].entrance_queue.append(agent_no)
                        g_agent_list[agent_no].veh_link_departure_time_in_simu_interval[g_agent_list[agent_no].current_link_seq_no_in_path]=t
                        g_agent_list[agent_no].veh_link_arrival_time_in_simu_interval[g_agent_list[agent_no].current_link_seq_no_in_path+1]=t
                        
                        actual_travel_time = t - g_agent_list[agent_no].veh_link_arrival_time_in_simu_interval[g_agent_list[agent_no].current_link_seq_no_in_path]
			#step 2.4.4 handling waiting vehicles in the queue 
                        waiting_time = actual_travel_time - g_link_list[link_seq_no].general_travel_time_in_min
                        temp_relative_time=g_A2R_simu_interval(g_agent_list[agent_no].veh_link_arrival_time_in_simu_interval[g_agent_list[agent_no].current_link_seq_no_in_path])
                        time_in_min=int(temp_relative_time/g_number_of_simu_intervals_per_min)
                  
                        g_link_list[link_seq_no].td_link_waiting_time[time_in_min] += waiting_time
                        g_link_list[link_seq_no].td_link_cumulative_departure[relative_t] += 1
                        g_link_list[next_link_seq].td_link_cumulative_arrival[relative_t] += 1
                    
                    g_agent_list[agent_no].current_link_seq_no_in_path+=1
                    #step 2.4.5 update time dependent link out flow capacity
                    g_link_list[link_seq_no].td_link_outflow_capacity[relative_t]-=1 
        
        
def g_output_files():
    with open ('link_performance.csv','w',newline='') as fp:
        writer=csv.writer(fp)
        line=["link_id","from_node_id","to_node_id","time_period","volume","CA","CD","density","queue","travel_time","waiting_time_in_sec","speed"]
        writer.writerow(line)
        
        for link in g_link_list:
            for relative_t in range(0,g_length_of_simulation_time_horizon_in_interval,1):
                abs_t = g_R2A_simu_interval(relative_t)
                
##                if(relative_t%(60/g_number_of_seconds_per_simu_interval)==0):
                time_in_min =  int(relative_t / (60 / g_number_of_seconds_per_simu_interval))
                abs_t_in_min= int(abs_t / (60 / g_number_of_seconds_per_simu_interval))
            
                volume = 0
                queue = 0
                waiting_time_in_sec = 0
                arrival_rate = 0
                avg_waiting_time_in_sec = 0
                avg_travel_time_in_min =float(link.general_travel_time_in_min + avg_waiting_time_in_sec/60.0)
                speed = link.length /(max(0.00001,avg_travel_time_in_min/60.0))
                virtual_arrival = 0
                if (time_in_min >= 1):
                
                    volume = link.td_link_cumulative_departure[relative_t] - link.td_link_cumulative_departure[relative_t - int(60 / g_number_of_seconds_per_simu_interval)]

                    if (relative_t- link.general_travel_time_in_simu_interval>= 0):  
                        virtual_arrival = link.td_link_cumulative_arrival[int(relative_t - link.general_travel_time_in_simu_interval)]
                                                                            
                    queue = virtual_arrival - link.td_link_cumulative_departure[relative_t]
                                                    

                    waiting_time_count = 0

                    waiting_time_in_sec = link.td_link_waiting_time[time_in_min] * g_number_of_seconds_per_simu_interval 
                    if(relative_t + int(60 / g_number_of_seconds_per_simu_interval)<g_length_of_simulation_time_horizon_in_interval):
                        arrival_rate = link.td_link_cumulative_arrival[relative_t + int(60 / g_number_of_seconds_per_simu_interval)] - link.td_link_cumulative_arrival[relative_t]
                    avg_waiting_time_in_sec = waiting_time_in_sec / max(1, arrival_rate)

                    avg_travel_time_in_min = link.general_travel_time_in_min+ avg_waiting_time_in_sec/60.0
                    speed = link.length / (max(0.00001,avg_travel_time_in_min) / 60)
                    
                density = (link.td_link_cumulative_arrival[relative_t] - link.td_link_cumulative_departure[relative_t]) / (link.length * link.lanes);

                line=[link.link_seq_no+1,link.external_from_node,link.external_to_node,time_stamp_to_HHMMSS(abs_t_in_min)+"_"+time_stamp_to_HHMMSS(abs_t_in_min+1),
                      volume,link.td_link_cumulative_arrival[relative_t],link.td_link_cumulative_departure[relative_t],
                      density,queue,avg_travel_time_in_min,waiting_time_in_sec,speed]
                
                writer.writerow(line)
    with open("agent.csv","w",newline='') as fp:
        writer=csv.writer(fp)
        line=["agent_id","agent_type","o_node_id","d_node_id","o_zone_id","d_zone_id","travel_time",
              "node_sequence","time_sequence"]
        writer.writerow(line)
        abs_satrt_time_in_min=int(g_start_simu_interval_no/g_number_of_simu_intervals_per_min )
        
        for agent in g_agent_list:
            if(agent.feasible_path_exist_flag==False):
                continue
            path_time_list=list()
            cost=(agent.veh_link_departure_time_in_simu_interval[-1]-agent.veh_link_arrival_time_in_simu_interval[0])\
                * g_number_of_seconds_per_simu_interval / 60.0
            for i in range(0,len(agent.path_link_seq_no_list),1):
                TA_in_min=agent.veh_link_arrival_time_in_simu_interval[i]* g_number_of_seconds_per_simu_interval / 60.0
                TD_in_min=agent.veh_link_departure_time_in_simu_interval[i]* g_number_of_seconds_per_simu_interval / 60.0
                if(i==0):    
                    path_time_list.extend([time_stamp_to_HHMMSS(TA_in_min),time_stamp_to_HHMMSS(TD_in_min)])
                else:
                    path_time_list.append(time_stamp_to_HHMMSS(TD_in_min))
                path_time_sequence=';'.join(path_time_list)
            path_node_str=';'.join([str(g_external_node_id_dict[i]) for i in agent.path_node_seq_no_list])
            line=[agent.agent_id,agent.agent_type,agent.o_node_id,
                  agent.d_node_id,agent.o_node_id,
                  agent.d_node_id,cost,
                  path_node_str,path_time_sequence]
            writer.writerow(line)
    with open("solution.csv",'w',newline='') as fp:
        writer=csv.writer(fp)
        line=["number_of_nodes","number_of_links","number_of_agents","CPU running time"]
        writer.writerow(line)
        line_value=[g_number_of_nodes,g_number_of_links,g_number_of_agents,end_time-begin_time]
        writer.writerow(line_value)
            
if __name__=="__main__":

    begin_time=time()
    g_read_input_data()  
    g_traffic_assignment()
    g_traffic_simulation()

    end_time=time()
    g_output_files()
    

    
