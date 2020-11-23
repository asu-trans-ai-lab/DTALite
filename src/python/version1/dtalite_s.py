import csv
from time import *
import operator

max_number_of_nodes=10000   
_MAX_LABEL_COST = 10000

g_number_of_assignment_iterations=1   #the number of assignment iterations
g_number_of_simulation_time=90    #the length of simulation time(min)
g_number_of_seconds_per_interval = 6;  #the number of seconds per simulation interval
g_number_of_interval_per_min=int(60/g_number_of_seconds_per_interval)  
g_number_of_simulation_intervals=int(g_number_of_simulation_time*60/g_number_of_seconds_per_interval)

g_Simulation_StartTimeInMin = 9999  # start time of the simulation
g_Simulation_EndTimeInMin = 0 # end time of the simulation
g_start_simu_interval_no=0  # the number of start simualation interval 
g_end_simu_interval_no=0 # the number of end simualation interval 


g_number_of_nodes=0
g_number_of_links=0
g_number_of_agents=0
g_Cumulative_Arrival_Count=0
g_Cumulative_Departure_Count = 0

#list definition 
g_node_list=list()
g_link_list=list()
g_agent_list=list()

#dictionary definition 
g_internal_node_seq_no_dict=dict() #key: external node id   value:internal node id
g_external_node_id_dict=dict()  #key: internal node id   value:external node id
g_agent_TD_list_dict=dict()     #TD:time-dependent    key:simulation time interval   value:agents(list) need to be activated



class Node:
    def __init__(self,external_node_id):  # the attribute of node 
        self.node_seq_no=0
        self.external_node_id=external_node_id
        self.m_outgoing_link_list=list()
        self.m_incoming_link_list=list()
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
        self.link_capacity = float(capacity)* int(lanes)
        self.length=float(length)
        self.free_flow_travel_time_in_min=self.length / max(0.001,int(free_speed)) * 60  # length:km   free_speed: km/h

        #the capacity of discharging for each link
        self.m_LinkOutFlowCapacity=[int(self.link_capacity/g_number_of_simulation_time)]*(g_number_of_simulation_intervals+1)

        #count the cumulative arrival and departure of links each simulation interval
        self.m_LinkCumulativeArrival=[0]*(g_number_of_simulation_intervals+1)
        self.m_LinkCumulativeDeparture=[0]*(g_number_of_simulation_intervals+1)

        #the time-dependent waiting time 
        self.m_LinkTDWaitingTime=[0]*(g_number_of_simulation_time+1)
        
        self.Entrance_queue=list() # link-in queue  of each link
        self.Exit_queue=list() # link-out queue  of each link
        

        self.Initialization()

    def Initialization(self): # update the number of links 
        global g_number_of_links
        self.link_seq_no=g_number_of_links
        g_number_of_links +=1
    def ResetMOE(self): # reset the measures of effectiveness
        self.m_CumulativeArrivalCount = 0
        self.m_CumulativeDepartureCount = 0
        self.m_CumulativeVirtualDelayCount = 0
    def CalculateBRPFunction(self): # update the cost of link
        self.general_travel_time_in_min = self.free_flow_travel_time_in_min*(1 + self.BPR_alpha*((self.flow_volume / max(0.00001, self.link_capacity))**self.BPR_beta))  #travel time in min
        self.general_travel_time_in_simu_interval=self.general_travel_time_in_min*g_number_of_interval_per_min  #travel time in simulation interval
        self.cost = self.general_travel_time_in_min
        

class Agent():
    def __init__(self,agent_id,agent_type,o_node_id,    # the attribute of agent
                 d_node_id,departure_time_in_min,PCE,
                 ):
        self.agent_id=agent_id
        self.agent_type=agent_type  #vehicle 
        self.from_origin_node_id=int(o_node_id)
        self.to_destination_node_id=int(d_node_id)
        self.path_node_seq_no_list=list()
        self.path_link_seq_no_list=list()
        self.departure_time_in_min=float(departure_time_in_min)
        self.PCE_factor=int(PCE)  #Passenger Car Equivalent (PCE) of the agent
        self.path_cost=0
        self.m_bGenerated=False
        self.m_bCompleteTrip=False
        self.m_current_link_seq_no=0
        self.departure_time_in_simu_interval = int(self.departure_time_in_min*60.0 / g_number_of_seconds_per_interval + 0.5)
        
    def Initialization(self): # update the number of agents 
        global g_number_of_agents
        self.agent_seq_no=g_number_of_agents
        g_number_of_agents+=1
        
    def Initialize_for_simulation(self): 
        if(len(self.path_node_seq_no_list)!=0):      # if the agent can find its path 
            self.m_Veh_LinkArrivalTime_in_simu_interval=[-1]*(len(self.path_node_seq_no_list)-1)  # the arrival time on each link in the agent's path
            self.m_Veh_LinkDepartureTime_in_simu_interval=[-1]*(len(self.path_node_seq_no_list)-1) # the depature time on each link in the agent's path
            self.m_Veh_LinkArrivalTime_in_simu_interval[0]=self.departure_time_in_simu_interval 
            self.find_path_flag=True
        else:
            self.find_path_flag=False

# convert absolute simulation interval to relative simulation interval
def g_A2R_simu_interval(abslute_time):
    return (abslute_time-g_start_simu_interval_no)
# convert relative simulation interval to absolute simulation interval
def g_R2A_simu_interval(relative_time):
    return (relative_time+g_start_simu_interval_no)

def g_ReadInputData():
    global g_node_list
    global g_link_list
    global g_agent_list
    global g_Simulation_StartTimeInMin
    global g_Simulation_EndTimeInMin
    global g_start_simu_interval_no
    global g_end_simu_interval_no


    #read input_node
    with open('node.csv','r',encoding='utf-8') as fp:
        reader=csv.DictReader(fp)
        for line in reader:
            node=Node(line['node_id'])
            g_node_list.append(node)
    print('the number of nodes is',g_number_of_nodes)

    #read input_link
    with open('link.csv','r',encoding='utf-8') as fp:
        reader=csv.DictReader(fp)
        for line in reader:
            link=Link(line['from_node_id'],line['to_node_id'],
                      line['length'],line['lanes'],line['free_speed'],
                      line['capacity'],line['link_type'],line['VDF_alpha1'],
                      line['VDF_beta1'])
            g_node_list[link.from_node_seq_no].m_outgoing_link_list.append(link)
            g_node_list[link.to_node_seq_no].m_incoming_link_list.append(link)

            g_link_list.append(link)
    print('the number of links is',g_number_of_links)

    #read input_agent
    with open('demand.csv','r',encoding='utf-8') as fp:
        reader=csv.DictReader(fp)
        for line in reader:
            agent=Agent(line['agent_id'],line["agent_type"],line['o_node_id'],
                        line['d_node_id'],line['departure_time_in_min'],
                        line['PCE']
                       )
            #if the from_origin_node_id or to_destination_node_id is not defined in nodes
            if (agent.from_origin_node_id not in g_internal_node_seq_no_dict.keys() or
                agent.to_destination_node_id not in g_internal_node_seq_no_dict.keys()):
                continue
            #if the from_origin_node_id equals to_destination_node_id 
            if(agent.from_origin_node_id==agent.to_destination_node_id ):
                continue
            agent.Initialization()

            #update the g_Simulation_StartTimeInMin and g_Simulation_EndTimeInMin 
            if agent.departure_time_in_min<g_Simulation_StartTimeInMin:
                g_Simulation_StartTimeInMin=agent.departure_time_in_min
            if agent.departure_time_in_min>g_Simulation_EndTimeInMin:
                g_Simulation_EndTimeInMin=agent.departure_time_in_min

            #Add the agent to the time dependent agent list     
            if(agent.departure_time_in_simu_interval not in g_agent_TD_list_dict.keys()):
                g_agent_TD_list_dict[agent.departure_time_in_simu_interval]=list()
                g_agent_TD_list_dict[agent.departure_time_in_simu_interval].append(agent.agent_seq_no)
            else:
                g_agent_TD_list_dict[agent.departure_time_in_simu_interval].append(agent.agent_seq_no)
            g_agent_list.append(agent)
            
    print('the number of agents is',g_number_of_agents)

    #sort agents by the departure time
    sort_fun = operator.attrgetter("departure_time_in_min")
    g_agent_list.sort(key=sort_fun)

    for agent_no in range(len(g_agent_list)):
        g_agent_list[agent_no].agent_seq_no=agent_no
        
    #start simulation interval and end simulation interval
    g_start_simu_interval_no =int(g_Simulation_StartTimeInMin * 60 / g_number_of_seconds_per_interval)
    g_end_simu_interval_no =g_start_simu_interval_no + g_number_of_simulation_intervals
    

class Network:
    global _MAX_LABEL_COST
    global g_number_of_links
    global g_agent_list 
    def allocate(self,number_of_nodes,number_of_links):
        self.node_predecessor = [-1]*number_of_nodes
        self.node_label_cost = [_MAX_LABEL_COST]*number_of_nodes
        self.link_predecessor = [-1]*number_of_nodes
        self.link_cost_array = [1.0]*number_of_links
        self.link_volume_array = [0.0]*number_of_links

        for j in range(number_of_links):
            self.link_volume_array[j] = 0.0
            self.link_cost_array[j] = 1.0
        
    def optimal_label_correcting(self,origin_node,destination_node,departure_time):
        origin_node=g_internal_node_seq_no_dict[origin_node]
        destination_node=g_internal_node_seq_no_dict[destination_node]
        if len(g_node_list[origin_node].m_outgoing_link_list) == 0:
            return 0
        
        for i in range(g_number_of_nodes): #Initialization for all nodes
            self.node_label_cost[i] = _MAX_LABEL_COST
            self.node_predecessor[i] = -1 #pointer to previous node index from the current label at current node and time
            self.link_predecessor[i] = -1 #pointer to previous node index from the current label at current node and time        
        
        self.node_label_cost[origin_node] = departure_time
        SEList = []
        SEList.append(origin_node)

        while len(SEList)>0:
            from_node = SEList[0]
            del SEList[0]
            for k in range(len(g_node_list[from_node].m_outgoing_link_list)):
                to_node = g_node_list[from_node].m_outgoing_link_list[k].to_node_seq_no
             
                new_to_node_cost = self.node_label_cost[from_node] + self.link_cost_array[g_node_list[from_node].m_outgoing_link_list[k].link_seq_no]
                
                if (new_to_node_cost < self.node_label_cost[to_node]):  #we only compare cost at the downstream node ToID at the new arrival time t
                    # update cost label and node/time predecessor
                    self.node_label_cost[to_node] = new_to_node_cost
                    self.node_predecessor[to_node] = from_node  #pointer to previous physical NODE INDEX from the current label at current node and time
                    self.link_predecessor[to_node] = g_node_list[from_node].m_outgoing_link_list[k].link_seq_no  #pointer to previous physical NODE INDEX from the current label at current node and time
                    SEList.append(to_node)
                                        

        if (destination_node >= 0 and self.node_label_cost[destination_node] < _MAX_LABEL_COST):
            return 1
        else: 
            return -1



    def find_path_for_agents(self,iteration_no):
            
        for s in range(g_number_of_links):
            self.link_volume_array[s] = 0

        #step 1: find shortest path if needed 
        for i in range(len(g_agent_list)):
            residual = i % (iteration_no + 1)
            if (residual != 0):  #no need to compute a new path at this iteration
                continue         #that is, it will reuse the path from the previous iteration, stored at p_agent->path_link_seq_no_list.
            #else move to the next line for finding the shortest path 
            g_agent_list[i].path_link_seq_no_list = []
            g_agent_list[i].path_node_seq_no_list = []
            #step 2 buil SP tree
            
            return_value = self.optimal_label_correcting(g_agent_list[i].from_origin_node_id, g_agent_list[i].to_destination_node_id, g_agent_list[i].departure_time_in_min)
            #step 3 update path

            if (return_value == -1):
                #print('agent ',i,'can not find destination node')
                continue

            current_node_seq_no = g_internal_node_seq_no_dict[g_agent_list[i].to_destination_node_id]
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
def g_TrafficAssignment():
    global g_link_list
    
    network = Network()
    network.allocate(g_number_of_nodes,g_number_of_links)

    print('Finding shortest path for all agents......')

    for i in range(g_number_of_assignment_iterations):
        print('iteration_no',i,'......')
        for l in range(g_number_of_links):
            g_link_list[l].CalculateBRPFunction()
            network.link_cost_array[l] = g_link_list[l].cost

            network.find_path_for_agents(i)     
                        
            for k in range(g_number_of_links):
                g_link_list[k].flow_volume = 0.0
                g_link_list[k].flow_volume += network.link_volume_array[k]
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
def g_TrafficSimulation():
    global g_Cumulative_Arrival_Count
    global g_Cumulative_Departure_Count

    #initialization for each link
    for li in range(len(g_link_list)):
        g_link_list[li].ResetMOE()
    #initialization for each agent
    for agent_no in range(len(g_agent_list)):
        g_agent_list[agent_no].Initialize_for_simulation()

    g_active_agent_list=list()
    current_active_agent_id = 0
    
    for t in range(g_start_simu_interval_no,g_end_simu_interval_no,1):
        number_of_simu_interval_per_min = 60 / g_number_of_seconds_per_interval;
        if(t%number_of_simu_interval_per_min==0):
            print("simu time= ",int(t / number_of_simu_interval_per_min),"min, ", "CA=",g_Cumulative_Arrival_Count,", CD=",g_Cumulative_Departure_Count)

        relative_t = g_A2R_simu_interval(t)
        for li in range(0,len(g_link_list)):
            link=g_link_list[li]
            if (relative_t >= 1):
                g_link_list[link.link_seq_no].m_LinkCumulativeDeparture[relative_t]=g_link_list[link.link_seq_no].m_LinkCumulativeDeparture[relative_t-1]
                g_link_list[link.link_seq_no].m_LinkCumulativeArrival[relative_t] =g_link_list[link.link_seq_no].m_LinkCumulativeArrival[relative_t-1]
	
      
  
        if(t in g_agent_TD_list_dict.keys()):  #activate the agent 
            for j in range(0,len(g_agent_TD_list_dict[t])):
                agent_no=g_agent_TD_list_dict[t][j]
                agent=g_agent_list[agent_no]
                if(agent.find_path_flag==True):
                    agent.m_bGenerated=True
                    first_link_seq=agent.path_link_seq_no_list[0]
                    g_link_list[first_link_seq].m_LinkCumulativeArrival[relative_t]+=1
                    g_link_list[first_link_seq].Entrance_queue.append(agent.agent_seq_no)
                    g_Cumulative_Arrival_Count+=1
                
        for li in range(0,len(g_link_list)):
            link=g_link_list[li]
            while(len(link.Entrance_queue)>0):  #there are agents in the Entrance_queue
                agent_seq=link.Entrance_queue[0]
                link.Entrance_queue.pop(0)
                link.Exit_queue.append(agent_seq)
                g_agent_list[agent_seq].m_Veh_LinkDepartureTime_in_simu_interval[g_agent_list[agent_seq].m_current_link_seq_no]=\
                    g_agent_list[agent_seq].m_Veh_LinkArrivalTime_in_simu_interval[g_agent_list[agent_seq].m_current_link_seq_no]+link.general_travel_time_in_simu_interval

        for no in range(0,len(g_node_list)):
            node=g_node_list[no]
            for l in range(0,len(node.m_incoming_link_list)):
                incoming_link_index=(l+t)%(len(node.m_incoming_link_list)) 
                #we will start with different first link from the incoming link list, 
                #equal change, regardless of # of lanes or main line vs. ramp, but one can use service arc, to control the effective capacity rates, e.g. through a metered ramp, to 
		        #allow mainline to use the remaining flow

                link_seq_no=node.m_incoming_link_list[incoming_link_index].link_seq_no

                #check if the current link has sufficient capacity
                #check if there are agents in exit queue
             
                while (g_link_list[link_seq_no].m_LinkOutFlowCapacity[relative_t]>=1 and len(g_link_list[link_seq_no].Exit_queue)>=1):
                    agent_no=node.m_incoming_link_list[incoming_link_index].Exit_queue[0]
                   
                    if(g_agent_list[agent_no].m_Veh_LinkDepartureTime_in_simu_interval[g_agent_list[agent_no].m_current_link_seq_no]>t):
                        break        #the future departure time on this link is later than the current time
               
                    if(g_agent_list[agent_no].m_current_link_seq_no==len(g_agent_list[agent_no].path_link_seq_no_list)-1): # end of the path
                        node.m_incoming_link_list[incoming_link_index].Exit_queue.pop(0)
                        g_agent_list[agent_no].m_bCompleteTrip = True
                        g_link_list[link_seq_no].m_LinkCumulativeDeparture[relative_t] += 1
                        g_Cumulative_Departure_Count += 1
                    else:  #  not complete the trip. move to the next link's entrance queue
                        next_link_seq=g_agent_list[agent_no].path_link_seq_no_list[g_agent_list[agent_no].m_current_link_seq_no+1]
                        node.m_incoming_link_list[incoming_link_index].Exit_queue.pop(0)
                        g_link_list[next_link_seq].Entrance_queue.append(agent_no)
                        g_agent_list[agent_no].m_Veh_LinkDepartureTime_in_simu_interval[g_agent_list[agent_no].m_current_link_seq_no]=t
                        g_agent_list[agent_no].m_Veh_LinkArrivalTime_in_simu_interval[g_agent_list[agent_no].m_current_link_seq_no+1]=t
                        
                        actual_travel_time = t - g_agent_list[agent_no].m_Veh_LinkArrivalTime_in_simu_interval[g_agent_list[agent_no].m_current_link_seq_no]
			#for each waited vehicle
                        waiting_time = actual_travel_time - g_link_list[link_seq_no].general_travel_time_in_min
                        temp_relative_time=g_A2R_simu_interval(g_agent_list[agent_no].m_Veh_LinkArrivalTime_in_simu_interval[g_agent_list[agent_no].m_current_link_seq_no])
                        time_in_min=int(temp_relative_time/g_number_of_interval_per_min)
                  
                        g_link_list[link_seq_no].m_LinkTDWaitingTime[time_in_min] += waiting_time
                        g_link_list[link_seq_no].m_LinkCumulativeDeparture[relative_t] += 1
                        g_link_list[next_link_seq].m_LinkCumulativeArrival[relative_t] += 1
                    
                    g_agent_list[agent_no].m_current_link_seq_no+=1
                    g_link_list[link_seq_no].m_LinkOutFlowCapacity[relative_t]-=1 
        
        
def g_OutputFiles():
    with open ('link_performance.csv','w',newline='') as fp:
        writer=csv.writer(fp)
        line=["link_id","from_node_id","to_node_id","time_period","volume","CA","CD","density","queue","travel_time","waiting_time_in_sec","speed"]
        writer.writerow(line)
        
        for link in g_link_list:
            for relative_t in range(0,g_number_of_simulation_intervals,1):
                abs_t = g_R2A_simu_interval(relative_t)
                
##                if(relative_t%(60/g_number_of_seconds_per_interval)==0):
                time_in_min =  int(relative_t / (60 / g_number_of_seconds_per_interval))
                abs_t_in_min= int(abs_t / (60 / g_number_of_seconds_per_interval))
            
                volume = 0
                queue = 0
                waiting_time_in_sec = 0
                arrival_rate = 0
                avg_waiting_time_in_sec = 0
                avg_travel_time_in_min =float(link.general_travel_time_in_min + avg_waiting_time_in_sec/60.0)
                speed = link.length /(max(0.00001,avg_travel_time_in_min/60.0))
                virtual_arrival = 0
                if (time_in_min >= 1):
                
                    volume = link.m_LinkCumulativeDeparture[relative_t] - link.m_LinkCumulativeDeparture[relative_t - int(60 / g_number_of_seconds_per_interval)]

                    if (relative_t- link.general_travel_time_in_simu_interval>= 0):  
                        virtual_arrival = link.m_LinkCumulativeArrival[int(relative_t - link.general_travel_time_in_simu_interval)]
                                                                            
                    queue = virtual_arrival - link.m_LinkCumulativeDeparture[relative_t]
                                                    

                    waiting_time_count = 0

                    waiting_time_in_sec = link.m_LinkTDWaitingTime[time_in_min] * g_number_of_seconds_per_interval 
                    if(relative_t + int(60 / g_number_of_seconds_per_interval)<g_number_of_simulation_intervals):
                        arrival_rate = link.m_LinkCumulativeArrival[relative_t + int(60 / g_number_of_seconds_per_interval)] - link.m_LinkCumulativeArrival[relative_t]
                    avg_waiting_time_in_sec = waiting_time_in_sec / max(1, arrival_rate)

                    avg_travel_time_in_min = link.general_travel_time_in_min+ avg_waiting_time_in_sec/60.0
                    speed = link.length / (max(0.00001,avg_travel_time_in_min) / 60)
                    
                density = (link.m_LinkCumulativeArrival[relative_t] - link.m_LinkCumulativeDeparture[relative_t]) / (link.length * link.lanes);

                line=[link.link_seq_no+1,link.external_from_node,link.external_to_node,time_stamp_to_HHMMSS(abs_t_in_min)+"_"+time_stamp_to_HHMMSS(abs_t_in_min+1),
                      volume,link.m_LinkCumulativeArrival[relative_t],link.m_LinkCumulativeDeparture[relative_t],
                      density,queue,avg_travel_time_in_min,waiting_time_in_sec,speed]
                
                writer.writerow(line)
    with open("agent.csv","w",newline='') as fp:
        writer=csv.writer(fp)
        line=["agent_id","agent_type","o_node_id","d_node_id","o_zone_id","d_zone_id","travel_time",
              "node_sequence","time_sequence"]
        writer.writerow(line)
        abs_satrt_time_in_min=int(g_start_simu_interval_no/g_number_of_interval_per_min )
        
        for agent in g_agent_list:
            if(agent.find_path_flag==False):
                continue
            path_time_list=list()
            cost=(agent.m_Veh_LinkDepartureTime_in_simu_interval[-1]-agent.m_Veh_LinkArrivalTime_in_simu_interval[0])\
                * g_number_of_seconds_per_interval / 60.0
            for i in range(0,len(agent.path_link_seq_no_list),1):
                TA_in_min=agent.m_Veh_LinkArrivalTime_in_simu_interval[i]* g_number_of_seconds_per_interval / 60.0
                TD_in_min=agent.m_Veh_LinkDepartureTime_in_simu_interval[i]* g_number_of_seconds_per_interval / 60.0
                if(i==0):
                    
                    path_time_list.extend([time_stamp_to_HHMMSS(TA_in_min),time_stamp_to_HHMMSS(TD_in_min)])
                else:
                    path_time_list.append(time_stamp_to_HHMMSS(TD_in_min))
                path_time_sequence=';'.join(path_time_list)
            path_node_str=';'.join([str(g_external_node_id_dict[i]) for i in agent.path_node_seq_no_list])
            line=[agent.agent_id,agent.agent_type,agent.from_origin_node_id,
                  agent.to_destination_node_id,agent.from_origin_node_id,
                  agent.to_destination_node_id,cost,
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
    g_ReadInputData()
    g_TrafficAssignment()
    g_TrafficSimulation()

    end_time=time()
    g_OutputFiles()
    

    
