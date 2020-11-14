import csv
from time import *
import operator
max_number_of_node=10000   
_MAX_LABEL_COST = 10000

g_number_of_assgnment_iterations=1
g_number_of_simulation_time=1440     #min
g_number_of_seconds_per_interval = 6;  #0.2 seconds for 300 intervals per min
g_number_of_simulation_intervals=int(g_number_of_simulation_time*60/g_number_of_seconds_per_interval)
g_Simulation_StartTimeInMin = 9999
g_Simulation_EndTimeInMin = 0
g_start_simu_interval_no=0
g_end_simu_interval_no=0

g_number_of_nodes=0
g_number_of_links=0
g_number_of_schedules=0
g_number_of_agents=0
g_number_of_overtime_agents=0
g_number_of_threads=4

g_node_list=list()
g_link_list=list()
g_schedule_list=list()
g_agent_list=list()
g_line_list=list()


g_internal_node_seq_no_dict=dict()
g_internal_node_seq_no_to_node_id_dict=dict()
g_link_key_to_seq_no_dict=dict()
g_map_agent_id_to_agent_seq_no=dict()
g_agent_TD_list_dict=dict()     #key:simulation time interval   value:agents(list) need to be activated



g_TotalCumulative_Arrival_Count=0
g_TotalCumulative_Departure_Count = 0

class Node:
    def __init__(self,external_node_id):
        self.node_seq_no=0
        self.external_node_id=external_node_id
        self.m_outgoing_link_list=list()
        self.m_incoming_link_list=list()
        self.Initialization()

    def Initialization(self):
        global g_number_of_nodes
        self.node_seq_no=g_number_of_nodes
        g_internal_node_seq_no_dict[int(self.external_node_id)]=self.node_seq_no
        g_internal_node_seq_no_to_node_id_dict[int(self.node_seq_no)]=self.external_node_id
        g_number_of_nodes+=1
        

class Link:
    def __init__(self,from_node,to_node,length,number_of_lanes,speed_limit,lane_cap,
                 link_type,BPR_alpha_term,BPR_beta_term):
        self.from_node_seq_no=g_internal_node_seq_no_dict[int(from_node)]
        self.to_node_seq_no=g_internal_node_seq_no_dict[int(to_node)]
        self.extern_from_node=int(from_node)
        self.extern_to_node=int(to_node)
        self.type=int(link_type) # 1:one direction 2:both way
        
        self.BRP_alpha=float(BPR_alpha_term)
        self.BRP_beta=float(BPR_beta_term)
        self.flow_volume=0
        self.link_capacity = float(lane_cap)* int(number_of_lanes)
        self.length=float(length)
        self.free_flow_travel_time_in_min=self.length / max(0.001,int(speed_limit)) * 60
        self.travel_time=0
    
        self.cost = self.length / int(speed_limit) * 60
      
        self.m_LinkOutFlowCapacity=[self.link_capacity]*(g_number_of_simulation_intervals+1)
        self.m_LinkCumulativeArrival=[0]*(g_number_of_simulation_intervals+1)
        self.m_LinkCumulativeDeparture=[0]*(g_number_of_simulation_intervals+1)
        self.m_LinkCumulativeVirtualDelay=[0]*(g_number_of_simulation_intervals+1)
        
        self.Entrance_queue=list() # link-in queue  of each link
        self.Exit_queue=list() # link-out queue  of each link
        

        self.Initialization()

    def Initialization(self):
        global g_number_of_links
        self.link_seq_no=g_number_of_links
        self.travel_time = self.free_flow_travel_time_in_min*(1 + self.BRP_alpha*pow(self.flow_volume / max(0.00001, self.link_capacity), self.BRP_beta))
        self.cost=self.travel_time
        self.free_flow_travel_time_in_simu_interval = int(self.free_flow_travel_time_in_min*60.0 / g_number_of_seconds_per_interval + 0.5)
        link_key=self.from_node_seq_no*max_number_of_node+self.to_node_seq_no
        g_link_key_to_seq_no_dict[link_key]=self.link_seq_no
        g_number_of_links +=1
    def ResetMOE(self):
        self.m_CumulativeArrivalCount = 0
        self.m_CumulativeDepartureCount = 0
        self.m_CumulativeVirtualDelayCount = 0
    def CalculateBRPFunction(self):
        self.travel_time = self.free_flow_travel_time_in_min*(1 + self.BRP_alpha*((self.flow_volume / max(0.00001, self.link_capacity))**self.BRP_beta))
        self.cost = self.travel_time
        

class Agent():
    def __init__(self,agent_id,agent_service_type,from_origin_node_id,
                 to_destination_node_id,departure_time_in_min,PCE,
                 ):
        self.agent_id=agent_id
       
        self.agent_service_type=agent_service_type
        self.from_origin_node_id=int(from_origin_node_id)
        self.to_destination_node_id=int(to_destination_node_id)
        self.path_node_seq_no_list=list()
        self.path_link_seq_no_list=list()
        self.departure_time_in_min=float(departure_time_in_min)

        self.PCE_factor=int(PCE)
        self.path_cost=0

        self.m_bMoveable=True
        self.m_bGenerated=False
        self.bActive = True
        self.m_bCompleteTrip=False
        
        self.m_current_link_seq_no=0

        

    def Initialization(self):
        global g_number_of_agents
        self.agent_seq_no=g_number_of_agents
        g_number_of_agents+=1
        self.departure_time_in_simu_interval = int(self.departure_time_in_min*60.0 / g_number_of_seconds_per_interval + 0.5);
    def Initialize_for_simulation(self):
        if(len(self.path_node_seq_no_list)!=0):
            self.m_Veh_LinkArrivalTime_in_simu_interval=[-1]*(len(self.path_node_seq_no_list)-1)
            self.m_Veh_LinkDepartureTime_in_simu_interval=[-1]*(len(self.path_node_seq_no_list)-1)
            self.m_Veh_LinkArrivalTime_in_simu_interval[0]=self.departure_time_in_simu_interval 
            self.find_path_flag=True
        else:
            self.find_path_flag=False

#g_convert_abs_simu_interval_to_relative_simu_interval
def g_A2R_simu_interval(abslute_time):
    return (abslute_time-g_start_simu_interval_no)

def g_ReadInputData():
    global g_node_list
    global g_link_list
    global g_agent_list
    global g_Simulation_StartTimeInMin
    global g_Simulation_EndTimeInMin
    global g_start_simu_interval_no
    global g_end_simu_interval_no
    global g_map_agent_id_to_agent_seq_no

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
                      line['capacity'],line['link_type'],line['BPR_alpha_term'],
                      line['BPR_beta_term'])
            g_node_list[link.from_node_seq_no].m_outgoing_link_list.append(link)
            g_node_list[link.to_node_seq_no].m_incoming_link_list.append(link)

            g_link_list.append(link)
    print('the number of links is',g_number_of_links)

    #read input_agent
    with open('input_agent.csv','r',encoding='utf-8') as fp:
        reader=csv.DictReader(fp)
        for line in reader:
            agent=Agent(line['agent_id'],line["agent_type"],line['o_node_id'],
                        line['d_node_id'],line['departure_time_in_min'],
                        line['PCE']
                       )

            if (agent.from_origin_node_id not in g_internal_node_seq_no_dict.keys() or
                agent.to_destination_node_id not in g_internal_node_seq_no_dict.keys()):
                continue
            if(agent.from_origin_node_id==agent.to_destination_node_id ):
                continue
            agent.Initialization()

            if agent.departure_time_in_min<g_Simulation_StartTimeInMin:
                g_Simulation_StartTimeInMin=agent.departure_time_in_min
            if agent.departure_time_in_min>g_Simulation_EndTimeInMin:
                g_Simulation_EndTimeInMin=agent.departure_time_in_min
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
        g_map_agent_id_to_agent_seq_no[g_agent_list[agent_no].agent_id]=agent_no
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
            self.node_predecessor[i] = -1 #ointer to previous node index from the current label at current node and time
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
                continue #that is, it will reuse the path from the previous iteration, stored at p_agent->path_link_seq_no_list.
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

    for i in range(g_number_of_assgnment_iterations):
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
    global g_TotalCumulative_Arrival_Count
    global g_TotalCumulative_Departure_Count

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
            print("simu time= ",int(t / number_of_simu_interval_per_min),"min, with TBL= ",len(g_active_agent_list),
                  ", CA=",g_TotalCumulative_Arrival_Count,", CD=",g_TotalCumulative_Departure_Count)

        relative_t = g_A2R_simu_interval(t);
  
        if(t in g_agent_TD_list_dict.keys()):  #activate the agent 
            for j in range(0,len(g_agent_TD_list_dict[t])):
                agent_no=g_agent_TD_list_dict[t][j]
                agent=g_agent_list[agent_no]
                if(agent.find_path_flag==True):
                    agent.m_bGenerated=True
                    first_link_seq=agent.path_link_seq_no_list[0]
                    g_link_list[first_link_seq].m_LinkCumulativeArrival[relative_t]+=1
                    g_link_list[first_link_seq].Entrance_queue.append(agent.agent_seq_no)
                    g_TotalCumulative_Arrival_Count+=1
                
        for li in range(0,len(g_link_list)):
            link=g_link_list[li]
            while(len(link.Entrance_queue)>0):  #there are agents in the Entrance_queue
                agent_seq=link.Entrance_queue[0]
                link.Entrance_queue.pop(0)
                link.Exit_queue.append(agent_seq)
                g_agent_list[agent_seq].m_Veh_LinkDepartureTime_in_simu_interval[g_agent_list[agent_seq].m_current_link_seq_no]=\
                    g_agent_list[agent_seq].m_Veh_LinkArrivalTime_in_simu_interval[g_agent_list[agent_seq].m_current_link_seq_no]+link.travel_time

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
                        g_TotalCumulative_Departure_Count += 1
                    else:  #  not complete the trip. move to the next link's entrance queue
                        next_link_seq=g_agent_list[agent_no].path_link_seq_no_list[g_agent_list[agent_no].m_current_link_seq_no+1]
                        node.m_incoming_link_list[incoming_link_index].Exit_queue.pop(0)
                        g_link_list[next_link_seq].Entrance_queue.append(agent_no)
                        g_agent_list[agent_no].m_Veh_LinkDepartureTime_in_simu_interval[g_agent_list[agent_no].m_current_link_seq_no]=t
                        g_agent_list[agent_no].m_Veh_LinkArrivalTime_in_simu_interval[g_agent_list[agent_no].m_current_link_seq_no+1]=t
                        g_link_list[link_seq_no].m_LinkCumulativeDeparture[relative_t] += 1
                        g_link_list[next_link_seq].m_LinkCumulativeArrival[relative_t] += 1
               
                    g_agent_list[agent_no].m_current_link_seq_no+=1
                    g_link_list[link_seq_no].m_LinkOutFlowCapacity[relative_t]-=1 
        
        
def g_OutputFiles():
    with open ('link_performance.csv','w',newline='') as fp:
        writer=csv.writer(fp)
        line=["o_node_id","d_node_id","cumulative_arrival_count","cumulative_departure_count",
              "travel_time_in_min"]
        writer.writerow(line)
        for t in range( int(g_Simulation_StartTimeInMin),int(g_Simulation_EndTimeInMin),1):
            time_in_simu_interval = int(t * 60 / g_number_of_seconds_per_interval)
            relative_time_interval=g_A2R_simu_interval(time_in_simu_interval)
            for link in g_link_list:
                line=[link.extern_from_node,link.extern_to_node,link.m_LinkCumulativeArrival[relative_time_interval],
                      link.m_LinkCumulativeDeparture[relative_time_interval],link.travel_time]
                writer.writerow(line)
    with open("agent.csv","w",newline='') as fp:
        writer=csv.writer(fp)
        line=["agent_id","agent_type","o_node_id","d_node_id","o_zone_id","d_zone_id","travel_time",
              "node_sequence","time_sequence"]
        writer.writerow(line)
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
            path_node_str=';'.join([str(i) for i in agent.path_node_seq_no_list])
            line=[agent.agent_id,agent.agent_service_type,agent.from_origin_node_id,
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
    

    
