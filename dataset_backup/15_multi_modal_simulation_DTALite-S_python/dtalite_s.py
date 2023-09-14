import csv
from time import *
import operator
max_number_of_node=10000   


g_number_of_simulation_time=120     #min
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
    def __init__(self,from_node,to_node,service_type,length,number_of_lanes,speed_limit,lane_cap,
                 link_type,BPR_alpha_term,BPR_beta_term,external_travel_time):
        self.from_node_seq_no=g_internal_node_seq_no_dict[int(from_node)]
        self.to_node_seq_no=g_internal_node_seq_no_dict[int(to_node)]
        self.extern_from_node=int(from_node)
        self.extern_to_node=int(to_node)
        self.type=int(link_type) # 1:one direction 2:both way
        self.service_type=int(service_type)  #0:moving link  1:drop-off link  2:pick-up link
        self.BRP_alpha=float(BPR_alpha_term)
        self.BRP_beta=float(BPR_beta_term)
        self.flow_volume=0
        self.link_capacity = int(lane_cap)* int(number_of_lanes)
        self.length=float(length)
        self.free_flow_travel_time_in_min=self.length / max(0.001,int(speed_limit)) * 60
        self.travel_time=0
    
        self.cost = self.length / int(speed_limit) * 60
      
        self.m_LinkOutFlowCapacity=[self.link_capacity]*g_number_of_simulation_intervals
        self.m_LinkCumulativeArrival=[0]*g_number_of_simulation_intervals
        self.m_LinkCumulativeDeparture=[0]*g_number_of_simulation_intervals
        self.m_LinkCumulativeVirtualDelay=[0]*g_number_of_simulation_intervals
        self.m_waiting_traveler_queue=list()

        if(float(external_travel_time)>0):
            self.free_flow_travel_time_in_min=float(external_travel_time)

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

class Agent():
    def __init__(self,agent_id,agent_service_type,vehicle_seat_capacity,from_origin_node_id,
                 to_destination_node_id,fixed_path_flag,PCE,path_node_sequence,
                 departure_time_in_min,path_schedule_time_sequence):
        self.agent_id=agent_id
        self.agent_service_type=int(agent_service_type)  #1:passenger   2:vehicle
        self.agent_seq_no=0
        self.from_origin_node_id=int(from_origin_node_id)
        self.to_destination_node_id=int(to_destination_node_id)
        self.path_node_str=path_node_sequence
        self.path_node_sequence=[int(node) for node in path_node_sequence.strip().split(';') if node!='']
        if(len(path_schedule_time_sequence)>0):
            self.path_schedule_time_sequence=[int(time) for time in path_schedule_time_sequence.strip().split(';')if time!='']
        self.path_link_seq_no_list=list()
        self.departure_time_in_min=float(departure_time_in_min)
        if(self.agent_service_type==2):
            self.vehicle_seat_capacity=int(vehicle_seat_capacity)

        self.path_cost=0

        self.m_bMoveable=True
        self.m_bGenerated=False
        self.bActive = True
        self.m_bCompleteTrip=False
        self.m_PassengerList=list()
        self.m_Veh_LinkArrivalTime_in_simu_interval=[-1]*(len(self.path_node_sequence)-1)
        self.m_Veh_LinkDepartureTime_in_simu_interval=[-1]*(len(self.path_node_sequence)-1)
        self.m_current_link_seq_no=0

        self.Initialization()

    def Initialization(self):
        global g_number_of_agents
        g_number_of_agents+=1
    def InitialStateSet(self):
        self.departure_time_in_simu_interval = int(self.departure_time_in_min*60.0 / g_number_of_seconds_per_interval + 0.5); 
        if(len(self.path_node_sequence)>0):
            self.m_Veh_LinkArrivalTime_in_simu_interval[0] = self.departure_time_in_simu_interval
            FirstLink =self.path_link_seq_no_list[0]
            self.m_Veh_LinkDepartureTime_in_simu_interval[0] = self.m_Veh_LinkArrivalTime_in_simu_interval[0] \
                + g_link_list[FirstLink].free_flow_travel_time_in_simu_interval;
            relative_simulation_interval = g_A2R_simu_interval(self.departure_time_in_simu_interval)
            g_link_list[FirstLink].m_LinkCumulativeArrival[relative_simulation_interval] += 1;
    def Pickup(self,p):
        if(len(self.m_PassengerList)<self.vehicle_seat_capacity):
            self.m_PassengerList.append(p)
    def GetRemainingCapacity(self):
        return (self.vehicle_seat_capacity-len(self.m_PassengerList))

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
    with open('input_node.csv','r',encoding='utf-8') as fp:
        reader=csv.DictReader(fp)
        for line in reader:
            node=Node(line['node_id'])
            g_node_list.append(node)
    print('the number of nodes is',g_number_of_nodes)

    #read input_link
    with open('input_link.csv','r',encoding='utf-8') as fp:
        reader=csv.DictReader(fp)
        for line in reader:
            link=Link(line['from_node_id'],line['to_node_id'],line['service_type'],
                      line['length'],line['number_of_lanes'],line['speed_limit'],
                      line['lane_cap'],line['link_type'],line['BPR_alpha_term'],
                      line['BPR_beta_term'],line['external_travel_time'])
            g_node_list[link.from_node_seq_no].m_outgoing_link_list.append(link)
            g_node_list[link.to_node_seq_no].m_incoming_link_list.append(link)

            g_link_list.append(link)
    print('the number of links is',g_number_of_links)

    #read input_agent
    with open('input_agent.csv','r',encoding='utf-8') as fp:
        reader=csv.DictReader(fp)
        for line in reader:
            agent=Agent(line['agent_id'],line['agent_service_type'],line['vehicle_seat_capacity'],
                        line['from_origin_node_id'],line['to_destination_node_id'],line['fixed_path_flag'],
                        line['PCE'],line['path_node_sequence'],line['departure_time_in_min'],
                        line['path_schedule_time_sequence'])

            if (agent.from_origin_node_id not in g_internal_node_seq_no_dict.keys() or
                agent.to_destination_node_id not in g_internal_node_seq_no_dict.keys()):
                continue

            if agent.departure_time_in_min<g_Simulation_StartTimeInMin:
                g_Simulation_StartTimeInMin=agent.departure_time_in_min
            if agent.departure_time_in_min>g_Simulation_EndTimeInMin:
                g_Simulation_EndTimeInMin=agent.departure_time_in_min
            
            if len(agent.path_node_sequence)>=2:
                for i in range(len(agent.path_node_sequence)-1):
                    from_node_seq=g_internal_node_seq_no_dict[agent.path_node_sequence[i]]
                    to_node_seq=g_internal_node_seq_no_dict[agent.path_node_sequence[i+1]]
                    link_key=from_node_seq*max_number_of_node+to_node_seq
                    link_seq_no=g_link_key_to_seq_no_dict[link_key]
                    agent.path_link_seq_no_list.append(link_seq_no)
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


def g_TrafficSimulation():
    global g_TotalCumulative_Arrival_Count
    global g_TotalCumulative_Departure_Count

    #initialization for each link
    for li in range(len(g_link_list)):
        g_link_list[li].ResetMOE()
    #initialization for each agent
    for agent_no in range(len(g_agent_list)):
        g_agent_list[agent_no].InitialStateSet()

    g_active_agent_list=list()
    current_active_agent_id = 0
    
    for t in range(g_start_simu_interval_no,g_end_simu_interval_no,1):
        number_of_simu_interval_per_min = 60 / g_number_of_seconds_per_interval;
        if(t%number_of_simu_interval_per_min==0):
            print("simu time= ",t / number_of_simu_interval_per_min,"min, with TBL= ",len(g_active_agent_list),
                  ", CA=",g_TotalCumulative_Arrival_Count,", CD=",g_TotalCumulative_Departure_Count)

        relative_t = g_A2R_simu_interval(t);
        if(t%number_of_simu_interval_per_min==0):  #unpdate the active agent list per min
            for agent_no in range(current_active_agent_id,len(g_agent_list),1):
                agent=g_agent_list[agent_no]
                if(t<=agent.departure_time_in_simu_interval<=t + number_of_simu_interval_per_min):
                    agent.m_bGenerated=True
                    agent.bActive = True
                    g_TotalCumulative_Arrival_Count+=1
                    g_active_agent_list.append(agent)
                    current_active_agent_id+=1
                else:
                    break
        for process in range(g_number_of_threads):  #1 Todo: parallel     2 deal with the same time movement(passenger first,vehicle second)
            for agent in g_active_agent_list:
                active_agent = g_agent_list[agent.agent_seq_no]
                if(active_agent.m_bGenerated == True and 
                   active_agent.m_bCompleteTrip == False and 
                   active_agent.m_bMoveable == True and 
                   active_agent.m_current_link_seq_no%g_number_of_threads==process and 
                   active_agent.m_current_link_seq_no < len(active_agent.path_link_seq_no_list)):
                    #ready to move
                    if(active_agent.m_Veh_LinkDepartureTime_in_simu_interval[active_agent.m_current_link_seq_no] == t):
                        link_seq_no = active_agent.path_link_seq_no_list[active_agent.m_current_link_seq_no]
                        #condition1:moving link
                        if(g_link_list[link_seq_no].service_type == 0):
                            #check check if the current link has sufficient
                            #capacity
                            if(g_link_list[link_seq_no].m_LinkOutFlowCapacity[relative_t] > 0):
                                if(active_agent.m_current_link_seq_no == len(active_agent.path_link_seq_no_list) - 1):
                                    #end of the path
                                    active_agent.m_bCompleteTrip = True
                                    active_agent.bActive = False
                                    g_link_list[link_seq_no].m_LinkCumulativeDeparture[relative_t] += 1
                                    g_TotalCumulative_Departure_Count += 1
                                else:   #not complete the trip
                                    next_link_seq_no = active_agent.path_link_seq_no_list[active_agent.m_current_link_seq_no + 1]
                                    active_agent.m_Veh_LinkArrivalTime_in_simu_interval[active_agent.m_current_link_seq_no + 1] = t
                                    active_agent.m_Veh_LinkDepartureTime_in_simu_interval[active_agent.m_current_link_seq_no + 1] = \
                                        t + g_link_list[next_link_seq_no].free_flow_travel_time_in_simu_interval
                                    for agent_id in active_agent.m_PassengerList:
                                        agent_seq_no = g_map_agent_id_to_agent_seq_no[agent_id]
                                        current_link_seq = g_agent_list[agent_seq_no].m_current_link_seq_no
                                        g_agent_list[agent_seq_no].m_Veh_LinkDepartureTime_in_simu_interval[current_link_seq] = t
                                        g_agent_list[agent_seq_no].m_Veh_LinkArrivalTime_in_simu_interval[current_link_seq + 1] = t
                                        g_agent_list[agent_seq_no].m_Veh_LinkDepartureTime_in_simu_interval[current_link_seq + 1] = \
                                            active_agent.m_Veh_LinkDepartureTime_in_simu_interval[active_agent.m_current_link_seq_no + 1]
                                        g_agent_list[agent_seq_no].m_current_link_seq_no+=1
                                    
                                    g_link_list[link_seq_no].m_CumulativeDepartureCount += 1
                                    g_link_list[link_seq_no].m_LinkCumulativeDeparture[relative_t] = g_link_list[link_seq_no].m_CumulativeDepartureCount
                                    g_link_list[next_link_seq_no].m_CumulativeArrivalCount += 1
                                    g_link_list[next_link_seq_no].m_LinkCumulativeArrival[relative_t] = g_link_list[next_link_seq_no].m_CumulativeArrivalCount
                                #move
                                active_agent.m_current_link_seq_no+=1
                                g_link_list[link_seq_no].m_LinkOutFlowCapacity[relative_t] -= 1
                            else:  #no outflow capacity ==0, cause delay
                                g_link_list[link_seq_no].m_CumulativeVirtualDelayCount += 1  #+1 means add one unit of simulation time interval of delay
                                g_link_list[link_seq_no].m_LinkCumulativeVirtualDelay[relative_t] = g_link_list[link_seq_no].m_CumulativeVirtualDelayCount
                                active_agent.m_Veh_LinkDepartureTime_in_simu_interval[active_agent.m_current_link_seq_no] = t + 1
                        #condition2: dropoff link
                        if(g_link_list[link_seq_no].service_type == -1):
                            next_vehicle_link_seq_no = active_agent.path_link_seq_no_list[active_agent.m_current_link_seq_no + 1]
                            active_agent.m_Veh_LinkArrivalTime_in_simu_interval[active_agent.m_current_link_seq_no + 1] = t
                            active_agent.m_Veh_LinkDepartureTime_in_simu_interval[active_agent.m_current_link_seq_no + 1] = t + \
                                g_link_list[next_vehicle_link_seq_no].free_flow_travel_time_in_simu_interval
                            for agent_id in active_agent.m_PassengerList:
                                agent_seq_no = g_map_agent_id_to_agent_seq_no[agent_id]
                                current_link_seq = g_agent_list[agent_seq_no].m_current_link_seq_no
                                next_pax_link_seq_no = current_link_seq + 1
                                g_agent_list[agent_seq_no].m_Veh_LinkDepartureTime_in_simu_interval[current_link_seq] = t
                                g_agent_list[agent_seq_no].m_Veh_LinkArrivalTime_in_simu_interval[current_link_seq + 1] = t
                                g_agent_list[agent_seq_no].m_Veh_LinkDepartureTime_in_simu_interval[current_link_seq + 1] = \
                                            active_agent.m_Veh_LinkDepartureTime_in_simu_interval[active_agent.m_current_link_seq_no + 1]
                                g_agent_list[agent_seq_no].m_current_link_seq_no+=1
                                
                                if(next_pax_link_seq_no != next_vehicle_link_seq_no):  # their next link seq no is not consistent
                                    active_agent.m_PassengerList.remove(agent_id)
                                    g_agent_list[agent_seq_no].m_bMoveable = True
                            #move
                            active_agent.m_current_link_seq_no+=1
                            g_link_list[link_seq_no].m_LinkOutFlowCapacity[relative_t] -= 1
                       #condition 3: pick up link
                        if(g_link_list[link_seq_no].service_type == 1):
                            if(active_agent.agent_id=='4'):
                                a=1
                            if(active_agent.agent_service_type == 1 and active_agent.m_bMoveable == True): #passenger: can enter link waiting queue
                                g_link_list[link_seq_no].m_waiting_traveler_queue.append(active_agent.agent_id)
                            if(active_agent.agent_service_type == 2): #vehicle: can pick up travelers from the link waiting queue
                                remaining_seat_capacity = active_agent.GetRemainingCapacity()
                                for p_ready in range(remaining_seat_capacity):
                                    if(len(g_link_list[link_seq_no].m_waiting_traveler_queue) == 0):
                                        break
                                    p = g_link_list[link_seq_no].m_waiting_traveler_queue[0]
                                    g_link_list[link_seq_no].m_waiting_traveler_queue.pop(0)
                                    active_agent.Pickup(p)

                                next_vehicle_link_seq_no = active_agent.path_link_seq_no_list[active_agent.m_current_link_seq_no + 1]
                                active_agent.m_Veh_LinkArrivalTime_in_simu_interval[active_agent.m_current_link_seq_no + 1] = t
                                active_agent.m_Veh_LinkDepartureTime_in_simu_interval[active_agent.m_current_link_seq_no + 1] = t + \
                                    g_link_list[next_vehicle_link_seq_no].free_flow_travel_time_in_simu_interval
                                for agent_id in active_agent.m_PassengerList:
                                    #move all passenger's timetable to the next
                                    #link along the vehicle trajectory
                                    agent_seq_no = g_map_agent_id_to_agent_seq_no[agent_id]
                                    current_link_seq = g_agent_list[agent_seq_no].m_current_link_seq_no
                                    next_pax_link_seq_no = current_link_seq + 1
                                    g_agent_list[agent_seq_no].m_Veh_LinkDepartureTime_in_simu_interval[current_link_seq] = t
                                    g_agent_list[agent_seq_no].m_Veh_LinkArrivalTime_in_simu_interval[current_link_seq + 1] = t
                                    g_agent_list[agent_seq_no].m_Veh_LinkDepartureTime_in_simu_interval[current_link_seq + 1] = \
                                            active_agent.m_Veh_LinkDepartureTime_in_simu_interval[active_agent.m_current_link_seq_no + 1]
                                    g_agent_list[agent_seq_no].m_current_link_seq_no+=1
                                    g_agent_list[agent_seq_no].m_bMoveable = False
                                #move
                                active_agent.m_current_link_seq_no += 1
                                g_link_list[link_seq_no].m_LinkOutFlowCapacity[relative_t] -= 1;
                                        
            for active_agent in g_active_agent_list[::-1]:  #dele the completed agent
                if(active_agent.bActive==False):
                    g_active_agent_list.remove(active_agent)
def g_OutputFiles():
    with open ('output_LinkMOE.csv','w',newline='') as fp:
        writer=csv.writer(fp)
        line=["from_node_id","to_node_id","cumulative_arrival_count","cumulative_departure_count",
              "travel_time_in_min"]
        writer.writerow(line)
        for t in range( int(g_Simulation_StartTimeInMin),int(g_Simulation_EndTimeInMin),1):
            time_in_simu_interval = int(t * 60 / g_number_of_seconds_per_interval)
            relative_time_interval=g_A2R_simu_interval(time_in_simu_interval)
            for link in g_link_list:
                line=[link.extern_from_node,link.extern_to_node,link.m_LinkCumulativeArrival[relative_time_interval],
                      link.m_LinkCumulativeDeparture[relative_time_interval],link.travel_time]
                writer.writerow(line)
    with open("output_agent.csv","w",newline='') as fp:
        writer=csv.writer(fp)
        line=["agent_id","agent_service_type","origin_node_id","destination_node_id","cost","departure_time_in_min",
              "path_node_sequence","path_time_sequence"]
        writer.writerow(line)
        for agent in g_agent_list:
            path_time_list=list()
            cost=(agent.m_Veh_LinkDepartureTime_in_simu_interval[-1]-agent.m_Veh_LinkArrivalTime_in_simu_interval[0])\
                * g_number_of_seconds_per_interval / 60.0
            for i in range(0,len(agent.path_link_seq_no_list),1):
                TA_in_min=agent.m_Veh_LinkArrivalTime_in_simu_interval[i]* g_number_of_seconds_per_interval / 60.0
                TD_in_min=agent.m_Veh_LinkDepartureTime_in_simu_interval[i]* g_number_of_seconds_per_interval / 60.0
                if(i==0):
                    path_time_list.extend([str(TA_in_min),str(TD_in_min)])
                else:
                    path_time_list.append(str(TD_in_min))
                path_time_sequence=';'.join(path_time_list)
            
            line=[agent.agent_id,agent.agent_service_type,agent.from_origin_node_id,
                  agent.to_destination_node_id,cost,agent.departure_time_in_min,
                  agent.path_node_str,path_time_sequence]
            writer.writerow(line)
    with open("output_solution.csv",'w',newline='') as fp:
        writer=csv.writer(fp)
        line=["number_of_nodes","number_of_links","number_of_agents","CPU running time"]
        writer.writerow(line)
        line_value=[g_number_of_nodes,g_number_of_links,g_number_of_agents,end_time-begin_time]
        writer.writerow(line_value)
            
if __name__=="__main__":

    begin_time=time()
    g_ReadInputData()

    g_TrafficSimulation()

    end_time=time()
    g_OutputFiles()
    

    
