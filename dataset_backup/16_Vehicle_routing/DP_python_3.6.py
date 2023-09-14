# coding=gb18030
'''
VRPPDTW in python
'''
import math
import csv
import copy

_MAX_LABEL_COST = 10000
g_node_list = [0]
g_agent_list = []
g_link_list = []
g_vehicle_list = []
g_passenger_list = []

g_internal_node_id_dict = {}
g_external_node_id_dict = {}
g_external_vehicle_id_dict = {}
g_internal_agent_no_dict ={}
g_external_passenger_id_dict={}

g_number_of_nodes = 0
g_number_of_links = 0
g_number_of_agents = 0
g_number_of_passengers=0
g_number_of_vehicles=0
g_number_of_time_intervals =50
internal_node_seq_no=1
g_passenger_vehicle_visit_flag =[]
g_passenger_number_of_visits=[0,0]


class Node:
    def __init__(self):
        self.node_id = 0
        self.x = 0.0
        self.y = 0.0
        self.outbound_node_list = []
        self.outbound_node_size=0
        self.node_type=0
        self.outbound_link_list = []
        self.g_node_passenger_id = 0
        self.g_activity_node_beginning_time = 0
        self.g_activity_node_ending_time=0


class Link:
    def __init__(self):
        self.length = 0.0
        self.link_id = 0
        self.from_node_seq_no = 0
        self.to_node_seq_no = 0
        self.direction=0
        self.number_of_lanes = 0
        self.lane_cap = 0
        self.speed_limit = 0.0
        self.jam_density= 0
        self.link_type=0
        self.free_flow_travel_time = 1.0
        self.node_type = 0



class Agent:
    def __init__(self):
        self.agent_id = 0
        self.from_node_id = 0
        self.to_node_id = 0
        self.departure_time_beginning = 0.0
        self.departure_time_ending = 0.0
        self.arrival_time_beginning = 0.0
        self.arrival_time_ending = 0.0
        self.capacity = 0
        self.base_profit = 0.0
        self.passenger_number_of_visits=0



def add_new_node_and_link_for_agent(agent_id,from_or_to_node_id,beginning_time,end_time,agent_type,\
                                    pickup_delivery_flag,origin_destination_flag):
    global g_number_of_nodes
    global g_number_of_links
    global internal_node_seq_no

    new_node = Node()
    new_link1 = Link()
    new_link2 = Link()
    if agent_type == 0:
        g_internal_node_id_dict[agent_id * 100] = internal_node_seq_no  # internal node seq
        g_external_node_id_dict[internal_node_seq_no] = agent_id * 100
        new_node.node_type = pickup_delivery_flag   #1 for pickup ,2 for delivery
        new_node.g_node_passenger_id = agent_id

        new_link1.link_type = agent_id*(2-pickup_delivery_flag)
        new_link2.link_type = agent_id*(2-pickup_delivery_flag)
        g_number_of_nodes+=1
        new_node.node_id = internal_node_seq_no
        new_node.outbound_node_size = 1
        g_node_list.append(new_node)
        new_node.outbound_node_list.append(from_or_to_node_id)
        g_node_list[from_or_to_node_id].outbound_node_list.append(g_number_of_nodes)
        g_node_list[from_or_to_node_id].outbound_node_size += 1
        new_node.g_activity_node_beginning_time = beginning_time
        new_node.g_activity_node_ending_time = end_time

        g_number_of_links += 1
        new_link1.link_id = g_number_of_links
        new_link1.from_node_seq_no = from_or_to_node_id
        new_link1.to_node_seq_no = internal_node_seq_no
        new_link1.direction = 1
        new_link1.number_of_lanes = 1
        new_link1.lane_cap = 0
        new_link1.speed_limit = 0.0
        new_link1.jam_density = 0

        new_link1.free_flow_travel_time = 1.0
        g_link_list.append(new_link1)
        g_node_list[from_or_to_node_id].outbound_link_list.append(new_link1)



        g_number_of_links += 1

        new_link2.link_id = g_number_of_links
        new_link2.from_node_seq_no = internal_node_seq_no
        new_link2.to_node_seq_no = from_or_to_node_id
        new_link2.direction = 1
        new_link2.number_of_lanes = 1
        new_link2.lane_cap = 0
        new_link2.speed_limit = 0.0
        new_link2.jam_density = 0

        new_link2.free_flow_travel_time = 1.0
        g_link_list.append(new_link2)
        g_node_list[internal_node_seq_no].outbound_link_list.append(new_link2)

    elif agent_type == 1:
        g_internal_node_id_dict[agent_id* 1000] = internal_node_seq_no
        g_external_node_id_dict[internal_node_seq_no] = agent_id * 1000
        new_link1.link_type = agent_id + 100
        new_link2.link_type = -(agent_id + 100)

        g_number_of_nodes+=1
        new_node.node_id = internal_node_seq_no
        new_node.node_type =3
        g_node_list.append(new_node)
        if origin_destination_flag==1:  #origin
            new_node.outbound_node_size = 1
            new_node.outbound_node_list.append(from_or_to_node_id)
            new_node.g_activity_node_beginning_time = beginning_time
            new_node.g_activity_node_ending_time = end_time

            g_number_of_links += 1
            new_link2.link_id = g_number_of_links
            new_link2.from_node_seq_no = internal_node_seq_no
            new_link2.to_node_seq_no = from_or_to_node_id
            new_link2.direction = 1
            new_link2.number_of_lanes = 1
            new_link2.lane_cap = 0
            new_link2.speed_limit = 0.0
            new_link2.jam_density = 0

            new_link2.free_flow_travel_time = 1.0
            g_link_list.append(new_link2)
            g_node_list[internal_node_seq_no].outbound_link_list.append(new_link2)




        if origin_destination_flag==2:
            g_node_list[from_or_to_node_id].outbound_node_list.append(g_number_of_nodes)
            g_node_list[from_or_to_node_id].outbound_node_size+=1
            new_node.g_activity_node_beginning_time = beginning_time
            new_node.g_activity_node_ending_time = end_time


            g_number_of_links += 1
            new_link1.link_id = g_number_of_links
            new_link1.from_node_seq_no = from_or_to_node_id
            new_link1.to_node_seq_no = internal_node_seq_no
            new_link1.direction = 1
            new_link1.number_of_lanes = 1
            new_link1.lane_cap = 0
            new_link1.speed_limit = 0.0
            new_link1.jam_density = 0

            new_link1.free_flow_travel_time = 1.0
            g_link_list.append(new_link1)
            g_node_list[from_or_to_node_id].outbound_link_list.append(new_link1)

    internal_node_seq_no += 1



def g_ReadInputData():
    #initialization
    global g_number_of_agents
    global g_number_of_vehicles
    global g_number_of_passengers
    global g_number_of_nodes
    global g_number_of_links
    global internal_node_seq_no
    global g_passenger_vehicle_visit_flag
    # read nodes information
    with open('node.csv', 'r') as fp:
        lines = fp.readlines()
        for l in lines[1:]:
            l = l.strip().split(',')
            if int(l[3]) ==1:
                continue
            try:
                node = Node()
                g_internal_node_id_dict[int(l[1])] = internal_node_seq_no  # internal node seq
                g_external_node_id_dict[internal_node_seq_no] = int(l[1])  # internal node seq--external_node_seq
                node.node_id = internal_node_seq_no #node id means internal node seq
                internal_node_seq_no += 1
                node.x = float(l[5])
                node.y = float(l[6])
                g_node_list.append(node)
                g_number_of_nodes += 1

                if g_number_of_nodes % 1000 == 0:
                    print('reading {} nodes..' \
                          .format(g_number_of_nodes))
            except:
                print('Bad read. Check file your self')
        print('nodes_number:{}'.format(g_number_of_nodes))

    with open('link.csv', 'r') as fl:
        linel = fl.readlines()
        for l in linel[1:]:
            l = l.strip().split(',')
            if int(l[9]) !=1:
                continue
            link = Link()
            link.link_id = g_number_of_links
            link.from_node_seq_no = g_internal_node_id_dict[int(l[2])]
            link.to_node_seq_no = g_internal_node_id_dict[int(l[3])]
            link.direction=1
            link.length = float(l[5])
            link.number_of_lanes = int(l[6])
            link.speed_limit = float(l[7])

            link.free_flow_travel_time = link.length / link.speed_limit * 60

            g_node_list[link.from_node_seq_no].outbound_node_list.append(link.to_node_seq_no)
            g_node_list[link.from_node_seq_no].outbound_node_size=len(g_node_list[link.from_node_seq_no].outbound_node_list)


            g_link_list.append(link)
            g_number_of_links += 1

            g_node_list[link.from_node_seq_no].outbound_link_list.append(link)

            if g_number_of_links % 8000 == 0:
                print('reading {} links..' \
                      .format(g_number_of_links))

        print('links_number:{}'.format(g_number_of_links))



    with open('input_agent.csv', 'r') as fa:
        linea = fa.readlines()
        for l in linea[1:]:
            try:
                l = l.strip().split(',')
                agent = Agent()
                agent.agent_id = int(l[0])
                agent.agent_type=int(l[1])
                agent.from_node_id = g_internal_node_id_dict[int(l[2])]
                agent.to_node_id = g_internal_node_id_dict[int(l[3])]
                agent.departure_time_beginning = int(l[4])
                agent.departure_time_ending = int(l[4])+int(l[5])
                agent.arrival_time_beginning = int(l[6])
                agent.arrival_time_ending = int(l[6])+int(l[7])
                g_agent_list.append(agent)

                if agent.agent_type==1:
                    g_number_of_vehicles += 1
                    vehicle_no=g_number_of_vehicles
                    g_internal_agent_no_dict[agent.agent_id] = vehicle_no
                    g_external_vehicle_id_dict[vehicle_no] = agent.agent_id
                    agent.capacity = int(l[8])

                    add_new_node_and_link_for_agent(agent.agent_id, agent.from_node_id,
                                                    agent.departure_time_beginning,agent.departure_time_ending,
                                                    agent.agent_type,1,1)
                    add_new_node_and_link_for_agent(agent.agent_id, agent.to_node_id,
                                                    agent.arrival_time_beginning,agent.arrival_time_ending,
                                                    agent.agent_type,2,2)
                    g_vehicle_list.append(agent)
                else:
                    g_number_of_passengers += 1
                    pax_no = g_number_of_passengers

                    g_internal_agent_no_dict[agent.agent_id] = pax_no
                    g_external_passenger_id_dict[pax_no] = agent.agent_id
                    agent.base_profit = float(l[9])

                    add_new_node_and_link_for_agent(agent.agent_id, agent.from_node_id,
                                                    agent.departure_time_beginning,agent.departure_time_ending,
                                                    agent.agent_type,1,1)
                    add_new_node_and_link_for_agent(agent.agent_id, agent.to_node_id,
                                                    agent.arrival_time_beginning,agent.arrival_time_ending,
                                                    agent.agent_type,2,2)
                    g_passenger_list.append(agent)
            except:
                print('Bad read. Check file your self')
            g_passenger_vehicle_visit_flag = [[0 for i in range(g_number_of_vehicles)]]*g_number_of_passengers

        print('agents_passenger:{}'.format(g_number_of_passengers))
        print('agents_vehicle:{}'.format(g_number_of_vehicles))

class CVSState:
    def __init__(self):
        self.current_node_id = 0
        self.passenger_service_state={i:0 for i in range(g_number_of_passengers)}
        self.passenger_service_time={}
        self.passenger_service_begin_time={}
        self.passenger_carrying_state={}
        self.m_visit_sequence=[]
        self.m_link_sequence =[]
        self.m_visit_time_sequence=[]
        self.m_vehicle_capacity = 0
        self.LabelCost = 9999     #with LR price
        self.PrimalLabelCost = 9999    #without LR price
        self.m_final_arrival_time =0

    def CVSState(self):
        self.m_final_arrival_time = 0
        self.LabelCost = _MAX_LABEL_COST
        self.m_vehicle_capacity = 0

    def mycopy(self,pElement):
        self.current_node_id = copy.copy(pElement.current_node_id)
        self.passenger_service_state = {}
        self.passenger_service_state = copy.copy(pElement.passenger_service_state)
        self.passenger_service_time = {}
        self.passenger_service_time = copy.copy(pElement.passenger_service_time)
        self.passenger_service_begin_time = {}
        self.passenger_service_begin_time = copy.copy(pElement.passenger_service_begin_time)
        self.passenger_carrying_state = {}
        self.passenger_carrying_state = copy.copy(pElement.passenger_carrying_state)
        self.m_visit_sequence = []
        self.m_visit_sequence = copy.copy(pElement.m_visit_sequence)
        self.m_link_sequence=[]
        self.m_link_sequence= copy.copy(pElement.m_link_sequence)
        self.m_visit_time_sequence = []
        self.m_visit_time_sequence = copy.copy(pElement.m_visit_time_sequence)
        self.LabelCost = copy.copy(pElement.LabelCost)


    def GetPassengerServiceState(self,passenger_id):
        if self.passenger_service_state[passenger_id] == 1 or 2:
            return self.passenger_service_state[passenger_id]
        else:
            return 0
    def StartCarryingService(self,passenger_id, service_time):
        self.passenger_carrying_state[passenger_id] = 1
        self.m_vehicle_capacity += 1
        self.passenger_service_begin_time[passenger_id] = service_time

    def CompleteCarryingService(self,passenger_id, service_time):
        if self.passenger_carrying_state[passenger_id] == 1:
            self.passenger_carrying_state[passenger_id] == 2


    def MarkCarryingService(self,passenger_id,node_type,service_time):
        if node_type == 1:
            self.StartCarryingService(self, passenger_id, service_time)
        if node_type == 2:
            self.CompleteCarryingService(self,passenger_id, service_time)



    def IsAllServiceComplete(self):
        pass
    def CalculateLabelCost(self,vehicle_id):
        self.LabelCost = 0
        self.PrimalLabelCost = 0
        for i in range(g_number_of_passengers):
            if  self.passenger_service_state[i]==2:
                self.LabelCost  = self.LabelCost - g_passenger_list[i].base_profit


        self.LabelCost = self.LabelCost + self.m_visit_time_sequence[-1] - g_vehicle_list[vehicle_id].departure_time_beginning
        self.PrimalLabelCost = self.PrimalLabelCost + self.m_visit_time_sequence[-1] - g_vehicle_list[vehicle_id].departure_time_beginning
        for i in range(1,len(self.m_visit_sequence)):
            if self.m_visit_sequence[i]==self.m_visit_sequence[i-1]:
                self.LabelCost -=(1 - 0.5)*(self.m_visit_time_sequence[i] - self.m_visit_time_sequence[i - 1])
                self.PrimalLabelCost -= (1 - 0.5) * (self.m_visit_time_sequence[i] - self.m_visit_time_sequence[i - 1])




    # class CVSState  record pax and vehicle's completed service
    def CountPassengerNumberOfVisits(self,vehicle_id):
        for i in range(g_number_of_passengers):
            if self.passenger_service_state[i] == 2:
                g_passenger_number_of_visits[i] += 1
                g_passenger_vehicle_visit_flag[i][vehicle_id] = 1

    def generate_string_key(self):
        str ='n'
        str = str + "%d"%(self.current_node_id)
        for i in range(g_number_of_passengers):
            if self.passenger_service_state[i] == 1 or self.passenger_service_state[i] == 2:
                str=str+ "_"+"%d"%(i)+"["+"%d"%(self.passenger_service_state[i])+"]"
        return str




class C_time_indexed_state_vector:
    def __init__(self):
        self.current_time=0
        self.m_VSStateVector=[]
        self.m_state_map={}
    def Reset(self):
        self.current_time = 0
        self.m_VSStateVector=[]
        self.m_state_map={}
    def m_find_state_index(self,string_key):
        if string_key in self.m_state_map.values():
            return list(self.m_state_map.keys())[list(self.m_state_map.values()).index(string_key)]


        else:
            return -1

    def update_state(self,new_element):

        string_key = new_element.generate_string_key()
        state_index = self.m_find_state_index(string_key)

        if state_index == -1:
            state_index = len(self.m_VSStateVector)
            self.m_VSStateVector.append(new_element)
            self.m_state_map[state_index] = string_key
        else:

            if new_element.LabelCost < self.m_VSStateVector[state_index].LabelCost:
                self.m_VSStateVector[state_index]=new_element


    def Sort(self):
        
        #self.m_VSStateVector=sorted(self.m_VSStateVector, cmp=lambda x, y: cmp(x.LabelCost, y.LabelCost))
        self.m_VSStateVector=sorted(self.m_VSStateVector,key=lambda x:x.LabelCost)
        #self.m_state_map = {}


    def SortAndCleanEndingState(self,BestKValue):
        if len(self.m_VSStateVector) > 2 * BestKValue:
            #self.m_VSStateVector[0].sort()
            #self.m_VSStateVector=sorted(self.m_VSStateVector, cmp=lambda x, y: cmp(x.LabelCost, y.LabelCost))
            self.m_VSStateVector=sorted(self.m_VSStateVector,key=lambda x:x.LabelCost)
            
            #self.m_state_map = {}
            self.m_VSStateVector=self.m_VSStateVector[0,BestKValue]

    def GetBestValue(self,DualPriceFlag,vehicle_id):
        if len(self.m_VSStateVector) >= 1:
            state_str = self.m_VSStateVector[-1].generate_string_key()
            self.m_VSStateVector[-1].CountPassengerNumberOfVisits(vehicle_id)
            if DualPriceFlag == 1:
                return self.m_VSStateVector[-1].LabelCost
            else:
                return self.m_VSStateVector[-1].PrimalLabelCost



def g_optimal_time_dependenet_dynamic_programming(
        vehicle_id,
        origin_node,
        departure_time_beginning,
        departure_time_ending,
        destination_node,
        arrival_time_beginning,
        arrival_time_ending,
        vehicle_capacity,
        BestKSize,
        DualCostFlag):

    global g_passenger_vehicle_visit_flag
    global g_time_dependent_state_vector
    global g_vehicle_passenger_visit_flag
    global g_vehicle_passenger_visit_allowed_flag
    global g_passenger_number_of_visits

    g_passenger_vehicle_visit_flag = [[0 for i in range(g_number_of_vehicles)]] * g_number_of_passengers
    g_passenger_vehicle_visit_allowed_flag = [[1 for i in range(g_number_of_vehicles)]] * g_number_of_passengers

    g_time_dependent_state_vector = [[None]*(arrival_time_ending-departure_time_beginning+2)]*g_number_of_vehicles

    if arrival_time_ending > g_number_of_time_intervals or g_node_list[origin_node].outbound_node_size == 0:

        return _MAX_LABEL_COST

    for p in range(g_number_of_passengers):
        g_passenger_vehicle_visit_flag[p][vehicle_id] = 0
     #step 2: Initialization  for origin node at the preferred departure time

    for t in range(departure_time_beginning,arrival_time_ending+1):

        g_time_dependent_state_vector[vehicle_id][t] = C_time_indexed_state_vector()
        g_time_dependent_state_vector[vehicle_id][t].Reset()
        g_time_dependent_state_vector[vehicle_id][t].current_time=t

    g_ending_state_vector=[None]*g_number_of_vehicles
    g_ending_state_vector[vehicle_id]= C_time_indexed_state_vector()
    g_ending_state_vector[vehicle_id].Reset()
    #origin_node
    element=CVSState()
    element.current_node_id = origin_node

    g_time_dependent_state_vector[vehicle_id][departure_time_beginning].update_state(element)

    # step 3:dynamic programming

    #1 sort m_VSStateVector by labelCost for scan best k elements in step2
    for t in range(departure_time_beginning,arrival_time_ending):
        g_time_dependent_state_vector[vehicle_id][t].Sort()
        #2 scan the best k elements
        for w_index in range(min(BestKSize,len(g_time_dependent_state_vector[vehicle_id][t].m_VSStateVector))):
            pElement=g_time_dependent_state_vector[vehicle_id][t].m_VSStateVector[w_index]
            from_node_id=pElement.current_node_id
            # step 2.1 link from_node to to_node
            from_node=g_node_list[from_node_id]

            for i in range(from_node.outbound_node_size):
                to_node_id=from_node.outbound_node_list[i]
                to_node = g_node_list[to_node_id]
                to_node_passenger_id = to_node.g_node_passenger_id
                to_node_type = to_node.node_type
                link_no = from_node.outbound_link_list[i]
                next_time = int(round(t + link_no.free_flow_travel_time ))


                #step 2.2 check feasibility of node type with the current element
                if next_time<= arrival_time_ending:
                    # 下一个点 activity_node
                    if to_node_passenger_id >= 1 and \
                        g_passenger_vehicle_visit_allowed_flag[to_node_passenger_id-1][vehicle_id] == 1:
                        # address the passengers' state transitions
                        # skip scanning when the origin/destination nodes arrival time is out of time window
                        if next_time > to_node.g_activity_node_ending_time:
                            continue
                        # feasible state transitions
                        if (to_node_type == 1 and pElement.GetPassengerServiceState(to_node_passenger_id-1) == 0)\
                               or (to_node_type == 2 and pElement.GetPassengerServiceState(to_node_passenger_id-1) == 1):
                            # waiting
                            # pickup process
                            if to_node_type == 1:
                                #skip pickup when the vehicle if on its capacity
                                if g_time_dependent_state_vector[vehicle_id][t].m_VSStateVector[w_index].m_vehicle_capacity >= vehicle_capacity:
                                    continue
                                #waiting
                                if next_time <to_node.g_activity_node_beginning_time:
                                    new_element = CVSState()
                                    new_element.mycopy(pElement)
                                    new_element.current_node_id = to_node_id
                                    new_element.passenger_service_state[to_node_passenger_id-1] = to_node_type
                                    #for arriving at activity node and begin wait
                                    new_element.m_visit_time_sequence.append(next_time)
                                    new_element.m_visit_sequence.append(to_node_id)

                                    #for wait until activity node's depature time
                                    new_element.m_visit_time_sequence.append(to_node.g_activity_node_beginning_time)
                                    new_element.m_visit_sequence.append(to_node_id)
                                    new_element.CalculateLabelCost(vehicle_id)

                                    g_time_dependent_state_vector[vehicle_id][to_node.g_activity_node_beginning_time].update_state(new_element)
                                    continue
                            # delivery process
                            if to_node_type == 2:
                                # waiting
                                if next_time < to_node.g_activity_node_beginning_time:
                                    new_element = CVSState()
                                    new_element.mycopy(pElement)
                                    new_element.current_node_id = to_node_id
                                    new_element.passenger_service_state[to_node_passenger_id-1] = to_node_type
                                    #for arriving at activity node and begin wait
                                    new_element.m_visit_time_sequence.append(next_time)
                                    new_element.m_visit_sequence.append(to_node_id)

                                    #for wait until activity node's depature time
                                    new_element.m_visit_time_sequence.append(to_node.g_activity_node_beginning_time)
                                    new_element.m_visit_sequence.append(to_node_id)
                                    new_element.CalculateLabelCost(vehicle_id)

                                    g_time_dependent_state_vector[vehicle_id][to_node.g_activity_node_beginning_time].update_state(new_element)
                                    continue
                            # donot need waiting
                            new_element = CVSState()
                            new_element.mycopy(pElement)
                            new_element.current_node_id = to_node_id
                            new_element.passenger_service_state[to_node_passenger_id-1] = to_node_type
                            new_element.m_visit_time_sequence.append(next_time)
                            new_element.m_visit_sequence.append(to_node_id)
                            new_element.m_link_sequence.append(link_no.link_id)
                            new_element.CalculateLabelCost(vehicle_id)

                            g_time_dependent_state_vector[vehicle_id][next_time].update_state(new_element)
                            continue

                    elif to_node_passenger_id < 1 and to_node_type == 0:
                        new_element = CVSState()
                        new_element.mycopy(pElement)
                        new_element.current_node_id = to_node_id
                        new_element.m_visit_time_sequence.append(next_time)
                        new_element.m_visit_sequence.append(to_node_id)
                        new_element.m_link_sequence.append(link_no.link_id)
                        new_element.CalculateLabelCost(vehicle_id)
                        g_time_dependent_state_vector[vehicle_id][next_time].update_state(new_element)
                        continue

                    elif to_node_id == destination_node:
                        new_element = CVSState()
                        new_element.mycopy(pElement)
                        #wait
                        if next_time < arrival_time_beginning:
                            new_element.m_visit_time_sequence.append(next_time)
                            new_element.m_visit_sequence.append(to_node_id)

                            #g_time_dependent_state_vector[vehicle_id][next_time].update_state(new_element)
                            new_element.m_visit_time_sequence.append(arrival_time_beginning)
                            new_element.m_visit_sequence.append(to_node_id)
                            new_element.m_link_sequence.append(link_no.link_id)
                            #g_time_dependent_state_vector[vehicle_id][next_time].update_state(new_element)
                            new_element.CalculateLabelCost(vehicle_id)

                        else:
                            new_element.m_visit_time_sequence.append(next_time)
                            new_element.m_visit_sequence.append(to_node_id)
                            new_element.m_link_sequence.append(link_no.link_id)
                            new_element.CalculateLabelCost(vehicle_id)
                            #g_time_dependent_state_vector[vehicle_id][next_time].update_state(new_element)

                        g_ending_state_vector[vehicle_id].update_state(new_element)
                        new_element.CalculateLabelCost(vehicle_id)
                        #print g_ending_state_vector[vehicle_id]
                        continue

    g_ending_state_vector[vehicle_id].Sort()
    print(g_ending_state_vector[vehicle_id].m_VSStateVector[0].m_visit_sequence)
    print(g_ending_state_vector[vehicle_id].m_VSStateVector[0].m_visit_time_sequence)
    print(g_ending_state_vector[vehicle_id].m_VSStateVector[0].m_link_sequence)
    return g_ending_state_vector[vehicle_id].m_VSStateVector[0].m_visit_sequence,g_ending_state_vector[vehicle_id].m_VSStateVector[0].m_visit_time_sequence,\
            g_ending_state_vector[vehicle_id].m_VSStateVector[0].m_link_sequence





if __name__=='__main__':
    print('Reading data......')
    g_ReadInputData()



    node_seq,time_seq,link_seq = g_optimal_time_dependenet_dynamic_programming(0, 11, 1, 3, 12, 14, 30, 3, 100, 1)
    for i in range(len(node_seq)):
        node_id = node_seq[i]
        if g_node_list[node_id].node_type == 1: #pickup
            passenger_id =  g_node_list[node_id].g_node_passenger_id
            g_passenger_list[passenger_id-1].start_index = i
        if g_node_list[node_id].node_type == 2: #dropoff
            passenger_id =  g_node_list[node_id].g_node_passenger_id
            g_passenger_list[passenger_id-1].end_index = i+1
    for p in g_passenger_list:
        p.node_sequence = node_seq[p.start_index:p.end_index]
        p.link_sequence = link_seq[(p.start_index+1):(p.end_index+1)]
        p.time_sequence = time_seq[p.start_index:p.end_index]

    f = open("agent.csv", "w")
    f.write("agent_id,o_zone_id,d_zone_id,path_id,o_node_id,to_node_id,agent_type,demand_period,\
    travel_time,distance,node_sequence,link_sequence,time_sequence,\n")

    for p in g_passenger_list:
        travel_time = p.time_sequence[-1]-p.time_sequence[0]
        distance = 0
        for i in p.link_sequence:
            distance+=g_link_list[i-1].length
        f.write(str(p.agent_id)+","+str(p.from_node_id)+","+str(p.to_node_id)+",0,"+str(p.from_node_id)+","+str(p.to_node_id)+","\
                +str(p.agent_type)+",AM,"+str(travel_time)+","+str(distance)+",")
        str1 = ""
        str2 = ""
        str3 = ""
        for s in range(len(p.time_sequence)):
            str1 = str1 + str(p.node_sequence[s]) + ";"
            str2 = str2 + str(p.link_sequence[s]) + ";"
            str3 = str3 +"07"+  str(p.time_sequence[s]) + ":00"+  ";"
        f.write((str1) + "," + (str2) + ","+(str3)+"\n")

    # vehicle
    for v in g_vehicle_list:
        travel_time = time_seq[-1] - time_seq[0]
        distance = 0
        for i in link_seq:

            distance += g_link_list[i-1].length
        f.write(
            str(3) + "," + str(v.from_node_id) + "," + str(v.to_node_id) + ",1," + str(v.from_node_id) + "," + str(
                v.to_node_id) + "," \
            + str(v.agent_type) + ",-," + str(travel_time) + "," + str(distance) + ",")
        str1 = ""
        str2 = ""
        str3 = ""
        for s in range(len(time_seq)):
            str1 = str1 + str(node_seq[s]) + ";"
            str2 = str2 + str(link_seq[s]) + ";"
            str3 = str3 + "07"+str(time_seq[s])+ ":00"+ ";"
        f.write((str1) + "," + (str2) + "," + (str3) + "\n")