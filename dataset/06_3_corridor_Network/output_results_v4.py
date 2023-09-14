# @author       Jiawei Lu (jiaweil9@asu.edu)
# @time         2022/7/27 12:45
# @desc         [script description]

import pandas as pd


route_assignment_filepath = 'route_assignment.csv'


route_assignment_df = pd.read_csv(route_assignment_filepath)
route_assignment_df['od'] = route_assignment_df['o_zone_id'].astype(str) + '_' + route_assignment_df['d_zone_id'].astype(str)

ods = route_assignment_df['od'].unique()
agent_types = route_assignment_df['agent_type'].unique()


od_res_list = []

for od in ods:
    ra_od = route_assignment_df[route_assignment_df['od'] == od]

    if (ra_od['OD_impact_flag'] == 0).all():
        pass
    elif (ra_od['OD_impact_flag'] == 1).all():
        total_volume = ra_od['volume'].sum()
        od_res_dict = {'od': od, 'volume':total_volume}
        od_res_dict['average_distance'] = round((ra_od['volume'] * ra_od['distance']).sum() / total_volume, 2)
        od_res_dict['average_travel_time'] = round((ra_od['volume'] * ra_od['travel_time']).sum() / total_volume, 2)

        for agent_type in agent_types:
            ra_od_agent = ra_od[ra_od['agent_type'] == agent_type]
            volume = ra_od_agent['volume'].sum() if len(ra_od_agent) != 0 else 0
            od_res_dict[f'volume_{agent_type}'] = volume
            od_res_dict[f'average_distance_{agent_type}'] = round((ra_od_agent['volume'] * ra_od_agent['distance']).sum() / volume, 2) if len(ra_od_agent) != 0 else -1
            od_res_dict[f'average_travel_time_{agent_type}'] = round((ra_od_agent['volume'] * ra_od_agent['travel_time']).sum() / volume, 2) if len(ra_od_agent) != 0 else -1

        ra_od_nd = ra_od[ra_od['vehicle_diverted_flag'] == 0]
        volume = ra_od_nd['volume'].sum() if len(ra_od_nd) != 0 else 0
        od_res_dict['volume_nondiverted'] = volume
        od_res_dict[f'average_distance_nondiverted'] = round((ra_od_nd['volume'] * ra_od_nd['distance']).sum() / volume, 2) if len(ra_od_nd) != 0 else -1
        od_res_dict[f'average_travel_time_nondiverted'] = round((ra_od_nd['volume'] * ra_od_nd['travel_time']).sum() / volume, 2) if len(ra_od_nd) != 0 else -1

        ra_od_d = ra_od[ra_od['vehicle_diverted_flag'] == 1]
        volume = ra_od_d['volume'].sum() if len(ra_od_d) != 0 else 0
        od_res_dict['volume_diverted'] = volume
        od_res_dict[f'average_distance_diverted'] = round((ra_od_d['volume'] * ra_od_d['distance']).sum() / volume, 2) if len(ra_od_d) != 0 else -1
        od_res_dict[f'average_travel_time_diverted'] = round((ra_od_d['volume'] * ra_od_d['travel_time']).sum() / volume, 2) if len(ra_od_d) != 0 else -1

        od_res_list.append(od_res_dict)

    else:
        print(f'warning: OD_impact_flag of od {od} is not consistent')

if len(od_res_list) == 0:
    print(f'error: no impacted od')
else:
    od_res_df = pd.DataFrame(od_res_list)

    average_dict = {'od':'average'}
    average_dict['volume'] = od_res_df['volume'].mean()
    average_dict['average_distance'] = round((od_res_df['volume'] * od_res_df['average_distance']).sum() / od_res_df['volume'].sum(), 2)
    average_dict['average_travel_time'] = round((od_res_df['volume'] * od_res_df['average_travel_time']).sum() / od_res_df['volume'].sum(), 2)

    for agent_type in agent_types:
        average_dict[f'volume_{agent_type}'] = od_res_df[f'volume_{agent_type}'].mean()
        average_dict[f'average_distance_{agent_type}'] = round((od_res_df[f'volume_{agent_type}'] * od_res_df[f'average_distance_{agent_type}']).sum() / od_res_df[f'volume_{agent_type}'].sum(), 2)
        average_dict[f'average_travel_time_{agent_type}'] = round((od_res_df[f'volume_{agent_type}'] * od_res_df[f'average_travel_time_{agent_type}']).sum() / od_res_df[f'volume_{agent_type}'].sum(), 2)

    average_dict['volume_nondiverted'] = od_res_df['volume_nondiverted'].mean()
    average_dict['average_distance_nondiverted'] = round((od_res_df['volume_nondiverted'] * od_res_df['average_distance_nondiverted']).sum() / od_res_df['volume_nondiverted'].sum(), 2)
    average_dict['average_travel_time_nondiverted'] = round((od_res_df['volume_nondiverted'] * od_res_df['average_travel_time_nondiverted']).sum() / od_res_df['volume_nondiverted'].sum(), 2)
    average_dict['volume_diverted'] = od_res_df['volume_diverted'].mean()
    average_dict['average_distance_diverted'] = round((od_res_df['volume_diverted'] * od_res_df['average_distance_diverted']).sum() / od_res_df['volume_diverted'].sum(), 2)
    average_dict['average_travel_time_diverted'] = round((od_res_df['volume_diverted'] * od_res_df['average_travel_time_diverted']).sum() / od_res_df['volume_diverted'].sum(), 2)

    od_res_list.append(average_dict)
    od_res_df2 = pd.DataFrame(od_res_list)

    od_res_df2.to_csv('results.csv', index=False)