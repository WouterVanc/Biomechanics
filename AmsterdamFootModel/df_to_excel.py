# Write data to excel

import pandas as pd
import os

def indiv_data_to_excel_R(wd, dynamicangles_right, staticangles_right, trialname, c):
    df_static_temp = pd.DataFrame(data=staticangles_right)
    planelist = ['x', 'y', 'z']
    romlist = ['Min', 'Max']
    pl = pd.DataFrame(data=planelist, columns=['Plane'])
    rl = pd.DataFrame(data=romlist, columns=['value'])

    df_static = pl.join(df_static_temp)
    df_dynamic = pd.DataFrame(data=dynamicangles_right['gaitcycle'])
    df_dynamic_norm = pd.DataFrame(data=dynamicangles_right['gaitcycle_normalized'])
    df_dynamic_ROM_temp = pd.DataFrame([dynamicangles_right['gaitcycle_ROM']])
    df_dynamic_ROM = rl.join(df_dynamic_ROM_temp)

    with pd.ExcelWriter(os.path.join(wd, 'AFM_kinematics', 'Excel_data', (c + '_' + trialname[0:-4] + '_right.xlsx'))) as writer:
        df_static.to_excel(writer, sheet_name='Static Angles')
        df_dynamic.to_excel(writer, sheet_name='Dynamic Angles')
        df_dynamic_norm.to_excel(writer, sheet_name='Dynamic Angles Normalized')
        df_dynamic_ROM.to_excel(writer, sheet_name='Dynamic Angles_ROM')

def indiv_data_to_excel_L(wd, dynamicangles_left, staticangles_left, trialname, c):
    df_static_temp = pd.DataFrame(data=staticangles_left)
    planelist = ['x', 'y', 'z']
    romlist = ['Min', 'Max']
    pl = pd.DataFrame(data=planelist, columns=['Plane'])
    rl = pd.DataFrame(data=romlist, columns=['value'])

    df_static = pl.join(df_static_temp)
    df_dynamic = pd.DataFrame(data=dynamicangles_left['gaitcycle'])
    df_dynamic_norm = pd.DataFrame(data=dynamicangles_left['gaitcycle_normalized'])
    df_dynamic_ROM_temp = pd.DataFrame([dynamicangles_left['gaitcycle_ROM']])
    df_dynamic_ROM = rl.join(df_dynamic_ROM_temp)

    with pd.ExcelWriter(os.path.join(wd, 'AFM_kinematics', 'Excel_data', (c + '_' + trialname[0:-4] + '_left' '.xlsx'))) as writer:
        df_static.to_excel(writer, sheet_name='Static Angles')
        df_dynamic.to_excel(writer, sheet_name='Dynamic Angles')
        df_dynamic_norm.to_excel(writer, sheet_name='Dynamic Angles Normalized')
        df_dynamic_ROM.to_excel(writer, sheet_name='Dynamic Angles_ROM')

def avg_data_to_excel_R(avgdata_right, conditions, wd):
    with pd.ExcelWriter(os.path.join(wd, 'AFM_kinematics', 'Excel_data', 'Condition_avg_right.xlsx')) as writer:
        for condition in conditions:
            df_avg = pd.DataFrame(data=avgdata_right[condition + '_avg'])
            df_std = pd.DataFrame(data=avgdata_right[condition + '_std'])
            rom_avg = pd.DataFrame(data=avgdata_right[condition + '_avg_ROM'])
            temp = avgdata_right[condition + '_std_ROM']

            df_avg.to_excel(writer, sheet_name=(condition + '_avg'))
            df_std.to_excel(writer, sheet_name=(condition + '_std'))
            rom_avg.to_excel(writer, sheet_name=(condition + '_avg_ROM'))
            temp.to_excel(writer, sheet_name=(condition + '_std_ROM'))

def avg_data_to_excel_L(avgdata_left, conditions, wd):
    with pd.ExcelWriter(os.path.join(wd, 'AFM_kinematics', 'Excel_data', 'Condition_avg_left.xlsx')) as writer:
        for condition in conditions:
            df_avg = pd.DataFrame(data=avgdata_left[condition + '_avg'])
            df_std = pd.DataFrame(data=avgdata_left[condition + '_std'])
            rom_avg = pd.DataFrame(data=avgdata_left[condition + '_avg_ROM'])
            temp = avgdata_left[condition + '_std_ROM']

            df_avg.to_excel(writer, sheet_name=(condition + '_avg'))
            df_std.to_excel(writer, sheet_name=(condition + '_std'))
            rom_avg.to_excel(writer, sheet_name=(condition + '_avg_ROM'))
            temp.to_excel(writer, sheet_name=(condition + '_std_ROM'))