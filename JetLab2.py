# -*- coding: utf-8 -*-
"""
Created on Thu Apr 21 19:23:09 2022

@author: deihi
"""

# Code Written in collaboration with Cameron White (Course 6).


# Despite the existence of libaries like pyCycle for the specific needs of this project, I had to write the entirety of thi code. 
# Identifyiing the point where all partial derivatives are 0 would be the ideal option without constraints. The original idea of optimization would have been better served with a Lagrangian field.
# Over this field, we would have used the technique of Lagrangian multipiers. Unfortunately, very few packages are useful at this, espeially with this many variables. 
# We had to settle for a by-hand computation where we adjusted each  variable; first one by one, then two at a time, then 3, then 4, so on until we were changing everything at once. 
# We would then sto pthe iteration process as soon as a local/global maximum was reached. This would functionally be identical to using an epsilon/input variation as characterized below. 
# You only record a value when the next iteration/incarement is less.


# Phase 1:
# Initializing and base parameters as fixed, although some will be incremented to variables and setting up the turbojet process

# Stations 1 to 20: Inlet
P_ambient = 101325;                                                                             #Unchangeable; can NOT be varied.
T_ambient = 288.15;                                                                             #Unchangeable; can NOT be varied.
R_s = 287.1;                                                                                    #Unchangeable; can NOT be varied.
rho = P_ambient/(R_s * T_ambient);                                                              #Formula 
diameter_inlet = 1.27;                                                                          #VARIABLE INPUT
A_inlet = 3.1415 * (0.5 * diameter_inlet)**2;                                                   #Formula
C_ambient = 137.2;                                                                              #Unchangeable; can NOT be varied.
m_dot_1_to_2 = rho * A_inlet * C_ambient;                                                       #Formula

# Stations 2 to 23: 87% Efficient Fan
m_dot_FAN = m_dot_1_to_2;                                                                       #Formula
Gamma_FAN = 1.4;                                                                                #Unchangeable; can NOT be varied.                                    
c_p_FAN = 1004.5;                                                                               #Unchangeable; can NOT be varied.
c_v_FAN = c_p_FAN / Gamma_FAN;                                                                  #Formula
Pressure_Ratio_FAN = 1.5;                                                                       #VARIABLE INPUT
Temperature_Ratio_FAN = ((Pressure_Ratio_FAN)**(Gamma_FAN - 1)/(Gamma_FAN));                    #Formula 
T_2 = T_ambient;                                                                                #Formula
T_23 = (Temperature_Ratio_FAN) * (T_2);                                                         #Formula
Total_Work_FAN = m_dot_FAN * (c_p_FAN + c_v_FAN) * (T_23 - T_2)                                 #Formula
Efficiency_FAN = 0.87;                                                                          #VARIABLE INPUT
Actual_Work_FAN = (1/Efficiency_FAN) * (Total_Work_FAN);                                        #Formula

# Stations 23 to 25: 87% Efficient Compressor 1 (LOW Pressure)
Gamma_COMP_1 = 1.4;                                                                             #Unchangeable; can NOT be varied.
c_p_COMP_1 = 1004.5;                                                                            #Unchangeable; can NOT be varied.
Bypass_Ratio = 6.0;                                                                             #VARIABLE INPUT
m_dot_COMP_1 = 1/(Bypass_Ratio + 1) * (m_dot_FAN);                                              #Formula
Pressure_Ratio_COMP_1 = 2.0;                                                                    #VARIABLE INPUT
Temperature_Ratio_COMP_1 = (Pressure_Ratio_COMP_1)**((Gamma_COMP_1 - 1)/(Gamma_COMP_1));        #Formula
T_25 = (T_23) * (Temperature_Ratio_COMP_1);                                                     #Formula
Total_Work_COMP_1 = (m_dot_COMP_1) * (c_p_COMP_1) * (T_25 - T_23);                              #Formula
Efficiency_COMP_1 = 0.87;                                                                       #VARIABLE INPUT
Actual_Work_COMP_1 = (1/Efficiency_COMP_1) * (Total_Work_COMP_1);                               #Formula

# Stations 25 to 3; 83% Efficient Compressor 2 (HIGH Pressure)
Gamma_COMP_2 = 1.4;                                                                             #Unchangeable; can NOT be varied.
c_p_COMP_2 = 1004.5;                                                                            #Unchangeable; can NOT be varied.
m_dot_COMP_2 = m_dot_COMP_1;                                                                    #Formula
Pressure_Ratio_COMP_2 = 8.0;                                                                    #VARIABLE INPUT
Temperature_Ratio_COMP_2 = (Pressure_Ratio_COMP_2)**((Gamma_COMP_2 - 1)/Gamma_COMP_2);          #Formula
T_3 = (T_25) * (Temperature_Ratio_COMP_2);                                                      #Formula
Total_Work_COMP_2 = (m_dot_COMP_2) * (c_p_COMP_2) * (T_3 - T_25);                               #Formula
Efficiency_COMP_2 = 0.83;                                                                       #VARIABLE INPUT
Actual_Work_COMP_2 = (1/Efficiency_COMP_2) * (Total_Work_COMP_2);                               #Formula

# Stations 3 to Station 4: Ignition at Combutor
c_p_IGNI = 1148.5;                                                                              #Unchangeable; can NOT be varied.
Min_Lift_Valid = 160135.978;                                                                    #Unchangeable; can NOT be varied.
Min_Fuel_Rate = 0.62;
Fuel_Energy_Release = 43.15e6;
Combustion_Energy_Release = (Fuel_Energy_Release) * (Min_Fuel_Rate);
Total_Mass_Flowing = (Min_Fuel_Rate) + (m_dot_COMP_2);
# Heat_Added = (Total_Mass_Flowing) * (c_p_IGNI) * (T_4 - T_3)
T_4 = (Combustion_Energy_Release)/((Total_Mass_Flowing)*(c_p_IGNI)) + (T_3);                                        

# Stations 4 to 45: 87% Efficient Turbine 1 (HIGH Pressure)
Gamma_TURB_1 = 1.33;   
c_p_TURB_1 = 1148.5;
m_dot_TURB_1 = Total_Mass_Flowing;
Efficiency_TURB_1 = 0.87;
T_45 = T_4 - (Actual_Work_COMP_2) * (m_dot_TURB_1) / ( (Efficiency_TURB_1) * (c_p_TURB_1));
Temperature_Ratio_TURB_1 = T_45 / T_4;                                                          #Formula
Pressure_Ratio_TURB_1 = (Temperature_Ratio_TURB_1)**(Gamma_TURB_1/(Gamma_TURB_1 - 1))           #VARIABLE INPUT/FORMULA;

# Stations 45 to 5: 89% Efficient Turbine 2 (LOW Pressure)
Gamma_TURB_2 =  1.33;                                                                                           #Unchangeable; can NOT be varied.
c_p_TURB_2 = 1148.5;                                                                                            #Unchangeable; can NOT be varied.
m_dot_TURB_2 = Total_Mass_Flowing;                                                                              #Formula
Efficiency_TURB_2 = 0.89;                                                                                       #VARIABLE INPUT
T_5 = T_45 - (Actual_Work_COMP_1 + Actual_Work_FAN) * (m_dot_TURB_1) / ( (Efficiency_TURB_1) * (c_p_TURB_1) );  #Formula
Temperature_Ratio_TURB_2 = T_5 / T_45;                                                                          #Formula
Pressure_Ratio_TURB_2 = (Temperature_Ratio_TURB_2)**(Gamma_TURB_2/(Gamma_TURB_2 - 1))                           #VARIABLE INPUT/FORMULA;

# Stations 5 to 9: Outlet
m_dot_Exit = m_dot_TURB_2                                                                      #Unchangeable; can NOT be varied.
c_p_Exit = 1148.5;                                                                             #Unchangeable; can NOT be varied.
T_Stag_Exit = T_5;                                                                             #Formula
Exit_Velocity = (2 * (c_p_Exit)*(T_Stag_Exit - T_ambient))**0.5;                               #Formula
# Thrust Produced
Thrust = (m_dot_Exit) * (Exit_Velocity);                                                       #Formula 


# Phase 2:
# Optimizing with Constraints over Stations and Global Constraitns for Entire System to minimize SFC
# Things we can change:
# (01) Efficiencies (May vary by up to 0.02 i.e. 2% closer or farther from 100%.)
# (02) Bypass Ratio (Bound by sense i.e. a realistic size of an engine.)
# (03) Compression Ratio: LOW PRESSURE Compressor 1 (Can NOT exceed 3.0.)
# (04) Temperature after LOW PRESSURE Turbine 2 (T5 Can NOT exceed 1700F.)
# (05) Temperature after HIGH PRESSURE Compressor 2 (T3 Can NOT exceed 1200F.)
# (06) Inlet Size (Diameter can NOT exceed 59in.)
# (07) Fan Compression Ratio (Keep it realistically small as fans induce little compression.)
# (08) Compression Ratio: HIGH PRESSURE Compressor (A temp. too high may melt the steel/Al.) 
# (09) Compression Ratio: HIGH PRESSURE Turbine
# (10) Compression Ratio: LOW PRESSURE Turbine 
# This means there could computationally be:
# 10C1 + 10C2 + 10C3 + ... + 10C10 changes possible. This would obviously take far too long. 
# We can instead opt to change only what actually matters. 
# Namely, this would be ...

    
#while diameter_inlet <= 1.4986:
#    diameter_inlet += 0.01;
#while Efficiency_FAN <= 0.89:
#    Efficiency_FAN += 0.005;
#Efficiencies = {Efficiency_FAN, Efficiency_COMP_1, Efficiency_COMP_2, Efficiency_TURB_1, Efficiency_TURB_2}
#for i in Efficencies;:
#    while i <= i + 0.02:
#        i += 
 
# A more efficient method would be to set a thrust target and then find the SFC needed to meet that target. Obviously, the omre thrust, the more fuel needed.
# All we need to do is let Python see the variables that can change, and change those, within the limits of our sliders/constraints and then see the SFC it returns.   
# The process is evaluated below as follows: we simply create one super-massive function to describe thrust and make the inputs known, then each of these inputs will be varied one by one and collectively.
# Once the thrust target is met, the input values of relevance are noted. 
# We simply create the function, set the output to have a fixed value, then proceed to vary the inputs, but keep the inputs within the bounds of the design constraints. 
    
Inputs = [ 
          ('Pressure_Ratio_FAN', (1.1, 1.75, 0.05)),
          ('Efficiency_FAN', (0.87, 0.89, 0.01)),
          ('Bypass_Ratio', (4.0, 20.0, 1.0)),
          ('Pressure_Ratio_COMP_1', (2.0, 3.0, 0.1)),
          ('Efficiency_COMP_1', (0.87,0.89, 0.01)),
          ('Pressure_Ratio_COMP_2', (5.0,15.0, 0.5)),
          ('Efficiency_COMP_2', (0.83,0.85, 0.01)),
          ('Efficiency_TURB_1', (0.87,0.89, 0.01)),
          ('Efficiency_TURB_2', (0.89,0.91, 0.01))]
          #('T_5', (373.15,1200,25)),
          #('T_3', (373.15,922,25))]

#Turbines 1 and 2 had pressure ratios undefined due to being realated to compressors 1 and 2 and that would mean changing too many things at once, preventing equality.

def find_vals(Inputs, Goal_Thrust, Epsilon, vals):
    if not Inputs:
        Thrust = compute_thrust(vals)
        if abs(Thrust - Goal_Thrust) <= Epsilon and Thrust >= Goal_Thrust:
            return vals
        else:
            return None
        
    var, constraints = Inputs[0]
    
    #Iterate for all values
    val = constraints[0]
    
    new_inputs = [item for item in Inputs if item[0] != var]
    while val <= constraints[1]:
        new_vals = vals.copy()
        new_vals[var] = val
        _vals = find_vals(new_inputs, Goal_Thrust, Epsilon, new_vals)
        if _vals is not None:
            return _vals
        val += constraints[2]
    return None

def compute_thrust(vals):
    #To compute the thrust
    #diameter_inlet = vals['diameter_inlet']
    Pressure_Ratio_FAN = vals['Pressure_Ratio_FAN']
    Efficiency_FAN = vals['Efficiency_FAN']
    Bypass_Ratio = vals['Bypass_Ratio']
    Pressure_Ratio_COMP_1 = vals['Pressure_Ratio_COMP_1']
    Efficiency_COMP_1 = vals['Efficiency_COMP_1']
    Pressure_Ratio_COMP_2 = vals['Pressure_Ratio_COMP_2']
    Efficiency_COMP_2 = vals['Efficiency_COMP_2']
    Efficiency_TURB_1 = vals['Efficiency_TURB_1']
    Efficiency_TURB_2 = vals['Efficiency_TURB_2']
    #T_3 = vals['T_3']
    #T_5 = vals['T_5']                             

    # Ambient to Stations 1 and 2
    P_ambient = 101325                                                                              #Unchangeable; can NOT be varied.
    T_ambient = 288.15;                                                                             #Unchangeable; can NOT be varied.
    R_s = 287.1;                                                                                    #Unchangeable; can NOT be varied.
    rho = P_ambient/(R_s * T_ambient);                                                              #Formula 
    #print(rho)
    diameter_inlet = max(1.27,1.27);                                                                #VARIABLE INPUT;
    A_inlet = 3.1415 * (0.5 * diameter_inlet)**2;                                                   #Formula
    #print(A_inlet)
    C_ambient = 137.2;                                                                              #Unchangeable; can NOT be varied.
    m_dot_1_to_2 = rho * A_inlet * C_ambient;                                                       #Formula
    #print(m_dot_1_to_2)

    # Stations 2 to 23: 87% Efficient Fan
    m_dot_FAN = m_dot_1_to_2;                                                                       #Formula
    Gamma_FAN = 1.4;                                                                                #Unchangeable; can NOT be varied.                                    
    c_p_FAN = 1004.5;                                                                               #Unchangeable; can NOT be varied.
    c_v_FAN = c_p_FAN / Gamma_FAN;                                                                  #Formula
    #print(c_v_FAN)
    #Pressure_Ratio_FAN = 1.5;                                                                      #VARIABLE INPUT
    Temperature_Ratio_FAN = (Pressure_Ratio_FAN)**((Gamma_FAN - 1)/(Gamma_FAN));                    #Formula 
    #print(Temperature_Ratio_FAN)
    T_2 = T_ambient;                                                                                #Formula
    T_23 = (Temperature_Ratio_FAN) * (T_2);                                                         #Formula
    Total_Work_FAN = m_dot_FAN * (7*(c_p_FAN)*(T_23 - T_2) + (c_v_FAN)*(T_23-T_2))                  #Formula
    #print(Total_Work_FAN) 
    #Efficiency_FAN = 0.87;                                                                         #VARIABLE INPUT
    Actual_Work_FAN = (1/Efficiency_FAN) * (Total_Work_FAN);                                        #Formula
    #print(Actual_Work_FAN)
    #BE CAREFUL TO NOT CONFUSE ACTUAL AND TOTAL FOR TURBINE CONNECTION

    # Stations 23 to 25: 87% Efficient Compressor 1 (LOW Pressure)
    Gamma_COMP_1 = 1.4;                                                                             #Unchangeable; can NOT be varied.
    c_p_COMP_1 = 1004.5;                                                                            #Unchangeable; can NOT be varied.
    # REMEMBER TO COMMENT BACK THE BYPASS RATIO
    #Bypass_Ratio = 6.0;                                                                            #VARIABLE INPUT
    m_dot_COMP_1 = 1/(Bypass_Ratio + 1) * (m_dot_FAN);                                              #Formula
    #print(m_dot_FAN)
    #print(m_dot_COMP_1)   
    # ERROR m_dot_COMP_1 = m_dot_FAN #(fan bypass corrected above)
    # REMEMBER TO COMMENT BACK THE PRESSURE RATIO TO FIX TEMP RATIO
    #Pressure_Ratio_COMP_1 = 2.0;                                                                   #VARIABLE INPUT
    Temperature_Ratio_COMP_1 = (Pressure_Ratio_COMP_1)**((Gamma_COMP_1 - 1)/(Gamma_COMP_1));        #Formula
    T_25 = (T_23) * (Temperature_Ratio_COMP_1);                                                     #Formula
    Total_Work_COMP_1 = (m_dot_COMP_1) * (c_p_COMP_1) * (T_25 - T_23);                              #Formula
    #print(Total_Work_COMP_1)
    #Efficiency_COMP_1 = 0.87;                                                                      #VARIABLE INPUT
    Actual_Work_COMP_1 = (1/Efficiency_COMP_1) * (Total_Work_COMP_1);                               #Formula
    #print(Actual_Work_COMP_1)

    # Stations 25 to 3; 83% Efficient Compressor 2 (HIGH Pressure)
    Gamma_COMP_2 = 1.4;                                                                             #Unchangeable; can NOT be varied.
    c_p_COMP_2 = 1004.5;                                                                            #Unchangeable; can NOT be varied.
    m_dot_COMP_2 = m_dot_COMP_1;                                                                    #Formula
    #REMEMBER TO COMMENT BACK PRESSURE RATIO
    #Pressure_Ratio_COMP_2 = 8.0;                                                                   #VARIABLE INPUT
    Temperature_Ratio_COMP_2 = (Pressure_Ratio_COMP_2)**((Gamma_COMP_2 - 1)/(Gamma_COMP_2));        #Formula
    T_3 = (T_25) * (Temperature_Ratio_COMP_2);                                                      #Formula
    Total_Work_COMP_2 = (m_dot_COMP_2) * (c_p_COMP_2) * (T_3 - T_25);                               #Formula
    #REMEMBER TO COMMENT BACK EFFICIENCY
    #Efficiency_COMP_2 = 0.83;                                                                      #VARIABLE INPUT
    Actual_Work_COMP_2 = (1/Efficiency_COMP_2) * (Total_Work_COMP_2);                               #Formula
    #print(Actual_Work_COMP_2)
    T_3 = min(922,T_3);

    # Stations 3 to Station 4: Ignition at Combutor
    c_p_IGNI = 1148.5;                                                                              #Unchangeable; can NOT be varied.
    Min_Lift_Valid = 160135.978;                                                                    #Unchangeable; can NOT be varied.
    Min_Fuel_Rate = 0.62;
    Fuel_Energy_Release = 43.15e6;
    Combustion_Energy_Release = (Fuel_Energy_Release) * (Min_Fuel_Rate);
    #print(Combustion_Energy_Release)
    Total_Mass_Flowing = (Min_Fuel_Rate) + (m_dot_COMP_2);
    #print(Total_Mass_Flowing)
    # Heat_Added = (Total_Mass_Flowing) * (c_p_IGNI) * (T_4 - T_3)
    T_4 = (Combustion_Energy_Release)/((Total_Mass_Flowing)*(c_p_IGNI)) + (T_3); 
    T_4 =  min(1475,T_4);                                     

    # Stations 4 to 45: 87% Efficient Turbine 1 (HIGH Pressure)
    Gamma_TURB_1 = 1.33;   
    c_p_TURB_1 = 1148.5;
    m_dot_TURB_1 = Total_Mass_Flowing;
    #print(m_dot_TURB_1)
    #Efficiency_TURB_1 = 0.87;
    T_45 = T_4 - (Actual_Work_COMP_2) *  (1/ ( (Efficiency_TURB_1) * (c_p_TURB_1) * (m_dot_TURB_1) ));
    #print(T_45)
    Temperature_Ratio_TURB_1 = T_45 / T_4;
    #print(Temperature_Ratio_TURB_1)                                                              #Formula
    Pressure_Ratio_TURB_1 = (Temperature_Ratio_TURB_1)**((Gamma_TURB_1)/(Gamma_TURB_1 - 1))       #VARIABLE INPUT/FORMULA;
    #print(Pressure_Ratio_TURB_1)

    # Stations 45 to 5: 89% Efficient Turbine 2 (LOW Pressure)
    Gamma_TURB_2 =  1.33;                                                                         #Unchangeable; can NOT be varied.
    c_p_TURB_2 = 1148.5;                                                                          #Unchangeable; can NOT be varied.
    m_dot_TURB_2 = Total_Mass_Flowing;                                                            #Formula
    #Efficiency_TURB_2 = 0.89;                                                                    #VARIABLE INPUT
    #T_45 = 1500
    T_5 = T_45 - (Actual_Work_COMP_1 + Actual_Work_FAN) * 1/( ( (Efficiency_TURB_2) * (c_p_TURB_2) * (m_dot_TURB_1) ));                 #Formula
    # (Efficiency_TURB_2)*(c_p_TURB_2)*(T_45 - T_5)*(m_dot_TURB_2) =  (1/Efficiency_TURB_1)*(Actual_Work_COMP_1 + Actual_Work_FAN)
    #print(T_5)
    T_5 = min(1200, T_5);
    Temperature_Ratio_TURB_2 = T_5 / T_45;                                                        #Formula
    Pressure_Ratio_TURB_2 = (Temperature_Ratio_TURB_2)**(Gamma_TURB_2/(Gamma_TURB_2 - 1))         #VARIABLE INPUT/FORMULA;

    # Stations 5 to 9: Outlet
    m_dot_Exit = m_dot_TURB_2                                                                     #Unchangeable; can NOT be varied.
    c_p_Exit = 1148.5;                                                                            #Unchangeable; can NOT be varied.
    T_Stag_Exit = T_5;                                                                            #Formula
    Exit_Velocity = (2 * (c_p_Exit)*(T_Stag_Exit - T_ambient))**0.5;                              #Formula; Modded by AL as ((T_Stag_Exit/T_ambient)-1);#
    Exit_Velocity = ((((T_5/T_ambient) - 1) * (2/(Gamma_TURB_2 - 1)))**(0.5))*(330)
    #print(T_Stag_Exit)
    # Thrust Produced
    Thrust = (m_dot_Exit) * (Exit_Velocity); 
    #print(Exit_Velocity)
    #print(Thrust)
    #sfc = Min_Fuel_Rate.index(len(Min_Fuel_Rate - 1))/Thrust                                     #Formula
    sfc = Thrust/Min_Fuel_Rate;
    
    return Thrust
    return sfc

# Print best result(s)

vals = find_vals(Inputs, 160136/2, 1000, {})

# Optmization in sets of incrementing 1 thing;

