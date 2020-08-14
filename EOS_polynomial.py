# -*- coding: utf-8 -*-
"""
Created on Sun Feb 23 17:08:00 2020

@author: Dell
"""
###############################################################################
###############################################################################
######                                                               ##########
######                NAME : JAYANT JHARKHANDE                       ##########
######                 ROLL NUMBER : PE19M002                        ##########
######   CODE FOR EQUATION OF STATE SINGLE AND MULTIPLE COMPONENT    ##########
######                                                               ##########
###############################################################################
###############################################################################


import numpy as np




###############################################################################
##############  components and its moleculer weight         ###################
###############################################################################

########################### components of gas  ################################
component = [ 0.45, 0.05, 0.05, 0.03, 0.01, 0.01, 0.40] 
print("\n component of gas is    =  \n")
print( component ) 

#####################  moleculer weight of component ########################## 
Moleculer_weight = [16.043, 30.070, 44.097, 58.123, 72.150, 84.00, 215] 
print("\n moleculer weight of gas is    =  \n")
print( Moleculer_weight ) 

############ total moleculer weight calculation ########### 
individual_moleculer_weight_list = [component[i] * Moleculer_weight[i] for i in range(len(component))]    
print ("\nThe individual moleculer weight list is \n: " + str(individual_moleculer_weight_list)) 

########################################################
########  total moleculer weight loop      ############# 
########################################################  

def Total_moleculer_weight(individual_moleculer_weight_list) : 
      
#start value for total moleculer weight is zero
    result = 0
    for x in individual_moleculer_weight_list: 
         result = result + x  
    return result  

###########################################################
#### calling function for total moleculer weight     ######
###########################################################    
Ma = Total_moleculer_weight(individual_moleculer_weight_list)
print("\n")
print(Ma)


##########################################################
######  Calling a function from vow family of EOS ########
##########################################################



print('\n In petroleum industry we use Van der waals euqation of family for compresibility factor \n') 

print(' for van der waal model for single component value  = 1 ' + ' \n for redlich kwong model for single component value  = 2 '+ '\n for Soave-Redlich-Kwong model for single component value  = 3 \n' + '\n for Redlich-Kwong model for mixture component value  = 4  \n' +  'for Soave redlich kwong and redlich kwong mixture component value = 5 \n'   )  

EOS_model = input('\n EOS model you want :  ' )  

EOS_model = int(EOS_model)
#
#if EOS_model == 1
#           van_der_waal()
#  elif EOS_model == 2
#                  redlich_kwong()
#   else EOS_model == 3
#                  Soave_redlich_kwong()
    
def van_der_waal():
    print("\n")
    print("\n")
    print("\n")
    print(   "\t ########################################################## \t"  )
    print(   "\t ####### Van Der waal EoS for single component         ####### \t"  )
    print(   "\t ########################################################## \t"  )
    
    MUa = 0.421875
    MUb = 0.125
    #print (MUa)
    #print (MUb)
    ################################### constant ##############################
    R = 10.73
    ########################## Critical Temperature ###########################
    Tc = 666 
    print ("\n critical temperature value in rankine : "+ str(Tc))
    print(" \n\t --------------------------------------------------- \t" )
    ########################### critical pressure #############################
    Pc = 616.3
    print ("\n critical pressure value in psi : "+ str(Pc))
    print(" \n\t --------------------------------------------------- \t" )
    ##############################  component #################################
    M = 44
    print ("\n moleculer weight of the single component : "+ str(M))
    print(" \n\t --------------------------------------------------- \t" )
    #################### container Pressure and temperature####################
    ############################## pressure value in psi ######################
    #P = input(" \npressure in the container in psi :  ")
    P = 616.3
    P = int(P)
    print ("\n pressure value in psi : "+ str(P))
    print(" \n\t --------------------------------------------------- \t" )
    ###################### temprerature value in fahrenheit ###################
    T = 206 
    T = 460 + T
    print ("\n temperature value rankine: "+ str(T))
    print(" \n\t --------------------------------------------------- \t" )
    ########################## calculation of parameter a #####################
    a = (MUa)*((R**2)*(Tc**2)/Pc)
    print( "\nvalue of a is  = " + str(a) + " . ")
    print(" \n\t --------------------------------------------------- \t" )
    ####################### calculation of parameter b ########################
    b = (MUb)*((R)*(Tc)/Pc)
    print("\nvalue of b is  = " + str(b) + " . ")
    print(" \n\t --------------------------------------------------- \t" )
    ############################### coefficient A #############################
    A = ( (a)*(P) )/ ((R**2)*(T**2))
    print("\nvalue of coefficient  A  is  = " + str(A) + " . ")
    print(" \n\t --------------------------------------------------- \t" )
    ############################### coefficient B #############################
    B = ((b)*(P))/ ((R)*(T))
    print("\nvalue of coefficient B is  = " + str(B) + " . ")
    print(" \n\t --------------------------------------------------- \t" )
    ############################ polynomial of order 3 ########################
    Z3 = 1
    Z2 = - (1 + B)
    Z1 = A
    Z0 = - A*B
    ########################### coefficient of polynomial #####################
    coeff = [Z3, Z2, Z1, Z0] 
    ############################## roots of polynomial ########################
    Z = np.roots(coeff)
    print("\n")
    print(Z)
    #########coefficient or z value for gas phase that is maximum one #########
    ZG = max(Z)
    print ("\n compressibility factor Gas phase is ZG   = " + str(ZG) + " . ")
    print(" \n\t --------------------------------------------------- \t" )
    ###### coefficient or z value for liquid phase that is minimum one  ######
    ZL = min(Z)
    print ("\n compressibility factor liquid phase is ZL   = " + str(ZL) + " . ")
    print(" \n\t --------------------------------------------------- \t" )
    ############################## density of gas phase #######################
    ROg = ((P)*(M))/((ZG)*(R)*(T))
    print ("\ndensity of Gas phase is    = " + str(ROg) + " . ")
    print(" \n\t --------------------------------------------------- \t" )
    ############################## density of liquid phase ####################
    ROl = ((P)*(M))/((ZL)*(R)*(T))
    print ("\ndensity of LIquid phase is    = " + str(ROl) + " . ")
    print(" \n\t --------------------------------------------------- \t" )


def redlich_kwong():
    print("\n")
    print("\n")
    print("\n")
    print(   "\t ########################################################## \t"  )
    print(   "\t ####### Redlich kwong EoS for single component        ####### \t"  )
    print(   "\t ########################################################## \t"  )
    
    MUa = 0.42747
    MUb = 0.08664
    #print (MUa)
    #print (MUb)
    ################################### constant ##############################
    R = 10.73
    ########################## Critical Temperature ###########################
    Tc = 666 
    print ("\n critical temperature value in rankine : "+ str(Tc))
    print(" \n\t --------------------------------------------------- \t" )
    ########################### critical pressure #############################
    Pc = 616.3
    print ("\n critical pressure value in psi : "+ str(Pc))
    print(" \n\t --------------------------------------------------- \t" )
    ##############################  component #################################
    M = 44
    print ("\n moleculer weight of the single component : "+ str(M))
    print(" \n\t --------------------------------------------------- \t" )
    #################### container Pressure and temperature####################
    ############################## pressure value in psi ######################
    #P = input(" \npressure in the container in psi :  ")
    P = 185
    P = int(P)
    print ("\n pressure value in psi : "+ str(P))
    print(" \n\t --------------------------------------------------- \t" )
    ###################### temprerature value in fahrenheit ###################
    T = 100 
    T = 460 + T
    print ("\n temperature value rankine: "+ str(T))
    print(" \n\t --------------------------------------------------- \t" )
    ####################### calculation of parameter a  #######################
    a = (MUa)*((R**2)*(Tc**2.5)/Pc)
    print( "\nvalue of a is  = " + str(a) + " . ")
    print(" \n\t --------------------------------------------------- \t" )
    ####################### calculation of parameter b ########################
    b = (MUb)*((R)*(Tc)/Pc)
    print("\nvalue of b is  = " + str(b) + " . ")
    print(" \n\t --------------------------------------------------- \t" )
    ########################### coefficient A  ################################
    A = ( (a)*(P) )/ ((R**2)*(T**2.5))
    print("\nvalue of coefficient  A  is  = " + str(A) + " . ")
    print(" \n\t --------------------------------------------------- \t" )
    ############################### coefficient B #############################
    B = ((b)*(P))/ ((R)*(T))
    print("\nvalue of coefficient B is  = " + str(B) + " . ")
    print(" \n\t --------------------------------------------------- \t" )
    ####################### polynomial of order 3 #############################
    Z3 = 1
    Z2 = -1
    Z1 = (A-B-(B**2))
    Z0 = - A*B
    ####################### coefficient of polynomial #########################
    coeff = [Z3, Z2, Z1, Z0] 
    ####################### roots of polynomial ###############################
    Z = np.roots(coeff)
    print("\n")
    print(Z)
    print(" \n\t --------------------------------------------------- \t" )
    ########## coefficient or z value for gas phase that is maximum one ####### 
    ZG = max(Z)
    print ("\n compressibility factor Gas phase is ZG   = " + str(ZG) + " . ")
    print(" \n\t --------------------------------------------------- \t" )
    ###### coefficient or z value for liquid phase that is minimum one ########
    ZL = min(Z)
    print ("\n compressibility factor liquid phase is ZL   = " + str(ZL) + " . ")
    print(" \n\t --------------------------------------------------- \t" )
    ###################### density of gas phase  ##############################
    ROg = ((P)*(M))/((ZG)*(R)*(T))
    print ("\ndensity of Gas phase is    = " + str(ROg) + " . ")
    print(" \n\t --------------------------------------------------- \t" )
    ####################### density of liquid phase ######################
    ROl = ((P)*(M))/((ZL)*(R)*(T))
    print ("\ndensity of LIquid phase is    = " + str(ROl) + " . ")
    print(" \n\t --------------------------------------------------- \t" )

#############################################################################################################
#############################################################################################################


def Soave_redlich_kwong():
    print("\n")
    print("\n")
    print("\n")
    print(   "\t ########################################################## \t"  )
    print(   "\t ####### Soave Redlich kwong EoS for single component ####### \t")
    print(   "\t ########################################################## \t"  )
    
    ######################### constant #######################################
    MUa = 0.42747
    MUb = 0.08664
    ############################## w acentricity value and m #################
    w = 1
    m = 0.480 + 1.574*w - 0.176*w**2       
    ######################### constant ####################################### 
    R = 10.73
    ########################## Critical Temperature ###########################
    Tc = 666  
    ########################## critical pressure ###########################
    Pc = 616.3
    ########################## critical pressure ###########################
    M = 44
    ################## container Pressure and temperature #####################
    ###################### pressure value in psi ##############################
    #P = input(" \npressure in the container in psi :  ")
    P = 185
    P = int(P)
    ##################### temprerature value in fahrenheit ####################
    T = 100 
    T = 460 + T
    print (T)
    ######################## Reduced pressure #################################
    Tr = (T / Tc)   
    ############################### alpha value ###############################
    alpha = (1 + m*(  1 - (Tr)**0.5 ) )**2
    ######################## calculation of parameter a #######################
    a = (MUa)*((R**2)*(Tc**2)/Pc)
    print( "\nvalue of a is  = " + str(a) + " . ")
    ######################## calculation of parameter b #######################
    b = (MUb)*((R)*(Tc)/Pc)
    print("\nvalue of b is  = " + str(b) + " . ")
    ############################ coefficient A  ###############################
    A = ( (a)*(P)*(alpha) ) / ((R**2)*(T**2.5))
    print("\nvalue of coefficient  A  is  = " + str(A) + " . ")
    ########################### coefficient B #################################
    B = ((b)*(P))/ ((R)*(T))
    print("\nvalue of coefficient B is  = " + str(B) + " . ")
    ######################## polynomial of order 3 ############################
    Z3 = 1
    Z2 = -1
    Z1 = (A-B-(B**2))
    Z0 = - A*B
    ########################## coefficient of polynomial  #####################
    coeff = [Z3, Z2, Z1, Z0] 
    ########################### roots of polynomial ########################### 
    Z = np.roots(coeff)
    print("\n")
    print(Z)
    ######## coefficient or z value for gas phase that is maximum one #########
    ZG = max(Z)
    print ("\n compressibility factor Gas phase is ZG   = " + str(ZG) + " . ")
    ##### coefficient or z value for liquid phase that is minimum one  ########
    ZL = min(Z)
    print ("\n compressibility factor liquid phase is ZL   = " + str(ZL) + " . ")
    ######################### density of gas phase ##########################
    ROg = ((P)*(M))/((ZG)*(R)*(T))
    print ("\ndensity of Gas phase is    = " + str(ROg) + " . ")
    ####################### density of liquid phase ##########################
    ROl = ((P)*(M))/((ZL)*(R)*(T))
    print ("\ndensity of LIquid phase is    = " + str(ROl) + " . ")

#############################################################################################################
#############################################################################################################
    
    
def redlich_kwong_mixture():
    print(   "\t ################################################################ \t"  )
    print(   "\t ######            Redlich Kwong EoS for liquid mixutre    ###### \t"  )
    print(   "\t ################################################################ \t"  )
    print(" \t --------------------------------------------------- \t" )
    ##########################################################
    ######  components and its moleculer weight         ######
    ##########################################################
    ############################ components of gas ############################
    component = [ 0.45, 0.05, 0.05, 0.03, 0.01, 0.01, 0.40] 
    print("\n component of gas is    =  \n")
    print( component ) 
    print(" \n\t --------------------------------------------------- \t" )
    ####################### moleculer weight of component #####################
    Moleculer_weight = [16.043, 30.070, 44.097, 58.123, 72.150, 84.00, 215] 
    print("\n moleculer weight of gas is    =  \n")
    print( Moleculer_weight ) 
    print(" \n\t --------------------------------------------------- \t" )
    ##########################################################
    ######  critical pressure and temperature           ######
    ##########################################################
    ####################### critical pressure of component ####################
    critical_pressure = [666.4, 706.5, 616.0, 527.9, 488.6, 453, 285]
    print("\n critical_pressure of gas is    =  \n")
    print( critical_pressure ) 
    print(" \n\t --------------------------------------------------- \t" )
    ####################  critical pressure of component ######################
    critical_temperature = [343.33, 549.92, 666.06, 765.62, 845.8, 923, 1287]
    print("\n critical_temperature of gas is    =  \n")
    print( critical_temperature ) 
    print(" \n\t --------------------------------------------------- \t" )
    ############################################################
    ######  am and bm coefficient for mixture             ######
    ############################################################
    ####################### induvidual component a value ##########################
    individual_a_component_list = [(0.42747*((10.73)**2)* (critical_temperature[i])**2.5)/ critical_pressure[i]  for i in range(len(component))]         
    print ("\nThe individual component a value in list is \n: " + str(individual_a_component_list))      
    print(" \n\t --------------------------------------------------- \t" )
    ##################### induvidual value of b coefficient ###################
    individual_b_component_list = [(0.08664*(10.73)* (critical_temperature[i]))/ critical_pressure[i]  for i in range(len(component))]         
    print ("\nThe individual component b coefficient value in list is \n: " + str(individual_b_component_list))      
    print(" \n\t --------------------------------------------------- \t" )
    ###################### total moleculer weight calculation ###############
    individual_moleculer_weight_list = [component[i] * Moleculer_weight[i] for i in range(len(component))]    
    print ("\nThe individual moleculer weight list is \n: " + str(individual_moleculer_weight_list)) 
    print(" \n\t --------------------------------------------------- \t" )
     ############################# a calculation ###########################
    multiplication_of_a_coefficient = [((individual_a_component_list[i])**0.5 * (component[i])) for i in range(len(component))]
    print ("\nThe multiplication of individual a to moleculer fraction list is \n: " + str(multiplication_of_a_coefficient)) 
    print(" \n\t --------------------------------------------------- \t" )
    ############################# b calculation ##############################
    multiplication_of_b_coefficient = [((individual_b_component_list[i]) * (component[i])) for i in range(len(component))]
    print ("\n \nThe multiplication of individual b to moleculer fraction list is \n: " + str(multiplication_of_b_coefficient)) 
    print(" \n\t --------------------------------------------------- \t" )
    ########################################################
    ########  summation of a coefficient       ############# 
    ########################################################
    final_a_coefficient = sum(multiplication_of_a_coefficient)
    final_a_coefficient = final_a_coefficient**2
    print("\n\n ## am coefficient value for mixture " + str(final_a_coefficient) + " ##")
    print(" \n\t --------------------------------------------------- \t" )   
    ########################################################
    ########  summation of b coefficient       ############# 
    ########################################################
    final_b_coefficient = sum(multiplication_of_b_coefficient)
    print("\n\n ## bm coefficient value for mixture" + str(final_b_coefficient) + "\n\n ##")
    print(" \n\t --------------------------------------------------- \t" )
    ########################################################
    ########  total moleculer weight loop      ############# 
    ########################################################  
    ###########################################################
    #### calling function for total moleculer weight     ######
    ########################################################### 
    Total_moleculer_weight = sum(individual_moleculer_weight_list)
    Ma = Total_moleculer_weight
    print("\n")
    print("total moleculer weight " + str(Ma) + ".")
    print(" \n\t --------------------------------------------------- \t" )
    ##################### container Pressure and temperature ##################
    ####################### pressure value in psi  ########################
    #P = input(" \npressure in the container in psi :  ")
    P = int(4000)
    print ("pressure value  :  "+ str(P) + ".")
    R =10.73 
    ####################### temprerature value in fahrenheit ########################
    T = 160 
    T = 460 + T
    print ("temperature value  : "+ str(T) + ".")
    ########################## coefficient A ########################
    A = ( (final_a_coefficient)*(P) )/ ((R**2)*(T**2.5))
    print("\nvalue of coefficient  A  is  = " + str(A) + " . ")   
    print(" \n\t --------------------------------------------------- \t" )
    ########################## coefficient B ########################
    B = ((final_b_coefficient)*(P))/ ((R)*(T))
    print("\nvalue of coefficient B is  = " + str(B) + " . ")   
    print(" \n\t --------------------------------------------------- \t" )
    ######################### polynomial of order 3 ########################
    Z3 = 1
    Z2 = -1
    Z1 = (A-B-(B**2))
    Z0 = - A*B
    ######################### coefficient of polynomial  ######################
    coeff = [Z3, Z2, Z1, Z0] 
    ######################### roots of polynomial ########################
    Z = np.roots(coeff)
    print("\n")
    print(Z)  
    print(" \n\t --------------------------------------------------- \t" )
    ######## coefficient or z value for gas phase that is maximum one #########
    ZG = max(Z)
    print ("\n compressibility factor Gas phase is ZG   = " + str(ZG) + " . ")   
    print(" \n\t --------------------------------------------------- \t" )
    ######## coefficient or z value for liquid phase that is minimum one  #####
    ZL = min(Z)
    print ("\n compressibility factor liquid phase is ZL   = " + str(ZL) + " . ")
    print(" \n\t --------------------------------------------------- \t" )
    ######################### density of gas phase ########################
    ROg = ((P)*(Ma))/((ZG)*(R)*(T))
    print ("\ndensity of Gas phase is    = " + str(ROg) + " . ")
    print(" \n\t --------------------------------------------------- \t" )
    ######################### density of liquid phase ########################
    ROl = ((P)*(Ma))/((ZL)*(R)*(T))
    print ("\ndensity of LIquid phase is    = " + str(ROl) + " . ")
    print(" \n\t --------------------------------------------------- \t" )
   
####################################################################################################################
    
def Soave_redlich_kwong_and_redlich_kwong_mixture():
    print(   "\t ########################################################## \t"  )
    print(   "\t ######      Mixing rule for SRK and RK EoS          ###### \t"  )
    print(   "\t ########################################################## \t"  )
    print(" \t --------------------------------------------------- \t" )
    ##########################################################
    ######           GAS and LIQUID components          ######
    ##########################################################
    ###################### components of mixture #####################
    ###################### components of liquid #####################
    liquid_component = [ 0.45, 0.05, 0.05, 0.03, 0.01, 0.01, 0.40] 
    print("\n liquid component of mixture is    =  \n")
    print( liquid_component ) 
    print(" \n\t --------------------------------------------------- \t" )
    ###################### components of gas #####################
    gas_component = [ 0.86, 0.05, 0.05, 0.02, 0.01, 0.005, 0.005]
    print("\n gas component of mixture is    =  \n")
    print( gas_component ) 
    print(" \n\t --------------------------------------------------- \t" )
    ###########################################################################
    ###############      Gas and liquid moleculer weight         ##############
    ###########################################################################
    #################### moleculer weight of liquid component #################
    liquid_Moleculer_weight = [16.043, 30.070, 44.097, 58.123, 72.150, 84.00, 215] 
    print("\n moleculer weight of liquid is    =  \n")
    print( liquid_Moleculer_weight ) 
    print(" \n\t --------------------------------------------------- \t" )
    #################### moleculer weight of liquid component #################
    gas_Moleculer_weight = [16.043, 30.070, 44.097, 58.123, 72.150, 84.00, 215] 
    print("\n moleculer weight of gas is    =  \n")
    print( gas_Moleculer_weight ) 
    print(" \n\t --------------------------------------------------- \t" )
    ###########################################################################
    ##############       critical pressure and temperature       ##############
    ###########################################################################
    ##################critical pressure of component ##################
    critical_pressure = [ 666.4, 706.5, 616.0, 527.9, 488.6, 453, 285]
    print("\n critical_pressure of gas in psi is    =  \n")
    print( critical_pressure ) 
    print(" \n\t --------------------------------------------------- \t" )
    #####################  critical pressure of component #####################
    critical_temperature = [343.33, 549.92, 666.06, 765.62, 845.8, 923, 1287]
    print("\n critical_temperature of gas is    =  \n")
    print( critical_temperature ) 
    print(" \n\t --------------------------------------------------- \t" )
    ############################################################
    ######  am and bm coefficient for mixture             ######
    ############################################################
    ############################################################
    #########                   Alpha value      ###############
    ######################### constant #######################################
    #MUa = 0.42747
    #MUb = 0.08664
    ######################### constant ####################################### 
    R = 10.73
    #################### container Pressure and temperature####################
    ############################## pressure value in psi ######################
    #P = input(" \npressure in the container in psi :  ")
    P = 4000
    P = int(P)
    print ("\n pressure value in psi : "+ str(P))
    print(" \n\t --------------------------------------------------- \t" )
    ###################### temprerature value in fahrenheit ###################
    T = 160 
    T = 460 + T
    print ("temperature of mixture is " + str(T) + " \n ")
    ############################## w acentricity value and m #################
    w = 0.52
    m = 0.480 + 1.574*w - 0.176*w**2    
    ########################### alpha value #################################
    ################# Alpha values of component ###############################
    ###############################################################################
    #####################  Alpha values of component  ############################
    ###############################################################################    
    alpha_value_list = [(1 + m*(  1 - (T / critical_temperature[i])**0.5 ) )**2  for i in range(len(liquid_component))]         
    print ("\n alpha value for individual component in the list is \n : " + str(alpha_value_list) + "\n")      
    print(" \n\t --------------------------------------------------- \t" )
    ####################### induvidual component a value ######################
    individual_a_component_list = [(0.42747*((10.73)**2)* (critical_temperature[i])**2)/ critical_pressure[i]  for i in range(len(liquid_component))]         
    print ("\nThe individual component a value in list is \n: " + str(individual_a_component_list))      
    print(" \n\t --------------------------------------------------- \t" )
    ####################### induvidual value of b coefficient ######################
    individual_b_component_list = [(0.08664*(10.73)* (critical_temperature[i]))/ critical_pressure[i]  for i in range(len(liquid_component))]         
    print ("\nThe individual component b coefficient value in list is \n: " + str(individual_b_component_list))      
    print(" \n\t --------------------------------------------------- \t" )
    ####################### total liquid moleculer weight calculation ####################### 
    liquid_individual_moleculer_weight_list = [liquid_component[i] * liquid_Moleculer_weight[i] for i in range(len(liquid_component))]    
    print ("\nThe individual moleculer weight list is \n: " + str(liquid_individual_moleculer_weight_list)) 
    print(" \n\t --------------------------------------------------- \t" )
    final_liquid_total_moleculer_weight = sum(liquid_individual_moleculer_weight_list)
    print("\n\n \t final_liquid_total_moleculer_weight : " + str(final_liquid_total_moleculer_weight) + "\t\n")
    print(" \n\t --------------------------------------------------- \t" )
    ####################### total gas moleculer weight calculation #######################
    gas_individual_moleculer_weight_list = [gas_component[i] * gas_Moleculer_weight[i] for i in range(len(gas_component))]    
    print ("\nThe individual moleculer weight list is \n: " + str(gas_individual_moleculer_weight_list))
    final_gas_total_moleculer_weight = sum(gas_individual_moleculer_weight_list)
    print("\n\n \t final_gas_total_moleculer_weight : " + str(final_gas_total_moleculer_weight) + "\t\n ")
    print(" \n\t --------------------------------------------------- \t" )
    #############################################################################################################
    #############                       am and bm                                      ##########################
    #############            coefficient am and bm for liquid calculation                      ##########################
    ############################################################################################################# 
    ####################### total am value ########################################
    X0 = 0
    for i in range(len(liquid_component)):
        #print(" value of i " +  str(i) + "." )
        for j in range(len(liquid_component)):
            #print(" value of j " +  str(j) + ".")
            X1 =  liquid_component[i]*liquid_component[j]*(individual_a_component_list[i]*individual_a_component_list[j]*alpha_value_list[i]*alpha_value_list[j])**0.5
            X0 = X0 + X1
    #print("\t value of (alpha)(am)  aml is " + str(X0) + ".")     
    ####################################################################################################################
    print(" \n\t --------------------------------------------------- \t" )        
    final_aml_coefficient = X0        
    print("\n final (alpha)(am)  aml coefficient value for liquid component in  mixture  :  " + str(final_aml_coefficient) + ("."))  
    print(" \n\t --------------------------------------------------- \t" )
    ############################## coefficient b for liquid calculation #################################
    ####################################################################################################################
    
    multiplication_of_bl_coefficient = [((individual_b_component_list[i]) * (liquid_component[i])) for i in range(len(liquid_component))]
    #print ("\n \nThe multiplication of individual b to moleculer fraction list is \n: " + str(multiplication_of_bl_coefficient)) 
    #print(" \n\t --------------------------------------------------- \t" )
    final_bl_coefficient = sum(multiplication_of_bl_coefficient)
    print("\n final bml coefficient value for liquid component in mixture  : " + str(final_bl_coefficient) + "\n")
    print(" \n\t --------------------------------------------------- \t" )

    #############################################################################################################
    #############                       am and bm                                      ##########################
    #############            coefficient am and bm for gas component calculation                      ##########################
    ############################################################################################################# 
    ####################### total am value ########################################
    Y0 = 0
    for i in range(len(gas_component)):
        #print(" value of i " +  str(i) + "." )
        for j in range(len(gas_component)):
            #print(" value of j " +  str(j) + ".")
            Y1 =  gas_component[i]*gas_component[j]*(individual_a_component_list[i]*individual_a_component_list[j]*alpha_value_list[i]*alpha_value_list[j])**0.5
            Y0 = Y0 + Y1
    #print("value of Y0 is " + str(Y0) + ".")     
    ####################################################################################################################
    print(" \n\t --------------------------------------------------- \t" )        
    final_amg_coefficient = Y0        
    print("\n final (alpha)(am)  amg coefficient value for gas component in mixture : " + str(final_amg_coefficient) + ("."))  
    print(" \n\t --------------------------------------------------- \t" )
    ############################## coefficient b for gas calculation #################################
    ####################################################################################################################
    
    multiplication_of_bg_coefficient = [((individual_b_component_list[i]) * (gas_component[i])) for i in range(len(gas_component))]
    #print ("\n \nThe multiplication of individual b to moleculer fraction list is \n: " + str(multiplication_of_bl_coefficient)) 
    #print(" \n\t --------------------------------------------------- \t" )
    final_bg_coefficient = sum(multiplication_of_bg_coefficient)
    print("\n final bmg coefficient value for gas component in mixture  :  " + str(final_bg_coefficient) + "\n")
    print(" \n\t --------------------------------------------------- \t" )

    #############################################################################################################
    #############                       A and B                                      ##########################
    #############    coefficient Al and Bl cofficient for liquid component calculation  ##########################
    ############################################################################################################# 
       ############################ coefficient A  ###############################
    Al = ( (final_aml_coefficient)*(P)) / ((R**2)*(T**2))
    print("\n value of coefficient  Al for liquid  is  = " + str(Al) + " . ")
    ########################### coefficient B #################################
    Bl = ((final_bl_coefficient)*(P))/ ((R)*(T))
    print("\n  value of coefficient Bl for liquid is  = " + str(Bl) + " . ")
    
    ######################## polynomial of order 3 ############################
    
    #############################################################################################################
    #############                       A and B                                      ##########################
    #############    coefficient Ag and Bg cofficient for Gas component calculation  ##########################
    ############################################################################################################# 
       ############################ coefficient A  ###############################
    Ag = ( (final_amg_coefficient)*(P)) / ((R**2)*(T**2))
    print("\n value of coefficient  Ag for liquid  is  = " + str(Ag) + " . ")
    ########################### coefficient B #################################
    Bg = ((final_bg_coefficient)*(P)) / ((R)*(T))
    print("\n  value of coefficient Bg for liquid is  = " + str(Bg) + " . ")
    
    ######################## polynomial of order 3 for liquid phase  ############################
    Z3 = 1
    Z2 = -1
    Z1 = (Al-Bl-(Bl**2))
    Z0 = - Al*Bl
    ########################## coefficient of polynomial  #####################
    coeff = [Z3, Z2, Z1, Z0] 
    ########################### roots of polynomial ########################### 
    Z = np.roots(coeff)
    print("\n")
    print(Z)
    
     ##### coefficient or z value for liquid phase that is minimum one  ########
    ZL = min(Z)
    print ("\n compressibility factor liquid phase is ZL   = " + str(ZL) + " . ")
   
    
    
    X = Z[:].real
    print("\n\t real part of roots "+str(X)) # outputs 3.6055
    
    print("modulus ")
    X = Z[0]
    x = abs(X)
    print(x)
    X = Z[1]
    x = abs(X)
    print(x)
    X = Z[2]
    x = abs(X)
    print(x)

    ######## coefficient or z value for gas phase that is maximum one #########
    #    ZL = min(Z)
    #    print ("\n compressibility factor liquid phase is ZL   = " + str(ZL) + " . ")
    ##    ####################### density of liquid phase ##########################
    #    ROl = ((P)*(M))/((ZL)*(R)*(T))
    #    print ("\ndensity of LIquid phase is    = " + str(ROl) + " . ")
    #  final_gas_total_moleculer_weight
    
    
        ######################## polynomial of order 3 for gas phase  ############################
    Z3 = 1
    Z2 = -1
    Z1 = (Ag-Bg-(Bg**2))
    Z0 =  Ag*Bg
    ########################## coefficient of polynomial  #####################
    coeff = [Z3, Z2, Z1, Z0] 
    ########################### roots of polynomial ########################### 
    Z = np.roots(coeff)
    print("\n")
    print(Z)
    X = Z[:].real
    
    print("\n\t modulus ")
    X = Z[0]
    x = abs(X)
    print(x)
    X = Z[1]
    x = abs(X)
    print(x)
    X = Z[2]
    x = abs(X)
    print(x)
    print("\n\t real part of roots "+str(X)) # outputs 3.6055
  
     ######## coefficient or z value for gas phase that is maximum one #########
    ZG = max(Z)
    print ("\n compressibility factor Gas phase is ZG   = " + str(ZG) + " . ")
   
    
#############################################################################################################
#############################################################################################################




    


if EOS_model == 1:
    van_der_waal()
elif EOS_model == 2:
    redlich_kwong()
elif EOS_model == 3:
    Soave_redlich_kwong()
elif EOS_model == 4:
    redlich_kwong_mixture() 
elif EOS_model == 5:    
    Soave_redlich_kwong_and_redlich_kwong_mixture()
else:    
     print(" \n Value is not defined \n Thanks Have a Good Day")
#else EOS_model == 3:
 #   print("jai bala ji")
    
 
 #
## constant values 
#MUa = 0.421875
#
#MUb = 0.125
#
#print (MUa)
#
## constant 
#R = 10.73
#
## Critical Temperature 
#
#Tc = 666 
#
## critical pressure 
#Pc = 616.3
#
#
#
#
##component 
##M= 44
#
#
#
#
#
## container Pressure and temperature
#
## pressure value in psi 
#
##P = input(" \npressure in the container in psi :  ")
#P = 185
#P = int(P)
#
## temprerature value in fahrenheit 
#T = 100 
#
#T = 460 + T
#
#print (T)
#
## calculation of parameter a
#a = (MUa)*((R**2)*(Tc**2)/Pc)
#
#print( "\nvalue of a is  = " + str(a) + " . ")
#  
#
## calculation of parameter b
#b = (MUb)*((R)*(Tc)/Pc)
#
#print("\nvalue of b is  = " + str(b) + " . ")
#
#
#
#
### coefficient A 
#
#A = ( (a)*(P) )/ ((R**2)*(T**2))
#
#print("\nvalue of coefficient  A  is  = " + str(A) + " . ")
#
### coefficient B
#
#B = ((b)*(P))/ ((R)*(T))
#
#print("\nvalue of coefficient B is  = " + str(B) + " . ")
#
## polynomial of order 3 
#
#Z3 = 1
#
#Z2 = - (1 + B)
#
#Z1 = A
#
#Z0 = - A*B
#
## coefficient of polynomial 
#coeff = [Z3, Z2, Z1, Z0]
# 
## roots of polynomial 
#Z = np.roots(coeff)
#
#print("\n")
#print(Z)
#
## coefficient or z value for gas phase that is maximum one 
#
#ZG = max(Z)
#
#print ("\n compressibility factor Gas phase is ZG   = " + str(ZG) + " . ")
#
#
## coefficient or z value for liquid phase that is minimum one 
#
#ZL = min(Z)
#
#
#print ("\n compressibility factor liquid phase is ZL   = " + str(ZL) + " . ")
#
#
## density of gas phase 
#
#ROg = ((P)*(M))/((ZG)*(R)*(T))
#
#print ("\ndensity of Gas phase is    = " + str(ROg) + " . ")
#    
#    
## density of liquid phase 
#
#ROl = ((P)*(M))/((ZL)*(R)*(T))
#
#print ("\ndensity of LIquid phase is    = " + str(ROl) + " . ")
#
#

    

    
    