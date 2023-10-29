Temp = 539.5
p_O = 6.67 * 10**(-5) #[mbar]
p_CO = 3 * 10**(-5) #[mbar]
time_to_run_over = 10 #[s]
dt = 0.001 #[s]
surface_coverage, surface_type = intilize_surface(500)

size_cat = len(surface_coverage) * len(surface_coverage) 
list_of_CO = []
list_of_O = []
list_of_1x1 = []
CO_counter = 0
O_counter = 0
surface_1x1_counter = 0
empty_coverage = size_cat
surface_1x2_counter = size_cat

for  i in np.arange( 0, time_to_run_over, dt ):
    surface_type = phaseshift (surface_type, dt, Temp, CO_counter, surface_1x1_counter, surface_1x2_counter, surface_coverage )
    surface_coverage = adsorbtion_CO(surface_coverage, p_CO, Temp, dt, CO_counter, empty_coverage,surface_type, p_O, O_counter)
    surface_coverage = adsorbtion_O(surface_coverage, p_CO, Temp, dt, CO_counter, empty_coverage, surface_type, p_O, O_counter, surface_1x1_counter, surface_1x2_counter)
    surface_coverage = desorption (surface_coverage, Temp, dt, CO_counter, empty_coverage)
    surface_coverage = reaction(surface_coverage, CO_counter, O_counter, dt, Temp)

    CO_counter = 0
    O_counter = 0
    empty_coverage = 0
    surface_1x1_counter = 0
    surface_1x2_counter=0
    
    for i in range (len(surface_coverage)):
        for j in range (len(surface_coverage)):
            if surface_coverage[i][j] == 1:
                CO_counter += 1
            elif surface_coverage[i][j]== 0:
                empty_coverage += 1
            elif surface_coverage[i][j] ==2:
                O_counter += 1

    for i in range(len(surface_type)):
        for j in range(len(surface_type)):
            if surface_type[i][j] == 0: 
                surface_1x1_counter += 1
            elif surface_type[i][j] == 1 : 
                surface_1x2_counter +=1
    

    list_of_CO.append(CO_counter/size_cat)
    list_of_O.append(O_counter/size_cat)
    list_of_1x1.append(surface_1x1_counter/size_cat)
    
time_run = np.arange( 0, time_to_run_over, dt )
plt.plot(time_run, list_of_CO, "b-" )
plt.plot(time_run, list_of_O, "r-" )
plt.plot(time_run, list_of_1x1, "g--" )
# plt.xscale("log")
plt.xlabel("Time [s]")
plt.ylabel( "red = O ; blue= CO; green = 1x1" )
plt.show()
