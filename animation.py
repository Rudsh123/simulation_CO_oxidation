import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import random
from matplotlib.animation import FuncAnimation, FFMpegWriter 
from matplotlib.text import Text
from matplotlib.patches import Rectangle

def intilize_surface(size):
    surface_coverage = np.zeros((size, size))
    surface_type = np.ones((size, size))
    return surface_coverage, surface_type

def adsorbtion_CO (surface_coverage, pressure_CO, temp,  dt, total_CO_surface, total_empty_surface, surface_type, pressure_O, total_O_surface):
    if total_empty_surface == 0: 
        return surface_coverage
    R = 8.31446261815324
    size_cat = len(surface_coverage) * len(surface_coverage) 
    adsorbtion_rate_of_CO = 3.135 * 10**(5) #[s^-1*mbar^-1]
    sticking_coefficient_for_CO = 1
    
    empty_space_counter = 0
    for i  in range (len(surface_coverage)):
        for j in range(len(surface_coverage[i])):
            if surface_coverage[i][j] == 0: 
                empty_space_counter += 1
    
    
    percent_change_of_CO_adsorbtion = pressure_CO * adsorbtion_rate_of_CO * sticking_coefficient_for_CO * dt *(1-(total_CO_surface/size_cat)**(3)) 
    absolut_CO_adsorbtion = percent_change_of_CO_adsorbtion * size_cat 
    chance_of_CO_adsorbtion = absolut_CO_adsorbtion / empty_space_counter 
    

    
    for i in range(len(surface_coverage)):
        for j in range(len(surface_coverage[i])):
            if surface_coverage[i][j] == 0 and  random.random() <= chance_of_CO_adsorbtion :
                surface_coverage[i][j] = 1
                continue 
    
    return surface_coverage
def adsorbtion_O (surface_coverage, pressure_CO, temp,  dt, total_CO_surface, total_empty_surface, surface_type, pressure_O, total_O_surface, surface_1x1_counter, surface_1x2_counter):
    adsorbtion_rate_of_O = 5.858 * 10**(5) #[s^-1*mbar^-1]
    R = 8.31446261815324
    size_cat = len(surface_coverage) * len(surface_coverage) 

    empty_space_counter_1x1 = 0
    empty_space_counter_1x2 = 0

    for i  in range (len(surface_coverage)):
        for j in range(len(surface_coverage[i])):
            if surface_coverage[i][j] == 0 and surface_type[i][j] == 0: 
                empty_space_counter_1x1 += 1
                
            elif surface_coverage[i][j] == 0 and surface_type[i][j] == 1:
                empty_space_counter_1x2 += 1
    
    sticking_coefficient_for_O_1x1 = 0.6*(surface_1x1_counter/size_cat)
    sticking_coefficient_for_O_1x2 = 0.4* (surface_1x2_counter/size_cat)
    
    chance_of_O_adsorbtion_1x1 = 0
    chance_of_O_adsorbtion_1x2 = 0
    
    absolut_O_adsorbtion_1x1=0
    absolut_O_adsorbtion_1x2=0
    
    if empty_space_counter_1x1+empty_space_counter_1x2 == 0:
        return surface_coverage
    
    percent_1x1_empty = empty_space_counter_1x1/(empty_space_counter_1x1+empty_space_counter_1x2)
    percent_1x2_empty = empty_space_counter_1x2/(empty_space_counter_1x1+empty_space_counter_1x2)
    
    
    if empty_space_counter_1x1 != 0:
        percent_change_of_O_adsorbtion_1x1 = pressure_O * adsorbtion_rate_of_O * sticking_coefficient_for_O_1x1* dt *(1-total_CO_surface/size_cat - total_O_surface/(size_cat * 0.8))**2
        absolut_O_adsorbtion_1x1 = percent_change_of_O_adsorbtion_1x1 * size_cat/2
        chance_of_O_adsorbtion_1x1 = absolut_O_adsorbtion_1x1/(empty_space_counter_1x1)#-(absolut_O_adsorbtion_1x1 + absolut_O_adsorbtion_1x2 )*percent_1x1_empty)
    
    if empty_space_counter_1x2 != 0:
        percent_change_of_O_adsorbtion_1x2 = pressure_O * adsorbtion_rate_of_O * sticking_coefficient_for_O_1x2* dt *(1-total_CO_surface/size_cat - total_O_surface/(size_cat * 0.8))**2
        absolut_O_adsorbtion_1x2 = percent_change_of_O_adsorbtion_1x2 * size_cat/2
        chance_of_O_adsorbtion_1x2 = absolut_O_adsorbtion_1x2/(empty_space_counter_1x2)#-((absolut_O_adsorbtion_1x1 + absolut_O_adsorbtion_1x2 )*percent_1x2_empty))

    for i in range(len(surface_coverage)):
        for j in range(len(surface_coverage[i])):
    
            if surface_coverage[i][j] == 0 and random.random() <= chance_of_O_adsorbtion_1x1 and surface_type[i][j] == 0:
                x_values = list(range(max(0, i - 1), min(len(surface_coverage), i + 2)))
                y_values = list(range(max(0, j - 1), min(len(surface_coverage[i]), j + 2)))

                random.shuffle(x_values)
                random.shuffle(y_values)
                neighbors = []
                for x in x_values:
                    if surface_coverage[i][j] != 0: 
                        break 
                    else:
                        for y in y_values:
                            if (x != i or y != j) and surface_coverage[x][y] == 0:
                                neighbors.append((x, y))
                                if len(neighbors) == 1:
                                    surface_coverage[x][y] = 2
                                    surface_coverage[i][j] = 2
                                    break
                                    
            elif surface_coverage[i][j] == 0 and random.random () <= chance_of_O_adsorbtion_1x2 and surface_type[i][j] == 1 :
                    x_values = list(range(max(0, i - 1), min(len(surface_coverage), i + 2)))
                    y_values = list(range(max(0, j - 1), min(len(surface_coverage[i]), j + 2)))

                    random.shuffle(x_values)
                    random.shuffle(y_values)
                    neighbors = []
                    for x in x_values:
                        if surface_coverage[i][j] != 0: 
                            break 
                        else:
                            for y in y_values:
                                if (x != i or y != j) and surface_coverage[x][y] == 0:
                                    neighbors.append((x, y))
                                    if len(neighbors) == 1:
                                        surface_coverage[x][y] = 2
                                        surface_coverage[i][j] = 2
                                        break
    
    
    return surface_coverage
def desorption (surface_coverage, temp , dt, total_CO_surface, total_empty_surface):
    
    R = 8.31446261815324
    desorption_rate = 2 * 10**(16) * np.exp(- ((38*4184)/(R*temp)))
    size_cat= len(surface_coverage) * len(surface_coverage) 
    
    if total_CO_surface == 0: 
        return surface_coverage
    new_CO_surface =0
    for i in range(len(surface_coverage)):
        for j in range(len(surface_coverage)):
            if surface_coverage[i][j] == 1 :
                new_CO_surface +=1
    
    change_due_desobtion = desorption_rate * dt * total_CO_surface / size_cat #/ total_CO_surface
    absolut_desobtion = change_due_desobtion * size_cat 
    percent_chance_desobtion = absolut_desobtion /(new_CO_surface)

    for i in range(len(surface_coverage)):
        for j in range(len(surface_coverage[i])):
            if surface_coverage[i][j] == 1:
                if random.random() <= percent_chance_desobtion:
                    surface_coverage[i][j]=0
        
    return surface_coverage
  
def reaction( surface_coverage, total_CO_surface, total_O_surface, dt, temp): 
    R=8.31446261815324 #[J/mol*K]
    k_r=3 * 10**(6) * np.exp(- ((10*4184)/(R*Temp)))  # [s^-1]
    size_cat= len(surface_coverage) * len(surface_coverage) 
    count_of_reactant = 0
    rows, cols = len(surface_coverage), len(surface_coverage[0])
    
    list_of_neighbor_CO = []
    
    
    
    for i in range(rows):
        for j in range(cols):
            if surface_coverage[i][j] == 1:
                neighbor_found=False
                for x in range(max(0, i - 1), min(rows, i + 2)):
                    if neighbor_found is True: 
                        break
                    for y in range(max(0, j - 1), min(cols, j + 2)):    
                        if (x != i or y != j) and surface_coverage[x][y] == 2:
                            count_of_reactant += 1
                            neighbor_found =True 
                            list_of_neighbor_CO.append([i,j])
                            break
                            
    
    if count_of_reactant == 0 or total_CO_surface == 0 or total_O_surface == 0:
        return surface_coverage

    
    
    change_due_reaction = k_r * dt* (total_CO_surface/size_cat) * (total_O_surface/size_cat)
    absolut_reaction = change_due_reaction * size_cat
    chance_of_reaction = absolut_reaction / count_of_reactant
    

    for i in list_of_neighbor_CO:
        if random.random() <= chance_of_reaction:

            x_values = list(range(max(0, i[0] - 1), min(len(surface_coverage), i[0] + 2)))
            y_values = list(range(max(0, i[1] - 1), min(len(surface_coverage), i[1] + 2)))

            random.shuffle(x_values)
            random.shuffle(y_values)

            for x in x_values:
                if surface_coverage[i[0]][i[1]] == 0:
                    break
                else:
                    for y in y_values:
                        if (x != i or y != j) and surface_coverage[x][y] == 2:
                            surface_coverage[i[0]][i[1]] = 0
                            surface_coverage[x][y] = 0
                            break
                        else:
                            continue
                                        
    return surface_coverage
  
def phaseshift (surface_type, dt, Temp, total_CO_surface, total_1x1_surface, total_1x2_surface, surface_coverage):
    size_cat= len(surface_type) * len(surface_type) 
    R=8.31446261815324 #[J/mol*K]
    rate_of_phaseshift = 10**2 * np.exp(-7*4184/(Temp*R))
    r3 = -(1/0.0135)
    r = [-0.026*r3, 0.3*r3, -1.05*r3, r3 ] #r0, r1, r2, r3
    poly_coefficent = 0 
    free_1x1_surf = 0
    free_1x2_surf = 0
    
    for i in range(len(surface_type)):
        for j in range(len(surface_type)):
            if surface_coverage[i][j] == 0 and surface_type[i][j] ==0:
                free_1x1_surf +=1
            elif surface_coverage[i][j] ==0 and surface_type[i][j] ==1:
                free_1x2_surf +=1 
    
    if total_CO_surface/size_cat<=0.2:
        percent_change_1x1=-total_1x1_surface/size_cat*rate_of_phaseshift
    elif total_CO_surface/size_cat>=0.5:
        percent_change_1x1=(1-total_1x1_surface/size_cat)*rate_of_phaseshift
    else:
        for i in range(len(r)):
            poly_coefficent += r[i]*(total_CO_surface/size_cat)**i
        percent_change_1x1 = (poly_coefficent - total_1x1_surface/size_cat)*rate_of_phaseshift
        
    absolut_change_1x1 = percent_change_1x1 * size_cat *dt
    
    if absolut_change_1x1 >=0: 
        
        percent_chance_1x1 = absolut_change_1x1 / total_1x2_surface 
        
        for i in range(len(surface_type)):
            for j in range (len(surface_type)):
                if surface_type[i][j] == 1 and random.random() <= percent_chance_1x1:
                    surface_type[i][j] = 0
        
    else: 
        percent_chance_1x1 = absolut_change_1x1 / total_1x1_surface  *(-1)
        for i in range(len(surface_type)):
            for j in range(len(surface_type)):
                if surface_type[i][j] == 0 and random.random() <= percent_chance_1x1:
                    surface_type[i][j] = 1
                    

    return surface_type


time_to_run_over = 11 #[s]
dt = 0.001 #[s]
surface_coverage, surface_type = intilize_surface(100)
size_cat = len(surface_coverage) * len(surface_coverage) 

fig, ax = plt.subplots(figsize=(8, 4)) 
custom_cmap = plt.cm.colors.ListedColormap(['white', 'blue', 'red'])
im = ax.imshow(surface_coverage, cmap=custom_cmap, vmin=0, vmax=2, animated=True)
plt.rcParams['animation.embed_limit'] = 10000

CO_counter = 0
O_counter = 0
empty_coverage = size_cat
surface_1x1_counter = 0
surface_1x2_counter= size_cat

Temp = 539.5
p_O = 6.67 * 10**(-5) #[mbar]
p_CO = 3 * 10**(-5) #[mbar]

CO_counter = 0
O_counter = 0
surface_1x1_counter = 0

bar_CO = Rectangle((415, 50 ), 0 , 10, color='blue')
bar_O = Rectangle((415, 50), 0 , 10, color='red')

fig.patches.extend([bar_CO, bar_O])

timer_text = ax.text(1.02, 1.02, '', transform=ax.transAxes, fontsize=12, color='black')
current_time=0
skip_frames = 10
update_counter =skip_frames-1 


def update(i):
    global surface_coverage, surface_type, CO_counter, O_counter, empty_coverage, surface_1x1_counter, surface_1x2_counter, current_time, p_CO, p_O, Temp, bar_CO, bar_O, update_counter, skip_frames
    update_counter +=1
    
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
    
    current_time += dt
    
    if  update_counter == skip_frames:
        bar_CO.set_x(415)
        bar_CO.set_y(80)
        bar_CO.set_width(CO_counter / size_cat * 100)

        bar_O.set_x(415)
        bar_O.set_y(60)


        bar_O.set_width(O_counter / size_cat * 100)


        timer_text.set_text(f'Time: {current_time:.2f} s')

        im.set_array(surface_coverage)
        update_counter = 0

        
    return im, bar_CO, bar_O, timer_text  

anim = FuncAnimation(fig, update, frames=range(int(time_to_run_over  / dt)), interval=1) # is time to run over
anim.save("animation.gif", writer='pillow', fps=60) #animation.gif   is the name the animation is saved under.
