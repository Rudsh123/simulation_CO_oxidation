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
