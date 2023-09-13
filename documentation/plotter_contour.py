
import numpy as np
import matplotlib.pyplot as plt
import os

# string magic
first_part = "/home/"
second_part = "/CLionProjects/NL-LatticeBoltzmann/cmake-build-relwithdebinfo/milestones/"
# second_part = "/CLionProjects/NL-LatticeBoltzmann/cmake-build-debug/milestones/"
pc = "christoph"
# pc = "sideproject"
basic = first_part + pc + second_part
milestone = "00"
basic = basic + milestone
string_surface = "xy_surfaces"
string_surface = basic + "/" + string_surface
node_dispersions = "node_type_file_"
node_dispersions = basic + "/" + node_dispersions
marker_surface = "markers"
marker_surface = basic + "/" + marker_surface
# the enum boundary_t gets converted into numbers this is the last relevant one
max_tag = 9

def plotter_contours():
    # plot
    print("Making Image")
    size = 400
    # Define a custom colormap for integer values (0 to 10)
    custom_cmap = plt.matplotlib.colors.ListedColormap(
        ['white', 'blue', 'green', 'red', 'purple', 'orange', 'yellow', 'gray', 'cyan', 'pink'])
    # Create a custom colormap based on the defined colors
    # plot surface
    plt.figure(figsize=(480 / 80, 480 / 80))
    # get in all the files for
    final_array = np.zeros((size, size))
    for i in range(0,max_tag + 1):
        filepath = node_dispersions + str(i)
        if(os.path.exists(filepath)):
            print(i)
            # read into a node array
            file = open(filepath, 'r')
            lines = file.readlines()
            node_array = np.zeros((size, size))
            k = 0
            for line in lines:
                node_array[k] = np.fromstring(line,dtype=int,sep = " ")
                k+= 1

            # plot
            node_array = np.flip(node_array, axis=(1))
            node_array = np.rot90(node_array,k = 1)
            final_array += i * node_array
            # close file
            file.close()

    # plot it
    print("ploting")
    # final_array = np.ones_like(final_array)*2
    # plt.imshow(final_array,cmap = custom_cmap,origin="lower",vmin = 0, vmax=max_tag)


    # surface
    if os.path.exists(string_surface):
        file = open(string_surface, 'r')
        lines = file.readlines()
        for line in lines :
            array = np.fromstring(line, dtype=float, sep=" ")
            # print(array)
            x1 = array[:2]
            y1 = array[2:]
            plt.plot(x1,y1, "k")

    # markers in
    if(1):
        if(os.path.exists(marker_surface)):
            file = open(marker_surface,'r')
            lines = file.readlines()
            for line in lines:
                array = np.fromstring(line,dtype=float,sep = " ")
                plt.scatter(array[0],array[1], color = 'red', marker='.')

    if(0):
        # plot some extra lines
        point_1 = np.array([185,185])
        point_2 = np.array([185,325])
        point_3 = np.array([250,25])
        plt.scatter(point_1[0],point_1[1],color = 'black', marker = '.')
        plt.scatter(point_2[0],point_2[1], color='black', marker='.')
        plt.scatter(point_3[0],point_3[1], color = 'black', marker= '.')
        xx = np.array([point_1[0],point_2[0]])
        yy = np.array([point_1[1],point_2[1]])
        plt.plot(xx,yy,"k")
        xx = [point_1[0],point_3[0]]
        yy = [point_1[1],point_3[1]]
        plt.plot(xx,yy,"k")

    ax = plt.gca()
    # standard limits
    if(1):
        ax.set_xlim([-5, size+5 ])
        ax.set_ylim([-5, size+5 ])
    if(0):
        ax.set_xlim([300, 500])
        ax.set_ylim([300, 500])
    if(0):
        ax.set_xlim([500, 700])
        ax.set_ylim([600, 800])
    if(0):
        ax.set_xlim([575, 625])
        ax.set_ylim([650, 750])
    if (0):
        ax.set_xlim([150, 200])
        ax.set_ylim([150, 200])
    if (0):
        ax.set_xlim([75, 125])
        ax.set_ylim([75, 125])
    if (1):
        ax.set_xlim([240, 260])
        ax.set_ylim([290, 310])
    titleString = "Zoom in on the node allocation"
    plt.title(titleString)
    plt.xlabel("x-Position")
    plt.ylabel("y-Position")
    plt.savefig('temp.png')
    plt.show()

plotter_contours()

