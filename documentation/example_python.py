# included libraries
import numpy as np
import matplotlib.pyplot as plt
import os
##########################
# Make sure the data is printed to file in the C++ code!
##########################
# folder magic
username = "user"
clion_project_folder = "/home/"+ username + "/CLionProjects/NL-LatticeBoltzmann/"
# debug or release
debug_or_release = "cmake-build-relwithdebinfo/"
# debug_or_release = "cmake-build-debug/"
# specific folder
#specific_folder = "test"
specific_folder = "milestones"
milestone = "15"                # Set accordingly
# basic file path
basic_file_path = clion_project_folder + debug_or_release + specific_folder +"/" + milestone
# files
string_ux = "ux_data_file"
string_uy = "uy_data_file"
string_surface = "xy_surfaces"
string_node_dispersions = "node_type_file_"
string_marker_surface = "markers"
# add the basic file path
string_ux = basic_file_path + "/" + string_ux
string_uy = basic_file_path + "/" + string_uy
string_surface = basic_file_path + "/" + string_surface
string_node_dispersions = basic_file_path + "/" + string_node_dispersions
string_marker_surface  = basic_file_path + "/" + string_marker_surface
####################
# size has to be set
size = 800
####################
def plotter_velocity_field():
    #
    print("Velocity field printer")
    # setup data
    ux = np.zeros((size, size))
    uy = np.zeros((size, size))
    # open and add ux data
    file = open(string_ux, 'r')
    lines = file.readlines()
    i = 0
    for line in lines:
        array = np.fromstring(line, dtype=float, sep=" ")
        ux[i] = array
        i += 1
    #
    file.close()
    # open and add uy data
    file = open(string_uy, 'r')
    lines = file.readlines()
    i = 0
    for line in lines:
        array = np.fromstring(line, dtype=float, sep=" ")
        uy[i] = array
        i += 1
    #
    file.close()
    # arrange data
    x = np.arange(0, size)
    y = np.arange(0, size)
    X, Y = np.meshgrid(x, y)
    speed = np.sqrt(ux.T ** 2 + uy.T ** 2)

    # plot surface
    if os.path.exists(string_surface):
        file = open(string_surface, 'r')
        lines = file.readlines()
        for line in lines:
            array = np.fromstring(line, dtype=float, sep=" ")
            # print(array)
            x1 = array[:2]
            y1 = array[2:]
            plt.plot(x1, y1, "k")

    # plot
    plt.streamplot(X, Y, ux.T, uy.T, color=speed, cmap=plt.cm.jet)
    ax = plt.gca()
    ax.set_xlim([0, size])
    ax.set_ylim([0, size])
    titleString = "Title"
    plt.title(titleString)
    plt.xlabel("x-Position")
    plt.ylabel("y-Position")
    fig = plt.colorbar()
    fig.set_label("Velocity u(x,y,t)", rotation=270, labelpad=15)
    plt.savefig('temp.png')
    plt.show()

def plotter_node_allocations():
    print("Node allocation plotter")
    # the enum boundary_t gets converted into numbers this is the last relevant one
    max_tag = 9
    # node allocation
    custom_cmap = plt.matplotlib.colors.ListedColormap(
        ['white', 'blue', 'green', 'red', 'purple', 'orange', 'yellow', 'gray', 'cyan', 'pink'])
    # Create a custom colormap based on the defined colors
    # plot surface
    plt.figure(figsize=(480 / 80, 480 / 80))
    # get in all the files for
    final_array = np.zeros((size, size))
    for i in range(0, max_tag + 1):
        filepath = string_node_dispersions + str(i)
        if (os.path.exists(filepath)):
            print(i)
            # read into a node array
            file = open(filepath, 'r')
            lines = file.readlines()
            node_array = np.zeros((size, size))
            k = 0
            for line in lines:
                node_array[k] = np.fromstring(line, dtype=int, sep=" ")
                k += 1

            # plot
            node_array = np.flip(node_array, axis=(1))
            node_array = np.rot90(node_array, k=1)
            final_array += i * node_array
            # close file
            file.close()

    # plot
    plt.imshow(final_array, cmap=custom_cmap, origin="lower", vmin=0, vmax=max_tag)

    # read and plot surfaces
    if os.path.exists(string_surface):
        file = open(string_surface, 'r')
        lines = file.readlines()
        for line in lines :
            array = np.fromstring(line, dtype=float, sep=" ")
            # print(array)
            x1 = array[:2]
            y1 = array[2:]
            plt.plot(x1,y1, "k")

    # plot and read markers
    if (os.path.exists(string_marker_surface)):
        file = open(string_marker_surface, 'r')
        lines = file.readlines()
        for line in lines:
            array = np.fromstring(line, dtype=float, sep=" ")
            plt.scatter(array[0], array[1], color='red', marker='.')

    # final plot
    ax = plt.gca()
    titleString = "Fancy Törö"
    plt.title(titleString)
    plt.xlabel("x-Position")
    plt.ylabel("y-Position")
    plt.savefig('temp.png')
    plt.show()

# choose what to do
plotter_velocity_field()
plotter_node_allocations()
