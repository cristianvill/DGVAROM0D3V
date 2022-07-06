########################
# 06/01/2020 A. Alekseenko
# This module contains functions and
# subroutines that create .vtk files
########################

#########################
# my_get_svdsols_file_names
# This subroutine returns a list
# of filenames that are located in
# directory given by path. Files will be
# located in subdirectories and listed with
# path, starting with the provided path.
#
###########################
def my_get_svdsols_file_names (path):
    import os

    files = []

    # r=root, d=directories, f = files
    for r, d, f in os.walk(path):
        for file in f:
            if 'svd_SolET_' in file:
                files.append(os.path.join(r, file))

    return files

###############################################
# make_vtk_files_trajectories
#
# This subroutine will process the pickled SVD encoded trajectories and
# create a vtk file to visualize those trajectories.
#
###############################################
#def make_vtk_files_trajectories (path)