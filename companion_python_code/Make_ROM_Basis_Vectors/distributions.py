##############################
###  12/30/2020 A. Alekseenko
###  This is a collection of subroutines to implement
###  velocity distribution functions.
###
###
##############################


####################################
# dimless_mxwln(nodes_u, nodes_v, nodes_w,n,ubar,vbar,wbar,T)
#
# This subroutine implements dimensionless Maxwellian
# nodes_u, nodes_v, nodes_w - numpy arrays containing nodes to evaluate the
#                             velocity distribution
# n,ubar,vbar,wbar, T --- macroparameters
#
###################################
def dimless_mxwln(nodes_u, nodes_v, nodes_w,n,ubar,vbar,wbar,T):
    #####
    import numpy as np
    pi25DT = 3.141592653589793238462643
    beta = np.sqrt(pi25DT * T) * (pi25DT * T)
    y = n * np.exp(-((nodes_u - ubar )**2 + (nodes_v - vbar)**2 + (nodes_w - wbar)**2) / max(T, 0.0000001)) / beta
    return y