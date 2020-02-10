"""
open master equation results
"""
import sys,os
import me_parser

####################################################################################################

#initialize the class in which to store the results
print('\n')
data = me_parser.paper()
data.reactions = []
# set some constants, depending upon whether the rate coefficients are to be used for CHEMKIN or something else.
chemkin = 1
if chemkin==1:
    data.T0 = 1.0
    data.R = 1.987 # cal/mol-K.  Note that PLOG formalism requires Ea in cal/mol-K!
    data.N_avo = 6.0221415E23 #convert bimolecular rate coefficients from cm^3/sec to cm^3/mol/s
elif chemkin==0:
    data.T0 = 298.0
    data.R = 1.0 # K^-1.
    data.N_avo = 1.0  #leave bimolecular rate coefficients in cm^3/sec

# set the minimum and maximum temperature
#data.Tmin = 1000.0
#data.Tmax = 2000.0

# read me.out file from the command line
command_line = sys.argv[1:]
me_dot_out = command_line[0]
results = open(me_dot_out,'r')
lines = results.readlines()
results.close()

# copy new plog executable to the path of the source file
path = os.path.abspath(os.path.dirname(me_dot_out))

# parse results for the temperature, pressure, and names of channels
me_parser.get_temp_pres(data,lines)
# parse results for the pdep rate constants
me_parser.get_pdep_k(data,lines)
# fit the results to PLOG expressions
me_parser.fit_pdep(data,nonlin_fit=True) #replace <True> with <False> if you don't want to use the nonlinear solver (not recommended)
# print the resulting PLOG expressions to file
me_parser.print_plog(data, me_dot_out)
# plot the results: dashed line = single PLOG, solid line = double PLOG
me_parser.plot_reactant(data, me_dot_out, show_plot=False, save_plot=True)
