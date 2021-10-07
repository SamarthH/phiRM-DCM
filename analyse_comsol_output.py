import numpy as np
import matplotlib.pyplot as plt
import subprocess as sp
from sys import argv

cleanup = False # Keep this as False if you need the txt output by the C program for later
plotmag = True # Keep this False to stop plotting magnitude
plotre = True # Keep this False to stop plotting real part
plotim = True # Keep this False to stop plotting imag part
savemag = True # Keep this False if you do not want magnitude plots saved
logmag = True # Make this True if you want logmag plots
showplots = False # Make this false if you do not want to show the plots and only save them

def analyseData(data, foldername= "Plots"): # data needs to be a numpy array with 2 rows: frequency and S21 (complex)
	if np.real(data[1][0]) < 0:
		data = -data
	input_to_prog = np.vstack((np.abs(data[0]), np.real(data[1]), np.imag(data[1])))
	input_to_prog = input_to_prog.T
	np.savetxt("infile.analyse", input_to_prog)
	try:	
		# Runs it assuming that the program is copied to bin
		params = sp.Popen(["resonanceFactor", "infile.analyse", "0", "1"], stdout=sp.PIPE).communicate()[0] # Performs linear detrend
	except:
		params = sp.Popen(["./resonanceFactor", "infile.analyse", "0", "1"], stdout=sp.PIPE).communicate()[0] # Performs linear detrend
	outfile = np.loadtxt("exp.txt", unpack=True)

	#sp.Popen(["rm", "infile.analyse"])
	sp.Popen(["mkdir", "-p", foldername])
	if cleanup:
		sp.Popen(["rm", "exp.txt"])
	else:
		sp.Popen(["mv", "exp.txt", foldername+"/exp_.txt"])

	fr, Qr, Qc, Qi, chi2 = (params.decode().split('\n'))[1].split('\t')

	sp.Popen(["mkdir", "-p", foldername + '/'])

	# Plotting Real Part
	if showplots and plotre:
		plt.scatter(outfile[0], outfile[1])
		plt.plot(outfile[0], outfile[4], color = 'red', label = ("f = {}\nQr = {}\nQc = {}\nQi = {}\nchi2 = {}".format(fr, Qr, Qc, Qi, chi2)))
		plt.title("Real Part")
		plt.xlabel("f (Comsol Units)")
		plt.ylabel("Re(S21)")
		plt.legend()
		plt.show()
		plt.clf()
		plt.cla()
	# Plotting Imag Part
	if showplots and plotim:
		plt.scatter(outfile[0], outfile[2])
		plt.plot(outfile[0], outfile[5], color = 'red', label = ("f = {}\nQr = {}\nQc = {}\nQi = {}\nchi2 = {}".format(fr, Qr, Qc, Qi, chi2)))
		plt.title("Imaginary Part")
		plt.xlabel("f (Comsol Units)")
		plt.ylabel("Im(S21)")
		plt.legend()
		plt.show()
		plt.clf()
		plt.cla()
	# Plotting Magnitude
	if plotmag == True and logmag == False:
		plt.scatter(outfile[0], outfile[3])
		plt.plot(outfile[0], outfile[6], color = 'red', label = ("f = {}\nQr = {}\nQc = {}\nQi = {}\nchi2 = {}".format(fr, Qr, Qc, Qi, chi2)))
		plt.title("Magnitude")
		plt.xlabel("f (Comsol Units)")
		plt.ylabel("|S21|")
		plt.legend()
		if savemag:
			plt.savefig(foldername + '/'+foldername+".png")
		if showplots:	
			plt.show()
		plt.clf()
		plt.cla()

	if plotmag == True and logmag == True:
		plt.scatter(outfile[0], 20*np.log10(outfile[3]))
		plt.plot(outfile[0], 20*np.log10(outfile[6]), color = 'red', label = ("f = {}\nQr = {}\nQc = {}\nQi = {}\nchi2 = {}".format(fr, Qr, Qc, Qi, chi2)))
		plt.title("Magnitude")
		plt.xlabel("f (Comsol Units)")
		plt.ylabel("|S21| (dB)")
		plt.legend()
		if savemag:
			plt.savefig(foldername + '/'+foldername+".png")
		if showplots:
			plt.show()
		plt.clf()
		plt.cla()

	#return (params, outfile)
	
infile_name = (str(argv[1])) # Name of input file
infile = open(infile_name)
infile_lines = infile.readlines()

# Getting the number of lines to skip
to_skip = 0
while infile_lines[to_skip][0][0] == '%':
	to_skip += 1

dat_raw = np.loadtxt(infile_name, skiprows=to_skip, unpack=True, dtype=np.complex64)
dat = dat_raw[:,:len(dat_raw[1])//2]
dat[1] = dat[1] + 1.0j*dat_raw[1,len(dat_raw[1])//2:]
analyseData(dat, infile_name[:-4])