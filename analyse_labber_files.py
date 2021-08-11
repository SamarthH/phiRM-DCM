import numpy as np
import matplotlib.pyplot as plt
import subprocess as sp
from sys import argv

cleanup = True # Keep this as False if you need the txt output by the C program for later
plotmag = True # Keep this False to stop plotting magnitude
plotre = True # Keep this False to stop plotting real part
plotim = True # Keep this False to stop plotting imag part
savemag = True # Keep this False if you do not want magnitude plots saved

def analyseData(data, power = ""): # data needs to be a numpy array with 2 rows: frequency and S21 (complex)
	if np.real(data[1][0]) < 0:
		data = -data
	input_to_prog = np.vstack((np.abs(data[0]), np.real(data[1]), np.imag(data[1])))
	input_to_prog = input_to_prog.T
	np.savetxt("infile.analyse", input_to_prog)
	params = sp.Popen(["./resonanceFactor", "infile.analyse", "0", "1"], stdout=sp.PIPE).communicate()[0] # Performs linear detren
	outfile = np.loadtxt("exp.txt", unpack=True)

	sp.Popen(["rm", "infile.analyse"])
	if cleanup:
		sp.Popen(["rm", "exp.txt"])
	else:
		sp.Popen(["mv", "exp.txt", "exp_"+power+".txt"])

	fr, Qr, Qc, Qi, chi2 = (params.decode().split('\n'))[1].split('\t')

	sp.Popen(["mkdir", "-p", "Plots"])

	# Plotting Real Part
	if plotre:
		plt.scatter(outfile[0], outfile[1])
		plt.plot(outfile[0], outfile[4], color = 'red', label = ("f = {}\nQr = {}\nQc = {}\nQi = {}\nchi2 = {}".format(fr, Qr, Qc, Qi, chi2)))
		plt.title("Power = "+power+" (Real Part)")
		plt.legend()
		plt.show()
	# Plotting Imag Part
	if plotim:
		plt.scatter(outfile[0], outfile[2])
		plt.plot(outfile[0], outfile[5], color = 'red', label = ("f = {}\nQr = {}\nQc = {}\nQi = {}\nchi2 = {}".format(fr, Qr, Qc, Qi, chi2)))
		plt.title("Power = "+power+" (Imaginary Part)")
		plt.legend()
		plt.show()
	# Plotting Magnitude
	if plotmag:
		plt.scatter(outfile[0], outfile[3])
		plt.plot(outfile[0], outfile[6], color = 'red', label = ("f = {}\nQr = {}\nQc = {}\nQi = {}\nchi2 = {}".format(fr, Qr, Qc, Qi, chi2)))
		plt.title("Power = "+power+" (Magnitude)")
		plt.legend()
		if savemag:
			plt.savefig("Plots/"+power+".png")
		plt.show()

	#return (params, outfile)

#print(str(argv[1]))
infile_name = (str(argv[1])) # Name of input file
infile = open(infile_name)
infile_lines = infile.readlines()

# Checking if this is a power sweep
checking_words = (infile_lines[2].split(' '))[-2] # Checking the last word of the third line
is_power_sweep = False
if checking_words == 'Output':
	is_power_sweep = True
	power_sweep = infile_lines[3].split('\t')

if is_power_sweep:
	dat = np.loadtxt(infile_name, skiprows=5, dtype = np.complex)
	for i in range(len(power_sweep)):
		analyseData(np.vstack((dat[0], dat[i+1])), power_sweep[i])
else:
	dat = np.loadtxt(infile_name, skiprows=3, dtype = np.complex)
	analyseData(np.vstack((dat[0], dat[1])), "Unknown")