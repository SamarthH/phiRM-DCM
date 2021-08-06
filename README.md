Dependencies of this program
You need to have the following libraries on your machine to compile and run this program:
- gcc
- GSL (The GNU Scientific Library)
- NLOPT (The Non-Linear Optimisation Library)

The developer for this code will try to decrease the requirements with time and try to remove the requirement of NLOPT.

Compiling the program:
Go into the main directory, and run `make`. It will generate `resonanceFactor`. This is the program that you need to run with your data.

How to use the program:

In the terminal (in Windows, open it by Shift+Right Click and then  select Open Linux shell here), go to where the `resonanceFactor` executable binary is.
Type ./resonanceFactor <filename for data file> <Cable delay in ns>
This will generate the output.

Eg. `./resonanceFactor data.txt 30.0 2 100 150`
This will take in input from the file data.txt ,correct for a cable delay of 30ns, and then perform quadratic detrending using the first 100 and the final 150 datapoints to find Qr, Qc from the phiRM and the DCM methods.

Format of the data file:
The data file should be a space separated, tab separated or a new-line separated file  with no headers, and data in the form freq(Arb Units) Re(S21) Im(S21). For example, see exampleData.txt.

Possible errors / Problems during compilation:
1. You do not have GCC. For most linux systems, you can use the command `sudo apt install build-essential` or the appropriate variant according to your distribution. For Windows, either use this command in WSL/WSL2, or use MinGW.
2. You do not have GSL. For most linux systems, you can use the command `sudo apt install libgsl-dev` or the appropriate variant according to your distribution. For Windows, either use this command in WSL/WSL2, or look up the appropriate method here https://www.gnu.org/software/gsl/
3. You do not have NLOPT. Look up the method here https://nlopt.readthedocs.io/en/latest/NLopt_Installation/
4. The GSL installation is not in the right location. To correct for this, go to phiRM_DCM.c and change all the headers of the form <gsl/xxxx.h> to just <xxxx.h>. That might help. Otherwise, look up how to do this.

Remember this:
This program is under development. In case if you find any issues, please let me know at my email address samarthh@iisc.ac.in. Please keep the Subject as "Bug Report for phiRM".
