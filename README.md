# Dependencies

- gcc
- GSL (The GNU Scientific Library)
- NLOPT (The Non-Linear Optimisation Library)
- Python 3 with numpy, matplotlib, sys, and subprocess modules (optional, Used for the helper script)

The developer for this code will try to decrease the requirements with time and try to remove the requirement of NLOPT.

# Compiling

1. Go into the main directory
2. run `make`
3. It will generate `resonanceFactor`. This is the program that you need to run with your data.

# Using the Program

## Using just the C program

In the terminal (in Windows, open it by Shift+Right Click and then  select Open Linux shell here), go to where the `resonanceFactor` executable binary is.
Type ./resonanceFactor \<filename for data file\> \<Cable delay in ns\> \<Detrend Mode\> \<Number of Initial Points fro Detrending\> \<Number of Final Points fro Detrending\>
This will generate the output.

Eg. `./resonanceFactor data.txt 30.0 2 100 150`
This will take in input from the file data.txt ,correct for a cable delay of 30ns, and then perform quadratic detrending using the first 100 and the final 150 datapoints to find Qr, Qc from the phiRM and the DCM methods.

### Input Format

The data file should be a space separated, tab separated or a new-line separated file  with no headers, and data in the form freq(Arb Units) Re(S21) Im(S21). For example, see exampleData.txt.

## Using the Python Script with Labber Outputs

As of now, the script handles only outputs from Labber for Power sweeps and single runs. To use it,
```
python3 analyse_labber_files.py <Input Labber File Name>
```

This will generate and show the plots of the sweeps. The magnitude plots would be saved.
The resonanceFactor executable binary should be in the same folder as the python script.

For example, to analyse data from `data.txt`, use

```
python3 analyse_labber_files.py data.txt
```

This would save the generated magnitude plots in a folder called Plots.

If you want logmag plots, just use `logmag = True` at the top of the script

# Possible errors / Problems during compilation

1. You do not have GCC. For most linux systems, you can use the command `sudo apt install build-essential` or the appropriate variant according to your distribution. For Windows, either use this command in WSL/WSL2, or use MinGW.
2. You do not have GSL. For most linux systems, you can use the command `sudo apt install libgsl-dev` or the appropriate variant according to your distribution. For Windows, either use this command in WSL/WSL2, or look up the appropriate method here <https://www.gnu.org/software/gsl/>
3. You do not have NLOPT. Look up the method here <https://nlopt.readthedocs.io/en/latest/NLopt_Installation/>
4. The GSL installation is not in the right location. To correct for this, go to phiRM_DCM.c and change all the headers of the form \<gsl/xxxx.h\> to just \<xxxx.h\>. That might help. Otherwise, look up how to do this.

# Reminder

This program is under development. In case if you find any issues, please let me know at my email address <samarthh@iisc.ac.in>. Please keep the Subject as \`Bug Report for phiRM'.