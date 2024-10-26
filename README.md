# Spin overlap calculator using Python
Firstly, you must obtained the calculated values of individual spin texture by the VASP as input files.

There are two forms of input files, one for calculating only the wave function (example1.txt) and the other for calculating both the wave function and spin overlap (example2.txt). Next, have a look at this example for calculating both the wave function and spin overlap (example2.txt):

The first line is fixed, that is, "Sx Sy Sz". The second line is the name of the type to be calculated. The name of this line will determine how many groups of data to calculate. Because the number of data is separated by spaces, please do not use spaces to separate the corresponding words for the same name, otherwise it will cause problems in reading data. For example, please write "MD 0.42 ps" as "MD0.42ps" or "MD_0.42_ps". Starting from the third line of data, the period of spin overlap is calculated for every three lines, which are the names of the calculation groups, the spin textures on the valence band (VB), and the conduction band (CB). When the input file is example1.txt, the data period becomes two lines, which are the group name and the data value at the valence band or conduction band.

Please prepare the input file for your calculation system strictly according to the format of the example input file.
After preparing all input files, please follow the following sequence:
1. Please prepare your input_file. In this example is 'spin texture with both CBM and VBM data.txt' (example2.txt) or 'spin texture with only CBM or VBM data.txt' (example1.txt).
2. Place your input_file and 'Spin overlap.py' in the same paths. Then, runing this Python script in your computer.
3. As a result, you will obtain the result file "Result.txt".
