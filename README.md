# TADClam

***Overlapped and nested topologically associating domains detection with community affiliation model

Requirements
enviorment : Python3.8 or above.

packages : including numpy, multiprocessing, and scipy.

Usage
python TADClam.py -in Hi_C_matrix

#Parameters

-in : the input file of a N*N Hi-C matrix separated by TAB for a chromosome i.

-threshold : the threshold for filtering candidate TADs. (Default value is 2.50)

-iter : the iteration number of table F. (Default value is 1500)

-out : the output file for the predicted TAD of chromosome i, where a line represents a TAD containing two columns that represent the start bin and the end bin of a TAD.

contact: zhaoling2-c@my.cityu.edu.hk
