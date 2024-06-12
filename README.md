# TADClam

***Overlapped and nested topologically associating domains detection with community affiliation model

Requirements
enviorment : Python3.7 or above.

packages : including numpy .

Usage
python TADClam.py in_Hi_C Ft_threshold

###Parameters

Hi-C matrix : the input file of a N*N Hi-C matrix separated by TAB for a chromosome i.

Ft_threshold : the threshold for filtering candidate TADs. (Default value is 2.50)

ou_TADs : the output file for the predicted TAD of chromosome i, where a line represents a TAD containing two columns that represent the start bin and the end bin of a TAD.



## Notes

This tool requires only a single user-defined parameter **CS_threthold**, which enable users to use easily. The default value of this parameter is 0.8.




