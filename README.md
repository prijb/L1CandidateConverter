# Candidate converter code 

## Function
This code was made during the DARTS hackathon (26/08-30/08) for the purpose of converting L1 objects (Jets, EGammas, Muons) to generic "Candidates". These candidates are then to be used for training a classifier for particle identification.

## Header Inclusion
The generic Candidate class header is included in the main function through ProcessLine statements. The code is compiled as follows:
```
g++ -o build/convert_candidates.exe src/convert_candidates.cpp -Iinclude $(root-config --cflags --libs)
```