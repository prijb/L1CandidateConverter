# Candidate converter code 

## Function
This code was made during the DARTS hackathon (26/08-30/08) for the purpose of converting L1 objects (Jets, EGammas, Muons) to generic "Candidates". These candidates are then to be used for training a classifier for particle identification

## Issues
The generic Candidate class header is not identified by the `convert_candidates.cpp` script. To repeat this error, perform the following steps
```
rootcling -f src/CandidatesDict.cxx -s build/CandidatesDict_rdict.pcm -Iinclude include/Candidates.h include/CandidatesLinkDef.h
g++ -shared -fPIC -o build/CandidatesDict.so src/CandidatesDict.cxx $(root-config --cflags --libs)
g++ -o build/convert_candidates.exe src/convert_candidates.cpp build/CandidatesDict.so -Iinclude $(root-config --cflags --libs)
build/convert_candidates.exe /vols/cms/pb4918/StoreNTuple/L1Scouting/PtRegression/QCD_15to7000 1.0 /vols/cms/pb4918/L1Scouting/Aug24/DARTS/pid_project/convert_candidates.root
```
I'm not sure if these compilation commands are the correct ones, so any changes are appreciated