
setradcoll 60 2

user Jim Low / Wei Pang/vkp
email vep@ansto.gov.au

title Cr2AlC run 1 full protocol
sampledescription Cr2AlC powder cyl
samplename Cr2AlC 
sampletitle Cr2AlC

ResetWest

# short hold points
sampletitle heat to 150 at 10deg/min
RadCollRampWest 150 600 1
sampletitle hold at 150 for 60min
RadCollTimed 1 60

sampletitle heat to 1100 at 10deg/min
RadCollRampWest 1100 600 1
sampletitle hold at 1100 for 50min
RadCollTimed 1 50

sampletitle heat to 1200 at 10deg/min
RadCollRampWest 1200 600 1
sampletitle hold at 1200 for 100min
RadCollTimed 1 100

# longer hold points

sampletitle heat to 1300 at 10deg/min
RadCollRampWest 1300 600 1
sampletitle hold at 1300 for 100min
RadCollTimed 1 100

sampletitle heat to 1400 at 10deg/min
RadCollRampWest 1400 600 1
sampletitle hold at 1400 for 100min
RadCollTimed 1 100

sampletitle heat to 1500 at 10deg/min
RadCollRampWest 1500 600 1
sampletitle hold at 1500 for 100min
RadCollTimed 1 100

#sampletitle heat to 1550 at 10deg/min
#RadCollRampWest 1550 600 1
#sampletitle hold at 1550 for 100min
#RadCollTimed 1 100
# cool down to rt

sampletitle cool to 500 at 10deg/min
RadCollRampWest 500 600 1

ResetWest

sampletitle cool to RT

RadCollTimed 1 100
RadCollTimed 1 100
RadCollTimed 1 100
RadCollTimed 1 100
RadCollTimed 1 100

