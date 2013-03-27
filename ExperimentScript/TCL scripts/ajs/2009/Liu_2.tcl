

title PLZT in cryofurnace 10-450K 
sampletitle PLZT
samplename mtth 90 Ge 115 1.54A
user Liu / ajs
email ajs@ansto.gov.au

proc RadCollRampCF1 {start step numsteps oscno} {

	RadCollOn $oscno
	histmem mode unlimited
	newfile HISTOGRAM_XY
	for {set i 0} {$i < $numsteps} {incr i} {
		set j [expr $i*$step+$start]
		tc1_asyncq send SETP 1,$j
		tc2_asyncq send SETP 1,$j
		histmem start block
		save $i
	}
	RadCollOff
}

SetRadColl 60 2
tc1 controlsensor sensorA
tc2 controlsensor sensorA
tc1 range 5
tc2 range 5

tc1 tolerance 50
tc2 tolerance 50
tc1 settle 0
tc2 settle 0

drive tc1 10 tc2 10


wait 3600

tc1_asyncq send RAMP 1,0,1.2
wait 3
tc1_asyncq send SETP 1,10
wait 3
tc1_asyncq send RAMP 1,1,1.2
wait 3
tc2_asyncq send RAMP 1,0,1.2
wait 3
tc2_asyncq send SETP 1,10
wait 3
tc2_asyncq send RAMP 1,1,1.2
wait 3

#
wait 3




sampletitle ramp to 450K
RadCollRampCF1 10 1 440 1

#turn ramp off
tc1_asyncq send RAMP 1,0,1.0
wait 3
tc2_asyncq send RAMP 1,0,1.0
wait 3
drive tc1 300 tc2 300
tc1_asyncq send RELAY 2,2,1
