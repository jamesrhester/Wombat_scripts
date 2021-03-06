# proc to set up rad coll

proc ::histogram_memory::pre_count {} {
#        global oct_cycles
#        oscmd start 1
}
proc CollAxisScan {motor start step numsteps} {

	histmem mode time
	histmem preset 15
	set spd [expr 4.0/60.0]
	broadcast $spd
	oct softlowerlim -1.5
	oct softupperlim 1.5
	oct speed $spd
	oct accel $spd
	oct maxretry 0
#	drive oct 1.0

	newfile HISTOGRAM_XY
	for {set i 0} {$i <$numsteps} {incr i} {
		drive $motor [expr $i*$step+$start]
		histmem start block
		save $i
	}
}	

#title BaFe1.74Co0.26As2 at T 25K tilt phi -2 chi +2
sampledescription BaFe1.74Co0.26As2 
samplename BaFe1.74Co0.26As2
user ajs
email ajs@ansto.gov.au

#drive sphi -2.46
#drive schi 91.72

#CollAxisScan som  0 0.2 100
#CollAxisScan som 20 0.2 100
#CollAxisScan som 40 0.2 100
#CollAxisScan som 60 0.2 100
#CollAxisScan som 80 0.2 100

#title BaFe1.74Co0.26As2 at T 25K tilt phi +2 chi -2
#drive sphi 2.46
#drive schi 88.28

#CollAxisScan som 00 0.2 50
#CollAxisScan som 20 0.2 100
#CollAxisScan som 40 0.2 100
#CollAxisScan som 60 0.2 100
#CollAxisScan som 80 0.2 100


title BaFe1.74Co0.26As2 at T30K 

tc1 settle 30
drive tc1 30.0

drive sphi 0
drive schi 90
CollAxisScan som  3.5 0.1 33
drive sphi 2.46
drive schi 88.28
CollAxisScan som  21.5 0.1 33 

title BaFe1.74Co0.26As2 at T40K 
tc1 settle 300
drive tc1 40.0
drive sphi 0
drive schi 90
CollAxisScan som  3.5 0.1 33
drive sphi 2.46
drive schi 88.28
CollAxisScan som  21.5 0.1 33 


title BaFe1.74Co0.26As2 at T60K 

drive tc1 60.0
drive sphi 0
drive schi 90
CollAxisScan som  3.5 0.1 33
drive sphi 2.46
drive schi 88.28
CollAxisScan som  21.5 0.1 33 

title BaFe1.74Co0.26As2 at T80K 

drive tc1 80.0
drive sphi 0
drive schi 90
CollAxisScan som  3.5 0.1 33
drive sphi 2.46
drive schi 88.28
CollAxisScan som  21.5 0.1 33 

title BaFe1.74Co0.26As2 at T100K 

drive tc1 100.0
drive sphi 0
drive schi 90
CollAxisScan som  3.5 0.1 33
drive sphi 2.46
drive schi 88.28
CollAxisScan som  21.5 0.1 33 

title BaFe1.74Co0.26As2 at T120K 

drive tc1 120.0
drive sphi 0
drive schi 90
CollAxisScan som  3.5 0.1 33
drive sphi 2.46
drive schi 88.28
CollAxisScan som  21.5 0.1 33 


# uncomment this to kill rad coll
proc ::histogram_memory::pre_count {} {}
