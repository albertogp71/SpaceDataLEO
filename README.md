# SpaceDataLEO
***

## Description
SpaceDataLEO is a simple satellite network simulator.
It builds a Walker Delta constellation where each satellite has four terminals for inter-satellite communication.
The simulator is a single-file Matlab program where the user specifies the following parameters:
 1. the satellite constellation altitude;
 1. the number of orbits;
 1. the number of satellites per orbit;
 1. the orbits inclination;
 1. the Walker Delta phasing parameter;
 1. the Inter-Satellite Links (ISLs) data rate;
 1. the ISLs optical exclusion zone half-cone angle;
 1. the Sun latitude and longitude. 

The simulator builds a Walker Delta constellation with circular orbits, where satellites have uniform angular separation along the orbit.
ISLs are assumed to be bidirectional. If one of the two links of a bidirectional pair is in outage at the physical layer, the link is assumed to be unavailable in both directions.

## Software structure
The simulator is a single-file Matlab program. It should work with plain Matlab - no toolboxes required.


Thanks for your interest in SpecaDataLEO!  
Alberto  
albertogp71.github@gmail.com


