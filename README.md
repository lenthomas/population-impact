# population-impact

R code to simulate the population impact of cumulative anthropogenic stressors on wildlife populations (operating via prey reduction).

Methods are documented in the manuscript: 
Williams, R., L. Thomas, E. Ashe, C.W. Clark and P.S. Hammond. Gauging allowable harm limits to cumulative, sub-lethal effects of human activities on wildlife populations.

File manifest:

1. killerSim.r - code to implement the killer whale example.  
Auxiliary input files:

  * CTCwcvi.xlsx
  * cord chinook cpue.xlsx
  * female_survival_restricted.csv
  * femaile_survival_unrestricted.csv
  * male_survival_restricted.csv
  * maile_survival_unrestricted.csv
  
2. humpbackSim.r - code to implement the humpback whale example.

The simulation .r files are very similar -- they differ in the details of the population biology in the simulation.  The intention is that these can act as templates for simulations on other biological systems of interest.

