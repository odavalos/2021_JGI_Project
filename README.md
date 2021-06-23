# 2021 JGI Internship Project




### gPPQ Calculator

The purpose of this script is to calculate the [PPQ and gPPQ scores from Pfeifer et.al. 2021](https://academic.oup.com/nar/article/49/5/2655/6137301).


What is considered a hit?

**Hit:**
- e-value < 10^-4
- identity >= 35%
- coverage <= 50%


**Phage-plasmid quotient (PPQ):**

<img src="https://render.githubusercontent.com/render/math?math=\large PPQ = \frac{H(phages)}{H(phages) %2B H(plasmids)}">


**gPPQ:**

Average of PPQ scores per P-P (Phage-Plasmid).
- P-P >= 10 protein sequences used

<img src="https://render.githubusercontent.com/render/math?math=\large gPPQ = \frac{\ %23 Phage \, PPQ \, Matches}{\ %23 \, Phage \, PPQ \, Matches \> %2B \ %23 \,Plasmid \, PPQ \, Matches}">


