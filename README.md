# 2021 JGI Internship Project

### gPPQ Calculator

gPPQ_calculator.py is a python script for calculating the [PPQ and gPPQ scores from Pfeifer et.al. 2021](https://academic.oup.com/nar/article/49/5/2655/6137301).




The pipeline can be summarized with 4 steps: 




1. Predict open reading frames with Prodigal
2. Create a reference database for diamond
3. Perform local alignments using Diamond to generate hits to Phages and or Plasmids
    - **Hits:**
        - e-value < 10^-4
        - identity >= 35%
        - coverage <= 50%
4. Calculate average phage-plasmid quotient (gPPQ) scores
    1. Calculate PPQ score for each genome with 10 <= protein sequences
        - How are the scores are calculated?
          - Phage Hits normalized by the length of the Phage database --> H(Phage) 
          - Plasmid Hits normalized by the length of the Plasmid database --> H(Plasmid)
          - Take normalized phage hits divide by the sum of normalized phage hits and normalized plasmid hits
    2. Calculate average PPQ score per genome.

**Phage-plasmid quotient (PPQ):**

<img src="https://render.githubusercontent.com/render/math?math=\large PPQ = \frac{\frac{H(Phages)}{length(Phage \, Database)}}{\frac{H(Phages)}{length(Phage \, Database)} %2B \frac{H(Plasmids)}{length(Plasmid \, Database)}}">


**gPPQ:**

<img src="https://render.githubusercontent.com/render/math?math=\large gPPQ = \frac{\sum Phage \, PPQs}{\sum Phage \, PPQs \> %2B \sum Plasmid \, PPQs}">


