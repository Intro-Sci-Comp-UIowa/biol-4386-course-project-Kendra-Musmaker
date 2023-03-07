# Intro
Telomeres are repetitive sequences bound by proteins at the ends of chromosomes. They shorten with each round of replication due to the “end replication problem,” the phenomenon in which Okazaki fragments cannot be placed at the very end of a replicating chromosome. This telomere shortening continues until the reach a critically short length. The shortening of the telomere to this critical point results in it appearing to the cell as a one-sided double strand break. The DNA damage response is then activated. This then usually results in cellular senescence. A subset of cells, however, will circumvent cell death by lengthening their telomeres. This often happens via the reactivation of telomerase. However, telomeres can also be lengthened via a mechanism known as alternative lengthening of telomeres (ALT) in which a critically short telomere uses the end of another chromosome as a template for synthesis. Cells that are able to bypass cell death are then referred to as “survivors”. Previous work has shown there to be two “types” of survivors in regard to chromosome end structure after undergoing ALT. Genes involved in ALT have historically been determined by deleting candidate genes to see if ALT can still occur, and if so, what type of survivor is formed. However, until recently researchers have not been able to determine a frequency for ALT. In this paper, a method for determining ALT frequency was developed. Once overall frequency was determined, genes suspected of being involved in ALT were deleted and changes in ALT frequency determined. Overall frequency was determined by deleting a component of telomerase. 
# Figure to Reproduce

The figure I want to recreate from Kockler et al 2022 is figure 1c, which shows the frequency of ALT in the context of various gene deletions. This figure shows the overall frequency of ALT to be 2x10-5. These frequencies were determined by deleting TLC1, the RNA component of telomerase. As the telomerase complex is nonfunctional without TLC1, passaging yeast to senescence in the absence of TLC1 therefore allows for the determination of ALT frequency. By deleting TLC1 alone, the overall frequency of ALT was determined to be 2x10^-5. Other genes believed to be involved in ALT were also deleted, and changes in frequency determined in these genetic backgrounds. It was shown that the deletion RAD51 or RAD59 decreases ALT frequency. This was expected, as it has been shown qualitatively that loss of Rad51 or Rad59 each eliminate a different type of ALT event. Cell cycle arrest is needed for ALT precursors to have sufficient time to become survivors, and so it is unsurprising  that loss of checkpoint proteins RAD9 and RAD24 all decrease ALT frequency. 
A major step of RAD51-dependent ALT is the invasion of a telomere end into another telomere end, facilitated by Rad51. Once the short telomere end invades into a template, however, Rad51 needs to be removed from the invaded DNA. Tis is done by SRS2. If SRS2 is deleted, Rad51 is never removed from the invading DNA, and something known as a toxic joint molecule is formed. This results in cell death and entirely eliminated. This can be partially rescued by reducing the amount of Rad51 loaded onto the short telomeric end, which is done by deleting a protein involved loading Rad51 onto DNA, Rad55. With Rad55 and Srs2 deleted in TLC`-deficient cells results in a frequency of ALT that is lower than the overall ALT frequency but higher than when SRS2 deletion alone (with TLC1 deletion). The decrease in ALT frequency in the absence of these proteins indicates that these are proteins that normally promote ALT. 
 The increase in ALT frequency when telomere capping proteins RIF1 and RIF2 are deleted indicates these proteins normally suppress ALT. This is a logical result, as the capping of telomeres prevents them from being recognized by ALT machinery. 
# Materials and Methods
- Strains
In this figure, the strain AM3692 and its derivatives are used. AM3692 is a diploid strain that is heterozygous for a deletion of the telomerase component TLC1. The yeast are maintained as diploids to prevent senescence.
- Protocols
To begin a frequency experiment the yeast undergo meiosis and the resulting spores (daughter cells) are separated from one another. TLC1-deleted cells are then select>
Rather than recreating this figure, I am interested in the effect of deleting other genes thought classically to be involved in the ALT process, such as RAD54, which assists in Rad51 function, and >
The next step is to use pandas within python to run anumber of calculations and add them to this csv. I first need to determine cell/ml(x10^6) (y.passage#) by dividing cells counted by area and con>
- Reagents
# Resullts
The file Day2_3data_counts_mlh1 contains data on cells counted, number of rectangles the cells were counted in and the concentration at which they were counted for day2 and day3 of the frequenciy experiments. There is also a column with the number of survivors counted. The data was acquired across 5 frequency experiments run late 2021 to early 2023. The day is formatted so that there is a column for culture name, 3 columns for day 2 cell, rectangle and conc (concentration), 3 columns for day 3 cell, rectangle and conc, and a final column for survivors. single mutant tlc1 and double mutant tlc1mlh1 data are in the same sheet and columns. May reformat the sheet at a later date.
# Discussion
# Conclusion
