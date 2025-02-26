# nutrients-fish-Pocillopora
R code and data to reproduce all analyses and plots found in the manuscript: Burgess SC, Johnston EC, Speare KE, McLachlan RH, Adam TC, Vega Thurber, R, Burkepile DE. Differential effects of nutrients and fish consumer pressure on sympatric cryptic coral species (Pocillopora spp.). Ecology, in press.

GENERAL INFORMATION

1. Title of Dataset: Differential effects of nutrients and fish consumer pressure on sympatric cryptic coral species (Pocillopora spp.)

2. Author Information
	A. Principal Investigator Contact Information
		Name: Scott Burgess
		Institution: Florida State University
		Address: 319 Stadium Drive, Tallahassee, FL, USA 32306
		Email: sburgess@bio.fsu.edu

3. Date of data collection (single date, range, approximate date): 2019-2021 

4. Geographic location of data collection: Moorea, French Polynesia

5. Information about funding sources that supported the collection of the data: National Science Foundation (NSF; OCE-1829867, OCE-2023701, OCE-2023424, OCE-1637396)



DATA & FILE OVERVIEW

Colony size over time.csv  
macroalgae_T10.csv  
Predation and Polyp data.csv  
Recruit mortality.csv  
Recruitment timepoint and size T10.csv  
Species IDs.csv  
Figure 2 and 3.R  
Figure 4 and 5.R  
Figure 6.R  
Figure S1.R  

2. Relationship between files:   
Figure 2 and 3.R uses Species IDs.csv, Recruitment timepoint and size T10.csv, and macroalgae_T10.csv  
Figure 4 and 5.R uses Species IDs.csv, Recruitment timepoint and size T10.csv, and Colony size over time.csv  
Figure 6.R uses Species IDs.csv and Predation and Polyp data.csv  

3. Metadata  

Colony size over time.csv  
Un_ID: Unique ID for each coral  
Timepoint: The timepoint of data collection. T0 = July/Aug 2018, T3 = July/Aug 2019, T6 = July/Aug 2020, T9=July/Aug 2021.   
Block_plot: The plot in which the coral is located. Plots are A3, A4, B1, B2, C2, C4, D1, D3. Each plot has 4 exclosures, each with a different level of herbivory  
Herbivory_trt: This is the consumer pressure treatment. Treatments are 1X1, 2x2, 3X3, and open, which correspond to very low, low, medium, and high levels of herbivory and corallivory  
Coral_Id: The correct coral ID's for all corals. Use these. These are the ID's used for the "Un_ID" column. These reflect all corrections. For most corals "Coral_ID" and "Coral_ID_past" are exactly the same. They only differ if we made corrections in the "Coral_ID_field" or "Coral_ID_corrected" columns.  
Coral_Genus: Genus of the coral	  
Enrichment_trt: This is the enrichment treatment. Plots are either Enriched or Ambient. Plots that are Enriched have nutrients added with fertilizer-filled nutrient diffusers.  
Size: The longest diameter of the coral (cm), measured in situ with a ruler.  
Recruitment: 1 = recruited, 0 = did not recruit.  
Alive: 1 = alive, 0 = dead  
Recruitment_timepoint: The first timepoint that we saw this coral. T0 = July/Aug 2018, T3 = July/Aug 2019, T6 = July/Aug 2020, T9=July/Aug 2021.   


macroalgae_T10.csv  
Block_Plot: The plot in which the coral is located. Plots are A3, A4, B1, B2, C2, C4, D1, D3. Each plot has 4 exclosures, each with a different level of herbivory  
Block_Plot_Herb_Trt: Combindation of Block_Plot and Herbivory_trt  
Herbivory_trt: This is the consumer pressure treatment. Treatments are 1X1, 2x2, 3X3, and open, which correspond to very low, low, medium, and high levels of herbivory and corallivory  
Enrichment_trt: This is the enrichment treatment. Plots are either Enriched or Ambient. Plots that are Enriched have nutrients added with fertilizer-filled nutrient diffusers.  
mean_macroalgae: The mean macroalgae percent cover between 2019 and 2021   
median_macroalgae: The median macroalgae percent cover between 2019 and 2021  

Predation and Polyp data.csv  
Un_ID: Unique ID for each coral  
percent_polyps_extended_Nov21_T10: Percentage of the colony surface area on which polyps were extended (in 5% increments ranging from 0% to 100%) in Nov 2021  
Bites_present_Nov21_T10: 1 = recent bite scars from excavating corallivores on the coral surface in Nov 2021, 0 = recent bite scars absent in Nov 2021  
percent_polyps_extended_Apr22_T11: Percentage of the colony surface area on which polyps were extended (in 5% increments ranging from 0% to 100%) in Apr 2022  
Bites_present_Apr22_T11: 1 = recent bite scars from excavating corallivores on the coral surface in Apr 2022, 0 = recent bite scars absent in Apr 2022  
percent_polyps_extended_Jul22_T12: Percentage of the colony surface area on which polyps were extended (in 5% increments ranging from 0% to 100%) in Jul 2022  
Bites_present_Jul22_T121 = recent bite scars from excavating corallivores on the coral surface in Jul 2022, 0 = recent bite scars absent in Jul 2022  

Recruit mortality.csv  
Un_ID: Unique ID for each coral  
Timepoint: The timepoint of data collection. T0 = July/Aug 2018, T3 = July/Aug 2019, T6 = July/Aug 2020, T9=July/Aug 2021.   
Un_id_tp: A column that combines Un_id and Timepoint  
Block_Plot: The plot in which the coral is located. Plots are A3, A4, B1, B2, C2, C4, D1, D3. Each plot has 4 exclosures, each with a different level of herbivory  
Herbivory_trt: This is the consumer pressure treatment. Treatments are 1X1, 2x2, 3X3, and open, which correspond to very low, low, medium, and high levels of herbivory and corallivory  
Coral_Id: The correct coral ID's for all corals. Use these. These are the ID's used for the "Un_ID" column. These reflect all corrections. For most corals "Coral_ID" and "Coral_ID_past" are exactly the same. They only differ if we made corrections in the "Coral_ID_field" or "Coral_ID_corrected" columns.  	
Coral_Genus: Genus of the coral	  
Enrichment_trt: This is the enrichment treatment. Plots are either Enriched or Ambient. Plots that are Enriched have nutrients added with fertilizer-filled nutrient diffusers.  
Recruitment  
Alive: 1 = alive, 0 = dead  
Notes1: Notes  
Recruitment_timepoint: The first timepoint that we saw this coral. T0 = July/Aug 2018, T3 = July/Aug 2019, T6 = July/Aug 2020, T9=July/Aug 2021.   

Recruitment timepoint and size T10.csv  
Un_ID: Unique ID for each coral  
Plot: The plot in which the coral is located. Plots are A3, A4, B1, B2, C2, C4, D1, D3. Each plot has 4 exclosures, each with a different level of herbivory  
Herbivory_trt: This is the consumer pressure treatment. Treatments are 1X1, 2x2, 3X3, and open, which correspond to very low, low, medium, and high levels of herbivory and corallivory  
Enrichment_trt: This is the enrichment treatment. Plots are either Enriched or Ambient. Plots that are Enriched have nutrients added with fertilizer-filled nutrient diffusers.  
Coral_Id: The correct coral ID's for all corals. Use these. These are the ID's used for the "Un_ID" column. These reflect all corrections. For most corals "Coral_ID" and "Coral_ID_past" are exactly the same. They only differ if we made corrections in the "Coral_ID_field" or "Coral_ID_corrected" columns.  
Size_T10: The longest diameter of the coral (cm), measured in situ with a ruler in Nov 2021  
Recruitment_timepoint: The first timepoint that we saw this coral. T0 = July/Aug 2018, T3 = July/Aug 2019, T6 = July/Aug 2020, T9=July/Aug 2021.   

Species IDs.csv  
Un_ID: Unique ID for each coral  
Plot: The plot in which the coral is located. Plots are A3, A4, B1, B2, C2, C4, D1, D3. Each plot has 4 exclosures, each with a different level of herbivory  
Herbivory_trt: This is the consumer pressure treatment. Treatments are 1X1, 2x2, 3X3, and open, which correspond to very low, low, medium, and high levels of herbivory and corallivory  
Enrichment_trt: This is the enrichment treatment. Plots are either Enriched or Ambient. Plots that are Enriched have nutrients added with fertilizer-filled nutrient diffusers.  
Coral_Id: The correct coral ID's for all corals. Use these. These are the ID's used for the "Un_ID" column. These reflect all corrections. For most corals "Coral_ID" and "Coral_ID_past" are exactly the same. They only differ if we made corrections in the "Coral_ID_field" or "Coral_ID_corrected" columns.  
Haplotype: the mtORF haplotype from Sanger sequences  
Species: The species name   
flag_uncertain_recruit_time: 1 = Recruitment_timepoint was uncertain. 0 = Recruitment_timepoint is known, can use  
flag_partial_mortality: 1 = Partial mortality. 0 = No partial mortality   
flag_unlabelled: 1 = "Poc in bags" samples  sampled from the new tiny recruits in plots C and D, but which  the datasheets were temporarily misplaced which allow us to figure out which samples are which.   
flag_bad_sequence: 1 = sample was sequenced but the mtORF sequence was poor quality, so no Species given; 0 = mtORF sequence was good and a Species name is provided  
flag_missing_sample: 1 = no tube or sample was provided; 0 = sample was provided    

