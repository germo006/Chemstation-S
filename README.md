# Chemstation-S
The codes corresponding to the field data used in my PhD thesis will live here. 
The title of the repository refers both to what this data is and where it came from. I sampled metabolites through benzoyl chloride derivatization on [BIOS-SCOPE](https://scope.bios.asu.edu/) cruise AE2123 (Nov. 2021) at Hydrostation S, which you can read about [here](https://bios.asu.edu/research/projects/hydrostation-s/). These data were taken over two days, with depth-resolved sampling via CTD rosette every six hours. As for the exact results of all this and what those results mean, you will eventually be able to read [the paper](endless.horse). 

## Who might use this repository, and how?
For now, if you have landed here I assume you are either (1) me, Noah Germolus or (2) interested in the analysis codes I used. If you are trying to process similar LC-MS data, starting from peak areas and using the same derivatization method that I used, feel free to contact me (germolus@mit.edu) with questions and/or start with the [SkyMat package](https://github.com/KujawinskiLaboratory/SkyMat), which the Kujawinski Lab uses to turn the instrument data into concentrations. If you do really want the meat of the processing codes post-calibration, this is the place for those. It's a bit complicated, though: some of the process requires data from my two other chapters. Once they are publicly available, you should be able to download spreadsheets that contain the relevant chemical data and metadata, but for now you'd have to follow along with the whole process, unless you *just* want metabolite profiles.
At the very least, you will need the data for this study, which can be found [here](https://www.windows93.net/).

*Photochemistry*
* [Code Repository](https://github.com/germo006/metabolitephotochem)
* [Data Repository](https://www.ebi.ac.uk/metabolights/editor/study/MTBLS7513/descriptors)

*Zooplankton Excretion*
* [Code Repository](https://github.com/germo006/zoopee)
* [Data Repository](https://www.ebi.ac.uk/metabolights/editor/study/MTBLS9061/descriptors)

In addition (I know, this is annoying, but trust me, it makes it feel more fun) the graph-generation codes rely upon the colormaps I made [here](https://github.com/germo006/NoahMaps). Use them however you like and/or listen to the albums they're based on. 

## You have all of the stuff downloaded. Now what?
I dunno, I guess you can redo all my graphs in different colors, or apply slightly different error bounds or something. The world is your oyster. Use this code wisely. 
