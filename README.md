# Alzheimer's PET analysis

Normalized microarray and RNA sequencing data (http://human.brain-map.org/static/download) from the Allen Human Brain Atlas project must first be downloaded to conduct analysis. These files are to be placed in /data/raw/allen_hba (folder for each donor).

The first step is to run the Combine MRI maps.ipynb notebook to create the combined sample annotation file. This needs the nilearn, numpy, nibabel, and pandas python (3.6) libraries.

Then the numbered R scripts can be ran in order. The first file will help load the required R libraries. However, this is not a complete listing of the packages needed and will potentally install newer versions than what we used. 

This repo contains or reuses data from several publications:

Grothe, Michel J., Jorge Sepulcre, Gabriel Gonzalez-Escamilla, Irina Jelistratova, Michael Schöll, Oskar Hansson, Stefan J. Teipel, and Alzheimer’s Disease Neuroimaging Initiative. 2018. “Molecular Properties Underlying Regional Vulnerability to Alzheimer’s Disease Pathology.” Brain: A Journal of Neurology 141 (9): 2755–71.

Crow, Megan, Nathaniel Lim, Sara Ballouz, Paul Pavlidis, and Jesse Gillis. 2019. “Predictability of Human Differential Gene Expression.” Proceedings of the National Academy of Sciences of the United States of America 116 (13): 6491–6500.

Darmanis, Spyros, Steven A. Sloan, Ye Zhang, Martin Enge, Christine Caneda, Lawrence M. Shuer, Melanie G. Hayden Gephart, Ben A. Barres, and Stephen R. Quake. 2015. “A Survey of Human Brain Transcriptome Diversity at the Single Cell Level.” Proceedings of the National Academy of Sciences of the United States of America 112 (23): 7285–90.

Hawrylycz, Michael J., Ed S. Lein, Angela L. Guillozet-Bongaarts, Elaine H. Shen, Lydia Ng, Jeremy A. Miller, Louie N. van de Lagemaat, et al. 2012. “An Anatomically Comprehensive Atlas of the Adult Human Brain Transcriptome.” Nature 489 (7416): 391–99.

Jansen, Iris E., Jeanne E. Savage, Kyoko Watanabe, Julien Bryois, Dylan M. Williams, Stacy Steinberg, Julia Sealock, et al. 2019. “Genome-Wide Meta-Analysis Identifies New Loci and Functional Pathways Influencing Alzheimer’s Disease Risk.” Nature Genetics 51 (3): 404–13.

Keren-Shaul, Hadas, Amit Spinrad, Assaf Weiner, Orit Matcovitch-Natan, Raz Dvir-Szternfeld, Tyler K. Ulland, Eyal David, et al. 2017. “A Unique Microglia Type Associated with Restricting Development of Alzheimer’s Disease.” Cell. https://doi.org/10.1016/j.cell.2017.05.018.

Mathys, Hansruedi, Jose Davila-Velderrain, Zhuyu Peng, Fan Gao, Shahin Mohammadi, Jennie Z. Young, Madhvi Menon, et al. 2019. “Single-Cell Transcriptomic Analysis of Alzheimer’s Disease.” Nature 570 (7761): 332–37.

Miller, Jeremy A., Vilas Menon, Jeff Goldy, Ajamete Kaykas, Chang-Kyu Lee, Kimberly A. Smith, Elaine H. Shen, John W. Phillips, Ed S. Lein, and Mike J. Hawrylycz. 2014. “Improving Reliability and Absolute Quantification of Human Brain Microarray Data by Filtering and Scaling Probes Using RNA-Seq.” BMC Genomics 15 (February): 154.

Ogan, M. B., Lilah Toker, Shreejoy J. Tripathy, Brenna Li, Brad Rocco, Etienne Sibille, and Paul Pavlidis. 2017. “Cross-Laboratory Analysis of Brain Cell Type Transcriptomes with Applications to Interpretation of Bulk Tissue Data.” eNeuro, November, ENEURO.0212–17.2017.

Zeisel, Amit, Hannah Hochgerner, Peter Lönnerberg, Anna Johnsson, Fatima Memic, Job van der Zwan, Martin Häring, et al. 2018. “Molecular Architecture of the Mouse Nervous System.” Cell 174 (4): 999–1014.e22.
