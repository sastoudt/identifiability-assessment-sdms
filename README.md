## identifiability in species distribution and abundance models README

check_reproducibility.R – reproducibility checks using fertile package

figure1.R – all code needed to produce Figure 1 in the manuscript

install_proj_packages.R – code to install all required packages needed to run the code in this repository, recommended to run as the first step, this file also notes package and R versions

MetadataS1.docx – description of files in supplementary code

modified_source_code
	svisitFormula.R – tweaked code for svisitFormula function found in the detect package for extra tracking
	svocc.fit.R - tweaked code for svocc.fit function found in the detect package for extra tracking
	svocc.R - tweaked code for svocc function found in the detect package for extra tracking


overall_code.R- high level script that calls helper scripts to execute all of the analysis in the manuscript and supplement and make the remaining plots, run after installing the requisite packages

presence_absence_only
	prevalence_large.R – code to fit presence-only and presence-absence models for large sample size, make corresponding plots
	
	prevalence_medium.R - code to fit presence-only and presence-absence models for medium sample size, make corresponding plots
	
	prevalence_small.R - code to fit presence-only and presence-absence models for small sample size, make corresponding plots
	
	prevalence_smallNoSpline.R- code to fit presence-only and presence-absence models for small sample size without spline, make corresponding plots
	
	prevalenceExampleDataSetup.R – get simulated presence-only and presence-absence data to fit models to
	
	raw_data
		
		popaData_large.RData – large presence-only and presence-absence data generated in prevalenceExampleDataSetup.R 
		
		popaData_medium.RData – medium presence-only and presence-absence data generated in prevalenceExampleDataSetup.R
		
		popaData_small.RData -  small presence-only and presence-absence data generated in prevalenceExampleDataSetup.R
	
	results
		
		paResults_small.RData –  presence-absence results for small data without splines
		
		paSplineResults_largeMoreKnots.RData –  presence-absence results for large data (not stored on GitHub because too large, generate yourself using the provided code)
		
		paSplineResults_mediumMoreKnots.RData - presence-absence results for medium data
		
		paSplineResults_smallMoreKnots.RData - presence-absence results for small data with splines
		
		poResults_small.RData - presence-only results for small data without splines
		
		poSplineResults_largeMoreKnots.RData - presence-only results for large data (not stored on GitHub because too large, generate yourself using the provided code)
		
		poSplineResults_mediumMoreKnots.RData - presence-only results for medium data
		
		poSplineResults_smallMoreKnots.RData - presence-only results for small data with splines
	
single_double_visit

	occurDetect_large.R - code to fit single-visit and double-visit models for large sample size, make corresponding plots
	
	occurDetect_medium.R - code to fit single-visit and double-visit models for medium sample size, make corresponding plots
	
	occurDetect_small.R - code to fit single-visit and double-visit models for small sample size, make corresponding plots
	
	occurDetect_small_noSpline.R - code to fit single-visit and double-visit models for small sample size without spline, make corresponding plots
	
	occurDetectDataSetup.R - get simulated single-visit and double-visit data to fit models to
	
	raw_data
		
		svdvOccurData_largeNoSaturateLine.RData - large single-visit and double-visit data generated in occurDetectDataSetup.R (not stored on GitHub because too large, generate yourself using the provided code)
		
		svdvOccurData_mediumNoSaturateLine.RData- medium single-visit and double-visit data generated in occurDetectDataSetup.R
		
		svdvOccurData_smallNoSaturateLine.RData- small single-visit and double-visit data generated in occurDetectDataSetup.R
	
	results
		
		dvOccurResults_smallNoSaturateLine.RData – double-visit results for small data without splines
		
		dvOccurSplineResults_largeNoSaturateLine.RData - double-visit results for large data  (not stored on GitHub because too large, generate yourself using the provided code)
		
		dvOccurSplineResults_mediumNoSaturateLine.RData - double-visit results for medium data
		
		dvOccurSplineResults_smallNoSaturateLine.RData – double-visit results for small data with splines
		
		svOccurResults_smallNoSaturateLine.RData - single-visit results for small data without splines
		
		svOccurSplineResults_largeNoSaturateLine.RData - single-visit results for large data  (not stored on GitHub because too large, generate yourself using the provided code)
		
		svOccurSplineResults_mediumNoSaturateLine.RData – single-visit results for medium data
		
		svOccurSplineResults_smallNoSaturateLine.RData - single-visit results for small data with splines

single_double_visit_abundance
	
	abunDetect_large.R - code to fit single-visit and double-visit abundance models for large sample size, make corresponding plots
	
	abunDetect_medium.R - code to fit single-visit and double-visit abundance models for medium sample size, make corresponding plots
	
	abunDetect_small.R - code to fit single-visit and double-visit abundance models for small sample size, make corresponding plots
	
	abunDetect_small_noSpline.R - code to fit single-visit and double-visit abundance models for small sample size without spline, make corresponding plots
	
	abunDetectDataSetup.R - get simulated single-visit and double-visit abundance data to fit models to

raw_data
	
	abunData_large.RData - large single-visit and double-visit abundance data generated in abunDetectDataSetup.R (not stored on GitHub because too large, generate yourself using the provided code)
	
	abunData_medium.RData - medium single-visit and double-visit abundance data generated in abunDetectDataSetup.R
	
	abunData_small.RData - small single-visit and double-visit abundance data generated in abunDetectDataSetup.R

results
	
	dvAbundanceLarge.RData – double-visit results for large data (not stored on GitHub because too large, generate yourself using the provided code)
	
	dvAbundanceMedium.RData - double-visit results for medium data
	
	dvAbundanceSmall_noSpline.RData - double-visit results for small data without splines
	
	dvAbundanceSmall.RData - double-visit results for small data with splines
	
	svAbundanceLarge.RData - single-visit results for large data (not stored on GitHub because too large, generate yourself using the provided code)
	
	svAbundanceMedium.RData - single-visit results for medium data
	
	svAbundanceSmall_noSpline.RData - single-visit results for small data without splines
	
	svAbundanceSmall.RData - single-visit results for small data with splines

