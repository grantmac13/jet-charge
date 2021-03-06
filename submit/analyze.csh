#!/bin/csh
#Isaac Mooney, WSU - June 2019
#This script will expand as I add in more of the analyses: pp, pA, AA
#It currently handles the QA, data (for all three systems), and the Pythia6 and Pythia6+Geant simulation on disk. It is being expanded to handle embedded pAu.

#args:
#1: the analysis type. Current options: QA, data, sim, embed
#2: the trigger. Current options for pp: JP2, HT2, VPDMB(-nobsmd); for pA: JP2, HT2, BBCMB. For MC analysis, only JP2 is allowed at the moment.
#3: the species. Current options: pp, pA, AA
#4: the jet radius. Options: any number > 0 & <= 9.9; passed in without the decimal, e.g. "06" for 0.6.
#5: an analysis-specific wildcard (not required). Currently being used in data to take either full (full) or charged (ch) jets.
#6: second analysis-specific wildcard (not required). Currently being used in sim to use data with or without decays at particle level

# changed
#6: second analysis-specific wildcard (not required). Currently being used in sim to either require matching (matched) - for responses - or not (unmatched) - for normal spectra

# e.g. csh submit/analyze.csh data JP2 pp 04 full
# or   csh submit/analyze.csh sim JP2 pp 04 full matched

# Arguments       
if ( $# < "4") then
    echo 'Error: illegal number of parameters'
    exit
endif
    
set ExecPath = `pwd`
set arg = '' 

set analysisType = $1
set trigger = $2
set species = $3
set radius = $4
set wc1 = $5
set wc2 = $6

#~~~input checking and directory / file setup~~~#

if (${analysisType} == 'QA') then
    make bin/QA || exit
    set execute = './bin/QA'
    # Create the folder name for output
    set outFile = QA
    echo 'Running in QA mode!'
else if (${analysisType} == 'data') then
    make bin/data || exit
    set execute = './bin/data'
    # Create the folder name for output
    set outFile = data
    echo 'Running in data mode!'
else if (${analysisType} == 'sim') then
    make bin/sim || exit
    set execute = './bin/sim'
    # Create the folder name for output
    set outFile = sim
    echo 'Running in sim mode!'
else if (${analysisType} == 'decay') then
    make bin/decay || exit
    set execute = './bin/decay'
    # Create the folder name for output
    set outFile = decay
    echo 'Running on decay data!'
else if (${analysisType} == 'embed') then
    make bin/embed || exit
    set execute = './bin/embed'
    # Create the folder name for output
    set outFile = embed
    echo 'Running in embed mode!'
else
    echo 'Error: unrecognized analysis type option'
    exit
endif

if (${analysisType} == 'QA' || ${analysisType} == 'data') then
    if (${trigger} == 'JP2' && ${species} == 'pp') then
	set trigger = "ppJP2"
	set base = /tier2/home/groups/rhi/STAR/Data/ppJP2Run12/sum	
	echo "Running on the ppJP2-triggered data!"
    else if (${trigger} == 'HT2' && ${species} == 'pp') then
	set trigger = "ppHT2"
	set base = /tier2/home/groups/rhi/STAR/Data/ppHT2Run12/pp12Pico
	echo "Running on the ppHT2-triggered data!"
    else if (${trigger} == 'VPDMB' && ${species} == 'pp') then
	set trigger = "ppVPDMB"
	set base = /tier2/home/groups/rhi/STAR/Data/ppMBRun12/sum
	echo "Running on the ppMB-triggered data!"
    else if (${trigger} == 'JP2' && ${species} == 'pA') then
	set trigger = "pAuJP2"
	set base = /tier2/home/groups/rhi/STAR/Data/P16id/production_pAu200_2015/HT/pAu_2015_200_HT_1
	echo "Running on the pAuJP2-triggered data!"
    else if (${trigger} == 'HT2' && ${species} == 'pA') then
	set trigger = "pAuHT2"
	set base = /tier2/home/groups/rhi/STAR/Data/P16id/production_pAu200_2015/HT/pAu_2015_200_HT_1
	echo "Running on the pAuHT2-triggered data!"	
    else if (${trigger} == 'BBCMB' && ${species} == 'pA') then
	set trigger = "pAuBBCMB"
	set base = /tier2/home/groups/rhi/STAR/Data/P16id/production_pAu200_2015_oldTrees/MB/pAu_2015_200_MB_1
	echo "Running on the pAuBBCMB-triggered data!"
    else if (${species} == 'AA' || ${species} == 'AuAu') then
	echo "AA is not ready yet - be patient!"
    else
	echo "Error: unrecognized trigger and/or species option!"
	exit
    endif
else if (${analysisType} == 'decay') then
#    set base = /tier2/home/groups/rhi/STAR/Data/AddedEmbedPythiaRun12pp200/Cleanpp12Pico_
#    set base = /tier2/home/groups/rhi/STAR/Data/RhicPythia6/PYTHIA6_decay_1_startune_1/pythia6_STARTune_
#    set base = /tier2/home/groups/rhi/STAR/Data/IsaacsPy6_hadds/decay_*_startune_*/decay_
    set base = /tier2/home/groups/rhi/STAR/Data/IsaacsPy6_hadds/decay_1_startune_1/decay_
    echo "Running on the Pythia6 and Pythia6+Geant (JP2-triggered) simulation with decays!"
else if (${analysisType} == 'sim') then
    set base = /tier2/home/groups/rhi/STAR/Data/AddedEmbedPythiaRun12pp200/Cleanpp12Pico_
    echo "Running on the Pythia6 and Pythia6+Geant (JP2-triggered) simulation!"
else
#    set base = /tier2/home/groups/rhi/STAR/Data/AddedEmbedPythiaRun12pp200/Cleanpp12Pico_pt11_15_g1
    #set base = /tier2/home/groups/rhi/STAR/Data/temp_embed_sample/pt-hat
    set base = /tier2/home/groups/rhi/STAR/Data/EmbedPythiaRun15pAu200_picos/P18ih/pAu_200_production_2015/out/HT2/pt-hat1525_
    echo "Running on the HT2-triggered pAu embedding!"
endif

# Make output directories if they don't already exist            
if ( ! -d out/${outFile} ) then
mkdir -p out/${outFile}
endif

if ( ! -d log/${outFile} ) then
mkdir -p log/${outFile}
endif

echo ${base} #for debugging

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Now Submit jobs for each data file
foreach input ( ${base}* )

    #synthesize the output file base name from the input file base name and the options passed
    set OutBase = `basename $input | sed 's/.root//g'`
    set uscore = "_" #useful for chaining variables together
    if (${analysisType} == 'QA' || ${analysisType} == 'data') then
	set OutBase = "$OutBase$uscore$trigger$uscore${wc1}${uscore}R${radius}"
    else if (${analysisType} == 'sim') then 
	set OutBase = "TEST_$OutBase$uscore${wc1}$uscore${wc2}${uscore}R${radius}"
    else
    set OutBase = "$OutBase${uscore}${wc2}${uscore}R${radius}"
#	set OutBase = "TESTSAMPLE_$OutBase${uscore}${wc2}${uscore}R${radius}"
    endif
    #make the output path and names
    set outLocation = out/${outFile}/
    set outName = ${OutBase}.root

    #input files
    set Files = ${input}
    
    #log files
    set LogFile = log/${outFile}/${OutBase}.log
    set ErrFile = log/${outFile}/${OutBase}.err

    echo "Logging output to " $LogFile
    echo "Logging errors to " $ErrFile

    set dummy = 'dummy' #this extra argument is for compatibility between QA and the other segments of the analysis (or for future expansion)

    #setting command-line args based on the analysis file being used:
    if (${analysisType} == 'QA') then
	set arg = "$outLocation $outName $trigger $dummy $Files"
    else if (${analysisType} == 'data') then
	set arg = "$outLocation $outName ${radius} $trigger ${wc1} $Files"
    else if (${analysisType} == 'sim' || ${analysisType} == 'embed' || ${analysisType} == 'decay') then
	set arg = "$outLocation $outName ${radius} ${wc1} ${wc2} $Files"
    endif

    echo "now submitting this script: "
#    echo qsub -V -q erhiq -l mem=4GB -o $LogFile -e $ErrFile -N ${analysisType} -- ${ExecPath}/submit/qwrap.sh ${ExecPath} $execute $arg

#    qsub -V -q erhiq -l mem=4GB -o $LogFile -e $ErrFile -N ${analysisType} -- ${ExecPath}/submit/qwrap.sh ${ExecPath} $execute $arg   
    echo sbatch -p erhip -q express --mem-per-cpu=32GB -o $LogFile e $ErrFile --job-name=${analysisType} -- ${ExecPath}/submit/qwrap.sh ${ExecPath} $execute $arg

    sbatch -p erhip -q express --mem-per-cpu=32GB -o $LogFile -e $ErrFile -t 120 --job-name=${analysisType} -- ${ExecPath}/submit/qwrap.sh ${ExecPath} $execute $arg   

end
