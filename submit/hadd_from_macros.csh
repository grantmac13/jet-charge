#!/bin/csh
# Grant McNamara with help from Isaac Mooney -- February 2022

# command line arguments
    #1: the executable, e.g. hists (the executable created by passing runroot hists --which calls hists.cxx)
    #2: input type, e.g. 'data', 'sim'
    #3: pp (for now this is the only system I use)
    #4: jet radius, passed as "04" for R = 0.4
    #5: wildcard, e.g. matched? or decayed? -- effectively this can be used as a pythia 6 vs pythia 8 differentiator
    #6: kappa (exponent for jet charge)

if($# < "2") then
    echo 'Error: Wrong number of parameters'
    exit
endif


cd macros # needed because macro_submit called from ppJCandBF, hists.cxx is in macros/
runroot $1 || exit # command not found?
cd ..

set base = ""

set analysisTag = $1 #options are 'hists', 'response', 'stat_err_scaling', 'unfold', 'closure', 'bin_drop', 'compare'
set inputType = $2 # data, sim ------------only hists can be called for data
set species = $3 # pp
set radius = $4 # e.g.'04' for 0.4
set tag = $5 # decayed/undecayed for pythia 8 or matched for pythia 6
set kappa = $6 # pass value of kappa to call hists so that I can create histogram for any kappa needed
# pass kappa the same way as radius

# so try:  csh submit/macro_submit_new.csh hists_new data pp 04 (anything?) 00
# if data, tag ($5) is arbitrary, never checked
# so try:  csh submit/macro_submit_new.csh hists_new sim pp 04 matched 00

echo 'pulling files from out/'${inputType}'!'
echo 'species is '${species}'.'



if(${species} == 'pp') then
    
    if(${inputType} == 'data') then
        set base = out/${inputType}/sum
#        set tag = 
    else if(${inputType} == 'sim') then
        if(${tag} == 'decayed' || ${tag} == 'undecayed') then
            set base = ../for_grant/out/star_mass_pythia8_IsaacsHepmcs_
	#	       ../for_grant/out/star_mass_pythia8_IsaacsHepmcs_25pthatbin30_R04_decayed.root
            # set to the filename from Rivet, move rivet output root files to this directory
        else if(${tag} == 'matched' || ${tag} == 'unmatched') then
            set base = out/${inputType}/TEST_Cleanpp12Pico_pt
    	endif
    endif
endif


set ExecPath = `pwd`

set uscore = "_"
set execute = "./macros/bin/${analysisTag}"


if(${tag} == 'decayed' || ${tag} == 'undecayed') then
    if(! -d out/${tag}/${analysisTag}) then
        mkdir -p out/${tag}/${analysisTag}
    endif
    if(! -d log/${tag}/${analysisTag}) then
        mkdir -p log/${tag}/${analysisTag}
    endif
else if(${analysisTag} == 'unfold' || ${analysisTag} == 'stat_err_scaling' || ${analysisTag} == 'closure' || ${analysisTag} == 'bin_drop' || ${analysisTag} == 'compare') then
    if(! -d out/${analysisTag}) then
        mkdir -p out/${analysisTag}
    endif
    if(! -d log/${analysisTag}) then
        mkdir -p log/${analysisTag}
    endif
else
    if(! -d out/${inputType}/${analysisTag}) then
        mkdir -p out/${inputType}/${analysisTag}
    endif
    if(! -d log/${inputType}/${analysisTag}) then
        mkdir -p log/${inputType}/${analysisTag}
    endif
endif


#echo ${tag}
#echo ${base}
#echo _R${radius}.root



if(${analysisTag} == 'stat_err_scaling' || ${analysisTag} == 'unfold' || ${analysisTag} == 'closure' || ${analysisTag} == 'bin_drop' || ${analysisTag} == 'compare') then
    set arg = "$radius $kappa"

    set LogFile = log/${analysisTag}/R${radius}_k${kappa}.log
    set ErrFile = log/${analysisTag}/R${radius}_k${kappa}.err

    echo "Logging output to " $LogFile
    echo "Logging errors to " $ErrFile

    echo qsub -V -l mem=4GB -o $LogFile -e $ErrFile -N ${analysisTag} -- ${ExecPath}/submit/qwrap.sh ${ExecPath} $execute $arg

    sbatch --mem-per-cpu=4GB -q express -p erhip -o $LogFile -e $ErrFile -t 180 --job-name=${analysisTag} -- ${ExecPath}/submit/qwrap.sh ${ExecPath} $execute $arg
else
    if(${tag} == 'decayed' || ${tag} == 'undecayed') then
        foreach input ( ${base}*_R${radius}_${tag}.root )
            set OutBase = `basename $input | sed 's/.root//g'`
            set OutBase = `basename $OutBase | sed 's:out/::g'`
            echo $OutBase
#            set OutBase = "$OutBase$uscore$analysisTag$uscore" + k + "$kappa"
            set OutBase = "$OutBase$uscore$analysisTag$uscore$kappa"

            set outLocation = out/${inputType}/${analysisTag}/
            set outName = ${OutBase}.root

            set inFiles = ${input}

            set LogFile = log/${tag}/${analysisTag}/${OutBase}.log
            set ErrFile = log/${tag}/${analysisTag}/${OutBase}.err

            echo "Logging output to " $LogFile
            echo "Logging errors to " $ErrFile

            set arg = "$outLocation $outName $inFiles $inputType $kappa"

            echo qsub -V -l mem=4GB -o $LogFile -e $ErrFile -N ${analysisTag} -- ${ExecPath}/submit/qwrap.sh ${ExecPath} $execute $arg

            sbatch --mem-per-cpu=4GB -q express -p erhip -o $LogFile -e $ErrFile -t 180 --job-name=${analysisTag} -- ${ExecPath}/submit/qwrap.sh ${ExecPath} $execute $arg
        end
    else
        if(${inputType} == 'data') then
            foreach input ( ${base}*_R${radius}.root )
                set OutBase = `basename $input | sed 's/.root//g'`
                set OutBase = `basename $OutBase | sed 's:out/::g'`
                echo $OutBase
                set OutBase = "$OutBase$uscore$analysisTag$uscore$kappa"

                set outLocation = out/${inputType}/${analysisTag}/
                set outName = ${OutBase}.root

                set inFiles = ${input}

                set LogFile = log/${inputType}/${analysisTag}/${OutBase}.log
                set ErrFile = log/${inputType}/${analysisTag}/${OutBase}.err

                echo "Logging output to " $LogFile
                echo "Logging errors to " $ErrFile

                set arg = "$outLocation $outName $inFiles $inputType $kappa"

                echo qsub -V -l mem=4GB -o $LogFile -e $ErrFile -N ${analysisTag} -- ${ExecPath}/submit/qwrap.sh ${ExecPath} $execute $arg

                sbatch --mem-per-cpu=4GB -q express -p erhip -o $LogFile -e $ErrFile -t 180 --job-name=${analysisTag} -- ${ExecPath}/submit/qwrap.sh ${ExecPath} $execute $arg
            end
        else
            foreach input ( ${base}*_${tag}_R${radius}.root )
                set OutBase = `basename $input | sed 's/.root//g'`
                set OutBase = `basename $OutBase | sed 's:out/::g'`
                echo $OutBase
                set OutBase = "$OutBase$uscore$analysisTag$uscore$kappa"

                set outLocation = out/${inputType}/${analysisTag}/
                set outName = ${OutBase}.root

                set inFiles = ${input}

                set LogFile = log/${inputType}/${analysisTag}/${OutBase}.log
                set ErrFile = log/${inputType}/${analysisTag}/${OutBase}.err

                echo "Logging output to " $LogFile
                echo "Logging errors to " $ErrFile

                set arg = "$outLocation $outName $inFiles $inputType $kappa"
            
                echo qsub -V -l mem=4GB -o $LogFile -e $ErrFile -N ${analysisTag} -- ${ExecPath}/submit/qwrap.sh ${ExecPath} $execute $arg

                sbatch --mem-per-cpu=4GB -q express -p erhip -o $LogFile -e $ErrFile -t 180 --job-name=${analysisTag} -- ${ExecPath}/submit/qwrap.sh ${ExecPath} $execute $arg
            end
        endif
    endif
endif



