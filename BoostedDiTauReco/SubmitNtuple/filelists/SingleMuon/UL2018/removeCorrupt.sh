#!/bin/tcsh

foreach d (`ls -d Ntuple_SingleMuon_Run2018*`)
    echo $d
    foreach f (`ls ${d}/*SingleMuon*txt`)
	echo $f
	foreach i (`cat corruptFile.txt`)
	    sed -i "\|${i}|d" $f
	end
    end
end

# foreach f (`ls *SingleMuon*txt`)
#     echo $f
#     foreach i (`cat corruptFile.txt`)
#     	sed -i "\|${i}|d" $f
#     end
# end
