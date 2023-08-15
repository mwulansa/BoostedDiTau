#!/bin/tcsh                                                                                                                                                   



foreach Sample (`ls *QCD*HT*py`)
    echo crab submit -c $Sample
    crab submit -c $Sample
end

# foreach Sample (`ls *DYJetsToLL_M-50_HT*py`)
#     echo crab submit -c $Sample
#     crab submit -c $Sample
# end

# foreach Sample (`ls *DY_M-4to50_HT*py`)
#     echo crab submit -c $Sample
#     crab submit -c $Sample
# end

# foreach Sample (`ls *WJetsToLNu_HT*py`)
#     echo crab submit -c $Sample
#     crab submit -c $Sample
# end
