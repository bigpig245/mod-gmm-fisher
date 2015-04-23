% compile
disp('Making mexGmmTrainSP...');
mex /home/ntrang/project/tools/gmm-fisher.git/matlab/mexGmmTrainSP.cxx /home/ntrang/project/tools/gmm-fisher.git/gmm.cxx /home/ntrang/project/tools/gmm-fisher.git/stat.cxx /home/ntrang/project/tools/gmm-fisher.git/simd_math.cxx -largeArrayDims CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
disp('Making mexGmmTrainDP...');
mex /home/ntrang/project/tools/gmm-fisher.git/matlab/mexGmmTrainDP.cxx /home/ntrang/project/tools/gmm-fisher.git/gmm.cxx /home/ntrang/project/tools/gmm-fisher.git/stat.cxx /home/ntrang/project/tools/gmm-fisher.git/simd_math.cxx -largeArrayDims CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"

disp('Making mexFisherEncodeHelperSP...');
mex /home/ntrang/project/tools/gmm-fisher.git/matlab/mexFisherEncodeHelperSP.cxx /home/ntrang/project/tools/gmm-fisher.git/fisher.cxx /home/ntrang/project/tools/gmm-fisher.git/gmm.cxx /home/ntrang/project/tools/gmm-fisher.git/stat.cxx /home/ntrang/project/tools/gmm-fisher.git/simd_math.cxx -largeArrayDims CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
