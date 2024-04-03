function mexCcode(debug,verbose)
  % This function helps to mex the c code on windows vs linux. Mac is not supported.

  arguments
    debug (1,1) logical = false;
    verbose (1,1) logical = false;
  end

MEX_API = "-R2018a";

if ispc
  blasLibrary = "-llibmwblas";
  lapackLibrary = "-llibmwlapack";
else
  blasLibrary = "-lmwblas";
  lapackLibrary = "-lmwlapack";
end

if debug
  debugFlag = "-g";
else
  debugFlag = "";
end

if verbose
  verboseFlag = "-v";
else
  verboseFlag = "";
end

% gaussianEliminationWithPartialPivotingC.c
mex("gaussianEliminationWithPartialPivotingC.c", MEX_API, debugFlag, verboseFlag);

% gaussianEliminationWithPartialPivotingCblas.c
mex("gaussianEliminationWithPartialPivotingCblas.c", MEX_API, blasLibrary, debugFlag, verboseFlag);

% dgesvLAPACK.c
mex("dgesvLAPACK.c", MEX_API, lapackLibrary, debugFlag, verboseFlag);

end