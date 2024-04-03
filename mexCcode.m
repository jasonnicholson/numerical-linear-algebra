function mexCcode(debug,verbose)
  % This function helps to mex the c code on windows vs linux. Mac is not supported.

  arguments
    debug (1,1) logical = false;
    verbose (1,1) logical = false;
  end

% the mex calls are different on Windows vs Linux
FILES_TO_MEX = ["gaussianEliminationWithPartialPivotingC.c";
    "gaussianEliminationWithPartialPivotingCblas.c"];

mexAPI = "-R2018a";

if ispc
  blasLibrary = "-libmwblas";
else
  blasLibrary = "-lmwblas";
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

for i=1:numel(FILES_TO_MEX)
  mex(FILES_TO_MEX(i),mexAPI,blasLibrary,debugFlag,verboseFlag);
end

end