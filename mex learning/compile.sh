#!/bin/bash

# Check if an argument is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 file.c"
    exit 1
fi

# Call gcc
/usr/bin/gcc -c -DMX_COMPAT_64  -DMATLAB_MEXCMD_RELEASE=R2018a  -DUSE_MEX_CMD   -D_GNU_SOURCE -DMATLAB_MEX_FILE  -I"/usr/local/MATLAB/R2024a/extern/include" -I"/usr/local/MATLAB/R2024a/simulink/include" -fexceptions -fPIC -fno-omit-frame-pointer -pthread -fwrapv -g $1
# /usr/bin/gcc -c -DMX_COMPAT_64  -DMATLAB_MEXCMD_RELEASE=R2018a  -DUSE_MEX_CMD   -D_GNU_SOURCE -DMATLAB_MEX_FILE  -I"/usr/local/MATLAB/R2024a/extern/include" -I"/usr/local/MATLAB/R2024a/simulink/include" -fexceptions -fPIC -fno-omit-frame-pointer -pthread -fwrapv -g "/usr/local/MATLAB/R2024a/extern/version/c_mexapi_version.c" -o /tmp/mex_7673166058528_11773/c_mexapi_version.o