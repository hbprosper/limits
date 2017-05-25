# need a robust way to find location of this file
setenv LIMITS_PATH $PWD
setenv PATH ${LIMITS_PATH}/bin:${PATH}
setenv LD_LIBRARY_PATH ${LIMITS_PATH}/lib:${LIMITS_PATH}/src:${LD_LIBRARY_PATH}
echo $LIMITS_PATH

