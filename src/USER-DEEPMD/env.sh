DEEPMD_ROOT=/home/dtisi/deepmd-1.3.3/deepmdlib
TENSORFLOW_INCLUDE_DIRS="/opt/sissa/libs/libtensorflow_cc/1.14.0/include;/opt/sissa/libs/libtensorflow_cc/1.14.0/include"
TENSORFLOW_LIBRARY_PATH="/opt/sissa/libs/libtensorflow_cc/1.14.0/lib;/opt/sissa/libs/libtensorflow_cc/1.14.0/lib"

TF_INCLUDE_DIRS=`echo $TENSORFLOW_INCLUDE_DIRS | sed "s/;/ -I/g"`
TF_LIBRARY_PATH=`echo $TENSORFLOW_LIBRARY_PATH | sed "s/;/ -L/g"`
TF_RPATH=`echo $TENSORFLOW_LIBRARY_PATH | sed "s/;/ -Wl,-rpath=/g"`

NNP_INC=" -std=c++11 -DHIGH_PREC   -I$TF_INCLUDE_DIRS -I$DEEPMD_ROOT/include/deepmd "
NNP_PATH=" -L$TF_LIBRARY_PATH -L$DEEPMD_ROOT/lib"
NNP_LIB=" -Wl,--no-as-needed -ldeepmd_op_cuda -ldeepmd_op -ldeepmd -ltensorflow_cc -ltensorflow_framework -Wl,-rpath=$TF_RPATH -Wl,-rpath=$DEEPMD_ROOT/lib"
