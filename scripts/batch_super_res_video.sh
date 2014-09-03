#!/bin/bash

if [ $# -ne 2 ]; then
    echo "$0 command absolute_path_to_working_directory"
    echo "example: $0 super_res_video /home/sun/sr"
    exit
fi

EXE_FILE=$1

CURRENT_DIR=`pwd`
FILES=`find $2 -name "*.cfg"`

for line in ${FILES[@]}; do
    START=$(date +%s)

    WORKING_PATH=${line%/*}
    FILE_NAME=${line##*/}

    cd $WORKING_PATH
    cmd="$EXE_FILE $line > sr.log"
    echo $cmd
    eval $cmd
    cd $CURRENT_DIR

    END=$(date +%s)
    DIFF=$(( $END-$START ))
    echo "$DIFF seconds"
done
