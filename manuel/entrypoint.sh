#!/bin/bash
##############################################################################

# Used in Dockerfile to executes pipeline according to MODE variable

##############################################################################


# Retrieve the execution mode from environment variable (MODE)
# If no value is provided, default execution is AUTO
MODE=${MODE:-AUTO}

# If MODE is "AUTO", execute control script in which while loop executes main python processing script
if [ "$MODE" = "AUTO" ]; then

    /opt/conda/bin/conda run -n dmpipeline python control_autorun.py > /proc/1/fd/1 2>&1

# IF MODE is "MANUAL"/"INDICES"/"DROUGHT", execute python main processing script (only once)
elif [ "$MODE" = "MANUAL" -o "$MODE" = "INDICES" -o "$MODE" = "DROUGHT" ]; then

    /opt/conda/bin/conda run -n dmpipeline ${SCRIPT_PATH} > /proc/1/fd/1 2>&1

# WRONG input MODE
else
    echo "WRONG INPUT FOR PROCESSING MODE : $MODE"

fi