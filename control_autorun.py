"""
##############################################################################

Control automatic running of processing chain

##############################################################################
"""

import os
import time
import subprocess
import logging
logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s', datefmt='%d/%m/%Y %H:%M:%S', level=logging.INFO)


SCRIPT_PATH = os.environ.get('SCRIPT_PATH')
TIME_SLEEP = 24*60*60 # => 24h time sleep between each run

while True:

  process = subprocess.Popen(f'{SCRIPT_PATH} > /proc/1/fd/1 2>&1', shell=True)
  process.wait()
  return_code = process.returncode

  if return_code==0:
    logging.info(f'\n\nScript {SCRIPT_PATH} completed successfully')
  else:
    logging.info(f'\n\nScript {SCRIPT_PATH} completed with return code : {return_code}')

  time.sleep(TIME_SLEEP)
