import os
import subprocess

month = './july'

ww3_utils = '../ww3_utils/post_processing/'

commands = [['ln', '-srf',  ww3_utils+'ww3_ounp.py',           month],
            ['ln', '-srf',  ww3_utils+'ww3_ounf.py',           month],
            ['ln', '-srf',  ww3_utils+'process_points.py',     month],
            ['ln', '-srf',  ww3_utils+'process_fields.py',     month],
            ['ln', '-srf',  ww3_utils+'plot_points.py',        month],
            ['ln', '-srf',  ww3_utils+'plot_fields.py',        month],
            ['cp',         ww3_utils+'ww3_ounp.config',        month],
            ['cp',         ww3_utils+'ww3_ounf.config',        month],
            ['cp',         ww3_utils+'process_points.config',  month],
            ['cp',         ww3_utils+'process_fields.config',  month],
            ['cp',         ww3_utils+'plot_points.config',     month],
            ['cp',         ww3_utils+'plot_fields.config',     month],
            ['cp',         './ww3_ounp',                       month],
            ['cp',         './ww3_ounf',                       month],
            ['ln', '-srf', './obs_data',                       month]]
            

for cmd in commands:
  print ' '.join(cmd)
  subprocess.call(cmd)
