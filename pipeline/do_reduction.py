# coding: utf-8
#
#    Project: BioCAT SAXS pipeline
#             https://github.com/biocatiit/saxs-pipeline
#
#
#    Principal author:       Jesse Hopkins
#
#    This is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This software is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this software.  If not, see <http://www.gnu.org/licenses/>.

import glob
import os.path

import bioxtasraw.RAWAPI as raw

"""
To do:
1) write functions simialr to do_analysis for simple reduction tasks
"""

def reduce_and_save_images(data_dir, out_dir, raw_settings, img_prefix='*'):
    files = glob.glob(os.path.join(data_dir, img_prefix))
    profiles = raw.load_and_integrate_images(files, raw_settings)

    for profile in profiles:
        raw.save_profile(profile, data_dir=out_dir, settings=raw_settings)
