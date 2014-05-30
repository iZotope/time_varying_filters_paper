# Copyright (c) 2014, Aaron S Wishnick
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in the
#   documentation and/or other materials provided with the distribution.
# * Neither the name of the Aaron S Wishnick, iZotope, Inc., nor the
#   names of its contributors may be used to endorse or promote products
#   derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""Generate all stimuli by running the "filter" tool. Requires this tool to be
built. Stimuli will be output in a "stimuli" folder inside the current working
directory.
"""
import os
import errno
import subprocess


def ensure_directory_exists(d):
    try:
        os.makedirs(d)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def run_filter(filter_type, freq, output_dir):
    filter_path = os.path.join(os.path.dirname(__file__), 'filter')
    p = subprocess.Popen([
        filter_path, str(filter_type), '{}'.format(freq), output_dir],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False)
    out, err = p.communicate()


def generate_trials(path, freq, filter_types):
    for filter_type in filter_types:
        output_dir = os.path.join(path, filter_type)
        ensure_directory_exists(output_dir)
        run_filter(filter_type, freq, output_dir)


def main():
    freq = 100
    filter_types = {
        'lowpass_freq',
        'lowpass_q',
        'peaking_gain',
        'peaking_freq',
        'peaking_q',
    }

    generate_trials('stimuli/MUSHRA', freq, filter_types)
    generate_trials('stimuli/dc', 'dc', filter_types)


if __name__ == '__main__':
    main()
