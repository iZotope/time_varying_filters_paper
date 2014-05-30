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
"""Anonymize data output by MUSHRATest, by renaming files, and removing
subjects' names.

usage: anonymize_data.py [input_path] [output_path]
"""
import os
import sys


def main(input_path, output_path):
    files = [f for f in os.listdir(input_path) if os.path.splitext(f)[1] ==
             '.txt']
    for idx, filename in enumerate(files):
        with open(os.path.join(input_path, filename), 'r') as f:
            lines = f.readlines()
        output_filename = '{}.txt'.format(idx+1)
        with open(os.path.join(output_path, output_filename), 'w') as f:
            f.writelines(lines[3:])


if __name__ == '__main__':
    try:
        sys.exit(main(sys.argv[1], sys.argv[2]))
    except:
        print(__doc__)
