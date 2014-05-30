# Time-Varying Filters for Musical Applications

This repository contains the code and data used by *Time-Varying Filters for
Musical Applications*, submitted to the DAFx-14 conference.

## Building

The filtering code requires [Eigen](http://eigen.tuxfamily.org/),
[libsndfile](http://www.mega-nerd.com/libsndfile/), and a modern C++ compiler
and standard library supporting C++11. [Clang](http://clang.llvm.org/) and
[libc++](http://libcxx.llvm.org/) are recommended, but modern GCC should work
as well, though it has not been tested.

Once these dependencies are installed, just run `make` in the `code` directory
to generate the `filter` tool.

## Generating stimuli

Note that the stimuli are included pre-generated in the `data/stimuli` folder, so it
is not necessary to build the filter tool, except for the purpose of
reproducibility.

The `filter` tool can be used to generate all the stimuli used in the subjective
and objective tests. Its usage is:

    ./filter [filter type] [sine_hz|dc] [output_dir]

The second argument can either be a number, to generate a sine wave, for the
listening tests, or "dc" to generate DC for the objective experiment. For the
listening tests in the paper, a 100 Hz sine wave was used.

There is a python script, `code/generate_stimuli.py`, which automatically generates
all the stimuli used in the paper, by running the `filter` tool.

## Data

The data collected from the subjective listening test is included in the
`data/results` folder. It was anonymized by running the
`code/anonymize_data.py` script, which removes dates and subjects' names. The
data was output by [MUSHRATest](http://mushra.kosobrodov.net/).

## License

The code, consisting of everything in the `code` folder, is released under the
BSD 3-Clause License. The data, consisting of everything in the `data` folder,
is released under the Creative Commons CC0 1.0 Universal license. See
`code/LICENSE.md` and `data/LICENSE.md` respectively for license details.
