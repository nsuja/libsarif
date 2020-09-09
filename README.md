# libsarif - Synthetic Aperture Radar Image Formation

Library to create a 2-d image from RAW SAR data

## Description

This library is able to do the image formation process from RAW SAR data using the Range-Doppler algorithm.

The library also includes the capability to do RCMC and multilooking. The Doppler centroid is estimated using the average phase method.

The input must be provided in patches and the output is a SLC converted patch per cycle.

## Build

A simple:

```shell
$ make
$ make install
```

The default install path is `/usr/lib/` and is not parameterized.

## Requirements

For the library:

* fftw3

For the tool:

* libmatio
* ers_raw_parser

## Usage

A usage example /can be found in: `./tools/`
