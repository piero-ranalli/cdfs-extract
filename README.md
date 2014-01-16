# NEWS

New version 3.1!  With a smarter detection of out-of-FOV sources and cleaner error messages.

# NAME

cdfs-extract -- A program to extract spectra and aperture photometry
from multiple sources in multiple XMM-Newton observations.

# ABOUT

This program can extract spectra and aperture photometry from
XMM-Newton observations, calling Ftools and SAS tasks if needed. The
observations may of may not overlap. This program checks if the source
and background are in the field of view of any observation, and
extract products accordingly.  Responses (RMFs and ARFs) are
calculated along with spectra.

The spectra and responses can be summed with a companion program
(cdfs-sumspectra.pl), available as part of this distribution. (See its
documentation with ./cdfs-sumspectra.pl --help).

A listing of all options, and a detailed discussion of the steps to
prepare the input files and run the program is in the following.

This program was originally written and used for the XMM-Newton deep
survey in the Chandra Deep Field South; see the paper: Ranalli et
al. 2013, A\\&A 555, A42.
 http://cdsads.u-strasbg.fr/abs/2013A%26A...555A..42R

cdfs-extract is written in the Perl programming language, using the
Perl Data Language modules.




# SYNOPSIS

./cdfs-extract.pl \[options\]

(for help: ./cdfs-extract.pl -h)

# OPTIONS

- \--check

    Only runs checks the source and background coverage and if the source
    and bkg fall on the same chips. Produces a screen log of the
    result. (The checks results are also written in
    cdfs-extract.log. Also, a useful unix command is tee:
     ./cdfs-extract.pl --check | tee checks.txt

    This will send the screen output also to the file checks.txt).

- \--dbupdate

    (experimental feature) Updates the mysql database 'spectralchecks'
    with the results of the checks.

- \--bkgresp

    Also calculates the RMFs and ARFs for the background spectra.

- \--exclude=XX

    Exclude a source/file combination from spectra extraction, if the exposed
    area is less than XX (default: 0.40).

- \--evtlist=filename

    Use filename instead of the default event.list.  (Especially useful to
    do photometry in different bands, or for tests on shorter event lists).

- \--photometry \[--energy=XX\]

    Runs (only) checks and aperture photometry. The energy (in keV;
    optional, default is 2 keV) is used to calculate the PSF fraction
    covered by the aperture.

- \--photoerror\[=XX\]

    Calculates the errors on net counts and net rate when doing
    photometry, using the external program BEHR. Also, switches the
    calculation of the net counts and rate from the "classical" definition
    to the median of the posterior distribution from BEHR.  If the latter
    program is not available, prints a warning and continues without
    calculating the errors.

    The optional parameter specifies the per cent probability level
    (default=90).

- \--noclobber

    Does not overwrite files. May be useful if you want to add one
    observation to the event.list and only create the relevant
    spectra/responses, without recreating any already existant
    file. Another use may be, together with --bkgresp, to extract only
    background responses for spectra which have already been extracted.

- \--mosix

    Enables the execution of selected commands (currently: backscale,
    rmfgen, arfgen) through the MOSIX clustering system.

- \--maxproc=X

    Uses a maximum number of X processes for the extraction of spectra
    (default: 8). For maximum speed, X should be the same of the number of
    CPU cores.

    However, during the check and photometry phases, the bottleneck is
    rather the disc I/O, and thus only for these phases the number of
    processes is limited to max(X,4).

- \--keep\_tempfiles

    Keeps the temporary files in /tmp (may be useful for debug).

- \--deadthreshks=X

    Sets a minimum threshold for exposure (in ks) when calculating the
    exposed fractions (default: 1 ks). May be useful in some cases
    where the chip gaps have tiny exposure times (e.g. 100 s) which fool
    the calculation.

- \--oldskytxt

    Uses old-format sky.txt for compatibility with versions up to 2.93.

- \--help

    Prints this help message.

# DESCRIPTION

## Prepare the files before you run:

1. SOURCES:

    prepare a file named sky.txt containing a table with the following info:
        ID\_src1  RA\_1   DEC\_1   RADIUS\_1   \[RADIUS2\_1\]   \[CAMERA\_1\]
        ID\_src2  RA\_2   DEC\_2   RADIUS\_2   \[RADIUS2\_2\]   \[CAMERA\_2\]
        ...      ...    ...     ...        \[...\]         \[...\]

    where RA and DEC are in degrees, and the RADII are in arcsec.

    Notes:

    - The columns marked with the square brackets above are optional.
    - If RADIUS2 is specified, then the program assumes that the
    source is an annulus; otherwise, it is a circle.
    - If CAMERA is specified, then the RA, DEC, and radii apply only
    to event files from that camera (this is especially useful to choose
    camera-dependent backgrounds: see below for bgsky.txt).
    - You can add comments in sky.txt by starting a line with
    the \# character.
    - There may be multiple lines with the same ID\_src. This means
    that the spectrum will that be extracted will have photons coming from
    all the regions: use this feature if you want average spectra.
    - If you want to exclude any area from the extraction of a
    spectrum (say, because there is a nearby source which should not be
    included), put the word "exclude" between ID\_src and RA\_src, as in the
    following:
        ID\_src1  exclude RA\_src99   DEC\_src99   RADIUS\_src99

        thus the area around RA\_src99,DEC\_src99 will be excluded from the
        spectrum of ID\_src1.  You can of course exclude multiple positions
        (use the same ID\_src for all).

2. BACKGROUND:

    prepare a file named bgsky.txt containing a table with the same
    structure of sky.txt. Use the same ID\_src to identify the sources,
    but give RA, DEC and radii relative to the background regions.
    The same notes for SOURCES apply.

3. EVENT FILES:

    prepare a file named event.list which contains the paths of all event
    files, exposure maps, images, ccf.cif and \*SUM.SAS, one per line as in
    the following scheme:
       evtfile1.fits  expmap1.fits img1.fits ccf1.cif xxxxxxSUM.SAS1
       evtfile2.fits  expmap2.fits img2.fits ccf2.cif xxxxxxSUM.SAS2
       ...            ...          ...       ...      ...

    Images are optional, since they are only used when doing aperture
    photometry. In this case, put NA as a placeholder.

    The ccf.cif and \*SUM.SAS files are needed to link the event files to
    the relevant calibration data.

    You may comment out events that you don't want to process, by
    prepending the \# character.  This program will take care of
    recognising the camera (MOS1, MOS2, PN) and use the right command
    templates accordingly.

4. COMMAND TEMPLATES:

    check that the command templates are ok, in terms of energy, pattern and
    flag filtering. The templates are:
     gsp.mos1.tmpl, gsp.mos2.tmpl, gsp.pn.tmpl  (extraction of spectra and calculations of backscale)
     response.tmpl   (generation of rmf, arf)

5. CLEAN SPACE:

    all the produced spectra and responses will be put in the Products/
    subdirectory. Please check that this does not contain older files
    which may get mixed with the new ones and make confusion. Take special
    care in this sense if you are going to use --noclobber.

6. RUN THE PROGRAM:

    start HEADAS (FTOOLS), SAS, and call ./cdfs-extract.pl with any needed
    option.

7. CHECK THE LOG:

    a log will be written in cdfs-extract.log. Please check that everything
    was processed correctly.

## Checks

The program checks that at least a significant fraction of the source
and background regions have been covered in each observation. Sources
may in fact not be exposed because they are outside of the FOV, or
behind chip gaps, etc.  Warnings are emitted (but spectra are still
extracted) if <50% of the PSF is exposed.  Spectra are not extracted
if <40% of the PSF is exposed (this threshold may be changed with
\--exclude=XX). The checks refer to "dead area" fractions which are
defined as 1 - exposed\_fraction.

The program also checks if the source and background spectra fall on
the same chip; a warning is shown if this is not the case.

If you just want the checks and need not to proceed further, you can
use the option --check .

## Spectra and responses

Source and background spectra, and source responses (rmf and arf) are
calculated by default.

If you also want background responses, add the option --bkgresp to the
command line. If you have already extracted the spectra and the source
responses, and you just need the background responses, you can add the
option --noclobber to avoid regenerating any already present file.

## Photometry

Prepare the files as above, and call ./cdfs-extract.pl --photometry .
Spectra and responses will not be extracted in this case. Instead, a
report with results from aperture photometry will be printed on
screen, and also saved to the file photometry.dat.

If you want the energy encircled fractions to be calculated at
energies different than the default (2 keV), there is the option
\--energy=XX , where XX is in keV.

Aperture photometry is calculated for each observation for each camera
in the following way:

1. Only the pixels inside the source region, and which are
active in the exposure map (i.e. are inside the FOV, do not fall in a
CCD gap, have actually been exposed, etc.) are considered. These are
referred to as "active pixels", in the following. The source counts
are extracted from the image, in the same active pixels of the expmap.
2. The source exposure is calculated, by averaging the values of
the expmap in the active pixels.
3. The source area is calculated, as the sum of the active pixel areas.
4. The same is done for the background region.
- 5a. Unless --photoerror is used: the net counts are

        netcts = src_cts - bkg_cts*src_area/bkg_area
       And the net rate is
        rate_nopsf = src_cts/src_exp - bkg_cts/bkg_exp * (src_area/bkg_area)
- 5b. If --photoerror is used: the net counts are taken as the
median of the (Bayesian) posterior distribution of source intensity,
given the observed source and bkg counts, and the relative areas and
exposures.  The 90% errors on the net counts and net rate are also
calculated. The external program BEHR is used for the calculations
(Park et al. 2006, ApJ 652, 610).
6. The encircled energy fraction is calculated based on the PSF,
as the fraction of flux actually falling on the (source) active pixels. The
corrected rate is thus

        rate = rate_nopsf / encircled_energy_fraction



# WARNINGS AND ERRORS

Some warnings emitted by SAS can usually be considered safe, such as the following:

- Could not read expmap (file) or source not specified for its
camera.  

        This signals one of two possibilities, which in the current version are
        not distiguished:
        1) a source region only applies to one or two cameras (eg MOS1 and
        MOS2), but an event file is specified for another too (eg PN); in
        this case it is safe to ignore it.
        2) the exposure map was not readable for some reason (wrong filename?).

- \* evselect: warning (NoFilteredEvents), No events have been selected - 
      filtering of the event list resulted in an empty table

        This warning signals that one source was completely out of the FOV
        in one observation.

Otherwise, if the following error is emitted, please report it along with
a copy of the log, of the sky.txt, bgsky.txt and event.list:

- .....: something wrong happened. Please check for errors in the log.

# VERSION HISTORY

- 3.1  -- 2013/12/23   internal refactor, smarter detection of out-of-FOV sources
- 3.0  -- 2013/11/01   first public version!
- (lots of development)
- 1.0  -- 2010/04/29   initial version

# LICENSE

Copyright (C) 2010-2013  Piero Ranalli

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU Affero General Public
License along with this program.  If not, see
<http://www.gnu.org/licenses/>.

# AUTHOR

Piero Ranalli

Post-doc researcher at IAASARS, National Observatory of Athens, Greece
Associate of INAF -- Osservatorio Astronomico di Bologna, Italy

pranalli.github@gmail.com
