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
(cdfs-sumspectra.pl), available on request (see "Author" at the end of
this document).

A listing of all options, and a detailed discussion of the steps to
prepare the input files and run the program is in the following.

This program was originally written and used for the XMM-Newton deep
survey in the Chandra Deep Field South; see the paper: Ranalli et
al. 2013, A\\&A 555, A42.
 http://cdsads.u-strasbg.fr/abs/2013A%26A...555A..42R





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

- \--help

    Prints this help message.

# DESCRIPTION

## Prepare the files before you run:

1. SOURCES:

    prepare a file named sky.txt containing a table with the following info:
        ID\_src1  RA\_1   DEC\_1   RADIUS\_1   \[RADIUS2\_1\]   \[CAMERA\_1\]
        ID\_src2  RA\_2   DEC\_2   RADIUS\_2   \[RADIUS2\_2\]   \[CAMERA\_2\]
        ...      ...    ...     ...        \[...\]         \[...\]

    where RA and DEC are in degrees, and the RADII are in "physical" units.

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

    prepare a file named bgsky.txt containing a table with the following info:
        ID\_src1  RA\_src1   DEC\_src1   RADIUS1\_src1  \[RADIUS2\_src1\]
        ID\_src2  RA\_src2   DEC\_src2   RADIUS1\_src1  \[RADIUS2\_src1\]
        ...      ...       ...        ...           \[...\]

    where RA and DEC are in degrees, and the RADII are in arcsec.
    (RADIUS2 is optional.  If present, the region is interpreted as an
    annulus, if missing, as a circle).  The same notes for SOURCES apply.

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

- \* evselect: warning (NoFilteredEvents), No events have been selected - 
      filtering of the event list resulted in an empty table
      (This warning signals that one source was completely out of the FOV
      in one observation).

Otherwise, if the following error is emitted, please report it along with
a copy of the log, of the sky.txt, bgsky.txt and event.list:

- .....: something wrong happened. Please check for errors in the log.

# VERSION HISTORY

3. 0  -- 2013/11/01 removed logic to find CDFS ccf.cif and sum.sas: these
                         files should now be stated in the evtlist; sky.txt
                         and bgsky.txt now use arcsecs instead of physical 
                         pixels;
                         evtlist logic refactored in EvtFiles package;
                         dbupdate is disabled;
                         exclusion threshold now set to exposed\_area>40% i.e.
                         deadarea<60% of PSF (and the meaning of the --exclude
                         option was reversed)
2. 93 -- 2013/9/27  started refactor and cleanup before releasing the code:
                         PSF and Image packages moved to separate files
2. 92 -- 2013/7/28  --photoerror=XXX now accepts probability level for BEHR
2. 91 -- 2012/4/27  skip photometry of sources without a bkg region
2. 9 -- 2011/9/23  reworked to use BEHR.pm; added timeout to BEHR;
corrected a bug in photometry code
2. 8 -- 2011/8/3  added checks of BACKSCAL
2. 71 -- 2011/6/21  corrected a bug in SQL code
2. 7  -- 2011/6/6   added results of spectral checks to mysql database
2. 61 -- 2011/6/6   check that SAS is initialized before starting
2. 6  -- 2011/5/12  adapted for MOSIX (--mosix option)
2. 52 -- 2011/5/11  due to bug appearing only on lafaro, modified
                         getSAS\_ODF\_CCF() to check the machine and behave
                         accordingly (see details in getSAS\_ODF\_CCF())
2. 51 -- 2011/4/14  added --evtlist=event.list option;
                         corrected a bug which prevent responses from
                         being calculated if bkg were not extracted
2. 5 -- 2011/4/6  added the choice of camera in bgsky.txt
2. 4.3 -- 2011/4/4  (internal use) made it possible to not specify any region in bgsky.txt
2. 4.2 -- 2010/11/25 sometimes it is not possible to open $$.cdfsextract.stderr in sas\_exec during the check phase. Let's try using File::Temp to avoid any possibile race condition..
2. 4.1 -- 2010/11/23 (internal use)   corrected bug in printfnlog
2. 4 -- 2010/10/21 --photoerror option added: calculates errors on net counts and rates, using BEHR (which must be available)
2. 31 -- 2010/09/29 added --noresp option
2. 3 -- 2010/08/30 in aperture photometry, psf corrections now
take into account which part of the PSF was exposed
2. 2 -- 2010/08/24 calculate, on request, rmf and arf for the bkg files;
also added options --noclobber and --maxproc
2. 1 -- 2010/07/21 can make aperture photometry
2. 0 -- 2010/06/18 added checks (whether src and bkg fall on
different chips), summary of checks, multiple positions for source or
bkg, exclusion of areas
1. 0 -- 2010/04/29  first public version

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
