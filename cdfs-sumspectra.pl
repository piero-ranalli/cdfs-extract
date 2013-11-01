#!/usr/bin/env perl
#

=head1 NAME

cdfs-sumspectra.pl -- Sum spectra extracted by cdfs-extract.pl.

=cut

# Copyright (C) 2010-2013  Piero Ranalli
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


=head1 SYNOPSIS

./cdfs-sumspectra.pl [--help] [--weights=exp|counts] src [cam [epoch]] [--keep_errors]
                     [--xspec11] [--nobkg] [--nospec] [--bkgresp]
                     [--productsdir=path/] [--joinmos]

=head1 OPTIONS

src is the src number;

cam is the camera (mos1/mos2/pn/all; default:all, meaning the sum is run
thrice and three sets of output files are produced);

and epoch (not yet available) will be used to select only a subset of the
spectra to be summed.

--weights: controls how the weights are calculated: using
the nominal (i.e. not vignetted) exposure times (--weights=exp; this
is the default) or the observed total counts (--weights=counts).

--help: prints this message and exits without summing any spectrum.

--keep_errors: do not delete the STAT_ERR and QUALITY columns from the summed spectra, nor
  change the POISSERR keyword.

--xspec11: same as --keep_errors

--nobkg: don't sum the backgrounds

--nospec: dont'sum the spectrum

--bkgresp: average the background responses

--productsdir: use specified path instead of Products/

The options names may be shortened (e.g. -k, --k, --ke, --keep and
--keep_errors all mean the same thing).

The spectra are expected to be in the Products/ directory and have names
as below.


=cut

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

{
    my $wmethod = 'exp'; #default
    my $help = 0;
    my $keep_staterr = 0; # delete STAT_ERR and QUALITY keywords in summed pha files
                          # unless explicitly stated
    my $nobkg = 0;
    my $nospec = 0;
    my $bkgresp = 0;
    my $joinmos = 0;
    my $dir = 'Products';
    GetOptions( 'weights=s'   => \$wmethod,
		'keep_errors' => \$keep_staterr,
		'xspec11'     => \$keep_staterr,
		'nobkg'       => \$nobkg,
		'nospec'      => \$nospec,
		'bkgresp'     => \$bkgresp,
		'productsdir=s' => \$dir,
		'joinmos'     => \$joinmos,
	        'help'        => \$help,
	      )
	or pod2usage(-verbose=>2);
    if ($help) {
	pod2usage(-verbose=>2);
	exit;
    }
    if ($wmethod ne 'exp' and $wmethod ne 'counts') {
	pod2usage(-msg=>'Unrecognized weighting method. Please check your --weights option.' );
	exit;
    }

    if ($#ARGV==-1) {
	print("Please provide the source name (hint: something like P34H_123).\n");
	print("If unsure, check the man page with: ./cdfs-sumspectra.pl --help\n");
	exit;
    }

    my $src = $ARGV[0];
    my $cam   = $#ARGV >= 1 ? $ARGV[1] : 'all';
    my $epoch = $#ARGV >= 2 ? $ARGV[2] : '';

    if ($cam eq 'all') {
	my @cams = $joinmos ?
	    qw/pn mos12/    :
	    qw/mos1 mos2 pn/ ;
	for my $singlecam (@cams) {
	    do_everything($dir,$src,$singlecam,$epoch,$wmethod,$keep_staterr,
		$nobkg,$nospec,$bkgresp);
	}
    } else {
	do_everything($dir,$src,$cam,$epoch,$wmethod,$keep_staterr,
		      $nobkg,$nospec,$bkgresp);
    }

}


sub do_everything {

    my ($dir,$src,$cam,$epoch,$wmethod,$keep_staterr,
	$nobkg,$nospec,$bkgresp) = @_;

    my ($spec,$bg,$rmf,$arf,$bgrmf,$bgarf) = 
	look4files($dir,$src,$cam,$epoch,$bkgresp);
    #$spec etc are array refs

    my $spec_backsc = getbackscale($spec); # xxx_backsc: array refs
    my $bg_backsc = getbackscale($bg);

    my $sumname = "${src}-${cam}-sum";
    unless( $epoch eq '' ) {
	$sumname .= "-${epoch}";
    }

    unless ($nospec) {
	sumspec($spec,$sumname.'.pha');
	setbackscale($sumname.'.pha', $spec_backsc);
	del_stat_cols($sumname.'.pha') unless $keep_staterr;
	checkbackscale($sumname.'.pha');
    }

    unless ($nobkg) {
	sumspec($bg,$sumname.'-bg.pha');
	setbackscale($sumname.'-bg.pha', $bg_backsc);
	del_stat_cols($sumname.'-bg.pha') unless $keep_staterr;
	checkbackscale($sumname.'.pha');
    }

    unless ($nospec) {
	my $weights = calc_weights($spec,$wmethod);  # a piddle
	calc_resp('addrmf',$rmf,$weights,$sumname.'.rmf');
	calc_resp('addarf',$arf,$weights,$sumname.'.arf');
    }

    if ($bkgresp) {
	my $weights = calc_weights($bg,$wmethod);  # a piddle
	calc_resp('addrmf',$bgrmf,$weights,$sumname.'-bg.rmf');
	calc_resp('addarf',$bgarf,$weights,$sumname.'-bg.arf');
    }

    print("Spectra summed for source $src epoch $epoch.\n");
}

=head1 DESCRIPTION

=head2 About filenames

Input filenames are assumed to have the following structure:

 src-obsid-cam.pha     for source spectra
 src-obsid-cam-bg.pha  for background spectra
 src-obsid-cam.rmf     for (source) responce matrices
 src-obsid-cam.arf     for (source) ancillary response files

while output filenames will be

 src-cam-sum[-epoch].pha
 src-cam-sum[-epoch]-bg.pha
 src-cam-sum[-epoch].rmf
 src-cam-sum[-epoch].arf

=cut

sub look4files {

    my ($dir,$src,$cam,$epoch,$bkgresp)=@_;

    # epoch will be selected here

    my $oldcam = $cam;
    if ($oldcam eq 'mos12') {
	$cam = 'mos1';
    }

    my @spectra = glob("$dir/${src}-*-${cam}.pha");
    my @bkgs = glob("$dir/${src}-*-${cam}-bg.pha");
    my @rmfs = glob("$dir/${src}-*-${cam}.rmf");
    my @arfs = glob("$dir/${src}-*-${cam}.arf");
    my (@bgrmfs,@bgarfs);
    if ($bkgresp) {
	@bgrmfs = glob("$dir/${src}-*-${cam}-bg.rmf");
	@bgarfs = glob("$dir/${src}-*-${cam}-bg.arf");
    }


    if ($oldcam eq 'mos12') {
	# add mos2 files
	$cam = 'mos2';

	push @spectra, glob("$dir/${src}-*-${cam}.pha");
	push @bkgs, glob("$dir/${src}-*-${cam}-bg.pha");
	push @rmfs, glob("$dir/${src}-*-${cam}.rmf");
	push @arfs, glob("$dir/${src}-*-${cam}.arf");
	if ($bkgresp) {
	    push @bgrmfs, glob("$dir/${src}-*-${cam}-bg.rmf");
	    push @bgarfs, glob("$dir/${src}-*-${cam}-bg.arf");
	}
    }


    if ($epoch) {
	my @obsid;
	if ($epoch eq '2001') {
	    @obsid = qw/0108060401 0108060501/;
	} elsif ($epoch eq '2002') {
	    @obsid = qw/0108060601 0108060701 0108061801 0108061901 0108062101 0108062301/;
	} elsif ($epoch eq 'sum2008') {
	    @obsid = qw/0555780101 0555780201 0555780301 0555780401/;
	} elsif ($epoch eq 'win2009') {
	    @obsid = qw/0555780501 0555780601 0555780701 0555780801 0555780901 0555781001 0555782301/;
	} elsif ($epoch eq 'sum2009') {
	    @obsid = qw/0604960101 0604960201 0604960301 0604960401/;
	} elsif ($epoch eq 'win2010') {
	    @obsid = qw/0604960501 0604960601 0604960701 0604960801 0604960901 0604961001 0604961101 0604961201 0604961301 0604961801/;
	}

	my @varlist = (\@spectra, \@bkgs, \@rmfs, \@arfs);
	if ($bkgresp) {
	    push(@varlist,(\@bgrmfs,\@bgarfs));
	}
	for my $var (@varlist) {
	    my @new = ();
	    for my $i (@obsid) {
		push(@new, grep(/$i/, @$var));
	    }
	    @$var = @new;
	}
    }

    if ($bkgresp) {
	return(\@spectra,\@bkgs,\@rmfs,\@arfs,\@bgrmfs,\@bgarfs);
    } else {
	return(\@spectra,\@bkgs,\@rmfs,\@arfs);
    }

}


sub getbackscale {
    my $spectra = shift;

    my @backscales;
    for my $spec (@$spectra) {

	open(my $fkey, "fkeyprint $spec+1 BACKSCAL |");
	while (my $line = <$fkey>) {
	    next unless ($line =~ m/BACKSCAL=/);

	    my (undef,$value) = split(/\s+/,$line);

	    if ($value < 1e-3) {
		die "BACKSCAL=0 in file $spec\n";
	    }
	    push(@backscales,$value);
	    last;
	}
	close($fkey);
    }

    return(\@backscales);
}


=head2 Sum of spectra: propagation of errors

There are different possibilities in mathpha to calculate the errors
on the counts. Here we use properr=no, which means NO propagation of
error calculation (actually, we want no errors at all in the
spectra). This makes sense as we regard the different observations in
a survey as discrete steps in accumulation of counts, and only as long
as the responses are not too different among the different
observations of a single source.  Also, this is probably the only way
to proceed in the regime of very low counts per bin, where any formula
is probably overestimating the final errors (in the spectra, we have
that most of the bins have no counts, then a sizable fraction has just
1 count per bin, and a small minority of bins may have a few counts).

=cut


sub sumspec {
    use File::Temp qw/tempfile/;

    my ($spec,$outfile) = @_;

    # create expression
    my $expr = "'".join("'+'",@$spec)."'";

    # for some unknown reason, mathpha has troubles with the filenames
    # when the expression is part of the command line:
    # 90-0555782301-mos1.pha is interpreted as an expression, trying
    # to subtract mos1.pha from 90 and from 0555782301. This does not
    # happen if mathpha is called with the @file sintax pointing to a
    # file containing the expression.

    # create temporary pha list file
    my ($sfileh,$sfile) = tempfile(
		'tmp-cdfslistXXXX', DIR=>'/tmp', UNLINK=>0 );
    print($sfileh "$expr\n");
    close($sfileh);


    # call mathpha
    my ($ok,$stdout,$stderr) = call(
	"mathpha expr=\@$sfile units='C' outfil=$outfile properr=no exposure=CALC backscal=NONE areascal='%' ncomments=0");

    die unless ($ok);
}




=pod

However, even with properr=no, mathpha calculates statistical errors
with Poissonian statistics. These need to be removed, by deleting the
STAT_ERR and QUALITY columns in the summed spectra. Finally, the POISSERR
keyword of the header is set to True.

=cut

sub del_stat_cols {

    my $name = shift;

    my ($ok,$out,$err) = call("fdelcol ${name}+1 STAT_ERR N Y");
    unless ($ok) {
	die("Could not delete STAT_ERR column in file $name.\n$err");
    }

    ($ok,$out,$err) = call("fdelcol ${name}+1 QUALITY N Y");
    unless ($ok) {
	die("Could not delete QUALITY column in file $name.\n$err");
    }

    ($ok,$out,$err) = call("fparkey T ${name}+1 POISSERR");
    unless ($ok) {
	die("Could not set POISSERR keyword in file $name.\n$err");
    }


}




sub call {
  my $cmd = shift;

  next if ($cmd =~ m/^\s*(\#.*)?$/);  # empty line or comment

  print("$$ COMMAND: $cmd\n");
  my $out = `$cmd   2> /tmp/$$.cdfsextr.stderr`;

  # check exit status and print output and/or errors
  if ($?) {
      print("$$: something wrong happened. Please check for errors.\n");
  }
  print("$$ STDOUT:\n$out");

  local $/;   #SLURP
  open(my $ERR, "< /tmp/$$.cdfsextr.stderr");
  my $err = <$ERR>;
  close($ERR);
  unlink("/tmp/$$.cdfsextr.stderr");

  # exit now if any error
  if ($?) {
      return(0,$out,$err);
  } else {
      return(1,$out,$err); # ok
  }
}


=head2 About BACKSCALes

The xspec manual says that the data to be fitted are the
background-subracted rates C, defined as:

C(i) = D(i) / a_d / t_d  -  B(i) * (b_d / b_b) / a_b / t_b

where a are the AREASCALes; t the EXPOSUREs; D(i) and B(i) the data
and background counts, respectively, in the bin i; and the _d and _b
pedices refer to the data and background.

(The a's may also depend on the bin, but in that case an AREASCAL
column should be present in the spectra, and the header keyword would
not be used. This, however, only happens in grating spectra, and will
be ignored here).

The backscales may thus be regarded as normalisations, and will be
treated as such: the summed source and background spectra will contain
in BACKSCAL the sum of all the backscales of the source_obsid and
background_obsid spectra, respectively.

=cut

sub setbackscale {
    my ($outfile, $scales) = @_;  # $scales is an array ref

    my $sum = 0;
    for my $x (@$scales) {
	$sum += $x;
    }

    my ($ok,$out,$err) = call("fparkey $sum $outfile+1 BACKSCAL");

    unless ($ok) {
	die("Something wrong in fparkey:\n$out");
    }
}


=head2 Averaging of RMFs and ARFs.

The RMFs are averaged by addrmf, using the nominal (unvignetted)
exposure times, or the actual counts in the spectra as the
weights. The same is made to the ARFs by addarf.

=cut


sub calc_resp {
    use File::Temp qw/tempfile/;

    my ($ftool,$infiles,$weights,$outfile) = @_;  # $rmf and $spec are array refs

    # open temporary RMF/ARF list file
    my ($sfileh,$sfile) = tempfile(
		'tmp-cdfslistXXXX', DIR=>'/tmp', UNLINK=>0 );
    # write list of responses and weights
    for my $i (0..$#$infiles) {
	printf($sfileh "%s %9.5e\n",$$infiles[$i],$$weights[$i]);
    }
    close($sfileh);

    # complete syntax of tool
    if ($ftool eq 'addrmf') {
	$outfile = "rmffile=$outfile";
    } elsif ($ftool eq 'addarf') {
	$outfile = "out_ARF=$outfile";
    } # otherwise do nothing

    # average the rmf/arf's
    my ($ok,$out,$err) = call("$ftool \@$sfile $outfile");

    unless($ok) {
	die("Error in $ftool:\n$out\n\n");
    }

}


sub calc_weights {  # in and out: array refs
    use PDL;
    # there is some going back and forth between lists and piddles,
    # but the way it works, PDL is only required here. And unless
    # there was a wcols which worked with strings, using PDL more
    # extensively would not help in this program.


    my $spectra = shift;  # array ref
    my $method = shift;
    my @weights;

    for my $spec (@$spectra) {

	my $s = rfits($spec.'[1]');

	if ($method eq 'counts') {
	    push(@weights,$s->{COUNTS}->sum);
	} else {
	    my $w = $s->{hdr}->{EXPOSURE};
	    push(@weights,$w);
	}
    }

    # normalize
    my $w = pdl @weights;
    $w /= $w->sum;
    @weights = $w->list;

    return(\@weights);
}




sub checkbackscale {
    my $spec = shift;

    open(my $fkey, "fkeyprint $spec+1 BACKSCAL |");
    while (my $line = <$fkey>) {
	next unless ($line =~ m/BACKSCAL=/);

	my (undef,$value) = split(/\s+/,$line);

	if ($value < 1e-3) {
	    die "BACKSCAL=0 in file $spec\n";
	}
	last;
    }
    close($fkey);
}


=head1 LICENSE

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



=head1 VERSION

=over 4

=item v. 2.71 -- 2013/11/1

public version

=item v. 2.7 -- 2012/6/6

can sum the MOSes together, either with --joinmos
or by choosing mos12 as camera

=item v. 2.6 -- 2011/8/3

added checks for BACKSCAL

=item v. 2.51 -- 2011/7/15

gives some help if called without source name

=item v. 2.5 -- 2011/5/13

added --productsdir

=item v. 2.4 -- 2011/4/26

added --nospec e --bkgresp options, to average bkg responses

=item v. 2.3 -- 2011/4/24

summing by epoch added (6 predefined epochs by now)

=item v. 2.2 -- 2011/4/13

added --nobkg option

=item v. 2.1 -- 2011/4/7

Changes POISSERR keyword to T to avoid warnings in XSPEC12, unless
--keep_errors is specified.  Also aliased --keep_errors to --xspec11.

=item v. 2 -- 2010/6/8

Added --weights option, weights responses on exposure times on
default, and deletes statistical error columns from the summed
spectra.

=item v. 1 -- 2010/6/4

First internal version.

=back

=head1 AUTHOR

Piero Ranalli

Post-doc researcher at IAASARS, National Observatory of Athens, Greece.
Associate of INAF -- Osservatorio Astronomico di Bologna, Italy.

pranalli.github@gmail.com

=cut

