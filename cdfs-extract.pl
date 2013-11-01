#!/usr/bin/env perl
#############################################################################
#
# cdfs-extract
#
#############################################################################
#
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
#
#############################################################################
#
# version 3.0   2013/11/01  PR
# version 2.93  2013/09/27  PR
# version 2.92  2013/07/28  PR
# version 2.91  2012/04/27  PR
# version 2.9   2011/09/22  PR
# version 2.8   2011/08/03  PR
# version 2.71  2011/06/21  PR
# version 2.7   2011/06/06  PR
# version 2.6   2011/05/12  PR
# version 2.52  2011/05/11  PR
# version 2.51  2011/04/14  PR
# version 2.5   2011/04/06  PR
# version 2.4   2010/10/21  PR
# version 2.31  2010/09/29  PR
# version 2.3   2010/08/30  PR
# version 2.2   2010/08/24  PR
# version 2.1   2010/07/21  PR
# version 2.0   2010/06/18  PR
# version 1.0   2010/04/29  PR
#
##############################################################################
#
# see documentation with "perldoc cdfs-extract.pl", or jump to the end of code
#
#
# TODO: refactor the photometry code in a Photometry object...
#       add sanity checks (look if event.list is there... if the sky.txt format
#                          is correct...)
#       rotate logs instead of overwrite
#       allow/bypass checks for srcs without bkgs


BEGIN {
    # purge @INC from sas stuff
    @INC = grep { ! /xmmsas/ } @INC;
}

use strict;
use warnings;
use Getopt::Long;
use Parallel::ForkManager;
use File::Temp qw/tempdir tempfile/;
use Pod::Usage;
use PDL;

use Image;
use PSF;
use EvtFiles;
use Astro::WCS2;




# this variables are visible in all the file
my $maindir = $ENV{PWD};  # save current dir to come back later
my $productdir = "$maindir/Products/";
my $LOG;
my $maxproc = 8;
my $keep_tempfiles = 0;
my $dontdo_resp = '';
my $do_bkgresp = '';
my $noclobber = '';
my $exclude = .40;  # meaning of exclude will be reversed after calling GetOptions()
my $mosix = 0;
my $DB_UPDATE = 0;

# signals
my $ALL_OK = 1;
my $NO_BKG = 2;
my $CHIP_NA = -99;
my $CHIP_OK = 0;
my $CHIP_DIFFERENT = 1;
my $BKG_NOT_SPECIFIED = -1;

main();


sub main {
    print <<LICENSE;
This is cdfs-extract, version 3.0.

Copyright (C) 2010-2013 Piero Ranalli
This program comes with ABSOLUTELY NO WARRANTY.
This is free software, and you are welcome to redistribute it
under certain conditions; for details type:   ./cdfs-extract.pl --help

LICENSE

    my $checkonly = '';
    my $photometry = '';
    my $photoerror = '';
    my $photoerrorlev = 90;
    my $energy = 2; # keV
    my $help = '';
    my $evtlist = 'event.list';
    GetOptions('check'           => \$checkonly,
	       'bkgresp'         => \$do_bkgresp,
	       'noresp'          => \$dontdo_resp,
	       'noclobber'       => \$noclobber,
	       'photometry'      => \$photometry,
	       'photoerror'      => \$photoerror,
	       'photoerrorlev=f' => \$photoerrorlev,
	       'energy=f'        => \$energy,
	       'keep_tempfiles'  => \$keep_tempfiles,
	       'exclude=f'       => \$exclude,
	       'help'            => \$help,
	       'evtlist=s'       => \$evtlist,
	       'maxproc=i'       => \$maxproc,
	       'mosix'           => \$mosix,
	       'dbupdate!'       => \$DB_UPDATE,
	      )
	or pod2usage(2);

    if ($help) {
	pod2usage( -verbose=>2 );
	exit;
    }

    if ($DB_UPDATE) {
	print_error("dbupdate is disabled in this version");
	$DB_UPDATE=0;
    }

    # override maxproc when running inside debugger
    $maxproc = 0 if defined &DB::DB;


    # meaning of --exclude option was reversed starting from version 3.0
    # but the internals still use the old convention
    $exclude = 1 - $exclude;


    unless( defined($ENV{SAS_DIR}) ) {
	die  "Please initialize the XMM SAS before extracting spectra.\n";
    }

    if ($mosix) {
	print("MOSIX clustering is enabled.\n");
    } else {
	print "\n";
    }

    if ($keep_tempfiles) {
	$keep_tempfiles = 1;
    }

#    rotatelog("cdfs-extract.log");
    open($LOG, "> cdfs-extract.log");

    my $behr;
    if ($photoerror) {
	eval { require Astro::BEHR };
	if ($@) {
	    printnlog('Cannot find the Astro::BEHR module.  Errors on net counts and rates will not\nbe calculated.\n\n');
	    $photoerror=0;
	} else {
	    $behr = Astro::BEHR->new({CHECK=>1});

	    unless(ref($behr) eq 'Astro::BEHR') { # correctly initialized
		$photoerror=0;
	    }
	}
    }

    # $events is an object of class EvtFiles
    my $events = EvtFiles->new({evtfile=>$evtlist});
    my %srcs = read_sources("sky.txt");
    my %bkgs = read_sources("bgsky.txt");

    # $srcs->{SRC} is a hash with the following structure:
    #    (  ID1 => [ 
    #                [TYPE, RA, DEC, RADIUS<, RADIUS2><, CAMERA>],
    #                [TYPE, RA, DEC, RADIUS<, RADIUS2><, CAMERA>],
    #                ...
    #              ],
    #   ,   ID2 => [ [...], [...] ] )
    # where radius2 is optional (expected only if TYPE==annulus), and
    # only the first type/ra/dec/radius element is needed. Other elements
    # allow to have sources spread in different places (i.e. stacked spectra).
    # However, photometry does not work with stacked spectra, as the PSF
    # routines only recognise the first element.
    # CAMERA for the moment applies only to bkgs
    #
    # same for %bkgs




    unless (-d $productdir) {
	system("mkdir -p $productdir");
    }

    my $good_ones = do_checks($photometry,$photoerror,$photoerrorlev,
			      \%srcs,\%bkgs,$events,
			      $energy);

    if ($checkonly or $photometry) {
	exit;
    }

    do_extractions($good_ones,\%srcs,\%bkgs,$events);

    close($LOG);
}




sub do_checks {

    if( $DB_UPDATE ) {
	require SQL;
    }

    my ($dophot,$doerrors,$errorlev,$srcs,$bkgs,$events,$energy) = @_;
                                  # hash  hash  object

    my (%deadsrc,%deadbkg,%chipstat,%src_cts,%bkg_cts,%src_avg_exp,
	%bkg_avg_exp,%src_area,%bkg_area,%psffrac,%net_cts,%net_lowerr,%net_higherr,
	%rate,%medianrate,%rate_lowerr,%rate_higherr);

    # using more than 4 cores here is worse, as the bottleneck is rather
    # the disc I/O than cpu time
    my $procnum = $maxproc > 4 && ! $doerrors ? 4 : $maxproc;
    my %talkfiles; my $pm = new Parallel::ForkManager($procnum);

    my $STOPPED_code = 99;

    $pm->run_on_finish( sub { my ($pid, $exit_code, $child_id) = @_;

			      if (defined($exit_code) and
				  $exit_code == $STOPPED_code) {
				  return;
			      }

			   # when child exits, read what it had to say
			   # and put in parent's memory
			   open( my $fileh, " < $talkfiles{$child_id}");
			   chomp($deadsrc{$child_id} = <$fileh>);
			   chomp($deadbkg{$child_id} = <$fileh>);
			   chomp($chipstat{$child_id} = <$fileh>);
			   if ($dophot) {
			       chomp($src_cts{$child_id} = <$fileh>);
			       chomp($bkg_cts{$child_id} = <$fileh>);
			       chomp($src_avg_exp{$child_id} = <$fileh>);
			       chomp($bkg_avg_exp{$child_id} = <$fileh>);
			       chomp($src_area{$child_id} = <$fileh>);
			       chomp($bkg_area{$child_id} = <$fileh>);
			       chomp($psffrac{$child_id} = <$fileh>);
			       chomp($net_cts{$child_id} = <$fileh>);
			       chomp($net_lowerr{$child_id} = <$fileh>);
			       chomp($net_higherr{$child_id} = <$fileh>);
			       chomp($rate{$child_id} = <$fileh>);
			       chomp($medianrate{$child_id} = <$fileh>);
			       chomp($rate_lowerr{$child_id} = <$fileh>);
			       chomp($rate_higherr{$child_id} = <$fileh>);
			   }
			   close($fileh);
			   system("unlink $talkfiles{$child_id}")
			       unless $keep_tempfiles;
		       }
		      );


    for my $evtfile ($events->evts) {
	for my $sourceID (keys %$srcs) {

	    my $child_id = unique_id($sourceID,$evtfile);

	    # open file for child->parent talk
	    # this file will be open for both parent and child,
	    # will be written and closed by child,
	    # re-opened and read and finally deleted by parent.
	    (my $fileh,$talkfiles{$child_id}) = tempfile(
				 'tmp-cdfsextrtalkXXXX',
				 DIR=>'/tmp', UNLINK=>0 );

	    ### begin parallel cycle ###
	    my $pid = $pm->start($child_id) and next;

	    # CHILD FROM HERE

	    unless (-e $events->sumsas($evtfile)) {
		die "Could not find *SUM.SAS or ccf.cif for event file $evtfile.\n";
	    }

	    # create and move to temporary directory
	    my $tmpdir = tempdir( "cdfs-extract-$pid-XXXX", DIR=>'/tmp', CLEANUP=>0 );
	    chdir($tmpdir);

	    # open DBI connection
	    my $sql;
	    if ($DB_UPDATE) {
		$sql = SQL->new;
	    }

	    #my ($src_dead_frac,$src_avg_exp,$src_mask,
	    #	$bkg_dead_frac,$bkg_avg_exp,$bkg_mask) = check_area(
            #    $srcs->{$sourceID},$bkgs->{$sourceID},$maps->{$evtfile} );

	    my $expmap = Image->new( $srcs->{$sourceID},$bkgs->{$sourceID},
				     $events->expmap($evtfile) );
	    unless ($expmap) {  # could not create object, probably because could
		# not read expmap/image, and this may happen if the cameras do not match
		$pm->finish($STOPPED_code);
		next; # in the case weŕe running with --maxproc=0
	    }
	    my ($src_dead_frac,$bkg_dead_frac) = $expmap->dead_fractions;

	    my ($chipstatus,$chipsrc,$chipbkg);
	    if (defined($bkg_dead_frac)) {

		if (($src_dead_frac<$exclude) and ($bkg_dead_frac<$exclude)) {
		    ($chipstatus,$chipsrc,$chipbkg) =
			check_samechips($srcs->{$sourceID},$bkgs->{$sourceID},
					$evtfile,$events->sumsas($evtfile) );
		} else {
		    ($chipstatus,$chipsrc,$chipbkg) =
			($CHIP_NA,$CHIP_NA,$CHIP_NA);
		}

	    } else {
		($chipstatus,$chipsrc,$chipbkg) =
		    ($CHIP_NA,$CHIP_NA,$CHIP_NA);
	        $bkg_dead_frac = $BKG_NOT_SPECIFIED;
	    }

	    # update database
	    if( $DB_UPDATE ) {
		my ($obsid,$instrume) = get_obsid_camera($evtfile);
		$sql->add('spectralchecks',
		          { SRCID=>$sourceID,
			    OBS_ID=>$obsid,
			    INSTRUME=>$instrume,
			    DEADAREA=>$src_dead_frac,
			    EXTRACTED=>($src_dead_frac<$exclude),
			    CHIP=>$chipsrc,
			    BKG_EXTRACTED=>
			    (($src_dead_frac<$exclude) and 
			     ($bkg_dead_frac<$exclude)),
			    BKG_CHIP=>$chipbkg,
			    BKG_DEADAREA=>$bkg_dead_frac,
			    DIFFERENTCHIPS=>$chipstatus
			  });
	    }


	    my $src_cts=my $bkg_cts=my $src_avg_exp=my $bkg_avg_exp=
		my $src_area=my $bkg_area=my $psffrac=my $net_cts=my $net_lowerr=my $net_higherr=
		    my $rate=my $medianrate=my $rate_lowerr=my $rate_higherr=0;
	    if ($dophot
		and $bkg_dead_frac != $BKG_NOT_SPECIFIED # don't do photometry for sources without bkg
		and $src_dead_frac<$exclude
		and $bkg_dead_frac<$exclude) {
		my $image = Image->new( $srcs->{$sourceID},$bkgs->{$sourceID},
					$events->img($evtfile) );
		unless ($image) {  # could not create object, probably because could
		    # not read expmap/image, and this may happen if the cameras do not match
		    $pm->finish($STOPPED_code);
		}
    		$image->masks( $expmap->exposed_masks );

		#print("$sourceID $evtfile\n");
		($src_cts,$bkg_cts) = $image->masked_sums;
		($src_avg_exp,$bkg_avg_exp) = $expmap->masked_averages;
		($src_area,$bkg_area) = $expmap->areas;

		# rate is always classical
		$rate = $src_cts/$src_avg_exp -
		    $bkg_cts/$bkg_avg_exp/$bkg_area*$src_area;

		if ($doerrors) {
		    # Bayesian posterior median of net counts and rate, with errors
		    ($net_cts,$net_lowerr,$net_higherr) = 
			get_errors_from_BEHR(
					     $src_cts, $bkg_cts, 
					     $bkg_area/$src_area*$bkg_avg_exp/$src_avg_exp,
					     $child_id,
					     $errorlev
					    );
		    $medianrate = $net_cts / $src_avg_exp;
		    $rate_lowerr = $net_lowerr / $src_avg_exp;
		    $rate_higherr = $net_higherr / $src_avg_exp;
		} else {
		    # 'classical' net counts and no rate errors
		    $net_cts = $src_cts - $bkg_cts/$bkg_area*$src_area;
		    $net_lowerr = $net_higherr = $rate_lowerr = $rate_higherr = $medianrate = 0;
		}

		my ($off_axis,$rotation) = $image->calc_offaxis;
		my $psf = PSF->new($off_axis,$energy);
		$psf->rotatepsf($rotation);
		$psffrac = $psf->eef($expmap);
		$rate /= $psffrac;
		$medianrate /= $psffrac;
		$rate_lowerr /= $psffrac;
		$rate_higherr /= $psffrac;

	    }

	    print($fileh $src_dead_frac."\n");
	    print($fileh $bkg_dead_frac."\n");
	    print($fileh $chipstatus."\n");
	    if ($dophot) {
		print($fileh $src_cts."\n");
		print($fileh $bkg_cts."\n");
		print($fileh $src_avg_exp."\n");
		print($fileh $bkg_avg_exp."\n");
		print($fileh $src_area."\n");
		print($fileh $bkg_area."\n");
		print($fileh $psffrac."\n");
		print($fileh $net_cts."\n");
		print($fileh $net_lowerr."\n");
		print($fileh $net_higherr."\n");
		print($fileh $rate."\n");
		print($fileh $medianrate."\n");
		print($fileh $rate_lowerr."\n");
		print($fileh $rate_higherr."\n");
	    }
	    close($fileh);

	    chdir($maindir);
	    sas_exec('',"rm -rf $tmpdir");

	    $pm->finish;
	}
    }

    # BACK TO PARENT
    # wait for all child processes before exit
    $pm->wait_all_children;

    my $good_ones = summary_of_checks(\%deadsrc,\%deadbkg,\%chipstat);
    if ($dophot) {
	photometry_report(\%src_cts,\%bkg_cts,\%src_avg_exp,\%bkg_avg_exp,
			  \%src_area,\%bkg_area,\%psffrac,
			  \%net_cts,\%net_lowerr,\%net_higherr,
			  \%rate,\%medianrate,\%rate_lowerr,\%rate_higherr,
			  $errorlev);

    }
    return($good_ones);
}



sub get_errors_from_BEHR {
    require Astro::BEHR;
    my ($ssrc,$sbkg,$sarea,$id,$errorlev) = @_;

    my $behr = Astro::BEHR->new;
    $behr->set($ssrc,$sbkg,$sarea,$ssrc,$sbkg,$sarea,90);

    $behr->{LEV} = $errorlev if ($errorlev);

    $behr->timeout(180);
    my $err = $behr->run;
    if ($err eq 'TIMEOUT') {
	printnlog("Timed out (quadr): $id\n");
        # try again:
        $behr->{ALGO} = 'gibbs';  # much faster
        my $err2 = $behr->run;

        if ($err2 eq 'TIMEOUT') {  # AGAIN?!?
	    printnlog("Timed out (gibbs): $id\n");
	    return(0,0,0);
        }
    }

    my $vals = $behr->get('S');
    return(@$vals[2..4]);
}




sub unique_id {
    my ($srcid,$evt) = @_;
    return( "${srcid}\@${evt}" );
}


sub do_extractions {
    my ($good,$srcs,$bkgs,$events,$maps) = @_;
      # hash  hash  hash  array hash (refs)

    my $pm = new Parallel::ForkManager($maxproc);

    for my $evtfile ($events->evts) {

	my $SAS = getSASenvironment( $events->ccfcif($evtfile),
				     $events->sumsas($evtfile) );

	my ($obsid,$camera) = get_obsid_camera( $evtfile );


	print("$evtfile\n");

	for my $sourceID (keys %$srcs) {

	    my $id = unique_id($sourceID,$evtfile);
	    next unless $$good{$id};

	    my $pid = $pm->start and next;

	    my $tmpdir = tempdir( "cdfs-extract-$pid-XXXX", DIR=>'/tmp', CLEANUP=>0 );
	    chdir($tmpdir);

	    # $good also contains info about bkg extraction
	    my $bkg = $$good{$id} == $ALL_OK ? $$bkgs{$sourceID} : undef;
	    extractor( $$srcs{$sourceID},$bkg,$evtfile,
		       $sourceID,$obsid,$camera,$SAS );


	    # go back to maindir and remove tmpdir
	    chdir($maindir);
	    sas_exec('',"rmdir $tmpdir");

	    $pm->finish;

	}
    }

    # wait for all child processes before exit
    $pm->wait_all_children;
}


sub extractor {

    my ($src,$bkg,$evt,$sourceID,$obsid,$camera,$SAS) = @_;

    my ($ok,@xy_rad) = radec2xy( $src, $evt, '+1' );
    unless ($ok) {
	chdir($maindir);
	return;
    }
    # it that was ok, we know that there is at list one source
    # in xy_rad fitting our camera

    # check if there is a background for our source
    my @bgxy_rad;
    my $withbackground = 1;

    if (defined( $bkg )) {
	($ok, @bgxy_rad) = radec2xy( $bkg, $evt, '+1' );
	unless( $ok ) {
	    $withbackground = 0;
	    print_error("Bkg spectra will not be extract for source $sourceID.\n");
	}
    } else {			# bkg not defined in bgsky.txt
	$withbackground = 0;
	print($LOG "No background regions for source $sourceID.\n Bkg spectra will not be extracted.\n");
    }

    ($ok, my $srcspec) = 
	extract_spectrum( $sourceID, $obsid, $SAS, $camera, $evt, \@xy_rad, '' );
    unless ($ok) {
	chdir($maindir);
	return;
    }

    my $bkgspec;
    if ($withbackground) {
	($ok, $bkgspec) = extract_spectrum( 
	    $sourceID, $obsid, $SAS, $camera, $evt, \@bgxy_rad, 'bg' );
	unless ($ok) {
	    chdir($maindir);
	    return;
	}
    }

    unless ($dontdo_resp) {
	calc_response( $SAS, $evt, $srcspec );

	if ($do_bkgresp and $withbackground) {
	    calc_response( $SAS, $evt, $bkgspec );
	}
    }
}









sub read_sources {
    # modified to allow multiple areas per source (see man perlol)
    # and camera-specific backgrounds

    use Scalar::Util qw/looks_like_number/;

    my $filename = shift;
    my %src;

    open(my $SRC, '<', $filename);
    while (my $f=<$SRC>) {

	next if ($f =~ m/^\s*\#/); # jump comments
	next if ($f =~ m/^\s*$/);  # and empty lines
	$f =~ s/\#.*$//;  # remove comments

	# look if region should be excluded or included
	# (and be smart enough not to look into comments)
	my $exclude = 0;
	my ($id,$ra,$dec,$radius1,$radius2,$camera);
	if ($f =~ m/exclude/) {
	    ($id,undef,$ra,$dec,$radius1,$radius2,$camera) = split(/\s+/, trim($f));
	    $exclude = 1;
	    # sources to be excluded have negative radii
	    $radius1 = -$radius1;
	    if (defined($radius2) and looks_like_number($radius2)) {
		$radius2 = -$radius2;
	    }
	} else {
	    ($id,$ra,$dec,$radius1,$radius2,$camera) = split(/\s+/, trim($f));
	}

	unless (defined($radius1)) {
	    die "radius not defined in $filename\n"
	}

	arcsec2physpix(\$radius1);
	arcsec2physpix(\$radius2) if looks_like_number($radius2);


	if (defined($camera)) { # sure it's an annulus, whatever $camera contains
	    push(@{ $src{$id} }, ['annulus',$ra,$dec,$radius1,$radius2,uc($camera)] );

	# or, it may or may not have a camera specified. We have to check.
	} elsif (defined($radius2) and not (
		     uc($radius2) eq 'EMOS1' or
		     uc($radius2) eq 'EMOS2' or
		     uc($radius2) eq 'EPN')) {
	    # then it's an annulus without camera
	    push(@{ $src{$id} }, ['annulus',$ra,$dec,$radius1,$radius2] );

	# otherwise, it's a circle; with or without a defined camera
	} else {
	    push(@{ $src{$id} }, ['circle',$ra,$dec,$radius1,
				  defined($radius2) ? uc($radius2) : undef
				 ] );
	}
    }
    close($SRC);
    return(%src);
}


sub arcsec2physpix {
    my $arcsec = shift;

    # the following is the XMM convention
    $$arcsec*=20;
}

sub getSASenvironment {
  my ($ccf,$sum) = @_;

  return "SAS_ODF=$sum SAS_CCF=$ccf ";
}


sub get_obsid_camera {

    my $evtfile = shift;
    my $cam;

    # now uses PDL instead of ftools because of (rare) race conditions when
    # different processes tried to do something at the same time

    my $hdr=rfits($evtfile.'[1]',{data=>0});
    $cam = $hdr->{INSTRUME};

    unless ($cam eq 'EMOS1'  or  $cam eq 'EMOS2'  or  $cam eq 'EPN') {
      print_error("Camera not regognised for event file $evtfile\n");
      print_error("Skipping processing for this event file.\n");
      return(0);
    }

    return($hdr->{OBS_ID},$cam);
}





sub print_error {
  my $s = shift;
  print($LOG $s);
  unless ($s =~ m/^\d+ STDERR:\n$/) {
      # only print to screen if there is actually any message
      print(STDERR $s);
  }
}


sub extract_spectrum { # ($SAS, $camera, $evtfile, \@xy_rad, '' or 'bg' );

  my ($src, $obsid, $SAS, $camera, $evtfile, $xy_rad, $bkgflag) = @_;

  # chooses and executes the template command


  # camera: from EMOS1 to mos1, ..., EPN to pn
  $camera =~ s/^E//;
  $camera = lc($camera);
  my $tmpl = "$maindir/gsp.$camera.tmpl";


  # decide spectrum name
  my $specname = "$productdir/$src-$obsid-$camera";
  unless ($bkgflag eq '') {
    $specname = $specname."-$bkgflag";
  }
  $specname = $specname.".pha";

  my $ok = 1;
  unless( $noclobber and -e $specname) {

    # source: from definition, to expression in SAS syntex
    my $expr = src2expr($xy_rad);

    # read template
    local $/;			#SLURP MODE!
    open(my $TMPL, "< $tmpl");
    my $cmd = <$TMPL>;
    close($TMPL);

    # modify
    $cmd =~ s/EVT/$evtfile/g;
    $cmd =~ s/NOME/$specname/g;
    $cmd =~ s/EXPR/$expr/g;

    my $backscal;
    my $i=0;
    do {
	$ok = sas_exec( $SAS, $cmd, {MOSIX=>1} );
	if (! $ok) { return($ok); }  # exit for error
	$i++;

	# now check for backscal value
	# and repeat extraction if backscal=0;
	# 1e-3 here is just a small value "indistinguishable from 0"
	my $h =rfitshdr($specname);
	$backscal = $h->{BACKSCAL};
	if ($backscal < 1e-3) {
	    printnlog("Repeating spectral extraction of $specname because of BACKSCAL=0.\n");
	    unlink($specname);
	}
    } until (($backscal > 1e-3) or ($i==3));
    
    if ($backscal <= 1e-3) {
	printnlog("ABORTED spectral extraction of $specname because of BACKSCAL=0 happening 3 times in a row.\n");
	$ok = 0;
    }

  }

  return($ok,$specname);
}


sub src2expr {
    # source: from definition, to expression in SAS syntax

    my $xy_rad = shift;   # is an array (of array) ref
    my $firstexpr = 1;
    my $expr = '(';
    my $firstexclude = 1;
    my $excludexpr = '(';

    for my $i (0..$#$xy_rad) {

	my $type = shift(@{$$xy_rad[$i]});
	my $x  = shift(@{$$xy_rad[$i]});
	my $y  = shift(@{$$xy_rad[$i]});
	my $r1 = shift(@{$$xy_rad[$i]});
	my $r2 = shift(@{$$xy_rad[$i]});

	if ($r1>0) {   # include
	    if($firstexpr) {
		$firstexpr = 0; # and do not add || to the expression
	    } else {
		$expr .= ' || ';
	    }
	} else { # exclude
	    if($firstexclude) {
		$firstexclude = 0; # and do not add &&
	    } else {
		$excludexpr .= ' && ';
	    }
	}




	if ($type eq 'circle') { 
	    if ($r1>0) {
		$expr .= "(X,Y) IN CIRCLE($x,$y,$r1)";
	    } else {
		$r1 = abs($r1);
		$excludexpr .= "! (X,Y) IN CIRCLE($x,$y,$r1)";
	    }
	} elsif ($type eq 'annulus') {
	    if ($r1>0) {
		$expr .= "(X,Y) IN ANNULUS($x,$y,$r1,$r2)";
	    } else {
		$r1 = abs($r1);
		$r2 = abs($r2);
		$excludexpr .= "! (X,Y) IN ANNULUS($x,$y,$r1,$r2)";
	    }
	} else {
	    print_error("Something wrong with source shape, which should not have happened.\nPlease report this error message along with sky.txt and bgsky.txt\n");
	    return(0);
	}

    }

    $expr .= ')';
    unless($firstexclude) { # only add exclude part if it was used
	$expr .= ' && '.$excludexpr.')';
    }
    return($expr);
}


sub calc_response { #( $SAS, $evtfile, $srcspec );
  my ( $SAS, $evtfile, $srcspec ) = @_;

  my $rmf = my $arf = $srcspec;
  $rmf =~ s/\.pha$/.rmf/;
  $arf =~ s/\.pha$/.arf/;

  unless( $noclobber and -e $rmf and -e $arf ) {

    local $/;			#SLURP!
    open(my $TMPL, "< $maindir/response.tmpl");
    my $cmd = <$TMPL>;
    close($TMPL);

    # replace keywords in template
    $cmd =~ s/SPEC/$srcspec/g;
    $cmd =~ s/RMF/$rmf/g;
    $cmd =~ s/ARF/$arf/g;
    $cmd =~ s/EVT/$evtfile/g;

    sas_exec($SAS,$cmd,{MOSIX=>1});
  }
}


sub sas_exec {
  my ($SAS,$script,$opt) = @_;

  # split template in commands
  my @cmds = split(/\n/,$script);

  # execute commands
  for my $cmd (@cmds) {


      my ($errfileh,$errfile) = tempfile(
	  "$$.cdfsextr.stderrXXXX",
	  DIR=>'/tmp', UNLINK=>0 );

      next if ($cmd =~ m/^\s*(\#.*)?$/);  # empty line or comment

      if ($mosix and
	  defined($opt) and exists($$opt{MOSIX}) and $$opt{MOSIX})
	  {
	  adapt4mosix(\$cmd);
      }

      print($LOG "$$ COMMAND: $SAS $cmd\n");
      my $out = `$SAS  $cmd   2> $errfile`;

      # check exit status and print output and/or errors
      if ($?) {
	print_error("$$: something wrong happened. Please check for errors in the log.\n");
      }
      print($LOG "$$ STDOUT:\n$out");

      close($errfileh);
      local $/;   #SLURP
      open(my $ERR, "< $errfile");
      my $err = <$ERR>;
      close($ERR);
      print_error("$$ STDERR:\n$err");
      unlink("$errfile");

      # exit now if any error
      if ($?) {
	return(0);
      }
  }

  return(1); # ok
}


sub adapt4mosix {
    my $comref = shift;
    my @mosix_allowed = qw/backscale rmfgen arfgen/;

    for my $ma (@mosix_allowed) {
	if ($$comref =~ m/^\s*$ma/) {
	    $$comref = 'mosenv -j1,3 -L '.$$comref;
	    return;
	}
    }
    return;
}


sub trim {
  my $s = shift;
  $s =~ s/^\s*//;
  $s =~ s/\s*$//;
  return($s);
}



sub phys2img {

    my ($phys,$cdelt) = @_;

    # phys/20 -> in arcsec   ( 20 is actually 1/abs($evt->hdr->{REFXCDLT}) )
    # /3600 -> in deg
    # /CDELT1 -> in image pixels
    return( $phys/20/3600/abs($cdelt) );
}




sub radec2xy {   # using sky2xy seems not to work in parallel
     # converts from ra,dec to x,y using PDL
     # now accepts multiple positions for same source (i.e. stacked
     # spectra)


    my ($src,$evt,$ext) = @_;	# $src is an array ref

    my ($wcs_sub,$header);
    if ($ext eq '+0') {
	$ext = '[0]';
	$header = rfits($evt.$ext,{data=>0});
	$wcs_sub = \&wcstransfinv;
    } elsif ($ext eq '+1') {
	$ext = '[1]';
	$header = rfits($evt.$ext,{data=>0});
	$wcs_sub = \&wcs_evt_transfinv;
    } else {
	print_error("Strange FITS extension \#$ext asked for file $evt -- could not read.\n");
	return(0);
    }

    my $cam = $header->{INSTRUME};

    my @new_lol;
    for my $i (0..$#$src) {
	my ($type,$ra,$dec,@radii) = @{$$src[$i]}  ;

	# for how read_sources works, the camera (if it's there) is
	# always the last item in the list.
	my $askedcam = $radii[$#radii];
	if (defined($askedcam) and (
		$askedcam eq 'EMOS1' or
		$askedcam eq 'EMOS2' or
		$askedcam eq 'EPN')) {

	    # a camera is specified, so let's cycle unless
	    # it is the good one
	    next unless ($askedcam eq $cam);

	    # remove from the radii
	    pop(@radii);
	}

	my ($x,$y);
	($x,$y) = &$wcs_sub($header,$ra,$dec);

	push(@new_lol, [ $type,$x,$y,@radii ]);
    }

    # check if the new lol actually contains any item, and return
    if ($#new_lol >= 0) {
	return(1,@new_lol);
    } else {
	return(0);
    }
}




sub printnlog {
    my $s = shift;
    print($s);
    print($LOG $s);
}

sub printfnlog {
    my ($f,@s) = @_;
    printf($f,@s);
    printf($LOG $f,@s);
}

sub summary_of_checks {

    my ($deadsrc,$deadbkg,$diffchip) = @_;  # hash refs
    my %goodness;

    my @k = sort(keys(%$deadsrc));

    my (@srcmsg1,@bkgmsg1,@srcmsg2,@bkgmsg2,@chipmsg);


    # SRC
    printnlog("Summary of checks.\n\n");

    for my $key (@k) {
	if ($$deadsrc{$key} >= $exclude) {
	    push(@srcmsg1,$key);
	} else {
	    $goodness{$key} = $ALL_OK;
	}
    }
    if ($#srcmsg1 >= 0) {
	printnlog("The following combinations of source and event file have dead areas >= 90%, thus\nno spectra will be extracted for them:\n");
	for my $key (@srcmsg1) { 
	    printfnlog("%s (%4.1f %%)\n",$key,100*$$deadsrc{$key}); 
	}
    } else {
	printnlog("There are no combinations of source and event file with dead areas >= 90%.\n");
    }



    for my $key (@k) {
	if ($$deadsrc{$key} >= .5 and $$deadsrc{$key} < $exclude) {
	    push(@srcmsg2,$key);
	}
    }
    if ($#srcmsg2 >= 0) {
	printnlog("\nHere are warnings: the following have source dead areas >= 50\%, but spectra\nwill still be extracted:\n");
	for my $key (@srcmsg2) { 
	    printfnlog("%s (%4.1f %%)\n",$key,100*$$deadsrc{$key}); 
	}
    } else {
	printnlog("There are no combinations of source and event file with dead areas >= 50%.\n");
    }

    # BKG

    for my $key (@k) {
	if ($$deadbkg{$key} != $BKG_NOT_SPECIFIED and
	    $$deadbkg{$key} >= $exclude and
	    $goodness{$key}) {

	    push(@bkgmsg1,$key);
	    $goodness{$key} = $NO_BKG;
	}
    }
    if ($#bkgmsg1 >= 0) {
	printnlog("\n\nThe following combination of background and event file have dead areas >= 90%,\n thus no bkg spectra will be extracted for them:\n");
	for my $key (@bkgmsg1) { 
	    printfnlog("%s (%4.1f %%)\n",$key,100*$$deadbkg{$key});
	}
    } else {
	printnlog("There are no combinations of background and event file with dead areas >= 90%.\n");
    }


    for my $key (@k) {
	if ($$deadbkg{$key} >= .5 and $$deadbkg{$key} < $exclude and $goodness{$key}) {
	    push(@bkgmsg2,$key);
	}
    }
    if ($#bkgmsg2 >= 0) {
	printnlog("\nHere are warnings: the following have background dead areas >= 50\%,\n but bkg spectra will still be extracted:\n");
	for my $key (@bkgmsg2) { 
	    printfnlog("%s (%4.1f %%)\n",$key,100*$$deadbkg{$key});
	}
    } else {
	printnlog("There are no combinations of background and event file with dead areas >= 50%.\n");
    }


    # CHIPS
    for my $key (@k) {
	if (($$diffchip{$key} != $CHIP_OK) and ($$diffchip{$key} != $CHIP_NA) 
	    and $goodness{$key}) {
	    push(@chipmsg,$key);
	}
    }
    if ($#chipmsg >= 0) {
	printnlog("\nThe following combinations of source and background fall on different chips:\n");
	for my $key (@chipmsg) { 
	    my $bc = int($$diffchip{$key} / 100);
	    my $sc = $$diffchip{$key}-100*$bc;
	    printnlog("$key (src: chip $sc; bkg: chip $bc)\n"); }
    } else {
	printnlog("There are no combinations of background and event file which fall on\ndifferent chips.\n");
    }



    return(\%goodness);
}



sub check_samechips {
    # checks if src and bkg fall on different chips

    my ($src,$bkg,$evt,$sumsas) = @_;

    my $SAS = "SAS_ODF=$sumsas";
    my (undef,$camera) = get_obsid_camera($evt);
    my $xmmea;
    if ($camera =~ m/EMOS/) {
	$xmmea = '#XMMEA_EM && (PATTERN<=12)';
    } elsif ($camera eq 'EPN') {
	$xmmea = '#XMMEA_EP && (PATTERN<=4) && (FLAG==0)';
    } else {
	die "Could not recognize camera $camera for event file $evt.\n";
    }

    # no need for tempfiles here -- we're already in a tempdir!
    my $histos = 'src-histogram.ds';
    my $histob = 'bkg-histogram.ds';

    # the SAS expression then gets recalculated in extract_spectrum
    # is it worthwhile to put all exprs in wholefile-scoped variables?
    my ($ok, @xy_rad) = radec2xy( $src, $evt, '+1' );
    my $expr = src2expr(\@xy_rad); # SAS expression

    # extract histogram of src
    my $cmd = "evselect table=$evt withhistogramset=yes histogramset=$histos histogramcolumn=CCDNR histogrambinsize=1 expression='$xmmea && $expr'";
    $ok = sas_exec( $SAS, $cmd );

    # now for bkg
    ($ok, @xy_rad) = radec2xy( $bkg, $evt, '+1' );
    $expr = src2expr(\@xy_rad);
    $cmd = "evselect table=$evt withhistogramset=yes histogramset=$histob histogramcolumn=CCDNR histogrambinsize=1 expression='$xmmea && $expr'";
    $ok = sas_exec( $SAS, $cmd );

    # look if most_populated_bin(src)==m_p_b_(bkg) in at least xx% of times
    my $sh = rfits($histos.'[1]');
    my $bh = rfits($histob.'[1]');

    # this could be made more clever, e.g. to look if regions lie on multiple
    # chips...
    my $maxs = $sh->{COUNTS}->maximum_ind;
    my $maxb = $bh->{COUNTS}->maximum_ind;
    if ($sh->{CCDNR}->at($maxs) != $bh->{CCDNR}->at($maxb)) {
	return($sh->{CCDNR}->at($maxs) + 100*$bh->{CCDNR}->at($maxb),
	       $sh->{CCDNR}->at($maxs),
	       $bh->{CCDNR}->at($maxb)
	      );
    } else {
	return ($CHIP_OK,
		$sh->{CCDNR}->at($maxs),
		$bh->{CCDNR}->at($maxb)
	       );
    }
}




sub photometry_report {
    my ($src_cts,$bkg_cts,$src_avg_exp,$bkg_avg_exp,   # all hash refs except $errorlev
	$src_area,$bkg_area,$psffrac,$net_cts,$net_lowerr,$net_higherr,
	$rate,$medianrate,$rate_lowerr,$rate_higherr,
        $errorlev) = @_;

    open(my $PHOT, "> $maindir/photometry.dat");

    $errorlev = 90 unless ($errorlev);

    pr_printf($PHOT,"%s","# Aperture photometry (per obsid per camera):\n(only the rates are corrected for the PSF fraction)\n(cameras: 1=M1, 2=M2, 3=PN)\nErrors are at ${errorlev}% level\n\nsrc    obsid   camera  src_cts bkg_cts src_exp bkg_exp src_area bkg_area psffrac net_cts net_lowerr net_higherr rate medianrate rate_lowerr rate_higherr\n");

    for my $k (sort(keys %$src_cts)) {

	# old, short format
# 	pr_printf($PHOT,"%5s %10s %2s %5i %5i %8.1f %8.1f %5.1f %5.1f %5.3f %8.1e %8.1e\n",
# 	      pr_getsrc($k),pr_getobsid($k),pr_getcamera($k),$$src_cts{$k},
# 	      $$bkg_cts{$k},
# 	      $$src_avg_exp{$k},$$bkg_avg_exp{$k},$$src_area{$k},$$bkg_area{$k},
# 	      $$psffrac{$k},$$rate{$k},$$rate_err{$k});

	# new, longer format with errors
	pr_printf($PHOT,
		  "%5s %10s %2s %5i %5i %8.1f %8.1f %5.1f %5.1f %5.3f %7.1f %7.1f %7.1f %8.1e %8.1e %8.1e %8.1e\n",
		  pr_getsrc($k),pr_getobsid($k),pr_getcamera($k),$$src_cts{$k},
		  $$bkg_cts{$k},
		  $$src_avg_exp{$k},$$bkg_avg_exp{$k},$$src_area{$k},$$bkg_area{$k},
		  $$psffrac{$k},$$net_cts{$k},$$net_lowerr{$k},$$net_higherr{$k},
		  $$rate{$k},$$medianrate{$k},$$rate_lowerr{$k},$$rate_higherr{$k});

    }

    close($PHOT);

    # TODO: per camera, and totals
#     print("Aperture photometry (per camera and total):\n\nsrc    obsid   camera  src_cts bkg_cts src_exp bkg_exp src_area bkg_area psffrac rate rate_err\n");

#     for my $k (pr_srcs(sort(keys %$src_cts))) {

# 	for my $cam (qw/M1 M2 PN/) {
# 	    print("%5s %10s %2s %5i %5i %5.1f %5.1f %5.1f %5.1f %4.1f %8.1e %8.1e\n",
# 	      pr_getsrc($k),pr_getobsid($k),pr_getcamera($k),$$src_cts($k),
# 	      $$bkg_cts($k),
# 	      $$src_avg_exp($k),$$bkg_avg_exp($k),$$src_area($k),$$bkg_area($k),
# 	      $$psffrac($k),$$rate($k),$$rate_err($k));
# 	}
#     }


    sub pr_printf {
	my $fh = shift;
	printf(@_);      # print to screen
	printf($fh @_);  # print to file
    }

    sub pr_getsrc {
	my $k = shift;
	(my $src) = ($k =~ m/^(.+)\@/);
	return($src);
    }
    sub pr_getobsid {
	my $k = shift;
	(my $obsid) = ($k =~ m|^.+\@.*/(\d+)/|);
	return($obsid);
    }
    sub pr_getcamera {
	my $k = shift;
	if ($k =~ m|/m1.?_|) {
	    return(1); #'M1');
	} elsif ($k =~ m|/m2.?_|) {
	    return(2); #'M2');
	} elsif ($k =~ m|/pn_|) {
	    return(3); #'PN');
	}
    }
}











__END__




=head1 NAME

cdfs-extract -- A program to extract spectra and aperture photometry
from multiple sources in multiple XMM-Newton observations.

=head1 ABOUT

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
al. 2013, A\&A 555, A42.
 http://cdsads.u-strasbg.fr/abs/2013A%26A...555A..42R

cdfs-extract is written in the Perl programming language, using the
Perl Data Language modules.


=head1 SYNOPSIS

./cdfs-extract.pl [options]

(for help: ./cdfs-extract.pl -h)

=head1 OPTIONS

=over 3

=item --check

Only runs checks the source and background coverage and if the source
and bkg fall on the same chips. Produces a screen log of the
result. (The checks results are also written in
cdfs-extract.log. Also, a useful unix command is tee:
 ./cdfs-extract.pl --check | tee checks.txt

This will send the screen output also to the file checks.txt).

=item --dbupdate

(experimental feature) Updates the mysql database 'spectralchecks'
with the results of the checks.

=item --bkgresp

Also calculates the RMFs and ARFs for the background spectra.

=item --exclude=XX

Exclude a source/file combination from spectra extraction, if the exposed
area is less than XX (default: 0.40).

=item --evtlist=filename

Use filename instead of the default event.list.  (Especially useful to
do photometry in different bands, or for tests on shorter event lists).

=item --photometry [--energy=XX]

Runs (only) checks and aperture photometry. The energy (in keV;
optional, default is 2 keV) is used to calculate the PSF fraction
covered by the aperture.

=item --photoerror[=XX]

Calculates the errors on net counts and net rate when doing
photometry, using the external program BEHR. Also, switches the
calculation of the net counts and rate from the "classical" definition
to the median of the posterior distribution from BEHR.  If the latter
program is not available, prints a warning and continues without
calculating the errors.

The optional parameter specifies the per cent probability level
(default=90).

=item --noclobber

Does not overwrite files. May be useful if you want to add one
observation to the event.list and only create the relevant
spectra/responses, without recreating any already existant
file. Another use may be, together with --bkgresp, to extract only
background responses for spectra which have already been extracted.

=item --mosix

Enables the execution of selected commands (currently: backscale,
rmfgen, arfgen) through the MOSIX clustering system.

=item --maxproc=X

Uses a maximum number of X processes for the extraction of spectra
(default: 8). For maximum speed, X should be the same of the number of
CPU cores.

However, during the check and photometry phases, the bottleneck is
rather the disc I/O, and thus only for these phases the number of
processes is limited to max(X,4).

=item --keep_tempfiles

Keeps the temporary files in /tmp (may be useful for debug).

=item --help

Prints this help message.

=back

=head1 DESCRIPTION

=head2 Prepare the files before you run:

=over 3

=item 1. SOURCES:

prepare a file named sky.txt containing a table with the following info:
    ID_src1  RA_1   DEC_1   RADIUS_1   [RADIUS2_1]   [CAMERA_1]
    ID_src2  RA_2   DEC_2   RADIUS_2   [RADIUS2_2]   [CAMERA_2]
    ...      ...    ...     ...        [...]         [...]

where RA and DEC are in degrees, and the RADII are in arcsec.

Notes:

=over 4

=item * 
The columns marked with the square brackets above are optional.

=item * 
If RADIUS2 is specified, then the program assumes that the
source is an annulus; otherwise, it is a circle.

=item * 
If CAMERA is specified, then the RA, DEC, and radii apply only
to event files from that camera (this is especially useful to choose
camera-dependent backgrounds: see below for bgsky.txt).

=item *
You can add comments in sky.txt by starting a line with
the # character.

=item *

There may be multiple lines with the same ID_src. This means
that the spectrum will that be extracted will have photons coming from
all the regions: use this feature if you want average spectra.

=item *

If you want to exclude any area from the extraction of a
spectrum (say, because there is a nearby source which should not be
included), put the word "exclude" between ID_src and RA_src, as in the
following:
    ID_src1  exclude RA_src99   DEC_src99   RADIUS_src99

thus the area around RA_src99,DEC_src99 will be excluded from the
spectrum of ID_src1.  You can of course exclude multiple positions
(use the same ID_src for all).

=back

=item 2. BACKGROUND:

prepare a file named bgsky.txt containing a table with the same
structure of sky.txt. Use the same ID_src to identify the sources,
but give RA, DEC and radii relative to the background regions.
The same notes for SOURCES apply.

=item 3. EVENT FILES:

prepare a file named event.list which contains the paths of all event
files, exposure maps, images, ccf.cif and *SUM.SAS, one per line as in
the following scheme:
   evtfile1.fits  expmap1.fits img1.fits ccf1.cif xxxxxxSUM.SAS1
   evtfile2.fits  expmap2.fits img2.fits ccf2.cif xxxxxxSUM.SAS2
   ...            ...          ...       ...      ...

Images are optional, since they are only used when doing aperture
photometry. In this case, put NA as a placeholder.

The ccf.cif and *SUM.SAS files are needed to link the event files to
the relevant calibration data.

You may comment out events that you don't want to process, by
prepending the # character.  This program will take care of
recognising the camera (MOS1, MOS2, PN) and use the right command
templates accordingly.

=item 4. COMMAND TEMPLATES:

check that the command templates are ok, in terms of energy, pattern and
flag filtering. The templates are:
 gsp.mos1.tmpl, gsp.mos2.tmpl, gsp.pn.tmpl  (extraction of spectra and calculations of backscale)
 response.tmpl   (generation of rmf, arf)

=item 5. CLEAN SPACE:

all the produced spectra and responses will be put in the Products/
subdirectory. Please check that this does not contain older files
which may get mixed with the new ones and make confusion. Take special
care in this sense if you are going to use --noclobber.

=item 6. RUN THE PROGRAM:

start HEADAS (FTOOLS), SAS, and call ./cdfs-extract.pl with any needed
option.

=item 7. CHECK THE LOG:

a log will be written in cdfs-extract.log. Please check that everything
was processed correctly.

=back

=head2 Checks

The program checks that at least a significant fraction of the source
and background regions have been covered in each observation. Sources
may in fact not be exposed because they are outside of the FOV, or
behind chip gaps, etc.  Warnings are emitted (but spectra are still
extracted) if <50% of the PSF is exposed.  Spectra are not extracted
if <40% of the PSF is exposed (this threshold may be changed with
--exclude=XX). The checks refer to "dead area" fractions which are
defined as 1 - exposed_fraction.

The program also checks if the source and background spectra fall on
the same chip; a warning is shown if this is not the case.

If you just want the checks and need not to proceed further, you can
use the option --check .

=head2 Spectra and responses

Source and background spectra, and source responses (rmf and arf) are
calculated by default.

If you also want background responses, add the option --bkgresp to the
command line. If you have already extracted the spectra and the source
responses, and you just need the background responses, you can add the
option --noclobber to avoid regenerating any already present file.

=head2 Photometry

Prepare the files as above, and call ./cdfs-extract.pl --photometry .
Spectra and responses will not be extracted in this case. Instead, a
report with results from aperture photometry will be printed on
screen, and also saved to the file photometry.dat.

If you want the energy encircled fractions to be calculated at
energies different than the default (2 keV), there is the option
--energy=XX , where XX is in keV.

Aperture photometry is calculated for each observation for each camera
in the following way:

=over 4

=item 1. Only the pixels inside the source region, and which are
active in the exposure map (i.e. are inside the FOV, do not fall in a
CCD gap, have actually been exposed, etc.) are considered. These are
referred to as "active pixels", in the following. The source counts
are extracted from the image, in the same active pixels of the expmap.

=item 2. The source exposure is calculated, by averaging the values of
the expmap in the active pixels.

=item 3. The source area is calculated, as the sum of the active pixel areas.

=item 4. The same is done for the background region.

=item 5a. Unless --photoerror is used: the net counts are

 netcts = src_cts - bkg_cts*src_area/bkg_area
And the net rate is
 rate_nopsf = src_cts/src_exp - bkg_cts/bkg_exp * (src_area/bkg_area)

=item 5b. If --photoerror is used: the net counts are taken as the
median of the (Bayesian) posterior distribution of source intensity,
given the observed source and bkg counts, and the relative areas and
exposures.  The 90% errors on the net counts and net rate are also
calculated. The external program BEHR is used for the calculations
(Park et al. 2006, ApJ 652, 610).

=item 6. The encircled energy fraction is calculated based on the PSF,
as the fraction of flux actually falling on the (source) active pixels. The
corrected rate is thus

 rate = rate_nopsf / encircled_energy_fraction


=back

=head1 WARNINGS AND ERRORS

Some warnings emitted by SAS can usually be considered safe, such as the following:

=over 4

=item ** evselect: warning (NoFilteredEvents), No events have been selected - 
      filtering of the event list resulted in an empty table
      (This warning signals that one source was completely out of the FOV
      in one observation).

=back

Otherwise, if the following error is emitted, please report it along with
a copy of the log, of the sky.txt, bgsky.txt and event.list:

=over 4

=item .....: something wrong happened. Please check for errors in the log.

=back

=head1 VERSION HISTORY

=over 4

=item 3.0  -- 2013/11/01 removed logic to find CDFS ccf.cif and sum.sas: these
                         files should now be stated in the evtlist; sky.txt
                         and bgsky.txt now use arcsecs instead of physical 
                         pixels;
                         evtlist logic refactored in EvtFiles package;
                         dbupdate is disabled;
                         exclusion threshold now set to exposed_area>40% i.e.
                         deadarea<60% of PSF (and the meaning of the --exclude
                         option was reversed)

=item 2.93 -- 2013/9/27  started refactor and cleanup before releasing the code:
                         PSF and Image packages moved to separate files

=item 2.92 -- 2013/7/28  --photoerror=XXX now accepts probability level for BEHR

=item 2.91 -- 2012/4/27  skip photometry of sources without a bkg region

=item 2.9 -- 2011/9/23  reworked to use BEHR.pm; added timeout to BEHR;
corrected a bug in photometry code

=item 2.8 -- 2011/8/3  added checks of BACKSCAL

=item 2.71 -- 2011/6/21  corrected a bug in SQL code

=item 2.7  -- 2011/6/6   added results of spectral checks to mysql database

=item 2.61 -- 2011/6/6   check that SAS is initialized before starting

=item 2.6  -- 2011/5/12  adapted for MOSIX (--mosix option)

=item 2.52 -- 2011/5/11  due to bug appearing only on lafaro, modified
                         getSAS_ODF_CCF() to check the machine and behave
                         accordingly (see details in getSAS_ODF_CCF())

=item 2.51 -- 2011/4/14  added --evtlist=event.list option;
                         corrected a bug which prevent responses from
                         being calculated if bkg were not extracted

=item 2.5 -- 2011/4/6  added the choice of camera in bgsky.txt

=item 2.4.3 -- 2011/4/4  (internal use) made it possible to not specify any region in bgsky.txt

=item 2.4.2 -- 2010/11/25 sometimes it is not possible to open $$.cdfsextract.stderr in sas_exec during the check phase. Let's try using File::Temp to avoid any possibile race condition..

=item 2.4.1 -- 2010/11/23 (internal use)   corrected bug in printfnlog

=item 2.4 -- 2010/10/21 --photoerror option added: calculates errors on net counts and rates, using BEHR (which must be available)

=item 2.31 -- 2010/09/29 added --noresp option

=item 2.3 -- 2010/08/30 in aperture photometry, psf corrections now
take into account which part of the PSF was exposed

=item 2.2 -- 2010/08/24 calculate, on request, rmf and arf for the bkg files;
also added options --noclobber and --maxproc

=item 2.1 -- 2010/07/21 can make aperture photometry

=item 2.0 -- 2010/06/18 added checks (whether src and bkg fall on
different chips), summary of checks, multiple positions for source or
bkg, exclusion of areas

=item 1.0 -- 2010/04/29  first public version

=back

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

=head1 AUTHOR

Piero Ranalli

Post-doc researcher at IAASARS, National Observatory of Athens, Greece
Associate of INAF -- Osservatorio Astronomico di Bologna, Italy

pranalli.github@gmail.com

=cut
