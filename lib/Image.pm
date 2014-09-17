package Image;

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

use PDL;


sub new {

    # $src and $bkg are Source objects

    my ($class, $src,$bkg,$img) = @_;
    my $self = {};

    $self->{SRC} = $src;
    $self->{BKG} = $bkg;
    $self->{IMG} = rfits($img.'[0]');
    $self->{IMGFILE} = $img;
    $self->{DEADTHRESH} = 1;  # in ks (or whatever the expmap units are)

    my @todo;
    if (defined($self->{BKG})) {
	@todo = qw/SRC BKG/;
    } else {
	@todo = qw/SRC/;
    }

    for my $what (@todo) {
	my ($ok,@xy_rad) = $self->{$what}->radec2xy( $self->{IMGFILE}, '+0' );
	unless ($ok) {
	    if ($what eq 'SRC') {
		main::print_error("Could not read expmap $self->{IMGFILE} or source not specified for its camera.\n");
	    } else {
		main::print_error("Could not read expmap $self->{IMGFILE} -- Or background not specified for its camera.\n");
	    }
	    return(0);
	}
	$self->{"${what}_XYRAD"} = \@xy_rad;
    }

    bless( $self, $class );
    return($self);
}

sub dead_fractions {
    my $self = shift;

    unless ( exists($self->{SRC_MASK}) ) {
	($self->{SRC_MASK},$self->{BKG_MASK}) = $self->calc_masks;
    }

    my @out;
    my @what;
    if (defined($self->{BKG_MASK})) {
	@what = (qw/SRC BKG/);
    } else {
	@what = (qw/SRC/);
    }

    for my $what (@what) {

	my $cut = $self->{IMG}->where($self->{"${what}_MASK"});
	my $pixtot = $cut->flat->dim(0);

	if ($pixtot==0) {  # src is completely outside image
	    push(@out,1);  # return 100% dead fraction
	    next;
	}

	my $pm1 = $cut->where($cut<($self->{DEADTHRESH}*1000))->dim(0);

	#  printf("pixtot=%5i pix<1=%5i pix<max/2=%5i mean=%10.3f sigma=%10.3f median=%10.3f\n",
	#	$pixtot,$pm1,$pm2,$mean,$sigma,$med);

	push(@out,$pm1/$pixtot);

    }
    return(@out);
}


sub calc_masks {
    my ($self) = @_;

    my @masks;
    my $firstpass = 1;
    my @list = ();
    if ($#{$self->{BKG_XYRAD}} >= 0) {
	@list = ($self->{SRC_XYRAD},$self->{BKG_XYRAD});
    } else {
	@list = ($self->{SRC_XYRAD});
    }

    for my $xy_rad_ref (@list) {

	my @xy_rad = @$xy_rad_ref;

	my $fullmask = zeroes $self->{IMG};

	# there may be more than 1 area for src/bkg
	for my $i (0..$#xy_rad) {

	    my $r1 = main::phys2img(int($xy_rad[$i][3]), $self->{IMG}->hdr->{CDELT1}) ;
	    if ($firstpass) {
		# ony set $self->{FIRST_SRC_AREA} on the first area
		# of the source
		$self->first_src_area( $xy_rad[$i][1],$xy_rad[$i][2],
				       $r1 );
		$firstpass = 0;
	    }
	    my $msk = rvals( $self->{IMG}, {CENTRE=>[$xy_rad[$i][1],$xy_rad[$i][2]],
					    SQUARED=>1} )
		<= $r1*$r1;

	    if ($xy_rad[$i][0] eq 'annulus') {
		my $r2 = main::phys2img(int($xy_rad[$i][4]), $self->{IMG}->hdr->{CDELT1}) ;
		my $msk2 = rvals( $self->{IMG}, {CENTRE=>[$xy_rad[$i][1],$xy_rad[$i][2]],
						 SQUARED=>1} ) <= $r2*$r2;
		$msk = $msk2 - $msk;
	    }

	    $fullmask += $msk;
	}
	push(@masks,$fullmask);
    }

    return(@masks);
}

sub first_src_area {
    my ($self,$x,$y,$r) = @_;

    if (defined($x) and defined($y) and defined($r)) {
	$self->{FIRST_SRC_AREA} = [$x,$y,$r];
    }

    return( @{$self->{FIRST_SRC_AREA}} );
}

sub masks {
    my ($self,$srcmask,$bkgmask) = @_;

    if (defined($srcmask)) {
	$self->{SRC_MASK} = $srcmask;
	$self->{BKG_MASK} = $bkgmask;
    }

    return($self->{SRC_MASK},$self->{BKG_MASK});
}

sub exposed_masks {
    # returns the mask (i.e. the extraction area) multiplied by the mask of
    # the 'on' pixels in the expmap)
    my $self = shift;
    my $on_msk = $self->{IMG} > 1;
    return( $self->{SRC_MASK} * $on_msk, $self->{BKG_MASK} * $on_msk );
}

sub masked_sums {
    my $self = shift;

    return( $self->{IMG}->where($self->{SRC_MASK})->sum,
	    $self->{IMG}->where($self->{BKG_MASK})->sum );
}
sub masked_averages {
    my $self = shift;

    return( $self->{IMG}->where($self->{SRC_MASK})->avg,
	    $self->{IMG}->where($self->{BKG_MASK})->avg );
}
sub areas {			# returns in arcsec
    my $self = shift;

    return( $self->{SRC_MASK}->sum * ($self->{IMG}->hdr->{CDELT1} * 3600)**2,
	    $self->{BKG_MASK}->sum * ($self->{IMG}->hdr->{CDELT2} * 3600)**2 );
}

sub src_radius {
    my $self = shift;
    #      my $r = shift;

    #       if (defined($r)) {
    # 	  $self->{RADIUS} = $r;
    #       }
    #       return($self->{RADIUS});

    # old line, used when sources where stored in hash of arrays of arrays:
    #return(${$self->{SRC}}[0][3] /20); # in arcsec!

    # new version using Source and Region objects:
    my $reg = $self->{SRC}->regions->[0]; # take first region of source
    return $reg->radius1;  # the /20 is no longer needed since we switched
                           # to arcsecs in the radii
}

sub calc_offaxis {
    my $self = shift;

    # take the position of the first source
    my $reg = $self->{SRC}->regions->[0];
    my $ra = $reg->ra;
    my $dec= $reg->dec;

    my $oa = sqrt( ($ra - $self->{IMG}->hdr->{RA_PNT}) ** 2
		   + ($dec - $self->{IMG}->hdr->{DEC_PNT}) ** 2 );

    # now get the angle to which the PSF will be rotated.
    # this angle starts counting on the right X axis, and is positive
    # counterclockwise (same as in the simulator)
    my $rotation = atan2(-$ra+$self->{IMG}->hdr->{RA_PNT},
			 $dec-$self->{IMG}->hdr->{DEC_PNT})*180/3.141592 +90;
    # convert to PDL::Transform convention (positive angle rotates
    # clockwise
    $rotation = -$rotation;

    return($oa,$rotation);
}

#   sub psf_fraction {  # (of the first region of the source)
#       use XMMPSF;

#       my $self = shift;

#       my $ra = $self->{SRC}[1];
#       my $dec= $self->{SRC}[2];
#       my $pointra = $self->{IMG}->hdr->{RA_PNT};
#       my $pointdec= $self->{IMG}->hdr->{DEC_PNT};

#       my $psf = XMMPSF->new($ra,$dec,$pointra,$pointdec);

#       return($psf->fraction($self->{SRC_XYRAD}[3])); # radius in physical units
#   }



1;
