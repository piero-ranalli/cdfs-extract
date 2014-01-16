package CDFSExtract::Region;

use Scalar::Util qw/looks_like_number/;
use Mo;
use Carp;

# isa's are ignored in Mo, but I'll specify them anyway for documentation purposes
has ra      => ( isa => 'Num' );
has dec     => ( isa => 'Num' );
has type    => ( isa => 'Str' );
has radius1 => ( isa => 'Num' );
has radius2 => ( isa => 'Num' );
has radius1xmm => ( isa => 'Num' );  # in XMM physical pixels
has radius2xmm => ( isa => 'Num' );
has camera  => ( isa => 'Str' );





sub setparams {
    my $self = shift;
    my @rest = @_;

    # look if region should be excluded or included
    my $exclude = 0;
    if ($rest[0] =~ m/exclude/) {
	$exclude = 1;
	shift @rest;
    }

    my ($ra,$dec,$radius1,$radius2,$camera);
    ($ra,$dec,$radius1,$radius2,$camera) = @rest;

    # sources to be excluded have negative radii
    if ($exclude) {
	$radius1 = -$radius1;
	if (defined($radius2) and looks_like_number($radius2)) {
	    $radius2 = -$radius2;
	}
    }

    unless (defined($radius1)) {
	die "radius not defined\n"
    }

    $self->ra( $ra );
    $self->dec( $dec );
    $self->radius1( $radius1 );
    $self->radius1xmm( $radius1 * 20 );


    if (defined($camera)) { # sure it's an annulus, whatever $camera contains
	$self->type('annulus');

	die "cannot make sense of radius2: $radius2\n" unless looks_like_number($radius2);

	$self->type('annulus');
	$self->radius2( $radius2 );
	$self->radius2xmm( $radius2 * 20 );
	$self->camera( uc($camera) );
	$self->checkcamera;


    # else, it may or may not have a camera specified. We have to check.
    } elsif (defined($radius2) and looks_like_number($radius2)) {

	# then it's an annulus without camera
	$self->type('annulus');
	$self->radius2( $radius2 );
	$self->radius2xmm( $radius2 * 20 );


    # otherwise, it's a circle; with or without a defined camera
    } else {
	$self->type('circle');

	if (defined($radius2)) { # it's the camera
	    $self->camera( uc($radius2) );
	    $self->checkcamera;
	}
    }

}






sub checkcamera {
    my $self = shift;

    # Let's warn if we don't recognise it.
    my $c = $self->camera;
    my $i = 0;
    for my $known (qw/EMOS1 EMOS2 EPN ACIS-I ACIS-S/) {
      $i++ if ($c eq $known);
    }
    carp "Carmera $c not recognised" unless $i;
}


sub arcsec2physpix {
    my $arcsec = shift;

    # the following is the XMM convention
    $$arcsec*=20;
}
