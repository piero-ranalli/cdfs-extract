package CDFSExtract::Source;

use CDFSExtract::Region;
use PDL;
use Mo;
use Astro::WCS2 qw/:img :evt/;

# isa's are ignored in Mo, but I'll specify them anyway for documentation purposes
has id      => ( isa => 'Str' );
has regions => ( isa => 'ArrayRef[Region]' );


sub add {
    my $self = shift;
    my @rest = @_;

    my $reg = CDFSExtract::Region->new( @rest );
    $reg->setparams( @rest );

    $self->regions( [] ) unless defined($self->regions);
    my $reglist = $self->regions;
    push @$reglist,$reg;
    $self->regions( $reglist );
}


sub radec2xy {   # using sky2xy seems not to work in parallel
     # converts from ra,dec to x,y using PDL
     # now accepts multiple positions for same source (i.e. stacked
     # spectra)

    my $self = shift;

    my ($evt,$ext) = @_;	# $src is an array ref

    # read fits header to get CDELTi
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
    my @reglist = @{ $self->regions };
    for my $i (0..$#reglist) {

	if (defined($reglist[$i]->camera)) {
	    # the region only applies to one camera, so let's cycle
	    # unless it is the good one
	    next unless ($reglist[$i]->camera eq $cam);
	}


	my ($x,$y);
	($x,$y) = &$wcs_sub($header,$reglist[$i]->ra,$reglist[$i]->dec);

	push(@new_lol, [ $reglist[$i]->type,
			 $x,
			 $y,
			 $reglist[$i]->radius1xmm,
			 $reglist[$i]->radius2xmm,
		       ]);
    }

    # check if the new lol actually contains any item, and return
    if ($#new_lol >= 0) {
	return(1,@new_lol);
    } else {
	return(0);
    }
}


sub old_style_radii {
    my $self = shift;

    for my $r ( @{ $self->regions } ) {
	$r->radius1xmm( $r->radius1xmm / 20 );
	$r->radius2xmm( $r->radius2xmm / 20 ) if ( $r->radius2xmm );
    }
}
