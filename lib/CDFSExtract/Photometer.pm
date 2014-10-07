package CDFSExtract::Photometer;

use PDL;
use Mo;

=head1 NAME

CDFSExtract::Photometer -- extract aperture photometry

=head1 SYNOPSYS

my $ph = CDFSExtract::Photometer->new;
$ph->extract($image,$expmap,$psf);
$ph->bayesian($child_id,$errorlevel);
$ph->print;

=head1 DESCRIPTION

Photometry code which is called on request during the check phase.
Performs aperture photometry, optionally with Bayesian estimate of count rate
and errors via the BEHR code. Finally prints the results in a format which can
be read by do_checks().

=cut


has src_cts      => (is => 'rw', isa => 'Num');
has bkg_cts      => (is => 'rw', isa => 'Num');
has src_avg_exp  => (is => 'rw', isa => 'Num');
has bkg_avg_exp  => (is => 'rw', isa => 'Num');
has src_area     => (is => 'rw', isa => 'Num');
has bkg_area     => (is => 'rw', isa => 'Num');
has psffrac      => (is => 'rw', isa => 'Num');
has net_cts      => (is => 'rw', isa => 'Num');
has net_lowerr   => (is => 'rw', isa => 'Num');
has net_higherr  => (is => 'rw', isa => 'Num');
has rate         => (is => 'rw', isa => 'Num');
has medianrate   => (is => 'rw', isa => 'Num');
has rate_lowerr  => (is => 'rw', isa => 'Num');
has rate_higherr => (is => 'rw', isa => 'Num');



sub extract {
    my $s = shift;
    my $image = shift;
    my $expmap = shift;
    my $psf = shift;


    #print("$sourceID $evtfile\n");
    my ($src_cts,$bkg_cts) = $image->masked_sums;
    my ($src_avg_exp,$bkg_avg_exp) = $expmap->masked_averages;
    my ($src_area,$bkg_area) = $expmap->areas;

    $s->src_cts($src_cts);
    $s->bkg_cts($bkg_cts);
    $s->src_avg_exp($src_avg_exp);
    $s->bkg_avg_exp($bkg_avg_exp);
    $s->src_area($src_area);
    $s->bkg_area($bkg_area);


    # rate is always classical
    $s->rate(
	     $src_cts/$src_avg_exp -
	     $bkg_cts/$bkg_avg_exp/$bkg_area*$src_area
	    );

    # 'classical' net counts
    $s->net_cts( $src_cts - $bkg_cts/$bkg_area*$src_area );

    #  and no rate errors (unless method bayesian is called)
    $s->net_lowerr( 0 );
    $s->net_higherr( 0 );
    $s->rate_lowerr( 0 );
    $s->rate_higherr( 0 );
    $s->medianrate( 0 );


    $s->psffrac( $psf->eef($expmap) );
    $s->rate( $s->rate / $s->psffrac );
}



sub bayesian {
    my $s = shift;
    my $child_id = shift;
    my $errorlev = shift;

    # Bayesian posterior median of net counts and rate, with errors
    my ($net_cts,$net_lowerr,$net_higherr) =
	$s->get_errors_from_BEHR(
			     $child_id,
			     $errorlev
			    );
    $s->net_cts( $net_cts );
    $s->net_lowerr( $net_lowerr );
    $s->net_higherr( $net_higherr );

    $s->medianrate(   $net_cts     / $s->src_avg_exp / $s->psffrac );
    $s->rate_lowerr(  $net_lowerr  / $s->src_avg_exp / $s->psffrac );
    $s->rate_higherr( $net_higherr / $s->src_avg_exp / $s->psffrac );
}


sub get_errors_from_BEHR {
    require Astro::BEHR;

    my $s = shift;
    my ($id,$errorlev) = @_;

    my $behr = Astro::BEHR->new;

    my $area = $s->bkg_area/$s->src_area*$s->bkg_avg_exp/$s->src_avg_exp;
    $behr->set($s->src_cts,$s->bkg_cts,$area,
	       $s->src_cts,$s->bkg_cts,$area,
	       90);

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


sub print {
    my $self = shift;
    my $fileh = shift;

    print($fileh $self->src_cts."\n");
    print($fileh $self->bkg_cts."\n");
    print($fileh $self->src_avg_exp."\n");
    print($fileh $self->bkg_avg_exp."\n");
    print($fileh $self->src_area."\n");
    print($fileh $self->bkg_area."\n");
    print($fileh $self->psffrac."\n");
    print($fileh $self->net_cts."\n");
    print($fileh $self->net_lowerr."\n");
    print($fileh $self->net_higherr."\n");
    print($fileh $self->rate."\n");
    print($fileh $self->medianrate."\n");
    print($fileh $self->rate_lowerr."\n");
    print($fileh $self->rate_higherr."\n");
}

sub printempty {
    my $self = shift;
    my $fileh = shift;

    print($fileh "0\n")  for 1..14;
}


1;
