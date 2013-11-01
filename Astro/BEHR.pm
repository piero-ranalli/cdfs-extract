package Astro::BEHR;

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


use strict;
use warnings;
#use PDL;

=head1 NAME

Astro::BEHR -- Bayesian Estimate of Hardness Ratios

A Perl interface to the program by Park et al. (2006) ApJ, 652, 610.

=head1 SYNOPSIS

 # check if the module can be used:
 # (the executable "BEHR" must be in your $PATH)
 if (BEHR::check()) {
     require BEHR;
 }

=cut

sub check {
    my $class = shift;

    # first, check that BEHR is actually present
    my $stat = system('BEHR > /dev/null');
    if ($stat % 256 == 0) { # the shell did return ok
	return(1);
    } else {
	print('BEHR could not be found in the path.\n\n');
	return(0);
    }
}


=pod

 my $behr = BEHR->new( [{CHECK=>1}] );  # build the object, with optional call to BEHR::check

=cut

sub new {
    my $class = shift;

    # if requested, check that BEHR is actually present
    my $opt = shift;
    if (defined($opt) and $$opt{CHECK}) {
	return 0 unless BEHR::check();
    }

    # now build the BEHR class
    my $self = { ALGO=>'quadrature' };
    bless($self,$class);
    return $self;
}


=pod

 $behr->{ALGO} = 'gibbs';  # 'quadrature' is the default

 $behr->mosix(1);   # turn on mosix
 $behr->mosix(0);   # turn off mosix

=cut

sub mosix {
    my $self = shift;
    my $flag = shift;
    $self->{MOSIX} = $flag;
}



=pod

 $behr->timeout($n);  # timeout after n seconds (timeout only works if mosix=0)
 $n = $behr->timeout; # returns the current timeout value (default: undef)

=cut

sub timeout {
    my $self = shift;
    my $t = shift;

    if (defined($t)) {
	$self->{TIMEOUT} = $t;

	if ($t>0) {
	    $self->{MOSIX} = 0;
	}

    } else {
	return($self->{TIMEOUT});
    }
}


=pod

 # set soft source counts, soft bkg counts, soft area,
hard source counts, hard bkg counts, hard area,
credible level:
 $behr->set($ssrc,$sbkg,$sarea,$hsrc,$hbkg,$harea,$lev);

=cut

sub set {
    my $self = shift;
    my ($ssrc,$sbkg,$sarea,$hsrc,$hbkg,$harea,$lev) = @_;

    $self->{SSRC} = $ssrc;
    $self->{SBKG} = $sbkg;
    $self->{SAREA} = $sarea;
    $self->{HSRC} = defined($hsrc) ? $hsrc : $ssrc;
    $self->{HBKG} = defined($hbkg) ? $hbkg : $sbkg;
    $self->{HAREA} = defined($harea) ? $harea : $sarea;
    $self->{LEV} = defined($lev) ? $lev : 90;
}


=pod

 $err = $self->run;   # does the calculation, returns any error

=cut

sub run {
    my $self = shift;

    my $s = "softsrc=".$self->{SSRC}.
	" softbkg=".$self->{SBKG}.
	" softarea=".$self->{SAREA}.
	" hardsrc=".$self->{HSRC}.
	" hardbkg=".$self->{HBKG}.
	" hardarea=".$self->{HAREA}.
	" level=".$self->{LEV}.
	" details=true algo=".$self->{ALGO};

    $self->{OUTPUT_HEADER} = "";

    my $BEHR;
    my $pid;
    if ($self->{MOSIX}) {
	open($BEHR, "mosrun BEHR $s |");
    } else {
	$pid = open($BEHR, "BEHR $s |") or die "Could not run BEHR\n";
    }

    eval {
	my $r;

	local $SIG{ALRM} = sub { die 'TIMEOUT' };
	if (defined($self->{TIMEOUT})) {
	    alarm($self->{TIMEOUT});
	}
	while($r = <$BEHR>) {

	    if ($r =~ m/^\#/) {
		$self->{OUTPUT_HEADER} .= $r;
		next;
	    } elsif ($r =~ m/^\s+/) {
		next;  # spazi vuoti e intestazioni
	    } else {  #il resto invece ci interessa

		chomp($r);
		my ($what,$mode,$mean,$median,$low,$high) = split(/\s+/,$r);
		$self->{$what} = [$mode,$mean,$median,$low,$high];
	    }
	}
	close($BEHR);
	if (defined($self->{TIMEOUT})) {
	    alarm(0);
	}
    };

    if ($@) {      # errors
	die $@ unless ($@ =~ m/TIMEOUT/);  # die in case of unexpected error
	                                          # also, match here is needed insteal of eq checking
	# otherwise, in case of timeout:
	kill 2, $pid;
	close($BEHR);
	return('TIMEOUT');
    } else {
	return(0); # no errors
    }


}


=pod

 $vals = $behr->get('(H-S)/(H+S)');  # put mode,mean,median,low,high in an arrayref

=cut

sub get {
    my $self=shift;
    my $what=shift;

    return($self->{$what});
}


=head1 DESCRIPTION

First, check that you have BEHR installed and reachable from your
$PATH. If you don't have it installed, you can get it from:
 http://hea-www.harvard.edu/astrostat/behr/

BEHR calculates posterior probabilities for hardness ratios, net
counts, and a few variations on these quantities. For more details,
see Park et al. 2006 (ApJ, 652, 610):
 http://cdsads.u-strasbg.fr/abs/2006ApJ...652..610P

=head1 TIMEOUTS

Sometimes BEHR using the quadrature algorithm goes on infinite
loop. You can use timeouts to kill the calculation if takes too long,
and switch to the gibbs algorithm. Here is how to do:

 $behr->{ALGO} = 'quadrature';  # precise but slow and prone to infinite loop

 $behr->timeout(60);  # 1 minute
 my $err = $behr->run;

 if ($err eq 'TIMEOUT') {

    # try again:
    $behr->{ALGO} = 'gibbs';  # much faster
    my $err2 = $behr->run;

    if ($err2 eq 'TIMEOUT') {  # AGAIN?!?
         die "Quitting because of double timeout.\n";
    }
 }

 my $vals = $behr->get($whatever);

 # and if you don't want timeouts any longer:
 $behr->timeout(0);

=cut


1;

