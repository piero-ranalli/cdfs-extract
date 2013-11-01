package EvtFiles;

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



# hold relationships between evtfiles and imgs, expmaps, ccf.cif and sum.sas
#
# (written in traditional OO style to avoid adding Moose as dependency,
#  but syntax is Moose-compatible)

# SYNOPSYS
#
# my $events = EvtFiles->new({evtfile=>'myfilename.txt'});
# for my $e ($events->evts) {
#      my $exp =  $events->expmap{$e};
#      my $img =  $events->img{$e};
#      my $ccf =  $events->ccfcif{$e};
#      my $sum =  $events->sumsas{$e};
# }

use strict;
use warnings;
use Carp;


sub new {
    my $class = shift;
    my $opts = shift;

    croak "event file not specified" unless 
	(defined($opts) and exists($$opts{evtfile}) and defined($$opts{evtfile}));

    croak "event file does not exist" unless (-e $$opts{evtfile});

    my $self = {};
    bless $self,$class;

    $self->{evtfile} = $$opts{evtfile};

    $self->readevt;

    return $self;
}



sub readevt {
    my $self = shift;

    #my (@evts,%maps,%imgs);
    open(my $EVT,'<', $self->{evtfile}) or croak
	'could not read event file '.$self->{evtfile};

    $self->{evts} = [];
    $self->{maps} = {};
    $self->{imgs} = {};
    $self->{ccfcif} = {};
    $self->{sumsas} = {};

    my $n = 0;
    while (my $f=<$EVT>) {

	$n++;

	next if ($f =~ /^\s*\#/);
	chomp($f);
	$f =~ s/\#.*$//;

	my @ff = split(' ',$f);

	if (@ff < 5) {
	    croak "missing elements at row $n in evtlist ".$self->evtfile;
	}

	croak "event file not found at $ff[0]" unless (-e $ff[0]);
	croak "expmap file not found at $ff[1]" unless (-e $ff[1]);
	if (uc($ff[2]) ne 'NA' and not -e $ff[2]) {
	    carp  "image file not found at $ff[2]";
	}
	croak "ccf.cif file not found at $ff[3]" unless (-e $ff[3]);
	croak "sum.sas file not found at $ff[4]" unless (-e $ff[4]);

	push(@{ $self->{evts} },$ff[0]);
	$self->{maps}->{$ff[0]} = $ff[1];
	$self->{imgs}->{$ff[0]} = $ff[2];
	$self->{ccfcif}->{$ff[0]} = $ff[3];
	$self->{sumsas}->{$ff[0]} = $ff[4];

    }
    close($EVT);
}


sub evts {
    my $self = shift;
    return @{ $self->{evts} };
}

sub expmap {
    my $self = shift;
    my $evt = shift;
    return $self->{maps}->{$evt};
}

sub img {
    my $self = shift;
    my $evt = shift;
    return $self->{imgs}->{$evt};
}

sub ccfcif {
    my $self = shift;
    my $evt = shift;
    return $self->{ccfcif}->{$evt};
}

sub sumsas {
    my $self = shift;
    my $evt = shift;
    return $self->{sumsas}->{$evt};
}


1;
