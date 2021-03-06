package Portable::minicpan;

use 5.008;
use strict;
use warnings;
use Carp ();

our $VERSION = '0.12';





#####################################################################
# Portable Driver API

sub new {
	my $class  = shift;
	my $parent = shift;
	unless ( Portable::_HASH($parent->portable_minicpan) ) {
		Carp::croak('Missing or invalid minicpan key in portable.perl');
	}

	# Create the object
	my $self = bless { }, $class;

	# Map paths to absolute paths
	my $minicpan = $parent->portable_minicpan;
	my $root     = $parent->dist_root;
	foreach my $key ( qw{ local } ) {
		unless (
			defined $minicpan->{$key}
			and
			length $minicpan->{$key}
		) {
			$self->{$key} = $minicpan->{$key};
			next;
		}
		$self->{$key} = File::Spec->catdir(
			$root, split /\//, $minicpan->{$key}
		);
	}

	# Add the literal params
	foreach my $key ( qw{
		remote
		quiet force offline
		also_mirror exact_mirror
		module_filters path_filters
		skip_cleanup skip_perl
	} ) {
		$self->{$key} = $minicpan->{$key};
	}

	return $self;
}

1;
