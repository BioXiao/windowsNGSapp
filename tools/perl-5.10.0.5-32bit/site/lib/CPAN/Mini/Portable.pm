package CPAN::Mini::Portable;

use 5.008;
use strict;
use warnings;
use Portable   ();
use CPAN::Mini ();

our $VERSION = '0.12';
our @ISA     = 'CPAN::Mini';

sub new {
	# Use the portable values as defaults,
	# completely ignoring any passed params
	my $minicpan = Portable->default->minicpan;

	# Hand off to the parent class
	return $_[0]->SUPER::new( %$minicpan );
}

1;
