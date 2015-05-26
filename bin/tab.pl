#!/usr/bin/perl

use strict;

my @lines;
# first read the file into a list of lists
while (<>)
{
    chomp; # remove the newline from the end of the line
    my @fields = split("\t");
    push @lines, \@fields;
}
my @lengths;
# calculate the maximum lengths of each field
foreach (@lines)
{
    for (my $i = 0; $i < scalar @$_; $i++)
    {
        $lengths[$i] = $lengths[$i] < length $$_[$i] ? length $$_[$i] : $lengths[$i];
    }
}
# now print the text aligned
foreach (@lines)
{
    for (my $i = 0; $i < scalar @$_; $i++)
    {
        print $$_[$i], " " x ($lengths[$i] - length ($$_[$i]) + 1);
    }
    print "\n";
}

