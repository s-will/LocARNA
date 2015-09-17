#!/usr/bin/env perl

## update index.html to include list release; move currently latest
## release to list of old releases

$RELEASE=shift @ARGV;
if ( ! defined($RELEASE) ) {
    print STDERR "USAGE: update_index.pl <release>\n";
    print STDERR int(@ARGV)."\n";
    exit -1;
}

my $latest_release="";
my $older_releases="";

while (<>) {
    if (/^<!--!LATEST RELEASE/) {
	while (<>) {
	    if (/^<!--!END LATEST RELEASE/) {
		print "<!--!LATEST RELEASE -- don't remove this line -- don't edit until END-->
<b>Latest Release: <a href=\"Releases/locarna-$RELEASE.tar.gz\"
	onClick=\"_gaq.push(['_trackEvent', 'Download', 'LocARNA', '$RELEASE']);\"
	>LocARNA $RELEASE</a></b>
<!--!END LATEST RELEASE --- don't remove this line-->
";
		last;
	    } else {
		$latest_release.=$_;
	    }
	}

	$latest_release=~s/^\s*<b>Latest Release: //;
	$latest_release=~s/<\b>\s*$//;
	
    } elsif (/^<!--!OLDER RELEASES/) {
	while (<>) {
	    if (/^<!--!END OLDER RELEASES/) {
		chomp $older_releases;
		print "<!--!OLDER RELEASES --- don't remove this line-->
<li>
  $latest_release
</li>
$older_releases
<!--!END OLDER RELEASES --- don't remove this line-->
";
		last;
		$latest_release="";
		$older_releases="";
	    } else {
		$older_releases.=$_;
	    }
	}
    } else {
	print $_;
    }
}
