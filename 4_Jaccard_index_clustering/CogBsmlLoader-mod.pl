#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1  NAME 

CogBsmlLoader.pl  -  Preprocess data stored in BSML pairwise alignment documents into BTAB
structure for COG analysis using best_hits.pl. 

=head1 SYNOPSIS

USAGE:  CogBsmlLoader.pl -b BsmlPairwiseAlignmentDirectory -o OutputFile

=head1 OPTIONS

=over 4

=item *

B<--help,-h> This help message

=back

=head1   DESCRIPTION

CogBsmlLoader.pl is designed to preprocess the data contained in a BSML pairwise alignment search 
for COGS analysis. Specifically it identifies the "best hit" per genome for each query gene. 
This data is packaged into the BTAB format for linkage analysis using best_hits.pl  

NOTE:  

Calling the script name with NO flags/options or --help will display the syntax requirement.

=cut

# Preprocess data stored in BSML pairwise alignment documents into BTAB
# structure for COG analysis.

############
# Arguments:
#
# outfile - btab output file
#
#

use strict;
use FileHandle;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use Ergatis::Logger;
use XML::Parser;

my $PROGRESS_INTERVAL = 100;
# dummy value for comp/ref genome
my $GENOME = 'META';

my %options = ();
my $results = GetOptions( \%options, 
              'bsmlSearchList|b=s', 
              'jaccardClusters|j=s', 
              'outfile|o=s', 
              'pvalcut|p=s', 
              'coverageCutoff|c=s',
              'identityCutoff|i=s',
              'similarityCutoff|s=s',
              'log|l=s',
              'debug=s',
              'help|h') || pod2usage();

my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
				  'LOG_LEVEL'=>$options{'debug'});
$logger = Ergatis::Logger::get_logger();

# display documentation
if( $options{'help'} ){
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDERR} );
}


&check_parameters(\%options);

#MAIN HERE

# If Jaccard data has been specified, use it for equivalence class filtering

# jaccardClusterHash maps sequence id to Jaccard Cluster id (jaccardClusterCount++)
my $jaccardClusterHash = {};
# jaccardRepSeqHash maps jaccardClusterCount to representative sequence (the first one)
my $jaccardRepSeqHash = {};
# used to generate Jaccard cluster ids, starting with 1
my $jaccardClusterCount = 0;

# parse Jaccard cluster .out file:
if( $options{'jaccardClusters'} && $options{'jaccardClusters'} ne "" )
  {
    if(-e $options{'jaccardClusters'}) {
      my $seqCount=0;
      my $jlnum = 0;
      my $cname = undef;
      my $csize = undef;
      my $cseqs = undef;
      
      my $jfh = FileHandle->new();
      $jfh->open($options{'jaccardClusters'}) || die "unable to read from $options{'jaccardClusters'}";
      
      # check size of previous cluster, if any
      my $check_size = sub {
        my($lnum) = @_;
        return unless (defined($cname));
        my $csize_observed = scalar(@$cseqs);
        if (defined($csize) && ($csize != $csize_observed)) {
          die "cluster size mismatch (reported=$csize, observed=$csize_observed) at line $lnum";
        }
      };
      
      while (my $line = <$jfh>) {
        chomp($line);
        ++$jlnum;
        # new cluster
        if ($line =~ /^COG = jaccard\_(\d+)\, size \d+\, connections  = \-1, perfect = \-1;$/) {
          my($cn, $cs) = ($1, $2);
          ++$jaccardClusterCount;
          print STDERR "$jaccardClusterCount Jaccard clusters read\n" if (($jaccardClusterCount % $PROGRESS_INTERVAL) == 0);
          &$check_size($jlnum);
          ($cname, $csize, $cseqs) = ($cn, $cs, []);
        } 
        # not a new cluster
        elsif ($line =~ /^\s+(\S+)$/) {
          my $prot = $1;
          $jaccardClusterHash->{$prot} = $jaccardClusterCount;
          # first seq will be the representative
          if (scalar(@$cseqs) == 0) {
            $jaccardRepSeqHash->{$jaccardClusterCount} = $cseqs;
          }
          push(@$cseqs, $prot);
          ++$seqCount;
        } 
        else {
          die "couldn't parse line $jlnum of $options{'jaccardClusters'}";
        }
      }
      print STDERR "done reading clusters\n";
      $jfh->close();
      &$check_size($jlnum);
      print STDERR "read $seqCount sequence(s) in $jaccardClusterCount cluster(s) from $options{'jaccardClusters'}\n";
    }
  }

#####################################

# structure for building the COGS input. For each query gene, the COGS analysis expects
# the single best scoring alignment for each reference genome. In BSML terms, COGS expects the
# highest scoring seqpair for each refseq compseq combination where all compseqs 
# are contained in the same genome. 


#  Genome A           Genome B            Genome C
#  compseqA1          compseqB1           compseqC1
#  compseqA2          compseqB2           compseqC2
#  compseqA3                              compseqC3


# If the above represent the sets of reference genes by genome. The following would 
# correspond to the expected output if A1, B2, and C1 were the best hits by genome. 

# refseq -> compseqA1
# refseq -> compseqB2
# refseq -> compseqC1

####################################
    
my $COGInput = {};

my %alnparams;
my $bestRunScore = 0;
my $isbestrun = 0;
my $bestSeqPairRun = undef;

my $alnhandlers = {'Seq-pair-alignment'=>
             sub {
             my ($expat,$elt,%params) = @_;
             #Process the previous alignment in file
             &process_alignment($alnparams{'compseq'},$GENOME,
                        $alnparams{'refseq'},$GENOME,
                        $bestRunScore,$bestSeqPairRun,
                        $COGInput,$jaccardRepSeqHash) if(keys %alnparams);

             my $compseq = $params{'compseq'};
             my $refseq = $params{'refseq'};
             $bestRunScore = 0;
             $isbestrun=0;
             $bestSeqPairRun = undef;
             %alnparams = ();

             # self-self alignments are not included 
             return if( $compseq eq $refseq );
             %alnparams = %params;
             },
           'Seq-pair-run'=>
               sub {
               my ($expat,$elt,%params) = @_;
               if(keys %alnparams){
                   my $runscore = $params{'runscore'};
                   my $runprob = $params{'runprob'};
                   if( defined $runscore && defined $runprob
                   && ($runscore > $bestRunScore) && ($runprob < $options{'pvalcut'}) ){
                   $logger->debug("$alnparams{'compseq'} $alnparams{'refseq'} using run with runscore $runscore $runprob. Previous bestrunscore $bestRunScore. pvalue cutoff $options{'pvalcut'}");
                   $bestRunScore = $runscore;
                   $logger->debug("bestrunscore $bestRunScore");
                   $bestSeqPairRun->{'reflength'} = $alnparams{'reflength'};
                   $bestSeqPairRun->{'method'} = $alnparams{'method'};
                   $bestSeqPairRun->{'compxref'} = $alnparams{'compxref'};
                   $bestSeqPairRun->{'refpos'} = $params{'refpos'};
                   $bestSeqPairRun->{'runlength'} = $params{'runlength'};
                   $bestSeqPairRun->{'comppos'} = $params{'comppos'};
                   $bestSeqPairRun->{'comprunlength'} = $params{'comprunlength'};
                   $bestSeqPairRun->{'runscore'} = $params{'runscore'};
                   $bestSeqPairRun->{'runprob'} = $runprob;
                   $isbestrun=1;
                   }
               }
               },
           'Attribute'=>
               sub {
               my ($expat,$elt,%params) = @_;
               my $index = scalar(@{$expat->{'Context'}}) - 1;
               if($isbestrun && $expat->{'Context'}->[$index] eq 'Seq-pair-run'){
                   $logger->debug("Dumping parameters for best run $bestRunScore");
                   if($params{'name'} eq 'p_value'){
                   $bestSeqPairRun->{'p_value'} = $params{'content'};
                   }
                   elsif($params{'name'} eq 'percent_identity'){
                   $bestSeqPairRun->{'percent_identity'} = $params{'content'};
                   }
                   elsif($params{'name'} eq 'percent_similarity'){
                   $bestSeqPairRun->{'percent_similarity'} = $params{'content'};
                   }
                   elsif($params{'name'} eq 'chain_number'){
                   $bestSeqPairRun->{'chain_number'} = $params{'content'};
                   }
                   elsif($params{'name'} eq 'segment_number'){
                   $bestSeqPairRun->{'segment_number'} = $params{'content'};
                   }
               }
               if($expat->{'Context'}->[$index] eq 'Seq-pair-alignment'){
                   if($params{'name'} eq 'percent_coverage_refseq'){
                       if($params{'content'} < $options{'coverageCutoff'}) {
                           %alnparams = ();
                       }
                   }
                   if($params{'name'} eq 'percent_identity'){
                       if($params{'content'} < $options{'identityCutoff'}) {
                           %alnparams = ();
                       }
                   }
                   if($params{'name'} eq 'percent_similarity'){
                       if($params{'content'} < $options{'similarityCutoff'}) {
                           %alnparams = ();
                       }
                   }
               }
           }
           };

my $alnparser = new XML::Parser(Handlers => 
                   {
                       Start =>
                       sub {
                        #$_[1] is the name of the element
                           if(exists $alnhandlers->{$_[1]}){
                           
                           $alnhandlers->{$_[1]}(@_);
                        }
                    }
                                }
                );


open( OUTFILE, ">$options{'outfile'}" ) or $logger->logdie("Can't open file $options{'outfile'}");

my $n_docs_processed = 0;
my $bsml_file_list = &get_list_from_file($options{'bsmlSearchList'});
my $n_docs = scalar(@$bsml_file_list);

foreach my $bsmlFile (@$bsml_file_list){

    # print timestamped progress message
    if (($n_docs_processed % $PROGRESS_INTERVAL) == 0) {
      my $date = `date`;
      chomp($date);
      my $progress_msg = sprintf("$date processed %d/%d BSML documents", $n_docs_processed, $n_docs);
      $logger->info($progress_msg);
      print STDERR "$progress_msg\n";
    }
    ++$n_docs_processed;
    
    if (!(-e $bsmlFile) && -e "$bsmlFile.gz") {
        $bsmlFile .= ".gz";
    }
    
    # builds the COGS input data structure

    $logger->debug("Parsing alignment file $bsmlFile") if($logger->is_debug());
    
    %alnparams = ();
    $bestRunScore = 0;
    $isbestrun = 0;
    $bestSeqPairRun = undef;

    my $ifh;
    if ($bsmlFile =~ /\.(gz|gzip)$/) {
        open ($ifh, "<:gzip", $bsmlFile) || die "can't read input file $bsmlFile: $!";
    } else {
        open ($ifh, "<$bsmlFile") || die "can't read input file $bsmlFile: $!";
    }
    $alnparser->parse( $ifh );
    close $ifh;
    
    #Process the last alignment in file
    &process_alignment($alnparams{'compseq'},$GENOME,
               $alnparams{'refseq'},$GENOME,
               $bestRunScore,$bestSeqPairRun,
               $COGInput,$jaccardRepSeqHash) if(keys %alnparams);;
    
    # print the results

    foreach my $k1 ( keys( %{$COGInput} ) )
    {
    foreach my $k2 (keys( %{$COGInput->{$k1}}))
    {
        my $member = $COGInput->{$k1}->{$k2}->[0];
        if(exists $jaccardClusterHash->{$member}){
        $COGInput->{$k1}->{$k2}->[21] = join(',',@{$jaccardRepSeqHash->{$jaccardClusterHash->{$member}}});
        }
        print OUTFILE join("\t", @{$COGInput->{$k1}->{$k2}});
        print OUTFILE "\n";
    }
    }

    $COGInput = {};
}


sub process_alignment{
    my($compseq,$compGenome,$refseq,$refGenome,$bestRunScore,$bestSeqPairRun,$COGInput,$jaccardRepSeqHash) = @_;
    # 
    if( ! defined $bestSeqPairRun ){
    $logger->warn("Best run not defined for $refseq");
    return;
    }
    else{
# If compseq (or refseq) is defined in a Jaccard equivalence class identify the class by
# its reference sequence. 
    
    if( defined( my $jId = $jaccardClusterHash->{$compseq} ) )
    {
        $logger->debug("Found jaccard cluster $jId for id $compseq. Using $jaccardRepSeqHash->{$jId}->[0] as cluster representative");
        $compseq = $jaccardRepSeqHash->{$jId}->[0];
    }
    
    if( defined( my $jId = $jaccardClusterHash->{$refseq} ) )
    {
        $logger->debug("Found jaccard cluster $jId for id $refseq. Using $jaccardRepSeqHash->{$jId}->[0] as cluster representative");
        $refseq = $jaccardRepSeqHash->{$jId}->[0];
    }
    
    my $lref = [];
    
    $lref->[0] = $refseq;  #query name
    $lref->[1] = '';       #date
    $lref->[2] = $bestSeqPairRun->{ 'reflength' }; #query length
    $lref->[3] = $bestSeqPairRun->{ 'method' }; #program
    $lref->[4] = $bestSeqPairRun->{ 'compxref' };
    $lref->[5] = $compseq;
    $lref->[6] = $bestSeqPairRun->{ 'refpos' };
    $lref->[7] = $bestSeqPairRun->{ 'refpos' } + $bestSeqPairRun->{'runlength'};
    $lref->[8] = $bestSeqPairRun->{ 'comppos' };
    $lref->[9] = $bestSeqPairRun->{ 'comppos' } + $bestSeqPairRun->{'comprunlength'};
    $lref->[10] = $bestSeqPairRun->{'percent_identity'};
    $lref->[11] = $bestSeqPairRun->{'percent_similarity'};
    $lref->[12] = $bestSeqPairRun->{ 'runscore' };
    $lref->[13] = $bestSeqPairRun->{'chain_number'};
    $lref->[14] = $bestSeqPairRun->{'segment_number'};
    $lref->[15] = '';
    $lref->[16] = '';
    $lref->[17] = '';
    $lref->[18] = $bestSeqPairRun->{'comprunlength'};
    $lref->[19] = $bestSeqPairRun->{'runprob' };
    $lref->[20] = $bestSeqPairRun->{'p_value'};
    
    if($refseq && $compseq){
        if(  $COGInput->{$refseq}->{$compGenome} )
        {
        if(  $COGInput->{$refseq}->{$compGenome}->[12] < $bestRunScore )
        {
            $logger->debug("$refseq match to $compGenome with score $bestRunScore is highest scoring match.  Previous high score is $COGInput->{$refseq}->{$compGenome}->[12]");
            $COGInput->{$refseq}->{$compGenome} = $lref;
        }
        }
        else
        {
        $logger->debug("$refseq match to $compGenome is first match found.");
        $COGInput->{$refseq}->{$compGenome} = $lref;
        }
    }
    }
}


sub get_list_from_file{
    my($file) = @_;
    my @lines;
    open( FH, $file ) or $logger->logdie("Could not open $file");
    while( my $line = <FH> ){
    chomp($line);
    push @lines, split(',',$line) if($line =~ /\S+/);
    }
    return \@lines;
}

sub check_parameters{
    my ($options) = @_;
    
    if(0){
    pod2usage({-exitval => 2,  -message => "error message", -verbose => 1, -output => \*STDERR});    
    }
}

# Returns only those BLAST HSPs that contributed to the best Sum(P) value
# reported for each subject/query sequence pair.
#
sub filterBlastpHsps {
    my($links) = @_;
    my $linksByQuery = &groupByMulti($links, ['from_tag', 'to_tag']);
    my $result = [];

    # Aggregate all HSPs with the same subject and query sequence
    foreach my $queryId (keys %$linksByQuery) {
	my $linksBySubject = $linksByQuery->{$queryId};

	foreach my $subjId (keys %$linksBySubject) {
	    my $slinks = $linksBySubject->{$subjId};
	    my @sortedLinks = sort { $a->{'p_value'} <=> $b->{'p_value'} } @$slinks;
	    # heuristic - assume that all HSPs with the same Sum(P) as the best are contributing to that Sum(p) score
	    my $bestScore = $sortedLinks[0]->{'p_value'};
	    
	    foreach my $sl (@sortedLinks) {
		last if ($sl->{'p_value'} > $bestScore);
		push(@$result, $sl);
	    }
	}
    }
    return $result;
}

# Generalized version of groupBy 
sub groupByMulti {
    my($arrayref, $keyFields) = @_;
    my $nKeys = scalar(@$keyFields);
    my $groups = {};

    foreach my $a (@$arrayref) {
	my @keyValues = map { $a->{$_} } @$keyFields;
	my $hash = $groups;

	for (my $i = 0;$i < $nKeys;++$i) {
	    my $kv = $keyValues[$i];

	    if ($i < ($nKeys-1)) {
		$hash->{$kv} = {} if (!defined($hash->{$kv}));
		$hash = $hash->{$kv};
	    } 
	    else {
		$hash->{$kv} = [] if (!defined($hash->{$kv}));
		push(@{$hash->{$kv}}, $a);
	    }
	}
    }
    return $groups;
}

sub getAvgBlastPPctCoverage {
    my($hsps) = @_;
    $hsps = &filterBlastpHsps($hsps);
    my $sum = 0;
    my $numHsps = 0;

    # Group by query and target id
    my $hspsByQuery = &groupByMulti($hsps, ['from_tag', 'to_tag']);

    foreach my $queryId (keys %$hspsByQuery) {
	my $hspsByTarget = $hspsByQuery->{$queryId};

	foreach my $subjId (keys %$hspsByTarget) {
	    ++$numHsps;
	    my $shsps = $hspsByTarget->{$subjId};
	    my $querySeqLen = $shsps->[0]->{'from_length'};
	    my $targetSeqLen = $shsps->[0]->{'to_length'};

	    my @queryIntervals = map { {'fmin' => $_->{'from_Nterm'}, 'fmax' => $_->{'from_Cterm'}, 'strand' => 1} } @$shsps;
	    my @targetIntervals = map { {'fmin' => $_->{'to_Nterm'}, 'fmax' => $_->{'to_Cterm'}, 'strand' => 1} } @$shsps;

	    my $mergedQueryIntervals = &mergeOverlappingIntervals(\@queryIntervals);
	    my $mergedTargetIntervals = &mergeOverlappingIntervals(\@targetIntervals);

	    my $queryHitLen = 0;
	    my $targetHitLen = 0;
	    
	    map { $queryHitLen += ($_->{'fmax'} - $_->{'fmin'}); } @$mergedQueryIntervals;
	    map { $targetHitLen += ($_->{'fmax'} - $_->{'fmin'}); } @$mergedTargetIntervals;

	    $sum += $queryHitLen / $querySeqLen;
	    $sum += $targetHitLen / $targetSeqLen;
	}
    }

    return ($numHsps > 0) ? ($sum/($numHsps * 2) * 100.0) : undef;
}

# Generate a new set of intervals by merging any that overlap in the original set.
#
sub mergeOverlappingIntervals {
    my($intervals) = @_;

    # result set of intervals
    my $merged = [];

    # sort all intervals by fmin
    my @sorted = sort { $a->{'fmin'} <=> $b->{'fmin'} } @$intervals;
    
    # current interval
    my $current = undef;

    foreach my $i (@sorted) {
	# case 1: no current interval
	if (!defined($current)) {
	    $current = $i;
	} 
	# case 2: compare current interval to interval $i
	else {
	    # case 2a: no overlap
	    if ($i->{'fmin'} > $current->{'fmax'}) {   
		push(@$merged, $current);
		$current = $i;
	    } 
	    # case 2b: overlap, with $i ending to the right of $current
	    elsif ($i->{'fmax'} > $current->{'fmax'}) {
		$current->{'fmax'} = $i->{'fmax'};
	    }
	}
    }
    push(@$merged, $current) if (defined($current));
    return $merged;
}
