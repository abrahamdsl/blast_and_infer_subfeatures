#!/usr/bin/perl
=doc
@task Move from band-aid OS-centric solutions ( ex. "wc -l") to Perl-ly solutions.
@task Move from print statements to utilization of handle_message() 
@task Dynamically named final output file.
@task Integrate/use library/simplify with C4HPC:/storage2/Allave/032_revisit-026-after-fishy-gene-success-ratio-discovery/./taskblastn_experiment/bulk_compare_gene_sequences.pl
  in mind.
@task No way to check yet if $fascistGeneEnd is > chromosome length.

a.llave@irri.org 12SEP2014-1040

A (quite messy but usable) script to infer gene and its subfeatures on a (longer, not sure with shorter) suspected corresponding gene sequence
in a new, but almost identical genome.


Requires:
  Benchmark
  Bio::DB::SeqFeature
  Bio::DB::SeqFeature::Store
  DateTime
  DBI
  DBD::mysql
  Getopt::Long
  
=cut
use Data::Dumper;
use warnings;
use strict;

# used for timing
use Benchmark;
# These modules are used for parsing, accessing, manipulating in memory, Bio objects
use Bio::DB::SeqFeature::Store;
use Bio::DB::SeqFeature;
# This module is used for parsing BLAST results
use Bio::SearchIO;
# Used in benchmarking/prediction of end time
use DateTime;
# Used for multiple BLAST matches analysis
use DBI;
use DBD::mysql;
# For parameters in running the program
use Getopt::Long;
use lib '.';
# Some subroutines repeatedly used
use InferAnnotationUtilities;

# debug levels
# 0 - no
# 1 - only those benchmarks
# 2 - all
my $debug=1;
my $blast6Table = 'blast6out';
my $database = '--insert-your-database-here--';
my $newSourceName = '';

my $communism;
my $fascism;
my $communismGFF;
my $fascismGFF;
my $geneListHandle;
my $communismHandle;
my $fascismHandle;
my $communismGFFHandle;
my $fascismGFFHandle;
my $communismOutHandle;
my $fascismOutHandle;
my $help;
my %genes;
# Number of extra bases to include from left and end of gene sequence in original. Default 100
my $endBuffer = 100;

my @communismGenes;
my @fascismGenes;
my $communismDB_handle;
my $fascismDB_handle;

# The ultimate file we want;
my $singleFinalGFF3 = 'sbicolor_2.1_gene_exons' . '.inferred-on-' ."--insert--genome--here" . '.gff3';
# Used for storing gene target when searching for appropriate BLAST match and storing
# that to central hash
my $geneTarget;
my %centralHash;

my $grandBenchMarkStart;
my $grandBenchMarkEnd;
my $preLoadBenchMark;
my $processProperBenchMark;
my $singleGFF3BenchMarkStart;
my $singleGFF3BenchMarkEnd;

# splash/guide info first
splash();
my $opt_success = GetOptions(
   'help' => \$help,
#   'genes' => \$geneList,
   'genome1=s' => \$communism,
   'genome2=s' => \$fascism,
   'gff1=s' => \$communismGFF,
   'gff2=s' => \$fascismGFF,
   'endbuffer=s' => \$endBuffer,
   'debug=s' => \$debug,
   'newsource=s' => \$newSourceName
);

die guide() if ( $help || ! $opt_success );
$grandBenchMarkStart = new Benchmark; 
handle_message( 'NOTICE', '', "Opening genome1 : $communism ...");
handle_message( 'FATAL', 'cannot_open_file', $communism ) unless open( $communismHandle, $communism );
handle_message( 'NOTICE', '', "Opening genome2 : $fascism ...");
handle_message( 'FATAL', 'cannot_open_file', $fascism ) unless  open( $fascismHandle, $fascism );
handle_message( 'NOTICE', '', "Opening gff1 : $communismGFF ...");
handle_message( 'FATAL', 'cannot_open_file', $communismGFF ) unless open( $communismGFFHandle, $communismGFF );
handle_message( 'NOTICE', '', "Opening gff2 : $fascismGFF ...");
handle_message( 'FATAL', 'cannot_open_file', $fascismGFF ) unless  open( $fascismGFFHandle, $fascismGFF );

handle_message( 'NOTICE', '', 'Loading genome 1 into memory ...');
eval {
   $communismDB_handle  = Bio::DB::SeqFeature::Store->new(
     -adaptor => 'memory',
     -fasta   => $communism,
     -gff     => $communismGFF
    );
};

if ($@) {
  print $@;
  print "[error] Error opening $communism or $communismGFF !\nProgram exiting.\n";
  exit 1;
}

handle_message( 'NOTICE', '', 'Loading genome 2 into memory ...');
eval {
   $fascismDB_handle  = Bio::DB::SeqFeature::Store->new(
     -adaptor => 'memory',
     -fasta   => $fascism,
     -gff     => $fascismGFF
    );
};

if ($@) {
  print $@;
  print "[error] Error opening $fascism or $fascismGFF !\nProgram exiting.\n";
  exit 1;
}
@communismGenes = $communismDB_handle->get_features_by_type('gene');
@fascismGenes   = $fascismDB_handle->get_features_by_type('gene');
print "[notice] " . localtime() . " Sorting and converting...\n";

my ($temp_arr_ref, $queryHashRef) = convert_BioDBSeqFeature_To_Hash( \@fascismGenes, 1 );
my ($temp_arr_ref2, $queryHashRef2) = convert_BioDBSeqFeature_To_Hash( \@communismGenes, 1 );
my @queryFeatureKeys = @$temp_arr_ref; # the gene names
my %queryFeaturesFast = %$queryHashRef;
my %queryFeaturesFast2 = %$queryHashRef2;
my $ax = 0;
my $bx = scalar(@queryFeatureKeys);
my $searchProperLoopSecs;
$preLoadBenchMark = new Benchmark;

handle_message(
  'NOTICE',
  'Benchmark, Pre-loading and processing took ' . timestr( timediff( $preLoadBenchMark, $grandBenchMarkStart ) ),
    ''
);

foreach my $hashKey (@queryFeatureKeys) {
  if ( ( $ax % 5 ) == 0 ) { 
    print "[notice] " . localtime() . " Progress $ax/$bx  " . ($ax/$bx) * 100 ; 
    if ($debug) { print "\n"; }else{ print "\r" };
  }  
  #
  # $communismGene/$fascismGene is object of type Bio::DB::SeqFeature
  #
  # Again, communist/left is the original, fascist/right is the new genome
  #
  my $fascistGeneSeqFile;
  my $communismGene = $queryFeaturesFast2{ $hashKey };
  my $fascismGene = $queryFeaturesFast{ $hashKey };
  my $fascismGeneSequence;  
  # For benchmarking performance
  my $start_time = new Benchmark;
  my $end_time; 
  my $difference;
  my $timeNotif;
  my $secsRemaining;
  my $averageTimeSoFar;
  my $remainingTimeSoFar;

  # Let's add a few base pairs to the start and end of the new genome gene area
  my $fascistGeneStart = $fascismGene->start - $endBuffer;
  my $fascistGeneEnd   = $fascismGene->end + $endBuffer;
  # This is for the case where the gene is ridiculously present near the start of the chromosome.
  if ( $fascistGeneStart < 0 ) {  $fascistGeneStart = 0; }
  # Now, how do we check for the length of the chromsome, and check if $fascistGeneEnd exceeds that?
  # Well, I don't think we'll arrive at that any time soon, so that's a bug for now.

  # Let's name it
  $fascistGeneSeqFile = "_tmp_targetseq_MNLICN_" .  $fascismGene->load_id . ".tmp.fa";

  $fascismGeneSequence = $fascismDB_handle->fetch_sequence( -seq_id=>$fascismGene->ref, -start=>$fascistGeneStart, -end=>$fascistGeneEnd );
# print to fasta, the target gene sequence (of new genome), with $endBuffer allowance on both ends.
  open( $fascismOutHandle, "> $fascistGeneSeqFile" ) or die ("Cannot open FASTA for target gene sequence.\n");
  print $fascismOutHandle  '>' . $fascismGene->ref . ":" . $fascistGeneStart  ."-" . $fascistGeneEnd . "\r\n";
  print $fascismOutHandle  $fascismGeneSequence . "\r\n";
  close $fascismOutHandle;
  setGeneTarget( $fascismGene->ref );
  searchBLASTProper( 0,  $fascistGeneSeqFile, $communismGene, $communismDB_handle, $fascistGeneStart, $fascistGeneEnd, \%centralHash, '' );
  clearGeneTarget();
  $ax++;
  $end_time = new Benchmark;
  $difference = timediff( $end_time, $start_time );
  $timeNotif =  "[debug] All processes for gene " . $communismGene->load_id . " took " . timestr( $difference ) . "\n"; 
  $timeNotif =~ /\s+(\d+)\s+wallclock\s+secs/;
  $searchProperLoopSecs += $1;
  $averageTimeSoFar = ( $searchProperLoopSecs / $ax );
  $remainingTimeSoFar = ( $averageTimeSoFar * ( $bx - $ax ) );
  if ( $debug > 0 ) { 
   my $dt = DateTime->now();
   my $expectedArrival = $dt->add( seconds => $remainingTimeSoFar );
   print $timeNotif;
   print "[debug] $ax features took $searchProperLoopSecs for an average of $averageTimeSoFar. Estimated arrival at " . $expectedArrival->datetime() . "\n" ;
  }
}
$processProperBenchMark = new Benchmark;
handle_message(
  'NOTICE',
  'Finished, Proper processing of features took ' . timestr( timediff( $processProperBenchMark, $preLoadBenchMark ) ),
  ''
 );

# Let's delete the older file if any
# Linux/Unix specific
if ( -e $singleFinalGFF3 ) {
  if ( execCommandEx(
	    "rm $singleFinalGFF3",
	    "Deleting any previous final GFF3...",
	    "Something went wrong with deleting previous final GFF3"
	  ) != 0
  ){
    die( handle_message( 'FATAL', 'Deletion of existing final GFF3 error, cannot proceed', '' ) );
  }
}

# Now, concatenate them to single file.
$singleGFF3BenchMarkStart = new Benchmark;
foreach my $hashKey (@queryFeatureKeys) {
  my $communismGene = $queryFeaturesFast2{ $hashKey };
#  outputToGFF3( 0, $communismGene, \%centralHash );
  writeGrandFinalOutputToGFF3(  $communismGene, $singleFinalGFF3 );
}
#=cut
$grandBenchMarkEnd = new Benchmark;
handle_message(
  'NOTICE',
  'Finished, Writing to single GFF3 took ' . timestr( timediff( $grandBenchMarkEnd, $singleGFF3BenchMarkStart ) ),
  ''
 );
handle_message( 
  'NOTICE', 
  'Finished, program took ' . timestr( timediff( $grandBenchMarkEnd, $grandBenchMarkStart ) ), 
  ''
 );

# ==end of main==

sub buildBLAST6LineHash {
=doc
  Accepts a line from BLAST result outfmt 6 and builds a hash accoridngly.
  
  Arguments:
    0 - string. The line/one result or match from BLAST results file.
	
  Returns:
    Reference to a hash. The hash then is accessed by BLAST6/Table labels
=cut

  my @temparr;
  my %thisHash;
 
  chomp $_[0];
  @temparr = split( /\s+/, $_[0] ); 
  $thisHash{ 'query_label' } = $temparr[ 0 ];
  $thisHash{ 'target' } = $temparr[ 1 ];
  $thisHash{ 'percent_identity' } = $temparr[ 2 ];
  $thisHash{ 'alignment_length' } = $temparr[ 3 ];
  $thisHash{ 'mismatch' } = $temparr[ 4 ];
  $thisHash{ 'gap' } = $temparr[ 5 ];
  $thisHash{ 'query_start' } = $temparr[ 6 ];
  $thisHash{ 'query_end' } = $temparr[ 7 ];
  $thisHash{ 'target_start' } = $temparr[ 8 ];
  $thisHash{ 'target_end' } = $temparr[ 9 ];
  $thisHash{ 'evalue' } = $temparr[ 10 ];
  $thisHash{ 'bitscore' } = $temparr[ 11 ];
  return \%thisHash;
} # sub 

sub ConnectToMySql {
=doc
  Used for connecting to MySQL database.
  Disclaimer: This function was sourced from somewhere else, I forgot where. All credits go to 
    original poster.

  Arguments
    0 - database name
  
  Returns
    DBI object for MySQL connection.
=cut
  my ($db) = @_;

  # assign the values in the accessDB file to the variables
  my $host = "localhost";
  my $userid = "--your-db--username-here--";
  my $passwd = "--your--db--password-here--";

  # assign the values to your connection variable
  my $connectionInfo="dbi:mysql:$db;$host";

  # make connection to database
  my $l_connection = DBI->connect($connectionInfo,$userid,$passwd);

  # the value of this connection is returned by the sub-routine
  return $l_connection;
} # Sub

sub execCommandEx {
=doc
  A subroutine to execute outside programs.

  Arguments:
  0 - string. Required. The command (and arguments to execute).
  1 - string. Required. Message to output upon start of execution.
  2 - string. Optional. Message to output upon execution error.
=cut
 my $pidx;
 my $executeCommand = $_[0];
 my $executeStartMsg = $_[1];
 my $executeErrMsg  =  ! defined $_[2]  ?  $_[2] : 'Execution error.\n';
 
 eval {
   $pidx = fork || exec $executeCommand;
   handle_message( 'NOTICE', $executeStartMsg . ' | PID ' . $pidx, '' );
 };
 if ($?){
   handle_message( 'ERROR', $executeErrMsg  . ' | PID ' . $pidx , '' );
   return $?;
 }
 waitpid( $pidx, 0 );
 
 return 0;
} #sub

sub insertBLAST6ToMySQL {
=doc
  Populates our main table in our own-devised MySQL database for BLAST 6 results manipulation.
 
  Arguments:
    0 - string. The table name.
    1 - reference to a HASH. The hash contains the BLAST matches, and whose keys are the line numbers
      from BLAST outfmt 6 results file.

  Returns:
    Nothing
=cut
  my $connection = ConnectToMySql( $database );
  my $query;
  my %resultsHash = %{ $_[1] };
  my $statement;

  $query = "INSERT INTO `$_[0]` (`id`,`query_label`,`target`,`percent_identity`,`alignment_length`,`mismatch`,`gap`,";
  $query .= "`query_start`,`query_end`,`target_start`,`target_end`,`evalue`,`bitscore`)";
  $query .= "  VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?);";
  $statement = $connection->prepare( "TRUNCATE TABLE `$_[0]`;" );
  $statement->execute();
  foreach my $lineNum ( sort keys %resultsHash ) {
    my %localHash = %{  $resultsHash{ $lineNum } };
    $statement = $connection->prepare( $query );
    $statement->execute(
      $lineNum,
      $localHash{ 'query_label' },
      $localHash{ 'target' },
      $localHash{ 'percent_identity' },
      $localHash{ 'alignment_length' },
      $localHash{ 'mismatch' },
      $localHash{ 'gap' },
      $localHash{ 'query_start' },
      $localHash{ 'query_end' },
      $localHash{ 'target_start' },
      $localHash{ 'target_end' },
      $localHash{ 'evalue' },
      $localHash{ 'bitscore' }
    )  || print "$DBI::errstr \n" ;
  }
} # sub


sub populateFeatureFinal {
=doc
  Adds features to one central hash and the final locations in new genome.

  Arguments:
    0 - reference to a hash. The central, multi-dimensional  hash.
    1 - string. The feature's parent's load_id.
    2 - string. The feature ID/name, or in Bio::DB::SeqFeature parlance, 'load_id'
    3 - Target locus ( i.e., chromosome )
    4 - new start
    5 - new end
  
  Returns:
    Nothing
=cut
  my %centralHashX = %{ $_[0] };
  my $target;

  # Fool proofing target
  if ( $_[3] =~ /:/ ) {
    #Ex: "Chr01_pilonV9:11023-15060"
    my @temparr = split( /:/, $_[3] );
    $target = $temparr[0];
  }else{
    $target = $_[3];
  }

  $centralHashX{ $_[2] } = {
    'parent' => $_[1],
    'target' => $target, 
    'start' => $_[4],
    'end' => $_[5]
  };
  print '[debug] [britney]' .  ${ $centralHashX{ $_[2] } }{ 'target' } . "\n";
} # sub


sub writeGrandFinalOutputToGFF3 {
=doc
  @devnote Might be currently limited to Linux systems

  Recursively traverses gene features, concatenates their GFF3 into the final, single GFF3.

  Arguments:
    0 - Bio::DB::SeqFeature object. The original gene object.
    1 - Grand, final, single GFF3 filename.
=cut
  my $communistGene = $_[0];
  my @subfeatures;
  my @sortedSubfeatures;

  # Linux/Unix specific!
  if ( ! ( checkIfRejected(  $communistGene->load_id ) or checkIfWaitlisted( $communistGene->load_id ) or checkIfNotImplemented( $communistGene->load_id )  ) ) {
    my $runResult =  appendLineToFile( 
      getSingleLineFromFile(
       getGFF3TmpFileName( $communistGene->load_id  )     
      ),
      $_[1]
    );
    if ( ! ( $runResult  == 0 or $runResult == 256 ) ){
      die( handle_message( 'FATAL', 'Something went wrong with writing ' . $communistGene->load_id, '' ) );
    }    
  }else{
    if( $debug > 1 ) { print '[debug] ' . $communistGene->load_id . 'is rejected!' . "\n"; }
  }
  @subfeatures = $communistGene->get_SeqFeatures();
  @sortedSubfeatures = sort_BioDBSeqFeature_by_PrimaryID( \@subfeatures );
  foreach my $subfeature (@sortedSubfeatures)
  {
    writeGrandFinalOutputToGFF3( $subfeature, $_[1] );
  }
} # sub

sub outputToGFF3 {
  # ( 0, $communismGene, \%centralHash );
=doc
  @devnote To be revisited later, when writeGrandFinalOutputToGFF3() is to be retired.

  Recursively traverses gene features, and outputs the final GFF3

  Arguments
    0 - int. Recursion level.
    1 - The feature in question
    2 - reference to a hash.

  Returns:
    nothing
=cut
  my $recursionLevel = $_[0];
  my $communistGene = $_[1];
  my %fascistData = %{ $_[2] };
  my %hashie = $fascistData{ $communistGene->load_id };
  my  @temparr;
  my $type;
  my $source = '.';  
  my $name;
  my @subfeatures;
  my @sortedSubfeatures; 
  
  if( $recursionLevel == 0 ) {
    $name = $communistGene->name; 
  }else{
    $name = $communistGene->load_id;
  }
  if( $communistGene->type =~ /:/ ) {
    @temparr = split( /:/,  $communistGene->type );
    $type = $temparr[0];
    $source = $temparr[1];
  }else{
    $type = $communistGene->type;
  }
    @subfeatures = $communistGene->get_SeqFeatures();
    @sortedSubfeatures = sort_BioDBSeqFeature_by_PrimaryID( \@subfeatures );
    foreach my $subfeature (@sortedSubfeatures)
    {
      outputToGFF3( $recursionLevel + 1, $subfeature, \%fascistData ); 
    }
} # sub

sub searchBLASTProper {
=doc
  Recursively finds gene features, and BLAST-s them on the target sequence.

  Arguments:
    0 - int. recursion level
    1 - string. The file name of the FASTA containing suspected gene sequence in new genome, prepended and appended with extra $endBuffer bases whenever possible.
    2 - Bio::DB::SeqFeature object. The original gene object (in original genome and annotation).
    3 - reference. reference to Bio::DB::SeqFeature handle
    4 - int. Start base pair of the target gene (at gene level)
    5 - int. End base pair of the target gene (at gene level)
    6 - reference to a hash. The main hash where the final coordinates are to be stored.
    7 - string. Load ID of the feature's parent. Not applicable to type gene.
=cut
 my $recursionLevel = shift;
 my $fascistGeneSeqFile = shift;
 my $communismGene =  shift ;
 my $communismDB_Handle = shift;
 my $geneLevelStart = shift;
 my $geneLevelEnd = shift;
 my $centralHashRef = shift;
 my $parentLoadID = shift;
 my @subfeatures;
 my @sortedSubfeatures;
 my %foundOut;
 my $pid;
 my %bestResult;
 my $orig_len =  $communismGene->end - $communismGene->start;
 my $new_len;
 my $samex;
 my $new_start;
 my $new_end;

 # Bullet-proofing first
 if ( $recursionLevel < 0 ) {
   print "[error][annoyed] recursionLevel < 0 is not meant to be!\n";
 }elsif 
 ( $recursionLevel == 0 ) {
   # Index first the target gene sequence in FASTA file
   my $blastDBIndexCMD = "makeblastdb -in $fascistGeneSeqFile -dbtype nucl";   
   if ( execCommandEx(
          $blastDBIndexCMD, 
          "Indexing FASTA for $fascistGeneSeqFile", 
          "Something went wrong with BLAST indexing" 
        ) != 0
    ){
      die( handle_message( 'FATAL', 'FASTA indexing for target sequence error, cannot proceed', '' ) ); 
    }
 } #endif

 # Get gene subfeatures
 @subfeatures = $communismGene->get_SeqFeatures();
 # so the tree has been fully traversed and reached a leaf
 if( scalar( @subfeatures ) < 1 ) {
    # No more subfeatures, therefore BLAST this one.
    if ($debug > 1 ) { print "[debug]\t\tNo more features under " . $communismGene->load_id .  ".\n"; 
      print "[debug] Will BLAST on $communismGene? ";
    }
    # MNLICN is supposed to be a unique, dynamically generated tag to identify from which run these
    # files were produced. But anyway, let's do that later.
    # BTW, that stands for MNL-ICN (Manila - Incheon), ahem ahem Daehan Hanggong and Asiana Hanggong,
    # trip to Korea please. :D 
    my $blastResultFile = "_tmp_blastresult_MNLICN_" . $communismGene->load_id . ".tmp";
    my $thisFeatureSequenceFile = "_tmp_featseq_MNLICN_" . $communismGene->load_id . ".fa";
    my $thisFeatureSequenceFileHandle;
    my $thisFeatureSequenceLen;
    my $thisFeatureSequence;
    my $shortTriggered = 0;
    my $blastCommand;
    # If BLAST sequence is too short, have to adjust, so these are added to the command to be executed
    my $blastAdditionForShort = "-evalue 1000 -word_size 7 -dust no ";
    my $pid;
   
 
    # Retrieve the original sequence    
    $thisFeatureSequence = $communismDB_Handle->fetch_sequence( 
      -seq_id=>$communismGene->ref,
      -start=>$communismGene->start,
      -end=>$communismGene->end
    );
    # Output to a file as BLAST requires file
    open( $thisFeatureSequenceFileHandle, "> $thisFeatureSequenceFile " ) or die ( "[fatal] Cannot open  $thisFeatureSequenceFile for writing.\n" );
    print $thisFeatureSequenceFileHandle '>' . $communismGene->ref . ':' . $communismGene->start . '-' . $communismGene->end . "\r\n";
    print $thisFeatureSequenceFileHandle $thisFeatureSequence;
    close(  $thisFeatureSequenceFileHandle );
    $blastCommand = "blastn -task blastn -num_threads 8 -outfmt 6 ";

    $thisFeatureSequenceLen = length( $thisFeatureSequence );
    if( $thisFeatureSequenceLen > 6 ) {
      if ( $debug > 1 ) { print "YES\n"; }
=doc
      If sequence to be BLAST-ed is less than 20 characters, there might not be any hits with the default parameters.
      So, let's adjust them. 

      For more info, visit: http://www.ncbi.nlm.nih.gov/blast/Why.shtml, under the heading,

	"Search for short and near exact matches" under Nucleotide BLAST is useful for primer or short nucleotide motif searches.  
=cut
      if( $thisFeatureSequenceLen < 20 ) {
 	 $blastCommand .= $blastAdditionForShort;
         $shortTriggered = 1;
         if ( $debug > 1 ) {  print "[debug] Short sequence parameters for BLAST activated .\n"; }
      }
      $blastCommand .= "-query $thisFeatureSequenceFile -db $fascistGeneSeqFile > $blastResultFile";
      if ( execCommandEx(  $blastCommand, "BLAST process launched ", " Something went wrong with BLAST!" ) != 0 ) {
        exit(2);
      }
      # Now check BLAST results
      %bestResult = parseBLAST6Results( $blastResultFile, $communismGene->load_id, 0 );
      if ( $bestResult{ 'target_start' } == -4 and $bestResult{ 'target_end' } == -4 ) {
      # NO BLAST result at first try? One More Chance please. ( <expletive> naman Bash, ganyan ka ba katigas? :D )
        my $dontReturn = 0;
        # 62 was chosen becaue of some feature 
        if ( $thisFeatureSequenceLen > 19 and $thisFeatureSequenceLen < 62  ) {
          if ( ! $shortTriggered ) {
            handle_message('WARNING','','No BLAST match on sequence with length ' . $thisFeatureSequenceLen . '? Calling again with short seq params...');
            # Insert the parameters for short sequences
            my $tmpquotemeta = quotemeta( $blastAdditionForShort . '-query ');
            $blastCommand =~ s/\-query\s+/$tmpquotemeta/;
            # Execute again
            $blastCommand =~ s/\\//g;
            if ( execCommandEx(  $blastCommand, "BLAST process launched ", " Something went wrong with BLAST!" ) != 0 ) {
              exit(2);
            }
            # Parse again
            %bestResult = parseBLAST6Results( $blastResultFile, $communismGene->load_id, 1 );
            # Now, if no result again, then really, this is rejected, sequence in genome really changed a lot
            if ( $bestResult{ 'target_start' } == -5 and $bestResult{ 'target_end' } == -5 ) {
              handle_message('ERROR','','No BLAST match at all for ' . $communismGene->load_id . ' even when parameters adjusted to for shorter seqs' );
            }else{
              $dontReturn = 1;
            }
          }
        }
        if( !$dontReturn ) {
          addToRejected( $communismGene->load_id, "BLASTMATCH404\t$thisFeatureSequenceLen\t" .  $bestResult{ 'target_start' }  );
	  return ( \%bestResult,  $bestResult{ 'target_start' }, $bestResult{ 'target_end' } );
        }
      }
      # The decrease by 1 is because of negating the effects of adding  $geneLevelStart in target sequence FA db  
      $new_start = ( $bestResult{ 'target_start' } + $geneLevelStart ) -1 ;         
      $new_end =  ( $bestResult{ 'target_end' } + $geneLevelStart ) - 1;
    }else{
       if ( $debug > 1 ) { print "NO\n"; }
       %bestResult = manualInferFeature( $communismGene, $communismDB_Handle );
       if(  $bestResult{ 'target_start' } == -2 and  $bestResult{ 'target_end' } == -2 ) {
         handle_message('WARNING','annotation for ' . $communismGene->load_id . ' does not make sense. skipping!' ,'');
         addToRejected( $communismGene->load_id, 'NONSENSICAL_ANNOT' );
         return ( \%bestResult,  $bestResult{ 'target_start' },  $bestResult{ 'target_end' } ) ;
       }elsif
       ( $bestResult{ 'target_start' } == -3 and  $bestResult{ 'target_end' } == -3 ){
	 addToWaitlisted( $communismGene->load_id );
         handle_message('WARNING','annotation for ' . $communismGene->load_id . ' on waitlist, pending parent completion...', '' );
         return ( \%bestResult,  $bestResult{ 'target_start' },  $bestResult{ 'target_end' } ) ;
       }else{
	 # No need for the addition of $geneLevelStart and decrease of 1 here since we're taking cues
	 # from the previous feature
	 $new_start = $bestResult{ 'target_start' };
	 $new_end = $bestResult{ 'target_end' };       
       }
    }

    # Check if hash is empty
    if ( !keys %bestResult ) {
      handle_message('ERROR', 'No keys for \%bestResult ' . $communismGene->load_id );
      addToRejected( $communismGene->load_id, 'HASHKEY404' );
      return ( \%bestResult, -1, -1 );
    }else{
      # output to its own GFF3 first
      writeGFF3Tmp( $communismGene, getGeneTarget(), $new_start, $new_end, $parentLoadID );
    }
#    return the best one, start and end
    return ( \%bestResult, $new_start, $new_end ) ;
 }else{
   # I guess there is not a genome that is larger than this value?
   my $lowerLimit = 9999999999999999999999;
   my $upperLimit = -1;
   my %returningHash;
   my %hashToBeReturned;
  
    
   if ( $debug > 1 ) { print "[notice] There are stil subfeatures under " .  $communismGene->load_id . "\n"; }
   # sort first
   @sortedSubfeatures = sort_BioDBSeqFeature_by_PrimaryID( \@subfeatures );
   foreach my $subfeature (@sortedSubfeatures)
   {   
      ( $foundOut{ $subfeature->load_id  }, my $tempX, my $tempY ) = searchBLASTProper( $recursionLevel + 1, $fascistGeneSeqFile, $subfeature, $communismDB_Handle, $geneLevelStart, $geneLevelEnd, $centralHashRef, $communismGene->load_id );
      if ( $tempX == -2 and $tempY == -2 ) {
        # This is met when manualInferFeature() is called and the annotation at that point does not make sense.
        handle_message('WARNING',$subfeature->load_id . ' will not be included due to nonsensical annotation','');
        next;
      }elsif
      ( $tempX == -3 and $tempY == -3 ) {
        # This is met when manualInferFeature() is called and the annotation for that feature depends on its parent feature, 
        # that, at this point, the output GFF3 is not yet written
        handle_message('WARNING',$subfeature->load_id . ' will not be written yet due to parent feature not analyzed yet','');
        next;
      }elsif
      ( ( $tempX == -4 and $tempY == -4 ) or ( $tempX == -5 and $tempY == -5 ) ) {
        # when no BLAST match really... :'(
        handle_message('WARNING',$subfeature->load_id . ' will not be included due to zero BLAST match', '' );
        next;
      }
      if( $tempX < $lowerLimit ) {  $lowerLimit = ( $tempX); }
      if( $tempY > $upperLimit ) {  $upperLimit = ( $tempY ); }
   }
#=test_start
   $new_len = $upperLimit - $lowerLimit;
   if ( $debug > 1 ) {
   print "[gff3]<target>\t" .  $communismGene->load_id .  "\t" .  $communismGene->type . "\t" .  "\t$lowerLimit\t$upperLimit";
   $samex = ( $orig_len == $new_len ) ? "same=yes" : "same=no";
   print "|" . "orig_len=" . $orig_len. "\t" . "new_len=" . $new_len  . "\t" . $samex . "\n";
   }
#=test_end
   $hashToBeReturned{ 'target_start' } = $lowerLimit;
   $hashToBeReturned{ 'target_end' } = $upperLimit;
#   populateFeatureFinal( $centralHashRef, $parentLoadID, $communismGene->load_id, getGeneTarget(),  $lowerLimit, $upperLimit );

   writeGFF3Tmp( $communismGene, getGeneTarget(),  $lowerLimit, $upperLimit, $parentLoadID );
   return ( \%hashToBeReturned, $lowerLimit, $upperLimit );
 }

} #sub

sub getGeneTarget {
  return $geneTarget;
}

sub setGeneTarget {
  $geneTarget = $_[0];  
}

sub clearGeneTarget {
   $geneTarget = '';
}
 
sub getFileNumLines {
=doc
  Arguments:
    0 - string. The filename of the file.
 
  Returns:
    int. The number of lines in a file. 
=cut
  my $lines = 0;
  my $fileHandle;
  open ( $fileHandle, $_[0] ) or die('ERROR',"Can't open $_[0] for line number checking!");
  $lines++ while( <$fileHandle> );
  close $fileHandle;
  return $lines;
} # getFileNumLines

sub guide{
  print "

  Synopsis:

  $0 [options]

  Options:
  --endbuffer Optional. Number of extra bases to include from left and end of gene sequence in original. Default 100.
  --genome1   Required. The original genome, in FASTA format. 
  --genome2   Required. The newly created genome (the one you're generating a new annotation for), in FASTA format.
  --gff1      Required. The GFF3 for --genome1. The subfeatures are required too. As in, complete.
  --gff2      Required. The GFF3 for --genome2. Only those type 'gene' are required.
  --newsource Optional. Specify if we have to have a new source to be put in the GFF3 column 2 (1-based indexing).
                If not specified, uses the original source as specified in --gff1
  ";
}

sub handle_message {
  my ($level, $code, $sentMsg) = @_;
  $sentMsg ||= "No message.";
  chomp $sentMsg;
  my $message .= ( localtime() . ' ' . $sentMsg . "\n" );
  if ($level eq 'FATAL') {
    die join ' : ', ($level, $code, $message);
  }else {
    print STDERR join ' : ', ($level, $code, $message);
  }
} # handle_message

sub parseBLAST6Results {
=doc
  Parses BLAST results where output format is 6 (tab-delimited, ala CSV), gets the best matching
  one as per our crteria and returns it to the caller.
  
  Arguments:
    0 - string. BLAST results output file
    1 - string. Feature (i.e. Gene) name
    2 - int. The zero-index, i-th time this function was called for the particular feature. 
          Useful for re-BLASTing short sequences where there were no results in the regular BLAST 
          (with default parameters).

  Returns:
    Hash accessed by BLAST6/Table labels

    Special note, the 'target_start' and 'target_end' are both assigned the following values as per
      the description:

    -4 - There is no BLAST match at all, for first time this was called
    -5 - There is no BLAST match at all, the second time this was called
=cut
  my @temparr;
  my %thatHash;
  my $lineCount = getFileNumLines( $_[0] );
  my $featureName = $_[1];
  my $timesCalled = $_[2];

# Is result just one line?
  if ( $lineCount == 1 ) {
    return %{
      buildBLAST6LineHash( getSingleLineFromFile( $_[0] ) )       
    };
  }elsif ( $lineCount == 0 ){   
      handle_message('ERROR', "Zero BLAST match for $featureName !", '');
      $thatHash{ -1 } = {        
	'target_start' => ( $timesCalled == 0 ) ? -4 : -5,
	'target_end' => ( $timesCalled == 0 ) ? -4 : -5
      };
      return %{ $thatHash{ -1 } };
  }else{
    if ( $debug > 1 ) { handle_message('NOTICE', 'More than one BLAST match detected for ' . $featureName, ''); }
    my %theResults;
    my $x = 0;
    my @arraytemp;

    # scan matches and put into hash
    open( my $fileHandle, $_[0] ) or die( handle_message('ERROR', "Can't open $_[0] for processing.", '') );
    while ( <$fileHandle> ) {
       $theResults{ $x++ } = buildBLAST6LineHash( $_ );
    }
    close( $fileHandle );
    # insert into database for easier selection
    insertBLAST6ToMySQL( $blast6Table , \%theResults );
    # We reuse $x !
    # Now eventually filter matches according to our criteria, progressively
    $x =  searchBestBLAST6Match( 0, $featureName, \@arraytemp, '', '' ); 
    if ( $x == -1 ) {
      handle_message('ERROR', 'No suitable BLAST match found for ' . $featureName, '' );
    } 
    if( ! exists( $theResults{ $x } ) ) {
       handle_message('ERROR', 'Incorrect BLAST match ' . $x . ' found for ' . $featureName, '' );
    }
    if ( $debug > 1 ) { print "[notice] Successfully found match $x for $featureName\n"; }
    return %{ $theResults{ $x } };
  }
} # parseBLAST6Results

sub searchBestBLAST6Match{
=doc
 @devnote Copied and adjusted from "our_new_script.pl"
 This is a recursive function, progressively filtering out rows. If at the upper level, there are still more than 
 one row, the topmost one will be returned.
 
=cut
  my $recursionLevel = $_[0] ;
  my $geneName = $_[1];  
  my @previousValues = @{ $_[2] };
  my $previousQuerySpecific = $_[3];
  my $previousTable = $_[4];   
  my $sqlQueryOrderBy; 
  my $sqlQuerySpecific;
  my $dispatchTable;
  my $connection = ConnectToMySql($database);
  my $statement;
  my $thisRecursionFeaturedValue;
  my $thisRecursionFeaturedColumn;
  my $thisTempTable;
  my $returnThis;

  if ( $debug > 1 ) { print "\t Recursion Level $recursionLevel \n"; }
  $dispatchTable = {
    0 => sub {
       $thisTempTable = "lessstrict_0";
       $thisRecursionFeaturedColumn = 'evalue';
       $previousTable = $blast6Table;
       $sqlQueryOrderBy = "SELECT DISTINCT($thisRecursionFeaturedColumn) FROM `$blast6Table` ORDER BY `$thisRecursionFeaturedColumn`;";
       $sqlQuerySpecific = "SELECT * FROM `specifictable` WHERE `$thisRecursionFeaturedColumn` = ? ";
      },
    1 => sub {
       $thisTempTable = "lessstrict_1";
       $thisRecursionFeaturedColumn = 'alignment_length';
       $sqlQueryOrderBy = "SELECT DISTINCT($thisRecursionFeaturedColumn) FROM $previousTable t ORDER BY t.$thisRecursionFeaturedColumn DESC; ";
       $sqlQuerySpecific = $previousQuerySpecific . " AND $thisRecursionFeaturedColumn = ? "; 
    },
    2 => sub {
       $thisTempTable = "lessstrict_2";
       $thisRecursionFeaturedColumn = 'percent_identity';
       $sqlQueryOrderBy = "SELECT DISTINCT($thisRecursionFeaturedColumn) FROM $previousTable u ORDER BY u.$thisRecursionFeaturedColumn DESC; ";
       $sqlQuerySpecific =  $previousQuerySpecific . " AND $thisRecursionFeaturedColumn = ? ";
    },
    3 => sub {
       $thisTempTable = "lessstrict_3";
       $thisRecursionFeaturedColumn = 'bitscore';
       $sqlQueryOrderBy = "SELECT DISTINCT($thisRecursionFeaturedColumn) FROM $previousTable v ORDER BY v.$thisRecursionFeaturedColumn DESC; ";
       $sqlQuerySpecific =  $previousQuerySpecific . " AND $thisRecursionFeaturedColumn = ? ";
    },
    4 => sub {
       $thisTempTable = "lessstrict_4";
       $thisRecursionFeaturedColumn = 'mismatch';
       $sqlQueryOrderBy =  "SELECT DISTINCT($thisRecursionFeaturedColumn) FROM $previousTable w ORDER BY w.$thisRecursionFeaturedColumn; ";
       $sqlQuerySpecific = $previousQuerySpecific . " AND $thisRecursionFeaturedColumn = ? ";
    },
    5 => sub {
       $thisTempTable = "lessstrict_5";
       $thisRecursionFeaturedColumn = 'gap';
       $sqlQueryOrderBy =  "SELECT DISTINCT($thisRecursionFeaturedColumn) FROM $previousTable x ORDER BY x.$thisRecursionFeaturedColumn; ";
       $sqlQuerySpecific = $previousQuerySpecific . " AND $thisRecursionFeaturedColumn  = ? ";
    },
    6 => sub {
      $thisTempTable = "lessstrict_6";
      print "[fatal] Exhausted all filters for $geneName at recursionLevel $recursionLevel ,still duplicated .. | $previousQuerySpecific \n";
      return "";
    },
   'default' => sub {
      print "[error] Unrecognized $recursionLevel \n";
      return "";
   }
  };
  # Execute dispatch table
  $dispatchTable->{ $recursionLevel }->();
  $statement = $connection->prepare(  $sqlQueryOrderBy  );
  $statement->execute() || print "$DBI::errstr \n";
  
  # Get the topmost row  
  my $rows =  $statement->rows;
  if( $rows < 1 ) {
    print "[fatal] At gene $geneName recursionLevel $recursionLevel No rows at ORDER by part \n";
    $returnThis = -1;
  } 
  # save row in a hash
  my $ref = $statement->fetchrow_hashref() ;
  # get the best value for this recursion
  $thisRecursionFeaturedValue = $ref->{ $thisRecursionFeaturedColumn };

  # some string changes so we can substitute the SQL query, and select only the appropriate data
  # to be progressively filtered
  my $localQuerySpecific =  $sqlQuerySpecific;
  my $quotemetaED = quotemeta( $previousTable );
  $localQuerySpecific =~ s/specifictable/$quotemetaED/;
  my $statement2 = $connection->prepare( "INSERT INTO `$thisTempTable` ( $localQuerySpecific  ) ;"  );
  push( @previousValues, $thisRecursionFeaturedValue );
  $statement2->execute( @previousValues )  || print "$DBI::errstr \n";;
  
  my $statement3 =  $connection->prepare("SELECT * FROM $thisTempTable ;");	
  $statement3->execute();	

  # zero condition not possible (this function won't be called if that's the case), it's either 1 or > 1
  my $s3rows =  $statement3->rows;
  if( $s3rows == 1 )
  {
    my $finalRef = $statement3->fetchrow_hashref();
    $returnThis = $finalRef->{ 'id' };
  }else{
    if ( $recursionLevel < 5 ) {
      if ( $debug > 1) { print "[notice] For $geneName $thisRecursionFeaturedColumn  $thisRecursionFeaturedValue is still not enough. Leveling up...\n"; }
      $returnThis = searchBestBLAST6Match( $recursionLevel + 1,  $geneName, \@previousValues, $sqlQuerySpecific, $thisTempTable);
    }else{
       # STill there are more results in the topmost criteria, return the top most
      my $finalRef = $statement3->fetchrow_hashref();
      if ( $debug > 1) { handle_message('WARNING', 'There are still more than one results in recursionLevel 5/max criterion. Returning ' . $finalRef->{ 'id' } . " for $geneName", ''); }
      $returnThis = $finalRef->{ 'id' };
    }
  }
  # TRUNCATE ANY TABLE TEMPORARY
  my $statement5 = $connection->prepare("TRUNCATE TABLE $thisTempTable ;");
  $statement5->execute()  || print "$DBI::errstr \n";;
  return $returnThis;
}#sub

sub sort_BioDBSeqFeature_by_PrimaryID
{
  #
  # Sorts an array that contains Bio::DB::SeqFeature objects via
  # their primary_ID,
  #
  # Arguments:
  #   0 - array reference.
  # Returns ARRAY. Of the sorted objects. Key starts at 0.
  #
#  if($debug){ print "[debug] " . localtime() . " Sorting some Bio::DB::SeqFeature object.\n"; }
  my @featuresx = @{ $_[0] };
  my %tempHash = ();
  my @sortedFeatures;
  my @sortedKeys;
  my $start_time = new Benchmark;
  my $arr_ref;
  my $hash_ref;

  ( $arr_ref, $hash_ref ) = convert_BioDBSeqFeature_To_Hash( \@featuresx, 0 );
  @sortedKeys = @$arr_ref;
  %tempHash = %$hash_ref;
  foreach my $hashKey( @sortedKeys )
  {
    my $feature =  $tempHash{ $hashKey };
    push( @sortedFeatures, $tempHash{ $hashKey } );
  }
  my $end_time = new Benchmark;
  my $difference = timediff( $end_time, $start_time );
  if($debug > 1){
    print "[debug] ". localtime() . " Sorting finished for sortedKeys " . scalar( @sortedKeys ) . "\n";
    print "[debug] Sorting took " . timestr( $difference ) . "\n";
  }
  return @sortedFeatures;
}


sub splash{
  print "\nBLAST and Infer Genes & Subfeatures Script \n";
  print "a.llave\@irri.org 01SEP2014-1514\n";
  print "-------------------------------------\n\n";
}

sub writeGFF3Tmp {
=doc
  Temporary function. Writes GFF3 line for the feature in a single file.

  Arguments:
    0 - Bio::DB::SeqFeature object. The original gene
    1 - string. Target locus (i.e. chromosome)
    2 - int. New start
    3 - int. New end
    4 - string. Parent load_id. If none, ''.
  Returns:
   Nothing
=cut
  my $communistGene = $_[0];
  my $thisFile = getGFF3TmpFileName( $communistGene->load_id );
  my $fileHandle;
  my $showChange = 1;

  my %col9Hash = ();
  my @col9HashKeys = ();
  my $col9HashKeysRef;
  my $col9HashRef;

  my $source;
  my $type;
  my $strand;
  my $phase = '.';
  my $parentLoadID = $_[4];

  # First and foremost, check the absolute values of start and end!!! Much time wasted because of this 
  # not being checked. Yep, at this point I don't trust the other parts of the code. :P
  if( $_[2] <= 0 or $_[2] <= 0 ) {
    handle_message('FATAL','Baloney start and end','For ' . $communistGene->load_id . ' invalid either start ' . $_[2] . ' or end ' . $_[3]  );
  }

  # Get type
  if ( $communistGene->type =~ /:/ ) {
    my @temparr = split ( /:/, $communistGene->type );
    $type = $temparr[0];
    $source = $temparr[1];
  }else{
    $type = $communistGene->type;
#   $source = '.'; 
    $source = $communistGene->source;
  }

  # Get strand
  $strand = getLiteralStrand( $communistGene );

  # How about phase?
  if ( defined( $communistGene->phase ) ){
     $phase = $communistGene->phase;
  }

  open( $fileHandle, ">$thisFile" ) or handle_message( 'FATAL',"Cannot write to $thisFile", '' );
  print $fileHandle $_[1] . "\t$source\t$type\t$_[2]\t$_[3]\t.\t$strand\t$phase\tID=" . $communistGene->load_id . ";Name=";
  if ( defined( $communistGene->name ) ){
    print $fileHandle $communistGene->name;
  }else{
    print $fileHandle $communistGene->load_id;
  }
  if( $parentLoadID ne '' ){
    print $fileHandle ";Parent=" . $parentLoadID ;
  }
  # Now, output the rest of column 9
  ($col9HashKeysRef, $col9HashRef)  = convertGFF3Column9_To_Hash( $communistGene, 0 );
  if( scalar( $col9HashKeysRef ) ){
    @col9HashKeys = @$col9HashKeysRef;
    %col9Hash = %$col9HashRef;
    foreach my $hashKey ( @col9HashKeys ){
      print $fileHandle ";" . $hashKey . "=" . $col9Hash{ $hashKey};
    }
  }
  if ( $showChange ) {
    my $communistGeneLen = $communistGene->end - $communistGene->start;
    my $fascistGeneLen = $_[3] - $_[2]; 
     if( $communistGeneLen != $fascistGeneLen ) {
       print $fileHandle ";change=";
       print $fileHandle ( $fascistGeneLen > $communistGeneLen  ) ? "larger" : "smaller";
       #print $fileHandle "$communistGeneLen..$fascistGeneLen";
     }
  }
  print $fileHandle "\n";
  close $fileHandle;
} # sub

sub getGFF3TmpFileName {
=doc
  Temporary function. Gets the appropriate temp GFF3 file name for the single feature.

  Arguments:
    0 - Feature name

  Returns:
    String. The file name.
=cut
 return "_tmp_gff3_MNLICN_" . $_[0] . ".tmp.gff3";
}

sub manualInferFeature {
=doc
  Manually determines the start and end of a subfeature that is not BLAST-able due
  to its short length.
  
  Arguments:
    0 - Bio::DB::SeqFeature object. The feature concerned.
    1 - Bio::DB::SeqFeature::Store reference/database handle. This is for the original genome/gene.	
 
  Returns:
    A Hash with only keys 'target_start' and 'target_end'
=cut

  my $communistFeature = $_[0];
  my $communistHandle = $_[1];
  my $communistFeatureType;
  my $communistFeatureStrand;
  my $previousFeatureObj;
  my $previousFeatureGFF3;
  my $fileHandle;  
  my $loneStr;
  my @temparr;
  my %returnThisHash;
  my %properOrdering;
  my $previousObjSupposedTypeData;
  my $previousObjSupposedType;

  # Decide later on what to do according to these combinations
  %properOrdering = %{
    {
      '-' => {
	# As seen in  Sobic.001G023600, PAC:28394962.five_prime_UTR.2 and PAC:28394962.CDS.1
	'five_prime_UTR' => { 'type' => 'CDS', 'area' => 'end', 'add' => 1 },
	# As seen in Sobic.001G023600, PAC:28394962.three_prime_UTR.1 and PAC:28394962.exon.5
	'three_prime_UTR' => { 'type' => 'exon', 'area' => 'start', 'add' => 0 },
	'CDS' => { 'type' => 'three_prime_UTR', 'area' => 'end', 'add' => 1 }
      },
      '+' => {
	# As seen in Sobic.001G091700, PAC:28394669.CDS.1 and PAC:28394669.five_prime_UTR.1
	'CDS' => { 'type' => 'five_prime_UTR', 'area' => 'end', 'add' => 1 },
	# ... PAC:28394669.five_prime_UTR.1 and PAC:28394669.exon.1
	'five_prime_UTR' => { 'type' => 'exon', 'area' => 'start', 'add' => 0 },
	# ... PAC:28394669.three_prime_UTR.1 and PAC:28394669.CDS.4
	'three_prime_UTR' => { 'type' => 'CDS', 'area' => 'end', 'add' => 1 },
        # As seen in Sobic.001G104200, PAC:28394984.exon.1 and PAC:28394984
        'exon' => { 'type' => 'mRNA', 'area' => 'start', 'add' => 0 },
      }
    }
  };  
  $communistFeatureStrand = getLiteralStrand( $communistFeature ); 

  # For some reason, Bio::DB::SeqFeature::Store makes it a colon-pairing, like
  # "three_prime_UTR:phytozome8_0" , so to avoid, this.
  if( $communistFeature->type =~ /:/ ) {
    $communistFeatureType = getLeftInColonPairing( $communistFeature->type );
  }else{
    $communistFeatureType = $communistFeature->type;
  }
  $previousFeatureObj = $communistHandle->get_feature_by_primary_id( $communistFeature->primary_id - 1 );
  if ( ! defined $previousFeatureObj->type ) {
    handle_message( 'FATAL', 'No previous feature to base on, in manualInferFeature()', '' ); 
  }
  if ( exists( $properOrdering{ $communistFeatureStrand }->{ $communistFeatureType } ) ) {
    if ( $debug > 0 ) { print "[debug] Fetching previous primary_id " .  ( $communistFeature->primary_id - 1) ."\n" ; }
    $previousObjSupposedType = quotemeta( $properOrdering{  $communistFeatureStrand  }->{ $communistFeatureType }->{ 'type'} );

    # This applies if previous feature is its parent feature    
    if ( ! -e  getGFF3TmpFileName( $previousFeatureObj->load_id ) ) {
      my $col9HashKeysRef;
      my $col9HashRef;
      my %col9Hash;

      ($col9HashKeysRef, $col9HashRef)  = convertGFF3Column9_To_Hash( $communistFeature, 1 );
      %col9Hash = %{ $col9HashRef };
      print Dumper %col9Hash;

      # Test 1, scenario between PAC:28394984.exon.1 and mother, with Name=Sobic.001G104200, ID=PAC:28394984
      if( $col9Hash{ 'parent_id' } eq $previousFeatureObj->load_id ) {
        $returnThisHash{ 'target_start' } = -3;
        $returnThisHash{ 'target_end' } = -3;
        return %returnThisHash;     
      }
      # Now, if load_id and parent_id are the same this time         
    } 
    if ( $previousFeatureObj->type !~ /$previousObjSupposedType/ ) {
      handle_message( 'ERROR', "Previous feature other than of type $previousObjSupposedType strand $communistFeatureStrand for "  . $communistFeatureType . ' not yet implemented in manualInferFeature(), requested ' .  $previousFeatureObj->type, '' );
      addToNotImplemented( $communistFeature->load_id, $communistFeature->type, $previousFeatureObj->load_id, $previousFeatureObj->type, $previousObjSupposedType );
      $returnThisHash{ 'target_start' } = -2;
      $returnThisHash{ 'target_end' } = -2;
      return %returnThisHash;
    }
    $previousFeatureGFF3 = getGFF3TmpFileName( $previousFeatureObj->load_id );
    open ( $fileHandle, $previousFeatureGFF3 ) or 
    die( handle_message('FATAL', 'Cannot open for reading - ' . $previousFeatureGFF3, '' ) );
    $loneStr = <$fileHandle>;
    close( $fileHandle );
#    print '[debug] Got previous feature ' . $previousFeatureGFF3 . "\n";
    chomp( $loneStr );
    @temparr = split( /\s+/, $loneStr );
    # GFF3 format, remember? This one's either start or end
    $returnThisHash{ 'target_start' } = ( $properOrdering{ $communistFeatureStrand }->{ $communistFeatureType }->{ 'area' } eq 'end' ) ? $temparr[4] : $temparr[3] ;
    $returnThisHash{ 'target_start' } += $properOrdering{ $communistFeatureStrand }->{ $communistFeatureType }->{ 'add' } ;
    $returnThisHash{ 'target_end' } = $returnThisHash{ 'target_start' }  + ( $communistFeature->end - $communistFeature->start );

    if ( $debug > 1 ) { print "[warning] Manually inferred location of " . $communistFeature->load_id . " via manualInferFeature() " . $returnThisHash{ 'target_start' } . "\t" .  $returnThisHash{ 'target_end' } . "\n"; }
    return %returnThisHash;
  }else{
    handle_message(
      'WARNING',
      'sub manualInferFeature() not sure how to be implemented for other feature types. Requested type and strand: ' . $communistFeature->type . "|" .  $communistFeatureStrand,
      ''
    );
    # well, reject
    $returnThisHash{ 'target_start' } = -2;
    $returnThisHash{ 'target_end' } = -2;
    return %returnThisHash;
  }
} # sub

sub checkIfRejected { 
  my $meow = "grep -cwP '^$_[0]\\s+' "  . getRejectionFile();
  if ( $debug> 1 ) {  print "[notice] checkIfRejectedPID $? \n"; }
  return ( `$meow` > 0 );
}


sub checkIfNotImplemented {
  my $meow = "grep -cwP '^$_[0]\\s+' "  . getNotImplementedFile();
  return ( `$meow` > 0 );
}


sub checkIfWaitlisted {
  my $meow = "grep -cwP '^$_[0]\$' "  . getWaitlistFile();
  return ( `$meow` > 0 );
}

sub getRejectionFile {
  return "_tmp_rejectionfile.tmp";
}

sub getWaitlistFile {
  return "_tmp_waitlistfile.tmp";
}

sub getNotImplementedFile {
  return "_tmp_notimplementedfile.tmp";
}

sub addToRejected {
=doc
 Arguments:
   0 - Feature name
   1 - Reason for rejection
=cut
   my $meow = "echo $_[0]\t$_[1]>> " .  getRejectionFile();
   `$meow`;
}

sub addToNotImplemented {
=doc
   0 - current feature load_id
   1 - 0's type
   2 - previous feature load_id
   3 - 1's type
   4 - supposed previous feature type
=cut
   my $meow = "echo $_[0]\t$_[1]\t$_[2]\t$_[3]\t$_[4]>> " . getNotImplementedFile();
   `$meow`;
}

sub addToWaitlisted {
   my $meow = "echo $_[0] >> " .  getWaitlistFile();
   `$meow`;
}


