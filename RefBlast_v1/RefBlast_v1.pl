#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
Getopt::Long::Configure qw(no_ignore_case);
use File::Basename;
use FindBin qw($Bin);
use lib "$Bin/.Modules";
use Parameter::BinList;

my ($HelpFlag,$BeginTime);
my $ThisScriptName = basename $0;
my $BinList = BinListGet();
my $FullMode = BinSearch("FullMode",$BinList,1);
my $Mode = BinSearch("Mode",$BinList,1);
my $FullNtDB = BinSearch("FullNtDB",$BinList,1);
my $NtDB = BinSearch("NtDBOfNt",$BinList,1);
my $FullBlastMode = BinSearch("FullBlastMode",$BinList,1);
my $BlastMode = BinSearch("BlastMode",$BinList,1);
my $FullResultFormat = BinSearch("FullResultFormat",$BinList,1);
my $ResultFormat = BinSearch("ResultFormat",$BinList,1);
my $MaxEvalue = BinSearch("MaxEvalue",$BinList,1);
my $MaxSeqNum = BinSearch("MaxSeqNum",$BinList,1);
my $CpuNum = BinSearch("CpuNum",$BinList,1);
my ($MaxEvalue4Short,$Makeblastdb,$Blastn,$File4Seq,$File4InitialResult);
my $HelpInfo = <<USAGE;

 $ThisScriptName
 Auther: zhangdong_xie\@foxmail.com

  This script was used to blast a seq to the reference.
  本脚本主要用于查找序列在基因组中上的比对位置。

 -i      ( Required ) Input file with fasta format;
 -o      ( Required ) Result logging file when nt blast or db prefix when nt db making;
                      在查询时指定结果文件，在创建数据库时指定数据库前缀（注意：包括路径 和 数据库名）。

 -mode   ( Optional ) Mode for software (Choices are: $FullMode, default: $Mode);
 For Blast:
 -db     ( Optional ) The name of database (Choices are: $FullNtDB,default: $NtDB);
 -s      ( Optional ) Mode for blast (Choices are: $FullBlastMode, default: $BlastMode);
 -f      ( Optional ) Format for result (Choices are: $FullResultFormat, default: $ResultFormat);
 -e      ( Optional ) Maximal E-value (default: $MaxEvalue);
 -n      ( Optional ) Maximal number of mapping items (default: $MaxSeqNum);
 -t      ( Optional ) The number of threads (default: $CpuNum);
 -bin    ( Optional ) List for searching of related bin or scripts; 
 -h      ( Optional ) Help infomation;

USAGE

GetOptions(
	'i=s' => \$File4Seq,
	'o=s' => \$File4InitialResult,
	'mode:s' => \$Mode,
	'db:s' => \$NtDB,
	's:s' => \$BlastMode,
	'f:s' => \$ResultFormat,
	'e:f' => \$MaxEvalue,
	'n:i' => \$MaxSeqNum,
	't:i' => \$CpuNum,
	'bin:s' => \$BinList,
	'h!' => \$HelpFlag
) or die $HelpInfo;

if($HelpFlag || !$File4Seq)
{
	die $HelpInfo;
}
else
{
	$BeginTime = ScriptBegin(0,$ThisScriptName);
	my $Dir = dirname $File4InitialResult;
	IfDirExist($Dir);
	
	$BinList = BinListGet() if(!$BinList);
	$Blastn = BinSearch("Blastn",$BinList);
	$Makeblastdb = BinSearch("Makeblastdb",$BinList);
	$NtDB = BinSearch("NtDB",$BinList,1) unless($NtDB);
	$FullBlastMode = BinSearch("FullBlastMode",$BinList,1) unless($FullBlastMode);
	$BlastMode = BinSearch("BlastMode",$BinList,1) unless($BlastMode);
	$Mode = BinSearch("Mode",$BinList,1) unless($Mode);
	$ResultFormat = BinSearch("ResultFormat",$BinList,1) unless($ResultFormat);
	$MaxEvalue = BinSearch("MaxEvalue",$BinList,1) unless($MaxEvalue);
	$MaxEvalue4Short = BinSearch("MaxEvalue4Short",$BinList,1) unless($MaxEvalue4Short);
	$MaxSeqNum = BinSearch("MaxSeqNum",$BinList,1) unless($MaxSeqNum);
	$CpuNum = BinSearch("CpuNum",$BinList,1) unless($CpuNum);
}

if(1)
{
	my $Return = "";
	my $CmdLine = "";
	if($Mode eq "Blast")
	{
		my $Seq = `cat $File4Seq | head -n2 | tail -n1`;
		chomp $Seq;
		my $SeqLen = length($Seq);
		
		$MaxEvalue = $MaxEvalue4Short if(($SeqLen < 25 || $BlastMode eq "blastn-short") && $MaxEvalue < $MaxEvalue4Short);
		$CmdLine = "$Blastn -db $NtDB -query $File4Seq -out $File4InitialResult -max_target_seqs $MaxSeqNum -num_threads $CpuNum -evalue $MaxEvalue -task $BlastMode";
		$CmdLine .= " -word_size 4" if($SeqLen < 25 || $BlastMode eq "blastn-short");
		my $FormatString = "7 qacc qstart qend sstart send length mismatch gapopen pident bitscore evalue stitle";
		$CmdLine .= " -outfmt \"$FormatString\"" if($ResultFormat eq "Simple");
		print "[ CmdLine ] $CmdLine\n";
		$Return = `$CmdLine`;
	}
	elsif($Mode eq "DBMake")
	{
		$CmdLine = "$Makeblastdb -dbtype nucl -in $File4Seq -input_type fasta -parse_seqids -out $File4InitialResult -blastdb_version 5";
		print "[ CmdLine ] $CmdLine\n";
		$Return = `$CmdLine`;
	}
	else
	{
		print "[ Error ] Unknown mode $Mode\.\n";
	}
	chomp $Return;
	print "\n$Return\n" if($Return);
}
printf "[ %s ] The end.\n",TimeString(time,$BeginTime);


######### Sub functions ##########
