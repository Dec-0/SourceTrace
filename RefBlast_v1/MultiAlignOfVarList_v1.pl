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
my ($VarList,$ExtendLen,$File4Log);
my ($Perl,$Reference,$RefBlastScript,$Bedtools);
my $BinList = BinListGet();
my $CpuNum = BinSearch("CpuNum",$BinList,1);
my $HelpInfo = <<USAGE;

 $ThisScriptName
 Auther: zhangdong_xie\@foxmail.com

  This script was used to check the multi alignment status of variants list.

 -i      ( Required ) List of variants, with at least 'chr,pos';
 -a      ( Required ) Extend length ,suggest half the length of single probe;
 -o      ( Required ) File for logging;

 -t      ( Optional ) Thread number (default: $CpuNum);
 -bin    ( Optional ) List for searching of related bin or scripts; 
 -h      ( Optional ) Help infomation;

USAGE

GetOptions(
	'i=s' => \$VarList,
	'a=i' => \$ExtendLen,
	'o=s' => \$File4Log,
	't:i' => \$CpuNum,
	'bin:s' => \$BinList,
	'h!' => \$HelpFlag
) or die $HelpInfo;

if($HelpFlag || !$VarList || !$ExtendLen || !$File4Log)
{
	die $HelpInfo;
}
else
{
	$BeginTime = ScriptBegin(0,$ThisScriptName);
	my $Dir = dirname $File4Log;
	IfDirExist($Dir);
	
	$Perl = BinSearch("Perl",$BinList);
	$Reference = BinSearch("Reference",$BinList);
	$RefBlastScript = BinSearch("RefBlast",$BinList);
	$Bedtools = BinSearch("Bedtools",$BinList);
	$CpuNum = BinSearch("CpuNum",$BinList,1) unless($CpuNum);
}

if(1)
{
	my $tFile4Fa = $File4Log . ".tmp.fa";
	my $tFile4Result = $File4Log . ".tmp.result";
	
	open(VL,"cat $VarList | grep -v ^# | cut -f 1,2 |") or die $! unless($VarList =~ /\.gz$/);
	open(VL,"zcat $VarList | grep -v ^# | cut -f 1,2 |") or die $! if($VarList =~ /\.gz$/);
	open(LOG,"> $File4Log") or die $! unless($File4Log =~ /\.gz$/);
	open(LOG,"| gzip > $File4Log") or die $! if($File4Log =~ /\.gz$/);
	while(my $Line = <VL>)
	{
		chomp $Line;
		my ($Chr,$Pos) = split /\t/, $Line;
		my $From = $Pos - $ExtendLen;
		my $To = $Pos + $ExtendLen;
		`echo -e '$Chr\\t$From\\t$To' | $Bedtools getfasta -fi $Reference -bed - 2>/dev/null > $tFile4Fa`;
		die "[ Error ] Fa file not exit ($tFile4Fa).\n" unless(-s $tFile4Fa);
		`$Perl $RefBlastScript -i $tFile4Fa -o $tFile4Result -s blastn -f Simple -t $CpuNum`;
		die "[ Error ] Result file not exit ($tFile4Result).\n" unless(-s $tFile4Result);
		my $MultiNum = `cat $tFile4Result | grep -v ^# | wc -l`;
		chomp $MultiNum;
		print LOG join("\t",$Line,$MultiNum),"\n";
		
		`rm $tFile4Fa $tFile4Result`;
	}
	close VL;
	close LOG;
}
printf "[ %s ] The end.\n",TimeString(time,$BeginTime);


######### Sub functions ##########
