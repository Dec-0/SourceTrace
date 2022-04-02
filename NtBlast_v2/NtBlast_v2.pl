#!/usr/bin/perl
use strict;
use Getopt::Long;
Getopt::Long::Configure qw(no_ignore_case);
use File::Basename;
use Sort::PureNum;
use Parameter::BinList;
use warnings;

my ($HelpFlag,$BinList,$BeginTime);
my ($FaFile,$Dir,$Prefix,$MaxEvalue,$CpuNum,$UniqFlag);
my ($Blastn,$NtDB,$Gi2TaxIdDB,$TaxId2SpeDB,$TaxNode);
my $HelpInfo = <<USAGE;

 NtBlast_v2.pl
 Auther: zhangdong_xie\@foxmail.com

  This script was used to blast fasta seq into NCBI nt database and calculate the species\' contents.
  这个脚本主要用来将fa格式的核苷酸序列比对到NCBI的nt库，并统计物种组成。

 -fa     ( Required ) A file with nucleotides in fasta format;
 -dir    ( Required ) The directory for result logging;

 -prefix ( Optional ) The prefix of the related files (default: the prefix of the fa file);
 -evalue ( Optional ) Maximal e-value (default: 0.0001);
 -cpu    ( Optional ) Cpu number for blastn (default: 5);
 -uniq   ( Optional ) If taken into account only the most best match (with -uinq only);
 -bin    ( Optional ) List for searching of related bin or scripts;
 -h      Help infomation;

USAGE

GetOptions(
	'fa=s' => \$FaFile,
	'dir=s' => \$Dir,
	'prefix:s' => \$Prefix,
	'evalue:f' => \$MaxEvalue,
	'cpu:i' => \$CpuNum,
	'uniq!' => \$UniqFlag,
	'bin:s' => \$BinList,
	'h!' => \$HelpFlag
) or die $HelpInfo;

if($HelpFlag || !$FaFile || !$Dir)
{
	die $HelpInfo;
}
else
{
	$BeginTime = ScriptBegin();
	
	die "[ Error ] File not exist ($FaFile).\n" unless(-e $FaFile);
	die "[ Error ] Directory not exist ($Dir).\n" unless(-d $Dir);
	unless($Prefix)
	{
		my $BaseName = basename $FaFile;
		my @Cols = split /\./, $BaseName;
		$Prefix = $Cols[0];
	}
	$MaxEvalue = 0.0001 unless($MaxEvalue);
	$CpuNum = 5 unless($CpuNum);
	
	$BinList = BinListGet() if(!$BinList);
	$Blastn = BinSearch("Blastn",$BinList);
	$NtDB = BinSearch("NtDb",$BinList,1);
	$Gi2TaxIdDB = BinSearch("Gi2TaxId",$BinList);
	$TaxId2SpeDB = BinSearch("TaxId2Spe",$BinList);
	$TaxNode= BinSearch("TaxNode",$BinList);
}

if(1)
{
	# 初步比对结果;
	my $InitialResult = $Dir . "/" . $Prefix . ".BlastnInitial.xls";
	if(1)
	{
		my $Return = `$Blastn -db $NtDB -query $FaFile -out $InitialResult -evalue $MaxEvalue -max_target_seqs 1 -num_threads $CpuNum -outfmt "7 qacc qstart qend sstart send length mismatch gapopen pident bitscore evalue stitle"`;
		chomp $Return;
		print "\n$Return\n" if($Return);
		printf "[ %s ] Initial result done.\n",TimeString(time,$BeginTime);
	}
	
	# 过滤及物种标记;
	my $FltResult = $Dir . "/" . $Prefix . ".BlastnFlt.xls";
	if(-e $InitialResult)
	{
		my (%SoloFlag,%TaxIdLog) = ();
		
		open(TRANS,"> $FltResult") or die $!;
		open(INI,"cat $InitialResult | grep -v ^# |") or die $!;
		print TRANS join("\t","#NameOfFasta","E-Value","MatchedLength","TaxID","SubjectName","NameNode","SubjectDescription"),"\n";
		while(my $Line = <INI>)
		{
			chomp $Line;
			my @Cols = split /\t/, $Line;
			my @SubCols = split /\|/, $Cols[11];
			my ($FaName,$GI,$EValue,$MatchLen,$Description) = ($Cols[0],$SubCols[1],$Cols[10],$Cols[5],$SubCols[-1]);
			
			# 每条fa序列只保存一条最优结果;
			if($UniqFlag)
			{
				next if($SoloFlag{$FaName});
				$SoloFlag{$FaName} = 1;
			}
			
			my ($TaxId,$NameNode,$SubjectName) = ();
			if($TaxIdLog{$GI})
			{
				($TaxId,$NameNode,$SubjectName) = split /,/, $TaxIdLog{$GI};
			}
			else
			{
				# 这两步里的grep太慢了，导致整体耗时达到blastn这一步的20倍，后面有时间考虑c重写grep语句代表的逻辑;
				($TaxId,$NameNode) = &TaxIdOfSpe($GI,$Gi2TaxIdDB,$TaxNode);
				$SubjectName = &TaxId2Name($TaxId,$TaxId2SpeDB);
				
				$TaxIdLog{$GI} = $TaxId . "," . $NameNode . "," . $SubjectName;
			}
			
			print TRANS join("\t",$FaName,$EValue,$MatchLen,$TaxId,$SubjectName,$NameNode,$Description),"\n";
		}
		close TRANS;
		close INI;
		printf "[ %s ] Filtered result done.\n",TimeString(time,$BeginTime);
	}
	
	# 数量统计;
	my $StaticResult = $Dir . "/" . $Prefix . ".BlastnStatics.xls";
	if(-e $FltResult)
	{
		my (@Num,@Other) = ();
		my %SpeCount = ();
		open(FLT,"cat $FltResult | grep -v ^# |") or die $!;
		while(my $Line = <FLT>)
		{
			chomp $Line;
			my @Cols = split /\t/, $Line;
			my $SpeId = "-\tUnKnownSpecies";
			$SpeId = $Cols[3] . "\t" . $Cols[4] if($Cols[3] && $Cols[4]);
			$SpeCount{$SpeId} = 0 unless($SpeCount{$SpeId});
			$SpeCount{$SpeId} ++;
		}
		close FLT;
		foreach my $Key (keys %SpeCount)
		{
			push @Num, $SpeCount{$Key};
			push @Other, $Key;
		}
		
		
		my @tRef = NumStringSort(\@Num,\@Other);
		@Num = @{$tRef[0]};
		@Other = @{$tRef[1]};
		open(STAT,"> $StaticResult") or die $!;
		print STAT join("\t","#TaxonID","Species","QueryNumber"),"\n";
		for(my $i = $#Num; $i >= 0; $i --)
		{
			print STAT join("\t",$Other[$i],$Num[$i]),"\n";
		}
		close STAT;
	}
}
printf "[ %s ] Done.\n",TimeString(time,$BeginTime);


######### Sub functions ##########
sub TaxIdOfSpe
{
	my ($GI,$Gi2TaxIdDB,$TaxNode) = @_;
	
	my $TaxId = `grep ^'$GI'\$'\\t' $Gi2TaxIdDB | cut -f 2 | head -n1`;
	chomp $TaxId;
	# 有些gi没有对应的taxid;
	return ("","") unless($TaxId);
	my $NameNode = `grep ^'$TaxId'\$'\\t' $TaxNode | cut -f 5 | head -n1`;
	chomp $NameNode;
	
	my $tTaxId = $TaxId;
	my $DefaultNode = "genus";
	my $LoopLimit = 20;
	while($LoopLimit > 0 && $NameNode ne $DefaultNode)
	{
		my $Return = `grep ^'$tTaxId'\$'\\t' $TaxNode | cut -f 3,5 | head -n1`;
		chomp $Return;
		last unless($Return);
		
		my @Cols = split /\t/, $Return;
		
		if($Cols[1] eq $DefaultNode)
		{
			$TaxId = $tTaxId;
			$NameNode = $Cols[1];
			last;
		}
		else
		{
			$tTaxId = $Cols[0];
		}
		$LoopLimit --;
	}
	
	return ($TaxId,$NameNode);
}

sub TaxId2Name
{
	my ($TaxId,$File) = @_;
	return "" unless($TaxId);
	
	my $Return = `grep ^'$TaxId'\$'\\t' $File | cut -f 3,7`;
	my @Array = split /\n/, $Return;
	my %NameH = ();
	for my $i (0 .. $#Array)
	{
		my @Cols = split /\t/, $Array[$i];
		$NameH{$Cols[1]} = $Cols[0];
	}
	
	my $Name = "";
	$Name = $NameH{"scientific name"} if(! $Name && $NameH{"scientific name"});
	$Name = $NameH{"common name"} if(! $Name && $NameH{"common name"});
	unless($Name)
	{
		my @Cols = split /\t/, $Array[0];
		$Name = $Cols[1];
	}
	
	return $Name;
}