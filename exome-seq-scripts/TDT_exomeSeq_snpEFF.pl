#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

## command line args ##
# mut_type,
# SNV_only,
# trio_dp,
# trio_GQ,
#


my ($file_in, 
	$comb,
	$fam_in,
	$dir_out,
	$file_out_name,
	$SNP_callRate,
	$trio_dp,
	$BG_AC,
	$BG_AF,
	$dbSNP, 
	$kg, 
	$func, 
	$EVS_Filter,	
	$EVS_EA_AC,	
	$EVS_AA_AC,	
	$EVS_TAC,	
	$EVS_DP,
	$EVS_GERP,
	$EVS_AA,	
	$EVS_AAC,	
	$EVS_PH,	
	$EVS_CG,	
	$EVS_GS,
	$Ignore_X,
	$BAF,
	$trio_GQ,
	$SNV_only,
	$Ignore_compHets,
	$mut_type,
	$EVS_only,
	$compHet,
	$freq,
	$BG_freq_new,
	$genoType_field,
	$gene_filter,
	$CNV_in,
	$ignore_multiAllele,
	$site_filter);
	
GetOptions ("file_in=s"=>\$file_in,
			"dir_out=s"=>\$dir_out,
			"file_out_name=s"=>\$file_out_name,
			"CNV_in=s"=>\$CNV_in,
			"fam_in=s"=>\$fam_in,
			"SNV_only=s"=>\$SNV_only,
			"comb=s"=>\$comb,
			"dbSNP=s"=>\$dbSNP,
			"mut_type=s"=>\$mut_type,
			"EVS_only=s"=>\$EVS_only,
			"BG_AF=s"=>\$BG_AF,
			"compHet=s"=>\$compHet,
			"EVS_Filter=s"=>\$EVS_Filter,	
			"EVS_GERP=s"=>\$EVS_GERP,
			"freq=s"=>\$freq,
			"EVS_AA_AC=s"=>\$EVS_AA_AC,	
			"EVS_TAC=s"=>\$EVS_TAC,	
			"EVS_DP=s"=>\$EVS_DP,	
			"EVS_AA=s"=>\$EVS_AA,	
			"BG_freq_new=s"=>\$BG_freq_new,	
			"EVS_AAC=s"=>\$EVS_AAC,	
			"EVS_PH=s"=>\$EVS_PH,	
			"EVS_CG=s"=>\$EVS_CG,	
			"ignore_multiAllele=s"=>\$ignore_multiAllele,
			"Ignore_X=s"=>\$Ignore_X,
			"Ignore_compHets=s"=>\$Ignore_compHets,
			"trio_dp=s"=>\$trio_dp,
			"site_filter=s"=>\$site_filter,
			"kg=s"=>\$kg,
			"func=s"=>\$func,
			"trio_GQ=s"=>\$trio_GQ,
			"SNP_callRate=s"=>\$SNP_callRate,
			"gene_filter=s"=>\$gene_filter,
			"BAF=s"=>\$BAF);

$dir_out .='/' unless ($dir_out =~/.*\/$/);

open (FILE1, "$file_in") || die "can't open annotated VCF file";
open (FAMFILE, "$fam_in") || die "can't open fam file";
open (FILE2, "$comb") || die "can't open genotype combinations file";

open (WRITE_GENE, ">$dir_out$file_out_name\_geneHits.txt") || die "can't write geneHits file\n$dir_out$file_out_name\_geneHits.txt\n";
open (WRITE_TDT, ">$dir_out$file_out_name\_TDT_Stats.txt") || die "can't write TDT_Stats file";
open (WRITE_T, ">$dir_out$file_out_name\_Transmitted_hits.txt") || die "can't write Transmitted file";
open (WRITE_NT, ">$dir_out$file_out_name\_NonTransmitted_hits.txt") || die "can't write Non_Transmitted file";
open (WRITE_DN, ">$dir_out$file_out_name\_deNovos.txt") || die "can't write DeNovo file";
open (WRITE_ME, ">$dir_out$file_out_name\_MendelErrors.txt") || die "can't write Menel file";

my (@Comb, @famFile, @fam2, @gene, @header);

#$Ignore_X = "True"
#	unless defined $Ignore_X;
$BAF = "True"
	unless defined $BAF;

my $chr_pos=0;
my $pos_pos=1;
my $ref_pos=3;
my $alt_pos=4;
my $info_pos=7;
my $filter_pos=6;


while (<FAMFILE>) { 
	my $line = $_;
	chomp $line;
	$line=~s/\r//g;
	my @temp = split(/\t/,$line);
	next if ($temp[5] != 2); ## only look at cases
	push (@famFile,[@temp]);
	}
close (FAMFILE);

##genotype matrix##
while (<FILE2>) {
	my $line = $_;
	chomp $line;
	$line=~s/\r//g;
	my @temp = split(/\t/,$line);
	push (@Comb,[@temp]);
	}
close (FILE2);
###################

my @EffAnoHeader;
## create headers from VCF file ##
open (FILE1, "$file_in") || die "can't open annotated VCF file";
	while (<FILE1>) {
		my $line = $_;
		chomp $line;
		$line=~s/\r//g;
		my @temp = split(/\t/,$line);
		if ($line =~m/^##INFO=<ID=ANN/) {
			my ($tempEffAnoHeader) = $line =~ m/ \'(.*)\' /;
			$tempEffAnoHeader =~s/\s*//g;
			@EffAnoHeader=split(/\|/,$tempEffAnoHeader);
			
			}
		if ($temp[0] eq "#CHROM"){
			push (@header,@temp);
			last;
			}
		}
close (FILE1);
###################################

## print headers ##
print WRITE_T "Chr\tPosition\tRef\tAlt\tFilter\tProband_ID\tP_Genotype\tFather_ID\tF_Genotype\tMother_ID\tM_Genotype\t",join("\t",@EffAnoHeader),"\n";
print WRITE_NT "Chr\tPosition\tRef\tAlt\tFilter\tProband_ID\tP_Genotype\tFather_ID\tF_Genotype\tMother_ID\tM_Genotype\t",join("\t",@EffAnoHeader),"\n";
print WRITE_DN "Chr\tPosition\tRef\tAlt\tFilter\tAlleleCount\tProband_ID\tP_Genotype\tFather_ID\tF_Genotype\tMother_ID\tM_Genotype\tP_altAlleles\tF_altAlleles\tM_altAlleles\t",join("\t",@EffAnoHeader),"\n";
#print WRITE_GENE "Proband_ID\tMutation\tGene\tChr\tPosition\tP_Genotype\tM_Genotype\tF_Genotype\tPlinkSeq_annotation\tPlinkSeq_annotation\tRef\tAlt\tProband_ID\tMutation\tGene\tChr\tPosition\tP_Genotype\tM_Genotype\tF_Genotype\tPlinkSeq_annotation\tPlinkSeq_annotation\tRef\tAlt\n";


## Count total number of loci, total transmitted/non-transmitted variants, mendelian errors
my $varCount=0;
my $T_count=0;
my $NT_count=0;
my $mend_error=0;
my $De_novo_count=0;
my $Trio_count=1;
my $Trio_N=($#famFile+1);
my $Trios_analysed=0;
my %parent_seen=();

my (@DN, @MD, @cnvSNP,@Trios_TDT);

##send each family to TDT sub-routine###
my %seen=();

foreach my $x (@famFile) {
	my $P_Pos;
	my $F_Pos;
	my $M_Pos;
	my $count = 0;
	my $trio_count_test=0;

	foreach my $y (@header) { ## get positions of trio members in VCF file
		if ($x->[1] eq $y) { ## proband position
			$P_Pos = $count;
			$trio_count_test++
			}
		if ($x->[2] eq $y) { ## father position
			$F_Pos = $count;
			$trio_count_test++
			}
		if ($x->[3] eq $y) { ## mother position
			$M_Pos = $count;
			$trio_count_test++
			}
		$count++;
		}

	if ($trio_count_test != 3) {
		print WRITE_TDT "$x->[0]-$x->[1] incomplete Trio\n\n";
		}
	else {
		print "Working on Trio $Trio_count out of $Trio_N\n";
		TDT ($P_Pos,$F_Pos,$M_Pos,$x->[0],$x->[1],$x->[2],$x->[3]); ## send trio to TDT subroutine
		$Trios_analysed++;
		}
	$Trio_count++;	 	
	}
########################################



##Print TDT Statistics##
$varCount = $varCount/$Trio_N; ##Average total number of variants for each trio##
$varCount=sprintf("%.0f",$varCount);

print WRITE_TDT "$Trios_analysed Parent-Offspring Trios analysed.\n\n";
print WRITE_TDT "Parameters Used:\n\nTrio DP = >$trio_dp\n";
#if ($segDup) {print WRITE_TDT "Segmental Duplications Filtered\n";}
#if ($dbSNP) {print WRITE_TDT "Exclude variants found in dbSNP\n";}
#if ($EVS_EA_AC) {print WRITE_TDT "Exclude variants with a freq in EVS Europeans above $EVS_EA_AC\n";}
#if ($EVS_AA_AC) {print WRITE_TDT "Exclude variants with a freq in EVS Africans above $EVS_AA_AC\n";}
#if ($EVS_TAC) {print WRITE_TDT "Exclude variants with a freq in EVS EA + AA above $EVS_TAC\n";}
#if ($kg) {print WRITE_TDT "Exclude variants with a freq in Thousand Genomes above $kg\n";}
#if ($BG_AF) {print WRITE_TDT "Exclude variants with a freq in Bulgarian Population above $BG_AF\n";}
#if ($exonicFunc) {print WRITE_TDT "Exclude all variants not predicted to be $exonicFunc\n";}
#if ($EVS_PH) {print WRITE_TDT "Exclude variants predicted by EVS_polyPhen to be $EVS_PH\n";}
#if ($Ignore_X) {print WRITE_TDT "Exclude variants on X Chromosome\n"};
#if ($SNP_callRate) {print WRITE_TDT "Exlude variants loci with a call rate < $SNP_callRate\n";}
#if ($SNV_only) {print WRITE_TDT "Inlcude only SNPs";}

print WRITE_TDT "\nSummary Statistics\n\n";
print WRITE_TDT "Total number of SNP loci counted = $varCount\n";
print WRITE_TDT "Total number of Transmitted variants = $T_count\nTotal number of Non-Transmitted variants  = $NT_count\nTotal number of mendelian errors = $mend_error\nTotal number of potential De Novos = $De_novo_count\n";

print WRITE_TDT "\nPotential De Novos\n\n";
foreach my $x (@DN) {
	print WRITE_TDT join("\t",@$x),"\n";
	}
	
print WRITE_TDT "\nMedelian Errors\n\n";
foreach my $x (@MD) {
	print WRITE_TDT join("\t",@$x),"\n";
	}
print WRITE_TDT "\n";

if (@cnvSNP) {
	print WRITE_CNV "\nCNV/SNP Overlap\n\n";
	print WRITE_CNV "Type\tVariant\tGene\tChr\tPos\tPGT\tMGT\tFGT\tVariant\tID\tChr\tBP1\tBP2\tType\n";
	foreach my $x (@cnvSNP) {
			print WRITE_CNV join("\t",@{$x}),"\n";
		}
	}

print WRITE_TDT "Trio\tTrio\tTransmissions\tNon-Transmissions\tMendel_Errors\tDe_novos\n";
foreach my $x (@Trios_TDT) {
	print WRITE_TDT join("\t",@{$x}),"\n";
	}
#############################################


##Sub-routines##

sub TDT {
	my ($P_POS, 
		$F_POS, 
		$M_POS,
		$TRIO,
		$FFI1,
		$FFI2,
		$FFI3) = @_;
		
	open (FILE1, "$file_in") || die "cant open variant file";

	my (@file_1, @CNV_vars);
	my $temp_T_count=0;
	my $temp_NT_count=0;
	my $temp_mend_error=0;
	my $temp_De_novo_count=0;

	while (<FILE1>) {
		my $line = $_;
		chomp $line;
		$line=~s/\r//g;
		$line=~s/"//g;
		next if ($line =~m/^#/);
		my @temp = split(/\t/,$line);
			
		
		
	if ($Ignore_X) {$temp[0]=~s/chr//g; if ($temp[0] eq "X") {next;}} ##Ignore X Chromosome
	if ($site_filter) {
		next if ($temp[$filter_pos] ne $site_filter);
		}
		
	##push trio genetype info into arrays
	my @p=split(":",$temp[$P_POS]); 
	my @m=split(":",$temp[$M_POS]);
	my @f=split(":",$temp[$F_POS]);
	my $temp_mCount=0;
	#####################################
	my @temp_p_altAlleles=split(",",$p[1]);
	my @temp_f_altAlleles=split(",",$f[1]);
	my @temp_m_altAlleles=split(",",$m[1]);
	my $p_altAlleles=$temp_p_altAlleles[1];
	my $f_altAlleles=$temp_f_altAlleles[1];
	my $m_altAlleles=$temp_m_altAlleles[1];
	
	## if data is missing for an individual skip variant##
		next if (($p[0] eq "./.") || ($m[0] eq "./.") || ($f[0] eq "./."));		
	##exclude if individual has a depth less than..##
		next if (($p[2] eq ".") || ($f[2] eq ".") || ($m[2] eq ".")); ## Think about if this is correct!!
		if ($trio_dp) {
			next if (($p[2] < $trio_dp) || ($f[2] < $trio_dp) || ($m[2] < $trio_dp));
			}
		
	##only look at positions in which a member of the trio has a variant##
		next if (($p[0] eq "0/0") && ($m[0] eq "0/0") && ($f[0] eq "0/0"));
		##exclude if individual has a GQ less than ..##
		if ($trio_GQ) {
			next if (($p[3] < $trio_GQ) || ($f[3] < $trio_GQ) || ($m[3] < $trio_GQ));
			}
		
	## exclude variants that fail BAF ##
	my @pGenoTemp=split(/\//,$p[0]);
		my @mGenoTemp=split(/\//,$m[0]);
		my @fGenoTemp=split(/\//,$f[0]);
		my @allGenoTemp= (@pGenoTemp,@mGenoTemp,@fGenoTemp);
		my @allGenoUnique=uniq(@allGenoTemp);
	if ($BAF) {
			
		
			################################		
			if ($pGenoTemp[0] == 0 && $pGenoTemp[1] == 0) {
				my $p_BAF = cal_BAF2($p[1],$p[0]);
				next if ($p_BAF <=0.9); 
				}
			elsif ($pGenoTemp[0] != $pGenoTemp[1]) {
				my $p_BAF = cal_BAF2($p[1],$p[0]);
				next if ($p_BAF <= 0.2 || $p_BAF >=0.8);
				}
			elsif (($pGenoTemp[0] != 0 && $pGenoTemp[1] != 0) && ($pGenoTemp[0] == $pGenoTemp[1])) {
				my $p_BAF = cal_BAF2($p[1],$p[0]);
				next if ($p_BAF <= 0.90);
				}			
			else {
				print "BAF error\n$pGenoTemp[0]\t$pGenoTemp[1]\nLine $.\n"; die;
				}
				
			if ($fGenoTemp[0] == 0 && $fGenoTemp[1] == 0) {
				my $f_BAF = cal_BAF2($f[1],$f[0]);
				next if ($f_BAF <=0.9); 
				}
			elsif ($fGenoTemp[0] != $fGenoTemp[1]) {
				my $f_BAF = cal_BAF2($f[1],$f[0]);
				next if ($f_BAF <= 0.2 || $f_BAF >=0.8);
				}
			elsif (($fGenoTemp[0] != 0 && $fGenoTemp[1] != 0) && ($fGenoTemp[0] == $fGenoTemp[1])) {
				my $f_BAF = cal_BAF2($f[1],$f[0]);
				next if ($f_BAF <= 0.90);
				}			
			else {
				print "BAF error\n$fGenoTemp[0]\t$fGenoTemp[1]\nLine $.\n"; die;
				}
				
			if ($mGenoTemp[0] == 0 && $mGenoTemp[1] == 0) {
				my $m_BAF = cal_BAF2($m[1],$m[0]);
				next if ($m_BAF <=0.9); 
				}
			elsif ($mGenoTemp[0] != $mGenoTemp[1]) {
				my $m_BAF = cal_BAF2($m[1],$m[0]);
				next if ($m_BAF <= 0.2 || $m_BAF >=0.8);
				}
			elsif (($mGenoTemp[0] != 0 && $mGenoTemp[1] != 0) && ($mGenoTemp[0] == $mGenoTemp[1])) {
				my $m_BAF = cal_BAF2($m[1],$m[0]);
				next if ($m_BAF <= 0.90);
				}			
			else {
				print "BAF error\n$mGenoTemp[0]\t$mGenoTemp[1]\nLine $.\n"; die;
				}
		}
	#############################################
	
	if ($#allGenoUnique > 2) {
#	print $.,"\t",$#allGenoUnique,"\t",$TRIO,,"\t@allGenoUnique","\n";		
	}
			
	### bi or multi allelic site ###	
	if ($temp[$alt_pos]=~m/,/g) { ## multi allelic
		my @multiVars=split(",",$temp[$alt_pos]);
		my $Nvars=$#multiVars+1;

		}

		
	else { ## bi-allelic
	###################################
	
	
	my %annotation = ();
	### get snpEFF annotations ###
	my @Info=split(";",$temp[$info_pos]);
	my @firstANN;
	my %splitInfo=();
	foreach my $x (@Info) {
		
		my @tempSplit=split("=",$x);
		$splitInfo{$tempSplit[0]} = $tempSplit[1];
		
		if ($x=~m/^ANN=/) {
			$x=~s/ANN=//g;
			my @ANN=split(",",$x);
			@firstANN=split(/\|/,$ANN[0]);
			@annotation {@EffAnoHeader} = @firstANN;					
			if ($temp[4] ne $annotation{"Allele"}) {
				die;
				}
			
			}
		}
	
	next if ($annotation{"Annotation_Impact"} eq "MODIFIER");
	
	##Determine Transmission Pattern##
		foreach my $c (@Comb) {
			if (($c->[0] eq $m[0]) and ($c->[1] eq $p[0]) and ($c->[2] eq $f[0])) {	##no variant found in trio (should not happen due to prior exclusion criteria)	
				if ($c->[3] eq "X") { 
					$temp_mCount++; ##temp_mCount > 0 means no Mendelian error
					next;
					}
				elsif ($c->[3] eq "T") { ##genetypes match that of a transmission
					$temp_mCount++;
					$temp_T_count+=$c->[4];
					$T_count+=$c->[4]; ## add to transmission 1 or 2 counts based on data
					T_forloop ($c->[4],$temp[$chr_pos],$temp[$pos_pos],$temp[$ref_pos],$temp[$alt_pos],$temp[$filter_pos],$FFI1,$FFI2,$FFI3,$temp[$P_POS],$temp[$F_POS],$temp[$M_POS],\@firstANN);##print transmitted variant N times
					}
				elsif ($c->[3] eq "NT") { ##genetypes match that of a non-transmission
					$temp_mCount++;
					$temp_NT_count+=$c->[4];
					$NT_count+=$c->[4];
					NT_forloop ($c->[4],$temp[$chr_pos],$temp[$pos_pos],$temp[$ref_pos],$temp[$alt_pos],$temp[$filter_pos],$FFI1,$FFI2,$FFI3,$temp[$P_POS],$temp[$F_POS],$temp[$M_POS],\@firstANN);##print non-transmitted variant N times
					}
				elsif ($c->[3] eq "DN")	{ ##genotypes match that of a de novo
							$De_novo_count++;
							$temp_De_novo_count++;
							$temp_mCount++;
							print WRITE_DN $temp[$chr_pos],"\t",$temp[$pos_pos],"\t",$temp[$ref_pos],"\t",$temp[$alt_pos],"\t",$temp[$filter_pos],"\t",$splitInfo{"AC"},"\t",$FFI1,"\t",$temp[$P_POS],"\t",$FFI2,"\t",$temp[$F_POS],"\t",$FFI3,"\t",$temp[$M_POS],"\t",$p_altAlleles,"\t",$f_altAlleles,"\t",$m_altAlleles,"\t",join("\t",@firstANN),"\n";
							#push (@DN,[$FFI1,$FFI2,$FFI3,$temp[$exonicFunc],$temp[$gene],$temp[$exonic],$temp[2],$temp[$chr],$temp[$start],$temp[$ref],$temp[$alt],$temp[$P_POS],$temp[$M_POS],$temp[$F_POS],$temp[59]]);					
						}
					}
				}
			if ($temp_mCount==0) { ##genotypes did not match any combination, mark as Medelian error
				print WRITE_ME $temp[$chr_pos],"\t",$temp[$pos_pos],"\t",$temp[$ref_pos],"\t",$temp[$alt_pos],"\t",$temp[$filter_pos],"\t",$splitInfo{"AC"},"\t",$FFI1,"\t",$temp[$P_POS],"\t",$FFI2,"\t",$temp[$F_POS],"\t",$FFI3,"\t",$temp[$M_POS],"\t",join("\t",@firstANN),"\n";
				$mend_error++;
				$temp_mend_error++;							
			}
	 ###############################################		
		
		
		}
	
	
	################################

	}
	
	close (FILE1);

	push (@Trios_TDT,[$TRIO,$FFI1,$temp_T_count,$temp_NT_count,$temp_mend_error,$temp_De_novo_count]);##push Trio TDT stats into array
}


sub cal_BAF {
	my ($AD) = @_;
	my @AD2=split(",",$AD);
	if ($#AD2 != 1) {print "AD more than 2 numbers\n"; die;}
		my $TAF=$AD2[0] + $AD2[1];
		my $BAF=$AD2[1]/$TAF;
		return $BAF;		
	}

sub cal_BAF2 {
	my ($AD,$geno) = @_;
	my @geno2=split(/\//,$geno); ## $geno2[1] is code for alt allele	
	my @AD2=split(",",$AD);
	my $TAF=0;
	foreach my $x (@AD2) {
		$TAF=$TAF + $x;
		}
	if ($TAF==0) {
		my $BAF=0;
		return $BAF;
		}
	else {		
		my $BAF=$AD2[$geno2[1]]/$TAF;
		return $BAF;	
		}	
	}

sub uniq {
    my %seen;
    grep !$seen{$_}++, @_;
}

sub T_forloop {  
	my ($N,$chr,$pos,$ref,$alt,$filter,$p,$f,$m,$pINFO,$fINFO,$mINFO,$annotations) = @_;
	for (my $i=1;$i<=$N;$i++) {
		print WRITE_T "$chr\t$pos\t$ref\t$alt\t$filter\t$p\t$pINFO\t$f\t$fINFO\t$m\t$mINFO\t",join("\t",@{$annotations}),"\n";
		}
	}
	
sub NT_forloop { 
	my ($N,$chr,$pos,$ref,$alt,$filter,$p,$f,$m,$pINFO,$fINFO,$mINFO,$annotations) = @_;
	for (my $i=1;$i<=$N;$i++) {
		print WRITE_NT "$chr\t$pos\t$ref\t$alt\t$filter\t$p\t$pINFO\t$f\t$fINFO\t$m\t$mINFO\t",join("\t",@{$annotations}),"\n";
		}
	}	
	
	
close (WRITE_GENE);
close (WRITE_TDT);
close (WRITE_T);
close (WRITE_NT);
close (WRITE_DN);
close (WRITE_ME);




=pod
#####conditions to exclude variants#####
	if ($mut_type) {##exclude variants if not S/NS
		if ($mut_type eq "nonsynonymous") {
			if ($temp[$exonicFunc] eq "silent" ||
				$temp[$exonicFunc] eq "synonymous SNV" || 
				$temp[$exonic] eq "intronic" ||
				$temp[$exonicFunc] eq "NA" ||
				$temp[$exonicFunc] eq "intronic" ||
				$temp[$exonicFunc] eq "UTR5" ||
				$temp[$exonic] eq "intergenic" ||
				$temp[$exonicFunc] eq "ncRNA_intronic" ||
				$temp[$exonicFunc] eq "unknown" ||
				$temp[$exonicFunc] eq "UTR3") {
					next;
					}
			}
		elsif ($mut_type eq "synonymous") {
			if ($temp[$exonicFunc] ne "synonymous SNV") {
				next;
				}
			}
		}


			

	
	next if ($temp[$filter] ne "PASS");
	
	########################
	## get allele count, frequency and number ONLY WORKS FOR SINGLE ALLELE##
	my @split_info=split(";",$temp[$info]);
	my @alleleFreqTemp=split("=",$split_info[1]);
	my $alleleFreq=$alleleFreqTemp[1];
	my @alleleCountTemp=split("=",$split_info[0]);
	my $alleleCount=$alleleCountTemp[1];
	my @alleleNumberTemp=split("=",$split_info[2]);
	my $alleleNumber=$alleleNumberTemp[1];
	my $alleleFreqTest=$alleleCount/$alleleNumber;

	if ($freq) {
		if ($temp[$evs_all] ne ".") {
			if ($temp[$evs_all] > $freq) {
		#		next;
				}
			}
		if ($temp[$kg_all] ne ".") {
			if ($temp[$kg_all] > $freq) {
		#		next;
				}
			}
		if ($alleleFreqTest > $freq) {
		#	next;
			}
		if ($alleleCount > 1) {
			next;
			}
		}
		


		$varCount++;
		my (@p, @m, @f);
 
	if ($temp[$format] ne "GT:AD:DP:GQ:PL" && $temp[$format] ne "GT:AD:DP:GQ:PGT:PID:PL") { ##Retrieve genotype fields, will differ for different datasets. 
		print "$temp[$format] Genotype field format does not match 'GT:AD:DP:GQ:PL' or 'GT:AD:DP:GQ:PGT:PID:PL'"; 
		next;
		}
		
		$temp[0]=~s/chr//g;
		@p=split(":",$temp[$P_POS]); ##push trio genetype info into arrays
		@m=split(":",$temp[$M_POS]);
		@f=split(":",$temp[$F_POS]);
		my $temp_mCount=0;
		## if data is missing for an individual skip variant##
		if (($p[0] eq "./.") || ($m[0] eq "./.") || ($f[0] eq "./.")) {
			next;
			}		
		##exclude if individual has a depth less than..##
		if ($trio_dp) {
			if (($p[2] < $trio_dp) || ($f[2] < $trio_dp) || ($m[2] < $trio_dp)) {
					next;
					}
				}	
		##only look at positions in which a member of the trio has a variant##
			if (($p[0] eq "0/0") && ($m[0] eq "0/0") && ($f[0] eq "0/0")) {
				next;
				}
		##exclude if individual has a GQ less than ..##
		if ($trio_GQ) {
			if (($p[3] < $trio_GQ) || ($f[3] < $trio_GQ) || ($m[3] < $trio_GQ)) {
					#print "@p\t@f\t@m\n";
					next;
					}
				}	
		##remove variants that fail BAF##
		if ($BAF) {
			if ($p[0] eq "0/1") {
				my $p_BAF = cal_BAF($p[1]);
				if ($p_BAF <= 0.2 || $p_BAF >=0.8) {
					#print "B allele fail $p_BAF\t$FFI1\t$temp[15]\t$temp[16]\n";
					next;
					}
				}
			if ($m[0] eq "0/1") {
				my $m_BAF = cal_BAF($m[1]);
				if ($m_BAF <= 0.2 || $m_BAF >=0.8) {
					#print "B allele fail $m_BAF\t$FFI1\t$temp[15]\t$temp[16]\n";
					next;
					}
				}
			if ($f[0] eq "0/1") {
				my $f_BAF = cal_BAF($f[1]);
				if ($f_BAF <= 0.2 || $f_BAF >=0.8) {
					#print "B allele fail $f_BAF\t$FFI1\t$temp[15]\t$temp[16]\n";
					next;
					}
				}
			}
		################################		
		if ($BAF) {
			if ($p[0] eq "0/0") {
				my $p_BAF = cal_BAF($p[1]);
				if ($p_BAF >=0.1) {
					#print "B allele fail $p_BAF\t$FFI1\t$temp[15]\t$temp[16]\n";
					next;
					}
				}
			if ($m[0] eq "0/0") {
				my $m_BAF = cal_BAF($m[1]);
				if ($m_BAF >=0.1) {
					#print "B allele fail $m_BAF\t$FFI1\t$temp[15]\t$temp[16]\n";
					next;
					}
				}
			if ($f[0] eq "0/0") {
				my $f_BAF = cal_BAF($f[1]);
				if ($f_BAF >=0.1) {
					#print "B allele fail $f_BAF\t$FFI1\t$temp[15]\t$temp[16]\n";
					next;
					}
				}
			}
		################################	
		if ($BAF) {
			if ($p[0] eq "1/1") {
				my $p_BAF = cal_BAF($p[1]);
				if ($p_BAF <= 0.90) {
					#print "B allele fail $p_BAF\t$FFI1\t$temp[15]\t$temp[16]\n";
					next;
					}
				}
			if ($m[0] eq "1/1") {
				my $m_BAF = cal_BAF($m[1]);
				if ($m_BAF <= 0.90) {
					#print "B allele fail $m_BAF\t$FFI1\t$temp[15]\t$temp[16]\n";
					next;
					}
				}
			if ($f[0] eq "1/1") {
				my $f_BAF = cal_BAF($f[1]);
				if ($f_BAF <= 0.90) {
					#print "B allele fail $f_BAF\t$FFI1\t$temp[15]\t$temp[16]\n";
					next;
					}
				}
			}
		################################	
				
		##push all variants that pass filtering into array for compound hets##
		push (@file_1, [$temp[$chr], $temp[$start], $temp[$P_POS], $temp[$M_POS], $temp[$F_POS], $temp[$gene],$temp[$exonicFunc],$temp[$exonic],$temp[$ref],$temp[$alt]]); 
	#	push (@file_1, [$temp[21], $temp[22], $temp[$P_POS], $temp[$M_POS], $temp[$F_POS], $temp[1],$temp[0],$temp[2],$temp[24],$temp[25]]);

			
		##Determine Transmission Pattern##
		foreach my $c (@Comb) {
			if (($c->[0] eq $m[0]) and ($c->[1] eq $p[0]) and ($c->[2] eq $f[0])) {	##no variant found in trio (should not happen due to prior exclusion criteria)	
				if ($c->[3] eq "X") { 
					$temp_mCount++; ##temp_mCount > 0 means no Mendelian error
					next;
					}
				elsif ($c->[3] eq "T") { ##genetypes match that of a transmission
					$temp_mCount++;
					$temp_T_count+=$c->[4];
					$T_count+=$c->[4]; ## add to transmission 1 or 2 counts based on data
					T_forloop ($c->[4],$temp[$chr],$temp[$start],$temp[$gene],$FFI1,$FFI2,$FFI3,$temp[$P_POS],$temp[$F_POS],$temp[$M_POS],$temp[$ref],$temp[$alt],$temp[$exonicFunc],$temp[$exonic],$alleleFreqTest);##print transmitted variant N times
					}
				elsif ($c->[3] eq "NT") { ##genetypes match that of a non-transmission
					$temp_mCount++;
					$temp_NT_count+=$c->[4];
					$NT_count+=$c->[4];
					NT_forloop ($c->[4],$temp[$chr],$temp[$start],$temp[$gene],$FFI1,$FFI2,$FFI3,$temp[$P_POS],$temp[$F_POS],$temp[$M_POS],$temp[$ref],$temp[$alt],$temp[$exonicFunc],$temp[$exonic],$alleleFreqTest);##print non-transmitted variant N times
					}
				elsif ($c->[3] eq "DN")	{ ##genotypes match that of a de novo
					if ($trio_dp) {
						if (($p[2] >= $trio_dp) && ($f[2] >= $trio_dp) && ($m[2] >= $trio_dp)) {
							$De_novo_count++;
							$temp_De_novo_count++;
							$temp_mCount++;
							push (@DN,[$FFI1,$FFI2,$FFI3,$temp[$exonicFunc],$temp[$gene],$temp[$exonic],$temp[2],$temp[$chr],$temp[$start],$temp[$ref],$temp[$alt],$temp[$P_POS],$temp[$M_POS],$temp[$F_POS],$temp[59]]);
							}
						}
					else {
						$De_novo_count++;
						$temp_De_novo_count++;
						$temp_mCount++;
						push (@DN,[$FFI1,$FFI2,$FFI3,$temp[$exonicFunc],$temp[$gene],$temp[$exonic],$temp[2],$temp[$chr],$temp[$start],$temp[$ref],$temp[$alt],$temp[$P_POS],$temp[$M_POS],$temp[$F_POS]]);
						}
					}
				}
			}
		if ($temp_mCount==0) { ##genotypes did not match any combination, mark as Medelian error
				push (@MD,[$FFI1,$FFI2,$FFI3,$temp[$exonicFunc],$temp[$gene],$temp[$exonic],$temp[2],$temp[$chr],$temp[$start],$temp[$ref],$temp[$alt],$temp[$P_POS],$temp[$M_POS],$temp[$F_POS]]);
				$mend_error++;
				$temp_mend_error++;							
			}
		
		}
	
	
close (FILE1);

push (@Trios_TDT,[$TRIO,$FFI1,$temp_T_count,$temp_NT_count,$temp_mend_error,$temp_De_novo_count]);##push Trio TDT stats into array
				
if ($Ignore_compHets) {
	##Bypass compHet and CNV hits subroutine##
	}
else {
##pass variants that survived filtering to compHet subroutine##
#compHet ($FFI1, $#file_1, \@file_1);
	}
	

}
###end sub TDT#####

##compound hets##
sub compHet {
my ($trio_ID,$N_Var, $file_2) = @_;
print "Analysing $trio_ID\n";
print "number of variants = $N_Var\n";
my@FAM_ID=$trio_ID;
	
	my $DN_compHet_count=0;
	my $compHet_count=0;
	my $recess_count=0;
	my %skipGene=();
	my %genes_inFam =();
	#my $fileCount=0;
	
	
	## find genes hit in trio ##
	foreach my $y (@{$file_2}) { 
		next if exists $skipGene {$y->[5]};		
		if ($y->[5] =~m/,/g) { ## Account for multiple genes hit by variant
			#$fileCount++;
			my @tempSplit = split(",",$y->[5]);
			$skipGene {$y->[5]} = 1;
			foreach my $x (@tempSplit) {
				$genes_inFam {$x} = 1;
				}
			}
		elsif ($y->[5] =~m/;/g) { ## this handles variants in different genes - not sure how this is coded in plinkseq
		#	print "$y->[5]\n";
			$skipGene {$y->[5]} =1;
			#$fileCount++;
			my @tempSplit = split(";",$y->[5]);
			if ($tempSplit[0] ne $tempSplit[1]) {
				$genes_inFam {$tempSplit[0]} = 1;
				$genes_inFam {$tempSplit[1]} = 1;
				}
			else {
				$genes_inFam {$tempSplit[0]} = 1;
				}
			}
		else {
			$skipGene {$y->[5]} = 1;
			$genes_inFam {$y->[5]} = 1;
			}
		
		}
	##############################
	
	## loop through each gene ##
	foreach my $key (keys %genes_inFam) {

		my $testCount=0;
		my $fileCount2=0;
		my (@geneHits,@compHits,@recess,@deNovo,@controlRecess,@controlCompHet, @CCH_m, @CCH_f, @parentalRecess);
		
		foreach my $y (@{$file_2}) { ## loop through variants in trio
			my $keyCount=0;		
			my @tempSplit;
			@tempSplit = split(",",$y->[5]);
			
			## check for variants hitting target gene ##
			if ($#tempSplit>0) {
				foreach my $tempKey (@tempSplit) {		
					if ($key eq $tempKey) {
						$keyCount++;					
						}				
					}				
				}
			else {
				$keyCount++ if $key eq $y->[5];
				}
		
			my @tempSplit2;
			@tempSplit2 = split(";",$y->[5]);
			if ($#tempSplit2>0) {
				foreach my $tempKey (@tempSplit2) {			
						if ($key eq $tempKey) {
							$keyCount++;					
							}
							
					}				
				}
			
			next if $keyCount==0; ## skip variants not hitting target gene		
			push (@geneHits, [$y->[2],$y->[3],$y->[4],$y->[0],$y->[1],$y->[6],$y->[7],$y->[8],$y->[9]]);
			$testCount++;
			}

		if ($testCount==0) {
			print "$key $trio_ID var not in hash\n";
			die;
			}
			
		## set transmission variables ##
		my $m_T=0;
		my $f_T=0;
		my $d_T=0;
		my $r_T=0;
		my $cm_NT=0;
		my $cf_NT=0;
		my $cr_NT=0;
		foreach my $h (@geneHits) {	
			my @p=split(":",$h->[0]);
			my @m=split(":",$h->[1]);
			my @f=split(":",$h->[2]);
	
				##print all homozygous alleles##
				if ($m[0] eq "1/1" && not exists $parent_seen {"$FAM_ID[0]"}) {
					print WRITE_GENE "$trio_ID\tmother_homoRecess\t$key\t$h->[3]\t$h->[4]\t$h->[0]\t$h->[1]\t$h->[2]\t$h->[5]\t$h->[6]\t$h->[7]\t$h->[8]\n";
					}
				if ($f[0] eq "1/1" && not exists $parent_seen {"$FAM_ID[0]"}) {
					print WRITE_GENE "$trio_ID\tfather_homoRecess\t$key\t$h->[3]\t$h->[4]\t$h->[0]\t$h->[1]\t$h->[2]\t$h->[5]\t$h->[6]\t$h->[7]\t$h->[8]\n";
					}	
				if ($p[0] eq "1/1") {
					print WRITE_GENE "$trio_ID\tproband_homoRecess\t$key\t$h->[3]\t$h->[4]\t$h->[0]\t$h->[1]\t$h->[2]\t$h->[5]\t$h->[6]\t$h->[7]\t$h->[8]\n";
					}
				################################
				
				## Determine transmission pattern for all variants in target gene ##
				foreach my $c (@Comb) {
				if (($c->[0] eq $m[0]) and ($c->[1] eq $p[0]) and ($c->[2] eq $f[0])) {
					if ($c->[3] eq "T") { ## transmissions
						if ($c->[5] eq "M") { ## mother transmits
							$m_T++;
							push (@compHits, [$h->[3], $h->[4], $h->[0],$h->[1],$h->[2],$h->[5],$h->[6],$h->[7],$h->[8]]);
							push (@CCH_m, [$h->[3], $h->[4], $h->[0],$h->[1],$h->[2],$h->[5],$h->[6],$h->[7],$h->[8]]);
							}
						elsif ($c->[5] eq "F") { ## father transmits
							$f_T++;
							push (@compHits, [$h->[3], $h->[4], $h->[0],$h->[1],$h->[2],$h->[5],$h->[6],$h->[7],$h->[8]]);
							push (@CCH_f, [$h->[3], $h->[4], $h->[0],$h->[1],$h->[2],$h->[5],$h->[6],$h->[7],$h->[8]]);
							}
						elsif ($c->[5] eq "R") { ## homozygous
							$r_T++;						
							push (@recess, [$h->[3], $h->[4], $h->[0],$h->[1],$h->[2],$h->[5],$h->[6],$h->[7],$h->[8]]);
							push (@CCH_m, [$h->[3], $h->[4], $h->[0],$h->[1],$h->[2],$h->[5],$h->[6],$h->[7],$h->[8]]);
							push (@CCH_f, [$h->[3], $h->[4], $h->[0],$h->[1],$h->[2],$h->[5],$h->[6],$h->[7],$h->[8]]);
							}
						}
					elsif ($c->[3] eq "DN")	{ ## de novo
						$d_T++;
						push (@deNovo, [$h->[3], $h->[4], $h->[0],$h->[1],$h->[2],$h->[5],$h->[6],$h->[7],$h->[8]]);
						}
					elsif ($c->[3] eq "NT")	{ ## non-transmissions
						if ($c->[5] eq "CR") { ## mother and father non-transmit
							$cr_NT++;
							push (@controlRecess, [$h->[3], $h->[4], $h->[0],$h->[1],$h->[2],$h->[5],$h->[6],$h->[7],$h->[8]]);
							push (@CCH_m, [$h->[3], $h->[4], $h->[0],$h->[1],$h->[2],$h->[5],$h->[6],$h->[7],$h->[8]]);
							push (@CCH_f, [$h->[3], $h->[4], $h->[0],$h->[1],$h->[2],$h->[5],$h->[6],$h->[7],$h->[8]]);
							}
						elsif ($c->[5] eq "M") { ## mother non-transmission
							$cm_NT++;
							push (@controlCompHet, [$h->[3], $h->[4], $h->[0],$h->[1],$h->[2],$h->[5],$h->[6],$h->[7],$h->[8]]);
							push (@CCH_m, [$h->[3], $h->[4], $h->[0],$h->[1],$h->[2],$h->[5],$h->[6],$h->[7],$h->[8]]);
							}
						elsif ($c->[5] eq "F") { ## father non-transmission
							$cf_NT++;
							push (@controlCompHet, [$h->[3], $h->[4], $h->[0],$h->[1],$h->[2],$h->[5],$h->[6],$h->[7],$h->[8]]);
							push (@CCH_f, [$h->[3], $h->[4], $h->[0],$h->[1],$h->[2],$h->[5],$h->[6],$h->[7],$h->[8]]);
							}
						}
					}
				}
				#######################################################################
			}
			
			
	## Identify compound heterozygous mutations ##
	
		if ($m_T >=1 && $f_T >=1) { ## both mother and father transmit - proband has a compound het
			$compHet_count++;
			foreach my $c (@compHits) {
				print WRITE_GENE "$trio_ID\tcompHet\t$key\t$c->[0]\t$c->[1]\t$c->[2]\t$c->[3]\t$c->[4]\t$c->[5]\t$c->[6]\t$c->[7]\t$c->[8]\t";
				}
			print WRITE_GENE "\n";
			}
		
		if (($m_T >=1 && $d_T >=1) || ($f_T >=1 && $d_T >=1)) { ## either mother or father transmits and proband has a de novo - potential proband compound het
			$DN_compHet_count++;
			foreach my $c (@compHits) {
				print WRITE_GENE "$trio_ID\tDN_compHet\t$key\t$c->[0]\t$c->[1]\t$c->[2]\t$c->[3]\t$c->[4]\t$c->[5]\t$c->[6]\t$c->[7]\t$c->[8]\t";
				}
			foreach my $c (@deNovo) {
				print WRITE_GENE "$trio_ID\tDN_compHet\t$key\t$c->[0]\t$c->[1]\t$c->[2]\t$c->[3]\t$c->[4]\t$c->[5]\t$c->[6]\t$c->[7]\t$c->[8]\t";
				}
			print WRITE_GENE "\n";
			}
		
		if ($r_T >=1) { ## print homozygous alleles ##	
			if ($#recess == 0 ){ ##if statment for 0 or multiple gene homo recess hits##
				foreach my $c (@recess) {
					$recess_count++;
					print WRITE_GENE "$trio_ID\thomoRecess\t$key\t$c->[0]\t$c->[1]\t$c->[2]\t$c->[3]\t$c->[4]\t$c->[5]\t$c->[6]\t$c->[7]\t$c->[8]\t";
					}
				print WRITE_GENE "\n";
				}
			else {
				foreach my $c (@recess) {
					$recess_count++;
					print WRITE_GENE "$trio_ID\thomoRecess_multiHit\t$key\t$c->[0]\t$c->[1]\t$c->[2]\t$c->[3]\t$c->[4]\t$c->[5]\t$c->[6]\t$c->[7]\t$c->[8]\t";
					}
				print WRITE_GENE "\n";
				}		
			}
			
		if ($cm_NT >=1 && $cf_NT >=1) { ## control compHet in TDT
			foreach my $c (@controlCompHet) {
				print WRITE_GENE "$trio_ID\tcontrol_compHet\t$key\t$c->[0]\t$c->[1]\t$c->[2]\t$c->[3]\t$c->[4]\t$c->[5]\t$c->[6]\t$c->[7]\t$c->[8]\t";
				}
			print WRITE_GENE "\n";
			}
		
		if ($cr_NT >=1) { ## control homozygous allele
			if ($#controlRecess == 0 ){ ##if statment for 0 or multiple gene homo recess hits##
				foreach my $c (@controlRecess) {
					print WRITE_GENE "$trio_ID\tcontrol_homoRecess\t$key\t$c->[0]\t$c->[1]\t$c->[2]\t$c->[3]\t$c->[4]\t$c->[5]\t$c->[6]\t$c->[7]\t$c->[8]\t";
					}
				print WRITE_GENE "\n";
				}
			else {
				foreach my $c (@controlRecess) {
					print WRITE_GENE "$trio_ID\tcontrol_homoRecess_multiHit\t$key\t$c->[0]\t$c->[1]\t$c->[2]\t$c->[3]\t$c->[4]\t$c->[5]\t$c->[6]\t$c->[7]\t$c->[8]\t";
					}
				print WRITE_GENE "\n";
				}		
			}
	
		## test mother compHet ##
		if (($m_T >=1 || $r_T>=1) && ($cm_NT >=1 || $cr_NT >=1) && not exists $parent_seen{"$FAM_ID[0]"}) { ## if mother transmits one variant and doesn't transmit another - mother has a compound het
			foreach my $c (@CCH_m) {
				print WRITE_GENE "$trio_ID\tmother_compHet\t$key\t$c->[0]\t$c->[1]\t$c->[2]\t$c->[3]\t$c->[4]\t$c->[5]\t$c->[6]\t$c->[7]\t$c->[8]\t";
				}
			print WRITE_GENE "\n";
			}
		
	## test father compHet ##
	if (($f_T >=1 || $r_T >=1) && ($cf_NT >=1 || $cr_NT >=1) && not exists $parent_seen{"$FAM_ID[0]"}) { ## if father transmits one variant and doesn't transmit another - father has a compound het
		foreach my $c (@CCH_f) {
			print WRITE_GENE "$trio_ID\tfather_compHet\t$key\t$c->[0]\t$c->[1]\t$c->[2]\t$c->[3]\t$c->[4]\t$c->[5]\t$c->[6]\t$c->[7]\t$c->[8]\t";
			}
		print WRITE_GENE "\n";
		}
		
	}

	$parent_seen {"$FAM_ID[0]"} = 1; ##hash for analysed parents - prevents the same parent being counted twice in multiplex families
}	
	

sub geneHits {
	my ($varChr,$varPos,@Gene) = @_;
		$varChr=~s/X/23/g;
		$varChr=~s/Y/24/g;
		$varChr=~s/M/25/g;
	my $count=0;
		foreach my $x (@Gene) {
				next if ($x->[0] != $varChr);
				if ($varPos >= $x->[1] && $varPos <= $x->[2]) {
					$count++;
					return $x->[3];
					last;
					}
				}
		if ($count==0) {
			return "NA";
			}
	}

sub T_forloop { 
	my ($N,$chr,$pos,$gene,$p,$f,$m,$pINFO,$fINFO,$mINFO,$ref,$alt,$func,$exonicFunc,$BG_freq) = @_;
	for (my $i=1;$i<=$N;$i++) {
		print WRITE_T "$chr\t$pos\t$ref\t$alt\t$func\t$exonicFunc\t$gene\t$p\t$pINFO\t$f\t$fINFO\t$m\t$mINFO\t$BG_freq\n";
		}
	}
	
sub NT_forloop { 
	my ($N,$chr,$pos,$gene,$p,$f,$m,$pINFO,$fINFO,$mINFO,$ref,$alt,$func,$exonicFunc,$BG_freq) = @_;
	for (my $i=1;$i<=$N;$i++) {
		print WRITE_NT "$chr\t$pos\t$ref\t$alt\t$func\t$exonicFunc\t$gene\t$p\t$pINFO\t$f\t$fINFO\t$m\t$mINFO\t$BG_freq\n";
		}
	}	
	
sub CNV_hits {
	my ($FID, @SNPs) = @_;
	my (@CNVs, @CNV_pHit);
	open (CNV, "$CNV_in") || die "can't open CNV file";
	while (<CNV>) { 
		my $line = $_;
		chomp $line;
		$line=~s/\r//g;
		my @temp = split(/\t/,$line);
		$temp[1]=~s/chr//g;
		$temp[1]=~s/X/23/g;
		$temp[1]=~s/Y/24/g;
		$temp[1]=~s/M/25/g;
		my @temp2=split("-",$temp[0]);
		if ($temp2[0] eq $FID) {
			push (@CNVs,[@temp2,$temp[1],$temp[2],$temp[3],$temp[4],$temp[0]]);
			}
		}
	close (CNV);

	foreach my $x (@SNPs) {
		$x->[0]=~s/chr//g;
		$x->[0]=~s/X/23/g;
		$x->[0]=~s/Y/24/g;
		$x->[0]=~s/M/25/g;
		foreach my $y (@CNVs) {
			next if ($x->[0] != $y->[2]);
			if (($x->[1] >= $y->[3]) && ($x->[1] <= $y->[4])) { ##for each T_SNP search for overlapping CNVs##
				print "CNV_hit\n";
				push(@cnvSNP,[@{$x},@{$y}]);					
			}		
		}
		
	}
	
}

sub cal_BAF {
	my ($AD) = @_;
	my @AD2=split(",",$AD);
	if ($#AD2 != 1) {print "AD more than 2 numbers\n"; die;}
		my $TAF=$AD2[0] + $AD2[1];
		my $BAF=$AD2[1]/$TAF;
		return $BAF;		
}		
	

	
close (WRITE_GENE);
close (WRITE_TDT);
close (WRITE_T);
close (WRITE_NT);









