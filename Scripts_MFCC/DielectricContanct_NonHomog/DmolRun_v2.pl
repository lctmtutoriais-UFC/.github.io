#!perl
# Author: Pablo Abreu de Morais
# Data: 18/04/2020
# Versao: 2.0
# Script para rodar o cálculo quântico do MFCC

#!!!!!!!! Na linha 16 deve ser informado as iniciais dos nomes dos arquivos .xsd


use strict;
use Getopt::Long;
use MaterialsScript qw(:all);
use List::Util qw[min max];
use Cwd;

my $init_file_name='[A][E]';    	         # Inicial dos arquivos que serão lidos; neste caso, arquivos iniciados com AB ou CD 
                                                 # Caso queira adicionar mais uma "inicial dos arquivos", basta adicionar |[chain1][chain2]
                                                 # ex: [A][B]|[C][D]|[E][F]; a "|" significa "ou", logo serao considerados arquivos que iniciam com AB ou CD ou EF  

my $dielectric_function = "no";                  # yes (tabela das constantes dieletricas calculadas) ou no (constante dieletrica homogenea).
my $constant =40.00;
my $table3 = "dielectric.std";                   # Inserir tabela com valores das constantes dieletricas nao homogeneas
if($dielectric_function eq "yes"){$table3 = $Documents{"$table3"}}


my $table2 = $Documents{"Residues.std"};         #Carregando tabela(existente) "Residues.std"
my $table1 = Documents->New("Energy.std");

$table1 -> ColumnHeading(0)="ID_1";
$table1 -> ColumnHeading(1)="Res_1";
$table1 -> ColumnHeading(2)="Res_1 Chain";
$table1 -> ColumnHeading(3)="ID_2";
$table1 -> ColumnHeading(4)="Res_2";
$table1 -> ColumnHeading(5)="Res_2 Chain";
$table1 -> ColumnHeading(6)="Distance (Angstrons)";
$table1 -> ColumnHeading(7)="Dielectric Constant";
$table1 -> ColumnHeading(8)="Energy A (kcal/mol)"; 
$table1 -> ColumnHeading(9)="Energy_B (kcal/mol)"; 
$table1 -> ColumnHeading(10)="Energy_C (kcal/mol)"; 
$table1 -> ColumnHeading(11)="Energy_D (kcal/mol)"; 
$table1 -> ColumnHeading(12)="Interaction Energy (kcal/mol)"; 
$table1 -> ColumnHeading(13)="Res_1 Type";
$table1 -> ColumnHeading(14)="Res_1 Class";
$table1 -> ColumnHeading(15)="Res_2 Type";
$table1 -> ColumnHeading(16)="Res_2 Class";
$table1 -> ColumnHeading(17)="Res_1 Charge";
$table1 -> ColumnHeading(18)="Res_2 Charge";
$table1 -> ColumnHeading(19)="Res_1 Atom";
$table1 -> ColumnHeading(20)="Res_2 Atom";
$table1 -> ColumnHeading(21)="Water";
$table1 -> ColumnHeading(22)="n. Water";
$table1 -> ColumnHeading(23)="W/hbond";
$table1 -> ColumnHeading(24)="no hbond";
$table1 -> ColumnHeading(25)="Interaction Energy (Ha)";
$table1 -> ColumnHeading(26)="File A";
$table1 -> ColumnHeading(27)="File B";
$table1 -> ColumnHeading(28)="File C";
$table1 -> ColumnHeading(29)="File D";
$table1->UpdateViews;

my $E_A;
my $E_B;
my $E_C;
my $E_D;

my $E=0;    #Energia de interacao: E = E_A - E_B - E_C + E_D

my $c1='i';

my $c2='Id1';
my $c3='Res1';
my $c4='Chain1';
my $c5='Id2';
my $c6='Res2';
my $c7='Chain2';
my $c8='Dist';
my $c9='Const';
my $c10='EA';
my $c11='EB';
my $c12='EC';
my $c13='ED';
my $c14='E';
my $c15='R1_type';
my $c16='R1_class';
my $c17='R2_type';
my $c18='R2_class';
my $c19='R1_charge';
my $c20='R2_charge';
my $c21='R1_atm';
my $c22='R2_atm';
my $c23='Water';
my $c24='n.Water';
my $c25='w/bond';
my $c26='no_hbond';
my $c27='File A';
my $c28='File B';
my $c29='File C';
my $c30='File D';

printf "%-10s %-10s %-10s %-10s %-10s %-10s %-10s %-10s %-10s %-20s %-20s %-20s %-20s %-20s %-10s %-10s %-10s %-10s %-10s %-10s %-10s %-10s %-30s %-30s %-30s %-30s %-10s %-10s %-s %-s \n",$c1,$c2,$c3,$c4,
$c5,$c6,$c7,$c8,$c9,$c10,$c11,$c12,$c13,$c14,$c15,$c16,$c17,$c18,$c19,$c20,$c21,$c22,$c27,$c28,$c29,$c30,$c23,$c24,$c25,$c26;

my $cwd="DmolRun_v2_Files//Documents";   # cwd="[nome_do_script]_Files//Documents";

opendir my $folder, $cwd or die "Couldn't open folder";
my @allfiles = grep {!/^\.\.?$/} readdir $folder;                #Armazenando o nome de todos os arquivos

############ Inserir inicio do nome dos arquivos como variavel ####################
@allfiles = grep(/^$init_file_name.*?\.xsd$/, @allfiles);	 #Selecionando só os nomes que comecam com "$init_file_name" com extensao .xsd	

my $row1=-1;
my @allfiles= sort @allfiles;                                    # ordena lista de nomes dos arquivos

my $FA;
my $FB;
my $FC;
my $FD;
foreach my $file (@allfiles) {
 
	my $item=substr($file, 0, index($file, '.'));    #extrai todos os caracteres anteriores ao ponto (retira a extensão .xsd)
	$item=substr($item,-1,1);                        # Identifica se é arquivo A, B, C ou D 
	
	if($item eq "A"){$row1=$row1+1;}
	   	  
	#------------------------------------------------------------------------
	#	   Dividindo nomes dos arquivos em campos delimitados por "_"
	#------------------------------------------------------------------------
	
	my @words = split /_/,$file;
	   
	my $residue1=$words[2];
	my $residue2=$words[3];
	
	my $id1=$residue1;
	$id1=~ s/\D//g;
	
	my $id2=$residue2;
	$id2=~ s/\D//g;
		  	 
	my $chain1=substr($words[0],0,1);
	my $chain2=substr($words[0],1,1);
	  	 
	#print "$file,$chain1,$chain2,$id1,$id2 \n";  	 
	  	   
	#------------------------------------------------------------------------
	#				lendo tabela do MFCC (Residues.std)
	#------------------------------------------------------------------------
	my $nline = $table2->RowCount;
		
	my $i;
	my $dist;
	
	my $charge1;
	my $charge2;
	
	my $atom1;
	my $atom2;
	
	my $residue1_class;
	my $residue2_class;
	
	my $residue1_type;
	my $residue2_type;
	
	my $water;
	my $nwater;
	my $hbond;
	my $no_hbond;
	
	for ($i=0;$i<$nline;$i++){   #i < numero de linhas da tabela
		if($residue1 eq $table2->cell($i,0) and $chain1 eq $table2->cell($i,3) and $residue2 eq $table2->cell($i,1) and $chain2 eq $table2->cell($i,6)){		
			
			$dist    = $table2->cell($i,2);
			$charge1 = $table2->cell($i,7);
			$charge2 = $table2->cell($i,8);
			
			if($table2->cell($i,3) ne "L"){ 
				$residue1_type=substr($residue1,0,3);
			}
			else{
				$residue1_type="Ligand";
			}
			
			if($table2->cell($i,6) ne "L"){ 
				$residue2_type=substr($residue2,0,3);
			}
			else{
				$residue2_type="Ligand";
			}	
			
			$residue1_class=$table2->cell($i,9);
			$residue2_class=$table2->cell($i,10);
			
			$atom1=$table2->cell($i,4);
			$atom2=$table2->cell($i,5);
			
			$water=$table2->cell($i,11);
			$nwater=$table2->cell($i,12);
			$hbond=$table2->cell($i,13);
			$no_hbond=$table2->cell($i,14); 
			
			last;	
		}
	}
	
	$table1->cell($row1,0)=$id1;
	$table1->cell($row1,1)=$residue1;  
	$table1->cell($row1,2)=$chain1;	
	
	$table1->cell($row1,3)=$id2;
	$table1->cell($row1,4)=$residue2;
	$table1->cell($row1,5)=$chain2;
	
	$table1->cell($row1,6)=$dist;
		
	$table1->cell($row1,13)=$residue1_type;
	$table1->cell($row1,14)=$residue1_class;
	
	$table1->cell($row1,15)=$residue2_type;
	$table1->cell($row1,16)=$residue2_class;
	
	$table1->cell($row1,17)=$charge1;
	$table1->cell($row1,18)=$charge2;
	
	$table1->cell($row1,19)=$atom1;	
	$table1->cell($row1,20)=$atom2;
	
	
	$table1->cell($row1,21)=$water;
	$table1->cell($row1,22)=$nwater;
	$table1->cell($row1,23)=$hbond;	
	$table1->cell($row1,24)=$no_hbond;
	$table1->UpdateViews;
		
	#------------------------------------------------------------------------
	#	               Atribuindo constante dieletrica
	#------------------------------------------------------------------------	
	   
	if($dielectric_function eq "yes"){
		$nline=$table3->RowCount;
		
		for ($i=0;$i<$nline;$i++){   #i < numero de linhas da tabela
			if($residue1 eq $table3->cell($i,0) and $chain1 eq $table3->cell($i,1) and $residue2 eq $table3->cell($i,2) and $chain2 eq $table3->cell($i,3)){		
			$constant=$table3->cell($i,5);
			last;
			}
		}
	}
	
	$table1->cell($row1,7)=$constant;
	$table1->UpdateViews;
		
	#---------------------------------------------------------------------
	#                     Carregando Arquivo
	#---------------------------------------------------------------------
	
	my $doc = $Documents{"$file"};       #Carregando o arquivo .xsd
	
	#---------------------------------------------------------------------
	#              Calculo da carga formal (Formal charge)
	#---------------------------------------------------------------------
	
	my $formalChargeSum = 0;
	foreach my $atom (@{$doc->Atoms}){
		my $formalCharge = $atom->FormalCharge->Value;
		#printf "Formal charge on atom %s = %f/%f\n", $atom->Name, $formalCharge;
		$formalChargeSum += $formalCharge;
	}
	
	#---------------------------------------------------------------------
	#                Configurando Dmol  (Set up Dmol)
	#---------------------------------------------------------------------
	
	my $dmol3 = Modules->DMol3;
			
		        $dmol3->ChangeSettings([Quality => "Fine",
		        			
		        	UseSymmetry => "No",
			                     
			        ElectronicQuality => "Fine",
			                       
			        Charge => "$formalChargeSum",
			                        
			        Multiplicity => "Auto",
			                        
			        SpinUnrestricted => "Yes",
			                   
			        TheoryLevel => "GGA",	     #LDA ou GGA.
						
							     # Usar LocalFunciontal p/ LDA
							     # Usar NonLocalFunctional p/ GGA
						
			       #LocalFunctional => "PWC",    #Comentar caso use GGA.
	
				NonLocalFunctional=> "PBE",  #Apenas para GGA.
				UseDFTD => "Yes",
						
				DFTDMethod => "TS",	     #OBS, TS ou Grimme.
	
			        Basis => "DNP+", 	     #MIN, DN, DND, DNP, MIXED.
			                        
			        BasisFile => "4.4",	     #"3.5" ou "4.4"
			                        
			        MaximumSCFCycles => "1000",
			                        
			        UseSmearing => "YES",
			                        
			        Smearing => "0.005",
			                      
			        UseDIIS => "YES",
			                        
			        CutoffType => "Custom",
			       					
			       	AtomCutoff => "5.5",
			                       
			        CoreTreatment => "All Electron", 	# "All Electron"
			                        									# "Effectove core potentials"
			                        									# "All electrons relativistic"
			                        									# "DFT semi-core pseudoposts"
			        UseCosmo => "YES", 					#Solvente
			                        
			        #CosmoSolvent => "Water",		# "Water", "Dimethyl sulfoxide"
			                        
			        SolventDielectric => "$constant",	# 1 a 10000 (Solvent dielectric constat)
			                        
			                                  
				#PopulationAnalysis => "Yes",

				]);
	
	
	#---------------------------------------------------------------------
	#                     Executando Dmol  (Run Dmol)
	#---------------------------------------------------------------------
	
		my $output=$dmol3->Energy->Run($doc);
		my $FinalEnergy=$output->TotalEnergy;
	 	 
	 	#my $FinalEnergy=20;
	 	my $mfcc_count=$row1+1;
	 	 	
	#---------------------------------------------------------------------
	#           Escrevendo as energias na tabela do Dmol (Energy.std)
	#---------------------------------------------------------------------
	
	
	if ($item eq "A"){
		$E_A = sprintf ("%20.10f",$FinalEnergy);
		my $E_A_aux=$E_A;
		$E_A_aux=~ s/\./,/;
		
		$table1->cell($row1,8)=$E_A_aux;
		$table1->cell($row1,26)=$doc;
		$FA=$file;
		$table1->UpdateViews;
	}
	
	if ($item eq "B"){
		$E_B = sprintf ("%20.10f",$FinalEnergy);
		my $E_B_aux=$E_B;
		$E_B_aux=~ s/\./,/;
		
		
		$table1->cell($row1,9)=$E_B_aux;
		$table1->cell($row1,27)=$doc;
		$FB=$file;
		$table1->UpdateViews;
	}
	
	if ($item eq "C"){
		$E_C = sprintf ("%20.10f",$FinalEnergy);
		my $E_C_aux=$E_C;
		$E_C_aux=~ s/\./,/;		
		
		$table1->cell($row1,10)=$E_C_aux;
		$table1->cell($row1,28)=$doc;
		$FC=$file;
		$table1->UpdateViews;
	}
	
	if ($item eq "D"){
		$E_D = sprintf ("%20.10f",$FinalEnergy);
		my $E_D_aux=$E_D;
		$E_D_aux=~ s/\./,/;
		
		$table1->cell($row1,11)=$E_D_aux;
		$E = $E_A-$E_B-$E_C+$E_D;
		$table1->cell($row1,12)=$E;
		$table1->cell($row1,29)=$doc;
		$FD=$file;
		$table1->cell($row1,25)=$E/627.509;	
		$table1->UpdateViews;

		printf "%-10d %-10d %-10s %-10s %-10d %-10s %-10s %-10.2f %-10.2f %-20.10f %-20.10f %-20.10f %-20.10f %-20.10f %-10s %-10s %-10s %-10s %-10d %-10d %-10s %-10s %-30s %-30s %-30s %-30s %-10s %-10s %-s %-s\n",
		$mfcc_count,$id1,$residue1,$chain1,$id2,$residue2,$chain2,$dist,$constant,$E_A,$E_B,$E_C,$E_D,$E,$residue1_type,$residue1_class,$residue2_type,
		$residue2_class,$charge1,$charge2,$atom1,$atom2,$FA,$FB,$FC,$FD,$water,$nwater,$hbond,$no_hbond;

	}
	
	$doc->Close;
}

