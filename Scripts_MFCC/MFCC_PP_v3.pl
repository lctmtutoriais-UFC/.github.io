#!perl
# Author: Pablo Abreu de Morais
# Data: 22/03/2020
# Versao: 3.0
# Script para fragmentar um sistema proteina/ligante ou proteina/proteina conforme MFCC

# !!!!!!!!!!!  As cadeias da proteina devem ter nomes diferentes (exceto "L") - ex: A, B, C, etc.
# !!!!!!!!!!!  A cadeia do(s) ligante(s) deve ter nome "L". NENHUMA cadeia da proteína pode se chamar "L".
# !!!!!!!!!!!  As moleculas de água devem ter sua(s) cadeia(s) nomeada(s) de "Water". Pode existir mais de uma cadeia de agua chamada "Water"
# !!!!!!!!!!!  CUIDADO: As ligacoes de H podem desaparecer nos arquivos do MFCC (A, B, C e D) devido a adicao de H para completar
# !!!!!!!!!!!  		as ligacoes peptidas, pois essa adicao gera uma torcao na estrutura do N.

use strict;
use Getopt::Long;
use MaterialsScript qw(:all);
use List::Util qw(min max);

#------------------------------------------------------------------------------------------------------------------------------------------
#                                                 Carregando input
#------------------------------------------------------------------------------------------------------------------------------------------

my $water="no";         # Atribua 'yes' caso queira incluir as moleculas de agua no MFCC
my $water_dist=2.5;      # Distancia maxima das moleculas de agua em relacao ao ligante, residuos e caps (deve ser <= a 5A)
my $R_cutoff=8.0;        # Raio de interacao. O valor default (1000000) considera todas as interacoes 
my $caps_number=1;       # Numero de caps considerado no MFCC. O valor default e 1.  
my $doc = $Documents{"input.xsd"};  #Fornecer o nome do arquivo .xsd ("nome_do_arquivo.xsd")

# !!!!!!!!!!!  O script calculara o MFCC para as interacoes das cadeias da @list_1 com as cadeias da @list_2
# !!!!!!!!!!!  Caso deseje calcular interacoes "intrachain", @list_1 e @list_2 devem ser iguais, ex:
# !!!!!!!!!!!  @list_1 = ('A','B') e @list_2 = ('A','B') -> Interacoes calculadas: AA, AB e BB
 
my @list_1=('B');
my @list_2=('L');
 
#------------------------------------------------------------------------------------------------------------------------------------------
#                                                 Calculando ligacoes
#------------------------------------------------------------------------------------------------------------------------------------------

Tools->BondCalculation->ChangeSettings(Settings(MinBondLength => 0.75, MaxBondLength => 1.20, ResonantBondRepresentation=>"Resonant")); # Calcula as ligações no arquivo xsd, pois essa informação se perde na conversão do pdb para o xsd
Tools->BondCalculation->Bonds->Calculate($doc);

#------------------------------------------------------------------------------------------------------------------------------------------
#                                             Criando tabela Residues.std
#------------------------------------------------------------------------------------------------------------------------------------------

my $table2 = Documents->New("Residues.std");
$table2 -> ColumnHeading(0) = "Residue_1"; 
$table2 -> ColumnHeading(1) = "Residue_2"; 
$table2 -> ColumnHeading(2) = "Distance"; 
$table2 -> ColumnHeading(3) = "Residue_1 chain"; 
$table2 -> ColumnHeading(4) = "Residue_1 atom"; 
$table2 -> ColumnHeading(5) = "Residue_2 atom"; 
$table2 -> ColumnHeading(6) = "Residue_2 chain"; 
$table2 -> ColumnHeading(7) = "Residue_1 charge";
$table2 -> ColumnHeading(8) = "Residue_2 charge";
$table2 -> ColumnHeading(9) = "Residue_1 classif";
$table2 -> ColumnHeading(10) = "Residue_2 classif";

$table2 -> ColumnHeading(11) = "Water";
$table2 -> ColumnHeading(12) = "n. water";
$table2 -> ColumnHeading(13) = "w/hbond";
$table2 -> ColumnHeading(14) = "no hbond";
$table2 -> ColumnHeading(15) = ".xsd";


#------------------------------------------------------------------------------------------------------------------------------------------
#                           Alterando nome dos residuos/Cadeia e transladando PDB para primeiro quadrante
#------------------------------------------------------------------------------------------------------------------------------------------

$doc=pdb_quadrante1($doc);
change_res_name($doc);
my $doc1=no_water($doc);

#------------------------------------------------------------------------------------------------------------------------------------------
#                     Criando listas com os id's dos residuos (list_chains)
#------------------------------------------------------------------------------------------------------------------------------------------

my $id_res=-1;   #id (global) do residuo

my $nchain_1=0;
my $nchain_2=0;

my $start_name_1="unknown";        
my $start_name_2="unknown";

my @list_chains_1;
my @list_chains_2;


foreach my $residue (@{$doc1->SubUnits}){

	$id_res=$id_res+1;
	my $chain=$residue->Ancestors->Chain->Name;

	foreach my $chain_name_aux (@list_1){
	
		if ($chain_name_aux eq $chain){
			
			if($chain_name_aux ne $start_name_1){
				$nchain_1=$nchain_1+1;
				$start_name_1=$chain_name_aux
			}		
			
			push(@{$list_chains_1[$nchain_1]},$id_res);
		}
	}

	
	foreach my $chain_name_aux (@list_2){
	
		if ($chain_name_aux eq $chain){
			
			if($chain_name_aux ne $start_name_2){
				$nchain_2=$nchain_2+1;
				$start_name_2=$chain_name_aux
			}		
			
			push(@{$list_chains_2[$nchain_2]},$id_res);
		}
	}
	
}

#------------------------------------------------------------------------------------------------------------------------------------------
#                                                          MFCC 
#------------------------------------------------------------------------------------------------------------------------------------------

my ($doc1a,$doc2a,$doc3a,$doc4a);
my $row=-1;
my $charge1;
my $charge2;
my $ok;
	
if($water eq "yes"){

	#------------------------------------------------------------------------------------------------------------------------------------------
        #                            Identificando pontes de hidrogenio dos residuos
        #------------------------------------------------------------------------------------------------------------------------------------------
	$doc->CalculateHBonds;	
	
	my @hbond_array;
	my $nbond=-1;
	foreach my $hbond (@{$doc->HydrogenBonds}){
		
		$nbond=$nbond+1;	
		my $atom1=$hbond->Acceptor;
		my $atom2=$hbond->Donor;
		
		my $res1=$atom1->Ancestors->SubUnit->Name;
		my $res2=$atom2->Ancestors->SubUnit->Name;
		
		my $chain1=$atom1->Ancestors->SubUnit->Ancestors->Chain->Name;
		my $chain2=$atom2->Ancestors->SubUnit->Ancestors->Chain->Name;
		
		my $id1=id_residue($doc,$res1);
		my $id2=id_residue($doc,$res2);	
		
		if($chain1 eq "Water" ^ $chain2 eq "Water"){
			#print "$nbond,$res1,$id1,$res2,$id2 \n";
			
					
			if($chain1 eq "Water"){
				push(@{$hbond_array[$id2]},$id1);
			}
			else{
				push(@{$hbond_array[$id1]},$id2);
			}
		}
	}
	
	@hbond_array=uniq_row_array(@hbond_array);
	
	
	# Parametros da rede cubica

	my $dmax=length_cube($doc);
	my $r_cube=12;                   #largura dos cubos pequenos
	my $L=int($dmax/$r_cube)+1;      #Quantidade de cubos por "linha"

	# Construindo rede cubica e mapeando moleculas de agua

	my @cube_array=cube_space($doc,$r_cube,$L);       

	# Associando moléculas de água aos residuos

	my @water_array=water_array($doc1,$doc,$r_cube,$L,\@cube_array,$water_dist);

	for (my $i=1;$i<=$nchain_1;$i++){
		
		my $cap_min_res1 = min(@{$list_chains_1[$i]});
		my $cap_max_res1 = max(@{$list_chains_1[$i]});
			
		foreach my $id1 (@{$list_chains_1[$i]}){	
			
			my $res1=$doc1->SubUnits->Item($id1)->Name;
			my $chain1=$doc1->SubUnits->Item($id1)->Ancestors->Chain->Name;
									
			for (my $j=1;$j<=$nchain_2;$j++){
			
				my $cap_min_res2 = min(@{$list_chains_2[$j]});
				my $cap_max_res2 = max(@{$list_chains_2[$j]});
				
				foreach my $id2 (@{$list_chains_2[$j]}){	
				
					my $res2=$doc1->SubUnits->Item($id2)->Name;
					my $chain2=$doc1->SubUnits->Item($id2)->Ancestors->Chain->Name;
					my($dist,$atom1,$atom2)=res_dist($doc1,$id1,$id2);
					
					$ok=bond($doc1,$id1,$id2);
					
					if($chain1 eq $chain2 and $id2<=$id1){$ok="yes"}
					if(@list_chains_1==@list_chains_2 and $j<$i){$ok="yes"}
					
					if($dist<=$R_cutoff and $ok eq "no"){				
						$row=$row+1;
						
						print "$res1,$res2,$dist \n";
			
						my $ok_water="no";  # status agua; yes -> participacao de H2O na interacao
						my $nw;             # quantidade total de moleculas de agua considerada
						my $n1;
						my $n2;
						
						$charge1=formal_charge($doc1,$id1);
						$charge2=formal_charge($doc1,$id2);
						
						my $res1_class=residue_class($res1);
						my $res2_class=residue_class($res2);
						
						my @caps1=res_caps($doc1,$id1,$cap_max_res1,$cap_min_res1,$caps_number);
						my @caps2=res_caps($doc1,$id2,$cap_max_res2,$cap_min_res2,$caps_number);
				
						my $c1=-1;
						my $c2=-1;
						
						if(substr($res1,0,3) eq "CYS"){$c1=cap_extra($doc1,$id1)}
						if(substr($res2,0,3) eq "CYS"){$c2=cap_extra($doc1,$id2)}
						
						if($c1>=0){
							push(@caps1,$c1);
							#print "$id1,$res1,@caps1 \n";
						}
						
						if($c2>=0){	
							push(@caps2,$c2);	
							#print "$id2,$res2,@caps2 \n";
						}
						
						my $aux=$chain1.$chain2;
						
						($doc1a,$doc2a,$doc3a,$doc4a)=create_mfcc_files_cysteine($doc1,$id1,$id2,$dist,$aux);
						mfcc_cysteine($doc1,$id1,$id2,\@caps1,\@caps2); 
						
						my @list_water=water_list($id1,$id2,\@caps1,\@caps2,\@water_array);
						copy_water($doc,\@list_water);
								
						my @list_hbond=hbond_index($id1,$id2,\@caps1,\@caps2,\@hbond_array); # lista com águas que fazem ponte de hidrogenio
						
						my @list_aux = grep!${{map{$_,1}@list_hbond}}{$_},@list_water;       #lista com aguas atendem ao criterio do cutoff e nao fazem ponte de hidrogenio  
						
						#print "Lista water -> @list_water \n";
						#print "Lista hbond -> @list_hbond \n";
						#print "Lista aux   -> @list_aux \n";
						
						$nw=$#list_water+1;
						
						
						#print "--------------- Lista water ---------------------\n";
						#for (my $i=0;$i<=$#list_water;$i++){
								#my $index=$list_water[$i];
								#$index=$doc->SubUnits->Item($index)->Name;
								#print "$i,$index \n";
						#}
						
						#print "--------------- Lista hbond ---------------------\n";
						#for (my $i=0;$i<=$#list_hbond;$i++){
								#my $index=$list_hbond[$i];
								#$index=$doc->SubUnits->Item($index)->Name;
								#print "$i,$index \n";
						#}
						
						#print "--------------- Lista aux ---------------------\n";
						#for (my $i=0;$i<=$#list_aux;$i++){
								#my $index=$list_aux[$i];
								#$index=$doc->SubUnits->Item($index)->Name;
								#print "$i,$index \n";
						#}
						
						
						if($nw>=1){
							
							
							$ok_water="yes";
							
							for (my $i=0;$i<=$#list_hbond;$i++){
								my $index=$list_hbond[$i];
								$list_hbond[$i]=$doc->SubUnits->Item($index)->Name;
								$n1=$n1.$list_hbond[$i].', ';			
								#print "------------- $i,$list_hbond[$i],$n1 -------------- \n";
							}
							
							for (my $i=0;$i<=$#list_aux;$i++){
								my $index=$list_aux[$i];
								$list_aux[$i]=$doc->SubUnits->Item($index)->Name;
								$n2=$n2.$list_aux[$i].', ';
								#print "############## $i,$list_aux[$i],$n2 ############ \n";
							}
						
						}
						else{
						$n1="-" ;
						$n2="-";
						}
						
										
						$table2->cell($row,0)  = substr($res1,0,index($res1, '_'));
						$table2->cell($row,1)  = substr($res2,0,index($res2, '_'));
						$table2->cell($row,2)  = $dist;
						$table2->cell($row,3)  = $chain1;
						$table2->cell($row,4)  = $atom1;
						$table2->cell($row,5)  = $atom2;
						$table2->cell($row,6)  = $chain2;
						$table2->cell($row,7)  = $charge1;
						$table2->cell($row,8)  = $charge2;
						$table2->cell($row,9)  = $res1_class;
						$table2->cell($row,10) = $res2_class;	
					
						$table2->cell($row,11) = $ok_water;      # status agua; yes -> participacao de H2O na interacao
						$table2->cell($row,12) = $nw;            # quantidade total de moleculas de agua consideradas
						$table2->cell($row,13) = $n1;            # moleculas de agua com ponte de hidrogenio      
						$table2->cell($row,14) = $n2;            # moleculas de agua sem ponte de hidrogenio
						$table2->cell($row,15) = $doc1a;         # moleculas de agua sem ponte de hidrogenio
						$table2->UpdateViews;
					}
				}
			}
		}
	}
}
else
{
	for (my $i=1;$i<=$nchain_1;$i++){
		
		my $cap_min_res1 = min(@{$list_chains_1[$i]});
		my $cap_max_res1 = max(@{$list_chains_1[$i]});
			
		foreach my $id1 (@{$list_chains_1[$i]}){	
			
			my $res1=$doc1->SubUnits->Item($id1)->Name;
			my $chain1=$doc1->SubUnits->Item($id1)->Ancestors->Chain->Name;
									
			for (my $j=1;$j<=$nchain_2;$j++){
			
				my $cap_min_res2 = min(@{$list_chains_2[$j]});
				my $cap_max_res2 = max(@{$list_chains_2[$j]});
				
				foreach my $id2 (@{$list_chains_2[$j]}){	
				
					my $res2=$doc1->SubUnits->Item($id2)->Name;
					my $chain2=$doc1->SubUnits->Item($id2)->Ancestors->Chain->Name;
					my($dist,$atom1,$atom2)=res_dist($doc1,$id1,$id2);
					
					$ok=bond($doc1,$id1,$id2);
					
					if($chain1 eq $chain2 and $id2<=$id1){$ok="yes"}
					if(@list_chains_1==@list_chains_2 and $j<$i){$ok="yes"}
					
					if($dist<=$R_cutoff and $ok eq "no"){				
						$row=$row+1;
						
						print "$res1,$res2,$dist \n";
						
						$charge1=formal_charge($doc1,$id1);
						$charge2=formal_charge($doc1,$id2);
						
						my $res1_class=residue_class($res1);
						my $res2_class=residue_class($res2);
						
						my @caps1=res_caps($doc1,$id1,$cap_max_res1,$cap_min_res1,$caps_number);
						my @caps2=res_caps($doc1,$id2,$cap_max_res2,$cap_min_res2,$caps_number);
				
						my $c1=-1;
						my $c2=-1;
						
						if(substr($res1,0,3) eq "CYS"){$c1=cap_extra($doc1,$id1)}
						if(substr($res2,0,3) eq "CYS"){$c2=cap_extra($doc1,$id2)}
						
						if($c1>=0){
							push(@caps1,$c1);
							#print "$id1,$res1,@caps1 \n";
						}
						
						if($c2>=0){	
							push(@caps2,$c2);	
							#print "$id2,$res2,@caps2 \n";
						}
						
						my $aux=$chain1.$chain2;
						
						($doc1a,$doc2a,$doc3a,$doc4a)=create_mfcc_files_cysteine($doc1,$id1,$id2,$dist,$aux);
						mfcc_cysteine($doc1,$id1,$id2,\@caps1,\@caps2);
													
						$table2->cell($row,0)  = substr($res1,0,index($res1, '_'));
						$table2->cell($row,1)  = substr($res2,0,index($res2, '_'));
						$table2->cell($row,2)  = $dist;
						$table2->cell($row,3)  = $chain1;
						$table2->cell($row,4)  = $atom1;
						$table2->cell($row,5)  = $atom2;
						$table2->cell($row,6)  = $chain2;
						$table2->cell($row,7)  = $charge1;
						$table2->cell($row,8)  = $charge2;
						$table2->cell($row,9)  = $res1_class;
						$table2->cell($row,10) = $res2_class;
						$table2->cell($row,11) = "NA";      
						$table2->cell($row,12) = "NA";            
						$table2->cell($row,13) = "NA";                  
						$table2->cell($row,14) = "NA";            
						$table2->cell($row,15) = $doc1a;         
						$table2->UpdateViews;
					}
				}
			}
		}
	}
	
}

#------------------------------------------------------------------------------------------------------------------------------------------
#                                                   Fim do programa principal/Inicio das funcoes
#------------------------------------------------------------------------------------------------------------------------------------------


#-------------------------------------------------------------------------
#                           sub no_water
#-------------------------------------------------------------------------

sub no_water{
	my $doc_aux=@_[0];

	my $doc1=Documents->New("PDB_no_water.xsd");
	$doc1->CopyFrom($doc_aux);
		
	foreach my $chain (@{$doc1->Chains}){
		my $name=$chain->Name;
		
		if($name eq "Water"){$chain->Delete}
	}

	return $doc1;
}

#-------------------------------------------------------------------------
#                           sub bond
#-------------------------------------------------------------------------

sub bond{
	my $doc_aux=@_[0];
	my $id1=@_[1];
	my $id2=@_[2];

	my $ok="no";

	my $res1=$doc_aux->SubUnits->Item($id1)->Name;
	my $res2=$doc_aux->SubUnits->Item($id2)->Name;

	my $chain1=$doc_aux->SubUnits->Item($id1)->Ancestors->Chain->Name;
	my $chain2=$doc_aux->SubUnits->Item($id2)->Ancestors->Chain->Name;
	
	if($chain1 eq "L" or $chain2 eq "L"){
		$ok="no";
	}
	else{
	
		if(substr($res1,0,3) eq "CYS" and substr($res2,0,3) eq "CYS"){
			if(cap_extra($doc1,$id1)==$id2){$ok eq "yes"}
		}
		else{
		
			if($id1==$id2){
			   $ok="yes"  	
			}		
			
			if($id2>$id1){
				
				my $C=$doc_aux->SubUnits->Item($id1)->Atoms->DisplayRange("C");
				my $nbonds=$C->AttachedAtoms->Count;
						
				for(my $j=0;$j<=$nbonds-1;$j++){
					my $neighbor=$C->AttachedAtoms->Item($j)->Ancestors->SubUnit->Name;
					if($neighbor eq $res2){$ok="yes"}	
				}
			}
				
			if($id2<$id1){
			
				my $C=$doc_aux->SubUnits->Item($id2)->Atoms->DisplayRange("C");
				my $nbonds=$C->AttachedAtoms->Count;
						
				for(my $j=0;$j<=$nbonds-1;$j++){
					my $neighbor=$C->AttachedAtoms->Item($j)->Ancestors->SubUnit->Name;
					if($neighbor eq $res1){$ok="yes"}	
				}
			}
		}	
	}
	return $ok;
}


#-------------------------------------------------------------------------
#                           sub res_caps -> FUNCAO ATUALIZADA
#-------------------------------------------------------------------------

sub res_caps{

	my $doc_aux=@_[0];
	my $id_residue=@_[1];
	my $id_max=@_[2];
	my $id_min=@_[3];
	my $caps_number=@_[4];
	
	my $i;
	my $ok;
	my @caps;

	my $chain=$doc_aux->SubUnits->Item($id_residue)->Ancestors->Chain->Name;
	
	if($chain ne "L"){
	
	#            -------------------------
	#                  Caps à direita
	#            -------------------------
	
		for ($i=1;$i<=$caps_number;$i++){
	
			my $id_cap_r = $id_residue + $i;
			my $id=$id_cap_r-1;
			
			if($id_cap_r <= $id_max){	
				$ok=bond($doc_aux,$id,$id_cap_r);
				#print "$id,$id_cap_r,$ok \n";
				
				if($ok eq "yes"){
					push(@caps,$id_cap_r);
					#print "$id_residue,$id_cap_r \n";
				}
				if($ok eq "no"){last}
			}
		}
	
	#            -------------------------
	#                  Caps à esquerda
	#            -------------------------
	
		for ($i=1;$i<=$caps_number;$i++){
	
			my $id_cap_l = $id_residue - $i;
			my $id=$id_cap_l+1;
		
			if($id_cap_l >= $id_min){	
				$ok=bond($doc_aux,$id,$id_cap_l);
				#print "$id,$id_cap_l,$ok \n";
				if($ok eq "yes"){
					push(@caps,$id_cap_l);
					#print "$id_residue,$id_cap_l \n";
				}
				if($ok eq "no"){last}
			}
		}
	}
	
	return @caps
}

#-------------------------------------------------------------------------
#                           sub cap_extra
#-------------------------------------------------------------------------

sub cap_extra{
	
	my $doc_aux=@_[0];
	my $id1=@_[1];
	
	my $cap_extra=-1;
	
	
	my $chain=$doc_aux->SubUnits->Item($id1)->Ancestors->Chain->Name;
		
	
	if($chain ne "L"){
	
	#            -------------------------
	#                     Cap extra
	#            -------------------------
			
		my $res= $doc_aux->SubUnits->Item($id1)->Name;		
		my $SG1=$doc_aux->SubUnits->Item($id1)->Atoms->DisplayRange("SG");
		my $nbonds=$SG1->AttachedAtoms->Count;
	
			for(my $i=0;$i<=$nbonds-1;$i++){
				my $res2=$SG1->AttachedAtoms->Item($i)->Ancestors->SubUnit->Name;
				my $SG2=$SG1->AttachedAtoms->Item($i)->Name;
			
				if($SG2 eq "SG"){
					$cap_extra=id_residue($doc_aux,$res2);
				}	
			}
	}

	return $cap_extra;
}

#-------------------------------------------------------------------------
#                           sub res_dist
#-------------------------------------------------------------------------

sub res_dist{

	my $doc_aux=@_[0];
	my $i=@_[1];
	my $j=@_[2];
	
	my $dist_min=1000;
	my $atom1_name;
	my $atom2_name;
	
	foreach my $atom1 (@{$doc_aux->SubUnits->Item($i)->Atoms}){
			foreach my $atom2 (@{$doc_aux->SubUnits->Item($j)->Atoms}){
				my $dist = sqrt ( ($atom1->X - $atom2->X)**2 + ($atom1->Y - $atom2->Y)**2 + ($atom1->Z - $atom2->Z)**2 );
				if ($dist<$dist_min){
						$dist_min=$dist;
						$atom1_name=$atom1->Name;
						$atom2_name=$atom2->Name;
				}	
			}
		}

	return ($dist_min,$atom1_name,$atom2_name);

}


#-------------------------------------------------------------------------
#                           sub res_dist_files
#-------------------------------------------------------------------------

sub res_dist_files{

	my $doc_aux=@_[0];
	my $doc2_aux=@_[1];
	my $i=@_[2];
	my $j=@_[3];
	
	my $dist_min=1000;
	my $atom1_name;
	my $atom2_name;
	
	foreach my $atom1 (@{$doc_aux->SubUnits->Item($i)->Atoms}){
			foreach my $atom2 (@{$doc2_aux->SubUnits->Item($j)->Atoms}){
				my $dist = sqrt ( ($atom1->X - $atom2->X)**2 + ($atom1->Y - $atom2->Y)**2 + ($atom1->Z - $atom2->Z)**2 );
				if ($dist<$dist_min){
						$dist_min=$dist;
						$atom1_name=$atom1->Name;
						$atom2_name=$atom2->Name;
				}	
			}
		}

	return ($dist_min,$atom1_name,$atom2_name);

}


#-------------------------------------------------------------------------
#                      sub create_mfcc_files_cysteine
#-------------------------------------------------------------------------

sub create_mfcc_files_cysteine{

		my $doc_aux=@_[0];
		my $id1=@_[1];
		my $id2=@_[2];
		my $dist=@_[3];
		my $aux=@_[4];
		
		my $dist=int($dist*1000);
						
		my $res1=$doc_aux->SubUnits->Item($id1)->Name;
		my $res2=$doc_aux->SubUnits->Item($id2)->Name;

		$res1=substr($res1,0,index($res1, '_'));
		$res2=substr($res2,0,index($res2, '_'));
		
		my $name;
		if($dist < 10000){$name="0"}
		if($dist > 9999){$name=""}
		
		my $new_name=$aux.'_'.$name.$dist .'_'.$res1.'_'.$res2;
	        
	        my $new_name2= $new_name . "_A";		
		my $doc1a=Documents->New("$new_name2.xsd");	
	
		$new_name2= $new_name . "_B";
		my $doc2a=Documents->New("$new_name2.xsd");
		
		$new_name2= $new_name . "_C";
		my $doc3a=Documents->New("$new_name2.xsd");
		
	        $new_name2= $new_name . "_D";
		my $doc4a=Documents->New("$new_name2.xsd");

		return($doc1a,$doc2a,$doc3a,$doc4a);
}

#-------------------------------------------------------------------------
#                      sub change_res_name
#-------------------------------------------------------------------------

sub change_res_name{

	my $doc_aux=@_[0];
	
	my $id=-1;
	foreach my $residue (@{$doc_aux->SubUnits}){
		$id=$id+1;
		my $chain=$residue->Ancestors->Chain->Name;
		$residue->Name=$residue->Name.'_'.$chain.'_'.$id;
	}
	
	return $doc_aux;
}

#-------------------------------------------------------------------------
#                           sub id_residue
#-------------------------------------------------------------------------

sub id_residue{
	
	my $doc_aux=@_[0];
	my $residue_name=@_[1];
	
	my $id_residue;
	my $id_aux=-1;
	foreach my $residue (@{$doc_aux->SubUnits}){
		$id_aux=$id_aux+1;
		my $name = $residue->Name;
		if($residue_name eq "$name"){
			$id_residue=$id_aux;
			last
		}
	}

	return $id_residue;	
}


#-------------------------------------------------------------------------
#                           sub mfcc_cysteine -> FUNCAO ATUALIZADA
#-------------------------------------------------------------------------

sub mfcc_cysteine{ 
	
	my $doc_aux=@_[0];
	my $id1=@_[1];
	my $id2=@_[2];
	my @caps1=@{$_[3]};
	my @caps2=@{$_[4]};

	my $res1=$doc_aux->SubUnits->Item($id1)->Name;
	my $res2=$doc_aux->SubUnits->Item($id2)->Name;
	
	my @id_list;
	
	push(@id_list,$id1);
	push(@id_list,$id2);
	push(@id_list,@caps1);
	push(@id_list,@caps2);
	
	#--------------------------------------------------------------------
	#                   Arquivo A  ->  k*-R1-k + C*-R2-C    
	#--------------------------------------------------------------------
	
	foreach my $id_aux (@id_list){
		$doc1a->CopyFrom($doc_aux->SubUnits->Item($id_aux));
	}	
	
	foreach my $atom (@{$doc1a->Atoms}){
		if (($atom->Name eq "N" and $atom->NumBonds<3) or ($atom->Name eq "SG" and $atom->NumBonds<2)){
			Tools->BondCalculation->Bonds->Calculate($atom);
		}
	}

	foreach my $atom (@{$doc1a->Atoms}){
		if ($atom->Name eq "N" and $atom->NumBonds<3){$atom->AdjustHydrogen}
		if ($atom->Name eq "C" and $atom->NumBonds<3){$atom->AdjustHydrogen}
		if ($atom->Name eq "SG" and $atom->NumBonds<2){$atom->AdjustHydrogen}
	}
	
	#--------------------------------------------------------------------
	#                  Arquivo B  ->   k*-R1-k +  C*- C
	#--------------------------------------------------------------------
	
	$doc2a->CopyFrom($doc1a);
	$doc2a->SubUnits->DisplayRange($res2)->Delete;

	foreach my $atom (@{$doc2a->Atoms}){
		if ($atom->Name eq "N" and $atom->NumBonds<3){$atom->AdjustHydrogen}
		if ($atom->Name eq "C" and $atom->NumBonds<3){$atom->AdjustHydrogen}
		if ($atom->Name eq "SG" and $atom->NumBonds<2){$atom->AdjustHydrogen}
	}

	#--------------------------------------------------------------------
	#                  Arquivo C  ->  k*-k + C*-R2-C  
	#--------------------------------------------------------------------
	
	$doc3a->CopyFrom($doc1a);
	$doc3a->SubUnits->DisplayRange($res1)->Delete;
	
	foreach my $atom (@{$doc3a->Atoms}){
		if ($atom->Name eq "N" and $atom->NumBonds<3){$atom->AdjustHydrogen}
		if ($atom->Name eq "C" and $atom->NumBonds<3){$atom->AdjustHydrogen}
		if ($atom->Name eq "SG" and $atom->NumBonds<2){$atom->AdjustHydrogen}
	}

	#--------------------------------------------------------------------
	#                       Arquivo D  ->   k*-k + C*-C   
	#--------------------------------------------------------------------
	
	$doc4a->CopyFrom($doc1a);
	$doc4a->SubUnits->DisplayRange($res1)->Delete;
	$doc4a->SubUnits->DisplayRange($res2)->Delete;

	foreach my $atom (@{$doc4a->Atoms}){
		if ($atom->Name eq "N" and $atom->NumBonds<3){$atom->AdjustHydrogen}
		if ($atom->Name eq "C" and $atom->NumBonds<3){$atom->AdjustHydrogen}
		if ($atom->Name eq "SG" and $atom->NumBonds<2){$atom->AdjustHydrogen}
	}
	
}

#-------------------------------------------------------------------------
#                           sub print_2d
#-------------------------------------------------------------------------

sub print_2d {
	my @array_2d=@_;
	for(my $i = 0; $i <= $#array_2d; $i++){
	   for(my $j = 0; $j <= $#{$array_2d[$i]} ; $j++){
	      print "$array_2d[$i][$j] ";
	   }
	   print "\n";
	}
}

#--------------------------------------------------------------------
#                         sub uniq_row_array   
#--------------------------------------------------------------------
	
sub uniq_row_array {
	my @array_2d=@_;
	
   	for(my $i = 0; $i <= $#array_2d; $i++){ @{$array_2d[$i]}=uniq(@{$array_2d[$i]})}
	
	return @array_2d;

}

#-------------------------------------------------------------------------
#                          sub length_cube
#-------------------------------------------------------------------------

sub length_cube{
	
	my $doc_aux=@_[0];
		
	my $dmax=-10000;
	foreach my $atom (@{$doc_aux->Atoms}){
		my $x=$atom->X;
		my $y=$atom->Y;
		my $z=$atom->Z;
	
		if($x>$dmax){$dmax=abs($x)}
		if($y>$dmax){$dmax=abs($y)}
		if($z>$dmax){$dmax=abs($z)}
	
	}

	return $dmax;	
}

#-------------------------------------------------------------------------
#                        sub id_cube
#-------------------------------------------------------------------------

sub id_cube{

	my $r=@_[0];
	my $L=@_[1];
	my $x=@_[2];
	my $y=@_[3];
	my $z=@_[4];
	
	my $idx=int($x/$r)+1;
	my $idy=int($y/$r)+1;
	my $idz=int($z/$r)+1;

	my $id = $idx + $L*($idy-1) + $L*$L*($idz-1);
	
	return $id;

}

#-------------------------------------------------------------------------
#                       sub cube_space
#-------------------------------------------------------------------------

sub cube_space{

	my $doc_aux=@_[0];
	my $r=@_[1];
	my $L=@_[2];
	
	my $id_cube;
	my $idx;
	my $idy;
	my $idz;
	
	my @cube_array;
	
	my $id_residue=-1;
	foreach my $residue (@{$doc_aux->SubUnits}){
		$id_residue=$id_residue+1;
		
		my $name=$residue->Name;
		my $chain=$residue->Ancestors->Chain->Name;
		
		if($chain eq "Water"){
		
			my $atom=$residue->Atoms->Item(0); 
								
			my $x=$atom->X;
			my $y=$atom->Y;
			my $z=$atom->Z;
					
			$idx=int($x/$r)+1;
			$idy=int($y/$r)+1;
			$idz=int($z/$r)+1;
			
			$id_cube = $idx + $L*($idy-1) + $L*$L*($idz-1);
			
			push(@{$cube_array[$id_cube]},$id_residue);
			
			#printf "%10d %10d %10d %10d %10.2f %10.2f %10.2f %10d %10d %10d \n",$r,$L,$id,$id_atom,$x,$y,$z,$idx,$idy,$idz;
			#print "$id_cube,$id_residue,$name \n";
		}
		
	}
	return @cube_array;
}

#-------------------------------------------------------------------------
#                       sub pdb_quadrante1
#-------------------------------------------------------------------------

sub pdb_quadrante1{
	
	my $doc_aux=@_[0];
	
	my $x;
	my $y;
	my $z;
	
	my $xmin=0;
	my $ymin=0;
	my $zmin=0;
	
	my @atom_set;
	foreach my $atom (@{$doc_aux->Atoms}){
		push(@atom_set,$atom);
		
		my $x=$atom->X;
		my $y=$atom->Y;
		my $z=$atom->Z;
			
		if($x<$xmin){$xmin=$x}
		if($y<$ymin){$ymin=$y}
		if($z<$zmin){$zmin=$z}
		
	}

	my $dx=abs($xmin);
	my $dy=abs($ymin);
	my $dz=abs($zmin);

	$doc_aux->CreateSet("all_atoms", \@atom_set); 
	$doc_aux->Sets("all_atoms")->Translate(Point(X => $dx, Y => $dy, Z => $dz));
	$doc_aux->Sets->Delete;

	return $doc_aux;
}

#-------------------------------------------------------------------------
#                          sub water_array
#-------------------------------------------------------------------------

sub water_array{
	
	my $doc_aux=@_[0];      #pdb sem água
	my $doc2_aux=@_[1];     #pdb completo
	my $r_cube=@_[2];
	my $L=@_[3];
	my @cube_array=@{$_[4]};
	my $water_dist=@_[5];
	
	my $N=$L*$L;
	my @water_array;
	
	my $id_residue=-1;
	foreach my $residue (@{$doc_aux->SubUnits}){
		$id_residue=$id_residue+1;
			
		my $center = $residue->CenterOfGeometry;
		my $x=$center->X;
		my $y=$center->Y;
		my $z=$center->Z;
		
		my $id=id_cube($r_cube,$L,$x,$y,$z);
	
		my @v;
		$v[0]=$id-1;
		$v[1]=$id+1;
		$v[2]=$id+$L;
		$v[3]=$id-$L;
		$v[4]=$id-$L-1;
		$v[5]=$id-$L+1;
		$v[6]=$id+$L-1;
		$v[7]=$id+$L+1;
		
		
		$v[8]=$id+$N;
		$v[9]=$id+$N+1;
		$v[10]=$id+$N-1;
		$v[11]=$id+$N+$L;
		$v[12]=$id+$N-$L;
		$v[13]=$id+$N+$L+1;
		$v[14]=$id+$N+$L-1;
		$v[15]=$id+$N-$L+1;
		$v[16]=$id+$N-$L-1;
		
	
		$v[17]=$id-$N;
		$v[18]=$id-$N+1;
		$v[19]=$id-$N-1;
		$v[20]=$id-$N+$L;
		$v[21]=$id-$N-$L;
		$v[22]=$id-$N+$L+1;
		$v[23]=$id-$N+$L-1;
		$v[24]=$id-$N-$L+1;
		$v[25]=$id-$N-$L-1;
	
		$v[26]=$id;
	
		foreach my $cube (@v){	
			if($cube >=1 and $cube <= ($L*$L*$L)){
				my $nwater=$#{$cube_array[$cube]};
				for(my $i = 0; $i <= $nwater; $i++){
					my $id_water=$cube_array[$cube][$i];
					my($dist,$atom1,$atom2)=res_dist_files($doc_aux,$doc2_aux,$id_residue,$id_water);								
					if($dist<=$water_dist){push(@{$water_array[$id_residue]},$id_water)}					
				}				
			}
		}
	}

	return @water_array;
}	
	
#--------------------------------------------------------------------
#              			     sub uniq   
#--------------------------------------------------------------------
	
sub uniq {
    my %seen;
    grep !$seen{$_}++, @_;
}

#--------------------------------------------------------------------
#              			   sub water_list   
#--------------------------------------------------------------------

sub water_list{
	  
	my $id1=$_[0];
	my $id2=$_[1];
	my @caps1=@{$_[2]};
	my @caps2=@{$_[3]};
	my @water_array=@{$_[4]};

	my @list_water;
		
	if($#{$water_array[$id1]} >= 0){push(@list_water,@{$water_array[$id1]})}					
	if($#{$water_array[$id2]} >= 0){push(@list_water,@{$water_array[$id2]})}	
	
	foreach my $id_cap1 (@caps1){
		if($#{$water_array[$id_cap1]} >= 0){push(@list_water,@{$water_array[$id_cap1]})};
	}
	
	foreach my $id_cap2 (@caps2){
		if($#{$water_array[$id_cap2]} >= 0){push(@list_water,@{$water_array[$id_cap2]})};
	}
	
	@list_water=uniq(@list_water);
	
	return @list_water;
}


#--------------------------------------------------------------------
#              			   sub hbond_index   
#--------------------------------------------------------------------

sub hbond_index{
	  
	my $id1=$_[0];
	my $id2=$_[1];
	my @caps1=@{$_[2]};
	my @caps2=@{$_[3]};
	my @hbond_array=@{$_[4]};

	my @list_hbond;
		
	if($#{$hbond_array[$id1]} >= 0){push(@list_hbond,@{$hbond_array[$id1]})}					
	if($#{$hbond_array[$id2]} >= 0){push(@list_hbond,@{$hbond_array[$id2]})}	
	
	foreach my $id_cap1 (@caps1){
		if($#{$hbond_array[$id_cap1]} >= 0){push(@list_hbond,@{$hbond_array[$id_cap1]})};
	}
	
	foreach my $id_cap2 (@caps2){
		if($#{$hbond_array[$id_cap2]} >= 0){push(@list_hbond,@{$hbond_array[$id_cap2]})};
	}
	
	@list_hbond=uniq(@list_hbond);
	
	return @list_hbond;
}

#--------------------------------------------------------------------
#              			   sub copy_water   
#--------------------------------------------------------------------

sub copy_water{
	my $doc_aux=@_[0];
	my @water_list=@{$_[1]};

	foreach my $id_water (@water_list){	
		$doc1a->CopyFrom($doc_aux->SubUnits->Item($id_water));
		$doc2a->CopyFrom($doc_aux->SubUnits->Item($id_water));
		$doc3a->CopyFrom($doc_aux->SubUnits->Item($id_water));
		$doc4a->CopyFrom($doc_aux->SubUnits->Item($id_water));
	}
}

#--------------------------------------------------------------------
#              			   sub formal_charge   
#--------------------------------------------------------------------

sub formal_charge{
	my $doc_aux=@_[0];
	my $id_residue=@_[1];

	my $formalChargeSum=0;
	foreach my $atom (@{$doc_aux->SubUnits->Item($id_residue)->Atoms}){
		my $formalCharge = $atom->FormalCharge->Value;
		$formalChargeSum += $formalCharge;
	}
	
	return $formalChargeSum;
}



sub residue_class{
	my $res_name=@_[0];

	my $res_name=substr($res_name,0,3);
	my $res_class="Ligand";
	
	if ($res_name eq 'ALA'){$res_class='Np'}
	if ($res_name eq 'GLY'){$res_class='Np'}
	if ($res_name eq 'VAL'){$res_class='Np'}
	if ($res_name eq 'LEU'){$res_class='Np'}
	if ($res_name eq 'ILE'){$res_class='Np'}
	if ($res_name eq 'PRO'){$res_class='Np'}
	if ($res_name eq 'MET'){$res_class='Np'}
	if ($res_name eq 'PHE'){$res_class='Np'}
	if ($res_name eq 'TRP'){$res_class='Np'}

	
	if ($res_name eq 'SER'){$res_class='Un'}
	if ($res_name eq 'THR'){$res_class='Un'}
	if ($res_name eq 'CYS'){$res_class='Un'}
	if ($res_name eq 'ASN'){$res_class='Un'}
	if ($res_name eq 'GLN'){$res_class='Un'}
	if ($res_name eq 'TYR'){$res_class='Un'}

	
	if ($res_name eq 'ASP'){$res_class='Ac'}
	if ($res_name eq 'GLU'){$res_class='Ac'}

	
	if ($res_name eq 'LYS'){$res_class='Ba'}
	if ($res_name eq 'ARG'){$res_class='Ba'}
	if ($res_name eq 'HIS'){$res_class='Ba'}
	

	return $res_class;

}


#-----------------------------------------------------------------------------------------------------------------------------
#                                    FUNCOES ANTIGAS - SUBSTUTUIDAS POR NOVAS VERSOES
#-----------------------------------------------------------------------------------------------------------------------------



#-------------------------------------------------------------------------
#                           sub res2_caps  - Funcao antiga
#-------------------------------------------------------------------------

sub res2_caps{

	my $doc_aux=@_[0];
	my $id_residue=@_[1];
	my $lim_down=@_[2];
	my $lim_up=@_[3];
	my $caps_number=@_[4];
	
	my $i;
	my @caps;
	
	for ($i=1;$i<=$caps_number;$i++){
		
#            -------------------------
#                   Caps à direita
#            -------------------------
		
		my $id_cap_r = $id_residue + $i;
		
		my $ok="false";
		
		my $cap_right="lim_up";
		if($id_cap_r <= $lim_up){$cap_right=$doc_aux->SubUnits->Item($id_cap_r)->Name}
		
		my $C=$doc_aux->SubUnits->Item($id_residue)->Atoms->DisplayRange("C");
		my $nbonds=$C->AttachedAtoms->Count;
		
		for(my $j=0;$j<=$nbonds-1;$j++){
			my $res2=$C->AttachedAtoms->Item($j)->Ancestors->SubUnit->Name;
			if($res2 eq $cap_right){$ok="true"}
		
		}
		
		
		if($id_cap_r <= $lim_up and $ok eq "true"){push(@caps,$id_cap_r)}
	
#            -------------------------
#                   Caps à esquerda
#            -------------------------
			
		my $id_cap_l = $id_residue - $i;

		my $ok="false";
		
		my $cap_left="lim_down";
		if($id_cap_l >= $lim_down){$cap_left=$doc_aux->SubUnits->Item($id_cap_l)->Name}
		
		my $N=$doc_aux->SubUnits->Item($id_residue)->Atoms->DisplayRange("N");
		my $nbonds=$N->AttachedAtoms->Count;
	
		for(my $j=0;$j<=$nbonds-1;$j++){
			my $res2=$N->AttachedAtoms->Item($j)->Ancestors->SubUnit->Name;
			if($res2 eq $cap_left){$ok="true"}
		}

		if($id_cap_l >= $lim_down and $ok eq "true"){push(@caps,$id_cap_l)}
	
	}
	
	return @caps;
}

#-------------------------------------------------------------------------
#                           sub mfcc2_cysteine -> Funcao antiga
#-------------------------------------------------------------------------

sub mfcc2_cysteine{ 
	
	my $doc_aux=@_[0];
	my $id1=@_[1];
	my $id2=@_[2];
	my @caps1=@{$_[3]};
	my @caps2=@{$_[4]};

	my $res1=$doc_aux->SubUnits->Item($id1)->Name;
	my $res2=$doc_aux->SubUnits->Item($id2)->Name;
	
	my @id_list;
	
	push(@id_list,$id1);
	push(@id_list,$id2);
	push(@id_list,@caps1);
	push(@id_list,@caps2);
	
	#--------------------------------------------------------------------
	#                   Arquivo A  ->  k*-R1-k + C*-R2-C    
	#--------------------------------------------------------------------

	$doc1a->CopyFrom($doc_aux);
	
	my $id=-1;
	foreach my $residue (@{$doc1a->SubUnits}){
		$id=$id+1;
		my $ok="true";
		
		foreach my $id_aux (@id_list){	
			if($id == $id_aux){
				$ok="false";
				last;
			}
		}
		if($ok eq "true"){$doc1a->SubUnits->DisplayRange($residue->Name)->Delete}
	}

	foreach my $atom (@{$doc1a->Atoms}){
		if ($atom->Name eq "N" and $atom->NumBonds<3){$atom->AdjustHydrogen}
		if ($atom->Name eq "C" and $atom->NumBonds<3){$atom->AdjustHydrogen}
		if ($atom->Name eq "SG" and $atom->NumBonds<2){$atom->AdjustHydrogen}
	}
	
	#--------------------------------------------------------------------
	#                  Arquivo B  ->   k*-R1-k +  C*- C
	#--------------------------------------------------------------------
	
	$doc2a->CopyFrom($doc1a);
	$doc2a->SubUnits->DisplayRange($res2)->Delete;

	foreach my $atom (@{$doc2a->Atoms}){
		if ($atom->Name eq "N" and $atom->NumBonds<3){$atom->AdjustHydrogen}
		if ($atom->Name eq "C" and $atom->NumBonds<3){$atom->AdjustHydrogen}
		if ($atom->Name eq "SG" and $atom->NumBonds<2){$atom->AdjustHydrogen}
	}

	#--------------------------------------------------------------------
	#                  Arquivo C  ->  k*-k + C*-R2-C  
	#--------------------------------------------------------------------
	
	$doc3a->CopyFrom($doc1a);
	$doc3a->SubUnits->DisplayRange($res1)->Delete;
	
	foreach my $atom (@{$doc3a->Atoms}){
		if ($atom->Name eq "N" and $atom->NumBonds<3){$atom->AdjustHydrogen}
		if ($atom->Name eq "C" and $atom->NumBonds<3){$atom->AdjustHydrogen}
		if ($atom->Name eq "SG" and $atom->NumBonds<2){$atom->AdjustHydrogen}
	}

	#--------------------------------------------------------------------
	#                       Arquivo D  ->   k*-k + C*-C   
	#--------------------------------------------------------------------
	
	$doc4a->CopyFrom($doc1a);
	$doc4a->SubUnits->DisplayRange($res1)->Delete;
	$doc4a->SubUnits->DisplayRange($res2)->Delete;

	foreach my $atom (@{$doc4a->Atoms}){
		if ($atom->Name eq "N" and $atom->NumBonds<3){$atom->AdjustHydrogen}
		if ($atom->Name eq "C" and $atom->NumBonds<3){$atom->AdjustHydrogen}
		if ($atom->Name eq "SG" and $atom->NumBonds<2){$atom->AdjustHydrogen}
	}
	
}


#my $id=-1;
#foreach my $residue (@{$doc->SubUnits}){
#	$id=$id+1;
#	my $name=$residue->Name;	
#	if(substr($name,0,3) eq "CSS"){
#		my $id_residue=substr($name,3,10);
#		$residue->Name='CYS'.$id_residue;
#	}			
#}




