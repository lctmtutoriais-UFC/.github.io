#!/bin/csh -f
# =====================================================================
# LCTMBio - GROMACS - Dinamica Molecular de Proteina / Proteina::Proteina
# / Proteina::Peptideo (Em agua)
# Script de automacao dos comandos do protocolo (ver pagina do site).
# Ajuste as variaveis abaixo antes de rodar.
# =====================================================================

source /usr/local/gromacs/bin/GMXRC

# ---- variaveis editaveis ----------------------------------------------
set PDBIN      = "rank_001_prep_charmmGui.pdb"
set FF_OPT     = "1"      # opcao do campo de forca (pdb2gmx)
set WATER_OPT  = "1"      # opcao do modelo de agua (pdb2gmx)
set ION_GROUP  = "13"     # grupo de solvente para genion (SOL)
set POT_GROUP  = "11 0"   # grupo p/ gmx energy (potencial) - ajustar se necessario
set TEMP_GROUP = "16 0"   # grupo p/ gmx energy (temperatura)
set PRES_GROUP = "18 0"   # grupo p/ gmx energy (pressao)
set DENS_GROUP = "24 0"   # grupo p/ gmx energy (densidade)
set CHAIN_A    = "1-255"  # faixa de residuos da cadeia A (ajustar)
set CHAIN_B    = "256-735" # faixa de residuos da cadeia B (ajustar)
set NTOMP      = "12"
# ------------------------------------------------------------------------

echo "== I) pdb2gmx =="
echo "$FF_OPT\n$WATER_OPT" | gmx pdb2gmx -ignh -f $PDBIN -o processed.gro

echo "== II) editconf - caixa de simulacao =="
gmx editconf -f processed.gro -o newbox.gro -c -d 1.0 -bt cubic

echo "== III) solvate =="
gmx solvate -cp newbox.gro -cs spc216.gro -o solv.gro -p topol.top

echo "== IV) grompp (ions) =="
gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr

echo "== V) genion =="
echo "$ION_GROUP" | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname SOD -nname CLA -neutral -conc 0.15

echo "== VI) grompp (minimizacao) =="
gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr

echo "== VII) mdrun (minimizacao) =="
gmx mdrun -v -deffnm em
echo "$POT_GROUP" | gmx energy -f em.edr -o potential.xvg
echo "0" | gmx trjconv -s em.tpr -f em.trr -o em.pdb -dump 0

echo "== VIII) grompp (NVT) =="
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr

echo "== IX) mdrun (NVT) =="
gmx mdrun -deffnm nvt
echo "$TEMP_GROUP" | gmx energy -f nvt.edr -o temperature.xvg
echo "0" | gmx trjconv -s nvt.tpr -f nvt.xtc -o nvt.pdb -dump 0

echo "== X) grompp (NPT) =="
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr

echo "== XI) mdrun (NPT) =="
gmx mdrun -deffnm npt
echo "$PRES_GROUP" | gmx energy -f npt.edr -o pressure.xvg
echo "$DENS_GROUP" | gmx energy -f npt.edr -o density.xvg
echo "0" | gmx trjconv -s npt.tpr -f npt.xtc -o npt.pdb -dump 0

echo "== XII) grompp (producao MD) =="
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr

echo "== XIII) mdrun (producao MD) =="
nohup gmx mdrun -deffnm md -nb gpu -pin on -ntmpi 1 -ntomp $NTOMP

# =====================================================================
# Pos-processamento (rodar depois que a producao md terminar)
# =====================================================================

echo "== Correcao de PBC =="
echo "0" | gmx trjconv -s md.tpr -f md.xtc -o md_whole.xtc -pbc whole
echo "0" | gmx trjconv -s md.tpr -f md_whole.xtc -o md_whole_nojump.xtc -pbc nojump
echo "1\n0" | gmx trjconv -s md.tpr -f md_whole_nojump.xtc -o md_whole_nojump_center.xtc -center -pbc mol -ur compact
echo "4\n1" | gmx trjconv -s md.tpr -f md_whole_nojump_center.xtc -o md_whole_nojump_center_fit.xtc -fit rot+trans
echo "1\n1" | gmx trjconv -s md.tpr -f md_whole_nojump_center_fit.xtc -o md_whole_nojump_center_fit_protOnly.xtc -center -pbc mol -ur compact

echo "== Extracao de frames especificas =="
echo "0" | gmx trjconv -s md.tpr -f md_whole_nojump_center_fit.xtc -o 100ns.pdb -dump 100000
echo "1" | gmx trjconv -s md.tpr -f md_whole_nojump_center_fit_protOnly.xtc -o 0ns_protOnly.pdb -dump 0

echo "== RMSD =="
echo "4\n1" | gmx rms -s md.tpr -f md_whole_nojump_center_fit.xtc -o rmsd_proteinALL.xvg -tu ns
echo "4\n3" | gmx rms -s md.tpr -f md_whole_nojump_center_fit.xtc -o rmsd_proteinCALPHA.xvg -tu ns

echo "== RMSF =="
echo "1" | gmx rmsf -s md.tpr -f md_whole_nojump_center_fit.xtc -o rmsf_proteinALL_res.xvg -res
echo "3" | gmx rmsf -s md.tpr -f md_whole_nojump_center_fit.xtc -o rmsf_proteinCALPHA_res.xvg -res

echo "== Rg (raio de giro) =="
echo "1" | gmx gyrate -s md.tpr -f md_whole_nojump_center_fit.xtc -o gyrate_complex.xvg -tu ns

echo "== SASA =="
echo "1" | gmx sasa -s md.tpr -f md_whole_nojump_center_fit.xtc -o sasa_protein.xvg -n index.ndx

echo "== PCA =="
echo "4\n1" | gmx covar -s md.tpr -f md_whole_nojump_center_fit.xtc -o eigenval.xvg -v eigenvec.trr -n index.ndx
echo "4\n1" | gmx anaeig -s md.tpr -f md_whole_nojump_center_fit.xtc -v eigenvec.trr -eig eigenval.xvg -2d 2dproj.xvg -first 1 -last 2 -n index.ndx
echo "4\n1" | gmx anaeig -s md.tpr -f md_whole_nojump_center_fit.xtc -v eigenvec.trr -eig eigenval.xvg -extr extreme.pdb -first 1 -last 1 -nframes 30 -n index.ndx

echo "== Index para separar cadeias (edite CHAIN_A / CHAIN_B acima) =="
echo "AVISO: gmx make_ndx e interativo. Rode manualmente:"
echo "  gmx make_ndx -f solv_ions.gro -n index.ndx"
echo "  ri $CHAIN_A"
echo "  ri $CHAIN_B"
echo "  q"

echo "== HBONDS (edite os numeros de grupo apos criar o index) =="
echo "AVISO: ajuste os grupos de cadeia A/B criados no make_ndx antes de rodar:"
echo "  gmx hbond -f md_whole_nojump_center_fit.xtc -s md.tpr -num hbnum.xvg -hbm hb_matrix -n index.ndx"

echo "== Contacts (edite os numeros de grupo apos criar o index) =="
echo "AVISO: ajuste os grupos de cadeia A/B criados no make_ndx antes de rodar:"
echo "  gmx mindist -s md.tpr -f md_whole_nojump_center_fit.xtc -n index.ndx -od mindist_protein1_protein2.xvg -on contacts_protein1_protein2.xvg"

echo "== Pairdist =="
gmx pairdist -f md_whole_nojump_center_fit.xtc -s md.tpr -n index.ndx -o dist.xvg -ref 'group "chain A"' -sel 'group "chain B"'

echo "== Script finalizado =="
