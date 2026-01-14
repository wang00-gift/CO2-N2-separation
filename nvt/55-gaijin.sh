#!/bash/bin
#source /home/opt/software/packmol/packmol-20.3.2/packmol
source ~/.bashrc

MAIN_DIR="20ns"
BEGIN_NUM="1"
END_NUM="100"

echo "正在创建主目录:$MAIN_DIR"
mkdir -p "$MAIN_DIR" 
cd "$MAIN_DIR"
for ((i="$BEGIN_NUM";i<="$END_NUM";i=i+1))
do
mkdir -p  $i/
cp -r ../start.inp    $i/
cp -r ../md.mdp      $i/
cp -r ../topol.top    $i/
cp -r ../CO2.pdb      $i/
cp -r ../N2.pdb       $i/
cp -r ../CGMD-3.py    $i/
cp -r ../forcefield   $i/
done
cd ../
echo "文件复制好了"

cp start.gro "$MAIN_DIR"/1/

cd "$MAIN_DIR"
for ((j="$BEGIN_NUM";j<="$END_NUM";j++))
do
cd $j/
python3 CGMD-3.py
cp 1.gro ../$((j+1))/start.gro
cd ..
done
cd ..

cd "$MAIN_DIR"
for ((z="$BEGIN_NUM";z<="$END_NUM";z++))
do
  cat $z/gas_rum.txt >> new-20il-gas.txt
done
echo "ok"
