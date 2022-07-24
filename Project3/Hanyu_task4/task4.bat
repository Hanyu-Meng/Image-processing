@ECHO OFF
ECHO:
ECHO Task 4 Erosion Boundary
ECHO ---------------------------------------------------
cd bin
ECHO Performing the following commands...
ECHO:
ECHO:Struture Element centered at (1,1) r = 1_symmetry extension
ECHO Project3_task4.exe ../data/bike_100.bmp ../data/bike_erosion_sym.bmp 1 1 1 0
Project3_task4.exe ../data/bike_100.bmp ../data/bike_erosion_sym.bmp 1 1 1 0
ECHO:
ECHO: Struture Element centered at (1,1) r = 1_zero extension
ECHO  Project3 task4.exe ../data/bike_100.bmp ../data/bike_erosion_0.bmp 1 1 1 1
Project3_task4.exe ../data/bike_100.bmp ../data/bike_erosion_0.bmp 1 1 1 1

cd ../
cd data
mi_viewer bike_100.bmp
mi_viewer bike_erosion_sym.bmp.bmp
mi_viewer bike_erosion_0.bmp
cd ../
ECHO Task complete
