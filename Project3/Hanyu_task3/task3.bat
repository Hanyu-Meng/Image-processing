@ECHO OFF
ECHO:
ECHO Task 3 Erosion (circle)
ECHO ---------------------------------------------------
cd bin
ECHO Performing the following commands...
ECHO:
ECHO:Struture Element centered at origin r = 1
ECHO Project3_task3.exe ../data/lenna_100.bmp ../data/lenna_erosion1.bmp 1 0 0
Project3_task3.exe ../data/lenna_100.bmp ../data/lenna_erosion1.bmp 1 0 0
ECHO:
ECHO: Struture Element centered at (1,1) r = 2
ECHO  Project3 task3.exe ../data/lenna_100.bmp ../data/lenna_erosion2.bmp 2 1 1
Project3_task3.exe ../data/lenna_100.bmp ../data/lenna_erosion2.bmp 2 1 1
ECHO:
ECHO:  Struture Element centered at (10,10) r = 1
ECHO Project3 task3.exe ../data/lenna_100.bmp ../data/lenna_erosion3.bmp 1 10 10
Project3_task3.exe ../data/lenna_100.bmp ../data/lenna_erosion3.bmp 1 10 10
ECHO:
ECHO Playing all the BMP files...
mi_pipe2 :: read_file -f ../data/lenna_100.bmp ../data/lenna_erosion1.bmp ../data/lenna_erosion2.bmp ../data/lenna_erosion3.bmp :: frame_repeat :: view -play -rate 1
ECHO:
cd ../
cd data
mi_viewer lenna_100.bmp
mi_viewer lenna_erosion1.bmp
mi_viewer lenna_erosion2.bmp
mi_viewer lenna_erosion3.bmp
cd ../
ECHO Task complete
