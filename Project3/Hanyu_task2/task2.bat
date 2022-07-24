@ECHO OFF
ECHO:
ECHO Task 1 Erosion
ECHO ---------------------------------------------------
cd bin
ECHO Performing the following commands...
ECHO:
ECHO: Struture Element centered at origin
ECHO Project3_task2.exe ../data/lenna_100.bmp ../data/lenna_erosion1.bmp 0 0 0 1 1 0 0 -1 -1 0
Project3_task2.exe ../data/lenna_100.bmp ../data/lenna_erosion1.bmp 0 0 0 1 1 0 0 -1 -1 0
ECHO:
ECHO: Struture Element with offset
ECHO  Project3 task2.exe ../data/lenna_100.bmp ../data/lenna_erosion2.bmp 0 1 1 0 1 1 1 2 2 1
Project3_task2.exe ../data/lenna_100.bmp ../data/lenna_erosion2.bmp 0 1 1 0 1 1 1 2 2 1
ECHO:
ECHO: Square
ECHO Project3 task2.exe ../data/lenna_100.bmp ../data/lenna_erosion3.bmp 0 0 0 1 0 -1 1 0 -1 0 1 1 1 -1 -1 1 -1 -1
Project3_task2.exe ../data/lenna_100.bmp ../data/lenna_erosion3.bmp 0 0 0 1 0 -1 1 0 -1 0 1 1 1 -1 -1 1 -1 -1
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
