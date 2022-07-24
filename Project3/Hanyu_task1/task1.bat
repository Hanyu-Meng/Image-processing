@ECHO OFF
ECHO:
ECHO Task 1 Black and White
ECHO ---------------------------------------------------
cd bin
ECHO Performing the following commands...
ECHO:
ECHO: pens_mono
ECHO Project3_task1 ../data/pens_mono.bmp ../data/pens_100.bmp 100
Project3_task1.exe ../data/pens_mono.bmp ../data/pens_100.bmp 100
ECHO:
ECHO: lenna_mono
ECHO  ../data/lenna_mono.bmp ../data/lenna_100.bmp 127
Project3_task1.exe ../data/lenna_mono.bmp ../data/lenna_100.bmp 127
ECHO:
ECHO: bike_mono
ECHO Project3_task1 ../data/bike_mono.bmp ../data/bike_100.bmp 100
Project3_task1.exe ../data/bike_mono.bmp ../data/bike_100.bmp 100
ECHO:
cd ../
cd data
mi_viewer lenna_100.bmp
mi_viewer pens_100.bmp
mi_viewer bike_100.bmp
cd ../
ECHO Task complete
