watermark: add_watermark verify_watermark
	 gimp cropped.xcf

watermark_whole_image: add_watermark verify_watermark
	gimp watermarked_lincoln_statue_whole_image.xcf

watermark_original_image: add_watermark verify_watermark
	gimp PXL_20221125_205752609.jpg

watermark_32: add_watermark verify_watermark
	gimp cropped_32.xcf

watermark_3024x4031: add_watermark verify_watermark
	gimp lincoln_statue_3024x4031.png

watermark_mandrill: add_watermark verify_watermark
	gimp mandrill.tiff

add_watermark: add_watermark_plugin.c
	LIBS="-lm -ljbig -L/home/isaac/Documents/C++/BLAKE3_install/lib/ -lblake3" CFLAGS="-g -I/home/isaac/Documents/C++/BLAKE3_install/include/" gimptool-2.0 --install add_watermark_plugin.c

verify_watermark: verify_watermark_plugin.c
	LIBS="-lm -ljbig -L/home/isaac/Documents/C++/BLAKE3_install/lib/ -lblake3" CFLAGS="-g -I/home/isaac/Documents/C++/BLAKE3_install/include/" gimptool-2.0 --install verify_watermark_plugin.c

blake3_test: blake3_test.c
	gcc blake3_test.c -o blake3_test -I /home/isaac/Documents/C++/BLAKE3_install/include/ -L /home/isaac/Documents/C++/BLAKE3_install/lib/ -lblake3

biggest_range: biggest_range.cpp
	g++ -O3 biggest_range.cpp -o biggest_range -Wall -Wextra
