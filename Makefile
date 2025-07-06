watermark: add_watermark verify_watermark
	 gimp cropped.xcf

add_watermark: add_watermark_plugin.c
	LIBS="-lm -ljbig -L/home/isaac/Documents/C++/BLAKE3_install/lib/ -lblake3" CFLAGS="-g -I/home/isaac/Documents/C++/BLAKE3_install/include/" gimptool-2.0 --install add_watermark_plugin.c

verify_watermark: verify_watermark_plugin.c
	LIBS="-lm -ljbig -L/home/isaac/Documents/C++/BLAKE3_install/lib/ -lblake3" CFLAGS="-g -I/home/isaac/Documents/C++/BLAKE3_install/include/" gimptool-2.0 --install verify_watermark_plugin.c

blake3_test: blake3_test.c
	gcc blake3_test.c -o blake3_test -I /home/isaac/Documents/C++/BLAKE3_install/include/ -L /home/isaac/Documents/C++/BLAKE3_install/lib/ -lblake3
