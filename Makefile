watermark: watermark_plugin.c
	LIBS="-lm -ljbig" CFLAGS=-g gimptool-2.0 --install watermark_plugin.c && gimp cropped.xcf

blake3_test: blake3_test.c
	gcc blake3_test.c -o blake3_test -I /home/isaac/Documents/C++/BLAKE3_install/include/ -L /home/isaac/Documents/C++/BLAKE3_install/lib/ -lblake3
