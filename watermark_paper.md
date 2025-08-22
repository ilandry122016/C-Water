# Watermarking #

## Introduction ##

Hyper-realistic but fake images have been a growing concern in recent years. Digital images are easy to forge, compared with analog images. While there are methods to watermark, it has been difficult to make them reversible. This project aims to make reversible watermarking. Unlike most methods of cryptography, this project aims to use fragile watermarking (cite). This is because watermarking will be used for authentication unlike identifying fake images.

## Method ##

The watermark we add is a cryprographic hash of the existing image. 

`add_watermark` is a function that adds a watermark. It does do by compressing some bits of the image and squeezing the hash into the space. `verify_watermark` is a function that removes a waermark and restores the original image.

If the image to be watermarked is too small, `add_watermark` will not watermark it, but will return an error.

   `add_watermark` will first get the bounds of the image, as well as the pixels in the image. Then it divides the image in blocks of 8x8 pixels. Assuming 3 channels (colors) per pixel, and 8 x 8 pixels per block, it aims to get the number of blocks in the image. computes  There are 2 bits per block stored in `watermark_bits`, a `guchar` pointer. It then iterates over each block using the `x_block` as the column index of each block, and the `y_block` as the row index of each block. The `watermark_bit_index` is computed to get the byte for the block.
   
     1 - α      -α       α  -1 + α  -1 + α       α      -α  1 - α
        -α   1 + α  -1 - α       α       α  -1 - α   1 + α     -α
         α  -1 + α   1 + α      -α      -α   1 + α  -1 - α      α
    -1 + α       α      -α   1 - α   1 - α      -α       α -1 + α
    -1 + α       α      -α   1 - α   1 - α      -α       α -1 + α
         α  -1 + α   1 + α      -α      -α   1 + α  -1 - α      α
        -α   1 + α  -1 - α       α       α  -1 - α   1 + α     -α
     1 - α      -α       α  -1 + α  -1 + α       α      -α  1 - α

The two bits to be computed are `G_6_6_central`, and `G_6_6_edge`. `G_6_6` is the DCT (discrete cosine transform) coefficient at coordinates (6, 6). `G_6_6_central` and `G_6_6_edge` are the components of `G_6_6`. G_6_6 was selected for simplicity. They were selected as part of a spread-spectrum, frequency-based robust watermarking algorithm with the following formula:

	D'_k = d_k + α S_k, k = 1, ... , N_m,
	
where `D'_k` denotes the modulated DCT coefficients, `α` is an image independent watermark strength, and `N_m` is the number of modified coefficients.

The 1 + α term is for `G_6_6_central` and the α term is for `G_6_6_edge`. They are then saved into the `original_bits` array. `add_watermark` then gets a byte to set with `original_bits[original_bit_index]`. The `orig_value` is shifted using  "`(orig_value << (sub_block_index * 2))`" so that the two bits to set are in the right place in the byte. They are then "|"  together to update the bit pattern for that byte.

The hash is finalized with the default output length, `BLAKE3_OUT_LEN`, set to 32 bytes.

The `original_bits` were copied into `copy_bits` because `jbg_enc_init` uses the pointer of `bitmaps`, a `char` pointer array. This would cause `original_bits` to be changed unexpectedly if it was used.

Then the image gets encoded and allocated resources used by JBIG are freed. The hash of all pixels in the image are stored in `new_bits`. The bits used for watermarking in a compressed form are then stored in `new_bits` starting from the index `BLAKE3_OUT_LEN`. But the rest of the `original_bits` gets copied into `new_bits` starting from index `BLAKE3_OUT_LEN + compressed_size`.

The blocks are then iterated over again from `i = y1` to `i < y1 + max_row_8` with `i` incremented by 8 each iteration. This time, it gets the byte with `new_bits[original_bit_index]`. The byte is then shifted  so that the two bits we want are all the way to the right with "`>> (sub_block_index * 2))`" (we multiply `sub_block_index` by 2 because there are 2 bits per block). "& 3" to select the lowest two bits because other bits are for other blocks. The result is then stored in `new_value`. A similar thing occurs for `original_value`, except that `original_bits` is in place of `new_bits`.

If `original_bit_1_p_alpha` is not equal to `bit_1_p_alpha`, this means `add_watermark` is adding a bit. It then iterates on `x = 1` and `x = 2`, and `y = 1` and `y = 2`. In each iteration it sets `row_arr[y][col_offset + x] = add_16(row_arr[y][col_offset + x], sgn)`, where `sgn = ((((x + y) % 2) == 0) ? 1 : -1) * add_subtract`.

If `original_bit_alpha` is not equal to `bit_alpha`, this means that `bit_alpha == 1`. It then iterates on `x = 1` and `x = 2`, and `y = 1` and `y = 2`. In each iteration it sets `row_arr[y][col_offset + x] = add_16(row_arr[y][col_offset + x], sgn)`, where `sgn = ((((x + y) % 2) == 0) ? 1 : -1) * add_subtract`.

`add_watermark` then finally flushes the drawable (gets the data transfered to the core), merges the `shadow_buffer` with drawable, and then updates the entirety of drawable.

`verify_watermark` will first get the bounds of the image, as well as the pixels in the image. Then it divides the image in blocks of 8x8 pixels. Assuming 3 channels (colors) per pixel, and 8 x 8 pixels per block, it aims to get the number of blocks in the image. computes  There are 2 bits per block stored in `watermark_bits`, a `guchar` pointer. It then iterates over each block using the `x_block` as the column index of each block, and the `y_block` as the row index of each block. The `watermark_bit_index` is computed to get the byte for the block.

The two bits to be computed are `G_6_6_central`, and `G_6_6_edge`. They are then saved into the `watermark_bits` array. `add_watermark` then gets a byte to set with `watermark_bits[watermark_bit_index]`. The `orig_value` is shifted using  "`(orig_value << (sub_block_index * 2))`" so that the two bits to set are in the right place in the byte. They are then "|"  together to update the bit pattern for that byte.

Then the watermark gets decompressed. The array of `recovered_bits`, which is pointed to by a `char` pointer, as well as the values of `result_size`, `result_width`, and `result_height`, are obtained.

The `BLAKE3_HASHER` is then initialized.

The blocks are then iterated over again from `i = y1` to `i < y1 + max_row_8` with `i` incremented by 8 each iteration. This time, it gets the byte with `recovered_bits[recovered_bit_index]`. The byte is then shifted so that the two bits we want are all the way to the right with "`>> (sub_block_index * 2))`" (we multiply `sub_block_index` by 2 because there are 2 bits per block). "`& 3`" to select the lowest two bits because other bits are for other blocks. The result is then stored in `new_value`. A similar thing occurs for `watermark_value`, except that `watermark_bits` is in place of `recovered_bits`.

If `watermark_bit_1_p_alpha` is not equal to `bit_1_p_alpha`, this means `add_watermark` is adding a bit. It then iterates on `x = 1` and `x = 2`, and `y = 1` and `y = 2`. In each iteration it sets `row_arr[y][col_offset + x] = add_16(row_arr[y][col_offset + x], sgn)`, where `sgn = ((((x + y) % 2) == 0) ? 1 : -1) * add_subtract`.

If `watermark_bit_alpha` is not equal to `bit_alpha`, this means that `bit_alpha == 1`. It then iterates on `x = 1` and `x = 2`, and `y = 1` and `y = 2`. In each iteration it sets `row_arr[y][col_offset + x] = add_16(row_arr[y][col_offset + x], sgn)`, where `sgn = ((((x + y) % 2) == 0) ? 1 : -1) * add_subtract`.

`verify_watermark` then finally flushes the `drawable` (get the data transfered to the core), merges the `shadow_buffer` with `drawable`, and then updates the entirety of `drawable`.

If the image to verify does not have a valid watermark, `verify_watermark` will not dewatermark the image, but will return an error. If the image's hash is too big, then `verify_watermark` will also not dewatermark the image. Instead, it will return an error message.

## Results ##

Here we have the two images, and it is very diffucult to see that it is there, but it is there and reversible.

Before
![Image](cropped.png "icon")

After
![Image](cropped_with_watermark.png "icon")

After zooming in to the top left corner of the same images, we get the following

Before
![Image]( "icon")

After
![Image]( "icon")
