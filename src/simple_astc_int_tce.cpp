#include <cstdint>

#include "platform.hpp"

#include <tceops.h>
#include <tce_vector.h>

/* Defines compile-time constants in an external repo using TCE */
#include "constants.hpp"

/** Read one block from the input image and calculate its min/max values
 */
inline void read_block_min_max(
    const uint8_t* __restrict__ inp,
    unsigned int block_id_x,
    unsigned int block_id_y,
    unsigned int img_w,
    uchar4* mincol,
    uchar4* maxcol,
    uchar4 block_pixels[BLOCK_PX_CNT]
) {
    unsigned int xb = block_id_x * BLOCK_X;
    unsigned int yb = block_id_y * BLOCK_Y;

    // Color endpoints - min and max colors in the current block
    uchar4 _mincol = { 255, 255, 255, 255 };
    uchar4 _maxcol = { 0, 0, 0, 0 };

    for (unsigned int y = 0; y < BLOCK_Y; ++y)
    {
        const unsigned int off = img_w * (yb + y) + xb;

        for (unsigned int x = 0; x < BLOCK_X; ++x)
        {
            const unsigned int i_inp = (off + x) << 2; // assuming NCH_RGB = 4
            const unsigned int i_out = y*BLOCK_X + x;

            uchar4 px;
            _load_4u8(inp+i_inp, &px);

            _min_4u8(px, _mincol, &_mincol);
            _max_4u8(px, _maxcol, &_maxcol);

            block_pixels[i_out] = px;
        }
    }

    *mincol = _mincol;
    *maxcol = _maxcol;
}

/** Downsample 12x12 block into 8x5 using bilinear filtering and quantize it
 *
 * The inp/out grid size is fixed and one channel per sample is assumed.
 * The output is also quantized.
 */
inline void downsample_12x12_to_8x5_u8_quant(
    const uint8_t inp[BLOCK_PX_CNT],
    uint8_t out[WGT_CNT]
){
    constexpr unsigned int w_inp = 12;
    constexpr unsigned int h_inp = 12;
    constexpr unsigned int w_out = 8;
    constexpr unsigned int h_out = 5;

    // size of the bilin. kernel
    constexpr unsigned int pixel_count_x = 3;
    constexpr unsigned int pixel_count_y = 6;

    // Buffer for holding intermediate results (columns stored as rows for easy
    // access in the second loop).
    // volatile uint8_t tmp[w_out*h_inp];
    volatile uint8_t tmp[h_inp*w_out];

    // First, interpolate rows.
    for (unsigned int y = 0; y < h_inp; ++y)
    {
        const unsigned int base_addr_i = y*w_inp;

        uchar4 inp0, inp4, inp8;
        _load_4u8(inp + base_addr_i + 0, &inp0);
        _load_4u8(inp + base_addr_i + 4, &inp4);
        _load_4u8(inp + base_addr_i + 8, &inp8);

        // m = 0
        const uchar4 wgt0 = BILIN_WEIGHTS_X_8_U8[0];

        uint32_t out0;
        _saturating_dot_acc_4u8(inp0, wgt0, 0, &out0);

        unsigned int addr_o = y;
        tmp[addr_o] = (uint8_t)(out0 >> 8);

        // m = 1
        const uchar4 wgt1 = BILIN_WEIGHTS_X_8_U8[1];

        uint32_t out1;
        _saturating_dot_acc_4u8(inp0, wgt1, 0, &out1);

        addr_o += h_inp;
        tmp[addr_o] = (uint8_t)(out1 >> 8);

        // m = 2
        const uchar4 wgt2a = BILIN_WEIGHTS_X_8_U8[2];
        const uchar4 wgt2b = BILIN_WEIGHTS_X_8_U8[3];

        uint32_t out2;
        _saturating_dot_acc_4u8(inp0, wgt2a, 0, &out2);
        _saturating_dot_acc_4u8(inp4, wgt2b, out2, &out2);

        addr_o += h_inp;
        tmp[addr_o] = (uint8_t)(out2 >> 8);

        // m = 3
        const uchar4 wgt3 = BILIN_WEIGHTS_X_8_U8[4];

        uint32_t out3;
        _saturating_dot_acc_4u8(inp4, wgt3, 0, &out3);

        addr_o += h_inp;
        tmp[addr_o] = (uint8_t)(out3 >> 8);

        // m = 4
        const uchar4 wgt4 = BILIN_WEIGHTS_X_8_U8[5];

        uint32_t out4;
        _saturating_dot_acc_4u8(inp4, wgt4, 0, &out4);

        addr_o += h_inp;
        tmp[addr_o] = (uint8_t)(out4 >> 8);

        // m = 5
        const uchar4 wgt5a = BILIN_WEIGHTS_X_8_U8[6];
        const uchar4 wgt5b = BILIN_WEIGHTS_X_8_U8[7];

        uint32_t out5;
        _saturating_dot_acc_4u8(inp4, wgt5a, 0, &out5);
        _saturating_dot_acc_4u8(inp8, wgt5b, out5, &out5);

        addr_o += h_inp;
        tmp[addr_o] = (uint8_t)(out5 >> 8);

        // m = 6
        const uchar4 wgt6 = BILIN_WEIGHTS_X_8_U8[8];

        uint32_t out6;
        _saturating_dot_acc_4u8(inp8, wgt6, 0, &out6);

        addr_o += h_inp;
        tmp[addr_o] = (uint8_t)(out6 >> 8);

        // m = 7
        const uchar4 wgt7 = BILIN_WEIGHTS_X_8_U8[9];

        uint32_t out7;
        _saturating_dot_acc_4u8(inp8, wgt7, 0, &out7);

        addr_o += h_inp;
        tmp[addr_o] = (uint8_t)(out7 >> 8);
    }

    // Next, columns
    constexpr unsigned int SHR_QUANT = 8 + 6; // Quantize to 2b while storing
    for (unsigned int m = 0; m < w_out; ++m)
    {
        const unsigned int base_addr_i = m*h_inp;

        uchar4 tmp0, tmp4, tmp8;
        _load_4u8(tmp + base_addr_i + 0, &tmp0);
        _load_4u8(tmp + base_addr_i + 4, &tmp4);
        _load_4u8(tmp + base_addr_i + 8, &tmp8);

        // n = 0
        const uchar4 wgt0 = BILIN_WEIGHTS_Y_5_U8[0];

        uint32_t out0;
        _saturating_dot_acc_4u8(tmp0, wgt0, 0, &out0);

        unsigned int addr_o = m;
        out[addr_o] = (uint8_t)(out0 >> SHR_QUANT);

        // n = 1
        const uchar4 wgt1a = BILIN_WEIGHTS_Y_5_U8[1];
        const uchar4 wgt1b = BILIN_WEIGHTS_Y_5_U8[2];

        uint32_t out1;
        _saturating_dot_acc_4u8(tmp0, wgt1a, 0, &out1);
        _saturating_dot_acc_4u8(tmp4, wgt1b, out1, &out1);

        addr_o += w_out;
        out[addr_o] = (uint8_t)(out1 >> SHR_QUANT);

        // n = 2
        const uchar4 wgt2a = BILIN_WEIGHTS_Y_5_U8[3];
        const uchar4 wgt2b = BILIN_WEIGHTS_Y_5_U8[4];
        const uchar4 wgt2c = BILIN_WEIGHTS_Y_5_U8[5];

        uint32_t out2;
        _saturating_dot_acc_4u8(tmp0, wgt2a, 0, &out2);
        _saturating_dot_acc_4u8(tmp4, wgt2b, out2, &out2);
        _saturating_dot_acc_4u8(tmp8, wgt2c, out2, &out2);

        addr_o += w_out;
        out[addr_o] = (uint8_t)(out2 >> SHR_QUANT);

        // n = 3
        const uchar4 wgt3a = BILIN_WEIGHTS_Y_5_U8[6];
        const uchar4 wgt3b = BILIN_WEIGHTS_Y_5_U8[7];

        uint32_t out3;
        _saturating_dot_acc_4u8(tmp4, wgt3a, 0, &out3);
        _saturating_dot_acc_4u8(tmp8, wgt3b, out3, &out3);

        addr_o += w_out;
        out[addr_o] = (uint8_t)(out3 >> SHR_QUANT);

        // n = 4
        const uchar4 wgt4 = BILIN_WEIGHTS_Y_5_U8[8];

        uint32_t out4;
        _saturating_dot_acc_4u8(tmp8, wgt4, 0, &out4);

        addr_o += w_out;
        out[addr_o] = (uint8_t)(out4 >> SHR_QUANT);
    }

    // for (int y = 0; y < h_out; ++y)
    // {
    //     for (int x = 0; x < w_out; ++x)
    //     {
    //         printf(" %#04x", out[y*w_out+x]);
    //     }
    //     printf("\n");
    // }
}

/** Encode a block of pixels
 *
 * The block_id_x/y could be substituted with get_global_id() if converted to
 * OpenCL.
 */
void encode_block_int(
    volatile const uint8_t* __restrict__ inp_img,
    unsigned int block_id_x,
    unsigned int block_id_y,
    unsigned int img_w,
    volatile uint8_t* __restrict__ out_data
){
    // Input pixel block as 1D array & min/max pixel values
    uchar4 block_pixels[BLOCK_PX_CNT / 4];
    uchar4 mincol, maxcol;

    // Read pixel block into a 1D array and select min/max at the same time
    read_block_min_max(
        inp_img,
        block_id_x,
        block_id_y,
        img_w,
        &mincol,
        &maxcol,
        block_pixels
    );

    // printf("block\n");
    // print_minmax_u8("   ", mincol, maxcol);
    // printf("mincol    : %#04x %#04x %#04x\n", mincol.x, mincol.y, mincol.z);

    // Move the endpoints line segment such that mincol is at zero
    uchar4 ep_vec;
    _saturating_sub_4u8(maxcol, mincol, &ep_vec);

    // Inset min/max (shrink the bounding box to reduce the effect of outliers)
    const uchar4 inset = ep_vec >> 5;
    mincol = mincol + inset;
    maxcol = maxcol - inset;

    // printf("mincol    : %#04x %#04x %#04x\n", mincol.x, mincol.y, mincol.z);

    // print_minmax_u8("ins", mincol, maxcol);

    // Quantize the endpoints
    const uchar4 mincol_quant = quantize_5b(&mincol);
    const uchar4 maxcol_quant = quantize_5b(&maxcol);

    // printf("mincol    : %#04x %#04x %#04x\n", mincol.x, mincol.y, mincol.z);

    // print_minmax_u8("mq ", mincol, maxcol);
    // print_minmax_u8("mqq", mincol_quant, maxcol_quant);
    // print_minmax_u8("ep ", ep_vec, ep_vec);
    // printf("e32: %#010x", *as_u32(&ep_vec));

    // Pixels quantized as 2b weights
    uint8_t quantized_weights[WGT_CNT] = { 0 };

    if (*as_u32(&mincol_quant) != *as_u32(&maxcol_quant))
    {
        // Projection of pixels onto ep_vec to get the ideal weights
        // First, we normalize the endpoint vector and scale it.
        uint32_t ep_dot;     // Q2.16
        _saturating_dot_acc_4u8(ep_vec, ep_vec, 0, &ep_dot);

        // Q10.22, max. value 1024 for 5b quant
        const uint32_t inv_ep_dot = approx_inv_u32(ep_dot);

        // printf("ep_dot: %#010x\n", ep_dot);
        // printf("inv   : %#010x\n", inv_ep_dot);

        // The scaled ep_vec can be max. 32 due to 5b quantization
        uint32_t ep_vec32[3];
        _unpack_rgb_u32(ep_vec, ep_vec32);

        uint32_t ep_sc32[3] = {
            ep_vec32[0] * (inv_ep_dot >> 8),  // Q0.8 * Q10.14 = Q10.22
            ep_vec32[1] * (inv_ep_dot >> 8),
            ep_vec32[2] * (inv_ep_dot >> 8),
        };

        // Scaling the endpoint vector to fit 8 bits. The Q representation
        // depends on the magnitude of the division result above: between Q5.3
        // and Q0.8.
        const uint32_t log2ep[3] = {
            log2(ep_sc32[0]),
            log2(ep_sc32[1]),
            log2(ep_sc32[2]),
        };

        uint32_t sc;
        _max_u32(log2ep[0], log2ep[1], &sc);
        _max_u32(log2ep[2], sc, &sc);
        const uint32_t q_int = sc >= 21 ? sc - 21 : 0; // no. integer bits
        const uint32_t shr = 14 + q_int;               // how much to rshift

        ep_sc32[0] >>= shr;
        ep_sc32[1] >>= shr;
        ep_sc32[2] >>= shr;

        uchar4 ep_sc8;
        pack_rgb_u8(ep_sc32, &ep_sc8);

        // One (almost), represented in the variable Q format.
        const uint32_t shr_res = 8 - q_int;
        const uint32_t one = (1 << shr_res) - 1;

        // printf("ep_sc32[0]: %#010x\n", ep_sc32[0]);
        // printf("ep_sc32[1]: %#010x\n", ep_sc32[1]);
        // printf("ep_sc32[2]: %#010x\n", ep_sc32[2]);
        // printf("ep_sc8    : %#010x\n", *as_u32(&ep_sc8));
        // printf("ep_sc8    : %#04x %#04x %#04x\n", ep_sc8.x, ep_sc8.y, ep_sc8.z);
        // printf("one       : %#010x\n", one);
        // printf("mincol    : %#04x %#04x %#04x\n", mincol.x, mincol.y, mincol.z);
        // print_minmax_u8("          :", mincol, maxcol);

        // Pixels mapped to the endpoint line without quantization and
        // downsampling
        uint8_t ideal_weights[BLOCK_PX_CNT];

        for (unsigned int i = 0; i < BLOCK_PX_CNT; ++i)
        {
            uchar4 diff;
            _saturating_sub_4u8(block_pixels[i], mincol, &diff);

            // dot product max: Q5.3 * Q0.8 + 3xADD = Q8.11 -> sat 5.11 (0.11)
            // dot product min: Q0.8 * Q0.8 + 3xADD = Q3.16 -> sat 0.16
            // recover lost precision by adding one
            uint32_t res;
            // TODO: this does overflow 16 bits:
            _saturating_dot_acc_4u8(diff, ep_sc8, one, &res);
            ideal_weights[i] = (uint8_t)(u32min(res >> shr_res, 255));
            // ideal_weights[i] = (uint8_t)(res >> shr_res);  // -> Q0.8

            // printf("px  [%3d] : %#04x %#04x %#04x\n", i, block_pixels[i].x, block_pixels[i].y, block_pixels[i].z);
            // printf("diff[%3d] : %#04x %#04x %#04x\n", i, diff.x, diff.y, diff.z);
            // printf("dot [%3d] : %6d\n", i, res);
            // printf("iwgt[%3d] : %3d\n", i, ideal_weights[i]);
        }

        downsample_12x12_to_8x5_u8_quant(
            ideal_weights,
            quantized_weights
        );
    }

    // Output buffer for quantized weights and output data
	uint8_t wgt_buf[BLOCK_SIZE] = { 0 };

    // weights ISE encoding
    constexpr uint8_t wgt_bits = 2;
    constexpr uint8_t wgt_byte_cnt = 10; // 40 * 2 / 8
    unsigned int j = 0;
    for (unsigned int i = 0; i < WGT_CNT; i += 4)
    {
        uint32_t wgt;
        _load_u32(quantized_weights + i , &wgt);

        uint8_t res = wgt & 0x03;
        res |= (wgt >> 6) & 0x0c;
        res |= (wgt >> 12) & 0x30;
        res |= (wgt >> 18) & 0xc0;

        //debug
        //_TCE_ST8(TEST+j, res);

        wgt_buf[j++] = res;
    }

    // Calculate output data address (16 bytes per block)
    out_data = out_data + ((block_id_y*NB_X + block_id_x) << BLOCK_SIZE_EXP);

    // write out weights
    for (unsigned int i = 0; i < BLOCK_SIZE; i += 4)
    {
        uint32_t wgt4 = *(uint32_t*)(&wgt_buf[12 - i]);
        uint32_t reflected;
        _reflect_u32(wgt4, &reflected);
        _store_u32(&out_data[i], reflected);
    }

    // write out mode, partition, CEM
    constexpr uint16_t block_mode = 102;
    _write_bits(block_mode, 11, 0, out_data);
    constexpr unsigned int partition_count = 1;
    _write_bits(partition_count - 1, 2, 11, out_data);
    constexpr unsigned int color_format = 8;
    _write_bits(color_format, 4, 13, out_data);

    // quantized endpoint output data (layout is R0 R1 G0 G1 B0 B1)
    uint8_t mincol_quant_unpacked[3];
    uint8_t maxcol_quant_unpacked[3];

    _unpack_rgb_u8(mincol_quant, (uint8_t*)(mincol_quant_unpacked));
    _unpack_rgb_u8(maxcol_quant, (uint8_t*)(maxcol_quant_unpacked));

    uint8_t endpoints_q[6] = {
        (uint8_t)(mincol_quant_unpacked[0]),
        (uint8_t)(maxcol_quant_unpacked[0]),
        (uint8_t)(mincol_quant_unpacked[1]),
        (uint8_t)(maxcol_quant_unpacked[1]),
        (uint8_t)(mincol_quant_unpacked[2]),
        (uint8_t)(maxcol_quant_unpacked[2]),
    };

    // write out endpoints
    unsigned int off = 17;  // starting bit position for endpoints data
    constexpr uint8_t ep_bits = 5;
    off = 17;  // starting bit position for endpoints data
    for (unsigned int i = 0; i < 6; ++i)
    {
        _write_bits(endpoints_q[i], ep_bits, off, out_data);
        off += ep_bits;
    }
}
