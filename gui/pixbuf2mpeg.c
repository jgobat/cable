/*
    This file is part of WHOI Cable, a program for the static and dynamic
    analysis of oceanographic cable structures.

    Copyright (C) 2008-2016 by Jason Gobat

    WHOI Cable is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    WHOI Cable is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with WHOI Cable.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef WINDOWS

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gdk/gdk.h>

#include "libavcodec/avcodec.h"

#define SCALEBITS 8
#define ONE_HALF  (1 << (SCALEBITS - 1))
#define FIX(x)    ((int) ((x) * (1L<<SCALEBITS) + 0.5))

static void 
pixbuf_to_yuv420p(uint8_t *lum, uint8_t *cb, uint8_t *cr, GdkPixbuf *buf)
{
    int width, height;
    int wrap, wrap3, x, y;
    int r, g, b, r1, g1, b1;
    uint8_t *p;
    int stride;
    int nchan;

    width = gdk_pixbuf_get_width(buf);
    height = gdk_pixbuf_get_height(buf);
    stride = gdk_pixbuf_get_rowstride(buf);
    nchan  = gdk_pixbuf_get_n_channels(buf);

    wrap = width;
    wrap3 = stride; // width * 3;

    p = gdk_pixbuf_get_pixels(buf);

    for(y=0;y<height;y+=2) {
        for(x=0;x<width;x+=2) {
            r = p[0];
            g = p[1];
            b = p[2];
            r1 = r;
            g1 = g;
            b1 = b;
            lum[0] = (FIX(0.29900) * r + FIX(0.58700) * g +
                      FIX(0.11400) * b + ONE_HALF) >> SCALEBITS;
            r = p[nchan];
            g = p[nchan + 1];
            b = p[nchan + 2];
            r1 += r;
            g1 += g;
            b1 += b;
            lum[1] = (FIX(0.29900) * r + FIX(0.58700) * g +
                      FIX(0.11400) * b + ONE_HALF) >> SCALEBITS;
            p += wrap3;
            lum += wrap;

            r = p[0];
            g = p[1];
            b = p[2];
            r1 += r;
            g1 += g;
            b1 += b;
            lum[0] = (FIX(0.29900) * r + FIX(0.58700) * g +
                      FIX(0.11400) * b + ONE_HALF) >> SCALEBITS;
            r = p[nchan];
            g = p[nchan + 1];
            b = p[nchan + 2];
            r1 += r;
            g1 += g;
            b1 += b;
            lum[1] = (FIX(0.29900) * r + FIX(0.58700) * g +
                      FIX(0.11400) * b + ONE_HALF) >> SCALEBITS;

            cb[0] = ((- FIX(0.16874) * r1 - FIX(0.33126) * g1 +
                      FIX(0.50000) * b1 + 4 * ONE_HALF - 1) >> (SCALEBITS + 2)) + 128;
            cr[0] = ((FIX(0.50000) * r1 - FIX(0.41869) * g1 -
                     FIX(0.08131) * b1 + 4 * ONE_HALF - 1) >> (SCALEBITS + 2)) + 128;

            cb++;
            cr++;
            p += -wrap3 + 2 * nchan;
            lum += -wrap + 2;
        }
        p += wrap3;
        lum += wrap;
    }
}

typedef struct {
    AVCodec        *codec;
    AVCodecContext *ctxt;
    FILE           *fp;
    AVFrame        *picture;
    uint8_t        *outbuf;
    int             outbuf_size;
} ChartMovie;
    
int 
InitializeMovie(ChartMovie *m, char *filename, int width, int height)
{
    int size;

    // avcodec_init();         // must be called before using avcodec lib 
    avcodec_register_all(); // register all codecs

    /* find the mpeg1 video encoder */
    m -> codec = avcodec_find_encoder(CODEC_ID_MPEG1VIDEO);
    if (!m -> codec) {
        fprintf(stderr, "codec not found\n");
        return 1;
    }

    m -> ctxt = avcodec_alloc_context3(m->codec);
    m -> picture = avcodec_alloc_frame();

    /* put sample parameters */
    m -> ctxt -> strict_std_compliance = FF_COMPLIANCE_EXPERIMENTAL;
    m -> ctxt -> bit_rate = 400000;

    /* resolution must be a multiple of two */
    m -> ctxt -> width = width;
    m -> ctxt -> height = height;

    m -> ctxt -> time_base= (AVRational){1,24}; // frames per second
    m -> ctxt -> gop_size = 10; // emit one intra frame every ten frames 
    m -> ctxt -> max_b_frames=1;
    m -> ctxt -> pix_fmt = PIX_FMT_YUV420P;

    /* open it */
    if (avcodec_open2(m -> ctxt, m -> codec, NULL) < 0) {
        fprintf(stderr, "could not open codec\n");
        return 1;
    }

    m -> fp = fopen(filename, "wb");
    if (!m -> fp) {
        fprintf(stderr, "could not open %s\n", filename);
        return 1;
    }

    /* alloc image and output buffer */
    size = 3 * m->ctxt->width * m->ctxt->height;
    m->outbuf_size = size;
    m->outbuf = malloc(m->outbuf_size);

    m->picture->data[0] = malloc((size * 3) / 2);
    m->picture->data[1] = m->picture->data[0] + size;
    m->picture->data[2] = m->picture->data[1] + size / 4;
    m->picture->linesize[0] = m -> ctxt -> width;
    m->picture->linesize[1] = m -> ctxt -> width / 2;
    m->picture->linesize[2] = m -> ctxt -> width / 2;

    return 0;
}

int
AddFrameToMovie(ChartMovie *m, GdkPixbuf *pixbuf, int copies)
{
    int out_size;
    int i;

    pixbuf_to_yuv420p(m -> picture->data[0], 
                      m -> picture->data[1], 
                      m -> picture->data[2], pixbuf);

    for (i = 0 ; i < copies ; i++) {
        out_size = avcodec_encode_video(m -> ctxt, m -> outbuf, 
                                        m -> outbuf_size, m -> picture);
        fwrite(m -> outbuf, 1, out_size, m -> fp);
    }

    return 0;
}

void
CloseMovie(ChartMovie *m)
{
    int out_size;

    /* get the delayed frames */
    while((out_size = avcodec_encode_video(m -> ctxt, 
                                        m -> outbuf, m -> outbuf_size, NULL))) {
         fwrite(m -> outbuf, 1, out_size, m -> fp);
    }

    /* add sequence end code to have a real mpeg file */
    m->outbuf[0] = 0x00;
    m->outbuf[1] = 0x00;
    m->outbuf[2] = 0x01;
    m->outbuf[3] = 0xb7;

    fwrite(m->outbuf, 1, 4, m -> fp);
    fclose(m -> fp);
    free(m->picture->data[0]);
    free(m->outbuf);

    avcodec_close(m -> ctxt);
    av_free(m -> ctxt);
    av_free(m -> picture);
}

#endif
