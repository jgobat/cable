/*
    This file is part of WHOI Cable, a program for the static and dynamic
    analysis of oceanographic cable structures.

    Copyright (C) 1997-2016 by Woods Hole Oceanographic Institution (WHOI)
    and Jason Gobat

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

# ifndef _COMPRESS_H
# define _COMPRESS_H

# ifdef ZLIB	/* use the zlib compressed file read/write routines */

# include <zlib.h>

# define ResFile gzFile
# define res_open gzopen
# define res_close gzclose
# define res_read(ptr, size, nmemb, fp) gzread((fp), (ptr), ((size)*(nmemb))) 
# define res_write(ptr, size, nmemb, fp) gzwrite((fp), (ptr), ((size)*(nmemb)))
# define res_flush(fp) gzflush(fp, Z_SYNC_FLUSH) 
# define res_rewind(fp) gzrewind((fp))
# define res_seek(fp, offset, whence) gzseek((fp), (offset), (whence))
# define res_eof gzeof
# define res_tell gztell

# else		/* use the regular old libc stream based routines */

# define ResFile FILE *
# define res_open fopen
# define res_close fclose
# define res_read fread
# define res_write fwrite
# define res_flush fflush
# define res_eof feof
# define res_seek fseek
# define res_rewind frewind
# define res_tell ftell

# endif

# endif /* _COMPRESS_H */
