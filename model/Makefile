#
#    This file is part of WHOI Cable, a program for the static and dynamic
#    analysis of oceanographic cable structures.
#
#    Copyright (C) 1997-2016 by Woods Hole Oceanographic Institution (WHOI)
#    and Jason Gobat
#
#    WHOI Cable is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    WHOI Cable is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with WHOI Cable.  If not, see <http://www.gnu.org/licenses/>.
#

CC	 = gcc
ifeq ($(OS),windows)
   CC = gcc
   CFLAGS   = -O6 -Wall -funroll-loops -I../include -DZLIB -DWINGUI \
              -DWINDOWS -mms-bitfields
else ifeq ($(OS),OSX-universal)
    CFLAGS = -g -O2 -Wall -I../include -DZLIB \
             -force_cpusubtype_ALL -mmacosx-version-min=10.4 \
             -arch i386 -arch ppc 
else ifeq ($(OS),OSX)
    CFLAGS = -g -O2 -Wall -I../include -DZLIB -DGUI -mmacosx-version-min=10.4
else
    CFLAGS = -g -Wall -I../include -DZLIB -DGUI
endif

AR	 = ar cq
RANLIB	 = ranlib
LEX	 = flex
YACC = bison
LFLAGS	 = -L -i -PCable_ 
YFLAGS   = -o parser.c -d -v -p Cable_ 

OBJS	= problem.o code.o objects.o \
          parser.o lexer.o Tree.o file.o \
          nodes.o segments.o setup.o check.o buoyancy.o tension.o

LIB	= libmodel.a


all:            parser.c lexer.c $(LIB)

clean:;         @rm -f $(OBJS) $(LIB) core parser.output

$(LIB):		$(OBJS)
		rm -f $(LIB)
		$(AR) $(LIB) $(OBJS)
		$(RANLIB) $(LIB)

parser.c:   parser.y
			$(YACC) -o parser.c -d -v -p Cable_ parser.y

lexer.c:    parser.y lexer.l

Tree.o: Tree.c ../include/Tree.h ../include/Item.h 

code.o: code.c ../include/code.h \
 ../include/allocate.h ../include/malloc.h

file.o: file.c ../include/problem.h ../include/Tree.h \
 ../include/Item.h ../include/model.h \
 ../include/code.h ../include/error.h

lexer.o: lexer.c ../include/error.h \
 ../include/problem.h ../include/Tree.h ../include/Item.h \
 ../include/model.h ../include/code.h 

objects.o: objects.c ../include/problem.h ../include/Tree.h \
 ../include/Item.h ../include/model.h \
 ../include/code.h ../include/error.h ../include/objects.h \
 ../include/allocate.h ../include/malloc.h

options.o: options.c ../include/allocate.h ../include/malloc.h \
 ../include/options.h ../include/error.h

parser.o: parser.c ../include/code.h \
 ../include/error.h ../include/objects.h ../include/model.h \
 ../include/Tree.h ../include/Item.h ../include/problem.h \
 ../include/allocate.h ../include/malloc.h

problem.o: problem.c ../include/error.h \
 ../include/problem.h ../include/Tree.h ../include/Item.h \
 ../include/model.h ../include/code.h ../include/allocate.h \
 ../include/malloc.h
