#
#    This file is part of WHOI Cable, a program for the static and dynamic
#    analysis of oceanographic cable structures.
#
#    Copyright (C) 1997-2012 by Woods Hole Oceanographic Institution (WHOI)
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

CC       = gcc
CFLAGS   = -O2 -g -I../include 
AR       = ar cq
RANLIB   = ranlib

ifeq ($(OS),OSX-universal)
    CFLAGS = -Wall -O2 -g -I../include \
             -force_cpusubtype_ALL -mmacosx-version-min=10.4 \
             -arch i386 -arch ppc
else ifeq ($(OS),OSX)
    CFLAGS   = -Wall -O2 -g -I../include -mmacosx-version-min=10.4
else
    CFLAGS   = -Wall -O2 -g -I../include 
endif

LIB	= libwidgets.a
OBJS	= post.o tabgroup.o scroll.o util.o outputdialog.o
LOBJS   = Layout.o laygram.o laylex.o
# LOBJS	= 

all:            $(LIB)

clean:;         @rm -f $(OBJS) $(LOBJS) $(LIB)

$(LIB):         $(OBJS) $(LOBJS)
		rm -f $(LIB)
		$(AR) $(LIB) $(OBJS) $(LOBJS)
		$(RANLIB) $(LIB)
